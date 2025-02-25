c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine gencsf3(nec,spin,mdon,mdoo,norb,ipt1,ismult,ncsf2,     12d12s19
     $     lprintx,ixw1,ixw2,ncsfx,iusewk,ibcb4,icsfxv,bc,ibc)          11d11s22
      implicit real*8 (a-h,o-z)
      external second                                                   7d3s19
      logical lprint,lprintx                                            12d12s19
      integer*1 icode(64)
      integer*8 i18,i28,i38,i48,i58,ipack8,icsfxv                       2d1s21
      dimension iost(3),ncomp(3),ndet(3),ncsf(3),ipt1(*),idata(6),      8d2s19
     $     ncsf2(4,*),ncsfx(*),ipack4(2)                                12d4s20
      equivalence (ipack8,ipack4)                                       12d4s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      include "common.store"
      include "common.print"                                            1d3s20
      if(lprintx)write(6,*)('hello, my name is gencsf3'),norb,ibcoff           12d12s19
      icsfxv=-1                                                         7d19s24
c
c     call cupo2 to set up cache storage
c
c
      is02=nint(spin*2d0)
      if(iprtr(2).ne.0)then                                             1d3s20
       lprint=.true.
      else                                                              1d3s20
       lprint=.false.                                                    12d12s19
      end if                                                            1d3s20
      ibctop0=ibcoff                                                    12d6s20
c
c     we have closed shell orbitals - call them c
c     we have open shell orbitals - call them o'
c     we have unoccupied orbitals - call them u".
c
c     in det space, singles are
c     co', converting one closed shell in to open, and one open to closed
c     cu", decreasing the number of closed shells by 1 and increasing the
c         number of open shells by two
c     o'u", changing one of the open shells, and
c     o'o', increasing the number of closed shells by 1 and decreasing the
c         number of open shells by two.
c
c     in det space, doubles are
      iptv=ibcoff                                                       6d18s19
      ibcoff=iptv+mdoo+1-mdon                                           6d18s19
      jptv=iptv-1                                                       6d18s19
      ismultm=ismult-1                                                  12d1s20
      ismultp=ismultm+2                                                 12d1s20
      ismultl=ismultm-2                                                 12d1s20
      do id=mdon,mdoo
       nopen=nec-id*2                                                   6d5s19
       idod=id+1-mdon                                                   6d11s19
       ipt1(idod)=ibcoff                                                6d11s19
       ibcoff=ibcoff+13                                                 6d25s19
       nopenmm=nopen-2                                                  12d1s20
       if(nopen.ge.2)then                                               12d1s20
        call fetchcsf2(nopenmm,ismultm,nsames,bc(ibcoff),ibc(ibcoff),   12d1s20
     $      ibc(ibcoff),ismultm,ndeto,bc,ibc)                           11d15s22
       else                                                             12d1s20
        nsames=0                                                        12d1s20
       end if                                                           12d1s20
       if(nopenmm.ge.ismultp)then                                       12d1s20
        call fetchcsf2(nopenmm,ismultp,nhighs,bc(ibcoff),ibc(ibcoff),    12d1s20
     $      ibc(ibcoff),ismultm,ndeto,bc,ibc)                           11d15s22
       else                                                             12d1s20
        nhighs=0                                                        12d1s20
       end if                                                           12d1s20
       if(ismultl.ge.0.and.nopenmm.ge.ismultl)then                      12d1s20
        call fetchcsf2(nopenmm,ismultl,nlows,bc(ibcoff),ibc(ibcoff),    12d1s20
     $      ibc(ibcoff),ismultl,ndeto,bc,ibc)                           11d15s22
       else                                                             12d1s20
        nlows=0                                                         12d1s20
       end if                                                           12d1s20
       call fetchcsf2(nopen,ismultm,nrooto,bc(ibcoff),ibc(ibcoff),      12d1s20
     $      ibc(ibcoff),ismultm,ndeto,bc,ibc)                           11d15s22
       ndet(1)=ndeto                                                    7d8s19
       if(lprint)                                                       12d12s19
     $   write(6,*)('ndeto,nrooto returned from fetchcsf '),ndeto,nrooto
     $      ,ibcoff
       incsf=ibcoff                                                     6d18s19
       ibc(jptv+idod)=incsf                                             6d18s19
       ibcoff=incsf+2                                                   6d18s19
       ibc(incsf)=ndeto                                                 7d7s19
       if(lprint)                                                       12d12s19
     $      write(6,*)('storing ndet '),ndet(1),(' at incsf = '),incsf, 12d12s19
     $      ipt1(idod),incsf-ipt1(idod),idod
       ias=ibcoff                                                       7d7s19
       ibcoff=ias+nopen*ndeto                                           7d7s19
       if(lprint)write(6,*)('ias address: '),ias
       if(lprint)write(6,*)('nopen ndeto for ias '),nopen,ndeto,ias     12d12s19
       call enough('gencsf3.  1',bc,ibc)
       ivecq=ibcoff
       ibcoff=ivecq+1                                                   12d1s20
       if(lprint)write(6,*)('ivecq address '),ivecq,ndeto,nrooto        12d18s19
       if(lprint)write(6,*)('top of bc: '),ibcoff
       call enough('gencsf3.  2',bc,ibc)
       call fetchcsf2(nopen,ismultm,nrooto,bc(ivecq),ibc(ias),          12d1s20
     $      ibc(ias),ismultm,ndeto,bc,ibc)                              11d15s22
       if(nopen.gt.2)then                                               7d30s19
        idata(3)=nsames                                                 12d1s20
        idata(4)=nlows                                                  12d1s20
        if(ismultm.gt.0)then                                            12d1s20
         idata(5)=nsames                                                12d1s20
        else                                                            12d1s20
         idata(5)=0                                                     12d1s20
        end if                                                          12d1s20
        idata(6)=nhighs                                                 12d1s20
        if(lprint)                                                      12d12s19
     $    write(6,*)('for doubles: nsp, low,mid,high '),(idata(j),j=3,6)12d12s19
       else if(nopen.eq.2)then                                          8d1s19
        if(lprint)write(6,*)('nopen = 2, ismult = '),ismult             12d12s19
        if(ismult.eq.1)then                                             8d1s19
         idata(3)=1                                                     8d1s19
         do j=4,6
          idata(j)=0                                                    8d1s19
         end do                                                         8d1s19
        else if(ismult.eq.3)then                                        8d1s19
         do j=3,6                                                       8d1s19
          idata(j)=0                                                    8d1s19
         end do                                                         8d1s19
         idata(4)=1                                                     11d3s19
        end if                                                          8d1s19
       else if(nopen.eq.1)then                                          11d17s19
        if(lprint)write(6,*)('nopen = 1, ismult = '),ismult             12d12s19
        idata(3)=1                                                      11d17s19
        do j=4,6                                                        11d17s19
         idata(j)=0                                                     11d17s19
        end do                                                          11d17s19
       else if(nopen.eq.0)then                                          10d28s20
        idata(3)=1                                                      11d17s19
        do j=4,6                                                        10d28s20
         idata(j)=0                                                     10d28s20
        end do                                                          10d28s20
       end if                                                           7d30s19
       ntest=0                                                          12d1s20
       do j=1,4                                                         8d2s19
        ncsf2(j,idod)=idata(2+j)                                        8d2s19
        ntest=ntest+ncsf2(j,idod)                                       12d1s20
       end do                                                           8d2s19
       if(lprint)                                                       12d12s19
     $      write(6,*)('ncsf2 of '),idod,(' is '),(ncsf2(j,idod),j=1,4) 12d12s19
       ncsf(1)=nrooto                                                   6d18s19
       iaddd=nint(bc(ivecq))
       if(id.eq.mdon.and.iusewk.ne.0.and.nopen.gt.2)then                6d6s21
        idistrib=1                                                      4d30s21
        i18=ncsf(1)                                                     2d1s21
        i28=ndet(1)
        call ddi_create(bc,ibc,i18,i28,icsfxv)                          11d15s22
        call ilimts(1,ndet(1),mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e) 2d1s21
        nhere=ih+1-il                                                   2d1s21
        i18=1                                                           2d1s21
        i28=ncsf(1)                                                     2d1s21
        i38=il                                                          2d1s21
        i48=ih                                                          2d1s21
        call ddi_put(bc,ibc,icsfxv,i18,i28,i38,i48,bc(iaddd))           11d15s22
       else                                                             2d1s21
        idistrib=0                                                      4d30s21
        nhere=ndet(1)                                                   2d1s21
       end if                                                           2d1s21
       if(nrooto.ne.ntest)then                                          12d1s20
        write(6,*)('nrooto = '),nrooto,(' but ntest = '),ntest          12d1s20
        call dws_sync                                                   12d1s20
        call dws_finalize                                               12d1s20
        stop                                                            12d1s20
       end if                                                           12d1s20
       ibc(incsf+1)=ncsf(1)                                             6d18s19
       if(nopen.gt.2)then                                               8d4s19
       else                                                             8d5s19
       end if                                                           8d4s19
       ivecout=ivecq                                                    6d18s19
       ibcoff=ivecout+1                                                 12d1s20
       if(lprint)write(6,*)('resetting top of bc to '),ibcoff,ndet(1),
     $      ncsf(1),nhere
       do i=1,6                                                         7d30s19
        ibc(ibcoff)=idata(i)                                            7d30s19
        ibcoff=ibcoff+1                                                 7d30s19
       end do                                                           7d30s19
       icmp1=ibcoff                                                     11d27s19
       icmp2=icmp1+nhere+1                                              2d1s21
       icmp3=icmp2+nhere*ncsf(1)                                        2d1s21
       ibcoff=icmp3+nhere*ncsf(1)                                       2d1s21
       call enough('gencsf3.  3',bc,ibc)
       ivecop=nint(bc(ivecout))                                         12d1s20
       call cmpvcsf(nhere,ncsf,ibc(icmp1),ibc(icmp2),bc(icmp3),         2d1s21
     $      bc(ivecop),nused)                                           12d1s20
       if(idistrib.ne.0)then                                            4d30s21
        ixused=ibcoff                                                   2d1s21
        ibcoff=ixused+mynprocg                                          2d1s21
        call enough('gencsf3.  4',bc,ibc)
        do i=0,mynprocg-1                                               2d1s21
         bc(ixused+i)=0d0                                               2d1s21
        end do                                                          2d1s21
        bc(ixused+mynowprog)=dfloat(nused)                              2d1s21
        call dws_gbor(bc(ixused),mynprocg)                              2d1s21
        nhall=0                                                         2d1s21
        do i=0,mynprocg-1                                               2d1s21
         nh=nint(bc(ixused+i))                                          2d1s21
         nhall=nhall+nh                                                 2d1s21
        end do                                                          2d1s21
        nb4=0                                                           2d1s21
        do i=0,mynowprog-1                                              2d1s21
         nb4=nb4+nint(bc(ixused+i))                                     2d1s21
        end do                                                          2d1s21
        nhallh=nhall/2                                                  2d1s21
        if(nhallh*2.ne.nhall)nhallh=nhallh+1                            2d1s21
c
c     make sure there is space under icmp1 for everything               3d17s21
c
        ibcneed=icmp1+ndet(1)+1+nhallh+nhall                            3d17s21
        if(ibcneed.gt.ibcoff)then                                       3d17s21
         ibcoff=ibcneed                                                 3d17s21
         call enough('gencsf3.  5',bc,ibc)
        end if                                                          3d17s21
        icmp1a=ibcoff                                                   2d1s21
        icmp2a=icmp1a+ndet(1)+1                                         2d1s21
        icmp3a=icmp2a+nhallh                                            2d1s21
        ibcoff=icmp3a+nhall                                             2d1s21
        call enough('gencsf3.  6',bc,ibc)
        do i=0,ndet(1)                                                  2d1s21
         ibc(icmp1a+i)=0                                                2d1s21
        end do                                                          2d1s21
        jcmp1a=icmp1a+il-1                                              2d1s21
        do i=0,nhere-1                                                  2d1s21
         ipack8=ibc(icmp1+i)                                            2d1s21
         ibc(jcmp1a+i)=ibc(icmp1+i)                                     2d1s21
        end do                                                          2d1s21
        do i=0,nhallh-1                                                 2d1s21
         ibc(icmp2a+i)=0                                                2d1s21
        end do                                                          2d1s21
        call i4to8copy(ibc(icmp2),ibc(icmp2a),nb4,nused)                2d1s21
        do i=0,nhall-1                                                  2d1s21
         bc(icmp3a+i)=0d0                                               2d1s21
        end do                                                          2d1s21
        jcmp3a=icmp3a+nb4                                               2d1s21
        do i=0,nused-1                                                  2d1s21
         bc(jcmp3a+i)=bc(icmp3+i)                                       2d1s21
        end do                                                          2d1s21
        nwds=ibcoff-icmp1a                                              2d1s21
        call dws_gbor(ibc(icmp1a),nwds)                                 2d1s21
        iadd=0                                                          2d1s21
        last=0                                                          2d1s21
        do ii=0,ndet(1)-1                                                   2d1s21
         ipack8=ibc(icmp1a+ii)                                          2d1s21
         ia=ipack4(1)                                                   2d1s21
         ib=ipack4(2)                                                   2d1s21
         if(ipack4(1).lt.last)then                                      2d1s21
          iadd=iadd+last                                                2d1s21
         end if                                                         2d1s21
         last=ipack4(2)                                                 2d1s21
         ipack4(1)=ipack4(1)+iadd                                       2d1s21
         ipack4(2)=ipack4(2)+iadd                                       2d1s21
         ibc(icmp1a+ii)=ipack8                                          2d1s21
        end do
        do i=0,ndet(1)-1                                                2d1s21
         ibc(icmp1+i)=ibc(icmp1a+i)                                     2d1s21
        end do                                                          2d1s21
        ibc(icmp1+ndet(1))=nhall                                        2d1s21
        icmp2=icmp1+ndet(1)+1                                           2d1s21
        do i=0,nhallh-1                                                 2d1s21
         ibc(icmp2+i)=ibc(icmp2a+i)                                     2d1s21
        end do                                                          2d1s21
        icmp3=icmp2+nhall                                               2d1s21
        do i=0,nhall                                                    2d1s21
         bc(icmp3+i)=bc(icmp3a+i)                                       2d1s21
        end do                                                          2d1s21
        ibcoff=icmp3+nhall                                              2d1s21
       else                                                             2d1s21
        ibc(icmp1+ndet(1))=nused                                         12d3s20
        imovto=icmp2+nused                                               12d3s20
        do ii=0,nused-1                                                  12d3s20
         bc(imovto+ii)=bc(icmp3+ii)                                      12d3s20
        end do                                                           12d3s20
        icmp3=imovto                                                     12d3s20
        ibcoff=icmp3+nused                                               12d1s20
       end if                                                           2d1s21
      end do                                                            6d18s19
      ioff=ibctop0-ibcb4                                                12d6s20
      nmove=ibcoff-ibctop0                                              12d6s20
      do i=0,nmove-1                                                    12d6s20
       bc(ibcb4+i)=bc(ibctop0+i)                                        12d6s20
      end do                                                            12d6s20
      ibcoff=ibcoff-ioff                                                12d6s20
      if(iusewk.ne.-132)then                                            5d18s21
c
c     addressing ixw1. We will be given io and icase and we need to
c     point to index specifying full or compacted, then address if full
c     or no. words, then address if compacted.
c     this is essential triangle storage, so ((io*(io-1))/2)+icase      12d3s20
c     would do it.
c
      ixw1=ibcoff                                                       5d8s20
      nopenx=nec-mdon*2                                                 5d8s20
      nopenn=nec-mdoo*2                                                 5d8s20
      nxx=nopenx+1                                                      5d8s20
      nxxa=(nxx*(nxx+1))/2                                              12d3s20
      ibcoff=ixw1+nxxa                                                  12d3s20
      ixw2=ibcoff                                                       5d11s20
      ibcoff=ixw2+nxxa                                                  12d3s20
      ibctop=ibcoff                                                     12d5s20
      do io=nopenn,nopenx,2                                             5d11s20
       nclo=(nec-io)/2                                                  5d8s20
       nclop=nclo+1                                                     5d8s20
       iarg=nclop-mdon                                                  5d8s20
       ncsfsum=ncsf2(1,iarg)+ncsf2(2,iarg)+ncsf2(3,iarg)+ncsf2(4,iarg)  10d26s20
       if(ncsfx(iarg).ne.ncsfsum)then                                   10d26s20
        ncsfx(iarg)=ncsfsum                                             10d26s20
       end if                                                           10d26s20
       incsf=ipt1(iarg)+13-ioff                                         12d6s20
       ndeti=ibc(incsf)                                                 11d4s19
       ias=incsf+2                                                      11d4s19
       ivec=ias+io*ndeti                                                5d8s20
       ip1=ivec+1+6                                                     12d1s20
       ip2=ip1+ndeti                                                    12d3s20
       nused=ibc(ip2)                                                   12d3s20
       ip2=ip2+1                                                        12d3s20
       ip3=ip2+nused                                                    12d3s20
       nx=io+1
       do i=1,io                                                        5d8s20
        icode(i)=3                                                      5d8s20
       end do                                                           5d8s20
       icode(nx)=2                                                      5d8s20
       nnxx=ncsfx(iarg)*ncsfx(iarg)*2                                   12d1s20
       i18=nnxx+1                                                       12d3s20
       i28=io                                                           12d1s20
       call ddi_create(bc,ibc,i18,i28,i58)                              11d15s22
       call ilimts(1,io,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)       12d1s20
       icasen=ibcoff                                                    12d3s20
       ibcoff=icasen+io                                                 12d3s20
       ioutbc=ibcoff                                                    12d1s20
       do icase=1,io                                                    5d8s20
        icode(icase)=1                                                  5d8s20
        iooa=ixw1+((io*(io-1))/2)+icase-1                               12d3s20
        ibc(iooa)=ibcoff                                                12d3s20
        call second(timeq)
        if(icase.ge.il.and.icase.le.ih)then                             12d1s20
         inused=ibcoff                                                  12d3s20
         ibcoff=ibcoff+1                                                12d3s20
         call cupo21(icode,nx,                                           5d11s20
     $       ibc(ip1),ibc(ip2),bc(ip3),ndeti,ncsfx(iarg),ibc(ias),io,    5d8s20
     $       ibc(ip1),ibc(ip2),bc(ip3),ndeti,ncsfx(iarg),ibc(ias),io,    5d11s20
     $       iout,nnot,bc,ibc)                                          11d15s22
         ioutt=ibcoff                                                    5d11s20
         ibcoff=ioutt+ncsfx(iarg)*ncsfx(iarg)                            5d11s20
         call enough('gencsf3.  7',bc,ibc)
         do i=0,ncsfx(iarg)-1                                            5d11s20
          do j=0,ncsfx(iarg)-1                                           5d11s20
           ji=iout+j+ncsfx(iarg)*i                                       5d11s20
           ij=ioutt+i+ncsfx(iarg)*j                                      5d11s20
           bc(ij)=bc(ji)                                                 5d11s20
          end do                                                         5d11s20
         end do                                                          5d11s20
         icmp1=ibcoff
         icmp2=icmp1+ncsfx(iarg)
         icmp3=icmp2+ncsfx(iarg)*ncsfx(iarg)
         ibcoff=icmp3+ncsfx(iarg)*ncsfx(iarg)
         call enough('gencsf3.  8',bc,ibc)
         call cmpvcsf(ncsfx(iarg),ncsfx(iarg),ibc(icmp1),ibc(icmp2),
     $       bc(icmp3),bc(iout),nused)                                  12d1s20
         nusedi=nused/2                                                 12d4s20
         if(nusedi*2.ne.nused)nusedi=nusedi+1                           12d4s20
         nused4=2*(nusedi+nused+ncsfx(iarg))+1                          12d4s20
         if(nused4.lt.nnxx)then                                         12d4s20
          ibc(inused)=nused                                             12d4s20
          imovto=icmp2+nusedi                                           12d4s20
          do ii=0,nused-1                                               12d4s20
           bc(imovto+ii)=bc(icmp3+ii)                                   12d4s20
          end do                                                        12d4s20
          icmp1t=imovto+nused                                            12d4s20
          ibc(icmp1t)=nused                                             12d4s20
          icmp1t=icmp1t+1                                               12d4s20
          icmp2t=icmp1t+ncsfx(iarg)
          icmp3t=icmp2t+ncsfx(iarg)*ncsfx(iarg)
          ibcoff=icmp3t+ncsfx(iarg)*ncsfx(iarg)
          call enough('gencsf3.  9',bc,ibc)
          call cmpvcsf(ncsfx(iarg),ncsfx(iarg),ibc(icmp1t),ibc(icmp2t),
     $       bc(icmp3t),bc(ioutt),nusedt)                                  12d1s20
          if(nusedt.ne.nused)then
           write(6,*)('nusedt = '),nusedt,(' is not equal to nused '),
     $          nused
           call dws_synca
           call dws_finalize
           stop
          end if                                                        12d4s20
          imovto=icmp2t+nusedi                                          12d4s20
          do ii=0,nused-1                                               12d4s20
           bc(imovto+ii)=bc(icmp3t+ii)                                  12d4s20
          end do                                                        12d4s20
          do ii=0,nused4-1                                              12d4s20
           bc(iout+ii)=bc(icmp1+ii)                                     12d4s20
           ipack8=ibc(icmp1+ii)                                         12d4s20
          end do                                                        12d4s20
         else                                                           12d4s20
          ibc(inused)=0                                                 12d4s20
         end if                                                         12d4s20
         ibcoff=icmp1                                                   12d4s20
         i18=1                                                          12d1s20
         i28=nnxx+1                                                     12d3s20
         i38=icase                                                      12d1s20
         call ddi_put(bc,ibc,i58,i18,i28,i38,i38,bc(inused))            11d15s22
        else                                                            12d1s20
         iout=ibcoff                                                    12d1s20
         ibc(iout)=0                                                    7d18s24
         ibcoff=iout+nnxx+1                                             12d3s20
         call enough('gencsf3. 10',bc,ibc)
        end if                                                          12d1s20
        icode(icase)=3                                                  5d8s20
       end do                                                           5d8s20
       call dws_synca                                                   12d1s20
       i18=1                                                            12d1s20
       i28=nnxx+1                                                       12d3s20
       i38=1                                                            12d1s20
       i48=io                                                           12d1s20
       call ddi_get(bc,ibc,i58,i18,i28,i38,i48,bc(ioutbc))              11d15s22
       call dws_synca
       call ddi_destroy(i58)                                            12d1s20
       call dws_synca
       joutbc=ioutbc
       do icase=1,io                                                    12d4s20
        iooa=ixw1+((io*(io-1))/2)+icase-1                               12d3s20
        itest=ibc(iooa)
        joutbc0=joutbc
        if(ibc(itest).eq.0)then                                         12d4s20
         do j=0,nnxx                                                    12d4s20
          bc(joutbc+j)=bc(itest+j)                                      12d4s20
         end do                                                         12d4s20
         ibc(iooa)=joutbc                                               12d4s20
         joutbc=joutbc+nnxx+1                                           12d4s20
        else
         nused=ibc(itest)
         nusedi=nused/2
         if(nusedi*2.ne.nused)nusedi=nusedi+1                           12d4s20
         need=2*(nused+nusedi+ncsfx(iarg))+2                            12d4s20
         do j=0,need                                                    12d4s20
          bc(joutbc+j)=bc(itest+j)                                      12d4s20
         end do                                                         12d4s20
         ibc(iooa)=joutbc                                               12d4s20
         joutbc=joutbc+need+1                                           12d4s20
        end if                                                          12d4s20
       end do                                                           12d4s20
       ibcoff=joutbc                                                    12d4s20
      end do                                                            5d8s20
c
c     addressing ixw2. We will be given io and icase and we need to
c     point to index specifying full or compacted, then address if full
c     or no. words, then address if compacted.
c     this is essential triangle storage, so ((iop*(iop-1))/2)+icase    12d3s20
c     with iop=io+1                                                     12d3s20
c
      call enough('gencsf3. 11',bc,ibc)
      do io=nopenn,nopenx-1,2                                           5d11s20
       iop=io+2
       nclo=(nec-io)/2                                                  5d8s20
       nclop=nclo+1                                                     5d8s20
       iarg=nclop-mdon                                                  5d8s20
       jarg=iarg-1
       incsf=ipt1(iarg)+13-ioff                                         12d6s20
       ndeti=ibc(incsf)                                                 11d4s19
       ias=incsf+2                                                      11d4s19
       ivec=ias+io*ndeti                                                5d8s20
       ip1=ivec+1+6                                                     12d1s20
       ip2=ip1+ndeti                                                    11d29s19
       nused=ibc(ip2)                                                   12d3s20
       ip2=ip2+1                                                        12d3s20
       ip3=ip2+nused                                                    12d3s20
       incsf=ipt1(jarg)+13-ioff                                         12d6s20
       ndetj=ibc(incsf)                                                 11d4s19
       jas=incsf+2                                                      11d4s19
       jvec=jas+iop*ndetj                                                5d8s20
       jp1=jvec+1+6                                                     12d1s20
       jp2=jp1+ndetj                                                    11d29s19
       nused=ibc(jp2)                                                   12d3s20
       jp2=jp2+1                                                        12d3s20
       jp3=jp2+nused                                                    12d3s20
       nx=iop                                                           5d11s20
       nxm=nx-1                                                         5d11s20
       do i=1,nx                                                        5d11s20
        icode(i)=3                                                      5d11s20
       end do                                                           5d11s20
       icode(nx)=1                                                      5d11s20
       nnxx=ncsfx(jarg)*ncsfx(iarg)*2                                   12d1s20
       i18=nnxx+1                                                       12d4s20
       i28=nxm                                                          12d1s20
       call ddi_create(bc,ibc,i18,i28,i58)                              11d15s22
       call ilimts(1,nxm,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)      12d1s20
       ioutbc=ibcoff                                                    12d1s20
       do icase=1,nxm                                                   5d11s20
        iooa=ixw2+((io*(io+1))/2)+icase-1                               12d3s20
        ibc(iooa)=ibcoff                                                12d3s20
        icode(icase)=7                                                  5d11s20
        call second(timeq)
        if(icase.ge.il.and.icase.le.ih)then                             12d1s20
         inused=ibcoff                                                  12d3s20
         ibcoff=ibcoff+1                                                12d3s20
         call cupo21(icode,nx,                                           5d11s20
     $       ibc(jp1),ibc(jp2),bc(jp3),ndetj,ncsfx(jarg),ibc(jas),iop,    5d8s20
     $       ibc(ip1),ibc(ip2),bc(ip3),ndeti,ncsfx(iarg),ibc(ias),io,    5d11s20
     $       iout,nnot,bc,ibc)                                          11d15s22
         ioutt=ibcoff                                                    5d11s20
         ibcoff=ioutt+ncsfx(jarg)*ncsfx(iarg)                            5d11s20
         call enough('gencsf3. 12',bc,ibc)
         do i=0,ncsfx(iarg)-1                                            5d11s20
          do j=0,ncsfx(jarg)-1                                           5d11s20
           ji=iout+j+ncsfx(jarg)*i                                       5d11s20
           ij=ioutt+i+ncsfx(iarg)*j                                      5d11s20
           bc(ij)=bc(ji)
          end do                                                         5d11s20
         end do                                                          5d11s20
         icmp1=ibcoff
         icmp2=icmp1+ncsfx(jarg)                                        12d4s20
         icmp3=icmp2+ncsfx(iarg)*ncsfx(jarg)
         ibcoff=icmp3+ncsfx(iarg)*ncsfx(jarg)
         call enough('gencsf3. 13',bc,ibc)
         call cmpvcsf(ncsfx(jarg),ncsfx(iarg),ibc(icmp1),ibc(icmp2),
     $       bc(icmp3),bc(iout),nused)                                  12d1s20
         nusedi=nused/2
         if(nusedi*2.ne.nused)nusedi=nusedi+1                           12d4s20
         nused4=ncsfx(jarg)+ncsfx(iarg)+2*(nused+nusedi)+1              12d4s20
         if(nused4.lt.nnxx)then                                         12d4s20
          ibc(inused)=nused                                             12d4s20
          imovto=icmp2+nusedi                                           12d4s20
          do ii=0,nused-1                                               12d4s20
           bc(imovto+ii)=bc(icmp3+ii)                                   12d4s20
          end do                                                        12d4s20
          icmp1t=imovto+nused                                            12d4s20
          ibc(icmp1t)=nused                                             12d4s20
          icmp1t=icmp1t+1                                               12d4s20
          icmp2t=icmp1t+ncsfx(iarg)                                     12d4s20
          icmp3t=icmp2t+ncsfx(iarg)*ncsfx(jarg)
          ibcoff=icmp3t+ncsfx(iarg)*ncsfx(jarg)
          call enough('gencsf3. 14',bc,ibc)
          call cmpvcsf(ncsfx(iarg),ncsfx(jarg),ibc(icmp1t),ibc(icmp2t),
     $       bc(icmp3t),bc(ioutt),nusedt)                                  12d1s20
          if(nusedt.ne.nused)then
           write(6,*)('nusedt = '),nusedt,(' is not equal to nused '),
     $          nused
           call dws_synca
           call dws_finalize
           stop
          end if                                                        12d4s20
          imovto=icmp2t+nusedi                                          12d4s20
          do ii=0,nused-1                                               12d4s20
           bc(imovto+ii)=bc(icmp3t+ii)                                  12d4s20
          end do                                                        12d4s20
          do ii=0,nused4-1                                              12d4s20
           bc(iout+ii)=bc(icmp1+ii)                                     12d4s20
           ipack8=ibc(icmp1+ii)                                         12d4s20
          end do                                                        12d4s20
         else                                                           12d4s20
          ibc(inused)=0                                                 12d4s20
         end if                                                         12d4s20
         ibcoff=icmp1
         i18=1                                                          12d1s20
         i28=nnxx+1                                                     12d4s20
         i38=icase                                                      12d1s20
         call ddi_put(bc,ibc,i58,i18,i28,i38,i38,bc(inused))            11d15s22
        else                                                            12d1s20
         iout=ibcoff                                                    12d1s20
         ibc(iout)=0                                                    7d18s24
         ibcoff=iout+nnxx+1                                             12d4s20
         call enough('gencsf3. 15',bc,ibc)
        end if                                                          12d1s20
        icode(icase)=3                                                  5d11s20
       end do                                                           5d11s20
       call dws_synca                                                   12d1s20
       i18=1                                                            12d1s20
       i28=nnxx+1                                                       12d4s20
       i38=1                                                            12d1s20
       i48=nxm                                                          12d1s20
       call ddi_get(bc,ibc,i58,i18,i28,i38,i48,bc(ioutbc))              11d15s22
       call dws_synca
       call ddi_destroy(i58)                                            12d1s20
       call dws_synca
       joutbc=ioutbc
       do icase=1,nxm                                                   12d4s20
        iooa=ixw2+((io*(io+1))/2)+icase-1                               12d3s20
        itest=ibc(iooa)
        joutbc0=joutbc
        if(ibc(itest).eq.0)then                                         12d4s20
         do j=0,nnxx                                                    12d4s20
          bc(joutbc+j)=bc(itest+j)                                      12d4s20
         end do                                                         12d4s20
         ibc(iooa)=joutbc                                               12d4s20
         joutbc=joutbc+nnxx+1                                           12d4s20
        else
         nused=ibc(itest)
         nusedi=nused/2
         if(nusedi*2.ne.nused)nusedi=nusedi+1                           12d4s20
         need=2*(nused+nusedi)+ncsfx(jarg)+ncsfx(iarg)+2                12d4s20
         do j=0,need                                                    12d4s20
          bc(joutbc+j)=bc(itest+j)                                      12d4s20
          ipack8=ibc(joutbc+j)
         end do                                                         12d4s20
         ibc(iooa)=joutbc                                               12d4s20
         joutbc=joutbc+need+1                                           12d4s20
        end if                                                          12d4s20
       end do                                                           12d4s20
       ibcoff=joutbc                                                    12d4s20
      end do                                                            5d11s20
      end if                                                            10d28s20
      ioff=ibcb4-ixw1                                                   12d6s20
      nmove=ibcoff-ixw1                                                 12d6s20
      do i=0,nmove-1                                                    12d6s20
       bc(ibcb4+i)=bc(ixw1+i)                                           12d6s20
      end do                                                            12d6s20
      ibcoff=ibcb4+nmove                                                12d6s20
      ixw1=ibcb4                                                        12d6s20
      ixw2=ixw1+nxxa                                                    12d6s20
      do io=nopenn,nopenx,2                                             5d11s20
       do icase=1,io                                                    5d8s20
        iooa=ixw1+((io*(io-1))/2)+icase-1                               12d3s20
        ibc(iooa)=ibc(iooa)+ioff                                        12d6s20
       end do                                                           12d6s20
      end do
      do io=nopenn,nopenx-1,2                                           5d11s20
       do icase=1,io+1                                                  12d6s20
        iooa=ixw2+((io*(io+1))/2)+icase-1                               12d3s20
        ibc(iooa)=ibc(iooa)+ioff                                        12d6s20
       end do                                                           12d6s20
      end do                                                            12d6s20
      ibcb4=ixw2+nxxa                                                   12d7s20
      return
      end
