c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine gencsf2(nek,spin,norb,idorb,isorb,iptr,nsymb,mdon,mdoo,7d15s19
     $     itmp,multh,ism,ncsf,nok,ntot,nod,nos,maxopens,nhole,ihole,   7d15s19
     $     nfill,ifill,idox,iflag,itestmrci,ismult,isymuse,nfcnr,idorbr,10d1s19
     $     isorbr,iptrr,ibasisr,icsfpd,nfcnx,interacts,nodu,nosu,nsymbl,1d5s20
     $     iuncd,nfcnt,iptrt,ibasist,irxinfo,isymmrci,ifcnpx,iptrbit,   10d19s20
     $     ifcnpxr,ncsf2,iusehole,bc,ibc)                               11d14s22
      implicit real*8 (a-h,o-z)
      real*16 fl                                                        7d15s19
      logical ldebug                                                    12d12s19
      integer*8 itesta,itestb                                           10d18s20
c
c     iflag = 0, do this for reference space
c     iflag = 1, do this for singles
c     iflag = 2, do this for doubles.
c     iflag = 3, do this for triples.
c     iflag = -1, do this for full reference space
c
      integer*1 idorb(*),isorb(*),iorb(64),idorbr(*),isorbr(*),         10d1s19
     $     iorbt(64,3),nab(2),iorbfill(64)                              10d25s20
      dimension iptr(4,mdoo+1,nsymb),multh(8,8),ism(norb),ncsf(*),      9d5s19
     $    itmp(*),iptrr(4,*),ibasisr(3,*),icsfpd(*),ifcnpx(2,4,nsymb,*),10d19s20
     $     ihole(idox,3),ifill(idox,*),nhole(3),nfill(4),nok(*),ntot(*),4d19s23
     $     nfcnt(*),iptrt(4,mdoo+1,nsymb),ibasist(*),irxinfo(*),        10d18s20
     $     iptrbit(2,mdoo+1,*),ifcnpxr(2,4,nsymb,*),ncsf2(4,*),nsum(4), 11d13s20
     $     nab4(2,2)                                                    11d13s20
      COMMON/FACT16/FL(922),NCALL                                       5d27s19
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      data icall/0/                                                     10d19s20
      include "common.store"                                            10d2s19
      include "common.print"                                            1d3s20
      save icall                                                        10d19s20
      icall=icall+1
      norbx=norb+1                                                      10d18s20
      isorbx=0                                                          10d3s24
      idorbx=0                                                          10d2s24
      nfillt=0                                                          9d3s24
      do i=1,4                                                          9d3s24
       if(nfill(i).ne.0)nfillt=i                                        9d3s24
      end do                                                            9d3s24
      if(iprtr(3).eq.0)then                                             1d3s20
       ldebug=.false.                                                    12d12s19
      else                                                              1d3s20
       ldebug=.true.
      end if                                                            1d3s20
      nodu=0                                                            11d24s19
      nosu=0                                                            11d24s19
      if(iflag.ge.0)then                                                10d1s19
       nec=nek-iflag                                                     7d15s19
      else                                                              10d1s19
       nec=nek                                                          10d1s19
      end if                                                            10d1s19
      if(ldebug)then                                                    12d12s19
       write(6,*)('in gencsf2 '),nec,mdon,mdoo
       write(6,*)('nsymb, nsymbl: '),nsymb,nsymbl                       12d13s19
       write(6,*)('iflag = '),iflag,icall
      end if                                                            12d12s19
      isn=nec-2*mdoo
      if(isn.lt.0)isn=mod(nec,2)                                        9d5s19
      isx=min(norb,maxopens,nec-2*mdon)                                 9d5s19
      if(mod(isx,2).ne.mod(nec,2))isx=isx-1                             9d5s19
      if(iflag.eq.0.or.iflag.eq.-1)then                                 10d1s19
       ismultm=ismult-1                                                 9d5s19
       if(ldebug)write(6,*)('isn, vs ismultm: '),isn,ismultm            12d12s19
       isn=max(isn,ismultm)                                             9d5s19
      else if(iflag.eq.1)then                                           9d5s19
       ispin2=ismult-1                                                  9d5s19
       ispin2l=iabs(ispin2-1)                                           9d5s19
       if(ldebug)write(6,*)('for singles, ispin2l = '),ispin2l,isn      12d12s19
       isn=max(isn,ispin2l)                                             9d5s19
      else if(iflag.eq.2)then                                           9d5s19
       ispin2=ismult-1                                                  9d5s19
       ispin2l=min(ispin2,iabs(ispin2-2))                               9d5s19
       if(ldebug)write(6,*)('for doubles, ispin2l = '),ispin2l,isn      12d12s19
       isn=max(isn,ispin2l)                                             9d5s19
      else if(iflag.eq.3)then                                           7d17s20
       ispin2=ismult-1                                                  9d5s19
       ispin2l=min(iabs(ispin2-1),iabs(ispin2-3))                       10d20s20
       if(ldebug)write(6,*)('for triples, ispin2l = '),ispin2l,isn      7d17s20
       isn=max(isn,ispin2l)                                             9d5s19
      end if                                                            9d5s19
      mxc=(nec-isn)/2                                                   9d5s19
      if(ldebug)then                                                    12d12s19
       write(6,*)('minimum no. of open shells'),isn
       write(6,*)('maximum no. of open shells'),isx,nec,mdon,norb
       write(6,*)('maximum no. of closed shells: '),mxc                  9d5s19
       write(6,*)('ism: '),ism
      end if                                                            12d12s19
      do i=1,nsymbl                                                     12d13s19
       do j=1,mdoo+1                                                    7d15s19
        do k=1,4                                                        7d15s19
         iptr(k,j,i)=0                                                  7d15s19
        end do                                                          7d15s19
       end do                                                           7d15s19
      end do                                                            7d15s19
      jorb=1                                                            6d12s19
      khere=1                                                           7d15s19
      ibasis=0                                                          7d15s19
      do isb=1,nsymbl                                                   12d13s19
       if(iflag.lt.0)then                                               10d1s19
        isbu=isymuse                                                    10d1s19
       else                                                             10d1s19
        isbu=isb                                                        10d1s19
       end if                                                           10d1s19
       loop=0                                                           11d19s20
       if(ldebug)write(6,*)('for symmetry '),isbu                       12d12s19
       if(mdon.eq.0)then
        if(ldebug)write(6,*)('storing data for no closed')              12d12s19
        iptr(1,1,isb)=1                                                 7d15s19
        iptr(2,1,isb)=1                                                 7d15s19
       end if
       if(isn.eq.0)then
        if(ldebug)write(6,*)('storing data for no open')                12d12s19
        iptr(4,mxc+1,isb)=1                                             9d5s19
        if(isbu.eq.1)then                                               11d3s19
         iptr(3,mxc+1,isb)=1                                            9d5s19
        else                                                            7d15s19
         iptr(3,mxc+1,isb)=0                                            9d5s19
        end if                                                          7d15s19
       end if
       do i=1,norb
        itmp(i+ibasis)=i                                                7d15s19
       end do
       n1c=norb-nfill(1)                                                10d25s20
       nthis=norb
       if(mdon.le.1.and.(nec.gt.2.or.(isb.eq.1.and.nec.eq.2)).and.      1d3s20
     $      mdoo.gt.0)then                                              1d3s20
        if(ldebug)write(6,*)('storing data for one closed'),jorb,nthis  12d12s19
        iptr(1,2,isb)=nthis                                             7d15s19
        iptr(2,2,isb)=jorb                                              7d15s19
        if(jorb+nthis.gt.2*nod*nsymb)then                               1d2s20
         write(6,*)('idorb dims exceeded! '),nthis,nod,jorb,jorb+nthis,
     $       nod*nsymb
         call dws_sync
         call dws_finalize
         stop 'gencsf2'
        end if                                                           6d12s19
        do i=1,nthis
         idorb(i+jorb-1)=itmp(i+ibasis)                                 7d15s19
         idorbx=max(idorbx,i+jorb-1)
        end do
        jorb=nthis+jorb                                                  7d15s19
       end if
       if(isn.eq.1)then
        iptr(4,mxc+1,isb)=khere                                         9d5s19
        nhere=0
        do i=1,norb
         if(ism(i).eq.isbu)then                                         10d1s19
          isorb(khere+nhere)=i                                           7d15s19
          isorbx=max(isorbx,khere+nhere)                                8d28s24
          nhere=nhere+1
          if(nhere+khere.gt.nos*nsymb)then                              7d15s19
           write(6,*)('isorb dimensions exceededa! '),nhere,nos,khere,  7d15s19
     $         nhere+khere,nos*nsymb                                    7d15s19
           call dws_sync                                                 6d12s19
           call dws_finalize                                             6d12s19
           stop 'gencsf2'                                                         6d12s19
          end if                                                         6d12s19
         end if
        end do
        iptr(3,mxc+1,isb)=nhere                                         9d5s19
        if(ldebug)write(6,*)('storing data for one open'),khere,nhere   12d12s19
        khere=nhere+khere                                               7d15s19
       end if
       do ido=2,max(mxc,isx)                                            9d5s19
        idol=ido-1
        if(ldebug)write(6,*)('generating strings for '),ido,            12d12s19
     $       (' orbitals')                                              12d12s19
     $      ,nthis
        itop=nthis*idol+1
        iadd=0
        do i=0,nthis-1
         ifirst=i*idol
         last=ifirst+idol
         last=itmp(last+ibasis)                                         7d15s19
         ngfill=0                                                       10d3s24
         if(nfillt.gt.0)then                                            9d3s24
          do k=1,idol                                                   9d3s24
           do j=1,nfill(nfillt)                                         9d3s24
            if(itmp(ifirst+k+ibasis).eq.ifill(j,nfillt))ngfill=ngfill+1 9d3s24
           end do                                                       9d3s24
          end do                                                        9d3s24
          if(ngfill.ne.nfillt)ngfill=0                                  9d3s24
         end if                                                         9d3s24
         do j=last+1,norb
          if(ngfill.ne.0)then                                           10d3s24
           do k=1,nfill(nfillt)                                         9d3s24
            if(ifill(k,nfillt).eq.j)then                                9d3s24
             go to 1066                                                 9d3s24
            end if                                                      9d3s24
           end do                                                       9d3s24
          end if                                                        9d3s24
          ibeg=itop
          do k=1,idol
           itmp(itop+ibasis)=itmp(ifirst+k+ibasis)                      7d15s19
           itop=itop+1
          end do
          iadd=iadd+1
          itmp(itop+ibasis)=j                                           7d15s19
          itop=itop+1
 1066     continue                                                      9d3s24
         end do
        end do
        itop=nthis*idol
        do i=1,iadd*ido
         itmp(i+ibasis)=itmp(itop+i+ibasis)                             7d15s19
        end do
        if(ido.ge.mdon.and.ido.le.mxc)then                              9d5s19
         if(ldebug)write(6,*)('store closed '),iadd,jorb,nod,ido+1      12d12s19
         iptr(2,ido+1,isb)=jorb                                         7d15s19
         jorbm=jorb-1
         jtop=jorbm+iadd*ido                                             6d12s19
         if(jtop.gt.nod*nsymb)then                                      7d15s19
          write(6,*)('idorb dims exceeded! '),jtop,nod,nsymb
          call dws_sync
          call dws_finalize
          stop 'gencsf2'
         end if                                                           6d12s19
         nuse=0                                                         7d19s19
         do i=0,iadd-1                                                  7d19s19
          iad1=ibasis+ido*i+1                                           7d19s19
          do j=1,norb                                                   7d19s19
           iorb(j)=0                                                    7d19s19
          end do                                                        7d19s19
          do j=0,ido-1                                                  8d5s19
           iom=itmp(iad1+j)                                             7d19s19
           iorb(iom)=2                                                  7d19s19
          end do                                                        7d19s19
          do j=1,4                                                      4d19s23
           jm=j-1                                                       7d19s19
           minthem=0                                                    7d19s19
           do k=1,nfill(j)                                              7d19s19
            iom=ifill(k,j)                                              7d19s19
            minthem=minthem+iorb(iom)                                   7d19s19
           end do                                                       7d19s19
           if(minthem.gt.j)then                                         8d5s19
            go to 1676                                                  7d19s19
           end if                                                       8d5s19
          end do                                                        7d19s19
          iad2=jorb+ido*nuse                                            7d19s19
          nuse=nuse+1                                                   7d19s19
          do j=0,ido-1                                                  7d19s19
           idorb(iad2+j)=itmp(iad1+j)                                   7d19s19
           idorbx=max(idorbx,iad2+j)
          end do                                                        7d19s19
 1676     continue                                                      7d19s19
         end do                                                         7d19s19
         iptr(1,ido+1,isb)=nuse                                         7d19s19
         jorb=jorb+nuse*ido                                             7d19s19
         if(ldebug)write(6,*)('jorb so far '),jorb                      12d12s19
        end if
        if(ido.ge.isn.and.mod(ido+isn,2).eq.0)then
         idou=(nec-ido)/2
         khere0=khere                                                   7d15s19
         iptr(4,idou+1,isb)=khere                                       7d15s19
         nhere=0
         do i=0,iadd-1
          ip=i+1
          ifirst=i*ido
          istry=ism(itmp(ifirst+1+ibasis))                              7d15s19
          do j=2,ido
           istry=multh(istry,ism(itmp(ifirst+j+ibasis)))                7d15s19
          end do
          if(istry.eq.isbu)then                                         10d1s19
           jtop=khere+ido-1                                              6d12s19
           if(jtop.gt.nos*nsymb)then                                    7d15s19
            write(6,*)('isorb dimensions exceededb! '),jtop,nos           6d12s19
            call dws_sync                                                 6d12s19
            call dws_finalize                                             6d12s19
            stop 'gencsf2'                                                          6d12s19
           end if                                                         6d12s19
           do j=1,norb                                                  7d19s19
            iorb(j)=0                                                   7d19s19
           end do                                                       7d19s19
           do j=1,ido                                                   7d19s19
            iom=itmp(ifirst+j+ibasis)                                   7d19s19
            iorb(iom)=1                                                 7d19s19
           end do                                                       7d19s19
           do j=1,4                                                     4d19s23
            jm=j-1                                                      7d19s19
            minthem=0                                                   7d19s19
            do k=1,nfill(j)                                             7d19s19
             iom=ifill(k,j)                                             7d19s19
             minthem=minthem+iorb(iom)                                  7d19s19
            end do                                                      7d19s19
            if(minthem.gt.j)go to 310                                   7d19s19
           end do                                                       7d19s19
           kherex=khere                                                 7d19s19
           do j=1,ido
            isorb(khere)=itmp(ifirst+j+ibasis)                          7d15s19
            isorbx=max(isorbx,khere)                                    8d28s24
            khere=khere+1
           end do
           nhere=nhere+1
  310      continue                                                     7d19s19
          end if
         end do
         iptr(3,idou+1,isb)=nhere                                       7d15s19
         if(ldebug)write(6,*)('store open shells under index '),idou+1, 12d12s19
     $        nhere,khere0                                              12d12s19
        end if
        nthis=iadd
       end do
       nok(isb)=0                                                       7d15s19
       isto=norb+ibasis+1                                               7d15s19
       do ido=mdon,mxc                                                  9d5s19
        nsing=nec-2*ido                                                  6d5s19
        if(nsing.le.isx)then                                            10d18s20
         if(ldebug)write(6,*)('for '),ido,(' closed shells orbitals '),  12d12s19
     $     ('we have '),iptr(1,ido+1,isb),('closed combinations, and '),7d15s19
     $      iptr(3,ido+1,isb),('open combinations'),ido+1                     7d15s19
         nsing=nec-2*ido                                                  6d5s19
         ios=iptr(4,ido+1,isb)                                           7d15s19
         do is=1,iptr(3,ido+1,isb)                                       7d15s19
          iod=iptr(2,ido+1,isb)                                          7d15s19
          do id=1,iptr(1,ido+1,isb)                                      7d15s19
           do i=1,norb                                                   7d16s19
            itmp(i+ibasis)=0                                             7d16s19
           end do                                                        7d16s19
           do i=0,nsing-1                                                7d16s19
            itmp(isorb(ios+i)+ibasis)=1                                  7d16s19
           end do                                                        7d16s19
           do i=0,ido-1
            if(itmp(idorb(iod+i)+ibasis).eq.1)then                       7d15s19
             go to 10
            else                                                         7d15s19
             itmp(idorb(iod+i)+ibasis)=2                                 7d15s19
            end if
           end do
c
c     nhole(1) no holes <- this will correspond to core orbs
c                          correlated, i.e. must always be in closed    7d15s19
c     nhole(2) one hole i.e. must always be in open or closed           7d15s19
c     nhole(3) two holes
c     ihole(j,i) is global orb no. for nhole(i)
c     nfill(1) max 1 e i.e. either empty or in open
c     nfill(2) max 2 e i.e. either empty, in open, or in double.
c     ifill(j,i) is global orb no. for nfile(i)
c
           if(iusehole.ne.0)then                                        10d25s20
            do i=1,3                                                      7d15s19
             im=i-1                                                       7d15s19
             minthem=0                                                    7d15s19
             do j=1,nhole(i)                                              7d15s19
              minthem=minthem+itmp(ihole(j,i)+ibasis)                     7d15s19
             end do                                                       7d15s19
             mgoal=nhole(i)*2-im                                          9d18s19
             if(iflag.ne.0)mgoal=mgoal-3                                 10d25s20
             if(minthem.lt.mgoal)go to 10                                 7d15s19
            end do                                                        7d15s19
            do i=1,4                                                    4d19s23
             minthem=0                                                    7d15s19
             do j=1,nfill(i)                                              7d15s19
              minthem=minthem+itmp(ifill(j,i)+ibasis)                     7d15s19
             end do                                                       7d15s19
             if(minthem.gt.i)then                                         8d5s19
              go to 10                                                    7d15s19
             end if                                                       8d5s19
            end do                                                        7d15s19
           end if                                                       10d25s20
           if(iflag.ne.0)then                                            10d1s19
            itestb=0                                                    10d18s20
            do i=1,nsing                                                 10d1s19
             im=i-1                                                      10d1s19
             iorbt(i,1)=isorb(ios+im)                                    10d1s19
             itestb=ibset(itestb,isorb(ios+im))                         10d18s20
            end do                                                       10d1s19
            notest=nsing                                                 10d1s19
            if(iflag.gt.0)then                                           10d1s19
             notest=notest+1                                             10d1s19
             iorbt(notest,1)=63                                          10d1s19
             if(iflag.eq.2)then                                          10d1s19
              notest=notest+1                                            10d1s19
              iorbt(notest,1)=64                                         10d1s19
             end if                                                      10d1s19
             itestb=ibset(itestb,norbx)                                 10d18s20
            end if                                                       10d1s19
            itesta=0                                                      10d18s20
            do i=0,ido-1                                                  10d18s20
             itesta=ibset(itesta,idorb(iod+i))                            10d18s20
            end do                                                        10d18s20
            iargt=ido+1-mdon                                             10d2s19
            incsf=icsfpd(iargt)+13                                       10d2s19
            ndett=ibc(incsf)                                             10d2s19
            iast=incsf+2                                                 10d2s19
            ivect=iast+notest*ndett                                      10d2s19
            if(interacts.eq.0)go to 1010                                10d18s20
            idnn=max(1,ido-1)                                           10d19s20
            idxx=min(ido+3,mdoo+1)                                      10d19s20
            do ifn=ifcnpxr(1,1,isymmrci,idnn),ifcnpxr(2,1,isymmrci,idxx)10d19s20
             nclor=ibasisr(1,ifn)                                        10d1s19
             nopenr=nek-2*nclor                                          10d1s19
             nclorp=nclor+1                                              10d1s19
             iargr=nclorp-mdon                                           10d1s19
             irc=iptrr(2,nclorp)+nclor*(ibasisr(2,ifn)-1)                10d1s19
             iro=iptrr(4,nclorp)+nopenr*(ibasisr(3,ifn)-1)               10d1s19
             iirc=iptrbit(1,nclorp,isymmrci)+ibasisr(2,ifn)-1           10d18s20
             iiro=iptrbit(2,nclorp,isymmrci)+ibasisr(3,ifn)-1           10d18s20
             if(interacts.ne.0)then                                      10d15s19
              icall=icall+1                                             10d19s20
              call gandc4(ibc(iirc),ibc(iiro),itesta,itestb,nopenr,      10d18s20
     $            notest,norbx,nnot,nab4,bc,ibc)                        11d14s22
             else                                                        10d15s19
              nnot=2                                                     10d15s19
             end if                                                      10d15s19
             if(nnot.ne.0)then                                           10d2s19
              go to 1010                                                10d2s19
             end if                                                      10d2s19
            end do                                                       10d1s19
            if(iflag.eq.2)then                                           10d1s19
             if(interacts.eq.0)then                                      1d23s19
              go to 1010                                                 1d23s19
             else                                                        1d23s19
              itesta=0                                                  10d18s20
              do i=1,ido                                                  10d1s19
               im=i-1                                                     10d1s19
               iorbt(i,1)=idorb(iod+im)                                   10d1s19
               itesta=ibset(itesta,idorb(iod+im))                       10d18s20
              end do                                                      10d1s19
              itestb=0                                                  10d18s20
              do i=0,nsing-1                                            10d18s20
               itestb=ibset(itestb,isorb(ios+i))                        10d18s20
              end do                                                    10d18s20
              itesta=ibset(itesta,norbx)                                10d18s20
              idop=ido+1                                                  10d1s19
              iorbt(idop,1)=64                                            10d1s19
              iargt=idoP+1-mdon                                             10d2s19
              incsf=icsfpd(iargt)+13                                       10d2s19
              ndett=ibc(incsf)                                             10d2s19
              iast=incsf+2                                                 10d2s19
              ivect=iast+nsing*ndett                                      10d2s19
              idnn=max(1,idop-2)                                        10d19s20
              idxx=min(idop+2,mdoo+1)                                   10d19s20
              do isbx=1,nsymb                                             1d23s19
c     when we use resolution of idenity, we need all symmetries.
               if(irxinfo(isbx).gt.0.or.isbx.eq.isymmrci)then              1d23s19
                ifcnl=ifcnpxr(1,1,isbx,idnn)                            10d19s20
                ifcnh=ifcnpxr(2,1,isbx,idxx)                            10d19s20
                call checkc(nfcnt(isbx),iptrt(1,1,isbx),                 1d23s19
     $              ibc(ibasist(isbx)),                                 1d23s19
     $             nek,mdon,idorbr,isorbr,iorbt,idop,isorb(ios),nsing,  1d23s19
     $             iorbt(1,2),iorbt(1,3),nnot,nx,nab,iptrbit(1,1,isbx), 10d18s20
     $              itesta,itestb,norbx,ifcnl,ifcnh)                    10d19s20
                if(nnot.ne.0)go to 1010                                  1d23s19
               end if                                                    1d23s19
              end do                                                     1d23s19
             end if                                                      1d23s19
            end if                                                       10d1s19
            if(ldebug)write(6,*)('skipping '),(idorb(iod+i),i=0,ido-1),  12d12s19
     $          (':'),(isorb(ios+i),i=0,nsing-1)                        12d12s19
            go to 10                                                     10d1s19
           end if                                                        10d1s19
 1010      continue                                                      10d1s19
           nok(isb)=nok(isb)+1                                           7d15s19
           if(iflag.ne.0)then                                            10d2s19
            if(nok(isb).gt.nfcnx)then                                    10d2s19
             write(6,*)('in gencsf2, nok = '),nok(isb),                  10d2s19
     $          (' exceeds nfcnx = '),nfcnx                             10d2s19
             call dws_sync                                               10d2s19
             call dws_finalize                                           10d2s19
             stop                                                        10d2s19
            end if                                                       10d2s19
           end if                                                        10d2s19
           itmp(isto)=ido
           itmp(isto+1)=id
           itmp(isto+2)=is
           isto=isto+3                                                    7d5s19
   10      continue
           iod=iod+ido
          loop=loop+1
          end do
          ios=ios+nsing
         end do
        end if                                                          10d18s20
       end do
       if(ldebug)write(6,*)('total no. of ok combinations '),nok(isb)   12d12s19
       isto=norb+1+ibasis                                               7d15s19
       ntot(isb)=0                                                      7d15s19
       if(ldebug)                                                       12d12s19
     $  write(6,*)('          #  #cs   ci   oi       #csfs     running')12d12s19
       do i=1,nok(isb)                                                  7d15s19
        icsf=itmp(isto)+1-mdon                                           6d5s19
        ntot(isb)=ntot(isb)+ncsf(icsf)                                  7d15s19
        if(ldebug)write(6,*)i,(itmp(isto+j),j=0,2),ncsf(icsf),ntot(isb) 12d12s19
        isto=isto+3                                                     7d15s19
       end do
       if(ldebug)write(6,*)('total no. of csfs: '),ntot(isb)            12d12s19
       do i=1,nok(isb)*3                                                7d15s19
        ip=i+norb
        itmp(i+ibasis)=itmp(ip+ibasis)                                  7d15s19
       end do
       do i=1,mdoo+1                                                    9d20s20
        ifcnpx(1,1,isb,i)=nok(isb)+1                                    10d19s20
        ifcnpx(2,1,isb,i)=0                                             10d19s20
        ifcnpx(1,3,isb,i)=0                                             10d19s20
        ifcnpx(2,3,isb,i)=0                                             10d19s20
        ifcnpx(1,4,isb,i)=0                                             10d19s20
        ifcnpx(2,4,isb,i)=0                                             10d19s20
       end do                                                           9d20s20
        do j=1,4                                                        10d19s20
         nsum(j)=0                                                      10d19s20
        end do                                                          10d19s20
        do i=0,nok(isb)-1                                               7d22s19
         jbasis=ibasis+i*3                                              7d22s19
         nclo=itmp(jbasis+1)                                            7d22s19
         nclop=nclo+1                                                   7d22s19
         if(ifcnpx(1,1,isb,nclop).gt.i+1)then                           10d19s20
          ifcnpx(1,3,isb,nclop)=nsum(1)                                 10d19s20
          ifcnpx(2,3,isb,nclop)=nsum(2)                                 10d19s20
          ifcnpx(1,4,isb,nclop)=nsum(3)                                 10d19s20
          ifcnpx(2,4,isb,nclop)=nsum(4)                                 10d19s20
         end if                                                         10d19s20
         iarg=nclop-mdon                                                10d19s20
         do j=1,4                                                       10d19s20
          nsum(j)=nsum(j)+ncsf2(j,iarg)                                 10d19s20
         end do                                                         10d19s20
         ifcnpx(1,1,isb,nclop)=min(ifcnpx(1,1,isb,nclop),i+1)           10d19s20
         ifcnpx(2,1,isb,nclop)=max(ifcnpx(2,1,isb,nclop),i+1)           10d19s20
         iod=iptr(2,nclop,isb)+nclo*(itmp(jbasis+2)-1)                  7d22s19
         nopen=nec-2*nclo                                               7d22s19
         ios=iptr(4,nclop,isb)+nopen*(itmp(jbasis+3)-1)                 7d22s19
         if(ldebug)write(6,*)i+1,(idorb(iod+j),j=0,nclo-1),(':'),       12d12s19
     $        (isorb(ios+j),j=0,nopen-1),('ifcn'),
     $        (ifcnpx(j,1,isb,nclop),j=1,2),('iod'),iod,nclop,
     $        itmp(jbasis+2)
         nodu=max(nodu,iod+nclo-1)                                      11d24s19
         nosu=max(nosu,ios+nopen-1)                                     11d24s19
        end do                                                          7d22s19
       do i=1,mdoo+1                                                    9d20s20
        ifcnpx(1,2,isb,i)=ifcnpx(1,1,isb,i)                             10d19s20
        ifcnpx(2,2,isb,i)=ifcnpx(2,1,isb,i)                             10d19s20
        if(ifcnpx(1,1,isb,i).eq.nok(isb)+1)ifcnpx(1,1,isb,i)=1          10d19s20
        if(ifcnpx(2,1,isb,i).eq.0)ifcnpx(2,1,isb,i)=nok(isb)            10d19s20
       end do                                                           9d20s20
       if(itestmrci.eq.3.and.iflag.gt.0.and.mynowprog.eq.0)then         10d1s19
        write(7,*)nok(isb)                                              7d22s19
        do i=0,nok(isb)-1                                               7d22s19
         jbasis=ibasis+i*3                                              7d22s19
         nclo=itmp(jbasis+1)                                            7d22s19
         nclop=nclo+1                                                   7d22s19
         iod=iptr(2,nclop,isb)+nclo*(itmp(jbasis+2)-1)                  7d22s19
         nopen=nec-2*nclo                                               7d22s19
         ios=iptr(4,nclop,isb)+nopen*(itmp(jbasis+3)-1)                 7d22s19
         write(7,*)('c: '),(idorb(iod+j),j=0,nclo-1)
         write(7,*)('o: '),(isorb(ios+j),j=0,nopen-1)
         if(ldebug)write(6,*)i+1,(idorb(iod+j),j=0,nclo-1),(':'),       12d12s19
     $        (isorb(ios+j),j=0,nopen-1),('@'),itmp(jbasis+1),
     $        itmp(jbasis+2),itmp(jbasis+3),iod,ios,isb,
     $        iptr(2,nclop,isb),iptr(4,nclop,isb)
        end do                                                          7d22s19
       end if                                                           7d22s19
       ibasis=ibasis+nok(isb)*3                                         7d15s19
       ibasis=ibasis+mod(nok(isb),2)                                    7d15s19
       if(ldebug)then                                                   12d12s19
        write(6,*)('iptr: '),loc(iptr(1,1,isb))
        do j=1,4                                                        12d28s19
         write(6,*)j,(iptr(j,i,isb),i=1,mxc+1)                           9d5s19
        end do
       end if                                                           12d12s19
      end do                                                            7d15s19
      if(mynowprog.eq.0)then                                            8d29s24
       write(6,*)('isorbx = '),isorbx                                   8d29s24
       write(6,*)('idorbx = '),idorbx                                   8d29s24
       write(6,*)('final isto = '),isto
      end if                                                            8d29s24
      return
      end
