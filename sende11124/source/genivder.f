c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine genivder(nder,idarot,nrxinfos,irxinfos,nlzzci,mdoo,    8d15s23
     $     iptr,ioooo,nsdlk,isblk,ionex,nsdlk1,isblk1,idoubo,irefo,noc, 8d17s23
     $     nvirt,nbasdws,nsymb,natom,ngaus,ibdat,isym,iapair,ibstor,    8d17s23
     $     isstor,idorel,ascale,multh,nbasisp,iorb,ih0mo,bc,ibc,isend1, 8d21s23
     $     isend2,isend3,isend4,iptoh,mdon,dorb,sorb,ibasis,icsf,nfcn,  8d23s23
     $     pthrs0,nec,ivintref,nrootci,isymmrci,nct,iptrbit,ixw1,ixw2,  8d23s23
     $     ih0ae,nfcnc,ibasisc,iptrcb,norb,irel,ism,ih0a,idvi,toldv,    9d16s23
     $     maxdiis,maxrest)                                             10d16s24
      implicit real*8 (a-h,o-z)                                         8d15s23
      external second                                                   10d7s24
      character*4 name                                                  8d15s23
      character*50 numbers                                              8d21s23
      include "common.store"                                            8d15s23
      include "common.basis"                                            8d28s23
      dimension ioooo(*),isblk(4,*),ionex(*),isblk1(4,*),irefo(*),      8d17s23
     $     nvirt(*),j4o(512),nbasdws(*),idas(8),npt(3),idatta(7,3),     8d17s23
     $     data(3),ihd(8),idoubo(*),noc(*),nbasisp(*),i4o(512),         8d21s23
     $     i1x(512),iptoh(8,8,8),idumi(512),i4odu(512),isou(8),iorb(*), 8d21s23
     $     iptoh2(8,8,8),multh(8,8),ihdu(8),irxinfos(5,*),nameciu(6),   8d23s23
     $     ibasis(*),nfcn(*),ivintref(*),jvintref(8),nct(*),ih0a(*),    8d28s23
     $     iapair(3,*),idvi(*),idvv(8),nxdata(2,8),idvit(8)             8d30s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      common/drsigncm/drsign                                            8d20s24
      if(mynowprog.eq.0)then
       write(6,*)('generate derivatives of internal vectors used '),     8d15s23
     $        ('in contraction')                                        8d15s23
       write(6,*)('i.e. genivder')
       write(6,*)('toldv = '),toldv                                     7d22s24
       write(6,*)('maxdiis = '),maxdiis,('maxrest = '),maxrest          9d16s23
       call second(time0)                                               10d7s24
      end if
      idwsdeb=0                                                         8d29s23
      do isb=1,nsymb                                                    8d29s23
       idvit(isb)=0                                                     8d30s23
       idvv(isb)=0                                                      8d29s23
      end do                                                            8d29s23
      do i=1,nrxinfos                                                   8d29s23
       isb=irxinfos(1,i)                                                4d25s21
       nrth=irxinfos(2,i)                                               4d25s21
       idvit(isb)=idvit(isb)+nrth                                       8d30s23
       idvv(isb)=idvv(isb)+nrth*nct(isb)                                8d29s23
      end do                                                            8d29s23
      nrb=0                                                             8d29s23
      do isb=1,nsymb                                                    8d29s23
       if(idvit(isb).gt.0)then                                           8d29s23
        nrb=nrb+1                                                       8d29s23
        nxdata(1,nrb)=idvit(isb)                                        8d30s23
        nxdata(2,nrb)=isb                                               8d29s23
       end if                                                           8d29s23
      end do                                                            8d29s23
      call enough('genivder.dvi',bc,ibc)                                8d29s23
      ibcoffo=ibcoff                                                    8d15s23
      ncomp=1                                                           8d17s23
      if(idorel.ne.0)ncomp=2                                            8d17s23
      ih0pp=ibcoff                                                      8d25s23
      jh0mo=ih0mo                                                       8d25s23
      jh0pp=ih0pp                                                       8d25s23
      do isb=1,nsymb                                                    8d25s23
       do ii=0,nbasdws(isb)*nbasdws(isb)-1                              8d25s23
        bc(jh0pp+ii)=bc(jh0mo+ii)                                       8d25s23
       end do                                                           8d25s23
       jh0pp=jh0pp+nbasdws(isb)**2                                      8d25s23
       jh0mo=jh0mo+(ncomp*nbasisp(isb))**2                              8d25s23
      end do                                                            8d25s23
      ibcoff=jh0pp                                                      8d25s23
      do isb=1,nsymb                                                    8d17s23
       ihd(isb)=ibcoff                                                  8d17s23
       ibcoff=ihd(isb)+(nbasisp(isb)*ncomp)**2                          8d17s23
      end do                                                            8d17s23
      nbb=ibcoff-ihd(1)                                                 8d17s23
      call enough('genivder.hd',bc,ibc)                                 8d17s23
      ibcb4=ibcoff                                                      8d21s23
      call paraeri(natom,ngaus,ibdat,idum,ihmat,iorb,noc,ipair,nhcolt,
     $     isym,iapair,ibstor,isstor,multh,iptoh,0,idwsdeb,idorel,
     $     ascale,nbasisp,1,0,idum,0,idum,idum,bc,ibc)                  5d20s24
      icol=ibcoff                                                       8d21s23
      imsg=icol+mynprocg                                                8d21s23
      ibcoff=imsg+mynnode                                               8d21s23
      call parajkfromh(ipair,ngaus,ibdat,nhcolt,ihmat,noc,icol,iorb,    8d21s23
     $     iapair,isstor,ibstor,multh,iptoh,imsg,idumi,idumi,i4o,i1x,     8d21s23
     $     noc,0,idwsdeb,ncomp,nvirt,idumi,0,ih0,idoubo,shift,nbasisp,   8d21s23
     $     idorel,1,idumi,idumi,ibcb4,bc,ibc,isend1,isend2,isend3,      8d21s23
     $     isend4,0,idum)                                               8d21s23
      ih0da=ibcoff                                                      8d21s23
      ibcoff=ih0da+nbb                                                  8d21s23
      call enough('genivder.h0da',bc,ibc)                               8d21s23
      nnp=0                                                             8d21s23
      do isb=1,nsymb                                                    8d21s23
       isou(isb)=nnp                                                    8d21s23
       nhp=nbasisp(isb)*ncomp                                           8d21s23
       nnp=nnp+nhp*nhp                                                  8d21s23
      end do                                                            8d21s23
      ih0dc=ibcoff                                                      8d21s23
      iovrc=ih0dc+nnp*2                                                 8d7s24
      iovrdc=iovrc+nnp*2                                                8d7s24
      ibcoff=iovrdc+nnp*2                                               8d7s24
      call enough('genivder.h0dc',bc,ibc)                               8d21s23
      jdarot=idarot                                                     8d15s23
      do i1=1,nsymb                                                     8d21s23
       do i2=1,nsymb                                                    8d21s23
        do i3=1,nsymb                                                   8d21s23
         iptoh2(i3,i2,i1)=0                                             8d21s23
        end do                                                          8d21s23
       end do                                                           8d21s23
      end do                                                            8d21s23
      ii=0                                                              8d21s23
      do isb=1,nsymb                                                    8d21s23
       if(nbasdws(isb).gt.0)then                                        8d21s23
        do isc=1,nsymb                                                  8d21s23
         if(nbasdws(isc).gt.0)then                                      8d21s23
          isbc=multh(isc,isb)                                            8d21s23
          do isd=1,nsymb                                                 8d21s23
           if(noc(isd).gt.0)then                                         8d21s23
            isa=multh(isbc,isd)                                         8d21s23
            if(nbasdws(isa).gt.0)then                                   8d21s23
             ii=ii+1                                                    8d21s23
             iptoh2(isd,isc,isb)=ii                                     8d21s23
            end if                                                      8d21s23
           end if                                                       8d21s23
          end do                                                        8d21s23
         end if                                                         8d21s23
        end do                                                          8d21s23
       end if                                                           8d21s23
      end do                                                            8d21s23
      do isb=1,nsymb                                                    8d29s23
       idvv(isb)=idvi(isb)                                              8d29s23
      end do                                                            8d29s23
      do ider=1,nder                                                    8d15s23
       ixyz=ibc(jdarot)                                                 7d13s23
       ia=ibc(jdarot+1)                                                 7d13s23
       ipass=ibc(jdarot+2)                                              7d13s23
       ia2=ibc(jdarot+3)                                                7d13s23
       dnuc=0d0                                                         8d28s23
       if(ipass.lt.0)then                                               7d13s23
        if(ixyz.le.3)then                                               8d3s23
         if(mynowprog.eq.0)write(6,*)('for '),char(ichar('w')+ixyz),    9d1s23
     $       (' component of dipole ')                                   7d13s23
         npt(1)=1                                                       3d31s23
         iosym=1                                                        8d18s23
         name='mu  '                                                    3d31s23
         name(3:3)=char(ichar('w')+ixyz)                                8d17s23
         do i=1,7                                                       8d30s22
          idatta(i,1)=0                                                 3d31s23
         end do                                                         8d30s22
         idatta(ixyz,1)=1                                               3d31s23
         data(1)=-1d0                                                   3d31s23
         do ia=1,natom                                                  8d28s23
          dnuc=dnuc+atnum(1,ia)*xcart(ixyz,ia)                          8d28s23
         end do                                                         8d28s23
        else if(ixyz.eq.4)then                                          8d3s23
         npt(1)=3                                                       3d31s23
         iosym=1                                                        8d18s23
         do j=1,3                                                       8d3s23
          do i=1,7                                                      8d3s23
           idatta(i,j)=0                                                8d3s23
          end do                                                        8d3s23
         end do                                                            8d30s22
         idatta(1,1)=2                                                  8d3s23
         data(1)=1d0/sqrt(6d0)                                          8d3s23
         idatta(2,2)=2                                                  8d3s23
         data(2)=1d0/sqrt(6d0)                                          8d3s23
         idatta(3,3)=2                                                  8d3s23
         data(3)=-2d0/sqrt(6d0)                                          8d3s23
         name='Q0  '                                                    8d3s23
         if(mynowprog.eq.0)write(6,*)('for Q0')                         9d1s23
         do ia=1,natom                                                  8d28s23
          dnuc=dnuc+atnum(1,ia)*(2d0*xcart(3,ia)*xcart(3,ia)            8d28s23
     $         -xcart(1,ia)*xcart(1,ia)-xcart(2,ia)*xcart(2,ia))        8d28s23
         end do                                                         8d28s23
         dnuc=dnuc/sqrt(6d0)                                            8d28s23
        else if(ixyz.eq.5)then                                          8d3s23
         npt(1)=2                                                       3d31s23
         iosym=1                                                        8d18s23
         do j=1,2
          do i=1,7                                                          8d30s22
           idatta(i,j)=0                                                   3d31s23
          end do                                                            8d30s22
         end do                                                         8d3s23
         idatta(1,1)=2                                                  8d3s23
         data(1)=-0.5d0
         idatta(2,2)=2                                                  8d3s23
         data(2)=0.5d0                                                     3d31s23
         name='ReQ2'                                                    8d3s23
         if(mynowprog.eq.0)write(6,*)('for ReQ2')                       9d1s23
         do ia=1,natom                                                  8d28s23
          dnuc=dnuc+atnum(1,ia)*(xcart(1,ia)*xcart(1,ia)                8d28s23
     $           -xcart(2,ia)*xcart(2,ia))                              8d28s23
         end do                                                         8d28s23
         dnuc=dnuc*0.5d0                                                8d28s23
        else if(ixyz.eq.6)then                                          8d3s23
         name='ImQ2'                                                    8d3s23
         npt(1)=1                                                         3d31s23
         iosym=1                                                        8d18s23
         do i=1,7                                                          8d30s22
          idatta(i,1)=0                                                 8d3s23
         end do                                                            8d30s22
         idatta(1,1)=1                                                  8d3s23
         idatta(2,1)=1                                                  8d3s23
         data(1)=-1d0
         if(mynowprog.eq.0)write(6,*)('for ImQ2')                       9d1s23
         do ia=1,natom                                                  8d28s23
          dnuc=dnuc+atnum(1,ia)*xcart(1,ia)*xcart(2,ia)                 8d28s23
         end do                                                         8d28s23
        else if(ixyz.eq.7)then                                          8d3s23
         name='ReQ1'                                                    8d3s23
         npt(1)=1                                                         3d31s23
         iosym=1                                                        8d18s23
         do i=1,7                                                          8d30s22
          idatta(i,1)=0                                                 8d3s23
         end do                                                            8d30s22
         idatta(1,1)=1                                                  8d3s23
         idatta(3,1)=1                                                  8d3s23
         data(1)=1d0
         if(mynowprog.eq.0)write(6,*)('for ReQ1')                       9d1s23
         do ia=1,natom                                                  8d28s23
          dnuc=dnuc+atnum(1,ia)*xcart(1,ia)*xcart(3,ia)                 8d28s23
         end do                                                         8d28s23
         dnuc=-dnuc                                                     8d28s23
        else                                                            8d3s23
         if(mynowprog.eq.0)write(6,*)('for ImQ1')                       9d1s23
         name='ImQ1'                                                    8d3s23
         npt(1)=1                                                         3d31s23
         iosym=1                                                        8d18s23
         do i=1,7                                                          8d30s22
          idatta(i,1)=0                                                 8d3s23
         end do                                                            8d30s22
         idatta(2,1)=1                                                  8d3s23
         idatta(3,1)=1                                                  8d3s23
         data(1)=1d0
         do ia=1,natom                                                  8d28s23
          dnuc=dnuc+atnum(1,ia)*xcart(2,ia)*xcart(3,ia)                 8d28s23
         end do                                                         8d28s23
         dnuc=-dnuc                                                     8d28s23
        end if                                                          8d3s23
        do iz=ihd(1),ihd(1)+nbb-1                                       8d25s23
         bc(iz)=0d0                                                     8d25s23
        end do                                                          8d25s23
        call parap(natom,ngaus,ibdat,ihd,isym,iapair,ibstor,isstor,     8d17s23
     $       idum,0,idorel,ascale,1,npt,data,idatta,iosym,1,multh,nbb,  8d17s23
     $       nbasisp,nbasdws,iorb,name,1,bc,ibc)                        8d17s23
        jhd=ihd(1)                                                      8d25s23
        dcore=0d0                                                       8d28s23
        do isb=1,nsymb
         if(nbasdws(isb).gt.0)then
          if(idwsdeb.gt.10)then                                         8d29s23
           write(6,*)('for symmetry block '),isb
           call prntm2(bc(ihd(isb)),nbasdws(isb),nbasdws(isb),
     $         nbasdws(isb))
          end if                                                        8d29s23
          do ii=0,nbasdws(isb)*nbasdws(isb)-1                           8d25s23
           bc(jhd+ii)=bc(ihd(isb)+ii)                                   8d25s23
          end do                                                        8d25s23
          jhd=jhd+nbasdws(isb)*nbasdws(isb)                             8d25s23
         end if
        end do
       else                                                             7d13s23
        dcore=0d0                                                       8d28s23
        idersign=1                                                      8d21s23
        if(ia2.eq.0)then                                                7d13s23
         if(mynowprog.eq.0)write(6,*)('for d/d'),char(ichar('w')+ixyz), 9d1s23
     $       ia                                                         9d1s23
        else                                                            7d13s23
         if(ipass.eq.1)then                                             7d13s23
          if(mynowprog.eq.0)write(6,*)('for d/d'),char(ichar('w')+ixyz),9d1s23
     $        ia,('+'),ia2                                              9d1s23
         else                                                           7d13s23
          if(mynowprog.eq.0)write(6,*)('for d/d'),char(ichar('w')+ixyz),9d1s23
     $         ia,('-'),ia2                                             9d1s23
          idersign=2                                                    8d21s23
         end if                                                         7d13s23
        end if                                                          7d13s23
        pnuc=0d0                                                        8d28s23
        do ixpass=1,ipass                                               8d28s23
         sig=1d0                                                        8d28s23
         if(ixpass.eq.1)then                                            8d28s23
          ja=ia                                                         8d28s23
         else                                                           8d28s23
          ja=iapair(1,ia)                                               8d28s23
          if(idersign.eq.2)sig=-1d0                                     8d28s23
         end if                                                         8d28s23
         do iai=1,natom                                                 8d28s23
          if(iai.ne.ja)then                                             8d28s23
           dist=0d0                                                     8d28s23
           do jxyz=1,3                                                  8d28s23
            dist=dist+(xcart(jxyz,iai)-xcart(jxyz,ja))**2               8d28s23
           end do                                                       8d28s23
           dist=1d0/dist                                                8d28s23
           dists=sqrt(dist)                                             8d28s23
           vij=atnum(1,iai)*atnum(1,ja)*dists                           8d28s23
           pnuc=pnuc+vij                                                8d28s23
           dvij=(xcart(ixyz,iai)-xcart(ixyz,ja))*vij*dist               8d28s23
           dnuc=dnuc-drsign*dvij*sig                                    8d20s24
          end if                                                        8d28s23
         end do                                                         8d28s23
        end do                                                          8d28s23
        if(idwsdeb.gt.10)write(6,*)('dnuc: '),dnuc                      8d29s23
        call parah0grad(natom,ngaus,ibdat,nbasis,bc(ih0dc),bc(iovrdc),  8d21s23
     $         bc(iovrdc),isym,iapair,ibstor,isstor,isou,nnp,idwsdeb,   8d29s23
     $         idorel,ascale,multh,ixyz,ia,1,idersign,nbasisp,.false.,  8d21s23
     $          dum,idum,idum,bc,ibc)                                   8d21s23
        jhd=ihd(1)                                                      8d25s23
        do isb=1,nsymb                                                  8d21s23
         nhp=nbasisp(isb)*ncomp                                         8d21s23
         if(nhp.gt.0)then                                               8d21s23
          iad=ih0dc+isou(isb)                                           8d21s23
          if(idwsdeb.gt.10)then                                         8d29s23
           write(6,*)('h0 after parah0grad for symmetry '),isb           8d21s23
           call mpprnt2(bc(iad),nhp)                                     8d21s23
           write(6,*)('orbitals ')
           call prntm2(bc(iorb(isb)),nhp,nbasdws(isb),nhp)
          end if                                                        8d29s23
          call square(bc(iad),nhp)                                      8d21s23
          nr=nhp                                                        8d21s23
          do ips=1,2                                                    8d21s23
           call dgemm('n','n',nr,nbasdws(isb),nhp,1d0,bc(iad),nr,       8d28s23
     $          bc(iorb(isb)),nhp,0d0,bc(iovrdc),nr,'genivder.h0')      8d21s23
           do i=0,nr-1                                                  8d21s23
            do j=0,nbasdws(isb)-1                                       8d21s23
             ji=iad+j+nbasdws(isb)*i                                    8d21s23
             ij=iovrdc+i+nr*j                                           8d21s23
             bc(ji)=bc(ij)                                              8d21s23
            end do                                                      8d21s23
           end do                                                       8d21s23
           nr=nbasdws(isb)                                              8d21s23
          end do                                                        8d21s23
          do i=0,nr*nr-1                                                8d21s23
           bc(jhd+i)=bc(iad+i)                                          8d25s23
          end do                                                        8d21s23
          if(idwsdeb.gt.10)then                                         8d29s23
           write(6,*)('in mo basis '),jhd-ihd(1)+1                                    8d21s23
           call prntm2(bc(jhd),nbasdws(isb),nbasdws(isb),                8d25s23
     $         nbasdws(isb))                                            8d21s23
          end if                                                        8d29s23
          jhd=jhd+nr*nr                                                 8d25s23
         end if                                                         8d21s23
        end do                                                          8d21s23
       end if                                                           7d13s23
       jdarot=jdarot+4                                                  7d13s23
       do isb=1,nsymb                                                   8d15s23
        if(nbasdws(isb).gt.0)then                                       8d15s23
         if(idwsdeb.gt.10)then                                          8d29s23
          write(6,*)('darot for sym '),isb                               8d15s23
          call prntm2(bc(jdarot),nbasdws(isb),nbasdws(isb),nbasdws(isb)) 8d15s23
         end if                                                         8d29s23
         idas(isb)=jdarot                                               8d15s23
         jdarot=jdarot+nbasdws(isb)*nbasdws(isb)                         8d15s23
         if(ipass.ge.0)then                                             8d21s23
          if(idwsdeb.gt.10)then                                         8d29s23
           write(6,*)('transder for sym '),isb                           8d21s23
           call prntm2(bc(jdarot),nbasdws(isb),nbasdws(isb),            8d29s23
     $          nbasdws(isb))                                           8d29s23
          end if                                                        8d29s23
          itmp=ibcoff                                                   8d21s23
          ibcoff=itmp+nbasdws(isb)*nbasdws(isb)                         8d21s23
          call enough('genivder.tmp',bc,ibc)                            8d21s23
          do i=0,nbasdws(isb)*nbasdws(isb)-1                            8d21s23
           bc(itmp+i)=bc(jdarot+i)+bc(idas(isb)+i)                      8d21s23
          end do                                                        8d21s23
          if(idwsdeb.gt.10)then                                         8d29s23
           write(6,*)('sum of darot and transder '),itmp
           call prntm2(bc(itmp),nbasdws(isb),nbasdws(isb),nbasdws(isb))  8d21s23
           if(isb.eq.2)itmpw=itmp
          end if                                                        8d29s23
          idas(isb)=itmp                                                8d21s23
          jdarot=jdarot+nbasdws(isb)*nbasdws(isb)                       8d21s23
         end if                                                         8d21s23
        end if                                                          8d15s23
       end do                                                           8d15s23
       call parajkfromhd0(noc,multh,idumi,idumi,i4o,i1x,itprt,          8d21s23
     $            idwsdeb,nvirt,1,idas,i4odu,ionexdu,nbasdws,           8d29s23
     $            isblkder,isblk1,nsblkder,nsdlk1,isblkkder,            8d21s23
     $            -10,jmatd,kmatd,i3x,i4od2b,ionexd2,isblkxder1,        5d16s22
     $            0,isblkder1,0,idorel,igoal,bc,ibc)                    11d10s22
       if(ipass.gt.0)then                                               8d21s23
        call paraerid(natom,ngaus,ibdat,nbasis,ihmat,iorb,noc,ipair,    8d21s23
     $      nhcolt,isym,iapair,ibstor,isstor,multh,iptoh2,0,idwsdeb,    8d29s23
     $      idorel,ascale,ia,ixyz,1,idersign,1,nbasisp,bc,ibc)          8d21s23
       idumix=-1
        call parajkfromhd1(ipair,ngaus,ibdat,nhcolt,ihmat,noc,icol,     8d21s23
     $          iorb,iapair,isstor,ibstor,multh,iptoh2,imsg,idumi,idumi,8d21s23
     $         idumi,idumi,noc,0,idwsdeb,ncomp,nvirt,1,idum,            8d29s23
     $            i4odu,idumix,idumi,idumi,idumi,idumi,der2e,nbasdws,   8d17s24
     $            isblk,isblk1,nsdlk,nsdlk1,                            8d28s23
     $            isblkxder1,0,isblkder1,0,isblkkder,                   8d21s23
     $            0,nbasisp,1,ihd,idorel,bc,ibc)                        8d21s23
       else
        call sym4o(i4odu,noc,isblk,isblk1,nsdlk,nsdlk1,bc,ibc)           8d21s23
       end if                                                           8d21s23
       call derh01b(bc(ih0pp),bc(ihd(1)),bc(ih0da),idas,nsymb,          8d25s23
     $          nbasdws,idwsdeb,multh,1,bc,ibc)                         8d29s23
       call foldr(ih0da,i4odu,noc,irefo,idoubo,nbasdws,nsymb,           8d21s23
     $      isblk,nsdlk,ihdu,bc,ibc,dcore,idwsdeb)                      8d29s23
       ikeep=ibcoff                                                     8d23s23
       ibcoff=ikeep+nfcnc                                               8d23s23
       call enough('genivder.keep',bc,ibc)                              8d23s23
       do isb=1,nsymb                                                   8d23s23
        jvintref(isb)=ivintref(isb)                                     8d23s23
       end do                                                           8d23s23
       dshift=dnuc+dcore                                                8d28s23
       if(idwsdeb.gt.10)write(6,*)('dshift: '),dshift,dnuc,dcore        8d29s23
       do i=1,nrb                                                       8d29s23
        nrth=nxdata(1,i)                                                8d29s23
        isb=nxdata(2,i)                                                 8d29s23
        if(idwsdeb.gt.10)then                                           8d29s23
         write(6,*)('for i = '),i
         write(6,*)('isb = '),isb
         write(6,*)('nfcn(isb) = '),nfcn(isb)
        end if                                                          8d29s23
        jptr=iptr+(mdoo+1)*2*(isb-1)                                     4d25s21
        idvec=idvv(isb)+nrth                                            8d29s23
        call intcsfder(dorb,sorb,mdon,mdoo,ibc(ibasis(isb)),ibc(jptr),  8d23s23
     $       icsf,nfcn(isb),ih0ae,ioooo,ihdu,i4odu,pthrs0,nec,multh,isb, 8d25s23
     $       nrth,mynowprog.eq.0,bc(jvintref(isb)),nct(isb),irefo,bc,   8d23s23
     $       ibc,ibc(iptrbit),ixw1,ixw2,ibc(ikeep),dorbf,sorbf,nfcnc,   8d23s23
     $       ibasisc,iptrcb,norb,irel,ism,ih0a,dshift,bc(idvv(isb)),    8d29s23
     $       bc(idvec),toldv,maxdiis,maxrest)                           10d16s24
        idvv(isb)=idvv(isb)+nrth*(1+nct(isb))                           8d29s23
       end do                                                           8d21s23
      end do                                                            8d15s23
      ibcoff=ibcoffo                                                    8d15s23
      if(mynowprog.eq.0)then                                            10d7s24
       call second(timen)                                               10d7s24
       telap=timen-time0                                                10d7s24
       write(6,*)('total cpu time on this proc for genivder is '),telap 10d7s24
      end if                                                            10d7s24
      return                                                            8d15s23
      end                                                               8d15s23
