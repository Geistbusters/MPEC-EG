c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine grado(natom,ngaus,ibdat,nbasis,isym,iapair,
     $     ibstor,isstor,iso,nbb,idwsdeb,idorel,ascale,multh,ipropsym,
     $     ivecs,icanon,ieigs,iorb,ih0,morb,idbly,ixinv,iooood,ionexd,  3d17s21
     $     jmatd,kmatd,ioooo,ionex,jmats,kmats,nbasdwsc,nvirtc,         3d28s16
     $     ih0mo,ipropmat,i3x,lprint,iact,nbasisp,nbasisc,ivecso,bc,ibc)11d10s22
      implicit real*8 (a-h,o-z)
      external second
c
c     compute gradients of orbitals and double gradients of orbitals    6d17s16
c     i.e. no mixed 2nd derivatives.                                    6d17s16
c
      integer*8 nn8                                                     2d19s10
      include "common.hf"
      include "common.store"
      include "common.spher"
      integer*8 ibstor,isstor                                           5d6s10
      include "common.basis"
      logical lprint                                                    3d2s17
      character*1 cphas
      character*4 olab(3)                                               4d15s22
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension isym(3,1),iapair(3,1),ibstor(1),isstor(1),iorb(8),
     $     iso(8),multh(8,8),ipropsym(6),isou(8),ivecs(1),ieigs(8),
     $     morb(8),icanon(8),noc(8),iptoh(8,8,8),ixinv(8),nvirtc(8),    2d23s16
     $     nbasdwsc(8),iamat(8),itrans(8),i4od(idbk),ionexd(idbk),      3d22s16
     $     igmat(8,8),isoud(8),kmatd(idbk),jmatd(idbk),i4od2(idbk),     7d27s16
     $     ionexd2(idbk),i3x(idbk),itran2(8),isblkder(4,idbk),          8d3s16
     $     isblkxder(4,idbk),isblkkder(4,idbk),isblkder1(4,idbk),       8d22s16
     $     isblkxder1(4,idbk),i4od2b(idbk),ihessd(8,8),iact(*),         3d17s21
     $     nbasisp(*),nbasisc(*),idbly(*),ivecso(*),isoub(8),           6d13s22
     $     ipxder(4,8,8,8)                                              11d28s22
      common/timerocm/tovr,telapo(15)                                   4d26s18
      common/drsigncm/drsign                                            8d20s24
      data olab/'d/dx','d/dy','d/dz'/                                   4d15s22
      do it=1,5                                                         3d2s17
       telapo(it)=0d0                                                   3d2s17
      end do                                                            3d2s17
      ncomp=1
      idwsdeb=0
      if(idorel.ne.0)ncomp=2
      if(lprint)write(6,*)('in grado ...'),ibcoff
      igoal=1
      do isb=1,nsymb
       noc(isb)=idbly(isb)+iact(isb)                                    3d17s21
      end do
      itprt=0                                                            3d3s17
      ipropmat=ibcoff                                                   3d29s16
      ispropmat=ipropmat+natom*4                                        6d20s16
      ieder=ispropmat+natom*3                                           4d7s16
      ibcoff=ieder+natom*3                                              4d7s16
      ibodc=ibcoff                                                      12d29s16
      ibcoff=ibodc+natom                                                12d29s16
      call enough('grado.  1',bc,ibc)
      call second(time1)                                                3d2s17
      call buildder3s(ionex,ider3,bc(ih0mo),idbly,nvirtc,nbasdwsc,      3d17s21
     $     idwsdeb,bc,ibc)                                              11d14s22
      call second(time2)                                                3d2s17
      telap=time2-time1-tovr                                            3d2s17
      telapo(4)=telap                                                   3d2s17
      do i=0,natom-1                                                    12d29s16
       bc(ibodc+i)=0d0                                                  12d29s16
      end do                                                            12d29s16
      n2ndsz=0                                                          6d21s16
      do isb=1,nsymb                                                    6d21s16
       n2ndsz=n2ndsz+nbasdwsc(isb)*nbasdwsc(isb)                        6d21s16
       n2ndsz=n2ndsz+noc(isb)*nvirtc(isb)                               12d2s16
      end do                                                            6d21s16
      do ipsym=nsymb,1,-1
       if(lprint)write(6,*)('for perturbations of symmetry '),ipsym
       nhit=0
       do ixyz=1,3                                                      3d28s16
        ixyzm=ixyz-1                                                    3d29s16
        ipuse=ipropsym(ixyz)                                            3d28s16
        do ia=1,natom                                                   3d28s16
         if(iapair(1,ia).eq.0)then
          if(ipuse.eq.ipsym)nhit=nhit+1
          iad=ispropmat+ixyzm+3*(ia-1)                                  3d29s16
          ibc(iad)=ipuse                                                3d29s16
         else if(iapair(1,ia).gt.0)then
          if(ipuse.eq.ipsym)nhit=nhit+1
          iad=ispropmat+ixyzm+3*(ia-1)                                  3d29s16
          ibc(iad)=ipuse                                                3d29s16
          itry=multh(ipuse,iapair(2,ia))
          if(itry.eq.ipsym)nhit=nhit+1
          iad=ispropmat+ixyzm+3*(iapair(1,ia)-1)                        3d29s16
          ibc(iad)=itry                                                 3d29s16
         end if
        end do                                                          3d28s16
       end do                                                           3d28s16
       nhit0=nhit                                                       6d17s16
       if(lprint)                                                       3d2s17
     $ write(6,*)('number of desired operators having this symmetry is')3d28s16
     $      ,nhit                                                       3d28s16
       if(nhit.gt.0)then                                                3d28s16
        nsz=0                                                           3d29s16
        do isb=1,nsymb                                                  3d29s16
         isk=multh(isb,ipsym)                                           3d29s16
         nsz=nsz+nbasdwsc(isb)*nbasdwsc(isk)                            3d29s16
        end do                                                          3d29s16
        do ixyzm=0,2                                                    3d29s16
         do ia=0,natom-1                                                3d29s16
          iad=ispropmat+ixyzm+3*ia                                      3d29s16
          if(ibc(iad).eq.ipsym)then                                     3d29s16
           iad=ipropmat+ixyzm+3*ia                                      3d29s16
           ibc(iad)=ibcoff                                              3d29s16
           ibcoff=ibcoff+nsz+n2ndsz                                     6d21s16
          end if                                                        3d29s16
         end do                                                         3d29s16
        end do                                                          3d29s16
        nsblkkder=-1                                                    8d4s16
        call gendertype(1,multh,noc,isblkder1,isblkxder1,nsblkder1,     8d4s16
     $        nsblkxder1,isblkkder,nsblkkder,ipxder)                    11d28s22
        nsblkkder=0                                                     8d4s16
        call gendertype(ipsym,multh,noc,isblkder,isblkxder,nsblkder,    8d3s16
     $        nsblkxder,isblkkder,nsblkkder,ipxder)                     6d13s22
        ibcoffo=ibcoff                                                    2d19s10
        icol=ibcoff
        ibcoff=icol+mynprocg
        imsg=ibcoff                                                       6d3s10
        ibcoff=imsg+mynnode                                               6d3s10
        ih0d=ibcoff                                                     3d29s16
        do isb=1,nsymb
         nh=nbasdwsc(isb)
         jvecs=ivecso(isb)                                              4d4s22
         if(idwsdeb.gt.10)then
          write(6,*)('eigenvectors of overlap matrix for symmetry '),
     $         isb,jvecs
          call prntm2(bc(jvecs),nh,nh,nh)
         end if
        end do
        call second(time1)                                              3d2s17
         call buildhess(ioooo,ionex,jmats,kmats,bc(ih0mo),noc,igmat,    3d28s16
     $       nvirtc,ipsym,multh,idwsdeb,bc,ibc)                         11d14s22
         ibcgmat=ibcoff                                                 4d1s22
         call second(time2)                                             3d2s17
         telap=time2-time1-tovr                                         3d2s17
         telapo(4)=telapo(4)+telap                                      3d2s17
        do ixyz=1,3
         do ia=1,natom
          if(iapair(1,ia).ge.0)then
           ipuse=ipropsym(ixyz)                                           2d3s16
           if(iapair(1,ia).eq.0)then                                      2d3s16
            npass=1                                                       2d3s16
           else                                                           2d3s16
            npass=2
           end if                                                         2d3s16
           bodcfact=1d0                                                 5d4s22
           do ipass=1,npass
            if(ipuse.eq.ipsym)then                                       3d28s16
             if(lprint)then                                             3d2s17
              write(6,*)('for nucleus no. '),ia                          4d19s23
              write(6,*)('atno, carts '),atnum(1,ia),(xcart(j,ia),j=1,3)
              if(iapair(1,ia).eq.0)then                                      2d3s16
               write(6,*)('symmetrically unique ')
              else                                                           2d3s16
               write(6,*)('symmetrically related to nucleus '),          4d19s23
     $              iapair(1,ia)                                        4d19s23
               write(6,*)('atno, carts '),atnum(1,ia),
     $            (xcart(j,iapair(1,ia)),j=1,3)
              end if
             end if                                                     3d2s17
             ixyzm=ixyz-1
             idersign=1                                                 4d21s22
             if(npass.eq.1)then
              cphas=' '                                                    2d3s16
              iad=ipropmat+ixyzm+3*(ia-1)                               3d29s16
              npropmat=ibc(iad)                                         3d29s16
             else if(ipass.eq.1)then                                       2d3s16
              cphas='+'                                                    2d3s16
              iad=ipropmat+ixyzm+3*(ia-1)                               3d29s16
              npropmat=ibc(iad)                                         3d29s16
              bodcfact=0.5d0                                            5d4s22
             else
              cphas='-'                                                    2d3s16
              idersign=2                                                4d21s22
              iad=ipropmat+ixyzm+3*(iapair(1,ia)-1)                     3d29s16
              npropmat=ibc(iad)                                         3d29s16
              bodcfact=0.5d0                                            5d4s22
             end if                                                        2d3s16
             if(lprint)then                                             4d28s22
              if(iapair(1,ia).eq.0)then                                 4d28s22
               write(6,1)olab(ixyz),ia,ipuse                            4d28s22
              else                                                      4d28s22
               write(6,11)olab(ixyz),ia,cphas,iapair(1,ia),ipuse        4d28s22
   11          format('for operator ',a4,i2,a1,i2,' with symmetry ',i1) 4d28s22
              end if                                                    4d28s22
             end if                                                     4d28s22
    1        format('for operator ',a4,i2,' with symmetry ',i1)         4d28s22
             if(ipuse.eq.1)then
              nnp=0                                                     4d8s22
              nn=0                                                      4d8s22
              do isb=1,nsymb
               nh=nbasdwsc(isb)
               nhp=nbasisp(isb)*ncomp                                   4d8s22
               isou(isb)=nnp                                            4d8s22
               isoub(isb)=nn                                            4d8s22
               isoud(isb)=nnp                                           4d8s22
               nn=nn+nh*nh
               nnp=nnp+nhp*nhp                                          4d8s22
              end do
              nn1=nnp                                                   4d8s22
             else
              nn=0
              nnp=0                                                     4d8s22
              nn1=0                                                     6d20s16
              do isb=1,nsymb
               isoud(isb)=nn1                                           6d20s16
               nn1=nn1+(nbasisp(isb)*ncomp)**2                          4d8s22
               isk=multh(isb,ipuse)                                        2d3s16
               if(isk.ge.isb)then
                isou(isb)=nnp                                           4d8s22
                isoub(isb)=nn                                           4d8s22
                nn=nn+nbasdwsc(isb)*nbasdwsc(isk)
                nnp=nnp+nbasisp(isb)*nbasisp(isk)*ncomp*ncomp           4d8s22
                isou(isk)=nnp                                           4d8s22
                isoub(isk)=nn                                           4d8s22
                nn=nn+nbasdwsc(isb)*nbasdwsc(isk)                       7d21s16
                nnp=nnp+nbasisp(isb)*nbasisp(isk)*ncomp*ncomp           4d8s22
               end if
              end do
             end if
             ih0d=ibcoff
             iovr=ih0d+nnp*2                                            4d8s22
             ibcoff=iovr+nnp*2                                          4d8s22
             iovrdk=ibcoff                                              4d4s16
             ibcoff=iovrdk+nnp*2                                        4d8s22
             ih0d2=ibcoff                                               6d20s16
             iovrd2=ih0d2+nn1*2                                         6d20s16
             iovrd2k=iovrd2+nn1*2                                       6d20s16
             ibcoff=iovrd2k+nn1*2                                       6d20s16
             call enough('grado.  2',bc,ibc)
             do iz=ih0d,ibcoff-1                                        4d14s22
              bc(iz)=0d0                                                4d14s22
             end do                                                     4d14s22
             derpotn=0d0
             do ixpass=1,npass                                          4d22s16
              sig=1d0                                                   4d22s16
              if(ixpass.eq.1)then                                       4d22s16
               ja=ia                                                    4d22s16
              else                                                      4d22s16
               ja=iapair(1,ia)                                          4d22s16
               if(idersign.eq.2)sig=-1d0                                4d21s22
              end if                                                    4d22s16
              do iai=1,natom                                            4d22s16
               if(iai.ne.ja)then                                        4d22s16
                dist=0d0                                                4d22s16
                do jxyz=1,3                                             4d22s16
                 dist=dist+(xcart(jxyz,iai)-xcart(jxyz,ja))**2          4d18s16
                end do                                                  4d22s16
                dist=1d0/dist                                           4d22s16
                dists=sqrt(dist)                                        4d22s16
                vij=atnum(1,iai)*atnum(1,ja)*dists                      4d22s16
                dvij=(xcart(ixyz,iai)-xcart(ixyz,ja))*vij*dist          4d22s16
                derpotn=derpotn+drsign*dvij*sig                         8d20s24
               end if                                                   4d22s16
              end do                                                    4d22s16
             end do                                                     4d22s16
             if(lprint)                                                 3d2s17
     $       write(6,*)('derivative of nuclear-nuclear repulsion '),       3d21s16
     $         derpotn                                                  3d21s16
             call second(time1)                                         3d2s17
             call parah0grad(natom,ngaus,ibdat,nbasis,bc(ih0d),bc(iovr),
     $         bc(iovrdk),isym,iapair,ibstor,isstor,isou,nnp,           8d4s22
     $            idwsdeb,                                              5d12s16
     $            idorel,
     $            ascale,multh,ixyz,ia,ipuse,idersign,nbasisp,.false.,  6d2s22
     $            dum,idum,idum,bc,ibc)                                 11d28s22
             call parah0d2(natom,ngaus,ibdat,nbasis,bc(ih0d2),
     $           bc(iovrd2),bc(iovrd2k),isym,iapair,ibstor,isstor,isoud,6d20s16
     $                nn1,idwsdeb,idorel,ascale,multh,ixyz,ia,ipuse,    3d7s23
     $            idersign,nbasisp,bc,ibc)                              11d10s22
             call second(time2)                                         3d2s17
             telap=time2-time1-tovr                                     3d2s17
             telapo(1)=telapo(1)+telap                                  3d2s17
             call derofxor(bc(iovr),isou,ivecs,ipuse,idorel,icanon,
     $       ieigs,ixor,multh,idwsdeb,iorb,bc(npropmat),bc(iovrdk),         6d20s16
     $            bc(iovrd2),bc(iovrd2k),ivecso,nbasisp,isoub,bc,ibc)   11d14s22
             call transder(noc,morb,ixor,ixinv,nbasdwsc,nvirtc,bc(ih0),    3d8s16
     $            iorb,itrans,ipuse,multh,idwsdeb,bc(npropmat),         10d27s16
     $            itran2,isoub,bc,ibc)                                  11d10s22
             call derh0(bc(ih0),bc(ih0d),ivecs,ixor,morb,idorel,
     $            ipuse,multh,ih0der,noc,iorb,idwsdeb,dere1,der2e1,     7d21s16
     $            bc(ih0d2),ih0der2,itrans,itran2,nbasisp,isou,bc,ibc)  11d14s22
             do i1=1,nsymb                                                     7d22s14
              do i2=1,nsymb                                                    7d22s14
               do i3=1,nsymb                                                   7d22s14
                iptoh(i3,i2,i1)=0                                              7d22s14
               end do                                                          7d22s14
              end do                                                           7d22s14
             end do                                                            7d22s14
             ii=0
             do isb=1,nsymb                                             9d26s16
              if(nbasdwsc(isb).gt.0)then                                9d26s16
               do isc=1,nsymb                                           9d27s16
                if(nbasdwsc(isc).gt.0)then                              10d3s16
                 isbc=multh(isc,isb)                                    9d27s16
                 do isd=1,nsymb                                         9d26s16
                  if(noc(isd).gt.0)then                                 9d26s16
                   isa=multh(ipuse,multh(isbc,isd))                     11d29s22
                   if(nbasdwsc(isa).gt.0)then                           9d26s16
                     ii=ii+1                                             9d26s16
                     iptoh(isd,isc,isb)=ii                               9d26s16
                   end if                                               9d26s16
                  end if                                                9d26s16
                 end do                                                 9d26s16
                end if                                                  9d26s16
               end do                                                   9d26s16
              end if                                                    9d26s16
             end do                                                     9d26s16
             ndum=-1                                                    8d4s16
             call second(time3)                                         3d2s17
             telap=time3-time2-tovr                                     3d2s17
             telapo(4)=telapo(4)+telap                                  3d2s17
             call parajkfromhd0(noc,multh,jmats,kmats,ioooo,ionex,itprt, 3d3s17
     $            idwsdeb,nvirtc,1,itran2,i4od2,ionexd2,nbasdwsc,       9d16s16
     $            isblkder1,isblkxder1,nsblkder1,nsblkxder1,            8d4s16
     $            isblkkder,ndum,idum,idum,idum,idum,idum,              11d28s22
     $            idum,0,idum,0,idorel,igoal,bc,ibc)                    11d28s22
             ilookat=81*9+ionexd2(1)
             call parajkfromhd0(noc,multh,jmats,kmats,ioooo,ionex,itprt, 3d3s17
     $            idwsdeb,nvirtc,ipuse,itrans,i4od,ionexd,nbasdwsc,     8d4s16
     $            isblkder,isblkxder,nsblkder,nsblkxder,isblkkder,      8d4s16
     $            nsblkkder,jmatd,kmatd,i3x,i4od2b,ionexd2,isblkxder1,  9d16s16
     $            nsblkxder1,isblkder1,nsblkder1,idorel,igoal,bc,ibc)   11d10s22
             call second(time2)                                         3d2s17
             telap=time2-time3-tovr                                     3d2s17
             telapo(3)=telapo(3)+telap                                  3d2s17
             ibcsav=ibcoff
             call paraerid(natom,ngaus,ibdat,nbasis,ihmat,iorb,noc,     3d28s16
     $            ipair,nhcolt,isym,iapair,ibstor,isstor,multh,iptoh,   4d21s22
     $            itprt,idwsdeb,idorel,ascale,ia,ixyz,ipuse,idersign,   4d21s22
     $            1,nbasisp,bc,ibc)                                     11d10s22
             call second(time1)                                         3d2s17
             telap=time1-time2-tovr                                     3d2s17
             telapo(2)=telapo(2)+telap                                  3d2s17
             call parajkfromhd1(ipair,ngaus,ibdat,nhcolt,ihmat,noc,icol,     2d22s16
     $           iorb,iapair,isstor,ibstor,multh,iptoh,imsg,jmats,kmats,      2d23s16
     $         ioooo,ionex,noc,itprt,idwsdeb,ncomp,nvirtc,ipuse,itrans,  3d3s17
     $            i4od,ionexd,kmatd,jmatd,i4od2b,ionexd2,der2e,nbasdwsc, 8d3s16
     $            isblkder,isblkxder,nsblkder,nsblkxder,                9d16s16
     $            isblkxder1,nsblkxder1,isblkder1,nsblkder1,isblkkder,  9d16s16
     $            nsblkkder,nbasisp,1,igoal,idorel,bc,ibc)              11d10s22
             ibcoff=ibcsav                                              9d20s16
             do i1=1,nsymb                                              9d26s16
              do i2=1,nsymb                                             9d26s16
               do i3=1,nsymb                                            9d26s16
                iptoh(i3,i2,i1)=0                                       9d26s16
               end do                                                   9d26s16
              end do                                                    9d26s16
             end do                                                     9d26s16
             ii=0
             do isb=1,nsymb                                             9d26s16
              if(nbasdwsc(isb).gt.0)then                                9d26s16
               do isc=1,nsymb                                           9d27s16
                if(noc(isc).gt.0)then                                   9d26s16
                 isbc=multh(isc,isb)                                    9d27s16
                 do isd=1,nsymb                                         9d26s16
                  if(noc(isd).gt.0)then                                 9d26s16
                   isad=multh(isbc,isd)                                 9d27s16
                   isa=multh(isad,ipuse)                                9d27s16
                   if(nbasdwsc(isa).gt.0)then                           9d26s16
                   if(isa.le.isb)then
                    ii=ii+1                                             9d26s16
                    iptoh(isd,isc,isb)=ii                               9d26s16
                   end if
                   end if                                               9d26s16
                  end if                                                9d26s16
                 end do                                                 9d26s16
                end if                                                  9d26s16
               end do                                                   9d26s16
              end if                                                    9d26s16
             end do                                                     9d26s16
             call second(time2)                                         3d2s17
             telap=time2-time1-tovr                                     3d2s17
             telapo(3)=telapo(3)+telap                                  3d2s17
             call paraeridj(natom,ngaus,ibdat,nbasis,ihmat,iorb,noc,    9d20s16
     $            ipair,nhcolt,isym,iapair,ibstor,isstor,multh,iptoh,
     $            itprt,idwsdeb,idorel,ascale,ia,ixyz,ipuse,idersign,   11d10s22
     $             nbasisp,bc,ibc)                                      11d10s22
             call second(time1)                                         3d2s17
             telap=time1-time2-tovr                                     3d2s17
             telapo(2)=telapo(2)+telap                                  3d2s17
             call parajkfromhd1j(ipair,ngaus,ibdat,nhcolt,ihmat,noc,icol     2d22s16
     $          ,iorb,iapair,isstor,ibstor,multh,iptoh,imsg,jmats,kmats,      2d23s16
     $         ioooo,ionex,noc,itprt,idwsdeb,ncomp,nvirtc,ipuse,itrans,  3d3s17
     $            i4od,ionexd,kmatd,jmatd,i4od2b,ionexd2,der2e,nbasdwsc, 8d3s16
     $            isblkder,isblkxder,nsblkder,nsblkxder,                9d16s16
     $            isblkxder1,nsblkxder1,isblkder1,nsblkder1,isblkkder,  9d16s16
     $            nsblkkder,nbasisp,idorel,bc,ibc)                      12d19s22
             ibcoff=ibcsav                                              9d20s16
             do i1=1,nsymb                                              9d26s16
              do i2=1,nsymb                                             9d26s16
               do i3=1,nsymb                                            9d26s16
                iptoh(i3,i2,i1)=0                                       9d26s16
               end do                                                   9d26s16
              end do                                                    9d26s16
             end do                                                     9d26s16
             ii=0
             do isb=1,nsymb                                             9d26s16
              if(nbasdwsc(isb).gt.0)then                                9d26s16
               do isc=1,nsymb                                           9d27s16
                if(nbasdwsc(isc).gt.0)then                              11d30s16
                 isbc=multh(isc,isb)                                    9d27s16
                 do isd=1,nsymb                                         9d26s16
                  if(noc(isd).gt.0)then                                 11d30s16
                   isa=multh(isbc,isd)                                  9d28s16
                   if(nbasdwsc(isa).gt.0)then                           9d26s16
                    if(isa.le.isb)then
                     ii=ii+1                                             9d26s16
                     iptoh(isd,isc,isb)=ii                               9d26s16
                    end if
                   end if                                               9d26s16
                  end if                                                9d26s16
                 end do                                                 9d26s16
                end if                                                  9d26s16
               end do                                                   9d26s16
              end if                                                    9d26s16
             end do                                                     9d26s16
             call second(time2)                                         3d2s17
             telap=time2-time1-tovr                                     3d2s17
             telapo(3)=telapo(3)+telap                                  3d2s17
             call paraerid2b(natom,ngaus,ibdat,nbasis,ihmat,iorb,noc,   9d28s16
     $            ipair,nhcolt,isym,iapair,ibstor,isstor,multh,iptoh,
     $            itprt,                                                 3d3s17
     $            idwsdeb,idorel,ascale,ia,ixyz,idersign,nbasisp,bc,ibc)11d10s22
             call second(time1)                                         3d2s17
             telap=time1-time2-tovr                                     3d2s17
             telapo(2)=telapo(2)+telap                                  3d2s17
             call parajkfromhd2b(ipair,ngaus,ibdat,nhcolt,ihmat,noc,icol9d28s16
     $          ,iorb,iapair,isstor,ibstor,multh,iptoh,imsg,itprt,       3d3s17
     $            idwsdeb,                                               3d3s17
     $            ncomp,nvirtc,i4od2b,ionexd2,der2e,nbasdwsc,           9d28s16
     $            isblkder,isblkxder,nsblkder,nsblkxder,                9d28s16
     $            isblkxder1,nsblkxder1,isblkder1,nsblkder1,isblkkder,  9d16s16
     $            nsblkkder,ipuse,nbasisp,bc(igoal),bc,ibc,idorel)      3d6s23
             ibcoff=ibcsav                                              9d20s16
             do i1=1,nsymb                                              9d26s16
              do i2=1,nsymb                                             9d26s16
               do i3=1,nsymb                                            9d26s16
                iptoh(i3,i2,i1)=0                                       9d26s16
               end do                                                   9d26s16
              end do                                                    9d26s16
             end do                                                     9d26s16
             ii=0
             do isb=1,nsymb                                             9d26s16
              if(nbasdwsc(isb).gt.0)then                                9d26s16
               do isc=1,nsymb                                           9d27s16
                if(noc(isc).gt.0)then                                   9d26s16
                 isbc=multh(isc,isb)                                    9d27s16
                 do isd=1,nsymb                                         9d26s16
                  if(noc(isd).gt.0)then                                 9d26s16
                   isa=multh(isbc,isd)                                  9d28s16
                   if(nbasdwsc(isa).gt.0)then                           9d26s16
                    ii=ii+1                                             9d26s16
                    iptoh(isd,isc,isb)=ii                               9d26s16
                   end if                                               9d26s16
                  end if                                                9d26s16
                 end do                                                 9d26s16
                end if                                                  9d26s16
               end do                                                   9d26s16
              end if                                                    9d26s16
             end do                                                     9d26s16
             call second(time2)                                         3d2s17
             telap=time2-time1-tovr                                     3d2s17
             telapo(3)=telapo(3)+telap                                  3d2s17
             call paraerid2c(natom,ngaus,ibdat,nbasis,ihmat,iorb,noc,   9d28s16
     $            ipair,nhcolt,isym,iapair,ibstor,isstor,multh,iptoh,
     $            itprt,                                                 3d3s17
     $            idwsdeb,idorel,ascale,ia,ixyz,idersign,nbasisp,bc,ibc)11d10s22
             call second(time1)                                         3d2s17
             telap=time1-time2-tovr                                     3d2s17
             telapo(2)=telapo(2)+telap                                  3d2s17
             call parajkfromhd2c(ipair,ngaus,ibdat,nhcolt,ihmat,noc,icol9d28s16
     $          ,iorb,iapair,isstor,ibstor,multh,iptoh,imsg,itprt,       3d3s17
     $            idwsdeb,                                               3d3s17
     $            ncomp,nvirtc,i4od2b,ionexd2,der2e,nbasdwsc,           9d28s16
     $            isblkder,isblkxder,nsblkder,nsblkxder,                9d28s16
     $            isblkxder1,nsblkxder1,isblkder1,nsblkder1,isblkkder,  9d16s16
     $            nsblkkder,i4od,nbasisp,bc(igoal),bc,ibc,idorel)       3d7s23
             ibcoff=ibcsav                                              9d20s16
             ii=0
             do isb=1,nsymb                                             9d26s16
              if(nbasdwsc(isb).gt.0)then                                9d26s16
               do isc=1,nsymb                                           9d27s16
                if(nbasdwsc(isc).gt.0)then                              10d3s16
                 isbc=multh(isc,isb)                                    9d27s16
                 do isd=1,nsymb                                         9d26s16
                  if(noc(isd).gt.0)then                                 9d26s16
                   isa=multh(isbc,isd)                                  9d28s16
                   if(nbasdwsc(isa).gt.0)then                           9d26s16
                     ii=ii+1                                             9d26s16
                     iptoh(isd,isc,isb)=ii                               9d26s16
                   end if                                               9d26s16
                  end if                                                9d26s16
                 end do                                                 9d26s16
                end if                                                  9d26s16
               end do                                                   9d26s16
              end if                                                    9d26s16
             end do                                                     9d26s16
             call second(time2)                                         3d2s17
             telap=time2-time1-tovr                                     3d2s17
             telapo(3)=telapo(3)+telap                                  3d2s17
             call paraerid(natom,ngaus,ibdat,nbasis,ihmat,iorb,noc,     3d28s16
     $            ipair,nhcolt,isym,iapair,ibstor,isstor,multh,iptoh,
     $            itprt,idwsdeb,idorel,ascale,ia,ixyz,1,idersign,2,     11d10s22
     $            nbasisp,bc,ibc)                                       11d10s22
             call second(time1)                                         3d2s17
             telap=time1-time2-tovr                                     3d2s17
             telapo(2)=telapo(2)+telap                                  3d2s17
             call parajkfromhd1(ipair,ngaus,ibdat,nhcolt,ihmat,noc,icol,     2d22s16
     $           iorb,iapair,isstor,ibstor,multh,iptoh,imsg,jmats,kmats,      2d23s16
     $         ioooo,ionex,noc,itprt,idwsdeb,ncomp,nvirtc,1,itrans,      3d3s17
     $            i4od2,ionexd2,kmatd,jmatd,i4od2b,ionexd2,der2e2,      12d28s16
     $            nbasdwsc,                                             10d3s16
     $            isblkder1,isblkxder1,nsblkder1,nsblkxder1,            10d6s16
     $            isblkxder1,-2,isblkder1,0,isblkkder,                  5d19s22
     $            0,nbasisp,2,igoal,idorel,bc,ibc)                      11d10s22
             derhf=derpotn+dere1+der2e
             if(lprint)write(6,*)('derivative of hf energy: '),derhf,
     $            derpotn,dere1,der2e
             if(derhf.ne.derhf)then
              call dws_sync
              call dws_finalize
              stop
             end if
             itry=iad+ieder-ipropmat                                    4d7s16
             bc(itry)=derhf                                             4d7s16
c
c     for 2nd der
c
            ibsv=ibcoff                                                 12d2s16
            call second(time2)                                          3d2s17
            telap=time2-time1-tovr                                      3d2s17
            telapo(3)=telapo(3)+telap                                   3d2s17
            call deramat(i4od2,ionexd2,iamat,noc,nvirtc,1,ih0der2,      12d2s16
     $          multh,idwsdeb,isblkder1,isblkxder1,nsblkder1,nsblkxder1,11d14s22
     $           bc,ibc)                                                11d14s22
            ioff=0                                                      12d2s16
            do isb=1,nsymb                                              12d2s16
             if(idwsdeb.gt.10)then                                      12d29s16
              write(6,*)('for symmetry block '),isb                      12d2s16
              write(6,*)('2nd der of amat '),iamat(isb)                             12d2s16
              call prntm2(bc(iamat(isb)),nvirtc(isb),noc(isb),           12d2s16
     $             nvirtc(isb))                                          12d2s16
              if(nsymb.eq.1)call printa(bc(iamat(isb)),nvirtc,noc,
     $             1,0,noc,0,1,0,bc(ibcoff))
              iamatp=iamat(isb)+nvirtc(isb)*noc(isb)
              write(6,*)('iamatp '),iamatp
              call mpprnt2(bc(iamatp),noc(isb))
             end if                                                     12d29s16
             isk=multh(isb,ipuse)                                       12d2s16
             ioff=ioff+nbasdwsc(isb)*(nbasdwsc(isb)+nbasdwsc(isk))      12d2s16
            end do                                                      12d2s16
            ioff=ioff+npropmat                                          12d2s16
            do isb=1,nsymb                                              12d2s16
             do i=0,noc(isb)*nvirtc(isb)-1                              12d2s16
              bc(ioff+i)=bc(iamat(isb)+i)                               12d2s16
             end do                                                     12d2s16
             ioff=ioff+noc(isb)*nvirtc(isb)                             12d2s16
            end do                                                      12d2s16
            ibcoff=ibsv                                                 12d2s16
c
c     for 1st der
c
             call deramat(i4od,ionexd,iamat,noc,nvirtc,ipuse,ih0der,       3d22s16
     $         multh,idwsdeb,isblkder,isblkxder,nsblkder,nsblkxder,     11d14s22
     $             bc,ibc)                                              11d14s22
             do isb=1,nsymb
              isk=multh(isb,ipuse)
              if(noc(isk)*nvirtc(isb).gt.0)then
              end if
             end do
c
c     first der of Hessian
c
             call buildhessdx(i4od,jmatd,kmatd,bc(ih0der),noc,ihessd,    12d2s16
     $            nvirtc,ipuse,multh,idwsdeb,isblkkder,nsblkkder,       12d2s16
     $            isblkder,nsblkder,bc,ibc)                             11d14s22
             noff=0
             mpropmat=npropmat
             do isb=1,nsymb
              isk=multh(ipuse,isb)
              if(isk.ge.isb)then
               mok=mpropmat                                             4d28s22
               mpropmat=mpropmat+nbasdwsc(isb)*nbasdwsc(isk)            4d28s22
               if(isk.ne.isb)then
                do i=0,nbasdwsc(isk)-1                                  4d28s22
                 do j=0,nbasdwsc(isb)-1                                 4d28s22
                  ji=mok+j+nbasdwsc(isb)*i                              4d28s22
                  ij=mpropmat+i+nbasdwsc(isk)*j                         4d28s22
                  bc(ij)=-bc(ji)                                        4d28s22
                 end do                                                 4d28s22
                end do                                                  4d28s22
                mpropmat=mpropmat+nbasdwsc(isb)*nbasdwsc(isk)            4d28s22
               end if
              end if
             end do                                                     4d28s22
             do isb=1,nsymb
              noff=noff+nbasdwsc(isb)*nbasdwsc(isb)
             end do
             do isb=1,nsymb                                             4d28s22
              noff=noff+nvirtc(isb)*noc(isb)
             end do                                                     4d28s22
             call second(time1)                                         3d2s17
             telap=time1-time2-tovr                                     3d2s17
             telapo(4)=telapo(4)+telap                                  3d2s17
              call orbder(iamat,igmat,noc,nvirtc,ipuse,multh,morb,      3d29s16
     $            bc(npropmat),nbasdwsc,dere1,bc(ih0mo),ihessd,1,       12d29s16
     $            bc(ibodc+ia-1),idwsdeb,ider3,ionex,i3x,bodcfact,bc,   11d10s22
     $            ibc,npass)                                            3d2s23
              call second(time2)                                        3d2s17
              telap=time2-time1-tovr                                    3d2s17
              telapo(5)=telapo(5)+telap                                 3d2s17
             noff=0
             mpropmat=npropmat
             do isb=1,nsymb
              isk=multh(ipuse,isb)
              if(isk.ge.isb)then
               mpropmat=mpropmat+nbasdwsc(isb)*nbasdwsc(isk)            4d28s22
               if(isk.ne.isb)then
                mpropmat=mpropmat+nbasdwsc(isb)*nbasdwsc(isk)            4d28s22
               end if
              end if
             end do                                                     4d28s22
             do isb=1,nsymb
              noff=noff+nbasdwsc(isb)*nbasdwsc(isb)
             end do
             do isb=1,nsymb                                             4d28s22
              noff=noff+nvirtc(isb)*noc(isb)
             end do                                                     4d28s22
            end if                                                      3d28s16
            if(npass.eq.2)then                                          2d17s23
             ipuse=multh(ipuse,iapair(2,ia))                               2d4s16
            end if                                                      2d17s23
           end do                                                         2d3s16
           if(ipsym.ne.1)then                                           4d4s22
            ibcoff=ibcgmat                                               4d3s22
           end if                                                       4d4s22
          end if                                                          2d3s16
         end do
        end do
        if(ipsym.eq.1)then                                              12d12s16
        if(lprint)                                                      3d2s17
     $   write(6,*)('now we have all the rhs for the 2nd der matrix '), 3d2s17
     $       ('elements')
         do ixyz=1,3                                                    12d12s16
          ixyzm=ixyz-1                                                  12d12s16
          if(lprint)write(6,*)('for ixyz = '),ixyz                      3d2s17
          do ia=1,natom                                                 12d12s16
           if(iapair(1,ia).ge.0)then                                    12d12s16
            if(lprint)write(6,*)('for nucleus no '),ia                   4d19s23
            ipuse=ipropsym(ixyz)                                        12d12s16
            if(iapair(1,ia).eq.0)then                                    12d12s16
             npass=1                                                     12d12s16
            else                                                         12d12s16
             npass=2
            end if                                                       12d12s16
            bodcfact=1d0                                                5d4s22
            do ipass=1,npass                                            12d12s16
             if(npass.eq.1)then
              cphas=' '                                                    2d3s16
              iad=ipropmat+ixyzm+3*(ia-1)                               3d29s16
              npropmat=ibc(iad)                                         3d29s16
             else if(ipass.eq.1)then                                       2d3s16
              cphas='+'                                                    2d3s16
              if(lprint)                                                4d19s23
     $             write(6,*)('symmetric combination with nucleus '),    4d19s23
     $             iapair(1,ia)
              iad=ipropmat+ixyzm+3*(ia-1)                               3d29s16
              npropmat=ibc(iad)                                         3d29s16
              bodcfact=0.5d0                                            5d4s22
             else
              if(lprint)                                                3d2s17
     $            write(6,*)('antisymmetric combination with nucleus '),4d19s23
     $             iapair(1,ia)
              cphas='-'                                                    2d3s16
              iad=ipropmat+ixyzm+3*(iapair(1,ia)-1)                     3d29s16
              npropmat=ibc(iad)                                         3d29s16
              bodcfact=0.5d0                                            5d4s22
             end if                                                        2d3s16
             if(lprint)then                                             4d28s22
              if(iapair(1,ia).eq.0)then                                 4d28s22
               write(6,1)olab(ixyz),ia,ipuse                            4d28s22
              else                                                      4d28s22
               write(6,11)olab(ixyz),ia,cphas,iapair(1,ia),ipuse        4d28s22
              end if                                                    4d28s22
             end if                                                     4d28s22
             noff=0
             do isb=1,nsymb                                             12d12s16
              isk=multh(ipuse,isb)                                      12d12s16
              npropmat=npropmat+nbasdwsc(isb)*nbasdwsc(isk)             12d12s16
              noff=noff+nbasdwsc(isb)*nbasdwsc(isb)                     12d12s16
             end do                                                     12d12s16
             do isb=1,nsymb                                             12d12s16
              iamat(isb)=npropmat+noff                                  12d12s16
              noff=noff+noc(isb)*nvirtc(isb)                            12d12s16
             end do                                                     12d12s16
              iflgu=2
              call second(time1)                                        3d2s17
             call orbder(iamat,igmat,noc,nvirtc,1,multh,morb,           12d12s16
     $            bc(npropmat),nbasdwsc,dere1,bc(ih0mo),ihessd,iflgu,       12d29s16
     $            bc(ibodc+ia-1),idwsdeb,ider3,ionex,i3x,bodcfact,bc,   11d10s22
     $             ibc,npass)                                           3d2s23
             call second(time2)                                         3d2s17
             telap=time2-time1-tovr                                     3d2s17
             telapo(5)=telapo(5)+telap                                  3d2s17
             if(npass.eq.2)then                                         2d17s23
              ipuse=multh(ipuse,iapair(2,ia))                            12d12s16
             end if                                                     2d17s23
            end do                                                      12d12s16
           end if                                                       12d12s16
          end do                                                        12d12s16
         end do                                                         12d12s16
        end if                                                          12d12s16
       end if
      end do
      if(lprint)then                                                    3d2s17
       do ia=1,natom                                                     12d29s16
        write(6,*)('>bodc for nucleus no. '),ia,bc(ibodc+ia-1)           4d19s23
       end do                                                            12d29s16
      end if                                                            3d2s17
      call dws_gsumf(telapo,5)                                          3d2s17
      if(lprint)write(6,3323)telapo                                     3d2s17
 3323 format('total time for 1e der ints: ',es10.3,/,                   3d2s17
     $       'total time for 2e der ints: ',es10.3,/,                   3d2s17
     $       'total time for trans      : ',es10.3,/,                   3d2s17
     $       'total time for manip      : ',es10.3,/,                   3d2s17
     $       'total time for orbder     : ',es10.3)                     3d2s17
      return
      end
