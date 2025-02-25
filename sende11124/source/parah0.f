c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine parah0(natom,ngaus,ibdat,nbasis,h0,ovr,isym,iapair,    5d3s10
     $                  ibstor,isstor,iso,nbb,idwsdeb,idorel,ascale,    8d29s22
     $     iftype,fstgth,ndfld,bc,ibc)                                  12d6s23
      implicit real*8 (a-h,o-z)
      external second
      integer*8 nn8                                                     2d19s10
      include "common.hf"
      include "common.store"
      include "common.spher"
      logical log1                                                      6d6s18
      integer*8 ibstor,isstor
      parameter (idf=9)                                                 8d29s22
      include "common.basis"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension h0(*),ovr(1),isym(3,1),iapair(3,1),ibstor(1),isstor(1), 5d3s10
     $     iso(8),cartb(3),cartk(3),ifdata(3,idf),iftype(*),fstgth(*)   12d6s23
      data ifdata/1,0,0, 0,1,0, 0,0,1, 1,1,0, 1,0,1, 0,1,1, 2,0,0,      8d29s22
     $     0,2,0, 0,0,2/                                                8d29s22
      ibcoffo=ibcoff                                                    2d19s10
      ib2=ibcoff                                                        1d11s23
      do ifi=1,ndfld                                                    12d6s23
       if(iftype(ifi).gt.idf)then                                       12d6s23
        write(6,*)('bad iftype in parah0: '),iftype(ifi)                12d6s23
        call dws_synca                                                   8d29s22
        call dws_finalize                                                8d29s22
        stop 'parah0'                                                   12d6s23
       end if                                                            8d29s22
      end do                                                            12d6s23
      if(idwsdeb.gt.100)then                                            5d25s18
       write(6,*)('in parah0 ')
      write(6,*)('ibdat '),ibdat
      call prntm2(bc(ibdat),ngaus,8,ngaus)
      end if                                                            5d25s18
      call second(time1)
      itry=13*25
c
c     build shell pair order
c
      srh=sqrt(0.5d0)                                                   5d3s10
      ascale2=ascale*2d0                                                8d20s15
      nn=(ngaus*(ngaus+1))/2                                            2d19s10
      nn2=nn*2                                                          2d19s10
      if(idorel.eq.0)then                                               8d20s15
c
c     overlap and kinetic energy together
c     and nuclear attraction for each atom
c
       ndo=1+natom                                                       2d19s10
      else                                                              8d20s15
c
c     overlap and kinetic energy together
c     and nuclear attraction
c     and pxVpx, pyVpy, pzVpz for each atom
c
       nbasall2=nbasall*2                                               8d20s15
       ndo=1+natom*4                                                    8d20s15
      end if                                                            8d20s15
      ndom=ndo-1                                                        8d29s22
      ipair=ibcoff                                                      2d19s10
      ibcoff=ipair+nn*ndo*3                                             2d19s10
      i12=ibcoff                                                        2d19s10
      ibcoff=i12+nn*4                                                   2d19s10
      call enough('parah0.  1',bc,ibc)
      j12=i12                                                           2d19s10
      jbdat=ibdat-1                                                     2d19s10
      ngaus7=ngaus*7                                                    5d3s10
      do i1=1,ngaus                                                     2d19s10
       do i2=1,i1                                                       2d19s10
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
      if(idwsdeb.gt.100)then                                            5d25s18
       write(6,*)('nn8: '),nn8
      end if                                                            5d25s18
      call idsortdws(ibc(i12),ibc(is12),nn8)                               1d18s23
      ii1=i12+nn                                                        2d19s10
      ii2=i12+nn2                                                       2d19s10
      jpair=ipair                                                       2d19s10
      do i=1,nn                                                         2d19s10
       j=ibc(is12+i-1)-1                                                2d19s10
    1  format(i5,i8,2i5)                                                2d19s10
       do k=1,ndo                                                       2d19s10
        ibc(jpair)=ibc(ii1+j)                                           2d19s10
        jpair=jpair+1                                                   2d19s10
        ibc(jpair)=ibc(ii2+j)                                           2d19s10
        jpair=jpair+1                                                   2d19s10
        ibc(jpair)=k-1                                                  2d22s10
        jpair=jpair+1                                                   2d19s10
       end do                                                           2d19s10
      end do                                                            2d19s10
      nneed=nn*ndo                                                      2d19s10
c
c     If we are doing a relativistic calculation, then parahf multiplied
c     nbb by 4, so we have room to store 2 component parts of ham and
c     overlap.
c
      if(idwsdeb.gt.100)then                                            5d25s18
       write(6,*)('nbb: '),nbb
      end if                                                            5d25s18
      do i=1,nbb                                                        12d28s19
       h0(i)=0d0
       ovr(i)=0d0                                                       2d24s10
      end do                                                            2d19s10
      call second(time2)
      telap=time2-time1
      loopit=0                                                          5d25s18
      do i=1+mynowprog,nneed,mynprocg                                   2d19s10
       loopit=loopit+1                                                  5d25s18
       jpair=ipair+3*(i-1)                                              2d19s10
       if(idwsdeb.gt.10)write(6,2)i,ibc(jpair),ibc(jpair+1),ibc(jpair+2),time1                 2d19s
    2  format('I am going to do ',i5,3i3,f18.5)                               2d19s10
       jbra=jbdat+ibc(jpair)                                            2d19s10
       jbra2=jbra+ngaus                                                 2d19s10
       jbra3=jbra2+ngaus                                                2d19s10
       jbra4=jbra3+ngaus                                                2d19s10
       jbra5=jbra4+ngaus                                                2d19s10
       jbra6=jbra5+ngaus                                                2d19s10
       jbra7=jbra6+ngaus                                                2d19s10
       jbra8=jbra7+ngaus                                                5d4s10
       jket=jbdat+ibc(jpair+1)                                          2d19s10
       jket2=jket+ngaus                                                 2d19s10
       jket3=jket2+ngaus                                                2d19s10
       jket4=jket3+ngaus                                                2d19s10
       jket5=jket4+ngaus                                                2d19s10
       jket6=jket5+ngaus                                                2d19s10
       jket7=jket6+ngaus                                                2d19s10
       jket8=jket7+ngaus                                                5d4s10
        nbra=2*ibc(jbra)+1                                              5d3s10
        nbraa=nbra                                                      5d3s10
        if(iapair(1,ibc(jbra8)).gt.0)nbraa=nbra*2                       5d3s10
        nket=2*ibc(jket)+1                                              5d3s10
        nketa=nket                                                      5d3s10
        if(iapair(1,ibc(jket8)).gt.0)nketa=nket*2                       5d3s10
        itmp1=ibcoff                                                    5d3s10
        itmp2=itmp1+nbraa*nketa                                         5d3s10
        ibcoff=itmp2+nbraa*nketa                                        5d3s10
        call enough('parah0.  2',bc,ibc)
       if(ibc(jpair+2).lt.1)then                                        2d22s10
        call onei(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),    2d19s10
     $      bc(jbra7),ibc(jbra4),ibc(jket),bc(jket2),bc(jket3),         2d19s10
     $      bc(jket5),bc(jket6),bc(jket7),ibc(jket4),ovr,h0,nsymb,ib1,  5d3s10
     $      ib2,bc,ibc)                                                 11d9s22
        npass=2                                                         5d7s10
      call second(time8)
       else                                                             5d7s10
        ia=ibc(jpair+2)                                                 2d22s10
        if(idorel.eq.0)then                                             8d20s15
         call nattrac(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),          8d20s15
     $        bc(jbra6),bc(jbra7),ibc(jbra4),ibc(jket),bc(jket2),       8d20s15
     $        bc(jket3),bc(jket5),bc(jket6),bc(jket7),ibc(jket4),       8d20s15
     $        atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),h0,nsymb, 8d20s15
     $        ib1,bc,ibc)                                               11d9s22
         if(ia.eq.1.and.ndfld.ne.0)then                                 12d6s23
          ibcoff=ib1+nbra*nket                                          12d23s22
          do ifi=1,ndfld                                                12d6s23
           call onep(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),    8d29s22
     $      bc(jbra7),ibc(jket),bc(jket2),bc(jket3),                    12d15s19
     $      bc(jket5),bc(jket6),bc(jket7),nsymb,kb1,0,0,0,              12d23s22
     $      ifdata(1,iftype(ifi)),ifdata(2,iftype(ifi)),                12d6s23
     $      ifdata(3,iftype(ifi)),0,0,0,bc,ibc)                         12d6s23
           do j=0,nbra*nket-1                                              8d29s22
            bc(ib1+j)=bc(ib1+j)-bc(kb1+j)*fstgth(ifi)                   12d6s23
           end do                                                        12d23s22
          end do                                                        12d6s23
          ibcoff=ib1                                                    12d23s22
         end if                                                         12d23s22
        else                                                            8d20s15
         if(ia.le.natom)then                                            8d20s15
          zm=-atnum(1,ia)                                               8d20s15
          call erir(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),   8d20s15
     $         bc(jbra7),ibc(jket),bc(jket2),bc(jket3),                 6d7s21
     $         bc(jket5),bc(jket6),bc(jket7),                           6d7s21
     $         0,atnum(2,ia),atnum(3,ia),xcart(1,ia),xcart(2,ia),       8d20s15
     $         xcart(3,ia),                                             6d7s21
     $         0,atnum(2,ia),atnum(3,ia),xcart(1,ia),xcart(2,ia),       8d20s15
     $         xcart(3,ia),ib1,0,zm,bc,ibc)                             11d14s22
          if(ia.eq.1.and.ndfld.ne.0)then                                12d6s23
           ibcoff=ib1+nbra*nket                                          12d23s22
           do ifi=1,ndfld                                               12d6s23
            call onep(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),    8d29s22
     $      bc(jbra7),ibc(jket),bc(jket2),bc(jket3),                    12d15s19
     $      bc(jket5),bc(jket6),bc(jket7),nsymb,kb1,0,0,0,              12d23s22
     $      ifdata(1,iftype(ifi)),ifdata(2,iftype(ifi)),                12d6s23
     $      ifdata(3,iftype(ifi)),0,0,0,bc,ibc)                         12d6s23
            do j=0,nbra*nket-1                                              8d29s22
             bc(ib1+j)=bc(ib1+j)-bc(kb1+j)*fstgth(ifi)                  12d6s23
            end do                                                        12d23s22
           end do                                                       12d6s23
           ibcoff=ib1                                                    12d23s22
          end if                                                         12d23s22
         else
          iderx=0                                                       8d20s15
          idery=0                                                       8d20s15
          iderz=0                                                       8d20s15
          if(ia.le.natom*2)then                                         8d20s15
           iderx=1
           iau=ia-natom                                                 8d20s15
          else if(ia.le.natom*3)then                                    8d20s15
           idery=1                                                      8d20s15
           iau=ia-natom*2                                               8d20s15
          else
           iderz=1
           iau=ia-natom*3
          end if                                                        8d20s15
          ioffa=ibc(jbra4)+nbasall                                      8d20s15
          ioffb=ibc(jket4)+nbasall                                      8d20s15
          zm=-atnum(1,iau)                                              8d20s15
          log1=.false.                                                  6d6s18
          call derid(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),  8d20s15
     $         bc(jbra7),ioffa,ibc(jket),bc(jket2),bc(jket3),           8d20s15
     $         bc(jket5),bc(jket6),bc(jket7),ioffb,                     8d20s15
     $         0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),   8d20s15
     $         xcart(3,iau),0,                                           8d20s15
     $         0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),       8d20s15
     $         xcart(3,iau),0,dum,nbasall2,nsymb,ib1,iderx,idery,       1d11s23
     $         iderz,iderx,idery,iderz,0,0,0,0,0,0,1,log1,zm,bc,ibc)    11d14s22
          if(iau.eq.1.and.ndfld.ne.0)then                               12d6s23
           ibcoff=ib1+nbra*nket                                          12d23s22
           do ifi=1,ndfld                                               12d6s23
            call onep(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),    8d29s22
     $      bc(jbra7),ibc(jket),bc(jket2),bc(jket3),                    12d15s19
     $      bc(jket5),bc(jket6),bc(jket7),nsymb,kb1,iderx,idery,iderz,              12d23s22
     $      ifdata(1,iftype(ifi)),ifdata(2,iftype(ifi)),                12d6s23
     $      ifdata(3,iftype(ifi)),iderx,idery,iderz,bc,ibc)             12d6s23
            do j=0,nbra*nket-1                                              8d29s22
             bc(ib1+j)=bc(ib1+j)-bc(kb1+j)*fstgth(ifi)                  12d6s23
            end do                                                        12d23s22
           end do                                                       12d6s23
           ibcoff=ib1                                                    12d23s22
          end if                                                         12d23s22
         end if
        end if                                                          8d20s15
        npass=1                                                         5d7s10
       end if                                                           5d7s10
       if(nsymb.ne.-1)then                                               5d3s10
        ibu=ib1                                                         5d7s10
        itmpu=itmp1                                                     5d7s10
        do ipass=1,npass                                                5d7s10
         do ik=1,nket                                                   5d7s10
          do ib=1,nbra                                                  5d7s10
           iad1=ik-1+nket*(ib-1)                                        5d7s10
           iad2=ib-1+nbraa*(ik-1)                                       5d7s10
           bc(itmpu+iad2)=bc(ibu+iad1)                                  8d20s15
          end do                                                        5d7s10
         end do                                                         5d7s10
         ibu=ib2                                                        5d7s10
         itmpu=itmp2                                                    5d7s10
        end do                                                          5d7s10
        if(nketa.ne.nket)then                                           5d3s10
         cartk(1)=bc(jket5)*dfloat(isym(1,iapair(2,ibc(jket8))))        5d5s10
         cartk(2)=bc(jket6)*dfloat(isym(2,iapair(2,ibc(jket8))))        5d5s10
         cartk(3)=bc(jket7)*dfloat(isym(3,iapair(2,ibc(jket8))))        5d5s10
         if(ibc(jpair+2).lt.1)then                                        2d22s10
          call onei(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),   5d3s10
     $     bc(jbra7),ibc(jbra4),ibc(jket),bc(jket2),bc(jket3),          2d19s10
     $         cartk,cartk(2),cartk(3),ibc(jket4),ovr,h0,nsymb,ib1,         5d3s10
     $         ib2,bc,ibc)                                              11d9s22
         else
          ia=ibc(jpair+2)                                                 2d22s10
          if(idorel.eq.0)then                                           8d20s15
           call nattrac(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),
     $    bc(jbra6),bc(jbra7),ibc(jbra4),ibc(jket),bc(jket2),bc(jket3),         2d19s10
     $     cartk,cartk(2),cartk(3),ibc(jket4),atnum(1,ia),              8d17s15
     $        xcart(1,ia),xcart(2,ia),xcart(3,ia),h0,nsymb,ib1,bc,ibc)  11d9s22
           if(ia.eq.1.and.ndfld.ne.0)then                               12d6s23
            ibcoff=ib1+nbra*nket                                        12d23s22
            do ifi=1,ndfld                                              12d6s23
             call onep(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),
     $            bc(jbra6),bc(jbra7),ibc(jket),bc(jket2),bc(jket3),    12d6s23
     $            cartk,cartk(2),cartk(3),nsymb,kb1,0,0,0,                    8d29s22
     $            ifdata(1,iftype(ifi)),ifdata(2,iftype(ifi)),                12d6s23
     $            ifdata(3,iftype(ifi)),0,0,0,bc,ibc)                         12d6s23
             do j=0,nbra*nket-1                                              8d29s22
              bc(ib1+j)=bc(ib1+j)-bc(kb1+j)*fstgth(ifi)                 12d6s23
             end do                                                          8d29s22
            end do                                                      12d6s23
            ibcoff=ib1                                                  12d23s22
           end if                                                       12d23s22
          else                                                          8d20s15
           if(ia.le.natom)then                                          8d20s15
            call erir(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),  8d20s15
     $        bc(jbra7),ibc(jket),bc(jket2),bc(jket3),                  6d7s21
     $        cartk,cartk(2),cartk(3),                                  6d7s21
     $        0,atnum(2,ia),atnum(3,ia),xcart(1,ia),xcart(2,ia),        8d20s15
     $        xcart(3,ia),                                              6d7s21
     $        0,atnum(2,ia),atnum(3,ia),xcart(1,ia),xcart(2,ia),        8d20s15
     $        xcart(3,ia),ib1,0,zm,bc,ibc)                              11d14s22
            if(ia.eq.1.and.ndfld.ne.0)then                              12d6s23
             ibcoff=ib1+nbra*nket                                        12d23s22
             do ifi=1,ndfld                                             12d6s23
              call onep(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),         12d23s22
     $            bc(jbra6),bc(jbra7),ibc(jket),bc(jket2),bc(jket3),    12d23s22
     $      cartk,cartk(2),cartk(3),nsymb,kb1,0,0,0,                    8d29s22
     $      ifdata(1,iftype(ifi)),ifdata(2,iftype(ifi)),                12d6s23
     $      ifdata(3,iftype(ifi)),0,0,0,bc,ibc)                         12d6s23
              do j=0,nbra*nket-1                                              8d29s22
               bc(ib1+j)=bc(ib1+j)-bc(kb1+j)*fstgth(ifi)                12d6s23
              end do                                                          8d29s22
             end do                                                     12d6s23
             ibcoff=ib1                                                  12d23s22
            end if                                                       12d23s22
           else                                                         8d20s15
            iderx=0                                                      8d20s15
            idery=0                                                      8d20s15
            iderz=0                                                      8d20s15
            if(ia.le.natom*2)then                                        8d20s15
             iderx=1
             iau=ia-natom                                                8d20s15
            else if(ia.le.natom*3)then                                   8d20s15
             idery=1                                                     8d20s15
             iau=ia-natom*2                                              8d20s15
            else
             iderz=1
             iau=ia-natom*3
            end if                                                       8d20s15
            call derid(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),         8d20s15
     $        bc(jbra6),bc(jbra7),ibc(jbra4),ibc(jket),bc(jket2),         2d19s10
     $        bc(jket3),cartk,cartk(2),cartk(3),ibc(jket4),             8d20s15
     $        0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),    8d20s15
     $        xcart(3,iau),0,                                           8d20s15
     $        0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),       8d20s15
     $        xcart(3,iau),0,dum,nbasall,nsymb,ib1,iderx,idery,         1d23s23
     $        iderz,iderx,idery,iderz,0,0,0,0,0,0,1,.false.,zm,bc,ibc)  11d14s22
            if(iau.eq.1.and.ndfld.ne.0)then                             12d6s23
             ibcoff=ib1+nbra*nket                                        12d23s22
             do ifi=1,ndfld                                             12d6s23
              call onep(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),         12d23s22
     $            bc(jbra6),bc(jbra7),ibc(jket),bc(jket2),bc(jket3),    12d23s22
     $      cartk,cartk(2),cartk(3),nsymb,kb1,iderx,idery,iderz,        12d23s22
     $      ifdata(1,iftype(ifi)),ifdata(2,iftype(ifi)),                12d6s23
     $      ifdata(3,iftype(ifi)),iderx,idery,iderz,bc,ibc)             12d6s23
              do j=0,nbra*nket-1                                              8d29s22
               bc(ib1+j)=bc(ib1+j)-bc(kb1+j)*fstgth(ifi)                12d6s23
              end do                                                          8d29s22
             end do                                                     12d6s23
             ibcoff=ib1                                                  12d23s22
            end if                                                       12d23s22
           end if                                                       8d20s15
          end if                                                        8d20s15
         end if                                                         5d7s10
         ibu=ib1                                                        5d7s10
         itmpu=itmp1                                                    5d7s10
         do ipass=1,npass                                               5d7s10
          do ik=1,nket                                                    5d3s10
           do ib=1,nbra                                                   5d3s10
            iad1=ik-1+nket*(ib-1)                                         5d3s10
            iad2=ib-1+nbraa*(ik-1+nket)                                  5d3s10
            bc(itmpu+iad2)=bc(ibu+iad1)                                 8d20s15
           end do                                                         5d3s10
          end do                                                          5d3s10
          ibu=ib2                                                       5d7s10
          itmpu=itmp2                                                   5d7s10
         end do                                                         5d7s10
        end if                                                          5d3s10
        if(nbraa.ne.nbra)then                                           5d3s10
         cartb(1)=bc(jbra5)*dfloat(isym(1,iapair(2,ibc(jbra8))))        5d5s10
         cartb(2)=bc(jbra6)*dfloat(isym(2,iapair(2,ibc(jbra8))))        5d5s10
         cartb(3)=bc(jbra7)*dfloat(isym(3,iapair(2,ibc(jbra8))))        5d5s10
         if(ibc(jpair+2).lt.1)then                                        2d22s10
          call onei(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),        5d3s10
     $        cartb(3),ibc(jbra4),ibc(jket),bc(jket2),bc(jket3),        5d3s10
     $     bc(jket5),bc(jket6),bc(jket7),ibc(jket4),ovr,h0,nsymb,ib1,   5d3s10
     $     ib2,bc,ibc)                                                  11d9s22
         else
          ia=ibc(jpair+2)                                                 2d22s10
          if(idorel.eq.0)then                                           8d20s15
           call nattrac(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),    5d7s10
     $     cartb(3),ibc(jbra4),ibc(jket),bc(jket2),bc(jket3),           5d7s10
     $     bc(jket5),bc(jket6),bc(jket7),ibc(jket4),atnum(1,ia),        8d17s15
     $        xcart(1,ia),xcart(2,ia),xcart(3,ia),h0,nsymb,ib1,bc,ibc)  11d9s22
           if(ia.eq.1.and.ndfld.ne.0)then                               12d6s23
            ibcoff=ib1+nbra*nket                                        12d23s22
            do ifi=1,ndfld                                              12d6s23
             call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),       8d29s22
     $      cartb(3),ibc(jket),bc(jket2),bc(jket3),                     8d29s22
     $      bc(jket5),bc(jket6),bc(jket7),nsymb,kb1,0,0,0,               5d27s21
     $      ifdata(1,iftype(ifi)),ifdata(2,iftype(ifi)),                12d6s23
     $      ifdata(3,iftype(ifi)),0,0,0,bc,ibc)                         12d6s23
             do j=0,nbra*nket-1                                              8d29s22
              bc(ib1+j)=bc(ib1+j)-bc(kb1+j)*fstgth(ifi)                 12d6s23
             end do                                                          8d29s22
            end do                                                      12d6s23
            ibcoff=ib1                                                  12d23s22
           end if                                                       12d23s22
          else                                                          8d20s15
           if(ia.le.natom)then                                           8d20s15
            call erir(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),      8d20s15
     $        cartb(3),ibc(jket),bc(jket2),bc(jket3),                   6d7s21
     $        bc(jket5),bc(jket6),bc(jket7),                            6d7s21
     $        0,atnum(2,ia),atnum(3,ia),xcart(1,ia),xcart(2,ia),        8d20s15
     $        xcart(3,ia),                                              6d7s21
     $        0,atnum(2,ia),atnum(3,ia),xcart(1,ia),xcart(2,ia),        8d20s15
     $        xcart(3,ia),ib1,0,zm,bc,ibc)                              11d14s22
            if(ia.eq.1.and.ndfld.ne.0)then                              12d6s23
             ibcoff=ib1+nbra*nket                                        12d23s22
             do ifi=1,ndfld                                             12d6s23
              call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),       8d29s22
     $            cartb(3),ibc(jket),bc(jket2),bc(jket3),                     8d29s22
     $            bc(jket5),bc(jket6),bc(jket7),nsymb,kb1,0,0,0,               5d27s21
     $            ifdata(1,iftype(ifi)),ifdata(2,iftype(ifi)),          12d6s23
     $            ifdata(3,iftype(ifi)),0,0,0,bc,ibc)                   12d6s23
              do j=0,nbra*nket-1                                              8d29s22
               bc(ib1+j)=bc(ib1+j)-bc(kb1+j)*fstgth(ifi)                12d6s23
              end do                                                          8d29s22
             end do                                                     12d6s23
             ibcoff=ib1                                                  12d23s22
            end if                                                       12d23s22
           else
            iderx=0                                                      8d20s15
            idery=0                                                      8d20s15
            iderz=0                                                      8d20s15
            if(ia.le.natom*2)then                                        8d20s15
             iderx=1
             iau=ia-natom                                                8d20s15
            else if(ia.le.natom*3)then                                   8d20s15
             idery=1                                                     8d20s15
             iau=ia-natom*2                                              8d20s15
            else
             iderz=1
             iau=ia-natom*3
            end if                                                       8d20s15
            call derid(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),    8d20s15
     $        cartb(3),ibc(jbra4),ibc(jket),bc(jket2),bc(jket3),        8d20s15
     $        bc(jket5),bc(jket6),bc(jket7),ibc(jket4),                 8d20s15
     $        0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),    8d20s15
     $        xcart(3,iau),0,                                           8d20s15
     $        0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),       8d20s15
     $        xcart(3,iau),0,dum,nbasall,nsymb,ib1,iderx,idery,         1d23s23
     $        iderz,iderx,idery,iderz,0,0,0,0,0,0,1,.false.,zm,bc,ibc)  11d14s22
            if(iau.eq.1.and.ndfld.ne.0)then                             12d6s23
             ibcoff=ib1+nbra*nket                                        12d23s22
             do ifi=1,4                                                 12d6s23
              call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),       8d29s22
     $            cartb(3),ibc(jket),bc(jket2),bc(jket3),                     8d29s22
     $            bc(jket5),bc(jket6),bc(jket7),nsymb,kb1,iderx,idery,  12d23s22
     $            iderz,ifdata(1,iftype(ifi)),ifdata(2,iftype(ifi)),    12d6s23
     $            ifdata(3,iftype(ifi)),iderx,idery,iderz,bc,ibc)       12d6s23
              do j=0,nbra*nket-1                                              8d29s22
               bc(ib1+j)=bc(ib1+j)-bc(kb1+j)*fstgth(ifi)                12d6s23
              end do                                                          8d29s22
             end do                                                     12d6s23
             ibcoff=ib1                                                  12d23s22
            end if                                                       12d23s22
           end if
          end if                                                         8d20s15
         end if                                                          5d7s10
         ibu=ib1                                                         5d7s10
         itmpu=itmp1                                                     5d7s10
         do ipass=1,npass                                                5d7s10
          do ik=1,nket                                                    5d3s10
           do ib=1,nbra                                                   5d3s10
            iad1=ik-1+nket*(ib-1)                                         5d3s10
            iad2=ib-1+nbra+nbraa*(ik-1)                                   5d3s10
            bc(itmpu+iad2)=bc(ibu+iad1)                                  8d20s15
           end do                                                        5d7s10
          end do                                                         5d3s10
          ibu=ib2                                                        5d7s10
          itmpu=itmp2                                                    5d7s10
         end do                                                          5d3s10
         if(nketa.ne.nket)then                                           5d3s10
          if(ibc(jpair+2).lt.1)then                                        2d22s10
           call onei(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),       5d3s10
     $        cartb(3),ibc(jbra4),ibc(jket),bc(jket2),bc(jket3),        5d3s10
     $     cartk,cartk(2),cartk(3),ibc(jket4),ovr,h0,nsymb,ib1,         5d3s10
     $     ib2,bc,ibc)                                                  11d9s22
          else
           ia=ibc(jpair+2)                                                 2d22s10
           if(idorel.eq.0)then                                          8d20s15
            call nattrac(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),    5d7s10
     $          cartb(3),ibc(jbra4),ibc(jket),bc(jket2),bc(jket3),          5d7s10
     $     cartk,cartk(2),cartk(3),ibc(jket4),atnum(1,ia),              8d17s15
     $          xcart(1,ia),xcart(2,ia),xcart(3,ia),h0,nsymb,ib1,bc,ibc)11d9s22
            if(ia.eq.1.and.ndfld.ne.0)then                              12d6s23
             ibcoff=ib1+nbra*nket                                       12d23s22
             do ifi=1,ndfld                                             12d6s23
              call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),      8d29s22
     $      cartb(3),ibc(jket),bc(jket2),bc(jket3),                     8d29s22
     $      cartk,cartk(2),cartk(3),nsymb,kb1,0,0,0,                    8d29s22
     $      ifdata(1,iftype(ifi)),ifdata(2,iftype(ifi)),                12d6s23
     $      ifdata(3,iftype(ifi)),0,0,0,bc,ibc)                         12d6s23
              do j=0,nbra*nket-1                                              8d29s22
               bc(ib1+j)=bc(ib1+j)-bc(kb1+j)*fstgth(ifi)                12d6s23
              end do                                                          8d29s22
             end do                                                     12d6s23
             ibcoff=ib1                                                 12d23s22
            end if                                                      12d23s22
           else                                                         8d20s15
            if(ia.le.natom)then                                           8d20s15
             call erir(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),     8d20s15
     $        cartb(3),ibc(jket),bc(jket2),bc(jket3),                   6d7s21
     $        cartk,cartk(2),cartk(3),                                  6d7s21
     $        0,atnum(2,ia),atnum(3,ia),xcart(1,ia),xcart(2,ia),        8d20s15
     $        xcart(3,ia),                                              6d7s21
     $        0,atnum(2,ia),atnum(3,ia),xcart(1,ia),xcart(2,ia),        8d20s15
     $        xcart(3,ia),ib1,0,zm,bc,ibc)                              11d14s22
             if(ia.eq.1.and.ndfld.ne.0)then                             12d6s23
              ibcoff=ib1+nbra*nket                                       12d23s22
              do ifi=1,ndfld                                            12d6s23
               call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),      8d29s22
     $      cartb(3),ibc(jket),bc(jket2),bc(jket3),                     8d29s22
     $      cartk,cartk(2),cartk(3),nsymb,kb1,0,0,0,                    8d29s22
     $      ifdata(1,iftype(ifi)),ifdata(2,iftype(ifi)),                12d6s23
     $      ifdata(3,iftype(ifi)),0,0,0,bc,ibc)                         12d6s23
               do j=0,nbra*nket-1                                              8d29s22
                bc(ib1+j)=bc(ib1+j)-bc(kb1+j)*fstgth(ifi)               12d6s23
               end do                                                          8d29s22
              end do                                                    12d6s23
              ibcoff=ib1                                                 12d23s22
             end if                                                      12d23s22
            else
             iderx=0                                                      8d20s15
             idery=0                                                      8d20s15
             iderz=0                                                      8d20s15
             if(ia.le.natom*2)then                                        8d20s15
              iderx=1
              iau=ia-natom                                                8d20s15
             else if(ia.le.natom*3)then                                   8d20s15
              idery=1                                                     8d20s15
              iau=ia-natom*2                                              8d20s15
             else
              iderz=1
              iau=ia-natom*3
             end if                                                       8d20s15
             call derid(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),   8d20s15
     $        cartb(3),ibc(jbra4),ibc(jket),bc(jket2),bc(jket3),         2d19s10
     $        cartk,cartk(2),cartk(3),ibc(jket4),                       8d20s15
     $        0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),    8d20s15
     $        xcart(3,iau),0,                                           8d20s15
     $        0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),       8d20s15
     $        xcart(3,iau),0,dum,nbasall,nsymb,ib1,iderx,idery,         1d23s23
     $        iderz,iderx,idery,iderz,0,0,0,0,0,0,1,.false.,zm,bc,ibc)  11d14s22
             if(iau.eq.1.and.ndfld.ne.0)then                            12d6s23
              ibcoff=ib1+nbra*nket                                       12d23s22
              do ifi=1,ndfld                                            12d6s23
               call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),      8d29s22
     $      cartb(3),ibc(jket),bc(jket2),bc(jket3),                     8d29s22
     $      cartk,cartk(2),cartk(3),nsymb,kb1,iderx,idery,iderz,        12d23s22
     $      ifdata(1,iftype(ifi)),ifdata(2,iftype(ifi)),                12d6s23
     $      ifdata(3,iftype(ifi)),iderx,idery,iderz,bc,ibc)             12d6s23
               do j=0,nbra*nket-1                                              8d29s22
                bc(ib1+j)=bc(ib1+j)-bc(kb1+j)*fstgth(ifi)               12d6s23
               end do                                                          8d29s22
              end do                                                    12d6s23
              ibcoff=ib1                                                 12d23s22
             end if                                                      12d23s22
            end if                                                      8d20s15
           end if                                                       8d20s15
          end if                                                         5d7s10
          ibu=ib1                                                        5d7s10
          itmpu=itmp1                                                    5d7s10
          do ipass=1,npass                                               5d7s10
           do ik=1,nket                                                  5d3s10
            do ib=1,nbra                                                 5d3s10
             iad1=ik-1+nket*(ib-1)                                       5d3s10
             iad2=ib-1+nbra+nbraa*(ik-1+nket)                            5d3s10
             bc(itmpu+iad2)=bc(ibu+iad1)                                 8d20s15
            end do                                                       5d7s10
           end do                                                        5d3s10
           ibu=ib2                                                       5d7s10
           itmpu=itmp2                                                   5d7s10
          end do                                                         5d3s10
         end if                                                          5d3s10
        end if                                                           5d3s10
        if(nketa.ne.nket)then
         itmpu=itmp1                                                     5d7s10
         do ipass=1,npass                                                5d7s10
          do ik=1,nket                                                   5d3s10
           do ib=1,nbraa                                                 5d3s10
            iad1=ib-1+nbraa*(ik-1)                                       5d3s10
            iad2=iad1+nbraa*nket                                         5d3s10
            sum=srh*(bc(itmpu+iad1)+bc(itmpu+iad2))                      5d3s10
            dif=srh*(-bc(itmpu+iad1)+bc(itmpu+iad2))                     5d3s10
            bc(itmpu+iad1)=sum                                           5d3s10
            bc(itmpu+iad2)=dif                                           5d3s10
           end do                                                        5d7s10
          end do                                                         5d3s10
          itmpu=itmp2                                                    5d7s10
         end do                                                          5d3s10
        end if                                                           5d3s10
        if(nbraa.ne.nbra)then
         itmpu=itmp1                                                     5d7s10
         do ipass=1,npass                                                5d7s10
          do ik=1,nketa                                                   5d3s10
           do ib=1,nbra                                                  5d3s10
            iad1=ib-1+nbraa*(ik-1)                                       5d3s10
            iad2=iad1+nbra                                               5d3s10
            sum=srh*(bc(itmpu+iad1)+bc(itmpu+iad2))                      5d3s10
            dif=srh*(-bc(itmpu+iad1)+bc(itmpu+iad2))                     5d3s10
            bc(itmpu+iad1)=sum                                           5d3s10
            bc(itmpu+iad2)=dif                                           5d3s10
           end do                                                        5d7s10
          end do                                                         5d3s10
          itmpu=itmp2                                                    5d7s10
         end do                                                          5d3s10
        end if                                                           5d3s10
        itmpu=itmp1                                                      5d7s10
        if(idorel.eq.0)then                                             8d20s15
         do ipass=1,npass                                                 5d7s10
          do ik=1,nketa                                                   5d3s10
           ikk=ik+ibc(jket4)                                              5d3s10
           isk=isstor(ikk)                                                5d3s10
           do ib=1,nbraa                                                  5d3s10
            ibb=ib+ibc(jbra4)                                             5d3s10
            iad1=ib-1+nbraa*(ik-1)                                        5d3s10
            if(isk.eq.isstor(ibb))then                                    5d3s10
             in=min(ibstor(ikk),ibstor(ibb))                              5d3s10
             ix=max(ibstor(ikk),ibstor(ibb))                              5d3s10
             iad=iso(isk)+((ix*(ix-1))/2)+in
             if(ipass.eq.npass)then                                       5d7s10
              if(ibstor(ikk).le.ibstor(ibb))then                          1d17s12
               h0(iad)=h0(iad)+bc(itmpu+iad1)                               5d3s10
              end if                                                      1d17s12
             else                                                         5d7s10
              ovr(iad)=bc(itmpu+iad1)                                     5d7s10
             end if                                                       5d7s10
            else                                                          5d3s10
             sz=abs(bc(itmpu+iad1))                                       5d7s10
             if(sz.gt.1d-10.and.npass.eq.2)then                                          5d3s10
              write(6,*)('symmetry transformation failure!!! ')           5d3s10
              write(6,*)ib,ik,ibb,ikk,sz                                  5d7s10
              write(6,*)('ipass,npass: '),ipass,npass                     5d7s10
              call prntm2(bc(itmpu),nbraa,nketa,nbraa)                    5d3s10
              stop 'parah0'                                                        5d3s10
             end if                                                       5d3s10
            end if                                                        5d3s10
           end do                                                         5d7s10
          end do                                                          5d3s10
          itmpu=itmp2                                                     5d7s10
         end do                                                           5d3s10
        else                                                            8d20s15
         ia=ibc(jpair+2)                                                 2d22s10
         do ipass=1,npass                                                 5d7s10
          do ik=1,nketa                                                   5d3s10
           ikk=ik+ibc(jket4)                                              5d3s10
           isk=isstor(ikk)                                                5d3s10
           do ib=1,nbraa                                                  5d3s10
            ibb=ib+ibc(jbra4)                                             5d3s10
            iad1=ib-1+nbraa*(ik-1)                                        5d3s10
            if(isk.eq.isstor(ibb))then                                    5d3s10
             inl=min(ibstor(ikk),ibstor(ibb))                              5d3s10
             ixl=max(ibstor(ikk),ibstor(ibb))                              5d3s10
             iadll=iso(isk)+((ixl*(ixl-1))/2)+inl
             ins=min(ibstor(ikk),ibstor(ibb))+nbasdws(isk)              8d20s15
             ixs=max(ibstor(ikk),ibstor(ibb))+nbasdws(isk)              8d20s15
             iadss=iso(isk)+((ixs*(ixs-1))/2)+ins
             ix=ibstor(ibb)+nbasdws(isk)
             in=ibstor(ikk)
             iadsl=iso(isk)+((ix*(ix-1))/2)+in
             ix=ibstor(ikk)+nbasdws(isk)
             in=ibstor(ibb)
             iadls=iso(isk)+((ix*(ix-1))/2)+in
c
c     npass=2: ipass=1 is overlap and ipass=2 is t
c     npass=1: Vnuc if ia le natom, dVnucd otherwise.
c     overlap goes into LL part of ovr,
c     t goes into LS part of h0, -t goes into SS part of h0,
c     and ascale*ascale*t goes into SS part of
c     ovr. Vnuc goes into LL part of h0 and dVnucd into SS part of h0.
c
             if(ipass.lt.npass)then
              if(ibstor(ikk).le.ibstor(ibb))then
               if(iadll.lt.1.or.iadll.gt.nbb)then
                write(6,*)('iadll is out of range: '),iadll,isk,iso(isk)
     $               ,ix,in
                stop
               end if
               ovr(iadll)=bc(itmpu+iad1)                                8d20s15
              end if                                                    8d20s15
             else if(npass.eq.2)then                                    8d20s15
              h0(iadsl)=bc(itmpu+iad1)                                  8d20s15
              h0(iadls)=bc(itmpu+iad1)                                  8d20s15
              if(ibstor(ikk).le.ibstor(ibb))then
               ovr(iadss)=bc(itmpu+iad1)*ascale2                        8d20s15
               h0(iadss)=h0(iadss)-bc(itmpu+iad1)                       8d20s15
              end if                                                    8d20s15
             else if(ia.le.natom)then                                   8d20s15
              if(ibstor(ikk).le.ibstor(ibb))then
               h0(iadll)=h0(iadll)+bc(itmpu+iad1)
              end if
             else                                                       8d20s15
              if(ibstor(ikk).le.ibstor(ibb))then                        8d20s15
               h0(iadss)=h0(iadss)+ascale*bc(itmpu+iad1)                8d20s15
              end if                                                    8d20s15
             end if                                                     8d20s15
            else                                                          5d3s10
             sz=abs(bc(itmpu+iad1))                                       5d7s10
             if(sz.gt.1d-10.and.npass.eq.2)then                                          5d3s10
              write(6,*)('symmetry transformation failure!!! ')           5d3s10
              write(6,*)ib,ik,ibb,ikk,sz                                  5d7s10
              write(6,*)('ipass,npass: '),ipass,npass                     5d7s10
              call prntm2(bc(itmpu),nbraa,nketa,nbraa)                    5d3s10
              stop                                                        5d3s10
             end if                                                       5d3s10
            end if                                                        5d3s10
           end do                                                         5d7s10
          end do                                                          5d3s10
          itmpu=itmp2                                                     5d7s10
         end do                                                           5d3s10
        end if                                                          8d20s15
        ibcoff=itmp1                                                     5d3s10
       end if                                                            5d3s10
      end do                                                            2d19s10
      call second(time3)
      call dws_sync                                                     2d22s10
      call dws_gsumf(ovr,nbb)                                           5d3s10
      call dws_gsumf(h0,nbb)                                            5d3s10
      call second(time4)
      nn=(nbasdws(1)*(nbasdws(1)+1))/2
      call dws_bcast(ovr,nn)
      if(idwsdeb.gt.10)then                                             3d16s12
      do isb=1,nsymb                                                    5d3s10
       write(6,*)('for symmetry block '),isb                            5d3s10
       write(6,*)('overlap: ')                                           2d22s10
       if(idorel.eq.0)then                                              8d20s15
        nbasz=nbasdws(isb)                                              8d20s15
       else                                                             8d20s15
        nbasz=nbasdws(isb)*2                                            8d20s15
       end if                                                           8d20s15
       call mpprnt2(ovr(iso(isb)+1),nbasz)                              8d20s15
       write(6,*)('h0: ')                                                2d22s10
       call mpprnt2(h0(iso(isb)+1),nbasz)                               8d20s15
      end do                                                            5d3s10
      end if                                                            3d16s12
      ibcoff=ibcoffo                                                    2d19s10
      return
      end                                                               2d19s10
