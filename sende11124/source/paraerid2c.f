c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine paraerid2c(natom,ngaus,ibdat,nbasis,ihmat,iorb,nocc,   9d28s16
     $     ipairab,nhcolt,isym,iapair,ibstor,isstor,multh,iptoh,iter,   3d2s16
     $     idwsdeb,idorel,ascale,iatom,ixyz,idersign,nbasisp,bc,ibc)    11d10s22
c
c     like paraerid, except form (ab'|cd')                              9d28s16
c     wrt ixyz coordinate of iatom.
c     output of this routine is ab ao and cd mo, with d only occ.
c
c
      implicit real*8 (a-h,o-z)
      character*5 lfn
      character*80 line(50)
      include "common.hf"
      include "common.store"
      include "common.basis"
      integer*8 ibstor,isstor                                           5d10s10
      logical ldeb,lpr                                                      4d21s16
      dimension iorb(1),isym(3,1),iapair(3,1),ibstor(1),isstor(1),      5d12s10
     $     carta(3),cartb(3),cartc(3),cartd(3),ieraw(8,8),itmpab(8,8),  5d12s10
     $     multh(8,8),iptoh(8,8,8),nocc(8),neraw(8,8),nueraw(8,8),      9d26s16
     $     nztype(8),nsumts(8,8),nbasisp(*)                             4d7s22
      dimension igsym(18,2),gsym(78,2)
      common/timerocm/tovr,telapo(15)                                       5d7s12
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      do isb=1,nsymb                                                    9d26s16
       nztype(isb)=0                                                    9d26s16
      end do                                                            9d26s16
      ibcoffo=ibcoff                                                    2d19s10
      if(idorel.eq.0)then                                               8d24s15
       ncomp=1                                                          8d24s15
      else                                                              8d24s15
       ncomp=2                                                          8d24s15
      end if                                                            8d24s15
      ntmpmat=2                                                         5d31s22
      if(idorel.eq.0.or.idorel.eq.1)then                                5d31s22
       nerimat=1                                                        5d31s22
       ntmpmat=1                                                        5d31s22
      else if(idorel.eq.2.or.idorel.eq.-2.or.idorel.eq.-4)then          5d31s22
       nerimat=3                                                        5d31s22
      else                                                              5d31s22
       nerimat=4                                                        5d31s22
      end if                                                            5d31s22
      idx=0                                                             2d1s16
      idy=0                                                             2d1s16
      idz=0                                                             2d1s16
      if(ixyz.eq.1)then                                                 2d1s16
       idx=1                                                            9d20s16
      else if(ixyz.eq.2)then                                            2d1s16
       idy=1                                                            9d20s16
      else                                                              2d1s16
       idz=1                                                            9d20s16
      end if                                                            2d1s16
      nbasallc=nbasall*ncomp                                            8d26s15
      srh=sqrt(0.5d0)                                                   5d12s10
      nbasx=nbasisp(1)*ncomp                                            4d7s22
      noccx=nocc(1)                                                     5d12s10
      do isb=2,nsymb                                                    5d12s10
       nbasx=max(nbasx,nbasisp(isb)*ncomp)                              4d7s22
       noccx=max(noccx,nocc(isb))                                       5d12s10
      end do                                                            5d12s10
c
c     build shell pair order
c
      nnab=ngaus*ngaus                                                  9d30s16
      nnab2=nnab*2                                                      3d2s16
      ipairab=ibcoff                                                    3d2s16
      ibcoff=ipairab+nnab2                                              3d2s16
      i12=ibcoff                                                        2d19s10
      ibcoff=i12+nnab2*2                                                3d2s16
      memneed(1:3)='i12'
      call enough('paraerid2c.  1',bc,ibc)
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
       do i2=1,ngaus                                                    9d28s16
        ibc(j12)=-((ibc(jbdat+i1)+ibc(jbdat+i2)+2)**2)                  2d19s10
        if(iapair(1,ibc(jbdat+i1+ngaus7)).gt.0)ibc(j12)=ibc(j12)*2      5d3s10
        if(iapair(1,ibc(jbdat+i2+ngaus7)).gt.0)ibc(j12)=ibc(j12)*2      5d3s10
        ibc(j12+nnab)=i1                                                  2d19s10
        ibc(j12+nnab2)=i2                                                 2d19s10
        j12=j12+1                                                       2d19s10
       end do                                                           2d19s10
      end do                                                            2d19s10
      is12=i12+nnab*3                                                     2d19s10
      nn8=nnab                                                            2d19s10
c
c     sort, so we can assign most expensive blocks first
c
      call idsortdws(ibc(i12),ibc(is12),nn8)                            1d18s23
      ii1=i12+nnab                                                        2d19s10
      ii2=i12+nnab2                                                       2d19s10
      jpairab=ipairab                                                       2d19s10
c
c     order shell pair indices
c
      do i=1,nnab                                                         2d19s10
       j=ibc(is12+i-1)-1                                                2d19s10
    1  format(i5,i8,2i5)                                                2d19s10
       ibc(jpairab)=ibc(ii1+j)                                            2d19s10
       jpairab=jpairab+1                                                    2d19s10
       ibc(jpairab)=ibc(ii2+j)                                            2d19s10
       jpairab=jpairab+1                                                    2d19s10
      end do                                                            2d19s10
      do i=1,nsdlkd                                                     3d1s16
       nhcol(i)=0                                                       5d12s10
       ihcol(i)=0                                                       5d12s10
      end do                                                            5d12s10
c
c     compute storage required on this proc
c
      idws=0                                                            1d10s11
      nsum1=0
      nsum2=0
      nsumt=0
      do i=1,nsymb
       do j=1,nsymb
        nsumts(j,i)=0
       end do
      end do
      do i=1+mynowprog,nnab,mynprocg                                      3d18s10
       jpairab=ipairab+2*(i-1)                                              2d19s10
       if(max(ibc(jpairab),ibc(jpairab+1)).gt.ngaus)then
        call dws_sync
        call dws_finalize
        stop
       end if
       la=ibc(jbdat+ibc(jpairab))                                         3d18s10
       na=ncomp*(2*la+1)                                                8d24s15
       lb=ibc(jbdat+ibc(jpairab+1))                                       3d18s10
       nb=ncomp*(2*lb+1)                                                8d24s15
       naa=na                                                           5d12s10
       nba=nb                                                           5d12s10
       if(iapair(1,ibc(jbdat+ibc(jpairab)+ngaus7)).gt.0)then              5d12s10
        naa=naa*2                                                       5d12s10
       end if                                                           5d12s10
       if(iapair(1,ibc(jbdat+ibc(jpairab+1)+ngaus7)).gt.0)then            5d12s10
        nba=nba*2                                                       5d12s10
       end if                                                           5d12s10
       nba0=nba/ncomp                                                   9d18s15
       naa0=naa/ncomp                                                   9d18s15
       do ib=1,nba0                                                     9d18s15
        ibb=ib+ibc(jbdat+ibc(jpairab+1)+ngaus3)                           6d4s10
        isb=isstor(ibb)                                                 5d12s10
        do ia=1,naa0                                                     9d18s15
         iaa=ia+ibc(jbdat+ibc(jpairab)+ngaus3)                            6d4s10
         isa=isstor(iaa)                                                 2d23s16
         nsumt=nsumt+1
         nsumts(isa,isb)=nsumts(isa,isb)+1
 3353    format(3i8,5x,3i8,5x,i8)
         isab=multh(isa,isb)                                            9d20s16
         do isd=1,nsymb                                                 5d12s10
          isc=multh(isd,isab)                                           5d12s10
          nadd=nbasdws(isc)*nocc(isd)                                   4d7s22
          if(nadd.gt.0)then                                             1d30s12
           nadd=nadd*ntmpmat                                            3d7s23
           if(iptoh(isd,isc,isb).gt.0)then                              9d20s16
            nhcol(iptoh(isd,isc,isb))=nhcol(iptoh(isd,isc,isb))+nadd    9d28s16
            nsum1=nsum1+1
           end if                                                       1d30s12
          end if                                                        1d30s12
         end do                                                         5d12s10
        end do                                                          5d12s10
       end do                                                           5d12s10
      end do                                                            3d18s10
      if(iter.eq.1)write(6,*)('this proc has ')                         3d15s12
      nhcolt=0                                                          5d12s10
      ii=0
      do isb=1,nsymb
       do isc=1,nsymb
        isbc=multh(isb,isc)
        do isd=1,nsymb
         isa=multh(isbc,isd)                                            9d28s16
         if(iptoh(isd,isc,isb).gt.0)then
          ii=ii+1
          nhcolt=nhcolt+nhcol(iptoh(isd,isc,isb))
         end if
        end do
       end do
      end do
      nhmt=ii                                                           9d26s16
      if(iter.eq.1)write(6,*)('for a total of '),nhcolt,                3d15s12
     $     (' elements of hmat')                                        3d15s12
      ihmat=ibcoff                                                      3d4s10
      ibcoff=ihmat+nhcolt                                               5d12s10
      xnan=0d0                                                          3d7s23
      do i=0,nhcolt-1
       bc(ihmat+i)=xnan
      end do
      ibcexit=ibcoff                                                    3d1s16
      memneed(1:4)='hmat'
      call enough('paraerid2c.  2',bc,ibc)
      ihcol(1)=ihmat                                                    5d12s10
      do i=2,nhmt                                                       9d26s16
       ihcol(i)=ihcol(i-1)+nhcol(i-1)                                   9d26s16
      end do                                                            9d26s16
      jhmat=ihmat                                                       3d4s10
      idoit=10
      ndoit=0
      nxtot=0
      iszz=ibcoff
      ibcoff=iszz+nsymb*nsymb
      call enough('paraerid2c.  3',bc,ibc)
      do i=0,nsymb*nsymb-1
       bc(iszz+i)=0d0
      end do
      ngothru=0                                                         9d26s16
      do i=1+mynowprog,nnab,mynprocg                                      3d4s10
       ngothru=ngothru+1
       jpair=ipairab+2*(i-1)                                              2d19s10
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
       nsza0=nsza                                                       5d31s22
       nsza=nsza*ncomp                                                  8d24s15
       nszaa=nsza                                                       5d12s10
       nszaa0=nsza0                                                     5d31s22
       if(iapair(1,ibc(ja8)).gt.0)then                                  5d11s10
        nszaa=nszaa*2                                                   5d12s10
        nszaa0=nszaa0*2                                                 5d31s22
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
       nszb0=nszb                                                       5d31s22
       nszb=nszb*ncomp                                                  8d24s15
       nszba=nszb                                                       5d12s10
       nszba0=nszb0                                                     5d31s22
 3352  format(6i5)
       if(iapair(1,ibc(jb8)).gt.0)then                                  5d11s10
        nszba=nszba*2                                                   5d12s10
        nszba0=nszba0*2                                                 5d31s22
       end if                                                           5d11s10
       ncsz=nsza*nszb                                                   3d2s16
       ncsz0=nsza0*nszb0                                                5d31s22
       ncsza0=nszaa0*nszba0                                             5d31s22
       ncsza=nszaa*nszba                                                4d21s16
       do isd=1,nsymb                                                   5d12s10
        do isc=1,nsymb                                                  5d12s10
         itmpab(isc,isd)=ibcoff                                         5d12s10
         ibcoff=itmpab(isc,isd)+ncsza0*nbasdws(isc)*nocc(isd)*ntmpmat   3d7s23
        end do                                                          5d12s10
       end do                                                           5d12s10
       memneed(1:5)='tmpab'
       call enough('paraerid2c.  4',bc,ibc)
       do ii=itmpab(1,1),ibcoff-1                                       5d31s22
        bc(ii)=0d0                                                      5d31s22
       end do                                                           5d31s22
       do isd=1,nsymb                                                   5d12s10
        do isc=1,nsymb                                                  5d12s10
         ieraw(isc,isd)=ibcoff                                          5d12s10
         ibcoff=ieraw(isc,isd)+ncsz0*nbasisp(isc)*nbasisp(isd)*nerimat  3d7s23
        end do                                                          5d12s10
       end do                                                           5d12s10
       memneed(1:4)='eraw'
       call enough('paraerid2c.  5',bc,ibc)
       do ipassb=1,nszba/nszb                                           5d12s10
        sigb=1d0                                                        4d18s16
        if(ipassb.eq.1)then
         cartb(1)=bc(jb5)                                               5d12s10
         cartb(2)=bc(jb6)                                               5d12s10
         cartb(3)=bc(jb7)                                               5d12s10
         iob=0                                                          5d12s10
        else                                                            5d12s10
         cartb(1)=bc(jb5)*dfloat(isym(1,iapair(2,ibc(jb8))))            5d11s10
         cartb(2)=bc(jb6)*dfloat(isym(2,iapair(2,ibc(jb8))))            5d11s10
         cartb(3)=bc(jb7)*dfloat(isym(3,iapair(2,ibc(jb8))))            5d11s10
         iob=nszb0                                                      3d7s23
        end if                                                          5d12s10
        sigbm=sigb                                                      3d2s23
        do ipassa=1,nszaa/nsza                                          5d12s10
         siga=1d0                                                       4d18s16
         if(ipassa.eq.1)then                                            5d12s10
          carta(1)=bc(ja5)                                              5d12s10
          carta(2)=bc(ja6)                                              8d15s11
          carta(3)=bc(ja7)                                              8d15s11
          ioa=0                                                         5d12s10
         else                                                           5d12s10
          carta(1)=bc(ja5)*dfloat(isym(1,iapair(2,ibc(ja8))))           5d12s10
          carta(2)=bc(ja6)*dfloat(isym(2,iapair(2,ibc(ja8))))           5d12s10
          carta(3)=bc(ja7)*dfloat(isym(3,iapair(2,ibc(ja8))))           5d12s10
          ioa=nsza0                                                     3d7s23
         end if                                                         5d12s10
         ioab=iob+nszba0*ioa                                            3d7s23
         if(idoit.eq.1)then
          write(6,*)('ipassa,ipassb: '),ipassa,ipassb
         end if
         do igc=1,ngaus                                                 3d2s16
          do igd=1,ngaus                                                9d28s16
           jc=jbdat+igc                                                  3d2s16
           jc2=jc+ngaus                                                    2d19s10
           jc3=jc2+ngaus                                                   2d19s10
           jc4=jc3+ngaus                                                   2d19s10
           jc5=jc4+ngaus                                                   2d19s10
           jc6=jc5+ngaus                                                   2d19s10
           jc7=jc6+ngaus                                                   2d19s10
           jc8=jc7+ngaus                                                   5d11s10
           jd=jbdat+igd                                                 3d2s16
           jd2=jd+ngaus                                                    2d19s10
           jd3=jd2+ngaus                                                   2d19s10
           jd4=jd3+ngaus                                                   2d19s10
           jd5=jd4+ngaus                                                   2d19s10
           jd6=jd5+ngaus                                                   2d19s10
           jd7=jd6+ngaus                                                   2d19s10
           jd8=jd7+ngaus                                                   5d11s10
           nszc=2*ibc(jc)+1                                              2d7s12
           nszc0=nszc                                                   5d31s22
           nszc=nszc*ncomp                                               8d24s15
           nszca=nszc                                                    2d7s12
           nszca0=nszc0                                                 5d31s22
           nszd=2*ibc(jd)+1                                              2d7s12
           nszd0=nszd                                                   5d31s22
           nszd=nszd*ncomp                                               8d24s15
           nszda=nszd                                                    2d7s12
           nszda0=nszd0                                                 5d31s22
           if(iapair(1,ibc(jc8)).gt.0)then                                5d12s10
            nszca=nszca*2                                                 5d12s10
            nszca0=nszca0*2                                             5d31s22
           end if                                                         5d12s10
           if(iapair(1,ibc(jd8)).gt.0)then                                5d12s10
            nszda=nszda*2                                                 5d12s10
            nszda0=nszda0*2                                             5d31s22
           end if                                                         5d12s10
           nwtmp=ncsz0*nszca0*nszda0                                    5d31s22
           nweri=ncsz0*nszc0*nszd0                                      5d31s22
           itmp=ibcoff                                                    5d12s10
           ibcoff=itmp+nwtmp*nerimat                                    3d7s23
           memneed(1:4)='itmp'
           call enough('paraerid2c.  6',bc,ibc)
           do iz=itmp,ibcoff-1                                          5d31s22
            bc(iz)=0d0                                                  5d31s22
           end do                                                       5d31s22
           nwds=nsza*nszb*nszc*nszd
           ldeb=idoit.eq.1                                              4d21s16
           if(ibc(jb8).eq.iatom.and.ibc(jd8).eq.iatom.                  3d3s23
     $          and.ipassb.eq.1)then                                    3d3s23
            if(idorel.eq.0)then                                         3d7s23
             call erird(ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3), 3d6s17
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       3d6s17
     $           ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),       3d6s17
     $           ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),       3d6s17
     $           ieri,0,sigb,0,0,0,idx,idy,idz,0,0,0,idx,idy,idz,bc,ibc)11d14s22
            else                                                        3d7s23
             call restind                                               3d7s23
     $            (ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),     3d7s23
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       3d6s17
     $           ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),       3d6s17
     $           ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),       3d6s17
     $           ieri,0,sigb,idorel,ascale,nmat,0,0,0,idx,idy,idz,0,0,0,3d7s23
     $            idx,idy,idz,bc,ibc)                                   3d7s23
            end if                                                      3d7s23
           else                                                          2d23s16
            ieri=ibcoff
            ibcoff=ieri+nwds                                            4d21s16
            call enough('paraerid2c.  7',bc,ibc)
            do k=0,nwds-1
             bc(ieri+k)=0d0
            end do
           end if
           if(idoit.eq.1)then
            write(6,*)('output from eri: '),ieri,nszd,nszc,ncsz
            call prntm2(bc(ieri),nszd*nszc,ncsz,nszd*nszc)
            call dws_sync
            call dws_finalize
            stop
           end if
c
            jtmp=itmp                                                     10d28s20
            jeri=ieri                                                     10d28s20
            do in=1,nerimat                                               10d28s20
             do id=0,nszd0-1                                              10d28s20
              do ic=0,nszc0-1                                             10d28s20
               iad1=jtmp+ncsz0*(ic+nszca0*id)                             10d28s20
               do iab=0,ncsz0-1                                           10d28s20
                iad2=jeri+id+nszd0*(ic+nszc0*iab)                         10d28s20
                bc(iad1+iab)=bc(iad2)                                       5d12s10
               end do                                                       5d12s10
              end do                                                        5d12s10
             end do                                                         5d12s10
             jtmp=jtmp+nwtmp                                              10d28s20
             jeri=jeri+nweri                                              10d28s20
            end do                                                        10d28s20
           ibcoff=ieri                                                   2d23s16
           if(idoit.eq.1)then
            write(6,*)('transposed ')
            call prntm2(bc(itmp),ncsz,nszca*nszd,ncsz)
           end if
           if(iapair(1,ibc(jc8)).gt.0)then                                5d12s10
            cartc(1)=bc(jc5)*dfloat(isym(1,iapair(2,ibc(jc8))))           5d12s10
            cartc(2)=bc(jc6)*dfloat(isym(2,iapair(2,ibc(jc8))))           5d12s10
            cartc(3)=bc(jc7)*dfloat(isym(3,iapair(2,ibc(jc8))))           5d12s10
            sigc=1d0                                                    4d18s16
            if(ibc(jb8).eq.iatom.and.ibc(jd8).eq.iatom                  3d3s23
     $           .and.ipassb.eq.1)then                                  3d3s23
             if(idorel.eq.0)then                                        3d7s23
              call erird(                                               3d7s23
     $            ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),      3d7s23
     $            ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),      3d7s23
     $         ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3),         3d7s23
     $         ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),         3d7s23
     $         ieri,0,sigb,0,0,0,idx,idy,idz,0,0,0,idx,idy,idz,bc,ibc)  3d7s23
             else                                                       3d7s23
              call restind(                                               3d7s23
     $            ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),      3d7s23
     $            ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),      3d7s23
     $         ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3),         3d7s23
     $         ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),         3d7s23
     $         ieri,0,sigb,idorel,ascale,nmat,0,0,0,idx,idy,idz,0,0,0,  3d7s23
     $             idx,idy,idz,bc,ibc)                                  3d7s23
             end if                                                     3d7s23
            else                                                        4d18s16
             ieri=ibcoff
             ibcoff=ieri+nwds                                           4d21s16
             call enough('paraerid2c.  8',bc,ibc)
             do k=0,nwds-1
              bc(ieri+k)=0d0
             end do
            end if
            if(idoit.eq.1)then
             write(6,*)('outputb from eri: '),ieri,nszd,nszc,ncsz
             call prntm2(bc(ieri),nszd*nszc,ncsz,nszd*nszc)
            end if
             jtmp=itmp                                                     10d28s20
             jeri=ieri                                                     10d28s20
             do in=1,nerimat                                               10d28s20
              do id=0,nszd0-1                                              10d28s20
               do ic=0,nszc0-1                                             10d28s20
                iad1=jtmp+ncsz0*(ic+nszc0+nszca0*id)                    5d31s22
                do iab=0,ncsz0-1                                           10d28s20
                 iad2=jeri+id+nszd0*(ic+nszc0*iab)                         10d28s20
                 bc(iad1+iab)=bc(iad2)                                       5d12s10
                end do                                                       5d12s10
               end do                                                        5d12s10
              end do                                                         5d12s10
              jtmp=jtmp+nwtmp                                              10d28s20
              jeri=jeri+nweri                                              10d28s20
             end do                                                        10d28s20
            ibcoff=ieri                                                 4d21s16
           end if                                                         5d12s10
           if(iapair(1,ibc(jd8)).gt.0)then                               1d10s11
            cartd(1)=bc(jd5)*dfloat(isym(1,iapair(2,ibc(jd8))))           5d12s10
            cartd(2)=bc(jd6)*dfloat(isym(2,iapair(2,ibc(jd8))))           5d12s10
            cartd(3)=bc(jd7)*dfloat(isym(3,iapair(2,ibc(jd8))))           5d12s10
            sigd=1d0                                                    4d18s16
            if(ibc(jb8).eq.iatom.and.ibc(jd8).eq.iatom                  3d3s23
     $           .and.ipassb.eq.2)then                                  3d3s23
             if(idorel.eq.0)then                                        3d7s23
              call erird(
     $            ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),      3d7s23
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       3d7s23
     $         ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),         3d7s23
     $         ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3),         3d7s23
     $         ieri,0,sigbm,0,0,0,idx,idy,idz,0,0,0,idx,idy,idz,bc,ibc) 3d7s23
             else                                                       3d7s23
              call restind(
     $            ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),      3d7s23
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       3d7s23
     $         ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),         3d7s23
     $         ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3),         3d7s23
     $         ieri,0,sigbm,idorel,ascale,nmat,0,0,0,idx,idy,idz,0,0,0, 3d7s23
     $             idx,idy,idz,bc,ibc)                                  3d7s23
             end if                                                     3d7s23
            else                                                        4d18s16
             ieri=ibcoff
             ibcoff=ieri+nwds                                           4d21s16
             call enough('paraerid2c.  9',bc,ibc)
             do k=0,nwds-1
              bc(ieri+k)=0d0
             end do
            end if
            if(idoit.eq.1)then
             write(6,*)('outputc from eri: '),ieri,nszd,nszc,ncsz
             call prntm2(bc(ieri),nszd*nszc,ncsz,nszd*nszc)
            end if
             jtmp=itmp                                                     10d28s20
             jeri=ieri                                                     10d28s20
             do in=1,nerimat                                               10d28s20
              do id=0,nszd0-1                                              10d28s20
               idp=id+nszd0                                             5d31s22
               do ic=0,nszc0-1                                             10d28s20
                iad1=jtmp+ncsz0*(ic+nszca0*idp)                         5d31s22
                do iab=0,ncsz0-1                                           10d28s20
                 iad2=jeri+id+nszd0*(ic+nszc0*iab)                         10d28s20
                 bc(iad1+iab)=bc(iad2)                                       5d12s10
                end do                                                       5d12s10
               end do                                                        5d12s10
              end do                                                         5d12s10
              jtmp=jtmp+nwtmp                                              10d28s20
              jeri=jeri+nweri                                              10d28s20
             end do                                                        10d28s20
            ibcoff=ieri                                                 4d21s16
            if(iapair(1,ibc(jc8)).gt.0)then                                5d12s10
             if(ibc(jb8).eq.iatom.and.ibc(jd8).eq.iatom                 3d3s23
     $           .and.ipassb.eq.2)then                                  3d3s23
              if(idorel.eq.0)then                                       3d7s23
               call erird                                               3d7s23
     $             (ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3)     3d7s23
     $            ,ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),     3d7s23
     $         ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3),         3d7s23
     $         ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3),         3d7s23
     $         ieri,0,sigbm,0,0,0,idx,idy,idz,0,0,0,idx,idy,idz,bc,ibc) 3d7s23
              else                                                      3d7s23
               call restind                                             3d7s23
     $             (ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3)     3d7s23
     $            ,ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),     3d7s23
     $         ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3),         3d7s23
     $         ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3),         3d7s23
     $         ieri,0,sigbm,idorel,ascale,nmat,0,0,0,idx,idy,idz,0,0,0, 3d7s23
     $              idx,idy,idz,bc,ibc)                                 3d7s23
              end if                                                    3d7s23
             else                                                        4d18s16
              ieri=ibcoff
             ibcoff=ieri+nwds                                           4d21s16
             call enough('paraerid2c. 10',bc,ibc)
              do k=0,nwds-1
               bc(ieri+k)=0d0
              end do
             end if
             if(idoit.eq.1)then
              write(6,*)('outputd from eri: '),ieri,nszd,nszc,ncsz
              call prntm2(bc(ieri),nszd*nszc,ncsz,nszd*nszc)
             end if
              jtmp=itmp                                                     10d28s20
              jeri=ieri                                                     10d28s20
              do in=1,nerimat                                               10d28s20
               do id=0,nszd0-1                                              10d28s20
                idp=id+nszd0                                             5d31s22
                do ic=0,nszc0-1                                             10d28s20
                 iad1=jtmp+ncsz0*(ic+nszc0+nszca0*idp)                  5d31s22
                 do iab=0,ncsz0-1                                           10d28s20
                  iad2=jeri+id+nszd0*(ic+nszc0*iab)                         10d28s20
                  bc(iad1+iab)=bc(iad2)                                       5d12s10
                 end do                                                       5d12s10
                end do                                                        5d12s10
               end do                                                         5d12s10
               jtmp=jtmp+nwtmp                                              10d28s20
               jeri=jeri+nweri                                              10d28s20
              end do                                                        10d28s20
             ibcoff=ieri                                                4d21s16
            end if                                                         5d12s10
           end if                                                         5d12s10
           nszabc=ncsz*nszca                                             5d12s10
           nszabc0=ncsz0*nszca0                                         5d31s22
           if(idoit.eq.1)then
            write(6,*)('before sum/dif: ')
            call prntm2(bc(itmp),nszabc,nszda,nszabc)
           end if
           if(nszd.ne.nszda)then                                          5d12s10
            jtmp=itmp                                                     10d28s20
            do in=1,nerimat                                               10d28s20
             if(idoit.eq.1)then
              write(6,*)('before sum/dif: '),in
              call prntm2(bc(jtmp),nszabc0,nszda0,nszabc0)
             end if
             do id=0,nszd0-1                                              10d28s20
              it1=jtmp+nszabc0*id                                         10d28s20
              it2=jtmp+nszabc0*(id+nszd0)                                 10d28s20
              do iabc=0,nszabc0-1                                         10d28s20
               sum=srh*(bc(it1+iabc)+bc(it2+iabc))                         5d12s10
               diff=srh*(-bc(it1+iabc)+bc(it2+iabc))                       5d12s10
               bc(it1+iabc)=sum                                            5d12s10
               bc(it2+iabc)=diff                                           5d12s10
              end do                                                       5d12s10
             end do                                                        5d12s10
             if(idoit.eq.1)then
              write(6,*)('after sum/dif: '),in
              call prntm2(bc(jtmp),nszabc0,nszda0,nszabc0)
             end if
             jtmp=jtmp+nwtmp                                              10d28s20
            end do                                                        10d28s20
           end if                                                         5d12s10
           if(nszc.ne.nszca)then                                          5d12s10
            jtmp=itmp                                                     10d28s20
            do in=1,nerimat                                               10d28s20
             do id=0,nszda0-1                                             10d28s20
              do ic=0,nszc0-1                                             10d28s20
               it1=jtmp+ncsz0*(ic+nszca0*id)                              10d28s20
               it2=jtmp+ncsz0*(ic+nszc0+nszca0*id)                        10d28s20
               do iab=0,ncsz0-1                                           10d28s20
                sum=srh*(bc(it1+iab)+bc(it2+iab))                          5d12s10
                diff=srh*(-bc(it1+iab)+bc(it2+iab))                        5d12s10
                bc(it1+iab)=sum                                            5d12s10
                bc(it2+iab)=diff                                           5d12s10
               end do                                                      5d12s10
              end do                                                       5d12s10
             end do                                                        5d12s10
             if(idoit.eq.1)then
              write(6,*)('after sum/dif: '),in
              call prntm2(bc(jtmp),nszabc0,nszda0,nszabc0)
             end if
             jtmp=jtmp+nwtmp                                              10d28s20
            end do                                                        10d28s20
           end if
           do in=0,nerimat-1                                              10d28s20
            do id=1,nszda0                                                8d26s15
             idd=id+ibc(jd4)                                               5d12s10
             isd=isstor(idd)                                               5d12s10
             do ic=1,nszca0                                               8d26s15
              icc=ic+ibc(jc4)                                              5d12s10
              isc=isstor(icc)                                              5d12s10
              iad1=ieraw(isc,isd)-1                                      9d20s16
     $         +ncsz0*(ibstor(icc)-1+nbasisp(isc)*(ibstor(idd)-1        3d7s23
     $             +nbasisp(isd)*in))                                   3d7s23
              iad2=itmp-1+ncsz0*(ic-1+nszca0*(id-1+nszda0*in))          3d7s23
              do iab=1,ncsz0                                            3d7s23
               bc(iad1+iab)=bc(iad2+iab)                                   5d12s10
              end do                                                       5d12s10
             end do                                                     3d7s23
            end do                                                        5d12s10
           end do                                                         5d12s10
           ibcoff=itmp                                                  3d7s23
          end do                                                           3d4s10
         end do
c
c     at this point, we have all of cd ints, and furthermore they have  1d10s11
c     been symmetrized (turned into so's) and stored at eraw            1d10s11
c     so transform cd from so's to mo's                                 1d10s11
c     1d10s11
         do isd=1,nsymb                                                  5d12s10
          if(nocc(isd).gt.0)then                                        6d2s10
           do isc=1,nsymb                                                 5d12s10
            if(nbasdws(isc).gt.0)then                                   9d28s16
             nabc=ncsz0*nbasisp(isc)                                    3d7s23
             ihalf=ibcoff                                                  5d12s10
             ibcoff=ihalf+nabc*nocc(isd)                                   5d12s10
             memneed(1:5)='ihalf'
             call enough('paraerid2c. 11',bc,ibc)
             if(idoit.eq.1)then
              write(6,*)('vectors: '),isc,isd
              call prntm2(bc(iorb(isd)),nbasisp(isd)*ncomp,nocc(isd),   4d7s22
     $             nbasisp(isd)*ncomp)                                  4d7s22
              write(6,*)('integrals ')
              nnz=0d0
              do ixd=1,nabc*nbasisp(isd)*ncomp                          4d7s22
               if(abs(bc(ieraw(isc,isd)+ixd-1)).gt.1d-6)nnz=nnz+1       9d20s16
              end do
              write(6,*)('nonzero: '),nnz,(' vs '),ndoit
              call prntm2(bc(ieraw(isc,isd)),nabc,nbasisp(isd)*ncomp,   4d7s22
     $             nabc)
             end if
             nbpn=nbasisp(isd)*ncomp                                    5d31s22
             if(min(nabc,nbpn,nocc(isd)).gt.0)then                      5d31s22
              nabc2=ncsz0*nocc(isd)                                     5d31s22
              nbpm=nbasisp(isc)*ncomp                                    4d30s18
              ncd=nbasdws(isc)*nocc(isd)                                5d31s22
              jeraw=ieraw(isc,isd)                                      10d28s20
              do in=1,nerimat                                           10d28s20
               if(in.eq.1.or.in.eq.3)then                               10d28s20
                jorb=iorb(isd)                                          10d28s20
               else                                                     10d28s20
                jorb=iorb(isd)+nbasisp(isd)                             10d28s20
               end if                                                   10d28s20
               call dgemm('n','n',nabc,nocc(isd),nbasisp(isd),1d0,      5d31s22
     $               bc(jeraw),nabc,bc(jorb),nbpn,0d0,bc(ihalf),nabc,    10d28s20
     d' paraeri.  1')
               if(idoit.eq.1)then
                write(6,*)('half transformed matrix '),in
                call prntm2(bc(ihalf),nabc,nocc(isd),nabc)              5d31s22
               end if
               do id=0,nocc(isd)-1                                      5d31s22
                do ic=0,nbasisp(isc)-1                                  10d28s20
                 iad1=ihalf+ncsz0*(ic+nbasisp(isc)*id)                  10d28s20
                 iad2=jeraw+ncsz0*(id+nocc(isd)*ic)                     5d31s22
                 do iab=0,ncsz0-1                                       10d28s20
                  bc(iad2+iab)=bc(iad1+iab)                             10d28s20
                 end do                                                 10d28s20
                end do                                                  10d28s20
               end do                                                   10d28s20
               if(in.eq.1.or.in.eq.3)then                               10d28s20
                jorb=iorb(isc)                                          10d28s20
               else                                                     10d28s20
                jorb=iorb(isc)+nbasisp(isc)                             10d28s20
               end if                                                   10d28s20
               call dgemm('n','n',nabc2,nbasdws(isc),nbasisp(isc),1d0,  10d28s20
     $               bc(jeraw),nabc2,bc(jorb),nbpm,0d0,                  10d28s20
     $               bc(ihalf),nabc2,                                    10d28s20
     d' paraeri.  2')
               if(idoit.eq.1)then
                write(6,*)('result '),in
                call prntm2(bc(ihalf),nabc2,nbasdws(isc),nabc2)         10d28s20
               end if
               if(in.le.2)then                                          10d28s20
                jtmpab=itmpab(isc,isd)                                    10d28s20
               else                                                     10d28s20
                jtmpab=itmpab(isc,isd)+ncsza0*ncd                       10d28s20
               end if
               do icd=0,ncd-1                                           10d28s20
                iad1=ihalf+ncsz0*icd                                    10d28s20
                do ia=0,nsza0-1                                         10d28s20
                 do ib=0,nszb0-1                                        10d28s20
                  iaba=ib+nszba0*ia                                     10d28s20
                  iab=ib+nszb0*ia                                       10d28s20
                  iad2=jtmpab+icd+ncd*(iaba+ioab)                       10d28s20
                  bc(iad2)=bc(iad2)+bc(iad1+iab)                        10d28s20
                 end do                                                      5d12s10
                end do                                                    2d9s12
               end do                                                       5d12s10
               jeraw=jeraw+nabc*nbasisp(isd)                            10d28s20
              end do                                                    10d28s20
             end if
             ibcoff=ihalf                                                5d12s10
            end if                                                      6d2s10
           end do                                                         5d12s10
          end if                                                        6d2s10
         end do                                                          5d12s10
c
c     at this point we have only 2 symmetry indices because ab have not 1d10s11
c     yet been transformed to so's                                      1d10s11
c
        end do                                                          5d12s10
       end do                                                           5d12s10
       ibcoff=ieraw(1,1)                                                5d12s10
c
c     complete ao to so transformation for ab                           1d10s11
c
       if(iapair(1,ibc(jb8)).gt.0)then                                  5d12s10
        do isd=1,nsymb                                                  5d12s10
         do isc=1,nsymb                                                 9d27s16
          jtmpab=itmpab(isc,isd)                                         10d28s20
          ncd=nbasdws(isc)*nocc(isd)                                    4d7s22
          do in=1,ntmpmat                                               3d7s23
           do ib=0,nszb0-1                                              3d7s23
            do ia=0,nszaa0-1                                            3d7s23
             it1=jtmpab+ncd*(ib+nszba0*ia)                              10d28s20
             it2=jtmpab+ncd*(ib+nszb0+nszba0*ia)                        10d28s20
             do icd=0,ncd-1                                             3d7s23
              sum=srh*(bc(it1+icd)+bc(it2+icd))                          5d12s10
              diff=srh*(-bc(it1+icd)+bc(it2+icd))                        5d12s10
              bc(it1+icd)=sum                                            5d12s10
              bc(it2+icd)=diff                                           5d12s10
             end do                                                       5d12s10
            end do                                                        5d12s10
           end do                                                        2d23s16
           jtmpab=jtmpab+ncd*nszaa0*nszba0                              10d28s20
          end do                                                        3d7s23
         end do                                                         5d12s10
        end do                                                          5d12s10
       end if                                                           5d12s10
       if(iapair(1,ibc(ja8)).gt.0)then                                  5d12s10
        do isd=1,nsymb                                                  5d12s10
         do isc=1,nsymb                                                 5d12s10
          ncdb=nbasdws(isc)*nocc(isd)*nszba0                            3d7s23
          jtmpab=itmpab(isc,isd)                                        10d28s20
          do in=1,ntmpmat                                               10d28s20
           do ia=0,nsza0-1                                              3d7s23
            it1=jtmpab+ncdb*ia                                          10d28s20
            it2=jtmpab+ncdb*(ia+nsza0)                                  10d28s20
            do iacd=0,ncdb-1                                            3d7s23
             sum=srh*(bc(it1+iacd)+bc(it2+iacd))                         5d12s10
             diff=srh*(-bc(it1+iacd)+bc(it2+iacd))                       5d12s10
             bc(it1+iacd)=sum                                            5d12s10
             bc(it2+iacd)=diff                                           5d12s10
            end do                                                       5d12s10
           end do                                                        5d12s10
           jtmpab=jtmpab+ncdb*nszaa0                                    10d28s20
          end do                                                        3d7s23
         end do                                                         5d12s10
        end do                                                          5d12s10
       end if                                                           5d12s10
c
c     final storage                                                     1d10s11
c
       do isd=1,nsymb                                                   5d12s10
        do isc=1,nsymb                                                  5d12s10
         ieraw(isc,isd)=nbasdws(isc)*nocc(isd)                          4d7s22
        end do                                                          5d12s10
       end do                                                           5d12s10
       do ia=1,nszaa0                                                   8d26s15
        iaa=ia+ibc(ja4)                                                 2d3s12
        isa=isstor(iaa)                                                 2d23s16
        do ib=1,nszba0                                                  8d26s15
         ibb=ib+ibc(jb4)                                                2d3s12
         isb=isstor(ibb)                                                2d3s12
          isab=multh(isa,isb)                                            9d20s16
          do isd=1,nsymb                                                 5d12s10
           isc=multh(isd,isab)                                          9d28s16
           ii=iptoh(isd,isc,isb)                                        9d26s16
           ncd=ieraw(isc,isd)                                           9d26s16
           if(ncd.gt.0)then                                             9d20s16
            if(ii.gt.0)then                                             9d20s16
             jj=itmpab(isc,isd)+ncd*(ib-1+nszba0*(ia-1))                3d7s23
             do in=0,ntmpmat-1
              do ix=0,ncd-1                                             3d7s23
               bc(ihcol(ii))=bc(jj+ix)                                  3d7s23
               bc(jj+ix)=0d0                                            3d7s23
               ihcol(ii)=ihcol(ii)+1                                       5d12s10
              end do
              jj=jj+ncd*nszba0*nszaa0                                     10d30s20
             end do                                                     3d7s23
            end if
           end if                                                        1d31s12
          end do                                                         5d12s10
        end do                                                          5d12s10
       end do                                                           5d12s10
       do isd=1,nsymb
        do isc=1,nsymb
         sz=0d0
         do j=0,ncsza*nbasdws(isc)*nocc(isd)-1                          4d7s22
          sz=sz+bc(itmpab(isc,isd)+j)**2
         end do
         iad=iszz+isc-1+nsymb*(isd-1)
         bc(iad)=bc(iad)+sz
        end do
       end do
       ibcoff=itmpab(1,1)                                               5d12s10
      end do                                                            2d19s10
      call dws_sync                                                     2d22s10
      ihcol(1)=ihmat                                                    1d10s11
      do i=2,nhmt                                                       9d26s16
       ihcol(i)=ihcol(i-1)+nhcol(i-1)                                   1d10s11
      end do                                                            5d12s10
      if(idwsdeb.gt.10)then                                             3d16s12
      write(6,*)('we have this procs part of hmat ')                    3d18s10
      write(6,*)('what we''ve got: ')                                   1d10s11
      memtop=0
      memtot=0
      do isb=1,nsymb                                                    1d10s11
       do isc=1,nsymb                                                   9d27s16
        isbc=multh(isc,isb)                                             9d26s16
        do isd=1,nsymb                                                  1d10s11
         isa=multh(isd,isbc)                                            9d28s16
         ii=iptoh(isd,isc,isb)                                          9d27s16
         ncd=nbasdws(isc)*nocc(isd)                                     4d7s22
         write(6,*)('ii,ncd: '),ii,ncd
         if(ii.gt.0.and.ncd.gt.0)then
          if(nhcol(ii).gt.0)then                                        1d10s11
           write(6,*)('integrals for symmetry '),isa,isb,isc,isd
           mcol=nhcol(ii)/ncd                                           1d10s11
           memtop=max(memtop,ihcol(ii)+ncd*mcol-1)                      3d19s12
           memtot=memtot+nhcol(ii)
           sz=0d0
           szx=-1d20
           do iq=0,ncd*mcol-1
            sz=sz+bc(ihcol(ii)+iq)**2
            if(abs(bc(ihcol(ii)+iq)).gt.szx)then
             szx=abs(bc(ihcol(ii)+iq))
             iszx=iq+1
            end if
           end do
           sz=sqrt(sz/dfloat(ncd*mcol))
           write(6,*)('rms size: '),sz,ihcol(ii),ihcol(ii)+nhcol(ii)
           write(6,*)('max '),szx,(' at '),iszx
           call prntm2(bc(ihcol(ii)),ncd,mcol,ncd)                      1d10s11
          end if                                                        1d10s11
         end if                                                         1d10s11
        end do                                                          1d10s11
       end do                                                           1d10s11
      end do                                                            1d10s11
      write(6,*)('memtop '),memtop
      write(6,*)('memtot '),memtot
      end if                                                            3d16s12
      ibcoff=ibcexit                                                    3d1s16
      return
      end                                                               2d19s10
