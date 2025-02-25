c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine paraeridd                                              6d11s24
     $     (natom,ngaus,ibdat,nbasis,ihmat,iorb,nocc,ipair,             6d11s24
     $     nhcolt,isym,iapair,ibstor,isstor,multh,iptoh,iter,idwsdeb,   8d24s15
     $     idorel,ascale,nbasisp,iaddr,naddr,ixyz,iatom,                6d11s24
     $     idersign,d4v,ifull4x,bc,ibc)                                 6d14s24
c
c     compute partial (rs|cd] with rs in ao basis and cd in mo basis.
c     wrt ixyz coordinate of iatom. We are only getting the derivative  6d11s24
c     of the ao integrals.                                              6d11s24
c
      implicit real*8 (a-h,o-z)
      character*4 ssss                                                  12d7s18
      include "common.hf"
      include "common.store"
      include "common.basis"
      include "common.print"                                            1d6s20
      logical lsame                                                     5d31s24
      integer*8 ibstor,isstor                                           5d10s10
      dimension iorb(*),isym(3,*),iapair(3,*),ibstor(*),isstor(*),      5d12s10
     $     carta(3),cartb(3),cartc(3),cartd(3),ieraw(8,8),itmpab(8,8),  5d12s10
     $     multh(8,8),iptoh(8,8,8),nocc(8),neraw(8,8),nueraw(8,8),      5d29s18
     $     nbasisp(8),fmulx(4),nduse(8),nduse2(8),nan(8),nbn(8),        4d9s24
     $     ieraw2(8,8,8),itmpab2(8,8,8),iptap(8),iptbp(8)               4d11s24
      common/timerocm/tovr,telapo(15)                                   4d26s18
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      common/drsigncm/drsign                                            8d20s24
      if(ifull4x.eq.1)then                                              6d11s21
       do isb=1,nsymb                                                   6d11s21
        nduse(isb)=nocc(isb)                                            6d11s21
        nduse2(isb)=nocc(isb)                                            6d11s21
       end do                                                           6d11s21
      else if(ifull4x.eq.2)then                                         8d18s23
       do isb=1,nsymb                                                    6d11s24
        nduse(isb)=nocc(isb)                                             6d11s24
        nduse2(isb)=nbasdws(isb)                                         6d11s24
       end do                                                            6d11s24
      else if(ifull4x.eq.3)then                                         8d18s23
       do isb=1,nsymb                                                   6d11s21
        nduse(isb)=nbasdws(isb)                                         6d11s21
        nduse2(isb)=nbasdws(isb)                                        8d18s23
       end do                                                           6d11s21
      else                                                              8d18s23
       write(6,*)('unknown ifull4x code in paraeri: '),ifull4x          8d18s23
       stop 'paraeri'                                                   8d18s23
      end if                                                            6d11s21
      idx=0                                                             2d1s16
      idy=0                                                             2d1s16
      idz=0                                                             2d1s16
      if(ixyz.eq.1)then                                                 2d1s16
       idx=1                                                            6d11s24
      else if(ixyz.eq.2)then                                            2d1s16
       idy=1                                                            6d11s24
      else                                                              2d1s16
       idz=1                                                            6d11s24
      end if                                                            2d1s16
      fmulx(1)=1d0                                                      3d13s20
      fmulx(2)=ascale                                                   3d13s20
      fmulx(3)=ascale                                                   3d13s20
      fmulx(4)=ascale*ascale                                            3d13s20
      ntmpmat=2                                                         10d28s20
      if(idorel.eq.0.or.idorel.eq.1)then                                10d28s20
       nerimat=1                                                        10d28s20
       ntmpmat=1                                                        10d28s20
      else if(idorel.eq.2.or.idorel.eq.-2.or.idorel.eq.-4)then          10d28s20
       nerimat=3                                                        10d28s20
      else                                                              10d28s20
       nerimat=4                                                        10d28s20
      end if                                                            10d28s20
      ibcoffo=ibcoff                                                    2d19s10
      if(idorel.eq.0)then                                               8d24s15
       ncomp=1                                                          8d24s15
      else                                                              8d24s15
       ncomp=2                                                          8d24s15
      end if                                                            8d24s15
      nbasallp=nbasisp(1)                                               1d30s18
      do isb=2,nsymb                                                    1d30s18
       nbasallp=nbasallp+nbasisp(isb)                                   1d30s18
      end do                                                            1d30s18
      nbasallc=nbasallp*ncomp                                           1d30s18
      if(idwsdeb.ne.0)then                                              6d11s24
      write(6,*)('in paraeridd '),ibcoff,iter,idorel,ixyz,iatom,idersign6d11s24
      write(6,*)('natom etc '),natom,ngaus,ibdat,nbasis,ihmat,iorb(1),
     $     nduse(1),ipair,nhcolt,(isym(j,1),j=1,3),iter,idwsdeb,idorel, 11d29s23
     $     ascale                                                       11d29s23
      end if
      srh=sqrt(0.5d0)                                                   5d12s10
      nbasx=nbasisp(1)*ncomp                                            8d24s15
      noccx=nduse(1)                                                    6d11s21
      do isb=2,nsymb                                                    5d12s10
       nbasx=max(nbasx,nbasisp(isb)*ncomp)                              8d24s15
       noccx=max(noccx,nduse(isb))                                      6d11s21
      end do                                                            5d12s10
c
c     build shell pair order
c
      nn=(ngaus*(ngaus+1))/2                                            2d19s10
      nn2=nn*2                                                          2d19s10
      ipair=ibcoff                                                      2d19s10
      if(iprtr(18).ne.0)then                                            5d4s20
       write(6,*)('nn,nn2,ibcoff '),nn,nn2,ibcoff
      end if                                                            7d27s16
      ibcoff=ipair+nn*2                                                 3d4s10
      i12=ibcoff                                                        2d19s10
      ibcoff=i12+nn*4                                                   2d19s10
      memneed(1:3)='i12'
      call enough('paraeridd.  1',bc,ibc)
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
c
c     sort, so we can assign most expensive blocks first
c
      if(iprtr(18).ne.0)then                                            5d4s20
       write(6,*)('calling ihpsort next ')
      end if
      call idsortdws(ibc(i12),ibc(is12),nn8)                            1d18s23
      ii1=i12+nn                                                        2d19s10
      ii2=i12+nn2                                                       2d19s10
      jpair=ipair                                                       2d19s10
c
c     order shell pair indices
c
      do i=1,nn                                                         2d19s10
       j=ibc(is12+i-1)-1                                                2d19s10
    1  format(i5,i8,2i5)                                                2d19s10
       ibc(jpair)=ibc(ii1+j)                                            2d19s10
       jpair=jpair+1                                                    2d19s10
       ibc(jpair)=ibc(ii2+j)                                            2d19s10
       jpair=jpair+1                                                    2d19s10
      end do                                                            2d19s10
      do i=1,nsdlkh                                                     5d12s10
       nhcol(i)=0                                                       5d12s10
       ihcol(i)=0                                                       5d12s10
      end do                                                            5d12s10
c
c     compute storage required on this proc
c
      idws=0                                                            1d10s11
      nsum1=0
      nsum2=0
      do i=1+mynowprog,nn,mynprocg                                      3d18s10
       jpair=ipair+2*(i-1)                                              2d19s10
       la=ibc(jbdat+ibc(jpair))                                         3d18s10
       na=ncomp*(2*la+1)                                                8d24s15
       lb=ibc(jbdat+ibc(jpair+1))                                       3d18s10
       nb=ncomp*(2*lb+1)                                                8d24s15
        naa=na                                                          5d12s10
        nba=nb                                                          5d12s10
        if(iapair(1,ibc(jbdat+ibc(jpair)+ngaus7)).gt.0)then             5d12s10
         naa=naa*2                                                      5d12s10
        end if                                                          5d12s10
        if(iapair(1,ibc(jbdat+ibc(jpair+1)+ngaus7)).gt.0)then           5d12s10
         nba=nba*2                                                      5d12s10
        end if                                                          5d12s10
        nba0=nba/ncomp                                                  9d18s15
        naa0=naa/ncomp                                                  9d18s15
        do ib=1,nba0                                                    9d18s15
         ibb=ib+ibc(jbdat+ibc(jpair+1)+ngaus3)                          6d4s10
         isb=isstor(ibb)                                                5d12s10
         do ia=1,naa0                                                   9d18s15
          iaa=ia+ibc(jbdat+ibc(jpair)+ngaus3)                           6d4s10
          isa=isstor(iaa)                                               5d12s10
          if(ibb.le.iaa)then                                            1d30s12
           isab=multh(isa,isb)                                           5d12s10
           do isd=1,nsymb                                                5d12s10
            isc=multh(isd,isab)                                          5d12s10
            nadd=nduse2(isc)*nduse(isd)                                 8d18s23
            nadd=nadd*ntmpmat                                           10d30s20
            if(nadd.gt.0)then                                           1d30s12
             if(iptoh(isd,isc,isb).gt.0)then                            1d30s12
              nhcol(iptoh(isd,isc,isb))=nhcol(iptoh(isd,isc,isb))+nadd  1d30s12
              nsum1=nsum1+1
             else if(iptoh(isd,isc,isa).gt.0)then                       1d30s12
              nhcol(iptoh(isd,isc,isa))=nhcol(iptoh(isd,isc,isa))+nadd  1d30s12
              nsum1=nsum1+1
             end if                                                     1d30s12
            end if                                                      1d30s12
           end do                                                        5d12s10
          end if                                                        5d12s10
         end do                                                         5d12s10
        end do                                                          5d12s10
      end do                                                            3d18s10
      if(iprtr(18).ne.0)then                                            5d4s20
       write(6,*)('nsum1,nsum2 '),nsum1,nsum2
      end if
      nhcolt=0                                                          5d12s10
      do i=1,nsdlkh                                                     5d12s10
 3315  format(i5,5x,4i2,5x,i8)                                          5d12s10
       nhcolt=nhcolt+nhcol(i)                                           5d12s10
      end do                                                            5d12s10
      ihmat=ibcoff                                                      3d4s10
      ibcoff=ihmat+nhcolt                                               5d12s10
      memneed(1:4)='hmat'
      call enough('paraeridd.  2',bc,ibc)
      ihcol(1)=ihmat                                                    5d12s10
      do i=2,nsdlkh                                                     5d12s10
       ihcol(i)=ihcol(i-1)+nhcol(i-1)                                   5d12s10
      end do                                                            5d12s10
      jhmat=ihmat                                                       3d4s10
      idoit=10
      jdoit=100
      ndoit=0
      nxtot=0
      do i=1+mynowprog,nn,mynprocg                                      3d4s10
       jpair=ipair+2*(i-1)                                              2d19s10
       if(iprtr(18).ne.0)then                                           5d4s20
        write(6,2)i,ibc(jpair),ibc(jpair+1),ibc(jpair+2),tcur            3d14s12
       end if
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
       nsza0=nsza
       nsza=nsza*ncomp                                                  8d24s15
       nszaa=nsza                                                       5d12s10
       nszaa0=nsza0                                                     5d4s20
       if(iapair(1,ibc(ja8)).gt.0)then                                  5d11s10
        nszaa=nszaa*2                                                   5d12s10
        nszaa0=nszaa0*2                                                 5d4s20
       end if                                                           5d11s10
       jb=jbdat+ibc(jpair+1)                                            2d19s10
       lsame=ja.eq.jb                                                   5d31s24
       jb2=jb+ngaus                                                     2d19s10
       jb3=jb2+ngaus                                                    2d19s10
       jb4=jb3+ngaus                                                    2d19s10
       jb5=jb4+ngaus                                                    2d19s10
       jb6=jb5+ngaus                                                    2d19s10
       jb7=jb6+ngaus                                                    2d19s10
       jb8=jb7+ngaus                                                    5d11s10
       nszb=2*ibc(jb)+1                                                 5d11s10
       nszb0=nszb
       nszb=nszb*ncomp                                                  8d24s15
       nszba=nszb                                                       5d12s10
       nszba0=nszb0                                                     5d4s20
 3352  format(6i5)
       if(iapair(1,ibc(jb8)).gt.0)then                                  5d11s10
        nszba=nszba*2                                                   5d12s10
        nszba0=nszba0*2                                                 5d4s20
       end if                                                           5d11s10
       do iss=1,nsymb                                                   4d10s24
        nan(iss)=0                                                      4d10s24
        nbn(iss)=0                                                      4d10s24
       end do                                                           4d9s24
       ipta=ibcoff                                                      4d9s24
       iptb=ipta+nszaa0                                                 4d9s24
       ibcoff=iptb+nszba0                                               4d9s24
       call enough('paraeridd.pta',bc,ibc)                                4d9s24
       jpta=ipta-1                                                      4d9s24
       jptb=iptb-1                                                      4d9s24
       do ia=1,nszaa0                                                   4d9s24
        iaa=ia+ibc(ja4)                                                 4d9s24
        isa=isstor(iaa)                                                 4d9s24
        ibc(jpta+ia)=nan(isa)                                           4d9s24
        nan(isa)=nan(isa)+1                                             4d9s24
       end do                                                           4d9s24
       do ib=1,nszba0                                                   4d9s24
        ibb=ib+ibc(jb4)                                                 4d9s24
        isb=isstor(ibb)                                                 4d9s24
        ibc(jptb+ib)=nbn(isb)                                           4d9s24
        nbn(isb)=nbn(isb)+1                                             4d9s24
       end do                                                           4d9s24
       do isa=1,nsymb                                                   4d10s24
        iptap(isa)=ibcoff                                               4d10s24
        iptbp(isa)=iptap(isa)+nan(isa)                                  4d10s24
        ibcoff=iptbp(isa)+nbn(isa)                                      4d10s24
       end do                                                           4d10s24
       do ia=1,nszaa0                                                   4d10s24
        iaa=ia+ibc(ja4)                                                 4d10s24
        isa=isstor(iaa)                                                 4d10s24
        ibc(iptap(isa)+ibc(jpta+ia))=ibstor(iaa)-1                      4d10s24
       end do                                                           4d10s24
       do ib=1,nszba0                                                   4d10s24
        ibb=ib+ibc(jb4)                                                 4d10s24
        isb=isstor(ibb)                                                 4d10s24
        ibc(iptbp(isb)+ibc(jptb+ib))=ibstor(ibb)-1                      4d10s24
       end do                                                           4d10s24
       iabstrt=ibcoff                                                   4d9s24
       do isa=1,nsymb                                                   4d9s24
        do isb=1,nsymb                                                  4d9s24
         isab=multh(isa,isb)                                            4d9s24
         if(min(nan(isa),nbn(isb)).gt.0)then                            4d9s24
          nab=nan(isa)*nbn(isb)                                          4d9s24
          do isc=1,nsymb                                                 4d9s24
           isd=multh(isab,isc)                                           4d9s24
           itmpab2(isa,isb,isc)=ibcoff                                    4d9s24
           ibcoff=itmpab2(isa,isb,isc)+nab*nduse2(isc)*nduse(isd)*      4d9s24
     $          ntmpmat                                                 4d9s24
          end do                                                        4d9s24
         end if                                                         4d9s24
        end do                                                          4d9s24
       end do                                                           4d9s24
       call enough('paraeridd.tmpab2',bc,ibc)                             4d9s24
       do iz=iabstrt,ibcoff-1                                           4d9s24
        bc(iz)=0d0                                                      4d9s24
       end do                                                           4d9s24
       do isa=1,nsymb                                                   4d9s24
        do isb=1,nsymb                                                  4d9s24
         isab=multh(isa,isb)                                            4d9s24
         if(min(nan(isa),nbn(isb)).gt.0)then                            4d9s24
          nab=nan(isa)*nbn(isb)                                          4d9s24
          do isc=1,nsymb                                                 4d9s24
           isd=multh(isab,isc)                                           4d9s24
           ieraw2(isa,isb,isc)=ibcoff                                    4d9s24
           ibcoff=ieraw2(isa,isb,isc)+nab*nbasisp(isc)*nbasisp(isd)*    4d9s24
     $          nerimat                                                 4d9s24
          end do                                                        4d9s24
         end if                                                         4d9s24
        end do                                                          4d9s24
       end do                                                           4d9s24
       ncsz=nsza*nszb                                                   5d11s10
       ncsz0=nsza0*nszb0                                                10d28s20
       ncsza0=nszaa0*nszba0                                             10d28s20
       ncsza=nszaa*nszba                                                1d10s11
       iptra=ibcoff                                                     6d11s24
       iptrb=iptra+nszaa                                                6d11s24
       ibcoff=iptrb+nszba                                               7d16s24
       memneed(1:5)='tmpab'
       call enough('paraeridd.  3',bc,ibc)
       npassb=nszba/nszb                                                2d29s24
       npassa=nszaa/nsza                                                2d29s24
       do j=1,nn                                                        3d4s10
        jpajr=ipair+2*(j-1)                                             3d4s10
        jc=jbdat+ibc(jpajr)                                             2d19s10
        jjjc=ibc(jpajr)
        jc2=jc+ngaus                                                    2d19s10
        jc3=jc2+ngaus                                                   2d19s10
        jc4=jc3+ngaus                                                   2d19s10
        jc5=jc4+ngaus                                                   2d19s10
        jc6=jc5+ngaus                                                   2d19s10
        jc7=jc6+ngaus                                                   2d19s10
        jc8=jc7+ngaus                                                   5d11s10
        jd=jbdat+ibc(jpajr+1)                                           2d19s10
        jjjd=ibc(jpajr+1)
        jd2=jd+ngaus                                                    2d19s10
        jd3=jd2+ngaus                                                   2d19s10
        jd4=jd3+ngaus                                                   2d19s10
        jd5=jd4+ngaus                                                   2d19s10
        jd6=jd5+ngaus                                                   2d19s10
        jd7=jd6+ngaus                                                   2d19s10
        jd8=jd7+ngaus                                                   5d11s10
        nszc=2*ibc(jc)+1                                                2d7s12
        nszc0=nszc                                                      3d13s20
        nszc=nszc*ncomp                                                 8d24s15
        nszca=nszc                                                      2d7s12
        nszca0=nszc0                                                    5d4s20
        nszd=2*ibc(jd)+1                                                2d7s12
        nszd0=nszd                                                      3d13s20
        nszd=nszd*ncomp                                                 8d24s15
        nszda=nszd                                                      2d7s12
        nszda0=nszd0                                                    5d4s20
        npassc=1                                                        4d9s24
        npassd=1                                                        4d9s24
        if(iapair(1,ibc(jc8)).gt.0)then                                 5d12s10
         nszca=nszca*2                                                  5d12s10
         nszca0=nszca0*2
         npassc=2                                                       4d9s24
        end if                                                          5d12s10
        if(iapair(1,ibc(jd8)).gt.0)then                                 5d12s10
         nszda=nszda*2                                                  5d12s10
         nszda0=nszda0*2                                                5d4s20
         npassd=2                                                       4d9s24
        end if                                                          5d12s10
        nwtmp=nszaa0*nszba0*nszca0*nszda0                               4d9s24
        nweri=nsza0*nszb0*nszc0*nszd0                                   4d9s24
        itmp=ibcoff                                                     4d9s24
        ibcoff=itmp+nwtmp*nerimat                                       4d9s24
        call enough('paraeridd.tmpn',bc,ibc)                              4d9s24
        do iz=itmp,ibcoff-1                                             6d11s24
         bc(iz)=0d0                                                     6d11s24
        end do                                                          6d11s24
        do ipassd=1,npassd                                              4d9s24
         sigd=1d0                                                       6d11s24
         if(ipassd.eq.1)then                                            4d9s24
          cartd(1)=bc(jd5)                                              4d9s24
          cartd(2)=bc(jd6)                                              4d9s24
          cartd(3)=bc(jd7)                                              4d9s24
          iod=0                                                         4d9s24
         else                                                           4d9s24
          cartd(1)=bc(jd5)*dfloat(isym(1,iapair(2,ibc(jd8))))
          cartd(2)=bc(jd6)*dfloat(isym(2,iapair(2,ibc(jd8))))
          cartd(3)=bc(jd7)*dfloat(isym(3,iapair(2,ibc(jd8))))
          if(idersign.eq.2)sigd=-1d0                                    6d11s24
          iod=nszd0                                                     4d9s24
         end if
         sigd=-drsign*sigd                                              8d21s24
         do ipassc=1,npassc                                              4d9s24
          sigc=1d0                                                      6d11s24
          if(ipassc.eq.1)then                                            4d9s24
           cartc(1)=bc(jc5)                                              4d9s24
           cartc(2)=bc(jc6)                                              4d9s24
           cartc(3)=bc(jc7)                                              4d9s24
           ioc=0                                                         4d9s24
          else                                                           4d9s24
           cartc(1)=bc(jc5)*dfloat(isym(1,iapair(2,ibc(jc8))))
           cartc(2)=bc(jc6)*dfloat(isym(2,iapair(2,ibc(jc8))))
           cartc(3)=bc(jc7)*dfloat(isym(3,iapair(2,ibc(jc8))))
           ioc=nszc0                                                     4d9s24
           if(idersign.eq.2)sigc=-1d0                                   6d11s24
          end if
          sigc=-drsign*sigc                                             8d21s24
          do ipassb=1,npassb                                               2d29s24
           sigb=1d0                                                     6d11s24
           if(ipassb.eq.1)then
            cartb(1)=bc(jb5)                                               5d12s10
            cartb(2)=bc(jb6)                                               5d12s10
            cartb(3)=bc(jb7)                                               5d12s10
            iob=0                                                          5d12s10
           else                                                            5d12s10
            cartb(1)=bc(jb5)*dfloat(isym(1,iapair(2,ibc(jb8))))            5d11s10
            cartb(2)=bc(jb6)*dfloat(isym(2,iapair(2,ibc(jb8))))            5d11s10
            cartb(3)=bc(jb7)*dfloat(isym(3,iapair(2,ibc(jb8))))            5d11s10
            iob=nszb0                                                      10d28s20
            if(idersign.eq.2)sigb=-1d0                                  6d11s24
           end if                                                          5d12s10
           sigb=-drsign*sigb                                            8d21s24
           do ipassa=1,npassa                                              2d29s24
            siga=1d0                                                    6d11s24
            if(ipassa.eq.1)then                                            5d12s10
             carta(1)=bc(ja5)                                              5d12s10
             carta(2)=bc(ja6)                                              8d15s11
             carta(3)=bc(ja7)                                              8d15s11
             ioa=0                                                         5d12s10
            else                                                           5d12s10
             carta(1)=bc(ja5)*dfloat(isym(1,iapair(2,ibc(ja8))))           5d12s10
             carta(2)=bc(ja6)*dfloat(isym(2,iapair(2,ibc(ja8))))           5d12s10
             carta(3)=bc(ja7)*dfloat(isym(3,iapair(2,ibc(ja8))))           5d12s10
             if(idersign.eq.2)siga=-1d0                                    4d26s22
             ioa=nsza0                                                     10d28s20
            end if                                                         5d12s10
            siga=-drsign*siga                                           8d21s24
            ibc00=ibcoff                                                6d11s24
            if(ibc(ja8).eq.iatom)then                                   6d11s24
             if(idorel.ne.0)then                                           3d13s20
              call restind(                                              6d11s24
     $           ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),       4d9s24
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       3d13s20
     $           ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3),       4d9s24
     $           ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3),       4d9s24
     $       ieri,jdoit,siga,idorel,ascale,nmat,idx,idy,idz,0,0,0,      6d11s24
     $            0,0,0, 0,0,0, bc,ibc)                                 6d11s24
             else                                                           3d13s20
              call erird                                                6d11s24
     $             (ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),    6d11s24
     $                 ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3), 4d9s24
     $                 ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3), 4d9s24
     $                 ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3), 4d9s24
     $         ieri,0,siga,idx,idy,idz,0,0,0,0,0,0,0,0,0,bc,ibc)        6d11s24
             end if                                                        3d13s20
             jtmp=itmp                                                     10d28s20
             jeri=ieri                                                     10d28s20
             do in=1,nerimat                                               10d28s20
              do id=0,nszd0-1                                              10d28s20
               idp=id+iod                                                4d9s24
               do ic=0,nszc0-1                                             10d28s20
                icp=ic+ioc                                               4d9s24
                do ib=0,nszb0-1                                          4d9s24
                 ibp=ib+iob                                              4d9s24
                 do ia=0,nsza0-1                                         4d9s24
                  iap=ia+ioa                                             4d9s24
                  idcba=jeri+id+nszd0*(ic+nszc0*(ib+nszb0*ia))           4d9s24
                  iabcd=jtmp+iap+nszaa0*(ibp+nszba0*(icp+nszca0*idp))    4d9s24
                  bc(iabcd)=bc(iabcd)+bc(idcba)                         6d20s24
                 end do                                                  4d9s24
                end do                                                       5d12s10
               end do                                                        5d12s10
              end do                                                         5d12s10
              jtmp=jtmp+nwtmp                                              10d28s20
              jeri=jeri+nweri                                              10d28s20
             end do                                                        10d28s20
             ibcoff=ibc00                                               6d11s24
            end if                                                      6d11s24
            if(ibc(jb8).eq.iatom)then                                   6d11s24
             if(idorel.ne.0)then                                           3d13s20
              call restind(                                              6d11s24
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       6d12s24
     $           ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),       6d12s24
     $           ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3),       4d9s24
     $           ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3),       4d9s24
     $       ieri,jdoit,sigb,idorel,ascale,nmat,idx,idy,idz,0,0,0,      6d12s24
     $            0,0,0, 0,0,0, bc,ibc)                                 6d11s24
             else                                                           3d13s20
              call erird                                                6d11s24
     $                (ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3), 4d9s24
     $                 ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),    6d11s24
     $                 ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3), 4d9s24
     $                 ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3), 4d9s24
     $         ieri,0,sigb,idx,idy,idz,0,0,0,0,0,0,0,0,0,bc,ibc)        6d12s24
             end if                                                        3d13s20
             jtmp=itmp                                                     10d28s20
             jeri=ieri                                                     10d28s20
             do in=1,nerimat                                               10d28s20
              do id=0,nszd0-1                                              10d28s20
               idp=id+iod                                                4d9s24
               do ic=0,nszc0-1                                             10d28s20
                icp=ic+ioc                                               4d9s24
                do ib=0,nszb0-1                                          4d9s24
                 ibp=ib+iob                                              4d9s24
                 do ia=0,nsza0-1                                         4d9s24
                  iap=ia+ioa                                             4d9s24
                  idcba=jeri+id+nszd0*(ic+nszc0*(ia+nsza0*ib))          6d12s24
                  iabcd=jtmp+iap+nszaa0*(ibp+nszba0*(icp+nszca0*idp))    4d9s24
                  bc(iabcd)=bc(iabcd)+bc(idcba)                         6d20s24
                 end do                                                  4d9s24
                end do                                                       5d12s10
               end do                                                        5d12s10
              end do                                                         5d12s10
              jtmp=jtmp+nwtmp                                              10d28s20
              jeri=jeri+nweri                                              10d28s20
             end do                                                        10d28s20
             ibcoff=ibc00                                               6d11s24
            end if                                                      6d11s24
            if(ibc(jc8).eq.iatom)then                                   6d11s24
             if(idorel.ne.0)then                                           3d13s20
              call restind(                                              6d11s24
     $           ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),       4d9s24
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       3d13s20
     $           ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3),       4d9s24
     $           ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3),       4d9s24
     $       ieri,jdoit,sigc,idorel,ascale,nmat,0,0,0,                  6d11s24
     $            0,0,0,idx,idy,idz,0,0,0,bc,ibc)                       6d11s24
             else                                                           3d13s20
              call erird                                                6d11s24
     $             (ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),    6d11s24
     $                 ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3), 4d9s24
     $                 ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3), 4d9s24
     $                 ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3), 4d9s24
     $         ieri,0,sigc,0,0,0,0,0,0,idx,idy,idz,0,0,0,bc,ibc)        6d11s24
             end if                                                        3d13s20
             jtmp=itmp                                                     10d28s20
             jeri=ieri                                                     10d28s20
             do in=1,nerimat                                               10d28s20
              do id=0,nszd0-1                                              10d28s20
               idp=id+iod                                                4d9s24
               do ic=0,nszc0-1                                             10d28s20
                icp=ic+ioc                                               4d9s24
                do ib=0,nszb0-1                                          4d9s24
                 ibp=ib+iob                                              4d9s24
                 do ia=0,nsza0-1                                         4d9s24
                  iap=ia+ioa                                             4d9s24
                  idcba=jeri+id+nszd0*(ic+nszc0*(ib+nszb0*ia))           4d9s24
                  iabcd=jtmp+iap+nszaa0*(ibp+nszba0*(icp+nszca0*idp))    4d9s24
                  bc(iabcd)=bc(iabcd)+bc(idcba)                         6d20s24
                 end do                                                  4d9s24
                end do                                                       5d12s10
               end do                                                        5d12s10
              end do                                                         5d12s10
              jtmp=jtmp+nwtmp                                              10d28s20
              jeri=jeri+nweri                                              10d28s20
             end do                                                        10d28s20
             ibcoff=ibc00                                               6d11s24
            end if                                                      6d11s24
            if(ibc(jd8).eq.iatom)then                                   6d11s24
             if(idorel.ne.0)then                                           3d13s20
              call restind(                                              6d11s24
     $           ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),       4d9s24
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       3d13s20
     $           ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3),       6d12s24
     $           ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3),       6d12s24
     $       ieri,jdoit,sigd,idorel,ascale,nmat,0,0,0,                  6d11s24
     $            0,0,0,idx,idy,idz, 0,0,0,bc,ibc)                      6d12s24
             else                                                           3d13s20
              call erird                                                6d11s24
     $             (ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),    6d11s24
     $                 ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3), 4d9s24
     $                 ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3), 6d12s24
     $                 ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3), 6d12s24
     $         ieri,0,sigd,0,0,0,0,0,0,idx,idy,idz,0,0,0,bc,ibc)        6d12s24
             end if                                                        3d13s20
             jtmp=itmp                                                     10d28s20
             jeri=ieri                                                     10d28s20
             do in=1,nerimat                                               10d28s20
              do id=0,nszd0-1                                              10d28s20
               idp=id+iod                                                4d9s24
               do ic=0,nszc0-1                                             10d28s20
                icp=ic+ioc                                               4d9s24
                do ib=0,nszb0-1                                          4d9s24
                 ibp=ib+iob                                              4d9s24
                 do ia=0,nsza0-1                                         4d9s24
                  iap=ia+ioa                                             4d9s24
                  idcba=jeri+ic+nszc0*(id+nszd0*(ib+nszb0*ia))          6d12s24
                  iabcd=jtmp+iap+nszaa0*(ibp+nszba0*(icp+nszca0*idp))    4d9s24
                  bc(iabcd)=bc(iabcd)+bc(idcba)                         6d20s24
                 end do                                                  4d9s24
                end do                                                       5d12s10
               end do                                                        5d12s10
              end do                                                         5d12s10
              jtmp=jtmp+nwtmp                                              10d28s20
              jeri=jeri+nweri                                              10d28s20
             end do                                                        10d28s20
             ibcoff=ibc00                                               6d11s24
            end if                                                      6d11s24
           end do                                                       4d9s24
          end do                                                        4d9s24
         end do
        end do
        nabc=nszaa0*nszba0*nszca0                                       4d9s24
        if(npassd.eq.2)then                                             4d9s24
         jtmp=itmp                                                      10d28s20
         do in=1,nerimat                                                10d28s20
          if(idoit.eq.1)then
           write(6,*)('before sum/dif: '),in
           call prntm2(bc(jtmp),nabc,nszda0,nabc)
          end if
          do id=0,nszd0-1                                               10d28s20
           it1=jtmp+nabc*id                                             4d9s24
           it2=jtmp+nabc*(id+nszd0)                                     4d9s24
           do iabc=0,nabc-1                                             4d9s24
            sum=srh*(bc(it1+iabc)+bc(it2+iabc))                         5d12s10
            diff=srh*(-bc(it1+iabc)+bc(it2+iabc))                       5d12s10
            bc(it1+iabc)=sum                                            5d12s10
            bc(it2+iabc)=diff                                           5d12s10
           end do                                                       5d12s10
          end do                                                        5d12s10
          if(idoit.eq.1)then                                            4d9s24
           write(6,*)('after sum/dif: '),in                             4d9s24
           call prntm2(bc(jtmp),nabc,nszda0,nabc)                       4d9s24
          end if                                                        4d9s24
          jtmp=jtmp+nwtmp                                               10d28s20
         end do                                                         10d28s20
        end if                                                          5d12s10
        nab=nszaa0*nszba0                                               4d9s24
        if(npassc.eq.2)then                                             4d9s24
         jtmp=itmp                                                      10d28s20
         do in=1,nerimat                                                10d28s20
          if(idoit.eq.1)then
           write(6,*)('before sum/dif: '),in
           call prntm2(bc(jtmp),nab,nszca0*nszda0,nab)
          end if
          do id=0,nszda0-1                                               10d28s20
           do ic=0,nszc0-1                                              4d9s24
            it1=jtmp+nab*(ic+nszca0*id)                                 4d9s24
            it2=jtmp+nab*(ic+nszc0+nszca0*id)                           4d9s24
            do iab=0,nab-1                                              4d9s24
             sum=srh*(bc(it1+iab)+bc(it2+iab))                          4d9s24
             diff=srh*(-bc(it1+iab)+bc(it2+iab))                        4d9s24
             bc(it1+iab)=sum                                            5d12s10
             bc(it2+iab)=diff                                           5d12s10
            end do                                                       5d12s10
           end do                                                        5d12s10
          end do                                                        4d9s24
          if(idoit.eq.1)then                                            4d9s24
           write(6,*)('after sum/dif: '),in                             4d9s24
           call prntm2(bc(jtmp),nab,nszda0*nszca0,nab)                  4d9s24
          end if                                                        4d9s24
          jtmp=jtmp+nwtmp                                               10d28s20
         end do                                                         10d28s20
        end if                                                          5d12s10
        ncd=nszca0*nszda0                                               4d9s24
        if(npassb.eq.2)then                                             4d9s24
         jtmp=itmp                                                      10d28s20
         do in=1,nerimat                                                10d28s20
          if(idoit.eq.1)then
           write(6,*)('before sum/dif: '),in
           call prntm2(bc(jtmp),nab,ncd,nab)
          end if
          do icd=0,ncd-1                                                4d9s24
           do ib=0,nszb0-1                                              4d9s24
            it1=jtmp+nszaa0*(ib+nszba0*icd)                             4d9s24
            it2=jtmp+nszaa0*(ib+nszb0+nszba0*icd)                       4d9s24
            do ia=0,nszaa0-1                                            4d9s24
             sum=srh*(bc(it1+ia)+bc(it2+ia))                            4d9s24
             diff=srh*(-bc(it1+ia)+bc(it2+ia))                          4d9s24
             bc(it1+ia)=sum                                             4d9s24
             bc(it2+ia)=diff                                            4d9s24
            end do                                                       5d12s10
           end do                                                        5d12s10
          end do                                                        4d9s24
          if(idoit.eq.1)then                                            4d9s24
           write(6,*)('after sum/dif: '),in                             4d9s24
           call prntm2(bc(jtmp),nab,ncd,nab)                            4d9s24
          end if                                                        4d9s24
          jtmp=jtmp+nwtmp                                               10d28s20
         end do                                                         10d28s20
        end if                                                          5d12s10
        nbcd=nszba0*nszca0*nszda0                                       4d9s24
        if(npassa.eq.2)then                                             4d9s24
         jtmp=itmp                                                      10d28s20
         do in=1,nerimat                                                10d28s20
          if(idoit.eq.1)then
           write(6,*)('before sum/dif: '),in
           call prntm2(bc(jtmp),nab,ncd,nab)
          end if
          do ibcd=0,nbcd-1                                                4d9s24
           it1=jtmp+nszaa0*ibcd                                         4d9s24
           it2=jtmp+nsza0+nszaa0*ibcd                                   4d9s24
           do ia=0,nsza0-1                                              4d9s24
            sum=srh*(bc(it1+ia)+bc(it2+ia))                             4d9s24
            diff=srh*(-bc(it1+ia)+bc(it2+ia))                           4d9s24
            bc(it1+ia)=sum                                              4d9s24
            bc(it2+ia)=diff                                             4d9s24
           end do                                                        5d12s10
          end do                                                        4d9s24
          if(idoit.eq.1)then                                            4d9s24
           write(6,*)('after sum/dif: '),in                             4d9s24
           call prntm2(bc(jtmp),nab,ncd,nab)                            4d9s24
          end if                                                        4d9s24
          jtmp=jtmp+nwtmp                                               10d28s20
         end do                                                         10d28s20
        end if                                                          5d12s10
        sz=0d0                                                          4d9s24
        do in=0,nerimat-1                                               4d9s24
         do id=1,nszda0                                                  4d9s24
          idm=id-1                                                      4d9s24
          idd=id+ibc(jd4)                                                4d9s24
          isd=isstor(idd)                                                4d9s24
          iddd=ibstor(idd)-1
          do ic=1,nszca0                                                4d9s24
           icm=ic-1                                                     4d9s24
           icc=ic+ibc(jc4)                                              4d9s24
           isc=isstor(icc)                                              4d9s24
           iccc=ibstor(icc)-1                                           4d9s24
           iscd=multh(isc,isd)                                          4d9s24
           do ib=1,nszba0                                               4d9s24
            ibm=ib-1                                                    4d9s24
            ibb=ib+ibc(jb4)                                             4d9s24
            isb=isstor(ibb)                                             4d9s24
            ibbb=ibc(jptb+ib)                                           4d9s24
            isbcd=multh(isb,iscd)                                       4d9s24
            iad1=itmp-1+nszaa0*(ibm+nszba0*(icm+nszca0*(idm+nszda0*in)))4d9s24
            iad2=ieraw2(isbcd,isb,isc)+nan(isbcd)*(ibbb                 4d9s24
     $           +nbn(isb)*(iccc+nbasisp(isc)*(iddd+nbasisp(isd)*in)))  4d9s24
            iad2b=ieraw2(isbcd,isb,isd)+nan(isbcd)*(ibbb                4d10s24
     $           +nbn(isb)*(iddd+nbasisp(isd)*(iccc+nbasisp(isc)*in)))  4d10s24
            do ia=1,nszaa0                                              4d9s24
             iaa=ia+ibc(ja4)                                            4d9s24
             isa=isstor(iaa)                                            4d9s24
             if(isa.eq.isbcd)then                                       4d9s24
              iaaa=ibc(jpta+ia)                                         4d9s24
              bc(iad2+iaaa)=bc(iad1+ia)                                 4d9s24
              bc(iad2b+iaaa)=bc(iad1+ia)                                4d10s24
             else                                                       4d9s24
              sz=sz+bc(iad1+ia)**2                                      4d9s24
             end if                                                      4d9s24
            end do                                                      4d9s24
           end do                                                       4d9s24
          end do                                                        4d9s24
         end do                                                         4d9s24
        end do                                                          4d9s24
        sz=sqrt(sz/dfloat(nab*ncd*nerimat))                             4d9s24
        if(sz.gt.1d-9.or.sz.ne.sz)then                                  6d18s24
         write(6,*)('we have symmetry breaking in paraeri!!!! '),sz     4d9s24
         stop 'paraeridd'                                                 4d9s24
        end if                                                          4d9s24
        ibcoff=itmp                                                     4d9s24
       end do
       do isd=1,nsymb                                                   4d9s24
        nbpd=nbasisp(isd)*ncomp                                         4d9s24
        ntest=nduse(isd)                                                6d11s24
        do isc=1,nsymb                                                  4d9s24
         ntestc=nduse2(isc)                                             6d11s24
         iscd=multh(isc,isd)                                            4d9s24
         nbpc=nbasisp(isc)*ncomp                                         4d9s24
         do isb=1,nsymb                                                 4d9s24
          isa=multh(isb,iscd)                                           4d9s24
          if(min(nan(isa),nbn(isb),nbasisp(isc),nbasisp(isd)).gt.0)then 6d20s24
           nab=nan(isa)*nbn(isb)                                        4d9s24
           nabc=nab*nbasisp(isc)                                        4d9s24
           nabd=nab*nduse(isd)                                          4d9s24
           nabcd=nabc*nbasisp(isd)                                      4d9s24
           itmp1=ibcoff                                                 4d9s24
           itmp2=itmp1+nabc*nduse(isd)                                  4d9s24
           ibcoff=itmp2+nabd*nduse2(isc)                                4d9s24
           call enough('paraeridd.tmp12',bc,ibc)                          4d9s24
           jeraw2=ieraw2(isa,isb,isc)                                   4d9s24
           do in=0,nerimat-1                                            4d9s24
            jorbd=iorb(isd)                                             4d9s24
            jorbc=iorb(isc)                                             4d9s24
            if(in.eq.1.or.in.eq.3)then                                  4d9s24
             jorbd=jorbd+nbasisp(isd)                                   4d9s24
             jorbc=jorbc+nbasisp(isc)                                   4d9s24
            end if                                                      4d9s24
            if(iaddr.gt.0)then                                          6d13s24
             call der4v(iaddr,naddr,d4v,isa,isb,isc,isd,nsymb,multh,    6d13s24
     $           nbasisp,nerimat,nan,iptap,nbn,iptbp,in,ncomp,lsame,    6d13s24
     $           bc(jeraw2),bc,ibc)                                     6d13s24
            end if                                                      6d13s24
            call dgemm('n','n',nabc,nbasdws(isd),nbasisp(isd),1d0,      6d11s24
     $           bc(jeraw2),nabc,bc(jorbd),nbpd,0d0,bc(itmp1),nabc,     4d9s24
     $           'paraeridd.eraw2')                                       4d9s24
            do id=0,nduse(isd)-1                                        4d9s24
             do ic=0,nbasisp(isc)-1                                     4d9s24
              icd=itmp1+nab*(ic+nbasisp(isc)*id)                        4d9s24
              idc=jeraw2+nab*(id+nduse(isd)*ic)                         4d9s24
              do iab=0,nab-1                                            4d9s24
               bc(idc+iab)=bc(icd+iab)                                  4d9s24
              end do                                                    4d9s24
             end do                                                     4d9s24
            end do                                                      4d9s24
            jtmp=itmpab2(isa,isb,isc)                                   4d10s24
            if(in.gt.1)then                                             4d9s24
             jtmp=jtmp+nabd*nduse2(isc)                                 4d9s24
            end if                                                      4d9s24
            if(min(nabd,nduse2(isc)).gt.0)then                          6d20s24
             call dgemm('n','n',nabd,nduse2(isc),nbasisp(isc),1d0,       4d9s24
     $           bc(jeraw2),nabd,bc(jorbc),nbpc,0d0,bc(itmp2),nabd,      4d9s24
     $           'paraeridd.eraw2c')                                      4d9s24
            end if                                                      6d20s24
            ndc=nduse(isd)*nduse2(isc)                                  4d9s24
            do iab=0,nab-1                                              4d9s24
             do idc=0,ndc-1                                             4d9s24
              iabdc=itmp2+iab+nab*idc                                   4d9s24
              idcab=jtmp+idc+ndc*iab                                    4d9s24
              bc(idcab)=bc(idcab)+bc(iabdc)                             4d9s24
             end do                                                     4d9s24
            end do                                                      4d9s24
            jeraw2=jeraw2+nabcd                                         4d9s24
           end do                                                       4d9s24
           ibcoff=itmp1                                                 4d9s24
          end if                                                        4d9s24
         end do                                                         4d9s24
        end do                                                          4d9s24
       end do                                                           4d9s24
c
c     final storage                                                     1d10s11
c
        do isd=1,nsymb                                                  5d12s10
         do isc=1,nsymb                                                 5d12s10
          ieraw(isc,isd)=nduse2(isc)*nduse(isd)                         8d18s23
         end do                                                         5d12s10
        end do                                                          5d12s10
        szq=0d0                                                         10d28s20
        iszq=0                                                          10d28s20
        do ia=1,nszaa0                                                  8d26s15
         iaa=ia+ibc(ja4)                                                2d3s12
         isa=isstor(iaa)                                                2d3s12
         iaaa=ibc(jpta+ia)                                              4d9s24
         do ib=1,nszba0                                                 8d26s15
          ibb=ib+ibc(jb4)                                               2d3s12
          isb=isstor(ibb)                                               2d3s12
          ibbb=ibc(jptb+ib)                                             4d9s24
          if(ibb.le.iaa)then                                            1d31s12
           isab=multh(isa,isb)                                           5d12s10
           do isd=1,nsymb                                                5d12s10
            isc=multh(isd,isab)                                          5d12s10
            ii=iptoh(isd,isc,isb)                                       1d10s11
            ncd=ieraw(isc,isd)                                          5d12s10
            if(ncd.gt.0)then                                            1d31s12
             if(ii.gt.0)then                                            1d31s12
              jj=itmpab(isc,isd)+ncd*(ib-1+nszba0*(ia-1))               10d30s20
              rms=0d0                                                   4d9s24
              nab=nan(isa)*nbn(isb)                                     4d10s24
              do in=0,ntmpmat-1
               kk=itmpab2(isa,isb,isc)+ncd*(iaaa                        4d9s24
     $              +nan(isa)*(ibbb+nbn(isb)*in))                         4d9s24
               do ix=0,ncd-1                                                1d10s11
                bc(ihcol(ii))=bc(kk+ix)                                 6d10s24
                ihcol(ii)=ihcol(ii)+1                                     5d12s10
               end do                                                     5d12s10
               jj=jj+ncd*nszba0*nszaa0                                  10d30s20
              end do                                                    10d30s20
             else if(iptoh(isd,isc,isa).gt.0)then                       1d31s12
              jj=itmpab(isc,isd)+ncd*(ib-1+nszba0*(ia-1))               10d31s20
              rms=0d0                                                   4d9s24
              ii=iptoh(isd,isc,isa)                                     1d31s12
              do in=0,ntmpmat-1                                         10d30s20
               kk=itmpab2(isa,isb,isc)+ncd*(iaaa                        4d9s24
     $              +nan(isa)*(ibbb+nbn(isb)*in))                         4d9s24
               do ix=0,ncd-1                                            10d30s20
                bc(ihcol(ii))=bc(kk+ix)                                 6d10s24
                ihcol(ii)=ihcol(ii)+1                                     5d12s10
               end do                                                     5d12s10
               jj=jj+ncd*nszba0*nszaa0                                  10d30s20
              end do                                                    10d30s20
             end if                                                     1d31s12
            end if                                                      1d31s12
           end do                                                        5d12s10
          end if                                                        5d12s10
         end do                                                         5d12s10
        end do                                                          5d12s10
       ibcoff=ipta                                                      4d9s24
      end do                                                            2d19s10
      call dws_sync                                                     2d22s10
      ihcol(1)=ihmat                                                    1d10s11
      do i=2,nsdlkh                                                     1d10s11
       ihcol(i)=ihcol(i-1)+nhcol(i-1)                                   1d10s11
      end do                                                            5d12s10
      if(iprtr(18).ne.0)then                                            5d4s20
      write(6,*)('we have this procs part of hmat ')                    3d18s10
      write(6,*)('what we''ve got: ')                                   1d10s11
      memtop=0d0
      do isb=1,nsymb                                                    1d10s11
       do isc=1,nsymb                                                   1d10s11
        isbc=multh(isc,isb)                                             1d10s11
        do isd=1,nsymb                                                  1d10s11
         isa=multh(isd,isbc)                                            1d10s11
         ii=iptoh(isd,isc,isb)                                          1d10s11
         ncd=nduse2(isc)*nduse(isd)                                     8d18s23
         if(ii.gt.0.and.ncd.gt.0)then
          if(nhcol(ii).gt.0)then                                        1d10s11
           write(6,*)('integrals for symmetry '),isa,isb,isc,isd,ii
           mcol=nhcol(ii)/ncd                                           1d10s11
           memtop=max(memtop,ihcol(ii)+ncd*mcol-1)                      3d19s12
           call prntm2(bc(ihcol(ii)),ncd,mcol,ncd)                      1d10s11
          end if                                                        1d10s11
         end if                                                         1d10s11
        end do                                                          1d10s11
       end do                                                           1d10s11
      end do                                                            1d10s11
      write(6,*)('memtop '),memtop
      end if                                                            3d16s12
      return
      end                                                               2d19s10
