c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine parap42(natom,ngaus,ibdat,ixmt,isym,iapair,            2d15s22
     $                  ibstor,isstor,iso,idwsdeb,ascale,               2d15s22
     $     ioppt,npt,opdata,iopdata,iosym,nop,multh,nbb,nbasisp,nbasdwx,2d15s22
     $     iorb,opname,isopt,nsopt,id,nopso,iopso,idoubo,bc,ibc)        11d9s22
      implicit real*8 (a-h,o-z)
      external second
      integer*8 nn8                                                     2d19s10
      include "common.hf"
      include "common.store"
      include "common.spher"
c
      logical log1,lprint                                               3d2s22
      character*10 opname(*)                                            2d17s22
      integer*8 ibstor,isstor                                           5d6s10
      include "common.basis"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension ixmt(8,id),isym(3,*),iapair(3,*),ibstor(*),isstor(*),   2d15s22
     $     iso(3),cartb(3),cartk(3),ioppt(*),npt(*),opdata(*),          12d16s19
     $     iopdata(7,*),iosym(*),multh(8,8),nbasisp(8),iorb(8),         5d27s21
     $     nbasdwx(8),ixyz(3,3),xcartb(3),xcartk(3),isopt(4,4),         2d15s22
     $     iopso(5,*),ider(3,2),idoubo(*)                               2d18s22
      data ixyz/1,0,0,0,1,0,0,0,1/                                      1d11s20
      lprint=mynowprog.eq.0                                             3d2s22
      ncomp=2                                                           1d11s20
      if(lprint)then                                                    3d2s22
       write(6,*)('Hi, my name is parap42')
       write(6,*)
     $('I''m going to generate matrices for'),
     $     (' spin-orbit transition moments')
      end if                                                            3d2s22
      nopso=0                                                           2d15s22
      do j=1,nop                                                        2d25s22
       if(opname(j)(1:1).eq.'m'.or.opname(j)(1:1).eq.'l'.or.            2d25s22
     $        opname(j)(1:1).eq.'Q'.or.opname(j)(1:1).eq.'R'.or.        2d15s22
     $        opname(j)(1:1).eq.'I')then                                2d15s22
        do i=1,nsopt
         if(isopt(2,i).ne.0.or.isopt(3,i).ne.0)then                     2d28s22
          if(isopt(2,i).eq.0)then                                         2d18s22
           id1=1                                                          2d15s22
           id2=3                                                          2d15s22
          else                                                            2d15s22
           if(isopt(3,i).eq.0)then                                        2d18s22
            id1=2                                                       2d28s22
            id2=1                                                       2d28s22
           else                                                           2d15s22
            id1=3                                                       2d28s22
            id2=2                                                       2d28s22
           end if                                                         2d15s22
          end if                                                          2d15s22
          isnew=multh(isopt(1,i),iosym(j))                              2d15s22
          nopso=nopso+1                                                 2d15s22
          iopso(1,nopso)=isnew                                          2d15s22
          iopso(2,nopso)=j                                              2d15s22
          iopso(3,nopso)=i                                              2d15s22
          iopso(4,nopso)=id1                                            2d15s22
          iopso(5,nopso)=id2                                            2d15s22
          do isa=1,nsymb                                                2d15s22
           isb=multh(isa,isnew)                                         2d15s22
           ixmt(isa,nopso)=ibcoff                                       2d15s22
           ibcoff=ixmt(isa,nopso)+nbasisp(isa)*nbasisp(isb)             2d16s22
          end do                                                        2d15s22
         end if                                                         2d15s22
        end do                                                          2d15s22
       end if                                                           2d15s22
      end do
      call enough('parap42.  1',bc,ibc)
      nbb=ibcoff-ixmt(1,1)                                              2d16s22
      do iz=ixmt(1,1),ibcoff-1                                          2d15s22
       bc(iz)=0d0                                                       2d15s22
      end do                                                            2d15s22
      ibcoffo=ibcoff                                                    2d19s10
      clight=0.5d0/sqrt(ascale)                                         2d21s20
      ascale2=ascale*2d0                                                8d20s15
c
c     build shell pair order
c
      srh=sqrt(0.5d0)                                                   5d3s10
      ascale2=ascale*2d0                                                8d20s15
      nn=ngaus*ngaus                                                    12d16s19
      nn2=nn*2                                                          2d19s10
c
c     properties:
c     there will be nopso of them, with rules given by ipt,npt,opdata,
c     iopdata, and iosym.
c
      ndo=nopso
      ipair=ibcoff                                                      2d19s10
      ibcoff=ipair+nn*ndo*3                                             2d19s10
      i12=ibcoff                                                        2d19s10
      ibcoff=i12+nn*4                                                   2d19s10
      call enough('parap42.  2',bc,ibc)
      j12=i12                                                           2d19s10
      jbdat=ibdat-1                                                     2d19s10
      ngaus7=ngaus*7                                                    5d3s10
      do i1=1,ngaus                                                     2d19s10
       do i2=1,ngaus                                                    12d16s19
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
      call idsortdws(ibc(i12),ibc(is12),nn8)                            1d18s23
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
        ibc(jpair)=k                                                    12d16s19
        jpair=jpair+1                                                   2d19s10
       end do                                                           2d19s10
      end do                                                            2d19s10
      nneed=nn*ndo                                                      2d19s10
      loopit=0                                                          5d25s18
      do i=1+mynowprog,nneed,mynprocg                                   2d19s10
       loopit=loopit+1                                                  5d25s18
       jpair=ipair+3*(i-1)                                              2d19s10
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
       ioprso=ibc(jpair+2)                                              2d15s22
       iopr=iopso(2,ioprso)                                             2d15s22
       nbra=2*ibc(jbra)+1                                               5d3s10
       nbraa=nbra                                                       5d3s10
       if(iapair(1,ibc(jbra8)).gt.0)nbraa=nbra*2                        5d3s10
       nket=2*ibc(jket)+1                                               5d3s10
       nketa=nket                                                       5d3s10
       if(iapair(1,ibc(jket8)).gt.0)nketa=nket*2                        5d3s10
       itmp1=ibcoff                                                     5d3s10
       itmp2=itmp1+nbraa*nketa                                          5d3s10
       ibcoff=itmp2+nbraa*nketa                                         5d3s10
       ltmp1=ibcoff                                                     1d11s20
       ltmp2=ltmp1+nbraa*nketa                                          1d11s20
       ibcoff=ltmp2+nbraa*nketa                                         1d11s20
       call enough('parap42.  3',bc,ibc)
       jopr=ioppt(iopr)                                                 12d16s19
       itmp3=ibcoff                                                     12d15s19
       ibcoff=itmp3+nbra*nket                                           12d15s19
       call enough('parap42.  4',bc,ibc)
       do ix=itmp3,ibcoff-1                                             12d15s19
        bc(ix)=0d0                                                      12d15s19
       end do                                                           12d15s19
       xcartb(1)=bc(jbra5)                                              5d27s21
       xcartb(2)=bc(jbra6)                                              5d27s21
       xcartb(3)=bc(jbra7)                                              5d27s21
       xcartk(1)=bc(jket5)                                              5d27s21
       xcartk(2)=bc(jket6)                                              5d27s21
       xcartk(3)=bc(jket7)                                              5d27s21
       do ipart=1,npt(iopr)                                             12d15s19
        ider(1,1)=0                                                     2d15s22
        ider(2,1)=0                                                     2d15s22
        ider(3,1)=0                                                     2d15s22
        ider(iopso(4,ioprso),1)=1
        ider(1,2)=iopdata(4,jopr)
        ider(2,2)=iopdata(5,jopr)
        ider(3,2)=iopdata(6,jopr)
        ider(iopso(5,ioprso),2)=ider(iopso(5,ioprso),2)+1               2d15s22
        call onep(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),    2d22s22
     $      bc(jbra7),ibc(jket),bc(jket2),bc(jket3),                    2d22s22
     $    bc(jket5),bc(jket6),bc(jket7),nsymb,ixmat,ider(1,1),ider(2,1),2d22s22
     $        ider(3,1),iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),2d15s22
     $      ider(1,2),ider(2,2),ider(3,2),bc,ibc)                       11d9s22
        do ix=0,nbra*nket-1                                             12d15s19
         bc(itmp3+ix)=bc(itmp3+ix)+opdata(jopr)*bc(ixmat+ix)             12d15s19
        end do                                                          12d15s19
        ider(1,1)=0                                                     2d15s22
        ider(2,1)=0                                                     2d15s22
        ider(3,1)=0                                                     2d15s22
        ider(iopso(5,ioprso),1)=1
        ider(1,2)=iopdata(4,jopr)
        ider(2,2)=iopdata(5,jopr)
        ider(3,2)=iopdata(6,jopr)
        ider(iopso(4,ioprso),2)=ider(iopso(4,ioprso),2)+1               2d18s22
        call onep(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),    2d22s22
     $      bc(jbra7),ibc(jket),bc(jket2),bc(jket3),                    2d22s22
     $   bc(jket5),bc(jket6),bc(jket7),nsymb,ixmat,ider(1,1),ider(2,1), 2d22s22
     $        ider(3,1),iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),2d15s22
     $      ider(1,2),ider(2,2),ider(3,2),bc,ibc)                       11d9s22
        do ix=0,nbra*nket-1                                             12d15s19
         bc(itmp3+ix)=bc(itmp3+ix)-opdata(jopr)*bc(ixmat+ix)             12d15s19
        end do                                                          12d15s19
        jopr=jopr+1                                                     12d15s19
       end do                                                           2d15s22
       ibcoff=itmp3                                                     12d15s19
       ib1=itmp3                                                        12d15s19
       ibu=ib1                                                          5d7s10
       itmpu=itmp1                                                      5d7s10
       do ik=1,nket                                                     12d15s19
        do ib=1,nbra                                                    12d15s19
         iad1=ik-1+nket*(ib-1)                                          12d15s19
         iad2=ib-1+nbraa*(ik-1)                                         12d15s19
         bc(itmpu+iad2)=bc(ibu+iad1)                                    12d15s19
        end do                                                          12d15s19
       end do                                                           12d15s19
       if(nketa.ne.nket)then                                            5d3s10
        cartk(1)=bc(jket5)*dfloat(isym(1,iapair(2,ibc(jket8))))         5d5s10
        cartk(2)=bc(jket6)*dfloat(isym(2,iapair(2,ibc(jket8))))         5d5s10
        cartk(3)=bc(jket7)*dfloat(isym(3,iapair(2,ibc(jket8))))         5d5s10
        jopr=ioppt(iopr)                                                12d16s19
        itmp3=ibcoff                                                    12d15s19
        ibcoff=itmp3+nbra*nket                                          12d15s19
        call enough('parap42.  5',bc,ibc)
        do ix=itmp3,ibcoff-1                                            12d15s19
         bc(ix)=0d0                                                     12d15s19
        end do                                                          12d15s19
        do ipart=1,npt(iopr)                                            12d15s19
         ider(1,1)=0                                                    2d15s22
         ider(2,1)=0                                                    2d15s22
         ider(3,1)=0                                                    2d15s22
         ider(iopso(4,ioprso),1)=1
         ider(1,2)=iopdata(4,jopr)
         ider(2,2)=iopdata(5,jopr)
         ider(3,2)=iopdata(6,jopr)
         ider(iopso(5,ioprso),2)=ider(iopso(5,ioprso),2)+1              2d15s22
         call onep(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),   2d15s22
     $      bc(jbra7),ibc(jket),bc(jket2),bc(jket3),                    12d15s19
     $      cartk,cartk(2),cartk(3),nsymb,ixmat,ider(1,1),ider(2,1),    2d15s22
     $        ider(3,1),iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),2d15s22
     $      ider(1,2),ider(2,2),ider(3,2),bc,ibc)                       11d9s22
         do ix=0,nbra*nket-1                                            12d15s19
          bc(itmp3+ix)=bc(itmp3+ix)+opdata(jopr)*bc(ixmat+ix)             12d15s19
         end do                                                          12d15s19
         ider(1,1)=0                                                    2d15s22
         ider(2,1)=0                                                    2d15s22
         ider(3,1)=0                                                    2d15s22
         ider(iopso(5,ioprso),1)=1
         ider(1,2)=iopdata(4,jopr)
         ider(2,2)=iopdata(5,jopr)
         ider(3,2)=iopdata(6,jopr)
         ider(iopso(4,ioprso),2)=ider(iopso(4,ioprso),2)+1              2d18s22
         call onep(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),   2d15s22
     $      bc(jbra7),ibc(jket),bc(jket2),bc(jket3),                    2d15s22
     $      cartk,cartk(2),cartk(3),nsymb,ixmat,ider(1,1),ider(2,1),    2d15s22
     $        ider(3,1),iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),2d15s22
     $      ider(1,2),ider(2,2),ider(3,2),bc,ibc)                       11d9s22
         do ix=0,nbra*nket-1                                            12d15s19
          bc(itmp3+ix)=bc(itmp3+ix)-opdata(jopr)*bc(ixmat+ix)             12d15s19
         end do                                                         12d15s19
         jopr=jopr+1                                                    12d15s19
        end do                                                          12d15s19
        ibcoff=itmp3                                                    12d15s19
        ib1=itmp3                                                       12d15s19
        ibu=ib1                                                         5d7s10
        itmpu=itmp1                                                     5d7s10
        do ik=1,nket                                                    5d3s10
         do ib=1,nbra                                                   5d3s10
          iad1=ik-1+nket*(ib-1)                                         5d3s10
          iad2=ib-1+nbraa*(ik-1+nket)                                   5d3s10
          bc(itmpu+iad2)=bc(ibu+iad1)                                   8d20s15
         end do                                                         5d3s10
        end do                                                          5d3s10
       end if                                                           2d15s22
       if(nbraa.ne.nbra)then                                            5d3s10
        cartb(1)=bc(jbra5)*dfloat(isym(1,iapair(2,ibc(jbra8))))         5d5s10
        cartb(2)=bc(jbra6)*dfloat(isym(2,iapair(2,ibc(jbra8))))         5d5s10
        cartb(3)=bc(jbra7)*dfloat(isym(3,iapair(2,ibc(jbra8))))         5d5s10
        jopr=ioppt(iopr)                                                12d16s19
        itmp3=ibcoff                                                    12d15s19
        ibcoff=itmp3+nbra*nket                                          12d15s19
        call enough('parap42.  6',bc,ibc)
        do ix=itmp3,ibcoff-1                                            12d15s19
         bc(ix)=0d0                                                     12d15s19
        end do                                                          12d15s19
        do ipart=1,npt(iopr)                                            12d15s19
         ider(1,1)=0                                                    2d15s22
         ider(2,1)=0                                                    2d15s22
         ider(3,1)=0                                                    2d15s22
         ider(iopso(4,ioprso),1)=1
         ider(1,2)=iopdata(4,jopr)
         ider(2,2)=iopdata(5,jopr)
         ider(3,2)=iopdata(6,jopr)
         ider(iopso(5,ioprso),2)=ider(iopso(5,ioprso),2)+1              2d15s22
         call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),        2d15s22
     $      cartb(3),ibc(jket),bc(jket2),bc(jket3),                     2d15s22
     $    bc(jket5),bc(jket6),bc(jket7),nsymb,ixmat,ider(1,1),ider(2,1),2d15s22
     $        ider(3,1),iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),2d15s22
     $      ider(1,2),ider(2,2),ider(3,2),bc,ibc)                       11d9s22
         do ix=0,nbra*nket-1                                            12d15s19
          bc(itmp3+ix)=bc(itmp3+ix)+opdata(jopr)*bc(ixmat+ix)             12d15s19
         end do                                                          12d15s19
         ider(1,1)=0                                                    2d15s22
         ider(2,1)=0                                                    2d15s22
         ider(3,1)=0                                                    2d15s22
         ider(iopso(5,ioprso),1)=1
         ider(1,2)=iopdata(4,jopr)
         ider(2,2)=iopdata(5,jopr)
         ider(3,2)=iopdata(6,jopr)
         ider(iopso(4,ioprso),2)=ider(iopso(4,ioprso),2)+1              2d18s22
         call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),        2d15s22
     $      cartb(3),ibc(jket),bc(jket2),bc(jket3),                     2d15s22
     $    bc(jket5),bc(jket6),bc(jket7),nsymb,ixmat,ider(1,1),ider(2,1),    2d15s22
     $        ider(3,1),iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),2d15s22
     $      ider(1,2),ider(2,2),ider(3,2),bc,ibc)                       11d9s22
         do ix=0,nbra*nket-1                                            12d15s19
          bc(itmp3+ix)=bc(itmp3+ix)-opdata(jopr)*bc(ixmat+ix)             12d15s19
         end do                                                         12d15s19
         jopr=jopr+1                                                    12d15s19
        end do                                                          12d15s19
        ibcoff=itmp3                                                    12d15s19
        ib1=itmp3                                                       12d15s19
        ibu=ib1                                                         5d7s10
        itmpu=itmp1                                                     5d7s10
        do ik=1,nket                                                    5d3s10
         do ib=1,nbra                                                   5d3s10
          iad1=ik-1+nket*(ib-1)                                         5d3s10
          iad2=ib-1+nbra+nbraa*(ik-1)                                   5d3s10
          bc(itmpu+iad2)=bc(ibu+iad1)                                   8d20s15
         end do                                                         5d7s10
        end do                                                          5d3s10
        if(nketa.ne.nket)then                                           5d3s10
         jopr=ioppt(iopr)                                               12d16s19
         itmp3=ibcoff                                                   12d15s19
         ibcoff=itmp3+nbra*nket                                         12d15s19
         call enough('parap42.  7',bc,ibc)
         do ix=itmp3,ibcoff-1                                            12d15s19
          bc(ix)=0d0                                                     12d15s19
         end do                                                         12d15s19
         do ipart=1,npt(iopr)                                           12d15s19
          ider(1,1)=0                                                   2d16s22
          ider(2,1)=0                                                   2d16s22
          ider(3,1)=0                                                   2d16s22
          ider(iopso(4,ioprso),1)=1                                     2d16s22
          ider(1,2)=iopdata(4,jopr)                                     2d16s22
          ider(2,2)=iopdata(5,jopr)                                     2d16s22
          ider(3,2)=iopdata(6,jopr)                                     2d16s22
          ider(iopso(5,ioprso),2)=ider(iopso(5,ioprso),2)+1             2d16s22
          call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),       2d16s22
     $          cartb(3),ibc(jket),bc(jket2),bc(jket3),                 2d16s22
     $          cartk,cartk(2),cartk(3),nsymb,ixmat,ider(1,1),ider(2,1),2d16s22
     $        ider(3,1),iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),2d16s22
     $          ider(1,2),ider(2,2),ider(3,2),bc,ibc)                   11d9s22
          do ix=0,nbra*nket-1                                           2d16s22
           bc(itmp3+ix)=bc(itmp3+ix)+opdata(jopr)*bc(ixmat+ix)          2d16s22
          end do                                                        2d16s22
          ider(1,1)=0                                                   2d16s22
          ider(2,1)=0                                                   2d16s22
          ider(3,1)=0                                                   2d16s22
          ider(iopso(5,ioprso),1)=1                                     2d16s22
          ider(1,2)=iopdata(4,jopr)                                     2d16s22
          ider(2,2)=iopdata(5,jopr)                                     2d16s22
          ider(3,2)=iopdata(6,jopr)                                     2d16s22
          ider(iopso(4,ioprso),2)=ider(iopso(4,ioprso),2)+1             2d18s22
          call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),       2d16s22
     $         cartb(3),ibc(jket),bc(jket2),bc(jket3),                  2d16s22
     $          cartk,cartk(2),cartk(3),nsymb,ixmat,ider(1,1),ider(2,1),2d16s22
     $        ider(3,1),iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),2d16s22
     $          ider(1,2),ider(2,2),ider(3,2),bc,ibc)                   11d9s22
          do ix=0,nbra*nket-1                                           2d16s22
           bc(itmp3+ix)=bc(itmp3+ix)-opdata(jopr)*bc(ixmat+ix)          2d16s22
          end do                                                         12d15s19
          jopr=jopr+1                                                   12d15s19
         end do                                                         12d15s19
         ibcoff=itmp3                                                   12d15s19
         ib1=itmp3                                                       12d15s19
         ibu=ib1                                                        5d7s10
         itmpu=itmp1                                                    5d7s10
         do ik=1,nket                                                   5d3s10
          do ib=1,nbra                                                  5d3s10
           iad1=ik-1+nket*(ib-1)                                        5d3s10
           iad2=ib-1+nbra+nbraa*(ik-1+nket)                             5d3s10
           bc(itmpu+iad2)=bc(ibu+iad1)                                  8d20s15
          end do                                                        5d7s10
         end do                                                         5d3s10
        end if                                                          5d3s10
       end if                                                           5d3s10
       if(nketa.ne.nket)then
        itmpu=itmp1                                                     5d7s10
        if(idwsdeb.gt.1)then
         write(6,*)('b4 folding ket ')
         call prntm2(bc(itmpu),nbraa,nketa,nbraa)
        end if
        do ik=1,nket                                                    5d3s10
         do ib=1,nbraa                                                  5d3s10
          iad1=ib-1+nbraa*(ik-1)                                        5d3s10
          iad2=iad1+nbraa*nket                                          5d3s10
          sum=srh*(bc(itmpu+iad1)+bc(itmpu+iad2))                       5d3s10
          dif=srh*(-bc(itmpu+iad1)+bc(itmpu+iad2))                      5d3s10
          bc(itmpu+iad1)=sum                                            5d3s10
          bc(itmpu+iad2)=dif                                            5d3s10
         end do                                                         5d7s10
        end do                                                          5d3s10
        if(idwsdeb.gt.1)then
         write(6,*)('after folding ket ')
         call prntm2(bc(itmpu),nbraa,nketa,nbraa)
        end if
       end if                                                           5d3s10
       if(nbraa.ne.nbra)then
        itmpu=itmp1                                                     5d7s10
        if(idwsdeb.ne.0)then
         write(6,*)('b4 folding bra ')
         call prntm2(bc(itmpu),nbraa,nketa,nbraa)
        end if
        do ik=1,nketa                                                   5d3s10
         do ib=1,nbra                                                   5d3s10
          iad1=ib-1+nbraa*(ik-1)                                        5d3s10
          iad2=iad1+nbra                                                5d3s10
          sum=srh*(bc(itmpu+iad1)+bc(itmpu+iad2))                       5d3s10
          dif=srh*(-bc(itmpu+iad1)+bc(itmpu+iad2))                      5d3s10
          bc(itmpu+iad1)=sum                                            5d3s10
          bc(itmpu+iad2)=dif                                            5d3s10
         end do                                                         5d7s10
        end do                                                          5d3s10
        if(idwsdeb.ne.0)then
         write(6,*)('after folding bra ')
         call prntm2(bc(itmpu),nbraa,nketa,nbraa)
        end if
       end if                                                           5d3s10
       itmpu=itmp1                                                      5d7s10
       do ik=1,nketa                                                    5d3s10
        ikk=ik+ibc(jket4)                                               5d3s10
        isk=multh(isstor(ikk),iopso(1,ioprso))                          2d18s22
        do ib=1,nbraa                                                   5d3s10
         ibb=ib+ibc(jbra4)                                              5d3s10
         iad1=ib-1+nbraa*(ik-1)                                         5d3s10
         if(isk.eq.isstor(ibb))then                                     5d3s10
          iad=ixmt(isstor(ibb),ioprso)+ibstor(ibb)-1                    2d22s22
     $          +nbasisp(isstor(ibb))*(ibstor(ikk)-1)                   2d16s22
          bc(iad)=bc(iad)+bc(itmpu+iad1)                                12d15s19
         end if                                                         12d16s19
        end do                                                          5d7s10
       end do                                                           5d3s10
      end do                                                            2d19s10
      call dws_sync                                                     2d22s10
      call dws_gsumf(bc(ixmt(1,1)),nbb)                                 12d15s19
c
c     now transform to mo basis
c
      do ioprso=1,nopso                                                 2d16s22
       iopr=iopso(2,ioprso)                                             2d16s22
       if(idwsdeb.gt.10)
     $      write(6,*)('for operator no. '),iopr,opname(iopr),('type '),
     $      iopso(3,ioprso)
       sum=0d0                                                          2d18s22
       do isb=1,nsymb                                                   12d15s19
        if(idwsdeb.gt.10)write(6,*)('for bra symmetry '),isb
        isk=multh(isb,iopso(1,ioprso))                                  2d16s22
        if(idwsdeb.gt.10)write(6,*)('for ket symmetry '),isk
        if(min(nbasisp(isb),nbasdwx(isb),nbasisp(isk),nbasdwx(isk))     12d16s19
     $       .gt.0)then                                                 12d16s19
         itmp=ibcoff                                                     12d15s19
         ibcoff=itmp+ncomp*nbasisp(isb)*nbasdwx(isk)                    2d18s22
         call enough('parap42.  8',bc,ibc)
         iorbk=iorb(isk)+nbasisp(isk)                                   2d16s22
         call dgemm('n','n',nbasisp(isb),nbasdwx(isk),nbasisp(isk),
     $       1d0,bc(ixmt(isb,ioprso)),nbasisp(isb),bc(iorbk),           2d18s22
     $       nbasisp(isk)*2,0d0,bc(itmp),nbasisp(isb),                  2d16s22
     d' parap42.  1')
         do i=0,nbasdwx(isk)-1                                           12d15s19
          do j=0,nbasisp(isb)-1                                         1d11s20
           ji=itmp+j+nbasisp(isb)*i                                     1d11s20
           ij=ixmt(isb,ioprso)+i+nbasdwx(isk)*j                            12d15s19
           bc(ij)=bc(ji)                                                 12d15s19
          end do                                                         12d15s19
         end do                                                          12d15s19
         iorbb=iorb(isb)+nbasisp(isb)                                   2d16s22
         call dgemm('n','n',nbasdwx(isk),nbasdwx(isb),nbasisp(isb),1d0,
     $        bc(ixmt(isb,ioprso)),nbasdwx(isk),bc(iorbb),              2d18s22
     $       nbasisp(isb)*ncomp,0d0,bc(itmp),nbasdwx(isk),              1d11s20
     d' parap42.  2')
         do i=0,nbasdwx(isb)-1                                           12d15s19
          do j=0,nbasdwx(isk)-1                                          12d15s19
           ji=itmp+j+nbasdwx(isk)*i                                      12d15s19
           ij=ixmt(isb,ioprso)+i+nbasdwx(isb)*j                            12d15s19
           bc(ij)=bc(ji)                                                 12d15s19
          end do                                                         12d15s19
         end do                                                          12d15s19
         if(isb.eq.isk)then                                             2d18s22
          do i=0,idoubo(isb)-1                                          2d18s22
           ii=ixmt(isb,ioprso)+i*(nbasdwx(isb)+1)                       2d18s22
           sum=sum+2d0*bc(ii)                                           2d18s22
          end do                                                        2d18s22
         end if                                                         2d18s22
         nnnb=nbasdwx(isb)-idoubo(isb)                                  2d18s22
         do i=idoubo(isk),nbasdwx(isk)-1                                2d18s22
          im=i-idoubo(isk)                                              2d18s22
          iad1=ixmt(isb,ioprso)+nbasdwx(isb)*i                          2d18s22
          iad2=itmp-idoubo(isb)+nnnb*im                                 2d18s22
          do j=idoubo(isb),nbasdwx(isb)-1                               2d18s22
           bc(iad2+j)=bc(iad1+j)                                        2d18s22
          end do                                                        2d18s22
         end do                                                         2d18s22
         nnnk=nbasdwx(isk)-idoubo(isk)                                  2d18s22
         do i=0,nnnb*nnnk-1                                             2d18s22
          bc(ixmt(isb,ioprso)+i)=bc(itmp+i)                             2d18s22
         end do                                                         2d18s22
         if(idwsdeb.gt.10)then                                          3d2s22
          write(6,*)('w/o doubles ')                                    3d2s22
          call prntm2(bc(ixmt(isb,ioprso)),nnnb,nnnk,nnnb)              3d2s22
         end if                                                         3d2s22
         ibcoff=itmp                                                     12d15s19
        end if                                                          12d16s19
       end do                                                           12d15s19
       if(abs(sum).gt.1d-10.and.lprint)then                             3d2s22
        write(6,*)('shift sum: '),sum,opname(iopr),('type '),           3d2s22
     $      iopso(3,ioprso)                                             3d2s22
        write(6,*)('but I''m throwing this away !!!')                   3d2s22
       end if                                                           3d2s22
      end do                                                            12d15s19
      ibcoff=ibcoffo                                                    2d19s10
      return
      end                                                               2d19s10
