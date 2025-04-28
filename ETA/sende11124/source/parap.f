c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine parap(natom,ngaus,ibdat,ixmt,isym,iapair,              12d15s19
     $                  ibstor,isstor,iso,idwsdeb,idorel,ascale,        12d15s19
     $     ioppt,npt,opdata,iopdata,iosym,nop,multh,nbb,nbasisp,nbasdwx,12d16s19
     $     iorb,opname,i2comp,bc,ibc)                                   11d9s22
      implicit real*8 (a-h,o-z)
      external second
      integer*8 nn8                                                     2d19s10
      include "common.hf"
      include "common.store"
      include "common.spher"
c
c     If we are working with the 2-component basis, lz etc should act
c     the same on the two components. In this case, i2comp=1.
c     If we are working with the 4-component basis, operators act on
c     the differentiated second function. In this case, i2comp=0.
c
      logical log1                                                      6d6s18
      character*(*) opname(*)                                           3d31s23
      integer*8 ibstor,isstor                                           5d6s10
      include "common.basis"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension ixmt(8,*),isym(3,*),iapair(3,*),ibstor(*),isstor(*),    12d16s19
     $     iso(3),cartb(3),cartk(3),ioppt(*),npt(*),opdata(*),          12d16s19
     $     iopdata(7,*),iosym(*),multh(8,8),nbasisp(8),iorb(8),         5d27s21
     $     nbasdwx(8),ixyz(3,3),xcartb(3),xcartk(3)                     5d27s21
      data ixyz/1,0,0,0,1,0,0,0,1/                                      1d11s20
      ibcoffo=ibcoff                                                    2d19s10
      ncomp=1                                                           1d11s20
      if(idorel.ne.0)then                                               1d8s19
       ncomp=2
      end if                                                            1d8s19
      if(idwsdeb.gt.100)then                                            5d25s18
       write(6,*)('in parap ')
       write(6,*)('i2comp = '),i2comp
       write(6,*)('ncomp '),ncomp,idorel
       write(6,*)('nbb = '),nbb
       write(6,*)('iapair ')
       do i=1,natom
        write(6,*)i,(iapair(j,i),j=1,3)
       end do
      write(6,*)('ibdat '),ibdat
      call prntm2(bc(ibdat),ngaus,8,ngaus)
      write(6,*)('bstor: '),(ibstor(i),i=1,20)
      write(6,*)('sstor: '),(ibstor(i),i=1,20)
      end if                                                            5d25s18
      call second(time1)
c
c     build shell pair order
c
      srh=sqrt(0.5d0)                                                   5d3s10
      ascale2=ascale*2d0                                                8d20s15
      nn=ngaus*ngaus                                                    12d16s19
      nn2=nn*2                                                          2d19s10
c
c     properties:
c     there will be nop of them, with rules given by ipt,npt,opdata,
c     iopdata, and iosym.
c
      if(idorel.eq.0)then                                               8d20s15
       ndo=nop                                                          12d15s19
      else                                                              8d20s15
c
c     px prop px, py prop py, pz prop pz
c
       nbasall2=nbasall*2                                               8d20s15
       ndo=nop                                                          1d11s20
      end if                                                            8d20s15
      ipair=ibcoff                                                      2d19s10
      ibcoff=ipair+nn*ndo*3                                             2d19s10
      i12=ibcoff                                                        2d19s10
      ibcoff=i12+nn*4                                                   2d19s10
      call enough('parap.  1',bc,ibc)
      j12=i12                                                           2d19s10
      jbdat=ibdat-1                                                     2d19s10
      ngaus7=ngaus*7                                                    5d3s10
      do i1=1,nsymb
 5515  format(i5,5x,3i5)
      end do
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
      if(idwsdeb.gt.100)then                                            5d25s18
       write(6,*)('nn8: '),nn8
      end if                                                            5d25s18
      call idsortdws(ibc(i12),ibc(is12),nn8)                            1d18s23
      ii1=i12+nn                                                        2d19s10
      ii2=i12+nn2                                                       2d19s10
      jpair=ipair                                                       2d19s10
      do i=1,nn                                                         2d19s10
       j=ibc(is12+i-1)-1                                                2d19s10
    1  format(i5,i8,2i5)                                                2d19s10
       do k=1,ndo                                                       2d19s10
        jpair0=jpair
        ibc(jpair)=ibc(ii1+j)                                           2d19s10
        jpair=jpair+1                                                   2d19s10
        ibc(jpair)=ibc(ii2+j)                                           2d19s10
        jpair=jpair+1                                                   2d19s10
        ibc(jpair)=k                                                    12d16s19
        jpair=jpair+1                                                   2d19s10
       end do                                                           2d19s10
      end do                                                            2d19s10
      nneed=nn*ndo                                                      2d19s10
      call second(time2)
      telap=time2-time1
      loopit=0                                                          5d25s18
      do i=1+mynowprog,nneed,mynprocg                                   2d19s10
       loopit=loopit+1                                                  5d25s18
       jpair=ipair+3*(i-1)                                              2d19s10
       if(idwsdeb.gt.10)
     $      write(6,2)i,ibc(jpair),ibc(jpair+1),ibc(jpair+2),time1                 2d19s10
    2  format('I am going to do ',i5,3i3,f18.5)                               2d19s10
       jbra=jbdat+ibc(jpair)                                            2d19s10
       jbra2=jbra+ngaus                                                 2d19s10
       jbra3=jbra2+ngaus                                                2d19s10
       jbra4=jbra3+ngaus                                                2d19s10
       jbra5=jbra4+ngaus                                                2d19s10
       jbra6=jbra5+ngaus                                                2d19s10
       jbra7=jbra6+ngaus                                                2d19s10
       jbra8=jbra7+ngaus                                                5d4s10
       do ia=1,natom                                                    5d27s21
        test=sqrt((bc(jbra5)-xcart(1,ia))**2+(bc(jbra6)-xcart(2,ia))**2 5d27s21
     $           +(bc(jbra7)-xcart(3,ia))**2)                           5d27s21
        if(test.lt.1d-10)then                                           5d27s21
         iabra=ia                                                       5d27s21
         go to 1001                                                     5d27s21
        end if                                                          5d27s21
       end do                                                           5d27s21
       write(6,*)('could not set iabra! ')
       write(6,*)('checking '),bc(jbra5),bc(jbra6),bc(jbra7)
       write(6,*)('against ')
       do ia=1,natom
        write(6,*)ia,xcart(1,ia),xcart(2,ia),xcart(3,ia)
       end do
       call dws_synca
       stop 'parap'
 1001  continue                                                         5d27s21
       jket=jbdat+ibc(jpair+1)                                          2d19s10
       jket2=jket+ngaus                                                 2d19s10
       jket3=jket2+ngaus                                                2d19s10
       jket4=jket3+ngaus                                                2d19s10
       jket5=jket4+ngaus                                                2d19s10
       jket6=jket5+ngaus                                                2d19s10
       jket7=jket6+ngaus                                                2d19s10
       jket8=jket7+ngaus                                                5d4s10
       do ia=1,natom                                                    5d27s21
        test=sqrt((bc(jket5)-xcart(1,ia))**2+(bc(jket6)-xcart(2,ia))**2 5d27s21
     $           +(bc(jket7)-xcart(3,ia))**2)                           5d27s21
        if(test.lt.1d-10)then                                           5d27s21
         iaket=ia                                                       5d27s21
         go to 1002                                                     5d27s21
        end if                                                          5d27s21
       end do                                                           5d27s21
       write(6,*)('could not set iaket! ')
       write(6,*)('checking '),bc(jket5),bc(jket6),bc(jket7)
       write(6,*)('against ')
       do ia=1,natom
        write(6,*)ia,xcart(1,ia),xcart(2,ia),xcart(3,ia)
       end do
       call dws_synca
       stop 'parap'
 1002  continue                                                         5d27s21
        iopr=ibc(jpair+2)                                               12d15s19
        nbra=2*ibc(jbra)+1                                              5d3s10
        nbraa=nbra                                                      5d3s10
        if(iapair(1,ibc(jbra8)).gt.0)nbraa=nbra*2                       5d3s10
        nket=2*ibc(jket)+1                                              5d3s10
        nketa=nket                                                      5d3s10
        if(iapair(1,ibc(jket8)).gt.0)nketa=nket*2                       5d3s10
        itmp1=ibcoff                                                    5d3s10
        itmp2=itmp1+nbraa*nketa                                         5d3s10
        ibcoff=itmp2+nbraa*nketa                                        5d3s10
        if(idorel.ne.0)then                                             1d13s20
         ltmp1=ibcoff                                                   1d11s20
         ltmp2=ltmp1+nbraa*nketa                                        1d11s20
         ibcoff=ltmp2+nbraa*nketa                                       1d11s20
        end if                                                          1d11s20
        call enough('parap.  2',bc,ibc)
        jopr=ioppt(iopr)                                                12d16s19
        if(npt(iopr).gt.1)then                                          12d15s19
         itmp3=ibcoff                                                   12d15s19
         ibcoff=itmp3+nbra*nket                                         12d15s19
         if(idorel.ne.0)then                                            1d13s20
          ltmp3=ibcoff                                                  1d11s20
          ibcoff=ltmp3+nbra*nket                                        1d11s20
         end if                                                         1d11s20
         call enough('parap.  3',bc,ibc)
         do ix=itmp3,ibcoff-1                                            12d15s19
          bc(ix)=0d0                                                     12d15s19
         end do                                                         12d15s19
         do ipart=1,npt(iopr)                                           12d15s19
          if(iopdata(7,jopr).eq.0.or.iabra.ne.iaket)then                2d1s22
           xcartb(1)=bc(jbra5)                                          5d27s21
           xcartb(2)=bc(jbra6)                                          5d27s21
           xcartb(3)=bc(jbra7)                                          5d27s21
           xcartk(1)=bc(jket5)                                          5d27s21
           xcartk(2)=bc(jket6)                                          5d27s21
           xcartk(3)=bc(jket7)                                          5d27s21
          else                                                          5d27s21
           do jxyz=1,3                                                  5d27s21
            xcartb(jxyz)=xcart(jxyz,iabra)-cmx(jxyz,iagrp(iabra))       5d27s21
            xcartk(jxyz)=xcart(jxyz,iaket)-cmx(jxyz,iagrp(iaket))       5d27s21
           end do                                                       5d27s21
          end if                                                        5d27s21
          call onep(ibc(jbra),bc(jbra2),bc(jbra3),xcartb,xcartb(2),     5d27s21
     $      xcartb(3),ibc(jket),bc(jket2),bc(jket3),                    12d15s19
     $      xcartk,xcartk(2),xcartk(3),nsymb,ixmat,0,0,0,               5d27s21
     $      iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),            12d15s19
     $      iopdata(4,jopr),iopdata(5,jopr),iopdata(6,jopr),bc,ibc)     11d9s22
          do ix=0,nbra*nket-1                                            12d15s19
           bc(itmp3+ix)=bc(itmp3+ix)+opdata(jopr)*bc(ixmat+ix)             12d15s19
          end do                                                        12d15s19
          if(idorel.ne.0)then                                           1d13s20
           if(i2comp.eq.0)then                                          1d13s20
            factx=opdata(jopr)*ascale
           else                                                         1d13s20
            factx=-opdata(jopr)*ascale
           end if                                                       1d13s20
           do jxyz=1,3                                                  1d11s20
            if(i2comp.eq.0)then                                         1d13s20
             jxp=ixyz(1,jxyz)                                           1d13s20
             jyp=ixyz(2,jxyz)                                           1d13s20
             jzp=ixyz(3,jxyz)                                           1d13s20
             ixp=iopdata(4,jopr)+ixyz(1,jxyz)                            1d11s20
             iyp=iopdata(5,jopr)+ixyz(2,jxyz)                            1d11s20
             izp=iopdata(6,jopr)+ixyz(3,jxyz)                            1d11s20
            else                                                        1d13s20
             jxp=2*ixyz(1,jxyz)                                           1d13s20
             jyp=2*ixyz(2,jxyz)                                           1d13s20
             jzp=2*ixyz(3,jxyz)                                           1d13s20
             ixp=iopdata(4,jopr)                                        1d13s20
             iyp=iopdata(5,jopr)                                        1d13s20
             izp=iopdata(6,jopr)                                        1d13s20
            end if                                                      1d13s20
            call onep(ibc(jbra),bc(jbra2),bc(jbra3),xcartb,xcartb(2),   5d27s21
     $           xcartb(3),ibc(jket),bc(jket2),bc(jket3),               5d27s21
     $      xcartk,xcartk(2),xcartk(3),nsymb,ixmat,jxp,jyp,jzp,         5d27s21
     $      iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),ixp,iyp,izp,11d9s22
     $           bc,ibc)                                                11d9s22
            do ix=0,nbra*nket-1                                            12d15s19
             bc(ltmp3+ix)=bc(ltmp3+ix)+factx*bc(ixmat+ix)               1d11s20
            end do                                                        12d15s19
           end do                                                       1d11s20
          end if                                                        1d11s20
          jopr=jopr+1                                                   12d15s19
         end do                                                         12d15s19
         ibcoff=itmp3                                                   12d15s19
         ib1=itmp3                                                       12d15s19
        else                                                             12d15s19
         call onep(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),    12d15s19
     $     bc(jbra7),ibc(jket),bc(jket2),bc(jket3),                     12d15s19
     $     bc(jket5),bc(jket6),bc(jket7),nsymb,ib1,0,0,0,               12d15s19
     $     iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),             12d15s19
     $     iopdata(4,jopr),iopdata(5,jopr),iopdata(6,jopr),bc,ibc)      11d9s22
         do ix=0,nbra*nket-1                                              12d15s19
          bc(ib1+ix)=bc(ib1+ix)*opdata(jopr)                               12d15s19
         end do                                                          12d15s19
         ibcoff=ib1+nbra*nket                                           3d3s20
         if(idorel.ne.0)then                                            1d13s20
          lb1=ib1+nbra*nket                                             1d11s20
          ibcoff=lb1+nbra*nket                                          1d11s20
          call enough('parap.  4',bc,ibc)
          do ii=lb1,ibcoff-1                                             3d3s20
           bc(ii)=0d0                                                    3d3s20
          end do                                                        3d3s20
          if(i2comp.eq.0)then                                           1d13s20
           factx=opdata(jopr)*ascale
          else                                                          1d13s20
           factx=-opdata(jopr)*ascale
          end if                                                        1d13s20
          do jxyz=1,3                                                   1d11s20
           if(i2comp.eq.0)then                                          1d13s20
            jxp=ixyz(1,jxyz)                                            1d13s20
            jyp=ixyz(2,jxyz)                                            1d13s20
            jzp=ixyz(3,jxyz)                                            1d13s20
            ixp=iopdata(4,jopr)+ixyz(1,jxyz)                             1d11s20
            iyp=iopdata(5,jopr)+ixyz(2,jxyz)                             1d11s20
            izp=iopdata(6,jopr)+ixyz(3,jxyz)                             1d11s20
           else                                                         1d13s20
            jxp=2*ixyz(1,jxyz)                                            1d13s20
            jyp=2*ixyz(2,jxyz)                                            1d13s20
            jzp=2*ixyz(3,jxyz)                                            1d13s20
            ixp=iopdata(4,jopr)                                         1d13s20
            iyp=iopdata(5,jopr)                                         1d13s20
            izp=iopdata(6,jopr)                                         1d13s20
           end if                                                       1d13s20
           call onep(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),  12d15s19
     $           bc(jbra7),ibc(jket),bc(jket2),bc(jket3),                    12d15s19
     $     bc(jket5),bc(jket6),bc(jket7),nsymb,ixmat,jxp,jyp,jzp,       1d13s20
     $     iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),ixp,iyp,izp, 11d9s22
     $          bc,ibc)                                                 11d9s22
           do ix=0,nbra*nket-1                                            12d15s19
            bc(lb1+ix)=bc(lb1+ix)+factx*bc(ixmat+ix)                    1d11s20
           end do                                                        12d15s19
          end do                                                        1d11s20
         ibcoff=ib1
         end if                                                         1d11s20
        end if                                                           12d15s19
       if(nsymb.ne.-1)then                                               5d3s10
        ibu=ib1                                                         5d7s10
        itmpu=itmp1                                                     5d7s10
        do ik=1,nket                                                    12d15s19
         do ib=1,nbra                                                   12d15s19
          iad1=ik-1+nket*(ib-1)                                         12d15s19
          iad2=ib-1+nbraa*(ik-1)                                        12d15s19
          bc(itmpu+iad2)=bc(ibu+iad1)                                   12d15s19
         end do                                                         12d15s19
        end do                                                          12d15s19
        if(idorel.ne.0)then                                             1d13s20
         ibu=ib1+nbra*nket                                              3d3s20
         itmpu=ltmp1                                                     5d7s10
         do ik=1,nket                                                    12d15s19
          do ib=1,nbra                                                   12d15s19
           iad1=ik-1+nket*(ib-1)                                         12d15s19
           iad2=ib-1+nbraa*(ik-1)                                        12d15s19
           bc(itmpu+iad2)=bc(ibu+iad1)                                   12d15s19
          end do                                                         12d15s19
         end do                                                          12d15s19
        end if                                                          1d11s20
        if(nketa.ne.nket)then                                           5d3s10
         cartk(1)=bc(jket5)*dfloat(isym(1,iapair(2,ibc(jket8))))        5d5s10
         cartk(2)=bc(jket6)*dfloat(isym(2,iapair(2,ibc(jket8))))        5d5s10
         cartk(3)=bc(jket7)*dfloat(isym(3,iapair(2,ibc(jket8))))        5d5s10
         do ia=1,natom                                                    5d27s21
          test=sqrt((cartk(1)-xcart(1,ia))**2+(cartk(2)-xcart(2,ia))**2 5d27s21
     $           +(cartk(3)-xcart(3,ia))**2)                            5d27s21
          if(test.lt.1d-10)then                                           5d27s21
           iakets=ia                                                       5d27s21
           go to 1003                                                     5d27s21
          end if                                                          5d27s21
         end do                                                           5d27s21
 1003    continue                                                         5d27s21
         iopr=ibc(jpair+2)                                               12d15s19
         jopr=ioppt(iopr)                                               12d16s19
         if(npt(iopr).gt.1)then                                          12d15s19
          itmp3=ibcoff                                                   12d15s19
          ibcoff=itmp3+nbra*nket                                         12d15s19
          if(idorel.ne.0)then                                           1d13s20
           ltmp3=ibcoff                                                   12d15s19
           ibcoff=ltmp3+nbra*nket                                         12d15s19
          end if                                                        1d11s20
          call enough('parap.  5',bc,ibc)
          do ix=itmp3,ibcoff-1                                            12d15s19
           bc(ix)=0d0                                                     12d15s19
          end do                                                         12d15s19
          do ipart=1,npt(iopr)                                           12d15s19
           if(iopdata(7,jopr).eq.0.or.iabra.ne.iakets)then              2d1s22
            call onep(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),  12d15s19
     $          bc(jbra7),ibc(jket),bc(jket2),bc(jket3),                    12d15s19
     $      cartk,cartk(2),cartk(3),nsymb,ixmat,0,0,0,                  12d15s19
     $      iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),            12d15s19
     $      iopdata(4,jopr),iopdata(5,jopr),iopdata(6,jopr),bc,ibc)     11d9s22
           else                                                         5d27s21
            do jxyz=1,3                                                  5d27s21
             xcartb(jxyz)=xcart(jxyz,iabra)-cmx(jxyz,iagrp(iabra))       5d27s21
             xcartk(jxyz)=xcart(jxyz,iakets)-cmx(jxyz,iagrp(iakets))    5d27s21
            end do                                                       5d27s21
            call onep(ibc(jbra),bc(jbra2),bc(jbra3),xcartb,xcartb(2),   5d27s21
     $          xcartb(3),ibc(jket),bc(jket2),bc(jket3),                5d27s21
     $      xcartk,xcartk(2),xcartk(3),nsymb,ixmat,0,0,0,               5d27s21
     $      iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),            12d15s19
     $      iopdata(4,jopr),iopdata(5,jopr),iopdata(6,jopr),bc,ibc)     11d9s22
           end if                                                       5d27s21
           do ix=0,nbra*nket-1                                            12d15s19
            bc(itmp3+ix)=bc(itmp3+ix)+opdata(jopr)*bc(ixmat+ix)             12d15s19
           end do                                                        12d15s19
           if(idorel.ne.0)then                                          1d13s20
            if(i2comp.eq.0)then                                         1d13s20
             factx=opdata(jopr)*ascale
            else                                                        1d13s20
             factx=-opdata(jopr)*ascale
            end if                                                      1d13s20
            do jxyz=1,3                                                  1d11s20
             if(i2comp.eq.0)then                                        1d13s20
              jxp=ixyz(1,jxyz)                                          1d13s20
              jyp=ixyz(2,jxyz)                                          1d13s20
              jzp=ixyz(3,jxyz)                                          1d13s20
              ixp=iopdata(4,jopr)+ixyz(1,jxyz)                            1d11s20
              iyp=iopdata(5,jopr)+ixyz(2,jxyz)                            1d11s20
              izp=iopdata(6,jopr)+ixyz(3,jxyz)                            1d11s20
             else                                                       1d13s20
              jxp=2*ixyz(1,jxyz)                                          1d13s20
              jyp=2*ixyz(2,jxyz)                                          1d13s20
              jzp=2*ixyz(3,jxyz)                                          1d13s20
              ixp=iopdata(4,jopr)                                       1d13s20
              iyp=iopdata(5,jopr)                                       1d13s20
              izp=iopdata(6,jopr)                                       1d13s20
             end if                                                     1d13s20
             if(iopdata(7,jopr).eq.0)then                               5d27s21
              call onep(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),         1d11s20
     $            bc(jbra6),bc(jbra7),ibc(jket),bc(jket2),bc(jket3),    1d11s20
     $      cartk,cartk(2),cartk(3),nsymb,ixmat,jxp,jyp,jzp,            10d31s20
     $      iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),ixp,iyp,izp,11d9s22
     $            bc,ibc)                                               11d9s22
             else                                                       5d27s21
              call onep(ibc(jbra),bc(jbra2),bc(jbra3),xcartb,           5d27s21
     $            xcartb(2),xcartb(3),ibc(jket),bc(jket2),bc(jket3),    5d27s21
     $      xcartk,xcartk(2),xcartk(3),nsymb,ixmat,jxp,jyp,jzp,            10d31s20
     $      iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),ixp,iyp,izp,11d9s22
     $             bc,ibc)                                              11d9s22
             end if                                                     5d27s21
             do ix=0,nbra*nket-1                                            12d15s19
              bc(ltmp3+ix)=bc(ltmp3+ix)+factx*bc(ixmat+ix)               1d11s20
             end do                                                        12d15s19
            end do                                                       1d11s20
           end if                                                        1d11s20
           jopr=jopr+1                                                   12d15s19
          end do                                                         12d15s19
          ibcoff=itmp3                                                   12d15s19
          ib1=itmp3                                                       12d15s19
         else                                                             12d15s19
          call onep(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),    12d15s19
     $     bc(jbra7),ibc(jket),bc(jket2),bc(jket3),                     12d15s19
     $     cartk,cartk(2),cartk(3),nsymb,ib1,0,0,0,                     12d15s19
     $     iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),             12d15s19
     $     iopdata(4,jopr),iopdata(5,jopr),iopdata(6,jopr),bc,ibc)      11d9s22
          do ix=0,nbra*nket-1                                              12d15s19
           bc(ib1+ix)=bc(ib1+ix)*opdata(jopr)                               12d15s19
          end do                                                          12d15s19
          if(idorel.ne.0)then                                           1d13s20
           lb1=ib1+nbra*nket                                             1d11s20
           ibcoff=lb1+nbra*nket                                          1d11s20
           call enough('parap.  6',bc,ibc)
           do ii=lb1,ibcoff-1                                           3d3s20
            bc(ii)=0d0                                                  3d3s20
           end do                                                       3d3s20
           factx=opdata(jopr)*ascale
           if(i2comp.ne.0)factx=-factx                                  1d13s20
           do jxyz=1,3                                                   1d11s20
            if(i2comp.eq.0)then                                         1d13s20
             jxp=ixyz(1,jxyz)                                           1d13s20
             jyp=ixyz(2,jxyz)                                           1d13s20
             jzp=ixyz(3,jxyz)                                           1d13s20
             ixp=iopdata(4,jopr)+ixyz(1,jxyz)                             1d11s20
             iyp=iopdata(5,jopr)+ixyz(2,jxyz)                             1d11s20
             izp=iopdata(6,jopr)+ixyz(3,jxyz)                             1d11s20
            else                                                        1d13s20
             jxp=2*ixyz(1,jxyz)                                           1d13s20
             jyp=2*ixyz(2,jxyz)                                           1d13s20
             jzp=2*ixyz(3,jxyz)                                           1d13s20
             ixp=iopdata(4,jopr)                                        1d13s20
             iyp=iopdata(5,jopr)                                        1d13s20
             izp=iopdata(6,jopr)                                        1d13s20
            end if                                                      1d13s20
            call onep(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),  12d15s19
     $           bc(jbra7),ibc(jket),bc(jket2),bc(jket3),                    12d15s19
     $     cartk,cartk(2),cartk(3),nsymb,ixmat,jxp,jyp,jzp,             10d31s20
     $     iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),ixp,iyp,izp, 11d9s22
     $           bc,ibc)                                                11d9s22
            do ix=0,nbra*nket-1                                            12d15s19
             bc(lb1+ix)=bc(lb1+ix)+factx*bc(ixmat+ix)                    1d11s20
            end do                                                        12d15s19
           end do                                                        1d11s20
           ibcoff=ib1                                                   3d3s20
          end if                                                         1d11s20
         end if                                                           12d15s19
         ibu=ib1                                                        5d7s10
         itmpu=itmp1                                                    5d7s10
         do ik=1,nket                                                    5d3s10
          do ib=1,nbra                                                   5d3s10
           iad1=ik-1+nket*(ib-1)                                         5d3s10
           iad2=ib-1+nbraa*(ik-1+nket)                                  5d3s10
           bc(itmpu+iad2)=bc(ibu+iad1)                                  8d20s15
          end do                                                         5d3s10
         end do                                                          5d3s10
         if(idorel.ne.0)then                                            1d13s20
          ibu=ib1+nket*nbra                                             1d11s20
          itmpu=ltmp1                                                    5d7s10
          do ik=1,nket                                                    5d3s10
           do ib=1,nbra                                                   5d3s10
            iad1=ik-1+nket*(ib-1)                                         5d3s10
            iad2=ib-1+nbraa*(ik-1+nket)                                  5d3s10
            bc(itmpu+iad2)=bc(ibu+iad1)                                  8d20s15
           end do                                                         5d3s10
          end do                                                          5d3s10
         end if                                                         1d11s20
        end if                                                          5d3s10
        if(nbraa.ne.nbra)then                                           5d3s10
         cartb(1)=bc(jbra5)*dfloat(isym(1,iapair(2,ibc(jbra8))))        5d5s10
         cartb(2)=bc(jbra6)*dfloat(isym(2,iapair(2,ibc(jbra8))))        5d5s10
         cartb(3)=bc(jbra7)*dfloat(isym(3,iapair(2,ibc(jbra8))))        5d5s10
         do ia=1,natom                                                    5d27s21
          test=sqrt((cartb(1)-xcart(1,ia))**2+(cartb(2)-xcart(2,ia))**2 5d27s21
     $           +(cartb(3)-xcart(3,ia))**2)                            5d27s21
          if(test.lt.1d-10)then                                           5d27s21
           iabras=ia                                                    5d27s21
           go to 1004                                                   5d27s21
          end if                                                          5d27s21
         end do                                                           5d27s21
 1004    continue                                                         5d27s21
         iopr=ibc(jpair+2)                                               12d15s19
         jopr=ioppt(iopr)                                               12d16s19
         if(npt(iopr).gt.1)then                                          12d15s19
          itmp3=ibcoff                                                   12d15s19
          ibcoff=itmp3+nbra*nket                                         12d15s19
          if(idorel.ne.0)then                                           1d13s20
           ltmp3=ibcoff                                                   12d15s19
           ibcoff=ltmp3+nbra*nket                                         12d15s19
          end if
          call enough('parap.  7',bc,ibc)
          do ix=itmp3,ibcoff-1                                            12d15s19
           bc(ix)=0d0                                                     12d15s19
          end do                                                         12d15s19
          do ipart=1,npt(iopr)                                           12d15s19
           if(iopdata(7,jopr).eq.0.or.iabras.ne.iaket)then              2d1s22
            call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),       12d15s19
     $      cartb(3),ibc(jket),bc(jket2),bc(jket3),                     12d15s19
     $      bc(jket5),bc(jket6),bc(jket7),nsymb,ixmat,0,0,0,
     $      iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),            12d15s19
     $      iopdata(4,jopr),iopdata(5,jopr),iopdata(6,jopr),bc,ibc)     11d9s22
           else                                                         5d27s21
            do jxyz=1,3                                                  5d27s21
             xcartb(jxyz)=xcart(jxyz,iabras)-cmx(jxyz,iagrp(iabras))    5d27s21
             xcartk(jxyz)=xcart(jxyz,iaket)-cmx(jxyz,iagrp(iaket))      5d27s21
            end do                                                       5d27s21
            call onep(ibc(jbra),bc(jbra2),bc(jbra3),xcartb,xcartb(2),   5d27s21
     $      xcartb(3),ibc(jket),bc(jket2),bc(jket3),                     12d15s19
     $      xcartk,xcartk(2),xcartk(3),nsymb,ixmat,0,0,0,               5d27s21
     $      iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),            12d15s19
     $      iopdata(4,jopr),iopdata(5,jopr),iopdata(6,jopr),bc,ibc)     11d9s22
           end if                                                       5d27s21
           do ix=0,nbra*nket-1                                            12d15s19
            bc(itmp3+ix)=bc(itmp3+ix)+opdata(jopr)*bc(ixmat+ix)             12d15s19
           end do                                                        12d15s19
           if(idorel.ne.0)then                                          1d13s20
            factx=opdata(jopr)*ascale
            if(i2comp.ne.0)factx=-factx                                 1d13s20
            do jxyz=1,3                                                  1d11s20
             if(i2comp.eq.0)then                                        1d13s20
              jxp=ixyz(1,jxyz)                                          1d13s20
              jyp=ixyz(2,jxyz)                                          1d13s20
              jzp=ixyz(3,jxyz)                                          1d13s20
              ixp=iopdata(4,jopr)+ixyz(1,jxyz)                            1d11s20
              iyp=iopdata(5,jopr)+ixyz(2,jxyz)                            1d11s20
              izp=iopdata(6,jopr)+ixyz(3,jxyz)                            1d11s20
             else                                                       1d13s20
              jxp=2*ixyz(1,jxyz)                                          1d13s20
              jyp=2*ixyz(2,jxyz)                                          1d13s20
              jzp=2*ixyz(3,jxyz)                                          1d13s20
              ixp=iopdata(4,jopr)                                       1d13s20
              iyp=iopdata(5,jopr)                                       1d13s20
              izp=iopdata(6,jopr)                                       1d13s20
             end if                                                     1d13s20
             if(iopdata(7,jopr).eq.0)then                               5d27s21
              call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),    10d31s20
     $      cartb(3),ibc(jket),bc(jket2),bc(jket3),                     10d31s20
     $      bc(jket5),bc(jket6),bc(jket7),nsymb,ixmat,jxp,jyp,jzp,      1d13s20
     $      iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),ixp,iyp,izp,11d9s22
     $            bc,ibc)                                               11d9s22
             else                                                       5d27s21
              call onep(ibc(jbra),bc(jbra2),bc(jbra3),xcartb,xcartb(2), 5d27s21
     $      xcartb(3),ibc(jket),bc(jket2),bc(jket3),                    5d27s21
     $      xcartk,xcartk(2),xcartk(3),nsymb,ixmat,jxp,jyp,jzp,         5d27s21
     $      iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),ixp,iyp,izp,11d9s22
     $             bc,ibc)                                              11d9s22
             end if                                                     5d27s21
             do ix=0,nbra*nket-1                                            12d15s19
              bc(ltmp3+ix)=bc(ltmp3+ix)+factx*bc(ixmat+ix)               1d11s20
             end do                                                        12d15s19
            end do                                                       1d11s20
           end if                                                        1d11s20
           jopr=jopr+1                                                   12d15s19
          end do                                                         12d15s19
          ibcoff=itmp3                                                   12d15s19
          ib1=itmp3                                                       12d15s19
         else                                                             12d15s19
          call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),        12d15s19
     $         cartb(3),ibc(jket),bc(jket2),bc(jket3),                      12d15s19
     $     bc(jket5),bc(jket6),bc(jket7),nsymb,ib1,0,0,0,               12d15s19
     $     iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),             12d15s19
     $     iopdata(4,jopr),iopdata(5,jopr),iopdata(6,jopr),bc,ibc)      11d9s22
          do ix=0,nbra*nket-1                                              12d15s19
           bc(ib1+ix)=bc(ib1+ix)*opdata(jopr)                               12d15s19
          end do                                                          12d15s19
          if(idorel.ne.0)then                                           1d13s20
           lb1=ib1+nbra*nket                                             1d11s20
           ibcoff=lb1+nbra*nket                                          1d11s20
           call enough('parap.  8',bc,ibc)
           do ii=lb1,ibcoff-1                                           3d3s20
            bc(ii)=0d0                                                  3d3s20
           end do                                                       3d3s20
           factx=opdata(jopr)*ascale
           if(i2comp.ne.0)factx=-factx                                  1d13s20
           do jxyz=1,3                                                   1d11s20
            if(i2comp.eq.0)then                                         1d13s20
             jxp=ixyz(1,jxyz)                                           1d13s20
             jyp=ixyz(2,jxyz)                                           1d13s20
             jzp=ixyz(3,jxyz)                                           1d13s20
             ixp=iopdata(4,jopr)+ixyz(1,jxyz)                             1d11s20
             iyp=iopdata(5,jopr)+ixyz(2,jxyz)                             1d11s20
             izp=iopdata(6,jopr)+ixyz(3,jxyz)                             1d11s20
            else                                                        1d13s20
             jxp=2*ixyz(1,jxyz)                                           1d13s20
             jyp=2*ixyz(2,jxyz)                                           1d13s20
             jzp=2*ixyz(3,jxyz)                                           1d13s20
             ixp=iopdata(4,jopr)                                        1d13s20
             iyp=iopdata(5,jopr)                                        1d13s20
             izp=iopdata(6,jopr)                                        1d13s20
            end if                                                      1d13s20
            call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),     10d31s20
     $      cartb(3),ibc(jket),bc(jket2),bc(jket3),                     10d31s20
     $     bc(jket5),bc(jket6),bc(jket7),nsymb,ixmat,jxp,jyp,jzp,       1d13s20
     $     iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),ixp,iyp,izp, 11d9s22
     $           bc,ibc)                                                11d9s22
            do ix=0,nbra*nket-1                                            12d15s19
             bc(lb1+ix)=bc(lb1+ix)+factx*bc(ixmat+ix)                    1d11s20
            end do                                                        12d15s19
           end do                                                        1d11s20
           ibcoff=ib1                                                   3d3s20
          end if                                                         1d11s20
         end if                                                           12d15s19
         ibu=ib1                                                         5d7s10
         itmpu=itmp1                                                     5d7s10
         do ik=1,nket                                                    5d3s10
          do ib=1,nbra                                                   5d3s10
           iad1=ik-1+nket*(ib-1)                                         5d3s10
           iad2=ib-1+nbra+nbraa*(ik-1)                                   5d3s10
           bc(itmpu+iad2)=bc(ibu+iad1)                                  8d20s15
          end do                                                        5d7s10
         end do                                                         5d3s10
         if(idorel.ne.0)then                                            1d13s20
          ibu=ib1+nbra*nket                                             1d11s20
          itmpu=ltmp1                                                     5d7s10
          do ik=1,nket                                                    5d3s10
           do ib=1,nbra                                                   5d3s10
            iad1=ik-1+nket*(ib-1)                                         5d3s10
            iad2=ib-1+nbra+nbraa*(ik-1)                                   5d3s10
            bc(itmpu+iad2)=bc(ibu+iad1)                                  8d20s15
           end do                                                        5d7s10
          end do                                                         5d3s10
         end if                                                         1d11s20
         if(nketa.ne.nket)then                                           5d3s10
          iopr=ibc(jpair+2)                                               12d15s19
          jopr=ioppt(iopr)                                              12d16s19
          if(npt(iopr).gt.1)then                                          12d15s19
           itmp3=ibcoff                                                   12d15s19
           ibcoff=itmp3+nbra*nket                                         12d15s19
           if(idorel.ne.0)then                                          1d13s20
            ltmp3=ibcoff                                                1d11s20
            ibcoff=ltmp3+nbra*nket                                      1d11s20
           end if                                                       1d11s20
           call enough('parap.  9',bc,ibc)
           do ix=itmp3,ibcoff-1                                            12d15s19
            bc(ix)=0d0                                                     12d15s19
           end do                                                         12d15s19
           do ipart=1,npt(iopr)                                           12d15s19
            if(iopdata(7,jopr).eq.0.or.iabras.ne.iakets)then            2d1s22
             call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),       12d15s19
     $      cartb(3),ibc(jket),bc(jket2),bc(jket3),                     12d15s19
     $      cartk,cartk(2),cartk(3),nsymb,ixmat,0,0,0,
     $      iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),            12d15s19
     $      iopdata(4,jopr),iopdata(5,jopr),iopdata(6,jopr),bc,ibc)     11d9s22
            else                                                        5d27s21
             do jxyz=1,3                                                  5d27s21
              xcartb(jxyz)=xcart(jxyz,iabras)-cmx(jxyz,iagrp(iabras))    5d27s21
              xcartk(jxyz)=xcart(jxyz,iakets)-cmx(jxyz,iagrp(iakets))      5d27s21
             end do                                                       5d27s21
             call onep(ibc(jbra),bc(jbra2),bc(jbra3),xcartb,xcartb(2),  5d27s21
     $      xcartb(3),ibc(jket),bc(jket2),bc(jket3),                    5d27s21
     $      xcartk,xcartk(2),xcartk(3),nsymb,ixmat,0,0,0,               5d27s21
     $      iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),            12d15s19
     $      iopdata(4,jopr),iopdata(5,jopr),iopdata(6,jopr),bc,ibc)     11d9s22
            end if                                                      5d27s21
            do ix=0,nbra*nket-1                                            12d15s19
             bc(itmp3+ix)=bc(itmp3+ix)+opdata(jopr)*bc(ixmat+ix)             12d15s19
            end do                                                        12d15s19
            if(idorel.ne.0)then                                         1d13s20
             factx=opdata(jopr)*ascale
             if(i2comp.ne.0)factx=-factx                                1d13s20
             do jxyz=1,3                                                  1d11s20
              if(i2comp.eq.0)then                                       1d13s20
               jxp=ixyz(1,jxyz)                                         1d13s20
               jyp=ixyz(2,jxyz)                                         1d13s20
               jzp=ixyz(3,jxyz)                                         1d13s20
               ixp=iopdata(4,jopr)+ixyz(1,jxyz)                            1d11s20
               iyp=iopdata(5,jopr)+ixyz(2,jxyz)                            1d11s20
               izp=iopdata(6,jopr)+ixyz(3,jxyz)                            1d11s20
              else                                                      1d13s20
               jxp=2*ixyz(1,jxyz)                                         1d13s20
               jyp=2*ixyz(2,jxyz)                                         1d13s20
               jzp=2*ixyz(3,jxyz)                                         1d13s20
               ixp=iopdata(4,jopr)                                      1d13s20
               iyp=iopdata(5,jopr)                                      1d13s20
               izp=iopdata(6,jopr)                                      1d13s20
              end if                                                    1d13s20
              if(iopdata(7,jopr).eq.0)then                              5d27s21
               call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),   10d31s20
     $      cartb(3),ibc(jket),bc(jket2),bc(jket3),                     10d31s20
     $      cartk,cartk(2),cartk(3),nsymb,ixmat,jxp,jyp,jzp,            10d31s20
     $      iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),ixp,iyp,izp,11d9s22
     $             bc,ibc)                                              11d9s22
              else                                                      5d27s21
               call onep(ibc(jbra),bc(jbra2),bc(jbra3),xcartb,xcartb(2),5d27s21
     $      xcartb(3),ibc(jket),bc(jket2),bc(jket3),                    5d27s21
     $      xcartk,xcartk(2),xcartk(3),nsymb,ixmat,jxp,jyp,jzp,         5d27s21
     $      iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),ixp,iyp,izp,11d9s22
     $              bc,ibc)                                             11d9s22
              end if                                                    5d27s21
              do ix=0,nbra*nket-1                                            12d15s19
               bc(ltmp3+ix)=bc(ltmp3+ix)+factx*bc(ixmat+ix)               1d11s20
              end do                                                        12d15s19
             end do                                                       1d11s20
            end if                                                        1d11s20
            jopr=jopr+1                                                   12d15s19
           end do                                                         12d15s19
           ibcoff=itmp3                                                   12d15s19
           ib1=itmp3                                                       12d15s19
          else                                                             12d15s19
           call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),        12d15s19
     $         cartb(3),ibc(jket),bc(jket2),bc(jket3),                      12d15s19
     $     cartk,cartk(2),cartk(3),nsymb,ib1,0,0,0,                     12d15s19
     $     iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),             12d15s19
     $     iopdata(4,jopr),iopdata(5,jopr),iopdata(6,jopr),bc,ibc)      11d9s22
           do ix=0,nbra*nket-1                                              12d15s19
            bc(ib1+ix)=bc(ib1+ix)*opdata(jopr)                               12d15s19
           end do                                                          12d15s19
           if(idorel.ne.0)then                                          1d13s20
            lb1=ib1+nbra*nket                                             1d11s20
            ibcoff=lb1+nbra*nket                                          1d11s20
            call enough('parap. 10',bc,ibc)
            do ii=lb1,ibcoff-1                                          3d3s20
             bc(ii)=0d0                                                 3d3s20
            end do                                                      3d3s20
            factx=opdata(jopr)*ascale
            if(i2comp.ne.0)factx=-factx                                 1d13s20
            do jxyz=1,3                                                   1d11s20
             if(i2comp.eq.0)then                                        1d13s20
              jxp=ixyz(1,jxyz)                                          1d13s20
              jyp=ixyz(2,jxyz)                                          1d13s20
              jzp=ixyz(3,jxyz)                                          1d13s20
              ixp=iopdata(4,jopr)+ixyz(1,jxyz)                             1d11s20
              iyp=iopdata(5,jopr)+ixyz(2,jxyz)                             1d11s20
              izp=iopdata(6,jopr)+ixyz(3,jxyz)                             1d11s20
             else                                                       1d13s20
              jxp=2*ixyz(1,jxyz)                                        1d13s20
              jyp=2*ixyz(2,jxyz)                                        1d13s20
              jzp=2*ixyz(3,jxyz)                                        1d13s20
              ixp=iopdata(4,jopr)                                       1d13s20
              iyp=iopdata(5,jopr)                                       1d13s20
              izp=iopdata(6,jopr)                                       1d13s20
             end if                                                     1d13s20
             call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),    10d31s20
     $      cartb(3),ibc(jket),bc(jket2),bc(jket3),                     10d31s20
     $      cartk,cartk(2),cartk(3),nsymb,ixmat,jxp,jyp,jzp,            10d31s20
     $     iopdata(1,jopr),iopdata(2,jopr),iopdata(3,jopr),ixp,iyp,izp, 11d9s22
     $            bc,ibc)                                               11d9s22
             do ix=0,nbra*nket-1                                            12d15s19
              bc(lb1+ix)=bc(lb1+ix)+factx*bc(ixmat+ix)                    1d11s20
             end do                                                        12d15s19
            end do                                                        1d11s20
           end if                                                         1d11s20
          end if                                                           12d15s19
          ibu=ib1                                                        5d7s10
          itmpu=itmp1                                                    5d7s10
          do ik=1,nket                                                  5d3s10
           do ib=1,nbra                                                 5d3s10
            iad1=ik-1+nket*(ib-1)                                       5d3s10
            iad2=ib-1+nbra+nbraa*(ik-1+nket)                            5d3s10
            bc(itmpu+iad2)=bc(ibu+iad1)                                 8d20s15
           end do                                                       5d7s10
          end do                                                        5d3s10
          if(idorel.ne.0)then                                           1d13s20
           itmpu=ltmp1                                                    5d7s10
           ibu=ib1+nbra*nket                                            1d11s20
           do ik=1,nket                                                  5d3s10
            do ib=1,nbra                                                 5d3s10
             iad1=ik-1+nket*(ib-1)                                       5d3s10
             iad2=ib-1+nbra+nbraa*(ik-1+nket)                            5d3s10
             bc(itmpu+iad2)=bc(ibu+iad1)                                 8d20s15
            end do                                                       5d7s10
           end do                                                        5d3s10
          end if                                                        1d11s20
         end if                                                          5d3s10
        end if                                                           5d3s10
        if(nketa.ne.nket)then
         itmpu=itmp1                                                     5d7s10
         if(idwsdeb.gt.1)then
          write(6,*)('b4 folding ket ')
          call prntm2(bc(itmpu),nbraa,nketa,nbraa)
         end if
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
         if(idwsdeb.gt.1)then
          write(6,*)('after folding ket ')
          call prntm2(bc(itmpu),nbraa,nketa,nbraa)
         end if
         if(idorel.ne.0)then                                            1d13s20
          itmpu=ltmp1                                                   1d11s20
          if(idwsdeb.gt.1)then
           write(6,*)('dorel part ')
           write(6,*)('b4 folding ket ')
           call prntm2(bc(itmpu),nbraa,nketa,nbraa)
          end if
          do ik=1,nket                                                  1d11s20
           do ib=1,nbraa                                                1d11s20
            iad1=ib-1+nbraa*(ik-1)                                      1d11s20
            iad2=iad1+nbraa*nket                                        1d11s20
            sum=srh*(bc(itmpu+iad1)+bc(itmpu+iad2))                     1d11s20
            dif=srh*(-bc(itmpu+iad1)+bc(itmpu+iad2))                    1d11s20
            bc(itmpu+iad1)=sum                                          1d11s20
            bc(itmpu+iad2)=dif                                          1d11s20
           end do                                                       1d11s20
          end do                                                        1d11s20
          if(idwsdeb.gt.1)then
           write(6,*)('after folding ket ')
           call prntm2(bc(itmpu),nbraa,nketa,nbraa)
          end if
         end if                                                         1d11s20
        end if                                                           5d3s10
        if(nbraa.ne.nbra)then
         itmpu=itmp1                                                     5d7s10
         if(idwsdeb.ne.0)then
          write(6,*)('b4 folding bra ')
          call prntm2(bc(itmpu),nbraa,nketa,nbraa)
         end if
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
         if(idwsdeb.ne.0)then
          write(6,*)('after folding bra ')
          call prntm2(bc(itmpu),nbraa,nketa,nbraa)
         end if
         if(idorel.ne.0)then                                            1d13s20
          itmpu=ltmp1                                                     5d7s10
         if(idwsdeb.ne.0)then
          write(6,*)('dorel part')
          write(6,*)('b4 folding bra ')
          call prntm2(bc(itmpu),nbraa,nketa,nbraa)
         end if
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
         if(idwsdeb.ne.0)then
          write(6,*)('dorel part')
          write(6,*)('after folding bra ')
          call prntm2(bc(itmpu),nbraa,nketa,nbraa)
         end if
         end if                                                         1d11s20
        end if                                                           5d3s10
        itmpu=itmp1                                                      5d7s10
         do ik=1,nketa                                                   5d3s10
          ikk=ik+ibc(jket4)                                              5d3s10
          isk=multh(isstor(ikk),iosym(iopr))                            12d15s19
          do ib=1,nbraa                                                  5d3s10
           ibb=ib+ibc(jbra4)                                             5d3s10
           iad1=ib-1+nbraa*(ik-1)                                        5d3s10
           if(isk.eq.isstor(ibb))then                                    5d3s10
            iad=ixmt(isstor(ibb),iopr)+ibstor(ibb)-1
     $          +ncomp*nbasisp(isstor(ibb))*(ibstor(ikk)-1)             1d11s20
            bc(iad)=bc(iad)+bc(itmpu+iad1)                              12d15s19
           else                                                          5d3s10
           end if                                                       12d16s19
          end do                                                         5d7s10
         end do                                                          5d3s10
        if(idorel.ne.0)then                                             8d20s15
         itmpu=ltmp1                                                    1d13s20
         do ik=1,nketa                                                   5d3s10
          ikk=ik+ibc(jket4)                                              5d3s10
          isk=multh(isstor(ikk),iosym(iopr))                            12d15s19
          do ib=1,nbraa                                                  5d3s10
           ibb=ib+ibc(jbra4)                                             5d3s10
           iad1=ib-1+nbraa*(ik-1)                                        5d3s10
           if(isk.eq.isstor(ibb))then                                    5d3s10
            iad=ixmt(isstor(ibb),iopr)+ibstor(ibb)+nbasisp(isstor(ibb)) 1d11s20
     $          -1+ncomp*nbasisp(isstor(ibb))*(ibstor(ikk)              1d11s20
     $          +nbasisp(isstor(ikk))-1)                                1d11s20
            bc(iad)=bc(iad)+bc(itmpu+iad1)                              12d15s19
           else                                                          5d3s10
           end if                                                       12d16s19
          end do                                                         5d7s10
         end do                                                          5d3s10
        end if                                                          8d20s15
        ibcoff=itmp1                                                     5d3s10
       end if                                                            5d3s10
      end do                                                            2d19s10
      call second(time3)
      call dws_sync                                                     2d22s10
      call dws_gsumf(bc(ixmt(1,1)),nbb)                                 12d15s19
c
c     now transform to mo basis
c
      do iopr=1,nop                                                     12d15s19
       do isb=1,nsymb                                                   12d15s19
        isk=multh(isb,iosym(iopr))                                      12d15s19
        if(isk.lt.1.or.isk.gt.nsymb)stop
        if(min(nbasisp(isb),nbasdwx(isb),nbasisp(isk),nbasdwx(isk))     12d16s19
     $       .gt.0)then                                                 12d16s19
         if(idwsdeb.ne.0)then                                           10d31s20
          write(6,*)('for operator '),opname(iopr)
          write(6,*)('for symmetry blocks '),isb,isk
          write(6,*)('in primitive basis: ')
          call prntm2(bc(ixmt(isb,iopr)),nbasisp(isb)*ncomp,             3d3s20
     $        nbasisp(isk)*ncomp,nbasisp(isb)*ncomp)                    3d3s20
          write(6,*)('using vectors '),iorb(isk),loc(bc(iorb(isk)))
          call prntm2(bc(iorb(isk)),nbasisp(isk)*ncomp,nbasdwx(isk),
     $         nbasisp(isk)*ncomp)
         end if                                                         10d31s20
         itmp=ibcoff                                                     12d15s19
         ibcoff=itmp+ncomp*nbasisp(isb)*nbasdwx(isk)                    1d11s20
         call enough('parap. 11',bc,ibc)
         call dgemm('n','n',ncomp*nbasisp(isb),nbasdwx(isk),
     $        ncomp*nbasisp(isk),                                       1d11s20
     $       1d0,bc(ixmt(isb,iopr)),nbasisp(isb)*ncomp,bc(iorb(isk)),   1d11s20
     $       nbasisp(isk)*ncomp,0d0,bc(itmp),nbasisp(isb)*ncomp,        1d11s20
     d' parap.  1')
         do i=0,nbasdwx(isk)-1                                           12d15s19
          do j=0,nbasisp(isb)*ncomp-1                                   1d11s20
           ji=itmp+j+ncomp*nbasisp(isb)*i                               1d11s20
           ij=ixmt(isb,iopr)+i+nbasdwx(isk)*j                            12d15s19
           bc(ij)=bc(ji)                                                 12d15s19
          end do                                                         12d15s19
         end do                                                          12d15s19
         call dgemm('n','n',nbasdwx(isk),nbasdwx(isb),
     $        nbasisp(isb)*ncomp,                                       1d11s20
     $        1d0,bc(ixmt(isb,iopr)),nbasdwx(isk),bc(iorb(isb)),         12d15s19
     $       nbasisp(isb)*ncomp,0d0,bc(itmp),nbasdwx(isk),              1d11s20
     d' parap.  2')
         do i=0,nbasdwx(isb)-1                                           12d15s19
          do j=0,nbasdwx(isk)-1                                          12d15s19
           ji=itmp+j+nbasdwx(isk)*i                                      12d15s19
           ij=ixmt(isb,iopr)+i+nbasdwx(isb)*j                            12d15s19
           bc(ij)=bc(ji)                                                 12d15s19
          end do                                                         12d15s19
         end do                                                          12d15s19
         if(idwsdeb.ne.0)then                                           10d31s20
         write(6,*)('in mo basis'),ixmt(isb,iopr)
         call prntm2(bc(ixmt(isb,iopr)),nbasdwx(isb),nbasdwx(isk),       12d16s19
     $       nbasdwx(isb))                                              12d16s19
         end if                                                         10d31s20
         ibcoff=itmp                                                     12d15s19
        end if                                                          12d16s19
       end do                                                           12d15s19
      end do                                                            12d15s19
      ibcoff=ibcoffo                                                    2d19s10
      return
      end                                                               2d19s10
