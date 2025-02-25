c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine parah042c(natom,ngaus,ibdat,ih0,ih0i,iovr,             3d19s20
     $                ibstor,isstor,idwsdeb,ascale,nbasisp,             3d4s20
     $                iorb,iapair,multh,ih0n,isopt,nsopt,srh,isym,bc,   11d9s22
     $                ibc,smsz)                                         11d16s23
c
c     compute one e hamitonian matrix and overlap matrix for 2 component3d4s20
c     spin-free spinors in the full 4 component space.                  3d4s20
c
      implicit real*8 (a-h,o-z)
      external second
      integer*8 nn8                                                     2d19s10
      include "common.hf"
      include "common.store"
      include "common.spher"
      logical log1,lprint                                               3d2s22
      integer*8 ibstor,isstor                                           3d4s20
      include "common.basis"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension ih0(8,8),iovr(8,8),ibstor(1),isstor(1),iorb(*),         3d4s20
     $     nbasisp(*),iapair(3,*),jh0(8,8),jovr(8,8),ih0i(8,8),         3d22s20
     $     multh(8,8),ih0n(8,4),isopt(4,4),cartb(3),cartk(3),isym(3,*)  3d11s22
      lprint=mynowprog.eq.0                                             3d2s22
      if(lprint)write(6,*)('hi, my name in parah042c '),ibcoff          3d2s22
      do i=1,4                                                          9d22s21
       do j=1,nsymb                                                     9d21s21
        ih0n(j,i)=0                                                     9d21s21
       end do                                                           9d21s21
      end do                                                            9d21s21
      ibcs=ibcoff                                                       10d13s21
      if(nsopt.eq.1.and.isopt(2,1).eq.0.and.isopt(3,1).eq.0)then        9d21s21
       if(lprint)                                                       3d2s22
     $    write(6,*)('there are no one-electron spin-orbital operators')3d2s22
       return                                                           9d21s21
      else                                                              9d21s21
       jj=0                                                             9d21s21
       do jsopt=1,nsopt                                                 9d21s21
         jj=jj+1                                                        9d21s21
         do isb=1,nsymb                                                 9d21s21
          jsb=multh(isb,isopt(1,jsopt))                                 9d21s21
          if(isb.le.jsb)then                                            9d22s21
           if(min(nbasdws(isb),nbasdws(jsb)).gt.0)then                   9d21s21
            ih0n(isb,jj)=ibcoff                                           9d21s21
            ibcoff=ih0n(isb,jj)+nbasdws(isb)*nbasdws(jsb)                 9d21s21
           end if                                                        9d21s21
          end if                                                        9d22s21
         end do                                                         9d21s21
       end do                                                           9d21s21
       call enough('parah042c.  1',bc,ibc)
      end if                                                            9d21s21
      ibcoffo=ibcoff                                                    2d19s10
      nwds=0                                                            3d4s20
      do isk=1,nsymb                                                    3d4s20
       do isb=1,nsymb                                                   3d4s20
        need=nbasisp(isb)*nbasisp(isk)*16                               3d4s20
        jh0(isb,isk)=ibcoff                                             3d4s20
        jovr(isb,isk)=jh0(isb,isk)+need*2                               3d4s20
        ibcoff=jovr(isb,isk)+need*2                                     3d4s20
        nwds=nwds+need*4                                                3d4s20
       end do                                                           3d4s20
      end do                                                            3d4s20
      call enough('parah042c.  2',bc,ibc)
      do i=ibcs,ibcoff-1                                                10d13s21
       bc(i)=0d0                                                        3d4s20
      end do                                                            3d4s20
      pi=acos(-1d0)                                                     2d21s20
      oneover4pi=0.25d0/pi                                              2d21s20
c     ascale=1/(4*c*c), thus c = 0.5/sqrt(ascale)
c
      clight=0.5d0/sqrt(ascale)                                         2d21s20
      if(idwsdeb.gt.100)then                                            5d25s18
       write(6,*)('in parah042c ')
       write(6,*)('speed of light: '),clight
       write(6,*)('ibdat '),ibdat
      do i=0,ngaus-1
       ip=i+1
       write(6,*)ip,ibc(ibdat+i),(bc(ibdat+i+ngaus*j),j=1,2),
     $      ibc(ibdat+i+3*ngaus),(bc(ibdat+i+ngaus*j),j=4,6),
     $      (ibc(ibdat+i+ngaus*j),j=7,8)
      end do
       call prntm2(bc(ibdat),ngaus,8,ngaus)
      end if                                                            5d25s18
      call second(time1)
      itry=13*25
c
c     build shell pair order
c
      ascale2=ascale*2d0                                                8d20s15
      nn=(ngaus*(ngaus+1))/2                                            2d19s10
      nn2=nn*2                                                          2d19s10
c
c     for large component,
c     overlap and nuclear attraction for each atom
c     for small component,
c     overlap and nuclear attraction for each atom
c     for ls coupling, dx,dy,dz
c     for sl coupling, dx,dy,dz
c
      ndo=2*(1+natom+3)                                                 2d21s20
      ipair=ibcoff                                                      2d19s10
      ibcoff=ipair+nn*ndo*3                                             2d19s10
      i12=ibcoff                                                        2d19s10
      ibcoff=i12+nn*4                                                   2d19s10
      call enough('parah042c.  3',bc,ibc)
      j12=i12                                                           2d19s10
      jbdat=ibdat-1                                                     2d19s10
      ngaus7=ngaus*7                                                    5d3s10
      do i1=1,ngaus                                                     2d19s10
       do i2=1,i1                                                       2d19s10
        ibc(j12)=-((ibc(jbdat+i1)+ibc(jbdat+i2)+2)**2)                  2d19s10
        ibc(j12+nn)=i1                                                  2d19s10
        ibc(j12+nn2)=i2                                                 2d19s10
        j12=j12+1                                                       2d19s10
       end do                                                           2d19s10
      end do                                                            2d19s10
      is12=i12+nn*3                                                     2d19s10
      nn8=nn                                                            2d19s10
      if(idwsdeb.gt.100)then                                            5d25s18
      end if                                                            5d25s18
      call idsortdws(ibc(i12),ibc(is12),nn8)                            1d18s23
      ii1=i12+nn                                                        2d19s10
      ii2=i12+nn2                                                       2d19s10
      jpair=ipair                                                       2d19s10
      iqy=0                                                             3d4s20
      do i=1,nn                                                         2d19s10
       j=ibc(is12+i-1)-1                                                2d19s10
    1  format(i5,i8,2i5)                                                2d19s10
       do k=1,ndo                                                       2d19s10
        iqy=iqy+1
        jpair0=jpair
        ibc(jpair)=ibc(ii1+j)                                           2d19s10
        jpair=jpair+1                                                   2d19s10
        ibc(jpair)=ibc(ii2+j)                                           2d19s10
        jpair=jpair+1                                                   2d19s10
        ibc(jpair)=k-1                                                  2d22s10
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
       if(idwsdeb.gt.10)write(6,2)i,ibc(jpair),ibc(jpair+1),
     $      ibc(jpair+2),time1
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
       nlb=2*ibc(jbra)+1                                                3d4s20
       nlk=2*ibc(jket)+1                                                3d4s20
c
c     the first two is for the two components and the second is for     3d4s20
c     alpha and beta spin.                                              3d4s20
c     we also have real and imaginary part ...                          3d4s20
c
       nlb2=nlb*2                                                       3d4s20
       nbsz=nlb2*2                                                      3d4s20
       nlk2=nlk*2                                                       3d4s20
       nksz=nlk2*2                                                      3d4s20
       nbsza=nbsz                                                       3d4s20
       nksza=nksz                                                       3d4s20
       ncb=1                                                            3d11s22
       if(iapair(1,ibc(jbra8)).gt.0)then                                3d11s22
        nbsza=nbsz*2                                                    3d11s22
        ncb=2                                                           3d11s22
       end if                                                           3d11s22
       nck=1                                                            3d11s22
       if(iapair(1,ibc(jket8)).gt.0)nck=2                               3d11s22
       nksza=nksz*nck                                                   3d11s22
       nlka=nlk*nck                                                     3d11s22
       nlba=nlb*ncb                                                     3d11s22
       nlka2=nlk2*nck                                                   3d11s22
       nlba2=nlb2*ncb                                                   3d11s22
       nnm=nbsza*nksza                                                  3d4s20
       itmpm=ibcoff                                                     3d4s20
       itmpi=itmpm+nnm                                                  3d4s20
       ibcoff=itmpi+nnm                                                 3d4s20
       call enough('parah042c.  4',bc,ibc)
       do ii=itmpm,ibcoff-1                                             3d4s20
        bc(ii)=0d0                                                      3d4s20
       end do                                                           3d4s20
       do iket=0,nck-1                                                  3d11s22
        if(iket.eq.0)then                                               3d11s22
         cartk(1)=bc(jket5)                                             3d11s22
         cartk(2)=bc(jket6)                                             3d11s22
         cartk(3)=bc(jket7)                                             3d11s22
        else                                                            3d11s22
         cartk(1)=bc(jket5)*dfloat(isym(1,iapair(2,ibc(jket8))))        3d11s22
         cartk(2)=bc(jket6)*dfloat(isym(2,iapair(2,ibc(jket8))))        3d11s22
         cartk(3)=bc(jket7)*dfloat(isym(3,iapair(2,ibc(jket8))))        3d11s22
        end if                                                          3d11s22
        do ibra=0,ncb-1                                                 3d11s22
         if(ibra.eq.0)then                                               3d11s22
          cartb(1)=bc(jbra5)                                             3d11s22
          cartb(2)=bc(jbra6)                                             3d11s22
          cartb(3)=bc(jbra7)                                             3d11s22
         else                                                            3d11s22
          cartb(1)=bc(jbra5)*dfloat(isym(1,iapair(2,ibc(jbra8))))        3d11s22
          cartb(2)=bc(jbra6)*dfloat(isym(2,iapair(2,ibc(jbra8))))        3d11s22
          cartb(3)=bc(jbra7)*dfloat(isym(3,iapair(2,ibc(jbra8))))        3d11s22
         end if                                                          3d11s22
         if(ibc(jpair+2).eq.0)then                                        2d20s20
c     large component overlap
          call onep(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),       3d11s22
     $      cartb(3),ibc(jket),bc(jket2),bc(jket3),                     3d11s22
     $      cartk,cartk(2),cartk(3),idum,ib1,0,0,0,0,0,0,0,0,0,bc,ibc)  11d9s22
          do i1=0,nlk-1                                                   3d4s20
           i1p=i1+iket*nlk                                              3d11s22
           do i2=0,nlb-1                                                3d11s22
            i2p=i2+ibra*nlb                                             3d11s22
            xint=bc(ib1+i1+nlk*i2)                                      3d11s22
            iaa=itmpm+i2p+nbsza*i1p                                     3d11s22
            ibb=itmpm+i2p+nlba2+nbsza*(i1p+nlka2)                       3d11s22
            bc(iaa)=xint                                                  3d4s20
            bc(ibb)=xint                                                  3d4s20
           end do                                                         3d4s20
          end do                                                          3d4s20
         else if(ibc(jpair+2).eq.1)then                                   2d20s20
c     small component overlap                                           2d20s20
          call ssdvrc(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),     3d11s22
     $      cartb(3),ibc(jket),bc(jket2),bc(jket3),                     3d11s22
     $      cartk,cartk(2),cartk(3),ib1,dum,0,clight,bc,ibc)            11d9s22
          iboff=nlb2*nlk2                                                 3d4s20
          do i1=0,nlk-1                                                   3d4s20
           i1p=i1+nlk                                                     3d4s20
           ii1p=i1+nlka+iket*nlk                                        3d11s22
           ii1ppp=ii1p+nlka2                                            3d11s22
           do i2=0,nlb-1                                                  3d4s20
            i2p=i2+nlb                                                    3d4s20
            ii2p=i2+nlba+ibra*nlb                                       3d11s22
            ii2ppp=ii2p+nlba2                                           3d11s22
c     alpha-alpha                                                       3d4s20
            iadaa=ib1+i1+nlk2*i2                                          3d4s20
            jadaa=itmpm+ii2p+nbsza*ii1p                                 3d11s22
c     real
            bc(jadaa)=bc(iadaa)                                           3d4s20
c     imag
            bc(jadaa+nnm)=bc(iadaa+iboff)                                 3d4s20
c     beta-alpha
            iadba=ib1+i1+nlk2*i2p
            jadba=itmpm+ii2ppp+nbsza*ii1p                               3d11s22
            bc(jadba)=bc(iadba)                                           3d4s20
            bc(jadba+nnm)=bc(iadba+iboff)                                 3d4s20
c     alpha-beta
            iadab=ib1+i1p+nlk2*i2                                         3d4s20
            jadab=itmpm+ii2p+nbsza*ii1ppp                               3d11s22
            bc(jadab)=bc(iadab)                                           3d4s20
            bc(jadab+nnm)=bc(iadab+iboff)                                 3d4s20
c     beta-beta
            iadbb=ib1+i1p+nlk2*i2p                                         3d4s20
            jadbb=itmpm+ii2ppp+nbsza*ii1ppp                             3d11s22
            bc(jadbb)=bc(iadbb)                                           3d4s20
            bc(jadbb+nnm)=bc(iadbb+iboff)                                 3d4s20
           end do                                                         3d4s20
          end do                                                          3d4s20
         else if(ibc(jpair+2).le.7)then                                   2d20s20
          idx=0                                                           2d20s20
          idy=0                                                           2d20s20
          idz=0                                                           2d20s20
          if(ibc(jpair+2).le.4)then
c     small-large coupling
           if(ibc(jpair+2).eq.2)idx=1                                      2d20s20
           if(ibc(jpair+2).eq.3)idy=1                                      2d20s20
           if(ibc(jpair+2).eq.4)idz=1                                      2d20s20
           call sldvrc(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),    3d11s22
     $        cartb(3),ibc(jket),bc(jket2),bc(jket3),cartk,             3d11s22
     $        cartk(2),cartk(3),ib1,idx,idy,idz,bc,ibc)                 11d9s22
           iboff=nlb2*nlk2                                                 3d4s20
           do i1=0,nlk-1                                                   3d4s20
            i1p=i1+nlk                                                     3d4s20
            ii1=i1+iket*nlk                                             3d11s22
            ii1p=ii1+nlka                                               3d11s22
            ii1pp=ii1p+nlka                                             3d11s22
            ii1ppp=ii1pp+nlka                                           3d11s22
            do i2=0,nlb-1                                                  3d4s20
             i2p=i2+nlb                                                    3d4s20
             ii2=i2+ibra*nlb                                            3d11s22
             ii2p=ii2+nlba                                              3d11s22
             ii2pp=ii2p+nlba                                            3d11s22
             ii2ppp=ii2pp+nlba                                          3d11s22
c     alpha-alpha                                                       3d4s20
             iadaa=ib1+i1+nlk2*i2                                          3d4s20
             jadaa=itmpm+ii2p+nbsza*ii1                                 3d11s22
c     real
             bc(jadaa)=bc(iadaa)                                           3d4s20
c     imag
             bc(jadaa+nnm)=bc(iadaa+iboff)                                 3d4s20
c     beta-alpha
             iadba=ib1+i1+nlk2*i2p
             jadba=itmpm+ii2ppp+nbsza*ii1                               3d11s22
             bc(jadba)=bc(iadba)                                           3d4s20
             bc(jadba+nnm)=bc(iadba+iboff)                                 3d4s20
c     alpha-beta
             iadab=ib1+i1p+nlk2*i2                                         3d4s20
             jadab=itmpm+ii2p+nbsza*ii1pp                               3d11s22
             bc(jadab)=bc(iadab)                                           3d4s20
             bc(jadab+nnm)=bc(iadab+iboff)                                 3d4s20
c     beta-beta
             iadbb=ib1+i1p+nlk2*i2p                                         3d4s20
             jadbb=itmpm+ii2ppp+nbsza*ii1pp                             3d11s22
             bc(jadbb)=bc(iadbb)                                           3d4s20
             bc(jadbb+nnm)=bc(iadbb+iboff)                                 3d4s20
            end do                                                         3d4s20
           end do                                                          3d4s20
          else if(ibc(jpair).ne.ibc(jpair+1))then                         3d5s20
c     large-small coupling
c     we will be storing transpose so this falls below the diagonal     2d21s20
           if(ibc(jpair+2).eq.5)idx=1                                      2d20s20
           if(ibc(jpair+2).eq.6)idy=1                                      2d20s20
           if(ibc(jpair+2).eq.7)idz=1                                      2d20s20
           call sldvrc(ibc(jket),bc(jket2),bc(jket3),cartk,cartk(2),    3d11s22
     $        cartk(3),ibc(jbra),bc(jbra2),bc(jbra3),cartb,             3d11s22
     $        cartb(2),cartb(3),ib1,idx,idy,idz,bc,ibc)                 11d9s22
           iboff=nlb2*nlk2                                                 3d4s20
           do i1=0,nlk-1                                                   3d4s20
            i1p=i1+nlk                                                     3d4s20
            ii1=i1+iket*nlk                                             3d11s22
            ii1p=ii1+nlka                                               3d11s22
            ii1pp=ii1p+nlka                                             3d11s22
            ii1ppp=ii1pp+nlka                                           3d11s22
            do i2=0,nlb-1                                                  3d4s20
             i2p=i2+nlb                                                    3d4s20
             ii2=i2+ibra*nlb                                            3d11s22
             ii2p=ii2+nlba                                               3d11s22
             ii2pp=ii2p+nlba                                             3d11s22
             ii2ppp=ii2pp+nlba                                           3d11s22
c     alpha-alpha                                                       3d4s20
             iadaa=ib1+i2+nlb2*i1                                          3d4s20
             jadaa=itmpm+ii2+nbsza*ii1p                                 3d11s22
c     real
             bc(jadaa)=bc(iadaa)                                           3d4s20
c     imag
             bc(jadaa+nnm)=-bc(iadaa+iboff)                                 3d4s20
c     beta-alpha
             iadba=ib1+i2p+nlb2*i1                                        3d4s20
             jadba=itmpm+ii2+nbsza*ii1ppp                               3d11s22
             bc(jadba)=bc(iadba)                                           3d4s20
             bc(jadba+nnm)=-bc(iadba+iboff)                               3d4s20
c     alpha-beta
             iadab=ib1+i2+nlb2*i1p                                        3d4s20
             jadab=itmpm+ii2pp+nbsza*ii1p                               3d11s22
             bc(jadab)=bc(iadab)                                           3d4s20
             bc(jadab+nnm)=bc(iadab+iboff)                                 3d4s20
c     beta-beta
             iadbb=ib1+i2p+nlb2*i1p                                         3d4s20
             jadbb=itmpm+ii2pp+nbsza*ii1ppp                             3d11s22
             bc(jadbb)=bc(iadbb)                                           3d4s20
             bc(jadbb+nnm)=bc(iadbb+iboff)                                 3d4s20
            end do                                                         3d4s20
           end do                                                          3d4s20
          end if                                                          3d4s20
         else if(ibc(jpair+2).le.(7+natom))then                             2d20s20
c     large nuclear attraction
          ia=ibc(jpair+2)-7                                               2d20s20
          zm=-atnum(1,ia)                                                 2d21s20
          call derid(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),      3d11s22
     $         cartb(3),ibc(jbra4),ibc(jket),bc(jket2),bc(jket3),       3d11s22
     $         cartk,cartk(2),cartk(3),ibc(jket4),                      3d11s22
     $         0,atnum(2,ia),atnum(3,ia),xcart(1,ia),xcart(2,ia),       8d20s15
     $         xcart(3,ia),0,                                           8d20s15
     $         0,atnum(2,ia),atnum(3,ia),xcart(1,ia),xcart(2,ia),       8d20s15
     $         xcart(3,ia),0,dum,1,idum,ib1,0,0,0,0,0,0,0,0,0,0,0,0,    2d21s20
     $         1,.false.,zm,bc,ibc)                                     11d14s22
          do i1=0,nlk-1                                                   3d4s20
           i1p=i1+iket*nlk                                              3d11s22
           do i2=0,nlb-1
            i2p=i2+ibra*nlb                                             3d11s22
            xint=bc(ib1+i1+nlk*i2)                                        3d4s20
            iaa=itmpm+i2p+nbsza*i1p                                     3d11s22
            ibb=itmpm+i2p+nlba2+nbsza*(i1p+nlka2)                         3d11s22
            bc(iaa)=xint                                                  3d4s20
            bc(ibb)=xint                                                  3d4s20
           end do                                                         3d4s20
          end do                                                          3d4s20
         else                                                             2d20s20
c     small nuclear attraction                                          2d20s20
          ia=ibc(jpair+2)-7-natom                                         2d20s20
          zm=-atnum(1,ia)                                                 2d21s20
          idata=ibcoff                                                    3d4s20
          ibcoff=idata+6                                                  3d4s20
          call enough('parah042c.  5',bc,ibc)
          bc(idata)=atnum(2,ia)                                           3d4s20
          bc(idata+1)=atnum(3,ia)                                         3d4s20
          bc(idata+2)=xcart(1,ia)                                         3d4s20
          bc(idata+3)=xcart(2,ia)                                         3d4s20
          bc(idata+4)=xcart(3,ia)                                         3d4s20
          bc(idata+5)=zm                                                  3d4s20
          call ssdvrc(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),     3d11s22
     $      cartb(3),ibc(jket),bc(jket2),bc(jket3),                     3d11s22
     $      cartk,cartk(2),cartk(3),ib1,bc(idata),2,clight,bc,ibc)      11d9s22
          iboff=nlb2*nlk2                                                 3d4s20
          do i1=0,nlk-1                                                   3d4s20
           i1p=i1+nlk                                                     3d4s20
           ii1p=i1+nlka+iket*nlk                                        3d11s22
           ii1ppp=ii1p+nlka2                                            3d11s22
           do i2=0,nlb-1
            i2p=i2+nlb                                                    3d4s20
            ii2p=i2+nlba+ibra*nlb                                       3d11s22
            ii2ppp=ii2p+nlba2                                           3d11s22
c     alpha-alpha                                                       3d4s20
            iadaa=ib1+i1+nlk2*i2                                          3d4s20
            jadaa=itmpm+ii2p+nbsza*ii1p                                 3d11s22
c     real
            bc(jadaa)=bc(iadaa)                                           3d4s20
c     imag
            bc(jadaa+nnm)=bc(iadaa+iboff)                                 3d4s20
c     beta-alpha
            iadba=ib1+i1+nlk2*i2p
            jadba=itmpm+ii2ppp+nbsza*ii1p                               3d11s22
            bc(jadba)=bc(iadba)                                           3d4s20
            bc(jadba+nnm)=bc(iadba+iboff)                                 3d4s20
c     alpha-beta
            iadab=ib1+i1p+nlk2*i2                                         3d4s20
            jadab=itmpm+ii2p+nbsza*ii1ppp                               3d11s22
            bc(jadab)=bc(iadab)                                           3d4s20
            bc(jadab+nnm)=bc(iadab+iboff)                                 3d4s20
c     beta-beta
            iadbb=ib1+i1p+nlk2*i2p                                         3d4s20
            jadbb=itmpm+ii2ppp+nbsza*ii1ppp                             3d11s22
            bc(jadbb)=bc(iadbb)                                           3d4s20
            bc(jadbb+nnm)=bc(iadbb+iboff)                                 3d4s20
           end do                                                         3d4s20
          end do                                                          3d4s20
         end if                                                           2d20s20
        end do                                                          3d11s22
       end do                                                           3d11s22
c
c     symmetrize                                                        3d11s22
c
       if(nck.eq.2)then                                                 3d11s22
        do ix=0,7                                                       3d11s22
         ioff=ix*nlk2                                                   3d11s22
         do ik=0,nlk-1                                                   3d11s22
          iad1=itmpm+nbsza*(ik+ioff)                                     3d11s22
          iad2=itmpm+nbsza*(ik+nlk+ioff)                                 3d11s22
          do ib=0,nbsza-1                                                3d11s22
           sum=srh*(bc(iad1+ib)+bc(iad2+ib))                            3d11s22
           dif=srh*(-bc(iad1+ib)+bc(iad2+ib))                           3d11s22
           bc(iad1+ib)=sum                                               3d11s22
           bc(iad2+ib)=dif                                               3d11s22
          end do                                                         3d11s22
         end do                                                          3d11s22
        end do                                                          3d11s22
       end if                                                           3d11s22
       if(ncb.eq.2)then                                                 3d11s22
        do iri=0,1                                                      3d11s22
         do ix=0,3                                                       3d11s22
          ioff=ix*nlb2+iri*nnm                                          3d11s22
          do ib=0,nlb-1                                                  3d11s22
           iad1=itmpm+ib+ioff                                            3d11s22
           iad2=itmpm+ib+nlb+ioff                                        3d11s22
           do ik=0,nksza-1                                                3d11s22
            kk=ik*nbsza                                                   3d11s22
            sum=srh*(bc(iad1+kk)+bc(iad2+kk))                            3d11s22
            dif=srh*(-bc(iad1+kk)+bc(iad2+kk))                           3d11s22
            bc(iad1+kk)=sum                                               3d11s22
            bc(iad2+kk)=dif                                               3d11s22
           end do                                                         3d11s22
          end do                                                          3d11s22
         end do                                                          3d11s22
        end do                                                          3d11s22
       end if                                                           3d11s22
c
c     store ...
c
       if(ibc(jpair+2).le.1)then                                        3d4s20
        io=1                                                            3d4s20
       else                                                             3d4s20
        io=0                                                            3d4s20
       end if                                                           3d4s20
       ih=1-io                                                          3d4s20
       do ik=0,nlka-1                                                   3d11s22
        ikk=ik+ibc(jket4)+1                                             3d4s20
        isk=isstor(ikk)                                                 3d4s20
        if(ibc(jpair).eq.ibc(jpair+1))then                              3d11s22
         ibtop=ik                                                       3d11s22
        else                                                            3d11s22
         ibtop=nlba-1                                                   3d11s22
        end if                                                          3d11s22
        do ib=0,ibtop                                                   3d11s22
         ibb=ib+ibc(jbra4)+1                                            3d4s20
         isb=isstor(ibb)                                                3d4s20
         iadlala=(io*jovr(isb,isk)+ih*jh0(isb,isk))+ibstor(ibb)-1+      3d4s20
     $        4*nbasisp(isb)*(ibstor(ikk)-1)                            3d4s20
         iadsala=iadlala+nbasisp(isb)                                   3d4s20
         iadlbla=iadsala+nbasisp(isb)                                   3d4s20
         iadsbla=iadlbla+nbasisp(isb)                                   3d4s20
         iadlasa=iadlala+4*nbasisp(isb)*nbasisp(isk)                    3d4s20
         iadsasa=iadlasa+nbasisp(isb)                                   3d4s20
         iadlbsa=iadsasa+nbasisp(isb)                                   3d4s20
         iadsbsa=iadlbsa+nbasisp(isb)                                   3d4s20
         iadlalb=iadlasa+4*nbasisp(isb)*nbasisp(isk)                    3d4s20
         iadsalb=iadlalb+nbasisp(isb)                                   3d4s20
         iadlblb=iadsalb+nbasisp(isb)                                   3d4s20
         iadsblb=iadlblb+nbasisp(isb)                                   3d4s20
         iadlasb=iadlalb+4*nbasisp(isb)*nbasisp(isk)                    3d4s20
         iadsasb=iadlasb+nbasisp(isb)                                   3d4s20
         iadlbsb=iadsasb+nbasisp(isb)                                   3d4s20
         iadsbsb=iadlbsb+nbasisp(isb)                                   3d4s20
         jadlala=itmpm+ib+nbsza*ik                                      3d4s20
         jadsala=jadlala+nlba                                           3d11s22
         jadlbla=jadsala+nlba                                           3d11s22
         jadsbla=jadlbla+nlba                                           3d11s22
         jadlasa=jadlala+nbsza*nlka                                     3d11s22
         jadsasa=jadlasa+nlba                                           3d11s22
         jadlbsa=jadsasa+nlba                                           3d11s22
         jadsbsa=jadlbsa+nlba                                           3d11s22
         jadlalb=jadlasa+nbsza*nlka                                     3d11s22
         jadsalb=jadlalb+nlba                                           3d11s22
         jadlblb=jadsalb+nlba                                           3d11s22
         jadsblb=jadlblb+nlba                                           3d11s22
         jadlasb=jadlalb+nbsza*nlka                                     3d11s22
         jadsasb=jadlasb+nlba                                           3d11s22
         jadlbsb=jadsasb+nlba                                           3d11s22
         jadsbsb=jadlbsb+nlba                                           3d11s22
         nns=nbasisp(isb)*nbasisp(isk)*16                               3d4s20
         bc(iadlala)=bc(iadlala)+bc(jadlala)                            3d4s20
         bc(iadlala+nns)=bc(iadlala+nns)+bc(jadlala+nnm)                3d4s20
         bc(iadsala)=bc(iadsala)+bc(jadsala)                                        3d4s20
         bc(iadsala+nns)=bc(iadsala+nns)+bc(jadsala+nnm)                                3d4s20
         bc(iadlbla)=bc(iadlbla)+bc(jadlbla)                                        3d4s20
         bc(iadlbla+nns)=bc(iadlbla+nns)+bc(jadlbla+nnm)                                3d4s20
         bc(iadsbla)=bc(iadsbla)+bc(jadsbla)                                        3d4s20
         bc(iadsbla+nns)=bc(iadsbla+nns)+bc(jadsbla+nnm)                                3d4s20
         bc(iadlasa)=bc(iadlasa)+bc(jadlasa)                                        3d4s20
         bc(iadlasa+nns)=bc(iadlasa+nns)+bc(jadlasa+nnm)                                3d4s20
         bc(iadsasa)=bc(iadsasa)+bc(jadsasa)                                        3d4s20
         bc(iadsasa+nns)=bc(iadsasa+nns)+bc(jadsasa+nnm)                                3d4s20
         bc(iadlbsa)=bc(iadlbsa)+bc(jadlbsa)                                        3d4s20
         bc(iadlbsa+nns)=bc(iadlbsa+nns)+bc(jadlbsa+nnm)                                3d4s20
         bc(iadsbsa)=bc(iadsbsa)+bc(jadsbsa)                                        3d4s20
         bc(iadsbsa+nns)=bc(iadsbsa+nns)+bc(jadsbsa+nnm)                                3d4s20
         bc(iadlalb)=bc(iadlalb)+bc(jadlalb)                                        3d4s20
         bc(nns+iadlalb)=bc(nns+iadlalb)+bc(nnm+jadlalb)                                3d4s20
         bc(iadsalb)=bc(iadsalb)+bc(jadsalb)                                        3d4s20
         bc(nns+iadsalb)=bc(nns+iadsalb)+bc(nnm+jadsalb)                                3d4s20
         bc(iadlblb)=bc(iadlblb)+bc(jadlblb)                                        3d4s20
         bc(nns+iadlblb)=bc(nns+iadlblb)+bc(nnm+jadlblb)                                3d4s20
         bc(iadsblb)=bc(iadsblb)+bc(jadsblb)                                        3d4s20
         bc(nns+iadsblb)=bc(nns+iadsblb)+bc(nnm+jadsblb)                                3d4s20
         bc(iadlasb)=bc(iadlasb)+bc(jadlasb)                                        3d4s20
         bc(nns+iadlasb)=bc(nns+iadlasb)+bc(nnm+jadlasb)                                3d4s20
         bc(iadsasb)=bc(iadsasb)+bc(jadsasb)                                        3d4s20
         bc(nns+iadsasb)=bc(nns+iadsasb)+bc(nnm+jadsasb)                                3d4s20
         bc(iadlbsb)=bc(iadlbsb)+bc(jadlbsb)                                        3d4s20
         bc(nns+iadlbsb)=bc(nns+iadlbsb)+bc(nnm+jadlbsb)                                3d4s20
         bc(iadsbsb)=bc(iadsbsb)+bc(jadsbsb)                                        3d4s20
         bc(nns+iadsbsb)=bc(nns+iadsbsb)+bc(nnm+jadsbsb)                                3d4s20
        end do                                                          3d4s20
       end do                                                           3d4s20
       ibcoff=itmpm                                                     3d4s20
      end do                                                            2d19s10
      call second(time3)
      call dws_sync                                                     2d22s10
      call dws_gsumf(bc(jh0(1,1)),nwds)                                 3d4s20
      twomc2=2d0*clight*clight                                          3d5s20
      do isb=1,nsymb                                                    3d5s20
       do jsb=1,nsymb                                                   3d5s20
        nb2=nbasisp(jsb)*2                                              3d5s20
        nb4=nb2*2                                                       3d5s20
        nk2=nbasisp(isb)*2                                              3d5s20
        nk4=nk2*2                                                       3d5s20
        nri=nb4*nk4                                                     3d5s20
        do ispink=0,1                                                   3d5s20
         iadko=jovr(jsb,isb)+nb4*(nbasisp(isb)+nk2*ispink)              3d5s20
         iadkh=jh0(jsb,isb)+nb4*(nbasisp(isb)+nk2*ispink)               3d5s20
         do ispinb=0,1                                                   3d5s20
          jadko=iadko+nbasisp(jsb)+nb2*ispinb                           3d5s20
          jadkh=iadkh+nbasisp(jsb)+nb2*ispinb                           3d5s20
          do i=0,nbasisp(isb)-1                                         3d5s20
           kadko=jadko+nb4*i                                            3d5s20
           kadkh=jadkh+nb4*i                                            3d5s20
           do j=0,nbasisp(jsb)-1                                        3d5s20
            bc(kadkh+j)=bc(kadkh+j)-twomc2*bc(kadko+j)                  3d5s20
            bc(kadkh+j+nri)=bc(kadkh+j+nri)-twomc2*bc(kadko+j+nri)      3d5s20
           end do                                                       3d5s20
          end do                                                        3d5s20
         end do                                                         3d5s20
        end do                                                          3d5s20
       end do                                                           3d5s20
      end do                                                            3d5s20
      do isb=1,nsymb                                                    3d5s20
       nb2=nbasisp(isb)*2                                               3d5s20
       nb4=nb2*2                                                        3d5s20
       nb16=nb4*nb4                                                     3d5s20
       do i=0,nb4-1                                                     3d5s20
        do j=0,i-1                                                      3d5s20
         ji=jovr(isb,isb)+j+nb4*i                                       3d5s20
         ij=jovr(isb,isb)+i+nb4*j                                       3d5s20
         sum=bc(ji)+bc(ij)                                              3d5s20
         bc(ji)=sum                                                     3d5s20
         bc(ij)=sum                                                     3d5s20
         diff=bc(ji+nb16)-bc(ij+nb16)                                    3d5s20
         bc(ji+nb16)=diff                                               3d5s20
         bc(ij+nb16)=-diff                                              3d5s20
         ji=jh0(isb,isb)+j+nb4*i                                        3d5s20
         ij=jh0(isb,isb)+i+nb4*j                                        3d5s20
         sum=bc(ji)+bc(ij)                                              3d5s20
         bc(ji)=sum                                                     3d5s20
         bc(ij)=sum                                                     3d5s20
         diff=bc(ji+nb16)-bc(ij+nb16)                                    3d5s20
         bc(ji+nb16)=diff                                               3d5s20
         bc(ij+nb16)=-diff                                              3d5s20
        end do                                                          3d5s20
       end do                                                           3d5s20
       do jsb=1,isb-1                                                   3d5s20
        nk2=nbasisp(jsb)*2                                              3d5s20
        nk4=nk2*2                                                       3d5s20
        nri=nb4*nk4                                                     3d5s20
        do i=0,nb4-1                                                    3d5s20
         do j=0,nk4-1                                                   3d5s20
          ji=jovr(jsb,isb)+j+nk4*i                                      3d5s20
          ij=jovr(isb,jsb)+i+nb4*j                                      3d5s20
          sum=bc(ji)+bc(ij)                                             3d5s20
          bc(ji)=sum                                                    3d5s20
          bc(ij)=sum                                                    3d5s20
          diff=bc(ji+nri)-bc(ij+nri)                                    3d5s20
          bc(ji+nri)=diff                                               3d5s20
          bc(ij+nri)=-diff                                              3d5s20
          ji=jh0(jsb,isb)+j+nk4*i                                       3d5s20
          ij=jh0(isb,jsb)+i+nb4*j                                       3d5s20
          sum=bc(ji)+bc(ij)                                             3d5s20
          bc(ji)=sum                                                    3d5s20
          bc(ij)=sum                                                    3d5s20
          diff=bc(ji+nri)-bc(ij+nri)                                    3d5s20
          bc(ji+nri)=diff                                               3d5s20
          bc(ij+nri)=-diff                                              3d5s20
         end do                                                         3d5s20
        end do                                                          3d5s20
       end do                                                           3d5s20
      end do                                                            3d5s20
      do isb=1,nsymb                                                    3d4s20
       if(nbasdws(isb).gt.0)then                                        3d4s20
        do itype=1,2
         if(itype.eq.1)then
          ibase=jovr(isb,isb)
         else
          ibase=jh0(isb,isb)
         end if                                                         3d4s20
         nb2=nbasisp(isb)*2                                              3d4s20
         nb4=nb2*2                                                      3d5s20
         nb16=nb4*nb4
         sz=0d0                                                         3d5s20
         do i=0,nb2-1                                                   3d5s20
          iad=ibase+nb16+nb4*i                                          3d5s20
          do j=0,nb2-1                                                  3d5s20
           sz=sz+bc(iad+j)**2                                           3d5s20
          end do                                                        3d5s20
         end do                                                         3d5s20
         sz=sqrt(sz/dfloat(nb2*nb2))                                    3d5s20
         jovrbb=ibase+nb2*(nb4+1)                                       3d4s20
         sz=0d0                                                         3d5s20
         do i=0,nb2-1                                                   3d5s20
          iad=jovrbb+nb16+nb4*i                                         3d5s20
          do j=0,nb2-1                                                  3d5s20
           sz=sz+bc(iad+j)**2                                           3d5s20
          end do
         end do                                                         3d5s20
         sz=sqrt(sz/dfloat(nb2*nb2))
         itmp=ibcoff
         itmpi=itmp+nb2*nb2
         ibcoff=itmpi+nb2*nb2
         call enough('parah042c.  6',bc,ibc)
         szaa=0d0                                                       3d5s20
         szbb=0d0                                                       3d5s20
         do i=0,nb2-1
          do j=0,nb2-1
           iad1=ibase+j+nb4*i
           iad2=jovrbb+j+nb4*i
           iad3=itmp+j+nb2*i
           bc(iad3)=0.5d0*(bc(iad1)+bc(iad2))
           bc(iad1)=bc(iad1)-bc(iad3)                                   3d5s20
           bc(iad2)=bc(iad2)-bc(iad3)                                   3d5s20
           szaa=szaa+bc(iad1)**2                                        3d5s20
           szbb=szbb+bc(iad2)**2                                        3d5s20
          end do
         end do
         szaa=sqrt(szaa/dfloat(nb2*nb2))                                3d5s20
         szbb=sqrt(szbb/dfloat(nb2*nb2))                                3d5s20
         ibcoff=itmp
         jovrab=ibase+nb4*nb2
         jovrba=ibase+nb2
        end do                                                          3d4s20
       end if                                                           3d4s20
      end do                                                            3d4s20
      isave=ibcoffo                                                     3d5s20
      iwtf=0                                                            3d10s22
      do isb=1,nsymb                                                    3d5s20
       nb2=nbasisp(isb)*2                                               3d5s20
       nb4=nb2*2                                                        3d5s20
       do jsb=1,nsymb                                                   3d5s20
        nk2=nbasisp(jsb)*2                                              3d5s20
        jisb=multh(jsb,isb)                                             3d22s20
        nk4=nk2*2                                                       3d5s20
        nbk=nb4*nk4                                                     3d5s20
        do itype=1,2                                                    3d5s20
         if(itype.eq.1)then
          ibase=jovr(jsb,isb)                                           3d5s20
          iovr(jsb,isb)=0                                               3d5s20
         else                                                           3d5s20
          ibase=jh0(jsb,isb)                                            3d5s20
          ih0(jsb,jisb)=0                                               3d22s20
          ih0i(jsb,jisb)=0                                              3d22s20
         end if                                                         3d5s20
         do ir=0,1
          sz=0d0                                                        3d5s20
          do i=0,nb4-1                                                  3d5s20
           iad=ibase+nk4*i                                              3d5s20
           do j=0,nk4-1                                                 3d5s20
            sz=sz+bc(iad+j)**2                                          3d5s20
           end do                                                       3d5s20
          end do                                                        3d5s20
          if(min(nbk,nk4).gt.0)then                                     1d13s23
           sz=sqrt(sz/dfloat(nb4*nk4))                                   3d5s20
          end if                                                        1d13s23
          if(sz.gt.smsz)then                                            6d12s23
           itmp=ibcoff                                                  3d5s20
           ibcoff=itmp+nk4*nb4                                          3d5s20
           call enough('parah042c.  7',bc,ibc)
           call dgemm('n','n',nk4,nbasdws(isb),nb2,1d0,bc(ibase),nk4,   3d5s20
     $          bc(iorb(isb)),nb2,0d0,bc(itmp),nk4,                     3d5s20
     d' parah042c.  1')
           itmpp=itmp+nk4*nbasdws(isb)                                  3d5s20
           ibasep=ibase+nk4*nb2                                         3d5s20
           call dgemm('n','n',nk4,nbasdws(isb),nb2,1d0,bc(ibasep),nk4,  3d5s20
     $          bc(iorb(isb)),nb2,0d0,bc(itmpp),nk4,                    3d5s20
     d' parah042c.  2')
           nbb2=nbasdws(isb)*2                                          3d5s20
           do i=0,nk4-1                                                 3d5s20
            do j=0,nbb2-1                                               3d5s20
             ji=ibase+j+nbb2*i                                          3d5s20
             ij=itmp+i+nk4*j                                            3d5s20
             bc(ji)=bc(ij)                                              3d5s20
            end do                                                      3d5s20
           end do                                                       3d5s20
           call dgemm('n','n',nbb2,nbasdws(jsb),nk2,1d0,bc(ibase),nbb2, 3d5s20
     $          bc(iorb(jsb)),nk2,0d0,bc(itmp),nbb2,                    3d5s20
     d' parah042c.  3')
           itmpp=itmp+nbb2*nbasdws(jsb)                                 3d5s20
           ibasep=ibase+nbb2*nk2                                        3d5s20
           call dgemm('n','n',nbb2,nbasdws(jsb),nk2,1d0,bc(ibasep),nbb2,3d5s20
     $          bc(iorb(jsb)),nk2,0d0,bc(itmpp),nbb2,                   3d5s20
     d' parah042c.  4')
           nkk2=nbasdws(jsb)*2                                          3d5s20
           do i=0,nbb2-1                                                 3d5s20
            do j=0,nkk2-1                                               3d5s20
             ji=ibase+j+nkk2*i                                          3d5s20
             ij=itmp+i+nbb2*j                                            3d5s20
             bc(ji)=bc(ij)                                              3d5s20
            end do                                                      3d5s20
           end do                                                       3d5s20
           ibcoff=itmp                                                  3d5s20
           sza=0d0                                                      9d2s21
           szb=0d0                                                      9d2s21
           rmsaa=0d0                                                    9d2s21
           rmsab=0d0                                                    9d2s21
           rmsba=0d0                                                    9d2s21
           rmsbb=0d0                                                    9d2s21
           do i=0,nbasdws(isb)-1                                        9d2s21
            ip=i+nbasdws(isb)                                           9d2s21
            do j=0,nbasdws(jsb)-1                                       9d2s21
             jp=j+nbasdws(jsb)                                          9d2s21
             ji=ibase+j+nkk2*i                                          9d2s21
             jip=ibase+jp+nkk2*ip                                       9d2s21
             sza=sza+bc(ji)**2                                          9d2s21
             rmsaa=rmsaa+(bc(ji)-bc(jip))**2                            9d2s21
             rmsbb=rmsbb+(bc(ji)+bc(jip))**2                            9d2s21
             ji=ibase+j+nkk2*ip                                         9d2s21
             jip=ibase+jp+nkk2*i                                        9d2s21
             szb=szb+bc(ji)**2                                          9d2s21
             rmsab=rmsab+(bc(ji)-bc(jip))**2                            9d2s21
             rmsba=rmsba+(bc(ji)+bc(jip))**2                            9d2s21
            end do                                                      9d2s21
           end do                                                       9d2s21
           sza=sqrt(sza/dfloat(nbasdws(isb)*nbasdws(jsb)))              9d2s21
           szb=sqrt(szb/dfloat(nbasdws(isb)*nbasdws(jsb)))              9d2s21
           rmsaa0=sqrt(rmsaa/dfloat(nbasdws(isb)*nbasdws(jsb)))          9d2s21
           rmsbb0=sqrt(rmsbb/dfloat(nbasdws(isb)*nbasdws(jsb)))          9d2s21
           rmsaa=rmsaa0/sza                                              6d7s23
           rmsbb=rmsbb0/sza                                              6d7s23
           rmsab0=sqrt(rmsab/dfloat(nbasdws(isb)*nbasdws(jsb)))          9d2s21
           rmsba0=sqrt(rmsba/dfloat(nbasdws(isb)*nbasdws(jsb)))          9d2s21
           rmsab=rmsab0/szb                                              6d7s23
           rmsba=rmsba0/szb                                              6d7s23
           jh0n=0                                                       9d21s21
           iiq=0                                                        9d21s21
           if(itype.eq.1)then                                           3d19s24
            iqtop=1                                                     3d19s24
           else                                                         3d19s24
            iqtop=nsopt                                                 3d19s24
           end if                                                       3d19s24
           do iq=1,iqtop                                                3d19s24
            iiq=iiq+1                                                   9d21s21
            if(isopt(2,iq).ne.0.or.isopt(3,iq).ne.0)then                9d21s21
             if(multh(isb,jsb).eq.isopt(1,iq).and.                      9d21s21
     $            isopt(2,iq).eq.ir)then                                9d21s21
              if(szb.gt.smsz.and.isopt(3,iq).eq.1)then                  6d12s23
               if(rmsab.lt.smsz.and.isopt(4,iq).eq.0)then               6d12s23
                jh0n=iiq                                                  9d21s21
                jjh0n=iq                                                9d21s21
               else if(rmsba.lt.smsz.and.isopt(4,iq).eq.1)then          6d12s23
                jh0n=iiq                                                  9d21s21
                jjh0n=iq                                                9d21s21
               end if                                                    9d21s21
              else if(sza.gt.smsz.and.isopt(3,iq).eq.0)then             6d12s23
               if(rmsaa.lt.smsz.and.isopt(4,iq).eq.0)then               6d12s23
                jh0n=iiq                                                  9d21s21
                jjh0n=iq                                                9d21s21
               else if(rmsbb.lt.smsz.and.isopt(4,iq).eq.1)then          6d12s23
                jh0n=iiq                                                  9d21s21
                jjh0n=iq                                                9d21s21
               end if                                                    9d21s21
              end if                                                     9d21s21
             end if                                                      9d21s21
            end if                                                      9d21s21
            if(jh0n.ne.0.or.itype.eq.1)then                             3d19s24
             if(szb.gt.smsz)then                                          6d12s23
              if(rmsab.lt.smsz)then                                       6d12s23
              else if(rmsba.lt.smsz)then                                  6d12s23
              else                                                        9d2s21
               write(6,*)('ab ? ba'),rmsab,rmsba,szb
               write(6,*)('threshold smsz = '),smsz                       11d16s23
               write(6,*)('in ms basis: ')
               call prntm2(bc(ibase),nkk2,nbb2,nkk2)
              end if                                                      9d2s21
             end if                                                       9d2s21
             if(sza.gt.smsz)then                                          6d12s23
              if(rmsaa.lt.smsz)then                                       6d12s23
              else if(rmsbb.lt.smsz)then                                  6d12s23
              else                                                        9d2s21
               write(6,*)('aa ? bb'),rmsaa,rmsbb,sza
               write(6,*)('threshold smsz = '),smsz                       11d16s23
               write(6,*)('in ms basis: ')
               call prntm2(bc(ibase),nkk2,nbb2,nkk2)
              end if                                                      9d2s21
             end if                                                       9d2s21
             if(jh0n.gt.0.or.itype.eq.1)then                              6d27s22
             else                                                         9d21s21
              write(6,*)('wtf? ')                                         9d21s21
              iwtf=iwtf+1
              go to 1066
             end if                                                       9d21s21
             if(ir.eq.0)then                                              3d5s20
              isaveu=isave                                                3d5s20
             else                                                         3d5s20
              isaveu=-isave                                               3d5s20
             end if                                                       3d5s20
             if(itype.eq.1)then                                           3d5s20
              iovr(jsb,isb)=isaveu                                        3d5s20
             else                                                         3d5s20
              if(ih0n(jsb,jh0n).gt.0)then                                 9d22s21
               do i=0,nbasdws(isb)-1                                       9d21s21
                ip=i+isopt(3,jjh0n)*nbasdws(isb)                           9d22s21
                iad2=ibase+nkk2*ip                                         9d21s21
                iad1=ih0n(jsb,jh0n)+nbasdws(jsb)*i                         9d21s21
                do j=0,nbasdws(jsb)-1                                      9d21s21
                 bc(iad1+j)=bc(iad2+j)                                     9d21s21
                end do                                                     9d21s21
               end do                                                      9d21s21
              end if                                                      9d22s21
              if(ir.eq.0)then                                             3d19s20
               ih0(jsb,jisb)=isave                                        3d22s20
              else                                                        3d19s20
               ih0i(jsb,jisb)=isave                                       3d22s20
              end if                                                      3d19s20
             end if                                                       3d5s20
             do i=0,nkk2*nbb2-1                                           3d5s20
              bc(isave+i)=bc(ibase+i)                                     3d5s20
             end do                                                       3d5s20
             isave=isave+nkk2*nbb2                                        3d5s20
            end if                                                        3d5s20
 1066       continue                                                      3d10s22
           end do                                                        3d19s24
          end if                                                        3d19s24
          ibase=ibase+nbk
c     ir
         end do
c     type
        end do                                                          3d5s20
c     jsb
       end do                                                           3d5s20
c     isb
      end do                                                            3d5s20
      if(iwtf.ne.0)then
       write(6,*)('final wtf = '),iwtf
       call dws_synca
       call dws_finalize
       stop
      end if
      ibcoff=isave                                                      3d5s20
      return
      end                                                               2d19s10
