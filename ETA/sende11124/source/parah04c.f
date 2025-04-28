c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine parah04c(natom,ngaus,ibdat,nbasis,h0,ovr,              2d24s20
     $                iblstor,ibsstor,nbb,idwsdeb,ascale,nbaslarge,     2d24s20
     $     nbassmall,ilsz,isnorm,nbtot,bc,ibc)                          11d9s22
      implicit real*8 (a-h,o-z)
      external second
      integer*8 nn8                                                     2d19s10
      include "common.hf"
      include "common.store"
      include "common.spher"
      logical log1                                                      6d6s18
      integer*8 iblstor,ibsstor,lbra8,lket8                             2d24s20
      include "common.basis"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension h0(*),ovr(1),iblstor(1),ibsstor(1)                      2d24s20
      ibcoffo=ibcoff                                                    2d19s10
      pi=acos(-1d0)                                                     2d21s20
      oneover4pi=0.25d0/pi                                              2d21s20
c     ascale=1/(4*c*c), thus c = 0.5/sqrt(ascale)
c
      clight=0.5d0/sqrt(ascale)                                         2d21s20
      if(idwsdeb.gt.100)then                                            5d25s18
       write(6,*)('in parah04c '),nbaslarge,nbassmall,nbb,ibcoff
       write(6,*)('speed of light: '),clight
       write(6,*)('ibdat '),ibdat
      do i=0,ngaus-1
       ip=i+1
       write(6,*)ip,ibc(ibdat+i),(bc(ibdat+i+ngaus*j),j=1,2),
     $      bc(isnorm+i),
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
      srh=sqrt(0.5d0)                                                   5d3s10
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
      call enough('parah04c.  1',bc,ibc)
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
c     overlap is real and block diagonal in alpha,beta,large,small.
c     h0 is complex hermition coupling everything.
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
       isnormb=isnorm+ibc(jpair)-1                                      2d20s20
       jbra2=jbra+ngaus                                                 2d19s10
       jbra3=jbra2+ngaus                                                2d19s10
       jbra4=jbra3+ngaus                                                2d19s10
       jbra5=jbra4+ngaus                                                2d19s10
       jbra6=jbra5+ngaus                                                2d19s10
       jbra7=jbra6+ngaus                                                2d19s10
       jbra8=jbra7+ngaus                                                5d4s10
       jket=jbdat+ibc(jpair+1)                                          2d19s10
       isnormk=isnorm+ibc(jpair+1)-1                                    2d20s20
       jket2=jket+ngaus                                                 2d19s10
       jket3=jket2+ngaus                                                2d19s10
       jket4=jket3+ngaus                                                2d19s10
       jket5=jket4+ngaus                                                2d19s10
       jket6=jket5+ngaus                                                2d19s10
       jket7=jket6+ngaus                                                2d19s10
       jket8=jket7+ngaus                                                5d4s10
       nbral=ibc(ilsz+ibc(jbra))                                        2d20s20
       nbras=ibc(ilsz+ibc(jbra)+1)                                      2d20s20
       nketl=ibc(ilsz+ibc(jket))                                        2d20s20
       nkets=ibc(ilsz+ibc(jket)+1)                                      2d20s20
       call enough('parah04c.  2',bc,ibc)
       lbra8=ibc(jbra)+1                                                2d20s20
       lket8=ibc(jket)+1                                                2d20s20
c                                                                       2d21s20
c     note: since cart to spher transformation is skipped, integrals    2d21s20
c     come out with indices swapped compared to nr or 2 component       2d21s20
c     calculation.                                                      2d21s20
c
       if(ibc(jpair+2).eq.0)then                                        2d20s20
c     large component overlap
        call onep(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),    2d19s10
     $      bc(jbra7),ibc(jket),bc(jket2),bc(jket3),                    2d20s20
     $      bc(jket5),bc(jket6),bc(jket7),idum,ib1,0,0,0,0,0,0,0,0,0,   11d9s22
     $      bc,ibc)                                                     11d9s22
        nlb=2*ibc(jbra)+1                                               2d21s20
        nlk=2*ibc(jket)+1
        xnorm=oneover4pi*sqrt(dfloat(nlb*nlk))
        if(ibc(jpair).eq.ibc(jpair+1))then                              2d21s20
         do i1=0,nketl-1                                                 2d21s20
          i1p=iblstor(ibc(jpair+1))+i1
          i1pp=i1p+nbaslarge                                             2d21s20
          jb1=ib1+nbral*i1                                               2d21s20
          do i2=0,nbral-1                                                2d21s20
           i2p=iblstor(ibc(jpair))+i2                                    2d21s20
           i2pp=i2p+nbaslarge                                            2d21s20
           lla=i2p+nbtot*i1p+1                                           2d21s20
           llb=lla+nbaslarge*(nbtot+1)                                   2d21s20
           ovr(lla)=bc(jb1+i2)*xnorm                                      2d21s20
           ovr(llb)=ovr(lla)                                               2d21s20
          end do                                                         2d21s20
         end do                                                          2d21s20
        else                                                            2d21s20
         do i1=0,nketl-1                                                 2d21s20
          i1p=iblstor(ibc(jpair+1))+i1
          i1pp=i1p+nbaslarge                                             2d21s20
          jb1=ib1+nbral*i1                                               2d21s20
          do i2=0,nbral-1                                                2d21s20
           i2p=iblstor(ibc(jpair))+i2                                    2d21s20
           i2pp=i2p+nbaslarge                                            2d21s20
           lla=i2p+nbtot*i1p+1                                           2d21s20
           llat=i1p+nbtot*i2p+1                                         2d21s20
           llb=lla+nbaslarge*(nbtot+1)                                   2d21s20
           llbt=llat+nbaslarge*(nbtot+1)                                2d21s20
           ovr(lla)=bc(jb1+i2)*xnorm                                      2d21s20
           ovr(llb)=ovr(lla)                                               2d21s20
           ovr(llat)=ovr(lla)                                           2d21s20
           ovr(llbt)=ovr(lla)                                           2d21s20
          end do                                                         2d21s20
         end do                                                          2d21s20
        end if                                                          2d21s20
       else if(ibc(jpair+2).eq.1)then                                   2d20s20
c     small component overlap                                           2d20s20
        call onep(lbra8,bc(jbra2),bc(isnormb),bc(jbra5),bc(jbra6),      2d20s20
     $      bc(jbra7),lket8,bc(jket2),bc(isnormk),bc(jket5),bc(jket6),
     $       bc(jket7),idum,ib1,0,0,0,0,0,0,0,0,0,bc,ibc)               11d9s22
        nlb=2*lbra8+1                                                   2d21s20
        nlk=2*lket8+1
        xnorm=oneover4pi*sqrt(dfloat(nlb*nlk))
        if(ibc(jpair).eq.ibc(jpair+1))then                              2d21s20
         do i1=0,nkets-1                                                 2d21s20
          i1p=ibsstor(ibc(jpair+1))+i1+2*nbaslarge                       2d21s20
          i1pp=i1p+nbassmall                                             2d21s20
          jb1=ib1+nbras*i1                                               2d21s20
          do i2=0,nbras-1                                                2d21s20
           i2p=ibsstor(ibc(jpair))+i2+2*nbaslarge                        2d21s20
           i2pp=i2p+nbassmall                                            2d21s20
           lssa=i2p+nbtot*i1p+1                                          2d21s20
           lssb=lssa+nbassmall*(nbtot+1)                                 2d21s20
           ovr(lssa)=bc(jb1+i2)*xnorm                                     2d21s20
           ovr(lssb)=ovr(lssa)                                           2d21s20
           term=-2d0*clight*clight*ovr(lssa)                            2d21s20
           h0(lssa)=h0(lssa)+term                                       2d21s20
           h0(lssb)=h0(lssb)+term                                       2d21s20
          end do                                                         2d21s20
         end do                                                          2d21s20
        else                                                            2d21s20
         do i1=0,nkets-1                                                 2d21s20
          i1p=ibsstor(ibc(jpair+1))+i1+2*nbaslarge                       2d21s20
          i1pp=i1p+nbassmall                                             2d21s20
          jb1=ib1+nbras*i1                                               2d21s20
          do i2=0,nbras-1                                                2d21s20
           i2p=ibsstor(ibc(jpair))+i2+2*nbaslarge                        2d21s20
           i2pp=i2p+nbassmall                                            2d21s20
           lssa=i2p+nbtot*i1p+1                                          2d21s20
           lssat=i1p+nbtot*i2p+1                                        2d21s20
           lssb=lssa+nbassmall*(nbtot+1)                                 2d21s20
           lssbt=lssat+nbassmall*(nbtot+1)                              2d21s20
           ovr(lssa)=bc(jb1+i2)*xnorm                                     2d21s20
           ovr(lssb)=ovr(lssa)                                           2d21s20
           ovr(lssat)=ovr(lssa)                                         2d21s20
           ovr(lssbt)=ovr(lssa)                                         2d21s20
           term=-2d0*clight*clight*ovr(lssa)                            2d21s20
           h0(lssa)=h0(lssa)+term                                       2d21s20
           h0(lssb)=h0(lssb)+term                                       2d21s20
          end do                                                         2d21s20
         end do                                                          2d21s20
        end if                                                          2d21s20
       else if(ibc(jpair+2).le.7)then                                   2d20s20
        idx=0                                                           2d20s20
        idy=0                                                           2d20s20
        idz=0                                                           2d20s20
        if(ibc(jpair+2).le.4)then
c     small-large coupling
         if(ibc(jpair+2).eq.2)idx=1                                      2d20s20
         if(ibc(jpair+2).eq.3)idy=1                                      2d20s20
         if(ibc(jpair+2).eq.4)idz=1                                      2d20s20
         nlb=2*lbra8+1                                                  2d21s20
         nlk=2*ibc(jket)+1                                              2d21s20
         xnorm=oneover4pi*sqrt(dfloat(nlb*nlk))*clight                  2d21s20
         call onep(lbra8,bc(jbra2),bc(isnormb),bc(jbra5),bc(jbra6),      2d20s20
     $       bc(jbra7),ibc(jket),bc(jket2),bc(jket3),bc(jket5),bc(jket6)2d20s20
     $       ,bc(jket7),idum,ib1,0,0,0,0,0,0,idx,idy,idz,bc,ibc)        11d9s22
         if(idx.ne.0)then                                               2d21s20
          do i1=0,nketl-1                                               2d21s20
           i1p=iblstor(ibc(jpair+1))+i1                                 2d21s20
           i1pp=i1p+nbaslarge                                           2d21s20
           jb1=ib1+nbras*i1                                             2d21s20
           do i2=0,nbras-1                                              2d21s20
            i2p=ibsstor(ibc(jpair))+i2+2*nbaslarge                      2d21s20
            i2pp=i2p+nbassmall                                          2d21s20
            lslab=i2p+nbtot*i1pp+1+nbtot*nbtot                          2d21s20
            lslba=i2pp+nbtot*i1p+1+nbtot*nbtot                          2d21s20
            term=-bc(jb1+i2)*xnorm                                      2d21s20
            h0(lslab)=h0(lslab)+term                                    2d21s20
            h0(lslba)=h0(lslba)+term                                    2d21s20
           end do                                                       2d21s20
          end do                                                        2d21s20
         else if(idy.ne.0)then                                          2d21s20
          do i1=0,nketl-1                                               2d21s20
           i1p=iblstor(ibc(jpair+1))+i1                                 2d21s20
           i1pp=i1p+nbaslarge                                           2d21s20
           jb1=ib1+nbras*i1                                             2d21s20
           do i2=0,nbras-1                                              2d21s20
            i2p=ibsstor(ibc(jpair))+i2+2*nbaslarge                      2d21s20
            i2pp=i2p+nbassmall                                          2d21s20
            lslab=i2p+nbtot*i1pp+1                                      2d21s20
            lslba=i2pp+nbtot*i1p+1                                      2d21s20
            term=bc(jb1+i2)*xnorm                                       2d21s20
            h0(lslab)=h0(lslab)-term                                    2d21s20
            h0(lslba)=h0(lslba)+term                                    2d21s20
           end do                                                       2d21s20
          end do                                                        2d21s20
         else                                                           2d21s20
          do i1=0,nketl-1                                               2d21s20
           i1p=iblstor(ibc(jpair+1))+i1                                 2d21s20
           i1pp=i1p+nbaslarge                                           2d21s20
           jb1=ib1+nbras*i1                                             2d21s20
           do i2=0,nbras-1                                              2d21s20
            i2p=ibsstor(ibc(jpair))+i2+2*nbaslarge                      2d21s20
            i2pp=i2p+nbassmall                                          2d21s20
            lslaa=i2p+nbtot*i1p+1+nbtot*nbtot                           2d21s20
            lslbb=i2pp+nbtot*i1pp+1+nbtot*nbtot                         2d21s20
            term=-bc(jb1+i2)*xnorm                                      2d21s20
            h0(lslaa)=h0(lslaa)+term                                    2d21s20
            h0(lslbb)=h0(lslbb)-term                                    2d21s20
           end do                                                       2d21s20
          end do                                                        2d21s20
         end if                                                         2d21s20
        else if(ibc(jpair).ne.ibc(jpair+1))then                         2d23s20
c     large-small coupling
c     we will be storing transpose so this falls below the diagonal     2d21s20
         if(ibc(jpair+2).eq.5)idx=1                                      2d20s20
         if(ibc(jpair+2).eq.6)idy=1                                      2d20s20
         if(ibc(jpair+2).eq.7)idz=1                                      2d20s20
         nlb=2*ibc(jbra)+1                                              2d21s20
         nlk=2*lket8+1                                                  2d21s20
         xnorm=oneover4pi*sqrt(dfloat(nlb*nlk))*clight                  2d21s20
         call onep(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),   2d21s20
     $       bc(jbra7),lket8,bc(jket2),bc(isnormk),bc(jket5),bc(jket6)  2d21s20
     $       ,bc(jket7),idum,ib1,0,0,0,0,0,0,idx,idy,idz,bc,ibc)        11d9s22
         if(idx.ne.0)then                                               2d21s20
          do i1=0,nkets-1                                               2d21s20
           i1p=ibsstor(ibc(jpair+1))+i1+2*nbaslarge                     2d21s20
           i1pp=i1p+nbassmall                                           2d21s20
           jb1=ib1+nbral*i1                                             2d21s20
           do i2=0,nbral-1                                              2d21s20
            i2p=iblstor(ibc(jpair))+i2                                  2d21s20
            i2pp=i2p+nbaslarge                                          2d21s20
            llsab=i1pp+nbtot*i2p+1+nbtot*nbtot                          2d21s20
            llsba=i1p+nbtot*i2pp+1+nbtot*nbtot                          2d21s20
            term=+bc(jb1+i2)*xnorm                                      2d21s20
            h0(llsab)=h0(llsab)+term                                    2d21s20
            h0(llsba)=h0(llsba)+term                                    2d21s20
           end do                                                       2d21s20
          end do                                                        2d21s20
         else if(idy.ne.0)then                                          2d21s20
          do i1=0,nkets-1                                               2d21s20
           i1p=ibsstor(ibc(jpair+1))+i1+2*nbaslarge                     2d21s20
           i1pp=i1p+nbassmall                                           2d21s20
           jb1=ib1+nbral*i1                                             2d21s20
           do i2=0,nbral-1                                              2d21s20
            i2p=iblstor(ibc(jpair))+i2
            i2pp=i2p+nbaslarge                                          2d21s20
            llsab=i1pp+nbtot*i2p+1                                      2d21s20
            llsba=i1p+nbtot*i2pp+1                                      2d21s20
            term=bc(jb1+i2)*xnorm                                       2d21s20
            h0(llsab)=h0(llsab)-term                                    2d21s20
            h0(llsba)=h0(llsba)+term                                    2d21s20
           end do                                                       2d21s20
          end do                                                        2d21s20
         else                                                           2d21s20
          do i1=0,nkets-1                                               2d21s20
           i1p=ibsstor(ibc(jpair+1))+i1+2*nbaslarge                     2d21s20
           i1pp=i1p+nbassmall                                           2d21s20
           jb1=ib1+nbral*i1                                             2d21s20
           do i2=0,nbral-1                                              2d21s20
            i2p=iblstor(ibc(jpair))+i2
            i2pp=i2p+nbaslarge                                          2d21s20
            llsaa=i1p+nbtot*i2p+1+nbtot*nbtot                           2d21s20
            llsbb=i1pp+nbtot*i2pp+1+nbtot*nbtot                         2d21s20
            term=+bc(jb1+i2)*xnorm                                      2d21s20
            h0(llsaa)=h0(llsaa)+term                                    2d21s20
            h0(llsbb)=h0(llsbb)-term                                    2d21s20
           end do                                                       2d21s20
          end do                                                        2d21s20
         end if                                                         2d21s20
        end if                                                          2d21s20
       else if(ibc(jpair+2).le.7+natom)then                             2d20s20
c     large nuclear attraction
        ia=ibc(jpair+2)-7                                               2d20s20
        zm=-atnum(1,ia)                                                 2d21s20
        nlb=2*ibc(jbra)+1                                               2d21s20
        nlk=2*ibc(jket)+1                                               2d21s20
        xnorm=oneover4pi*oneover4pi*sqrt(dfloat(nlb*nlk))
        zm=zm*xnorm                                                     2d21s20
        call derid(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),   8d20s15
     $         bc(jbra7),ibc(jbra4),ibc(jket),bc(jket2),bc(jket3),      2d21s20
     $         bc(jket5),bc(jket6),bc(jket7),ibc(jket4),                8d20s15
     $         0,atnum(2,ia),atnum(3,ia),xcart(1,ia),xcart(2,ia),       8d20s15
     $         xcart(3,ia),0,                                           8d20s15
     $         0,atnum(2,ia),atnum(3,ia),xcart(1,ia),xcart(2,ia),       8d20s15
     $         xcart(3,ia),0,dum,1,idum,ib1,0,0,0,0,0,0,0,0,0,0,0,0,    2d21s20
     $         1,.false.,zm,bc,ibc)                                     11d14s22
        do i1=0,nketl-1                                                 2d21s20
         i1p=iblstor(ibc(jpair+1))+i1
         i1pp=i1p+nbaslarge                                             2d21s20
         jb1=ib1+nbral*i1                                               2d21s20
         do i2=0,nbral-1                                                2d21s20
          i2p=iblstor(ibc(jpair))+i2                                    2d21s20
          i2pp=i2p+nbaslarge                                            2d21s20
          lla=i2p+nbtot*i1p+1                                           2d21s20
          llb=lla+nbaslarge*(nbtot+1)                                   2d21s20
          h0(lla)=h0(lla)+bc(jb1+i2)                                    2d21s20
          h0(llb)=h0(llb)+bc(jb1+i2)                                    2d21s20
         end do                                                         2d21s20
        end do                                                          2d21s20
       else                                                             2d20s20
c     small nuclear attraction                                          2d20s20
        ia=ibc(jpair+2)-7-natom                                         2d20s20
        zm=-atnum(1,ia)                                                 2d21s20
        nlb=2*lbra8+1                                                   2d21s20
        nlk=2*lket8+1                                                   2d21s20
        xnorm=oneover4pi*oneover4pi*sqrt(dfloat(nlb*nlk))
        zm=zm*xnorm                                                     2d21s20
        call derid(lbra8,bc(jbra2),bc(isnormb),bc(jbra5),bc(jbra6),     2d21s20
     $         bc(jbra7),ibc(jbra4),lket8,bc(jket2),bc(isnormk),        2d21s20
     $         bc(jket5),bc(jket6),bc(jket7),ibc(jket4),                8d20s15
     $         0,atnum(2,ia),atnum(3,ia),xcart(1,ia),xcart(2,ia),       8d20s15
     $         xcart(3,ia),0,                                           8d20s15
     $         0,atnum(2,ia),atnum(3,ia),xcart(1,ia),xcart(2,ia),       8d20s15
     $         xcart(3,ia),0,dum,1,idum,ib1,0,0,0,0,0,0,0,0,0,0,0,0,    2d21s20
     $         1,.false.,zm,bc,ibc)                                     11d14s22
        do i1=0,nkets-1                                                 2d21s20
         i1p=ibsstor(ibc(jpair+1))+i1+2*nbaslarge                       2d21s20
         i1pp=i1p+nbassmall                                             2d21s20
         jb1=ib1+nbras*i1                                               2d21s20
         do i2=0,nbras-1                                                2d21s20
          i2p=ibsstor(ibc(jpair))+i2+2*nbaslarge                        2d21s20
          i2pp=i2p+nbassmall                                            2d21s20
          lssa=i2p+nbtot*i1p+1                                          2d21s20
          lssb=lssa+nbassmall*(nbtot+1)                                 2d21s20
          h0(lssa)=h0(lssa)+bc(jb1+i2)                                  2d21s20
          h0(lssb)=h0(lssb)+bc(jb1+i2)                                  2d21s20
         end do                                                         2d21s20
        end do                                                          2d21s20
       end if                                                           2d20s20
      end do                                                            2d19s10
      call second(time3)
      call dws_sync                                                     2d22s10
      call dws_gsumf(ovr,nbb)                                           5d3s10
      call dws_gsumf(h0,nbb)                                            5d3s10
      ibcoff=ibcoffo                                                    2d19s10
      return
      end                                                               2d19s10
