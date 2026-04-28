c mpec2.1 version theta. For copyright and Disclaimers, see start.f
      subroutine catorbs(nfname,nsymb,ngaus,natom,nwcont,nbasdws,       1d22s25
     $     nlorcont,nbasisc,nbasisp,isym,idoub,iacto,iapair,ibdat,      1d22s25
     $     ibstor,isstor,ibcode,iorbov,iorbao,nbasallp,bc,ibc)          1d22s25
      implicit real*8 (a-h,o-z)                                         1d21s25
      logical ldebug                                                    1d23s25
      include 'common.store'                                            1d21s25
      include 'common.basis'                                            1d22s25
      include 'common.print'                                            1d23s25
      dimension nsymb(*),ngaus(*),natom(*),nwcont(*),                   1d21s25
      integer*8 ibcode,ibdat,ibdatt,ibstor,ibstort,iextradata,iorbao
     $iorbov
     $ivecn
     $          ipropsym(6),isinfo(11),multh(8,8),ipropsymw(6),         1d21s25
     $          multhw(8,8),nbasdws(8,*),nlorcont(8,*),nbasisc(8,*),    1d21s25
     $          nbasisp(8,*),isym(3,8,*),idoub(8,*),iacto(8,*),         1d22s25
     $          iapair(3,*),ibdat(*),ibstor(*),isstor(*),ibcode(8,*),   1d22s25
     $          iorbov(8,*),iorbao(8,*),nbasallp(*),nadd(8),ivecn(8),   1d22s25
     $          ibcoden(8)                                              1d22s25
      ldebug=iprtr(35).ne.0                                             1d23s25
      lmaxx=0                                                           1d21s25
      numeminust=0                                                      1d21s25
      naa=1                                                             1d22s25
      eet=0d0                                                           1d22s25
      ngaust=0                                                          1d22s25
      nwcontt=0                                                         1d22s25
      nbasallpt=0                                                       1d22s25
      do i=1,nfname                                                     1d21s25
       iunit=6+i                                                        1d21s25
       read(iunit)nsymb(i),idorel,ngaus(i),natom(i),nwcont(i),          1d21s25
     $      numeminus,lmax,nbasallp(i),multh,ascale,potdws,ipropsym,    1d22s25
     $      isinfo,nextradatatag                                        1d21s25
       write(6,*)('from orbital file no. '),i
       ngaust=ngaust+ngaus(i)                                           1d22s25
       nwcontt=nwcontt+nwcont(i)                                        1d22s25
       nbasallpt=nbasallpt+nbasallp(i)                                  1d22s25
       if(i.eq.1)then                                                   1d21s25
        idorelw=idorel                                                  1d21s25
        ncomp=1                                                         1d22s25
        if(idorelw.ne.0)ncomp=2                                         1d22s25
        ascalew=ascale                                                  1d21s25
        nsymbw=nsymb(i)                                                 1d22s25
        do j=1,6                                                        1d21s25
         ipropsymw(j)=ipropsym(j)                                       1d22s25
        end do                                                          1d21s25
        do j=1,nsymbw                                                   1d21s25
         do k=1,nsymbw                                                  1d21s25
          multhw(k,j)=multh(k,j)                                        1d21s25
         end do                                                         1d21s25
        end do                                                          1d21s25
       end if                                                           1d21s25
       write(6,*)('nsymb etc '),nsymb(i),idorel,ngaus(i),natom(i),      1d21s25
     $      nwcont(i),numeminus,lmax,nbasallp(i),ascale                 1d22s25
       lmaxx=max(lmaxx,lmax)                                            1d21s25
       numeminust=numeminust+numeminus                                  1d21s25
       if(nsymb(i).ne.nsymbw)then                                       1d21s25
        write(6,*)('point group order changed!')                        1d21s25
        stop                                                             1d21s25
       end if                                                           1d21s25
       if(idorel.ne.idorelw)then                                        1d21s25
        write(6,*)('idorel is not consistent!')                         1d21s25
        stop                                                            1d21s25
       end if                                                           1d21s25
       if(abs(ascale-ascalew).gt.1d-10)then                             1d21s25
        write(6,*)('ascale is not consistent!')                         1d21s25
        stop                                                            1d21s25
       end if                                                           1d21s25
       nbad=0                                                           1d21s25
       do j=1,nsymbw                                                    1d21s25
        do k=1,nsymbw                                                   1d21s25
         if(multh(k,j).ne.multhw(k,j))nbad=nbad+1                       1d21s25
        end do                                                          1d21s25
       end do                                                           1d21s25
       if(nbad.ne.0)then                                                1d21s25
        write(6,*)('multiplication table not consistent1!')             1d21s25
        stop                                                            1d21s25
       end if                                                           1d21s25
       read(iunit)(nbasdws(isb,i),isb=1,nsymbw),                            1d21s25
     $      (nlorcont(isb,i),isb=1,nsymbw)                              1d21s25
       read(iunit)(nbasisc(isb,i),isb=1,nsymbw)                             1d21s25
       read(iunit)(nbasisp(isb,i),isb=1,nsymbw)                         1d21s25
       iextradata=ibcoff                                                1d22s25
       ibcoff=iextradata+nextradatatag                                  1d22s25
       read(iunit)((isym(j,k,i),j=1,3),k=1,8),                          1d22s25
     $      ((xcart(ixyz,ia),atnum(ixyz,ia),ixyz=1,3),                  1d22s25
     $    ia=naa,naa+natom(i)-1)
     $      ,(bc(iextradata+j),j=0,nextradatatag-1)                     1d22s25
     $      ,ee                                                         1d22s25
     $      ,(idoub(isb,i),iacto(isb,i),isb=1,nsymbw),nlzz              1d22s25
       ibcoff=iextradata                                                1d22s25
       write(6,*)('ee = '),ee
       write(6,*)('nlzz = '),nlzz
       read(iunit)((iapair(j,ia),j=1,3),ia=naa,naa+natom(i)-1)          1d22s25
       nbdat=ngaus(i)*9+nwcont(i)                                       1d22s25
       ibdat(i)=ibcoff                                                  1d22s25
       ibcoff=ibdat(i)+nbdat                                            1d22s25
       read(iunit)(bc(ibdat(i)+j),j=0,nbdat-1)                          1d22s25
       if(ldebug)write(6,*)('what we have for ibdat ')                  1d23s25
       do j=1,ngaus(i)                                                  1d22s25
        iad1=ibdat(i)+j-1
        iad2=iad1+ngaus(i)
        iad3=iad2+ngaus(i)
        iad4=iad3+ngaus(i)
        iad5=iad4+ngaus(i)
        iad6=iad5+ngaus(i)
        iad7=iad6+ngaus(i)
        iad8=iad7+ngaus(i)
        iad9=iad8+ngaus(i)
        if(ldebug)write(6,34)j,ibc(iad1),bc(iad2),bc(iad3),ibc(iad4),   1d23s25
     $       bc(iad5),bc(iad6),bc(iad7),ibc(iad8),ibc(iad9)             1d23s25
   34   format(i5,5x,i8,2es15.7,i8,3es15.7,2i8)
       end do
       if(ldebug)then                                                   1d23s25
        write(6,*)('start of nwcont '),iad9+1-ibdat(i)
        call prntm2(bc(iad9+1),1,nwcont(i),1)
       end if                                                           1d23s25
       ibstor(i)=ibcoff                                                 1d22s25
       isstor(i)=ibstor(i)+nbasallp(i)*3                                1d22s25
       ibcoff=isstor(i)+nbasallp(i)*3                                   1d22s25
       read(iunit)(ibc(ibstor(i)+j),ibc(isstor(i)+j),                   1d22s25
     $      j=0,nbasallp(i)*3-1)                                        1d22s25
       do isb=1,nsymbw                                                  1d22s25
        ibcode(isb,i)=ibcoff                                                 1d22s25
        ibcoff=ibcode(isb,i)+nbasdws(isb,i)                             1d22s25
        read(iunit)(ibc(ibcode(isb,i)+j),j=0,nbasdws(isb,i)-1)          1d22s25
        if(ldebug)then                                                  1d23s25
         write(6,*)('bcode for symmetry block '),isb
         write(6,*)(ibc(ibcode(isb,i)+j),j=0,nbasdws(isb,i)-1)           1d22s25
        end if                                                          1d23s25
       end do                                                           1d22s25
       do isb=1,nsymbw                                                  1d22s25
        nrow=nbasisp(isb,i)*ncomp                                       1d22s25
        if(nrow.gt.0)then                                               1d22s25
         iorbao(isb,i)=ibcoff                                           1d22s25
         ibcoff=iorbao(isb,i)+nbasdws(isb,i)*nrow                       1d22s25
         read(iunit)(bc(iorbao(isb,i)+j),j=0,nbasdws(isb,i)*nrow-1)     1d22s25
         iorbov(isb,i)=ibcoff                                           1d22s25
         ibcoff=iorbov(isb,i)+nbasisc(isb,i)*nbasdws(isb,i)             1d22s25
         read(iunit)(bc(iorbov(isb,i)+j),j=0,                           1d22s25
     $        nbasisc(isb,i)*nbasdws(isb,i)-1)                          1d22s25
        else                                                            1d22s25
         read(iunit)                                                    1d22s25
         read(iunit)                                                    1d22s25
        end if                                                          1d22s25
       end do                                                           1d22s25
       eet=eet+ee                                                       1d22s25
       naa=naa+natom(i)                                                 1d22s25
      end do                                                            1d21s25
      naa=naa-1                                                         1d22s25
      write(6,*)('sum of energies '),eet                                1d22s25
      write(6,*)('total no. of atoms: '),naa
      write(6,*)('total no. of electrons '),numeminust
      if(ldebug)then                                                    1d23s25
       write(6,*)('cartesian coordinates ')
       call prntm2(xcart,3,naa,3)
       write(6,*)('atnum ')
       call prntm2(atnum,3,naa,3)
       write(6,*)('building full ibdat ')                                1d22s25
      end if                                                            1d23s25
      ibdatt=ibcoff                                                     1d22s25
      ibcoef=ibdatt+9*ngaust                                            1d22s25
      nbdatt=9*ngaust+nwcontt                                           1d22s25
      ibcoff=ibcoef+nwcontt                                             1d22s25
      jbdatt=ibdatt                                                     1d22s25
      jbcoef=ibcoef                                                     1d22s25
      nas=0                                                             1d22s25
      lhere=-1                                                          1d22s25
      iahere=-1                                                         1d22s25
      ntype=0                                                           1d22s25
      nbb=0                                                             1d23s25
      do i=1,nfname                                                     1d22s25
       if(ldebug)write(6,*)('from file '),i                             1d23s25
       do j=1,ngaus(i)                                                  1d22s25
        iad1=ibdat(i)+j-1                                               1d22s25
        iad2=iad1+ngaus(i)
        iad3=iad2+ngaus(i)
        iad4=iad3+ngaus(i)
        iad5=iad4+ngaus(i)
        iad6=iad5+ngaus(i)
        iad7=iad6+ngaus(i)
        iad8=iad7+ngaus(i)
        iad9=iad8+ngaus(i)
        iat1=jbdatt+j-1                                                 1d22s25
        iat2=iat1+ngaust                                                1d22s25
        iat3=iat2+ngaust                                                1d22s25
        iat4=iat3+ngaust                                                1d22s25
        iat5=iat4+ngaust                                                1d22s25
        iat6=iat5+ngaust                                                1d22s25
        iat7=iat6+ngaust                                                1d22s25
        iat8=iat7+ngaust                                                1d22s25
        iat9=iat8+ngaust                                                1d22s25
        ibc(iat1)=ibc(iad1)                                             1d22s25
        bc(iat2)=bc(iad2)                                               1d22s25
        bc(iat3)=bc(iad3)                                               1d22s25
        ibc(iat4)=ibc(iad4)+nbb                                         1d23s25
        bc(iat5)=bc(iad5)                                               1d22s25
        bc(iat6)=bc(iad6)                                               1d22s25
        bc(iat7)=bc(iad7)                                               1d22s25
        ibc(iat8)=ibc(iad8)+nas                                         1d22s25
        if(ibc(iad1).eq.lhere.and.ibc(iat8).eq.iahere)then              1d22s25
         ntype=ntype+1                                                  1d22s25
        else                                                            1d22s25
         if(ntype.gt.0)then                                             1d22s25
          istart=ibdat(i)+ibc(itt2)                                     1d22s25
          if(bc(istart).ne.-132d0)then                                  1d22s25
           inew=jbcoef-ibdatt                                           1d22s25
           ntypep=ntype*ncomp                                           1d24s25
           if(ldebug)then                                               1d23s25
            write(6,*)('we have a new type ...')                           1d22s25
            write(6,*)('l = '),lhere
            write(6,*)('iatom no. '),iahere
            write(6,*)('number of primitive fcns: '),ntype
            write(6,*)('number of contracted fcns: '),ncf                  1d22s25
            write(6,*)('contraction coefficients at '),ibc(itt2)           1d22s25
            write(6,*)('new address: '),inew,jbcoef
            call prntm2(bc(istart),ntypep,ncf,ntypep)                   1d24s25
           end if                                                       1d23s25
           do k=0,ntypep*ncf-1                                          1d24s25
            bc(k+jbcoef)=bc(istart+k)                                   1d22s25
           end do                                                       1d22s25
           jbcoef=jbcoef+ntypep*ncf                                     1d24s25
          end if                                                        1d22s25
          ibc(iat9)=inew                                                1d22s25
          bc(istart)=-132d0                                             1d22s25
         end if                                                         1d22s25
        end if                                                          1d22s25
        if(ibc(iad9).gt.0)then                                          1d22s25
         if(ldebug)write(6,*)('we have contract offset '),ibc(iad9)     1d23s25
         itt1=iat9                                                      1d22s25
         itt2=iad9                                                      1d22s25
         lhere=ibc(iad1)                                                1d22s25
         iahere=ibc(iat8)                                               1d22s25
         ntype=1                                                        1d22s25
        else                                                            1d22s25
         if(ibc(iad9).lt.-1)then                                         1d22s25
          ncf=-ibc(iad9)                                                1d22s25
          if(ldebug)write(6,*)('number of contracted functions '),ncf   1d23s25
         end if                                                         1d22s25
         ibc(iat9)=ibc(iad9)                                            1d22s25
        end if                                                          1d22s25
       end do                                                           1d22s25
       if(ldebug)then                                                   1d23s25
        write(6,*)('we have reached the end of gausian loop ')
        write(6,*)('ntype = '),ntype
       end if                                                           1d23s25
       if(ntype.gt.0)then                                               1d22s25
        istart=ibdat(i)+ibc(itt2)                                       1d22s25
        if(bc(istart).ne.-132d0)then
         inew=jbcoef-ibdatt                                             1d22s25
         ntypep=ntype*ncomp                                             1d24s25
         if(ldebug)then                                                 1d23s25
          write(6,*)('we have a new type ...')                            1d22s25
          write(6,*)('l = '),lhere
          write(6,*)('iatom no. '),iahere
          write(6,*)('number of primitive fcns: '),ntype
          write(6,*)('number of contracted fcns: '),ncf                   1d22s25
          write(6,*)('contraction coefficients at '),ibc(itt2)            1d22s25
          call prntm2(bc(istart),ntypep,ncf,ntypep)                     1d24s25
          write(6,*)('new address: '),inew,jbcoef
         end if                                                         1d23s25
         do k=0,ntypep*ncf-1                                            1d24s25
          bc(k+jbcoef)=bc(istart+k)                                     1d22s25
         end do                                                         1d22s25
         jbcoef=jbcoef+ntypep*ncf                                       1d24s25
        end if
        ibc(iat9)=inew                                                  1d22s25
       end if                                                           1d22s25
       lhere=-1                                                          1d22s25
       iahere=-1                                                         1d22s25
       ntype=0                                                           1d22s25
       nas=nas+natom(i)                                                 1d22s25
       jbdatt=jbdatt+ngaus(i)                                           1d22s25
       nbb=nbb+nbasallp(i)                                              1d23s25
      end do                                                            1d22s25
      if(ldebug)then                                                    1d23s25
       write(6,*)('altogether now ')
       do i=1,ngaust
        iat1=ibdatt+i-1                                                  1d22s25
        iat2=iat1+ngaust                                                 1d22s25
        iat3=iat2+ngaust                                                 1d22s25
        iat4=iat3+ngaust                                                 1d22s25
        iat5=iat4+ngaust                                                 1d22s25
        iat6=iat5+ngaust                                                 1d22s25
        iat7=iat6+ngaust                                                 1d22s25
        iat8=iat7+ngaust                                                 1d22s25
        iat9=iat8+ngaust                                                 1d22s25
        write(6,34)i,ibc(iat1),bc(iat2),bc(iat3),ibc(iat4),bc(iat5),
     $       bc(iat6),bc(iat7),ibc(iat8),ibc(iat9)
       end do
       write(6,*)('now update bstor and sstor ')                         1d22s25
      end if                                                            1d23s25
      ibstort=ibcoff                                                    1d22s25
      isstort=ibstort+nbasallpt*3                                       1d22s25
      ibcoff=isstort+nbasallpt*3                                        1d22s25
      jbstort=ibstort                                                   1d22s25
      jsstort=isstort                                                   1d22s25
      do isb=1,nsymbw                                                   1d22s25
       nadd(isb)=0                                                      1d22s25
      end do                                                            1d22s25
      do i=1,nfname                                                     1d22s25
       if(ldebug)then                                                   1d23s25
        write(6,*)('from file '),i
        write(6,*)(nbasisp(isb,i),isb=1,nsymbw)
       end if                                                           1d23s25
       do j=1,nbasallp(i)                                               1d22s25
        jm=j-1                                                          1d22s25
        isb=ibc(isstor(i)+jm)                                           1d22s25
        ibc(jbstort)=ibc(ibstor(i)+jm)+nadd(isb)                        1d22s25
        ibc(jsstort)=isb                                                1d22s25
        if(ldebug)                                                      1d23s25
     $      write(6,*)j,ibc(ibstor(i)+jm),ibc(isstor(i)+jm),ibc(jbstort)1d23s25
        jbstort=jbstort+1                                               1d22s25
        jsstort=jsstort+1                                               1d22s25
       end do                                                           1d22s25
       do isb=1,nsymbw                                                  1d22s25
        nadd(isb)=nadd(isb)+nbasisp(isb,i)                              1d22s25
       end do                                                           1d22s25
      end do                                                            1d22s25
      if(ldebug)then                                                    1d23s25
       write(6,*)('final bstor,sstor ')
       do j=1,nbasallpt
        jm=j-1
        write(6,*)j,ibc(ibstort+jm),ibc(isstort+jm)
       end do
      end if                                                            1d23s25
      iat9=iat9+1
      if(ldebug)then                                                    1d23s25
       write(6,*)('starting '),iat9-ibdatt,iat9
       call prntm2(bc(iat9),1,nwcontt,1)
       write(6,*)('now for orbs in ao basis ...')
      end if                                                            1d23s25
      do isb=1,nsymbw                                                   1d22s25
       if(ldebug)write(6,*)('for symmetry '),isb                        1d23s25
       nptot=0                                                          1d22s25
       nctot=0                                                          1d22s25
       nbtot=0                                                          1d22s25
       do i=1,nfname                                                    1d22s25
        nptot=nptot+nbasisp(isb,i)                                      1d22s25
        nctot=nctot+nbasisc(isb,i)                                      1d22s25
        nbtot=nbtot+nbasdws(isb,i)                                      1d22s25
       end do                                                           1d22s25
       nptotp=nptot*ncomp                                               1d22s25
       if(ldebug)then                                                   1d23s25
        write(6,*)('total no. primitive fcns: '),nptot                   1d22s25
        write(6,*)('total no. contracted fcns: '),nctot                  1d22s25
        write(6,*)('total no. of vectors : '),nbtot                      1d22s25
        write(6,*)('first dim of ao vector matrix '),nptotp              1d22s25
       end if                                                           1d23s25
       ivecn(isb)=ibcoff                                                1d22s25
       ibcoden(isb)=ivecn(isb)+nptotp*nbtot                             1d22s25
       ibcoff=ibcoden(isb)+nbtot                                        1d22s25
       jvecn=ivecn(isb)                                                 1d22s25
       jbcoden=ibcoden(isb)                                             1d22s25
       do iz=ivecn(isb),ibcoff-1                                        1d22s25
        bc(iz)=0d0                                                      1d22s25
       end do                                                           1d22s25
       do ipass=1,4                                                     1d22s25
        iwant=ipass                                                     1d22s25
        if(ipass.eq.4)then                                              1d22s25
         iwant=0                                                        1d22s25
        end if                                                          1d22s25
        if(ldebug)write(6,*)('for bcode = '),iwant                      1d23s25
        ioff=0                                                          1d22s25
        do i=1,nfname                                                   1d22s25
         if(ldebug)write(6,*)('from file '),i                           1d23s25
         nnn=nbasisp(isb,i)*ncomp                                       1d22s25
         do j=0,nbasdws(isb,i)-1                                        1d22s25
          itry=ibc(ibcode(isb,i)+j)                                     1d22s25
          if(itry.eq.iwant)then                                         1d22s25
           if(ldebug)write(6,*)('fcn '),j+1,('matches bcode ')          1d23s25
           ibc(jbcoden)=iwant                                           1d22s25
           jbcoden=jbcoden+1                                            1d22s25
           iad1=iorbao(isb,i)+nnn*j                                     1d22s25
           iad2=jvecn+ioff                                              1d22s25
           do k=0,nbasisp(isb,i)-1                                      1d22s25
            bc(iad2+k)=bc(iad1+k)                                       1d22s25
           end do                                                       1d22s25
           if(ncomp.eq.2)then                                           1d22s25
            iad2=iad2+nptot                                             1d22s25
            iad1=iad1+nbasisp(isb,i)                                    1d22s25
            do k=0,nbasisp(isb,i)-1                                     1d22s25
             bc(iad2+k)=bc(iad1+k)                                      1d22s25
            end do                                                      1d22s25
           end if                                                       1d22s25
           jvecn=jvecn+nptotp                                           1d22s25
          end if                                                        1d22s25
         end do                                                         1d22s25
         ioff=ioff+nbasisp(isb,i)                                       1d22s25
        end do                                                          1d22s25
       end do                                                           1d22s25
       if(ldebug)then                                                   1d23s25
        write(6,*)('new vectors ')
        nbtotp=nbtot*ncomp                                              1d24s25
        call prntm2(bc(ivecn(isb)),nptot,nbtotp,nptot)                  1d22s25
       end if                                                           1d23s25
      end do                                                            1d22s25
      nextradatatag=0                                                   1d22s25
      write(6,*)('save to forbs')                                       1d22s25
      write(4)nsymbw,idorelw,ngaust,natomt,nwcontt,numeminust,lmaxx,    1d22s25
     $     nbasallpt,multh,ascale,potdws,ipropsym,isinfo,nextradatatag  1d22s25
      do i=2,nfname                                                     1d22s25
       do isb=1,nsymbw                                                  1d22s25
        nbasdws(isb,1)=nbasdws(isb,1)+nbasdws(isb,i)                    1d22s25
        nbasisc(isb,1)=nbasisc(isb,1)+nbasisc(isb,i)                    1d22s25
        nbasisp(isb,1)=nbasisp(isb,1)+nbasisp(isb,i)                    1d22s25
        idoub(isb,1)=idoub(isb,1)+idoub(isb,i)                          1d22s25
        iacto(isb,1)=iacto(isb,1)+iacto(isb,i)                          1d22s25
       end do                                                           1d22s25
      end do                                                            1d22s25
      write(6,*)('concatenated idoub: '),(idoub(isb,1),isb=1,nsymbw)    1d23s25
      write(6,*)('concatenated iacto: '),(iacto(isb,1),isb=1,nsymbw)    1d23s25
      write(4)(nbasdws(isb,1),isb=1,nsymbw),                            1d21s25
     $      (nlorcont(isb,1),isb=1,nsymbw)                              1d21s25
      write(4)(nbasisc(isb,1),isb=1,nsymbw)                             1d21s25
      write(4)(nbasisp(isb,1),isb=1,nsymbw)                             1d22s25
      write(4)((isym(j,k,1),j=1,3),k=1,8),                              1d22s25
     $     ((xcart(ixyz,ia),atnum(ixyz,ia),ixyz=1,3),ia=1,naa),eet,
     $     (idoub(isb,1),iacto(isb,1),isb=1,nsymbw),nlzz                1d22s25
      write(4)((iapair(j,ia),j=1,3),ia=1,naa)                           1d22s25
      write(4)(bc(ibdatt+j),j=0,nbdatt-1)                               1d22s25
      write(4)(ibc(ibstort+j),ibc(isstort+j),j=0,nbasallpt*3-1)         1d22s25
      if(ldebug)then                                                    1d23s25
       write(6,*)('bstor saved ')
       do j=1,nbasallpt
        jm=j-1
        write(6,*)j,ibc(ibstort+jm),ibc(isstort+jm)
       end do
      end if                                                            1d23s25
      do isb=1,nsymbw                                                   1d22s25
       write(4)(ibc(ibcoden(isb)+j),j=0,nbasdws(isb,1)-1)               1d22s25
      end do                                                            1d22s25
      do isb=1,nsymbw                                                   1d22s25
       nrow=nbasisp(isb,1)*ncomp                                        1d22s25
       if(nrow.gt.0)then                                                1d22s25
        write(4)(bc(ivecn(isb)+j),j=0,nbasdws(isb,1)*nrow-1)            1d22s25
        write(4)(0d0,j=0,nbasisc(isb,1)*nbasdws(isb,1)-1)               1d22s25
       else                                                             1d22s25
        write(4)                                                        1d22s25
        write(4)                                                        1d22s25
       end if                                                           1d22s25
      end do                                                            1d22s25
      return                                                            1d22s25
      end                                                               1d22s25
