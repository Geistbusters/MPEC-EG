c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine move2min(bc,ibc,idata,nunique,cfile,ehf,i3x,iapairg,   3d13s23
     $     ibdat,ibstor,idenergy,ieigr,ih0ao,ionex,ioooo,ispt1,isstor,  3d13s23
     $     isymdws,ivecr,jmats,kmats,maxdws,mgmat,ncfile,potdws,        3d13s23
     $     nder)                                                        3d13s23
      implicit real*8 (a-h,o-z)                                         2d14s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      logical lpowell,lgrad,ldebug,lprint                               3d14s23
      integer*1 idogrado1(4)                                            6d18s22
      include "common.basis"                                            2d14s23
      include "common.input"                                            2d14s23
      include "common.mrci"
      include "common.store"                                            2d14s23
      include "common.print"                                            3d21s23
      equivalence (idogrado4,idogrado1(1))                              3d14s23
      parameter(ide=100)                                                3d23s23
      character*3 element(ide)                                          3d23s23
      data element/'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ',     3d23s23
     $     'F  ','Ne ','Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ', 3d23s23
     $     'K  ','Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ', 3d23s23
     $     'Cu ','Zn ','Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ', 3d23s23
     $     'Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ', 3d23s23
     $     'In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ', 3d23s23
     $     'Pr ','Nd ','Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ', 3d23s23
     $     'Tm ','Yb ','Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ', 3d23s23
     $     'Au ','Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ', 3d23s23
     $     'Ac ','Th ','Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ', 3d23s23
     $     'Es ','Fm '/                                                 3d23s23
      common/drsigncm/drsign                                            8d20s24
      idogrado4=idogrado                                                3d14s23
      lprint=mynowprog.eq.0                                             3d14s23
      ldebug=.false.                                                    3d14s23
      if(lprint)write(6,*)('hi, my name is move2min '),idogrado1        3d14s23
      natom3=natom*3                                                    3d13s23
      if(ldebug)then                                                    3d14s23
       write(6,*)('what is under ehf '),ehf
       write(6,*)('what is under idenergy: '),
     $      (bc(idenergy+i),i=0,nder-1)
       write(6,*)('coordinate transformation matrix:')
       call prntm2(bc(idata),natom3,nder,natom3)
      end if                                                            3d14s23
      itrans=idata+natom3*natom                                         3d13s23
      if(ldebug)then                                                    3d14s23
       write(6,*)('unique directions '),itrans
       call prntm2(bc(itrans),nder,nunique,nder)
       write(6,*)('product of the two: ')
      end if                                                            3d14s23
      iprod=ibcoff                                                      3d13s23
      ibcoff=iprod+natom3*nunique                                       3d13s23
      call enough('move2min.prod',bc,ibc)                               3d13s23
      call dgemm('n','n',natom3,nunique,nder,1d0,bc(idata),natom3,      3d13s23
     $     bc(itrans),nder,0d0,bc(iprod),natom3)                        3d13s23
      if(ldebug)then                                                    3d14s23
       call prntm2(bc(iprod),natom3,nunique,natom3)                      3d13s23
      end if                                                            3d14s23
      lpowell=.false.                                                    3d14s23
      lgrad=.false.
      if(idogrado1(4).eq.2)then                                         3d14s23
       lpowell=.true.                                                   3d14s23
       lgrad=.true.                                                     3d14s23
       if(lprint)write(6,*)('we will be minimizing the gradient ')      3d14s23
      else                                                              3d14s23
       if(lprint)write(6,*)('we will be minimizing the energy ')        3d14s23
      end if                                                            3d14s23
      if(lprint)write(6,*)('lpowell '),lpowell,('lgrad '),lgrad         3d14s23
      if(natom.gt.3)pi=acos(-1d0)                                       3d23s23
      if(lprint)then                                                    3d22s23
       sz=0d0
       do ib=1,natom                                                    3d22s23
        do ia=1,natom                                                   3d22s23
         sz=sz+coord(ia,ib)**2                                          3d22s23
        end do                                                          3d22s23
       end do                                                           3d22s23
       sz=sqrt(sz)/dfloat(natom)                                        3d22s23
       if(ldebug)write(6,*)('size of coord matrix '),sz                 4d12s23
       if(sz.lt.1d-2)then                                               3d22s23
        icpy=-1                                                         3d22s23
       else                                                             3d22s23
        if(ldebug)then                                                  4d12s23
         write(6,*)('coord: ')                                             3d22s23
         call prntm2(coord,natom,natom,ida)                                3d22s23
        end if                                                          4d12s23
        icpy=ibcoff                                                       3d22s23
        imve=icpy+natom*max(3,natom)                                      3d21s23
        ipvt=imve+natom*natom                                             3d22s23
        ibcoff=ipvt+natom                                                 3d22s23
        call enough('move2min.cpy',bc,ibc)                                3d22s23
        do i=0,natom*natom-1                                              3d22s23
         bc(imve+i)=0d0                                                   3d22s23
        end do                                                            3d22s23
        jcpy=icpy                                                         3d22s23
        do i=1,natom                                                      3d22s23
         ii=imve+(i-1)*(natom+1)                                          3d22s23
         bc(ii)=1d0                                                       3d22s23
         do j=1,natom                                                     3d22s23
          bc(jcpy)=coord(j,i)                                             3d22s23
          jcpy=jcpy+1                                                     3d22s23
         end do                                                           3d22s23
        end do                                                            3d22s23
        call lusolv(bc(icpy),natom,natom,bc(imve),natom,natom,ibc(ipvt),  3d22s23
     $       ierr,3)                                                    3d21s23
        if(ldebug)then                                                  4d12s23
         write(6,*)('coordi ')
         call prntm2(bc(imve),natom,natom,natom)
         write(6,*)('transpose dirs ')
        end if                                                          4d12s23
        n3=nunique*3
        itmp1=ibcoff
        itmp2=itmp1+n3*natom
        ibcoff=itmp2+n3*natom
        call enough('move2min.tmp1',bc,ibc)                             3d23s23
        do iu=0,nunique-1
         do ixyz=0,2                                                    3d23s23
          do ia=0,natom-1                                               3d23s23
           ifrom=iprod+ia+natom*(ixyz+3*iu)
           ito=itmp1+ixyz+3*(iu+nunique*ia)
           bc(ito)=bc(ifrom)
          end do
         end do
        end do
        if(ldebug)call prntm2(bc(itmp1),n3,natom,n3)                    4d12s23
        call dgemm('n','n',n3,natom,natom,1d0,bc(itmp1),n3,bc(imve),
     $       natom,0d0,bc(itmp2),n3,'move2min.tmp2')
        if(ldebug)write(6,*)('times coordi ')
        icm=natom-1
        do ixyz=0,2
         do iu=0,nunique-1
          if(ldebug)write(6,*)('for ixyz '),ixyz+1,(' unique '),iu+1
          sum=0d0
          do ia=0,natom-1
           iadxy=imve+ia+natom*icm
           iadyx=itmp1+ixyz+3*(iu+nunique*ia)
           term=bc(iadxy)*bc(iadyx)
           sum=sum+term
          end do
          if(ldebug)write(6,*)('contribution to cm motion: '),sum       4d12s23
         end do
        end do
        if(ldebug)call prntm2(bc(itmp2),n3,natom,n3)                    4d12s23
       end if                                                           3d22s23
      end if                                                            3d22s23
      idirs=ibcoff                                                      3d13s23
      icoord0=idirs+nunique*nunique                                      3d13s23
      icuse=icoord0+nunique                                             3d13s23
      icoord=icuse+nunique                                              3d13s23
      iccuse=icoord+nunique                                             3d13s23
      icstart=iccuse+natom3                                             3d13s23
      ibcoff=icstart+natom3                                             3d13s23
      call enough('move2min.idirs',bc,ibc)                              3d13s23
      jcstart=icstart-1                                                 3d13s23
      do ixyz=1,3                                                       3d13s23
       do i=1,natom                                                      3d13s23
        bc(jcstart+i)=xcart(ixyz,i)                                     3d13s23
       end do                                                           3d13s23
       jcstart=jcstart+natom                                            3d13s23
      end do                                                            3d13s23
      do i=0,nunique-1                                                  3d13s23
       iad=idirs+nunique*i                                              3d13s23
       do j=0,nunique-1                                                 3d13s23
        bc(iad+j)=0d0                                                   3d13s23
       end do                                                           3d13s23
       bc(iad+i)=1d0                                                    3d13s23
       bc(icoord+i)=0d0                                                 3d13s23
      end do                                                            3d13s23
      if(.not.lpowell)then                                              3d14s23
       igdir=ibcoff                                                     3d14s23
       idenergy0=igdir+nunique                                          3d14s23
       idxi=idenergy0+nunique                                           3d14s23
       idfi=idxi+nunique                                                3d14s23
       ibcoff=idfi+nunique                                              3d14s23
       call enough('move2min.gdir',bc,ibc)                              3d14s23
      end if
      step=1d-1                                                         3d13s23
      thrsx=1d-3                                                        3d13s23
      thrse=1d-6
      if(lgrad)then
       v=0d0
       do i=0,nder-1
        v=v+bc(idenergy+i)**2
       end do
       v=sqrt(v/dfloat(nder))
      else                                                              3d14s23
       v=ehf                                                             3d13s23
      end if                                                            3d14s23
      icasguess=3                                                       3d13s23
      iparahfcall=1
      maxcall=100
      ihessupdate=0                                                     3d23s23
      vbest=v                                                           4d12s23
    1 continue                                                          3d13s23
       if(iparahfcall.gt.maxcall)then
        if(lprint)write(6,*)('tooo many iterations for minimum!')
        call dws_synca
        call dws_finalize
        stop 'move2min'
       end if
       do i=0,nunique-1                                                 3d13s23
        bc(icoord0+i)=bc(icoord+i)                                      3d13s23
       end do                                                           3d13s23
       dex=0d0                                                          3d13s23
       if(lpowell)then
        itop=nunique
       else
        itop=1
        do i=0,nunique-1                                                3d14s23
         bc(idenergy0+i)=bc(idenergy+i)                                 3d14s23
        end do                                                          3d14s23
       end if
       do i=1,itop                                                      3d14s23
        v0=v                                                            3d13s23
        if(lpowell.and.lprint)write(6,*)('for direction '),i,           3d14s23
     $       (' starting gradient norm is '),v0                         4d12s23
        idone=1                                                         3d13s23
        copt=0d0                                                        3d13s23
        step=0.1d0
        if(.not.lpowell)then
         if(ldebug)then                                                 4d12s23
          write(6,*)('multiply dirs ')
          call prntm2(bc(idirs),nunique,nunique,nunique)
          write(6,*)('denergy ')
          call prntm2(bc(idenergy),nunique,1,nunique)
         end if                                                         4d12s23
         fff=-drsign                                                    8d21s24
         call dgemm('n','n',nunique,1,nunique,fff,bc(idirs),nunique,    8d21s24
     $        bc(idenergy),nunique,0d0,bc(igdir),nunique)
         step=1d0                                                       3d14s23
         if(ldebug)then                                                 4d12s23
          write(6,*)('igdir: ')
          call prntm2(bc(igdir),nunique,1,nunique)
          write(6,*)('coord ')
          call prntm2(bc(icoord),nunique,1,nunique)
         end if                                                         4d12s23
        end if                                                          3d14s23
    2   continue                                                        3d13s23
         call minimize(copt,v,step,idone,thrsx,thrse,lpowell,ldebug)    3d23s23
         if(lpowell)then
          iad=idirs+nunique*(i-1)                                        3d13s23
         else
          iad=igdir
         end if
         do j=0,nunique-1                                               3d13s23
          bc(icuse+j)=bc(icoord+j)+bc(iad+j)*copt                       3d13s23
         end do                                                         3d13s23
         if(ldebug)then                                                 3d14s23
          write(6,*)('direction set with this value of copt: '),copt
          call prntm2(bc(icuse),nunique,1,nunique)
         end if                                                         3d14s23
         call dgemm('n','n',natom3,1,nunique,1d0,bc(iprod),natom3,      3d13s23
     $         bc(icuse),nunique,0d0,bc(iccuse),natom3)                 3d13s23
         if(ldebug)then                                                 3d14s23
          write(6,*)('times prod ')
          call prntm2(bc(iccuse),3,natom,3)
         end if
         do j=0,natom3-1                                                3d13s23
          bc(iccuse+j)=bc(icstart+j)+bc(iccuse+j)                       3d13s23
         end do                                                         3d13s23
         if(ldebug)then                                                 3d14s23
          write(6,*)('with starting coordinates ')
          call prntm2(bc(iccuse),natom,3,natom)
         end if                                                         3d14s23
         jccuse=iccuse-1                                                3d13s23
         do ixyz=1,3                                                    3d13s23
          do j=1,natom                                                   3d13s23
           xcart(ixyz,j)=bc(jccuse+j)                                   3d13s23
          end do                                                        3d13s23
          jccuse=jccuse+natom                                           3d13s23
         end do                                                         3d13s23
         if(ldebug)then
          write(6,*)('moved to xcart ')
          call prntm2(xcart,3,natom,3)
         end if                                                         3d14s23
         call dws_bcast(xcart,natom*3)                                  3d22s23
         potn=0d0
         do ia=1,natom
          do ib=ia+1,natom                                                 2d14s23
           dist=sqrt((xcart(1,ia)-xcart(1,ib))**2                          2d14s23
     $           +(xcart(2,ia)-xcart(2,ib))**2                          2d14s23
     $           +(xcart(3,ia)-xcart(3,ib))**2)                         2d14s23
           potn=potn+((atnum(1,ia)*atnum(1,ib))/dist)                      2d14s23
          end do                                                           2d14s23
         end do
         if(ldebug)write(6,*)('computed potn: '),potn,potdws-potn       3d14s23
         potdws=potn                                                    3d13s23
         ngaus2=ngaus*2
         ngaus3=ngaus2+ngaus
         ngaus4=ngaus3+ngaus
         ngaus5=ngaus4+ngaus
         ngaus6=ngaus5+ngaus
         ngaus7=ngaus6+ngaus
         jbdat=ibdat-1
         rms=0d0
         do ig=1,ngaus
          iad=jbdat+ig                                                     2d14s23
          l=ibc(iad)
          iad=iad+ngaus
          xp=bc(iad)
          iad=iad+ngaus
          xnorm=bc(iad)
          iad=iad+ngaus
          ioff=ibc(iad)
          iad=iad+ngaus
          iadx=iad+ngaus3
          nhere=ibc(iadx)
          bc(iad)=xcart(1,nhere)
          iad=iad+ngaus
          bc(iad)=xcart(2,nhere)
          iad=iad+ngaus
          bc(iad)=xcart(3,nhere)
         end do
         if(lprint)then                                                 3d14s23
          write(6,*)(' ')                                               3d14s23
          write(6,*)('calling parahf with coordinates ')                3d14s23
          if(icpy.lt.1)then                                             3d23s23
           do ia=1,natom                                                 3d23s23
            nuc=nint(atnum(1,ia))                                        3d23s23
            write(6,10)ia,atnum(1,ia),element(nuc),(xcart(j,ia),j=1,3)   3d23s23
   10       format('>',i3,f5.0,x,a3,2x,3f15.8,2x,f10.4,2es11.3)         3d23s23
           end do                                                        3d23s23
          else                                                          3d23s23
           call dgemm('n','n',3,natom,natom,1d0,xcart,3,bc(imve),natom,    3d22s23
     $         0d0,bc(icpy),3,'move2min.fint1')                          3d22s23
           if(ldebug)then                                               4d12s23
            write(6,*)('internals ')
            call prntm2(bc(icpy),3,natom,3)
           end if                                                       4d12s23
           nvec=natom-1                                                    3d22s23
           irs=ibcoff                                                      3d22s23
           ithetas=irs+nvec                                                3d22s23
           iphs=ithetas+nvec-1                                             3d22s23
           ibcoff=iphs+max(0,nvec-2)                                       3d22s23
           call enough('move2min.rs',bc,ibc)                               3d22s23
           do iv=0,nvec-1                                                  3d22s23
            jcart=icpy+3*iv                                                3d22s23
            rme=sqrt(bc(jcart)**2+bc(jcart+1)**2+bc(jcart+2)**2)           3d22s23
            bc(irs+iv)=rme                                                 3d22s23
            if(iv.gt.0)then                                                3d22s23
             rme=1d0/rme                                                    3d22s23
             theme=acos(bc(jcart+2)*rme)                                   3d22s23
             bc(ithetas+iv-1)=theme                                        3d22s23
             if(iv.gt.1)then                                               3d22s23
              ff=sin(theme)*bc(irs+iv)                                     3d22s23
              phme=atan2(bc(jcart+1),bc(jcart))                            3d22s23
              tryy=ff*sin(phme)                                            3d22s23
              tryx=ff*cos(phme)                                            3d22s23
              if(abs(tryy-bc(jcart+1)).gt.1d-10.or.
     $         abs(tryx-bc(jcart)).gt.1d-10)then                         3d22s23
               phme=phme+pi
               tryy=ff*sin(phme)                                            3d22s23
               tryx=ff*cos(phme)                                            3d22s23
               if(abs(tryy-bc(jcart+1)).gt.1d-10.or.                       3d22s23
     $          abs(tryx-bc(jcart)).gt.1d-10)then                         3d22s23
                write(6,*)('unable to match phi! '),bc(irs+iv),theme,      3d22s23
     $           tryx,bc(jcart),tryy,bc(jcart+1)                        3d22s23
                stop 'move2min'                                            3d22s23
               end if                                                      3d22s23
              end if                                                       3d22s23
              bc(iphs+iv-2)=phme                                           3d22s23
             end if                                                        3d22s23
            end if                                                         3d22s23
           end do                                                          3d22s23
           write(6,*)('>rs '),(bc(irs+iv),iv=0,nvec-1)                      3d22s23
           if(nvec.gt.1)write(6,*)('>theta '),                          3d23s23
     $          (bc(ithetas+iv),iv=0,nvec-2)                            3d23s23
           if(nvec.gt.2)write(6,*)('>phi '),(bc(iphs+iv),iv=0,nvec-3)   3d23s23
           ibcoff=irs                                                   3d23s23
          end if                                                        3d23s23
          write(6,*)(' ')                                               3d14s23
         end if                                                         3d14s23
         call parahf(mynowprog,mynprocg,jmats,kmats,mgmat,              3d13s23
     $     nocc,i3x,itmax,ieigr,ivecr,v,ih0ao,maxdws,ispt1,icore,       3d13s23
     $     nbasis,isymdws,myguess,natom,ngaus,ibdat,isym,iapair,        5d3s10
     $     multh,ibc(ibstor),ibc(isstor),potdws,ioooo,ionex,iwerner,    7d9s13
     $     idavopt,1,nlzz,idorel,ascale,idogrado,ipropsym,inocan,       10d27s17
     $     iapairg,dynw,nbasb,nbasc,nwcont,lmax,xcart,numeminus,atnum,  1d10s19
     $     bc(iextradatad),nextradatx,idorel4c,smallest,bc(idenergy),   8d29s22
     $      iftype,fstgth,ndfld,bc,ibc,1)                               12d6s23
         iparahfcall=iparahfcall+1
         if(lgrad.and.lpowell)then                                                  3d14s23
          v=0d0
          do j=0,nder-1
           v=v+bc(idenergy+j)**2
          end do
          v=sqrt(v/dfloat(nder))
          if(v.le.vbest)then                                            4d12s23
           vbest=v                                                      4d12s23
           if(lprint)write(6,*)('>new gradient norm: '),v,              4d12s23
     $          (' < best so far')                                      4d12s23
          else                                                          4d12s23
           if(lprint)write(6,*)('>new gradient norm: '),v
          end if                                                        4d12s23
         else
          if(ldebug)then                                                3d14s23
           write(6,*)('new energy '),v
           write(6,*)('copt: '),copt,v-v0
          end if                                                        3d14s23
          if(.not.lpowell.and.abs(copt-1d0).lt.1d-8.and.v.lt.v0)idone=1 3d14s23
          if(ldebug)write(6,*)('idone is now '),idone                   4d12s23
         end if
         if(idone.eq.0)go to 2                                                       3d13s23
        if(lpowell)then
         iad=idirs+nunique*(i-1)                                         3d13s23
        else
         iad=igdir
        end if
        do j=0,nunique-1                                                3d13s23
         bc(icoord+j)=bc(icoord+j)+bc(iad+j)*copt                       3d13s23
        end do                                                          3d13s23
        de=v0-v                                                         3d13s23
        if(lprint.and.lpowell)write(6,*)('copt, de: '),copt,de          3d14s23
        if(de.gt.dex)then                                               3d13s23
         dex=de                                                         3d13s23
         iex=i                                                          3d13s23
        end if                                                          3d13s23
       end do                                                           3d13s23
       if(lpowell)then
        if(lprint)write(6,*)('direction with greatest lowering: '),iex, 3d14s23
     $       dex                                                        3d14s23
        iad=idirs+nunique*(iex-1)                                        3d13s23
        do i=0,nunique-1                                                 3d13s23
         bc(iad+i)=bc(icoord+i)-bc(icoord0+i)                            3d13s23
        end do
       else
        ihessupdate=1                                                   3d23s23
        dotxf=0d0
        do i=0,nunique-1                                                3d14s23
         bc(idxi+i)=bc(icoord+i)-bc(icoord0+i)                          3d14s23
         bc(idfi+i)=drsign*(bc(idenergy+i)-bc(idenergy0+i))             8d21s24
         dotxf=dotxf+bc(idxi+i)*bc(idfi+i)
         ip=i+1
         if(ldebug)write(6,22)ip,bc(icoord+i),bc(icoord0+i),bc(idxi+i),
     $        bc(idenergy+i),bc(idenergy0+i),bc(idfi+i)
   22    format(i5,3es15.7,5x,3es15.7)
        end do
        dotxfi=1d0/dotxf
        if(ldebug)write(6,*)('dotxf '),dotxf,dotxfi
        call dgemm('n','n',nunique,1,nunique,1d0,bc(idirs),nunique,
     $       bc(idfi),nunique,0d0,bc(igdir),nunique)                    3d14s23
        if(ldebug)then
         write(6,*)('H*df ')
         call prntm2(bc(igdir),nunique,1,nunique)
        end if
        dotfhf=0d0
        do i=0,nunique-1
         dotfhf=dotfhf+bc(igdir+i)*bc(idfi+i)
        end do
        dotfhfi=1d0/dotfhf
        if(ldebug)write(6,*)('dotfhf = '),dotfhf,dotfhfi
        do i=0,nunique-1                                                3d14s23
         iad=idirs+nunique*i                                            3d14s23
         f1=bc(idxi+i)*dotxfi                                           3d14s23
         f2=-bc(igdir+i)*dotfhfi                                        3d14s23
         do j=0,nunique-1                                               3d14s23
          bc(iad+j)=bc(iad+j)+bc(idxi+j)*f1+bc(igdir+j)*f2              3d14s23
         end do                                                         3d14s23
        end do                                                          3d14s23
        if(ldebug)then
         write(6,*)('after BFP update: ')
         call prntm2(bc(idirs),nunique,nunique,nunique)
         do i=0,nunique-1
          bc(igdir+i)=dotxfi*bc(idxi+i)-dotfhfi*bc(igdir+i)
         end do
        end if                                                          3d14s23
        do i=0,nunique-1                                                3d14s23
         iad=idirs+nunique*i                                            3d14s23
         f1=bc(igdir+i)*dotfhf                                          3d14s23
         do j=0,nunique-1                                               3d14s23
          bc(iad+j)=bc(iad+j)+f1*bc(igdir+j)                            3d14s23
         end do                                                         3d14s23
        end do                                                          3d14s23
        if(ldebug)then                                                  3d14s23
         write(6,*)('after BFGS update: ')
         call prntm2(bc(idirs),nunique,nunique,nunique)
        end if                                                          3d14s23
       end if
       if(.not.lpowell)dex=dex*1d1                                      3d14s23
       if(dex.gt.thrse)go to 1
       if(lprint)then
        write(6,*)('iterations for minimum converged ')
        write(6,*)('no. of calls = '),iparahfcall
        gnorm=0d0                                                       4d12s23
        do i=0,nder-1                                                   4d12s23
         gnorm=gnorm+bc(idenergy+i)**2                                  4d12s23
        end do                                                          4d12s23
        gnorm=sqrt(gnorm/dfloat(nder))                                  4d12s23
        write(6,*)('>rms gradient norm '),gnorm                         4d12s23
        write(6,*)('>coordinates at minimum: ')
        if(icpy.gt.0)then                                               3d22s23
         call dgemm('n','n',3,natom,natom,1d0,xcart,3,bc(imve),natom,    3d22s23
     $         0d0,bc(icpy),3,'move2min.fint')                          3d22s23
         nvec=natom-1                                                    3d22s23
         irs=ibcoff                                                      3d22s23
         ithetas=irs+nvec                                                3d22s23
         iphs=ithetas+nvec-1                                             3d22s23
         ibcoff=iphs+max(0,nvec-2)                                       3d22s23
         call enough('move2min.rs',bc,ibc)                               3d22s23
         if(natom.gt.3)pi=acos(-1d0)                                     3d22s23
         do iv=0,nvec-1                                                  3d22s23
          jcart=icpy+3*iv                                                3d22s23
          rme=sqrt(bc(jcart)**2+bc(jcart+1)**2+bc(jcart+2)**2)           3d22s23
          bc(irs+iv)=rme                                                 3d22s23
          if(iv.gt.0)then                                                3d22s23
           rme=1d0/rme                                                    3d22s23
           theme=acos(bc(jcart+2)*rme)                                   3d22s23
           bc(ithetas+iv-1)=theme                                        3d22s23
           if(iv.gt.1)then                                               3d22s23
            ff=sin(theme)*bc(irs+iv)                                     3d22s23
            phme=atan2(bc(jcart+1),bc(jcart))                            3d22s23
            tryy=ff*sin(phme)                                            3d22s23
            tryx=ff*cos(phme)                                            3d22s23
            if(abs(tryy-bc(jcart+1)).gt.1d-10.or.
     $         abs(tryx-bc(jcart)).gt.1d-10)then                         3d22s23
             phme=phme+pi
             tryy=ff*sin(phme)                                            3d22s23
             tryx=ff*cos(phme)                                            3d22s23
             if(abs(tryy-bc(jcart+1)).gt.1d-10.or.                       3d22s23
     $          abs(tryx-bc(jcart)).gt.1d-10)then                         3d22s23
              write(6,*)('unable to match phi! '),bc(irs+iv),theme,      3d22s23
     $           tryx,bc(jcart),tryy,bc(jcart+1)                        3d22s23
              stop 'move2min'                                            3d22s23
             end if                                                      3d22s23
            end if                                                       3d22s23
            bc(iphs+iv-2)=phme                                           3d22s23
           end if                                                        3d22s23
          end if                                                         3d22s23
         end do                                                          3d22s23
         write(6,*)('>rs '),(bc(irs+iv),iv=0,nvec-1)                      3d22s23
         if(nvec.gt.1)write(6,*)('>theta '),(bc(ithetas+iv),iv=0,nvec-2)  3d22s23
         if(nvec.gt.2)write(6,*)('>phi '),(bc(iphs+iv),iv=0,nvec-3)     3d23s23
        else
         do ia=1,natom                                                  3d23s23
          nuc=nint(atnum(1,ia))                                         3d23s23
          write(6,110)ia,atnum(1,ia),element(nuc),(xcart(j,ia),j=1,3)   3d23s23
  110     format('>',i3,f5.0,x,a3,2x,3f21.14,2x,f10.4,2es11.3)          3d23s23
         end do                                                         3d23s23
        end if                                                          3d22s23
       end if
      call dws_synca
      return
      end
