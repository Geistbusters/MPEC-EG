c mpec2.1 version theta. For copyright and Disclaimers, see start.f
      subroutine prebodcfd(ifname,nfname,natom,nuniq,idatal,ndatal,
     $     ibdatf,nbdat,bc,ibc)                                         2d3s25
      implicit real*8 (a-h,o-z)                                         2d3s25
      character*200 fname                                               11d22s19
      integer*8 ipack8                                                  2d3s25
      integer*4 ipack4(2)                                               2d3s25
      equivalence (ipack8,ipack4)                                       2d3s25
      include "common.store"                                            2d3s25
      dimension nameci(6),idumx18(18),multhx(8,8),isymx(3,8),ibdatf(*)  2d3s25
      integer*8 iawgt,ibcoff0,icart,icartl,idatal,itmp1
      write(6,*)('Hi, my name is prebodcfd ')                           2d3s25
      write(6,*)('no. of atoms = '),natom
      write(6,*)('no. of files = '),nfname
      write(6,*)('let us collect coordinates ')
      icart=ibcoff                                                      2d3s25
      ibcoff=icart+3*natom*nfname                                       2d3s25
      call enough('prebodcfd.cart',bc,ibc)                              2d3s25
      jfname=ifname                                                     2d3s25
      do i=1,nfname                                                     2d3s25
       ibcoff0=ibcoff                                                   2d3s25
       im=i-1                                                           2d3s25
       write(6,*)('for file no. '),i
       nhere=ibc(jfname)
       jfname=jfname+1                                                  2d13s25
       if(nhere.gt.200)then
        write(6,*)('name is toooo long!')                               2d3s25
        stop 'prebodcfd'                                                2d3s25
       end if                                                           2d3s25
       do j=1,nhere
        fname(j:j)=char(ibc(jfname))                                    11d22s19
        jfname=jfname+1
       end do                                                           11d22s19
       write(6,*)('file name: '),fname(1:nhere)                         2d3s25
       open(unit=1,file=fname(1:nhere),form='unformatted')              2d3s25
       read(1)isumlt,isymmrci,nroot,norb,nsymb,nameci                   2d3s25
       read(1)idum                                                      2d3s25
       if(norb.gt.0)then                                                2d3s25
        read(1)idum                                                      2d3s25
       end if                                                           2d3s25
       read(1)nsymbi,idorelx,ngaus,natom,nwcont,numeminus,lmax,nbasallp,8d10s22
     $     multhx,ascale,potdws,idumx18                                 5d25s21
       read(1)
c     nbasc
       read(1)
       read(1)
       write(6,*)('nsymb,idorel: '),nsymbi,idorelx
       write(6,*)('ngauss,natom: '),ngaus,natom
       itmp1=ibcoff                                                     2d3s25
       itmp2=itmp1+natom*3                                              2d3s25
       ibcoff=itmp2+natom*2                                             2d3s25
       call enough('prebodcfd.tmp1',bc,ibc)                             2d3s25
       read(1)isymx,(bc(itmp1+ii),bc(itmp2+ii),ii=0,natom*3-1)          2d12s25
       jtmp1=itmp1                                                      2d3s25
       do ii=0,natom-1                                                  2d12s25
        jcart=icart+3*(ii+natom*im)                                     2d12s25
        do ixyz=0,2                                                     2d3s25
         bc(jcart+ixyz)=bc(jtmp1+ixyz)                                  2d3s25
        end do                                                          2d3s25
        jtmp1=jtmp1+3                                                   2d3s25
       end do                                                           2d3s25
       read(1)dum
       ibcoff=itmp1
       nbdat=ngaus*9+nwcont
       ibcoff=itmp1+nbdat
       call enough('getwf.  2',bc,ibc)
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)(bc(itmp1+ii),ii=0,nbdat-1)                               2d12s25
       ibdatf(i)=itmp1                                                  2d3s25
      end do                                                            2d3s25
      ibcoff0=ibcoff                                                    2d3s25
      write(6,*)('cartesians read in ')                                 2d3s25
      natom3=natom*3                                                    2d3s25
      call prntm2(bc(icart),natom3,nfname,natom3)                       2d3s25
      write(6,*)('wrt reference coordinates ')                          2d3s25
      nfnamem=nfname-1                                                  2d3s25
      do i=1,nfnamem                                                    2d3s25
       jcart=icart+natom3*i                                             2d3s25
       do j=0,natom3-1                                                  2d3s25
        bc(jcart+j)=bc(jcart+j)-bc(icart+j)                             2d3s25
       end do                                                           2d3s25
      end do                                                            2d3s25
      icartp=icart+natom3                                               2d3s25
      call prntm2(bc(icartp),natom3,nfnamem,natom3)                       2d3s25
      icartl=ibcoff                                                     2d3s25
      ibcoff=icartl+natom3*nfnamem                                      2d3s25
      call enough('prebodcfd.cartl',bc,ibc)                             2d3s25
c
c     orthogonalize to determine dispacements
c
      nuniq=0
      do i=0,nfnamem-1                                                  2d3s25
       jcartl=icartl+natom3*nuniq                                       2d3s25
       jcartp=icartp+natom3*i                                           2d3s25
       do k=0,natom3-1                                                  2d3s25
        bc(jcartl+k)=bc(jcartp+k)                                       2d3s25
       end do                                                           2d3s25
       do j=0,nuniq-1                                                   2d3s25
        jjcartl=icartl+natom3*j                                         2d3s25
        dot=0d0                                                         2d3s25
        do k=0,natom3-1                                                 2d3s25
         dot=dot+bc(jcartl+k)*bc(jjcartl+k)                             2d3s25
        end do                                                          2d3s25
        do k=0,natom3-1                                                 2d3s25
         bc(jcartl+k)=bc(jcartl+k)-dot*bc(jjcartl+k)                    2d3s25
        end do                                                          2d3s25
       end do                                                           2d3s25
       sz=0d0                                                           2d3s25
       do k=0,natom3-1                                                  2d3s25
        sz=sz+bc(jcartl+k)**2                                           2d3s25
       end do                                                           2d3s25
       write(6,*)('sz = '),sz
       if(sz.gt.1d-8)then
        write(6,*)('unique! ')
        sz=1d0/sqrt(sz)                                                 2d3s25
        do k=0,natom3-1                                                 2d3s25
         bc(jcartl+k)=bc(jcartl+k)*sz                                   2d3s25
        end do                                                          2d3s25
        nuniq=nuniq+1                                                   2d3s25
       end if
      end do                                                            2d3s25
      write(6,*)('unique dispacement directions ')                      2d3s25
      call prntm2(bc(icartl),natom3,nuniq,natom3)                       2d3s25
      iawgt=ibcoff                                                      2d3s25
      ibcoff=iawgt+natom                                                2d3s25
      call enough('prebodcfd.awgt',bc,ibc)                              2d3s25
      do i=0,nuniq-1                                                    2d3s25
       write(6,*)('displacement direction '),i+1                        2d3s25
       do iz=iawgt,ibcoff-1                                             2d3s25
        bc(iz)=0d0                                                      2d3s25
       end do                                                           2d3s25
       ax=0d0                                                           2d3s25
       do ia=0,natom-1                                                  2d3s25
        jcartl=icartl+3*(ia+natom*i)                                    2d3s25
        do ixyz=0,2                                                     2d3s25
         bc(iawgt+ia)=bc(iawgt+iz)+bc(jcartl+ixyz)**2                   2d3s25
        end do                                                          2d3s25
        ax=max(ax,bc(iawgt+ia))                                         2d3s25
       end do                                                           2d3s25
       ax9=0.39d0*ax                                                    2d3s25
       do ia=0,natom-1                                                  2d3s25
        if(bc(iawgt+ia).gt.ax9)then                                     2d3s25
         write(6,*)('has contribution from atom '),ia+1                 2d3s25
        end if                                                          2d3s25
       end do                                                           2d3s25
      end do                                                            2d3s25
      ibcoff=iawgt                                                      2d3s25
      idatal=ibcoff                                                     2d3s25
      jdatal=idatal                                                     2d3s25
      do i=0,nuniq-1                                                    2d3s25
       jcartl=icartl+natom3*i                                           2d3s25
       write(6,*)('for dispacement no. '),i+1
       nstep=0                                                          2d3s25
       jdatal0=jdatal                                                   2d3s25
       jdatal=jdatal+1                                                  2d3s25
       do j=0,nfnamem-1                                                 2d3s25
        jcartp=icartp+natom3*j                                          2d3s25
        dot=0d0                                                         2d3s25
        do k=0,natom3-1                                                 2d3s25
         dot=dot+bc(jcartl+k)*bc(jcartp+k)                              2d3s25
        end do                                                          2d3s25
        if(abs(dot).gt.1d-8)then
         write(6,*)('filename '),j+1,('has stepsize '),dot
         bc(jdatal)=dot                                                 2d3s25
         jdatal=jdatal+1                                                2d3s25
         ibc(jdatal)=j+2                                                2d13s25
         jdatal=jdatal+1                                                2d3s25
         nstep=nstep+1                                                  2d3s25
        end if
       end do
       write(6,*)('no. of steps = '),nstep                              2d3s25
       nstepp=nstep+1                                                   2d3s25
       ifdstp=jdatal+2                                                  2d3s25
       ibcoff=ifdstp+nstepp                                             2d3s25
       call enough('prebodcfd.fdstp',bc,ibc)                            2d3s25
       jj=jdatal0+1                                                     2d3s25
       jfdstp=ifdstp-1                                                  2d3s25
       do j=1,nstep                                                     2d3s25
        bc(jfdstp+j)=bc(jj)                                             2d3s25
        jj=jj+2                                                         2d3s25
       end do                                                           2d3s25
       bc(ifdstp+nstep)=0d0                                             2d3s25
       call step(nstepp,bc(ifdstp),bc(ibcoff),1)                        2d3s25
       write(6,*)('fd coefficients: ')                                  2d3s25
       call prntm2(bc(ibcoff),nstepp,1,nstepp)                          2d3s25
       write(6,*)(' testing on '),bc(ibcoff+nstep)
       if(abs(bc(ibcoff+nstep)).lt.1d-10)then                           2d3s25
        nsteppu=nstep                                                   2d3s25
       else                                                             2d3s25
        nsteppu=nstepp                                                  2d3s25
       end if                                                           2d3s25
       jj=jdatal0+1                                                     2d3s25
       jbcoff=ibcoff-1                                                  2d3s25
       do j=1,nstepp                                                    2d3s25
        bc(jj)=bc(jbcoff+j)                                             2d3s25
        jj=jj+2                                                         2d3s25
       end do                                                           2d3s25
       if(abs(bc(ibcoff+nstep)).gt.1d-10)then                           2d3s25
        ibc(jj+1)=1                                                      2d3s25
       end if                                                           2d3s25
       ibc(jdatal0)=nsteppu                                             2d13s25
       write(6,*)('storing nsteppu '),nsteppu,('at '),jdatal0
      end do                                                            2d3s25
      ibcoff=jdatal                                                     2d3s25
      ndatal=ibcoff-idatal                                              2d3s25
      return                                                            2d3s25
      end                                                               2d3s25
