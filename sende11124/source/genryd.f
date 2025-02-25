c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine genryd(jmats,kmats,h0,iaorb,                           5d7s18
     $     nalpha,iborb,nbeta,nsbeta,nconf,iveca,nvirtc,m12sym,idata1,  5d7s18
     $     idatb1,eigcation,nroot,iorbn,morb,nbasisp,ihryd,istate,      1d18s20
     $     nstatex,lprint,nryd,ncomp,nlzz,iorbsym,iorbsymz,wwww,wsum,   10d13s22
     $     shift,bc,ibc,icanon)                                         5d5s23
      implicit real*8 (a-h,o-z)
      include "common.hf"
      include "common.store"
      include "common.cas"                                              5d7s18
      include "common.print"                                            5d6s20
      logical lwrite,lprint                                             1d18s20
      integer*1 iaorb(nalpha,*),iborb(nbeta,*)
      dimension jmats(*),kmats(*),h0(*),iveca(*),                       5d7s18
     $     nsbeta(*),noc(8),nvirtc(*),jtype(8),ktype(8),nrj(8),nrk(8),  5d7s18
     $     idata1(*),idatb1(*),m12sym(8,2,2),eigcation(nroot),          12d31s19
     $     iorbn(*),morb(*),nbasisp(*),ihryd(8),nryd(8),iorbsym(*),     5d3s21
     $     iorbsymz(*),wwww(*),icanon(*)                                5d5s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      data ecation/1d10/                                                5d12s23
      lwrite=lprint.and.istate.eq.nstatex                               1d18s20
      if(mynowprog.eq.0)then
      end if
      if(lwrite)write(6,*)                                              1d18s20
     $     ('in genryd for state '),istate,(' of '),nstatex             1d18s20
      do nsrydb=1,nsymb
       nryd(nsrydb)=0                                                   1d18s20
       nvrt=nbasdws(nsrydb)-idoub(nsrydb)-iacto(nsrydb)
       if(nvrt.gt.0)then                                                1d18s20
        if(lwrite)then                                                  1d18s20
         write(6,*)('for symmetry '),nsrydb                               1d18s20
         write(6,*)('idoub of this sym: '),idoub(nsrydb)
         write(6,*)('iacto of this sym: '),iacto(nsrydb)
         write(6,*)('nbasdws of this sym: '),nbasdws(nsrydb)
         write(6,*)('mxryd of this sym: '),mxryd(nsrydb)                9d10s21
         write(6,*)('energies will be given wrt N-1 core state')
         write(6,*)('size of ci matrix = '),nvrt
        end if                                                          1d18s20
        call ilimts(nvirtc(nsrydb),nvirtc(nsrydb),mynprocg,mynowprog,     5d7s18
     $     il,ih,i1s,i1e,i2s,i2e)                                       5d7s18
        nact=0                                                            5d7s18
        do isb=1,nsymb                                                    5d7s18
         nact=nact+iacto(isb)                                             5d7s18
         noc(isb)=idoub(isb)+iacto(isb)                                   5d7s18
         nrj(isb)=(noc(isb)*(noc(isb)+1))/2                               5d7s18
         nrk(isb)=noc(isb)*noc(isb)                                       5d7s18
         jtype(isb)=1                                                     5d7s18
         ktype(isb)=1                                                     5d7s18
         do i=1,nsdlk                                                     5d7s18
          if(isblk(1,i).eq.isb.and.isblk(2,i).eq.isb.and.                 5d7s18
     $       isblk(3,i).eq.nsrydb)then                                    5d7s18
           jtype(isb)=i                                                   5d7s18
           go to 1                                                        5d7s18
          end if                                                          5d7s18
         end do                                                           5d7s18
         if(nrj(isb).gt.0)write(6,*)('we have a miss for j'),isb,nsrydb, 1d18s20
     $       nrj(isb)                                                   1d18s20
    1    continue
         do i=1,nsdlkk                                                     5d7s18
          if(isblkk(1,i).eq.isb.and.isblkk(2,i).eq.isb.and.               5d7s18
     $       isblkk(3,i).eq.nsrydb)then                                 5d7s18
           ktype(isb)=i                                                   5d7s18
           go to 2                                                        5d7s18
          end if                                                          5d7s18
         end do                                                           5d7s18
         if(nrk(isb).gt.0)write(6,*)('we have a miss for k'),isb,nsrydb, 1d18s20
     $       nrk(isb)                                                   1d18s20
    2    continue
        end do                                                            5d7s18
        ism=ibcoff                                                        5d7s18
        irelo=ism+nact                                                    5d7s18
        ibcoff=irelo+nact                                                 5d7s18
        call enough('genryd.  2',bc,ibc)
        jsm=ism                                                           5d7s18
        jrelo=irelo                                                       5d7s18
        do isb=1,nsymb                                                    5d7s18
         do i=1,iacto(isb)                                                5d7s18
          ibc(jsm)=isb                                                    5d7s18
          ibc(jrelo)=i+idoub(isb)                                         5d7s18
          jsm=jsm+1                                                       5d7s18
          jrelo=jrelo+1                                                   5d7s18
         end do                                                           5d7s18
        end do                                                            5d7s18
        jrelo=irelo-1                                                     5d7s18
        jsm=ism-1                                                         5d7s18
c
c     fold doubly occupied orbital contribution into h0.
c
        ioffh=1                                                           5d7s18
        do isb=1,nsrydb-1                                                 5d7s18
         nb=noc(isb)+nvirtc(isb)                                          5d7s18
         ioffh=ioffh+nb*nb                                                5d7s18
        end do                                                            5d7s18
        ih0f=ibcoff                                                       5d7s18
        nb=noc(nsrydb)+nvirtc(nsrydb)                                     5d7s18
        ibcoff=ih0f+nvirtc(nsrydb)*nvirtc(nsrydb)                         5d7s18
        call enough('genryd.  3',bc,ibc)
        if(mynowprog.eq.0)then                                            5d7s18
         do i=0,nvirtc(nsrydb)-1                                          5d7s18
          iad1=ih0f+nvirtc(nsrydb)*i                                      5d7s18
          iad2=ioffh+noc(nsrydb)+nb*(i+noc(nsrydb))                       5d7s18
          do j=0,nvirtc(nsrydb)-1                                         5d7s18
           bc(iad1+j)=h0(iad2+j)                                          5d7s18
          end do                                                          5d7s18
         end do                                                           5d7s18
        else                                                              5d7s18
         do i=0,nvirtc(nsrydb)*nvirtc(nsrydb)-1                           5d7s18
          bc(ih0f+i)=0d0                                                  5d7s18
         end do                                                           5d7s18
        end if                                                            5d7s18
        if(iprtr(19).ne.0)then
         write(6,*)('for symmetry block '),nsrydb                          5d7s18
         write(6,*)('vv part of h0 ')                                      5d7s18
         call prntm2(bc(ih0f),nvirtc(nsrydb),nvirtc(nsrydb),            5d6s20
     $        nvirtc(nsrydb))                                           5d6s20
        end if
        do isd=1,nsymb                                                    5d7s18
         if(idoub(isd).gt.0)then                                          5d7s18
          i10=i1s                                                          5d7s18
          i1n=nvirtc(nsrydb)                                               5d7s18
          jj=jmats(jtype(isd))-1                                           5d7s18
          np=noc(isd)+1                                                    5d7s18
          kk=kmats(ktype(isd))-np                                          5d7s18
          do i2=i2s,i2e                                                    5d7s18
           i2m=i2-1                                                        5d7s18
           if(i2.eq.i2e)i1n=i1e                                            5d7s18
           do i1=i10,i1n                                                   5d7s18
            i1m=i1-1                                                       5d7s18
            iad1=ih0f+i1m+nvirtc(nsrydb)*i2m                               5d7s18
            do i34=1,idoub(isd)                                            5d7s18
             iad2=jj+(i34*(i34+1))/2                                       5d7s18
             iad3=kk+i34*np                                                5d7s18
             orig=bc(iad1)
             bc(iad1)=bc(iad1)+2d0*bc(iad2)-bc(iad3)                      5d8s18
            end do                                                         5d7s18
            jj=jj+nrj(isd)                                                 5d7s18
            kk=kk+nrk(isd)                                                 5d7s18
           end do                                                          5d7s18
           i10=1                                                           5d7s18
          end do                                                           5d7s18
         end if                                                           5d7s18
        end do                                                            5d7s18
        if(iprtr(19).ne.0)then                                          5d6s20
         write(6,*)('vv part of h0 with doubles part folded in')
         call prntm2(bc(ih0f),nvirtc(nsrydb),nvirtc(nsrydb),            5d6s20
     $        nvirtc(nsrydb))                                           5d6s20
         nnvc=nvirtc(nsrydb)*nvirtc(nsrydb)                             5d6s21
         ih0fgs=ibcoff                                                  5d6s21
         ibcoff=ih0fgs+nnvc                                             5d6s21
         call enough('genryd.  4',bc,ibc)
         do ii=0,nnvc-1                                                 5d6s21
          bc(ih0fgs+ii)=bc(ih0f+ii)                                     5d6s21
         end do                                                         5d6s21
         call dws_gsumf(bc(ih0fgs),nnvc)                                5d6s21
         write(6,*)('global summed')                                      5d7s18
         call prntm2(bc(ih0fgs),nvirtc(nsrydb),nvirtc(nsrydb),            5d6s20
     $        nvirtc(nsrydb))                                           5d6s20
         ibcoff=ih0fgs                                                  5d6s21
        end if                                                          5d6s20
c
c     diagonal correction
c
        nn=(nvrt*(nvrt+1))/2                                              5d7s18
        iham=ibcoff                                                       5d7s18
        ibcoff=iham+nn                                                    5d7s18
        call enough('genryd.  5',bc,ibc)
        ibctop=ibcoff                                                     5d10s18
        wsumc=0d0                                                       1d31s22
        do iroot=1,nroot                                                  5d10s18
         ecation=min(ecation,eigcation(iroot))                          5d12s23
         wsumc=wsumc+wwww(iroot)                                        1d31s22
         irotm=iroot-1                                                    5d10s18
         do i=0,nn-1                                                       5d7s18
          bc(iham+i)=0d0                                                   5d7s18
         end do                                                            5d7s18
         ioff=0                                                            5d7s18
         do isb=1,nsymb                                                    5d7s18
          jsb=nsbeta(isb)                                                  5d7s18
          if(min(nadet(isb),nbdet(jsb)).gt.0)then                        1d18s20
           if(iprtr(19).ne.0)write(6,*)('for symmetry block '),isb,jsb  5d7s21
           ioffa=0                                                          5d7s18
           do i=1,isb-1                                                     5d7s18
            ioffa=ioffa+nadet(i)                                            5d7s18
           end do                                                           5d7s18
           ioffb=0                                                          5d7s18
           do i=1,jsb-1                                                     5d7s18
            ioffb=ioffb+nbdet(i)                                            5d7s18
           end do                                                           5d7s18
           itmpa=ibcoff                                                     5d7s18
           itmpb=itmpa+nadet(isb)*nn                                        5d7s18
           ibcoff=itmpb+nbdet(jsb)*nn                                       5d7s18
           call enough('genryd.  6',bc,ibc)
           do i=0,nn*(nadet(isb)+nbdet(jsb))-1                              5d7s18
            bc(itmpa+i)=0d0                                                 5d7s18
           end do                                                           5d7s18
           do ia=1,nadet(isb)                                               5d7s18
            iap=ia+ioffa                                                    5d7s18
            iad=itmpa+nn*(ia-1)                                             5d8s18
            do ie=1,nalpha                                                  5d7s18
             is=ibc(jsm+iaorb(ie,iap))                                      5d7s18
             io=ibc(jrelo+iaorb(ie,iap))                                    5d7s18
             ioj=((io*(io+1))/2)-1                                          5d7s18
             np=noc(is)+1                                                   5d7s18
             iok=(io-1)*np                                                  5d8s18
             do iva=0,nvrt-1                                                5d7s18
              do ivb=0,iva                                                  5d7s18
               icol=1+ivb+nvirtc(nsrydb)*iva                                5d7s18
               iiv=((iva*(iva+1))/2)+ivb                                     5d7s18
               if(icol.ge.il.and.icol.le.ih)then                             5d7s18
                iadj=jmats(jtype(is))+ioj+nrj(is)*(icol-il)                  5d7s18
                iadk=kmats(ktype(is))+iok+nrk(is)*(icol-il)                  5d7s18
                bc(iad+iiv)=bc(iad+iiv)+bc(iadj)-bc(iadk)                    5d7s18
               end if                                                        5d7s18
              end do                                                        5d7s18
             end do                                                         5d7s18
            end do                                                          5d7s18
           end do                                                           5d7s18
           do ib=1,nbdet(jsb)                                               5d7s18
            ibp=ib+ioffb                                                    5d7s18
            iad=itmpb+nn*(ib-1)                                             5d8s18
            do ie=1,nbeta                                                   5d7s18
             is=ibc(jsm+iborb(ie,ibp))                                      5d7s18
             io=ibc(jrelo+iborb(ie,ibp))                                    5d7s18
             ioj=((io*(io+1))/2)-1                                          5d7s18
             do iva=0,nvrt-1                                                5d7s18
              do ivb=0,iva                                                  5d7s18
               ivv=((iva*(iva+1))/2)+ivb                                    5d7s18
               icol=1+ivb+nvirtc(nsrydb)*iva                                5d7s18
               if(icol.ge.il.and.icol.le.ih)then                             5d7s18
                iadj=jmats(jtype(is))+ioj+nrj(is)*(icol-il)                  5d7s18
                bc(iad+ivv)=bc(iad+ivv)+bc(iadj)                            5d8s18
               end if                                                        5d7s18
              end do                                                        5d7s18
             end do                                                         5d7s18
            end do                                                          5d7s18
           end do                                                           5d7s18
           if(iprtr(19).ne.0)then                                       5d7s21
            write(6,*)('my part of tmpa ')
            call prntm2(bc(itmpa),nn,nadet(isb),nn)
            write(6,*)('my part of tmpb ')
            call prntm2(bc(itmpb),nn,nbdet(jsb),nn)
            nntotnn=nn*(nadet(isb)+nbdet(jsb))
            do i=0,nntotnn-1
             bc(ibcoff+i)=bc(itmpa+i)
            end do
            call dws_gsumf(bc(ibcoff),nntotnn)
            write(6,*)('global summed: ')
            call prntm2(bc(ibcoff),nn,nadet(isb),nn)
            jtmpbxx=ibcoff+nn*nadet(isb)
            call prntm2(bc(jtmpbxx),nn,nbdet(jsb),nn)
           end if                                                       5d7s21
           do ib=0,nbdet(jsb)-1                                             5d7s18
            iadb=itmpb+nn*ib                                                5d7s18
            do ia=0,nadet(isb)-1                                            5d7s18
             iada=itmpa+nn*ia                                               5d7s18
             iadv=iveca(isb)+irotm+nroot*(ia+nadet(isb)*ib)                5d10s18
             ww=bc(iadv)**2                                                 5d7s18
             ww=ww*wwww(iroot)                                          1d31s22
             do iv=0,nn-1                                                   5d7s18
              bc(iham+iv)=bc(iham+iv)+ww*(bc(iada+iv)+bc(iadb+iv))          5d8s18
             end do                                                         5d7s18
            end do                                                          5d7s18
           end do                                                           5d7s18
           ibcoff=itmpa                                                     5d7s18
          end if                                                           1d18s20
         end do                                                            5d7s18
         if(iprtr(19).ne.0)then                                         5d6s20
          write(6,*)('my part of ham b4 h0f')
          do i=0,nn-1
           bc(ibcoff+i)=bc(iham+i)
          end do
          call dws_gsumf(bc(ibcoff),nn)
          write(6,*)('global summed ')
          call mpprnt2(bc(ibcoff),nvrt)
          call prntm2(bc(ibcoff),nn,1,nn)
         end if                                                         5d6s20
         jham=iham                                                         5d7s18
         do iva=0,nvrt-1                                                   5d7s18
          do ivb=0,iva                                                     5d7s18
           iadh=ih0f+ivb+iva*nvirtc(nsrydb)                                5d7s18
           bc(jham)=bc(jham)+bc(iadh)                                      5d7s18
     $          *wwww(iroot)                                            1d31s22
           jham=jham+1                                                     5d7s18
          end do                                                           5d7s18
         end do                                                            5d7s18
         if(iprtr(19).ne.0)then                                         5d6s20
          write(6,*)('my part of ham after diagonal N-1 fcns:')             5d7s18
          call mpprnt2(bc(iham),nvrt)
          do i=0,nn-1
           bc(ibcoff+i)=bc(iham+i)
          end do
          call dws_gsumf(bc(ibcoff),nn)
          write(6,*)('global summed ')
          call mpprnt2(bc(ibcoff),nvrt)
         end if                                                         5d6s20
c
c     single differences of N-1 fcns
c
         do isb=1,nsymb                                                    5d7s18
          jsb=nsbeta(isb)                                                  5d7s18
          call singnm1(bc(iham),nvrt,bc(iveca(isb)),nadet(isb),          1d18s20
     $        nbdet(jsb),                                               1d18s20
     $      m12sym(isb,1,1),ibc(idata1(isb)),m12sym(jsb,1,2),           5d7s18
     $      ibc(idatb1(jsb)),ibc(ism),ibc(irelo),jmats,jtype,kmats,     5d7s18
     $      ktype,il,ih,nvirtc(nsrydb),nrj,nrk,noc,iroot,nroot,wwww,bc, 11d10s22
     $         ibc)                                                     11d10s22
         end do                                                            5d7s18
         if(iprtr(19).ne.0)then                                         5d6s20
          write(6,*)('my part of ham after single dif N-1 fcns:')             5d7s18
          call mpprnt2(bc(iham),nvrt)
         end if                                                         5d6s20
         do i=0,nn-1                                                     1d18s20
          bc(ihryd(nsrydb)+i)=bc(ihryd(nsrydb)+i)+bc(iham+i)             1d18s20
         end do                                                          1d18s20
         if(istate.eq.nstatex.and.iroot.eq.nroot)then                   5d12s23
          fff=1d0/wsum                                                  1d31s22
          do i=0,nn-1                                                    1d18s20
           bc(iham+i)=bc(ihryd(nsrydb)+i)*fff                            1d18s20
          end do                                                         1d18s20
          call dws_gsumf(bc(iham),nn)
          if(iprtr(19).ne.0)then                                        5d6s20
           write(6,*)('global summed ...')
           call mpprnt2(bc(iham),nvrt)
          end if                                                        5d6s20
          ieig=ibcoff                                                   5d3s21
          ivec=ieig+nvrt                                                5d3s21
          isymv=ivec+nvrt*nvrt                                          5d3s21
          ibcoff=isymv+nvrt                                             5d3s21
          call enough('genryd.  7',bc,ibc)
          ihsqr=ibcoff                                                  5d3s21
          ibcoff=ihsqr+nvrt*nvrt                                        5d3s21
          call enough('genryd.  8',bc,ibc)
          do i=0,nn-1                                                   5d3s21
           bc(ihsqr+i)=bc(iham+i)                                       5d3s21
          end do                                                        5d3s21
          call square(bc(ihsqr),nvrt)                                   5d3s21
          if(nlzz.ne.0)then                                             5d3s21
           jorbs=iorbsym(nsrydb)+idoub(nsrydb)+iacto(nsrydb)            5d3s21
           if(nlzz.eq.2)then                                            5d3s21
            do i=0,nvrt-1                                               5d3s21
             ibc(isymv+i)=ibc(jorbs+i)                                  5d3s21
            end do                                                      5d3s21
           else                                                         5d3s21
            jorbz=iorbsymz(nsrydb)+idoub(nsrydb)+iacto(nsrydb)          5d3s21
            do i=0,nvrt-1                                               5d3s21
             ibc(isymv+i)=ibc(jorbs+i)+100*ibc(jorbz+i)                 5d3s21
            end do                                                      5d3s21
           end if                                                       5d3s21
           call diagy(nvrt,bc(ihsqr),bc(ieig),bc(ivec),ibc(isymv),bc,   11d14s22
     $          ibc,0,idum,dum)                                         9d1s23
           if(nlzz.eq.2)then                                            5d3s21
            do i=0,nvrt-1                                               5d3s21
             ibc(jorbs+i)=ibc(isymv+i)                                  5d3s21
            end do                                                      5d3s21
           else                                                         5d3s21
            do i=0,nvrt-1                                               5d3s21
             ibc(jorbz+i)=ibc(isymv+i)/100                              5d3s21
             ibc(jorbs+i)=ibc(isymv+i)-100*ibc(jorbz+i)                 5d3s21
            end do                                                      5d3s21
           end if                                                       5d3s21
          else                                                          5d3s21
           call diagx(nvrt,bc(ihsqr),bc(ieig),bc(ivec),ibc(isymv),bc,   11d14s22
     $          ibc)                                                    11d14s22
          end if                                                        5d3s21
          ibcoff=ihsqr                                                  5d3s21
          nrydb=0                                                       1d18s20
          do i=0,nvrt-1                                                   1d18s20
           if(bc(ieig+i).lt.0d0)then                                      1d18s20
            if(lwrite)then                                              5d3s21
             if(nlzz.eq.0)then                                          5d3s21
              if(i.eq.0)write(6,*)('High-spin Rydberg eigenvalues: ')    5d3s21
              write(6,*)i,bc(ieig+i)                                    5d3s21
             else                                                       5d3s21
              if(nlzz.eq.2)then                                         5d3s21
               if(i.eq.0)
     $             write(6,*)('High-spin Rydberg eigenvalues with Lz: ')5d3s21
               write(6,*)i,bc(ieig+i),ibc(isymv+i)                      5d3s21
              else                                                      5d3s21
               lz=ibc(isymv+i)/100                                      5d3s21
               ll=ibc(isymv+i)-100*lz                                   5d3s21
               if(i.eq.0)write(6,*)                                     5d3s21
     $              ('High-spin Rydberg eigenvalues with L,Lz: ')       5d3s21
               write(6,*)i,bc(ieig+i),ll,lz                             5d3s21
              end if                                                    5d3s21
             end if                                                     5d3s21
            end if                                                      5d3s21
            nrydb=i+1                                                     1d18s20
           end if                                                         1d18s20
          end do                                                          1d18s20
          if(mxryd(nsrydb).gt.0)then                                    9d10s21
           nryd(nsrydb)=mxryd(nsrydb)                                   9d10s21
          else                                                          9d10s21
           nryd(nsrydb)=nrydb                                            1d18s20
          end if                                                        9d10s21
          if(lwrite)then                                                1d18s20
           ecations=ecation+shift                                       6d7s23
           if(nrydb.gt.0)write(6,*)('with lowest cation energy: ')      6d7s23
           do i=0,nrydb-1
            bc(ieig+i)=bc(ieig+i)+ecations                              6d7s23
            write(6,*)i,bc(ieig+i)
           end do
          end if                                                        1d18s20
          iovt=iorbn(nsrydb)                                             1d18s20
     $        +ncomp*nbasisp(nsrydb)*(idoub(nsrydb)+iacto(nsrydb))      3d3s20
          itmp=ibcoff                                                      12d31s19
          ibcoff=itmp+ncomp*nbasisp(nsrydb)*nvrt                        3d3s20
          call enough('genryd.  9',bc,ibc)
          call dgemm('n','n', ncomp*nbasisp(nsrydb),nvrt,nvrt,1d0,      3d3s20
     $         bc(iovt),nbasisp(nsrydb)*ncomp,bc(ivec),nvrt,0d0,        3d3s20
     $         bc(itmp),nbasisp(nsrydb)*ncomp,                          3d3s20
     d' genryd.  1')
          do i=0,nbasisp(nsrydb)*ncomp*nvrt-1                           3d3s20
           bc(iovt+i)=bc(itmp+i)                                           12d31s19
          end do                                                           12d31s19
          ibcoff=ibctop                                                    5d10s18
         end if                                                          1d18s20
c
c     end of loop over roots!
        end do                                                            5d10s18
        ibcoff=ism
       end if                                                           1d18s20
      end do                                                            1d18s20
      return
      end
