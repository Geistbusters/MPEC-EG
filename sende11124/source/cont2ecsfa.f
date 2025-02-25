c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine cont2ecsfa(ibasis,iptr,idorb,isorb,nfcn,nec,nrootu,    9d10s19
     $     icsfpd,vint,ndim,mdon,                                       11d19s20
     $     ncsf,ncsf2,isbc,isbr,vc,nct2,ndimv,ivcc,nct2c,multh,irooto,  9d28s20
     $     fact3jst,fact3jht,fact3jlt,fact3jst2,fact3jht2,srhh,ilcol,   10d17s20
     $     ihcol,iff2,iptrto,jptrto,nff2,mdoo,nsymb,ndetc,icsfxv,nrtot, 11d14s22
     $     bc,ibc)                                                      11d14s22
      implicit real*8 (a-h,o-z)
      integer*1 idorb(*),isorb(*),iorb(64)                              11d19s20
      integer*8 iff2(*),iptrto(*),itesta,itestb,icsfxv,i18,i28,i38,i48  2d1s21
      logical ldebug                                                    12d12s19
      dimension ibasis(3,*),icsfpd(*),iptr(4,*),                        11d19s20
     $     vint(ndim,*),ncsf(*),multh(8,8),iotmp(3),                    11d19s20
     $     vc(nct2,ndimv,*),ntmp(2),ncsf2(4,*),ivcc(4),nct2c(4),        8d14s19
     $     joffg(4),inst(2),ncsf2x(4),nff2(mdoo+1,nsymb,*),ndetc(*)     1d12s23
      include "common.store"
      include "common.mrci"                                             9d10s19
      include "common.print"                                            1d5s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      if(iprtr(9).eq.0)then                                             1d5s20
       ldebug=.false.                                                    12d12s19
      else                                                              1d5s20
       ldebug=.true.                                                    1d5s20
      end if                                                            1d5s20
      ismultm=ismult-1                                                  9d16s19
      if(ldebug)then                                                    12d12s19
       write(6,*)('hi, my name is cont2ecsfa')
       write(6,*)('starting jptrto: '),jptrto
      end if                                                            12d12s19
      srh=1d0
      ndimh=nct2*ndimv                                                  8d1s19
      ioff=1                                                            7d26s19
      do if=1,nfcn                                                      7d26s19
       nclo=ibasis(1,if)
       nclom=nclo-1                                                     9d10s19
       nclomm=nclom-1                                                   9d10s19
       nclop=nclo+1
       nopen=nec-2*nclo
       nopenpp=nopen+2                                                  9d10s19
       nopenmm=nopen-2                                                  9d10s19
       ic=iptr(2,nclop)+nclo*(ibasis(2,if)-1)
       io=iptr(4,nclop)+nopen*(ibasis(3,if)-1)
       iarg=nclop-mdon
       if(ldebug)then                                                   12d12s19
        write(6,*)('for function ')
        write(6,*)if,(idorb(ic+j),j=0,nclo-1),(':'),
     $      (isorb(io+j),j=0,nopen-1)
        write(6,*)('vectors in csf basis: ')
        call prntm2(vint(ioff,1),ncsf(iarg),nrootu,ndim)                 9d10s19
       end if                                                           12d12s19
       call fetchcsf2(nopen,ismultm,nrooto,vecadd,ibc(ibcoff),
     $      ibc(ibcoff),ismultm,ndet,bc,ibc)                            11d15s22
       ias=ibcoff                                                       12d7s20
       ibcoff=ias+ndet*nopen                                            12d7s20
       call enough('cont2ecsfa.  1',bc,ibc)
       call fetchcsf2(nopen,ismultm,nrooto,vecadd,ibc(ias),             12d7s20
     $      ibc(ias),ismultm,ndet,bc,ibc)                               11d15s22
       ivecp=nint(vecadd)                                               12d7s20
       if(ivecp.lt.0.or.ivecp.gt.ibcoff)then                            12d1s20
        write(6,*)('bad vectors address: '),ivecp                       12d1s20
        call dws_sync                                                   12d1s20
        call dws_finalize                                               12d1s20
        stop                                                            12d1s20
       end if                                                           12d1s20
       nn=ncsf(iarg)*ncsf(iarg)
       ivnorm=ibcoff                                                    1d22s20
       ivdet=ibcoff                                                     9d10s19
       ibcoff=ivdet+max(1,ndet)*nrootu                                  11d4s19
       call enough('cont2ecsfa.  2',bc,ibc)
       if(ldebug.and.mdon.ne.nclo)then                                  2d1s21
        write(6,*)('transformation to use: '),ivecp
        call prntm2(bc(ivecp),ndet,ncsf(iarg),ndet)
       end if                                                           12d12s19
       if(ndet.eq.0)then                                                11d4s19
        do i=1,nrootu                                                   11d4s19
         im=i-1                                                         11d4s19
         bc(ivdet+im)=vint(ioff,i)                                      11d4s19
        end do                                                          11d4s19
       else                                                             11d4s19
        if(nclo.ne.mdon.or.nopen.eq.2)then                              6d6s21
         call dgemm('n','n',ndet,nrootu,ncsf(iarg),1d0,bc(ivecp),ndet,    9d10s19
     $      vint(ioff,1),ndim,0d0,bc(ivdet),ndet,                       9d10s19
     d' cont2ecsfa.  1')
        else                                                            2d1s21
         i18=1                                                          2d1s21
         i28=ncsf(iarg)                                                 2d1s21
         do iproc=0,mynprocg-1                                          2d1s21
          call ilimts(1,ndet,mynprocg,iproc,il,ih,i1s,i1e,i2s,i2e)      3d30s21
          nhere=ih+1-il                                                 2d1s21
          if(nhere.gt.0)then                                            4d14s21
           itmpt=ibcoff                                                  2d1s21
           ibcoff=itmpt+nhere*ncsf(iarg)                                 2d1s21
           call enough('cont2ecsfa.  3',bc,ibc)
           i38=il                                                        2d1s21
           i48=ih                                                        2d1s21
           if(icsfxv.ge.0)then                                          8d22s24
            call ddi_get(bc,ibc,icsfxv,i18,i28,i38,i48,bc(itmpt))        11d15s22
            jvdet=ivdet+il-1                                              2d1s21
            call dgemm('n','n',nhere,nrootu,ncsf(iarg),1d0,bc(itmpt),
     $         nhere,vint(ioff,1),ndim,0d0,bc(jvdet),ndet,              2d1s21
     d' cont2ecsfa.  2')
           else                                                         7d19s24
            bc(ivdet)=vint(ioff,1)                                      7d19s24
           end if                                                       7d19s24
           ibcoff=itmpt                                                  2d1s21
          end if                                                        4d14s21
         end do                                                         2d1s21
        end if                                                          2d1s21
       end if                                                           11d4s19
       if(ldebug)then                                                   12d12s19
        write(6,*)('in det basis ')                                      9d10s19
        call prntm2(bc(ivdet),max(1,ndet),nrootu,max(1,ndet))            11d4s19
       end if                                                           12d12s19
c
c     there are 4 scenarios:
c     1: excite both electrons from a single cs orbital
c     2: excite electrons from two distinct cs orbitals
c     3: excite one electron from a cs orbital and one from an os orbital
c     4: excite electrons from two distinct os orbitals.
c
       if(nclo.gt.0)then                                                11d16s20
        if(isbc.eq.isbr)then                                            11d16s20
         do i=1,nclo                                                     9d10s19
          im=i-1                                                         9d10s19
          jj=0                                                           9d10s19
          itesta=0                                                      11d16s20
          do j=1,nclo                                                    9d10s19
           if(j.ne.i)then                                                9d10s19
            jm=j-1                                                       9d10s19
            jj=jj+1                                                       9d10s19
            iorb(jj)=idorb(ic+jm)                                         9d10s19
            itesta=ibset(itesta,idorb(ic+jm))                           11d16s20
           end if                                                        9d10s19
          end do                                                         9d10s19
          do j=0,nopen-1                                                11d16s20
           itesta=ibset(itesta,isorb(io+j)+32)                          11d18s20
          end do                                                        11d16s20
          icol=((idorb(ic+im)*(idorb(ic+im)+1))/2)-1                     9d12s19
          nfcnrel=iptrto(jptrto)-nff2(nclo,isbc,2)                      11d16s20
          if(icol.ge.ilcol.and.icol.le.ihcol)then
           if(ldebug)then                                                 12d12s19
            write(6,*)('scenario 1')
            write(6,*)('icol = '),icol
           end if                                                         12d12s19
           icol=icol-ilcol                                               10d17s20
           if(ldebug)write(6,*)('store to col '),icol
          jarg=nclo-mdon
          ndetj=ndetc(jarg+1)                                           11d19s20
          do l=1,4
           lll=nff2(nclo,isbc,2+l)+ncsf2(l,jarg)*nfcnrel
           joffg(l)=ivcc(l)+lll                                         11d19s20
          end do
c
c     get neg sign since my convention on doubly occupied orbs is
c     beta first.
c
           do j=1,nrootu                                                  9d12s19
            iadd=joffg(1)+nct2c(1)*(j-1+irooto+nrtot*icol)              4d28s21
            do k=0,ncsf(iarg)-1                                           9d12s19
             bc(iadd+k)=bc(iadd+k)-vint(ioff+k,j)                         9d16s19
            end do                                                        9d12s19
           end do                                                         9d12s19
          end if
          jptrto=jptrto+1                                                11d16s20
         end do                                                          9d10s19
        else                                                            11d16s20
         jptrto=jptrto+nclo                                             11d16s20
        end if                                                          11d16s20
       end if                                                           9d10s19
       if(nclo.gt.1)then                                                9d10s19
        do i=1,nclo                                                     9d10s19
         im=i-1                                                         9d10s19
         do j=1,i-1                                                     9d10s19
          jm=j-1                                                        9d10s19
          isji=multh(ism(idorb(ic+im)),ism(idorb(ic+jm)))               9d10s19
          if(multh(isji,isbc).eq.isbr)then                              9d10s19
           if(ldebug)then                                               12d12s19
            write(6,*)('scenario 2')
            write(6,*)('passed symmetry test ')
           end if                                                       12d12s19
           if(idorb(ic+im).lt.idorb(ic+jm))then                          9d10s19
            in=idorb(ic+im)                                             9d10s19
            ix=idorb(ic+jm)                                             9d10s19
           else                                                         9d10s19
            in=idorb(ic+jm)                                             9d10s19
            ix=idorb(ic+im)                                             9d10s19
           end if                                                       9d10s19
           if(ldebug)write(6,*)('in,ix: '),in,ix                        12d12s19
           icol=((ix*(ix-1))/2)+in-1                                    9d12s19
           if(icol.ge.ilcol.and.icol.le.ihcol)then                      10d17s20
            if(ldebug)write(6,*)('raw col '),icol
            icol=icol-ilcol                                             10d17s20
            if(ldebug)write(6,*)('store to col '),icol
            jj=0                                                          9d10s19
            itesta=0                                                    11d16s20
            do k=1,nclo                                                   9d10s19
             if(k.ne.i.and.k.ne.j)then                                    9d10s19
              km=k-1                                                     9d10s19
              jj=jj+1
              iorb(jj)=idorb(ic+km)                                      9d10s19
              itesta=ibset(itesta,idorb(ic+km))                         11d16s20
             end if                                                      9d10s19
            end do                                                        9d10s19
            iput=0                                                       9d10s19
            jjo=jj+1                                                     9d10s19
            do k=0,nopen-1                                               9d10s19
             if(isorb(io+k).lt.in)then                                   9d10s19
              jj=jj+1                                                    9d10s19
              iorb(jj)=isorb(io+k)                                       9d10s19
              itesta=ibset(itesta,isorb(io+k)+32)                       11d18s20
             else if(isorb(io+k).lt.ix)then                              9d10s19
              if(iput.eq.0)then                                          9d10s19
               jj=jj+1                                                   9d10s19
               iorb(jj)=in                                               9d10s19
               itesta=ibset(itesta,in+32)                               11d18s20
               inst(1)=k+1                                               9d10s19
               iput=1                                                    9d10s19
              end if                                                     9d10s19
              jj=jj+1                                                    9d10s19
              iorb(jj)=isorb(io+k)                                       9d10s19
             else if(isorb(io+k).gt.ix)then                              9d10s19
              if(iput.eq.0)then                                          9d10s19
               jj=jj+1                                                   9d10s19
               iorb(jj)=in                                               9d10s19
               itesta=ibset(itesta,in+32)                               11d18s20
               inst(1)=k+1                                               9d10s19
               iput=1                                                    9d10s19
              end if                                                     9d10s19
              if(iput.eq.1)then                                          9d10s19
               jj=jj+1                                                   9d10s19
               iorb(jj)=ix                                               9d10s19
               itesta=ibset(itesta,ix+32)                               11d18s20
               inst(2)=k+1                                               9d10s19
               iput=2                                                    9d10s19
              end if                                                     9d10s19
              jj=jj+1                                                    9d10s19
              iorb(jj)=isorb(io+k)                                       9d10s19
              itesta=ibset(itesta,isorb(io+k)+32)                       11d18s20
             end if                                                      9d10s19
            end do                                                       9d10s19
            if(iput.eq.0)then                                            9d10s19
             jj=jj+1                                                     9d10s19
             iorb(jj)=in                                                 9d10s19
             itesta=ibset(itesta,in+32)                                 11d18s20
             inst(1)=nopen+1                                             9d10s19
             iput=1                                                      9d10s19
            end if                                                       9d10s19
            if(iput.eq.1)then                                            9d10s19
             jj=jj+1                                                     9d10s19
             iorb(jj)=ix                                                 9d10s19
             inst(2)=nopen+1                                             9d10s19
             itesta=ibset(itesta,ix+32)                                 11d18s20
            end if                                                       9d10s19
            if(ldebug)write(6,*)('looking for '),(iorb(ij),ij=1,nclomm), 12d12s19
     $          (':'),(iorb(ij),ij=jjo,jj)                              12d12s19
            nfcnrel=iptrto(jptrto)-nff2(nclom,isbc,2)                     11d16s20
            jarg=nclom-mdon
            ndetj=ndetc(jarg+1)                                           11d19s20
            do l=1,4
             lll=nff2(nclom,isbc,2+l)+ncsf2(l,jarg)*nfcnrel
             joffg(l)=ivcc(l)+lll                                       11d19s20
            end do
            ivsp=ibcoff                                                  9d12s19
            ivtp=ivsp+ndetj*nrootu                                       9d12s19
            ibcoff=ivtp+ndetj*nrootu                                     9d12s19
            call enough('cont2ecsfa.  4',bc,ibc)
            iaorb=ibcoff                                                 9d16s19
            ivtmp=iaorb+ndetj*nopenpp                                    3d13s20
            ibcoff=ivtmp+ndetj*ncsf2(1,jarg)                                  9d16s19
            call enough('cont2ecsfa.  5',bc,ibc)
            call fetchcsf2(nopenpp,ismultm,nfcnx,bc(ivtmp),idum,         9d16s19
     $          ibc(iaorb),ismultm,idum,bc,ibc)                         11d15s22
            ibcoff=ivtmp                                                 9d16s19
            call matchc(inst,ibc(ias),ndet,nopen,2,ndetj,ibc(iaorb),     9d16s19
     $          nopenpp,bc(ivdet),bc(ivsp),bc(ivtp),nrootu,ismult,      9d12s19
     $          ncsf2(1,jarg),joffg,nct2c,icol,ndimv,0,vint(ioff,1),    11d4s19
     $          ndim,lcrenorm,bc(ivnorm),irooto,fact3jst,fact3jht,      9d28s20
     $          fact3jlt,fact3jst2,fact3jht2,srhh,nrtot,bc,ibc)         11d14s22
            ibcoff=ivsp                                                  9d16s19
           end if                                                       10d17s20
          end if                                                        9d10s19
          jptrto=jptrto+1                                                11d16s20
         end do                                                         9d10s19
        end do                                                          9d10s19
       end if                                                           9d10s19
       if(nclo.gt.0.and.nopen.gt.0)then                                 9d10s19
        do i=0,nclo-1                                                   9d10s19
         do j=0,nopen-1                                                 9d10s19
          if(ldebug)                                                    12d12s19
     $        write(6,*)('try exciting out of(2) '),
     $         idorb(ic+i),isorb(io+j)
          inst(1)=j+1                                                   9d10s19
          isij=multh(ism(idorb(ic+i)),ism(isorb(io+j)))
          if(multh(isij,isbc).eq.isbr)then                              9d10s19
           if(ldebug)write(6,*)('scenario 3')                            12d12s19
           if(idorb(ic+i).lt.isorb(io+j))then                           9d12s19
            in=idorb(ic+i)                                              9d12s19
            ix=isorb(io+j)                                              9d12s19
           else                                                         9d12s19
            ix=idorb(ic+i)                                              9d12s19
            in=isorb(io+j)                                              9d12s19
           end if                                                       9d12s19
           icol=((ix*(ix-1))/2)+in-1                                    9d12s19
           if(icol.ge.ilcol.and.icol.le.ihcol)then                      10d17s20
            if(ldebug)write(6,*)('raw col '),icol
            icol=icol-ilcol                                             10d17s20
            if(ldebug)write(6,*)('store to col '),icol
            jj=0                                                         9d10s19
            itesta=0                                                    11d16s20
            do k=0,nclo-1                                                9d10s19
             if(k.ne.i)then                                              9d10s19
              jj=jj+1                                                    9d10s19
              iorb(jj)=idorb(ic+k)                                       9d10s19
              itesta=ibset(itesta,idorb(ic+k))                          11d16s20
             end if                                                      9d10s19
            end do                                                       9d10s19
            jjo=jj+1                                                     9d10s19
            iput=0                                                       9d10s19
            do k=0,nopen-1                                               9d10s19
             if(k.ne.j)then                                              9d10s19
              if(isorb(io+k).lt.idorb(ic+i))then                         9d10s19
               jj=jj+1                                                   9d10s19
               iorb(jj)=isorb(io+k)                                      9d10s19
               itesta=ibset(itesta,isorb(io+k)+32)                      11d18s20
              else                                                       9d10s19
               if(iput.eq.0)then                                         9d10s19
                jj=jj+1                                                  9d10s19
                iorb(jj)=idorb(ic+i)                                     9d10s19
                itesta=ibset(itesta,idorb(ic+i)+32)                     11d18s20
                inst(2)=k+1                                              9d13s19
                iput=1                                                   9d10s19
               end if                                                    9d10s19
               jj=jj+1                                                   9d10s19
               iorb(jj)=isorb(io+k)                                      9d10s19
               itesta=ibset(itesta,isorb(io+k)+32)                      11d18s20
              end if                                                     9d10s19
             else if(iput.eq.0)then                                      11d18s19
              if(idorb(ic+i).lt.isorb(io+k))then
               jj=jj+1                                                    11d18s19
               iorb(jj)=idorb(ic+i)                                       11d18s19
               itesta=ibset(itesta,idorb(ic+i)+32)                      11d18s20
               inst(2)=k+1                                                11d18s19
               iput=1                                                     11d18s19
              end if                                                     11d18s19
             end if                                                      9d10s19
            end do                                                       9d10s19
            if(iput.eq.0)then                                            9d10s19
             jj=jj+1                                                     9d10s19
             iorb(jj)=idorb(ic+i)                                        9d10s19
             itesta=ibset(itesta,idorb(ic+i)+32)                        11d18s20
             inst(2)=nopen+1                                             11d18s19
            end if                                                       9d10s19
            if(ldebug)write(6,*)('inst(1),2: '),inst(1),inst(2)          12d12s19
            if(idorb(ic+i).lt.isorb(io+j))then                           11d18s19
             iextra=1                                                    9d16s19
            else                                                         9d16s19
             iextra=0                                                    9d16s19
            end if                                                       9d16s19
            if(ldebug)write(6,*)('looking for '),(iorb(ij),ij=1,nclom),  12d12s19
     $          (':'),(iorb(ij),ij=jjo,jj)                              12d12s19
            nfcnrel=iptrto(jptrto)-nff2(nclo,isbc,2)                    11d16s20
            jarg=nclo-mdon
            ndetj=ndetc(jarg+1)                                           11d19s20
            do l=1,4
             lll=nff2(nclo,isbc,2+l)+ncsf2(l,jarg)*nfcnrel
             joffg(l)=ivcc(l)+lll                                       11d19s20
            end do
            ivsp=ibcoff                                                  9d12s19
            ivtp=ivsp+ndetj*nrootu                                       9d12s19
            ibcoff=ivtp+ndetj*nrootu                                     9d12s19
            call enough('cont2ecsfa.  6',bc,ibc)
            iaorb=ibcoff                                                 9d16s19
            ivtmp=iaorb+ndetj*nopen                                      3d13s20
            ibcoff=ivtmp+ndetj*ncsf2(1,jarg)                                  9d16s19
            call enough('cont2ecsfa.  7',bc,ibc)
            call fetchcsf2(nopen,ismultm,nfcnx,bc(ivtmp),idum,           9d16s19
     $          ibc(iaorb),ismultm,idum,bc,ibc)                         11d15s22
            ibcoff=ivtmp                                                 9d16s19
            call matchc(inst,ibc(ias),ndet,nopen,3,ndetj,ibc(iaorb),     9d16s19
     $          nopen,bc(ivdet),bc(ivsp),bc(ivtp),nrootu,ismult,
     $          ncsf2(1,jarg),joffg,nct2c,icol,ndimv,iextra,vint(ioff,1)11d4s19
     $          ,ndim,lcrenorm,bc(ivnorm),irooto,fact3jst,fact3jht,     9d28s20
     $          fact3jlt,fact3jst2,fact3jht2,srhh,nrtot,bc,ibc)         11d14s22
            ibcoff=ivsp                                                  9d16s19
           end if                                                       10d17s20
          end if                                                        9d10s19
          jptrto=jptrto+1                                                11d16s20
         end do
        end do                                                          9d10s19
       end if                                                           9d10s19
       if(nopen.gt.1)then                                               9d10s19
        do i=0,nopen-1                                                  9d10s19
         do j=0,i-1                                                     9d10s19
          isij=multh(ism(isorb(io+j)),ism(isorb(io+i)))                  9d10s19
          if(multh(isij,isbr).eq.isbc)then                              9d10s19
           if(ldebug)write(6,*)('scenario 4')                            12d12s19
           icol=((isorb(io+i)*(isorb(io+i)-1))/2)+isorb(io+j)-1         9d12s19
           if(icol.ge.ilcol.and.icol.le.ihcol)then                      10d17s20
            if(ldebug)write(6,*)('icol = '),icol                         12d12s19
            icol=icol-ilcol                                             10d17s20
            if(ldebug)write(6,*)('store to col '),icol
            inst(1)=j+1                                                  9d10s19
            inst(2)=i+1                                                  9d10s19
            jj=0                                                         9d10s19
            itesta=0                                                    11d16s20
            do k=0,nclo-1                                               11d16s20
             itesta=ibset(itesta,idorb(ic+k))                           11d16s20
            end do                                                      11d16s20
            do k=0,nopen-1                                               9d10s19
             if(k.ne.i.and.k.ne.j)then                                   9d10s19
              jj=jj+1                                                    9d10s19
              iorb(jj)=isorb(io+k)                                       9d10s19
              itesta=ibset(itesta,isorb(io+k)+32)                       11d18s20
             end if                                                      9d10s19
            end do                                                       9d10s19
            nfcnrel=iptrto(jptrto)-nff2(nclop,isbc,2)                   11d16s20
            jarg=nclop-mdon
            ndetj=ndetc(jarg+1)                                           11d19s20
            do l=1,4
             lll=nff2(nclop,isbc,2+l)+ncsf2(l,jarg)*nfcnrel
             joffg(l)=ivcc(l)+lll                                       11d19s20
            end do
            if(nopenmm.lt.ismultm)then                                   10d2s20
             isss=ismultm-2                                              10d2s20
             call fetchcsf2(nopenmm,isss,nfcnx,bc(ibcoff),idum,
     $           ibc(ibcoff),isss,ndetj,bc,ibc)                         11d15s22
             ncsf2x(1)=0                                                 10d2s20
             ncsf2x(2)=nfcnx                                             10d2s20
             ncsf2x(3)=0                                                 10d2s20
             ncsf2x(4)=0                                                 10d2s20
            else                                                         10d2s20
             isss=ismultm                                                10d2s20
            end if
            ivsp=ibcoff                                                  9d12s19
            ivtp=ivsp+ndetj*nrootu                                       9d12s19
            ibcoff=ivtp+ndetj*nrootu                                     9d12s19
            call enough('cont2ecsfa.  8',bc,ibc)
            iaorb=ibcoff                                                 9d16s19
            ivtmp=iaorb+ndetj*nopenmm                                    3d13s20
            ibcoff=ivtmp+ndetj*ndetj                                    11d17s20
            call enough('cont2ecsfa.  9',bc,ibc)
            call fetchcsf2(nopenmm,isss,nfcnx,bc(ivtmp),idum,            10d2s20
     $          ibc(iaorb),isss,idum,bc,ibc)                            11d15s22
            ibcoff=ivtmp                                                 9d16s19
            if(nopenmm.lt.ismultm)then                                   10d2s20
             call matchc(inst,ibc(ias),ndet,nopen,4,ndetj,ibc(iaorb),     9d16s19
     $          nopenmm,bc(ivdet),bc(ivsp),bc(ivtp),nrootu,ismult,      9d12s19
     $          ncsf2x,joffg,nct2c,icol,ndimv,0,vint(ioff,1),           10d2s20
     $          ndim,lcrenorm,bc(ivnorm),irooto,fact3jst,fact3jht,      9d28s20
     $          fact3jlt,fact3jst2,fact3jht2,srhh,nrtot,bc,ibc)         11d14s22
            else                                                         10d2s20
             call matchc(inst,ibc(ias),ndet,nopen,4,ndetj,ibc(iaorb),     9d16s19
     $          nopenmm,bc(ivdet),bc(ivsp),bc(ivtp),nrootu,ismult,      9d12s19
     $          ncsf2(1,jarg),joffg,nct2c,icol,ndimv,0,vint(ioff,1),    11d4s19
     $          ndim,lcrenorm,bc(ivnorm),irooto,fact3jst,fact3jht,      9d28s20
     $          fact3jlt,fact3jst2,fact3jht2,srhh,nrtot,bc,ibc)         11d14s22
            end if                                                       10d2s20
            ibcoff=ivsp                                                  9d16s19
           end if                                                       10d17s20
          end if                                                        9d10s19
          jptrto=jptrto+1                                                11d16s20
         end do                                                         9d10s19
        end do                                                          9d10s19
       end if                                                           9d10s19
       ioff=ioff+ncsf(iarg)
      end do
      return
      end
