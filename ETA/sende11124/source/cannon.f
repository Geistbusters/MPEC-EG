c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine cannon(ih0mo,numpro,nowpro,kmats,jmats,ioooo,noc,iorbn,3d19s12
     $     idwsdebx,ieigv,ionex,nbasdwsc,nvirtc,idorel,idoub,iacto,     3d2s22
     $     iden1,noocc,nbasisp,iovr,ifockeig,nlzz,iorbsym,iorbsymz,     7d21s21
     $     iorbsymc,bc,ibc)                                             11d9s22
      implicit real*8 (a-h,o-z)                                         3d16s12
c
c     cannonicalize orbitals
c
      logical lpr,motion                                                5d26s22
      include "common.hf"                                               3d16s12
      include "common.store"                                            3d16s12
      include "common.print"
      dimension kmats(1),jmats(1),ioooo(1),noc(1),iorbn(8),ieigv(8),
     $     ionex(1),nbasdwsc(8),nvirtc(8),idoub(*),iacto(*),iden1(*),   4d5s18
     $     noocc(*),nbasisp(8),iovr(8),nsyo(8),nxf(8),ifockeig(*),      4d15s21
     $     iorbsym(*),iorbsymz(*),iorbsymc(*)                           7d21s21
      motion=.false.
      if(motion)write(6,*)                                              5d26s22
     $     ('we are in cannon, but just going through the motions')     5d26s22
      if(iprtr(25).eq.0)then                                            7d21s21
       lpr=.false.                                                      7d21s21
       idwsdeb=0                                                        3d2s22
      else                                                              7d21s21
       lpr=.true.                                                       7d21s21
       idwsdeb=1000
      end if                                                            7d21s21
      noff=0                                                            5d14s19
      do isb=1,nsymb
       nsyo(isb)=noff                                                   5d14s19
       noff=noff+nbasdwsc(isb)*nbasdwsc(isb)                            5d14s19
      end do                                                            5d14s19
      efromf=0d0
      if(lpr)then                                                       7d21s21
       write(6,*)('in cannon: ')
       write(6,*)('idoub: '),(idoub(isb),isb=1,nsymb)
       write(6,*)('noc: '),(noc(isb),isb=1,nsymb)
       write(6,*)('iacto: '),(iacto(isb),isb=1,nsymb)
       write(6,*)('starting orbitals ')
       do isb=1,nsymb                                                   3d2s22
        if(nbasdwsc(isb).gt.0)then                                      3d2s22
         write(6,*)('for symmetry '),isb
         call prntm2(bc(iorbn(isb)),nbasdwsc(isb),nbasdwsc(isb),
     $        nbasdwsc(isb))
        end if                                                          3d2s22
       end do                                                           3d2s22
      end if                                                            7d21s21
      if(idwsdeb.gt.10)then                                             3d19s12
       write(6,*)('cannon'),nsymb,ih0mo                                           3d16s12
      end if                                                            3d19s12
      do isb=1,nsymb                                                    9d30s20
       nxf(isb)=0                                                       9d30s20
       if(idoub(isb).gt.0)then                                          9d30s20
        if(iacto(isb).gt.0)then                                         9d30s20
         do ia=0,iacto(isb)-1                                           9d30s20
          ii=ia+iacto(isb)*ia                                           9d30s20
          delta1=abs(bc(iden1(isb)+ii)-2d0)                             9d30s20
          delta2=0d0                                                    9d30s20
          idelta2=0
          do ja=ia+1,iacto(isb)-1                                       9d30s20
           ji=ja+iacto(isb)*ia                                          9d30s20
           delta2=delta2+bc(iden1(isb)+ji)**2                           9d30s20
           idelta2=idelta2+1                                            9d30s20
          end do                                                        9d30s20
          if(idelta2.gt.0)delta2=sqrt(delta2/dfloat(idelta2))           9d30s20
         end do                                                         9d30s20
        end if                                                          9d30s20
       end if                                                           9d30s20
      end do                                                            9d30s20
      do isb=1,nsymb                                                    3d16s12
       if(nbasdwsc(isb).gt.0)then                                        3d16s12
         nmrow=nbasisp(isb)                                              2d1s19
         if(idorel.ne.0)nmrow=nmrow*2                                     2d1s19
        if(idoub(isb).gt.0)then                                         4d5s18
         if(idwsdeb.gt.10)write(6,*)('for dd part')                     4d5s18
         idoubb=idoub(isb)+nxf(isb)                                     9d30s20
         ifock=ibcoff                                                     3d16s12
         ibcoff=ifock+idoubb*idoubb                                     9d30s20
         call enough('cannon.  1',bc,ibc)
         if(nowpro.eq.0)then                                              3d16s12
          jfock=ifock                                                   3d19s12
          do i=1,idoubb                                                 9d30s20
           iad=ih0mo+nsyo(isb)+i-1+nbasdwsc(isb)*(i-1)                  5d14s19
           efromf=efromf+bc(iad)                                        2d19s14
           do j=1,idoubb                                                9d30s20
            iad=ih0mo+nsyo(isb)+j-1+nbasdwsc(isb)*(i-1)                 5d14s19
            bc(jfock)=bc(iad)                                           3d19s12
            jfock=jfock+1                                               3d19s12
           end do                                                       3d19s12
          end do                                                          3d16s12
         else                                                             3d16s12
          do i=0,idoubb*idoubb-1                                        9d30s20
           bc(ifock+i)=0d0                                                3d16s12
          end do                                                          3d16s12
         end if                                                           3d16s12
         do isa=1,nsymb                                                   3d16s12
          if(noc(isa).gt.0)then                                         5d10s19
           do i=1,nsdlk                                                   3d16s12
            if(isblk(1,i).eq.isb.and.isblk(2,i).eq.isa.and.              3d16s12
     $         isblk(3,i).eq.isa.and.isblk(4,i).eq.isb.and.              3d16s12
     $           noc(isb).gt.0)then                                      3d16s12
             call ilimts(noc(isa),noc(isb),numpro,nowpro,il,ih,i1s,i1e, 3d16s12
     $          i2s,i2e)                                                3d16s12
             i10=i1s                                                     3d16s12
             i1n=noc(isa)                                                3d16s12
             if(isa.eq.isb)then                                          3d16s12
              nrow=(noc(isa)*(noc(isa)+1))/2                             3d16s12
              ii=ioooo(i)-1                                              3d16s12
             else                                                        3d16s12
              nrow=noc(isa)*noc(isb)                                     3d16s12
              ii=ioooo(i)-1-noc(isb)                                      3d16s12
             end if                                                      3d16s12
             do i2=i2s,i2e                                               3d16s12
              if(i2.eq.i2e)i1n=i1e                                       3d16s12
              do i1=i10,i1n                                              3d16s12
               if(i1.le.idoub(isa).and.i2.le.idoubb)then                9d30s20
                if(isa.eq.isb)then
                 do i3=1,idoubb                                         9d30s20
                  in=min(i3,i1)                                           3d16s12
                  ix=max(i3,i1)                                           3d16s12
                  ixn=((ix*(ix-1))/2)+in+ii                               3d16s12
                  jfock=ifock+i3-1+idoubb*(i2-1)                        9d30s20
                  bc(jfock)=bc(jfock)-bc(ixn)                             3d16s12
                 end do                                                   3d16s12
                else                                                      3d16s12
                 do i3=1,idoubb                                         9d30s20
                  jfock=ifock+i3-1+idoubb*(i2-1)                        9d30s20
                  ixn=ii+i3+i1*noc(isb)                                    3d16s12
                  bc(jfock)=bc(jfock)-bc(ixn)                             3d16s12
                 end do                                                    3d16s12
                end if                                                    3d16s12
               else if(i1.gt.idoub(isa).and.i2.le.idoub(isb))then       5d9s19
c     (i3i1|i1i2)->Fi3i2 i.e. k part
                i1a=i1-idoub(isa)-1                                     5d9s19
                if(isa.eq.isb)then                                      5d9s19
                 do i1b=0,iacto(isa)-1                                  5d9s19
                  i1bp=i1b+idoub(isa)+1                                 5d9s19
                  do i3=1,idoubb                                        9d30s20
                   in=min(i3,i1bp)                                      5d9s19
                   ix=max(i3,i1bp)                                      5d9s19
                   ixn=((ix*(ix-1))/2)+in+ii                               3d16s12
                   jfock=ifock+i3-1+idoubb*(i2-1)                       9d30s20
                   iad3=iden1(isa)+i1b+iacto(isa)*i1a                    5d9s19
                   bc(jfock)=bc(jfock)-bc(ixn)*0.5d0*bc(iad3)            5d9s19
                  end do                                                5d9s19
                 end do                                                   3d16s12
                else                                                    5d9s19
                 do i1b=0,iacto(isa)-1                                  5d9s19
                  i1bp=i1b+idoub(isa)+1                                 5d9s19
                  do i3=1,idoubb                                        9d30s20
                   jfock=ifock+i3-1+idoubb*(i2-1)                       9d30s20
                   ixn=ii+i3+i1bp*noc(isb)                                    3d16s12
                   iad3=iden1(isa)+i1b+iacto(isa)*i1a                    5d9s19
                   bc(jfock)=bc(jfock)-bc(ixn)*0.5d0*bc(iad3)            5d9s19
                  end do                                                5d9s19
                 end do                                                    3d16s12
                end if                                                  5d9s19
               end if                                                   4d5s18
               ii=ii+nrow                                                3d16s12
              end do                                                     3d16s12
              i10=1                                                      3d16s12
             end do                                                      3d16s12
            else if(isblk(1,i).eq.isa.and.isblk(2,i).eq.isb.and.        3d19s12
     $         isblk(3,i).eq.isa.and.isblk(4,i).eq.isb.and.             3d19s12
     $           noc(isb).gt.0)then                                     5d10s19
             call ilimts(noc(isa),noc(isb),numpro,nowpro,il,ih,i1s,i1e, 3d16s12
     $          i2s,i2e)                                                3d16s12
             i10=i1s                                                     3d16s12
             i1n=noc(isa)                                                3d16s12
             nrow=noc(isa)*noc(isb)                                     3d16s12
             ii=ioooo(i)-1-noc(isa)                                      3d16s12
             do i2=i2s,i2e                                               3d16s12
              if(i2.eq.i2e)i1n=i1e                                       3d16s12
              do i1=i10,i1n                                              3d16s12
               if(i1.le.idoub(isa).and.i2.le.idoubb)then                9d30s20
                do i3=1,idoubb                                          9d30s20
                 jfock=ifock+i3-1+idoubb*(i2-1)                         9d30s20
                 ixn=ii+i1+i3*noc(isa)                                    3d16s12
                 bc(jfock)=bc(jfock)-bc(ixn)                             3d16s12
                end do                                                    3d16s12
               else if(i1.gt.idoub(isa).and.i2.le.idoubb)then           9d30s20
                i1a=i1-idoub(isa)-1
                do i1b=0,iacto(isa)-1                                   5d9s19
                 i1bp=i1b+1+idoub(isa)                                  5d9s19
                 do i3=1,idoubb                                         9d30s20
                  jfock=ifock+i3-1+idoubb*(i2-1)                        9d30s20
                  ixn=ii+i1bp+i3*noc(isa)                               5d9s19
                  iad3=iden1(isa)+i1b+iacto(isa)*i1a                    5d9s19
                  bc(jfock)=bc(jfock)-bc(ixn)*0.5d0*bc(iad3)             5d9s19
                 end do                                                 5d9s19
                end do                                                    3d16s12
               end if                                                   4d5s18
               ii=ii+nrow                                                3d16s12
              end do                                                     3d16s12
              i10=1                                                      3d16s12
             end do                                                      3d16s12
            else if(isblk(1,i).eq.isb.and.isblk(2,i).eq.isa.and.        3d19s12
     $         isblk(3,i).eq.isb.and.isblk(4,i).eq.isa.and.             3d19s12
     $           noc(isb).gt.0)then                                     5d10s19
             call ilimts(noc(isb),noc(isa),numpro,nowpro,il,ih,i1s,i1e, 3d16s12
     $          i2s,i2e)                                                3d16s12
             i10=i1s                                                     3d16s12
             i1n=noc(isb)                                                3d16s12
             nrow=noc(isa)*noc(isb)                                     3d16s12
             ii=ioooo(i)-1-noc(isb)                                     2d15s13
             do i2=i2s,i2e                                               3d16s12
              if(i2.eq.i2e)i1n=i1e                                       3d16s12
              do i1=i10,i1n                                              3d16s12
               if(i1.le.idoubb.and.i2.le.idoub(isa))then                9d30s20
                do i3=1,idoubb                                          9d30s20
                 jfock=ifock+i3-1+idoubb*(i1-1)                         9d30s20
                 ixn=ii+i3+i2*noc(isb)                                   2d15s13
                 bc(jfock)=bc(jfock)-bc(ixn)                             3d16s12
                end do                                                    3d16s12
               else if(i1.le.idoubb.and.i2.gt.idoub(isa))then           9d30s20
                i2a=i2-idoub(isa)-1                                     5d9s19
                do i2b=0,iacto(isa)-1                                   5d9s19
                 i2bp=i2b+1+idoub(isa)                                  5d9s19
                 do i3=1,idoubb                                         9d30s20
                  jfock=ifock+i3-1+idoubb*(i1-1)                        9d30s20
                  ixn=ii+i3+i2bp*noc(isb)                               5d9s19
                  iad3=iden1(isa)+i2b+iacto(isa)*i2a                     5d9s19
                  bc(jfock)=bc(jfock)-bc(ixn)*0.5d0*bc(iad3)             5d9s19
                 end do                                                 5d9s19
                end do                                                    3d16s12
               end if                                                   4d5s18
               ii=ii+nrow                                                3d16s12
              end do                                                     3d16s12
              i10=1                                                      3d16s12
             end do                                                      3d16s12
            end if                                                       3d16s12
            if(isblk(1,i).eq.isa.and.isblk(2,i).eq.isa.and.               3d16s12
     $       isblk(3,i).eq.isb.and.isblk(4,i).eq.isb.and.               3d19s12
     $           noc(isb).gt.0)then                                     5d10s19
             call ilimts(noc(isb),noc(isb),numpro,nowpro,il,ih,i1s,i1e,  3d16s12
     $          i2s,i2e)                                                3d16s12
             i10=i1s                                                     3d16s12
             i1n=noc(isb)                                                3d16s12
             ii=ioooo(i)-1                                               3d16s12
             nrow=(noc(isa)*(noc(isa)+1))/2                              3d16s12
             do i2=i2s,i2e                                               3d16s12
              if(i2.eq.i2e)i1n=i1e                                       3d16s12
              do i1=i10,i1n                                              3d16s12
               if(i1.le.idoubb.and.i2.le.idoubb)then                    9d30s20
                jfock=ifock+i1-1+idoubb*(i2-1)                          9d30s20
                do i34=1,idoub(isa)                                     4d5s18
                 ixn=ii+((i34*(i34+1))/2)                                 3d16s12
                 bc(jfock)=bc(jfock)+2d0*bc(ixn)                          3d16s12
                end do                                                    3d16s12
                do i3=0,iacto(isa)-1                                    5d9s19
                 i3p=i3+idoub(isa)+1                                    5d10s19
                 do i4=0,iacto(isa)-1                                   5d9s19
                  i4p=i4+idoub(isa)+1                                   5d10s19
                  ix=max(i3p,i4p)                                       5d9s19
                  in=min(i3p,i4p)                                       5d9s19
                  ixn=ii+((ix*(ix-1))/2)+in                             5d10s19
                  iad3=iden1(isa)+i4+iacto(isa)*i3                      5d9s19
                  bc(jfock)=bc(jfock)+bc(ixn)*bc(iad3)                  5d9s19
                 end do                                                 5d9s19
                end do                                                  5d9s19
               end if                                                   4d5s18
               ii=ii+nrow                                                3d16s12
              end do                                                     3d16s12
              i10=1                                                      3d16s12
             end do                                                      3d16s12
            end if                                                       3d16s12
           end do                                                       3d19s12
          end if                                                        3d19s12
         end do                                                         3d19s12
         iarg=idoubb*idoubb                                             9d30s20
         call dws_gsumf(bc(ifock),iarg)                                 3d19s12
         if(lpr)then                                                    7d21s21
          write(6,*)('dd fock matrix ')
          call prntm2(bc(ifock),idoubb,idoubb,idoubb)                   9d30s20
         end if                                                         3d19s12
         ivec=ibcoff                                                    3d19s12
         ieig=ivec+idoubb*idoubb                                        9d30s20
         isym=ieig+idoubb                                               9d30s20
         ibcoff=isym+idoubb                                             9d30s20
         call enough('cannon.  2',bc,ibc)
         nsend=idoubb*(idoubb+1)                                        9d30s20
         if(nowpro.eq.0)then
          if(nlzz.ne.0)then                                             4d19s21
           if(nlzz.eq.2)then                                            4d19s21
            do i=0,idoubb-1                                             4d19s21
             ibc(isym+i)=ibc(iorbsym(isb)+i)                            4d19s21
            end do                                                      4d19s21
           else                                                         4d19s21
            do i=0,idoubb-1                                             4d19s21
             ibc(isym+i)=ibc(iorbsymz(isb)+i)+100*ibc(iorbsym(isb)+i)   8d23s21
            end do                                                      4d19s21
           end if                                                       4d19s21
           call diagy(idoubb,bc(ifock),bc(ieig),bc(ivec),ibc(isym),bc,  11d14s22
     $          ibc,0,idum,dum)                                         9d1s23
           do i=0,idoubb-1                                              7d21s21
            ibc(iorbsymc(isb)+i)=ibc(isym+i)                            7d21s21
           end do                                                       7d21s21
          else                                                          4d19s21
           call diagx(idoubb,bc(ifock),bc(ieig),bc(ivec),bc(isym),bc,   11d14s22
     $          ibc)                                                    11d14s22
          end if                                                        4d19s21
         end if
         call dws_bcast(bc(ivec),nsend)
         do i=1,idoub(isb)                                              4d5s18
          im=i-1                                                        2d28s21
          ii=ieig+im                                                    2d28s21
          bc(ifockeig(isb)+im)=bc(ii)                                   2d28s21
          efromf=efromf+bc(ii)
         end do
         if(idwsdeb.gt.10)then                                          3d19s12
          write(6,*)('eigenvalues ')
          call prntm2(bc(ieig),1,idoubb,1)                              9d30s20
          write(6,*)('eigenvectors ')
          call prntm2(bc(ivec),idoubb,idoubb,idoubb)                    9d30s20
         end if                                                         3d19s12
         itmp=ibcoff                                                    3d19s12
         ibcoff=itmp+nmrow*idoubb                                       9d30s20
         if(nmrow.gt.0.and.idoubb.gt.0)then                             9d30s20
         call dgemm('n','n',nmrow,idoubb,idoubb,1d0,                    9d30s20
     $        bc(iorbn(isb)),nmrow,bc(ivec),idoubb,0d0,                 9d30s20
     $        bc(itmp),nmrow,                                            2d1s19
     d' cannon.  1')
         end if                                                         2d22s19
         if(idwsdeb.gt.10)then                                          3d19s12
          write(6,*)('new occ orbs '),isb
          call prntm2(bc(itmp),nmrow,idoubb,nmrow)                      9d30s20
         end if                                                         3d19s12
         do i=1,idoubb                                                  9d30s20
          bc(ieigv(isb)+i-1)=bc(ieig+i-1)                               3d19s12
          if(.not.motion)then                                           5d26s22
           do j=1,nmrow                                                   2d1s19
            iad1=iorbn(isb)+j-1+nmrow*(i-1)                               2d1s19
            iad2=itmp+j-1+nmrow*(i-1)                                     2d1s19
            bc(iad1)=bc(iad2)                                            3d19s12
           end do                                                        3d19s12
          end if                                                        5d26s22
         end do                                                         3d19s12
         if(idwsdeb.gt.10)then                                          3d19s12
        write(6,*)('orbs after storing ')
        call prntm2(bc(iorbn(isb)),nmrow,nbasdwsc(isb),nmrow)             2d1s19
         end if                                                         3d19s12
         ibcoff=ifock                                                   3d19s12
        end if                                                          3d19s12
        if(iacto(isb).gt.nxf(isb))then                                  9d30s20
c
c     to order nos so biggest occupation numbers come first,            4d5s18
c     diagonalize negative of density                                   4d5s18
c                                                                       4d5s18
        iactou=iacto(isb)-nxf(isb)                                      9d30s20
        idencp=ibcoff                                                   9d30s20
        ibcoff=idencp+iactou*iactou                                     9d30s20
        call enough('cannon.  3',bc,ibc)
        do ia=0,iactou-1                                                9d30s20
         iap=ia+nxf(isb)                                                9d30s20
         do ib=0,iactou-1                                               9d30s20
          ibp=ib+nxf(isb)                                               9d30s20
          iad1=idencp+ib+iactou*ia                                      9d30s20
          iad2=iden1(isb)+ibp+iacto(isb)*iap                            9d30s20
          bc(iad1)=-bc(iad2)                                            9d30s20
         end do                                                         9d30s20
        end do                                                          9d30s20
        i=iden1(isb)+nxf(isb)*(iacto(isb)+1)                            9d30s20
         ieig=ibcoff                                                    4d5s18
         ivec=ieig+iactou                                               9d30s20
         isym=ivec+iactou*iactou                                        9d30s20
         ibcoff=isym+iactou                                             9d30s20
         call enough('cannon.  4',bc,ibc)
         nsend=iactou*(iactou+1)                                        9d30s20
         if(nowpro.eq.0)then                                            4d5s18
          if(nlzz.ne.0)then                                             4d19s21
           iorbsymu=iorbsym(isb)+idoub(isb)+nxf(isb)                    4d19s21
           iorbsymcu=iorbsymc(isb)+idoub(isb)+nxf(isb)                  7d21s21
           if(nlzz.eq.2)then                                            4d19s21
            do i=0,iactou-1                                             4d19s21
             ibc(isym+i)=ibc(iorbsymu+i)                                4d19s21
            end do                                                      4d19s21
           else                                                         4d19s21
            iorbsymzu=iorbsymz(isb)+idoub(isb)+nxf(isb)                 4d19s21
            do i=0,iactou-1                                             4d19s21
             ibc(isym+i)=ibc(iorbsymzu+i)+100*ibc(iorbsymu+i)           8d23s21
            end do                                                      4d19s21
           end if                                                       4d19s21
           call diagy(iactou,bc(idencp),bc(ieig),bc(ivec),ibc(isym),    11d14s22
     $          bc,ibc,0,idum,dum)                                      9d1s23
           do i=0,iactou-1                                              7d21s21
            ibc(iorbsymcu+i)=ibc(isym+i)                                7d21s21
           end do                                                       7d21s21
          else                                                          4d19s21
           call diagx(iactou,bc(idencp),bc(ieig),bc(ivec),bc(isym),bc,  11d14s22
     $          ibc)                                                    11d14s22
          end if                                                        4d19s21
         end if                                                         4d5s18
         call dws_bcast(bc(ieig),nsend)                                 4d5s18
         do i=0,nxf(isb)-1                                              9d30s20
          bc(noocc(isb)+i)=2d0                                          9d30s20
         end do                                                         9d30s20
         do i=0,iactou-1                                                9d30s20
          ip=i+nxf(isb)                                                 9d30s20
          bc(noocc(isb)+ip)=-bc(ieig+i)                                 9d30s20
         end do                                                         9d30s20
c                                                                       12d24s19
c     now switch density back so its correct for futher fock matrix     12d24s19
c     constructions                                                     12d24s19
c                                                                       12d24s19
         if(idwsdeb.gt.10)then                                          3d19s12
          write(6,*)('occupation nos ')
          do i=0,iactou-1                                               9d30s20
           bc(ieig+i)=-bc(ieig+i)                                       4d5s18
          end do                                                        4d5s18
          call prntm2(bc(ieig),1,iactou,1)                              9d30s20
          write(6,*)('eigenvectors ')
          call prntm2(bc(ivec),iactou,iactou,iactou)                    9d30s20
         end if                                                         3d19s12
         itmp=ibcoff                                                    3d19s12
         ibcoff=itmp+nmrow*iactou                                       9d30s20
         jorbn=iorbn(isb)+(idoub(isb)+nxf(isb))*nmrow                   9d30s20
         if(nmrow.gt.0.and.iactou.gt.0)then                             9d30s20
         call dgemm('n','n',nmrow,iactou,iactou,1d0,                    9d30s20
     $        bc(jorbn),nmrow,bc(ivec),iactou,0d0,                      9d30s20
     $        bc(itmp),nmrow,                                            2d1s19
     d' cannon.  2')
         end if                                                         2d22s19
         if(idwsdeb.gt.10)then                                          3d19s12
          write(6,*)('new active orbs '),isb
          call prntm2(bc(itmp),nmrow,iactou,nmrow)                      9d30s20
         end if                                                         3d19s12
         if(.not.motion)then                                            5d26s22
          do i=1,iactou                                                  9d30s20
           do j=1,nmrow                                                   2d1s19
            iad1=jorbn+j-1+nmrow*(i-1)                                    2d1s19
            iad2=itmp+j-1+nmrow*(i-1)                                     2d1s19
            bc(iad1)=bc(iad2)                                            3d19s12
           end do                                                        3d19s12
          end do                                                         3d19s12
         end if                                                         5d26s22
         if(idwsdeb.gt.10)then                                          3d19s12
        write(6,*)('orbs after storing ')
        call prntm2(bc(iorbn(isb)),nmrow,nbasdwsc(isb),nmrow)
         end if                                                         3d19s12
         ibcoff=ieig                                                    4d5s18
         ibcoff=idencp                                                  9d30s20
        end if                                                          4d5s18
        if(nvirtc(isb).gt.0)then                                         3d19s12
         if(idwsdeb.gt.10)write(6,*)('for vv part'),isb,noc(isb)                     3d19s12
         ifock=ibcoff                                                     3d16s12
         ibcoff=ifock+nvirtc(isb)*nvirtc(isb)                             3d19s12
         jorbsym=iorbsym(isb)+noc(isb)-1                                4d20s21
         call enough('cannon.  5',bc,ibc)
         if(nowpro.eq.0)then                                              3d16s12
          jfock=ifock                                                   3d19s12
          do i=1,nvirtc(isb)                                             3d19s12
           ip=i+noc(isb)                                                3d19s12
           do j=1,nvirtc(isb)                                            3d19s12
            jp=j+noc(isb)                                               3d19s12
            iad=ih0mo+nsyo(isb)+jp-1+nbasdwsc(isb)*(ip-1)               5d14s19
            bc(jfock)=bc(iad)                                           3d19s12
            jfock=jfock+1                                               3d19s12
           end do                                                       3d19s12
          end do                                                          3d16s12
         else                                                             3d16s12
          do i=0,nvirtc(isb)*nvirtc(isb)-1                                3d19s12
           bc(ifock+i)=0d0                                                3d16s12
          end do                                                          3d16s12
         end if                                                           3d16s12
         igoal=ifock+1
         do isa=1,nsymb                                                   3d16s12
          if(noc(isa).gt.0)then                                           3d16s12
           do i=1,nsdlk                                                   3d16s12
            if(isblk(1,i).eq.isa.and.isblk(2,i).eq.isa.and.               3d16s12
     $       isblk(3,i).eq.isb.and.isblk(4,i).eq.isb)then               3d16s12
             call ilimts(nvirtc(isb),nvirtc(isb),numpro,nowpro,il,ih,i1s8d31s15
     $           ,i1e,i2s,i2e)                                          8d31s15
             i10=i1s                                                     3d16s12
             i1n=nvirtc(isb)                                              3d16s12
             ii=jmats(i)-1                                               3d16s12
             nrow=(noc(isa)*(noc(isa)+1))/2                              3d16s12
             do i2=i2s,i2e                                               3d16s12
              if(i2.eq.i2e)i1n=i1e                                       3d16s12
              do i1=i10,i1n                                              3d16s12
               jfock=ifock+i1-1+nvirtc(isb)*(i2-1)                        3d19s12
               do i34=1,idoub(isa)                                      4d5s18
                ixn=ii+((i34*(i34+1))/2)                                 3d16s12
                bc(jfock)=bc(jfock)+2d0*bc(ixn)                          3d16s12
               end do                                                    3d16s12
               do i3=0,iacto(isa)-1                                     5d10s19
                i3p=i3+idoub(isa)+1                                     5d10s19
                do i4=0,iacto(isa)-1                                    5d10s19
                 i4p=i4+idoub(isa)+1                                    5d10s19
                 ix=max(i3p,i4p)                                        5d10s19
                 in=min(i3p,i4p)                                        5d10s19
                 ixn=ii+((ix*(ix-1))/2)+in                              5d10s19
                 iad3=iden1(isa)+i4+iacto(isa)*i3                       5d10s19
                 bc(jfock)=bc(jfock)+bc(ixn)*bc(iad3)                   5d10s19
                end do                                                  5d10s19
               end do                                                   5d10s19
               ii=ii+nrow                                                3d16s12
              end do                                                     3d16s12
              i10=1                                                      3d16s12
             end do                                                      3d16s12
            end if                                                        3d16s12
           end do                                                         3d16s12
           do i=1,nsdlkk                                                  3d16s12
            if(isblkk(1,i).eq.isa.and.isblkk(2,i).eq.isa.and.             3d16s12
     $       isblkk(3,i).eq.isb.and.isblkk(4,i).eq.isb)then             3d16s12
             call ilimts(nvirtc(isb),nvirtc(isb),numpro,nowpro,il,ih,i1s8d31s15
     $           ,i1e,i2s,i2e)                                          8d31s15
             i10=i1s                                                     3d16s12
             i1n=nvirtc(isb)                                              3d16s12
             nocp=noc(isa)+1                                             3d16s12
             ii=kmats(i)-nocp                                            3d16s12
             nrow=noc(isa)*noc(isa)                                      3d16s12
             do i2=i2s,i2e                                               3d16s12
              if(i2.eq.i2e)i1n=i1e                                       3d16s12
              do i1=i10,i1n                                              3d16s12
               jfock=ifock+i1-1+nvirtc(isb)*(i2-1)                        3d19s12
               do i34=1,idoub(isa)                                      4d5s18
                ixn=ii+i34*nocp                                          3d16s12
                bc(jfock)=bc(jfock)-bc(ixn)                              3d16s12
               end do                                                    3d16s12
               do i3=0,iacto(isa)-1                                     5d10s19
                i3p=i3+idoub(isa)+1                                     5d10s19
                do i4=0,iacto(isa)-1                                    5d10s19
                 i4p=i4+idoub(isa)+1                                    5d10s19
                 iad3=iden1(isa)+i4+iacto(isa)*i3                       5d10s19
                 ixn=ii+i3p+noc(isa)*i4p                                5d10s19
                 bc(jfock)=bc(jfock)-bc(ixn)*0.5d0*bc(iad3)             5d10s19
                end do                                                  5d10s19
               end do                                                   5d10s19
               ii=ii+nrow                                                3d16s12
              end do                                                     3d16s12
              i10=1                                                      3d16s12
             end do                                                      3d16s12
            end if                                                        3d16s12
           end do                                                         3d16s12
          end if                                                          3d16s12
         end do                                                           3d16s12
         iarg=nvirtc(isb)*nvirtc(isb)                                     3d19s12
         call dws_gsumf(bc(ifock),iarg)                                 3d19s12
         if(idwsdeb.gt.10)then                                          3d19s12
          write(6,*)('fock matrix ')
          call prntm2(bc(ifock),nvirtc(isb),nvirtc(isb),nvirtc(isb))        3d19s12
         end if                                                         3d19s12
         ivec=ibcoff                                                    3d19s12
         ieig=ivec+nvirtc(isb)*nvirtc(isb)                                3d19s12
         isym=ieig+nvirtc(isb)
         ibcoff=isym+nvirtc(isb)
         call enough('cannon.  6',bc,ibc)
         nsend=nvirtc(isb)*(nvirtc(isb)+1)
         if(nowpro.eq.0)then
          if(nlzz.ne.0)then                                             4d19s21
           iorbsymu=iorbsym(isb)+idoub(isb)+iacto(isb)                  4d19s21
           iorbsymuc=iorbsymc(isb)+idoub(isb)+iacto(isb)                7d21s21
           if(nlzz.eq.2)then                                            4d19s21
            do i=0,nvirtc(isb)-1                                        4d19s21
             ibc(isym+i)=ibc(iorbsymu+i)                                4d19s21
            end do                                                      4d19s21
           else                                                         4d19s21
            iorbsymzu=iorbsymz(isb)+idoub(isb)+iacto(isb)               4d19s21
            do i=0,nvirtc(isb)-1                                        4d19s21
             ibc(isym+i)=ibc(iorbsymzu+i)+100*ibc(iorbsymu+i)           8d23s21
            end do                                                      4d19s21
           end if                                                       4d19s21
           call diagy(nvirtc(isb),bc(ifock),bc(ieig),bc(ivec),          4d19s21
     $         ibc(isym),bc,ibc,0,idum,dum)                             9d1s23
           do i=0,nvirtc(isb)-1                                         7d21s21
            ibc(iorbsymuc+i)=ibc(isym+i)                                7d21s21
           end do                                                       7d21s21
          else                                                          4d19s21
           call diagx(nvirtc(isb),bc(ifock),bc(ieig),bc(ivec),             3d19s12
     $         bc(isym),bc,ibc)                                         11d14s22
          end if                                                        4d19s21
         end if
         call dws_bcast(bc(ivec),nsend)
         if(idwsdeb.gt.10)then                                          3d19s12
          write(6,*)('eigenvalues ')
          call prntm2(bc(ieig),1,nvirtc(isb),1)
         end if                                                         3d19s12
         itmp=ibcoff                                                    3d19s12
         ibcoff=itmp+nmrow*nvirtc(isb)                                   2d1s19
         jorbn=iorbn(isb)+nmrow*noc(isb)                                 2d1s19
         if(nmrow.gt.0.and.nvirtc(isb).gt.0)then                        2d22s19
         call dgemm('n','n',nmrow,nvirtc(isb),nvirtc(isb),1d0,           2d1s19
     $        bc(jorbn),nmrow,bc(ivec),nvirtc(isb),0d0,                  2d1s19
     $        bc(itmp),nmrow,                                            2d1s19
     d' cannon.  3')
         end if                                                         2d22s19
         do i=1,nvirtc(isb)                                              3d19s12
          bc(ieigv(isb)+i+noc(isb)-1)=bc(ieig+i-1)                      3d19s12
          if(.not.motion)then                                           5d26s22
           do j=1,nmrow                                                   2d1s19
            iad1=jorbn+j-1+nmrow*(i-1)                                    2d1s19
            iad2=itmp+j-1+nmrow*(i-1)                                     2d1s19
            bc(iad1)=bc(iad2)                                            3d19s12
           end do                                                        3d19s12
          end if                                                        5d26s22
         end do                                                         3d19s12
         ibcoff=ifock                                                   3d19s12
        end if                                                           3d16s12
        if(idwsdeb.gt.10)then                                           3d19s12
        write(6,*)('brand new orbitals ')
        call prntm2(bc(iorbn(isb)),nmrow,nbasdwsc(isb),nmrow)
        write(6,*)('overlap matrix? '),iovr(isb)
        call prntm2(bc(iovr(isb)),nmrow,nmrow,nmrow)
        itmpo1=ibcoff
        itmpo2=itmpo1+nmrow*nbasdwsc(isb)
        ibcoff=itmpo2+nbasdwsc(isb)*nbasdwsc(isb)
        call enough('cannon.  7',bc,ibc)
        write(6,*)('ortho test')
        call dgemm('n','n',nmrow,nbasdwsc(isb),nmrow,1d0,bc(iovr(isb)),
     $       nmrow,bc(iorbn(isb)),nmrow,0d0,bc(itmpo1),nmrow,
     d' cannon.  4')
        call dgemm('t','n',nbasdwsc(isb),nbasdwsc(isb),nmrow,1d0,
     $       bc(iorbn(isb)),nmrow,bc(itmpo1),nmrow,0d0,bc(itmpo2),
     $       nbasdwsc(isb),
     d' cannon.  5')
        call prntm2(bc(itmpo2),nbasdwsc(isb),nbasdwsc(isb),nbasdwsc(isb)
     $       )
        ibcoff=itmpo1
        end if                                                          3d19s12
       end if                                                           3d19s12
      end do
      if(lpr)then                                                       3d2s22
       write(6,*)('ending orbitals ')
       do isb=1,nsymb                                                   3d2s22
        if(nbasdwsc(isb).gt.0)then                                      3d2s22
         write(6,*)('for symmetry '),isb
         call prntm2(bc(iorbn(isb)),nbasdwsc(isb),nbasdwsc(isb),
     $        nbasdwsc(isb))
        end if                                                          3d2s22
       end do                                                           3d2s22
      end if                                                            3d2s22
      return
      end
