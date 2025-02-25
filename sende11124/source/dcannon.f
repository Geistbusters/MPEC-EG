c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dcannon(ih0mo,idh0,numpro,nowpro,ioooo,id4o,kmats,     5d20s22
     $     kmatd,jmats,jmatd,noc,morb,idarot,idwsdeb,idoub,iacto,       7d14s22
     $     nvirtx,iden1,iden1d,nsblkkder,isblkkder,bc,ibc)              11d14s22
      implicit real*8 (a-h,o-z)                                         3d16s12
c
c     derivative of cannonicalization of orbitals
c
      logical lpr                                                       7d21s21
      include "common.hf"                                               3d16s12
      include "common.store"                                            3d16s12
      include "common.print"
      dimension kmats(1),jmats(1),ioooo(1),noc(1),morb(8),idarot(8),    5d19s22
     $     idoub(*),iacto(*),iden1(*),nsyo(8),id4o(*),iden1d(*),        5d20s22
     $     isblkkder(4,*),nvirtx(*),kmatd(*),jmatd(*)                   5d20s22
      lpr=idwsdeb.ne.0                                                  6d24s24
      noff=0                                                            5d14s19
      do isb=1,nsymb
       nsyo(isb)=noff                                                   5d14s19
       noff=noff+nbasdws(isb)*nbasdws(isb)                              5d19s22
      end do                                                            5d14s19
      efromf=0d0
      if(lpr)then                                                       7d21s21
       write(6,*)('in dcannon: ')
       write(6,*)('idoub: '),(idoub(isb),isb=1,nsymb)
       write(6,*)('noc: '),(noc(isb),isb=1,nsymb)
       write(6,*)('iacto: '),(iacto(isb),isb=1,nsymb)
       write(6,*)('starting orbitals ')
       do isb=1,nsymb                                                   3d2s22
        if(nbasdws(isb).gt.0)then                                       5d19s22
         write(6,*)('for symmetry '),isb
         call prntm2(bc(morb(isb)),nbasdws(isb),nbasdws(isb),           5d19s22
     $        nbasdws(isb))                                             5d19s22
        end if                                                          3d2s22
       end do                                                           3d2s22
      end if                                                            7d21s21
      if(idwsdeb.gt.10)then                                             3d19s12
       write(6,*)('dcannon'),nsymb,ih0mo                                           3d16s12
      end if                                                            3d19s12
      do isb=1,nsymb                                                    3d16s12
       if(nbasdws(isb).gt.0)then                                        3d16s12
        ix2=ibcoff                                                      5d19s22
        ibcoff=ix2+nbasdws(isb)*nbasdws(isb)                            5d19s22
        call enough('dcannon.  1',bc,ibc)
        call dgemm('n','n',nbasdws(isb),nbasdws(isb),nbasdws(isb),1d0,  5d19s22
     $       bc(morb(isb)),nbasdws(isb),bc(idarot(isb)),nbasdws(isb),   5d19s22
     $       0d0,bc(ix2),nbasdws(isb),                                  5d19s22
     d'dcannon.  1')
        if(lpr)then                                                     7d25s22
         write(6,*)('starting ix2 matrix ')
         call prntm2(bc(ix2),nbasdws(isb),nbasdws(isb),nbasdws(isb))
        end if                                                          7d25s22
        if(idoub(isb).gt.1)then                                         4d5s18
         if(idwsdeb.gt.10)write(6,*)('for dd part')                     4d5s18
         idoubb=idoub(isb)                                              5d19s22
         ifock=ibcoff                                                     3d16s12
         idfock=ifock+idoubb*idoubb                                     5d19s22
         ibcoff=idfock+idoubb*idoubb                                    5d19s22
         call enough('dcannon.  2',bc,ibc)
         if(nowpro.eq.0)then                                              3d16s12
          jfock=ifock                                                   3d19s12
          jdfock=idfock                                                 5d19s22
          do i=1,idoubb                                                 9d30s20
           do j=1,idoubb                                                9d30s20
            idad=idh0+nsyo(isb)+j-1+nbasdws(isb)*(i-1)                  5d19s22
            iad=ih0mo+nsyo(isb)+j-1+nbasdws(isb)*(i-1)                  5d19s22
            bc(jfock)=bc(iad)                                           3d19s12
            jfock=jfock+1                                               3d19s12
            bc(jdfock)=bc(idad)                                         5d19s22
            jdfock=jdfock+1                                             5d19s22
           end do                                                       3d19s12
          end do                                                          3d16s12
         else                                                             3d16s12
          do i=0,idoubb*idoubb-1                                        9d30s20
           bc(ifock+i)=0d0                                                3d16s12
           bc(idfock+i)=0d0                                             5d19s22
          end do                                                          3d16s12
         end if                                                           3d16s12
         if(lpr)then
          write(6,*)('after h0 part ')
          call prntm2(bc(ifock),idoub(isb),idoub(isb),idoub(isb))
          call prntm2(bc(idfock),idoub(isb),idoub(isb),idoub(isb))
         end if
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
              ii=id4o(i)-1                                              5d19s22
              iiu=ioooo(i)-1                                            5d19s22
             else                                                        3d16s12
              nrow=noc(isa)*noc(isb)                                     3d16s12
              ii=id4o(i)-1-noc(isb)                                     5d19s22
              iiu=ioooo(i)-1-noc(isb)                                   5d19s22
             end if                                                      3d16s12
             do i2=i2s,i2e                                               3d16s12
              if(i2.eq.i2e)i1n=i1e                                       3d16s12
              do i1=i10,i1n                                              3d16s12
               if(i1.le.idoub(isa).and.i2.le.idoubb)then                9d30s20
                if(isa.eq.isb)then
                 do i3=1,idoubb                                         9d30s20
                  in=min(i3,i1)                                           3d16s12
                  ix=max(i3,i1)                                           3d16s12
                  ixn=((ix*(ix-1))/2)+in                                5d19s22
                  iixn=ixn
                  ixnu=ixn+iiu                                          5d19s22
                  ixn=ixn+ii                                            5d19s22
                  jfock=ifock+i3-1+idoubb*(i2-1)                        9d30s20
                  jdfock=idfock+i3-1+idoubb*(i2-1)                      5d19s22
                  bc(jfock)=bc(jfock)-bc(ixnu)                          5d19s22
                  bc(jdfock)=bc(jdfock)-bc(ixn)                             3d16s12
                 end do                                                   3d16s12
                else                                                      3d16s12
                 do i3=1,idoubb                                         9d30s20
                  jfock=ifock+i3-1+idoubb*(i2-1)                        9d30s20
                  jdfock=idfock+i3-1+idoubb*(i2-1)                      5d19s22
                  ixn=i3+i1*noc(isb)                                    5d19s22
                  ixnu=iiu+ixn                                          5d19s22
                  ixn=ii+ixn                                            5d19s22
                  bc(jdfock)=bc(jdfock)-bc(ixn)                         5d19s22
                  bc(jfock)=bc(jfock)-bc(ixnu)                          5d19s22
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
                   ixn=((ix*(ix-1))/2)+in                               5d19s22
                   ixnu=iiu+ixn                                         5d19s22
                   ixn=ii+ixn                                           5d19s22
                   jfock=ifock+i3-1+idoubb*(i2-1)                       9d30s20
                   jdfock=idfock+i3-1+idoubb*(i2-1)                     5d19s22
                   iad3=iden1(isa)+i1b+iacto(isa)*i1a                    5d9s19
                   idad3=iden1d(isa)+i1b+iacto(isa)*i1a                 5d19s22
                   bc(jfock)=bc(jfock)-bc(ixnu)*0.5d0*bc(iad3)          5d19s22
                   bc(jdfock)=bc(jdfock)                                5d19s22
     $                  -0.5d0*(bc(ixn)*bc(iad3)+bc(ixnu)*bc(idad3))    5d19s22
                  end do                                                5d9s19
                 end do                                                   3d16s12
                else                                                    5d9s19
                 do i1b=0,iacto(isa)-1                                  5d9s19
                  i1bp=i1b+idoub(isa)+1                                 5d9s19
                  do i3=1,idoubb                                        9d30s20
                   jfock=ifock+i3-1+idoubb*(i2-1)                       9d30s20
                   jdfock=idfock+i3-1+idoubb*(i2-1)                     5d19s22
                   ixn=i3+i1bp*noc(isb)                                 5d19s22
                   ixnu=iiu+ixn                                         5d19s22
                   ixn=ii+ixn                                           5d19s22
                   iad3=iden1(isa)+i1b+iacto(isa)*i1a                    5d9s19
                   idad3=iden1d(isa)+i1b+iacto(isa)*i1a                 5d19s22
                   bc(jfock)=bc(jfock)-0.5d0*bc(ixnu)*bc(iad3)          5d19s22
                   bc(jdfock)=bc(jdfock)                                5d19s22
     $                  -0.5d0*(bc(ixn)*bc(iad3)+bc(ixnu)*bc(idad3))    5d19s22
                  end do                                                5d9s19
                 end do                                                    3d16s12
                end if                                                  5d9s19
               end if                                                   4d5s18
               ii=ii+nrow                                                3d16s12
               iiu=iiu+nrow                                                3d16s12
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
             ii=id4o(i)-1-noc(isa)                                      5d19s22
             iiu=ioooo(i)-1-noc(isa)                                    5d19s22
             do i2=i2s,i2e                                               3d16s12
              if(i2.eq.i2e)i1n=i1e                                       3d16s12
              do i1=i10,i1n                                              3d16s12
               if(i1.le.idoub(isa).and.i2.le.idoubb)then                9d30s20
                do i3=1,idoubb                                          9d30s20
                 jfock=ifock+i3-1+idoubb*(i2-1)                         9d30s20
                 jdfock=idfock+i3-1+idoubb*(i2-1)                       5d19s22
                 ixn=i1+i3*noc(isa)                                     5d19s22
                 ixnu=iiu+ixn                                           5d19s22
                 ixn=ii+ixn                                             5d19s22
                 bc(jfock)=bc(jfock)-bc(ixnu)                           5d19s22
                 bc(jdfock)=bc(jdfock)-bc(ixn)                          5d19s22
                end do                                                    3d16s12
               else if(i1.gt.idoub(isa).and.i2.le.idoubb)then           9d30s20
                i1a=i1-idoub(isa)-1
                do i1b=0,iacto(isa)-1                                   5d9s19
                 i1bp=i1b+1+idoub(isa)                                  5d9s19
                 do i3=1,idoubb                                         9d30s20
                  jfock=ifock+i3-1+idoubb*(i2-1)                        9d30s20
                  jdfock=idfock+i3-1+idoubb*(i2-1)                      5d19s22
                  ixn=i1bp+i3*noc(isa)                                  5d19s22
                  ixnu=iiu+ixn                                          5d19s22
                  ixn=ii+ixn                                            5d19s22
                  iad3=iden1(isa)+i1b+iacto(isa)*i1a                    5d9s19
                  idad3=iden1d(isa)+i1b+iacto(isa)*i1a                  5d19s22
                  bc(jfock)=bc(jfock)-0.5d0*bc(ixnu)*bc(iad3)           5d19s22
                  bc(jdfock)=bc(jdfock)                                 5d19s22
     $                 -0.5d0*(bc(ixn)*bc(iad3)+bc(ixnu)*bc(idad3))     5d19s22
                 end do                                                 5d9s19
                end do                                                    3d16s12
               end if                                                   4d5s18
               ii=ii+nrow                                                3d16s12
               iiu=iiu+nrow                                             5d19s22
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
             ii=id4o(i)-1-noc(isb)                                      5d19s22
             iiu=ioooo(i)-1-noc(isb)                                    5d19s22
             do i2=i2s,i2e                                               3d16s12
              if(i2.eq.i2e)i1n=i1e                                       3d16s12
              do i1=i10,i1n                                              3d16s12
               if(i1.le.idoubb.and.i2.le.idoub(isa))then                9d30s20
                do i3=1,idoubb                                          9d30s20
                 jfock=ifock+i3-1+idoubb*(i1-1)                         9d30s20
                 jdfock=idfock+i3-1+idoubb*(i1-1)                       5d19s22
                 ixn=i3+i2*noc(isb)                                     5d26s22
                 ixnu=iiu+ixn                                           5d26s22
                 ixn=ii+ixn                                             5d26s22
                 bc(jfock)=bc(jfock)-bc(ixnu)                           5d19s22
                 bc(jdfock)=bc(jdfock)-bc(ixn)                          5d19s22
                end do                                                    3d16s12
               else if(i1.le.idoubb.and.i2.gt.idoub(isa))then           9d30s20
                i2a=i2-idoub(isa)-1                                     5d9s19
                do i2b=0,iacto(isa)-1                                   5d9s19
                 i2bp=i2b+1+idoub(isa)                                  5d9s19
                 do i3=1,idoubb                                         9d30s20
                  jfock=ifock+i3-1+idoubb*(i1-1)                        9d30s20
                  jdfock=idfock+i3-1+idoubb*(i1-1)                      5d19s22
                  ixn=i3+i2bp*noc(isb)                                  5d19s22
                  ixnu=iiu+ixn                                          5d19s22
                  ixn=ii+ixn                                            5d19s22
                  iad3=iden1(isa)+i2b+iacto(isa)*i2a                     5d9s19
                  idad3=iden1d(isa)+i2b+iacto(isa)*i2a                  5d19s22
                  bc(jfock)=bc(jfock)-0.5d0*bc(ixnu)*bc(iad3)           5d19s22
                  bc(jdfock)=bc(jdfock)                                 5d19s22
     $                 -0.5d0*(bc(ixn)*bc(iad3)+bc(ixnu)*bc(idad3))     5d19s22
                 end do                                                 5d9s19
                end do                                                    3d16s12
               end if                                                   4d5s18
               ii=ii+nrow                                                3d16s12
               iiu=iiu+nrow                                             5d19s22
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
             ii=id4o(i)-1                                               5d19s22
             iiu=ioooo(i)-1                                             5d19s22
             nrow=(noc(isa)*(noc(isa)+1))/2                              3d16s12
             do i2=i2s,i2e                                               3d16s12
              if(i2.eq.i2e)i1n=i1e                                       3d16s12
              do i1=i10,i1n                                              3d16s12
               if(i1.le.idoubb.and.i2.le.idoubb)then                    9d30s20
                jfock=ifock+i1-1+idoubb*(i2-1)                          9d30s20
                jdfock=idfock+i1-1+idoubb*(i2-1)                        5d19s22
                do i34=1,idoub(isa)                                     4d5s18
                 ixn=((i34*(i34+1))/2)                                  5d19s22
                 ixnu=iiu+ixn                                           5d19s22
                 ixn=ii+ixn                                             5d19s22
                 bc(jfock)=bc(jfock)+2d0*bc(ixnu)                       5d19s22
                 bc(jdfock)=bc(jdfock)+2d0*bc(ixn)                      5d19s22
                end do                                                    3d16s12
                do i3=0,iacto(isa)-1                                    5d9s19
                 i3p=i3+idoub(isa)+1                                    5d10s19
                 do i4=0,iacto(isa)-1                                   5d9s19
                  i4p=i4+idoub(isa)+1                                   5d10s19
                  ix=max(i3p,i4p)                                       5d9s19
                  in=min(i3p,i4p)                                       5d9s19
                  ixn=((ix*(ix-1))/2)+in                                5d19s22
                  ixnu=iiu+ixn                                          5d19s22
                  ixn=ii+ixn                                            5d19s22
                  iad3=iden1(isa)+i4+iacto(isa)*i3                      5d9s19
                  idad3=iden1d(isa)+i4+iacto(isa)*i3                    5d19s22
                  bc(jfock)=bc(jfock)+bc(ixnu)*bc(iad3)                 5d19s22
                  bc(jdfock)=bc(jdfock)+bc(ixn)*bc(iad3)                5d19s22
     $                 +bc(ixnu)*bc(idad3)                              5d19s22
                 end do                                                 5d9s19
                end do                                                  5d9s19
               end if                                                   4d5s18
               ii=ii+nrow                                                3d16s12
               iiu=iiu+nrow                                             5d19s22
              end do                                                     3d16s12
              i10=1                                                      3d16s12
             end do                                                      3d16s12
            end if                                                       3d16s12
           end do                                                       3d19s12
          end if                                                        3d19s12
         end do                                                         3d19s12
         iarg=idoubb*idoubb*2                                           5d19s22
         call dws_gsumf(bc(ifock),iarg)                                 3d19s12
         if(lpr)then                                                    7d21s21
          write(6,*)('dd fock matrix ')
          call prntm2(bc(ifock),idoubb,idoubb,idoubb)                   9d30s20
          write(6,*)('derivative of dd fock matrix ')
          call prntm2(bc(idfock),idoubb,idoubb,idoubb)                   9d30s20
         end if                                                         3d19s12
         idvec=ibcoff                                                   5d19s22
         ibcoff=idvec+idoubb*idoubb                                     5d19s22
         call enough('dcannon.  3',bc,ibc)
         do i=0,idoubb-1                                                5d19s22
          iad=ifock+i*(idoubb+1)                                        5d19s22
          do j=0,i-1                                                    5d19s22
           jad=ifock+j*(idoubb+1)                                       5d19s22
           ji=j+idoubb*i                                                5d19s22
           ij=i+idoubb*j                                                5d19s22
           val=bc(idfock+ij)/(bc(iad)-bc(jad))                          5d19s22
           bc(idvec+ij)=-val                                             5d19s22
           bc(idvec+ji)=val                                             5d19s22
          end do                                                        5d19s22
          ii=idvec+i*(idoubb+1)                                         5d19s22
          bc(ii)=0d0                                                    5d19s22
         end do                                                         5d19s22
         if(lpr)then                                                    7d25s22
          write(6,*)('derivative of fock vectors: ')                     5d19s22
          call prntm2(bc(idvec),idoubb,idoubb,idoubb)                    5d19s22
          write(6,*)('derivative matrix goes from ')
          call prntm2(bc(ix2),nbasdws(isb),idoubb,nbasdws(isb))          5d19s22
          write(6,*)('to ')
          call dgemm('n','n',nbasdws(isb),idoubb,idoubb,1d0,             5d19s22
     $        bc(morb(isb)),nbasdws(isb),bc(idvec),idoubb,0d0,          5d19s22
     $        bc(ibcoff),nbasdws(isb),                                     5d19s22
     d'dcannon.  2')
          call prntm2(bc(ibcoff),nbasdws(isb),idoubb,nbasdws(isb))
         end if                                                         7d25s22
         call dgemm('n','n',nbasdws(isb),idoubb,idoubb,1d0,             5d19s22
     $        bc(morb(isb)),nbasdws(isb),bc(idvec),idoubb,1d0,          5d19s22
     $        bc(ix2),nbasdws(isb),                                     5d19s22
     d'dcannon.  3')
         if(lpr)then                                                    7d25s22
          call prntm2(bc(ix2),nbasdws(isb),idoubb,nbasdws(isb))          5d19s22
         end if                                                         7d25s22
         ibcoff=ifock                                                   3d19s12
        end if                                                          3d19s12
        if(iacto(isb).gt.1)then                                         4d5s18
         if(lpr)then                                                    7d25s22
          write(6,*)('density matrix: ')
          call prntm2(bc(iden1(isb)),iacto(isb),iacto(isb),iacto(isb))
          write(6,*)('derivative of density matrix: ')
          call prntm2(bc(iden1d(isb)),iacto(isb),iacto(isb),iacto(isb))
         end if
         idvec=ibcoff                                                   5d19s22
         ibcoff=idvec+iacto(isb)*iacto(isb)                             5d19s22
         call enough('dcannon.  4',bc,ibc)
         do i=0,iacto(isb)-1                                                5d19s22
          iad=iden1(isb)+i*(iacto(isb)+1)                               5d19s22
          do j=0,i-1                                                    5d19s22
           jad=iden1(isb)+j*(iacto(isb)+1)                              5d19s22
           ji=j+iacto(isb)*i                                                5d19s22
           ij=i+iacto(isb)*j                                            5d19s22
           bot=bc(iad)-bc(jad)                                          11d30s22
           if(abs(bot).lt.1d-10)then                                    11d30s22
            boti=0d0                                                    11d30s22
           else                                                         11d30s22
            boti=1d0/bot                                                11d30s22
           end if                                                       11d30s22
           val=bc(iden1d(isb)+ij)*boti                                  11d30s22
           bc(idvec+ij)=-val                                             5d19s22
           bc(idvec+ji)=val                                             5d19s22
          end do                                                        5d19s22
          ii=idvec+i*(iacto(isb)+1)                                     5d19s22
          bc(ii)=0d0                                                    5d19s22
         end do                                                         5d19s22
         jx2=ix2+nbasdws(isb)*idoub(isb)                                5d19s22
         mm=morb(isb)+nbasdws(isb)*idoub(isb)                           5d19s22
         if(lpr)then                                                    7d25s22
          write(6,*)('derivative of density vectors: ')                     5d19s22
          call prntm2(bc(idvec),iacto(isb),iacto(isb),iacto(isb))        5d19s22
          write(6,*)('derivative matrix goes from ')
          call prntm2(bc(jx2),nbasdws(isb),iacto(isb),nbasdws(isb))      5d19s22
          write(6,*)('to ')
          call dgemm('n','n',nbasdws(isb),iacto(isb),iacto(isb),1d0,     5d19s22
     $        bc(mm),nbasdws(isb),bc(idvec),iacto(isb),0d0,             5d19s22
     $        bc(ibcoff),nbasdws(isb),                                     5d19s22
     d'dcannon.  4')
          call prntm2(bc(ibcoff),nbasdws(isb),iacto(isb),nbasdws(isb))   5d19s22
         end if                                                         7d25s22
         call dgemm('n','n',nbasdws(isb),iacto(isb),iacto(isb),1d0,     5d19s22
     $        bc(mm),nbasdws(isb),bc(idvec),iacto(isb),1d0,             5d19s22
     $        bc(jx2),nbasdws(isb),                                     5d19s22
     d'dcannon.  5')
         if(lpr)then                                                    7d25s22
          call prntm2(bc(jx2),nbasdws(isb),iacto(isb),nbasdws(isb))      5d19s22
         end if                                                         7d25s22
         ibcoff=idvec                                                   5d19s22
        end if                                                          5d19s22
        if(nvirt(isb).gt.0)then                                         5d20s22
         if(idwsdeb.gt.10)write(6,*)('for vv part'),isb,noc(isb)                     3d19s12
         ifock=ibcoff                                                     3d16s12
         idfock=ifock+nvirt(isb)*nvirt(isb)                             5d20s22
         ibcoff=idfock+nvirt(isb)*nvirt(isb)                            5d20s22
         call enough('dcannon.  5',bc,ibc)
         if(nowpro.eq.0)then                                              3d16s12
          jfock=ifock                                                   3d19s12
          jdfock=idfock                                                 5d20s22
          do i=1,nvirt(isb)                                             3d19s12
           ip=i+noc(isb)                                                3d19s12
           do j=1,nvirt(isb)                                            3d19s12
            jp=j+noc(isb)                                               3d19s12
            iad=ih0mo+nsyo(isb)+jp-1+nbasdws(isb)*(ip-1)                5d20s22
            bc(jfock)=bc(iad)                                           3d19s12
            jfock=jfock+1                                               3d19s12
            iad=idh0+nsyo(isb)+jp-1+nbasdws(isb)*(ip-1)                 5d20s22
            bc(jdfock)=bc(iad)                                          5d20s22
            jdfock=jdfock+1                                              5d20s22
           end do                                                       3d19s12
          end do                                                          3d16s12
         else                                                             3d16s12
          do i=0,nvirt(isb)*nvirt(isb)-1                                3d19s12
           bc(ifock+i)=0d0                                                3d16s12
           bc(idfock+i)=0d0                                             5d20s22
          end do                                                          3d16s12
         end if                                                           3d16s12
         do isa=1,nsymb                                                   3d16s12
          if(noc(isa).gt.0)then                                           3d16s12
           do i=1,nsdlk                                                   3d16s12
            if(isblk(1,i).eq.isa.and.isblk(2,i).eq.isa.and.               3d16s12
     $       isblk(3,i).eq.isb.and.isblk(4,i).eq.isb)then               3d16s12
             call ilimts(nvirt(isb),nvirt(isb),numpro,nowpro,il,ih,i1s, 5d20s22
     $           i1e,i2s,i2e)                                           5d20s22
             i10=i1s                                                     3d16s12
             i1n=nvirt(isb)                                              3d16s12
             ii=jmats(i)-1                                               3d16s12
             nrow=(noc(isa)*(noc(isa)+1))/2                              3d16s12
             do i2=i2s,i2e                                               3d16s12
              if(i2.eq.i2e)i1n=i1e                                       3d16s12
              do i1=i10,i1n                                              3d16s12
               jfock=ifock+i1-1+nvirt(isb)*(i2-1)                        3d19s12
               jdfock=idfock+i1-1+nvirt(isb)*(i2-1)                     5d20s22
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
                 idad3=iden1d(isa)+i4+iacto(isa)*i3                     5d20s22
                 bc(jfock)=bc(jfock)+bc(ixn)*bc(iad3)                   5d10s19
                 bc(jdfock)=bc(jdfock)+bc(ixn)*bc(idad3)                5d20s22
                end do                                                  5d10s19
               end do                                                   5d10s19
               ii=ii+nrow                                                3d16s12
              end do                                                     3d16s12
              i10=1                                                      3d16s12
             end do                                                      3d16s12
            end if                                                        3d16s12
           end do                                                         3d16s12
           do i=1,nsblkkder                                             5d20s22
            if(isblkkder(1,i).eq.isa.and.isblkkder(2,i).eq.isa.and.     5d20s22
     $       isblkkder(3,i).eq.isb.and.isblkkder(4,i).eq.isb)then       5d20s22
             call ilimts(nvirt(isb),nvirt(isb),numpro,nowpro,il,ih,i1s, 5d20s22
     $           i1e,i2s,i2e)                                           5d20s22
             i10=i1s                                                     3d16s12
             i1n=nvirt(isb)                                              3d16s12
             ii=jmatd(i)                                                5d20s22
             nrow=noc(isa)*noc(isa)                                     5d20s22
             do i2=i2s,i2e                                               3d16s12
              if(i2.eq.i2e)i1n=i1e                                       3d16s12
              do i1=i10,i1n                                              3d16s12
               jdfock=idfock+i1-1+nvirt(isb)*(i2-1)                     5d20s22
               do i34=0,idoub(isa)-1                                    5d20s22
                ixn=ii+i34*(noc(isa)+1)                                 5d20s22
                bc(jdfock)=bc(jdfock)+2d0*bc(ixn)                       5d20s22
               end do                                                    3d16s12
               do i3=0,iacto(isa)-1                                     5d10s19
                i3p=i3+idoub(isa)                                       5d20s22
                do i4=0,iacto(isa)-1                                    5d10s19
                 i4p=i4+idoub(isa)                                      5d20s22
                 ixn=ii+i4p+noc(isa)*i3p                                5d20s22
                 iad3=iden1(isa)+i4+iacto(isa)*i3                       5d10s19
                 bc(jdfock)=bc(jdfock)+bc(ixn)*bc(iad3)                 5d20s22
                end do                                                  5d10s19
               end do                                                   5d10s19
               ii=ii+nrow                                                3d16s12
              end do                                                     3d16s12
              i10=1                                                      3d16s12
             end do                                                      3d16s12
            end if                                                        3d16s12
            if(isblkkder(1,i).eq.isa.and.isblkkder(2,i).eq.isa.and.     5d20s22
     $       isblkkder(3,i).eq.isb.and.isblkkder(4,i).eq.isb)then       5d20s22
             call ilimts(nvirt(isb),nvirt(isb),numpro,nowpro,il,ih,i1s, 5d20s22
     $           i1e,i2s,i2e)                                           5d20s22
             i10=i1s                                                     3d16s12
             i1n=nvirt(isb)                                              3d16s12
             nocp=noc(isa)+1                                             3d16s12
             ii=kmatd(i)                                                5d20s22
             nrow=noc(isa)*noc(isa)                                      3d16s12
             do i2=i2s,i2e                                               3d16s12
              if(i2.eq.i2e)i1n=i1e                                       3d16s12
              do i1=i10,i1n                                              3d16s12
               jdfock=idfock+i1-1+nvirt(isb)*(i2-1)                     5d20s22
               do i34=0,idoub(isa)-1                                    5d20s22
                ixn=ii+i34*nocp                                          3d16s12
                bc(jdfock)=bc(jdfock)-bc(ixn)                           5d20s22
               end do                                                    3d16s12
               do i3=0,iacto(isa)-1                                     5d10s19
                i3p=i3+idoub(isa)                                       5d20s22
                do i4=0,iacto(isa)-1                                    5d10s19
                 i4p=i4+idoub(isa)                                      5d20s22
                 iad3=iden1(isa)+i4+iacto(isa)*i3                       5d10s19
                 ixn=ii+i3p+noc(isa)*i4p                                5d10s19
                 bc(jdfock)=bc(jdfock)-bc(ixn)*0.5d0*bc(iad3)           5d20s22
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
             call ilimts(nvirt(isb),nvirt(isb),numpro,nowpro,il,ih,i1s, 5d20s22
     $           i1e,i2s,i2e)                                           5d20s22
             i10=i1s                                                     3d16s12
             i1n=nvirt(isb)                                              3d16s12
             nocp=noc(isa)+1                                             3d16s12
             ii=kmats(i)-nocp                                            3d16s12
             nrow=noc(isa)*noc(isa)                                      3d16s12
             do i2=i2s,i2e                                               3d16s12
              if(i2.eq.i2e)i1n=i1e                                       3d16s12
              do i1=i10,i1n                                              3d16s12
               jfock=ifock+i1-1+nvirt(isb)*(i2-1)                        3d19s12
               jdfock=idfock+i1-1+nvirt(isb)*(i2-1)                     5d20s22
               do i34=1,idoub(isa)                                      4d5s18
                ixn=ii+i34*nocp                                          3d16s12
                bc(jfock)=bc(jfock)-bc(ixn)                              3d16s12
               end do                                                    3d16s12
               do i3=0,iacto(isa)-1                                     5d10s19
                i3p=i3+idoub(isa)+1                                     5d10s19
                do i4=0,iacto(isa)-1                                    5d10s19
                 i4p=i4+idoub(isa)+1                                    5d10s19
                 iad3=iden1(isa)+i4+iacto(isa)*i3                       5d10s19
                 idad3=iden1d(isa)+i4+iacto(isa)*i3                     5d20s22
                 ixn=ii+i3p+noc(isa)*i4p                                5d10s19
                 bc(jfock)=bc(jfock)-bc(ixn)*0.5d0*bc(iad3)             5d10s19
                 bc(jdfock)=bc(jdfock)-bc(ixn)*0.5d0*bc(idad3)          5d20s22
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
         iarg=2*nvirt(isb)*nvirt(isb)                                   5d20s22
         call dws_gsumf(bc(ifock),iarg)                                 3d19s12
         if(idwsdeb.gt.10)then                                          3d19s12
          write(6,*)('fock matrix ')
          call prntm2(bc(ifock),nvirt(isb),nvirt(isb),nvirt(isb))        3d19s12
          write(6,*)('derivative of fock matrix ')
          call prntm2(bc(idfock),nvirt(isb),nvirt(isb),nvirt(isb))      5d20s22
         end if                                                         3d19s12
         idvec=ibcoff                                                   5d19s22
         ibcoff=idvec+nvirt(isb)*nvirt(isb)                             5d19s22
         call enough('dcannon.  6',bc,ibc)
         do i=0,nvirt(isb)-1                                                5d19s22
          iad=ifock+i*(nvirt(isb)+1)                                    5d20s22
          do j=0,i-1                                                    5d19s22
           jad=ifock+j*(nvirt(isb)+1)                                   5d20s22
           ji=j+nvirt(isb)*i                                                5d19s22
           ij=i+nvirt(isb)*j                                            5d19s22
           val=bc(idfock+ij)/(bc(iad)-bc(jad))                          5d20s22
           bc(idvec+ij)=-val                                             5d19s22
           bc(idvec+ji)=val                                             5d19s22
          end do                                                        5d19s22
          ii=idvec+i*(nvirt(isb)+1)                                     5d19s22
          bc(ii)=0d0                                                    5d19s22
         end do                                                         5d19s22
         jx2=ix2+nbasdws(isb)*noc(isb)                                  5d20s22
         mm=morb(isb)+nbasdws(isb)*noc(isb)                             5d20s22
         if(lpr)then                                                    7d25s22
          write(6,*)('derivative of fock vectors: ')                     5d19s22
          call prntm2(bc(idvec),nvirt(isb),nvirt(isb),nvirt(isb))        5d19s22
          write(6,*)('derivative matrix goes from ')
          call prntm2(bc(jx2),nbasdws(isb),nvirt(isb),nbasdws(isb))      5d20s22
          write(6,*)('to ')
         end if                                                         7d25s22
         call dgemm('n','n',nbasdws(isb),nvirt(isb),nvirt(isb),1d0,     5d20s22
     $        bc(mm),nbasdws(isb),bc(idvec),nvirt(isb),1d0,             5d20s22
     $        bc(jx2),nbasdws(isb),                                     5d19s22
     d'dcannon.  6')
         if(lpr)then                                                    7d25s22
          call prntm2(bc(jx2),nbasdws(isb),nvirt(isb),nbasdws(isb))      5d20s22
         end if                                                         7d25s22
         ibcoff=idvec                                                   5d19s22
        end if
        if(lpr)then                                                     7d25s22
         write(6,*)('ending ix2 matrix ')
         call prntm2(bc(ix2),nbasdws(isb),nbasdws(isb),nbasdws(isb))
        end if                                                          7d25s22
        call dgemm('t','n',nbasdws(isb),nbasdws(isb),nbasdws(isb),1d0,  6d24s24
     $       bc(morb(isb)),nbasdws(isb),bc(ix2),nbasdws(isb),           6d24s24
     $       0d0,bc(idarot(isb)),nbasdws(isb),                          6d24s24
     d'dcannon.  2')                                                    6d24s24
        if(lpr)then                                                     6d24s24
         write(6,*)('ending darot matrix ')                             6d24s24
         call prntm2(bc(idarot(isb)),nbasdws(isb),nbasdws(isb),         6d24s24
     $        nbasdws(isb))                                             6d24s24
        end if                                                          6d24s24
       end if                                                           3d19s12
      end do
      return
      end
