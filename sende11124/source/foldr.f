c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine foldr(ih0,i4od,noc,irefo,idoubo,nbasdws,nsymb,isblk,    8d21s23
     $     nsdlk,ihdu,bc,ibc,dcore,idwsdeb)                             8d29s23
      implicit real*8 (a-h,o-z)                                         8d21s23
      include "common.store"                                            8d21s23
      dimension i4od(*),noc(*),irefo(*),idoubo(*),nbasdws(*),           8d21s23
     $     isblk(4,*),ihdu(*)                                           8d21s23
      if(idwsdeb.gt.10)                                                 8d29s23
     $  write(6,*)('fold ints and h0 down to what we need for intcsf')  8d29s23
      jh0=ih0                                                           8d21s23
      kh0=ih0                                                           8d23s23
      do isb=1,nsymb                                                    8d21s23
       if(nbasdws(isb).gt.0)then                                        8d21s23
        if(idwsdeb.gt.10)then                                           8d29s23
         write(6,*)('dh0 for symmetry '),isb                             8d21s23
         call prntm2(bc(jh0),nbasdws(isb),nbasdws(isb),nbasdws(isb))     8d21s23
        end if                                                          8d29s23
        do i=0,idoubo(isb)-1                                            8d28s23
         iad=jh0+i*(nbasdws(isb)+1)                                     8d28s23
         dcore=dcore+2d0*bc(iad)                                        8d28s23
        end do                                                          8d28s23
        itmp=ibcoff                                                     8d21s23
        ibcoff=itmp+irefo(isb)*irefo(isb)                               8d21s23
        call enough('foldr.tmp',bc,ibc)                                 8d21s23
        do i=0,irefo(isb)-1                                             8d21s23
         ip=i+idoubo(isb)                                               8d21s23
         iad1=itmp+irefo(isb)*i                                         8d21s23
         iad2=jh0+idoubo(isb)+nbasdws(isb)*ip                           8d21s23
         do j=0,irefo(isb)-1                                            8d21s23
          bc(iad1+j)=bc(iad2+j)                                         8d21s23
         end do                                                         8d21s23
        end do                                                          8d21s23
        if(idwsdeb.gt.10)then                                           8d29s23
         write(6,*)('irefo part ')
         call prntm2(bc(itmp),irefo(isb),irefo(isb),irefo(isb))
        end if                                                          8d29s23
        do i=0,irefo(isb)*irefo(isb)-1                                  8d21s23
         bc(kh0+i)=bc(itmp+i)                                           8d23s23
        end do                                                          8d21s23
        ibcoff=itmp                                                     8d21s23
        ihdu(isb)=kh0                                                   8d21s23
        kh0=kh0+irefo(isb)*irefo(isb)                                   8d23s23
        jh0=jh0+nbasdws(isb)*nbasdws(isb)                               8d21s23
       end if                                                           8d21s23
      end do                                                            8d21s23
      dcore1=dcore
      if(idwsdeb.gt.10)write(6,*)('dcore after 1 e part '),dcore1       8d29s23
      dcore=0d0
      do is=1,nsdlk                                                     8d21s23
       nrow=noc(isblk(1,is))*noc(isblk(2,is))                           8d21s23
       ncol=noc(isblk(3,is))*noc(isblk(4,is))                           8d21s23
       isw=1                                                            8d21s23
       if(min(nrow,ncol).gt.0)then                                      8d21s23
        if(idwsdeb.gt.10)then                                           8d29s23
         write(6,*)('for integral type '),(isblk(j,is),j=1,4),
     $      i4od(is)
         call prntm2(bc(i4od(is)),nrow,ncol,nrow)                        8d21s23
        end if                                                          8d29s23
        if(isblk(1,is).eq.isblk(2,is))then                              8d21s23
         do id=0,idoubo(isblk(3,is))-1                                  8d28s23
          icol=id*(noc(isblk(3,is))+1)                                  8d28s23
          do idp=0,idoubo(isblk(1,is))-1                                8d28s23
           irow=idp*(noc(isblk(1,is))+1)                                8d28s23
           iad=i4od(is)+irow+nrow*icol                                  8d28s23
           dcore=dcore+2d0*bc(iad)                                      8d28s23
          end do                                                        8d28s23
         end do                                                         8d28s23
         if(idoubo(isblk(1,is)).gt.0)then                               8d21s23
          if(idwsdeb.gt.10)write(6,*)('contribution for h12:')
          do ia=0,irefo(isblk(3,is))-1                                   8d21s23
           iap=ia+idoubo(isblk(3,is))                                   8d21s23
           do ja=0,irefo(isblk(3,is))-1                                 8d21s23
            jap=ja+idoubo(isblk(3,is))                                  8d21s23
            icol=jap+noc(isblk(3,is))*iap                               8d21s23
            do id=0,idoubo(isblk(1,is))-1                               8d21s23
             irow=id+noc(isblk(1,is))*id                                8d21s23
             ii=i4od(is)+irow+nrow*icol                                 8d21s23
             iad=ihdu(isblk(3,is))+ja+irefo(isblk(3,is))*ia             8d21s23
             bc(iad)=bc(iad)+2d0*bc(ii)                                 8d21s23
            end do                                                      8d21s23
           end do                                                       8d21s23
          end do                                                        8d21s23
         else if(isblk(1,is).ne.isblk(3,is).and.                        8d21s23
     $         idoubo(isblk(3,is)).gt.0)then                            8d21s23
          do id=0,idoubo(isblk(3,is))-1                                  8d28s23
           icol=id*(noc(isblk(3,is))+1)                                  8d28s23
           do idp=0,idoubo(isblk(1,is))-1                                8d28s23
            irow=idp*(noc(isblk(1,is))+1)                                8d28s23
            iad=i4od(is)+irow+nrow*icol                                  8d28s23
            dcore=dcore+2d0*bc(iad)                                      8d28s23
           end do                                                        8d28s23
          end do                                                         8d28s23
         end if                                                         8d21s23
        end if                                                          8d21s23
        if(isblk(1,is).eq.isblk(3,is).and.                              8d21s23
     $       isblk(2,is).eq.isblk(4,is))then                            8d21s23
         if(idwsdeb.gt.10)write(6,*)('K type ...'),(isblk(j,is),j=1,4)
         do id=0,idoubo(isblk(2,is))-1                                  8d28s23
          do idp=0,idoubo(isblk(1,is))-1                                8d28s23
           irow=idp+noc(isblk(1,is))*id                                 8d28s23
           iad=i4od(is)+irow*(nrow+1)                                   8d28s23
           dcore=dcore-bc(iad)                                          8d28s23
          end do                                                        8d28s23
         end do                                                         8d28s23
          if(idwsdeb.gt.10)write(6,*)('contribution to h13: ')
          do ia=0,irefo(isblk(2,is))-1                                     8d21s23
           iap=ia+idoubo(isblk(2,is))                                      8d21s23
           do ja=0,irefo(isblk(2,is))-1                                    8d21s23
            jap=ja+idoubo(isblk(2,is))                                     8d21s23
            do id=0,idoubo(isblk(1,is))-1                               8d21s23
             irow=id+noc(isblk(1,is))*jap                               8d21s23
             icol=id+noc(isblk(1,is))*iap                               8d21s23
             ii=i4od(is)+irow+nrow*icol                                 8d21s23
             iad=ihdu(isblk(2,is))+ja+irefo(isblk(2,is))*ia             8d21s23
             bc(iad)=bc(iad)-bc(ii)                                     8d21s23
            end do                                                      8d28s23
            do id=0,irefo(isblk(1,is))-1                                8d25s23
             idp=id+idoubo(isblk(1,is))                                 8d25s23
             irow=idp+noc(isblk(1,is))*jap                               8d21s23
             icol=idp+noc(isblk(1,is))*iap                               8d21s23
             ii=i4od(is)+irow+nrow*icol                                 8d21s23
             iad=ihdu(isblk(2,is))+ja+irefo(isblk(2,is))*ia             8d21s23
             orig=bc(iad)
             bc(iad)=bc(iad)-bc(ii)*0.5d0                               8d25s23
            end do                                                      8d21s23
           end do                                                       8d21s23
          end do                                                        8d21s23
         if(isblk(1,is).ne.isblk(2,is))then                             8d25s23
          if(idwsdeb.gt.10)write(6,*)('contribution to h24: ')
          do id=0,idoubo(isblk(2,is))-1                                  8d28s23
           do idp=0,idoubo(isblk(1,is))-1                                8d28s23
            irow=idp+noc(isblk(1,is))*id                                 8d28s23
            iad=i4od(is)+irow*(nrow+1)                                   8d28s23
            dcore=dcore-bc(iad)                                          8d28s23
           end do                                                        8d28s23
          end do                                                         8d28s23
          do ia=0,irefo(isblk(1,is))-1                                     8d21s23
           iap=ia+idoubo(isblk(1,is))                                      8d21s23
           do ja=0,irefo(isblk(1,is))-1                                    8d21s23
            jap=ja+idoubo(isblk(1,is))                                  8d21s23
            do id=0,idoubo(isblk(2,is))-1                               8d21s23
             irow=jap+noc(isblk(1,is))*id                               8d28s23
             icol=iap+noc(isblk(1,is))*id                               8d28s23
             ii=i4od(is)+irow+nrow*icol                                 8d21s23
             iad=ihdu(isblk(1,is))+ja+irefo(isblk(1,is))*ia             8d21s23
             bc(iad)=bc(iad)-bc(ii)                                     8d21s23
            end do                                                      8d28s23
            do id=0,irefo(isblk(2,is))-1                                8d25s23
             idp=id+idoubo(isblk(2,is))                                 8d25s23
             irow=jap+noc(isblk(1,is))*idp                              8d25s23
             icol=iap+noc(isblk(1,is))*idp                              8d25s23
             ii=i4od(is)+irow+nrow*icol                                 8d21s23
             iad=ihdu(isblk(1,is))+ja+irefo(isblk(1,is))*ia             8d21s23
             orig=bc(iad)
             bc(iad)=bc(iad)-bc(ii)*0.5d0                               8d25s23
            end do                                                      8d21s23
           end do                                                       8d21s23
          end do                                                        8d21s23
         end if                                                         8d21s23
        end if                                                          8d21s23
        if(isblk(1,is).eq.isblk(2,is))then                              8d21s23
         nrown=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2            8d21s23
         ncoln=(irefo(isblk(3,is))*(irefo(isblk(3,is))+1))/2            8d21s23
         itmp=ibcoff                                                    8d21s23
         ibcoff=itmp+nrown*ncoln                                        8d21s23
         call enough('foldr.tmp',bc,ibc)                                8d21s23
         do i4=0,irefo(isblk(4,is))-1                                   8d21s23
          i4p=i4+idoubo(isblk(4,is))                                    8d21s23
          do i3=0,i4                                                    8d21s23
           i3p=i3+idoubo(isblk(4,is))                                   8d21s23
           icol=i3p+noc(isblk(3,is))*i4p                                8d21s23
           jcol=((i4*(i4+1))/2)+i3                                      8d21s23
           do i2=0,irefo(isblk(2,is))-1                                 8d21s23
            i2p=i2+idoubo(isblk(2,is))                                  8d21s23
            do i1=0,i2                                                  8d21s23
             i1p=i1+idoubo(isblk(2,is))                                 8d21s23
             irow=i1p+noc(isblk(1,is))*i2p                              8d21s23
             jrow=((i2*(i2+1))/2)+i1                                    8d21s23
             ii=i4od(is)+irow+nrow*icol                                 8d21s23
             jj=itmp+jrow+nrown*jcol                                    8d21s23
             bc(jj)=bc(ii)                                              8d21s23
            end do                                                      8d21s23
           end do                                                       8d21s23
          end do                                                        8d21s23
         end do                                                         8d21s23
         do i=0,nrown*ncoln-1                                           8d21s23
          bc(i4od(is)+i)=bc(itmp+i)                                     8d21s23
         end do                                                         8d21s23
         ibcoff=itmp                                                    8d21s23
        else                                                            8d21s23
         nrown=irefo(isblk(1,is))*irefo(isblk(2,is))                    8d21s23
         ncoln=irefo(isblk(3,is))*irefo(isblk(4,is))                    8d21s23
         itmp=ibcoff                                                    8d21s23
         ibcoff=itmp+nrown*ncoln                                        8d21s23
         call enough('foldr.tmpb',bc,ibc)                               8d21s23
         do i4=0,irefo(isblk(4,is))-1                                   8d21s23
          i4p=i4+idoubo(isblk(4,is))                                    8d21s23
          do i3=0,irefo(isblk(3,is))-1                                  8d21s23
           i3p=i3+idoubo(isblk(3,is))                                   8d21s23
           icol=i3p+noc(isblk(3,is))*i4p                                8d21s23
           jcol=i3+irefo(isblk(3,is))*i4                                8d21s23
           do i2=0,irefo(isblk(2,is))-1                                 8d21s23
            i2p=i2+idoubo(isblk(2,is))                                  8d21s23
            do i1=0,irefo(isblk(1,is))                                  8d21s23
             i1p=i1+idoubo(isblk(1,is))                                 8d21s23
             irow=i1p+noc(isblk(1,is))*i2p                              8d21s23
             jrow=i1+irefo(isblk(1,is))*i2                              8d21s23
             ii=i4od(is)+irow+nrow*icol                                 8d21s23
             jj=itmp+jrow+nrown*jcol                                    8d21s23
             bc(jj)=bc(ii)                                              8d21s23
            end do                                                      8d21s23
           end do                                                       8d21s23
          end do                                                        8d21s23
         end do                                                         8d21s23
         do i=0,nrown*ncoln-1                                           8d21s23
          bc(i4od(is)+i)=bc(itmp+i)                                     8d21s23
         end do                                                         8d21s23
         ibcoff=itmp                                                    8d21s23
        end if                                                          8d21s23
        if(idwsdeb.gt.10)then                                           8d29s23
         write(6,*)('compacted matrix ')                                 8d21s23
         call prntm2(bc(i4od(is)),nrown,ncoln,nrown)                     8d21s23
        end if                                                          8d29s23
       end if                                                           8d21s23
      end do                                                            8d21s23
      if(idwsdeb.gt.10)then                                             8d29s23
       write(6,*)('modified h0 part ')
       do isb=1,nsymb
        if(irefo(isb).gt.0)then
         write(6,*)('for symmetry = '),isb,ihdu(isb)
         call prntm2(bc(ihdu(isb)),irefo(isb),irefo(isb),irefo(isb))
        end if
       end do
       write(6,*)('dcore after 2e part '),dcore                                      8d28s23
      end if                                                            8d29s23
      dcore=dcore+dcore1                                                8d28s23
      if(idwsdeb.gt.10)write(6,*)('total dcore '),dcore                 8d29s23
      return                                                            8d21s23
      end                                                               8d21s23
