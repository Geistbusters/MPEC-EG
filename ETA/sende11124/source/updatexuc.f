c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine updatexuc(hdd,vdold,nhddiag,nhsdiag,                   6d9s21
     $     vsold,hss,nff1,nff2,ncsf,ncsf2,nsymb,mdon,mdoo,nvirt,        6d9s21
     $     isymmrci,eig,                                                6d9s21
     $     nroot,isto,vstry,nsing,lprt,ngots,ngotd,                     6d9s21
     $     iamconverged,nlocald,shift0,pairs,bc,ibc)                    11d10s22
      implicit real*8 (a-h,o-z)                                         6d22s18
c                                                                       6d22s18
c     davidson update for uncontracted doubles.                           11d10s20
c
      logical lprint,lprt                                               1d25s21
      integer*8 nhsdiag(mdoo+1,nsymb,2),irow,icol,il8,ih8,i1,           4d21s21
     $     iamconverged(*),nhddiag(mdoo+1,nsymb,2)                      6d9s21
      dimension hdd(*),vsold(*),hss(*),vdold(*),                        6d9s21
     $     nff1(mdoo+1,nsymb),ncsf(*),ncsf2(4,*),nvirt(*),eig(*),       6d9s21
     $     nff2(mdoo+1,nsymb),vstry(*),pairs(4,*)                       8d19s21
      include "common.basis"                                            7d10s19
      include "common.input"                                            7d10s19
      include "common.store"
      include "common.print"                                            1d13s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data icall/0/
      save
      i1=1                                                              11d10s20
      icall=icall+1
      ibcoffo=ibcoff
      if(iprtr(11).ne.0)then                                            1d13s20
       lprint=.true.                                                     6d26s18
      else                                                              1d13s20
       lprint=.false.                                                    6d26s18
      end if                                                            1d13s20
      if(lprint)then                                                    6d26s18
       write(6,*)('in updatexuc ')
       write(6,*)('input eigenvalues ')
       call prntm2(eig,1,nroot,1)
       ivstest=ibcoff                                                   1d24s21
       igstest=ivstest+nsing*nroot                                      1d24s21
       ibcoff=igstest+nsing*nroot                                       1d24s21
       call enough('updatexuc.  1',bc,ibc)
       if(nsing.gt.0)then                                               3d19s21
        jvstest=ivstest                                                  1d24s21
        jgstest=igstest                                                  1d24s21
        do isb=1,nsymb                                                   11d10s20
         isbv=multh(isb,isymmrci)                                        11d10s20
         irow=nvirt(isbv)*nroot                                          11d10s20
         irow4=irow                                                      11d10s20
         do ii=mdon+1,mdoo+1                                             11d10s20
          iarg=ii-mdon                                                   11d10s20
          ndblock=ncsf(iarg)*nff1(ii,isb)
          if(ndblock.gt.0)then                                           11d10s20
           icol=ndblock                                                  11d10s20
           icol4=icol                                                     11d10s20
           call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)  11d10s20
           nhere=ih+1-il                                                  11d10s20
           il8=il                                                        11d10s20
           ih8=ih                                                        11d10s20
           ivec=ibcoff                                                   11d10s20
           ibcoff=ivec+irow4*nhere                                       11d10s20
           call ddi_get(bc,ibc,nhsdiag(ii,isb,1),i1,irow,il8,ih8,       11d15s22
     $          bc(ivec))                                               11d15s22
           jvec=ivec                                                     1d24s21
           do if1=1,nff1(ii,isb)                                         1d24s21
            do i=0,ncsf(iarg)-1                                          1d24s21
             do iv=0,nvirt(isbv)-1                                       1d24s21
              do ir=0,nroot-1                                            1d24s21
               iadg=jvstest+nsing*ir                                     1d24s21
               bc(iadg)=bc(jvec+ir)                                      1d24s21
              end do                                                     1d24s21
              jvstest=jvstest+1                                          1d24s21
              jvec=jvec+nroot                                            1d24s21
             end do                                                      1d24s21
            end do                                                       1d24s21
           end do                                                        1d24s21
           call ddi_get(bc,ibc,nhsdiag(ii,isb,2),i1,irow,il8,ih8,       11d15s22
     $          bc(ivec))                                               11d15s22
           jvec=ivec                                                     1d24s21
           do if1=1,nff1(ii,isb)                                         1d24s21
            do i=0,ncsf(iarg)-1                                          1d24s21
             do iv=0,nvirt(isbv)-1                                       1d24s21
              do ir=0,nroot-1                                            1d24s21
               iadg=jgstest+nsing*ir                                     1d24s21
               bc(iadg)=bc(jvec+ir)                                      1d24s21
              end do                                                     1d24s21
              jgstest=jgstest+1                                          1d24s21
              jvec=jvec+nroot                                            1d24s21
             end do                                                      1d24s21
            end do                                                       1d24s21
           end do                                                        1d24s21
           ibcoff=ivec                                                   11d10s20
          end if                                                         11d10s20
         end do                                                          11d10s20
        end do                                                           11d10s20
       end if                                                           3d19s21
       ibcoff=ivstest                                                   10d7s21
      end if                                                            6d26s18
      do ir=1,nroot                                                     8d19s21
       do l=1,4                                                         8d19s21
        pairs(l,ir)=0d0                                                 8d19s21
       end do                                                           8d19s21
      end do                                                            8d19s21
      ioff=1
      ioffd=1                                                           12d11s20
      jhdd=1                                                            6d9s21
      ivdtry=ibcoff                                                     6d9s21
      ibcoff=ivdtry+nlocald                                             6d9s21
      call enough('updatexuc.  2',bc,ibc)
      jvdtry=ivdtry                                                     6d9s21
      do jsb=1,nsymb                                                    12d8s20
       jsbv12=multh(jsb,isymmrci)                                       12d8s20
       nvisv=0                                                          6d9s21
       nvnotv=0                                                         6d9s21
       do jsbv1=1,nsymb                                                 6d9s21
        jsbv2=multh(jsbv1,jsbv12)                                       6d9s21
        if(jsbv2.ge.jsbv1)then                                          6d9s21
         if(jsbv1.eq.jsbv2)then                                         6d9s21
          nvisv=nvisv+nvirt(jsbv1)                                      6d9s21
          nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                         6d9s21
         else                                                           6d9s21
          nvv=nvirt(jsbv1)*nvirt(jsbv2)                                 6d9s21
         end if                                                         6d9s21
         nvnotv=nvnotv+nvv                                              6d9s21
        end if                                                          6d9s21
       end do                                                           6d9s21
       do ii=mdon+1,mdoo+1                                              6d9s21
        if(nff2(ii,jsb).gt.0)then                                       6d9s21
         iarg=ii-mdon                                                    6d9s21
         icol=nff2(ii,jsb)                                              6d9s21
         call ilimts(1,nff2(ii,jsb),mynprocg,mynowprog,il,ih,i1s,i1e,   6d9s21
     $        i2s,i2e)                                                  6d9s21
         nhere=ih+1-il                                                  6d9s21
         nrow=nroot*(nvisv*ncsf2(1,iarg)+nvnotv*ncsf(iarg))             6d9s21
         il8=il                                                         6d9s21
         ih8=ih                                                         6d9s21
         ivecd=ibcoff                                                   6d9s21
         igd=ivecd+nrow*nhere                                           6d9s21
         ibcoff=igd+nrow*nhere                                          6d9s21
         call enough('updatexuc.  3',bc,ibc)
         irow=nrow                                                      6d9s21
         call ddi_get(bc,ibc,nhddiag(ii,jsb,1),i1,irow,il8,ih8,         11d15s22
     $        bc(ivecd))                                                11d15s22
         call ddi_get(bc,ibc,nhddiag(ii,jsb,2),i1,irow,il8,ih8,bc(igd)) 11d15s22
         jvecd=ivecd                                                    6d9s21
         jgd=igd                                                        6d9s21
         do iff=il,ih                                                   6d9s21
          do jsbv1=1,nsymb                                               6d13s21
           jsbv2=multh(jsbv1,jsbv12)                                     6d9s21
           if(jsbv2.ge.jsbv1)then                                        6d13s21
            if(jsbv1.eq.jsbv2)then                                      6d9s21
             do iv=0,nvirt(jsbv1)-1                                     6d9s21
              do i=0,ncsf2(1,iarg)-1                                    6d9s21
               nadd=-1+i*nroot                                          6d9s21
               kgd=jgd+nadd                                             6d9s21
               kvecd=jvecd+nadd                                         6d9s21
               kvdtry=jvdtry+nadd                                       6d9s21
               do ir=1,nroot                                            6d9s21
                resid=bc(kgd+ir)-(eig(ir)+shift0)*bc(kvecd+ir)                   6d9s21
                bot=hdd(jhdd)-(eig(ir)+shift0)                                   6d9s21
                pairs(1,ir)=pairs(1,ir)+bc(kgd+ir)*bc(kvecd+ir)         8d19s21
                bot=bot+1d-14                                           6d9s21
                update=-resid/bot                                       6d9s21
                bc(kvdtry+ir)=update                                    6d9s21
  333           format(4es15.7,7i8)
               end do                                                   6d9s21
              end do                                                    6d9s21
              jhdd=jhdd+1                                               6d9s21
              nadd=nroot*ncsf2(1,iarg)                                  6d9s21
              jgd=jgd+nadd                                              6d9s21
              jvecd=jvecd+nadd                                          6d9s21
              jvdtry=jvdtry+nadd                                        6d9s21
             end do                                                     6d9s21
             nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                      6d9s21
            else                                                        6d9s21
             nvv=nvirt(jsbv1)*nvirt(jsbv2)                              6d9s21
            end if                                                      6d9s21
            do ivv=0,nvv-1                                              6d9s21
             l=1                                                        8d19s21
             lsum=ncsf2(1,iarg)                                         8d19s21
             do i=0,ncsf(iarg)-1                                        6d9s21
              if(i.gt.lsum)then                                         8d19s21
               l=l+1                                                    8d19s21
               lsum=lsum+ncsf2(l,iarg)                                  8d19s21
              end if                                                    8d19s21
              nadd=-1+i*nroot                                           6d9s21
              kgd=jgd+nadd                                              6d9s21
              kvecd=jvecd+nadd                                          6d9s21
              kvdtry=jvdtry+nadd                                        6d9s21
              do ir=1,nroot                                              6d9s21
               resid=bc(kgd+ir)-(eig(ir)+shift0)*bc(kvecd+ir)                    6d9s21
               pairs(l,ir)=pairs(l,ir)+bc(kgd+ir)*bc(kvecd+ir)          8d19s21
               bot=hdd(jhdd)-(eig(ir)+shift0)                                    6d9s21
               bot=bot+1d-14                                            6d9s21
               update=-resid/bot                                        6d9s21
               bc(kvdtry+ir)=update                                     6d9s21
              end do                                                     6d9s21
             end do                                                     6d9s21
             jhdd=jhdd+1                                                6d9s21
             nadd=nroot*ncsf(iarg)                                      6d9s21
             jgd=jgd+nadd                                               6d9s21
             jvecd=jvecd+nadd                                           6d9s21
             jvdtry=jvdtry+nadd                                         6d9s21
            end do                                                      6d9s21
           end if                                                        6d13s21
          end do                                                         6d9s21
         end do                                                         6d13s21
         ibcoff=ivecd                                                   6d9s21
        end if                                                          6d9s21
       end do                                                           6d9s21
      end do                                                            6d9s21
      jvsold=1                                                          11d10s20
      jhss=1                                                            11d10s20
      if(nsing.gt.0)then                                                3d19s21
       do isb=1,nsymb                                                    11d10s20
        isbv=multh(isb,isymmrci)                                         11d10s20
        irow=nvirt(isbv)*nroot                                           11d10s20
        irow4=irow                                                       11d10s20
        do ii=mdon+1,mdoo+1                                              11d10s20
         iarg=ii-mdon                                                    11d10s20
         ndblock=ncsf(iarg)*nff1(ii,isb)
         if(ndblock.gt.0)then                                            11d10s20
          icol=ndblock                                                   11d10s20
          icol4=icol                                                      11d10s20
          call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)   11d10s20
          nhere=ih+1-il                                                   11d10s20
          il8=il                                                         11d10s20
          ih8=ih                                                         11d10s20
          ivecs=ibcoff                                                    11d10s20
          igs=ivecs+irow4*nhere                                          11d10s20
          ibcoff=igs+nvirt(isbv)*nhere                                   11d10s20
          call enough('updatexuc.  4',bc,ibc)
          call ddi_get(bc,ibc,nhsdiag(ii,isb,1),i1,irow,il8,ih8,        11d15s22
     $         bc(ivecs))                                               11d15s22
          call ddi_get(bc,ibc,nhsdiag(ii,isb,2),i1,irow,il8,ih8,bc(igs))11d15s22
          ic=0                                                           11d10s20
          jvecs=ivecs                                                    11d10s20
          jgs=igs                                                        11d10s20
          jvsold0=jvsold
          do iff=1,nff1(ii,isb)                                          11d10s20
           khss=jhss                                                     11d10s20
           do m=1,ncsf(iarg)                                             11d10s20
            ic=ic+1                                                      11d10s20
            if(ic.ge.il.and.ic.le.ih)then                                11d10s20
             khss=jhss                                                   12d14s20
             do iv=1,nvirt(isbv)                                         11d10s20
              do ir=1,nroot                                                 11d10s20
               resid=bc(jgs)-(eig(ir)+shift0)*bc(jvecs)                           11d10s20
               bot=hss(khss)-(eig(ir)+shift0)                                     11d10s20
               bot=bot+1d-14                                              11d10s20
               update=-resid/bot                                         11d10s20
               vstry(jvsold)=update                                      11d10s20
               jgs=jgs+1                                                 11d10s20
               jvecs=jvecs+1                                             11d10s20
               jvsold=jvsold+1                                           11d10s20
              end do                                                     11d10s20
              khss=khss+1                                                12d14s20
             end do                                                       11d10s20
            end if                                                       11d10s20
           end do                                                        11d10s20
           jhss=khss                                                     11d10s20
          end do                                                         11d10s20
          ibcoff=ivecs                                                   11d10s20
         end if                                                          11d10s20
        end do                                                           11d10s20
       end do                                                            11d10s20
      end if                                                            3d19s21
c
c     print out what we have so far
c
c
c     modified gram-schmidtz orthogonalization of previous (old) vectors
c
      ngot=max(ngots,ngotd)                                             4d14s21
      jvdtry=ivdtry-1                                                   6d9s21
      do ir=1,nroot                                                     11d10s20
       if(lprint)write(6,*)('for new root no.: '),ir
       do iro=1,ngot                                                    1d27s21
        if(lprint)write(6,*)('for old root no.: '),iro
        sum2=0d0                                                        11d10s20
        if(iro.le.ngotd)then                                            6d9s21
         ioff=0                                                          11d10s20
         do jsb=1,nsymb                                                    12d8s20
          jsbv12=multh(jsb,isymmrci)                                       12d8s20
          do ii=mdon+1,mdoo+1                                           6d9s21
           if(nff2(ii,jsb).gt.0)then                                    6d9s21
            iarg=ii-mdon                                                6d9s21
            call ilimts(1,nff2(ii,jsb),mynprocg,mynowprog,il,ih,i1s,i1e,6d9s21
     $        i2s,i2e)                                                  6d9s21
            do iff=il,ih                                                6d9s21
             do jsbv1=1,nsymb                                                 12d8s20
              jsbv2=multh(jsbv1,jsbv12)
              if(jsbv2.ge.jsbv1)then
               if(jsbv1.eq.jsbv2)then                                    6d9s21
                do iv=0,nvirt(jsbv1)-1                                   6d9s21
                 koff=ioff+iro                                           6d9s21
                 joff=jvdtry+ioff+ir                                     6d9s21
                 do i=0,ncsf2(1,iarg)-1                                  6d9s21
                  sum2=sum2+bc(joff+i*nroot)*vdold(koff+i*nroot)         6d9s21
                 end do                                                  6d9s21
                 ioff=ioff+nroot*ncsf2(1,iarg)                           6d9s21
                end do                                                   6d9s21
                nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                    6d9s21
               else                                                      6d9s21
                nvv=nvirt(jsbv1)*nvirt(jsbv2)                            6d9s21
               end if                                                   6d9s21
               do ivv=0,nvv-1                                            6d9s21
                koff=ioff+iro                                            6d9s21
                joff=jvdtry+ioff+ir                                      6d9s21
                do i=0,ncsf(iarg)-1                                      6d9s21
                 sum2=sum2+bc(joff+i*nroot)*vdold(koff+i*nroot)          6d9s21
                end do                                                   6d9s21
                ioff=ioff+nroot*ncsf(iarg)                               6d9s21
               end do                                                    6d9s21
              end if                                                     6d9s21
             end do                                                      6d9s21
            end do                                                      6d9s21
           end if                                                       6d9s21
          end do                                                        6d9s21
         end do                                                         6d9s21
        end if                                                          6d9s21
        sum1=0d0                                                        11d10s20
        if(nsing.gt.0)then                                              3d19s21
         jvsold=0                                                        11d10s20
         jvstry=0                                                        11d10s20
         do isb=1,nsymb                                                  11d10s20
          isbv=multh(isb,isymmrci)                                       11d10s20
          do ii=mdon+1,mdoo+1                                            11d10s20
           iarg=ii-mdon                                                  11d10s20
           ndblock=ncsf(iarg)*nff1(ii,isb)                               11d10s20
           if(ndblock.gt.0)then                                          11d10s20
            icol4=ndblock                                                11d10s20
            call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,   3d19s21
     $           i2e)                                                   3d19s21
            ic=0                                                           11d10s20
            do iff=1,nff1(ii,isb)                                        11d10s20
             do m=1,ncsf(iarg)                                            11d10s20
              ic=ic+1                                                    11d10s20
              if(ic.ge.il.and.ic.le.ih)then                              11d10s20
               jvsoldu=jvsold+iro                                        12d14s20
               jvstryu=jvstry+ir                                         12d14s20
               do iv=1,nvirt(isbv)                                       12d14s20
                sum1=sum1+vsold(jvsoldu)*vstry(jvstryu)                  12d14s20
                jvsoldu=jvsoldu+nroot                                    12d14s20
                jvstryu=jvstryu+nroot                                    12d14s20
               end do                                                    12d14s20
               jvsold=jvsold+nvirt(isbv)*nroot                           12d14s20
               jvstry=jvstry+nvirt(isbv)*nroot                           12d14s20
              end if                                                     11d10s20
             end do                                                      11d10s20
            end do                                                       11d10s20
           end if                                                        11d10s20
          end do                                                         11d10s20
         end do                                                          11d10s20
         if(lprint)write(6,*)('my part of dot for singles: '),sum1
        end if                                                          3d19s21
        dot=sum1+sum2                                                   11d10s20
        call dws_gsumf(dot,1)                                           6d9s21
        if(lprint)then                                                    1d25s21
         write(6,*)('total dot: '),dot                                   11d10s20
         write(6,*)('subtracting off ...')
        end if                                                          1d25s21
        if(iro.le.ngotd)then                                            6d9s21
         ioff=0                                                          11d10s20
         do jsb=1,nsymb                                                    12d8s20
          jsbv12=multh(jsb,isymmrci)                                       12d8s20
          do ii=mdon+1,mdoo+1                                           6d9s21
           if(nff2(ii,jsb).gt.0)then                                    6d9s21
            iarg=ii-mdon                                                6d9s21
            call ilimts(1,nff2(ii,jsb),mynprocg,mynowprog,il,ih,i1s,i1e,6d9s21
     $        i2s,i2e)                                                  6d9s21
            do iff=il,ih                                                6d9s21
             do jsbv1=1,nsymb                                                 12d8s20
              jsbv2=multh(jsbv1,jsbv12)
              if(jsbv2.ge.jsbv1)then
               if(jsbv1.eq.jsbv2)then                                    6d9s21
                do iv=0,nvirt(jsbv1)-1                                   6d9s21
                 koff=ioff+iro                                           6d9s21
                 joff=jvdtry+ioff+ir                                     6d9s21
                 do i=0,ncsf2(1,iarg)-1                                  6d9s21
                  bc(joff+i*nroot)=bc(joff+i*nroot)                     6d9s21
     $                 -dot*vdold(koff+i*nroot)                         6d9s21
                 end do                                                  6d9s21
                 ioff=ioff+nroot*ncsf2(1,iarg)                           6d9s21
                end do                                                   6d9s21
                nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                    6d9s21
               else                                                      6d9s21
                nvv=nvirt(jsbv1)*nvirt(jsbv2)                            6d9s21
               end if                                                   6d9s21
               do ivv=0,nvv-1                                            6d9s21
                koff=ioff+iro                                            6d9s21
                joff=jvdtry+ioff+ir                                      6d9s21
                do i=0,ncsf(iarg)-1                                      6d9s21
                 bc(joff+i*nroot)=bc(joff+i*nroot)                      6d9s21
     $                -dot*vdold(koff+i*nroot)                          6d9s21
                end do                                                   6d9s21
                ioff=ioff+nroot*ncsf(iarg)                               6d9s21
               end do                                                    6d9s21
              end if                                                     6d9s21
             end do                                                      6d9s21
            end do                                                      6d9s21
           end if                                                       6d9s21
          end do                                                        6d9s21
         end do                                                          11d10s20
        end if                                                          2d17s21
        if(nsing.gt.0)then                                              3d19s21
         jvsold=0                                                        11d10s20
         jvstry=0                                                        11d10s20
         do isb=1,nsymb                                                  11d10s20
          isbv=multh(isb,isymmrci)                                       11d10s20
          do ii=mdon+1,mdoo+1                                            11d10s20
           iarg=ii-mdon                                                  11d10s20
           ndblock=ncsf(iarg)*nff1(ii,isb)                               11d10s20
           if(ndblock.gt.0)then                                          11d10s20
            icol4=ndblock                                                11d10s20
            call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,   3d19s21
     $           i2e)                                                   3d19s21
            ic=0                                                           11d10s20
            do iff=1,nff1(ii,isb)                                        11d10s20
             do m=1,ncsf(iarg)                                            11d10s20
              ic=ic+1                                                    11d10s20
              if(ic.ge.il.and.ic.le.ih)then                              11d10s20
               jvsoldu=jvsold+iro                                        12d14s20
               jvstryu=jvstry+ir                                         12d14s20
               do iv=1,nvirt(isbv)                                       12d14s20
                vstry(jvstryu)=vstry(jvstryu)-dot*vsold(jvsoldu)         12d14s20
                jvstryu=jvstryu+nroot                                    12d14s20
                jvsoldu=jvsoldu+nroot                                    12d14s20
               end do                                                    12d14s20
               jvsold=jvsold+nvirt(isbv)*nroot                           12d14s20
               jvstry=jvstry+nvirt(isbv)*nroot                           12d14s20
              end if                                                     11d10s20
             end do                                                      11d10s20
            end do                                                       11d10s20
           end if                                                        11d10s20
          end do                                                         11d10s20
         end do                                                          11d10s20
        end if                                                          3d19s21
       end do                                                           11d10s20
      end do                                                            11d10s20
c
c     now modified gram-schmidt among the new vectors
c
      if(lprint)write(6,*)('now among new vectors ')                      1d25s21
      isto=0                                                             3d30s21
      do ir=1,nroot                                                     11d10s20
       if(lprint)write(6,*)('for new root no.: '),ir                      1d25s21
       do iro=1,isto                                                     3d30s21
        if(lprint)write(6,*)('for new root no.: '),iro                    1d25s21
        sum2=0d0                                                        11d10s20
        ioff=0                                                          11d10s20
        do jsb=1,nsymb                                                    12d8s20
         jsbv12=multh(jsb,isymmrci)                                       12d8s20
         do ii=mdon+1,mdoo+1                                            6d9s21
          if(nff2(ii,jsb).gt.0)then                                     6d9s21
           iarg=ii-mdon                                                 6d9s21
           call ilimts(1,nff2(ii,jsb),mynprocg,mynowprog,il,ih,i1s,i1e, 6d9s21
     $       i2s,i2e)                                                   6d9s21
           do iff=il,ih                                                 6d9s21
            do jsbv1=1,nsymb                                                 12d8s20
             jsbv2=multh(jsbv1,jsbv12)
             if(jsbv2.ge.jsbv1)then
              if(jsbv1.eq.jsbv2)then                                    6d9s21
               do iv=0,nvirt(jsbv1)-1                                   6d9s21
                koff=jvdtry+ioff+iro                                           6d9s21
                joff=jvdtry+ioff+ir                                     6d9s21
                do i=0,ncsf2(1,iarg)-1                                  6d9s21
                 sum2=sum2+bc(joff+i*nroot)*bc(koff+i*nroot)            6d9s21
                end do                                                  6d9s21
                ioff=ioff+nroot*ncsf2(1,iarg)                           6d9s21
               end do                                                   6d9s21
               nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                    6d9s21
              else                                                      6d9s21
               nvv=nvirt(jsbv1)*nvirt(jsbv2)                            6d9s21
              end if                                                    6d9s21
              do ivv=0,nvv-1                                            6d9s21
               koff=jvdtry+ioff+iro                                     6d9s21
               joff=jvdtry+ioff+ir                                      6d9s21
               do i=0,ncsf(iarg)-1                                      6d9s21
                sum2=sum2+bc(joff+i*nroot)*bc(koff+i*nroot)             6d9s21
               end do                                                   6d9s21
               ioff=ioff+nroot*ncsf(iarg)                               6d9s21
              end do                                                    6d9s21
             end if                                                     6d9s21
            end do                                                      6d9s21
           end do                                                       6d9s21
          end if                                                        6d9s21
         end do                                                         6d9s21
        end do                                                          6d9s21
        sum1=0d0                                                        11d10s20
        if(nsing.gt.0)then                                              3d19s21
         jvsold=0                                                        11d10s20
         jvstry=0                                                        11d10s20
         do isb=1,nsymb                                                  11d10s20
          isbv=multh(isb,isymmrci)                                       11d10s20
          do ii=mdon+1,mdoo+1                                            11d10s20
           iarg=ii-mdon                                                  11d10s20
           ndblock=ncsf(iarg)*nff1(ii,isb)                               11d10s20
           if(ndblock.gt.0)then                                          11d10s20
            icol4=ndblock                                                11d10s20
            call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,   3d19s21
     $           i2e)                                                   3d19s21
            ic=0                                                           11d10s20
            do iff=1,nff1(ii,isb)                                        11d10s20
             do m=1,ncsf(iarg)                                            11d10s20
              ic=ic+1                                                    11d10s20
              if(ic.ge.il.and.ic.le.ih)then                              11d10s20
               jvsoldu=jvsold+iro                                         12d14s20
               jvstryu=jvstry+ir                                          12d14s20
               do iv=1,nvirt(isbv)                                        12d14s20
                sum1=sum1+vstry(jvsoldu)*vstry(jvstryu)                   12d14s20
                jvsoldu=jvsoldu+nroot                                     12d14s20
                jvstryu=jvstryu+nroot                                     12d14s20
               end do                                                     12d14s20
               jvsold=jvsold+nvirt(isbv)*nroot                           12d14s20
               jvstry=jvstry+nvirt(isbv)*nroot                           12d14s20
              end if                                                     11d10s20
             end do                                                      11d10s20
            end do                                                       11d10s20
           end if                                                        11d10s20
          end do                                                         11d10s20
         end do                                                          11d10s20
         if(lprint)write(6,*)('my part of dot for singles: '),sum1         1d25s21
        end if                                                          3d19s21
        dot=sum1+sum2                                                   11d10s20
        call dws_gsumf(dot,1)
        if(lprint)then                                                    1d25s21
         write(6,*)('total dot: '),dot                                   11d10s20
         write(6,*)('project out ...')
        end if                                                          1d25s21
        ioff=0                                                          11d10s20
        do jsb=1,nsymb                                                    12d8s20
         jsbv12=multh(jsb,isymmrci)                                       12d8s20
         do ii=mdon+1,mdoo+1                                            6d9s21
          if(nff2(ii,jsb).gt.0)then                                     6d9s21
           iarg=ii-mdon                                                 6d9s21
           call ilimts(1,nff2(ii,jsb),mynprocg,mynowprog,il,ih,i1s,i1e, 6d9s21
     $       i2s,i2e)                                                   6d9s21
           do iff=il,ih                                                 6d9s21
            do jsbv1=1,nsymb                                                 12d8s20
             jsbv2=multh(jsbv1,jsbv12)
             if(jsbv2.ge.jsbv1)then
              if(jsbv1.eq.jsbv2)then                                    6d9s21
               do iv=0,nvirt(jsbv1)-1                                   6d9s21
                koff=jvdtry+ioff+iro                                           6d9s21
                joff=jvdtry+ioff+ir                                     6d9s21
                do i=0,ncsf2(1,iarg)-1                                  6d9s21
                 bc(joff+i*nroot)=bc(joff+i*nroot)-dot*bc(koff+i*nroot) 6d9s21
                end do                                                  6d9s21
                ioff=ioff+nroot*ncsf2(1,iarg)                           6d9s21
               end do                                                   6d9s21
               nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                    6d9s21
              else                                                      6d9s21
               nvv=nvirt(jsbv1)*nvirt(jsbv2)                            6d9s21
              end if                                                    6d9s21
              do ivv=0,nvv-1                                            6d9s21
               koff=jvdtry+ioff+iro                                     6d9s21
               joff=jvdtry+ioff+ir                                      6d9s21
               do i=0,ncsf(iarg)-1                                      6d9s21
                bc(joff+i*nroot)=bc(joff+i*nroot)-dot*bc(koff+i*nroot)  6d9s21
               end do                                                   6d9s21
               ioff=ioff+nroot*ncsf(iarg)                               6d9s21
              end do                                                    6d9s21
             end if                                                     6d9s21
            end do                                                      6d9s21
           end do                                                       6d9s21
          end if                                                        6d9s21
         end do                                                         6d9s21
        end do                                                          6d9s21
        if(nsing.gt.0)then                                              3d19s21
         jvsold=0                                                        11d10s20
         jvstry=0                                                        11d10s20
         do isb=1,nsymb                                                  11d10s20
          isbv=multh(isb,isymmrci)                                       11d10s20
          do ii=mdon+1,mdoo+1                                            11d10s20
           iarg=ii-mdon                                                  11d10s20
           ndblock=ncsf(iarg)*nff1(ii,isb)                               11d10s20
           if(ndblock.gt.0)then                                          11d10s20
            icol4=ndblock                                                11d10s20
            call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,   3d19s21
     $           i2e)                                                   3d19s21
            ic=0                                                           11d10s20
            do iff=1,nff1(ii,isb)                                        11d10s20
             do m=1,ncsf(iarg)                                            11d10s20
              ic=ic+1                                                    11d10s20
              if(ic.ge.il.and.ic.le.ih)then                              11d10s20
               jvsoldu=jvsold+iro                                        12d14s20
               jvstryu=jvstry+ir                                         12d14s20
               do iv=1,nvirt(isbv)                                       11d10s20
                vstry(jvstryu)=vstry(jvstryu)-dot*vstry(jvsoldu)         12d14s20
                jvstryu=jvstryu+nroot                                    12d14s20
                jvsoldu=jvsoldu+nroot                                    12d14s20
               end do                                                    11d10s20
               jvsold=jvsold+nvirt(isbv)*nroot                           12d14s20
               jvstry=jvstry+nvirt(isbv)*nroot                           12d14s20
              end if                                                     11d10s20
             end do                                                      11d10s20
            end do                                                       11d10s20
           end if                                                        11d10s20
          end do                                                         11d10s20
         end do                                                          11d10s20
        end if                                                          3d19s21
       end do                                                           11d10s20
c
c     normalize                                                         11d10s20
c
       if(lprint)write(6,*)('normalize')                                  1d25s21
       sum2=0d0                                                         11d10s20
       ioff=0                                                           11d10s20
       do jsb=1,nsymb                                                     12d8s20
        jsbv12=multh(jsb,isymmrci)                                        12d8s20
        do ii=mdon+1,mdoo+1                                             6d9s21
         if(nff2(ii,jsb).gt.0)then                                      6d9s21
          iarg=ii-mdon                                                  6d9s21
          call ilimts(1,nff2(ii,jsb),mynprocg,mynowprog,il,ih,i1s,i1e,  6d9s21
     $      i2s,i2e)                                                    6d9s21
          do iff=il,ih                                                  6d9s21
           do jsbv1=1,nsymb                                                  12d8s20
            jsbv2=multh(jsbv1,jsbv12)
            if(jsbv2.ge.jsbv1)then
             if(jsbv1.eq.jsbv2)then                                     6d9s21
              do iv=0,nvirt(jsbv1)-1                                    6d9s21
                      joff=jvdtry+ioff+ir                                      6d9s21
               do i=0,ncsf2(1,iarg)-1                                   6d9s21
                sum2=sum2+bc(joff+i*nroot)*bc(joff+i*nroot)             6d9s21
               end do                                                   6d9s21
               ioff=ioff+nroot*ncsf2(1,iarg)                            6d9s21
              end do                                                    6d9s21
              nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                     6d9s21
             else                                                       6d9s21
              nvv=nvirt(jsbv1)*nvirt(jsbv2)                             6d9s21
             end if                                                     6d9s21
             do ivv=0,nvv-1                                             6d9s21
              joff=jvdtry+ioff+ir                                       6d9s21
              do i=0,ncsf(iarg)-1                                       6d9s21
               sum2=sum2+bc(joff+i*nroot)*bc(joff+i*nroot)              6d9s21
              end do                                                    6d9s21
              ioff=ioff+nroot*ncsf(iarg)                                6d9s21
             end do                                                     6d9s21
            end if                                                      6d9s21
           end do                                                       6d9s21
          end do                                                        6d9s21
         end if                                                         6d9s21
        end do                                                          6d9s21
       end do                                                           6d9s21
       sum1=0d0                                                         11d10s20
       if(nsing.gt.0)then                                               3d19s21
        jvstry=0                                                         11d10s20
        do isb=1,nsymb                                                   11d10s20
         isbv=multh(isb,isymmrci)                                        11d10s20
         do ii=mdon+1,mdoo+1                                             11d10s20
          iarg=ii-mdon                                                   11d10s20
          ndblock=ncsf(iarg)*nff1(ii,isb)                                11d10s20
          if(ndblock.gt.0)then                                           11d10s20
           icol4=ndblock                                                 11d10s20
           call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)   11d10s20
           ic=0                                                           11d10s20
           do iff=1,nff1(ii,isb)                                         11d10s20
            do m=1,ncsf(iarg)                                            11d10s20
             ic=ic+1                                                     11d10s20
             if(ic.ge.il.and.ic.le.ih)then                               11d10s20
              jvstryu=jvstry+ir                                          12d14s20
              do iv=1,nvirt(isbv)                                        11d10s20
               sum1=sum1+vstry(jvstryu)**2                               12d14s20
               jvstryu=jvstryu+nroot                                     12d14s20
              end do                                                     11d10s20
              jvstry=jvstry+nvirt(isbv)*nroot                            12d14s20
             end if                                                      11d10s20
            end do                                                       11d10s20
           end do                                                        11d10s20
          end if                                                         11d10s20
         end do                                                          11d10s20
        end do                                                           11d10s20
        if(lprint)write(6,*)('my part of dot for singles: '),sum1          1d25s21
       end if                                                           3d19s21
       dot=sum1+sum2                                                    11d10s20
       call dws_gsumf(dot,1)                                            6d9s21
       xnorm=1d0/sqrt(dot)                                              11d10s20
       if(lprint)then                                                     1d25s21
        write(6,*)('for root '),ir
        write(6,*)('total dot: '),dot                                    11d10s20
        write(6,*)('xnorm = '),xnorm
       end if                                                           1d25s21
       if(iamconverged(ir).eq.0)then                                    4d21s21
        isto=isto+1                                                       3d30s21
        sum2=0d0                                                        11d10s20
        ioff=0                                                          11d10s20
        do jsb=1,nsymb                                                    12d8s20
         jsbv12=multh(jsb,isymmrci)                                       12d8s20
         do ii=mdon+1,mdoo+1                                            6d9s21
          if(nff2(ii,jsb).gt.0)then                                     6d9s21
           iarg=ii-mdon                                                 6d9s21
           call ilimts(1,nff2(ii,jsb),mynprocg,mynowprog,il,ih,i1s,i1e, 6d9s21
     $       i2s,i2e)                                                   6d9s21
           do iff=il,ih                                                 6d9s21
            do jsbv1=1,nsymb                                                 12d8s20
             jsbv2=multh(jsbv1,jsbv12)
             if(jsbv2.ge.jsbv1)then
              if(jsbv1.eq.jsbv2)then                                    6d9s21
               do iv=0,nvirt(jsbv1)-1                                   6d9s21
                koff=jvdtry+ioff+isto                                   6d9s21
                joff=jvdtry+ioff+ir                                     6d9s21
                do i=0,ncsf2(1,iarg)-1                                  6d9s21
                 bc(koff+i*nroot)=bc(joff+i*nroot)*xnorm                6d9s21
                end do                                                  6d9s21
                ioff=ioff+nroot*ncsf2(1,iarg)                           6d9s21
               end do                                                   6d9s21
               nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                    6d9s21
              else                                                      6d9s21
               nvv=nvirt(jsbv1)*nvirt(jsbv2)                            6d9s21
              end if                                                    6d9s21
              do ivv=0,nvv-1                                            6d9s21
               koff=jvdtry+ioff+isto                                    6d9s21
               joff=jvdtry+ioff+ir                                      6d9s21
               do i=0,ncsf(iarg)-1                                      6d9s21
                bc(koff+i*nroot)=bc(joff+i*nroot)*xnorm                 6d9s21
               end do                                                   6d9s21
               ioff=ioff+nroot*ncsf(iarg)                               6d9s21
              end do                                                    6d9s21
             end if                                                     6d9s21
            end do                                                      6d9s21
           end do                                                       6d9s21
          end if                                                        6d9s21
         end do                                                         6d9s21
        end do                                                          6d9s21
        if(nsing.gt.0)then                                               3d19s21
         jvstry=0                                                         11d10s20
         do isb=1,nsymb                                                   11d10s20
          isbv=multh(isb,isymmrci)                                        11d10s20
          irow4=nvirt(isbv)*nroot
          do ii=mdon+1,mdoo+1                                             11d10s20
           iarg=ii-mdon                                                   11d10s20
           ndblock=ncsf(iarg)*nff1(ii,isb)                                11d10s20
           if(ndblock.gt.0)then                                           11d10s20
            icol4=ndblock                                                 11d10s20
            call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,   3d30s21
     $           i2e)                                                   3d30s21
            ic=0                                                           11d10s20
            ncol=ih+1-il
            jvstry0=jvstry+1                                              11d10s20
            do iff=1,nff1(ii,isb)                                         11d10s20
             do m=1,ncsf(iarg)                                            11d10s20
              ic=ic+1                                                     11d10s20
              if(ic.ge.il.and.ic.le.ih)then                               11d10s20
               jvstryu=jvstry+ir                                          12d14s20
               jvstrysu=jvstry+isto                                      3d30s21
               do iv=1,nvirt(isbv)                                        11d10s20
                vstry(jvstrysu)=vstry(jvstryu)*xnorm                       12d14s20
                jvstryu=jvstryu+nroot                                     12d14s20
                jvstrysu=jvstrysu+nroot                                 3d30s21
               end do                                                     11d10s20
               jvstry=jvstry+nvirt(isbv)*nroot                            12d14s20
              end if                                                      11d10s20
             end do                                                       11d10s20
            end do                                                        11d10s20
           end if                                                         11d10s20
          end do                                                          11d10s20
         end do                                                           11d10s20
        end if                                                           3d19s21
       end if                                                           3d30s21
      end do                                                            11d10s20
      if(isto.ne.nroot)then                                             3d30s21
       if(lprint)then                                                     7d7s21
        write(6,*)('total no. of vectors kept = '),isto
        write(6,*)('we need to compress vector storage')                 3d30s21
       end if                                                           7d7s21
       ioff=0                                                           11d10s20
       ioffs=0                                                          6d9s21
       do jsb=1,nsymb                                                     12d8s20
        jsbv12=multh(jsb,isymmrci)                                        12d8s20
        do ii=mdon+1,mdoo+1                                             6d9s21
         if(nff2(ii,jsb).gt.0)then                                      6d9s21
          iarg=ii-mdon                                                  6d9s21
          call ilimts(1,nff2(ii,jsb),mynprocg,mynowprog,il,ih,i1s,i1e,  6d9s21
     $      i2s,i2e)                                                    6d9s21
          do iff=il,ih                                                  6d9s21
           do jsbv1=1,nsymb                                                  12d8s20
            jsbv2=multh(jsbv1,jsbv12)
            if(jsbv2.ge.jsbv1)then
             if(jsbv1.eq.jsbv2)then                                     6d9s21
              do iv=0,nvirt(jsbv1)-1                                    6d9s21
               do i=0,ncsf2(1,iarg)-1                                   6d9s21
                koff=jvdtry+ioffs+isto*i                                6d9s21
                joff=jvdtry+ioff+nroot*i                                6d9s21
                do ir=1,isto                                            6d9s21
                 bc(koff+ir)=bc(joff+ir)                                6d9s21
                end do                                                  6d9s21
               end do                                                   6d9s21
               ioff=ioff+nroot*ncsf2(1,iarg)                            6d9s21
               ioffs=ioffs+isto*ncsf2(1,iarg)                           6d9s21
              end do                                                    6d9s21
              nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                     6d9s21
             else                                                       6d9s21
              nvv=nvirt(jsbv1)*nvirt(jsbv2)                             6d9s21
             end if                                                     6d9s21
             do ivv=0,nvv-1                                             6d9s21
              do i=0,ncsf(iarg)-1                                       6d9s21
               koff=jvdtry+ioffs+isto*i                                 6d9s21
               joff=jvdtry+ioff+nroot*i                                 6d9s21
               do ir=1,isto                                             6d9s21
                bc(koff+ir)=bc(joff+ir)                                 6d9s21
               end do                                                   6d9s21
              end do                                                    6d9s21
              ioff=ioff+nroot*ncsf(iarg)                                6d9s21
              ioffs=ioffs+isto*ncsf(iarg)                               6d9s21
             end do                                                     6d9s21
            end if                                                      6d9s21
           end do                                                       6d9s21
          end do                                                        6d9s21
         end if                                                         6d9s21
        end do                                                          6d9s21
       end do                                                           6d9s21
      end if                                                            6d9s21
      ioff=0                                                            11d10s20
      jhdd=1                                                            6d9s21
      do jsb=1,nsymb                                                      12d8s20
       jsbv12=multh(jsb,isymmrci)                                         12d8s20
       nvisv=0                                                          6d9s21
       nvnotv=0                                                         6d9s21
       do jsbv1=1,nsymb                                                 6d9s21
        jsbv2=multh(jsbv1,jsbv12)                                       6d9s21
        if(jsbv2.ge.jsbv1)then                                          6d9s21
         if(jsbv1.eq.jsbv2)then                                         6d9s21
          nvisv=nvisv+nvirt(jsbv1)                                      6d9s21
          nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                         6d9s21
         else                                                           6d9s21
          nvv=nvirt(jsbv1)*nvirt(jsbv2)                                 6d9s21
         end if                                                         6d9s21
         nvnotv=nvnotv+nvv                                              6d9s21
        end if                                                          6d9s21
       end do                                                           6d9s21
       do ii=mdon+1,mdoo+1                                              6d9s21
        if(nff2(ii,jsb).gt.0)then                                       6d9s21
         iarg=ii-mdon                                                   6d9s21
         call ilimts(1,nff2(ii,jsb),mynprocg,mynowprog,il,ih,i1s,i1e,   6d9s21
     $     i2s,i2e)                                                     6d9s21
         il8=il                                                         6d9s21
         ih8=ih                                                         6d9s21
         irow=isto*(nvisv*ncsf2(1,iarg)+nvnotv*ncsf(iarg))              7d8s21
         ioffs=ioff                                                     6d9s21
         kvdtry=ivdtry+ioff                                             6d9s21
         call ddi_put(bc,ibc,nhddiag(ii,jsb,1),i1,irow,il8,ih8,         11d15s22
     $        bc(kvdtry))                                               11d15s22
         do iff=il,ih                                                   6d9s21
          do jsbv1=1,nsymb                                                   12d8s20
           jsbv2=multh(jsbv1,jsbv12)
           if(jsbv2.ge.jsbv1)then
            if(jsbv1.eq.jsbv2)then                                      6d9s21
             do iv=0,nvirt(jsbv1)-1                                     6d9s21
              do i=0,ncsf2(1,iarg)*isto-1                               6d9s21
               bc(kvdtry+i)=bc(kvdtry+i)*hdd(jhdd)                      6d9s21
              end do                                                    6d9s21
              jhdd=jhdd+1                                               6d9s21
              kvdtry=kvdtry+ncsf2(1,iarg)*isto                          6d9s21
             end do                                                     6d9s21
             ioff=ioff+ncsf2(1,iarg)*isto*nvirt(jsbv1)                  6d9s21
             nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                      6d9s21
            else                                                        6d9s21
             nvv=nvirt(jsbv1)*nvirt(jsbv2)                              6d9s21
            end if                                                      6d9s21
            do ivv=0,nvv-1                                              6d9s21
             do i=0,ncsf(iarg)*isto-1                                   6d9s21
              bc(kvdtry+i)=bc(kvdtry+i)*hdd(jhdd)                       6d9s21
             end do                                                     6d9s21
             jhdd=jhdd+1                                                6d9s21
             kvdtry=kvdtry+ncsf(iarg)*isto                              6d9s21
            end do                                                      6d9s21
            ioff=ioff+ncsf(iarg)*isto*nvv                               6d9s21
           end if                                                       6d9s21
          end do                                                        6d9s21
         end do                                                         6d9s21
         kvdtry=ivdtry+ioffs                                             6d9s21
         call ddi_put(bc,ibc,nhddiag(ii,jsb,2),i1,irow,il8,ih8,         11d15s22
     $        bc(kvdtry))                                               11d15s22
        end if                                                          6d9s21
       end do                                                           6d9s21
      end do                                                            6d9s21
      if(nsing.gt.0)then                                                3d19s21
       jvstry=1                                                          11d10s20
       jvstryn=1                                                        3d30s21
       jhss=1                                                            11d10s20
       do isb=1,nsymb                                                    11d10s20
        isbv=multh(isb,isymmrci)                                         11d10s20
        irow4=nvirt(isbv)*nroot                                         3d30s21
        irown=nvirt(isbv)*isto                                          3d30s21
        irow=irown                                                      3d30s21
        do ii=mdon+1,mdoo+1                                              11d10s20
         iarg=ii-mdon                                                    11d10s20
         ndblock=ncsf(iarg)*nff1(ii,isb)                                 11d10s20
         if(ndblock.gt.0)then                                            11d10s20
          icol4=ndblock                                                  11d10s20
          call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)   11d10s20
          il8=il                                                         11d10s20
          ih8=ih                                                         11d10s20
          nhere=ih+1-il                                                  11d10s20
          if(isto.ne.nroot)then                                         3d30s21
           jo=jvstry-1                                                  3d30s21
           jn=jvstryn-1                                                 3d30s21
           do if=1,nhere                                                3d30s21
            do iv=1,nvirt(isbv)                                         3d30s21
             do ir=1,isto                                               3d30s21
              vstry(jn+ir)=vstry(jo+ir)                                 3d30s21
             end do                                                     3d30s21
             jn=jn+isto                                                 3d30s21
             jo=jo+nroot                                                3d30s21
            end do                                                      3d30s21
           end do                                                       3d30s21
          end if                                                        3d30s21
          call ddi_put(bc,ibc,nhsdiag(ii,isb,1),i1,irow,il8,ih8,        11d15s22
     $         vstry(jvstryn))                                          11d15s22
          igtmp=ibcoff                                                   12d15s20
          ibcoff=igtmp+irown*nhere                                      3d30s21
          call enough('updatexuc.  5',bc,ibc)
          jjvstry=jvstryn                                               3d30s21
          jjgtmp=igtmp                                                   12d15s20
          ic=0                                                           12d15s20
          do iff=1,nff1(ii,isb)                                          11d10s20
           khss=jhss                                                     12d17s20
           do m=1,ncsf(iarg)                                             11d10s20
            ic=ic+1                                                      11d10s20
            if(ic.ge.il.and.ic.le.ih)then                                11d10s20
             khss=jhss                                                   12d14s20
             do iv=1,nvirt(isbv)                                         11d10s20
              do ir=1,isto                                              3d30s21
               bc(jjgtmp)=vstry(jjvstry)*hss(khss)                       12d15s20
               jjgtmp=jjgtmp+1                                           12d15s20
               jjvstry=jjvstry+1                                         12d15s20
              end do                                                     12d15s20
              khss=khss+1                                                12d15s20
             end do                                                      12d15s20
            end if                                                       12d15s20
           end do                                                        12d15s20
           jhss=khss                                                     12d15s20
          end do                                                         12d15s20
          call ddi_put(bc,ibc,nhsdiag(ii,isb,2),i1,irow,il8,ih8,        11d15s22
     $         bc(igtmp))                                               11d15s22
          ibcoff=igtmp                                                   12d15s20
          jvstry=jvstry+irow4*nhere                                      3d30s21
          jvstryn=jvstryn+irown*nhere                                   3d30s21
         end if                                                          11d10s20
        end do                                                           11d10s20
       end do                                                            11d10s20
      end if                                                            3d19s21
      ibcoff=ibcoffo                                                    10d7s21
      call dws_synca                                                    11d10s20
      return
      end
