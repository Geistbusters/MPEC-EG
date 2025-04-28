c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hdoubpartuc(hxy,oxy,nvdim,nrootu,ihddiag,nff2,ncsf,
     $     ncsf2,nvirt,mdon,mdoo,nsymb,multh,isymmrci,vdold,gdold,ngot, 7d12s21
     $     epart,idep,bc,ibc)                                           11d10s22
      implicit real*8 (a-h,o-z)
c
c     compute my processors part of singles contribution to h matrix.
c
      integer*8 ihddiag(mdoo+1,nsymb,2),i18,i28,i38,i48                 1d21s21
      dimension hxy(nvdim,*),oxy(nvdim,*),nff2(mdoo+1,nsymb,*),ncsf(*), 1d21s21
     $     multh(8,8),nvirt(*),vdold(*),gdold(*),ncsf2(4,*),            7d12s21
     $     epart(idep,4,nsymb)                                          1d17s22
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      do isb=1,nsymb                                                    1d17s22
       do l=1,4                                                          7d12s21
        do i=1,idep                                                      7d12s21
         epart(i,l,isb)=0d0                                             1d17s22
        end do                                                           7d12s21
       end do                                                            7d12s21
      end do                                                            1d17s22
      joff=0                                                            1d23s21
      do isb=1,nsymb
       isbv12=multh(isb,isymmrci)                                       7d7s21
       nvisv=0                                                          7d7s21
       nvnotv=0                                                         7d7s21
       do isbv1=1,nsymb                                                 7d7s21
        isbv2=multh(isbv1,isbv12)                                       7d7s21
        if(isbv2.ge.isbv1)then                                          7d7s21
         if(isbv1.eq.isbv2)then                                         7d7s21
          nvisv=nvisv+nvirt(isbv1)                                      7d7s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         7d7s21
         else                                                           7d7s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 7d7s21
         end if                                                         7d7s21
         nvnotv=nvnotv+nvv                                              7d7s21
        end if                                                          7d7s21
       end do                                                           7d7s21
       do nclo=mdon,mdoo                                                1d21s21
        nclop=nclo+1                                                    1d21s21
        if(min(nvisv+nvnotv,nff2(nclop,isb,1)).gt.0)then                7d7s21
         iarg=nclop-mdon                                                 1d21s21
         nxx=nff2(nclop,isb,1)                                          7d7s21
         call ilimts(1,nxx,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)    1d21s21
         nhere=ih+1-il                                                  1d21s21
         if(nhere.gt.0)then                                             1d21s21
          nrow=(nvisv*ncsf2(1,iarg)+nvnotv*ncsf(iarg))*nrootu           7d7s21
          ntt=nrow*nhere                                                1d21s21
          ivec=ibcoff                                                   1d21s21
          ig=ivec+ntt                                                   1d21s21
          ibcoff=ig+ntt                                                 1d21s21
          call enough('hdoubpartuc.  1',bc,ibc)
          i18=1                                                         1d21s21
          i28=nrow                                                      1d21s21
          i38=il                                                        1d21s21
          i48=ih                                                        1d21s21
          call ddi_get(bc,ibc,ihddiag(nclop,isb,1),i18,i28,i38,i48,     11d15s22
     $         bc(ivec))                                                11d15s22
          call ddi_get(bc,ibc,ihddiag(nclop,isb,2),i18,i28,i38,i48,     11d15s22
     $         bc(ig))                                                  11d15s22
          jvec=ivec-1                                                   7d7s21
          jg=ig-1                                                       7d7s21
          do icol=0,nhere-1                                             1d21s21
           do isbv1=1,nsymb                                             7d7s21
            isbv2=multh(isbv1,isbv12)                                   7d7s21
            if(isbv2.ge.isbv1)then                                      7d7s21
             if(isbv1.eq.isbv2)then                                     7d7s21
              do iv=0,nvirt(isbv1)*ncsf2(1,iarg)-1                      7d7s21
               do ir=1,nrootu                                           7d7s21
                do jr=1,ir                                              7d7s21
                 oxy(jr,ir)=oxy(jr,ir)+bc(jvec+jr)*bc(jvec+ir)          7d7s21
                 hxy(jr,ir)=hxy(jr,ir)+bc(jvec+jr)*bc(jg+ir)            7d7s21
                end do                                                  7d7s21
               end do                                                   7d7s21
               do ir=1,ngot                                                1d23s21
                irp=ir+nrootu                                              1d23s21
                do jr=1,nrootu                                             1d23s21
                 oxy(jr,irp)=oxy(jr,irp)+bc(jvec+jr)*vdold(joff+ir)        1d23s21
                 hxy(jr,irp)=hxy(jr,irp)+bc(jvec+jr)*gdold(joff+ir)        1d23s21
                end do                                                  7d7s21
                do jr=1,ir                                                 1d23s21
                 jrp=jr+nrootu                                             1d23s21
                 oxy(jrp,irp)=oxy(jrp,irp)+vdold(joff+jr)*vdold(joff+ir)   1d23s21
                 hxy(jrp,irp)=hxy(jrp,irp)+vdold(joff+jr)*gdold(joff+ir)   1d23s21
                end do                                                     1d23s21
                epart(ir,1,isb)=epart(ir,1,isb)                         1d17s22
     $               +vdold(joff+ir)*gdold(joff+ir)                     1d15s22
               end do                                                   7d7s21
               jvec=jvec+nrootu                                         7d7s21
               jg=jg+nrootu                                             7d7s21
               joff=joff+ngot                                           7d7s21
              end do                                                    7d7s21
              nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                     7d7s21
             else                                                       7d7s21
              nvv=nvirt(isbv1)*nvirt(isbv2)                             7d7s21
             end if                                                     7d7s21
             joffp=joff                                                 7d12s21
             do l=1,4                                                   7d12s21
              do i=0,ncsf2(l,iarg)-1                                    7d12s21
               do ivv=0,nvv-1                                           7d12s21
                do ir=1,ngot                                            7d12s21
                 epart(ir,l,isb)=epart(ir,l,isb)                        1d17s22
     $                 +vdold(joffp+ir)*gdold(joffp+ir)                 7d12s21
                end do                                                  7d12s21
                joffp=joffp+ngot                                        7d12s21
               end do                                                   7d12s21
              end do                                                    7d12s21
             end do                                                     7d12s21
             nvv=nvv*ncsf(iarg)-1                                       7d7s21
             do iv=0,nvv                                                7d7s21
              do ir=1,nrootu                                            7d7s21
               do jr=1,ir                                               7d7s21
                oxy(jr,ir)=oxy(jr,ir)+bc(jvec+jr)*bc(jvec+ir)           7d7s21
                hxy(jr,ir)=hxy(jr,ir)+bc(jvec+jr)*bc(jg+ir)             7d7s21
               end do                                                   7d7s21
              end do                                                    7d7s21
              do ir=1,ngot                                                1d23s21
               irp=ir+nrootu                                              1d23s21
               do jr=1,nrootu                                             1d23s21
                oxy(jr,irp)=oxy(jr,irp)+bc(jvec+jr)*vdold(joff+ir)        1d23s21
                hxy(jr,irp)=hxy(jr,irp)+bc(jvec+jr)*gdold(joff+ir)        1d23s21
               end do                                                   7d7s21
               do jr=1,ir                                                 1d23s21
                jrp=jr+nrootu                                             1d23s21
                oxy(jrp,irp)=oxy(jrp,irp)+vdold(joff+jr)*vdold(joff+ir)   1d23s21
                hxy(jrp,irp)=hxy(jrp,irp)+vdold(joff+jr)*gdold(joff+ir)   1d23s21
               end do                                                     1d23s21
              end do                                                    7d7s21
              jvec=jvec+nrootu                                          7d7s21
              jg=jg+nrootu                                              7d7s21
              joff=joff+ngot                                            7d7s21
             end do                                                     7d7s21
            end if                                                      7d7s21
           end do                                                       7d7s21
          end do                                                        1d21s21
          ibcoff=ivec                                                   1d21s21
         end if                                                         1d21s21
        end if
       end do                                                           1d21s21
      end do
      return
      end
