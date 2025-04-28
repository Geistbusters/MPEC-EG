c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hsingpart(hxy,oxy,nvdim,nrootu,ihsdiag,nff1,ncsf,nvirt,1d21s21
     $     mdon,mdoo,nsymb,multh,isymmrci,vsold,gsold,ngot,bc,ibc)      11d10s22
      implicit real*8 (a-h,o-z)
c
c     compute my processors part of singles contribution to h matrix.
c
      integer*8 ihsdiag(mdoo+1,nsymb,2),i18,i28,i38,i48                 1d21s21
      dimension hxy(nvdim,*),oxy(nvdim,*),nff1(mdoo+1,nsymb,*),ncsf(*), 1d21s21
     $     multh(8,8),nvirt(*),vsold(*),gsold(*)                        1d23s21
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      joff=0                                                            1d23s21
      do isb=1,nsymb
       isbv=multh(isb,isymmrci)
       do nclo=mdon,mdoo                                                1d21s21
        nclop=nclo+1                                                    1d21s21
        if(min(nvirt(isbv),nff1(nclop,isb,1)).gt.0)then                 1d21s21
         iarg=nclop-mdon                                                 1d21s21
         nxx=nff1(nclop,isb,1)*ncsf(iarg)                               1d21s21
         call ilimts(1,nxx,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)    1d21s21
         nhere=ih+1-il                                                  1d21s21
         if(nhere.gt.0)then                                             1d21s21
          nrow=nvirt(isbv)*nrootu                                       1d21s21
          ntt=nrow*nhere                                                1d21s21
          ivec=ibcoff                                                   1d21s21
          ig=ivec+ntt                                                   1d21s21
          ibcoff=ig+ntt                                                 1d21s21
          call enough('hsingpart.  1',bc,ibc)
          i18=1                                                         1d21s21
          i28=nrow                                                      1d21s21
          i38=il                                                        1d21s21
          i48=ih                                                        1d21s21
          call ddi_get(bc,ibc,ihsdiag(nclop,isb,1),i18,i28,i38,i48,     11d15s22
     $         bc(ivec))                                                11d15s22
          call ddi_get(bc,ibc,ihsdiag(nclop,isb,2),i18,i28,i38,i48,     11d15s22
     $         bc(ig))                                                  11d15s22
          do icol=0,nhere-1                                             1d21s21
           do iv=0,nvirt(isbv)-1                                        1d21s21
            jvec=ivec+nrootu*(iv+nvirt(isbv)*icol)-1                    1d21s21
            jg=ig+nrootu*(iv+nvirt(isbv)*icol)-1                        1d21s21
            do ir=1,nrootu                                              1d21s21
             do jr=1,ir                                                 1d21s21
              oxy(jr,ir)=oxy(jr,ir)+bc(jvec+jr)*bc(jvec+ir)             1d21s21
              hxy(jr,ir)=hxy(jr,ir)+bc(jvec+jr)*bc(jg+ir)               1d21s21
             end do                                                     1d21s21
            end do                                                      1d21s21
            do ir=1,ngot                                                1d23s21
             irp=ir+nrootu                                              1d23s21
             do jr=1,nrootu                                             1d23s21
              jrp=jr+nrootu                                             2d11s21
              oxy(jr,irp)=oxy(jr,irp)+bc(jvec+jr)*vsold(joff+ir)        1d23s21
              hxy(jr,irp)=hxy(jr,irp)+bc(jvec+jr)*gsold(joff+ir)        1d23s21
             end do                                                     1d23s21
             do jr=1,ir                                                 1d23s21
              jrp=jr+nrootu                                             1d23s21
              oxy(jrp,irp)=oxy(jrp,irp)+vsold(joff+jr)*vsold(joff+ir)   1d23s21
              hxy(jrp,irp)=hxy(jrp,irp)+vsold(joff+jr)*gsold(joff+ir)   1d23s21
             end do                                                     1d23s21
            end do                                                      1d23s21
            joff=joff+ngot                                              1d23s21
           end do                                                       1d21s21
          end do                                                        1d21s21
          ibcoff=ivec                                                   1d21s21
         end if                                                         1d21s21
        end if
       end do                                                           1d21s21
      end do
      return
      end
