c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine reloadvs(vsold,vstry,ihsdiag,nroot,nff1,ncsf,nvirt,    4d14s21
     $     mdon,mdoo,nsymb,multh,isymmrci,sdig,bc,ibc)                  11d9s22
      implicit real*8 (a-h,o-z)                                         4d14s21
      integer*8 ihsdiag(mdoo+1,nsymb,2),i18,i28,i38,i48                 4d14s21
      dimension vsold(*),nff1(mdoo+1,nsymb),ncsf(*),nvirt(*),           4d15s21
     $     multh(8,8),vstry(*),sdig(*)                                  4d14s21
      include "common.store"                                            4d14s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
c
c     take current best vectors (in vsold) and overwrite guess vectors  4d14s21
c     in hsdiag
c
      ioff=1                                                            4d14s21
      ioffd=1                                                           4d14s21
      sumt=0d0
      sumt2=0d0                                                         4d14s21
      do isb=1,nsymb                                                    4d14s21
       isbv=multh(isb,isymmrci)                                         4d14s21
       do nclo=mdon,mdoo                                                 4d14s21
        nclop=nclo+1                                                     4d14s21
        iarg=nclop-mdon                                                  4d14s21
        if(min(nff1(nclop,isb),nvirt(isbv)).gt.0)then                   4d15s21
         nnc=nff1(nclop,isb)*ncsf(iarg)                                 4d15s21
         call ilimts(ncsf(iarg),nff1(nclop,isb),mynprocg,mynowprog,il,  4d15s21
     $        ih,i1s,i1e,i2s,i2e)                                       4d15s21
         nhere=ih+1-il                                                  4d14s21
         if(nhere.gt.0)then                                             4d14s21
          i18=1                                                         4d14s21
          i28=nvirt(isbv)*nroot                                         4d14s21
          i38=il                                                        4d14s21
          i48=ih                                                        4d14s21
          call ddi_put(bc,ibc,ihsdiag(nclop,isb,1),i18,i28,i38,i48,     11d15s22
     $         vsold(ioff))                                             11d15s22
          igtmp=ibcoff                                                  4d14s21
          ibcoff=igtmp+nhere*nvirt(isbv)*nroot                          4d14s21
          call enough('reloadvs.  1',bc,ibc)
          jgtmp=igtmp                                                   4d14s21
          joff=ioff                                                     4d14s21
          i10=i1s                                                       4d15s21
          i1n=ncsf(iarg)                                                4d15s21
          do i2=i2s,i2e                                                 4d15s21
           if(i2.eq.i2e)i1n=i1e                                         4d15s21
           do i1=i10,i1n                                                4d15s21
            do iv=0,nvirt(isbv)-1                                       4d15s21
             do ir=0,nroot-1                                            4d15s21
              bc(jgtmp+ir)=vsold(joff+ir)*sdig(ioffd+iv)                4d15s21
             end do                                                     4d15s21
             jgtmp=jgtmp+nroot                                          4d15s21
             joff=joff+nroot                                            4d15s21
            end do                                                      4d15s21
           end do                                                       4d15s21
           ioffd=ioffd+nvirt(isbv)                                      4d15s21
           i10=1                                                        4d15s21
          end do                                                        4d15s21
          call ddi_put(bc,ibc,ihsdiag(nclop,isb,2),i18,i28,i38,i48,     11d15s22
     $         bc(igtmp))                                               11d15s22
          ibcoff=igtmp                                                  4d14s21
          do i=0,nhere*nvirt(isbv)*nroot-1
           sumt=sumt+vsold(ioff+i)**2
           sumt2=sumt2+vstry(ioff+i)**2                                    4d14s21
           vstry(ioff+i)=vsold(ioff+i)                                  4d14s21
          end do
          ioff=ioff+nhere*nvirt(isbv)*nroot                             4d14s21
         end if                                                         4d14s21
        end if                                                          4d14s21
       end do                                                            4d14s21
      end do                                                            4d14s21
      call dws_synca                                                    4d16s21
      return                                                            4d14s21
      end                                                               4d14s21
