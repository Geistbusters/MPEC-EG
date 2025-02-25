c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine getgss(igss,nsymb,isymmrci,mdon,mdoo,nff1,nvirt,       1d27s23
     $     multh,ncsf,nrootu,ihsdiag,bc,ibc)                            1d27s23
      implicit real*8 (a-h,o-z)
c
c     copy gs for use in davidson correction.
c
      integer*8 ihsdiag(mdoo+1,nsymb,2),i18,i28,i38,i48                 1d27s23
      dimension nff1(mdoo+1,nsymb,2),ncsf(*),multh(8,8),nvirt(*)        1d27s23
      include "common.store"                                            1d27s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      jgss=igss                                                         1d27s23
      do isb=1,nsymb                                                    2d17s21
       isbv=multh(isb,isymmrci)                                         2d17s21
       do nclo=mdon,mdoo                                                2d17s21
        nclop=nclo+1                                                    2d17s21
        if(min(nff1(nclop,isb,1),nvirt(isbv)).gt.0)then                 2d17s21
         iarg=nclop-mdon                                                2d17s21
         call ilimts(ncsf(iarg),nff1(nclop,isb,1),mynprocg,mynowprog,   2d17s21
     $          il,ih,i1s,i1e,i2s,i2e)                                   2d17s21
         nhere=ih+1-il                                                  2d17s21
         if(nhere.gt.0)then                                             1d27s23
          ngg=nvirt(isbv)*nrootu                                         2d17s21
          i18=1                                                          2d17s21
          i28=ngg                                                        2d17s21
          i38=il                                                         2d17s21
          i48=ih                                                         2d17s21
          call ddi_get(bc,ibc,ihsdiag(nclop,isb,2),i18,i28,i38,i48,      11d15s22
     $         bc(jgss))                                                 11d15s22
          jgss=jgss+ngg*nhere                                            2d17s21
         end if                                                         1d27s23
        end if                                                          2d17s21
       end do                                                           2d17s21
      end do                                                            2d17s21
      return                                                            1d27s23
      end                                                               1d27s23
