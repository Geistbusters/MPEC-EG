c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine sotest(ihsdiag,nff1,ncsf,isymmrci,nvirt,test,nroot,    7d9s21
     $     mdon,mdoop,nsymb,multh,nsing,vx,bc,ibc)                      11d10s22
      implicit real*8 (a-h,o-z)                                         7d9s21
      integer*8 ihsdiag(mdoop,*),i18,i28,i38,i48                        7d9s21
      dimension nff1(mdoop,*),ncsf(*),nvirt(*),test(nroot,*),multh(8,8),7d30s21
     $     vx(*)                                                        7d30s21
      include "common.store"                                            7d9s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      nn=(nroot*(nroot+1))/2                                            7d9s21
      itest=ibcoff                                                      7d9s21
      ibcoff=itest+nn                                                   7d9s21
      call enough('sotest.  1',bc,ibc)
      do i=itest,ibcoff-1                                               7d9s21
       bc(i)=0d0                                                        7d9s21
      end do                                                            7d9s21
      do i=1,nroot
       vx(i)=0d0
      end do
      nsing=0
      i18=1                                                             7d9s21
      do isb=1,nsymb                                                    7d9s21
       isbv=multh(isb,isymmrci)
       do ii=mdon+1,mdoop                                               7d9s21
        if(min(nvirt(isbv),nff1(ii,isb)).gt.0)then                      7d9s21
         iarg=ii-mdon                                                    7d9s21
         nrow=nvirt(isbv)*nroot                                         7d9s21
         ncol=nff1(ii,isb)*ncsf(iarg)                                   7d9s21
         call ilimts(1,ncol,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)   7d9s21
         nhere=ih+1-il                                                  7d9s21
         if(nhere.gt.0)then                                             7d9s21
          itmp=ibcoff                                                   7d9s21
          ibcoff=itmp+nrow*nhere                                        7d9s21
          call enough('sotest.  2',bc,ibc)
          i28=nrow                                                      7d9s21
          i38=il                                                        7d9s21
          i48=ih                                                        7d9s21
          call ddi_get(bc,ibc,ihsdiag(ii,isb),i18,i28,i38,i48,bc(itmp)) 11d15s22
          jtmp=itmp                                                     7d9s21
          do icol=0,nhere-1                                             7d9s21
           do iv=0,nvirt(isbv)-1                                        7d9s21
            jtest=itest                                                 7d9s21
            do jr=0,nroot-1                                             7d9s21
             jrp=jr+1
             vx(jrp)=max(vx(jrp),abs(bc(jtmp+jr)))                      7d30s21
             do kr=0,jr                                                 7d9s21
              bc(jtest)=bc(jtest)+bc(jtmp+kr)*bc(jtmp+jr)               7d9s21
              jtest=jtest+1                                             7d9s21
             end do                                                     7d9s21
            end do                                                      7d9s21
            jtmp=jtmp+nroot                                             7d9s21
           end do                                                       7d9s21
          end do                                                        7d9s21
          ibcoff=itmp                                                   7d9s21
         end if                                                         7d9s21
         nsing=nsing+ncol*nvirt(isbv)                                   7d9s21
        end if                                                          7d9s21
       end do                                                           7d9s21
      end do                                                            7d9s21
      call dws_gsumf(bc(itest),nn)                                      7d9s21
      jtest=itest                                                       7d9s21
      do jr=1,nroot                                                     7d9s21
       do kr=1,jr-1                                                     7d9s21
        test(kr,jr)=test(kr,jr)+bc(jtest)                               7d9s21
        test(jr,kr)=test(jr,kr)+bc(jtest)                               7d9s21
        jtest=jtest+1                                                   7d9s21
       end do                                                           7d9s21
       test(jr,jr)=test(jr,jr)+bc(jtest)                                7d9s21
       jtest=jtest+1                                                    7d9s21
      end do                                                            7d9s21
      ibcoff=itest                                                      7d9s21
      return                                                            7d9s21
      end                                                               7d9s21
