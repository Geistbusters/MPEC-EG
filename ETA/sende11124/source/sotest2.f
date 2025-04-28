c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine sotest2(ihddiag,nff2,ncsf,ncsf2,isymmrci,nvirt,test,   7d21s21
     $     nroot,mdon,mdoop,nsymb,multh,ndoub,vx,bc,ibc)                11d10s22
      implicit real*8 (a-h,o-z)                                         7d9s21
      integer*8 ihddiag(mdoop,*),i18,i28,i38,i48                        7d9s21
      dimension nff2(mdoop,*),ncsf(*),nvirt(*),test(nroot,*),multh(8,8),7d21s21
     $     ncsf2(*),vx(*)                                               7d22s21
      include "common.store"                                            7d9s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      write(6,*)('in sotest2'),isymmrci
      write(6,*)('ncsf: '),(ncsf(i),i=1,mdoop-mdon)
      write(6,*)('ncsf2: '),(ncsf2(i),i=1,mdoop-mdon)
      nn=(nroot*(nroot+1))/2                                            7d9s21
      itest=ibcoff                                                      7d9s21
      ibcoff=itest+nn                                                   7d9s21
      call enough('sotest2.  1',bc,ibc)
      do i=itest,ibcoff-1                                               7d9s21
       bc(i)=0d0                                                        7d9s21
      end do                                                            7d9s21
      do ir=1,nroot                                                     7d22s21
       vx(ir)=0d0                                                       7d22s21
      end do                                                            7d22s21
      ndoub=0
      i18=1                                                             7d9s21
      do isb=1,nsymb                                                    7d9s21
       isbv12=multh(isb,isymmrci)                                       7d21s21
       nvisv=0                                                          7d21s21
       nvnotv=0                                                         7d21s21
       do isbv1=1,nsymb                                                 7d21s21
        isbv2=multh(isbv1,isbv12)                                       7d21s21
        if(isbv2.ge.isbv1)then                                          7d21s21
         if(isbv1.eq.isbv2)then                                         7d21s21
          nvisv=nvisv+nvirt(isbv1)                                      7d21s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         7d21s21
         else                                                           7d21s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 7d21s21
         end if                                                         7d21s21
         nvnotv=nvnotv+nvv                                              7d21s21
        end if                                                          7d21s21
       end do                                                           7d21s21
       do ii=mdon+1,mdoop                                               7d9s21
        if(min(nvisv+nvnotv,nff2(ii,isb)).gt.0)then                     7d21s21
         iarg=ii-mdon                                                    7d9s21
         write(6,*)('for ii,isb: '),ii,isb,nff2(ii,isb),iarg,ncsf(iarg)
         ncol=nff2(ii,isb)                                              7d21s21
         nvv=nvisv*ncsf2(iarg)+nvnotv*ncsf(iarg)                        7d21s21
         nrow=nvv*nroot                                                 7d21s21
         call ilimts(1,ncol,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)   7d9s21
         nhere=ih+1-il                                                  7d9s21
         if(nhere.gt.0)then                                             7d9s21
          itmp=ibcoff                                                   7d9s21
          ibcoff=itmp+nrow*nhere                                        7d9s21
          call enough('sotest2.  2',bc,ibc)
          i28=nrow                                                      7d9s21
          i38=il                                                        7d9s21
          i48=ih                                                        7d9s21
          call ddi_get(bc,ibc,ihddiag(ii,isb),i18,i28,i38,i48,bc(itmp)) 11d15s22
          jtmp=itmp                                                     7d9s21
          do icol=0,nhere-1                                             7d9s21
           do iv=0,nvv-1                                                7d21s21
            jtest=itest                                                 7d9s21
            do jr=0,nroot-1                                             7d9s21
             if(abs(bc(jtmp+jr)).gt.vx(jr+1))vx(jr+1)=abs(bc(jtmp+jr))
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
         ndoub=ndoub+ncol*nvv                                           7d21s21
        end if                                                          7d9s21
       end do                                                           7d9s21
      end do                                                            7d9s21
      write(6,*)('checking ... ndoub = '),ndoub                         7d9s21
      call dws_gsumf(bc(itest),nn)                                      7d9s21
      call mpprnt2(bc(itest),nroot)
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
      write(6,*)('overlap with doubles contribution ')
      call prntm2(test,nroot,nroot,nroot)
      ibcoff=itest                                                      7d9s21
      return                                                            7d9s21
      end                                                               7d9s21
