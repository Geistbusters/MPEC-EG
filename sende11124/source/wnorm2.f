c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine wnorm2(iddiag,nff2,ncsf,ncsf2,mdon,mdoop,isymmrci,     8d19s21
     $     nvirt,nsymb,multh,ovr,nroot,xnorm,bc,ibc)                    11d10s22
      implicit real*8 (a-h,o-z)
      integer*8 iddiag(mdoop,nsymb),i18,i28,i38,i48
      dimension nff2(mdoop,nsymb),ncsf(*),nvirt(*),multh(8,8),
     $     ovr(nroot,*),ncsf2(*)                                        8d19s21
      data loopx/5120/
      include "common.store"                                            8d19s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      i18=1                                                             8d19s21
      do ii=mdon+1,mdoop                                                8d19s21
       iarg=ii-mdon                                                     8d19s21
       do isb=1,nsymb                                                   8d19s21
        isbv12=multh(isb,isymmrci)                                      8d19s21
        nvisv=0                                                         8d19s21
        nvnotv=0                                                        8d19s21
        do isbv1=1,nsymb                                                8d19s21
         isbv2=multh(isbv1,isbv12)                                      8d19s21
         if(isbv2.ge.isbv1)then                                         8d19s21
          if(isbv1.eq.isbv2)then                                        8d19s21
           nvisv=nvisv+nvirt(isbv1)                                     8d19s21
           nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        8d19s21
          else                                                          8d19s21
           nvv=nvirt(isbv1)*nvirt(isbv2)                                8d19s21
          end if                                                        8d19s21
          nvnotv=nvnotv+nvv                                             8d19s21
         end if                                                         8d19s21
        end do                                                          8d19s21
        if(min(nff2(ii,isb),nvisv+nvnotv).gt.0)then                     8d20s21
         nrow=nroot*(ncsf2(iarg)*nvisv+ncsf(iarg)*nvnotv)               8d19s21
         ncol=nff2(ii,isb)                                              8d19s21
         call ilimts(1,ncol,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)   8d21s21
         nhere=ih+1-il                                                  8d21s21
         if(nhere.gt.0)then                                             8d21s21
          itmp=ibcoff                                                    8d19s21
          ibcoff=itmp+nrow*nhere                                        8d21s21
          call enough('wnorm2.  1',bc,ibc)
          i28=nrow                                                       8d19s21
          i38=il                                                        8d21s21
          i48=ih                                                        8d21s21
          call ddi_get(bc,ibc,iddiag(ii,isb),i18,i28,i38,i48,bc(itmp))  11d15s22
          do i=itmp,ibcoff-1                                            8d21s21
           bc(i)=bc(i)*xnorm                                            8d21s21
          end do                                                        8d21s21
          jtmp=itmp-1                                                    8d19s21
          do icol=1,nhere                                               8d21s21
           do isbv1=1,nsymb                                              8d19s21
            isbv2=multh(isbv1,isbv12)                                    8d19s21
            if(isbv2.ge.isbv1)then                                       8d19s21
             if(isbv1.eq.isbv2)then                                      8d19s21
              do iv=1,nvirt(isbv1)                                       8d19s21
               do i=1,ncsf2(iarg)                                        8d19s21
                do ir=1,nroot                                                8d19s21
                 do jr=1,nroot                                               8d19s21
                  ovr(jr,ir)=ovr(jr,ir)+bc(jtmp+jr)*bc(jtmp+ir)              8d19s21
                 end do                                                      8d19s21
                end do                                                       8d19s21
                jtmp=jtmp+nroot                                              8d19s21
               end do                                                    8d19s21
              end do                                                     8d19s21
              nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                      8d19s21
             else                                                        8d19s21
              nvv=nvirt(isbv1)*nvirt(isbv2)                              8d19s21
             end if                                                      8d19s21
             do ivv=1,nvv                                                8d19s21
              do i=1,ncsf(iarg)                                          8d19s21
               do ir=1,nroot                                                8d19s21
                do jr=1,nroot                                               8d19s21
                 ovr(jr,ir)=ovr(jr,ir)+bc(jtmp+jr)*bc(jtmp+ir)              8d19s21
                end do                                                      8d19s21
               end do                                                       8d19s21
               jtmp=jtmp+nroot                                              8d19s21
              end do                                                     8d19s21
             end do                                                      8d19s21
            end if                                                       8d19s21
           end do                                                        8d19s21
          end do                                                         8d19s21
          call ddi_put(bc,ibc,iddiag(ii,isb),i18,i28,i38,i48,bc(itmp))  11d15s22
          ibcoff=itmp                                                    8d19s21
         end if                                                         8d21s21
        end if                                                          8d19s21
       end do                                                           8d19s21
      end do                                                            8d19s21
      return                                                            8d18s21
      end                                                               8d18s21
