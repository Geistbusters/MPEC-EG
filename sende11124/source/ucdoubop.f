c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine ucdoubop(nsymb,isymmrci,multh,nvirt,mdon,mdoo,nff2,    8d24s22
     $     ncsf,ncsf2,ihddiag,jrm,krm,dot,icode,nroot,bc,ibc)           11d10s22
      implicit real*8 (a-h,o-z)
      integer*8 i18,i28,i38,i48,ihddiag(mdoo+1,nsymb,*)
      dimension multh(8,8),nvirt(*),nff2(mdoo+1,*),ncsf(*),ncsf2(4,*)
c
c     icode=1
c     dot between root jrm+1 and krm+1, return in dot
c     icode=2
c     subtract dot times root krm+1 from root jrm+1
c     icode=3
c     scale root jrm by dot
c
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      include "common.store"                                            8d24s22
      do isb=1,nsymb                                                    8d18s22
       isbv12=multh(isb,isymmrci)                                       8d18s22
       nvisv=0                                                             6d9s21
       nvnotv=0                                                            6d9s21
       do isbv1=1,nsymb                                                 8d18s22
        isbv2=multh(isbv1,isbv12)                                       8d18s22
        if(isbv2.ge.isbv1)then                                          8d18s22
         if(isbv1.eq.isbv2)then                                         8d18s22
          nvisv=nvisv+nvirt(isbv1)                                      8d19s22
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         8d18s22
         else                                                           8d18s22
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 8d18s22
         end if                                                         8d18s22
         nvnotv=nvnotv+nvv                                              8d19s22
        end if                                                          8d19s22
       end do                                                           8d19s22
       do nclop=mdon+1,mdoo+1                                           8d19s22
        if(nff2(nclop,isb).gt.0)then                                    8d19s22
         iarg=nclop-mdon                                                8d19s22
         call ilimts(1,nff2(nclop,isb),mynprocg,mynowprog,il,ih,i1s,    8d19s22
     $           i1e,i2s,i2e)                                              8d19s22
         nhere=ih+1-il                                                  8d19s22
         if(nhere.gt.0)then                                             8d19s22
          nrow=nroot*(nvisv*ncsf2(1,iarg)+nvnotv*ncsf(iarg))            8d19s22
          itmp=ibcoff                                                   8d19s22
          ibcoff=itmp+nrow*nhere                                        8d19s22
          i18=1                                                         8d19s22
          i28=nrow                                                      8d19s22
          i38=il                                                        8d19s22
          i48=ih                                                        8d19s22
          call ddi_get(bc,ibc,ihddiag(nclop,isb,1),i18,i28,i38,i48,     11d15s22
     $         bc(itmp))                                                11d15s22
          jtmp=itmp+jrm                                                 8d19s22
          ktmp=itmp+krm
          do iff=il,ih                                                  8d19s22
           do isbv1=1,nsymb                                             8d19s22
            isbv2=multh(isbv1,isbv12)                                   8d19s22
            if(isbv2.ge.isbv1)then                                      8d19s22
             if(isbv1.eq.isbv2)then                                     8d19s22
              if(icode.eq.1)then
               do iv=0,nvirt(isbv1)-1                                   8d19s22
                do i=0,ncsf2(1,iarg)-1                                  8d19s22
                 dot=dot+bc(jtmp)*bc(ktmp)                                8d24s22
                 jtmp=jtmp+nroot                                        8d19s22
                 ktmp=ktmp+nroot                                          8d24s22
                end do                                                  8d19s22
               end do                                                   8d19s22
              else if(icode.eq.2)then                                     8d24s22
               do iv=0,nvirt(isbv1)-1                                   8d19s22
                do i=0,ncsf2(1,iarg)-1                                  8d19s22
                 bc(jtmp)=bc(jtmp)-dot*bc(ktmp)
                 jtmp=jtmp+nroot                                        8d19s22
                 ktmp=ktmp+nroot                                          8d24s22
                end do                                                  8d19s22
               end do                                                   8d19s22
              else if(icode.eq.3)then                                     8d24s22
               do iv=0,nvirt(isbv1)-1                                   8d19s22
                do i=0,ncsf2(1,iarg)-1                                  8d19s22
                 bc(jtmp)=bc(jtmp)*dot                                    8d24s22
                 jtmp=jtmp+nroot                                        8d19s22
                end do                                                  8d19s22
               end do                                                   8d19s22
              end if                                                    8d24s22
              nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                     8d19s22
             else                                                       8d19s22
              nvv=nvirt(isbv1)*nvirt(isbv2)                             8d19s22
             end if                                                     8d19s22
             if(icode.eq.1)then                                           8d24s22
              do ivv=0,nvv-1                                            8d19s22
               do i=0,ncsf(iarg)-1                                      8d19s22
                dot=dot+bc(jtmp)*bc(ktmp)                                 8d24s22
                jtmp=jtmp+nroot                                         8d19s22
                ktmp=ktmp+nroot                                         8d19s22
               end do                                                   8d19s22
              end do                                                    8d19s22
             else if(icode.eq.2)then                                           8d24s22
              do ivv=0,nvv-1                                            8d19s22
               do i=0,ncsf(iarg)-1                                      8d19s22
                bc(jtmp)=bc(jtmp)-dot*bc(ktmp)                            8d24s22
                jtmp=jtmp+nroot                                         8d19s22
                ktmp=ktmp+nroot                                         8d19s22
               end do                                                   8d19s22
              end do                                                    8d19s22
             else if(icode.eq.3)then                                           8d24s22
              do ivv=0,nvv-1                                            8d19s22
               do i=0,ncsf(iarg)-1                                      8d19s22
                bc(jtmp)=bc(jtmp)*dot                                     8d24s22
                jtmp=jtmp+nroot                                         8d19s22
               end do                                                   8d19s22
              end do                                                    8d19s22
             end if                                                     8d24s22
            end if                                                      8d19s22
           end do                                                       8d19s22
          end do                                                        8d19s22
          if(icode.eq.2.or.icode.eq.3)then                              8d24s22
           call ddi_put(bc,ibc,ihddiag(nclop,isb,1),i18,i28,i38,i48,    11d15s22
     $         bc(itmp))                                                11d15s22
          end if                                                        8d24s22
          ibcoff=itmp                                                   8d19s22
         end if                                                         8d19s22
        end if                                                          8d19s22
       end do                                                           8d19s22
      end do                                                            8d19s22
      return
      end
