c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hdoubpart(hxy,oxy,nvdim,nrootu,vec,gd,nfdat,nvirt,     1d21s21
     $     nsymb,multh,isymmrci,vdold,gdold,ngot,epart,idep,bc,ibc)     11d10s22
      implicit real*8 (a-h,o-z)
c
c     doubles part of h matrix
c
      dimension hxy(nvdim,nvdim),oxy(nvdim,nvdim),vec(*),gd(*),
     $     nfdat(5,4,*),nvirt(*),multh(8,8),vdold(*),gdold(*),
     $     epart(idep,4,nsymb)                                          1d17s22
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      include "common.store"
      ioff=0                                                            1d24s21
      joff=0                                                            1d24s21
      do isb=1,nsymb                                                    1d17s22
       do l=1,4                                                         1d17s22
        do i=1,idep                                                     1d17s22
         epart(i,l,isb)=0d0                                             1d17s22
        end do                                                           7d12s21
       end do                                                            7d12s21
      end do                                                            1d17s22
      do isb=1,nsymb                                                    1d21s21
       isbv12=multh(isb,isymmrci)                                       1d21s21
       do isbv1=1,nsymb                                                 1d21s21
        isbv2=multh(isbv1,isbv12)                                       1d21s21
        if(isbv2.ge.isbv1)then                                          1d21s21
         if(isbv1.eq.isbv2)then                                         1d21s21
          do i=0,nfdat(3,1,isb)-1                                       1d21s21
           do ir=1,nrootu                                               1d21s21
            iadi=ioff+nvirt(isbv1)*(ir-1+nrootu*i)                      1d21s21
            do jr=1,ir                                                  1d21s21
             iadj=ioff+nvirt(isbv1)*(jr-1+nrootu*i)                     1d21s21
             do iv=1,nvirt(isbv1)                                       1d24s21
              orig=oxy(jr,ir)
              oxy(jr,ir)=oxy(jr,ir)+vec(iv+iadj)*vec(iv+iadi)           1d21s21
              hxy(jr,ir)=hxy(jr,ir)+vec(iv+iadj)*gd(iv+iadi)            1d21s21
             end do                                                     1d21s21
            end do                                                      1d21s21
           end do                                                       1d21s21
           do ir=1,ngot                                                 1d23s21
            irp=ir+nrootu                                               1d23s21
            iadi=joff+nvirt(isbv1)*(ir-1+ngot*i)                        1d24s21
            do jr=1,nrootu                                              1d23s21
             jrp=jr+nrootu                                              2d11s21
             iadj=ioff+nvirt(isbv1)*(jr-1+nrootu*i)                     1d21s21
             do iv=1,nvirt(isbv1)                                       1d24s21
              oxy(jr,irp)=oxy(jr,irp)+vec(iv+iadj)*vdold(iv+iadi)       1d23s21
              hxy(jr,irp)=hxy(jr,irp)+vec(iv+iadj)*gdold(iv+iadi)       1d23s21
             end do                                                     1d21s21
            end do                                                      1d21s21
           end do                                                       1d21s21
           do ir=1,ngot                                                 1d23s21
            irp=ir+nrootu                                               1d23s21
            iadi=joff+nvirt(isbv1)*(ir-1+ngot*i)                        1d24s21
            do jr=1,ir                                                  1d23s21
             jrp=jr+nrootu                                              1d23s21
             iadj=joff+nvirt(isbv1)*(jr-1+ngot*i)                       1d24s21
             do iv=1,nvirt(isbv1)                                       1d24s21
              orig=oxy(jrp,irp)
              oxy(jrp,irp)=oxy(jrp,irp)+vdold(iv+iadj)*vdold(iv+iadi)   1d23s21
              hxy(jrp,irp)=hxy(jrp,irp)+vdold(iv+iadj)*gdold(iv+iadi)   1d23s21
             end do                                                     1d21s21
            end do                                                      1d21s21
            do iv=1,nvirt(isbv1)                                        7d12s21
             epart(ir,1,isb)=epart(ir,1,isb)                            1d17s22
     $            +vdold(iv+iadi)*gdold(iv+iadi)                        1d15s22
            end do                                                      7d12s21
           end do                                                       1d21s21
          end do                                                        1d21s21
          ioff=ioff+nvirt(isbv1)*nrootu*nfdat(3,1,isb)                  1d21s21
          joff=joff+nvirt(isbv1)*ngot*nfdat(3,1,isb)                    1d24s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         1d24s21
         else                                                           1d21s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 1d21s21
         end if                                                         1d21s21
         do l=1,4                                                       1d21s21
          if(nfdat(3,l,isb).gt.0)then                                   1d21s21
           do i=0,nfdat(3,l,isb)-1                                      1d21s21
            do ir=1,nrootu                                               1d21s21
             iadi=ioff+nvv*(ir-1+nrootu*i)                              1d21s21
             do jr=1,ir                                                  1d21s21
              iadj=ioff+nvv*(jr-1+nrootu*i)                             1d21s21
              do iv=1,nvv                                               1d24s21
               orig=oxy(jr,ir)
               oxy(jr,ir)=oxy(jr,ir)+vec(iv+iadj)*vec(iv+iadi)          1d21s21
               hxy(jr,ir)=hxy(jr,ir)+vec(iv+iadj)*gd(iv+iadi)           1d21s21
              end do                                                    1d21s21
             end do                                                     1d21s21
            end do                                                      1d21s21
            do ir=1,ngot                                                1d23s21
             irp=ir+nrootu                                              1d27s21
             iadi=joff+nvv*(ir-1+ngot*i)                                1d24s21
             do jr=1,nrootu                                             1d23s21
              jrp=jr+nrootu                                             2d11s21
              iadj=ioff+nvv*(jr-1+nrootu*i)                             1d21s21
              do iv=1,nvv                                                1d24s21
               oxy(jr,irp)=oxy(jr,irp)+vec(iv+iadj)*vdold(iv+iadi)      1d23s21
               hxy(jr,irp)=hxy(jr,irp)+vec(iv+iadj)*gdold(iv+iadi)      1d23s21
              end do                                                     1d21s21
             end do                                                      1d21s21
             do jr=1,ir                                                 1d23s21
              jrp=jr+nrootu                                             1d23s21
              iadj=joff+nvv*(jr-1+ngot*i)                               1d24s21
              do iv=1,nvv                                                1d24s21
               orig=oxy(jrp,irp)
               oxy(jrp,irp)=oxy(jrp,irp)+vdold(iv+iadj)*vdold(iv+iadi)  1d23s21
               hxy(jrp,irp)=hxy(jrp,irp)+vdold(iv+iadj)*gdold(iv+iadi)  1d23s21
              end do                                                    1d21s21
             end do                                                     1d21s21
             do iv=1,nvv                                                7d12s21
              epart(ir,l,isb)=epart(ir,l,isb)                           1d17s22
     $             +vdold(iv+iadi)*gdold(iv+iadi)                       1d15s22
             end do                                                     7d12s21
            end do                                                      1d21s21
           end do                                                       1d21s21
           ioff=ioff+nvv*nrootu*nfdat(3,l,isb)                          1d21s21
           joff=joff+nvv*ngot*nfdat(3,l,isb)                            1d24s21
          end if                                                        1d21s21
         end do                                                         1d21s21
        end if                                                          1d21s21
       end do                                                           1d21s21
      end do                                                            1d21s21
      return                                                            1d21s21
      end                                                               1d21s21
