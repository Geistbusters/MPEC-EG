c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine wnorm2c(vd,nfdat,isymmrci,nvirt,nsymb,multh,ovr,nroot, 8d21s21
     $     xnorm)                                                       8d21s21
      implicit real*8 (a-h,o-z)                                         8d21s21
      dimension vd(*),nfdat(5,4,*),nvirt(*),multh(8,8),ovr(nroot,*)     8d21s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      ioff=0                                                            8d21s21
      do isb=1,nsymb                                                    8d21s21
       isbv12=multh(isb,isymmrci)                                       8d21s21
       do isbv1=1,nsymb                                                 8d21s21
        isbv2=multh(isbv1,isbv12)                                       8d21s21
        if(isbv2.ge.isbv1)then                                          8d21s21
         if(isbv1.eq.isbv2)then                                         8d21s21
          sz=0d0
          do i=1,nvirt(isbv1)*nroot*nfdat(3,1,isb)                      8d21s21
           vd(ioff+i)=vd(ioff+i)*xnorm                                  8d21s21
           sz=sz+vd(ioff+i)**2
          end do                                                        8d21s21
          do i=0,nfdat(3,1,isb)-1                                       8d21s21
           do ik=1,nroot                                                8d21s21
            iadk=ioff+nvirt(isbv1)*(ik-1+nroot*i)                       8d21s21
            do ib=1,nroot                                               8d21s21
             iadb=ioff+nvirt(isbv1)*(ib-1+nroot*i)                      8d21s21
             do iv=1+mynowprog,nvirt(isbv1),mynprocg                    8d21s21
              ovr(ib,ik)=ovr(ib,ik)+vd(iadk+iv)*vd(iadb+iv)             8d21s21
             end do                                                     8d21s21
            end do                                                      8d21s21
           end do                                                       8d21s21
          end do                                                        8d21s21
          ioff=ioff+nvirt(isbv1)*nroot*nfdat(3,1,isb)                   8d21s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         8d21s21
         else                                                           8d21s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 8d21s21
         end if                                                         8d21s21
         do l=1,4                                                       8d21s21
          if(nfdat(3,l,isb).gt.0)then
           sz=0d0
           do i=1,nvv*nroot*nfdat(3,l,isb)                              8d21s21
            vd(ioff+i)=vd(ioff+i)*xnorm                                  8d21s21
            sz=sz+vd(ioff+i)**2
           end do                                                        8d21s21
           nrow=nvv*nroot
           do i=0,nfdat(3,l,isb)-1                                      8d21s21
            do ik=1,nroot                                               8d21s21
             iadk=ioff+nvv*(ik-1+nroot*i)                               8d21s21
             do ib=1,nroot                                              8d21s21
              iadb=ioff+nvv*(ib-1+nroot*i)                              8d21s21
              do ivv=1+mynowprog,nvv,mynprocg                           8d21s21
               ovr(ib,ik)=ovr(ib,ik)+vd(iadk+ivv)*vd(iadb+ivv)          8d21s21
              end do                                                    8d21s21
             end do                                                     8d21s21
            end do                                                      8d21s21
           end do                                                       8d21s21
           ioff=ioff+nfdat(3,l,isb)*nvv*nroot                           8d21s21
          end if                                                        8d21s21
         end do                                                         8d21s21
        end if                                                          8d21s21
       end do                                                           8d21s21
      end do                                                            8d21s21
      return                                                            8d21s21
      end                                                               8d21s21
