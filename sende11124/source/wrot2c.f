c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine wrot2c(vd,vdn,nfdat,isymmrci,nvirt,nsymb,multh,ovr,    2d13s22
     $     nroot,vecon,bc,ibc)                                          11d9s22
      implicit real*8 (a-h,o-z)                                         8d21s21
      dimension vd(*),nfdat(5,4,*),nvirt(*),multh(8,8),ovr(nroot,*),    2d13s22
     $     vdn(*),vecon(*)                                              2d13s22
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      include "common.store"                                            2d13s22
      ioff=0                                                            8d21s21
      ioffn=0                                                           2d13s22
      do isb=1,nsymb                                                    8d21s21
       isbv12=multh(isb,isymmrci)                                       8d21s21
       do isbv1=1,nsymb                                                 8d21s21
        isbv2=multh(isbv1,isbv12)                                       8d21s21
        if(isbv2.ge.isbv1)then                                          8d21s21
         if(isbv1.eq.isbv2)then                                         8d21s21
          ivdn=ioffn+1                                                  2d13s22
          nn=nvirt(isbv1)*nroot                                         2d13s22
          nnm=nn-1                                                      2d13s22
          itmp=ibcoff
          ibcoff=itmp+nn                                                2d13s22
          call enough('wrot2c.  1',bc,ibc)
          do i=1,nfdat(2,1,isb)                                         2d13s22
           call dgemm('n','n',nvirt(isbv1),nroot,nroot,1d0,vdn(ivdn),   2d13s22
     $          nvirt(isbv1),vecon,nroot,0d0,bc(itmp),nvirt(isbv1),     2d13s22
     d' wrot2c.  1')
           do j=0,nnm                                                   2d13s22
            vdn(ivdn+j)=bc(itmp+j)                                      2d13s22
           end do                                                       2d13s22
           ivdn=ivdn+nn                                                 2d13s22
          end do                                                        2d13s22
          ivd=ioff+1                                                    2d13s22
          do i=1,nfdat(3,1,isb)                                         2d13s22
           call dgemm('n','n',nvirt(isbv1),nroot,nroot,1d0,vd(ivd),     2d13s22
     $          nvirt(isbv1),vecon,nroot,0d0,bc(itmp),nvirt(isbv1),     2d13s22
     d' wrot2c.  2')
           do j=0,nnm                                                   2d13s22
            vd(ivd+j)=bc(itmp+j)                                        2d13s22
           end do                                                       2d13s22
           ivd=ivd+nn                                                   2d13s22
          end do                                                        2d13s22
          ibcoff=itmp                                                   2d13s22
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
          ioffn=ioffn+nvirt(isbv1)*nroot*nfdat(2,1,isb)                 2d13s22
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         8d21s21
         else                                                           8d21s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 8d21s21
         end if                                                         8d21s21
         nn=nvv*nroot                                                   2d13s22
         nnm=nn-1                                                       2d13s22
         itmp=ibcoff                                                    2d13s22
         ibcoff=itmp+nn                                                 2d13s22
         call enough('wrot2c.  2',bc,ibc)
         do l=1,4                                                       8d21s21
          if(nfdat(3,l,isb).gt.0)then
           ivdn=ioffn+1                                                  2d13s22
           do i=1,nfdat(2,l,isb)                                         2d13s22
            call dgemm('n','n',nvv,nroot,nroot,1d0,vdn(ivdn),           2d13s22
     $          nvv,vecon,nroot,0d0,bc(itmp),nvv,                       2d13s22
     d' wrot2c.  3')
            do j=0,nnm                                                   2d13s22
             vdn(ivdn+j)=bc(itmp+j)                                      2d13s22
            end do                                                       2d13s22
            ivdn=ivdn+nn                                                 2d13s22
           end do                                                        2d13s22
           ivd=ioff+1                                                    2d13s22
           do i=1,nfdat(3,l,isb)                                         2d13s22
            call dgemm('n','n',nvv,nroot,nroot,1d0,vd(ivd),             2d13s22
     $          nvv,vecon,nroot,0d0,bc(itmp),nvv,                       2d13s22
     d' wrot2c.  4')
            do j=0,nnm                                                   2d13s22
             vd(ivd+j)=bc(itmp+j)                                        2d13s22
            end do                                                       2d13s22
            ivd=ivd+nn                                                   2d13s22
           end do                                                        2d13s22
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
           ioffn=ioffn+nfdat(2,l,isb)*nvv*nroot                         2d13s22
          end if                                                        8d21s21
         end do                                                         8d21s21
         ibcoff=itmp                                                    2d13s22
        end if                                                          8d21s21
       end do                                                           8d21s21
      end do                                                            8d21s21
      return                                                            8d21s21
      end                                                               8d21s21
