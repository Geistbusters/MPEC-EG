c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine mtimesi(q,p,nun,nbk,multh,isw,istoh,icbk,nvirt,ihcis,  6d26s23
     $     noca,nthism,ntot,ibk,bc,ibc)                                 6d26s23
      implicit real*8 (a-h,o-z)                                         6d26s23
      dimension q(*),multh(8,8),icbk(5,*),nvirt(*),ihcis(*),noca(*),    6d26s23
     $     ibk(*),p(*)                                                  6d26s23
      include "common.store"                                            6d26s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      do i=1,nun                                                        6d26s23
       q(i)=0d0                                                         6d26s23
      end do                                                            3d3s17
      ioff1=0
      do ib1=1,nbk
       iv1=ibk(ib1)
       io1=multh(iv1,isw)
       ioff2=0
       do ib2=1,ib1
        iv2=ibk(ib2)
        io2=multh(iv2,isw)
        do i=1,istoh
         if(icbk(5,i).eq.1.and.icbk(1,i).eq.io1.and.
     $         icbk(2,i).eq.io2.and.icbk(3,i).eq.iv1.and.
     $         icbk(4,i).eq.iv2)then
          nv1=nvirt(iv1)
          nv2=nvirt(iv2)
          call ilimts(nv1,nv2,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,
     $            i2e)
          i10=i1s
          i1n=nv1
          iih=ihcis(i)
          nrow=noca(io1)*noca(io2)
          do i2=i2s,i2e
           if(i2.eq.i2e)i1n=i1e
           do i1=i10,i1n
            do i4=1,noca(io2)                                           6d26s23
             i4p=i4+noca(io2)*(i2-1)+ioff2
             do i3=1,noca(io1)                                          6d26s23
              i3p=i3+noca(io1)*(i1-1)+ioff1
              do icol=0,nthism                                          3d3s17
               iadp=i3p+ntot*icol                                       6d26s23
               iadq=i4p+ntot*icol                                       6d26s23
               q(iadq)=q(iadq)+bc(iih)*p(iadp)                          6d26s23
              end do                                                    3d3s17
              if(icbk(1,i).ne.icbk(2,i))then                            3d6s17
               do icol=0,nthism                                          3d3s17
                iadp=i4p+ntot*icol                                      6d26s23
                iadq=i3p+ntot*icol                                      6d26s23
                q(iadq)=q(iadq)+bc(iih)*p(iadp)                         6d26s23
               end do                                                    3d3s17
              end if                                                    3d3s17
              iih=iih+1
             end do
            end do
           end do
           i10=1
          end do
         else if(icbk(5,i).ne.1.and.icbk(1,i).eq.io1.and.
     $          icbk(2,i).eq.io2.and.icbk(3,i).eq.iv2.and.
     $          icbk(4,i).eq.iv1)then
          nv1=nvirt(iv2)
          nv2=nvirt(iv1)
          call ilimts(nv1,nv2,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,
     $            i2e)
          i10=i1s
          i1n=nv1
          iih=ihcis(i)
          nrow=noca(io1)*noca(io2)
          do i2=i2s,i2e
           if(i2.eq.i2e)i1n=i1e
           do i1=i10,i1n
            do i4=1,noca(io2)                                           6d26s23
             i4p=i4+noca(io2)*(i1-1)+ioff2
             do i3=1,noca(io1)                                          6d26s23
              i3p=i3+noca(io1)*(i2-1)+ioff1
              do icol=0,nthism                                          3d3s17
               iadp=i3p+ntot*icol                                       6d26s23
               iadq=i4p+ntot*icol                                       6d26s23
               q(iadq)=q(iadq)+bc(iih)*p(iadp)                          6d26s23
              end do                                                    3d3s17
              iih=iih+1
             end do
            end do
           end do
           i10=1
          end do
         end if
        end do
        ioff2=ioff2+noca(io2)*nvirt(iv2)
       end do
       ioff1=ioff1+noca(io1)*nvirt(iv1)
      end do
      call dws_gsumf(q,nun)                                             6d26s23
      return
      end                                                               6d26s23
