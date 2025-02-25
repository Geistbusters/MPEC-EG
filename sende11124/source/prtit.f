c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine prtit(hx,nrootu,nfdat,nvirt,nsymb,multh,isymmrci,ndoub,11d10s22
     $     bc,ibc)                                                      11d10s22
      implicit real*8 (a-h,o-z)
      dimension hx(*),nfdat(5,4,*),nvirt(*),multh(8,8)
      include "common.store"
      write(6,*)('hx in orthogonal basis ')
      xnan=-1d0
      itmpv=ibcoff
      ibcoff=itmpv+ndoub*nrootu
      call enough('prtit.  1',bc,ibc)
      do i=itmpv,ibcoff-1
       bc(i)=xnan
      end do
      ioff=1
      jtmpv=itmpv
      do isb=1,nsymb
       isbv12=multh(isb,isymmrci)
       do isbv1=1,nsymb
        isbv2=multh(isbv1,isbv12)
        if(isbv2.ge.isbv1)then
         if(isbv1.eq.isbv2)then                                         1d18s21
          write(6,*)('visv for symmetry '),isb,isbv1
          call prntm2(hx(ioff),nvirt(isbv1)*nrootu,nfdat(3,1,isb),
     $         nvirt(isbv1)*nrootu)
          do i=0,nfdat(3,1,isb)-1
           do ir=0,nrootu-1
            do iv=0,nvirt(isbv1)-1
             ivv=((iv*(iv+1))/2)+iv
             iad=jtmpv+i+nfdat(3,1,isb)*ivv+ndoub*ir
             bc(iad)=hx(ioff+iv)
            end do
            ioff=ioff+nvirt(isbv1)
           end do
          end do
          nvvs=(nvirt(isbv1)*(nvirt(isbv1)+1))/2
          nvvt=(nvirt(isbv1)*(nvirt(isbv1)-1))/2
          isw=0
         else
          nvvs=nvirt(isbv1)*nvirt(isbv2)
          nvvt=nvvs
          isw=1
         end if
         do l=1,4
          if(nfdat(3,l,isb).gt.0)then
           if(nvvt.gt.0)then
            write(6,*)('vnotv for l = '),l,ioff
            nn=nvvt*nrootu
            call prntm2(hx(ioff),nn,nfdat(3,l,isb),nn)                   12d23s20
           end if                                                       1d18s21
           if(isbv1.eq.isbv2)then                                       12d23s20
            if(l.eq.1)then
             nvvu=nvvs                                                  1d18s21
             ip=+1
            else
             nvvu=nvvt
             ip=-1
            end if
            do i=0,nfdat(3,l,isb)-1
             do ir=0,nrootu-1
              do iv2=0,nvirt(isbv2)-1
               do iv1=0,iv2-1
                iv12=((iv2*(iv2+ip))/2)+iv1
                iadto=jtmpv+i+nfdat(3,l,isb)*iv12+ndoub*ir
                bc(iadto)=hx(ioff+iv1)
               end do
               ioff=ioff+iv2
              end do
             end do
            end do
            jtmpv=jtmpv+nvvu*nfdat(3,l,isb)
           else
            do i=0,nfdat(3,l,isb)-1
             do ir=0,nrootu-1
              do iv12=0,nvvt-1
               iadto=jtmpv+i+nfdat(3,l,isb)*iv12+ndoub*ir
               bc(iadto)=hx(ioff+iv12)
              end do
              ioff=ioff+nvvt
             end do
            end do
            jtmpv=jtmpv+nvvt*nfdat(3,l,isb)
           end if
          end if                                                        1d18s21
         end do
        end if
       end do
      end do
      call prntm2(bc(itmpv),ndoub,nrootu,ndoub)
      ibcoff=itmpv
      return
      end
