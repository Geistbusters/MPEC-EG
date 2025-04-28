c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine reordergv(hx,nrootu,nfdat,nvirt,nsymb,multh,isymmrci,  1d20s21
     $     ndoub,iflag,bc,ibc)                                          11d9s22
c                                                                       1d20s21
c     iflag gt 0: hx is in visv,vnotv order, return it in ntot order.   1d20s21
c     iflag lt 0: hx is in ntot order, return it in visv,vnotv order.   1d20s21
c                                                                       1d20s21
      implicit real*8 (a-h,o-z)                                         1d20s21
      dimension hx(*),nfdat(5,4,*),nvirt(*),multh(8,8)                  1d20s21
      include "common.store"                                            1d20s21
      itmpv=ibcoff
      ibcoff=itmpv+ndoub*nrootu
      call enough('reordergv.  1',bc,ibc)
      ioff=1
      jtmpv=itmpv
      do isb=1,nsymb
       isbv12=multh(isb,isymmrci)
       do isbv1=1,nsymb
        isbv2=multh(isbv1,isbv12)
        if(isbv2.ge.isbv1)then
         if(isbv1.eq.isbv2)then                                         1d18s21
          do i=0,nfdat(3,1,isb)-1
           do ir=0,nrootu-1
            if(iflag.gt.0)then                                          1d20s21
             do iv=0,nvirt(isbv1)-1
              ivv=((iv*(iv+1))/2)+iv
              iad=jtmpv+i+nfdat(3,1,isb)*ivv+ndoub*ir
              bc(iad)=hx(ioff+iv)
             end do
             ioff=ioff+nvirt(isbv1)
            else                                                        1d20s21
             do iv=0,nvirt(isbv1)-1
              ivv=((iv*(iv+1))/2)+iv
              iad=ioff+i+nfdat(3,1,isb)*ivv+ndoub*ir
              bc(jtmpv+iv)=hx(iad)
             end do
             jtmpv=jtmpv+nvirt(isbv1)                                   1d20s21
            end if                                                      1d20s21
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
              if(iflag.gt.0)then                                        1d20s21
               do iv2=0,nvirt(isbv2)-1
                do iv1=0,iv2-1
                 iv12=((iv2*(iv2+ip))/2)+iv1
                 iadto=jtmpv+i+nfdat(3,l,isb)*iv12+ndoub*ir
                 bc(iadto)=hx(ioff+iv1)
                end do
                ioff=ioff+iv2
               end do
              else                                                      1d20s21
               do iv2=0,nvirt(isbv2)-1                                  1d20s21
                do iv1=0,iv2-1                                          1d20s21
                 iv12=((iv2*(iv2+ip))/2)+iv1                            1d20s21
                 iadf=ioff+i+nfdat(3,l,isb)*iv12+ndoub*ir               1d20s21
                 bc(jtmpv+iv1)=hx(iadf)                                 1d20s21
                end do                                                  1d20s21
                jtmpv=jtmpv+iv2                                         1d20s21
               end do                                                   1d20s21
              end if                                                    1d20s21
             end do                                                     1d20s21
            end do                                                      1d20s21
            if(iflag.gt.0)then                                          1d20s21
             jtmpv=jtmpv+nvvu*nfdat(3,l,isb)
            else                                                        1d20s21
             ioff=ioff+nvvu*nfdat(3,l,isb)                              1d20s21
            end if                                                      1d20s21
           else
            do i=0,nfdat(3,l,isb)-1
             do ir=0,nrootu-1
              if(iflag.gt.0)then                                        1d20s21
               do iv12=0,nvvt-1
                iadto=jtmpv+i+nfdat(3,l,isb)*iv12+ndoub*ir
                bc(iadto)=hx(ioff+iv12)
               end do
               ioff=ioff+nvvt
              else                                                      1d20s21
               do iv12=0,nvvt-1                                         1d20s21
                iadf=ioff+i+nfdat(3,l,isb)*iv12+ndoub*ir                1d20s21
                bc(jtmpv+iv12)=hx(iadf)                                 1d20s21
               end do
               jtmpv=jtmpv+nvvt
              end if                                                    1d20s21
             end do
            end do
            if(iflag.gt.0)then                                          1d20s21
             jtmpv=jtmpv+nvvt*nfdat(3,l,isb)
            else                                                        1d20s21
             ioff=ioff+nvvt*nfdat(3,l,isb)                              1d20s21
            end if                                                      1d20s21
           end if
          end if                                                        1d18s21
         end do
        end if
       end do
      end do
      jtmpv=itmpv-1                                                     1d20s21
      do i=1,ndoub*nrootu                                               1d20s21
       hx(i)=bc(jtmpv+i)                                                1d20s21
      end do                                                            1d20s21
      ibcoff=itmpv
      return
      end
