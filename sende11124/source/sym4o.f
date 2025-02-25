c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine sym4o(i4od,nocc,isblkder,isblkxder,nsblkder,nsblkxder, 11d10s22
     $     bc,ibc)                                                      11d10s22
      implicit real*8 (a-h,o-z)                                         5d16s22
      dimension i4od(*),nocc(*),isblkder(4,*),isblkxder(4,*)             5d16s22
      include "common.store"                                            5d16s22
      im0=ibcoff
      do is=1,nsblkxder
       nrow=nocc(isblkxder(1,is))*nocc(isblkxder(2,is))
       ncol=nocc(isblkxder(3,is))*nocc(isblkxder(4,is))
       if(min(nrow,ncol).gt.0)then
        call dws_gsumf(bc(i4od(is)),nrow*ncol)                          5d31s22
       end if
      end do
      do is=1,nsblkder
       nrow=nocc(isblkder(1,is))*nocc(isblkder(2,is))
       ncol=nocc(isblkder(3,is))*nocc(isblkder(4,is))
       if(nrow*ncol.gt.0)then
        imt=ibcoff
        ibcoff=imt+nrow*ncol
        call enough('sym4o.  1',bc,ibc)
        do i=0,nrow*ncol-1
         bc(imt+i)=0d0
        end do
        do is2=1,nsblkxder
         if(isblkder(1,is).eq.isblkxder(4,is2).and.
     $        isblkder(2,is).eq.isblkxder(3,is2))then
          if(isblkder(3,is).eq.isblkxder(1,is2))then
           iad1=i4od(is2)
           do i4=0,nocc(isblkxder(4,is2))-1
            do i3=0,nocc(isblkxder(3,is2))-1
             do i2=0,nocc(isblkxder(2,is2))-1
              do i1=0,nocc(isblkxder(1,is2))-1
               iad2=imt+i4+nocc(isblkxder(4,is2))*(i3
     $              +nocc(isblkxder(3,is2))*(i1
     $              +nocc(isblkxder(1,is2))*i2))
               bc(iad2)=bc(iad2)+bc(iad1)
               iad1=iad1+1
              end do
             end do
            end do
           end do
          else if(isblkder(3,is).eq.isblkxder(2,is2))then               6d29s22
           iad1=i4od(is2)
           do i4=0,nocc(isblkxder(4,is2))-1
            do i3=0,nocc(isblkxder(3,is2))-1
             do i2=0,nocc(isblkxder(2,is2))-1
              do i1=0,nocc(isblkxder(1,is2))-1
               iad2=imt+i4+nocc(isblkxder(4,is2))*(i3
     $              +nocc(isblkxder(3,is2))*(i2
     $              +nocc(isblkxder(2,is2))*i1))
               bc(iad2)=bc(iad2)+bc(iad1)
               iad1=iad1+1
              end do
             end do
            end do
           end do
          end if
         end if
         if(isblkder(2,is).eq.isblkxder(4,is2).and.
     $        isblkder(1,is).eq.isblkxder(3,is2))then
          if(isblkder(3,is).eq.isblkxder(1,is2))then
           iad1=i4od(is2)
           do i4=0,nocc(isblkxder(4,is2))-1
            do i3=0,nocc(isblkxder(3,is2))-1
             do i2=0,nocc(isblkxder(2,is2))-1
              do i1=0,nocc(isblkxder(1,is2))-1
               iad2=imt+i3+nocc(isblkxder(3,is2))*(i4
     $              +nocc(isblkxder(4,is2))*(i1
     $              +nocc(isblkxder(1,is2))*i2))
               bc(iad2)=bc(iad2)+bc(iad1)
               iad1=iad1+1
              end do
             end do
            end do
           end do
          else if(isblkder(3,is).eq.isblkxder(2,is2))then               6d29s22
           iad1=i4od(is2)
           do i4=0,nocc(isblkxder(4,is2))-1
            do i3=0,nocc(isblkxder(3,is2))-1
             do i2=0,nocc(isblkxder(2,is2))-1
              do i1=0,nocc(isblkxder(1,is2))-1
               iad2=imt+i3+nocc(isblkxder(3,is2))*(i4
     $              +nocc(isblkxder(4,is2))*(i2                         6d29s22
     $              +nocc(isblkxder(2,is2))*i1))
               bc(iad2)=bc(iad2)+bc(iad1)
               iad1=iad1+1
              end do
             end do
            end do
           end do
          end if
         end if
         if(isblkder(3,is).eq.isblkxder(4,is2).and.
     $        isblkder(4,is).eq.isblkxder(3,is2))then
          if(isblkder(1,is).eq.isblkxder(1,is2))then
           iad1=i4od(is2)
           do i4=0,nocc(isblkxder(4,is2))-1
            do i3=0,nocc(isblkxder(3,is2))-1
             do i2=0,nocc(isblkxder(2,is2))-1
              do i1=0,nocc(isblkxder(1,is2))-1
               iad2=imt+i1+nocc(isblkxder(1,is2))*(i2
     $              +nocc(isblkxder(2,is2))*(i4
     $              +nocc(isblkxder(4,is2))*i3))
               bc(iad2)=bc(iad2)+bc(iad1)
               iad1=iad1+1
              end do
             end do
            end do
           end do
          else if(isblkder(1,is).eq.isblkxder(2,is2))then
           iad1=i4od(is2)
           do i4=0,nocc(isblkxder(4,is2))-1
            do i3=0,nocc(isblkxder(3,is2))-1
             do i2=0,nocc(isblkxder(2,is2))-1
              do i1=0,nocc(isblkxder(1,is2))-1
               iad2=imt+i2+nocc(isblkxder(2,is2))*(i1
     $              +nocc(isblkxder(1,is2))*(i4
     $              +nocc(isblkxder(4,is2))*i3))
               bc(iad2)=bc(iad2)+bc(iad1)
               iad1=iad1+1
              end do
             end do
            end do
           end do
          end if
         end if
         if(isblkder(4,is).eq.isblkxder(4,is2).and.
     $        isblkder(3,is).eq.isblkxder(3,is2))then
          if(isblkder(1,is).eq.isblkxder(1,is2))then
           iad1=i4od(is2)
           do i4=0,nocc(isblkxder(4,is2))-1
            do i3=0,nocc(isblkxder(3,is2))-1
             do i2=0,nocc(isblkxder(2,is2))-1
              do i1=0,nocc(isblkxder(1,is2))-1
               iad2=imt+i1+nocc(isblkxder(1,is2))*(i2
     $              +nocc(isblkxder(2,is2))*(i3
     $              +nocc(isblkxder(3,is2))*i4))
               bc(iad2)=bc(iad2)+bc(iad1)
               iad1=iad1+1
              end do
             end do
            end do
           end do
          else if(isblkder(1,is).eq.isblkxder(2,is2))then
           iad1=i4od(is2)
           do i4=0,nocc(isblkxder(4,is2))-1
            do i3=0,nocc(isblkxder(3,is2))-1
             do i2=0,nocc(isblkxder(2,is2))-1
              do i1=0,nocc(isblkxder(1,is2))-1
               iad2=imt+i2+nocc(isblkxder(2,is2))*(i1
     $              +nocc(isblkxder(1,is2))*(i3
     $              +nocc(isblkxder(3,is2))*i4))
               bc(iad2)=bc(iad2)+bc(iad1)
               iad1=iad1+1
              end do
             end do
            end do
           end do
          end if
         end if
        end do
       end if
      end do
      imt=im0
c
c     i want to copy imt to i4od, however i4od is pointed to by
c     nsblkxder and isblkxder while imt is point to by nsblkder and
c     isblkder. now the storage initially allocated for i4od will
c     not be smaller than that allocated for imt, but the symmetries
c     (and sizes) pointed to by isblkxder could be different from
c     those used in imt. however since memory for both were allocated
c     contiguously, we just need to reset the i4od(i),i gt 1.
c
      i4od0=i4od(1)                                                     12d8s16
      do is=1,nsblkder                                                  12d8s16
       nrow=nocc(isblkder(1,is))*nocc(isblkder(2,is))                   12d8s16
       ncol=nocc(isblkder(3,is))*nocc(isblkder(4,is))                   12d8s16
       nwds=nrow*ncol                                                   12d8s16
       if(nwds.gt.0)then                                                12d8s16
        i4od(is)=i4od0                                                  12d8s16
        do i=0,nwds-1                                                   12d8s16
         bc(i4od(is)+i)=bc(imt+i)                                       12d8s16
        end do                                                          12d8s16
        imt=imt+nwds                                                    12d8s16
        i4od0=i4od0+nwds                                                12d8s16
       end if                                                           12d8s16
      end do                                                            12d8s16
      ibcoff=im0                                                        12d8s16
      return
      end
