c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine der3part(propmat,ider3,ionex,i3x,da,noc,nvirtc,ipuse,
     $     multh,isbu,itmp,bc,ibc)                                      11d14s22
      implicit real*8 (a-h,o-z)
c
c     compute sum jk Cijk daj dak contribution to lhs for d2a.
c
      include "common.store"
      include "common.hf"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension propmat(*),ionex(*),i3x(*),da(*),noc(8),nvirtc(8),
     $     ipt(8),multh(8,8),ipt2(8)
      sixteenthree=16d0/3d0
      ibcoffo=ibcoff
      ioff=0                                                            1d24s17
      ioff2=0
      do isb=1,nsymb                                                    1d24s17
       ipt(isb)=ioff                                                    1d24s17
       ipt2(isb)=ioff2
       isk=multh(isb,ipuse)                                             1d24s17
       ioff=ioff+noc(isk)*nvirtc(isb)                                   1d24s17
       ioff2=ioff2+noc(isb)*nvirtc(isb)
      end do                                                            1d24s17
      isb=isbu
       isk=multh(isb,ipuse)                                             1d25s17
       itmp=ibcoff                                                      1d24s17
       nwds=nvirtc(isb)*noc(isb)                                        1d24s17
       ibcoff=itmp+nwds                                                 1d24s17
       call enough('der3part.  1',bc,ibc)
       do i=0,nwds-1                                                    1d24s17
        bc(itmp+i)=0d0                                                  1d24s17
       end do                                                           1d24s17
       do is=1,nsdlk1                                                   1d24s17
        is1=isblk1(1,is)
        is2=isblk1(2,is)
        is3=isblk1(3,is)
        is4=isblk1(4,is)
        call ilimts(noc(is3),nvirtc(is4),mynprocg,mynowprog,
     $        il,ih,i3s,i3e,i4s,i4e)
        if(is3.eq.isb.and.multh(is2,ipuse).eq.is4)then                  1d24s17
c
c     do lqo (qo|ml) ok
c
         if(is1.eq.is2)then
          nrow=(noc(is1)*(noc(is1)+1))/2
         else
          nrow=noc(is1)*noc(is2)
         end if
         i30=i3s
         i3n=noc(isb)
         ii=ionex(is)
         do i4=i4s,i4e
          if(i4.eq.i4e)i3n=i3e
          do i3=i30,i3n
           i3m=i3-1
           if(is1.eq.is2)then
            do io=0,noc(is2)-1
             do iq=0,noc(is1)-1
              ix=max(io,iq)
              in=min(io,iq)
              iad=ii+((ix*(ix+1))/2)+in
              iadlo=ipt(is4)+i4+nvirtc(is4)*io
              val=8d0*bc(iad)*da(iadlo)
              iadiq=ipt(isb)+nvirtc(isb)*iq+1
              iadr=itmp+nvirtc(isb)*i3m
              do i=0,nvirtc(isb)-1
               bc(iadr+i)=bc(iadr+i)+val*da(iadiq+i)
              end do
             end do
            end do
           else
            do io=0,noc(is2)-1
             do iq=0,noc(is1)-1
              iad=ii+iq+noc(is1)*io
              iadlo=ipt(is4)+i4+nvirtc(is4)*io
              val=8d0*bc(iad)*da(iadlo)
              iadiq=ipt(isb)+nvirtc(isb)*iq+1
              iadr=itmp+nvirtc(isb)*i3m                                 2d1s17
              do i=0,nvirtc(isb)-1
               bc(iadr+i)=bc(iadr+i)+val*da(iadiq+i)
              end do
             end do
            end do
           end if
           ii=ii+nrow
          end do
          i30=1
         end do
        else if(is3.eq.isb.and.multh(is1,ipuse).eq.is4)then             1d25s17
         nrow=noc(is1)*noc(is2)
         i30=i3s
         i3n=noc(isb)
         ii=ionex(is)
         do i4=i4s,i4e
          if(i4.eq.i4e)i3n=i3e
          do i3=i30,i3n
           i3m=i3-1
           do iq=0,noc(is2)-1
            do io=0,noc(is1)-1
             iad=ii+io+noc(is1)*iq
             iadlo=ipt(is4)+i4+nvirtc(is4)*io
             val=8d0*bc(iad)*da(iadlo)
             iadiq=ipt(isb)+nvirtc(isb)*iq+1
             iadr=itmp+nvirtc(isb)*i3m
             do i=0,nvirtc(isb)-1
              bc(iadr+i)=bc(iadr+i)+val*da(iadiq+i)
             end do
            end do
           end do
           ii=ii+nrow
          end do
          i30=1
         end do
        end if
        if(is1.eq.isb.and.multh(is2,ipuse).eq.is4)then                  1d25s17
c
c     do lqo (mo|ql) ok
c
         if(is1.eq.is2)then
          nrow=(noc(is1)*(noc(is1)+1))/2
         else
          nrow=noc(is1)*noc(is2)
         end if
         i30=i3s
         i3n=noc(is3)
         ii=ionex(is)
         do l=i4s,i4e
          if(l.eq.i4e)i3n=i3e
          do i3=i30,i3n
           iq=i3-1
           if(is1.eq.is2)then
            do io=0,noc(is2)-1
             do m=0,noc(is1)-1
              ix=max(io,m)
              in=min(io,m)
              iad=ii+((ix*(ix+1))/2)+in
              iadlo=ipt(is4)+l+nvirtc(is4)*io
              val=8d0*bc(iad)*da(iadlo)
              iadiq=ipt(isb)+nvirtc(isb)*iq+1
              iadr=itmp+nvirtc(isb)*m
              do i=0,nvirtc(isb)-1
               bc(iadr+i)=bc(iadr+i)+val*da(iadiq+i)
              end do
             end do
            end do
           else
            do io=0,noc(is2)-1
             do m=0,noc(is1)-1
              iad=ii+m+noc(is1)*io
              iadlo=ipt(is4)+l+nvirtc(is4)*io
              val=8d0*bc(iad)*da(iadlo)
              iadiq=ipt(isb)+nvirtc(isb)*iq+1
              iadr=itmp+nvirtc(isb)*m
              do i=0,nvirtc(isb)-1
               bc(iadr+i)=bc(iadr+i)+val*da(iadiq+i)
              end do
             end do
            end do
           end if
           ii=ii+nrow
          end do
          i30=1
         end do
        else if(is2.eq.isb.and.multh(is1,ipuse).eq.is4)then             1d25s17
         nrow=noc(is1)*noc(is2)
         i30=i3s
         i3n=noc(is3)
         ii=ionex(is)
         do l=i4s,i4e
          if(l.eq.i4e)i3n=i3e
          do i3=i30,i3n
           iq=i3-1
           do m=0,noc(is2)-1
            do io=0,noc(is1)-1
             iad=ii+io+noc(is1)*m
             iadlo=ipt(is4)+l+nvirtc(is4)*io
             val=8d0*bc(iad)*da(iadlo)
             iadiq=ipt(isb)+nvirtc(isb)*iq+1
             iadr=itmp+nvirtc(isb)*m
             do i=0,nvirtc(isb)-1
              bc(iadr+i)=bc(iadr+i)+val*da(iadiq+i)
             end do
            end do
           end do
           ii=ii+nrow
          end do
          i30=1
         end do
        end if
        if(is1.eq.isb.and.multh(is3,ipuse).eq.is4)then                  1d25s17
c
c     do lqo (mq|ol) ok
c
         if(is1.eq.is2)then
          nrow=(noc(is1)*(noc(is1)+1))/2
         else
          nrow=noc(is1)*noc(is2)
         end if
         i30=i3s
         i3n=noc(is3)
         ii=ionex(is)
         do l=i4s,i4e
          if(l.eq.i4e)i3n=i3e
          do i3=i30,i3n
           io=i3-1
           if(is1.eq.is2)then
            do iq=0,noc(is2)-1
             do m=0,noc(is1)-1
              ix=max(iq,m)
              in=min(iq,m)
              iad=ii+((ix*(ix+1))/2)+in
              iadlo=ipt(is4)+l+nvirtc(is4)*io
              val=-32d0*bc(iad)*da(iadlo)
              iadiq=ipt(isb)+nvirtc(isb)*iq+1
              iadr=itmp+nvirtc(isb)*m
              do i=0,nvirtc(isb)-1
               bc(iadr+i)=bc(iadr+i)+val*da(iadiq+i)
              end do
             end do
            end do
           else
            do iq=0,noc(is2)-1
             do m=0,noc(is1)-1
              iad=ii+m+noc(is1)*iq
              iadlo=ipt(is4)+l+nvirtc(is4)*io
              val=-32d0*bc(iad)*da(iadlo)
              iadiq=ipt(isb)+nvirtc(isb)*iq+1
              iadr=itmp+nvirtc(isb)*m
              do i=0,nvirtc(isb)-1
               bc(iadr+i)=bc(iadr+i)+val*da(iadiq+i)
              end do
             end do
            end do
           end if
           ii=ii+nrow
          end do
          i30=1
         end do
        else if(is2.eq.isb.and.multh(is3,ipuse).eq.is4)then             1d25s17
         nrow=noc(is1)*noc(is2)
         i30=i3s
         i3n=noc(is3)
         ii=ionex(is)
         do l=i4s,i4e
          if(l.eq.i4e)i3n=i3e
          do i3=i30,i3n
           io=i3-1
           do m=0,noc(is2)-1
            do iq=0,noc(is1)-1
             iad=ii+iq+noc(is1)*m
             iadlo=ipt(is4)+l+nvirtc(is4)*io
             val=-32d0*bc(iad)*da(iadlo)
             iadiq=ipt(isb)+nvirtc(isb)*iq+1
             iadr=itmp+nvirtc(isb)*m
             do i=0,nvirtc(isb)-1
              bc(iadr+i)=bc(iadr+i)+val*da(iadiq+i)
             end do
            end do
           end do
           ii=ii+nrow
          end do
          i30=1
         end do
        end if
        if(is1.eq.isb.and.is4.eq.isb)then                               1d26s17
c
c     do lqo (mo|qi) dalq dalo ok
c
         if(is1.eq.is2)then
          nrow=(noc(is1)*(noc(is1)+1))/2
         else
          nrow=noc(is1)*noc(is2)
         end if
         i30=i3s
         i3n=noc(is3)
         ii=ionex(is)
         isl=multh(ipuse,is2)                                           1d26s17
         do i4=i4s,i4e
          if(i4.eq.i4e)i3n=i3e
          i=i4-1
          do i3=i30,i3n
           iq=i3-1
           if(is1.eq.is2)then
            do io=0,noc(is2)-1
             do m=0,noc(is1)-1
              ix=max(io,m)
              in=min(io,m)
              iad=ii+((ix*(ix+1))/2)+in
              iadlo=ipt(isl)+nvirtc(isl)*io+1
              iadlq=ipt(isl)+nvirtc(isl)*iq+1
              val=8d0*bc(iad)
              iadr=itmp+i+nvirtc(isb)*m
              do l=0,nvirtc(isl)-1
               bc(iadr)=bc(iadr)+val*da(iadlq+l)*da(iadlo+l)
              end do
             end do
            end do
           else
            do io=0,noc(is2)-1
             do m=0,noc(is1)-1
              iad=ii+m+noc(is1)*io
              iadlo=ipt(isl)+nvirtc(isl)*io+1
              val=8d0*bc(iad)
              iadlq=ipt(isl)+nvirtc(isl)*iq+1
              iadr=itmp+i+nvirtc(isb)*m
              do l=0,nvirtc(isl)-1
               bc(iadr)=bc(iadr)+val*da(iadlq+l)*da(iadlo+l)
              end do
             end do
            end do
           end if
           ii=ii+nrow
          end do
          i30=1
         end do
        else if(is2.eq.isb.and.is4.eq.isb)then                          1d26s17
         nrow=noc(is1)*noc(is2)
         i30=i3s
         i3n=noc(is3)
         ii=ionex(is)
         isl=multh(is1,ipuse)                                           1d26s17
         do i4=i4s,i4e
          if(i4.eq.i4e)i3n=i3e
          i=i4-1                                                        1d26s17
          do i3=i30,i3n
           iq=i3-1
           do m=0,noc(is2)-1
            do io=0,noc(is1)-1
             iad=ii+io+noc(is1)*m
             iadlo=ipt(isl)+nvirtc(isl)*io+1
             val=8d0*bc(iad)
             iadlq=ipt(isl)+nvirtc(isl)*iq+1
             iadr=itmp+i+nvirtc(isb)*m
             do l=0,nvirtc(isl)-1
              bc(iadr)=bc(iadr)+val*da(iadlq+l)*da(iadlo+l)
             end do
            end do
           end do
           ii=ii+nrow
          end do
          i30=1
         end do
        end if
        if(is3.eq.isb.and.is4.eq.isb)then                               1d26s17
c
c     do lqo (qo|mi) dalq dalo ok
c
         nrow=(noc(is1)*(noc(is1)+1))/2
         i30=i3s
         i3n=noc(is3)
         ii=ionex(is)
         isl=multh(ipuse,is2)                                           1d26s17
         do i4=i4s,i4e
          if(i4.eq.i4e)i3n=i3e
          i=i4-1
          do i3=i30,i3n
           m=i3-1
           do io=0,noc(is2)-1
            do iq=0,noc(is1)-1
             ix=max(iq,io)
             in=min(iq,io)
             iad=ii+((ix*(ix+1))/2)+in
             iadlo=ipt(isl)+nvirtc(isl)*io+1
             iadlq=ipt(isl)+nvirtc(isl)*iq+1
             val=-16d0*bc(iad)
             iadr=itmp+i+nvirtc(isb)*m
             do l=0,nvirtc(isl)-1
              bc(iadr)=bc(iadr)+val*da(iadlq+l)*da(iadlo+l)
             end do
            end do
           end do
           ii=ii+nrow
          end do
          i30=1
         end do
        end if
        if(is3.eq.isb.and.is4.eq.isb)then                               1d26s17
c
c     do jql (lj|mi) dajq dalq ok
c
         nrow=(nvirtc(is1)*(nvirtc(is1)+1))/2
         i30=i3s
         i3n=noc(is3)
         ii=i3x(is)
         isq=multh(ipuse,is2)                                           1d26s17
         do i4=i4s,i4e
          if(i4.eq.i4e)i3n=i3e
          i=i4-1
          do i3=i30,i3n
           m=i3-1
           do j=0,nvirtc(is2)-1
            do l=0,nvirtc(is1)-1
             ix=max(j,l)
             in=min(j,l)
             iad=ii+((ix*(ix+1))/2)+in
             iadlq=ipt(is1)+l+1
             iadjq=ipt(is1)+j+1
             val=16d0*bc(iad)
             iadr=itmp+i+nvirtc(isb)*m
             do iq=0,noc(isq)-1
              iqn=iq*nvirtc(is1)
              bc(iadr)=bc(iadr)+val*da(iadlq+iqn)*da(iadjq+iqn)
             end do
            end do
           end do
           ii=ii+nrow
          end do
          i30=1
         end do
        end if
        if(is1.eq.isb.and.is3.eq.isb)then                               1d26s17
c
c     do jql (ij|ml) dajq dalq ok
c
         if(is1.eq.is2)then
          nrow=(nvirtc(is1)*(nvirtc(is1)+1))/2
         else
          nrow=nvirtc(is1)*nvirtc(is2)
         end if
         i30=i3s
         i3n=noc(is3)
         ii=i3x(is)
         isq=multh(ipuse,is4)                                           1d26s17
         do i4=i4s,i4e
          if(i4.eq.i4e)i3n=i3e
          l=i4-1
          do i3=i30,i3n
           m=i3-1
           if(is1.eq.is2)then
            do j=0,nvirtc(is2)-1
             do i=0,nvirtc(is1)-1
              ix=max(j,i)
              in=min(j,i)
              iad=ii+((ix*(ix+1))/2)+in
              iadlq=ipt(is4)+l+1
              iadjq=ipt(is4)+j+1
              val=-8d0*bc(iad)
              iadr=itmp+i+nvirtc(isb)*m
              do iq=0,noc(isq)-1
               iqn=iq*nvirtc(is4)
               bc(iadr)=bc(iadr)+val*da(iadlq+iqn)*da(iadjq+iqn)
              end do
             end do
            end do
           else
            do j=0,nvirtc(is2)-1
             do i=0,nvirtc(is1)-1
              iad=ii+i+nvirtc(is1)*j
              iadlq=ipt(is4)+l+1
              iadjq=ipt(is4)+j+1
              val=-8d0*bc(iad)
              iadr=itmp+i+nvirtc(isb)*m
              do iq=0,noc(isq)-1
               iqn=iq*nvirtc(is4)
               bc(iadr)=bc(iadr)+val*da(iadlq+iqn)*da(iadjq+iqn)
              end do
             end do
            end do
           end if
           ii=ii+nrow
          end do
          i30=1
         end do
        else if(is2.eq.isb.and.is3.eq.isb)then
         nrow=nvirtc(is1)*nvirtc(is2)
         i30=i3s
         i3n=noc(is3)
         ii=i3x(is)
         isq=multh(ipuse,is4)                                           1d26s17
         do i4=i4s,i4e
          if(i4.eq.i4e)i3n=i3e
          l=i4-1
          do i3=i30,i3n
           m=i3-1
           do i=0,nvirtc(is2)-1
            do j=0,nvirtc(is1)-1
             iad=ii+j+nvirtc(is1)*i
             iadlq=ipt(is4)+l+1
             iadjq=ipt(is4)+j+1
             val=-8d0*bc(iad)
             iadr=itmp+i+nvirtc(isb)*m
             do iq=0,noc(isq)-1
              iqn=iq*nvirtc(is4)
              bc(iadr)=bc(iadr)+val*da(iadlq+iqn)*da(iadjq+iqn)
             end do
            end do
           end do
           ii=ii+nrow
          end do
          i30=1
         end do
        end if
        if(is1.eq.isb.and.multh(is3,is4).eq.ipuse)then                  1d26s17
c
c     do jql (il|qj) dajq dalm ok
c
         if(is1.eq.is2)then
          nrow=(nvirtc(is1)*(nvirtc(is1)+1))/2
         else
          nrow=nvirtc(is1)*nvirtc(is2)
         end if
         i30=i3s
         i3n=noc(is3)
         ii=i3x(is)
         do i4=i4s,i4e
          if(i4.eq.i4e)i3n=i3e
          j=i4-1
          do i3=i30,i3n
           iq=i3-1
           if(is1.eq.is2)then
            do l=0,nvirtc(is2)-1
             do i=0,nvirtc(is1)-1
              ix=max(l,i)
              in=min(l,i)
              iad=ii+((ix*(ix+1))/2)+in
              iadlm=ipt(is2)+l+1
              iadjq=ipt(is4)+j+1+nvirtc(is4)*iq
              val=32d0*bc(iad)*da(iadjq)
              iadr=itmp+i
              do m=0,noc(isb)-1
               im=m*nvirtc(isb)
               il=m*nvirtc(is2)
               bc(iadr+im)=bc(iadr+im)+val*da(iadlm+il)
              end do
             end do
            end do
           else
            do l=0,nvirtc(is2)-1
             do i=0,nvirtc(is1)-1
              iad=ii+i+nvirtc(is1)*l
              iadlm=ipt(is2)+l+1
              iadjq=ipt(is4)+j+1+nvirtc(is4)*iq
              val=32d0*bc(iad)*da(iadjq)
              iadr=itmp+i
              do m=0,noc(isb)-1
               im=m*nvirtc(isb)
               il=m*nvirtc(is2)
               bc(iadr+im)=bc(iadr+im)+val*da(iadlm+il)
              end do
             end do
            end do
           end if
           ii=ii+nrow
          end do
          i30=1
         end do
        else if(is2.eq.isb.and.multh(is3,is4).eq.ipuse)then
         nrow=nvirtc(is1)*nvirtc(is2)
         i30=i3s
         i3n=noc(is3)
         ii=i3x(is)
         do i4=i4s,i4e
          if(i4.eq.i4e)i3n=i3e
          j=i4-1
          do i3=i30,i3n
           iq=i3-1
           do i=0,nvirtc(is2)-1
            do l=0,nvirtc(is1)-1
             iad=ii+l+nvirtc(is1)*i
             iadlm=ipt(is1)+l+1
             iadjq=ipt(is4)+j+1+nvirtc(is4)*iq
             val=32d0*bc(iad)*da(iadjq)
             iadr=itmp+i
             do m=0,noc(isb)-1
              im=m*nvirtc(isb)
              il=m*nvirtc(is1)
              bc(iadr+im)=bc(iadr+im)+val*da(iadlm+il)
             end do
            end do
           end do
           ii=ii+nrow
          end do
          i30=1
         end do
        end if
        if(is1.eq.isb.and.multh(is2,is3).eq.ipuse)then                  1d26s17
c
c     do jql (ij|ql) dajq dalm ok
c
         if(is1.eq.is2)then
          nrow=(nvirtc(is1)*(nvirtc(is1)+1))/2
         else
          nrow=nvirtc(is1)*nvirtc(is2)
         end if
         i30=i3s
         i3n=noc(is3)
         ii=i3x(is)
         do i4=i4s,i4e
          if(i4.eq.i4e)i3n=i3e
          l=i4-1
          do i3=i30,i3n
           iq=i3-1
           if(is1.eq.is2)then
            do j=0,nvirtc(is2)-1
             do i=0,nvirtc(is1)-1
              ix=max(j,i)
              in=min(j,i)
              iad=ii+((ix*(ix+1))/2)+in
              iadlm=ipt(is4)+l+1
              iadjq=ipt(is2)+j+1+nvirtc(is2)*iq
              val=-8d0*bc(iad)*da(iadjq)
              iadr=itmp+i
              do m=0,noc(isb)-1
               im=m*nvirtc(isb)
               il=m*nvirtc(is4)
               bc(iadr+im)=bc(iadr+im)+val*da(iadlm+il)
              end do
             end do
            end do
           else
            do j=0,nvirtc(is2)-1
             do i=0,nvirtc(is1)-1
              iad=ii+i+nvirtc(is1)*j
              iadlm=ipt(is4)+l+1
              iadjq=ipt(is2)+j+1+nvirtc(is2)*iq
              val=-8d0*bc(iad)*da(iadjq)
              iadr=itmp+i
              do m=0,noc(isb)-1
               im=m*nvirtc(isb)
               il=m*nvirtc(is4)
               bc(iadr+im)=bc(iadr+im)+val*da(iadlm+il)
              end do
             end do
            end do
           end if
           ii=ii+nrow
          end do
          i30=1
         end do
        else if(is2.eq.isb.and.multh(is3,is1).eq.ipuse)then
         nrow=nvirtc(is1)*nvirtc(is2)
         i30=i3s
         i3n=noc(is3)
         ii=i3x(is)
         do i4=i4s,i4e
          if(i4.eq.i4e)i3n=i3e
          l=i4-1
          do i3=i30,i3n
           iq=i3-1
           do i=0,nvirtc(is2)-1
            do j=0,nvirtc(is1)-1
             iad=ii+j+nvirtc(is1)*i
             iadlm=ipt(is4)+l+1
             iadjq=ipt(is1)+j+1+nvirtc(is1)*iq
             val=-8d0*bc(iad)*da(iadjq)
             iadr=itmp+i
             do m=0,noc(isb)-1
              im=m*nvirtc(isb)
              il=m*nvirtc(is4)
              bc(iadr+im)=bc(iadr+im)+val*da(iadlm+il)
             end do
            end do
           end do
           ii=ii+nrow
          end do
          i30=1
         end do
        end if
        if(is4.eq.isb.and.multh(is2,is3).eq.ipuse)then                  1d26s17
c
c     do jql (lj|qi) dajq dalm ok
c
         if(is1.eq.is2)then
          nrow=(nvirtc(is1)*(nvirtc(is1)+1))/2
         else
          nrow=nvirtc(is1)*nvirtc(is2)
         end if
         i30=i3s
         i3n=noc(is3)
         ii=i3x(is)
         do i4=i4s,i4e
          if(i4.eq.i4e)i3n=i3e
          i=i4-1
          do i3=i30,i3n
           iq=i3-1
           if(is1.eq.is2)then
            do j=0,nvirtc(is2)-1
             do l=0,nvirtc(is1)-1
              ix=max(j,l)
              in=min(j,l)
              iad=ii+((ix*(ix+1))/2)+in
              iadlm=ipt(is1)+l+1
              iadjq=ipt(is2)+j+1+nvirtc(is2)*iq
              val=-8d0*bc(iad)*da(iadjq)
              iadr=itmp+i
              do m=0,noc(isb)-1
               im=m*nvirtc(isb)
               il=m*nvirtc(is1)
               bc(iadr+im)=bc(iadr+im)+val*da(iadlm+il)
              end do
             end do
            end do
           else
            do j=0,nvirtc(is2)-1
             do l=0,nvirtc(is1)-1
              iad=ii+l+nvirtc(is1)*j
              iadlm=ipt(is1)+l+1
              iadjq=ipt(is2)+j+1+nvirtc(is2)*iq
              val=-8d0*bc(iad)*da(iadjq)
              iadr=itmp+i
              do m=0,noc(isb)-1
               im=m*nvirtc(isb)
               il=m*nvirtc(is1)
               bc(iadr+im)=bc(iadr+im)+val*da(iadlm+il)
              end do
             end do
            end do
           end if
           ii=ii+nrow
          end do
          i30=1
         end do
        else if(is4.eq.isb.and.multh(is3,is1).eq.ipuse)then
         nrow=nvirtc(is1)*nvirtc(is2)
         i30=i3s
         i3n=noc(is3)
         ii=i3x(is)
         do i4=i4s,i4e
          if(i4.eq.i4e)i3n=i3e
          i=i4-1
          do i3=i30,i3n
           iq=i3-1
           do l=0,nvirtc(is2)-1
            do j=0,nvirtc(is1)-1
             iad=ii+j+nvirtc(is1)*l
             iadlm=ipt(is2)+l+1
             iadjq=ipt(is1)+j+1+nvirtc(is1)*iq
             val=-8d0*bc(iad)*da(iadjq)
             iadr=itmp+i
             do m=0,noc(isb)-1
              im=m*nvirtc(isb)
              il=m*nvirtc(is2)
              bc(iadr+im)=bc(iadr+im)+val*da(iadlm+il)
             end do
            end do
           end do
           ii=ii+nrow
          end do
          i30=1
         end do
        end if
       end do
       call dws_gsumf(bc(itmp),nwds)
c
c     sum lq Slq da iq da lm ok
c
       islq=multh(isb,ipuse)
       jder3=ider3+ipt2(islq)
       do iq=0,noc(islq)-1
        do l=0,nvirtc(islq)-1
         iads=jder3+l+nvirtc(islq)*iq
         do m=0,noc(isb)-1
          iadlm=ipt(islq)+l+nvirtc(islq)*m+1
          val=sixteenthree*bc(iads)*da(iadlm)
          iadiq=ipt(isb)+1+nvirtc(isb)*iq
          iadr=itmp+nvirtc(isb)*m
          do i=0,nvirtc(isb)-1
           bc(iadr+i)=bc(iadr+i)+val*da(iadiq+i)
          end do
         end do
        end do
       end do
c
c     sum lq Slm da iq da lq ok
c
       isq=multh(isb,ipuse)
       isl=multh(isq,ipuse)
       jder3=ider3+ipt2(isl)
       do iq=0,noc(isq)-1
        iadiq=ipt(isb)+1+nvirtc(isb)*iq
        do l=0,nvirtc(isl)-1
         iadlq=ipt(isl)+l+nvirtc(isl)*iq+1
         do m=0,noc(isb)-1
          iads=jder3+l+nvirtc(isl)*m
          val=sixteenthree*bc(iads)*da(iadlq)
          iadr=itmp+nvirtc(isb)*m
          do i=0,nvirtc(isb)-1
           bc(iadr+i)=bc(iadr+i)+val*da(iadiq+i)
          end do
         end do
        end do
       end do
c
c     sum lq Siq da lq da lm ok
c
       isl=multh(isb,ipuse)
       isq=isb
       jder3=ider3+ipt2(isb)
       do iq=0,noc(isq)-1
        iads=jder3+nvirtc(isb)*iq
        do l=0,nvirtc(isl)-1
         iadlq=ipt(isl)+l+nvirtc(isl)*iq+1
         do m=0,noc(isb)-1
          iadlm=ipt(isl)+l+1+nvirtc(isl)*m
          val=sixteenthree*da(iadlq)*da(iadlm)
          iadr=itmp+nvirtc(isb)*m
          do i=0,nvirtc(isb)-1
           bc(iadr+i)=bc(iadr+i)+val*bc(iads+i)
          end do
         end do
        end do
       end do
      return
      end
