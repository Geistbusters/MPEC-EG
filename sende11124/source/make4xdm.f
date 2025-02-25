c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine make4xdm(nsymb,nvirt,isblk4x,nx,multh)                 7d5s21
      implicit real*8 (a-h,o-z)
c
c     storage for 4 virt integrals. Since they are used fairly simply,
c     store as compact as possible.
c     abcd: i1+na*(i2+nb*(i3+nc*i4)
c     abab: ((ix*(ix-1))/2)+in, ix=i1+na*i2 ge in=i3+na*i4
c     aabb: ((i1*(i1-1))/2)+i2+naa*[((i3*(i3-1))/2)+i4]
c     aaaa: ((ix*(ix-1))/2)+in, ix=((i1*(i1-1))/2)+i2, in=((i3*(i3-1))/2+i4.
c
      common/fnd4xcm/inv4x(2,8,8,8)                                     7d4s21
      dimension nvirt(*),isblk4x(5,*),multh(8,8)                        7d4s21
      nx=0                                                              7d4s21
      do isa=1,nsymb                                                    7d4s21
       do isb=1,isa                                                     7d4s21
        itriab=((isa*(isa-1))/2)+isb
        isab=multh(isa,isb)                                             7d4s21
        if(isa.eq.isb)then                                              7d4s21
         nab=(nvirt(isa)*(nvirt(isa)+1))/2                              7d4s21
        else                                                            7d4s21
         nab=nvirt(isa)*nvirt(isb)                                      7d4s21
        end if                                                          7d4s21
        do isc=1,nsymb                                                  7d4s21
         isd=multh(isab,isc)                                            7d4s21
         if(isd.le.isc)then                                             7d4s21
          itricd=((isc*(isc-1))/2)+isd                                  7d4s21
          if(itricd.le.itriab.and.nab.gt.0)then                         7d4s21
           if(isd.eq.isc)then                                           7d4s21
            ncd=(nvirt(isc)*(nvirt(isc)+1))/2                           7d4s21
           else                                                         7d4s21
            ncd=nvirt(isc)*nvirt(isd)                                   7d4s21
           end if                                                       7d4s21
           if(itricd.eq.itriab)then                                     7d4s21
            ntot=(nab*(nab+1))/2                                        7d4s21
           else                                                         7d4s21
            ntot=nab*ncd                                                7d4s21
           end if                                                       7d4s21
           if(ntot.gt.0)then                                            7d4s21
            nx=nx+1                                                     7d4s21
            isblk4x(1,nx)=isa                                           7d4s21
            isblk4x(2,nx)=isb                                           7d4s21
            isblk4x(3,nx)=isc                                           7d4s21
            isblk4x(4,nx)=isd                                           7d4s21
    1       format(i4,4x,4i2,i4)
            inv4x(1,isa,isb,isc)=nx                                     7d4s21
            inv4x(2,isa,isb,isc)=1                                      7d4s21
            if(isa.ne.isb)then                                          7d6s21
c                      2   1   3   4
             inv4x(1,isb,isa,isc)=nx                                     7d4s21
             inv4x(2,isb,isa,isc)=2                                      7d4s21
            end if                                                      7d6s21
            if(isc.ne.isd)then                                          7d6s21
c                     1   2    4   3
             inv4x(1,isa,isb,isd)=nx                                     7d4s21
             inv4x(2,isa,isb,isd)=3                                      7d4s21
             if(isa.ne.isb)then                                         7d6s21
c                      2   1   4   3
              inv4x(1,isb,isa,isd)=nx                                     7d4s21
              inv4x(2,isb,isa,isd)=4                                      7d4s21
             end if                                                     7d6s21
            end if                                                      7d6s21
            if(isa.ne.isc.and.isd.ne.isb)then                           7d6s21
c                     3   4   1   2
             inv4x(1,isc,isd,isa)=nx                                     7d4s21
             inv4x(2,isc,isd,isa)=5                                      7d4s21
             if(isc.ne.isd)then
c                      4   3   1   2
              inv4x(1,isd,isc,isa)=nx                                     7d4s21
              inv4x(2,isd,isc,isa)=6                                      7d4s21
             end if
             if(isa.ne.isb)then
c                      3   4   2   1
              inv4x(1,isc,isd,isb)=nx                                     7d4s21
              inv4x(2,isc,isd,isb)=7                                      7d4s21
              if(isd.ne.isc)then
c                       4   3   2   1
               inv4x(1,isd,isc,isb)=nx                                     7d4s21
               inv4x(2,isd,isc,isb)=8                                      7d4s21
              end if                                                    7d6s21
             end if                                                     7d6s21
            end if                                                      7d6s21
           end if                                                       7d4s21
          end if                                                        7d4s21
         end if                                                         7d4s21
        end do                                                          7d4s21
       end do                                                           7d4s21
      end do                                                            7d4s21
      return
      end                                                               7d4s21
