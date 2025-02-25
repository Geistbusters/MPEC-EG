c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine gandc4(idb,iob,idk,iok,nopenb,nopenk,norbx,nnot,nab4,  11d14s22
     $     bc,ibc)                                                      11d14s22
      implicit real*8 (a-h,o-z)                                         11d1s22
c
c     like gandc, but just return in nnot 0 for no coupling, 1 for same,
c     2 for single and 4 for double.
c
      integer*8 idb,iob,idk,iok,itempc,itempd                           5d21s20
      dimension ioccd(4,2),nab4(2,3)                                    11d13s20
      include "common.store"                                            5d21s20
      data icall/0/
      save icall
      icall=icall+1
      nnot=0                                                            5d21s20
      nab4=0                                                            11d13s20
      itempc=ieor(idb,idk)                                              5d21s20
      ndifd=popcnt(itempc)                                              5d21s20
      ncsfmid=-1                                                        5d21s20
      if(ndifd.le.4)then                                                5d21s20
       itempd=ieor(iob,iok)                                             5d21s20
       ndifs=popcnt(itempd)                                             5d21s20
       if(ndifs.le.4)then                                               5d21s20
        if(ndifd.gt.0)then                                              5d21s20
         iq=0                                                           5d21s20
         do i=1,norbx                                                   5d21s20
          if(btest(itempc,i))then                                       5d21s20
           iq=iq+1                                                      5d21s20
           ioccd(iq,1)=i                                                5d21s20
           if(btest(idb,i))ioccd(iq,1)=-i                               5d21s20
           if(iq.eq.ndifd)go to 3011                                    5d21s20
          end if                                                        5d21s20
         end do                                                         5d21s20
 3011    continue                                                       5d21s20
        end if                                                          5d21s20
        if(ndifs.gt.0)then                                              5d21s20
         iq=0                                                           5d21s20
         do i=1,norbx                                                   5d21s20
          if(btest(itempd,i))then                                       5d21s20
           iq=iq+1                                                      5d21s20
           ioccd(iq,2)=i                                                5d21s20
           if(btest(iob,i))ioccd(iq,2)=-i                               5d21s20
           if(iq.eq.ndifs)go to 3012                                    5d21s20
          end if                                                        5d21s20
         end do                                                         5d21s20
 3012    continue                                                       5d21s20
        end if                                                          5d21s20
        if(ndifd.eq.0)then                                              10d19s20
         if(ndifs.eq.0)then                                             10d19s20
          nnot=1                                                        10d19s20
         else if(ndifs.eq.2)then                                        10d19s20
          nnot=2                                                         5d21s20
          if(ioccd(1,2).gt.0)then                                       11d13s20
           nab4(2,1)=ioccd(1,2)                                         11d27s20
           nab4(1,1)=-ioccd(2,2)                                        11d27s20
          else                                                          11d13s20
           nab4(2,1)=ioccd(2,2)                                         11d27s20
           nab4(1,1)=-ioccd(1,2)                                        11d27s20
          end if                                                        11d13s20
          nab4(2,3)=1
         else if(ndifs.eq.3.or.ndifs.eq.4)then                          10d19s20
          nnot=4                                                         10d18s20
          nab4(2,3)=2
          iq1=1                                                         11d13s20
          iq2=1                                                         11d13s20
          do is=1,ndifs                                                 11d13s20
           if(ioccd(is,2).gt.0)then                                     11d13s20
            nab4(2,iq1)=ioccd(is,2)                                     11d27s20
            iq1=iq1+1                                                   11d13s20
           else                                                         11d13s20
            nab4(1,iq2)=-ioccd(is,2)                                    11d27s20
            iq2=iq2+1                                                   11d13s20
           end if                                                       11d13s20
          end do                                                        11d13s20
         end if                                                         10d19s20
        else if(ndifd.eq.1)then                                         10d19s20
         if(ndifs.eq.1)then                                             10d19s20
          nnot=3                                                        10d19s20
          nab4(2,3)=3
         else if(ndifs.eq.2)then                                        10d19s20
          if(ioccd(1,1).eq.-ioccd(1,2).or.ioccd(1,1).eq.-ioccd(2,2))     5d21s20
     $            then
           nnot=2                                                        5d21s20
           if(ioccd(1,1).gt.0)then                                      11d13s20
            nab4(2,1)=ioccd(1,1)                                        11d27s20
            if(ioccd(1,1).eq.-ioccd(1,2))then                           11d13s20
             nab4(1,1)=-ioccd(2,2)                                      11d27s20
            else                                                        11d13s20
             nab4(1,1)=-ioccd(1,2)                                      11d27s20
            end if                                                      11d13s20
           else                                                         11d13s20
            nab4(1,1)=-ioccd(1,1)                                       11d27s20
            if(ioccd(1,1).eq.-ioccd(1,2))then                           11d13s20
             nab4(2,1)=ioccd(2,2)                                       11d27s20
            else                                                        11d13s20
             nab4(2,1)=ioccd(1,2)                                       11d27s20
            end if                                                      11d13s20
           end if                                                       11d13s20
           nab4(2,3)=4
          else                                                          10d19s20
           nnot=3                                                       10d19s20
           if(ioccd(1,1).gt.0)then                                      11d13s20
            nab4(2,1)=ioccd(1,1)                                        11d27s20
            nab4(2,2)=ioccd(1,1)                                        11d27s20
            nab4(1,1)=-ioccd(1,2)                                       11d27s20
            nab4(1,2)=-ioccd(2,2)                                       11d27s20
           else                                                         11d13s20
            nab4(1,1)=-ioccd(1,1)                                       11d27s20
            nab4(1,2)=-ioccd(1,1)                                       11d27s20
            nab4(2,1)=ioccd(1,2)                                        11d27s20
            nab4(2,2)=ioccd(2,2)                                        11d27s20
           end if                                                       11d13s20
           nab4(2,3)=5
          end if                                                        10d19s20
         else if(ndifs.eq.3)then                                        10d19s20
          if(-ioccd(1,1).eq.ioccd(1,2).or.-ioccd(1,1).eq.ioccd(2,2).or.  10d19s20
     $      -ioccd(1,1).eq.ioccd(3,2))nnot=4                            10d19s20
         else if(ndifs.eq.4)then                                        10d19s20
          if(-ioccd(1,1).eq.ioccd(1,2).or.                              10d19s20
     $      -ioccd(1,1).eq.ioccd(2,2).or.                               10d19s20
     $      -ioccd(1,1).eq.ioccd(3,2).or.                               10d19s20
     $         -ioccd(1,1).eq.ioccd(4,2))then
           nnot=4                                                       11d13s20
           if(ioccd(1,1).gt.0)then                                      11d13s20
            nab4(2,1)=ioccd(1,1)                                        11d27s20
            iq1=2                                                       11d13s20
            iq2=1                                                       11d13s20
           else                                                         11d13s20
            nab4(1,1)=-ioccd(1,1)                                       11d27s20
            iq1=1                                                       11d13s20
            iq2=2                                                       11d13s20
           end if                                                       11d13s20
           do is=1,4                                                    11d13s20
            if(ioccd(1,1).ne.-ioccd(is,2))then                          11d13s20
             if(ioccd(is,2).gt.0)then                                   11d13s20
              nab4(2,iq1)=ioccd(is,2)                                   11d27s20
              iq1=iq1+1                                                 11d13s20
             else                                                       11d13s20
              nab4(1,iq2)=-ioccd(is,2)                                  11d27s20
              iq2=iq2+1                                                 11d13s20
             end if                                                     11d13s20
            end if                                                      11d13s20
           end do                                                       11d13s20
           nab4(2,3)=6
          end if
         end if                                                         10d19s20
        else if(ndifd.eq.2)then                                         10d19s20
         if(ndifs.eq.0)then                                             10d19s20
          nnot=3                                                        10d25s20
          if(ioccd(1,1).gt.0)then                                       12d15s20
           nab4(2,1)=ioccd(1,1)                                         12d15s20
           nab4(1,1)=-ioccd(2,1)                                         12d15s20
          else                                                          12d15s20
           nab4(2,1)=ioccd(2,1)                                         12d15s20
           nab4(1,1)=-ioccd(1,1)                                         12d15s20
          end if                                                        12d15s20
          nab4(2,2)=nab4(2,1)                                           12d15s20
          nab4(1,2)=nab4(1,1)                                           12d15s20
          nab4(2,3)=7
         else if(ndifs.eq.2)then                                        10d19s20
          if((ioccd(1,1).eq.-ioccd(1,2).and.ioccd(2,1).eq.               5d21s20
     $        -ioccd(2,2)).or.(ioccd(2,1).eq.-ioccd(1,2).and.           5d21s20
     $            ioccd(1,1).eq.-ioccd(2,2)))                           5d21s20
     $           then                                                   5d21s20
           nnot=2                                                        5d21s20
           if(ioccd(1,1).gt.0)then                                      11d13s20
            nab4(2,1)=ioccd(1,1)                                        11d27s20
            nab4(1,1)=-ioccd(2,1)                                       11d27s20
           else                                                         11d13s20
            nab4(2,1)=ioccd(2,1)                                        11d27s20
            nab4(1,1)=-ioccd(1,1)                                       11d27s20
           end if                                                       11d13s20
           nab4(2,3)=8
          else                                                           10d18s20
           if(-ioccd(1,1).eq.ioccd(1,2).or.-ioccd(1,1).eq.ioccd(2,2).or.10d19s20
     $        -ioccd(2,1).eq.ioccd(1,2).or.-ioccd(2,1).eq.ioccd(2,2))   10d19s20
     $          then
            if(-ioccd(1,1).eq.ioccd(1,2).or.-ioccd(1,1).eq.ioccd(2,2))  11d13s20
     $           then                                                   11d13s20
             if(ioccd(1,1).gt.0)then                                    11d13s20
              nab4(2,1)=ioccd(1,1)                                      11d27s20
              iq1=2                                                     11d13s20
              iq2=1                                                     11d13s20
             else                                                       11d13s20
              nab4(1,1)=-ioccd(1,1)                                     11d27s20
              iq1=1                                                     11d13s20
              iq2=2                                                     11d13s20
             end if                                                     11d13s20
             if(-ioccd(1,1).eq.ioccd(1,2))then                          11d13s20
              if(ioccd(2,2).gt.0)then                                   11d13s20
               nab4(2,iq1)=ioccd(2,2)                                   11d27s20
               iq1=iq1+1                                                11d13s20
              else                                                      11d13s20
               nab4(1,iq2)=-ioccd(2,2)                                  11d27s20
               iq2=iq2+1                                                11d13s20
              end if                                                    11d13s20
             else                                                       11d13s20
              if(ioccd(1,2).gt.0)then                                   11d13s20
               nab4(2,iq1)=ioccd(1,2)                                   11d27s20
               iq1=iq1+1                                                11d13s20
              else                                                      11d13s20
               nab4(1,iq2)=-ioccd(1,2)                                  11d27s20
               iq2=iq2+1                                                11d13s20
              end if                                                    11d13s20
             end if                                                     11d13s20
             if(ioccd(2,1).gt.0)then                                    11d13s20
              nab4(2,iq1)=ioccd(2,1)                                    11d27s20
              nab4(2,iq1+1)=ioccd(2,1)                                  11d27s20
             else                                                       11d13s20
              nab4(1,iq2)=-ioccd(2,1)                                   11d27s20
              nab4(1,iq2+1)=-ioccd(2,1)                                 11d27s20
             end if                                                     11d13s20
            else                                                        11d13s20
             if(ioccd(2,1).gt.0)then                                    11d13s20
              nab4(2,1)=ioccd(2,1)                                      11d27s20
              iq1=2                                                     11d13s20
              iq2=1                                                     11d13s20
             else                                                       11d13s20
              nab4(1,1)=-ioccd(2,1)                                     11d27s20
              iq1=1                                                     11d13s20
              iq2=2                                                     11d13s20
             end if                                                     11d13s20
             if(-ioccd(2,1).eq.ioccd(1,2))then                          11d13s20
              if(ioccd(2,2).gt.0)then                                   11d13s20
               nab4(2,iq1)=ioccd(2,2)                                   11d27s20
               iq1=iq1+1                                                11d13s20
              else                                                      11d13s20
               nab4(1,iq2)=-ioccd(2,2)                                  11d27s20
               iq2=iq2+1                                                11d13s20
              end if                                                    11d13s20
             else                                                       11d13s20
              if(ioccd(1,2).gt.0)then                                   11d13s20
               nab4(2,iq1)=ioccd(1,2)                                   11d27s20
               iq1=iq1+1                                                11d13s20
              else                                                      11d13s20
               nab4(1,iq2)=-ioccd(1,2)                                  11d27s20
               iq2=iq2+1                                                11d13s20
              end if                                                    11d13s20
             end if                                                     11d13s20
             if(ioccd(1,1).gt.0)then                                    11d13s20
              nab4(2,iq1)=ioccd(1,1)                                    11d27s20
              nab4(2,iq1+1)=ioccd(1,1)                                  11d27s20
             else                                                       11d13s20
              nab4(1,iq2)=-ioccd(1,1)                                   11d27s20
              nab4(1,iq2+1)=-ioccd(1,1)                                 11d27s20
             end if                                                     11d13s20
            end if                                                      11d13s20
            nnot=3                                                      10d19s20
            nab4(2,3)=9
           end if                                                       11d13s20
          end if                                                        10d19s20
         else if(ndifs.eq.3)then                                        10d19s20
          if(ioccd(1,1).eq.-ioccd(1,2))then
           if(ioccd(2,1).eq.-ioccd(2,2).or.
     $          ioccd(2,1).eq.-ioccd(3,2))then
            nnot=4
            nab4(2,3)=10
           end if
          else if(ioccd(1,1).eq.-ioccd(2,2).and.
     $         ioccd(2,1).eq.-ioccd(3,2))then
           nab4(2,3)=11
           nnot=4
          end if                                                         10d18s20
         else if(ndifs.eq.4)then                                        10d19s20
          if(-ioccd(1,1).eq.ioccd(1,2))then                             10d19s20
           if(ioccd(1,1).gt.0)then                                      11d13s20
            nab4(2,1)=ioccd(1,1)                                        11d27s20
            iq1=2                                                       11d13s20
            iq2=1                                                       11d13s20
           else                                                         11d13s20
            nab4(1,1)=-ioccd(1,1)                                       11d27s20
            iq1=1                                                       11d13s20
            iq2=2                                                       11d13s20
           end if                                                       11d13s20
           if(-ioccd(2,1).eq.ioccd(2,2).or.                             10d19s20
     $       -ioccd(2,1).eq.ioccd(3,2).or.                              10d19s20
     $         -ioccd(2,1).eq.ioccd(4,2))then
            if(ioccd(2,1).gt.0)then                                     11d13s20
             nab4(2,iq1)=ioccd(2,1)                                     11d27s20
             iq1=iq1+1                                                  11d13s20
            else                                                        11d13s20
             nab4(1,iq2)=-ioccd(2,1)                                    11d27s20
             iq2=iq2+1                                                  11d13s20
            end if                                                      11d13s20
            if(-ioccd(2,1).eq.ioccd(2,2))then                           11d13s20
             if(ioccd(3,2).gt.0)then                                    11d13s20
              nab4(2,iq1)=ioccd(3,2)                                    11d27s20
              iq1=iq1+1                                                 11d13s20
             else                                                       11d13s20
              nab4(1,iq2)=-ioccd(3,2)                                   11d27s20
              iq2=iq2+1                                                 11d13s20
             end if                                                     11d13s20
             if(ioccd(4,2).gt.0)then                                    11d13s20
              nab4(2,iq1)=ioccd(4,2)                                    11d27s20
             else                                                       11d13s20
              nab4(1,iq2)=-ioccd(4,2)                                   11d27s20
             end if                                                     11d13s20
            else if(-ioccd(2,1).eq.ioccd(3,2))then                      11d13s20
             if(ioccd(2,2).gt.0)then                                    11d13s20
              nab4(2,iq1)=ioccd(2,2)                                    11d27s20
              iq1=iq1+1                                                 11d13s20
             else                                                       11d13s20
              nab4(1,iq2)=-ioccd(2,2)                                   11d27s20
              iq2=iq2+1                                                 11d13s20
             end if                                                     11d13s20
             if(ioccd(4,2).gt.0)then                                    11d13s20
              nab4(2,iq1)=ioccd(4,2)                                    11d27s20
             else                                                       11d13s20
              nab4(1,iq2)=-ioccd(4,2)                                   11d27s20
             end if                                                     11d13s20
            else if(-ioccd(2,1).eq.ioccd(4,2))then                      11d13s20
             if(ioccd(2,2).gt.0)then                                    11d13s20
              nab4(2,iq1)=ioccd(2,2)                                    11d27s20
              iq1=iq1+1                                                 11d13s20
             else                                                       11d13s20
              nab4(1,iq2)=-ioccd(2,2)                                   11d27s20
              iq2=iq2+1                                                 11d13s20
             end if                                                     11d13s20
             if(ioccd(3,2).gt.0)then                                    11d13s20
              nab4(2,iq1)=ioccd(3,2)                                    11d27s20
             else                                                       11d13s20
              nab4(1,iq2)=-ioccd(3,2)                                   11d27s20
             end if                                                     11d13s20
            end if                                                      11d13s20
            nnot=4                                                      11d13s20
            nab4(2,3)=12
           end if
          else if(-ioccd(1,1).eq.ioccd(2,2))then                        10d19s20
           if(ioccd(1,1).gt.0)then                                      11d13s20
            nab4(2,1)=ioccd(1,1)                                        11d27s20
            iq1=2                                                       11d13s20
            iq2=1                                                       11d13s20
           else                                                         11d13s20
            nab4(1,1)=-ioccd(1,1)                                       11d27s20
            iq1=1                                                       11d13s20
            iq2=2                                                       11d13s20
           end if                                                       11d13s20
           if(ioccd(1,2).gt.0)then                                      11d13s20
            nab4(2,iq1)=ioccd(1,2)                                      11d27s20
            iq1=iq1+1                                                   11d13s20
           else                                                         11d13s20
            nab4(1,iq2)=-ioccd(1,2)                                     11d27s20
            iq2=iq2+1                                                   11d13s20
           end if                                                       11d13s20
           if(-ioccd(2,1).eq.ioccd(3,2).or.                             10d19s20
     $          -ioccd(2,1).eq.ioccd(4,2))then                          11d13s20
            if(ioccd(2,1).gt.0)then                                     11d13s20
             nab4(2,iq1)=ioccd(2,1)                                     11d27s20
             iq1=iq1+1                                                  11d13s20
            else                                                        11d13s20
             nab4(1,iq2)=-ioccd(2,1)                                    11d27s20
             iq2=iq2+1                                                  11d13s20
            end if                                                      11d13s20
            if(-ioccd(2,1).eq.ioccd(3,2))then                           11d13s20
             if(ioccd(4,2).gt.0)then                                    11d13s20
              nab4(2,iq1)=ioccd(4,2)                                    11d27s20
             else                                                       11d13s20
              nab4(1,iq2)=-ioccd(4,2)                                   11d27s20
             end if                                                     11d13s20
             nnot=4                                                      11d13s20
             nab4(2,3)=22
            else if(-ioccd(2,1).eq.ioccd(4,2))then                           11d13s20
             if(ioccd(3,2).gt.0)then                                    11d13s20
              nab4(2,iq1)=ioccd(3,2)                                    11d27s20
             else                                                       11d13s20
              nab4(1,iq2)=-ioccd(3,2)                                   11d27s20
             end if                                                     11d13s20
             nnot=4                                                      11d13s20
             nab4(2,3)=13
            end if                                                      11d13s20
           end if                                                       11d13s20
          else if(-ioccd(1,1).eq.ioccd(3,2).and.                        10d19s20
     $         -ioccd(2,1).eq.ioccd(4,2))then                           10d19s20
           if(ioccd(1,1).gt.0)then                                      11d13s20
            nab4(2,1)=ioccd(1,1)                                        11d27s20
            iq1=2                                                       11d13s20
            iq2=1                                                       11d13s20
           else                                                         11d13s20
            nab4(1,1)=-ioccd(1,1)                                       11d27s20
            iq1=1                                                       11d13s20
            iq2=2                                                       11d13s20
           end if                                                       11d13s20
           if(ioccd(1,2).gt.0)then                                      11d13s20
            nab4(2,iq1)=ioccd(1,2)                                      11d27s20
            iq1=iq1+1                                                   11d13s20
           else                                                         11d13s20
            nab4(1,iq2)=-ioccd(1,2)                                     11d27s20
            iq2=iq2+1                                                   11d13s20
           end if                                                       11d13s20
           if(ioccd(2,1).gt.0)then                                      11d13s20
            nab4(2,iq1)=ioccd(2,1)                                      11d27s20
            iq1=iq1+1                                                   11d13s20
           else                                                         11d13s20
            nab4(1,iq2)=-ioccd(2,1)                                     11d27s20
            iq2=iq2+1                                                   11d13s20
           end if                                                       11d13s20
           if(ioccd(2,2).gt.0)then                                      2d25s21
            nab4(2,iq1)=ioccd(2,2)                                      2d25s21
           else                                                         2d25s21
            nab4(1,iq2)=-ioccd(2,2 )                                    2d25s21
           end if                                                       2d25s21
           nnot=4                                                       10d19s20
           nab4(2,3)=14
          end if                                                        10d19s20
         end if                                                         10d19s20
        else if(ndifd.eq.3)then                                         10d19s20
         if(ndifs.eq.0)then                                             10d19s20
          nnot=4                                                        10d19s20
          nab4(2,3)=15
         else if(ndifs.eq.2)then                                        10d19s20
          if(-ioccd(1,2).eq.ioccd(1,1))then                             10d19s20
           if(ioccd(1,1).gt.0)then                                      11d13s20
            nab4(2,1)=ioccd(1,1)                                        11d27s20
            iq1=2                                                       11d13s20
            iq2=1                                                       11d13s20
           else                                                         11d13s20
            nab4(1,1)=-ioccd(1,1)                                       11d27s20
            iq1=1                                                       11d13s20
            iq2=2                                                       11d13s20
           end if                                                       11d13s20
           if(-ioccd(2,2).eq.ioccd(2,1).or.-ioccd(2,2).eq.ioccd(3,1))
     $         then                                                     11d13s20
            if(-ioccd(2,2).eq.ioccd(2,1))then                           11d13s20
             if(ioccd(2,1).gt.0)then                                    11d13s20
              nab4(2,iq1)=ioccd(2,1)                                    11d27s20
              iq1=iq1+1                                                 11d13s20
             else                                                       11d13s20
              nab4(1,iq2)=-ioccd(2,1)                                   11d27s20
              iq2=iq2+1                                                 11d13s20
             end if                                                     11d13s20
             if(ioccd(3,1).gt.0)then                                    11d13s20
              nab4(2,iq1)=ioccd(3,1)                                    11d27s20
              nab4(2,iq1+1)=ioccd(3,1)                                  11d27s20
             else                                                       11d13s20
              nab4(1,iq2)=-ioccd(3,1)                                   11d27s20
              nab4(1,iq2+1)=-ioccd(3,1)                                 11d27s20
             end if                                                     11d13s20
            else                                                        11d13s20
             if(ioccd(3,1).gt.0)then                                    11d13s20
              nab4(2,iq1)=ioccd(3,1)                                    11d27s20
              iq1=iq1+1                                                 11d13s20
             else                                                       11d13s20
              nab4(1,iq2)=-ioccd(3,1)                                   11d27s20
              iq2=iq2+1                                                 11d13s20
             end if                                                     11d13s20
             if(ioccd(2,1).gt.0)then                                    11d13s20
              nab4(2,iq1)=ioccd(2,1)                                    11d27s20
              nab4(2,iq1+1)=ioccd(2,1)                                  11d27s20
             else                                                       11d13s20
              nab4(1,iq2)=-ioccd(2,1)                                   11d27s20
              nab4(1,iq2+1)=-ioccd(2,1)                                 11d27s20
             end if                                                     11d13s20
            end if                                                      11d13s20
            nnot=3                                                      11d13s20
            nab4(2,3)=16
           end if                                                       11d13s20
          else if(-ioccd(1,2).eq.ioccd(2,1).and.                        10d19s20
     $         -ioccd(2,2).eq.ioccd(3,1))then                           10d19s20
           nnot=3                                                       10d19s20
           if(ioccd(1,1).gt.0)then                                      11d13s20
            nab4(2,1)=ioccd(1,1)                                        11d27s20
            nab4(2,2)=ioccd(1,1)                                        11d27s20
           else                                                         11d13s20
            nab4(1,1)=-ioccd(1,1)                                       11d27s20
            nab4(1,2)=-ioccd(1,1)                                       11d27s20
           end if                                                       11d13s20
           if(ioccd(2,1).gt.0)then                                      11d13s20
            nab4(2,1)=ioccd(2,1)                                        11d27s20
            nab4(2,2)=ioccd(3,1)                                        11d27s20
           else                                                         11d13s20
            nab4(1,1)=-ioccd(2,1)                                       11d27s20
            nab4(1,2)=-ioccd(3,1)                                       11d27s20
           end if                                                       11d13s20
           nab4(2,3)=21
          end if                                                        10d19s20
         else if(ndifs.eq.4)then                                        10d19s20
          if(-ioccd(1,1).eq.ioccd(1,2))then                             10d19s20
           if(ioccd(1,1).gt.0)then                                      11d13s20
            nab4(2,1)=ioccd(1,1)                                        11d27s20
            iq1=2                                                       11d13s20
            iq2=1
           else                                                         11d13s20
            nab4(1,1)=-ioccd(1,1)                                       11d27s20
            iq2=2                                                       11d13s20
            iq1=1
           end if                                                       11d13s20
           if(-ioccd(2,1).eq.ioccd(2,2))then                            10d19s20
            if(ioccd(2,1).gt.0)then                                     11d13s20
             nab4(2,iq1)=ioccd(2,1)                                     11d27s20
             iq1=iq1+1                                                  11d13s20
            else                                                        11d13s20
             nab4(1,iq2)=-ioccd(2,1)                                    11d27s20
             iq2=iq2+1                                                  11d13s20
            end if                                                      11d13s20
            if(-ioccd(3,1).eq.ioccd(3,2).or.                            10d19s20
     $          -ioccd(3,1).eq.ioccd(4,2))then                          11d13s20
             if(ioccd(3,1).gt.0)then                                    11d13s20
              nab4(2,iq1)=ioccd(3,1)                                    11d27s20
              iq1=iq1+1                                                 2d25s21
              if(-ioccd(3,1).eq.ioccd(3,2))then                         11d13s20
               if(ioccd(4,2).gt.0)then                                  2d25s21
                nab4(2,iq1)=ioccd(4,2)                                  2d25s21
               else                                                     2d25s21
                nab4(1,iq2)=-ioccd(4,2)                                  11d27s20
               end if                                                   2d25s21
              else                                                      11d13s20
               if(ioccd(3,2).gt.0)then                                  2d25s21
                nab4(2,iq1)=ioccd(3,2)                                  11d13s20
               else                                                     2d25s21
                nab4(1,iq2)=ioccd(3,2)                                  2d25s21
               end if                                                   2d25s21
              end if                                                    11d13s20
             else                                                       11d13s20
              nab4(1,iq2)=-ioccd(3,1)                                   11d27s20
              iq2=iq2+1                                                 2d25s21
              if(-ioccd(3,1).eq.ioccd(3,2))then                         11d13s20
               nab4(2,iq1)=ioccd(4,2)                                   11d27s20
              else                                                      11d13s20
               nab4(2,iq1)=ioccd(3,2)                                   11d27s20
              end if                                                    11d13s20
             end if                                                     11d13s20
             nab4(2,3)=17                                               11d13s20
             nnot=4                                                     11d13s20
            end if                                                      11d13s20
           else if(-ioccd(2,1).eq.ioccd(3,2).and.                       10d19s20
     $           -ioccd(3,1).eq.ioccd(4,2))then                         10d19s20
            if(ioccd(2,1).gt.0)then                                     11d13s20
             nab4(2,iq1)=ioccd(2,1)                                     11d27s20
             iq1=iq1+1                                                  11d13s20
            else                                                        11d13s20
             nab4(1,iq2)=-ioccd(2,1)                                    11d27s20
             iq2=iq2+1                                                  11d13s20
            end if                                                      11d13s20
            if(ioccd(3,1).gt.0)then                                     11d13s20
             nab4(2,iq1)=ioccd(3,1)                                     11d27s20
             iq1=iq1+1                                                  2d25s21
             nab4(2,iq1)=ioccd(2,2)                                     2d25s21
            else                                                        11d13s20
             nab4(1,iq2)=-ioccd(3,1)                                    11d27s20
             nab4(2,iq1)=ioccd(2,2)                                     11d27s20
            end if                                                      11d13s20
            nnot=4                                                      10d18s20
            nab4(2,3)=18                                                11d13s20
           end if                                                       10d18s20
          else if(-ioccd(1,1).eq.ioccd(2,2).and.                        10d19s20
     $          -ioccd(2,1).eq.ioccd(3,2).and.                          10d19s20
     $          -ioccd(3,1).eq.ioccd(4,2))then                          10d19s20
           nnot=4                                                       10d18s20
           nab4(2,3)=19                                                 11d13s20
           if(ioccd(1,1).gt.0)then                                      11d13s20
            nab4(2,1)=ioccd(1,1)                                        11d27s20
            iq1=2                                                       11d13s20
            iq2=1                                                       11d13s20
           else                                                         11d13s20
            nab4(1,1)=-ioccd(1,1)                                       11d27s20
            iq1=1                                                       11d13s20
            iq2=2                                                       11d13s20
           end if                                                       11d13s20
           if(ioccd(2,1).gt.0)then                                      11d13s20
            nab4(2,iq1)=ioccd(2,1)                                      11d27s20
            iq1=iq1+1                                                   11d13s20
           else                                                         11d13s20
            nab4(1,iq2)=-ioccd(2,1)                                     11d27s20
            iq2=iq2+1                                                   11d13s20
           end if                                                       11d13s20
           if(ioccd(3,1).gt.0)then                                      11d13s20
            nab4(2,iq1)=ioccd(3,1)                                      11d27s20
            iq1=iq1+1                                                   11d13s20
           else                                                         11d13s20
            nab4(1,iq2)=-ioccd(3,1)                                     11d27s20
            iq2=iq2+1                                                   11d13s20
           end if                                                       11d13s20
           if(ioccd(1,2).gt.0)then                                      11d13s20
            nab4(2,iq1)=ioccd(1,2)                                      11d27s20
           else                                                         11d13s20
            nab4(1,iq2)=-ioccd(1,2)                                     11d27s20
           end if                                                       11d13s20
          end if                                                        10d18s20
         end if                                                         10d19s20
        else if(ndifd.eq.4.and.ndifs.eq.4)then                          11d13s20
         if(ioccd(1,1).eq.-ioccd(1,2).and.ioccd(2,1).eq.-ioccd(2,2).and.11d13s20
     $      ioccd(3,1).eq.-ioccd(3,2).and.ioccd(4,1).eq.-ioccd(4,2))then11d13s20
          nnot=4                                                        11d13s20
          iq1=1                                                         11d13s20
          iq2=1                                                         11d13s20
          do ii=1,4                                                     11d13s20
           if(ioccd(ii,1).gt.0)then                                     11d13s20
            nab4(2,iq1)=ioccd(ii,1)                                     11d27s20
            iq1=iq1+1                                                   11d13s20
           else                                                          11d13s20
            nab4(1,iq2)=-ioccd(ii,1)                                    11d27s20
            iq2=iq2+1                                                   11d13s20
           end if                                                        11d13s20
          end do                                                        11d13s20
          nab4(2,3)=20                                                  11d13s20
         end if                                                         11d13s20
        end if                                                          10d19s20
       end if                                                           5d21s20
      end if                                                            5d21s20
      return
      end
