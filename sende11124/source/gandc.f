c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine gandc(idb,iob,idk,iok,nopenb,nopenk,barg,karg,ncsf,    5d21s20
     $     norbx,ixw1,ixw2,nnot,nab,iwpb,iwpk,ncsfmid,bc,ibc)           11d14s22
      implicit real*8 (a-h,o-z)                                         12d3s20
c
c     no, not guidance and control but rather a combinaton of gcode1
c     and cupo21.
c
      integer*8 idb,iob,idk,iok,itempc,itempd                           5d21s20
      integer barg                                                      5d21s20
      integer*1 nab(2)                                                  5d21s20
      logical lprint                                                    6d23s21
      dimension ioccd(2,2),ncsf(*)                                      5d21s20
      data icall/0/
      include "common.store"                                            5d21s20
      save icall
      icall=icall+1
      nnot=0                                                            5d21s20
      itempc=ieor(idb,idk)                                              5d21s20
      ndifd=popcnt(itempc)                                              5d21s20
      ncsfmid=-1                                                        5d21s20
      if(ndifd.le.2)then                                                5d21s20
       itempd=ieor(iob,iok)                                             5d21s20
       ndifs=popcnt(itempd)                                             5d21s20
       if(ndifs.le.2)then                                               5d21s20
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
        if(ndifs.eq.0.and.ndifd.eq.0)then                               5d21s20
         nnot=1                                                         5d21s20
        else if(ndifd.eq.1.and.ndifs.eq.2)then                          5d21s20
         if(ioccd(1,1).eq.-ioccd(1,2).or.ioccd(1,1).eq.-ioccd(2,2))     5d21s20
     $            then
          itempd=ior(iob,iok)                                           5d21s20
          nnot=2                                                        5d21s20
          if(ioccd(1,1).lt.0)then                                       5d21s20
           nab(1)=-ioccd(1,1)                                           5d21s20
           if(ioccd(1,2).ne.nab(1))then                                 5d21s20
            nab(2)=ioccd(1,2)                                           5d21s20
           else                                                         5d21s20
            nab(2)=ioccd(2,2)                                           5d21s20
           end if                                                       5d21s20
          else                                                          5d21s20
           nab(2)=ioccd(1,1)                                            5d21s20
           if(ioccd(1,2).ne.-nab(2))then                                5d21s20
            nab(1)=-ioccd(1,2)                                          5d21s20
           else                                                         5d21s20
            nab(1)=-ioccd(2,2)                                          5d21s20
           end if                                                       5d21s20
          end if                                                        5d21s20
          nx=0                                                          5d21s20
          i1o=0                                                         5d21s20
          n1o=0                                                         5d21s20
          i2o=0                                                         5d21s20
          n2o=0                                                         5d21s20
          do i=1,norbx                                                  5d21s20
           if(btest(itempd,i))then                                      5d21s20
            nx=nx+1                                                     5d21s20
            icodex=3                                                    5d21s20
            if(i.eq.iabs(ioccd(1,1)))then                               5d21s20
             if(ioccd(1,1).gt.0)then                                    5d21s20
              icodex=7                                                  5d21s20
             else                                                       5d21s20
              icodex=5                                                  5d21s20
             end if                                                     5d21s20
             i2o=icodex                                                 5d21s20
            else if(i.eq.iabs(ioccd(1,2)))then                          5d21s20
             if(ioccd(1,2).lt.0)then                                    5d21s20
              icodex=1                                                  5d21s20
             else                                                       5d21s20
              icodex=2                                                  5d21s20
             end if                                                     5d21s20
             i1o=icodex                                                 5d21s20
            else if(i.eq.iabs(ioccd(2,2)))then                          5d21s20
             if(ioccd(2,2).lt.0)then                                    5d21s20
              icodex=1                                                  5d21s20
             else                                                       5d21s20
              icodex=2                                                  5d21s20
             end if                                                     5d21s20
             i1o=icodex                                                 5d21s20
            end if                                                      5d21s20
            if(i1o.eq.0.and.icodex.ge.3)n1o=n1o+1                       5d21s20
            if(i2o.eq.0.and.icodex.ge.3)n2o=n2o+1                       5d21s20
           end if                                                       5d21s20
          end do                                                        5d21s20
          if(i1o.eq.1)then                                              5d21s20
           iicase=1
           iooa=ixw1+((nopenb*(nopenb-1))/2)+n1o                        12d3s20
           itest=ibc(iooa)                                              12d4s20
           if(ibc(itest).ne.0)then                                      12d4s20
            nusedi=ibc(itest)/2                                         12d4s20
            if(2*nusedi.ne.ibc(itest))nusedi=nusedi+1                   12d4s20
            iuse=itest+ibc(itest)+nusedi+ncsf(barg)+1                   12d4s20
            iwpb=-iuse
           else
            iwpb=ibc(iooa)+1                                             12d3s20
           end if
           iooa=ixw2+((nopenk*(nopenk+1))/2)+n2o                        12d3s20
           itest=ibc(iooa)                                              12d4s20
           if(ibc(itest).ne.0)then                                      12d4s20
            nusedi=ibc(itest)/2                                         12d4s20
            if(2*nusedi.ne.ibc(itest))nusedi=nusedi+1                   12d4s20
            iuse=itest+ibc(itest)+nusedi+ncsf(barg)+1                   12d4s20
            iwpk=-iuse
           else                                                         12d4s20
            iwpk=ibc(iooa)+1                                             12d4s20
           end if
           ncsfmid=ncsf(barg)                                           5d21s20
          else                                                          5d21s20
           iicase=2
           iooa=ixw2+((nopenb*(nopenb+1))/2)+n2o                        12d3s20
           itest=ibc(iooa)                                              12d4s20
           if(ibc(itest).ne.0)then                                      12d4s20
            iwpb=-itest                                                 12d5s20
           else                                                         12d4s20
            iwpb=ibc(iooa)+ncsf(karg)*ncsf(barg)+1                       12d4s20
           end if                                                       12d4s20
           iooa=ixw1+((nopenk*(nopenk-1))/2)+n1o                        12d3s20
           itest=ibc(iooa)                                              12d4s20
           if(ibc(itest).ne.0)then                                      12d4s20
            iwpk=-itest                                                 12d5s20
           else
            iwpk=ibc(iooa)+ncsf(karg)*ncsf(karg)+1                       12d3s20
           end if
           ncsfmid=ncsf(karg)                                           5d21s20
          end if                                                        5d21s20
         end if                                                         5d21s20
        else if(ndifd.eq.2.and.ndifs.eq.2)then                          5d21s20
         if((ioccd(1,1).eq.-ioccd(1,2).and.ioccd(2,1).eq.               5d21s20
     $        -ioccd(2,2)).or.(ioccd(2,1).eq.-ioccd(1,2).and.           5d21s20
     $            ioccd(1,1).eq.-ioccd(2,2)))                           5d21s20
     $           then                                                   5d21s20
          nnot=2                                                        5d21s20
          itempd=ior(iob,iok)                                           5d21s20
          if(ioccd(1,1).lt.0)then                                       5d21s20
           nab(1)=-ioccd(1,1)                                           5d21s20
           nab(2)=ioccd(2,1)                                            5d21s20
          else                                                          5d21s20
           nab(1)=-ioccd(2,1)                                           5d21s20
           nab(2)=ioccd(1,1)                                            5d21s20
          end if                                                        5d21s20
          nx=0                                                          5d21s20
          do i=1,norbx                                                  5d21s20
           if(btest(itempd,i))then                                      5d21s20
            nx=nx+1                                                     5d21s20
            icodex=3                                                    5d21s20
            if(i.eq.nab(1))then                                         5d21s20
             icodex=5                                                   5d21s20
             i1o=nx                                                     5d21s20
            else if(i.eq.nab(2))then                                    5d21s20
             i2o=nx                                                     5d21s20
             icodex=7                                                   5d21s20
            end if                                                      5d21s20
           end if                                                       5d21s20
          end do                                                        5d21s20
          karga=karg-1                                                  5d21s20
          iicase=3
          iooa=ixw2+((nopenb*(nopenb+1))/2)+i1o-1                       12d3s20
          itest=ibc(iooa)                                               12d4s20
          if(ibc(itest).ne.0)then                                       12d4s20
           iwpb=-itest                                                  12d5s20
          else
           iwpb=ibc(iooa)+ncsf(barg)*ncsf(karga)+1                       12d4s20
          end if                                                        12d4s20
          iooa=ixw2+((nopenk*(nopenk+1))/2)+i2o-1                       12d3s20
          itest=ibc(iooa)                                               12d4s20
          if(ibc(itest).ne.0)then                                       12d4s20
           nusedi=ibc(itest)/2                                          12d5s20
           if(2*nusedi.ne.ibc(itest))nusedi=nusedi+1                    12d5s20
           iuse=itest+ibc(itest)+nusedi+ncsf(karga)+1                   12d5s20
           iwpk=-iuse                                                   12d5s20
          else
           iwpk=ibc(iooa)+1                                              12d4s20
          end if                                                        12d4s20
          ncsfmid=ncsf(karga)                                           5d21s20
         end if                                                         5d21s20
        else if(ndifd.eq.0.and.ndifs.eq.2)then                          5d21s20
         nnot=2                                                         5d21s20
         itempd=ior(iob,iok)                                            5d21s20
         nab=0                                                          5d21s20
         if(ioccd(1,2).lt.0)then                                        5d21s20
          nab(1)=-ioccd(1,2)                                            5d21s20
          nab(2)=ioccd(2,2)                                             5d21s20
         else                                                           5d21s20
          nab(1)=-ioccd(2,2)                                            5d21s20
          nab(2)=ioccd(1,2)                                             5d21s20
         end if                                                         5d21s20
         nx=0                                                           5d21s20
         n1o=0                                                          5d21s20
         i1=0                                                           5d21s20
         n2o=0                                                          5d21s20
         i2=0                                                           5d21s20
         do i=1,norbx                                                   5d21s20
          if(btest(itempd,i))then                                       5d21s20
           nx=nx+1                                                      5d21s20
           icodex=3                                                     5d21s20
           if(i.eq.nab(1))then                                          5d21s20
            i1=1                                                        5d21s20
            icodex=1                                                    5d21s20
           else if(i.eq.nab(2))then                                     5d21s20
            i2=1                                                        5d21s20
            icodex=2                                                    5d21s20
           end if                                                       5d21s20
           if(i1.eq.0.and.icodex.eq.3)n1o=n1o+1                         5d21s20
           if(i2.eq.0.and.icodex.eq.3)n2o=n2o+1                         5d21s20
          end if                                                        5d21s20
         end do                                                         5d21s20
         iooa=ixw1+((nopenk*(nopenk-1))/2)+n1o                          12d3s20
         itest=ibc(iooa)                                                12d4s20
         if(ibc(itest).ne.0)then                                        12d4s20
          nusedi=ibc(itest)/2                                           12d5s20
          if(2*nusedi.ne.ibc(itest))nusedi=nusedi+1                     12d5s20
          iuse=itest+ibc(itest)+nusedi+ncsf(barg)+1                     12d5s20
          iwpb=-iuse                                                    12d5s20
         else
          iwpb=ibc(iooa)+1                                               12d3s20
         end if
         iooa=ixw1+((nopenk*(nopenk-1))/2)+n2o                          12d3s20
         itest=ibc(iooa)                                                12d4s20
         if(ibc(itest).ne.0)then                                         12d4s20
          iwpk=-itest                                                   12d5s20
         else
          iwpk=ibc(iooa)+ncsf(barg)*ncsf(barg)+1                         12d3s20
         end if                                                         12d4s20
         iicase=4
         ncsfmid=ncsf(barg)                                             5d21s20
        end if                                                          5d21s20
       end if                                                           5d21s20
      end if                                                            5d21s20
      return
      end
