c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine gandcr(idb,iob,idk,iok,nopenb,nopenk,                  9d6s21
     $     norbx,nnot,nab,icode,imap,nx,irw1,irw2,iwpb,iwpk,bc,ibc)     11d14s22
      implicit real*8 (a-h,o-z)                                         12d3s20
c
c     try to use gandc methods in spin-orbit case.
c
      integer*8 idb,iob,idk,iok,itempc,itempd                           5d21s20
      integer barg                                                      5d21s20
      integer*1 nab(2),icode(*),imap(*)                                 9d6s21
      logical lprint                                                    6d23s21
      dimension ioccd(2,2),iwpb(*),iwpk(*)                              7d29s22
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
        irwoffb=(nopenb-mod(nopenb,2))/2                                9d17s21
        irwoffk=(nopenk-mod(nopenb,2))/2                                9d17s21
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
            imap(nx)=i                                                  9d6s21
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
            icode(nx)=icodex                                            9d6s21
           end if                                                       5d21s20
          end do                                                        5d21s20
          if(i1o.eq.1)then                                              5d21s20
           iwpb(1)=irw1+irwoffb                                         9d17s21
           iwpb(2)=n1o                                                  9d16s21
           iwpb(3)=0                                                    9d16s21
           iwpb(4)=1                                                    9d17s21
           iwpk(1)=irw2+irwoffk                                         9d17s21
           iwpk(2)=n2o                                                  9d16s21
           iwpk(3)=0                                                    9d16s21
           iwpk(4)=2                                                    9d17s21
          else                                                          5d21s20
           iwpb(1)=irw2+irwoffb                                         9d17s21
           iwpb(2)=n2o                                                  9d16s21
           iwpb(3)=1                                                    9d16s21
           iwpb(4)=2                                                    9d17s21
           iwpk(1)=irw1+irwoffk                                         9d17s21
           iwpk(2)=n1o                                                  9d16s21
           iwpk(3)=1                                                    9d16s21
           iwpk(4)=1                                                    9d17s21
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
            imap(nx)=i                                                  9d6s21
            icode(nx)=icodex                                            9d6s21
           end if                                                       5d21s20
          end do                                                        5d21s20
          iwpb(1)=irw2+irwoffb                                          9d17s21
          iwpb(2)=i1o-1                                                 9d16s21
          iwpb(3)=1                                                     9d16s21
          iwpb(4)=2                                                     9d17s21
          iwpk(1)=irw2+irwoffk                                          9d17s21
          iwpk(2)=i2o-1                                                 9d16s21
          iwpk(3)=0                                                     9d16s21
          iwpk(4)=2                                                     9d17s21
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
           imap(nx)=i                                                   9d6s21
           icode(nx)=icodex                                             9d6s21
          end if                                                        5d21s20
         end do                                                         5d21s20
         iwpb(1)=irw1+irwoffk                                           9d17s21
         iwpb(2)=n1o                                                    9d16s21
         iwpb(3)=0                                                      9d16s21
         iwpb(4)=1                                                      9d17s21
         iwpk(1)=irw1+irwoffk                                           9d17s21
         iwpk(2)=n2o
         iwpk(3)=1                                                      9d16s21
         iwpk(4)=1                                                      9d17s21
        end if                                                          5d21s20
       end if                                                           5d21s20
      end if                                                            5d21s20
      return
      end
