c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine paraerikbk4(natom,ngaus,ibdat,isym,iapair,ibstor,      12d20s21
     $     isstor,multh,ascale,nbasisp,iaddr,nfcn,nbasispc,iorb,nrootu, 12d20s21
     $     nfdat,isymw,isymwk,gd,idorel,itu,ntype,srh,bc,ibc)           11d14s22
c
c     form 4v Hdd*Vd calculation in so basis.                           11d21s18
c     input
c
      implicit real*8 (a-h,o-z)
      include "common.hf"
      include "common.store"
      include "common.basis"
      include "common.mrci"                                             12d12s18
      integer*8 ibstor,isstor                                           5d10s10
      dimension isym(3,1),iapair(3,1),ibstor(1),isstor(1),              11d23s18
     $     carta(3),cartb(3),cartc(3),cartd(3),ieraw(8,8),itmpab(8),    12d7s18
     $     multh(8,8),iptoh(8,8),nbasispc(8),ntmpab(8),                 12d7s18
     $     nbasisp(8),iaddr(36,*),nfcn(36,3),ncd(3,8,2),iorb(8),        1d13s23
     $     gd(*),ispino(4,2,4),igaunt(4,3,4,4),spino(4,4),icoul(4,2,4), 12d20s21
     $     ibrt(4,11,4,4),iorbt(8),nfdat(5,4,*),fmulx(4)                12d29s21
      dimension nvprt(8),nsz(2),nfcnb(8)                                1d25s22
      common/timerocm/tovr,telapo(15)                                   4d26s18
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      data fmulx/4*1d0/                                                 12d29s21
      data loopx/1000000/
      data icall/0/
      data ((icoul(j,i,2),j=1,4),i=1,2)/2,-2,0,0, 2,2,0,0/              12d20s21
      data ((icoul(j,i,3),j=1,4),i=1,2)/2,-2,0,0, 0,0,2,2/              12d20s21
      data ((icoul(j,i,4),j=1,4),i=1,2)/2, 2,0,0, 0,0,2,2/              12d20s21
c     vtilde spins:
c
c     integral spins:
c iifmx for operator            1         j                k
c           1   0   1   1   0 = +--+ = table col 4 ++-- = --++ col 2
c           2   1   0   1   0 = -+-+ = table col 3 -+-+  col 3
c           3   1   1   0   0 = --++ = table col 2 -++- = +--+ col 4
c           4   0   0   0   0 = ++++ = table col 1 ++++ col 1
c iifmx for operator            2         j                k
c           1   0   0   0   0 = ++++ = table col 1 ++++  col 1
c           2   1   1   0   0 = --++ = table col 2 -++- = +--+ col 4 -
c           3   0   1   1   0 = +--+ = table col 4 ++-- = --++ col 2 -
c           4   1   0   1   0 = -+-+ = table col 3 -+-+ col 3
c iifmx for operator            3         j                k
c           1   0   0   1   0 = ++-+ = table col 3 ++-+ col 3
c           2   1   1   1   0 = ---+ = table col 4 -+-- = +-++ col 2 -
c           3   0   1   0   0 = +-++ = table col 2 +++- = ---+ col 4 -
c           4   1   0   0   0 = -+++ = table col 1 -+++ col 1
c iifmx for operator            4         j                k
c           1   0   0   1   0 = ++-+ = table col 3  ++-+ = table col 3
c           2   1   1   1   0 = ---+ = table col 4  -+-- = +-++ table col 2
c           3   0   1   0   0 = +-++ = table col 2  +++- = ---+ table col 4
c           4   1   0   0   0 = -+++ = table col 1  -+++ = table col 1
      data ispino/4,3,2,1, 2,3,4,1, 1,2,4,3, 1,4,2,3,                   12d20s21
     $     3,4,2,1, 3,2,4,1, 3,4,2,1, 3,2,4,1/                          1d24s22
      data spino/4*1d0, 1d0,2*-1d0,1d0, 1d0,2*-1d0,1d0, 4*1d0/          1d24s22
      data (((igaunt(j,i,k,1),j=1,4),i=1,3),k=1,4)                      12d20s21
     $     /0,4,2,-2, 0,4,-2,-2, 2,2,0,-4,
     $     -4,0,2,-2, 4,0,2,2, 2,2,0,4,
     $     -4,0,2,-2, 4,0,2,2, 2,2,0,4,
     $     0,4,2,-2, 0,4,-2,-2, 2,2,0,-4/
      data (((igaunt(j,i,k,2),j=1,4),i=1,2),k=1,4)                      12d20s21
     $     /0,-4,2,-2, 0,4,2,2, 4,0,2,-2, -4,0,2,2,                     12d20s21
     $     -4,0,2,-2, 4,0,2,2, 0,4,2,-2, 0,-4,2,2/                      12d20s21
      data (((igaunt(j,i,k,3),j=1,4),i=1,2),k=1,4)                      12d20s21
     $     /2,-2,0,-4, 0,4,2,2, 2,-2,4,0, 0,4,-2,-2,                    12d20s21
     $     -2,2,0,-4, 4,0,2,2, -2,2,4,0, 4,0,-2,-2/                     12d20s21
      data (((igaunt(j,i,k,4),j=1,4),i=1,2),k=1,4)                      12d20s21
     $     /2,2,0,-4, 0,-4,2,2, 2,2,4,0, 0,4,2,2,                       12d20s21
     $     -2,-2,0,-4, 4,0,2,2, -2,-2,4,0, 4,0,-2,-2/                   12d20s21
      data (((ibrt(j,i,k,1),j=1,4),i=1,7),k=1,4)                        12d20s21
     $     /0,0,1,-1, 0,0,-1,-1, 1,1,0,0, -1,1,0,0, 2,0,0,0, 1,1,1,1,   12d20s21
     $     1,1,-1,1, 0,0,1,-1, 0,0,-1,-1, -1,-1,0,0, -1,1,0,0, 0,-2,0,0,12d20s21
     $     -1,-1,-1,1, -1,-1,1,1, 0,0,1,-1, 0,0,-1,-1, -1,-1,0,0,       12d20s21
     $     -1,1,0,0, 0,-2,0,0, -1,-1,-1,1, -1,-1,1,1, 0,0,1,-1,         12d20s21
     $     0,0,-1,-1, 1,1,0,0, -1,1,0,0, 2,0,0,0, 1,1,1,1, 1,1,-1,1/    1d31s22
      data (((ibrt(j,i,k,2),j=1,4),i=1,10),k=1,4)                       12d20s21
     $     /0,0,1,1, 0,0,1,-1, 1,1,0,0, 1,-1,0,0, 2,0,0,0, 0,0,2,0,     12d20s21
     $     1,-1,-1,1, -1,1,-1,-1, 1,1,-1,-1, -1,-1,-1,1, 0,0,1,1,       12d20s21
     $     0,0,1,-1, 1,1,0,0, 1,-1,0,0, 0,2,0,0, 0,0,2,0, 1,-1,-1,1,    12d20s21
     $     -1,1,-1,-1, -1,-1,-1,-1, 1,1,-1,1, 0,0,1,1, 0,0,1,-1,        12d20s21
     $     1,1,0,0, 1,-1,0,0, 0,2,0,0, 0,0,2,0, -1,1,-1,1, 1,-1,-1,-1,  12d20s21
     $     1,1,-1,-1, -1,-1,-1,1, 0,0,1,1, 0,0,1,-1, 1,1,0,0,           12d20s21
     $     1,-1,0,0, 2,0,0,0, 0,0,2,0, -1,1,-1,1, 1,-1,-1,-1,           12d20s21
     $     -1,-1,-1,-1, 1,1,-1,1/                                       12d20s21
      data (((ibrt(j,i,k,3),j=1,4),i=1,11),k=1,4)                       12d20s21
     $     /1,1,0,0, 1,-1,0,0, 0,0,1,1, 0,0,1,-1, 0,0,-2,0, -2,0,0,0,   12d20s21
     $     1,-1,1,1, 1,-1,-1,1, -1,-1,1,1, 1,1,1,-1, -1,1,-1,-1,        12d20s21
     $     0,0,1,1, 0,0,1,-1, 1,1,0,0, 1,-1,0,0, -2,0,0,0, 0,0,0,2,     12d20s21
     $     1,-1,-1,-1, 1,-1,-1,1, -1,-1,-1,-1, 1,1,1,-1, -1,1,1,1,      12d20s21
     $     0,0,1,1, 0,0,1,-1, 1,1,0,0, 1,-1,0,0, 0,-2,0,0, 0,0,-2,0,    12d20s21
     $     -1,1,1,1, -1,1,-1,1, -1,-1,1,1, 1,1,1,-1, 1,-1,-1,-1,        12d20s21
     $     0,0,1,1, 0,0,1,-1, 1,1,0,0, 1,-1,0,0, 0,0,0,2, 0,-2,0,0,     12d20s21
     $     -1,1,-1,-1, -1,1,-1,1, -1,-1,-1,-1, 1,1,1,-1, 1,-1,1,1/      12d20s21
      data (((ibrt(j,i,k,4),j=1,4),i=1,11),k=1,4)                       12d20s21
     $     /0,0,1,1, 0,0,1,-1, 1,1,0,0, 1,-1,0,0, 0,0,-2,0, -2,0,0,0,   12d20s21
     $     1,-1,1,-1, -1,1,1,1, 1,1,-1,1, 1,1,1,1, -1,-1,-1,-1,         12d20s21
     $     0,0,1,1, 0,0,1,-1, 1,1,0,0, 1,-1,0,0, -2,0,0,0, 0,0,0,2,     12d20s21
     $     1,-1,1,-1, -1,1,-1,-1, 1,1,-1,1, 1,1,-1,-1, -1,-1,1,1,       12d20s21
     $     0,0,1,1, 0,0,1,-1, 1,1,0,0, 1,-1,0,0, 0,2,0,0, 0,0,-2,0,     12d20s21
     $     1,-1,1,-1, -1,1,1,1, -1,-1,-1,1, -1,-1,1,1, 1,1,-1,-1,       12d20s21
     $     0,0,1,1, 0,0,1,-1, 1,1,0,0, 1,-1,0,0, 0,2,0,0, 0,0,0,2,      12d20s21
     $     1,-1,1,-1, -1,1,-1,-1, -1,-1,-1,1, -1,-1,-1,-1, 1,1,1,1/     12d20s21
c
c itu=1
c nn'mm'   nmn'm'      nm'n'm
c +--+ 4   +--+  4     ++--  2
c -+-+ 3   --++  2     -++-  4
c --++ 2   -+-+  3     -+-+  3
c ++++ 1   ++++  1     ++++  1
c itu=2
c nn'mm'   nmn'm'      nm'n'm
c ++++ 1   ++++  1     ++++  1
c --++ 2   -+-+  3     -+-+  3
c +--+ 4   +--+  4     ++--  -2
c -+-+ 3   --++  2     -++-  -4
c itu=3
c nn'mm'   nmn'm'      nm'n'm
c ++-+ 3   +-++  2     +++-  -4
c ---+ 4   ---+  4     -+--  -2
c +-++ 2   ++-+  3     ++-+  3
c -+++ 1   -+++  1     -+++  1
c itu=4
c nn'mm'   nmn'm'      nm'n'm
c ++-+ 3   +-++  2     +++- 4
c ---+ 4   ---+  4     -+-- 2
c +-++ 2   ++-+  3     ++-+ 3
c -+++ 1   -+++  1     -+++ 1
c      k"k"'      k"kk"'k'     k"k'k"'k       kk'
c     G           (ac|bd)      (ad|bc)       V
c       11        1212gb4      1212gb4       22
c       21        2211c1       2112gb3       21
c                 2112gb3      2211c1        12
c       12        1221gb2      1122c2        21
c                 1122c2       1221gb2       12
c       22        2121gb1      2121gb1       11
c
      icall=icall+1
      loop=0
      idwsdeb=0
      igoal=1+5*2
      igoal1=2040128
      igoal2=2040105
      igoal3=2036113
      igoal4=1
      igoal5=1912564
      igoal6=1933267
      igoal7=2
      igoal8=1
      ituc=itu                                                          1d31s22
      if(idorel.eq.2)then                                               12d20s21
       nerimat=2                                                        12d20s21
       ituc=itu+1                                                       1d31s22
       mycode=10+ituc                                                   1d31s22
      else if(idorel.eq.-2)then
       if(itu.eq.1)then                                                 12d20s21
        nerimat=40                                                      12d20s21
       else if(itu.eq.2)then                                            12d20s21
        nerimat=50                                                      12d20s21
       else if(itu.ge.3)then                                            12d20s21
        nerimat=54                                                      12d20s21
       end if                                                           12d20s21
       mycode=18+itu                                                    12d20s21
      else if(idorel.eq.-4)then                                              12d20s21
       if(itu.eq.1)then                                                 12d20s21
        nerimat=12                                                      12d20s21
       else                                                             12d20s21
        nerimat=10                                                      12d20s21
       end if                                                           12d20s21
       mycode=14+itu                                                    12d20s21
      else                                                              12d20s21
       write(6,*)('unknown idorel code: '),idorel
       stop 'paraerikbk4'                                               12d20s21
      end if
      clight=0.5d0/sqrt(ascale)                                         2d21s20
      dgemmf=1d0/(8d0*clight*clight)                                    12d20s21
      do isb=1,nsymb                                                    12d23s21
       if(min(nvirt(isb),nbasisp(isb)).gt.0)then                        12d23s21
        jorb=iorb(isb)+nbasispc(isb)*(irefo(isb)+idoubo(isb))           12d23s21
        iorbt(isb)=ibcoff                                               12d23s21
        ibcoff=iorbt(isb)+nbasispc(isb)*nvirt(isb)                      12d23s21
        call enough('paraerikbk4.  1',bc,ibc)
        do iv=0,nvirt(isb)-1                                            12d23s21
         do i=0,nbasispc(isb)-1                                         12d23s21
          jorb1=jorb+i+nbasispc(isb)*iv                                 12d23s21
          jorbt1=iorbt(isb)+iv+nvirt(isb)*i                             12d23s21
          bc(jorbt1)=bc(jorb1)                                          12d23s21
         end do                                                         12d23s21
        end do                                                          12d23s21
       end if                                                           12d23s21
      end do                                                            12d23s21
c
c     build shell pair order
c
      nn=ngaus*ngaus                                                    12d20s21
      nn2=nn*2                                                          2d19s10
      ipair=ibcoff                                                      2d19s10
      ibcoff=ipair+nn*2                                                 3d4s10
      i12=ibcoff                                                        2d19s10
      ibcoff=i12+nn*4                                                   2d19s10
      call enough('paraerikbk4.  2',bc,ibc)
      j12=i12                                                           2d19s10
      jbdat=ibdat-1                                                     2d19s10
      ngaus3=ngaus*3                                                    5d12s10
      ngaus7=ngaus*7                                                    5d3s10
c
c     what's under ibdat
c     for each gaussian...
c     1: l
c     2: exponential parameter
c     3: normalization coefficient
c     4: offset in total basis function list
c     5: x of center
c     6: y of center
c     7: z of center
c     8: nhere -> address in iapair of possible symmetry equivalent
c        center
c
c     what's under iapair:
c     1: if nonzero, iabs is symmetry equivalent basis fcn
c     2: group operation relating symmetry equivalent basis fcns
c     3: phase?
c
c
c     estimate work for shell pair
c
      do i1=1,ngaus                                                     2d19s10
       do i2=1,ngaus                                                    12d20s21
        ibc(j12)=-((ibc(jbdat+i1)+ibc(jbdat+i2)+2)**2)                  2d19s10
        if(iapair(1,ibc(jbdat+i1+ngaus7)).gt.0)ibc(j12)=ibc(j12)*2      5d3s10
        if(iapair(1,ibc(jbdat+i2+ngaus7)).gt.0)ibc(j12)=ibc(j12)*2      5d3s10
        ibc(j12+nn)=i1                                                  2d19s10
        ibc(j12+nn2)=i2                                                 2d19s10
        j12=j12+1                                                       2d19s10
       end do                                                           2d19s10
      end do                                                            2d19s10
      is12=i12+nn*3                                                     2d19s10
      nn8=nn                                                            2d19s10
c
c     sort, so we can assign most expensive blocks first
c
      call idsortdws(ibc(i12),ibc(is12),nn8)                            1d18s23
      ii1=i12+nn                                                        2d19s10
      ii2=i12+nn2                                                       2d19s10
      jpair=ipair                                                       2d19s10
c
c     order shell pair indices
c
      do i=1,nn                                                         2d19s10
       j=ibc(is12+i-1)-1                                                2d19s10
       ibc(jpair)=ibc(ii1+j)                                            2d19s10
       jpair=jpair+1                                                    2d19s10
       ibc(jpair)=ibc(ii2+j)                                            2d19s10
       jpair=jpair+1                                                    2d19s10
      end do                                                            2d19s10
      needt=0                                                           12d12s18
      ibptoh=ibcoff                                                     12d12s18
      do isb=1,nsymb                                                    12d12s18
       isbv12=multh(isb,isymw)                                          12d23s21
       nfcnb(isb)=nrootu*(                                              1d26s22
     $      nfdat(2,1,isb)+nfdat(2,2,isb)+nfdat(2,3,isb)+nfdat(2,4,isb))1d26s22
       do isbv1=1,nsymb                                                 12d12s18
        isbv2=multh(isbv1,isbv12)                                       12d12s18
        if(isbv2.ge.isbv1)then                                          12d12s18
         iptoh(isbv1,isbv2)=ibcoff                                      12d12s18
         need=nvirt(isbv1)*nbasisp(isbv2)                               11d2s20
         iiii=((isbv2*(isbv2-1))/2)+isbv1
         need=need*nfcnb(isb)                                           1d25s22
         ibcoff=ibcoff+need*2                                           12d23s21
         needt=needt+need*2                                             12d23s21
        end if                                                          12d12s18
       end do                                                           12d12s18
      end do                                                            12d12s18
      call enough('paraerikbk4.  3',bc,ibc)
      do i=0,needt-1                                                    12d12s18
       bc(ibptoh+i)=0d0                                                 12d12s18
      end do                                                            12d12s18
      idoit=10
      do i=1+mynowprog,nn,mynprocg                                      3d4s10
       jpair=ipair+2*(i-1)                                              2d19s10
    2  format('I am going to do ',i5,3i3,f10.4)                               2d19s10
       idoit=idoit+1
       ja=jbdat+ibc(jpair)                                              2d19s10
       ja2=ja+ngaus                                                     2d19s10
       ja3=ja2+ngaus                                                    2d19s10
       ja4=ja3+ngaus                                                    2d19s10
       ja5=ja4+ngaus                                                    2d19s10
       ja6=ja5+ngaus                                                    2d19s10
       ja7=ja6+ngaus                                                    2d19s10
       ja8=ja7+ngaus                                                    5d11s10
       nsza=2*ibc(ja)+1                                                 5d11s10
       nsza0=nsza                                                       11d26s18
       nszaa=nsza                                                       5d12s10
       nszaa0=nsza0                                                     10d31s20
       if(iapair(1,ibc(ja8)).gt.0)then                                  5d11s10
        nszaa=nszaa*2                                                   5d12s10
        nszaa0=nszaa0*2                                                 10d31s20
       end if                                                           5d11s10
       jb=jbdat+ibc(jpair+1)                                            2d19s10
       jb2=jb+ngaus                                                     2d19s10
       jb3=jb2+ngaus                                                    2d19s10
       jb4=jb3+ngaus                                                    2d19s10
       jb5=jb4+ngaus                                                    2d19s10
       jb6=jb5+ngaus                                                    2d19s10
       jb7=jb6+ngaus                                                    2d19s10
       jb8=jb7+ngaus                                                    5d11s10
       nszb=2*ibc(jb)+1                                                 5d11s10
       nszb0=nszb                                                       11d26s18
       nszba=nszb                                                       5d12s10
       nszba0=nszb0                                                     11d2s20
 3352  format(6i5)
       if(iapair(1,ibc(jb8)).gt.0)then                                  5d11s10
        nszba=nszba*2                                                   5d12s10
        nszba0=nszba0*2                                                 10d31s20
       end if                                                           5d11s10
       ncsz=nsza*nszb                                                   5d11s10
       ncsz0=nsza0*nszb0                                                10d31s20
       ncsza=nszaa*nszba                                                1d10s11
       ncsza0=nszaa0*nszba0                                             10d31s20
c
c     index of of itmpab is symmetry of contracted n-2 electron fcns
c
       do isd=1,nsymb                                                   12d7s18
        itmpab(isd)=ibcoff                                              12d7s18
        nhere=nfcnb(isd)                                                1d25s22
        ibcoff=itmpab(isd)+nhere*ncsza0*4                               12d23s21
        ntmpab(isd)=nhere                                               12d7s18
        call enough('paraerikbk4.  4',bc,ibc)
        do iz=itmpab(isd),ibcoff-1                                      10d31s20
         bc(iz)=0d0                                                     10d31s20
        end do                                                          12d7s18
       end do                                                           5d12s10
       do isd=1,nsymb                                                   5d12s10
        do isc=1,nsymb                                                  5d12s10
         ieraw(isc,isd)=ibcoff                                          5d12s10
         ibcoff=ieraw(isc,isd)+ncsz0*nbasisp(isc)*nbasisp(isd)*nerimat  10d31s20
         call enough('paraerikbk4.  5',bc,ibc)
         do iz=ieraw(isc,isd),ibcoff-1                                  10d31s20
          bc(iz)=xnan                                                   10d31s20
         end do
        end do                                                          5d12s10
       end do                                                           5d12s10
       memneed(1:5)='tmpab'
       call enough('paraerikbk4.  6',bc,ibc)
       do ipassb=1,nszba/nszb                                           5d12s10
        if(ipassb.eq.1)then
         cartb(1)=bc(jb5)                                               5d12s10
         cartb(2)=bc(jb6)                                               5d12s10
         cartb(3)=bc(jb7)                                               5d12s10
         iob=0                                                          5d12s10
        else                                                            5d12s10
         cartb(1)=bc(jb5)*dfloat(isym(1,iapair(2,ibc(jb8))))            5d11s10
         cartb(2)=bc(jb6)*dfloat(isym(2,iapair(2,ibc(jb8))))            5d11s10
         cartb(3)=bc(jb7)*dfloat(isym(3,iapair(2,ibc(jb8))))            5d11s10
         iob=nszb0                                                      11d2s20
        end if                                                          5d12s10
        do ipassa=1,nszaa/nsza                                          5d12s10
         if(ipassa.eq.1)then                                            5d12s10
          carta(1)=bc(ja5)                                              5d12s10
          carta(2)=bc(ja6)                                              8d15s11
          carta(3)=bc(ja7)                                              8d15s11
          ioa=0                                                         5d12s10
         else                                                           5d12s10
          carta(1)=bc(ja5)*dfloat(isym(1,iapair(2,ibc(ja8))))           5d12s10
          carta(2)=bc(ja6)*dfloat(isym(2,iapair(2,ibc(ja8))))           5d12s10
          carta(3)=bc(ja7)*dfloat(isym(3,iapair(2,ibc(ja8))))           5d12s10
          ioa=nsza0                                                     11d2s20
         end if                                                         5d12s10
         if(idoit.eq.1)then
          write(6,*)('ipassa,ipassb: '),ipassa,ipassb
         end if
         do jjc=1,ngaus                                                 11d28s18
          jc=jbdat+jjc
          jc2=jc+ngaus                                                    2d19s10
          jc3=jc2+ngaus                                                   2d19s10
          jc4=jc3+ngaus                                                   2d19s10
          jc5=jc4+ngaus                                                   2d19s10
          jc6=jc5+ngaus                                                   2d19s10
          jc7=jc6+ngaus                                                   2d19s10
          jc8=jc7+ngaus                                                   5d11s10
          do jjd=1,ngaus                                                11d28s18
           jd=jbdat+jjd                                                 11d23s18
           jd2=jd+ngaus                                                    2d19s10
           jd3=jd2+ngaus                                                   2d19s10
           jd4=jd3+ngaus                                                   2d19s10
           jd5=jd4+ngaus                                                   2d19s10
           jd6=jd5+ngaus                                                   2d19s10
           jd7=jd6+ngaus                                                   2d19s10
           jd8=jd7+ngaus                                                   5d11s10
           nszc=2*ibc(jc)+1                                              2d7s12
           nszc0=nszc                                                   10d28s20
           nszca=nszc                                                    2d7s12
           nszca0=nszc0                                                 10d31s20
           nszd=2*ibc(jd)+1                                              2d7s12
           nszd0=nszd                                                   10d28s20
           nszda=nszd                                                    2d7s12
           nszda0=nszd0                                                 10d31s20
           if(iapair(1,ibc(jc8)).gt.0)then                                5d12s10
            nszca=nszca*2                                                 5d12s10
            nszca0=nszca0*2                                             10d31s20
           end if                                                         5d12s10
           if(iapair(1,ibc(jd8)).gt.0)then                                5d12s10
            nszda=nszda*2                                                 5d12s10
            nszda0=nszda0*2                                             10d31s20
           end if                                                         5d12s10
           itmp=ibcoff                                                    5d12s10
           jtmp=itmp-1                                                  11d26s18
           ibcoff=itmp+ncsz0*nszca0*nszda0*nerimat                      10d31s20
           memneed(1:4)='itmp'
           call enough('paraerikbk4.  7',bc,ibc)
           do iz=itmp,ibcoff-1                                          10d31s20
            bc(iz)=0d0                                                  10d31s20
           end do                                                       12d10s18
           jdoit=10
c     (ac|bd)=(ab||cd), inner loop over shells is c,d
           call dgerid(
     $            ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),      9d24s21
     $            ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),      12d20s21
     $            ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),      9d24s21
     $            ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),      12d20s21
     $            0,0,0, 0,0,0, 0,0,0, 0,0,0,ieri,.false.,mycode,       9d24s21
     $            fmulx,bc,ibc)                                         11d14s22
           if(idoit.eq.1.or.idwsdeb.gt.10)then
            write(6,*)('output from eri: '),ieri                        10d31s20
            call prntm2(bc(ieri),nszd0*nszc0,ncsz0*nerimat,nszd0*nszc0) 10d31s20
           end if
           do in=0,nerimat-1                                            10d31s20
            do id=0,nszd0-1                                               12d10s18
             do ic=0,nszc0-1                                              12d10s18
              do ib=0,nszb0-1                                             12d10s18
               do ia=0,nsza0-1                                            12d10s18
                iad1=ieri+id+nszd0*(ib+nszb0*(ic+nszc0*(ia+nsza0*in)))  10d31s20
                iad2=itmp+ia+nsza0*(ib+nszb0*(ic+nszca0*(id+nszda0*in))) 10d31s20
                bc(iad2)=bc(iad1)                                        12d10s18
               end do                                                    12d10s18
              end do                                                     12d10s18
             end do                                                      12d10s18
            end do                                                       12d10s18
           end do                                                       11d2s20
           if(idoit.eq.1.or.idwsdeb.gt.10)then
            write(6,*)('transposed ')
            call prntm2(bc(itmp),ncsz0,nszca0*nszd0*nerimat,ncsz0)      10d31s20
           end if
           if(iapair(1,ibc(jc8)).gt.0)then                                5d12s10
            cartc(1)=bc(jc5)*dfloat(isym(1,iapair(2,ibc(jc8))))           5d12s10
            cartc(2)=bc(jc6)*dfloat(isym(2,iapair(2,ibc(jc8))))           5d12s10
            cartc(3)=bc(jc7)*dfloat(isym(3,iapair(2,ibc(jc8))))           5d12s10
            call dgerid(ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),  4d23s18
     $         ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3),         4d23s18
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       4d23s18
     $         ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),         4d23s18
     $            0,0,0, 0,0,0, 0,0,0, 0,0,0,ieri,.false.,mycode,       9d24s21
     $            fmulx,bc,ibc)                                         11d14s22
            if(idoit.eq.1.or.idwsdeb.gt.10)then
             write(6,*)('outputb from eri: '),ieri
             call prntm2(bc(ieri),nszd0*nszc0,ncsz0*nerimat,nszd0*nszc0)
            end if
            do in=0,nerimat-1                                           10d31s20
             do id=0,nszd0-1                                               12d10s18
              do ic=0,nszc0-1                                              12d10s18
               do ib=0,nszb0-1                                             12d10s18
                do ia=0,nsza0-1                                            12d10s18
                 iad1=ieri+id+nszd0*(ib+nszb0*(ic+nszc0*(ia+nsza0*in))) 10d31s20
                 iad2=itmp+ia+nsza0*(ib+nszb0*(ic+nszc0+nszca0*(id      10d31s20
     $                +nszda0*in)))                                     10d31s20
                 bc(iad2)=bc(iad1)                                        12d10s18
                end do                                                    12d10s18
               end do                                                     12d10s18
              end do                                                      12d10s18
             end do                                                       12d10s18
            end do                                                      10d31s20
            if(idoit.eq.1.or.idwsdeb.gt.10)then
             write(6,*)('transposed ')
             call prntm2(bc(itmp),ncsz0,nszca0*nszd0*nerimat,ncsz0)      10d31s20
            end if
           end if                                                         5d12s10
           if(iapair(1,ibc(jd8)).gt.0)then                              1d10s11
            cartd(1)=bc(jd5)*dfloat(isym(1,iapair(2,ibc(jd8))))           5d12s10
            cartd(2)=bc(jd6)*dfloat(isym(2,iapair(2,ibc(jd8))))           5d12s10
            cartd(3)=bc(jd7)*dfloat(isym(3,iapair(2,ibc(jd8))))           5d12s10
            call dgerid(ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),  4d23s18
     $         ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),         4d23s18
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       4d23s18
     $         ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3),         4d23s18
     $            0,0,0, 0,0,0, 0,0,0, 0,0,0,ieri,.false.,mycode,       9d24s21
     $            fmulx,bc,ibc)                                         11d14s22
            if(idoit.eq.1.or.idwsdeb.gt.10)then
             write(6,*)('outputc from eri: '),ieri,nszd,nszc,ncsz
             call prntm2(bc(ieri),nszd0*nszc0,ncsz0*nerimat,nszd0*nszc0)10d31s20
            end if
            do in=0,nerimat-1                                           10d31s20
             do id=0,nszd0-1                                               12d10s18
              do ic=0,nszc0-1                                              12d10s18
               do ib=0,nszb0-1                                             12d10s18
                do ia=0,nsza0-1                                            12d10s18
                 iad1=ieri+id+nszd0*(ib+nszb0*(ic+nszc0*(ia+ncsza0*in)))10d31s20
                 iad2=itmp+ia+nsza0*(ib+nszb0*(ic+nszca0*(id+nszd0      10d31s20
     $                +nszda0*in)))                                     10d31s20
                 bc(iad2)=bc(iad1)                                        12d10s18
                end do                                                  10d31s20
               end do                                                    12d10s18
              end do                                                     12d10s18
             end do                                                      12d10s18
            end do                                                       12d10s18
            if(idoit.eq.1.or.idwsdeb.gt.10)then
             write(6,*)('transposed ')
             call prntm2(bc(itmp),ncsz0,nszca0*nszd0*nerimat,ncsz0)      10d31s20
            end if
            if(iapair(1,ibc(jc8)).gt.0)then                                5d12s10
             call dgerid(                                               12d20s21
     $           ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),       12d20s21
     $         ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3),         4d23s18
     $           ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),       4d23s18
     $         ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3),         4d23s18
     $            0,0,0, 0,0,0, 0,0,0, 0,0,0,ieri,.false.,mycode,       9d24s21
     $            fmulx,bc,ibc)                                         11d14s22
             if(idoit.eq.1.or.idwsdeb.gt.10)then
              write(6,*)('outputd from eri: '),ieri,nszd,nszc,ncsz
              call prntm2(bc(ieri),nszd0*nszc0,ncsz0*nerimat,
     $             nszd0*nszc0)
             end if
             do in=0,nerimat-1                                          10d31s20
              do id=0,nszd0-1                                               12d10s18
               do ic=0,nszc0-1                                              12d10s18
                do ib=0,nszb0-1                                             12d10s18
                 do ia=0,nsza0-1                                            12d10s18
                  iad1=ieri+id+nszd0*(ib+nszb0*(ic+nszc0*(ia+nsza0*in)))10d31s20
                  iad2=itmp+ia+nsza0*(ib+nszb0*(ic+nszc0+nszca*(id      10d31s20
     $                 +nszd0+nszda0*in)))                              10d31s20
                  bc(iad2)=bc(iad1)                                        12d10s18
                 end do                                                 11d2s20
                end do                                                    12d10s18
               end do                                                     12d10s18
              end do                                                      12d10s18
             end do                                                       12d10s18
             if(idoit.eq.1.or.idwsdeb.gt.10)then
              write(6,*)('transposed ')
              call prntm2(bc(itmp),ncsz0,nszca0*nszd0*nerimat,ncsz0)      10d31s20
             end if
            end if                                                         5d12s10
           end if                                                         5d12s10
           nszabc0=ncsz0*nszca0                                             5d12s10
           if(nszd.ne.nszda)then                                          5d12s10
            do in=0,nerimat-1                                           10d31s20
             jtmp=itmp+nszabc0*nszda0*in                                10d31s20
             if(idoit.eq.1.or.idwsdeb.gt.10)then
              write(6,*)('before sum/dif: matrix '),in
              call prntm2(bc(jtmp),nszabc0,nszda0,nszabc0)              10d31s20
             end if
             do id=0,nszd0-1                                            10d31s20
              it1=jtmp+nszabc0*id                                       10d31s20
              it2=jtmp+nszabc0*(id+nszd0)                               10d31s20
              do iabc=0,nszabc0-1                                       1d21s21
               sum=srh*(bc(it1+iabc)+bc(it2+iabc))                         5d12s10
               diff=srh*(-bc(it1+iabc)+bc(it2+iabc))                       5d12s10
               bc(it1+iabc)=sum                                            5d12s10
               bc(it2+iabc)=diff                                           5d12s10
              end do                                                       5d12s10
             end do                                                        5d12s10
             if(idoit.eq.1.or.idwsdeb.gt.10)then
              write(6,*)('after sum/dif: ')
              call prntm2(bc(jtmp),nszabc0,nszda0,nszabc0)
             end if
            end do                                                      10d31s20
           end if                                                         5d12s10
           if(nszc.ne.nszca)then                                          5d12s10
            do in=0,nerimat-1                                           10d31s20
             jtmp=itmp+nszc0*nszda0*in                                  10d31s20
             do id=0,nszda0-1                                           10d31s20
              do ic=0,nszc0-1                                           10d31s20
               it1=jtmp+ncsz0*(ic+nszca0*id)                            10d31s20
               it2=jtmp+ncsz0*(ic+nszc0+nszca0*id)                      10d31s20
               do iab=0,ncsz0-1                                         10d31s20
                sum=srh*(bc(it1+iab)+bc(it2+iab))                          5d12s10
                diff=srh*(-bc(it1+iab)+bc(it2+iab))                        5d12s10
                bc(it1+iab)=sum                                            5d12s10
                bc(it2+iab)=diff                                           5d12s10
               end do                                                      5d12s10
              end do                                                       5d12s10
             end do                                                        5d12s10
             if(idoit.eq.1.or.idwsdeb.gt.10)then
              write(6,*)('after sum/dif: for matrix '),in
              call prntm2(bc(jtmp),nszabc0,nszda0,nszabc0)
             end if
            end do                                                      10d31s20
           end if                                                         5d12s10
c
           do id=1,nszda0                                               8d26s15
            idm=id-1                                                    10d31s20
            idd=id+ibc(jd4)                                               5d12s10
            isd=isstor(idd)                                               5d12s10
            iddd=ibstor(idd)-1                                          10d31s20
            do ic=1,nszca0                                              8d26s15
             icm=ic-1                                                   10d31s20
             icc=ic+ibc(jc4)                                              5d12s10
             isc=isstor(icc)                                              5d12s10
             iccc=ibstor(icc)-1                                         10d31s20
             do in=0,nerimat-1                                          10d31s20
              iad1=ieraw(isc,isd)+ncsz0*(iccc+nbasisp(isc)*(iddd        10d31s20
     $             +nbasisp(isd)*in))                                   10d31s20
              iad2=itmp+ncsz0*(icm+nszca0*(idm+nszda0*in))               10d31s20
              do iab=0,ncsz0-1                                          11d2s20
               bc(iad1+iab)=bc(iad2+iab)                                   5d12s10
              end do                                                       5d12s10
             end do                                                     10d31s20
            end do                                                        5d12s10
           end do                                                         5d12s10
           ibcoff=itmp                                                    5d12s10
          end do                                                        11d26s18
c
c     end of jjc loop
c
         end do                                                           3d4s10
c
c     at this point, we have all of cd ints,
c     and furthermore they have been symmetrized (turned into so's) and 11d27s18
c     stored under eraw.                                                11d27s18
c     so do the b<d restriction and then multiply by ci vector in ao    11d27s18
c     basis.
c                                                                       1d10s11
         do isd=1,nsymb                                                 11d27s18
          do isc=1,isd                                                  11d27s18
           iiii=((isd*(isd-1))/2)+isc                                   11d27s18
c           k"k"' spin kk'      k"kk"'k'f
c     form E           =   sum c (ac|bd) - c (ad|bc)
c           abcd             f  spin
c
c
c     iaddr(iiii,in) with in=ispin,k,k'
c     2 coulomb, -4 gaunt, -2 breit
c
           nvv=nbasisp(isc)*nbasisp(isd)                                12d23s21
           if(min(nvv,nfcn(iiii,1)).gt.0)then                           12d23s21
            itmpg21=ibcoff                                               12d23s21
            itmpg12=itmpg21+ncsz0*nfcn(iiii,1)                           12d23s21
            ibcoff=itmpg12+ncsz0*nfcn(iiii,1)                            12d23s21
            if(idorel.lt.0)then                                          12d23s21
             itmpg11=ibcoff                                              12d23s21
             itmpg22=itmpg11+ncsz0*nfcn(iiii,1)                           12d23s21
             ibcoff=itmpg22+ncsz0*nfcn(iiii,1)                            12d23s21
            end if                                                       12d23s21
            call enough('paraerikbk4.  8',bc,ibc)
            do iz=itmpg21,ibcoff-1                                      12d30s21
             bc(iz)=0d0                                                  12d23s21
            end do                                                       12d23s21
            itmpe2211=ibcoff                                                 12d20s21
            nnx=ncsz0*nvv                                                12d23s21
            itmpe1122=itmpe2211+nnx                                      12d23s21
            itmpe2121=itmpe1122+nnx                                      12d23s21
            itmpe2112=itmpe2121+nnx                                        12d20s21
            itmpe1221=itmpe2112+nnx                                        12d20s21
            itmpe1212=itmpe1221+nnx                                        12d20s21
            ibcoff=itmpe1212+nnx                                          12d20s21
            call enough('paraerikbk4.  9',bc,ibc)
            nnxm=nnx-1                                                   12d20s21
            nnx6m=nnx*6-1                                                12d20s21
c
c      k"k"'      k"kk"'k'     k"k'k"'k       kk'
c     G           (ac|bd)  -   (ad|bc)       V
c       22        2121gb1      2121gb1       11
c       11        1212gb4      1212gb4       22
c       21        2211c1       2112gb3       21
c                 2112gb3      2211c1        12
c       12        1221gb2      1122c2        21
c                 1122c2       1221gb2       12
c
            do ispn=1,ntype                                             2d1s22
             jspn1=ispino(ispn,1,ituc)                                  3d6s22
             jspn2=ispino(ispn,2,ituc)                                  3d6s22
             do iz=0,nnx6m                                               12d20s21
              bc(itmpe2211+iz)=0d0                                       12d23s21
             end do                                                      12d20s21
             if0=0                                                       12d20s21
             if(itu.ne.1.or.idorel.eq.2)then                            1d31s22
              if0=2                                                      12d20s21
              if(icoul(jspn1,1,ituc).ne.0)then                          1d31s22
               ff=dfloat(icoul(jspn1,1,ituc))*dgemmf                    1d31s22
               jj=ieraw(isc,isd)                                         12d20s21
               do ii=0,nnxm                                               12d20s21
                bc(itmpe2211+ii)=bc(itmpe2211+ii)+ff*bc(jj+ii)          12d30s21
               end do                                                    12d20s21
              end if                                                     12d20s21
              if(icoul(jspn2,1,ituc).ne.0)then                          1d31s22
               ff=-dfloat(icoul(jspn2,1,ituc))*dgemmf                   1d31s22
               ff=ff*spino(ispn,ituc)                                   1d31s22
               jj=ieraw(isd,isc)                                         12d20s21
               do id=0,nbasisp(isd)-1                                   12d20s21
                do ic=0,nbasisp(isc)-1                                  12d20s21
                 idc=jj+ncsz0*(id+nbasisp(isd)*ic)                      12d20s21
                 icd=itmpe2112+ncsz0*(ic+nbasisp(isc)*id)                12d20s21
                 do ii=0,ncsz0-1                                         12d20s21
                  bc(icd+ii)=bc(icd+ii)+ff*bc(idc+ii)                   12d23s21
                 end do                                                 12d20s21
                end do                                                  12d20s21
               end do                                                   12d20s21
              end if                                                     12d20s21
              if(icoul(jspn1,2,ituc).ne.0)then                          1d31s22
               ff=dfloat(icoul(jspn1,2,ituc))*dgemmf                    1d31s22
               jj=ieraw(isc,isd)+nnx                                     12d20s21
               do ii=0,nnxm                                               12d20s21
                bc(itmpe1122+ii)=bc(itmpe1122+ii)+ff*bc(jj+ii)          12d30s21
               end do                                                    12d20s21
              end if                                                     12d20s21
              if(icoul(jspn2,2,ituc).ne.0)then                          1d31s22
               ff=-dfloat(icoul(jspn2,2,ituc))*dgemmf                   1d31s22
               ff=ff*spino(ispn,ituc)                                   1d31s22
               jj=ieraw(isd,isc)+nnx                                     12d20s21
               do id=0,nbasisp(isd)-1                                   12d20s21
                do ic=0,nbasisp(isc)-1                                  12d20s21
                 idc=jj+ncsz0*(id+nbasisp(isd)*ic)                      12d20s21
                 icd=itmpe1221+ncsz0*(ic+nbasisp(isc)*id)                12d20s21
                 do ii=0,ncsz0-1                                         12d20s21
                  bc(icd+ii)=bc(icd+ii)+ff*bc(idc+ii)                   12d23s21
                 end do                                                 12d20s21
                end do                                                  12d20s21
               end do                                                   12d20s21
              end if                                                     12d20s21
             end if                                                      12d20s21
             if(idorel.lt.0)then                                         12d20s21
              if0=if0-1                                                  12d20s21
              if(itu.eq.1)then                                           12d20s21
               nf=3                                                      12d20s21
              else                                                       12d20s21
               nf=2                                                      12d20s21
              end if                                                     12d20s21
              nf2=nf*2                                                   12d20s21
              nf3=nf*3                                                   12d20s21
c
c         k
c        1234
c     1: 2121
c     2: 1221
c     3: 2112
c     4: 1212
c
              do if=1,nf                                                 12d20s21
               if(igaunt(jspn1,if,1,itu).ne.0)then                       12d20s21
                ff=dfloat(igaunt(jspn1,if,1,itu))*dgemmf                 12d20s21
                in=if+if0                                                12d20s21
                jj=ieraw(isc,isd)+nnx*in                                12d30s21
                do ii=0,nnxm                                              12d20s21
                 bc(itmpe2121+ii)=bc(itmpe2121+ii)+ff*bc(jj+ii)         12d30s21
                end do                                                   12d20s21
               end if                                                    12d20s21
               if(igaunt(jspn1,if,2,itu).ne.0)then                       12d20s21
                ff=dfloat(igaunt(jspn1,if,2,itu))*dgemmf                 12d20s21
                in=if+if0+nf                                                12d20s21
                jj=ieraw(isc,isd)+nnx*in                                12d30s21
                do ii=0,nnxm                                            12d23s21
                 bc(itmpe1221+ii)=bc(itmpe1221+ii)+ff*bc(jj+ii)         12d23s21
                end do                                                   12d20s21
               end if                                                    12d20s21
               if(igaunt(jspn1,if,3,itu).ne.0)then                       12d20s21
                ff=dfloat(igaunt(jspn1,if,3,itu))*dgemmf                 12d20s21
                in=if+if0+nf2                                                12d20s21
                jj=ieraw(isc,isd)+nnx*in                                12d30s21
                do ii=0,nnxm                                            12d23s21
                 bc(itmpe2112+ii)=bc(itmpe2112+ii)+ff*bc(jj+ii)         12d23s21
                end do                                                   12d20s21
               end if                                                    12d20s21
               if(igaunt(jspn1,if,4,itu).ne.0)then                       12d20s21
                ff=dfloat(igaunt(jspn1,if,4,itu))*dgemmf                 12d20s21
                in=if+if0+nf3                                                12d20s21
                jj=ieraw(isc,isd)+nnx*in                                12d30s21
                do ii=0,nnxm                                              12d20s21
                 bc(itmpe1212+ii)=bc(itmpe1212+ii)+ff*bc(jj+ii)         12d30s21
                end do                                                   12d20s21
               end if                                                    12d20s21
               if(igaunt(jspn2,if,1,itu).ne.0)then                       12d20s21
                ff=-dfloat(igaunt(jspn2,if,1,itu))*dgemmf                12d20s21
                ff=ff*spino(ispn,itu)                                   1d25s22
                in=if+if0                                                12d20s21
                jj=ieraw(isd,isc)+nnx*in                                 12d20s21
                do id=0,nbasisp(isd)-1                                   12d20s21
                 do ic=0,nbasisp(isc)-1                                  12d20s21
                  idc=jj+ncsz0*(id+nbasisp(isd)*ic)                      12d20s21
                  icd=itmpe2121+ncsz0*(ic+nbasisp(isc)*id)              1d28s22
                  do ii=0,ncsz0-1                                         12d20s21
                   bc(icd+ii)=bc(icd+ii)+ff*bc(idc+ii)                  12d23s21
                  end do                                                 12d20s21
                 end do                                                  12d20s21
                end do                                                   12d20s21
               end if                                                    12d20s21
               if(igaunt(jspn2,if,2,itu).ne.0)then                       12d20s21
                ff=-dfloat(igaunt(jspn2,if,2,itu))*dgemmf                12d20s21
                ff=ff*spino(ispn,itu)                                   1d25s22
                in=if+if0+nf                                                12d20s21
                jj=ieraw(isd,isc)+nnx*in                                 12d20s21
                do id=0,nbasisp(isd)-1                                   12d20s21
                 do ic=0,nbasisp(isc)-1                                  12d20s21
                  idc=jj+ncsz0*(id+nbasisp(isd)*ic)                      12d20s21
                  icd=itmpe1122+ncsz0*(ic+nbasisp(isc)*id)              1d28s22
                  do ii=0,ncsz0-1                                       12d23s21
                   bc(icd+ii)=bc(icd+ii)+ff*bc(idc+ii)                  12d23s21
                  end do                                                 12d20s21
                 end do                                                  12d20s21
                end do                                                   12d20s21
               end if                                                    12d20s21
               if(igaunt(jspn2,if,3,itu).ne.0)then                       12d20s21
                ff=-dfloat(igaunt(jspn2,if,3,itu))*dgemmf                12d20s21
                ff=ff*spino(ispn,itu)                                   1d25s22
                in=if+if0+nf2                                            12d20s21
                jj=ieraw(isd,isc)+nnx*in                                 12d20s21
                do id=0,nbasisp(isd)-1                                   12d20s21
                 do ic=0,nbasisp(isc)-1                                  12d20s21
                  idc=jj+ncsz0*(id+nbasisp(isd)*ic)                      12d20s21
                  icd=itmpe2211+ncsz0*(ic+nbasisp(isc)*id)              1d28s22
                  do ii=0,ncsz0-1                                       12d23s21
                   bc(icd+ii)=bc(icd+ii)+ff*bc(idc+ii)                  12d23s21
                  end do                                                 12d20s21
                 end do                                                  12d20s21
                end do                                                   12d20s21
               end if                                                    12d20s21
               if(igaunt(jspn2,if,4,itu).ne.0)then                       12d20s21
                ff=-dfloat(igaunt(jspn2,if,4,itu))*dgemmf                12d20s21
                ff=ff*spino(ispn,itu)                                   1d25s22
                in=if+if0+nf3                                            12d20s21
                jj=ieraw(isd,isc)+nnx*in                                 12d20s21
                do id=0,nbasisp(isd)-1                                   12d20s21
                 do ic=0,nbasisp(isc)-1                                  12d20s21
                  idc=jj+ncsz0*(id+nbasisp(isd)*ic)                      12d20s21
                  icd=itmpe1212+ncsz0*(ic+nbasisp(isc)*id)              1d28s22
                  do ii=0,ncsz0-1                                       12d23s21
                   bc(icd+ii)=bc(icd+ii)+ff*bc(idc+ii)                  12d23s21
                  end do                                                 12d20s21
                 end do                                                  12d20s21
                end do                                                   12d20s21
               end if                                                    12d20s21
              end do                                                     12d20s21
              if0=if0+nf*4                                               12d20s21
              if(idorel.eq.-2)then                                       12d20s21
               if(itu.eq.1)then                                          12d20s21
                nf=7                                                     12d20s21
               else if(itu.eq.2)then                                     12d20s21
                nf=10
               else if(itu.ge.3)then                                     12d20s21
                nf=11                                                    12d20s21
               end if                                                    12d20s21
               nf2=nf*2                                                   12d20s21
               nf3=nf*3                                                   12d20s21
               do if=1,nf                                                 12d20s21
                if(ibrt(jspn1,if,1,itu).ne.0)then                        12d20s21
                 ff=dfloat(ibrt(jspn1,if,1,itu))*dgemmf                  12d20s21
                 in=if+if0                                                12d20s21
                 jj=ieraw(isc,isd)+nnx*in                               12d30s21
                 do ii=0,nnxm                                              12d20s21
                  bc(itmpe2121+ii)=bc(itmpe2121+ii)+ff*bc(jj+ii)        12d30s21
                 end do                                                   12d20s21
                end if                                                    12d20s21
                if(ibrt(jspn1,if,2,itu).ne.0)then                        12d20s21
                 ff=dfloat(ibrt(jspn1,if,2,itu))*dgemmf                  12d20s21
                 in=if+if0+nf                                                12d20s21
                 jj=ieraw(isc,isd)+nnx*in                               12d30s21
                 do ii=0,nnxm                                              12d20s21
                  bc(itmpe1221+ii)=bc(itmpe1221+ii)+ff*bc(jj+ii)        12d23s21
                 end do                                                   12d20s21
                end if                                                    12d20s21
                if(ibrt(jspn1,if,3,itu).ne.0)then                        12d20s21
                 ff=dfloat(ibrt(jspn1,if,3,itu))*dgemmf                  12d20s21
                 in=if+if0+nf2                                                12d20s21
                 jj=ieraw(isc,isd)+nnx*in                               12d30s21
                 do ii=0,nnxm                                           12d23s21
                  bc(itmpe2112+ii)=bc(itmpe2112+ii)+ff*bc(jj+ii)        12d23s21
                 end do                                                   12d20s21
                end if                                                    12d20s21
                if(ibrt(jspn1,if,4,itu).ne.0)then                        12d20s21
                 ff=dfloat(ibrt(jspn1,if,4,itu))*dgemmf                  12d20s21
                 in=if+if0+nf3                                                12d20s21
                 jj=ieraw(isc,isd)+nnx*in                               12d30s21
                 do ii=0,nnxm                                           12d23s21
                  bc(itmpe1212+ii)=bc(itmpe1212+ii)+ff*bc(jj+ii)        12d30s21
                 end do                                                   12d20s21
                end if                                                    12d20s21
                if(ibrt(jspn2,if,1,itu).ne.0)then                        12d20s21
                 ff=-dfloat(ibrt(jspn2,if,1,itu))*dgemmf                 12d20s21
                 ff=ff*spino(ispn,itu)                                  1d25s22
                 in=if+if0                                                12d20s21
                 jj=ieraw(isd,isc)+nnx*in                               12d30s21
                 do id=0,nbasisp(isd)-1                                   12d20s21
                  do ic=0,nbasisp(isc)-1                                  12d20s21
                   idc=jj+ncsz0*(id+nbasisp(isd)*ic)                      12d20s21
                   icd=itmpe2121+ncsz0*(ic+nbasisp(isc)*id)             1d28s22
                   do ii=0,ncsz0-1                                      12d23s21
                    bc(icd+ii)=bc(icd+ii)+ff*bc(idc+ii)                 12d23s21
                   end do                                                 12d20s21
                  end do                                                  12d20s21
                 end do                                                   12d20s21
                end if                                                    12d20s21
                if(ibrt(jspn2,if,2,itu).ne.0)then                        12d20s21
                 ff=-dfloat(ibrt(jspn2,if,2,itu))*dgemmf                 12d20s21
                 ff=ff*spino(ispn,itu)                                  1d25s22
                 in=if+if0+nf                                                12d20s21
                 jj=ieraw(isd,isc)+nnx*in                               12d30s21
                 do id=0,nbasisp(isd)-1                                   12d20s21
                  do ic=0,nbasisp(isc)-1                                  12d20s21
                   idc=jj+ncsz0*(id+nbasisp(isd)*ic)                      12d20s21
                   icd=itmpe1122+ncsz0*(ic+nbasisp(isc)*id)             1d28s22
                   do ii=0,ncsz0-1                                      12d23s21
                    bc(icd+ii)=bc(icd+ii)+ff*bc(idc+ii)                 12d23s21
                   end do                                                 12d20s21
                  end do                                                  12d20s21
                 end do                                                   12d20s21
                end if                                                    12d20s21
                if(ibrt(jspn2,if,3,itu).ne.0)then                        12d20s21
                 ff=-dfloat(ibrt(jspn2,if,3,itu))*dgemmf                 12d20s21
                 ff=ff*spino(ispn,itu)                                  1d25s22
                 in=if+if0+nf2                                            12d20s21
                 jj=ieraw(isd,isc)+nnx*in                               12d30s21
                 do id=0,nbasisp(isd)-1                                   12d20s21
                  do ic=0,nbasisp(isc)-1                                  12d20s21
                   idc=jj+ncsz0*(id+nbasisp(isd)*ic)                      12d20s21
                   icd=itmpe2211+ncsz0*(ic+nbasisp(isc)*id)             1d28s22
                   do ii=0,ncsz0-1                                         12d20s21
                    bc(icd+ii)=bc(icd+ii)+ff*bc(idc+ii)                 12d23s21
                   end do                                                 12d20s21
                  end do                                                  12d20s21
                 end do                                                   12d20s21
                end if                                                    12d20s21
                if(ibrt(jspn2,if,4,itu).ne.0)then                        12d20s21
                 ff=-dfloat(ibrt(jspn2,if,4,itu))*dgemmf                 12d20s21
                 ff=ff*spino(ispn,itu)                                  1d25s22
                 in=if+if0+nf3                                            12d20s21
                 jj=ieraw(isd,isc)+nnx*in                               12d30s21
                 do id=0,nbasisp(isd)-1                                   12d20s21
                  do ic=0,nbasisp(isc)-1                                  12d20s21
                   idc=jj+ncsz0*(id+nbasisp(isd)*ic)                      12d20s21
                   icd=itmpe1212+ncsz0*(ic+nbasisp(isc)*id)             1d28s22
                   do ii=0,ncsz0-1                                         12d20s21
                    bc(icd+ii)=bc(icd+ii)+ff*bc(idc+ii)                 12d23s21
                   end do                                                 12d20s21
                  end do                                                  12d20s21
                 end do                                                   12d20s21
                end if                                                    12d20s21
               end do                                                     12d20s21
              end if                                                     12d20s21
             end if                                                      12d20s21
c
c      k"k"'      k"kk"'k'     k"k'k"'k       kk'
c     G           (ac|bd)  -   (ad|bc)       V
c       22        2121gb1      2121gb1       11
c       11        1212gb4      1212gb4       22
c       21        2211c1       2112gb3       21
c                 2112gb3      2211c1        12
c       12        1221gb2      1122c2        21
c                 1122c2       1221gb2       12
c
c
             inv11=ispn                                                  12d23s21
             inv21=ispn+ntype                                            12d23s21
             inv12=ispn+ntype*2                                          12d23s21
             inv22=ispn+ntype*3                                          12d23s21
             if(idorel.eq.2)then                                        1d31s22
              inv21=ispn                                                1d31s22
              inv12=ispn+ntype                                          1d31s22
              inv11=inv21                                               1d31s22
              inv22=inv21                                               1d31s22
             end if                                                     1d31s22
             call dgemm('n','n',ncsz0,nfcn(iiii,1),nvv,1d0,             12d23s21
     $            bc(itmpe2211),ncsz0,bc(iaddr(iiii,inv21)),nvv,1d0,    12d23s21
     $            bc(itmpg21),ncsz0,                                    12d23s21
     d' paraerikbk4.  1')
             call dgemm('n','n',ncsz0,nfcn(iiii,1),nvv,1d0,             12d23s21
     $            bc(itmpe2112),ncsz0,bc(iaddr(iiii,inv12)),nvv,1d0,    12d23s21
     $            bc(itmpg21),ncsz0,                                    12d23s21
     d' paraerikbk4.  2')
             call dgemm('n','n',ncsz0,nfcn(iiii,1),nvv,1d0,             12d23s21
     $            bc(itmpe1221),ncsz0,bc(iaddr(iiii,inv21)),nvv,1d0,    12d23s21
     $            bc(itmpg12),ncsz0,                                    12d23s21
     d' paraerikbk4.  3')
             call dgemm('n','n',ncsz0,nfcn(iiii,1),nvv,1d0,             12d23s21
     $            bc(itmpe1122),ncsz0,bc(iaddr(iiii,inv12)),nvv,1d0,    12d23s21
     $            bc(itmpg12),ncsz0,                                    12d23s21
     d' paraerikbk4.  4')
             if(idorel.lt.0)then                                         12d23s21
              sz=0d0
              do ixi=0,nvv*nfcn(iiii,1)-1
               sz=sz+bc(iaddr(iiii,inv11)+ixi)**2
              end do
              sz=sqrt(sz/dfloat(nvv*nfcn(iiii,1)))
              call dgemm('n','n',ncsz0,nfcn(iiii,1),nvv,1d0,             12d23s21
     $           bc(itmpe2121),ncsz0,bc(iaddr(iiii,inv11)),nvv,1d0,     12d23s21
     $           bc(itmpg22),ncsz0,                                     12d23s21
     d' paraerikbk4.  5')
              sz=0d0
              do ixi=0,nvv*nfcn(iiii,1)-1
               sz=sz+bc(iaddr(iiii,inv22)+ixi)**2
              end do
              sz=sqrt(sz/dfloat(nvv*nfcn(iiii,1)))
              call dgemm('n','n',ncsz0,nfcn(iiii,1),nvv,1d0,             12d23s21
     $           bc(itmpe1212),ncsz0,bc(iaddr(iiii,inv22)),nvv,1d0,     12d23s21
     $           bc(itmpg11),ncsz0,                                     12d23s21
     d' paraerikbk4.  6')
             end if                                                      12d23s21
            end do                                                       12d23s21
            nnn=nfcn(iiii,1)*nszaa0*nszba0                              12d23s21
            nnn3=nnn*3
            sz=0d0
            do ii=0,nfcn(iiii,1)-1                                      12d23s21
             do ib=0,nszb0-1                                            12d23s21
              ibp=ib+iob                                                12d23s21
              do ia=0,nsza0-1                                           12d23s21
               iap=ia+ioa                                               12d23s21
               iad21=itmpab(nfcn(iiii,2))+ii+nfcn(iiii,1)               12d23s21
     $              *(iap+nszaa0*(ibp+nszba0))                          12d23s21
               iad12=iad21+nnn                                          12d23s21
               ig21=itmpg21+ia+nsza0*(ib+nszb0*ii)                      12d23s21
               ig12=itmpg12+ia+nsza0*(ib+nszb0*ii)                      12d23s21
               bc(iad21)=bc(iad21)+bc(ig21)                             12d23s21
               bc(iad12)=bc(iad12)+bc(ig12)                             12d23s21
              end do                                                    12d23s21
             end do                                                      11d27s18
            end do                                                       11d27s18
            if(idorel.lt.0)then                                         12d23s21
             do ii=0,nfcn(iiii,1)-1                                      12d23s21
              do ib=0,nszb0-1                                            12d23s21
               ibp=ib+iob                                                12d23s21
               do ia=0,nsza0-1                                           12d23s21
                iap=ia+ioa                                               12d23s21
                iad11=itmpab(nfcn(iiii,2))+ii+nfcn(iiii,1)              12d23s21
     $              *(iap+nszaa0*ibp)                                   12d23s21
                iad22=iad11+nnn3                                        12d23s21
                ig11=itmpg11+ia+nsza0*(ib+nszb0*ii)                      12d23s21
                ig22=itmpg22+ia+nsza0*(ib+nszb0*ii)                      12d23s21
                bc(iad11)=bc(iad11)+bc(ig11)                             12d23s21
                bc(iad22)=bc(iad22)+bc(ig22)                             12d23s21
               end do                                                    12d23s21
              end do                                                      11d27s18
             end do                                                       11d27s18
            end if                                                      12d23s21
            ibcoff=itmpg21                                              12d23s21
           end if                                                       11d28s18
          end do                                                        11d27s18
         end do                                                         11d27s18
c     end of ipassa loop
        end do                                                          5d12s10
c     end of ipassb loop
       end do                                                           5d12s10
c
c     complete ao to so transformation for ab                           1d10s11
c
       if(iapair(1,ibc(jb8)).gt.0)then                                  5d12s10
        do isd=1,nsymb                                                  5d12s10
         mcd=nfcnb(isd)                                                 1d25s22
         do in=0,3                                                      12d23s21
          do ib=0,nszb0-1                                               11d2s20
           do ia=0,nszaa0-1                                             11d2s20
            it1=itmpab(isd)+mcd*(ia+nszaa0*(ib+nszba0*in))              11d2s20
            it2=itmpab(isd)+mcd*(ia+nszaa0*(ib+nszb0+nszba0*in))        11d2s20
            do icd=0,mcd-1                                              11d2s20
             sum=srh*(bc(it1+icd)+bc(it2+icd))                           12d7s18
             diff=srh*(-bc(it1+icd)+bc(it2+icd))                         12d7s18
             bc(it1+icd)=sum                                             12d7s18
             bc(it2+icd)=diff                                            12d7s18
            end do                                                       12d7s18
           end do                                                        12d7s18
          end do                                                         12d7s18
         end do                                                          5d12s10
        end do                                                          11d2s20
       end if                                                           5d12s10
       if(iapair(1,ibc(ja8)).gt.0)then                                  5d12s10
        do isd=1,nsymb                                                  5d12s10
         mcd=nfcnb(isd)                                                 1d25s22
         do in=0,3                                                      12d23s21
          do ia=0,nsza0-1                                               11d2s20
           do ib=0,nszba0-1                                             11d2s20
            it1=itmpab(isd)+mcd*(ia+nszaa0*(ib+nszba0*in))              11d2s20
            it2=itmpab(isd)+mcd*(ia+nsza0+nszaa0*(ib+nszba0*in))        11d2s20
            do icd=0,mcd-1                                              11d2s20
             sum=srh*(bc(it1+icd)+bc(it2+icd))                           12d11s18
             diff=srh*(-bc(it1+icd)+bc(it2+icd))                         12d11s18
             bc(it1+icd)=sum                                             12d11s18
             bc(it2+icd)=diff                                            12d11s18
            end do                                                       12d11s18
           end do                                                        12d7s18
          end do                                                         5d12s10
         end do                                                          5d12s10
        end do                                                          11d2s20
       end if                                                           5d12s10
c
c     final storage                                                     1d10s11
c
       nmei=0                                                           12d11s18
       do ia=1,nszaa0                                                   8d26s15
        iaa=ia+ibc(ja4)                                                 2d3s12
        isa=isstor(iaa)                                                 2d3s12
        iaaa=ibstor(iaa)                                                12d11s18
        iaaap=iaaa+nbasdws(isa)                                         12d11s18
        do ib=1,nszba0                                                  8d26s15
         ibb=ib+ibc(jb4)                                                2d3s12
         isb=isstor(ibb)                                                2d3s12
         ibbb=ibstor(ibb)                                               12d12s18
         if(isb.ge.isa)then                                             12d23s21
          isc=multh(isymw,multh(isa,isb))
          ii=iptoh(isa,isb)                                             12d23s21
          mcd=nfcnb(isc)                                                1d25s22
          mcdn=mcd*ncsza0                                               12d27s21
          mcdv=mcd*nbasisp(isb)                                         12d27s21
          nvn=nvirt(isa)*nbasisp(isb)                                   12d27s21
c     thus tmpab is 1 LL, 2 SL, 3 LS, 4 SS
          if(mcdn.gt.0)then                                             12d27s21
           jorb1=iorbt(isa)+nvirt(isa)*(iaaa-1)                         12d23s21
           jorb2=iorbt(isa)+nvirt(isa)*(iaaa-1+nbasisp(isa))            12d23s21
           jtmp11=itmpab(isc)+mcd*(ia-1+nszaa0*(ib-1))                  12d30s21
           jtmp21=jtmp11+mcdn                                           12d27s21
           jtmp12=jtmp21+mcdn                                           12d27s21
           jtmp22=jtmp12+mcdn                                           12d27s21
           iip=ii+nvirt(isa)*(ibbb-1)                                   12d27s21
           iiq=iip+nvn*mcd                                              12d30s21
           do ixi=0,mcd-1                                                12d27s21
            jjp=iiq+nvn*ixi                                             12d30s21
            do jj=0,nvirt(isa)-1                                        12d27s21
             bc(jjp+jj)=bc(jjp+jj)+bc(jorb1+jj)*bc(jtmp12+ixi)          12d27s21
            end do                                                      12d27s21
           end do                                                       12d27s21
           do ixi=0,mcd-1                                                12d27s21
            jjp=iip+nvn*ixi                                             12d30s21
            do jj=0,nvirt(isa)-1                                        12d27s21
             bc(jjp+jj)=bc(jjp+jj)+bc(jorb2+jj)*bc(jtmp21+ixi)          12d30s21
            end do                                                      12d27s21
           end do                                                       12d27s21
           if(idorel.lt.0)then                                          3d2s22
            do ixi=0,mcd-1                                                12d27s21
             jjp=iip+nvn*ixi                                            12d30s21
             do jj=0,nvirt(isa)-1                                        12d27s21
              bc(jjp+jj)=bc(jjp+jj)+bc(jorb1+jj)*bc(jtmp11+ixi)         12d30s21
             end do                                                      12d27s21
            end do                                                       12d27s21
            do ixi=0,mcd-1                                                12d27s21
             jjp=iiq+nvn*ixi                                            12d30s21
             do jj=0,nvirt(isa)-1                                        12d27s21
              bc(jjp+jj)=bc(jjp+jj)+bc(jorb2+jj)*bc(jtmp22+ixi)           12d27s21
             end do                                                      12d27s21
            end do                                                       12d27s21
           end if                                                       12d23s21
          end if                                                        12d23s21
         end if                                                         12d12s18
        end do                                                          5d12s10
       end do                                                           5d12s10
       ibcoff=itmpab(1)                                                 12d7s18
c
c     end of do i loop                                                  11d27s18
c
      end do                                                            2d19s10
      call dws_gsumf(bc(ibptoh),needt)                                  12d12s18
      ioff=1                                                            12d23s21
      frac=1d0/dfloat(mynprocg)                                         2d1s22
      do isb=1,nsymb                                                    12d12s18
       isbv12=multh(isb,isymw)                                          12d23s21
       do isbv1=1,nsymb                                                 12d12s18
        isbv2=multh(isbv1,isbv12)                                       12d12s18
        if(isbv2.ge.isbv1)then                                          12d12s18
         nvvh=nvirt(isbv1)*nbasisp(isbv2)                                11d2s20
         mcd=nfcnb(isb)                                                 1d25s22
         if(mcd*min(nvirt(isbv2),nvvh).gt.0)then                        11d5s20
c     we have v1,b,ll,k
          itmpt=ibcoff                                                  12d23s21
          ibcoff=itmpt+nvvh*mcd*2                                       12d23s21
          call enough('paraerikbk4. 10',bc,ibc)
          sz=0d0
          do k=0,1                                                      12d23s21
           do i=0,mcd-1                                                 12d23s21
            do ib=0,nbasisp(isbv2)-1                                     12d23s21
             do ia=0,nvirt(isbv1)-1                                      12d23s21
              iad1=iptoh(isbv1,isbv2)+ia+nvirt(isbv1)*(ib+nbasisp(isbv2) 12d23s21
     $            *(i+mcd*k))                                           12d23s21
              iad2=itmpt+ib+nbasisp(isbv2)*(k+2*(ia+nvirt(isbv1)*i))     12d23s21
              bc(iad2)=bc(iad1)                                          12d23s21
             end do                                                     12d23s21
            end do                                                      12d23s21
           end do                                                       12d23s21
          end do                                                        12d23s21
          itmpa=ibcoff                                                  12d12s18
          ibcoff=itmpa+mcd*nvirt(isbv1)*nvirt(isbv2)                    12d12s18
          call enough('paraerikbk4. 11',bc,ibc)
          ncol=mcd*nvirt(isbv1)                                         12d23s21
          call dgemm('n','n',nvirt(isbv2),ncol,nbasispc(isbv2),frac,    2d1s22
     $         bc(iorbt(isbv2)),nvirt(isbv2),bc(itmpt),nbasispc(isbv2), 12d23s21
     $         0d0,bc(itmpa),nvirt(isbv2),                              12d23s21
     d' paraerikbk4.  7')
          if(isbv1.eq.isbv2)then                                        12d23s21
           do k=0,nfdat(2,1,isb)-1                                      12d23s21
            do ir=0,nrootu-1
             do iv=0,nvirt(isbv1)-1                                     12d23s21
              iad1=itmpa+iv+nvirt(isbv1)*(iv+nvirt(isbv1)*(ir+nrootu*k))12d23s21
              iad2=ioff+iv+nvirt(isbv1)*(ir+nrootu*k)                   12d23s21
              gd(iad2)=gd(iad2)+bc(iad1)*srh                            1d27s22
             end do                                                     12d23s21
            end do                                                      12d23s21
           end do                                                       12d23s21
           ioff=ioff+nvirt(isbv1)*nrootu*nfdat(2,1,isb)                 12d23s21
           isw=0                                                        12d23s21
           nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        12d23s21
          else                                                          12d23s21
           isw=1                                                        12d23s21
           nvv=nvirt(isbv1)*nvirt(isbv2)                                12d23s21
          end if                                                        12d23s21
          jtmpa=itmpa                                                   12d23s21
          do l=1,4                                                      12d23s21
           if(nfdat(2,l,isb).gt.0)then                                  12d23s21
            do k=0,nfdat(2,l,isb)-1                                     12d23s21
             do ir=0,nrootu-1                                           12d23s21
              do iv2=0,nvirt(isbv2)-1                                   12d23s21
               itop=(iv2+isw*(nvirt(isbv1)-iv2))-1                      12d23s21
               do iv1=0,itop                                            12d23s21
                iad1=jtmpa+iv2+nvirt(isbv2)*(iv1+nvirt(isbv1)*(ir       12d23s21
     $               +nrootu*k))                                        12d23s21
                itri=((iv2*(iv2-1))/2)+iv1                              12d23s21
                irec=iv1+nvirt(isbv1)*iv2                               12d23s21
                iad2=ioff+itri+isw*(irec-itri)+nvv*(ir+nrootu*k)        12d27s21
                gd(iad2)=gd(iad2)+bc(iad1)                              12d23s21
               end do                                                   12d23s21
              end do                                                    12d23s21
             end do                                                     12d23s21
            end do                                                      12d23s21
            ioff=ioff+nvv*nrootu*nfdat(2,l,isb)                         12d23s21
            iorig=jtmpa
            jtmpa=jtmpa+nvirt(isbv2)*nvirt(isbv1)*nrootu*nfdat(2,l,isb) 12d23s21
           end if                                                       12d23s21
          end do                                                        12d23s21
          ibcoff=itmpt                                                  12d23s21
         end if                                                         12d12s18
        end if                                                          12d12s18
       end do                                                           12d12s18
      end do                                                            12d12s18
      return
      end                                                               2d19s10
