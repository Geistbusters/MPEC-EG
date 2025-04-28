c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine paraeri42cf(natom,ngaus,ibdat,                          3d5s20
     $                ibstor,isstor,idwsdeb,ascale,nbasisp,             9d23s21
     $     iorb,iapair,idoubo,irefo,idorel,isnd,nsnd,ircv,nrcv,ihrow,   9d22s21
     $     multh,i4o,ionex,jmats,kmats,kmatsb,i3x,shift,ih0n,isopt,     11d15s21
     $     nsopt,scals,srh,isym,iifmx,ntype,nh0,iall,bc,ibc)            11d9s22
c
c     compute two e integrals for 2 component                           3d5s20
c     spin-free spinors in the full 4 component space.                  3d4s20
c
      implicit real*8 (a-h,o-z)
      external second
      integer*8 nn8                                                     2d19s10
      include "common.hf"
      include "common.store"
      include "common.spher"
      include "common.print"                                            3d2s22
      parameter (nfmx=1344,nfmx2=32,nfmxl=6,nfmxt=265)                            9d23s21
      logical logab,logcd,lprint                                        3d2s22
      integer*8 ibstor,isstor,i08
      integer*1 ifmx(10,nfmx),ifmxd(9,nfmx),iltype(4,nfmxl)             9d23s21
c     17*8 = 80+56=136 +3 = 139
      character*139 line
      character*2 xcode(5)
      character*1 alphab(2),rori(2)                                     3d8s20
      character*2 aaorbb(4)                                             9d21s21
      include "common.basis"
      data alphab/'+','-'/                                              3d8s20
      data aaorbb/'aa','ab','ba','bb'/                                  9d21s21
      data rori/'R','I'/
      data loopx/10000000/
      data xcode/'4o','1x','2x','3x','4x'/
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension irtyp(5),ibstor(*),isstor(*),iorb(*),idoubo(*),irefo(*),12d9s22
     $     nbasisp(*),iapair(3,*),nbasisp4(8),nbasisp2(8),              9d24s21
     $     nc(8),ihalf(8,8,8,4),jhalf(8,8,8,4),isnd(*),                 9d24s21
     $     nsnd(*),nhrow(8,8),ircv(*),nrcv(*),ihrow(*),multh(8,8),      10d7s21
     $     i4o(8,8,8,*),ionex(8,8,8,*),jmats(8,8,8,*),kmats(8,8,8,*),   10d5s21
     $     i3x(8,8,8,*),fmulx(4),ih0n(8,*),nc2(8),irefo2(8),            1d13s23
     $     iltyp(2),isopt(4,4),jfmx(3,nfmx),carta(3),cartb(3),cartc(3), 9d23s21
     $     cartd(3),ieraw(8,8,8,4),itmpc(8,8),iifmx(32,4),isym(3,8),    10d8s21
     $     nltype(nfmxl),neraw(8,8),nerimat(4),nerityp(2,nfmxt),        9d24s21
     $     meraw(8,8),ntype(4),nhalf(8,8,8,4),nd1(8),nd2(8),nd12(8),    10d12s21
     $     iqptr(32),ih0t(8,4),nh0(*),iall(8,8,8,*),kmatsb(8,8,8,*)     12d27s21
      dimension nkmats(8,8,8,4,4),njmats(8,8,8,4,4),n3x(8,8,8,4,4),
     $     igoalq(50),nallx(8,8,8,4,4)
      dimension ncode(4,3),icode(72,4,4,3),ilcode(4,6),iabcode(4,16)    11d5s21
      data ilcode/0,0,1,1, 1,1,0,0, 0,1,0,1, 0,1,1,0, 1,0,0,1, 1,0,1,0/  11d3s21
      data iabcode/0,0,0,0, 1,1,0,0, 0,0,1,0, 1,1,1,0, 0,1,0,0, 1,0,0,0,
     $    0,1,1,0, 1,0,1,0, 0,0,1,1, 1,1,1,1, 0,0,0,1, 1,1,0,1, 0,1,1,1,
     $    1,0,1,1, 0,1,0,1, 1,0,0,1/
      data (ncode(j,1),j=1,4)/8,8,8,0/
      data (icode(j,1,1,1),j=1,8)
     $     /0, 0, 0, 0, 1, 1, 1, 1/
      data (icode(j,2,1,1),j=1,8)
     $     /1, 1, 0, 0, 1, 1,-1, -1/
      data (icode(j,3,1,1),j=1,8)
     $     /1,1,1,1,2,2,2,2/
      data (icode(j,4,1,1),j=1,8)
     $     /1, 2, 7, 8, 1, 9, 2, 10/                                    3d2s22
      data (icode(j,1,2,1),j=1,8)
     $     /2, 2, 2, 2, 3, 3, 3, 3/
      data (icode(j,2,2,1),j=1,8)
     $     /1,-1, 1,-1, 1, 1,-1,-1/
      data (icode(j,3,2,1),j=1,8)
     $     /1,1,1,1,2,2,2,2/
      data (icode(j,4,2,1),j=1,8)
     $     /11, 3,12, 4, 5,13, 6,14/
      data (icode(j,1,3,1),j=1,8)
     $     /4, 4, 4, 4, 5, 5, 5, 5/
      data (icode(j,2,3,1),j=1,8)
     $     /1, 1, 1, 1, 1, 1, 1, 1/
      data (icode(j,3,3,1),j=1,8)
     $     /1,1,1,1,2,2,2,2/
      data (icode(j,4,3,1),j=1,8)
     $     /11, 3,12, 4, 5,13, 6,14/
      data (ncode(j,2),j=1,4)/72,3*64/
      data (icode(j,1,1,2),j=1,72)
     $     /12,12,12,12,12,12,16,16,16,16,16,16,20,20,20,20,20,20,21,21,
     $21,21,21,21,25,25,25,25,25,25,29,29,29,29,29,29,30,30,30,30,
     $30,30,34,34,34,34,34,34,38,38,38,38,38,38,39,39,39,39,39,39,
     $43,43,43,43,43,43,47,47,47,47,47,47/
      data (icode(j,2,1,2),j=1,72)
     $     /4, 2,-2,-2, 2, 4, 4,-2,-2,-2,-2, 4, 2, 2,-4,-4, 2, 2,-4, 2,
     $-2,-2, 2,-4,-4,-2,-2,-2,-2,-4,-2,-2,-4,-4,-2,-2,-4, 2,-2,-2,
     $ 2,-4,-4,-2,-2,-2,-2,-4,-2,-2,-4,-4,-2,-2, 4, 2,-2,-2, 2, 4,
     $ 4,-2,-2,-2,-2, 4, 2, 2,-4,-4, 2, 2/
      data (icode(j,3,1,2),j=1,72)
     $     /3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,
     $4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,
     $6,6,6,6,6,6,6,6,6,6,6,6/
      data (icode(j,4,1,2),j=1,72)
     $     /9,15, 7,16, 8, 2, 9,15, 7,16, 8, 2, 1, 9, 7,16, 2,10, 1,15,
     $ 7,16, 8,10, 1,15, 7,16, 8,10, 1, 9, 7,16, 2,10, 1,15, 7,16,
     $ 8,10, 1,15, 7,16, 8,10, 1, 9, 7,16, 2,10, 9,15, 7,16, 8, 2,
     $ 9,15, 7,16, 8, 2, 1, 9, 7,16, 2,10/
      data (icode(j,1,2,2),j=1,64)
     $     /0, 0, 0, 0, 2, 2, 2, 2, 6, 6, 6, 6, 8, 8, 8, 8,13,13,13,13,
     $13,13,15,15,15,15,15,15,22,22,22,22,22,22,24,24,24,24,24,24,
     $31,31,31,31,31,31,33,33,33,33,33,33,40,40,40,40,40,40,42,42,
     $42,42,42,42/
      data (icode(j,2,2,2),j=1,64)
     $     /2,-2, 2,-2,-2, 2,-2, 2, 2, 2,-2,-2,-2,-2, 2, 2,-4,-2,-2, 2,
     $ 2, 4, 4,-2, 2,-2, 2,-4, 4,-2,-2, 2, 2,-4,-4,-2, 2,-2, 2, 4,
     $-4,-2,-2, 2, 2, 4, 4,-2, 2,-2, 2,-4, 4,-2,-2, 2, 2,-4,-4,-2,
     $ 2,-2, 2, 4/
      data (icode(j,3,2,2),j=1,64)
     $     /1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,4,4,
     $4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,
     $6,6,6,6/
      data (icode(j,4,2,2),j=1,64)
     $     /1, 9, 2,10, 1, 9, 2,10, 1, 9, 2,10, 1, 9, 2,10, 9,15, 7,16,
     $ 8, 2, 9,15, 7,16, 8, 2, 1,15, 7,16, 8,10, 1,15, 7,16, 8,10,
     $ 1,15, 7,16, 8,10, 1,15, 7,16, 8,10, 9,15, 7,16, 8, 2, 9,15,
     $ 7,16, 8, 2/
      data (icode(j,1,3,2),j=1,64)
     $     /1, 1, 1, 1, 4, 4, 4, 4, 7, 7, 7, 7,10,10,10,10,14,14,14,14,
     $14,14,18,18,18,18,18,18,23,23,23,23,23,23,27,27,27,27,27,27,
     $32,32,32,32,32,32,36,36,36,36,36,36,41,41,41,41,41,41,45,45,
     $45,45,45,45/
      data (icode(j,2,3,2),j=1,64)
     $     /-2, 2,-2, 2, 2,-2, 2,-2,-2,-2, 2, 2, 2, 2,-2,-2, 4, 2, 2,-2,
     $-2,-4, 2,-2,-4, 4, 2,-2, 4,-2,-2, 2, 2,-4, 2,-2, 4,-4, 2,-2,
     $ 4, 2, 2,-2,-2,-4,-2, 2,-4, 4,-2, 2, 4,-2,-2, 2, 2,-4,-2, 2,
     $ 4,-4,-2, 2/
      data (icode(j,3,3,2),j=1,64)
     $     /1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,4,4,
     $4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,
     $6,6,6,6/
      data (icode(j,4,3,2),j=1,64)
     $     /11, 3,12, 4,11, 3,12, 4, 5,13, 6,14, 5,13, 6,14, 3, 5,13, 6,
     $14,12,11, 3,13, 6,12, 4, 3, 5,13, 6,14,12,11, 3, 5,14,12, 4,
     $11, 5,13, 6,14, 4,11, 3,13, 6,12, 4,11, 5,13, 6,14, 4,11, 3,
     $ 5,14,12, 4/
      data (icode(j,1,4,2),j=1,64)
     $     /3, 3, 3, 3, 5, 5, 5, 5, 9, 9, 9, 9,11,11,11,11,17,17,17,17,
     $17,17,19,19,19,19,19,19,26,26,26,26,26,26,28,28,28,28,28,28,
     $35,35,35,35,35,35,37,37,37,37,37,37,44,44,44,44,44,44,46,46,
     $46,46,46,46/
      data (icode(j,2,4,2),j=1,64)
     $     /2, 2, 2, 2,-2,-2,-2,-2, 2, 2, 2, 2,-2,-2,-2,-2, 4,-2,-2,-2,
     $-2, 4,-2,-2, 4, 4,-2,-2, 4, 2, 2, 2, 2, 4,-2,-2,-4,-4,-2,-2,
     $-4,-2,-2,-2,-2,-4, 2, 2, 4, 4, 2, 2,-4, 2, 2, 2, 2,-4, 2, 2,
     $-4,-4, 2, 2/
      data (icode(j,3,4,2),j=1,64)
     $     /1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,4,4,
     $4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,
     $6,6,6,6/
      data (icode(j,4,4,2),j=1,64)
     $     /11, 3,12, 4,11, 3,12, 4, 5,13, 6,14, 5,13, 6,14, 3, 5,13, 6,
     $14,12,11, 3,13, 6,12, 4, 3, 5,13, 6,14,12,11, 3, 5,14,12, 4,
     $11, 5,13, 6,14, 4,11, 3,13, 6,12, 4,11, 5,13, 6,14, 4,11, 3,
     $ 5,14,12, 4/
      data (ncode(j,3),j=1,4)/32,40,40,40/
      data (icode(j,1,1,3),j=1,32)
     $     /0, 1, 2, 2, 1, 0, 3, 3, 4, 5, 6, 6, 5, 4, 7, 7, 8, 9,10,10,
     $      9, 8,11,11,12,13,14,14,13,12,15,15/
      data (icode(j,2,1,3),j=1,32)
     $     /32*1/
      data (icode(j,3,1,3),j=1,32)
     $     /3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,6,6,6,6,6,6,
     $     6,6/
      data (icode(j,4,1,3),j=1,32)
     $     /9,15, 7,16, 8, 2, 1,10, 1,15, 7,16, 8,10, 9, 2, 1,15, 7,16,
     $      8,10, 9, 2, 9,15, 7,16, 8, 2, 1,10/
      data (icode(j,1,2,3),j=1,40)
     $     /16,16,16,16,17,17,17,17,18,18,19,20,21,21,20,19,22,22,23,24,
     $      25,25,24,23,26,26,27,28,29,29,28,27,30,30,31,32,33,33,32,31/
      data (icode(j,2,2,3),j=1,40)
     $     /1,-1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1, 1,
     $      1,-1,-1,-1, 1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1, 1, 1,-1,-1,-1/
      data (icode(j,3,2,3),j=1,40)
     $     /1,1,1,1,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,
     $      5,5,6,6,6,6,6,6,6,6/
      data (icode(j,4,2,3),j=1,40)
     $     /1, 9, 2,10, 1, 9, 2,10, 1,10, 9,15, 7,16, 8, 2, 9, 2, 1,15,
     $      7,16, 8,10, 9, 2, 1,15, 7,16, 8,10, 1,10, 9,15, 7,16, 8, 2/
      data (icode(j,1,3,3),j=1,40)
     $     /34,34,34,34,35,35,35,35,36,37,38,39,39,38,37,36,40,41,42,43,
     $      43,42,41,40,44,45,46,47,47,46,45,44,48,49,50,51,51,50,49,48/
      data (icode(j,2,3,3),j=1,40)
     $     /1,-1, 1,-1, 1, 1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1, 1, 1, 1, 1,
     $     -1,-1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1/
      data (icode(j,3,3,3),j=1,40)
     $     /1,1,1,1,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,
     $      5,5,6,6,6,6,6,6,6,6/
      data (icode(j,4,3,3),j=1,40)
     $     /11, 3,12, 4, 5,13, 6,14,11, 3, 5,13, 6,14,12, 4,11, 3, 5,13,
     $       6,14,12, 4,11, 3, 5,13, 6,14,12, 4,11, 3, 5,13, 6,14,12, 4/
      data (icode(j,1,4,3),j=1,40)
     $     /52,52,52,52,53,53,53,53,54,55,56,57,57,56,55,54,58,59,60,61,
     $      61,60,59,58,62,63,64,65,65,64,63,62,66,67,68,69,69,68,67,66/
      data (icode(j,2,4,3),j=1,40)
     $     /40*1/
      data (icode(j,3,4,3),j=1,40)
     $     /1,1,1,1,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,
     $      5,5,6,6,6,6,6,6,6,6/
      data (icode(j,4,4,3),j=1,40)
     $     /11, 3,12, 4, 5,13, 6,14,11, 3, 5,13, 6,14,12, 4,11, 3, 5,13,
     $      6,14,12, 4,11, 3, 5,13, 6,14,12, 4,11, 3, 5,13, 6,14,12, 4/
      data icall/0/
      save
      clight=0.5d0/sqrt(ascale)                                         2d21s20
      icall=icall+1
      lprint=mynowprog.eq.0                                             3d2s22
      fmulx(1)=1d0                                                      3d13s20
      fmulx(2)=1d0
      fmulx(3)=1d0
      fmulx(4)=1d0
      if(icall.eq.1)then                                                9d21s21
c
c     overall factor to multiple ints is ifmx(10,i)/(8*clight*clight)   9d27s21
c
       dgemmf=1d0/(8d0*clight*clight)                                   9d27s21
       iboost=1
       dgemmf=dgemmf/dfloat(iboost)
       nadd=0
       iuse=1
       nfoff=0                                                          11d4s21
       mycode=7                                                         11d4s21
       if(idorel.eq.-4)then
        iuse=2
        if(nsopt.eq.1)then                                              11d4s21
         mycode=9                                                       11d5s21
        else                                                            11d4s21
         mycode=8                                                       11d5s21
        end if                                                          11d4s21
       else if(idorel.lt.0)then                                         11d3s21
        iuse=3                                                          11d3s21
        if(nsopt.eq.1)then                                              11d4s21
         mycode=10                                                      11d4s21
        else                                                            11d4s21
         mycode=11                                                      11d4s21
        end if                                                          11d4s21
       end if                                                           11d3s21
       do it=1,nsopt                                                    11d3s21
        nfx=0
        nerimat(it)=0                                                   11d5s21
        do i=1,ncode(it,iuse)                                           11d3s21
         ip=i+nadd
         i3=icode(i,3,it,iuse)
         if(ilcode(1,i3).eq.
     $      ilcode(2,i3))icode(i,2,it,iuse)=
     $        icode(i,2,it,iuse)*iboost
         jfmx(1,ip)=icode(i,1,it,iuse)                                  11d4s21
         do j=1,nerimat(it)                                             11d5s21
          if(jfmx(1,ip).eq.jfmx(3,j))go to 2022                         11d5s21
         end do                                                         11d5s21
         nfoff=jfmx(1,ip)+1                                             11d5s21
         nerityp(1,nfoff)=it                                            11d5s21
         nerityp(2,nfoff)=nerimat(it)                                   11d5s21
         nerimat(it)=nerimat(it)+1                                      11d5s21
         jfmx(3,nerimat(it))=jfmx(1,ip)                                 11d5s21
 2022    continue                                                       11d5s21
         ifmx(1,ip)=it                                                  11d4s21
         do j=1,4                                                       11d4s21
          jp=j+1                                                        11d4s21
          ifmx(jp,ip)=ilcode(j,icode(i,3,it,iuse))                      11d4s21
          jp=jp+4                                                       11d4s21
          ifmx(jp,ip)=iabcode(j,icode(i,4,it,iuse))                     11d4s21
         end do                                                         11d4s21
         ifmx(10,ip)=icode(i,2,it,iuse)                                 11d4s21
 2020    format(i5,5x,2i3,5x,4i2,5x,4i2,5x,2i8)
        end do                                                          11d3s21
        nadd=nadd+ncode(it,iuse)
       end do                                                           11d3s21
       mfmx=nadd                                                        11d4s21
c
c     nerimat=265, ntype=152
c
       do i=1,nsopt                                                     9d24s21
        ntype(i)=0                                                      9d24s21
       end do                                                           9d24s21
       do jpass=0,1                                                     12d3s21
        do i=1,mfmx                                                      9d23s21
         if(ifmx(9,i).eq.jpass)then                                     12d3s21
          itest=ifmx(6,i)+2*(ifmx(7,i)+2*(ifmx(8,i)+2*ifmx(9,i)))         9d24s21
          it=ifmx(1,i)                                                    9d24s21
          do j=1,ntype(it)                                                9d24s21
           if(itest.eq.iifmx(j,it))then                                   9d24s21
            jfmx(2,i)=j                                                   9d23s21
            go to 1517                                                    9d23s21
           end if                                                         9d23s21
          end do                                                          9d23s21
          ntype(it)=ntype(it)+1                                           9d24s21
          iifmx(ntype(it),it)=itest                                       9d24s21
          jfmx(2,i)=ntype(it)                                             9d24s21
 1517     continue                                                        9d23s21
         end if                                                         12d3s21
        end do                                                          12d3s21
       end do                                                           9d23s21
       ltype=0                                                          9d23s21
       nds=0                                                            9d23s21
       ncs=0
       nbs=0                                                            9d23s21
       nas=0                                                            9d23s21
       do lsd=0,1                                                       9d23s21
        do lsc=0,1                                                      9d23s21
         do lsb=0,1                                                     9d23s21
          do lsa=0,1                                                    9d23s21
           ihit=0                                                       9d23s21
           do i=1,mfmx                                                  9d23s21
            if(ifmx(2,i).eq.lsa.and.ifmx(3,i).eq.lsb.and.ifmx(4,i).     9d23s21
     $           eq.lsc.and.ifmx(5,i).eq.lsd)then                       9d23s21
             jfmx(3,i)=ltype+1                                          9d24s21
             ihit=ihit+1                                                9d23s21
            end if                                                      9d23s21
           end do                                                       9d23s21
           if(ihit.gt.0)then                                            9d23s21
            if(lsd.eq.0)then                                            9d23s21
             nds=nds+1                                                  9d23s21
             if(lsc.eq.0)then                                           9d23s21
              ncs=ncs+1                                                 9d23s21
              if(lsb.eq.0)then                                          9d23s21
               nbs=nbs+1                                                9d23s21
               if(lsa.eq.0)nas=nas+1                                    9d23s21
              end if                                                    9d23s21
             end if                                                     9d23s21
            end if                                                      9d23s21
            ltype=ltype+1                                               9d23s21
            if(ltype.gt.nfmxl)then                                      9d23s21
             call dws_synca                                             9d23s21
             call dws_finalize                                          9d23s21
             write(6,*)('too many ltypes!')                             9d23s21
             stop 'paraeri42cf'                                         9d23s21
            end if                                                      9d23s21
            nltype(ltype)=ihit                                          9d23s21
            iltype(1,ltype)=lsa                                         9d23s21
            iltype(2,ltype)=lsb                                         9d23s21
            iltype(3,ltype)=lsc                                         9d23s21
            iltype(4,ltype)=lsd                                         9d23s21
           end if                                                       9d23s21
          end do                                                        9d23s21
         end do                                                         9d23s21
        end do                                                          9d23s21
       end do                                                           9d23s21
      end if                                                            9d21s21
      nerimatt=0                                                        9d24s21
      do i=1,nsopt                                                      9d24s21
       nerimatt=nerimatt+nerimat(i)                                     9d24s21
      end do                                                            9d24s21
      val=-1d0                                                          3d17s20
c     ascale=1/(4*c*c), thus c = 0.5/sqrt(ascale)
c
      do isb=1,nsymb                                                    3d5s20
c
       nc(isb)=idoubo(isb)+irefo(isb)                                   3d6s20
       nc2(isb)=nc(isb)*2                                               3d24s20
       irefo2(isb)=irefo(isb)*2                                         3d24s20
       nbasisp2(isb)=nbasisp(isb)*2
       nbasisp4(isb)=nbasisp2(isb)*2
       nd1(isb)=nbasdws(isb)                                            10d20s21
       nd12(isb)=nd1(isb)*2                                             10d4s21
       nd2(isb)=nbasdws(isb)                                            10d4s21
       nh0(isb)=nbasdws(isb)-idoubo(isb)                                10d13s21
      end do                                                            3d5s20
      ibcoffo=ibcoff                                                    2d19s10
      if(idwsdeb.gt.100)then                                            5d25s18
       write(6,*)('in paraeri42c ')
       write(6,*)('speed of light: '),clight
       write(6,*)('ibdat '),ibdat
       write(6,*)('isstor: '),(isstor(i),i=1,28)
       do i=0,ngaus-1
        ip=i+1
        write(6,*)ip,ibc(ibdat+i),(bc(ibdat+i+ngaus*j),j=1,2),
     $      ibc(ibdat+i+3*ngaus),(bc(ibdat+i+ngaus*j),j=4,6),
     $      (ibc(ibdat+i+ngaus*j),j=7,8)
       end do
      end if                                                            5d25s18
      call second(time1)
c
c     build shell pair order
c
      ascale2=ascale*2d0                                                8d20s15
      nn=ngaus*ngaus                                                    9d23s21
      nn2=nn*2                                                          2d19s10
c
c     for large component,
c     overlap and nuclear attraction for each atom
c     for small component,
c     overlap and nuclear attraction for each atom
c     for ls coupling, dx,dy,dz
c     for sl coupling, dx,dy,dz
c
      ipair=ibcoff                                                      2d19s10
      ibcoff=ipair+nn*2                                                 3d6s20
      i12=ibcoff                                                        2d19s10
      ibcoff=i12+nn*4                                                   2d19s10
      call enough('paraeri42cf.  1',bc,ibc)
      j12=i12                                                           2d19s10
      jbdat=ibdat-1                                                     2d19s10
      ngaus7=ngaus*7                                                    5d3s10
      ngaus3=ngaus*3                                                    9d24s21
      do i1=1,ngaus                                                     2d19s10
       do i2=1,ngaus                                                    9d23s21
        ibc(j12)=-((ibc(jbdat+i1)+ibc(jbdat+i2)+2)**2)                  2d19s10
        ibc(j12+nn)=i1                                                  2d19s10
        ibc(j12+nn2)=i2                                                 2d19s10
        j12=j12+1                                                       2d19s10
       end do                                                           2d19s10
      end do                                                            2d19s10
      is12=i12+nn*3                                                     2d19s10
      nn8=nn                                                            2d19s10
      call idsortdws(ibc(i12),ibc(is12),nn8)                            1d18s23
      ii1=i12+nn                                                        2d19s10
      ii2=i12+nn2                                                       2d19s10
      jpair=ipair                                                       2d19s10
      iqy=0                                                             3d4s20
      do i=1,nn                                                         2d19s10
       j=ibc(is12+i-1)-1                                                2d19s10
    1  format(i5,i8,2i5)                                                2d19s10
       iqy=iqy+1
       jpair0=jpair
       ibc(jpair)=ibc(ii1+j)                                            2d19s10
       jpair=jpair+1                                                    2d19s10
       ibc(jpair)=ibc(ii2+j)                                            2d19s10
       jpair=jpair+1                                                    2d19s10
      end do                                                            2d19s10
      nneed=nn                                                          3d6s20
      call second(time2)
      telap=time2-time1
      do isb=1,nsymb                                                    9d24s21
       do isa=1,nsymb                                                   9d24s21
        nhrow(isa,isb)=0                                                9d24s21
       end do                                                           9d24s21
      end do                                                            9d24s21
      do ip=1,mynprocg*nsymb*nsymb                                      9d24s21
       ihrow(ip)=0                                                      3d7s20
      end do                                                            3d7s20
      do i=0,nneed-1                                                    3d6s20
       ip=mod(i,mynprocg)                                               3d6s20
       ip=ip+1                                                          3d6s20
       jpair=ipair+2*i                                                  3d6s20
       ja=jbdat+ibc(jpair)                                              3d6s20
       ja4=ja+ngaus3                                                    9d24s21
       ja8=ja+ngaus7                                                    9d23s21
       nla=ibc(ja)*2+1                                                  3d6s20
       jb=jbdat+ibc(jpair+1)                                              3d6s20
       jb4=jb+ngaus3                                                    9d24s21
       jb8=jb+ngaus7                                                    9d23s21
       nlb=ibc(jb)*2+1                                                  3d6s20
       if(iapair(1,ibc(ja8)).gt.0)nla=nla*2                             9d23s21
       if(iapair(1,ibc(jb8)).gt.0)nlb=nlb*2                             9d23s21
       do ib=1,nlb                                                      9d24s21
        ibb=ib+ibc(jb4)                                                 9d24s21
        isb=isstor(ibb)                                                 9d24s21
        do ia=1,nla                                                     9d24s21
         iaa=ia+ibc(ja4)                                                9d24s21
         isa=isstor(iaa)                                                9d24s21
         nhrow(isa,isb)=nhrow(isa,isb)+1                                9d24s21
         iddr=isa+nsymb*(isb-1+nsymb*(ip-1))                            9d24s21
         ihrow(iddr)=ihrow(iddr)+1                                      9d24s21
        end do                                                          9d24s21
       end do                                                           9d24s21
      end do                                                            3d6s20
      do it=1,nsopt                                                     9d24s21
       do isd=1,nsymb                                                   9d24s21
        do isc=1,nsymb                                                  9d24s21
         do isb=1,nsymb                                                 9d24s21
          nhalf(isb,isc,isd,it)=0                                       9d24s21
         end do                                                         9d24s21
        end do                                                          9d24s21
       end do                                                           9d24s21
      end do                                                            9d24s21
      nabn=0                                                            3d6s20
      do i=1+mynowprog,nneed,mynprocg                                   3d6s20
       jpair=ipair+2*(i-1)                                              3d6s20
       ja=jbdat+ibc(jpair)                                              3d5s20
       ja8=ja+ngaus7                                                    9d23s21
       ja4=ja+ngaus3                                                    9d24s21
       jb=jbdat+ibc(jpair+1)                                            3d6s20
       jb4=jb+ngaus3                                                    9d24s21
       jb8=jb+ngaus7                                                    9d23s21
       nla=ibc(ja)*2+1                                                  3d6s20
       nlb=ibc(jb)*2+1                                                  3d6s20
       if(iapair(1,ibc(ja8)).gt.0)nla=nla*2                             9d23s21
       if(iapair(1,ibc(jb8)).gt.0)nlb=nlb*2                             9d23s21
       do isb=1,nsymb                                                   9d24s21
        do isa=1,nsymb                                                  9d24s21
         neraw(isa,isb)=0                                               9d24s21
        end do                                                          9d24s21
       end do                                                           9d24s21
       do ib=1,nlb                                                      9d24s21
        ibb=ib+ibc(jb4)                                                 9d24s21
        isb=isstor(ibb)                                                 9d24s21
        do ia=1,nla                                                     9d24s21
         iaa=ia+ibc(ja4)                                                9d24s21
         isa=isstor(iaa)                                                9d24s21
         neraw(isa,isb)=neraw(isa,isb)+1                                9d24s21
        end do                                                          9d24s21
       end do                                                           9d24s21
       do it=1,nsopt                                                    9d24s21
        do isd=1,nsymb                                                  9d24s21
         istd=multh(isopt(1,it),isd)                                    9d24s21
         do isc=1,nsymb                                                 9d24s21
          istdc=multh(istd,isc)                                         9d24s21
          do isb=1,nsymb                                                9d24s21
           isa=multh(istdc,isb)                                         9d24s21
           nhalf(isb,isc,isd,it)=nhalf(isb,isc,isd,it)+neraw(isa,isb)   9d24s21
          end do                                                        9d24s21
         end do                                                         9d24s21
        end do                                                          9d24s21
       end do                                                           9d24s21
      end do                                                            3d6s20
      ibchalf=ibcoff                                                    10d5s21
      do it=1,nsopt                                                     9d24s21
       ntype4=ntype(it)*4                                               9d24s21
       do isd=1,nsymb
        istd=multh(isopt(1,it),isd)                                     9d24s21
        do isc=1,nsymb                                                   3d6s20
         istdc=multh(istd,isc)                                          9d24s21
         do isb=1,nsymb                                                 9d24s21
          jhalf(isb,isc,isd,it)=ibcoff                                  9d24s21
          ibcoff=jhalf(isb,isc,isd,it)+nhalf(isb,isc,isd,it)            9d24s21
     $         *ntype4*nd1(isd)*nd2(isc)                                10d4s21
          ihalf(isb,isc,isd,it)=jhalf(isb,isc,isd,it)                   9d24s21
         end do                                                         9d24s21
        end do                                                          9d24s21
       end do                                                           3d6s20
      end do                                                            3d6s20
      nwhalf=ibcoff-ihalf(1,1,1,1)
      call enough('paraeri42cf.  2',bc,ibc)
      do iz=ihalf(1,1,1,1),ibcoff-1                                     9d24s21
       bc(iz)=0d0                                                       9d23s21
      end do                                                            9d23s21
      loopit=0                                                          5d25s18
      ngaus3=ngaus*3                                                    3d6s20
      tgeri=0d0                                                         3d16s20
      tssdvr=0d0                                                        3d16s20
      call second(time1)                                                3d16s20
      call second(time2)                                                3d16s20
      tovr=time2-time1                                                  3d16s20
      loop=0
            igoalq(1)=2215676
            igoalq(2)=2215661
            nkq=1
      do i=1+mynowprog,nneed,mynprocg                                   2d19s10
       loopit=loopit+1                                                  5d25s18
       jpair=ipair+2*(i-1)                                              3d6s20
       call second(time1)
       if(idwsdeb.gt.10)write(6,2)i,ibc(jpair),ibc(jpair+1),loop,time1             3d6s20
    2  format('I am going to do ',i5,3i3,f18.5)                               2d19s10
       logab=.false.
       ja=jbdat+ibc(jpair)                                              3d5s20
       ja2=ja+ngaus                                                     3d5s20
       ja3=ja2+ngaus                                                    3d5s20
       ja4=ja3+ngaus                                                    3d5s20
       ja5=ja4+ngaus                                                    3d5s20
       ja6=ja5+ngaus                                                    3d5s20
       ja7=ja6+ngaus                                                    3d5s20
       ja8=ja7+ngaus                                                    3d5s20
       jb=jbdat+ibc(jpair+1)                                            3d5s20
       jb2=jb+ngaus                                                     3d5s20
       jb3=jb2+ngaus                                                    3d5s20
       jb4=jb3+ngaus                                                    3d5s20
       jb5=jb4+ngaus                                                    3d5s20
       jb6=jb5+ngaus                                                    3d5s20
       jb7=jb6+ngaus                                                    3d5s20
       jb8=jb7+ngaus                                                    3d5s20
       nsza=2*ibc(ja)+1                                                  3d5s20
       nsza0=nsza                                                       9d22s21
       if(iapair(1,ibc(ja8)).gt.0)then                                  9d22s21
        nsza=nsza*2                                                       9d22s21
       end if                                                           9d22s21
       nszb=2*ibc(jb)+1                                                 9d22s21
       nszb0=nszb                                                       9d22s21
       if(iapair(1,ibc(jb8)).gt.0)then                                  9d22s21
        nszb=nszb*2                                                     9d22s21
       end if                                                           9d22s21
       do isb=1,nsymb                                                   9d24s21
        do isa=1,nsymb                                                  9d24s21
         neraw(isa,isb)=0                                               9d24s21
        end do                                                          9d24s21
       end do                                                           9d24s21
       do ib=1,nszb                                                     9d24s21
        ibb=ib+ibc(jb4)                                                 9d24s21
        isb=isstor(ibb)                                                 9d24s21
        do ia=1,nsza                                                    9d24s21
         iaa=ia+ibc(ja4)                                                9d24s21
         isa=isstor(iaa)                                                9d24s21
         neraw(isa,isb)=neraw(isa,isb)+1                                9d24s21
        end do                                                          9d24s21
       end do                                                           9d24s21
       ncsz=nsza*nszb                                                   5d11s10
       ncsz0=nsza0*nszb0                                                10d28s20
       iberaw=ibcoff                                                    10d28s21
       do it=1,nsopt                                                    10d28s21
        do isd=1,nsymb                                                  10d28s21
         istd=multh(isopt(1,it),isd)                                    10d28s21
         do isc=1,nsymb                                                 10d28s21
          istdc=multh(istd,isc)                                         10d28s21
          do isb=1,nsymb                                                10d28s21
           isa=multh(istdc,isb)                                         10d28s21
           ieraw(isb,isc,isd,it)=ibcoff                                 10d28s21
           ibcoff=ieraw(isb,isc,isd,it)+neraw(isa,isb)*nbasisp(isc)     10d28s21
     $          *nbasisp(isd)*nerimat(it)                               10d28s21
          end do                                                        10d28s21
         end do                                                         10d28s21
        end do                                                          10d28s21
       end do                                                           10d28s21
       call enough('paraeri42cf.  3',bc,ibc)
       do iz=iberaw,ibcoff-1                                            10d28s21
        bc(iz)=val                                                      10d28s21
       end do                                                           10d28s21
       do jjc=1,ngaus                                                   9d24s21
        jc=jbdat+jjc                                                    9d24s21
        jc2=jc+ngaus                                                     3d5s20
        jc3=jc2+ngaus                                                    3d5s20
        jc4=jc3+ngaus                                                    3d5s20
        jc5=jc4+ngaus                                                    3d5s20
        jc6=jc5+ngaus                                                    3d5s20
        jc7=jc6+ngaus                                                    3d5s20
        jc8=jc7+ngaus                                                    3d5s20
        nszc=2*ibc(jc)+1                                                  3d5s20
        nszc0=nszc                                                      9d24s21
        if(iapair(1,ibc(jc8)).gt.0)then                                 9d24s21
         nszc=nszc*2                                                    9d24s21
        end if                                                          9d24s21
        do jjd=1,ngaus                                                  9d24s21
         jd=jbdat+jjd                                                   9d24s21
         jd2=jd+ngaus                                                     3d5s20
         jd3=jd2+ngaus                                                    3d5s20
         jd4=jd3+ngaus                                                    3d5s20
         jd5=jd4+ngaus                                                    3d5s20
         jd6=jd5+ngaus                                                    3d5s20
         jd7=jd6+ngaus                                                    3d5s20
         jd8=jd7+ngaus                                                    3d5s20
         nszd=2*ibc(jd)+1                                                  3d5s20
         nszd0=nszd                                                     9d24s21
         if(iapair(1,ibc(jd8)).gt.0)then                                9d24s21
          nszd=nszd*2                                                   9d24s21
         end if                                                         9d24s21
         nwtmp=ncsz*nszc*nszd                                           9d24s21
         nweri=ncsz0*nszc0*nszd0                                        9d24s21
         itmp=ibcoff                                                    9d24s21
         ibcoff=itmp+nwtmp*nerimatt                                     9d24s21
         call enough('paraeri42cf.  4',bc,ibc)
         do ipassd=1,nszd/nszd0                                         9d24s21
          if(ipassd.eq.1)then                                           9d24s21
           cartd(1)=bc(jd5)                                             9d24s21
           cartd(2)=bc(jd6)                                             9d24s21
           cartd(3)=bc(jd7)                                             9d24s21
           iod=0                                                          5d12s10
          else                                                          9d24s21
           cartd(1)=bc(jd5)*dfloat(isym(1,iapair(2,ibc(jd8))))            5d11s10
           cartd(2)=bc(jd6)*dfloat(isym(2,iapair(2,ibc(jd8))))            5d11s10
           cartd(3)=bc(jd7)*dfloat(isym(3,iapair(2,ibc(jd8))))            5d11s10
           iod=nszd0                                                      10d28s20
          end if                                                        9d24s21
          do ipassc=1,nszc/nszc0                                         9d24s21
           if(ipassc.eq.1)then                                           9d24s21
            cartc(1)=bc(jc5)                                             9d24s21
            cartc(2)=bc(jc6)                                             9d24s21
            cartc(3)=bc(jc7)                                             9d24s21
            ioc=0                                                          5d12s10
           else                                                          9d24s21
            cartc(1)=bc(jc5)*dfloat(isym(1,iapair(2,ibc(jc8))))            5d11s10
            cartc(2)=bc(jc6)*dfloat(isym(2,iapair(2,ibc(jc8))))            5d11s10
            cartc(3)=bc(jc7)*dfloat(isym(3,iapair(2,ibc(jc8))))            5d11s10
            ioc=nszc0                                                      10d28s20
           end if                                                        9d24s21
           do ipassb=1,nszb/nszb0                                           9d22s21Q
            if(ipassb.eq.1)then
             cartb(1)=bc(jb5)                                               5d12s10
             cartb(2)=bc(jb6)                                               5d12s10
             cartb(3)=bc(jb7)                                               5d12s10
             iob=0                                                          5d12s10
            else                                                            5d12s10
             cartb(1)=bc(jb5)*dfloat(isym(1,iapair(2,ibc(jb8))))            5d11s10
             cartb(2)=bc(jb6)*dfloat(isym(2,iapair(2,ibc(jb8))))            5d11s10
             cartb(3)=bc(jb7)*dfloat(isym(3,iapair(2,ibc(jb8))))            5d11s10
             iob=nszb0                                                      10d28s20
            end if                                                          9d22s21
            do ipassa=1,nsza/nsza0                                          9d22s21
             if(ipassa.eq.1)then                                            5d12s10
              carta(1)=bc(ja5)                                              5d12s10
              carta(2)=bc(ja6)                                              8d15s11
              carta(3)=bc(ja7)                                              8d15s11
              ioa=0                                                         5d12s10
             else                                                           5d12s10
              carta(1)=bc(ja5)*dfloat(isym(1,iapair(2,ibc(ja8))))           5d12s10
              carta(2)=bc(ja6)*dfloat(isym(2,iapair(2,ibc(ja8))))           5d12s10
              carta(3)=bc(ja7)*dfloat(isym(3,iapair(2,ibc(ja8))))           5d12s10
              ioa=nsza0                                                     10d28s20
             end if                                                         9d22s21
             call dgerid(
     $            ibc(ja),bc(ja2),bc(ja3),carta,carta(2),carta(3),      9d24s21
     $            ibc(jb),bc(jb2),bc(jb3),cartb,cartb(2),cartb(3),      9d24s21
     $            ibc(jc),bc(jc2),bc(jc3),cartc,cartc(2),cartc(3),      9d24s21
     $            ibc(jd),bc(jd2),bc(jd3),cartd,cartd(2),cartd(3),      9d24s21
     $            0,0,0, 0,0,0, 0,0,0, 0,0,0,ieri,.false.,mycode,       9d24s21
     $            fmulx,bc,ibc)                                         11d14s22
             do in=0,nerimatt-1                                         9d24s21
              do id=0,nszd0-1                                              10d28s20
               idp=id+iod                                               9d24s21
               do ic=0,nszc0-1                                             10d28s20
                icp=ic+ioc                                              9d24s21
                do ib=0,nszb0-1                                           9d24s21
                 ibp=ib+iob                                               9d24s21
                 do ia=0,nsza0-1                                          9d24s21
                  iab=ib+nszb0*ia                                         9d24s21
                  iad1=itmp+in+nerimatt*(ia+ioa+nsza*(ibp+nszb*(icp     9d24s21
     $               +nszc*idp)))                                        9d24s21
                  iad2=ieri+id+nszd0*(ic+nszc0*(iab+ncsz0*in))             9d23s21
                  bc(iad1)=bc(iad2)                                       9d24s21
                 end do                                                       5d12s10
                end do                                                        5d12s10
               end do                                                         5d12s10
              end do                                                        10d28s20
             end do                                                       9d24s21
            end do                                                      9d24s21
           end do                                                       9d24s21
          end do                                                        9d24s21
         end do                                                         9d24s21
         if(nszd.ne.nszd0)then                                          9d24s21
          nrow=ncsz*nszc*nerimatt                                       9d24s21
          do id=0,nszd0-1                                               9d24s21
           it1=itmp+nrow*id                                             9d24s21
           it2=itmp+nrow*(id+nszd0)                                     9d24s21
           do iabc=0,nrow-1                                             9d24s21
            sum=srh*(bc(it1+iabc)+bc(it2+iabc))                         9d24s21
            diff=srh*(-bc(it1+iabc)+bc(it2+iabc))                       9d24s21
            bc(it1+iabc)=sum                                            9d24s21
            bc(it2+iabc)=diff                                           9d24s21
           end do                                                       9d24s21
          end do                                                        9d24s21
         end if                                                         9d24s21
         if(nszc.ne.nszc0)then                                          9d24s21
          nrow=ncsz*nerimatt                                            9d24s21
          do id=0,nszd-1                                                9d24s21
           do ic=0,nszc0-1                                              9d24s21
            it1=itmp+nrow*(ic+nszc*id)                                  9d24s21
            it2=itmp+nrow*(ic+nszc0+nszc*id)                            9d24s21
            do iab=0,nrow-1                                             9d24s21
             sum=srh*(bc(it1+iab)+bc(it2+iab))                          5d12s10
             diff=srh*(-bc(it1+iab)+bc(it2+iab))                        5d12s10
             bc(it1+iab)=sum                                            5d12s10
             bc(it2+iab)=diff                                           5d12s10
            end do                                                      5d12s10
           end do                                                       5d12s10
          end do                                                        5d12s10
         end if                                                         5d12s10
         if(nszb.ne.nszb0)then                                          9d24s21
          nrow=nsza*nerimatt                                            9d24s21
          ncol=nszc*nszd                                                9d24s21
          do icol=0,ncol-1                                              9d24s21
           do ib=0,nszb0-1                                              9d24s21
            it1=itmp+nrow*(ib+nszb*icol)                                9d24s21
            it2=itmp+nrow*(ib+nszb0+nszb*icol)                            9d24s21
            do irow=0,nrow-1                                             9d24s21
             sum=srh*(bc(it1+irow)+bc(it2+irow))                          5d12s10
             diff=srh*(-bc(it1+irow)+bc(it2+irow))                        5d12s10
             bc(it1+irow)=sum                                            5d12s10
             bc(it2+irow)=diff                                           5d12s10
            end do                                                      5d12s10
           end do                                                       5d12s10
          end do                                                        5d12s10
         end if                                                         5d12s10
         if(nsza.ne.nsza0)then                                          9d24s21
          nrow=nerimatt                                                 9d24s21
          ncol=nszb*nszc*nszd                                                9d24s21
          do icol=0,ncol-1                                              9d24s21
           do ia=0,nsza0-1                                              9d24s21
            it1=itmp+nrow*(ia+nsza*icol)                                9d24s21
            it2=itmp+nrow*(ia+nsza0+nsza*icol)                            9d24s21
            do irow=0,nrow-1                                             9d24s21
             sum=srh*(bc(it1+irow)+bc(it2+irow))                          5d12s10
             diff=srh*(-bc(it1+irow)+bc(it2+irow))                        5d12s10
             bc(it1+irow)=sum                                            5d12s10
             bc(it2+irow)=diff                                           5d12s10
            end do                                                      5d12s10
           end do                                                       5d12s10
          end do                                                        5d12s10
         end if                                                         5d12s10
         do id=1,nszd                                                   9d24s21
          idd=id+ibc(jd4)                                               9d24s21
          iddd=ibstor(idd)-1                                            9d24s21
          idm=id-1                                                      9d24s21
          isd=isstor(idd)                                               9d24s21
          do ic=1,nszc                                                  9d24s21
           icc=ic+ibc(jc4)                                              9d24s21
           iccc=ibstor(icc)-1                                           9d24s21
           icm=ic-1                                                     9d24s21
           isc=isstor(icc)                                              9d24s21
           iscd=multh(isc,isd)                                          9d24s21
           do isb=1,nsymb                                               9d24s21
            do isa=1,nsymb                                              9d24s21
             meraw(isa,isb)=0                                           9d24s21
            end do                                                      9d24s21
           end do                                                       9d24s21
           do ib=1,nszb                                                 9d24s21
            ibm=ib-1                                                    9d24s21
            ibb=ib+ibc(jb4)                                             9d24s21
            isb=isstor(ibb)                                             9d24s21
            isbcd=multh(iscd,isb)                                       9d24s21
            do ia=1,nsza                                                9d24s21
             iam=ia-1                                                   9d24s21
             iaa=ia+ibc(ja4)                                            9d24s21
             isa=isstor(iaa)                                            9d24s21
             isabcd=multh(isa,isbcd)                                    9d24s21
             jtmp=itmp-1+nerimatt*(iam+nsza*(ibm+nszb*(icm+nszc*idm)))  9d24s21
             do in=1,nerimatt                                           9d24s21
              if(isopt(1,nerityp(1,in)).eq.isabcd)then                  9d24s21
               iad1=ieraw(isb,isc,isd,nerityp(1,in))+                   9d24s21
     $             nerityp(2,in)+nerimat(nerityp(1,in))                 9d27s21
     $             *(meraw(isa,isb)+neraw(isa,isb)                      9d24s21
     $             *(iccc+nbasisp(isc)*iddd))                           9d24s21
               bc(iad1)=bc(jtmp+in)                                     9d24s21
               bc(jtmp+in)=0d0                                          10d4s21
              end if                                                    9d24s21
             end do                                                     9d24s21
             meraw(isa,isb)=meraw(isa,isb)+1                            9d24s21
            end do                                                         5d12s10
           end do                                                         10d28s20
          end do                                                        9d24s21
         end do                                                         9d24s21
         sz=0d0                                                         10d4s21
         do is=0,nwtmp*nerimatt-1                                       10d4s21
          sz=sz+bc(itmp+is)**2                                          10d4s21
         end do                                                         10d4s21
         sz=sqrt(sz/dfloat(nwtmp*nerimatt))                             10d4s21
         if(sz.gt.1d-8)then                                             12d19s22
          write(6,*)
     $        ('after tmp to eraw copy, rms size of what''s left is'),sz
          write(6,*)('nwtmp,nerimatt '),nwtmp,nerimatt
          call prntm2(bc(itmp),nwtmp,nerimatt,nwtmp)
          do it=1,4
           write(6,*)('size for type: '),nerimat(it)
          end do
          do in=1,nerimatt
           write(6,*)('for '),in,('type is '),nerityp(1,in),
     $          ('loc is '),nerityp(2,in)
          end do
          stop 'paraeri42cf'
         end if                                                         10d4s21
         ibcoff=itmp                                                    5d12s10
c     jjd
        end do                                                          9d24s21
c     jjc
       end do                                                           9d24s21c
c     under ieraw, we have all of cd in order nerimat,ab,cd.
c     now contract cd into mo basis. for mycode=11 case, we are allready
c     in ntype,ltype order.
c
       do it=1,nsopt                                                    9d24s21
        ntype4=ntype(it)*4                                              9d24s21
        do isb=1,nsymb                                                  9d24s21
         istb=multh(isb,isopt(1,it))                                    9d24s21
         do isd=1,nsymb                                                 9d23s21
          istbd=multh(istb,isd)                                         9d24s21
          do isc=1,nsymb                                                9d23s21
           isa=multh(istbd,isc)                                         9d24s21
           if(min(neraw(isa,isb),nd2(isc),nd1(isd)).gt.0)then           10d4s21
            ncd=nd1(isd)*nd2(isc)                                       10d4s21
            nalls=neraw(isa,isb)*nerimat(it)                            9d24s21
            nrow=nalls*nbasisp(isc)                                     9d24s21
            nrow2=nrow*2                                                9d23s21
            itmp=ibcoff                                                 9d23s21
            ibcoff=itmp+nrow2*nd1(isd)                                  10d4s21
            call enough('paraeri42cf.  5',bc,ibc)
            call dgemm('n','n',nrow,nd1(isd),nbasisp(isd),dgemmf,       10d4s21
     $       bc(ieraw(isb,isc,isd,it)),nrow,bc(iorb(isd)),nbasisp2(isd),9d24s21
     $           0d0,bc(itmp),nrow,                                     9d23s21
     d' paraeri42cf.  1')
            jtmp=itmp+nrow*nd1(isd)                                     10d4s21
            jorb=iorb(isd)+nbasisp(isd)                                 9d23s21
            call dgemm('n','n',nrow,nd1(isd),nbasisp(isd),dgemmf,       10d4s21
     $           bc(ieraw(isb,isc,isd,it)),nrow,bc(jorb),nbasisp2(isd), 9d24s21
     $           0d0,bc(jtmp),nrow,                                     9d23s21
     d' paraeri42cf.  2')
            mrow=2*nalls*nd1(isd)                                       10d4s21
            mrow2=mrow*2                                                9d23s21
            ittmp=ibcoff                                                9d23s21
            ibcoff=ittmp+mrow*nbasisp(isc)                              9d24s21
            call enough('paraeri42cf.  6',bc,ibc)
            do id=0,nd12(isd)-1                                         10d4s21
             do ic=0,nbasisp(isc)-1                                     9d23s21
              icd=itmp+nalls*(ic+nbasisp(isc)*id)                       9d24s21
              idc=ittmp+nalls*(id+nd12(isd)*ic)                         10d4s21
              do ii=0,nalls-1                                           9d24s21
               bc(idc+ii)=bc(icd+ii)                                    9d23s21
              end do                                                    9d23s21
             end do                                                     9d23s21
            end do                                                      9d23s21
            ltmpc=ibcoff                                                9d23s21
            ibcoff=ltmpc+mrow2*nd2(isc)                                 10d4s21
            call enough('paraeri42cf.  7',bc,ibc)
            call dgemm('n','n',mrow,nd2(isc),nbasisp(isc),1d0,          10d4s21
     $           bc(ittmp),mrow,bc(iorb(isc)),nbasisp2(isc),0d0,        9d23s21
     $           bc(ltmpc),mrow,                                        9d23s21
     d' paraeri42cf.  3')
            jtmpc=ltmpc+mrow*nd2(isc)                                   10d4s21
            jorb=iorb(isc)+nbasisp(isc)                                 9d23s21
            call dgemm('n','n',mrow,nd2(isc),nbasisp(isc),1d0,          10d4s21
     $            bc(ittmp),mrow,bc(jorb),nbasisp2(isc),0d0,             9d23s21
     $           bc(jtmpc),mrow,                                        9d23s21
     d' paraeri42cf.  4')
            do ic=0,nd2(isc)-1                                          10d4s21
             do id=0,nd1(isd)-1                                         10d4s21
              do ii=1,mfmx
               if(ifmx(1,ii).eq.it)then                                 9d24s21
                imul=ifmx(10,ii)                                          9d23s21
                xmul=dfloat(imul)                                         9d23s21
                ifrm=nerityp(2,jfmx(1,ii)+1)                            9d27s21
                in=jfmx(2,ii)-1+ntype(it)*(iltype(1,jfmx(3,ii))         9d24s21
     $               +2*iltype(2,jfmx(3,ii)))                           9d23s21
                jdc=id+nd1(isd)*ic                                      10d4s21
                idc=id+nd1(isd)*(iltype(4,jfmx(3,ii))                   10d4s21
     $               +2*(ic+nd2(isc)*iltype(3,jfmx(3,ii))))             10d4s21
c
c     cases:
c     1 - Coulomb SO. then lx,ly,lz can have only ab=ss or ab=ll.
c     2 - Breit symmetric operator. The ab=sl or ab=ls.
c     3 - Breit-Coulomb SO. then lx,ly,lz have ab=ss,sl,ls, or ll.
                do iab=0,neraw(isa,isb)-1                               9d24s21
                 iadc=ltmpc+ifrm+nerimat(it)*(iab+neraw(isa,isb)*idc)
                 iadh=jhalf(isb,isc,isd,it)+jdc+ncd*(in+ntype4*iab)     9d24s21
                 orig=bc(iadh)
                 bc(iadh)=bc(iadh)+xmul*bc(iadc)                        9d24s21
                end do                                                  9d24s21
               end if                                                   9d24s21
              end do                                                    9d24s21
             end do                                                     9d24s21
            end do                                                      9d24s21
            ibcoff=itmp                                                 9d24s21
           end if                                                       9d23s21
          end do                                                        9d23s21
         end do                                                         9d23s21
        end do
       end do                                                           9d22s21
       ibcoff=ieraw(1,1,1,1)                                            9d24s21
       do it=1,nsopt                                                    9d24s21
        ntype4=ntype(it)*4                                              9d24s21
        do isd=1,nsymb                                                   9d23s21
         istd=multh(isopt(1,it),isd)                                    9d24s21
         do isc=1,nsymb                                                  9d23s21
          istdc=multh(istd,isc)                                         9d24s21
          do isb=1,nsymb                                                  9d24s21
           isa=multh(istdc,isb)                                         9d24s21
           ngot=ntype4*neraw(isa,isb)                                   9d24s21
           jhalf(isb,isc,isd,it)=jhalf(isb,isc,isd,it)
     $          +ngot*nd2(isc)*nd1(isd)                                 10d4s21
          end do                                                          9d23s21
         end do                                                           9d23s21
        end do                                                          9d24s21
       end do                                                           9d24s21
      end do                                                            2d19s10
      call second(time3)
c
c     now transpose this across processors so we have all of ab here.
c
      do ipass=1,2                                                      3d6s20
c
c     ipass=1 count, ipass=2 pass about.
c
c     for sending
c
       if(ipass.eq.1)then                                               3d6s20
        do ip=1,mynprocg                                                 3d6s20
         nsnd(ip)=0                                                      3d6s20
        end do                                                           3d6s20
       else                                                             3d6s20
        isnd(1)=ibuffs                                                  3d6s20
        do ip=2,mynprocg                                                3d6s20
         im=ip-1                                                        3d6s20
         isnd(ip)=isnd(im)+nsnd(im)                                     3d6s20
        end do                                                          3d6s20
       end if                                                           3d6s20
       do it=1,nsopt                                                    9d24s21
        ntype4=ntype(it)*4                                              9d24s21
        do isd=1,nsymb                                                   3d6s20
         istd=multh(isd,isopt(1,it))                                    9d24s21
         do isc=1,nsymb                                                  3d6s20
          istdc=multh(istd,isc)                                         9d24s21
          do isb=1,nsymb                                                9d24s21
           if(min(nhalf(isb,isc,isd,it),nd1(isd),nd2(isc)).gt.0)then    10d4s21
            isa=multh(istdc,isb)                                        9d24s21
            nrow=nd1(isd)*nd2(isc)                                      10d4s21
            nhere=nrow/mynprocg                                            3d6s20
            nleft=nrow-nhere*mynprocg                                      3d6s20
            nherep=nhere+1                                                 3d6s20
            irow=0                                                         3d6s20
            do ip=1,mynprocg                                               3d6s20
             if(ip.le.nleft)then                                           3d6s20
              nhereu=nherep                                                3d6s20
             else                                                          3d6s20
              nhereu=nhere                                                 3d6s20
             end if                                                        3d6s20
             if(ipass.eq.1)then                                            3d6s20
              nsnd(ip)=nsnd(ip)+nhereu*nhalf(isb,isc,isd,it)*ntype4     9d24s21
             else                                                          3d6s20
c     ihalf is stored id,ic,ntype4,iab
c     so snd is stored ntype4,iab,id,ic
              do i=1,nhereu                                                3d6s20
               do icol=0,nhalf(isb,isc,isd,it)*ntype4-1                 9d24s21
                iad=ihalf(isb,isc,isd,it)+nrow*icol+irow                9d24s21
                bc(isnd(ip))=bc(iad)                                       3d6s20
                isnd(ip)=isnd(ip)+1                                        3d6s20
               end do                                                      3d6s20
               irow=irow+1                                                 3d6s20
              end do                                                       3d6s20
             end if                                                        3d6s20
            end do                                                         3d6s20
           end if                                                         3d6s20
          end do                                                          3d6s20
         end do                                                           3d6s20
        end do                                                          9d24s21
       end do                                                           9d24s21
       if(ipass.eq.1)then                                               3d6s20
        ntots=nsnd(1)                                                   3d6s20
        do ip=2,mynprocg                                                3d6s20
         ntots=ntots+nsnd(ip)                                           3d6s20
        end do                                                          3d6s20
       else                                                             3d6s20
        isnd(1)=0                                                       3d6s20
        do ip=2,mynprocg                                                3d6s20
         im=ip-1                                                        3d6s20
         isnd(ip)=isnd(im)+nsnd(im)                                     3d6s20
        end do                                                          3d6s20
        call dws_all2allvb(bc(ibuffs),nsnd,isnd,bc(ibuffr),nrcv,ircv)   3d6s20
        ibcoff=ibuffs                                                   3d6s20
        isto1=ibcoff                                                    3d6s20
        ibcoff=isto1+nrtot                                              3d6s20
        ibcsto=ibcoff                                                   10d5s21
        call enough('paraeri42cf.  8',bc,ibc)
        do i=isto1,ibcoff-1
         bc(i)=val
        end do
       end if                                                           3d6s20
       dgemmf=1d0
c
c     for recieving
c
       if(ipass.eq.1)then                                               3d6s20
        do ip=1,mynprocg                                                 3d6s20
         nrcv(ip)=0                                                      3d6s20
        end do                                                           3d6s20
       else                                                             3d6s20
        ircv(1)=ibuffr                                                  3d6s20
        do ip=2,mynprocg                                                3d6s20
         im=ip-1                                                        3d6s20
         ircv(ip)=ircv(im)+nrcv(im)                                     3d6s20
        end do                                                          3d6s20
       end if                                                           3d6s20
       do it=1,nsopt                                                    9d24s21
        ntype4=ntype(it)*4                                              9d24s21
        do isd=1,nsymb                                                   3d6s20
         istd=multh(isopt(1,it),isd)                                    9d24s21
         do isc=1,nsymb                                                  3d6s20
          istdc=multh(istd,isc)                                         9d24s21
          do isb=1,nsymb                                                9d24s21
           if(ipass.eq.2)ihalf(isb,isc,isd,it)=isto1                    9d24s21
           if(min(nd1(isd),nd2(isc)).gt.0)then                          2d13s22
            isa=multh(istdc,isb)                                        9d24s21
            ncol=nd1(isd)*nd2(isc)                                      10d4s21
            nhere=ncol/mynprocg                                            3d6s20
            nleft=ncol-nhere*mynprocg                                      3d6s20
            if(mynowprog.lt.nleft)nhere=nhere+1                            3d6s20
            irow=0
            do ip=1,mynprocg                                               3d6s20
             iddr=isa+nsymb*(isb-1+nsymb*(ip-1))                        9d24s21
             if(ipass.eq.1)then                                            3d6s20
              nrcv(ip)=nrcv(ip)+ntype4*ihrow(iddr)*nhere                9d24s21
             else                                                          3d6s20
c     rcv is stored ntype4,iab,id,ic: so is half.
              do icol=0,nhere-1                                            3d6s20
               jsto1=isto1+ntype4*nhrow(isa,isb)*icol+irow              9d24s21
               do i=0,ntype4*ihrow(iddr)-1                              9d24s21
                bc(jsto1+i)=bc(ircv(ip))                                   3d6s20
                ircv(ip)=ircv(ip)+1                                        3d6s20
               end do                                                      3d6s20
              end do                                                       3d6s20
              irow=irow+ntype4*ihrow(iddr)                              9d24s21
             end if                                                        3d6s20
            end do                                                         3d6s20
            if(ipass.eq.2)isto1=isto1+ntype4*nhrow(isa,isb)*nhere       9d24s21
           end if                                                         3d6s20
          end do                                                          3d6s20
         end do                                                           3d6s20
        end do                                                          9d24s21
       end do                                                           9d24s21
       if(ipass.eq.1)then                                               3d6s20
        nrtot=nrcv(1)                                                   3d6s20
        ircv(1)=0                                                       3d6s20
        do ip=2,mynprocg                                                3d6s20
         nrtot=nrtot+nrcv(ip)                                           3d6s20
         im=ip-1                                                        3d6s20
         ircv(ip)=ircv(im)+nrcv(im)                                     3d6s20
        end do                                                          3d6s20
c
c     we want ibuffr to occupy space that was under ihalf...
c
        ibuffs=ibchalf+max(nwhalf,nrtot)                                10d5s21
        ibcoff=ibuffs+ntots                                             10d5s21
        ibuffr=ibchalf                                                  10d5s21
        call enough('paraeri42cf.  9',bc,ibc)
       end if                                                           3d6s20
      end do                                                            3d6s20
      ibcoff=ibcsto                                                     10d5s21
      nww=0                                                             10d5s21
      do it=1,nsopt                                                     10d5s21
       ntype4=ntype(it)*4                                               10d5s21
       do isd=1,nsymb                                                   10d5s21
        istd=multh(isopt(1,it),isd)                                     10d5s21
        do isc=1,nsymb                                                  10d5s21
         istdc=multh(isc,istd)                                          10d5s21
         do isb=1,nsymb                                                 10d5s21
          isa=multh(isb,istdc)                                          10d5s21
          if(min(nhrow(isa,isb),nd2(isc),nd1(isd)).gt.0)then            10d4s21
           call ilimts(nd1(isd),nd2(isc),                               10d4s21
     $         mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)                9d24s21
           nhere=ih+1-il                                                  9d23s21
           nww=nww+nhere*nbasdws(isa)*nbasdws(isb)*ntype(it)            10d5s21
          end if                                                        10d5s21
         end do                                                         10d5s21
        end do                                                          10d5s21
       end do                                                           10d5s21
      end do                                                            10d5s21
      if(ibuffs-ibchalf.gt.nww)then                                     10d5s21
       err=0d0                                                          12d3s21
      else
       write(6,*)('no, we don''t have enough space :( ')
       err=1d0                                                          12d3s21
      end if
      call dws_gsumf(err,1)
      if(err.ne.0d0)then
       call dws_synca
       call dws_finalize
       stop
      end if
      isto2=ibchalf                                                     10d5s21
c
c     ok, now irows are ab and columns are part of cd, but need to      3d6s20
c     take ab shell order and turn it into normal order.                3d6s20
c
      do it=1,nsopt                                                     9d24s21
       ntype4=ntype(it)*4                                               9d24s21
       do isd=1,nsymb                                                    9d23s21
        istd=multh(isopt(1,it),isd)                                     9d24s21
        do isc=1,nsymb                                                   9d23s21
         istdc=multh(isc,istd)                                          9d24s21
         do isb=1,nsymb                                                 9d24s21
          isa=multh(isb,istdc)                                          9d24s21
          if(min(nhrow(isa,isb),nd2(isc),nd1(isd)).gt.0)then            10d4s21
           call ilimts(nd1(isd),nd2(isc),                               10d4s21
     $         mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)                9d24s21
           nhere=ih+1-il                                                  9d23s21
           if(nhere.gt.0)then                                           10d8s21
            itmp=ibcoff                                                  9d24s21
            ibcoff=itmp+nhere*ntype4*nbasisp(isa)*nbasisp(isb)           9d24s21
            call enough('paraeri42cf. 10',bc,ibc)
            do iz=itmp,ibcoff-1
             bc(iz)=val
            end do
c     half is stored ntype4,iab,id,ic
            irow=0                                                         3d6s20
            do ip=0,mynprocg-1                                             3d6s20
             do i=1+ip,nneed,mynprocg                                      3d6s20
              jpair=ipair+2*(i-1)                                          3d6s20
              ja=jbdat+ibc(jpair)                                          3d6s20
              ja4=ja+ngaus3                                                3d6s20
              ja8=ja+ngaus7                                              9d24s21
              nla=2*ibc(ja)+1                                              3d6s20
              if(iapair(1,ibc(ja8)).gt.0)nla=nla*2                       9d24s21
              jb=jbdat+ibc(jpair+1)                                        3d6s20
              jb4=jb+ngaus3                                                3d6s20
              jb8=jb+ngaus7                                              9d24s21
              nlb=2*ibc(jb)+1                                              3d6s20
              if(iapair(1,ibc(jb8)).gt.0)nlb=nlb*2                       9d24s21
               do ib=0,nlb-1                                               3d6s20
                ibb=ib+ibc(jb4)+1                                        9d24s21
                isbt=isstor(ibb)                                            3d6s20
                ibp=ibstor(ibb)-1                                        9d24s21
                do ia=0,nla-1                                                3d6s20
                 iaa=ia+ibc(ja4)+1                                         3d6s20
                 isat=isstor(iaa)                                             3d6s20
                 iap=ibstor(iaa)-1                                         9d24s21
                 if(isa.eq.isat.and.isb.eq.isbt)then                      9d24s21
                  do icol=0,nhere-1                                       9d24s21
                   do iq=0,ntype4-1                                       9d24s21
                    iadto=itmp+icol+nhere*(iq+ntype4*(ibp                  9d24s21
     $                 +nbasisp(isb)*iap))                              9d24s21
                    iadfrom=ihalf(isb,isc,isd,it)+iq+ntype4*(irow         9d24s21
     $                 +nhrow(isa,isb)*icol)                            9d24s21
                    bc(iadto)=bc(iadfrom)                                  9d24s21
                   end do                                                  9d24s21
                  end do                                                  9d24s21
                  irow=irow+1                                             9d24s21
                 end if                                                   9d24s21
                end do                                                    9d24s21
               end do                                                     9d24s21
              end do                                                      9d24s21
             end do                                                       9d24s21
            nrow=nhere*ntype(it)*nbasisp(isb)                            9d24s21
            nrow4=nrow*4
            nrow2=nrow*2                                                 9d24s21
            ibtmp=ibcoff                                                 9d24s21
            ibcoff=ibtmp+nrow2*nd2(isa)                                  10d4s21
            ltmp=ibcoff                                                  9d24s21
            ibcoff=ltmp+nrow2*nbasisp(isa)                               9d24s21
            call enough('paraeri42cf. 11',bc,ibc)
            do iz=ibtmp,ibcoff-1
             bc(iz)=val
            end do
            do lsb=0,1                                                   9d24s21
             do lsa=0,1                                                  9d24s21
              do ia=0,nbasisp(isa)-1                                     9d24s21
               do ib=0,nbasisp(isb)-1                                    9d24s21
                do iq=0,ntype(it)-1                                      9d24s21
                 iadi=itmp+nhere*(iq+ntype(it)*(lsa+2*(lsb+2*(ib         9d24s21
     $                +nbasisp(isb)*ia))))                              9d24s21
                 iadl=ltmp+nhere*(iq+ntype(it)*(ib+nbasisp(isb)          9d24s21
     $                *(ia+nbasisp(isa)*lsa)))                          9d24s21
                 do i=0,nhere-1                                          9d24s21
                  bc(iadl+i)=bc(iadi+i)                                  9d24s21
                 end do                                                  9d24s21
                end do                                                   9d24s21
               end do                                                    9d24s21
              end do                                                     9d24s21
             end do                                                      9d24s21
             itmpm=ibcoff                                                9d24s21
             ibcoff=itmpm+nrow*nd2(isa)                                  10d4s21
             call enough('paraeri42cf. 12',bc,ibc)
             call dgemm('n','n',nrow,nd2(isa),nbasisp2(isa),dgemmf,      7d6s23
     $           bc(ltmp),nrow,bc(iorb(isa)),nbasisp2(isa),0d0,         9d24s21
     $           bc(itmpm),nrow,                                        9d24s21
     d' paraeri42cf.  5')
             do ia=0,nd2(isa)-1                                          10d4s21
              do ib=0,nbasisp(isb)-1                                     9d24s21
               do iq=0,ntype(it)-1                                       9d24s21
                iadm=itmpm+nhere*(iq+ntype(it)*(ib+nbasisp(isb)*ia))     9d24s21
                iadb=ibtmp+nhere*(iq+ntype(it)*(ia+nd2(isa)*             10d4s21
     $              (ib+nbasisp(isb)*lsb)))                             9d24s21
                do i=0,nhere-1                                           9d24s21
                 bc(iadb+i)=bc(iadm+i)                                   9d24s21
                end do                                                   9d24s21
               end do                                                    9d24s21
              end do                                                     9d24s21
             end do                                                      9d24s21
            end do                                                       9d24s21
            nrow=nhere*ntype(it)*nd2(isa)                                10d4s21
            ihalf(isb,isc,isd,it)=isto2                                  10d6s21
            itmpf=isto2                                                  10d5s21
            isto2=isto2+nrow*nd2(isb)                                    10d5s21
            call enough('paraeri42cf. 13',bc,ibc)
            call dgemm('n','n',nrow,nd2(isb),nbasisp2(isb),1d0,          10d4s21
     $          bc(ibtmp),nrow,bc(iorb(isb)),nbasisp2(isb),0d0,         9d24s21
     $          bc(itmpf),nrow,                                         9d24s21
     d' paraeri42cf.  6')
            ibcoff=itmp                                                 2d14s23
           end if                                                       10d6s21
          end if                                                          9d23s21
         end do                                                           9d23s21
        end do                                                            9d23s21
       end do                                                           9d24s21
      end do                                                            9d24s21
      ibcoff=isto2                                                      10d5s21
c
c     now it's time to fan out to the various procs.
c     there are 2 types of integrals:
c     a) 4o and onex, which are replicated on each proc
c     b) J,K,3x, which are distributed across procs.
c     Let us stow 4o and onex first, finish off by global sum,
c     then do the others.
c
      do it=1,nsopt                                                     10d5s21
       do isd=1,nsymb                                                   10d5s21
        istd=multh(isopt(1,it),isd)                                     10d5s21
        do isc=1,nsymb                                                  10d5s21
         istdc=multh(isc,istd)                                          10d5s21
         do isb=1,nsymb                                                 10d5s21
          isa=multh(isb,istdc)                                          10d5s21
          nrow=irefo(isa)*irefo(isb)                                    10d5s21
          i4o(isb,isc,isd,it)=ibcoff                                    10d5s21
          ionex(isb,isa,isd,it)=i4o(isb,isc,isd,it)                     10d22s21
     $         +nrow*ntype(it)*irefo(isc)*irefo(isd)                    10d5s21
          ibcoff=ionex(isb,isa,isd,it)                                  10d22s21
     $         +nrow*ntype(it)*irefo(isd)*nvirt(isc)                    10d5s21
          if(iprtr(26).ne.0)then                                        3d2s22
           iall(isb,isc,isd,it)=ibcoff                                   10d20s21
           ibcoff=iall(isb,isc,isd,it)+nh0(isa)*nh0(isb)                 10d20s21
     $         *nh0(isc)*nh0(isd)*ntype(it)                             10d20s21
          end if                                                        3d2s22
          call enough('paraeri42cf. 14',bc,ibc)
          do iz=i4o(isb,isc,isd,it),ibcoff-1                            10d5s21
           bc(iz)=0d0                                                   10d5s21
          end do                                                        10d5s21
          call ilimts(nd1(isd),nd2(isc),                                10d5s21
     $         mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)                9d24s21
          nhere=ih+1-il                                                  9d23s21
          if(isopt(2,it).ne.0)then                                      10d5s21
           ssig=-1d0                                                    10d5s21
          else                                                          10d5s21
           ssig=+1d0                                                    10d5s21
          end if                                                        10d5s21
          do iq=0,ntype(it)-1                                           9d24s21
           itest=iifmx(iq+1,it)
           ids=itest/8                                                  9d24s21
           itest=itest-8*ids                                            9d24s21
           ics=itest/4                                                  9d24s21
           itest=itest-4*ics                                            9d24s21
           ibs=itest/2                                                  9d24s21
           ias=itest-2*ibs                                              9d24s21
           do ib=0,nd2(isb)-1                                           10d5s21
            ibr=ib-idoubo(isb)                                          10d1s21
            ibv=ib-idoubo(isb)-irefo(isb)                               10d5s21
            do ia=0,nd2(isa)-1                                          10d5s21
             iar=ia-idoubo(isa)                                         10d4s21
             iav=ia-idoubo(isa)-irefo(isa)                              10d5s21
             i=0                                                        9d24s21
             i10=i1s                                                    9d24s21
             i1n=nd1(isd)                                               10d4s21
             do i2=i2s,i2e                                              9d24s21
              icr=i2-1-idoubo(isc)                                       10d1s21
              icv=icr-irefo(isc)                                        10d8s21
              if(i2.eq.i2e)i1n=i1e                                      9d24s21
              do i1=i10,i1n                                             9d24s21
               idr=i1-1-idoubo(isd)                                      10d1s21
               iadf=ihalf(isb,isc,isd,it)+i+nhere*(iq+ntype(it)*(ia                     10d5s21
     $               +nd2(isa)*ib))                                      10d5s21
               if(min(iar,ibr,icr,idr).ge.0.and.iprtr(26).ne.0)then     3d2s22
                iad2=iall(isb,isc,isd,it)+iar+nh0(isa)*(ibr+nh0(isb)    10d20s21
     $              *(icr+nh0(isc)*(idr+nh0(isd)*iq)))                  10d20s21
                bc(iad2)=bc(iadf)                                       10d20s21
               end if                                                   10d20s21
               if(min(iar,ibr).ge.0.and.iar.lt.irefo(isa).and.          10d5s21
     $               ibr.lt.irefo(isb).and.idr.lt.irefo(isd))then       10d20s21
                if(min(icr,idr).ge.0.and.icr.lt.irefo(isc))then         10d5s21
                 iadm=i4o(isb,isc,isd,it)+iar+irefo(isa)*(ibr
     $               +irefo(isb)*(icr+irefo(isc)*(idr+irefo(isd)*iq)))  10d6s21
                 bc(iadm)=bc(iadf)                                       9d24s21
                else if(icv.ge.0.and.idr.ge.0)then                      10d5s21
                 iadx=ionex(isb,isa,isd,it)+ibr+irefo(isb)*(iar         10d22s21
     $                +irefo(isa)*(idr+irefo(isd)*(icv+nvirt(isc)*iq))) 10d6s21
                 bc(iadx)=bc(iadf)
                end if                                                  10d5s21
               end if                                                    10d5s21
               i=i+1                                                    9d24s21
              end do                                                    9d24s21
              i10=1                                                     9d24s21
             end do                                                     9d24s21
            end do                                                      9d24s21
           end do                                                       9d24s21
          end do                                                        10d5s21
         end do                                                         10d5s21
        end do                                                          10d5s21
       end do                                                           10d5s21
      end do                                                            10d5s21
c
c     fold Breit closed shell part into h0t:                            10d13s21
c
      ibcs=ibcoff                                                       10d13s21
      shift=0d0                                                         10d13s21
      do it=1,nsopt                                                     10d13s21
       ibc0=ibcoff                                                      10d13s21
       do isb=1,nsymb                                                   10d13s21
        jsb=multh(isb,isopt(1,it))                                      10d13s21
        ih0t(isb,it)=ibcoff                                             10d13s21
        ibcoff=ih0t(isb,it)+nh0(jsb)*nh0(isb)                           10d13s21
       end do                                                           10d13s21
       call enough('paraeri42cf. 15',bc,ibc)
       do iz=ibc0,ibcoff-1                                              10d13s21
        bc(iz)=0d0                                                      10d13s21
       end do                                                           10d13s21
       do isd=1,nsymb                                                   10d13s21
        isdt=multh(isd,isopt(1,it))                                     10d13s21
        do isc=1,nsymb                                                  10d13s21
         iscdt=multh(isc,isdt)                                          10d13s21
         call ilimts(nd1(isd),nd2(isc),                                 10d5s21
     $         mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)                 9d24s21
         nhere=ih+1-il                                                  10d13s21
         do isb=1,nsymb                                                 10d13s21
          isa=multh(isb,iscdt)                                          10d13s21
          if(isc.eq.isd.and.min(idoubo(isc),nvirt(isb),nvirt(isa)).gt.0)10d13s21
     $         then                                                     10d13s21
           do iq=1,ntype(it)                                            10d13s21
            ids=iifmx(iq,it)/8                                          10d13s21
            ltest=iifmx(iq,it)-8*ids                                    10d13s21
            ics=ltest/4                                                 10d13s21
            ltest=ltest-4*ics                                           10d13s21
            ibs=ltest/2                                                 10d13s21
            ias=ltest-2*ibs                                             10d13s21
            if(ics.eq.ids.and.ias.eq.0)then                             10d13s21
             iqm=iq-1                                                   10d13s21
             do ib=0,nbasdws(isb)-1                                     10d13s21
              ibv=ib-idoubo(isb)                                        10d13s21
              do ia=0,nbasdws(isa)-1                                    10d13s21
               iav=ia-idoubo(isa)                                       10d13s21
               i=0                                                        9d24s21
               i10=i1s                                                    9d24s21
               i1n=nd1(isd)                                               10d4s21
               do i2=i2s,i2e                                              9d24s21
                if(i2.eq.i2e)i1n=i1e                                      9d24s21
                do i1=i10,i1n                                             9d24s21
                 iadf=ihalf(isb,isc,isd,it)+i+nhere*(iqm+ntype(it)*(ia                     10d5s21
     $               +nd2(isa)*ib))                                      10d5s21
                 if(i1.le.idoubo(isd).and.i2.eq.i1)then                 10d13s21
                  if(min(iav,ibv).ge.0)then                             10d13s21
                   iad=ih0t(isa,it)+iav+nh0(isa)*ibv                    10d13s21
                   bc(iad)=bc(iad)+bc(iadf)                             10d13s21
                  else if(isa.eq.isb.and.ia.eq.ib.and.ib.lt.idoubo(isb))10d13s21
     $                  then                                            10d13s21
                   shift=shift+bc(iadf)                                 10d13s21
                  end if                                                10d13s21
                 end if                                                 10d13s21
                 i=i+1
                end do                                                  10d13s21
                i10=1                                                   10d13s21
               end do                                                   10d13s21
              end do
             end do
            end if                                                      10d13s21
           end do                                                       10d13s21
          end if                                                        10d13s21
          if(isa.eq.isd.and.min(idoubo(isd),nvirt(isb),nvirt(isc)).gt.0)10d13s21
     $         then                                                     10d13s21
            phs=1d0                                                     10d13s21
           phs=-phs                                                     10d13s21
           do iq=1,ntype(it)                                            10d13s21
            ids=iifmx(iq,it)/8                                          10d13s21
            ltest=iifmx(iq,it)-8*ids                                    10d13s21
            ics=ltest/4                                                 10d13s21
            ltest=ltest-4*ics                                           10d13s21
            ibs=ltest/2                                                 10d13s21
            ias=ltest-2*ibs                                             10d13s21
            if(ias.eq.ids.and.ics.eq.0)then                             10d14s21
             iqm=iq-1                                                   10d13s21
             do ib=0,nbasdws(isb)-1                                     10d13s21
              ibv=ib-idoubo(isb)                                        10d13s21
              do ia=0,idoubo(isa)-1                                     10d13s21
               i=0                                                        9d24s21
               i10=i1s                                                    9d24s21
               i1n=nd1(isd)                                               10d4s21
               do i2=i2s,i2e                                              9d24s21
                icv=i2-idoubo(isc)-1                                    10d13s21
                if(i2.eq.i2e)i1n=i1e                                      9d24s21
                do i1=i10,i1n                                             9d24s21
                 i1m=i1-1                                               10d13s21
                 iadf=ihalf(isb,isc,isd,it)+i+nhere*(iqm+ntype(it)*(ia                     10d5s21
     $               +nd2(isa)*ib))                                      10d5s21
                 if(i1.le.idoubo(isd).and.ia.eq.i1m)then                10d13s21
                  if(min(icv,ibv).ge.0)then                             10d13s21
                   iad=ih0t(isc,it)+icv+nh0(isc)*ibv                    10d13s21
                   bc(iad)=bc(iad)+phs*bc(iadf)                         10d13s21
                  else if(isc.eq.isb.and.icv.eq.ibv.and.                10d13s21
     $                  ib.lt.idoubo(isb))then                          10d13s21
                   shift=shift+phs*bc(iadf)                             10d13s21
                  end if                                                10d13s21
                 end if                                                 10d13s21
                 i=i+1
                end do                                                  10d13s21
                i10=1                                                   10d13s21
               end do                                                   10d13s21
              end do
             end do
            end if                                                      10d13s21
           end do                                                       10d13s21
          end if                                                        10d13s21
         end do                                                         10d13s21
        end do                                                          10d13s21
       end do                                                           10d13s21
      end do                                                            10d13s21
      bc(ibcoff)=shift                                                  10d13s21
      nwds=ibcoff-isto2+1                                               10d13s21
      call dws_gsumf(bc(isto2),nwds)                                    10d5s21
      shift=bc(ibcoff)                                                  10d13s21
      do it=1,nsopt                                                     10d13s21
       do isb=1,nsymb                                                   10d13s21
        jsb=multh(isb,isopt(1,it))                                      10d13s21
        if(min(nh0(jsb),nh0(isb),ih0n(isb,it)).gt.0)then                10d13s21
         do ij=0,nh0(jsb)-1                                             10d13s21
          ijp=ij+idoubo(jsb)                                            10d13s21
          do ii=0,nh0(isb)-1                                            10d13s21
           iip=ii+idoubo(isb)                                           10d13s21
           iadn=ih0n(isb,it)+iip+nbasdws(isb)*ijp                       10d13s21
           jtmp=ih0t(isb,it)+ii+nh0(isb)*ij                             7d10s23
           bc(jtmp)=bc(jtmp)+bc(iadn)                                   10d13s21
          end do                                                        10d13s21
         end do                                                         10d13s21
         do ii=0,nh0(jsb)*nh0(isb)-1                                    10d13s21
          bc(ih0n(isb,it)+ii)=bc(ih0t(isb,it)+ii)                               10d13s21
         end do                                                         10d13s21
        end if                                                          10d13s21
       end do                                                           10d13s21
      end do                                                            10d13s21
      ibcoff=ibcs                                                       10d13s21
c
c     current memory map:
c     ibchalf to isto2-1: my part of fully transformed.
c     isto2 to ibcoff-1: 4o and onex.
c
c     now fling about!
c
      do ipass=1,2                                                      10d5s21
c
c     for sending
c
       if(ipass.eq.1)then                                               3d6s20
        do ip=1,mynprocg                                                 3d6s20
         nsnd(ip)=0                                                      3d6s20
        end do                                                           3d6s20
       else                                                             3d6s20
        isnd(1)=ibuffs                                                  3d6s20
        do ip=2,mynprocg                                                3d6s20
         im=ip-1                                                        3d6s20
         isnd(ip)=isnd(im)+nsnd(im)                                     3d6s20
        end do                                                          3d6s20
       end if                                                           3d6s20
       do it=1,nsopt                                                    10d5s21
        do isd=1,nsymb                                                   3d6s20
         istd=multh(isd,isopt(1,it))                                    9d24s21
         do isc=1,nsymb                                                  3d6s20
          call ilimts(nd1(isd),nd2(isc),                                10d5s21
     $         mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)                9d24s21
          nheree=ih+1-il                                                10d8s21
          istdc=multh(istd,isc)                                         9d24s21
          do isb=1,nsymb                                                9d24s21
           isa=multh(istdc,isb)                                         10d5s21
           ncol=nvirt(isb)*nvirt(isc)                                   10d5s21
           nhere=ncol/mynprocg                                          10d5s21
           nleft=ncol-mynprocg*nhere                                    10d5s21
           nhere0=nhere                                                 10d5s21
           if(nleft.gt.0)nhere=nhere+1                                  10d5s21
           ncolac=nvirt(isa)*nvirt(isc)                                 11d15s21
           nhereac=ncolac/mynprocg                                      11d15s21
           nleftac=ncolac-mynprocg*nhereac                              2d14s22
           nhereac0=nhereac                                             11d15s21
           if(nleftac.gt.0)nhereac=nhereac+1                            11d15s21
           ncolab=nvirt(isa)*nvirt(isb)                                   10d5s21
           nhereab=ncolab/mynprocg                                          10d5s21
           nleftab=ncolab-mynprocg*nhereab                                    10d5s21
           nhereab0=nhereab                                             10d5s21
           if(nleftab.gt.0)nhereab=nhereab+1                            10d5s21
           ncoldc=irefo(isd)*nvirt(isc)                                 10d5s21
           nheredc=ncoldc/mynprocg                                      10d5s21
           nleftdc=ncoldc-mynprocg*nheredc                              10d5s21
           nheredc0=nheredc                                             10d5s21
           if(nleftdc.gt.0)nheredc=nheredc+1                            10d5s21
           do ib=0,nd2(isb)-1                                           10d5s21
            ibr=ib-idoubo(isb)                                          10d1s21
            ibv=ib-idoubo(isb)-irefo(isb)                               10d5s21
            do ia=0,nd2(isa)-1                                          10d5s21
             iar=ia-idoubo(isa)                                         10d4s21
             iav=ia-idoubo(isa)-irefo(isa)                              10d5s21
             icol=iav+nvirt(isa)*ibv                                    10d5s21
             if(nhereab.gt.0)then                                       11d19s21
              ipab=icol/nhereab                                          10d5s21
              if(ipab.ge.nleftab)then                                    10d5s21
               ipab=(icol-nleftab*nhereab)/nhereab0                      10d5s21
               ipab=ipab+nleftab                                         10d5s21
              end if
             else                                                       10d28s21
              ipab=0
             end if
             icolab=icol
             ippab=ipab+1                                               10d5s21
             i=0                                                        9d24s21
             i10=i1s                                                    9d24s21
             i1n=nd1(isd)                                               10d4s21
             do i2=i2s,i2e                                              9d24s21
              icr=i2-1-idoubo(isc)                                       10d1s21
              icv=i2-1-idoubo(isc)-irefo(isc)                           10d5s21
              icolac=iav+nvirt(isa)*icv                                 11d15s21
              if(nhereac.gt.0)then                                        11d15s21
               ipac=icolac/nhereac                                      11d15s21
               if(ipac.ge.nleftac)then                                  11d15s21
                ipac=(icolac-nleftac*nhereac)/nhereac0                  11d15s21
                ipac=ipac+nleftac                                       11d15s21
               end if                                                   11d15s21
              else                                                      11d15s21
               ipac=0                                                   11d15s21
              end if                                                    11d15s21
              ipac=ipac+1                                               11d15s21
              icol=ibv+nvirt(isb)*icv                                    10d5s21
              if(nhere.gt.0)then                                        10d28s21
               ip=icol/nhere                                              10d5s21
               if(ip.ge.nleft)then                                       10d6s21
                ip=(icol-nleft*nhere)/nhere0                              10d5s21
                ip=ip+nleft                                               10d5s21
               end if
              else                                                      10d28s21
               ip=0                                                     10d28s21
              end if                                                    10d28s21
              ippx=ip+1                                                 10d7s21
              icolx=icol
              if(i2.eq.i2e)i1n=i1e                                      9d24s21
              do i1=i10,i1n                                             9d24s21
               idr=i1-1-idoubo(isd)                                     10d1s21
               iadf=ihalf(isb,isc,isd,it)+i+nheree*ntype(it)*(ia        10d8s21
     $               +nd2(isa)*ib)                                      10d5s21
c
c     for always occupied spinors ... ???
c
               if(idr.lt.irefo(isd))then                                10d20s21
                if(min(iav,icv).ge.0.and.ibr.lt.irefo(isb).and.ibr.ge.0 12d3s21
     $              .and.idr.ge.0)then                                  12d3s21
                 if(ipass.eq.1)then                                     11d15s21
                  nsnd(ipac)=nsnd(ipac)+ntype(it)                       11d15s21
                 else                                                   11d15s21
                  do iq=0,ntype(it)-1                                   11d15s21
                   bc(isnd(ipac))=bc(iadf+iq*nheree)                    11d15s21
                   isnd(ipac)=isnd(ipac)+1                              11d15s21
                  end do                                                11d15s21
                 end if                                                 11d15s21
                end if                                                  11d15s21
                if(min(iar,ibr).ge.0.and.iar.lt.irefo(isa).and.          10d8s21
     $               ibr.lt.irefo(isb))then                             10d8s21
                else if                                                  10d8s21
     $              (iar.ge.0.and.iar.lt.irefo(isa).and.ibv.ge.0)then   10d8s21
                 if(min(idr,icv).ge.0)then                               10d5s21
                  if(ipass.eq.1)then                                     10d5s21
                   nsnd(ippx)=nsnd(ippx)+ntype(it)                         10d5s21
                  else                                                   10d5s21
                   do iq=0,ntype(it)-1                                   10d5s21
                    bc(isnd(ippx))=bc(iadf+iq*nheree)                   10d20s21
                    isnd(ippx)=isnd(ippx)+1                              10d7s21
                   end do                                                10d5s21
                  end if                                                 10d5s21
                 end if                                                  10d5s21
                else if(min(iav,ibv).ge.0)then                           10d5s21
                 if(min(icr,idr).ge.0.and.icr.lt.irefo(isc))then         10d5s21
                  if(ipass.eq.1)then                                     10d5s21
                   nsnd(ippab)=nsnd(ippab)+ntype(it)                     10d5s21
                  else                                                   10d5s21
                   do iq=0,ntype(it)-1                                   10d5s21
                    bc(isnd(ippab))=bc(iadf+iq*nheree)                  10d20s21
                    isnd(ippab)=isnd(ippab)+1                                10d5s21
                   end do                                                10d5s21
                  end if                                                 10d5s21
                 else if(idr.ge.0.and.icv.ge.0)then                      10d5s21
                  icol=idr+irefo(isd)*icv                                10d5s21
                  if(nheredc.gt.0)then                                  11d19s21
                   ipdc=icol/nheredc                                      10d5s21
                   if(ipdc.ge.nleftdc)then                                10d5s21
                    ipdc=(icol-nleftdc*nheredc)/nheredc0                  10d5s21
                    ipdc=ipdc+nleftdc                                     10d5s21
                   end if                                                 10d5s21
                  else                                                  10d28s21
                   ipdc=0                                               10d28s21
                  end if                                                10d28s21
                  ippdc=ipdc+1                                           10d5s21
                  if(ipass.eq.1)then                                     10d5s21
                   nsnd(ippdc)=nsnd(ippdc)+ntype(it)                     10d5s21
                  else                                                   10d5s21
                   do iq=0,ntype(it)-1                                   10d5s21
                    bc(isnd(ippdc))=bc(iadf+iq*nheree)                  10d20s21
                    isnd(ippdc)=isnd(ippdc)+1                                10d5s21
                   end do                                                10d5s21
                  end if                                                 10d5s21
                 end if                                                  10d5s21
                end if                                                   9d24s21
               end if                                                   10d20s21
               i=i+1                                                    9d24s21
              end do                                                    9d24s21
              i10=1                                                     9d24s21
             end do                                                     9d24s21
            end do                                                      9d24s21
           end do                                                       9d24s21
          end do                                                        10d5s21
         end do                                                         10d5s21
        end do                                                          10d5s21
       end do                                                           10d5s21
c
       if(ipass.eq.1)then                                               10d6s21
        ntots=nsnd(1)                                                   3d6s20
        do ip=2,mynprocg                                                3d6s20
         ntots=ntots+nsnd(ip)                                           3d6s20
        end do                                                          3d6s20
       else                                                             3d6s20
        isnd(1)=0                                                       3d6s20
        do ip=2,mynprocg                                                3d6s20
         im=ip-1                                                        3d6s20
         isnd(ip)=isnd(im)+nsnd(im)                                     3d6s20
        end do                                                          3d6s20
        call dws_all2allvb(bc(ibuffs),nsnd,isnd,bc(ibuffr),nrcv,ircv)   3d6s20
c
c     current memory map:
c     ibchalf=ibuffr to ibuffr+nrtot < isto2:
c     isto2 to ibuffs-1: 4o and onex.
c
        ibcoff=ibuffs                                                   3d6s20
        kq=0
        jq=0
        i3xq=0
        do it=1,nsopt                                                   10d6s21
         do isd=1,nsymb                                                 10d6s21
          istd=multh(isd,isopt(1,it))                                   10d6s21
          do isc=1,nsymb                                                10d6s21
           call ilimts(irefo(isd),nvirt(isc),mynprocg,mynowprog,ildc,   10d6s21
     $          ihdc,j1s,j1e,j2s,j2e)                                   10d6s21
           nheredc=ihdc+1-ildc                                          10d6s21
           istdc=multh(istd,isc)                                         9d24s21
           do isb=1,nsymb                                                9d24s21
            isa=multh(istdc,isb)                                         10d5s21
            nnheredc=nheredc*nvirt(isb)*nvirt(isa)                      10d8s21
            call ilimts(nvirt(isb),nvirt(isc),mynprocg,mynowprog,ilbc,  10d5s21
     $           ihbc,j1s,j1e,j2s,j2e)                                  10d5s21
            nherebc=ihbc+1-ilbc                                         10d5s21
            nnherebc=nherebc*irefo(isa)*irefo(isd)                      10d6s21
            call ilimts(nvirt(isa),nvirt(isb),mynprocg,mynowprog,ilab,  10d6s21
     $           ihab,j1s,j1e,j2s,j2e)                                  10d6s21
            nhereab=ihab+1-ilab                                         10d6s21
            nnhereab=nhereab*irefo(isc)*irefo(isd)                      10d6s21
            call ilimts(nvirt(isa),nvirt(isc),mynprocg,mynowprog,ilac,  11d15s21
     $           ihac,j1s,j1e,j2s,j2e)                                  11d15s21
            nhereac=ihac+1-ilac                                         11d15s21
            nnhereac=nhereac*irefo(isb)*irefo(isd)                      11d15s21
            kmats(isb,isc,isd,it)=ibcoff                                10d6s21
            kmatsb(isb,isc,isd,it)=kmats(isb,isc,isd,it)                 10d6s21
     $           +nnherebc*ntype(it)                                    10d6s21
            jmats(isb,isc,isd,it)=kmatsb(isb,isc,isd,it)                 10d6s21
     $           +nnhereac*ntype(it)                                    11d15s21
            i3x(isb,isc,isd,it)=jmats(isb,isc,isd,it)                   10d6s21
     $           +nnhereab*ntype(it)                                    10d6s21
            ibcoff=i3x(isb,isc,isd,it)+nnheredc*ntype(it)               10d6s21
           end do                                                       10d6s21
          end do                                                        10d6s21
         end do                                                         10d6s21
        end do                                                          10d6s21
        call enough('paraeri42cf. 16',bc,ibc)
        do iz=ibuffs,ibcoff-1
         bc(iz)=val
        end do
       end if                                                           3d6s20
c
c     for receiving
c
       if(ipass.eq.1)then                                               3d6s20
        do ip=1,mynprocg                                                 3d6s20
         nrcv(ip)=0                                                      3d6s20
        end do                                                           3d6s20
       else                                                             3d6s20
        ircv(1)=ibuffr                                                  3d6s20
        do ip=2,mynprocg                                                3d6s20
         im=ip-1                                                        3d6s20
         ircv(ip)=ircv(im)+nrcv(im)                                     3d6s20
        end do                                                          3d6s20
       end if                                                           3d6s20
       do ip=0,mynprocg-1                                               10d5s21
        ippx=ip+1                                                       10d7s21
        do it=1,nsopt                                                    9d24s21
         if(isopt(2,it).ne.0)then                                       10d5s21
          ssig=-1d0                                                     10d5s21
         else                                                           10d5s21
          ssig=+1d0                                                     10d5s21
         end if                                                         10d5s21
         do isd=1,nsymb                                                   3d6s20
          istd=multh(isd,isopt(1,it))                                    9d24s21
          do isc=1,nsymb                                                  3d6s20
           call ilimts(irefo(isd),nvirt(isc),mynprocg,mynowprog,ildc,   10d6s21
     $          ihdc,j1s,j1e,j2s,j2e)                                   10d6s21
           nheredc=ihdc+1-ildc                                          10d6s21
           call ilimts(nd1(isd),nd2(isc),                                10d5s21
     $         mynprocg,ip,il,ih,i1s,i1e,i2s,i2e)                       10d5s21
           nhere=ih+1-il                                                  9d23s21
           istdc=multh(istd,isc)                                         9d24s21
           do isb=1,nsymb                                                9d24s21
            isa=multh(istdc,isb)                                         10d5s21
            nnheredc=nheredc*nvirt(isb)*nvirt(isa)                      12d7s21
            call ilimts(nvirt(isb),nvirt(isc),mynprocg,mynowprog,ilbc,  10d5s21
     $           ihbc,j1s,j1e,j2s,j2e)                                  10d5s21
            nherebc=ihbc+1-ilbc                                         10d5s21
            nnherebc=nherebc*irefo(isa)*irefo(isd)                      10d6s21
            call ilimts(nvirt(isa),nvirt(isb),mynprocg,mynowprog,ilab,  10d6s21
     $           ihab,j1s,j1e,j2s,j2e)                                  10d8s21
            nhereab=ihab+1-ilab                                         10d6s21
            nnhereab=nhereab*irefo(isc)*irefo(isd)                      10d6s21
            call ilimts(nvirt(isa),nvirt(isc),mynprocg,mynowprog,ilac,  11d15s21
     $           ihac,j1s,j1e,j2s,j2e)                                  11d15s21
            nhereac=ihac+1-ilac                                         11d15s21
            nnhereac=nhereac*irefo(isb)*irefo(isd)                      11d15s21
            do ib=0,nd2(isb)-1                                           10d5s21
             ibr=ib-idoubo(isb)                                          10d1s21
             ibv=ib-idoubo(isb)-irefo(isb)                               10d5s21
             do ia=0,nd2(isa)-1                                          10d5s21
              iar=ia-idoubo(isa)                                         10d4s21
              iav=ia-idoubo(isa)-irefo(isa)                              10d5s21
              icolab=iav+nvirt(isa)*ibv+1-ilab                          10d6s21
              i=0                                                        9d24s21
              i10=i1s                                                    9d24s21
              i1n=nd1(isd)                                               10d4s21
              do i2=i2s,i2e                                              9d24s21
               icr=i2-1-idoubo(isc)                                       10d1s21
               icv=i2-1-idoubo(isc)-irefo(isc)                           10d5s21
               icolbc=ibv+nvirt(isb)*icv+1-ilbc                         10d6s21
               icolac=iav+nvirt(isa)*icv+1-ilac                         11d15s21
               if(i2.eq.i2e)i1n=i1e                                      9d24s21
               do i1=i10,i1n                                             9d24s21
                idr=i1-1-idoubo(isd)                                     10d1s21
                icoldc=idr+irefo(isd)*icv+1-ildc                        10d6s21
                if(idr.lt.irefo(isd))then                               10d20s21
                 if(min(iav,icv).ge.0.and.ibr.lt.irefo(isb))then        11d15s21
                  if(icolac.ge.0.and.icolac.lt.nhereac.and.ibr.ge.0.and.12d3s21
     $                idr.ge.0)then                                     12d3s21
                   if(ipass.eq.1)then                                    11d15s21
                    nrcv(ippx)=nrcv(ippx)+ntype(it)                      11d15s21
                   else                                                  11d15s21
                    ito=kmatsb(isb,isc,isd,it)+icolac+nhereac*(          11d15s21
     $                  ibr+irefo(isb)*idr)                             11d15s21
                    do iq=0,ntype(it)-1                                   10d5s21
                     iad=ito+iq*nnhereac                                 11d15s21
                     bc(iad)=bc(ircv(ippx))                              10d7s21
                     ircv(ippx)=ircv(ippx)+1                             11d15s21
                    end do                                               11d15s21
                   end if                                                11d15s21
                  end if                                                12d3s21
                 end if                                                 11d15s21
                 if(min(iar,ibr).ge.0.and.iar.lt.irefo(isa).and.         10d8s21
     $               ibr.lt.irefo(isb))then                             10d8s21
                 else if                                                 10d8s21
     $               (iar.ge.0.and.iar.lt.irefo(isa).and.ibv.ge.0)then  10d8s21
                  if(min(idr,icv).ge.0.and.icolbc.ge.0.and.              10d6s21
     $               icolbc.lt.nherebc)then                             10d6s21
                   if(ipass.eq.1)then                                     10d5s21
                    nrcv(ippx)=nrcv(ippx)+ntype(it)                      10d7s21
                   else                                                   10d5s21
                    ito=kmats(isb,isc,isd,it)+icolbc+nherebc*(           10d6s21
     $                  iar+irefo(isa)*idr)                             10d6s21
                    do iq=0,ntype(it)-1                                   10d5s21
                     iad=ito+iq*nnherebc                                 10d6s21
                     bc(iad)=bc(ircv(ippx))                              10d7s21
                     ircv(ippx)=ircv(ippx)+1                             10d7s21
                    end do                                                10d5s21
                   end if                                                 10d5s21
                  end if                                                  10d5s21
                 else if(min(iav,ibv).ge.0)then                           10d5s21
                  if(min(icr,idr).ge.0.and.icr.lt.irefo(isc).and.         10d6s21
     $                 icolab.ge.0.and.icolab.lt.nhereab)then            10d6s21
                   if(ipass.eq.1)then                                     10d5s21
                    nrcv(ippx)=nrcv(ippx)+ntype(it)                      10d7s21
                   else                                                   10d5s21
                    ito=jmats(isb,isc,isd,it)+icolab+nhereab*(icr        10d6s21
     $                  +irefo(isc)*idr)                                10d6s21
                    do iq=0,ntype(it)-1                                   10d5s21
                     bc(ito+nnhereab*iq)=bc(ircv(ippx))                  10d7s21
                     ircv(ippx)=ircv(ippx)+1                             10d7s21
                    end do                                                10d5s21
                   end if                                                 10d5s21
                  else if(idr.ge.0.and.icv.ge.0.and.icoldc.ge.0.and.     10d6s21
     $                 icoldc.lt.nheredc)then                           10d6s21
                   if(ipass.eq.1)then                                     10d5s21
                    nrcv(ippx)=nrcv(ippx)+ntype(it)                      10d7s21
                   else                                                   10d5s21
                    ito=i3x(isb,isc,isd,it)+ibv+nvirt(isb)*(iav          10d6s21
     $                  +nvirt(isa)*icoldc)                             10d6s21
                    do iq=0,ntype(it)-1                                   10d5s21
                     bc(ito+nnheredc*iq)=bc(ircv(ippx))
                     ircv(ippx)=ircv(ippx)+1                             10d7s21
                    end do                                                10d5s21
                   end if                                                 10d5s21
                  end if                                                  10d5s21
                 end if                                                   9d24s21
                end if                                                  10d20s21
                i=i+1                                                    9d24s21
               end do                                                    9d24s21
               i10=1                                                     9d24s21
              end do                                                     9d24s21
             end do                                                      9d24s21
            end do                                                       9d24s21
           end do                                                        10d5s21
          end do                                                         10d5s21
         end do                                                          10d5s21
        end do                                                           10d5s21
       end do                                                           10d6s21
       if(ipass.eq.1)then                                               3d6s20
        nrtot=nrcv(1)                                                   3d6s20
        ircv(1)=0                                                       3d6s20
        do ip=2,mynprocg                                                3d6s20
         nrtot=nrtot+nrcv(ip)                                           3d6s20
         im=ip-1                                                        3d6s20
         ircv(ip)=ircv(im)+nrcv(im)                                     3d6s20
        end do                                                          3d6s20
c
c     current memory map:
c     ibchalf to isto2-1: my part of fully transformed.
c     isto2 to ibcoff-1: 4o and onex.
c
        ibuffs=ibcoff                                                   10d6s21
        ibcoff=ibuffs+ntots                                             10d6s21
        call enough('paraeri42cf. 17',bc,ibc)
        ibuffr=ibchalf                                                  10d6s21
        itest=ibuffr+nrtot                                              10d6s21
        if(itest.gt.isto2)then                                          10d6s21
         write(6,*)('itest exceeds isto2!! '),itest,isto2               10d6s21
         stop 'paraeri42cf'                                             10d6s21
        end if                                                          10d6s21
       end if                                                           3d6s20
      end do                                                            10d5s21
c
c     current memory map:
c     ibchalf=ibuffr to ibuffr+nrtot < isto2:
c     isto2 to ibuffs-1: 4o and onex etc.
c
      nmov=ibcoff-isto2                                                 10d6s21
      idelta=isto2-ibcoffo                                              10d6s21
      do i=0,nmov-1                                                     10d6s21
       bc(ibcoffo+i)=bc(isto2+i)                                        10d6s21
      end do                                                            10d6s21
      ibcoff=ibcoffo+nmov
      do it=1,nsopt                                                     10d6s21
       do isd=1,nsymb                                                   10d6s21
        istd=multh(isd,isopt(1,it))                                     10d6s21
        do isc=1,nsymb                                                  3d6s20
         ncola=nh0(isc)*nh0(isd)                                        10d20s21
         call ilimts(irefo(isd),nvirt(isc),mynprocg,mynowprog,ildc,     10d6s21
     $          ihdc,j1s,j1e,j2s,j2e)                                   10d6s21
         nheredc=ihdc+1-ildc                                            10d6s21
         istdc=multh(istd,isc)                                          10d6s21
         do isb=1,nsymb                                                 10d6s21
          isa=multh(istdc,isb)                                          10d6s21
          nrowa=nh0(isa)*nh0(isb)                                       10d20s21
          nnheredc=nheredc*nvirt(isb)*nvirt(isa)                         10d6s21
          i4o(isb,isc,isd,it)=i4o(isb,isc,isd,it)-idelta                10d6s21
          ionex(isb,isa,isd,it)=ionex(isb,isa,isd,it)-idelta            10d22s21
          jmats(isb,isc,isd,it)=jmats(isb,isc,isd,it)-idelta            10d6s21
          kmats(isb,isc,isd,it)=kmats(isb,isc,isd,it)-idelta            10d6s21
          kmatsb(isb,isc,isd,it)=kmatsb(isb,isc,isd,it)-idelta          11d15s21
          i3x(isb,isc,isd,it)=i3x(isb,isc,isd,it)-idelta                10d6s21
          if(iprtr(26).ne.0)then                                        3d2s22
           iall(isb,isc,isd,it)=iall(isb,isc,isd,it)-idelta              10d20s21
          end if                                                        3d2s22
          call ilimts(nvirt(isb),nvirt(isc),mynprocg,mynowprog,ilbc,    10d6s21
     $           ihbc,j1s,j1e,j2s,j2e)                                  10d6s21
          nherebc=ihbc+1-ilbc                                           10d6s21
          nnherebc=nherebc*irefo(isa)*irefo(isd)                        10d6s21
          call ilimts(nvirt(isa),nvirt(isb),mynprocg,mynowprog,ilab,    10d6s21
     $           ihab,j1s,j1e,j2s,j2e)                                  10d6s21
          nhereab=ihab+1-ilab                                           10d6s21
          nnhereab=nhereab*irefo(isc)*irefo(isd)                        10d6s21
          call ilimts(nvirt(isa),nvirt(isc),mynprocg,mynowprog,ilac,    11d15s21
     $           ihac,j1s,j1e,j2s,j2e)                                  11d15s21
          nhereac=ihac+1-ilac                                           11d15s21
          nnhereac=nhereac*irefo(isb)*irefo(isd)                        11d15s21
          nrow4=irefo(isa)*irefo(isb)                                   10d6s21
          ncol4=irefo(isc)*irefo(isd)                                   10d6s21
          ncol1=irefo(isd)*nvirt(isc)                                   10d6s21
         end do                                                         10d6s21
        end do                                                          10d6s21
       end do                                                           10d6s21
      end do                                                            10d6s21
      do it=1,nsopt                                                     10d12s21
       do i=1,ntype(it)                                                  10d12s21
        jtest=iifmx(i,it)                                                10d12s21
        ids=jtest/8
        jtest=jtest-8*ids
        ics=jtest/4
        jtest=jtest-4*ics
        ibs=jtest/2
        ias=jtest-2*ibs
       end do
       do i4=0,1                                                         10d8s21
        do i3=0,1                                                        10d8s21
         do i2=0,1                                                       10d8s21
          do i1=0,1                                                      10d8s21
           ltest=i1+2*(i2+2*(i3+2*i4))                                   10d8s21
           itestp=ltest+1                                                10d8s21
           iqptr(itestp)=-1                                              10d8s21
           do i=1,ntype(it)                                              10d12s21
            if(iifmx(i,it).eq.ltest.and.i4.eq.0)then                    11d12s21
             iqptr(itestp)=i-1                                           10d8s21
            end if                                                       10d8s21
           end do                                                        10d8s21
          end do                                                         10d8s21
         end do                                                          10d8s21
        end do                                                           10d8s21
       end do                                                            10d8s21
       nx=-1                                                            12d3s21
       do i=1,32                                                        10d12s21
        iifmx(i,it)=iqptr(i)                                            10d12s21
        nx=max(nx,iqptr(i))                                             12d3s21
       end do                                                           10d12s21
       nx=nx+1                                                          12d3s21
       ntype(it)=nx                                                     12d3s21
      end do                                                            10d12s21
      return
      end                                                               2d19s10
