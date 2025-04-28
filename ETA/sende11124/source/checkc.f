c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine checkc(nfcnt,iptrt,ibasis,nek,mdon,idorbr,isorbr,      1d23s19
     $     iorb1,idop,isorb,nsing,iorb2,iorb3,nnot,nx,nab,              10d18s20
     $     iptrbit,itesta,itestb,norbx,ifcnl,ifcnh,bc,ibc)              11d14s22
      implicit real*8 (a-h,o-z)                                         10d31s22
      integer*1 idorbr(*),isorbr(*),iorb1(*),isorb(*),iorb2(*),iorb3(*),1d23s19
     $     nab(*)                                                       1d23s19
      integer*8 itesta,itestb                                           10d18s20
      dimension iptrt(4,*),ibasis(3,*),iptrbit(2,*),nab4(2,2)           11d13s20
      include "common.store"                                            10d18s20
      data icall/0/                                                     10d19s20
      save icall                                                        10d19s20
      write(6,*)('hi, my name is checkc ...')
      if(bc(132).ne.-132d0)then
       stop 'checkc'
      end if
      do ifcn=ifcnl,ifcnh                                               10d19s20
       icall=icall+1                                                    10d19s20
       nclor=ibasis(1,ifcn)                                              1d23s19
       nopenr=nek-2*nclor                                               1d23s19
       nclorp=nclor+1                                                   1d23s19
       irc=iptrt(2,nclorp)+nclor*(ibasis(2,ifcn)-1)                     1d23s19
       iro=iptrt(4,nclorp)+nopenr*(ibasis(3,ifcn)-1)                    1d23s19
       iirc=iptrbit(1,nclorp)+ibasis(2,ifcn)-1                          10d18s20
       iiro=iptrbit(2,nclorp)+ibasis(3,ifcn)-1                          10d18s20
       call gandc4(ibc(iirc),ibc(iiro),itesta,itestb,nopenr,nsing,      10d18s20
     $      norbx,nnot,nab4,bc,ibc)                                     11d14s22
       if(nnot.ne.0)return                                              1d23s19
      end do                                                            1d23s19
      return                                                            1d23s19
      end                                                               1d23s19
