c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine restin(la,zza,xna,xa,ya,za,                            3d13s20
     $               lb,zzb,xnb,xb,yb,zb,                               3d13s20
     $               lc,zzc,xnc,xc,yc,zc,                               3d13s20
     $               ld,zzd,xnd,xd,yd,zd,                               3d13s20
     $               icartt,idoit,fmulq,idorel,ascale,nmat,bc,ibc)      11d9s22
      implicit real*8 (a-h,o-z)
      logical ldeb,lab,lcd                                              10d8s15
c
c     take electron repulsion integrals and store them properly
c     for 2 component calculation, and optionally compute additional
c     2 electron integrals.
c     return nszd*nszc*nszb*nsza block of LLLL ints,
c     and if idorel ne 1,
c     next the nszd*nszc*nszb*nsza block of LLSS ints,
c     next the nszd*nszc*nszb*nsza block of SSLL ints,
c     and if idorel indicates so,
c     next the nszd*nszc*nszb*nsza block of SSSS ints.
c
      include "common.store"
      include "common.spher"
      include "common.rys"                                              7d5s12
      dimension fmulx(4)                                                3d13s20
      common/timerocm/tovr,telapo(15)                                    5d7s12
      data icall/0/                                                     7d22s14
      save                                                              7d22s14
      ldeb=.false.                                                      10d8s15
      ibcoffo=ibcoff                                                    8d26s15
      nsza=2*la+1
      nszac=nsza*2
      nszb=2*lb+1
      nszbc=nszb*2
      nszc=2*lc+1
      nszcc=nszc*2
      nszd=2*ld+1
      nszdc=nszd*2
      ntemp=nsza*nszb*nszc*nszd                                         10d8s15
      if(idorel.eq.1)then                                               3d13s20
       call erir(la,zza,xna,xa,ya,za,                                   3d13s20
     $            lb,zzb,xnb,xb,yb,zb,                                  3d13s20
     $            lc,zzc,xnc,xc,yc,zc,                                  3d13s20
     $            ld,zzd,xnd,xd,yd,zd,icartt,idoit,1d0,bc,ibc)          11d14s22
       nmat=1                                                           3d13s20
      else                                                              3d13s20
       fmulx(1)=1d0                                                      3d13s20
       fmulx(2)=ascale                                                   3d13s20
       fmulx(3)=ascale                                                   3d13s20
       fmulx(4)=ascale*ascale                                            3d13s20
       if(idorel.eq.2.or.idorel.eq.-2.or.idorel.eq.-4)then              10d24s20
        icode=5
        nmat=2                                                          3d18s20
       else
        icode=4                                                         3d13s20
        nmat=4                                                          3d13s20
       end if                                                           3d13s20
       call dgerid(la,zza,xna,xa,ya,za,                                  3d13s20
     $            lb,zzb,xnb,xb,yb,zb,                                  3d13s20
     $            lc,zzc,xnc,xc,yc,zc,                                  3d13s20
     $            ld,zzd,xnd,xd,yd,zd,                                  3d13s20
     $         0,0,0, 0,0,0, 0,0,0, 0,0,0,icartt,.false.,icode,fmulx,   11d14s22
     $      bc,ibc)                                                     11d14s22
      end if                                                            3d13s20
      ibcoff=icartt+ntemp*nmat                                          3d13s20
      if(icode.eq.5)then                                                10d28s20
       call dgerid(la,zza,xna,xa,ya,za,                                  3d13s20
     $            lb,zzb,xnb,xb,yb,zb,                                  3d13s20
     $            lc,zzc,xnc,xc,yc,zc,                                  3d13s20
     $            ld,zzd,xnd,xd,yd,zd,                                  3d13s20
     $         0,0,0, 0,0,0, 0,0,0, 0,0,0,icartt2,.false.,6,fmulx(2),   11d14s22
     $     bc,ibc)                                                      11d14s22
       ibcoff=icartt2+ntemp                                              10d28s20
      end if                                                            10d28s20
      return
      end
