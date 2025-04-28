c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dovr(iwfb,iwfk,ovr,mdon,mdoop,ixw1,ixw2,ism,irel,irefo,7d27s21
     $      norb,nbasdws,idoubo,nvirt,maxbx,maxbxd,srh,sr2,multh,nsymb, 7d27s21
     $      ncsf,npadddi,bc,ibc)                                        11d10s22
      implicit real*8 (a-h,o-z)                                         2d11s22
      dimension iwfb(*),iwfk(*),ovr(*),ism(*),irel(*),irefo(*),         2d11s22
     $     nbasdws(*),idoubo(*),nvirt(*),multh(8,8),ncsf(*),ixmt(8)     2d11s22
      include "common.store"                                            2d11s22
      ibcoffo=ibcoff                                                    2d11s22
      do isb=1,nsymb                                                    2d11s22
       ixmt(isb)=ibcoff                                                 2d11s22
       ibcoff=ixmt(isb)+nbasdws(isb)*nbasdws(isb)                       2d11s22
       call enough('dovr.  1',bc,ibc)
       do i=0,nbasdws(isb)-1                                            2d11s22
        ix=ixmt(isb)+nbasdws(isb)*i                                     2d11s22
        do j=0,nbasdws(isb)-1                                           2d11s22
         bc(ix+j)=0d0                                                   2d11s22
        end do                                                          2d11s22
        if(i.ge.idoubo(isb))then                                        2d11s22
         bc(ix+i)=1d0                                                    2d11s22
        end if                                                          2d11s22
       end do                                                           2d11s22
      end do                                                            2d11s22
      nroot=iwfb(3)                                                     2d11s22
      call psioppsi(iwfb,nroot,1,1d0,ixmt,1,0,idum,dum,iwfk,nroot,      2d11s22
     $     ovr,mdon,mdoop,ixw1,ixw2,ism,irel,irefo,norb,nbasdws,        2d11s22
     $     idoubo,nvirt,maxbx,maxbxd,srh,sr2,multh,nsymb,iloob,ncsf,    2d11s22
     $     0,npadddi,bc,ibc)                                            11d10s22
      ibcoff=ibcoffo                                                    2d11s22
      return                                                            2d11s22
      end                                                               2d11s22
