c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine sldvrc(lb,b1,b2,b3,b4,b5,lk,k1,k2,k3,k4,k5,iout,idx,   3d25s20
     $     idy,idz,bc,ibc)                                              11d9s22
      implicit real*8 (a-h,o-z)
      real*8 k1,k2,k3,k4,k5
c
c     get sigma p contrabution to small-large coupling.                 3d4s20
c
      dimension ider(3,3),nspin(2,2),ispin(3,2,2,2),phs(2),             3d4s20
     $     ioffb(2),ioffk(2),cvec(6,2,2),pref(2)                        3d25s20
      data ider/1,0,0, 0,1,0, 0,0,1/                                    3d4s20
      data phs/1d0,-1d0/                                                3d4s20
      data nspin/1,2,2,1/
      include "common.store"                                            3d4s20
      ibcoffo=ibcoff                                                    3d4s20
      nlb=2*lb+1                                                        3d4s20
      nlb2=nlb*2                                                        3d4s20
      nlk=2*lk+1                                                        3d4s20
      nlk2=nlk*2                                                        3d4s20
      iout=ibcoff                                                       3d4s20
      iouti=iout+nlb2*nlk2                                              3d4s20
      ibcoff=iouti+nlb2*nlk2                                            3d4s20
      iprim=ibcoff                                                      3d25s20
      ibcoff=iprim+nlb2*nlk2*9                                          3d25s20
      call enough('sldvrc.  1',bc,ibc)
      do i=iout,ibcoff-1                                                3d4s20
       bc(i)=0d0                                                        3d4s20
      end do                                                            3d4s20
      fact=0.5d0
c
c     contraction from real 4 component basis functions to complex
c     2 component basis functions
c
      do i3=1,2                                                         1d13s23
       do i2=1,2                                                        1d13s23
        do i1=1,6                                                       1d13s23
         cvec(i1,i2,i3)=0d0                                             1d13s23
        end do                                                          1d13s23
       end do                                                           1d13s23
      end do                                                            1d13s23
c     what dyall used
      cvec(5,1,1)=fact
      cvec(3,1,2)=-fact
      cvec(4,1,2)=-fact
      cvec(2,2,1)=-fact
      cvec(1,2,2)=-fact
      cvec(6,2,2)=fact
      nlb6=nlb2*3                                                       3d25s20
      do ixyz=1,3                                                       3d25s20
       iof=nlb2*(ixyz-1)                                                3d25s20
       do ispinlk=0,1                                                   3d25s20
        if(idz.ne.0)then                                                3d25s20
         ispinsb=ispinlk                                                3d25s20
        else                                                            3d25s20
         ispinsb=1-ispinlk                                              3d25s20
        end if                                                          3d25s20
        sigma=1d0                                                       3d25s20
        if(idz.ne.0.and.ispinlk.eq.1)sigma=-sigma                       3d25s20
        if(idy.ne.0.and.ispinlk.eq.1)sigma=-sigma                       3d25s20
        jprim=iprim+nlb*(ispinsb+2*(nlk*(ispinlk+2*(ixyz-1))))          3d28s20
        call onep(lb,b1,b2,b3,b4,b5,lk,k1,k2,k3,k4,k5,idum,ibo,         3d25s20
     $        ider(1,ixyz),ider(2,ixyz),ider(3,ixyz),0,0,0,             3d25s20
     $        idx,idy,idz,bc,ibc)                                       11d9s22
c
c     recall the transpose comes out of onep
c
        do ib=0,nlb-1                                                   3d25s20
         do ik=0,nlk-1                                                   3d25s20
          iad=jprim+nlb2*ik                                              3d28s20
          bc(iad+ib)=bc(ibo+ik)*sigma                                   3d25s20
         end do                                                         3d25s20
         ibo=ibo+nlk                                                    3d31s20
        end do                                                          3d25s20
       end do                                                           3d25s20
      end do                                                            3d25s20
      pref(1)=0d0                                                       3d25s20
      pref(2)=0d0                                                       3d25s20
      if(idy.ne.0)then                                                  3d25s20
       pref(1)=1d0                                                      3d25s20
      else                                                              3d25s20
       pref(2)=-1d0                                                     3d25s20
      end if                                                            3d25s20
      do ispinlb=0,1                                                    3d25s20
       jspinlb=ispinlb+1                                                3d25s20
       iaddb=ispinlb*nlb                                                3d25s20
       do ibx=0,2                                                       3d25s20
        do ispinsb=0,1                                                  3d25s20
         kb=1+ibx+3*ispinsb                                             3d28s20
         ffr=cvec(kb,jspinlb,1)*pref(1)                                 3d25s20
     $      +cvec(kb,jspinlb,2)*pref(2)                                 3d25s20
         ffi=cvec(kb,jspinlb,1)*pref(2)                                 3d25s20
     $      -cvec(kb,jspinlb,2)*pref(1)                                 3d25s20
         if(abs(ffr).gt.1d-10)then                                      3d25s20
          do ispinlk=0,1                                                    3d25s20
           iaddk=ispinlk*nlk                                                3d25s20
           do ik=0,nlk-1                                                 3d25s20
            do ib=0,nlb-1                                                3d25s20
c
c     we are expecting transpose on output
c
             ito=iout+ik+iaddk+nlk2*(ib+iaddb)                          3d25s20
             ifrom=iprim+ib+nlb*(ispinsb+2*(ik+nlk*(ispinlk+2*ibx)))
             bc(ito)=bc(ito)+bc(ifrom)*ffr                              3d25s20
            end do
           end do
          end do
         end if                                                         3d25s20
         if(abs(ffi).gt.1d-10)then                                      3d25s20
          do ispinlk=0,1                                                    3d25s20
           iaddk=ispinlk*nlk                                                3d25s20
           do ik=0,nlk-1                                                 3d25s20
            do ib=0,nlb-1                                                3d25s20
             ito=iouti+ik+iaddk+nlk2*(ib+iaddb)                          3d25s20
             ifrom=iprim+ib+nlb*(ispinsb+2*(ik+nlk*(ispinlk+2*ibx)))
             bc(ito)=bc(ito)+bc(ifrom)*ffi                              3d25s20
            end do
           end do
          end do
         end if                                                         3d25s20
        end do                                                          3d25s20
       end do                                                           3d25s20
      end do                                                            3d25s20
      ibcoff=iout                                                       3d4s20
      return                                                            3d4s20
      end                                                               3d4s20
