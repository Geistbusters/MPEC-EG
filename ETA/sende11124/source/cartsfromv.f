c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine cartsfromv(xm,ends,iends,rjac,idax,mya,nvec,nblba,     1d10s19
     $     extradata)                                                   5d26s21
      implicit real*8 (a-h,o-z)
c
c     turn internal coordinates into cartesians
c
      include "common.basis"
      dimension xm(*),iends(idax,2),rjac(idax,3),mya(*),                5d26s21
     $     cart(3,ida),nends(4,2,ida),iradau(ida)
      dimension coordi(ida,ida),ipvt(ida),tmp(3,ida),extradata(*)       1d10s19
      character*(*) ends(idax,2)                                        4d12s19
      write(6,*)(' ')                                                   12d20s19
      iok=0                                                             5d26s21
      natom=nvec+1                                                      4d12s19
      nradau=0                                                          4d12s19
      do i=1,nvec                                                       4d12s19
       do j=1,2                                                         4d12s19
        do k=1,iends(i,j)                                               4d12s19
         itry1=ichar(ends(i,j)(k:k))-ichar('a')                         4d12s19
         itry2=ichar(ends(i,j)(k:k))-ichar('A')                         4d12s19
         if(itry1.ge.0.and.itry1.le.25)then                             4d12s19
          iatom=itry1+1                                                 4d12s19
         else if(itry2.ge.0.and.itry.le.25)then                         4d12s19
          iatom=itry2+1                                                 4d12s19
         else
          write(6,*)('can not decode name: '),ends(i,j)(k:k),itry1,itry2
          stop 'decoder'
         end if                                                         4d12s19
         nends(k,j,i)=iatom                                             4d12s19
        end do                                                          4d12s19
       end do                                                           4d12s19
c
c     radau coordinates are used if one end has 3 or 4 atoms,
c     the other has 1, and it is one of the atoms at the other end.
c
       iradau(i)=0
       if(iends(i,1).gt.iends(i,2))then
        imx=1
        imn=2
       else
        imx=2
        imn=1
       end if
       if(iends(i,imx).ge.3.and.iends(i,imn).eq.1)then
        do j=1,iends(i,imx)
         if(nends(1,imn,i).eq.nends(j,imx,i))iradau(i)=1
        end do
       end if
       nradau=nradau+iradau(i)                                          6d9s94
      end do                                                            4d12s19
      if(nvec.eq.2.and.nblba.eq.2)then                                  9d10s19
       if(nends(1,1,1).eq.nends(1,1,2).or.                              9d10s19
     $    nends(1,1,1).eq.nends(1,2,2))then                             9d10s19
        icommon=nends(1,1,1)                                            9d10s19
       else if(nends(1,2,1).eq.nends(1,1,2).or.                         9d10s19
     $    nends(1,2,1).eq.nends(1,2,2))then                             9d10s19
        icommon=nends(1,2,1)                                            9d10s19
       end if                                                           9d10s19
      end if                                                            9d10s19
      do i=1,natom
       do j=1,natom
        coord(j,i)=0d0
       end do
      end do
      do i=1,natom
       coord(natom,i)=1d0
      end do
      nradaus=nradau                                                    1d2s20
      nradau=0                                                          4d12s19
      do iv=1,nvec
       if(iradau(iv).eq.0)then
        if(nvec.eq.2.and.nblba.eq.2)then                                9d10s19
         write(6,*)('vector '),iv,(' is a blba vector')                 9d10s19
         ft=0d0                                                         3d31s21
         fh=1d0                                                         3d31s21
         do i=1,iends(iv,1)                                             9d10s19
          coord(iv,nends(i,1,iv))=ft                                    9d10s19
         end do                                                         9d10s19
         do i=1,iends(iv,2)                                             9d10s19
          coord(iv,nends(i,2,iv))=fh                                    9d10s19
         end do                                                         9d10s19
        else                                                            9d10s19
         write(6,*)('vector '),iv,(' is a jacobi ')
         if(iv.eq.1)iok=1                                               5d26s21
c
c     for jacobi vectors
c
         xmh=0d0
         do i=1,iends(iv,2)
          xmh=xmh+xm(nends(i,2,iv))
         end do
         xmt=0d0
         do i=1,iends(iv,1)
          xmt=xmt+xm(nends(i,1,iv))
         end do
         fh=xmt/(xmt+xmh)
         ft=-xmh/(xmt+xmh)
         do i=1,iends(iv,2)
          coord(iv,nends(i,2,iv))=fh
         end do
         do i=1,iends(iv,1)
          coord(iv,nends(i,1,iv))=ft
         end do
        end if                                                          9d10s19
       else
c
c     for radau vectors
c
        nradau=nradau+1                                                 4d12s19
        write(6,*)('radau '),nradau
        if(iends(iv,1).gt.iends(iv,2))then
         imx=1
         imn=2
        else
         imx=2
         imn=1
        end if
         xmhvy=0d0                                                       3d25s94
         nlight=iends(iv,imx)-iends(iv,imn)                             9d6s21
         do 997 i=nlight+1,iends(iv,imx)                                9d7s21
          xmhvy=xmhvy+xm(nends(i,imx,iv))                                3d25s94
  997    continue                                                        3d25s94
         xmall=0d0                                                       3d25s94
         do 998 i=1,iends(iv,imx)                                       4d12s19
          xmall=xmall+xm(nends(i,imx,iv))                                3d25s94
  998    continue                                                        3d25s94
         xmo=xmall-xmhvy                                                 3d25s94
         gamma=(1d0-sqrt(xmall/xmhvy))/xmo
         beta=-(xmhvy*gamma+1d0)/xmall
         coord(iv,nends(1,imn,iv))=1d0
         bb=beta*xm(nends(1,imn,iv))
         do 68 i=iends(iv,imn),iends(iv,imx)                            2d26s21
          coord(iv,nends(i,imx,iv))=coord(iv,nends(i,imx,iv))+bb
   68    continue
         do 999 i=nlight+1,iends(iv,imx)                                9d6s21
          coord(iv,nends(i,imx,iv))=(gamma+beta)*xm(nends(1,imn,iv))     3d25s94
  999    continue                                                        3d25s94
         if(imx.eq.2)then                                                3d25s94
          do 996 i=1,iends(iv,imx)                                      4d12s19
           coord(iv,nends(i,imx,iv))=-coord(iv,nends(i,imx,iv))          3d25s94
  996     continue                                                       3d25s94
         end if                                                          3d25s94
        end if
      end do
      write(6,*)('we have coord matrix: ')
      call prntm2(coord,natom,natom,ida)                                4d12s19
c
c     carts of internal vectors in standard orientation
c
      cart(1,1)=0d0
      cart(2,1)=0d0
      cart(3,1)=rjac(1,1)                                               4d12s19
      if(nvec.gt.1)then                                                 4d12s19
       ca=cos(rjac(1,2))                                                4d12s19
       sa=sin(rjac(1,2))                                                4d12s19
       cart(1,2)=rjac(2,1)*sa                                           4d12s19
       cart(2,2)=0d0                                                    4d12s19
       cart(3,2)=rjac(2,1)*ca                                           4d12s19
       do iv=3,nvec                                                     4d12s19
        ivm=iv-1                                                        4d12s19
        ivmm=ivm-1                                                      4d12s19
        ca=cos(rjac(ivm,2))                                              4d12s19
        sa=sin(rjac(ivm,2))                                             4d12s19
        sp=sin(rjac(ivmm,3))                                            4d12s19
        cp=cos(rjac(ivmm,3))                                            4d12s19
        cart(1,iv)=rjac(iv,1)*sa*cp                                     4d12s19
        cart(2,iv)=rjac(iv,1)*sa*sp                                     4d12s19
        cart(3,iv)=rjac(iv,1)*ca                                        4d12s19
       end do                                                           4d12s19
      end if                                                            4d12s19
      write(6,*)('atomic positions computed from internal coordinates:')12d20s19
      write(6,*)('>rs: '),(rjac(iv,1),iv=1,nvec)                         12d20s19
      do iv=1,nvec                                                      1d10s19
       extradata(iv)=rjac(iv,1)                                         1d10s19
      end do                                                            1d10s19
      if(nvec.gt.1)then                                                 12d20s19
       write(6,*)('>thetas: '),(rjac(iv,2),iv=1,nvec-1)                  12d20s19
       do iv=1,nvec-1                                                   1d10s19
        extradata(iv+nvec)=rjac(iv,2)                                   1d10s19
       end do                                                           1d10s19
       if(nvec.gt.2)then                                                12d20s19
        write(6,*)('>phis: '),(rjac(iv,3),iv=1,nvec-2)                   12d20s19
        noff=nvec*2-1                                                   1d10s19
        do iv=1,nvec-2                                                  1d10s19
         extradata(iv+noff)=rjac(iv,3)                                  1d10s19
        end do                                                          1d10s19
       end if                                                           12d20s19
      end if                                                            12d20s19
      write(6,*)(' ')
      do i=1,natom                                                      3d21s23
       do j=1,natom                                                     3d21s23
        tmp(1,j)=coord(i,mya(j))                                        3d21s23
       end do
       do j=1,natom
        coord(i,j)=tmp(1,j)
       end do
      end do
      call prntm2(coord,natom,natom,ida)
      call dgemm('n','n',3,natom,nvec,1d0,cart,3,coord,ida,0d0,tmp,3,   4d12s19
     d' cartsfromv.  1')
      do i=1,natom
       do ixyz=1,3                                                      4d12s19
        xcart(ixyz,i)=tmp(ixyz,i)                                       3d21s23
       end do                                                           4d12s19
      end do
      if(iok.ne.0)then                                                  5d26s21
       write(6,*)('vector 1 is a jacobi vector')                        5d26s21
       ntarg=iends(1,1)                                                 5d26s21
       nproj=iends(1,2)                                                 5d26s21
       ntargp=ntarg+nproj                                               5d26s21
       if(ntargp.eq.natom)then                                          5d26s21
        write(6,*)('we have target-projectile situation ')              4d27s23
        do i=1,2                                                        5d26s21
         write(6,*)('for group '),i
         xmh=0d0                                                        5d26s21
         do j=1,iends(1,i)                                              5d26s21
          iagrp(nends(j,i,1))=i                                         5d26s21
          xmh=xmh+xm(nends(j,i,1))                                      5d26s21
         end do                                                         5d26s21
         write(6,*)('mass of group = '),xmh
         xmhi=1d0/xmh                                                   5d26s21
         do ixyz=1,3                                                    5d26s21
          cmx(ixyz,i)=0d0                                                  5d26s21
          do j=1,iends(1,i)                                             5d26s21
           cmx(ixyz,i)=cmx(ixyz,i)+xm(nends(j,i,1))*xmhi*xcart(ixyz,    5d26s21
     $          nends(j,i,1))                                           5d26s21
          end do                                                        5d26s21
         end do                                                         5d26s21
        end do                                                          5d26s21
        write(6,*)('centers of mass of groups: ')                       2d4s23
        call prntm2(cmx,3,2,3)                                          2d4s23
       end if                                                           5d26s21
      end if                                                            5d26s21
      return                                                            4d12s19
      end
