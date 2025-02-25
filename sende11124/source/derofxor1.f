c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine derofxor1(ovr,iso,ivecs,ipropsym,idorel,ieigs,         5d4s22
     $     ixor0,multh,idwsdeb,iorb,ivecso,nbasisp,isob,bc,ibc)         11d14s22
c
c     differentiate ao to orthogonal basis transformation.
c
      implicit real*8 (a-h,o-z)
      logical lbodc                                                     6d3s22
      include "common.hf"
      include "common.store"
      include "common.spher"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension ovr(1),iso(1),ivecs(1),ieigs(1),multh(8,8),isob(*),     5d4s22
     $     iorb(1),propmat(1),iptb(8),ivecso(*),nbasisp(*)              6d7s22
      if(idwsdeb.ne.0)write(6,*)('Hi, my name is derofxor1! ')
      ibcoffo=ibcoff
      ncomp=1
      if(idorel.ne.0)ncomp=2
      if(ipropsym.eq.1)then
       ioff=1
       ioffb=1                                                          4d15s22
       do isb=1,nsymb
        if(idwsdeb.ne.0)write(6,*)('for symmetry block '),isb
        nh=nbasdws(isb)                                                 4d4s22
        nhp=nbasisp(isb)*ncomp                                          4d4s22
        ixor=ibcoff
        if(isb.eq.1)ixor0=ixor
        ibcoff=ixor+nhp*nhp                                             5d4s22
        ibcoffo=ibcoff
        itmp1=ibcoff
        itmp2=itmp1+nhp*nhp                                             4d4s22
        itmp3=itmp2+nhp*nhp                                             4d4s22
        ibcoff=itmp3+nhp*nhp                                            4d4s22
        call enough('derofxor1.  1',bc,ibc)
        if(nhp.gt.0)then                                                5d20s22
         if(idwsdeb.ne.0)write(6,*)('calling square .... '),ioff,nhp
         call square(ovr(ioff),nhp)                                      4d4s22
         iad=ivecso(isb)+nh*nh                                           4d6s22
         if(idwsdeb.gt.10)then
          write(6,*)('der of overlap ')
          call prntm2(ovr(ioff),nhp,nhp,nhp)                              4d4s22
          call prntm2(bc(iad),nhp,nh,nhp)                                 4d4s22
         end if
         call dgemm('n','n',nhp,nh,nhp,1d0,ovr(ioff),nhp,bc(iad),nhp,
     $        0d0,bc(itmp1),nhp,                                        5d20s22
     d' derofxor.  5')
         if(idwsdeb.gt.10)then
          write(6,*)('transformed from the right ')
          call prntm2(bc(itmp1),nhp,nh,nhp)                               4d4s22
         end if
         do i=0,nhp-1                                                    4d4s22
          do j=0,nh-1                                                    4d4s22
           ij=j+nh*i
           ji=i+nhp*j                                                    4d4s22
           bc(itmp2+ij)=bc(itmp1+ji)
          end do
         end do
         call dgemm('n','n',nh,nh,nhp,1d0,bc(itmp2),nh,bc(iad),nhp,0d0,  4d4s22
     $       bc(itmp1),nh,
     d' derofxor.  6')
         if(idwsdeb.gt.10)then
          write(6,*)('dS in diagonal basis: '),sum
          call prntm2(bc(itmp1),nh,nh,nh)
         end if
         ids=ibcoff
         ibcoff=ids+nh
         ise=ibcoff
         ibcoff=ise+nh
         call enough('derofxor1.  2',bc,ibc)
         n0=0
         do i=0,nh-1
          iad=itmp1+i*(nh+1)
          if(bc(ieigs(isb)+i).ne.0d0)then
           bc(ise+i)=1d0/(bc(ieigs(isb)+i)**2)
           bc(ids+i)=bc(iad)
          else
           bc(ise+i)=0d0
           bc(ids+i)=0d0
           n0=n0+1
          end if
         end do
         if(idwsdeb.gt.10)then
          write(6,*)('number of vectors set to zero = '),n0
          write(6,*)('eigenvalues of overlap matrix: ')
          call prntm2(bc(ise),1,nh,1)
          write(6,*)('derivatives of eigenvalues: ')
          call prntm2(bc(ids),1,nh,1)
         end if
         do i=0,n0-1
          do j=0,n0-1
           iad=itmp2+j+nh*i
           bc(iad)=0d0
          end do
         end do
         do i=n0,nh-1
          iad=itmp2+i*(nh+1)
          bc(iad)=0d0
          do j=i+1,nh-1
           ji=itmp2+i+nh*j
           ji2=ji+itmp1-itmp2
           ij=itmp2+j+nh*i
           bot=bc(ise+i)-bc(ise+j)
           if(abs(bot).lt.1d-14)then
            if(abs(bc(ji2)).gt.1d-10)then
             write(6,*)('small denominator ...'),i,j,bot,
     $            ('but not so small numberator! '),bc(ji2),
     $            bc(ise+i),bc(ji2)
            end if
            bot=1d0
           end if
           bc(ji)=-bc(ji2)/bot                                           1d21s16
           bc(ij)=-bc(ji)
          end do
         end do
         if(idwsdeb.gt.10)then
          write(6,*)('Vt dV ')
          call prntm2(bc(itmp2),nh,nh,nh)
         end if
         iad=ivecso(isb)+nh*nh                                           4d6s22
         iad=ivecso(isb)+nh*n0                                           4d6s22
         call dgemm('n','n',nh,nh,nh-n0,1d0,bc(iad),nh,bc(itmp2+n0),nh,  4d6s22
     $       0d0,bc(itmp1),nh,                                          4d6s22
     d' derofxor. 11')
         if(idwsdeb.gt.10)then
          write(6,*)('dV ')
          call prntm2(bc(itmp1),nh,nh,nh)                                 4d6s22
         end if
         call enough('derofxor1.  3',bc,ibc)
         if(idwsdeb.gt.10)then
          write(6,*)('symmetric orthogonalization: ')
          write(6,*)('tmp1: = dV')
          call prntm2(bc(itmp1),nh,nh,nh)
          write(6,*)('ivecso')
          call prntm2(bc(ivecso(isb)),nh,nh,nh)
         end if
         itmp4=ibcoff                                                    6d22s16
         ibcoff=itmp4+nhp*nh                                             4d4s22
         itmp5=ibcoff
         ibcoff=itmp5+nhp*nh                                             4d4s22
         call enough('derofxor1.  4',bc,ibc)
         do i=0,nh-1
          smh=bc(ieigs(isb)+i)
          part1=-0.5d0*smh*smh*smh
          dsmh=part1*bc(ids+i)
          bc(ise+i)=dsmh
          do j=0,nh-1                                                    4d6s22
           ij=i+j*nh
           ji=j+i*nh                                                     4d6s22
           iad1=itmp2+ij
           iad2=itmp1+ji
           iad3=ivecso(isb)+ji                                           4d4s22
           iad4=itmp3+ij
           iad5=itmp4+ij
c     s^{-1/2}VsT: dVs times goes to dX, d2Vs times goes to d2X
           bc(iad4)=smh*bc(iad3)
c     s^{-1/2}dVsT-.5 ds s^{-3/2} VsT, VS times goes to dX, 2dVs times goes to d2X
           bc(iad1)=smh*bc(iad2)+dsmh*bc(iad3)
          end do
         end do
         iad=ivecso(isb)                                                 4d1s22
         call dgemm('n','n',nh,nh,nh,1d0,bc(iad),nh,bc(itmp2),nh,0d0,    4d9s22
     $        bc(ixor),nh,                                              4d6s22
     d' derofxor. 15')
         if(idwsdeb.gt.10)then
          write(6,*)('xor of 2nd bit = Vsmh dV')
          call prntm2(bc(ixor),nh,nh,nh)
          if(nsymb.eq.1)call printa(bc(ixor),nh,0,1,0,nh,0,1,0,
     $         bc(ibcoff))
         end if
         call dgemm('n','n',nh,nh,nh,1d0,bc(itmp1),nh,bc(itmp3),nh,     4d6s22
     $        1d0,bc(ixor),nh,                                          4d6s22
     d' derofxor. 19')
         if(idwsdeb.gt.10)then
          write(6,*)('derivative of orthogonalization matrix '),ixor
          call prntm2(bc(ixor),nh,nh,nh)                                  4d6s22
         end if
         ibcoff=ixor+nh*nh
         ioff=ioff+nhp*nhp                                               4d4s22
         ioffb=ioffb+nh*nh*2                                            6d3s22
        end if                                                          5d20s22
       end do
      else
       ioff=1
       ioffp=1                                                          3d29s16
       n2ndsz=0                                                          6d21s16
       nx=0                                                             7d5s16
       do isb=1,nsymb                                                    6d21s16
        isk=multh(isb,ipropsym)                                         7d5s16
        nx=nx+nbasdws(isb)**2                                           4d14s22
        n2ndsz=n2ndsz+nbasdws(isk)*nbasdws(isb)                         4d4s22
       end do                                                            6d21s16
       ixor0=ibcoff                                                     7d5s16
       ibcoff=ixor0+n2ndsz                                              5d4s22
       ibcoffo=ibcoff                                                   5d4s22
       call enough('derofxor1.  5',bc,ibc)
       ixor=ixor0                                                       7d5s16
       do isb=1,nsymb
        isk=multh(isb,ipropsym)
        ioffp=ioffp+nbasdws(isk)*nbasdws(isb)                           4d4s22
        if(isk.gt.isb)then
         nk=nbasdws(isk)                                                4d4s22
         nb=nbasdws(isb)                                                4d4s22
         nkp=nbasisp(isk)*ncomp                                         4d4s22
         nbp=nbasisp(isb)*ncomp                                         4d4s22
         itmp1=ibcoff
         itmp2=itmp1+nbp*nkp                                            4d4s22
         ibcoff=itmp2+nbp*nkp                                           4d4s22
         call enough('derofxor1.  6',bc,ibc)
         iad=ivecso(isk)+nk*nk                                          4d6s22
         if(min(nbp,nkp).gt.0)then                                      7d11s22
          ioff=iso(isb)+1                                                6d27s22
          call dgemm('n','n',nbp,nk,nkp,1d0,ovr(ioff),nbp,bc(iad),nkp,   4d4s22
     $        0d0,bc(itmp1),nbp,                                        4d4s22
     d' derofxor. 24')
          ioff=ioff+nbp*nkp                                              4d4s22
          do i=0,nbp-1                                                   4d4s22
           do j=0,nk-1
            ji=j+nk*i+itmp2
            ij=i+nbp*j+itmp1                                             4d4s22
            bc(ji)=bc(ij)
           end do
          end do
          iad=ivecso(isb)+nb*nb                                          4d6s22
          call dgemm('n','n',nk,nb,nbp,1d0,bc(itmp2),nk,bc(iad),nbp,     4d4s22
     $        0d0,bc(itmp1),nk,
     d' derofxor. 25')
         end if                                                         7d11s22
         iptb(isb)=ioffp-1                                              7d6s16
         do i=0,nb-1
          do j=0,nk-1
           iad2=itmp2+i+nb*j
           iad1=itmp1+j+nk*i
           bc(iad2)=bc(iad1)
          end do
         end do
         if(idwsdeb.gt.10)then                                          4d26s18
          write(6,*)('dS in diagonal basis ')
          call prntm2(bc(itmp2),nb,nk,nb)
         end if                                                         4d26s18
         ise=ibcoff
         ibcoff=ise+nk
         call enough('derofxor1.  7',bc,ibc)
         do j=0,nk-1
          bc(ise+j)=1d0/(bc(ieigs(isk)+j)**2)
         end do
         do i=0,nb-1
          val=1d0/(bc(ieigs(isb)+i)**2)
          do j=0,nk-1
           ji=itmp1+j+nk*i
           ji2=itmp2+j+nk*i
           bot=bc(ise+j)-val
           if(abs(bot).lt.1d-12)then
            if(abs(bc(ji)).gt.1d-10)then
             write(6,*)('small denominator '),j,i,bot,
     $            ('but not so small numerator '),bc(ji),bc(ise+j),val
            end if                                                      4d6s23
            bot=1d0                                                     4d4s23
           end if
           bc(ji2)=-bc(ji)/bot
          end do
         end do
         ibcoff=ise
c
c     itmp2 now contains VtdV with isk labels on left and isb labels on right.
c
         if(min(nk,nb).gt.0)then                                        7d11s22
          if(idwsdeb.ne.0)then                                          7d14s22
           write(6,*)('VtdV ')
           call prntm2(bc(itmp2),nk,nb,nk)
          end if                                                        7d14s22
          iad=ivecso(isk)                                                4d1s22
          call dgemm('n','n',nk,nb,nk,1d0,bc(iad),nk,bc(itmp2),nk,       4d6s22
     $        0d0,bc(itmp1),nk,                                         4d6s22
     d' derofxor. 26')
          itmp3=ibcoff
          ibcoff=itmp3+nb*nk
          call enough('derofxor1.  8',bc,ibc)
          do i=0,nb-1
           do j=0,nk-1
            iad2=itmp2+j+nk*i
            iad3=itmp3+i+nb*j
            bc(iad3)=-bc(iad2)
           end do
          end do
          iad=ivecso(isb)                                                4d1s22
          call dgemm('n','n',nb,nk,nb,1d0,bc(iad),nb,bc(itmp3),nb,       4d6s22
     $        0d0,bc(ixor),nb,                                          4d6s22
     d' derofxor. 27')
         end if                                                         7d11s22
c
c     itmp1 now contains dV with isk on left and isb on right, and
c     ixor now contains  dV with isb on the left and isk on the right.
c     this perturbation does not effect the eigenvalues.
c     for our final result, we want bk.
c
c
         iad=ivecso(isb)                                                4d1s22
c
c     compute smh dVt and smh Vt
c
c
         do i=0,nb-1                                                    4d6s22
          do j=0,nb-1                                                   4d4s22
           ji=j+nb*i                                                    4d4s22
           ij=i+nb*j                                                    4d6s22
           iad1=itmp3+ji
           iad2=iad+ij
           bc(iad1)=bc(ieigs(isb)+j)*bc(iad2)
          end do
          do j=0,nk-1
           ji=j+nk*i
           ij=i+nb*j                                                    4d6s22
           iad3=itmp2+ji
           iad4=ixor+ij
           bc(iad3)=bc(ieigs(isk)+j)*bc(iad4)
          end do
         end do
         iad=ivecso(isk)                                                4d1s22
         if(min(nk,nb).gt.0)then                                        7d11s22
          call dgemm('n','n',nk,nb,nk,1d0,bc(iad),nk,bc(itmp2),nk,       4d6s22
     $        0d0,bc(ixor),nk,                                          4d6s22
     d' derofxor. 28')
c
c     for this mxm, we need dv kb, ie tmp1
c
          call dgemm('n','n',nk,nb,nb,+1d0,bc(itmp1),nk,bc(itmp3),nb,    4d6s22
     $        1d0,bc(ixor),nk,                                          4d6s22
     d' derofxor. 29')
          do i=0,nb-1                                                    4d6s22
           do j=0,nk-1                                                   4d6s22
            ji=j+nk*i+ixor                                               4d6s22
            ij=i+nb*j+itmp2                                              4d6s22
            bc(ij)=bc(ji)                                                4d4s22
           end do
          end do
          if(idwsdeb.gt.10)then
           write(6,*)('dXor: '),ixor
           call prntm2(bc(itmp2),nb,nk,nb)                                4d6s22
          end if
          do i=0,nb*nk-1                                                 4d6s22
           bc(ixor+i)=bc(itmp2+i)
          end do
         end if                                                         7d11s22
         ibcoff=itmp1
         ioffp=ioffp+nb*nk                                              3d29s16
         ixor=ixor+nb*nk                                                4d6s22
        end if
       end do
       ibcoff=ibcoffo
      end if                                                            5d4s22
      return
      end
