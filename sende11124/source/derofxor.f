c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine derofxor(ovr,iso,ivecs,ipropsym,idorel,icanon,ieigs,
     $     ixor0,multh,idwsdeb,iorb,propmat,ovrdk,ovrd2,ovrdkd2,        4d8s22
     $     ivecso,nbasisp,isob,bc,ibc)                                  11d14s22
c
c     differentiate ao to orthogonal basis transformation.
c
      implicit real*8 (a-h,o-z)
      include "common.hf"
      include "common.store"
      include "common.spher"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension ovr(1),iso(1),ivecs(1),ieigs(1),multh(8,8),
     $     icanon(8),iorb(1),propmat(1),ovrdk(1),ovrd2(1),              4d8s22
     $     ovrdkd2(1),iptb(8),ivecso(*),nbasisp(*),isob(*)              4d8s22
      ibcoffo=ibcoff
      ncomp=1
      if(idorel.ne.0)ncomp=2
      if(ipropsym.eq.1)then
       n2ndsz=0                                                          6d21s16
       do isb=1,nsymb                                                    6d21s16
        n2ndsz=n2ndsz+nbasdws(isb)**2                                   4d4s22
       end do                                                            6d21s16
       ioff=1
       ioffb=1                                                          4d15s22
       do isb=1,nsymb
        nh=nbasdws(isb)                                                 4d4s22
        nhp=nbasisp(isb)*ncomp                                          4d4s22
        ixor=ibcoff
        if(isb.eq.1)ixor0=ixor
        ibcoff=ixor+nhp*nhp*2                                           4d4s22
        ibcoffo=ibcoff
        itmp1=ibcoff
        itmp2=itmp1+nhp*nhp                                             4d4s22
        itmp3=itmp2+nhp*nhp                                             4d4s22
        ibcoff=itmp3+nhp*nhp                                            4d4s22
        call enough('derofxor.  1',bc,ibc)
        call square(ovr(ioff),nhp)                                      4d4s22
        call square(ovrd2(ioff),nhp)                                    4d4s22
c                                                                       3d29s16
c     first part of (n|d/dq|m) matrix is der of overlap transformed     3d29s16
c              and (n|d2/dq2|m)                                         6d20s16
c     to mo basis.                                                      3d29s16
c                                                                       3d29s16
        if(idwsdeb.gt.10)then
         write(6,*)('ovrdk: '),ioff
         call prntm2(ovrdk(ioff),nhp,nhp,nhp)
         write(6,*)('let me add transpose: this should be over')
         itrantmp=ibcoff
         ibcoff=itrantmp+nhp*nhp
         call enough('derofxor.  2',bc,ibc)
         jtrantmp=itrantmp
         do i=0,nhp-1
          do j=0,i
           ji=j+nhp*i
           ij=i+nhp*j
           bc(jtrantmp)=ovrdk(ioff+ij)+ovrdk(ioff+ji)
           jtrantmp=jtrantmp+1
          end do
         end do
         call mpprnt2(bc(itrantmp),nhp)
         call square(bc(itrantmp),nhp)
         call prntm2(bc(itrantmp),nhp,nhp,nhp)
         ibcoff=itrantmp
         write(6,*)('times orbs ')
         call prntm2(bc(iorb(isb)),nhp,nh,nhp)
        end if
        call dgemm('n','n',nhp,nh,nhp,1d0,ovrdk(ioff),nhp,              4d4s22
     $       bc(iorb(isb)),nhp,0d0,bc(itmp1),nhp,                       4d4s22
     d' derofxor.  1')
        do i=0,nhp-1                                                    4d4s22
         do j=0,nh-1                                                    3d29s16
          ij=itmp1+i+nhp*j                                              4d4s22
          ji=itmp2+j+nh*i                                               3d29s16
          bc(ji)=bc(ij)                                                 3d29s16
         end do                                                         3d29s16
        end do                                                          3d29s16
        call dgemm('n','n',nh,nh,nhp,1d0,bc(itmp2),nh,bc(iorb(isb)),nhp,4d4s22
     $       0d0,bc(itmp1),nh,                                          3d29s16
     d' derofxor.  2')
        if(idwsdeb.gt.10)then
         xmul=1d0
         if(nsymb.eq.1)xmul=2d0
         sum=0d0
         do i=0,nhp-1
          iad1=itmp2+nh*i
          val=bc(iad1)*xmul
          iad2=iorb(isb)+i
          term=val*bc(iad2)
          sum=sum+term
          if(abs(term).gt.1d-10)write(6,*)i,val,bc(iad2),sum
         end do
         write(6,*)('yields '),ioff,sum
         call prntm2(bc(itmp1),nh,nh,nh)
        end if
        do i=0,nh-1                                                     3d29s16
         do j=0,nh-1                                                    3d29s16
          ij=itmp1+i+nh*j                                               3d29s16
          ji=ioffb+j+nh*i                                               4d15s22
          propmat(ji)=bc(ij)                                            3d29s16
         end do                                                         3d29s16
        end do                                                          3d29s16
        if(idwsdeb.gt.10)then
         write(6,*)('starting ovrdk2: ')
         call prntm2(ovrdkd2(ioff),nhp,nhp,nhp)                         4d4s22
        end if
        call dgemm('n','n',nhp,nh,nhp,1d0,ovrdkd2(ioff),nhp,            4d4s22
     $       bc(iorb(isb)),nhp,0d0,bc(itmp1),nhp,                       4d4s22
     d' derofxor.  3')
        do i=0,nhp-1                                                    4d4s22
         do j=0,nh-1                                                    3d29s16
          ij=itmp1+i+nhp*j                                              4d4s22
          ji=itmp2+j+nh*i                                               3d29s16
          bc(ji)=bc(ij)                                                 3d29s16
         end do                                                         3d29s16
        end do                                                          3d29s16
        call dgemm('n','n',nh,nh,nhp,1d0,bc(itmp2),nh,bc(iorb(isb)),nhp,4d4s22
     $       0d0,bc(itmp1),nh,                                          3d29s16
     d' derofxor.  4')
        do i=0,nh-1                                                     3d29s16
         do j=0,nh-1                                                    3d29s16
          ij=itmp1+i+nh*j                                               3d29s16
          ji=ioffb+n2ndsz+j+nh*i                                        4d15s22
          propmat(ji)=bc(ij)                                            3d29s16
         end do                                                         3d29s16
        end do                                                          3d29s16
        if(idwsdeb.gt.10)then                                           12d23s16
         write(6,*)('first part of propmat: '),ioffb,loc(propmat(ioffb))
         itry=(loc(propmat(ioffb))-loc(bc))/8
         itry=itry+1
         write(6,*)('itry? '),itry,bc(itry),loc(bc(itry))
         write(6,*)('first der')
         call prntm2(propmat(ioffb),nh,nh,nh)
         if(nsymb.eq.1)then
          call printa(propmat(ioffb),nh,0,1,0,nh,0,1,0,bc(ibcoff))
         end if
         ioff2=ioffb+n2ndsz                                             4d15s22
         write(6,*)('2nd der')
         call prntm2(propmat(ioff2),nh,nh,nh)
         if(nsymb.eq.1)call printa(propmat(ioff2),nh,0,1,0,nh,0,1,0,
     $        bc(ibcoff))
        end if
        iad=ivecso(isb)+nh*nh                                           4d6s22
        call dgemm('n','n',nhp,nh,nhp,1d0,ovr(ioff),nhp,bc(iad),nhp,0d0,4d8s22
     $       bc(itmp1),nhp,                                             4d4s22
     d' derofxor.  5')
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
        ids=ibcoff
        ibcoff=ids+nh
        ise=ibcoff
        ibcoff=ise+nh
        call enough('derofxor.  5',bc,ibc)
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
           write(6,*)('small denominator ...'),i,j,bot,bc(ise+i),
     $          bc(ise+j),bc(ji2)
           bot=1d0
          end if
          bc(ji)=-bc(ji2)/bot                                           1d21s16
          bc(ij)=-bc(ji)
         end do
        end do
        if(idwsdeb.gt.10)then
        write(6,*)('Vt dV ')
        call prntm2(bc(itmp2),nh,nh,nh)
        if(nsymb.eq.1)call printb(bc(itmp2),nh,0,1,0,nh,0,1,0,
     $       bc(ibcoff))
        end if
        ids2=ibcoff                                                     6d21s16
        ibcoff=ids2+nh*nh                                               6d21s16
        call enough('derofxor.  6',bc,ibc)
        call dgemm('n','n',nh,nh,nh,1d0,bc(itmp1),nh,bc(itmp2),nh,0d0,  6d21s16
     $       bc(ids2),nh,                                               6d21s16
     d' derofxor.  7')
        do i=0,nh-1                                                     6d21s16
         do j=0,i                                                       6d21s16
          ij=i+nh*j+ids2                                                6d21s16
          ji=j+nh*i+ids2                                                6d21s16
          sum=bc(ij)+bc(ji)                                             6d21s16
          bc(ij)=sum                                                    6d21s16
          bc(ji)=sum                                                    6d21s16
         end do                                                         6d21s16
        end do                                                          6d21s16
        iad=ivecso(isb)+nh*nh                                           4d6s22
        call dgemm('n','n',nhp,nh,nhp,1d0,ovrd2(ioff),nhp,bc(iad),nhp,  4d4s22
     $       0d0,bc(itmp1),nhp,                                         4d4s22
     d' derofxor.  8')
        do i=0,nh-1                                                     6d21s16
         do j=0,nhp-1                                                    6d21s16
          ji=itmp1+j+nhp*i                                              4d4s22
          ij=itmp3+i+nh*j                                               6d21s16
          bc(ij)=bc(ji)                                                 6d21s16
         end do                                                         6d21s16
        end do                                                          6d21s16
        call dgemm('n','n',nh,nh,nhp,1d0,bc(itmp3),nh,bc(iad),nhp,0d0,  4d4s22
     $       bc(itmp1),nh,                                               6d21s16
     d' derofxor.  9')
        call dgemm('n','n',nh,nh,nhp,1d0,bc(itmp3),nh,bc(iad),nhp,1d0,  4d4s22
     $       bc(ids2),nh,                                               6d21s16
     d' derofxor. 10')
        iders2=ibcoff                                                   6d21s16
        ibcoff=iders2+nh                                                6d21s16
        call enough('derofxor.  7',bc,ibc)
        jders2=iders2-1                                                 6d21s16
        do i=1,nh                                                       6d21s16
         iad=ids2+(nh+1)*(i-1)                                          6d21s16
         bc(jders2+i)=bc(iad)                                           6d21s16
        end do                                                          6d21s16
        do i=0,nh-1                                                     6d21s16
         iad=ids2+(nh+1)*i                                              6d21s16
         bc(iad)=0d0                                                    6d21s16
         do j=0,i-1                                                     6d21s16
          bot=bc(ise+j)-bc(ise+i)                                       6d21s16
          ds=bc(ids+j)-bc(ids+i)                                        6d22s16
          itop=ids2+j+nh*i                                              6d21s16
          if(abs(bot).lt.1d-14)then                                     6d21s16
           write(6,*)('small denominator ...'),j,i,bot,bc(itop)         6d21s16
           bot=1d0                                                      6d21s16
          end if                                                        6d21s16
          boti=1d0/bot                                                  6d21s16
          itop2=itmp2+j+nh*i                                            6d21s16
          bc(itop)=-(bc(itop2)*ds+bc(itop))*boti                         6d21s16
          itopt=ids2+i+nh*j                                             6d21s16
          bc(itopt)=-bc(itop)                                           6d21s16
         end do                                                         6d21s16
        end do                                                          6d21s16
        if(idwsdeb.gt.10)then
         write(6,*)('dVT*dV+VT*d2 V')
         call prntm2(bc(ids2),nh,nh,nh)
         if(nsymb.eq.1)then
          call printb(bc(ids2),nh,0,1,0,nh,0,1,0,bc(ibcoff))
         end if
        end if
        iad=ivecso(isb)+nh*n0                                           4d6s22
        call dgemm('n','n',nh,nh,nh-n0,1d0,bc(iad),nh,bc(itmp2+n0),nh,  4d6s22
     $       0d0,bc(itmp1),nh,                                          4d6s22
     d' derofxor. 11')
        if(idwsdeb.gt.10)then
        write(6,*)('dV ')
        call prntm2(bc(itmp1),nh,nh,nh)                                 4d6s22
        end if
        do i=0,nh-1                                                     6d21s16
         do j=0,nh-1                                                    4d6s22
          ji=itmp1+j+nh*i                                               4d6s22
          ij=itmp3+i+nh*j                                               6d21s16
          bc(ij)=-bc(ji)                                                6d21s16
         end do                                                         6d21s16
        end do                                                          6d21s16
        if(idwsdeb.gt.10)then
         call dgemm('n','n',nh,nh,nh,1d0,bc(itmp3),nh,bc(itmp1),nh,     4d6s22
     $        0d0,bc(ibcoff),nh,                                        4d4s22
     d' derofxor. 12')
         write(6,*)('dVT dV')
         call prntm2(bc(ibcoff),nh,nh,nh)
         if(nsymb.eq.1)then
          call printb(bc(ibcoff),nh,0,1,0,nh,0,1,0,bc(ibcoff+nh*nh))
         end if
        end if
        call dgemm('n','n',nh,nh,nh,1d0,bc(itmp3),nh,bc(itmp1),nh,1d0,  4d6s22
     $       bc(ids2),nh,                                               6d21s16
     d' derofxor. 13')
        if(idwsdeb.gt.10)then                                           6d21s16
         write(6,*)('VT*d2V')                                           6d21s16
         call prntm2(bc(ids2),nh,nh,nh)                                 6d21s16
         if(nsymb.eq.1)call printb(bc(ids2),nh,0,1,0,nh,0,1,0,
     $        bc(ibcoff))
        end if                                                          6d21s16
        iad=ivecso(isb)                                                 4d6s22
        call dgemm('n','n',nh,nh,nh,1d0,bc(iad),nh,bc(ids2),nh,0d0,     4d6s22
     $       bc(itmp3),nh,                                              4d6s22
     d' derofxor. 14')
        do ij=0,nh*nh-1                                                 4d6s22
         bc(ids2+ij)=bc(itmp3+ij)                                       6d21s16
        end do                                                          6d21s16
        if(idwsdeb.gt.10)then                                           6d21s16
         write(6,*)('d2V ')                                             6d21s16
         call prntm2(bc(ids2),nh,nh,nh)                                 4d6s22
        end if                                                          6d21s16
        call enough('derofxor.  8',bc,ibc)
        if(icanon(isb).eq.0)then                                        1d29s16
         if(idwsdeb.gt.10)write(6,*)('symmetric orthogonalization: ')
         itmp4=ibcoff                                                   6d22s16
         ibcoff=itmp4+nhp*nh                                            4d4s22
         itmp5=ibcoff
         ibcoff=itmp5+nhp*nh                                            4d4s22
         call enough('derofxor.  9',bc,ibc)
         do i=0,nh-1
          smh=bc(ieigs(isb)+i)
          part1=-0.5d0*smh*smh*smh
          dsmh=part1*bc(ids+i)
          d2smh=part1*(bc(iders2+i)-1.5d0*((smh*bc(ids+i))**2))
          d2smhx=part1*bc(iders2+i)
          bc(ise+i)=dsmh
          do j=0,nh-1                                                   4d6s22
           ij=i+j*nh
           ji=j+i*nh                                                    4d6s22
           iad1=itmp2+ij
           iad2=itmp1+ji
           iad3=ivecso(isb)+ji                                          4d4s22
           iad4=itmp3+ij
           iad5=itmp4+ij
           iad6=ids2+ji
c     s^{-1/2}VsT: dVs times goes to dX, d2Vs times goes to d2X
           bc(iad4)=smh*bc(iad3)
c     s^{-1/2}dVsT-.5 ds s^{-3/2} VsT, VS times goes to dX, 2dVs times goes to d2X
           bc(iad1)=smh*bc(iad2)+dsmh*bc(iad3)
c    [-.5 d2s s^{-3/2}+1.5 ds*ds s^{-5/2}] VsT-ds s^{-3/2}dVsT+s^{-1/2}d2VsT: Vs times goes to d2X
           bc(iad5)=d2smh*bc(iad3)+2d0*dsmh*bc(iad2)+smh*bc(iad6)
          end do
         end do
         iad=ivecso(isb)                                                4d1s22
         call dgemm('n','n',nh,nh,nh,1d0,bc(iad),nh,bc(itmp2),nh,0d0,   4d9s22
     $        bc(ixor),nh,                                              4d6s22
     d' derofxor. 15')
        ixort=ibcoff
        ibcoff=ixort+nh*nh
        call enough('derofxor. 10',bc,ibc)
        do i=0,nh-1
         do j=0,nh-1
          ji=ixor+j+nh*i
          ij=ixort+i+nh*j
          bc(ij)=bc(ji)
         end do
        end do
        do i=0,nh*nh-1
         bc(ixort+i)=bc(ixort+i)+bc(ixor+i)
        end do
        ibcoff=ixort
         ixor2=ixor+nh*nh                                               4d6s22
         call dgemm('n','n',nh,nh,nh,1d0,bc(ids2),nh,bc(itmp3),nh,      4d9s22
     $        0d0,bc(ixor2),nh,                                         4d6s22
     d' derofxor. 16')
         call dgemm('n','n',nh,nh,nh,2d0,bc(itmp1),nh,bc(itmp2),nh,     4d6s22
     $        1d0,bc(ixor2),nh,                                         4d6s22
     d' derofxor. 17')
         call dgemm('n','n',nh,nh,nh,1d0,bc(iad),nh,bc(itmp4),nh,1d0,   4d6s22
     $        bc(ixor2),nh,                                             4d6s22
     d' derofxor. 18')
         call dgemm('n','n',nh,nh,nh,1d0,bc(itmp1),nh,bc(itmp3),nh,     4d6s22
     $        1d0,bc(ixor),nh,                                          4d6s22
     d' derofxor. 19')
        else
         write(6,*)('cannonical orthogonalization: ')
         write(6,*)('no code for this yet in derofxor ')
         call dws_sync
         call dws_finalize
         stop
        end if
        if(idwsdeb.gt.10)then
        write(6,*)('derivative of orthogonalization matrix '),ixor
        call prntm2(bc(ixor),nh,nh,nh)                                  4d6s22
        if(nsymb.eq.1)call printc(bc(ixor),bc(ibcoff))
        write(6,*)('2nd derivative of orthogonalization matrix ')
        call prntm2(bc(ixor2),nh,nh,nh)                                 4d6s22
        if(nsymb.eq.1)call printc(bc(ixor2),bc(ibcoff))
        end if
        ibcoff=ixor2+nh*nh                                              4d6s22
        ioff=ioff+nhp*nhp                                               4d4s22
        ioffb=ioffb+nh*nh                                               4d15s22
       end do
      else
       ioff=1
       ioffp=1                                                          3d29s16
       n2ndsz=0                                                          6d21s16
       nx=0                                                             7d5s16
       do isb=1,nsymb                                                    6d21s16
        isk=multh(isb,ipropsym)                                         7d5s16
        nx=nx+nbasdws(isb)**2                                           4d14s22
         n2ndsz=n2ndsz+nbasdws(isk)*nbasdws(isb)                        4d4s22
       end do                                                            6d21s16
       ixor0=ibcoff                                                     7d5s16
       ixor2=ixor0+n2ndsz                                               4d6s22
       idssav=ixor0+n2ndsz+nx                                           4d6s22
       ivdvsav=idssav+n2ndsz                                            7d5s16
       idvsav=ivdvsav+n2ndsz                                            7d5s16
       ibcoff=idvsav+n2ndsz                                             7d5s16
       ibcoffo=idssav                                                   7d19s16
       call enough('derofxor. 11',bc,ibc)
       ixor=ixor0                                                       7d5s16
       do isb=1,nsymb
        isk=multh(isb,ipropsym)
         ioffp=ioffp+nbasdws(isk)*nbasdws(isb)                          4d4s22
        if(isk.gt.isb)then
         if(icanon(isk).ne.0.or.icanon(isb).ne.0)then
          write(6,*)('cannonical orthogonalization: ')
          write(6,*)('no code for this yet in derofxor ')
          call dws_sync
          call dws_finalize
          stop
         end if
         nk=nbasdws(isk)                                                4d4s22
         nb=nbasdws(isb)                                                4d4s22
         if(min(nk,nb).gt.0)then                                        11d28s22
          nkp=nbasisp(isk)*ncomp                                         4d4s22
          nbp=nbasisp(isb)*ncomp                                         4d4s22
          itmp1=ibcoff
          itmp2=itmp1+nbp*nkp                                            4d4s22
          ibcoff=itmp2+nbp*nkp                                           4d4s22
          call enough('derofxor. 12',bc,ibc)
c                                                                       3d29s16
c     first part of (n|d/dq|m) matrix is der of overlap transformed     3d29s16
c     to mo basis.                                                      3d29s16
c                                                                       3d29s16
          ioff=iso(isb)+1                                                4d15s22
          if(idwsdeb.gt.10)then
           write(6,*)('ioff vs. isos: '),ioff,iso(isb),iso(isk),isb,isk
           write(6,*)('bk '),ioff
           call prntm2(ovrdk(ioff),nbp,nkp,nbp)                           4d4s22
           call prntm2(ovrdk(iso(isb)+1),nbp,nkp,nbp)                           4d4s22
           write(6,*)('kb: '),iso(isk)+1
           call prntm2(ovrdk(iso(isk)+1),nkp,nbp,nkp)                     4d4s22
          end if
          call dgemm('n','n',nbp,nk,nkp,1d0,ovrdk(ioff),nbp,             4d4s22
     $        bc(iorb(isk)),nkp,0d0,bc(itmp1),nbp,                      4d4s22
     d' derofxor. 20')
          do i=0,nbp-1                                                   4d4s22
           do j=0,nk-1                                                    3d29s16
            ij=itmp1+i+nbp*j                                             4d4s22
            ji=itmp2+j+nk*i                                               3d29s16
            bc(ji)=bc(ij)                                                 3d29s16
           end do                                                         3d29s16
          end do                                                          3d29s16
          call dgemm('n','n',nk,nb,nbp,1d0,bc(itmp2),nk,bc(iorb(isb)),   4d4s22
     $        nbp,0d0,bc(itmp1),nk,                                     4d4s22
     d' derofxor. 21')
          ioffp=isob(isb)+1                                              4d15s22
          do i=0,nk-1                                                     3d29s16
           do j=0,nb-1                                                    3d29s16
            ij=itmp1+i+nk*j                                               3d29s16
            ji=ioffp+j+nb*i                                               3d29s16
            propmat(ji)=bc(ij)                                            3d29s16
           end do                                                         3d29s16
          end do                                                          3d29s16
          if(idwsdeb.gt.10)then
           write(6,*)('propmat: der of overlap part of (n|d/dq|m) '),
     $          isb,ioffp,isob(isb)
           call prntm2(propmat(ioffp),nb,nk,nb)                             3d29s16
          end if
          call dgemm('n','n',nkp,nb,nbp,1d0,ovrdk(iso(isk)+1),nkp,       4d4s22
     $        bc(iorb(isb)),nbp,0d0,bc(itmp1),nkp,                      4d4s22
     d' derofxor. 22')
          do i=0,nkp-1                                                   4d4s22
           do j=0,nb-1                                                   7d21s16
            ij=itmp1+i+nkp*j                                             4d4s22
            ji=itmp2+j+nb*i                                              7d21s16
            bc(ji)=bc(ij)                                                7d21s16
           end do                                                        7d21s16
          end do                                                         7d21s16
          call dgemm('n','n',nb,nk,nkp,1d0,bc(itmp2),nb,bc(iorb(isk)),   4d4s22
     $        nkp,0d0,bc(itmp1),nb,                                     4d4s22
     d' derofxor. 23')
          do i=0,nb-1                                                    7d21s16
           do j=0,nk-1                                                   7d21s16
            ij=itmp1+i+nb*j                                              7d21s16
            ji=isob(isk)+1+j+nk*i                                        4d8s22
            propmat(ji)=bc(ij)                                           7d21s16
           end do                                                        7d21s16
          end do                                                         7d21s16
          if(idwsdeb.gt.10)then
           write(6,*)('propmat: reverse symmetry block of (n|d/dq|m)'),   4d15s22
     $         isk,isob(isk)+1
           call prntm2(propmat(isob(isk)+1),nk,nb,nk)                     4d8s22
          end if
          iad=ivecso(isk)+nk*nk                                          4d6s22
          ioff=iso(isb)+1                                                8d4s22
          if(min(nbp,nkp).gt.0)then                                      7d11s22
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
            iad3=idssav+i+nb*j+ioffp-1                                   7d6s16
            bc(iad3)=bc(iad2)                                            7d5s16
           end do
          end do
          ise=ibcoff
          ibcoff=ise+nk
          call enough('derofxor. 13',bc,ibc)
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
             write(6,*)('small denominator '),j,i,bot,bc(ise+j),val,
     $           bc(ji)
            end if
            bc(ji2)=-bc(ji)/bot
           end do
          end do
          ibcoff=ise
c
c     itmp2 now contains VtdV with isk labels on left and isb labels on right.
c
          iad3=ivdvsav+isob(isk)                                         4d15s22
          do i=0,nk*nb-1                                                 7d5s16
           bc(iad3+i)=bc(itmp2+i)                                        7d5s16
          end do                                                         7d5s16
          iad=ivecso(isk)                                                4d1s22
          call dgemm('n','n',nk,nb,nk,1d0,bc(iad),nk,bc(itmp2),nk,       4d6s22
     $        0d0,bc(itmp1),nk,                                         4d6s22
     d' derofxor. 26')
          iad3=idvsav+isob(isk)                                          4d15s22
          do i=0,nk*nb-1                                                 4d6s22
           bc(iad3+i)=bc(itmp1+i)                                        7d5s16
          end do                                                         7d5s16
          itmp3=ibcoff
          ibcoff=itmp3+nb*nk
          call enough('derofxor. 14',bc,ibc)
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
          ibcoff=itmp1
         end if                                                         11d28s22
         ioffp=ioffp+nb*nk                                              3d29s16
         ixor=ixor+nb*nk                                                4d6s22
        end if
       end do
c                                                                       7d5s16
c     now for d2 terms, which are all symmetric                         7d5s16
c                                                                       7d5s16
       ioff=1                                                           7d5s16
       ioffp=1+n2ndsz                                                   4d18s22
       jxor2=ixor2                                                      7d8s16
       do isb=1,nsymb                                                   7d5s16
        nh=nbasdws(isb)                                                 4d4s22
        nhp=nbasisp(isb)*ncomp                                          4d4s22
        isk=multh(isb,ipropsym)                                         7d5s16
        nk=nbasdws(isk)                                                 4d4s22
        nkp=nbasisp(isk)*ncomp                                          4d4s22
        nhk=max(nhp,nkp)                                                4d4s22
        itmp1=ibcoff                                                    7d5s16
        itmp2=itmp1+nhk*nhk                                             7d11s16
        itmp3=itmp2+nhk*nhk                                             7d11s16
        itmp4=itmp3+nhk*nhk                                             7d11s16
        itmp5=itmp4+nhk*nhk                                             7d12s16
        ibcoff=itmp5+nhk*nhk                                            7d12s16
        call enough('derofxor. 15',bc,ibc)
        if(idwsdeb.gt.10)then
        write(6,*)('starting ovrdkd2 for symmetry block '),isb,ioff
        call prntm2(ovrdkd2(ioff),nhp,nhp,nhp)                          4d4s22
        end if
        if(min(nhp,nh).gt.0)then                                        11d28s22
         call dgemm('n','n',nhp,nh,nhp,1d0,ovrdkd2(ioff),nhp,            4d4s22
     $       bc(iorb(isb)),nhp,0d0,bc(itmp1),nhp,                       4d4s22
     d' derofxor. 30')
         do i=0,nhp-1                                                    4d4s22
          do j=0,nh-1                                                    7d5s16
           ij=itmp1+i+nhp*j                                              4d4s22
           ji=itmp2+j+nh*i                                               7d5s16
           bc(ji)=bc(ij)                                                 7d5s16
          end do                                                         7d5s16
         end do                                                          7d5s16
         call dgemm('n','n',nh,nh,nhp,1d0,bc(itmp2),nh,bc(iorb(isb)),   11d28s22
     $        nhp,0d0,bc(itmp1),nh,                                     11d28s22
     d' derofxor. 31')
        end if                                                          11d28s22
        do i=0,nh-1                                                     7d5s16
         do j=0,nh-1                                                    7d5s16
          ij=itmp1+i+nh*j                                               7d5s16
          ji=ioffp+j+nh*i                                               7d5s16
          propmat(ji)=bc(ij)                                            7d5s16
         end do                                                         7d5s16
        end do                                                          7d5s16
        if(idwsdeb.gt.10)then
        write(6,*)('propmat2: '),ioffp
        call prntm2(propmat(ioffp),nh,nh,nh)                            7d21s16
        end if
c     dVT dS ... dS is saved b*k, and dV is saved k*b. dS is symmetric
c     while dV is antisymmetic.
c     dV*dS is k*k and dVT*dST is b*b. but one has minus sign ...
c
        if(isk.gt.isb)then                                              7d5s16
         if(nk.gt.0)then                                                11d28s22
          call dgemm('n','n',nh,nh,nk,1d0,bc(idssav+isob(isb)),nh,       4d15s22
     $                   bc(ivdvsav+isob(isk)),nk,0d0,bc(itmp1),nh,     4d15s22
     d' derofxor. 32')
         else                                                           11d28s22
          do iz=itmp1,itmp1+nh*nh-1                                     11d28s22
           bc(iz)=0d0                                                   11d28s22
          end do                                                        11d28s22
         end if                                                         11d28s22
        else                                                            7d5s16
         if(min(nh,nk).gt.0)then                                        11d28s22
          call dgemm('n','n',nh,nh,nk,-1d0,bc(ivdvsav+isob(isb)),nh,     4d15s22
     $        bc(idssav+isob(isk)),nk,                                  4d15s22
     $        0d0,bc(itmp1),nh,                                         7d5s16
     d' derofxor. 33')
         else                                                           11d28s22
          do iz=itmp1,itmp1+nh*nh-1                                     11d28s22
           bc(iz)=0d0                                                   11d28s22
          end do                                                        11d28s22
         end if                                                         11d28s22
        end if                                                          7d5s16
        do i=0,nh-1                                                     7d6s16
         do j=i,nh-1                                                    7d6s16
          ji=itmp1+j+nh*i                                               7d6s16
          ij=itmp1+i+nh*j                                               7d6s16
          sum=bc(ji)+bc(ij)                                             7d6s16
          bc(ji)=sum                                                    7d6s16
          bc(ij)=sum                                                    7d6s16
         end do                                                         7d6s16
        end do                                                          7d6s16
        call square(ovrd2(ioff),nhp)                                    4d4s22
        iad=ivecso(isb)+nh*nh                                           4d6s22
        if(min(nhp,nh).gt.0)then                                        11d28s22
         call dgemm('n','n',nhp,nh,nhp,1d0,ovrd2(ioff),nhp,bc(iad),nhp,  4d4s22
     $       0d0,bc(itmp2),nhp,                                         4d4s22
     d' derofxor. 34')
         do i=0,nh-1                                                     6d21s16
          do j=0,nhp-1                                                   4d4s22
           ji=itmp2+j+nhp*i                                              4d4s22
           ij=itmp3+i+nh*j                                               6d21s16
           bc(ij)=bc(ji)                                                 6d21s16
          end do                                                         6d21s16
         end do                                                          6d21s16
         call dgemm('n','n',nh,nh,nhp,1d0,bc(itmp3),nh,bc(iad),nhp,0d0,  4d4s22
     $       bc(itmp2),nh,                                               6d21s16
     d' derofxor. 35')
         call dgemm('n','n',nh,nh,nhp,1d0,bc(itmp3),nh,bc(iad),nhp,1d0,  4d14s22
     $       bc(itmp1),nh,                                               6d21s16
     d' derofxor. 36')
        else                                                            11d28s22
         do iz=itmp2,itmp2+nh*nh-1                                      11d28s22
          bc(iz)=0d0                                                    11d28s22
         end do                                                         11d28s22
        end if                                                          11d28s22
        iad=ivecso(isb)                                                 4d6s22
        id2e=ibcoff                                                     7d6s16
        ibcoff=id2e+nh                                                  7d6s16
        call enough('derofxor. 16',bc,ibc)
        do j=0,nh-1
         bc(id2e+j)=1d0/(bc(ieigs(isb)+j)**2)
        end do
        do i=0,nh-1                                                     7d6s16
         do j=0,nh-1                                                    7d6s16
          if(j.ne.i)then                                                7d6s16
           bot=1d0/(bc(id2e+j)-bc(id2e+i))                              7d6s16
           ji=itmp1+j+nh*i                                              7d6s16
           bc(ji)=-bc(ji)*bot                                            7d6s16
          end if                                                        7d6s16
         end do                                                         7d6s16
        end do                                                          7d6s16
        do i=0,nh-1                                                     7d6s16
         ii=itmp1+i*(nh+1)                                              7d6s16
         bc(id2e+i)=bc(ii)                                              7d6s16
         bc(ii)=0d0                                                     7d6s16
        end do                                                          7d6s16
        if(isk.gt.isb)then
         do i=0,nh-1
          do j=0,nk-1
           ji=ivdvsav+isob(isk)+j+nk*i                                  4d15s22
           ij=itmp4+i+nh*j                                              7d12s16
           bc(ij)=-bc(ji)
          end do
         end do
         if(nk.gt.0)then                                                11d28s22
          call dgemm('n','n',nh,nh,nk,1d0,bc(itmp4),nh,                  7d13s16
     $        bc(ivdvsav+isob(isk)),nk,0d0,bc(itmp3),nh,                4d15s22
     d' derofxor. 37')
         else                                                           11d28s22
          do iz=itmp3,itmp3+nh*nh-1                                     11d28s22
           bc(iz)=0d0                                                   11d28s22
          end do                                                        11d28s22
         end if                                                         11d28s22
        else
         do i=0,nk-1
          do j=0,nh-1
           ji=ivdvsav+isob(isb)+j+nh*i                                  4d15s22
           ij=itmp2+i+nk*j
           bc(ij)=-bc(ji)
          end do
         end do
         do i=0,nk*nh-1                                                 7d12s16
          bc(itmp4+i)=bc(ivdvsav+isob(isb)+i)                           4d15s22
         end do                                                         7d12s16
         if(min(nh,nk).gt.0)then                                        11d28s22
          call dgemm('n','n',nh,nh,nk,1d0,                               7d6s16
     $        bc(ivdvsav+isob(isb)),nh,bc(itmp2),nk,0d0,bc(itmp3),nh,   4d15s22
     d' derofxor. 38')
         else                                                           11d28s22
          do iz=itmp3,itmp3+nh*nh-1                                     11d28s22
           bc(iz)=0d0                                                   11d28s22
          end do                                                        11d28s22
         end if                                                         11d28s22
        end if                                                          7d6s16
        do i=0,nh*nh-1
         bc(itmp2+i)=bc(itmp1+i)+bc(itmp3+i)
        end do
        do i=0,nh-1
         do j=0,nh-1
          ji=itmp2+j+nh*i
          bc(ji)=bc(ji)*bc(ieigs(isb)+i)
         end do
        end do
        if(nh.gt.0)then                                                 11d28s22
         call dgemm('n','n',nh,nh,nh,1d0,bc(iad),nh,bc(itmp2),nh,0d0,    4d6s22
     $       bc(jxor2),nh,                                              4d6s22
     d' derofxor. 39')
         do i=0,nh-1                                                     7d8s16
          do j=0,nh-1                                                    7d8s16
           ji=jxor2+j+nh*i                                               4d6s22
           ij=itmp2+i+nh*j                                               7d8s16
           bc(ij)=bc(ji)                                                 7d8s16
          end do                                                         7d8s16
         end do                                                          7d8s16
         call dgemm('n','n',nh,nh,nh,1d0,bc(iad),nh,bc(itmp2),nh,0d0,    4d6s22
     $       bc(jxor2),nh,                                              4d6s22
     d' derofxor. 40')
        end if                                                          11d28s22
        do i=0,nh-1                                                     4d6s22
         do j=0,i                                                       7d8s16
          ji=jxor2+j+nh*i                                               4d6s22
          ij=jxor2+i+nh*j                                               4d6s22
          sum=bc(ji)+bc(ij)                                             7d8s16
          bc(ji)=sum                                                    7d8s16
          bc(ij)=sum                                                    7d8s16
         end do                                                         7d8s16
        end do                                                          7d8s16
         do i=0,nh-1                                                    4d4s22
          do j=0,nk-1
           ji=itmp2+j+nk*i
           ij=itmp4+i+nh*j                                              4d4s22
           bc(ji)=bc(ij)*bc(ieigs(isk)+j)
          end do
         end do
        if(min(nh,nk).gt.0)then                                         11d28s22
         call dgemm('n','n',nh,nh,nk,1d0,bc(itmp4),nh,bc(itmp2),nk,      4d4s22
     $        0d0,bc(itmp5),nh,                                         4d4s22
     d' derofxor. 41')
        else                                                            11d28s22
         do iz=itmp5,itmp5+nh*nh-1                                      11d28s22
          bc(iz)=0d0                                                    11d28s22
         end do                                                         11d28s22
        end if                                                          11d28s22
        if(nh.gt.0)then                                                 11d28s22
         call dgemm('n','n',nh,nh,nh,1d0,bc(iad),nh,bc(itmp5),nh,0d0,    4d6s22
     $       bc(itmp2),nh,                                              4d6s22
     d' derofxor. 42')
         do i=0,nh-1                                                     7d12s16
          do j=0,nh-1                                                    4d6s22
           ji=itmp2+j+nh*i                                               4d6s22
           ij=itmp5+i+nh*j                                               7d12s16
           bc(ij)=bc(ji)                                                 7d12s16
          end do                                                         7d12s16
         end do                                                          7d12s16
         call dgemm('n','n',nh,nh,nh,2d0,bc(iad),nh,bc(itmp5),nh,1d0,    4d6s22
     $       bc(jxor2),nh,                                              4d6s22
     d' derofxor. 43')
        end if                                                          11d28s22
        do i=0,nh-1                                                     7d8s16
         xmult=-0.5d0*(bc(ieigs(isb)+i)**3)*bc(id2e+i)                  7d8s16
         do j=0,nh-1                                                    4d6s22
          ij=itmp4+i+nh*j
          ji=iad+j+nh*i
          bc(ij)=bc(ji)*xmult
         end do                                                         7d8s16
        end do                                                          7d8s16
        if(nh.gt.0)then                                                 11d28s22
         call dgemm('n','n',nh,nh,nh,1d0,bc(iad),nh,bc(itmp4),nh,1d0,    4d6s22
     $       bc(jxor2),nh,                                              4d6s22
     d' derofxor. 44')
        end if                                                          11d28s22
        if(idwsdeb.gt.10)then
        write(6,*)('d2X: '),jxor2
        call prntm2(bc(jxor2),nh,nh,nh)                                 4d6s22
        end if
        ibcoff=itmp1                                                    7d8s16
        jxor2=jxor2+nh*nh                                               4d6s22
        ioffp=ioffp+nh*nh                                               7d5s16
        ioff=ioff+nhp*nhp                                               4d4s22
       end do                                                           7d5s16
      end if
      ibcoff=ibcoffo
      return
      end
