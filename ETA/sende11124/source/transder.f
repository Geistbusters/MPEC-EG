c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine transder(noc,morb,ixor,ixinv,nbasdwsc,nvertc,h0,iorb,
     $     itrans,ipuse,multh,idwsdeb,propmat,itran2,isob,bc,ibc)       11d10s22
      implicit real*8 (a-h,o-z)
c
c     compute transformation that takes integrals in mo basis to
c     integrals in the mo basis differentiating the orthogonality
c     transformation.
c
      include "common.store"
      include "common.hf"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension noc(*),morb(*),ixinv(*),itrans(8),ipt2(8),iptx(8),      7d21s16
     $     nbasdwsc(*),h0(*),iorb(*),multh(8,8),propmat(*),isob(*)      4d15s22
      dimension igap(60,2,2),sgap(60,2),itran2(8)                       7d27s16
      data (igap(i,1,1),i=1,60)/   1,   2,   1,   3,   2,   4,   3,   5,
     $     4,  6,   5,   7,   8,   9,   6,  10,   7,  11,  12,  13,   8,
     $     9,  14,  10, 11,  15,  16,  17,  12,  13,  18,  19,  14,  20,
     $     15,  21,  16,  22,  17,  23,  18,  24,  19,  20,  25,  26,
     $     27,  21,  28,  29,  22,  30,  31,  23,  32,  33,  24,  34,
     $     35,  25/
      data (igap(i,2,1),i=1,60) /  1,   1,   3,   1,   3,   1,   3,   1,
     $     3,  1,   3,   1,   1,   1,   3,   1,   3,   1,   1,   1,   3,
     $     3,  1,   3,   3,   1,   1,   1,   3,   3,   1,   1,   3,   1,
     $     3,  1,   3,   1,   3,   1,   3,   1,   3,   3,   1,   1,   1,
     $     3,  1,   1,   3,   1,   1,   3,   1,   1,   3,   1,   1,   3/
      data (sgap(i,1),i=1,60)/1d0, 1d0,-1d0, 1d0, 1d0, 1d0, 1d0, 1d0,
     $     -1d0, 1d0,-1d0,
     $    -1d0, 1d0,-1d0, 1d0, 1d0,-1d0, 1d0, 1d0, 1d0, 1d0, 1d0,-1d0,
     $     1d0, 1d0, 1d0,-1d0, 1d0, 1d0, 1d0, 1d0,-1d0,-1d0, 1d0,-1d0,
     $     1d0,-1d0, 1d0,-1d0, 1d0, 1d0, 1d0,-1d0, 1d0, 1d0, 1d0, 1d0,
     $    -1d0, 1d0,-1d0,-1d0, 1d0,-1d0,-1d0, 1d0,-1d0,-1d0, 1d0,-1d0,
     $    -1d0/
      data (igap(i,1,2),i=1,18)/   1,   1,   2,   3,   2,   4,   5,   3,
     $     6,   4,    5,   7,   8,   9,   6,  10,  11,   7/
      data (igap(i,2,2),i=1,18)/   2,   4,   2,   2,   4,   2,   2,   4,
     $     2,   4,    4,   2,   2,   2,   4,   2,   2,   4/
      data (sgap(i,2),i=1,18)/ 1d0, 1d0,-1d0,-1d0, 1d0, 1d0, 1d0, 1d0,
     $     -1d0, 1d0,
     $     -1d0, 1d0, 1d0, 1d0,-1d0, 1d0, 1d0,-1d0/
      dimension gsym(78),igsym(18)
      data gsym/78*1d0/
      data igsym/6,10,13,16,19,23,26,30,32,35,37,40,43,44,48,56,57,60/
      do i=1,18
       gsym(igsym(i))=-1d0
      end do
      if(idwsdeb.gt.10)write(6,*)('in transder ')
      ibcoffo=ibcoff
      if(idwsdeb.gt.10)write(6,*)('generating transformation ...')
      if(ipuse.eq.1)then
       if(idwsdeb.gt.10)write(6,*)('totally symmetric operator ')
       jxor=ixor
       ioff=1
       ioffh=1                                                          12d23s16
       nwds=0                                                           7d27s16
       do isb=1,nsymb                                                   7d27s16
        nn=nbasdwsc(isb)**2                                             7d27s16
        itrans(isb)=ibcoff                                              7d27s16
        ibcoff=ibcoff+nn                                                7d27s16
        nwds=nwds+nn                                                    7d27s16
       end do                                                           7d27s16
       jtran2=ibcoff                                                    7d27s16
       ibcoff=ibcoff+nwds                                               7d27s16
       call enough('transder.  1',bc,ibc)
       do isb=1,nsymb
        itran2(isb)=jtran2                                              7d27s16
        itmp1=ibcoff                                                    7d27s16
        ibcoffo=itmp1
        itmp2=itmp1+nbasdwsc(isb)**2
        ibcoff=itmp2+nbasdwsc(isb)**2
        itmp3=ibcoff
        ibcoff=itmp3+nbasdwsc(isb)**2
        itmp4=ibcoff                                                    7d20s16
        ibcoff=itmp4+nbasdwsc(isb)**2                                   7d20s16
        call enough('transder.  2',bc,ibc)
        call dgemm('n','n',nbasdwsc(isb),nbasdwsc(isb),nbasdwsc(isb),
     $      1d0,bc(jxor),nbasdwsc(isb),bc(morb(isb)),nbasdwsc(isb),0d0,       2d23s16
     $       bc(itmp1),nbasdwsc(isb),                                     2d23s16
     d' transder.  3')
        call dgemm('n','n',nbasdwsc(isb),nbasdwsc(isb),nbasdwsc(isb),   2d23s16
     $     1d0,bc(ixinv(isb)),nbasdwsc(isb),bc(itmp1),nbasdwsc(isb),0d0,     2d23s16
     $       bc(itmp2),nbasdwsc(isb),                                     2d23s16
     d' transder.  4')
        do i=0,nbasdwsc(isb)-1                                            2d23s16
         do j=0,nbasdwsc(isb)-1
          ij=i+nbasdwsc(isb)*j
          ji=j+nbasdwsc(isb)*i
          bc(itmp1+ij)=bc(morb(isb)+ji)
         end do
        end do
        call dgemm('n','n',nbasdwsc(isb),nbasdwsc(isb),nbasdwsc(isb),
     $      1d0,bc(itmp1),nbasdwsc(isb),bc(itmp2),nbasdwsc(isb),0d0,
     $       bc(itrans(isb)),nbasdwsc(isb),
     d' transder.  5')
        jxor2=jxor+nbasdwsc(isb)**2                                     7d20s16
          call dgemm('n','n',nbasdwsc(isb),nbasdwsc(isb),nbasdwsc(isb),
     $       1d0,bc(ixinv(isb)),nbasdwsc(isb),
     $        bc(jxor2),nbasdwsc(isb),0d0,bc(ibcoff),nbasdwsc(isb))
        call dgemm('n','n',nbasdwsc(isb),nbasdwsc(isb),nbasdwsc(isb),   7d20s16
     $      1d0,bc(jxor2),nbasdwsc(isb),bc(morb(isb)),nbasdwsc(isb),0d0,7d20s16
     $      bc(itmp4),nbasdwsc(isb),                                    7d20s16
     d' transder.  6')
        call dgemm('n','n',nbasdwsc(isb),nbasdwsc(isb),nbasdwsc(isb),   2d23s16
     $     1d0,bc(ixinv(isb)),nbasdwsc(isb),bc(itmp4),nbasdwsc(isb),0d0,7d20s16
     $      bc(itmp2),nbasdwsc(isb),                                    7d20s16
     d' transder.  7')
        call dgemm('n','n',nbasdwsc(isb),nbasdwsc(isb),nbasdwsc(isb),   7d20s16
     $      1d0,bc(itmp1),nbasdwsc(isb),bc(itmp2),nbasdwsc(isb),0d0,    7d20s16
     $       bc(jtran2),nbasdwsc(isb),                                  7d20s16
     d' transder.  8')
        if(idwsdeb.gt.10)then
         write(6,*)('trans: ')
         call prntm2(bc(itrans(isb)),nbasdwsc(isb),nbasdwsc(isb),       12d23s16
     $        nbasdwsc(isb))
         if(nsymb.eq.1)call printa(bc(itrans(isb)),nbasdwsc,0,1,0,
     $        nbasdwsc,0,1,0,bc(ibcoff))
         write(6,*)('tran2:a '),jtran2
         call prntm2(bc(jtran2),nbasdwsc(isb),nbasdwsc(isb),
     $        nbasdwsc(isb))
         if(nsymb.eq.1)call printa(bc(jtran2),nbasdwsc,0,1,0,
     $        nbasdwsc,0,1,0,bc(ibcoff))
        end if
        ioffp=ioff+nwds                                                 12d23s16
        itry=ioffp+8*(nbasdwsc(isb)+1)
        if(idwsdeb.gt.10)then
         write(6,*)('starting propmat for 1st der: '),ioff,
     $        loc(propmat(ioff))
         call prntm2(propmat(ioff),nbasdwsc(isb),nbasdwsc(isb),
     $        nbasdwsc(isb))
         if(nsymb.eq.1)call printa(propmat(ioff),nbasdwsc,0,1,0,
     $        nbasdwsc,0,1,0,bc(ibcoff))
         write(6,*)('starting propmat for 2nd der: ')
         call prntm2(propmat(ioffp),nbasdwsc(isb),nbasdwsc(isb),
     $        nbasdwsc(isb))
         if(nsymb.eq.1)call printa(propmat(ioffp),nbasdwsc,0,1,0,
     $        nbasdwsc,0,1,0,bc(ibcoff))
        end if                                                          7d21s16
        itest=ibcoff
        ibcoff=itest+nbasdwsc(isb)**2
        call enough('transder.  3',bc,ibc)
        call dgemm('n','n',nbasdwsc(isb),nbasdwsc(isb),nbasdwsc(isb),
     $        2d0,propmat(ioff),nbasdwsc(isb),bc(itrans(isb)),
     $       nbasdwsc(isb),0d0,bc(itest),nbasdwsc(isb),
     d' transder.  9')
        ibcoff=itest                                                    5d4s22
        call dgemm('n','n',nbasdwsc(isb),nbasdwsc(isb),nbasdwsc(isb),
     $        2d0,propmat(ioff),nbasdwsc(isb),bc(itrans(isb)),
     $       nbasdwsc(isb),1d0,propmat(ioffp),nbasdwsc(isb),
     d' transder.  9')
        if(idwsdeb.gt.10)then
         write(6,*)('after D(1)ao * trans contribution to propmat2: ')
         call prntm2(propmat(ioffp),nbasdwsc(isb),nbasdwsc(isb),
     $        nbasdwsc(isb))
         if(nsymb.eq.1)call printa(propmat(ioffp),nbasdwsc,0,1,0,
     $        nbasdwsc,0,1,0,bc(ibcoff))
        end if
        ibcoff=itmp1
c
c     second contribution to (n|d/dq|m): the der of the orthogonality matrix
c     but since (n|m) is the unit matrix, this is just the transformation
c     matrix.
c
        do i=0,nbasdwsc(isb)-1                                          3d29s16
         do j=0,nbasdwsc(isb)-1                                         3d29s16
          ji=j+nbasdwsc(isb)*i                                          3d29s16
          propmat(ioff+ji)=propmat(ioff+ji)+bc(itrans(isb)+ji)          3d29s16
          propmat(ioffp+ji)=propmat(ioffp+ji)+bc(jtran2+ji)             7d21s16
         end do                                                         3d29s16
        end do                                                          3d29s16
        if(idwsdeb.gt.10)then
         write(6,*)('final propmat after 2nd step '),ioff,
     $        loc(propmat(ioff))
         call prntm2(propmat(ioff),nbasdwsc(isb),nbasdwsc(isb),
     $        nbasdwsc(isb))
         if(nsymb.eq.1)then
          write(6,*)('times 2 ...')
          i2x=ibcoff
          ibcoff=i2x+nbasdwsc(isb)**2
          call enough('transder.  5',bc,ibc)
          do i=0,nbasdwsc(isb)**2-1
           bc(i2x+i)=propmat(ioff+i)*2d0
          end do
          ibcoff=i2x
         end if
         write(6,*)('final propmat2 after 2nd step '),ioffp,
     $        loc(propmat(ioffp))
         call prntm2(propmat(ioffp),nbasdwsc(isb),nbasdwsc(isb),
     $        nbasdwsc(isb))
         if(nsymb.eq.1)then
          call printa(propmat(ioffp),nbasdwsc,0,1,0,nbasdwsc,0,1,0,
     $         bc(ibcoff))
         end if
        end if                                                          7d21s16
        jxor=jxor+nbasdwsc(isb)*nbasdwsc(isb)*2                         6d22s16
        ioff=ioff+nbasdwsc(isb)*nbasdwsc(isb)                           12d23s16
        ioffh=ioffh+nbasdwsc(isb)*nbasdwsc(isb)                         12d23s16
        ibcoff=ibcoffo
        jtran2=jtran2+nbasdwsc(isb)*nbasdwsc(isb)                       7d27s16
       end do
      else                                                              3d15s16
       if(idwsdeb.gt.10)
     $      write(6,*)('derivative operator has symmetry '),ipuse
       n2ndz=0                                                          7d21s16
       nx=0                                                             7d21s16
       do isb=1,nsymb                                                   7d21s16
        isk=multh(isb,ipuse)
        nb=nbasdwsc(isb)                                                7d21s16
        ipt2(isb)=nx                                                    7d21s16
        nx=nx+nb*nb                                                     7d21s16
        if(isk.gt.isb)then                                              7d21s16
         nk=nbasdwsc(isk)                                               7d21s16
         iptx(isb)=n2ndz                                                7d21s16
         n2ndz=n2ndz+nb*nk                                              7d21s16
         iptx(isk)=n2ndz                                                7d21s16
         n2ndz=n2ndz+nb*nk                                              7d21s16
        end if                                                          7d21s16
       end do                                                           7d21s16
       itran1=ibcoff                                                    7d21s16
       ktran2=itran1+n2ndz                                              7d27s16
       ibcoff=ktran2+nx                                                 7d27s16
       n2ndzh=n2ndz/2                                                   7d21s16
       call enough('transder.  7',bc,ibc)
       jxor=ixor
       ioffp=1                                                          3d29s16
       do isb=1,nsymb
        if(idwsdeb.gt.10)write(6,*)('for isb = '),isb,iptx(isb)
        nb=nbasdwsc(isb)
        isk=multh(isb,ipuse)
        if(idwsdeb.gt.10)write(6,*)('isk '),isk
        if(isk.lt.isb)then                                              3d29s16
         ioffp=ioffp+nbasdwsc(isk)*nbasdwsc(isb)                        3d29s16
        end if                                                          3d29s16
        if(isk.gt.isb)then
         nk=nbasdwsc(isk)                                               3d15s16
         itrans(isb)=itran1+iptx(isb)                                   7d21s16
         itmp1=ibcoff                                                   7d21s16
         itmp2=itmp1+max(nb,nk)*nb                                      3d15s16
         ibcoff=itmp2+nb*nk                                             3d15s16
         call enough('transder.  8',bc,ibc)
         if(min(nb,nk).gt.0)then                                        11d28s22
          call dgemm('n','n',nb,nk,nk,1d0,bc(jxor),nb,bc(morb(isk)),nk,  3d15s16
     $        0d0,bc(itmp1),nb,                                         3d15s16
     d' transder. 13')
          call dgemm('n','n',nb,nk,nb,1d0,bc(ixinv(isb)),nb,bc(itmp1),  11d28s22
     $         nb,0d0,bc(itmp2),nb,                                     11d28s22
     d' transder. 14')
          do i=0,nb-1                                                    3d15s16
           do j=0,nb-1                                                   3d15s16
            ij=i+nb*j                                                    3d15s16
            ji=j+nb*i                                                    3d15s16
            bc(itmp1+ij)=bc(morb(isb)+ji)                                3d15s16
           end do                                                        3d15s16
          end do                                                         3d15s16
          call dgemm('n','n',nb,nk,nb,1d0,bc(itmp1),nb,bc(itmp2),nb,0d0, 3d15s16
     $        bc(itrans(isb)),nb,                                       3d15s16
     d' transder. 15')
         else                                                           11d28s22
          do iz=itrans(isb),itrans(isb)+nb*nk-1                         11d28s22
           bc(iz)=0d0                                                   11d28s22
          end do                                                        11d28s22
         end if                                                         11d28s22
         ioffp=isob(isb)+1                                              4d15s22
         if(idwsdeb.gt.10)then
          write(6,*)('starting propmat: '),ioffp,isb,isk
          call prntm2(propmat(ioffp),nb,nk,nb)
          write(6,*)('trans for bra: ')
          call prntm2(bc(itrans(isb)),nb,nk,nb)
         end if
         ibcoff=itmp1                                                   7d21s16
         itrans(isk)=itran1+iptx(isk)                                   7d21s16
         itmp1=ibcoff                                                   7d21s16
         itmp2=itmp1+max(nb,nk)*nk                                      3d15s16
         ibcoff=itmp2+nb*nk                                             3d15s16
         call enough('transder.  9',bc,ibc)
         do i=0,nb-1
          do j=0,nk-1
           ji=itmp1+j+nk*i
           ij=jxor+i+nb*j
           bc(ji)=bc(ij)
          end do
         end do
         if(min(nk,nb).gt.0)then                                        11d28s22
          call dgemm('n','n',nk,nb,nb,1d0,bc(itmp1),nk,bc(morb(isb)),nb,  3d15s16
     $        0d0,bc(itrans(isk)),nk,                                   3d15s16
     d' transder. 16')
          call dgemm('n','n',nk,nb,nk,1d0,bc(ixinv(isk)),nk,             3d15s16
     $        bc(itrans(isk)),nk,0d0,bc(itmp2),nk,                      3d15s16
     d' transder. 17')
          do i=0,nk-1                                                    3d15s16
           do j=0,nk-1                                                   3d15s16
            ij=i+nk*j                                                    3d15s16
            ji=j+nk*i                                                    3d15s16
            bc(itmp1+ij)=bc(morb(isk)+ji)                                3d15s16
           end do                                                        3d15s16
          end do                                                         3d15s16
          call dgemm('n','n',nk,nb,nk,1d0,bc(itmp1),nk,bc(itmp2),nk,0d0, 3d15s16
     $        bc(itrans(isk)),nk,                                       3d15s16
     d' transder. 18')
         end if                                                         11d28s22
         ibcoff=itmp1                                                   3d15s16
         ioffb=ipt2(isb)+n2ndz+1                                        7d21s16
         if(idwsdeb.gt.10)then
          write(6,*)('transformation matrix for ket: ')                          3d15s16
          call prntm2(bc(itrans(isk)),nk,nb,nk)                          3d15s16
          write(6,*)('propmat for bra ')
          call prntm2(propmat(ioffp),nb,nk,nb)
         end if
c
         if(idwsdeb.gt.10)then
          write(6,*)('propmat2 for symmetry '),isb,ioffb
          call prntm2(propmat(ioffb),nb,nb,nb)
         end if
         if(min(nb,nk).gt.0)then                                        11d28s22
          call dgemm('n','n',nb,nb,nk,2d0,propmat(ioffp),nb,             7d21s16
     $        bc(itrans(isk)),nk,1d0,propmat(ioffb),nb,                 7d21s16
     d' transder. 19')
         end if                                                         11d28s22
         ioffk=ipt2(isk)+n2ndz+1
         if(idwsdeb.gt.10)then
          write(6,*)('with propmat contribution: 2*Tk*Pb ')
          write(6,*)('propmat2 for symmetry '),isk,ioffk,
     $        loc(propmat(ioffk))
          call prntm2(propmat(ioffk),nk,nk,nk)
          write(6,*)('propmat for other symmetry (ket): '),
     $         loc(propmat(isob(isk)+1)),isob(isk)+1                     4d15s22
          call prntm2(propmat(isob(isk)+1),nk,nb,nk)                     4d15s22
         end if                                                         7d21s16
         if(idwsdeb.gt.10)then                                          7d21s16
          write(6,*)('trans again (bra)')
          call prntm2(bc(itrans(isb)),nb,nk,nb)
         end if                                                         7d21s16
         if(min(nk,nb).gt.0)then                                        11d28s22
          call dgemm('n','n',nk,nk,nb,2d0,propmat(isob(isk)+1),nk,       4d15s22
     $        bc(itrans(isb)),
     $        nb,1d0,propmat(ioffk),nk,                                 7d21s16
     d' transder. 20')
         end if                                                         11d28s22
         if(idwsdeb.gt.10)then
          write(6,*)('with other propmat contribution: 2Pk*Tb ')
          call prntm2(propmat(ioffk),nk,nk,nk)
         end if
         itmp1=ibcoff                                                   7d21s16
         itmp2=itmp1+max(nb,nk)**2                                      7d21s16
         ibcoff=itmp2+max(nb,nk)**2                                     7d21s16
         call enough('transder. 10',bc,ibc)
         ixor2=ixor+n2ndz+ipt2(isb)                                     4d18s22
         if(idwsdeb.gt.10)then
          write(6,*)('d2X for symmetry '),isb,ixor2,ipt2(isb)
          call prntm2(bc(ixor2),nb,nb,nb)
         end if
         if(nb.gt.0)then                                                11d28s22
          call dgemm('n','n',nb,nb,nb,1d0,bc(ixinv(isb)),nb,
     $        bc(ixor2),nb,0d0,bc(ibcoff),nb)
          call dgemm('n','n',nb,nb,nb,1d0,bc(ixor2),nb,bc(morb(isb)),nb, 7d21s16
     $        0d0,bc(itmp1),nb,                                         7d21s16
     d' transder. 21')
          call dgemm('n','n',nb,nb,nb,1d0,bc(ixinv(isb)),nb,bc(itmp1),  11d28s22
     $         nb,0d0,bc(itmp2),nb,                                     11d28s22
     d' transder. 22')
          do i=0,nb-1                                                    7d21s16
           do j=0,nb-1                                                   7d21s16
            ji=morb(isb)+j+nb*i                                          7d21s16
            ij=itmp1+i+nb*j                                              7d21s16
            bc(ij)=bc(ji)                                                7d21s16
           end do                                                        7d21s16
          end do                                                         7d21s16
          jtran2=ktran2+ipt2(isb)                                        7d21s16
          itran2(isb)=jtran2                                             7d27s16
          call dgemm('n','n',nb,nb,nb,1d0,bc(itmp1),nb,bc(itmp2),nb,0d0, 7d21s16
     $        bc(jtran2),nb,                                            7d21s16
     d' transder. 23')
         end if                                                         11d28s22
         if(idwsdeb.gt.10)then
          write(6,*)('tran2 bra:MT*Xinv*xor2*M ')
          call prntm2(bc(jtran2),nb,nb,nb)                               7d21s16
         end if
         do i=0,nb*nb-1                                                 7d21s16
          propmat(ioffb+i)=propmat(ioffb+i)+bc(jtran2+i)                7d21s16
         end do                                                         7d21s16
         ixor2=ixor+n2ndz+ipt2(isk)                                     4d19s22
         if(idwsdeb.gt.10)then
          write(6,*)('bra propmat2 with tran2 contribution: ')               7d21s16
          call prntm2(propmat(ioffb),nb,nb,nb)                           7d21s16
          write(6,*)('d2X for symmetry '),isk
          call prntm2(bc(ixor2),nk,nk,nk)
         end if
         if(nk.gt.0)then                                                11d28s22
          call dgemm('n','n',nk,nk,nk,1d0,bc(ixor2),nk,bc(morb(isk)),nk, 7d21s16
     $        0d0,bc(itmp1),nk,                                         7d21s16
     d' transder. 24')
          call dgemm('n','n',nk,nk,nk,1d0,bc(ixinv(isk)),nk,bc(itmp1),  11d28s22
     $         nk,0d0,bc(itmp2),nk,                                     11d28s22
     d' transder. 25')
          do i=0,nk-1                                                    7d21s16
           do j=0,nk-1                                                   7d21s16
            ji=morb(isk)+j+nk*i                                          7d21s16
            ij=itmp1+i+nk*j                                              7d21s16
            bc(ij)=bc(ji)                                                7d21s16
           end do                                                        7d21s16
          end do                                                         7d21s16
          jtran2=ktran2+ipt2(isk)                                        7d27s16
          itran2(isk)=jtran2                                             11d1s16
          call dgemm('n','n',nk,nk,nk,1d0,bc(itmp1),nk,bc(itmp2),nk,0d0, 7d21s16
     $        bc(jtran2),nk,                                            7d21s16
     d' transder. 26')
         end if                                                         11d28s22
         do i=0,nk*nk-1                                                 7d21s16
          propmat(ioffk+i)=propmat(ioffk+i)+bc(jtran2+i)                7d21s16
         end do                                                         7d21s16
         if(idwsdeb.gt.10)then
          write(6,*)('tran2 ket: '),jtran2                                          7d21s16
          call prntm2(bc(jtran2),nk,nk,nk)                               7d21s16
          write(6,*)('final ket propmat2 with tran2 contribution: ')
          call prntm2(propmat(ioffk),nk,nk,nk)                           7d21s16
         end if
c
c     second part of (n|d/dq|m) is der of orthogonality matrix, ie      3d29s16
c     the transformation matrix since (n|m) is the unit matrix.
c
         do i=0,nk-1                                                    3d29s16
          do j=0,nb-1                                                   3d29s16
           ji=j+nb*i                                                    3d29s16
           propmat(ioffp+ji)=propmat(ioffp+ji)+bc(itrans(isb)+ji)       3d29s16
          end do                                                        3d29s16
         end do                                                         3d29s16
         if(idwsdeb.gt.10)then
          write(6,*)('final propmat for sym '),isb,isk,ioffp,
     $         loc(propmat(ioffp))
          call prntm2(propmat(ioffp),nb,nk,nb)
         end if
         ioffp=ioffp+nb*nk                                              3d29s16
         jxor=jxor+nb*nk                                                3d15s16
        end if                                                          3d15s16
       end do                                                           3d15s16
      end if                                                            3d15s16
      return
      end
