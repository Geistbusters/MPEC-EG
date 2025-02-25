c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine transder1(morb,ixor,ixinv,nbasdwsc,iorb,               5d4s22
     $     itrans,ipuse,multh,idwsdeb,bc,ibc)                           11d10s22
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
      dimension morb(*),ixinv(*),itrans(8),ipt2(8),iptx(8),             5d4s22
     $     nbasdwsc(*),iorb(*),multh(8,8)                               6d8s22
      if(idwsdeb.gt.10)write(6,*)('Hi, my name is transder1!')
      ibcoffo=ibcoff
      if(idwsdeb.gt.10)write(6,*)('generating transformation ...')
      if(ipuse.eq.1)then
       if(idwsdeb.gt.10)write(6,*)('totally symmetric operator ')
       jxor=ixor
       ioff=1
       ioffh=1                                                          12d23s16
       do isb=1,nsymb                                                   7d27s16
        nn=nbasdwsc(isb)**2                                             7d27s16
        itrans(isb)=ibcoff                                              7d27s16
        ibcoff=ibcoff+nn                                                7d27s16
       end do                                                           7d27s16
       call enough('transder1.  1',bc,ibc)
       do isb=1,nsymb
        if(nbasdwsc(isb).gt.0)then                                      5d20s22
         itmp1=ibcoff                                                    7d27s16
         ibcoffo=itmp1
         itmp2=itmp1+nbasdwsc(isb)**2
         ibcoff=itmp2+nbasdwsc(isb)**2
         itmp3=ibcoff
         ibcoff=itmp3+nbasdwsc(isb)**2
         itmp4=ibcoff                                                    7d20s16
         ibcoff=itmp4+nbasdwsc(isb)**2                                   7d20s16
         call enough('transder1.  2',bc,ibc)
         call dgemm('n','n',nbasdwsc(isb),nbasdwsc(isb),nbasdwsc(isb),
     $      1d0,bc(jxor),nbasdwsc(isb),bc(morb(isb)),nbasdwsc(isb),0d0,       2d23s16
     $       bc(itmp1),nbasdwsc(isb),                                     2d23s16
     d' transder.  3')
         if(idwsdeb.gt.10)then
          write(6,*)('jxor ')
          call prntm2(bc(jxor),nbasdwsc(isb),nbasdwsc(isb),
     $         nbasdwsc(isb))
          write(6,*)('jxor ')
          call prntm2(bc(morb(isb)),nbasdwsc(isb),nbasdwsc(isb),
     $         nbasdwsc(isb))
          write(6,*)('jxor * morb ')
          call prntm2(bc(itmp1),nbasdwsc,nbasdwsc,nbasdwsc)
         end if
         call dgemm('n','n',nbasdwsc(isb),nbasdwsc(isb),nbasdwsc(isb),   2d23s16
     $     1d0,bc(ixinv(isb)),nbasdwsc(isb),bc(itmp1),nbasdwsc(isb),0d0,     2d23s16
     $       bc(itmp2),nbasdwsc(isb),                                     2d23s16
     d' transder.  4')
         if(idwsdeb.gt.10)then
          write(6,*)('X-1 jxor * morb')
          call prntm2(bc(itmp2),nbasdwsc,nbasdwsc,nbasdwsc)
         end if
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
         if(idwsdeb.gt.10)then
          write(6,*)('trans: ')
          call prntm2(bc(itrans(isb)),nbasdwsc(isb),nbasdwsc(isb),       12d23s16
     $        nbasdwsc(isb))
          if(nsymb.eq.1)call printa(bc(itrans(isb)),nbasdwsc(isb),0,
     $         1,0,nbasdwsc(isb),0,1,0,bc(ibcoff))
         end if
         ibcoff=itmp1
        end if                                                          5d20s22
        jxor=jxor+nbasdwsc(isb)*nbasdwsc(isb)                           5d9s22
        ioff=ioff+nbasdwsc(isb)*nbasdwsc(isb)*2                         6d3s22
        ioffh=ioffh+nbasdwsc(isb)*nbasdwsc(isb)                         12d23s16
        ibcoff=ibcoffo
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
       call enough('transder1.  3',bc,ibc)
       jxor=ixor
       do isb=1,nsymb
        if(idwsdeb.gt.10)write(6,*)('for isb = '),isb,iptx(isb)
        nb=nbasdwsc(isb)
        isk=multh(isb,ipuse)
        if(idwsdeb.gt.10)write(6,*)('isk '),isk
        if(isk.gt.isb)then
         nk=nbasdwsc(isk)                                               3d15s16
         itrans(isb)=itran1+iptx(isb)                                   7d21s16
         itmp1=ibcoff                                                   7d21s16
         itmp2=itmp1+max(nb,nk)*nb                                      3d15s16
         ibcoff=itmp2+nb*nk                                             3d15s16
         call enough('transder1.  4',bc,ibc)
         if(min(nb,nk).gt.0)then                                        7d11s22
          if(idwsdeb.ne.0)then                                          7d14s22
           write(6,*)('jxor ')
           call prntm2(bc(jxor),nb,nk,nb)
           write(6,*)('morb '),isk
           call prntm2(bc(morb(isk)),nk,nk,nk)
          end if                                                        7d14s22
          call dgemm('n','n',nb,nk,nk,1d0,bc(jxor),nb,bc(morb(isk)),nk,  3d15s16
     $        0d0,bc(itmp1),nb,                                         3d15s16
     d' transder. 13')
          if(idwsdeb.ne.0)then                                          7d14s22
           write(6,*)('xinv '),isb
           call prntm2(bc(ixinv(isb)),nb,nb,nb)
           write(6,*)('itmp1 ')
           call prntm2(bc(itmp1),nb,nk,nb)
          end if                                                        7d14s22
          call dgemm('n','n',nb,nk,nb,1d0,bc(ixinv(isb)),nb,bc(itmp1),  7d11s22
     $         nb,0d0,bc(itmp2),nb,                                     7d11s22
     d' transder. 14')
          do i=0,nb-1                                                    3d15s16
           do j=0,nb-1                                                   3d15s16
            ij=i+nb*j                                                    3d15s16
            ji=j+nb*i                                                    3d15s16
            bc(itmp1+ij)=bc(morb(isb)+ji)                                3d15s16
           end do                                                        3d15s16
          end do                                                         3d15s16
          if(idwsdeb.ne.0)then                                          7d14s22
           write(6,*)('tmp2')
           call prntm2(bc(itmp2),nb,nk,nb)
           call prntm2(bc(itmp1),nb,nb,nb)
          end if                                                        7d14s22
          call dgemm('n','n',nb,nk,nb,1d0,bc(itmp1),nb,bc(itmp2),nb,0d0, 3d15s16
     $        bc(itrans(isb)),nb,                                       3d15s16
     d' transder. 15')
          ibcoff=itmp1                                                   7d21s16
          itrans(isk)=itran1+iptx(isk)                                   7d21s16
          itmp1=ibcoff                                                   7d21s16
          itmp2=itmp1+max(nb,nk)*nk                                      3d15s16
          ibcoff=itmp2+nb*nk                                             3d15s16
          call enough('transder1.  5',bc,ibc)
          do i=0,nb-1
           do j=0,nk-1
            ji=itmp1+j+nk*i
            ij=jxor+i+nb*j
            bc(ji)=bc(ij)
           end do
          end do
          call dgemm('n','n',nk,nb,nb,1d0,bc(itmp1),nk,bc(morb(isb)),nb,  3d15s16
     $        0d0,bc(itrans(isk)),nk,                                   3d15s16
     d' transder. 16')
          call dgemm('n','n',nk,nb,nk,1d0,bc(ixinv(isk)),nk,             3d15s16
     $        bc(itrans(isk)),nk,0d0,bc(itmp2),nk,                      3d15s16
     d' transder. 17')
          if(idwsdeb.ne.0)then                                          7d14s22
           write(6,*)('xinv*trans ')
           call prntm2(bc(itmp2),nk,nb,nk)
          end if                                                        7d14s22
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
          ibcoff=itmp1                                                   3d15s16
          if(idwsdeb.gt.10)then
           write(6,*)('transformation matrix for ket: ')                          3d15s16
           call prntm2(bc(itrans(isk)),nk,nb,nk)                          3d15s16
          end if
         end if                                                         7d11s22
         itmp1=ibcoff                                                   7d21s16
         itmp2=itmp1+max(nb,nk)**2                                      7d21s16
         ibcoff=itmp2+max(nb,nk)**2                                     7d21s16
         call enough('transder1.  6',bc,ibc)
c
c     second part of (n|d/dq|m) is der of orthogonality matrix, ie      3d29s16
c     the transformation matrix since (n|m) is the unit matrix.
c
         jxor=jxor+nb*nk                                                3d15s16
        end if                                                          3d15s16
       end do                                                           3d15s16
      end if                                                            3d15s16
      return
      end
