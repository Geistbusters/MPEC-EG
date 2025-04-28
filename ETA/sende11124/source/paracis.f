c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine paracis(idorr,ioooo,jmats,kmats,h0mo,noc,nvirtc,
     $     nbasdwsc,icore,multh,ieigv,ehf,ipropmat,natom,iapair,bc,ibc) 11d9s22
      implicit real*8 (a-h,o-z)
      integer*8 ibk8,nproc8,ileft8,itrans8,ipropmat
c
c     ci singles. we can either diagonalize to get natural orbitals
c     or do response for non-adiabatic effects.
c     idorr=+1 means diagonalize
c     idorr=-1 means response
c
      include "common.store"
      include "common.hf"
      character*1 cart(3)                                               4d28s22
      dimension ioooo(1),jmats(1),kmats(1),h0mo(1),noca(8),nvirtc(1),
     $ nbasdwsc(1),icore(1),noc(8),i4o(idbk),iph(8),iham(idbk),
     $     multh(8,8),ieigv(8),ibk(8),ihcis(idbk),icbk(5,idbk),
     $     ipropmat(1),iapair(3,*),got(2),igot(2),ipt2(8,8)             4d29s22
      equivalence (itrans8,trans8)
      data cart/'x','y','z'/                                            4d28s22
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      common/drsigncm/drsign                                            8d20s24
      if(mynowprog.eq.0)write(6,*)('in paracis'),idorr
      iawgt=ibcoff                                                      2d22s23
      ibcoff=iawgt+natom                                                2d22s23
      call enough('paracis.awgt',bc,ibc)                                2d22s23
      do ia=1,natom                                                     2d22s23
       iam=ia-1                                                         2d22s23
       bc(iawgt+iam)=1d0                                                2d22s23
       if(iapair(1,ia).ne.0)then                                        2d22s23
        bc(iawgt+iam)=0.25d0                                            2d22s23
       end if                                                           2d22s23
      end do                                                            2d22s23
      if(nsymb.gt.1)then                                                4d28s22
       ic2mat=ibcoff                                                     4d28s22
       ic1mat=ic2mat+(3*natom)**2                                        4d28s22
       ic0mat=ic1mat+3*(natom**2)                                       3d3s23
       ibcoff=ic0mat+natom*natom                                        3d3s23
       call enough('paracis.  1',bc,ibc)
       do iz=ic2mat,ibcoff-1                                            4d28s22
        bc(iz)=0d0                                                      4d28s22
       end do                                                           4d28s22
      end if                                                            4d28s22
      kgoal=2471824
      n2sz=0                                                            12d12s16
      do isb=1,nsymb
       noca(isb)=noc(isb)-icore(isb)
       n2sz=n2sz+nbasdwsc(isb)*nbasdwsc(isb)                            12d12s16
      end do
      if(mynowprog.eq.0)write(6,1)(icore(isb),isb=1,nsymb)               3d3s17
    1 format('no. of core orbitals    ',8i5)
      if(mynowprog.eq.0)write(6,2)(noc(isb),isb=1,nsymb)
    2 format('no. of valence orbitals ',8i5)
      if(mynowprog.eq.0)write(6,3)(noca(isb),isb=1,nsymb)
    3 format('no. of active orbitals  ',8i5)
      if(mynowprog.eq.0)write(6,4)(nvirt(isb),isb=1,nsymb)
    4 format('no. of virtual orbs     ',8i5)
      ibcoffo=ibcoff
      if(idorr.lt.0)then
       natom3=natom*3
       natom4=natom3+natom                                              6d20s16
       do ixyz=1,3                                                      4d28s22
        igot(1)=0                                                       4d28s22
        igot(2)=0                                                       4d28s22
        ngot=0                                                          4d28s22
        do ia=1,natom                                                   4d28s22
         ift=ixyz+3*(ia-1)
         ift1=ift+natom4                                                6d20s16
         ift2=ift1+natom3
         itrans8=ipropmat(ift2)
         if(mynowprog.eq.0)then                                         4d28s22
          if(iapair(1,ia).gt.0)then                                      4d28s22
           igot(1)=ift                                                  4d28s22
           got(1)=trans8                                                 4d28s22
           ngot=ngot+1                                                   4d28s22
           iause=ia                                                     4d28s22
          else if(iapair(1,ia).lt.0)then                                 4d28s22
           igot(2)=ift                                                  4d28s22
           got(2)=trans8                                                 4d28s22
           ngot=ngot+1                                                   4d28s22
          else                                                          4d28s22
           write(6,253)ift,cart(ixyz),cart(ixyz),                        4d28s22
     $        ia,trans8                                                 4d28s22
          end if                                                         4d28s22
          if(ngot.eq.2)then                                              4d28s22
           got1=0.5d0*(got(1)+got(2))                                   4d28s22
           got2=0.5d0*(got(1)-got(2))                                   4d28s22
           write(6,254)igot,cart(ixyz),cart(ixyz),iause,got1            4d28s22
  254      format('fcn ',i2,'+',i2,x,2a1,' nucleus ',i2,' hf der ',      4d19s23
     $          es15.7)                                                 4d28s22
           write(6,255)igot,cart(ixyz),cart(ixyz),iapair(1,iause),got2  4d28s22
  255      format('fcn ',i2,'-',i2,x,2a1,' nucleus ',i2,' hf der ',      4d19s23
     $          es15.7)                                                 4d28s22
          end if                                                         4d28s22
         end if                                                         4d28s22
  253    format('fcn ',i3,': ',2a1,' nucleus ',i2,' hf der ',es15.7)     4d19s23
         ioff=ipropmat(ift)
         do isb=1,nsymb
          isbo=multh(isb,ipropmat(ift1))
          if(isbo.ge.isb)then                                           4d28s22
           ipt2(isb,isbo)=ioff-ipropmat(ift)                            4d29s22
           ioff=ioff+nbasdwsc(isb)*nbasdwsc(isbo)                       4d28s22
           if(isb.ne.isbo)then                                          4d28s22
            ipt2(isbo,isb)=ioff-ipropmat(ift)                           4d29s22
            ioff=ioff+nbasdwsc(isb)*nbasdwsc(isbo)                       4d28s22
           end if                                                       4d28s22
          end if
         end do
         do isb=1,nsymb
          ioff=ioff+nbasdwsc(isb)*nbasdwsc(isb)
         end do
        end do
       end do
      end if
c
c     (na|H|mb)=delta_ab delta_nm (Ehf+Fbb-Fmm)-J_nm^ab+K_nm^ab+K_mn^ab ?
c     (na|H|mb)=delta_ab delta_nm (Ehf+Fbb-Fmm)-J_nm^ab+2K_mn^ab ?
c     On my proc, compute my contribution to H from my J and K.
c     Then do H*guess.
c
      do isw=1,nsymb
       if(mynowprog.eq.0)write(6,*)('functions of symmetry '),isw
       igoal=0
       igoalc=0
       kgoal=2652328
       if(nsymb.eq.1)then
        igoal=3
        igoalc=11
       else
        if(isw.eq.2)then
         igoal=1
         igoalc=10
        end if
       end if
       nbk=0
       ntot=0
       do isv=1,nsymb
        iso=multh(isv,isw)
        nsz=nvirt(isv)*noca(iso)
        if(nsz.gt.0)then
         ntot=ntot+nsz
         nbk=nbk+1
         ibk(nbk)=isv
         if(mynowprog.eq.0)write(6,5)nbk,iso,isv,nsz
    5    format('block ',i2,' with active symmetry ',i1,
     $        ' and virtual symmetry ',i1,' has size ',3i5)
        end if
       end do
       if(mynowprog.eq.0)write(6,*)('total number of fcns = '),ntot
       istoh=0
       do ib1=1,nbk
        iv1=ibk(ib1)
        io1=multh(iv1,isw)
        do ib2=1,ib1
         iv2=ibk(ib2)
         io2=multh(iv2,isw)
         do is=1,nsdlk
          if(isblk(1,is).eq.io1.and.isblk(2,is).eq.io2.and.
     $         isblk(3,is).eq.iv1.and.isblk(4,is).eq.iv2)then
           jcase=1
           js=is
           go to 8
          else if(isblk(1,is).eq.io2.and.isblk(2,is).eq.io1.and.
     $          isblk(3,is).eq.iv1.and.isblk(4,is).eq.iv2)then
           jcase=2
           js=is
           go to 8
          end if
         end do
         write(6,*)('jmat not found to match '),io1,io2,iv1,iv2
         do is=1,nsdlk
          write(6,7)is,(isblk(j,is),j=1,4)
    7     format(i5,5x,4i1)
         end do
         call dws_sync
         call dws_finalize
         stop
    8    continue
         istoh=istoh+1
         if(istoh.gt.idbk)then
          write(6,*)('istoh is tooooo large!! '),istoh
          write(6,*)('in paracis ')
          call dws_sync
          call dws_finalize
          stop
         end if
         nv1=nvirt(iv1)
         nv2=nvirt(iv2)
         call ilimts(nv1,nv2,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)
         nhere=ih+1-il
         nrow=noca(io1)*noca(io2)
         ihcis(istoh)=ibcoff
         ibcoff=ihcis(istoh)+nrow*nhere
         icbk(1,istoh)=io1
         icbk(2,istoh)=io2
         icbk(3,istoh)=iv1
         icbk(4,istoh)=iv2
         icbk(5,istoh)=1
  393    format('new istoh ',4i1,i2)
         call enough('paracis.  3',bc,ibc)
c
c     -J part and diagonal fock matrix part
c
         i10=i1s
         i1n=nv1
         iih=ihcis(istoh)
         iij=jmats(js)
         if(io1.eq.io2)then
          nrowj=(noc(io1)*(noc(io1)+1))/2
         else
          nrowj=noc(io1)*noc(io2)
         end if
         do i2=i2s,i2e
          if(i2.eq.i2e)i1n=i1e
          do i1=i10,i1n
           iih0=iih
           if(jcase.eq.1)then
            if(io1.eq.io2)then
             do i4=0,noca(io1)-1
              do i3=0,noca(io1)-1
               in=min(i3,i4)+icore(io1)
               ix=max(i3,i4)+icore(io1)
               iadj=iij+((ix*(ix+1))/2)+in
               bc(iih)=-bc(iadj)
               iih=iih+1
              end do
             end do
            else
             do i4=0,noca(io2)-1
              iadj=iij+icore(io1)+noc(io1)*(i4+icore(io2))
              do i3=0,noca(io1)-1
               bc(iih)=-bc(iadj+i3)
               iih=iih+1
              end do
             end do
            end if
           else
            do i4=0,noca(io2)-1
             do i3=0,noca(io1)-1
              iadj=iij+i4+icore(io2)+noc(io2)*(i3+icore(io1))
              bc(iih)=-bc(iadj)
              iih=iih+1
             end do
            end do
           end if
           if(i1.eq.i2.and.iv1.eq.iv2)then
            iadv=ieigv(iv1)+noc(iv1)+i1-1
            eigv=bc(iadv)                                               5d3s16
            do i34=0,noca(io1)-1
             iadh=iih0+i34*(noca(io1)+1)
             iado=ieigv(io1)+i34+icore(io1)
             bc(iadh)=bc(iadh)+eigv-bc(iado)
            end do
           end if
           iij=iij+nrowj
          end do
          i10=1
         end do
c
c     K_{nm}^{ab} part
c
         do is=1,nsdlkk
          if(isblkk(1,is).eq.io1.and.isblkk(2,is).eq.io2.and.
     $         isblkk(3,is).eq.iv1.and.isblkk(4,is).eq.iv2)then
           kcase=1
           ks=is
           go to 6
          else if(isblkk(1,is).eq.io2.and.isblkk(2,is).eq.io1.and.
     $          isblkk(3,is).eq.iv2.and.isblkk(4,is).eq.iv1)then
           kcase=2
           ks=is
           go to 6
          end if
         end do
         write(6,*)('k_{nm}^{ab} not found to match '),io1,io2,iv1,iv2
         do is=1,nsdlkk
          write(6,7)is,(isblkk(j,is),j=1,4)
         end do
         call dws_sync
         call dws_finalize
         stop
    6    continue
         nrowk=noc(io1)*noc(io2)
         istoh0=0
c
c     K_{mm}^{ab} part
c
         do is=1,nsdlkk
          if(isblkk(1,is).eq.io2.and.isblkk(2,is).eq.io1.and.
     $         isblkk(3,is).eq.iv1.and.isblkk(4,is).eq.iv2)then
           lcase=1
           ls=is
           go to 16
          else if(isblkk(1,is).eq.io1.and.isblkk(2,is).eq.io2.and.
     $          isblkk(3,is).eq.iv2.and.isblkk(4,is).eq.iv1)then
           lcase=2
           ls=is
           go to 16
          end if
         end do
         write(6,*)('k_{mn}^{ab} not found to match '),io1,io2,iv1,iv2
         do is=1,nsdlkk
          write(6,7)is,(isblkk(j,is),j=1,4)
         end do
         call dws_sync
         call dws_finalize
         stop
   16    continue
         if(lcase.eq.1)then
          if(istoh0.eq.0)then
           jstoh=istoh
          else
           jstoh=istoh0
          end if
          i10=i1s
          i1n=nv1
          iik=kmats(ls)
          iih=ihcis(jstoh)
          do i2=i2s,i2e
           if(i2.eq.i2e)i1n=i1e
           do i1=i10,i1n
            do i4=0,noca(io2)-1
             i4p=i4+icore(io2)
             do i3=0,noca(io1)-1
              i3p=i3+icore(io1)
              iadk=iik+i4p+noc(io2)*i3p
              bc(iih)=bc(iih)+bc(iadk)*2d0
              iih=iih+1
             end do
            end do
            iik=iik+nrowk
           end do
           i10=1
          end do
         else
          if(istoh0.eq.0)then
           istoh0=istoh
           istoh=istoh+1
           if(istoh.gt.idbk)then
            write(6,*)('istoh is tooooo large for lcase=2!! '),istoh
            write(6,*)('in paracis ')
            call dws_sync
            call dws_finalize
            stop
           end if
           call ilimts(nv2,nv1,mynprocg,mynowprog,ilo,iho,i1so,i1eo,i2so,
     $         i2eo)
           nhereo=iho+1-ilo
           ihcis(istoh)=ibcoff
           ibcoff=ihcis(istoh)+nrow*nhereo
           icbk(1,istoh)=io1
           icbk(2,istoh)=io2
           icbk(3,istoh)=iv2
           icbk(4,istoh)=iv1
           icbk(5,istoh)=3
           write(6,393)(icbk(j,istoh),j=1,5)
           call enough('paracis.  4',bc,ibc)
          end if
          i10=i1so
          i1n=nv2
          iih=ihcis(istoh)
          iik=kmats(ls)
          do i2=i2so,i2eo
           if(i2.eq.i2eo)i1n=i1eo
           do i1=i10,i1n
            do i4=0,noca(io2)-1
             i4p=i4+icore(io2)
             do i3=0,noca(io1)-1
              i3p=i3+icore(io1)
              iadk=iik+i3p+noc(io1)*i4p
              bc(iih)=bc(iadk)*2d0
              iih=iih+1
             end do
            end do
            iik=iik+nrowk
           end do
           i10=1
          end do
         end if
        end do
       end do
       if(mynowprog.eq.0)write(6,*)('we have all of cis matrix ')
       do i=1,istoh
        if(icbk(5,i).eq.1)then
        else if(icbk(5,i).eq.2)then
        else
        end if
   12   format('integral type: ',4i1)
        nv1=nvirt(icbk(3,i))
        nv2=nvirt(icbk(4,i))
        call ilimts(nv1,nv2,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)
        nrow=noca(icbk(1,i))*noca(icbk(2,i))
        nhere=ih+1-il
       end do
        if(mynowprog.eq.0)write(6,*)('setting up rhs for response ')    3d6s17
        ioffp=0                                                         5d9s16
        do isb=1,nsymb                                                  5d9s16
         iph(isb)=ioffp                                                 5d9s16
         ioffp=ioffp+nbasdws(isb)*nbasdws(isb)                          5d9s16
        end do                                                          5d9s16
        natom3=natom*3
        natom4=natom3+natom                                             6d20s16
        nthis=0                                                         5d3s16
        do irhs=1,natom3                                                5d3s16
         irhsp=irhs+natom4                                              6d20s16
         if(ipropmat(irhsp).ne.-isw)nthis=nthis+1                        5d3s16
        end do                                                          5d3s16
        if(mynowprog.eq.0)write(6,*)                                    3d6s17
     $       ('number of right hand sides of this symmetry = '),
     $       nthis
        nthis0=nthis                                                    5d9s16
        if(isw.ne.-1)then                                                5d9s16
         if(mynowprog.eq.0)write(6,*)('adding squared contributions ')  3d6s17
         nthis=nthis+natom                                              5d9s16
         if(mynowprog.eq.0)                                             3d6s17
     $        write(6,*)('total number of right hand sides for this '), 3d6s17
     $        ('symmetry = '),nthis                                     5d9s16
         i2nd=ibcoff                                                    12d12s16
         ibcoff=i2nd+natom*n2sz                                         12d12s16
         call enough('paracis.  5',bc,ibc)
         do i=0,natom*n2sz-1                                            12d12s16
          bc(i2nd+i)=0d0                                                12d12s16
         end do                                                         12d12s16
         do iat=1,natom
          do ixyz=1,3
           j2nd=i2nd+n2sz*(iat-1)                                        12d12s16
           iad=ixyz+3*(iat-1)
           ipad=iad+natom4                                              6d20s16
           ioffp=ipropmat(iad)                                          5d9s16
           do isb=1,nsymb                                                  5d9s16
            isk=multh(isb,ipropmat(ipad))
            iph(isb)=ioffp                                                 5d9s16
            ioffp=ioffp+nbasdwsc(isb)*nbasdwsc(isk)                          5d9s16
           end do                                                          5d9s16
           jsq=ioffp                                                    12d12s16
           do isb=1,nsymb
c
c     the -1/2 comes from definition of b0                               3d2s17
c
            do i=0,nbasdwsc(isb)*nbasdwsc(isb)-1                        12d12s16
             bc(j2nd+i)=bc(j2nd+i)-0.5d0*bc(jsq+i)                      3d2s17
            end do                                                      12d12s16
            jsq=jsq+nbasdwsc(isb)*nbasdwsc(isb)                         12d12s16
            j2nd=j2nd+nbasdwsc(isb)*nbasdwsc(isb)                       12d12s16
           end do
          end do
          j2nd=i2nd+n2sz*(iat-1)                                        12d12s16
          do isb=1,nsymb                                                12d12s16
           j2nd=j2nd+nbasdwsc(isb)*nbasdwsc(isb)                        2d1s17
          end do                                                        12d12s16
         end do                                                         12d12s16
        end if                                                          5d9s16
        if(nthis.gt.0)then
         irhs=ibcoff
         ibcoff=irhs+ntot*nthis
         igoal=irhs+ntot*9
         ipvt=ibcoff                                                    5d3s16
         ibcoff=ipvt+ntot                                               5d3s16
         call enough('paracis.  6',bc,ibc)
         do i=0,ntot*nthis-1                                            2d27s17
          bc(irhs+i)=0d0                                                2d27s17
         end do                                                         2d27s17
         ioff=0
         do ib=1,nbk                                                    5d3s16
          iv=ibk(ib)
          io=multh(iv,isw)
          ii=0
          do i=1,natom3
           ip=i+natom4                                                  6d20s16
           if(ipropmat(ip).ne.-isw)then
            ioffp=ipropmat(i)
            do isb=1,nsymb
             isbo=multh(isb,ipropmat(ip))                               5d4s22
             ioffp=ipropmat(i)+ipt2(iv,io)                              5d2s22
c
c     propmat stored under ioffp is nbasdwsc(isb),nbasdwsc(isbo)        3d1s17
c     we want der to be on occupied orbital, io - this is ket,          3d1s17
c     so need io=isbo and iv=isb.                                       3d1s17
c
             if(io.eq.isbo.and.iv.eq.isb)then                           3d1s17
              iad=irhs+ioff+ii
              do i2=1,nvirtc(iv)                                        3d1s17
               i2p=i2+noc(iv)-1+ioffp+(icore(io)-1)*nbasdwsc(iv)        3d1s17
               do i1=1,noca(io)                                         3d1s17
                bc(iad)=-bc(i2p+i1*nbasdwsc(iv))                        3d2s17
                iad=iad+1
               end do
              end do
             end if
             ioffp=ioffp+nbasdwsc(isb)*nbasdwsc(isbo)                   5d3s16
            end do
            ii=ii+ntot
           end if
          end do
          if(isw.ne.-1)then                                              5d9s16
           do iat=1,natom                                               5d9s16
            ip=iat+natom4                                                 12d12s16
            j2nd=i2nd+n2sz*(iat-1)                                      12d12s16
c
c     one electron contribution
c
            do isb=1,nsymb                                              5d9s16
             if(isb.eq.io.and.iv.eq.isb)then                            5d4s22
              iad=irhs+ioff+ii                                          5d9s16
c
c     2nd der matrix is not symmetric ... so we need ket to
c     be the hf orbital that is excited out of.
c
              do i2=1,nvirtc(isb)                                       5d9s16
               i2p=i2+noc(isb)-1+j2nd+(icore(isb)-1)*nbasdwsc(isb)      2d28s17
               do i1=1,noca(isb)                                        5d9s16
                bc(iad)=bc(iad)+bc(i2p+i1*nbasdwsc(isb))                2d28s17
                iad=iad+1                                               5d9s16
               end do                                                   5d9s16
              end do                                                    5d9s16
             end if                                                     5d9s16
             j2nd=j2nd+nbasdwsc(isb)*nbasdwsc(isb)                      12d12s16
            end do                                                      5d9s16
c
c     two electron contribution                                         12d29s16
c     [mp][nn] only contributes when perturbation is totally symmetric
c     [mn][np] only contributes when sym m* sym n=sym n * sym p = sym of
c     perturbation
c     in this part of code, isb points to n
c
            do ixyz=1,3                                                 12d29s16
             ift=ixyz+3*(iat-1)                                         12d29s16
             ift1=ift+natom4                                            12d29s16
c
c     [mp][nn]
c
c
c     -[mn][np]
c
             isb=multh(io,ipropmat(ift1))
              joff=ipropmat(ift)
              joffn=ipt2(isb,io)+joff                                   4d29s22
              joffmn=ipt2(io,isb)+joff                                  4d29s22
              iad=irhs+ioff+ii                                          12d29s16
              if(isw.eq.1)then                                          2d24s23
c     joffn iv,isb
c     joffmn io,isb
               do ip=0,nvirtc(iv)-1                                      12d29s16
                ipp=(ip+noc(iv))*nbasdwsc(isb)+joffn                     12d29s16
                do im=0,noca(io)-1                                       12d29s16
                 imm=im+icore(io)+joffmn                                 12d29s16
                 do in=0,noc(isb)-1                                      12d29s16
c                                                                       3d2s17
c     factor of 2 comes from 2 multiplying mixed der                    3d2s17
c     but then b0 is divided by -2, so no net factor.                   3d2s17
c                                                                       3d2s17
                  bc(iad)=bc(iad)+bc(ipp+in)*bc(imm+nbasdwsc(io)*in)     3d2s17
                 end do                                                  12d29s16
                 iad=iad+1                                               12d29s16
                end do                                                   12d29s16
               end do
              end if                                                    5d4s22
             end do                                                      12d29s16
            ii=ii+ntot                                                  5d9s16
           end do                                                       5d9s16
          end if
          ioff=ioff+noca(io)*nvirtc(iv)
         end do
         nsymx=0                                                        3d3s23
         do ia=1,natom                                                  3d3s23
          if(iapair(1,ia).gt.0)then                                     3d3s23
           nsymx=nsymx+1                                                3d3s23
           iad1=irhs+ntot*(natom*3+ia-1)-1                              3d3s23
           iad2=irhs+ntot*(natom*3+iapair(1,ia)-1)-1                    3d3s23
           do j=1,ntot                                                  3d3s23
            avg=0.25d0*(bc(iad1+j)+bc(iad2+j))                          3d3s23
            bc(iad1+j)=avg                                              3d3s23
            bc(iad2+j)=avg                                              3d3s23
           end do                                                       3d3s23
          end if                                                        3d3s23
         end do                                                         3d3s23
         itmp=ibcoff
         ibcoff=itmp+ntot*nthis
         call enough('paracis.  7',bc,ibc)
         do i=0,nthis*ntot-1
          bc(itmp+i)=bc(irhs+i)
         end do
c
c     solve linear equations via conjugate gradient
c
         ntotm=ntot-1                                                    3d3s17
         nthism=nthis-1                                                  3d3s17
         nun=ntot*nthis                                                 4d4s23
         nunm=nun-1                                                     4d4s23
         ihdiag=ibcoff                                                  4d4s23
         ix=ihdiag+ntot                                                 4d4s23
         ir=ix+nun                                                      4d4s23
         iz=ir+nun                                                      4d4s23
         iq=iz+nun                                                      4d4s23
         ip=iq+nun                                                      4d4s23
         ibeta=ip+nun                                                   4d4s23
         iscale=ibeta+nthis
         idotrz=iscale+nthis
         idotrznew=idotrz+nthis
         ialpha=idotrznew+nthis                                         4d4s23
         ibcoff=ialpha+nthis                                            4d4s23
         call enough('paracis.x',bc,ibc)
         do i=0,ntot-1                                                   3d3s17
          bc(ihdiag+i)=0d0                                               3d3s17
         end do                                                          3d3s17
         jx=ix                                                          4d4s23
         jrhs=irhs                                                      4d4s23
         jr=ir                                                          4d4s23
         ss=0d0                                                         6d26s23
         do icol=0,nthis-1                                               3d3s17
          scale=0d0                                                     4d4s23
          do irow=0,ntot-1                                               3d3s17
           bc(jx+irow)=0d0                                              4d4s23
           bc(jr+irow)=bc(jrhs+irow)                                    4d4s23
           scale=scale+bc(jrhs+irow)**2                                 4d4s23
          end do                                                         3d3s17
          ss=ss+scale                                                   6d26s23
          if(abs(scale).ne.0d0)then                                     4d12s23
           bc(iscale+icol)=1d0/scale                                     4d4s23
          else                                                          4d4s23
           bc(iscale+icol)=1d0                                          4d4s23
          end if                                                        4d4s23
          jx=jx+ntot                                                    4d4s23
          jrhs=jrhs+ntot                                                 3d3s17
          jr=jr+ntot                                                    4d4s23
         end do                                                          3d3s17
         ioff1=0
         do ib1=1,nbk
          iv1=ibk(ib1)
          io1=multh(iv1,isw)
          ioff2=0
          do ib2=1,ib1
           iv2=ibk(ib2)
           io2=multh(iv2,isw)
           do i=1,istoh
            if(icbk(5,i).eq.1.and.icbk(1,i).eq.io1.and.
     $         icbk(2,i).eq.io2.and.icbk(3,i).eq.iv1.and.
     $         icbk(4,i).eq.iv2)then
             nv1=nvirt(iv1)
             nv2=nvirt(iv2)
             call ilimts(nv1,nv2,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,
     $            i2e)
             i10=i1s
             i1n=nv1
             iih=ihcis(i)
             nrow=noca(io1)*noca(io2)
             do i2=i2s,i2e
              if(i2.eq.i2e)i1n=i1e
              do i1=i10,i1n
               do i4=0,noca(io2)-1
                i4p=i4+noca(io2)*(i2-1)+ioff2
                do i3=0,noca(io1)-1
                 i3p=i3+noca(io1)*(i1-1)+ioff1
                 if(i3p.eq.i4p)then                                        3d3s17
                  iad=ihdiag+i3p                                           3d3s17
                  bc(iad)=bc(iad)+bc(iih)                                  3d3s17
                 end if                                                    3d3s17
                 iih=iih+1
                end do
               end do
              end do
              i10=1
             end do
            else if(icbk(5,i).ne.1.and.icbk(1,i).eq.io1.and.
     $          icbk(2,i).eq.io2.and.icbk(3,i).eq.iv2.and.
     $          icbk(4,i).eq.iv1)then
             nv1=nvirt(iv2)
             nv2=nvirt(iv1)
             call ilimts(nv1,nv2,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,
     $            i2e)                                                  4d4s23
             i10=i1s
             i1n=nv1
             iih=ihcis(i)
             nrow=noca(io1)*noca(io2)
             do i2=i2s,i2e
              if(i2.eq.i2e)i1n=i1e
              do i1=i10,i1n
               do i4=0,noca(io2)-1
                i4p=i4+noca(io2)*(i1-1)+ioff2
                do i3=0,noca(io1)-1
                 i3p=i3+noca(io1)*(i2-1)+ioff1
                 if(i3p.eq.i4p)then                                        3d3s17
                  iad=ihdiag+i3p                                           3d3s17
                  bc(iad)=bc(iad)+bc(iih)                                  3d3s17
                 end if                                                    3d3s17
                 iih=iih+1
                end do
               end do
              end do
              i10=1
             end do
            end if
           end do
           ioff2=ioff2+noca(io2)*nvirt(iv2)
          end do
          ioff1=ioff1+noca(io1)*nvirt(iv1)
         end do
         call dws_gsumf(bc(ihdiag),ntot)                                   3d3s17
         do i=0,ntotm                                                      3d3s17
          bc(ihdiag+i)=1d0/bc(ihdiag+i)                                    3d3s17
         end do                                                            3d3s17
         jz=iz                                                             3d3s17
         jr=ir                                                          4d4s23
         xtest=0d0                                                      4d4s23
         jx=ix                                                          4d4s23
         jp=ip                                                          4d4s23
         do icol=0,nthis-1                                                 3d3s17
          xdot=0d0                                                      4d4s23
          dotrz=0d0                                                     4d4s23
          do irow=0,ntotm                                                  3d3s17
           bc(jz+irow)=bc(jr+irow)*bc(ihdiag+irow)                      4d4s23
           bc(jx+irow)=0d0                                              4d4s23
           xdot=xdot+bc(jr+irow)**2                                     4d4s23
           dotrz=dotrz+bc(jz+irow)*bc(jr+irow)                          4d4s23
          end do                                                           3d3s17
          bc(idotrz+icol)=dotrz                                         4d4s23
          xtest=xtest+xdot*bc(iscale+icol)                              4d4s23
          jz=jz+ntot                                                       3d3s17
          jr=jr+ntot                                                    4d4s23
          jp=jp+ntot                                                    4d4s23
          jx=jx+ntot                                                    4d4s23
         end do                                                            3d3s17
         xtest=sqrt(xtest/dfloat(nthis))                                4d4s23
         ss=sqrt(ss/dfloat(nun))                                        6d26s23
         if(ss.gt.1d-10)then                                            6d26s23
          macit=0                                                        6d26s23
          tol=1d-14                                                         3d3s17
          itmax=10                                                       6d26s23
          ierror=ibcoff                                                  6d26s23
          isoln=ierror+itmax                                             6d26s23
          ibcoff=isoln+ntot*nthis                                       6d26s23
          call enough('paracis.err',bc,ibc)                              6d26s23
          do izz=isoln,ibcoff-1                                           6d26s23
           bc(izz)=0d0                                                    6d26s23
          end do                                                         6d26s23
  101     continue                                                       6d26s23
           macit=macit+1                                                 6d26s23
           if(macit.gt.100)then                                          6d26s23
            write(6,*)('too many restarts!!!')
            call dws_synca
            call dws_finalize
            stop 'paracisa'
           end if
           do i=0,nunm                                                   6d26s23
            bc(ix+i)=0d0                                                 6d26s23
            bc(ip+i)=bc(iz+i)                                            6d26s23
           end do                                                        6d26s23
           iter=0                                                        6d26s23
  103      continue                                                      6d26s23
            iter=iter+1                                                  6d26s23
            if(iter.gt.itmax)then                                        6d26s23
             if(mynowprog.eq.0)then                                      6d26s23
              write(6,52)(bc(ierror+izz),izz=0,iter-2)                   6d26s23
             end if                                                      6d26s23
             do i=0,nunm                                                  6d26s23
              bc(isoln+i)=bc(isoln+i)+bc(ix+i)                           6d26s23
             end do                                                       6d26s23
             call mtimesi(bc(iq),bc(isoln),nun,nbk,multh,isw,istoh,icbk,  6d26s23
     $          nvirt,ihcis,noca,nthism,ntot,ibk,bc,ibc)                6d26s23
             jr=ir                                                        6d26s23
             jq=iq                                                        6d26s23
             jz=iz                                                        6d26s23
             jrhs=irhs                                                  6d26s23
             xtest=0d0                                                    6d26s23
             do icol=0,nthism                                           6d26s23
              xdot=0d0                                                     4d4s23
              dotrz=0d0                                                     6d26s23
              do irow=0,ntot-1                                               3d3s17
               bc(jr+irow)=bc(jrhs+irow)-bc(jq+irow)                     6d26s23
               bc(jz+irow)=bc(jr+irow)*bc(ihdiag+irow)                      6d26s23
               xdot=xdot+bc(jr+irow)**2                                   6d26s23
               dotrz=dotrz+bc(jr+irow)*bc(jz+irow)                          6d26s23
              end do                                                         3d3s17
              xtest=xtest+xdot*bc(iscale+icol)                            6d26s23
              bc(idotrz+icol)=dotrz                                       6d26s23
              jr=jr+ntot                                                 6d26s23
              jrhs=jrhs+ntot                                            6d26s23
              jq=jq+ntot                                                 6d26s23
              jz=jz+ntot                                                 6d26s23
             end do                                                       6d26s23
             xtest=sqrt(xtest/dfloat(nthis))                            6d26s23
             go to 101                                                  6d26s23
            end if                                                        6d26s23
            if(iter.ne.1)then                                             6d26s23
             jp=ip                                                        6d26s23
             jz=iz                                                        6d26s23
             do icol=0,nthism                                           6d26s23
              if(bc(idotrz+icol).ne.0d0)then                            6d26s23
               beta=bc(idotrznew+icol)/bc(idotrz+icol)                     6d26s23
               do irow=0,ntot-1                                            6d26s23
                bc(jp+irow)=bc(jz+irow)+beta*bc(jp+irow)                   6d26s23
               end do                                                      6d26s23
              end if                                                    6d26s23
              bc(idotrz+icol)=bc(idotrznew+icol)                          6d26s23
              jp=jp+ntot                                                  6d26s23
              jz=jz+ntot                                                  6d26s23
             end do                                                       6d26s23
            end if                                                        6d26s23
            call mtimesi(bc(iq),bc(ip),nun,nbk,multh,isw,istoh,icbk,      6d26s23
     $          nvirt,ihcis,noca,nthism,ntot,ibk,bc,ibc)                6d26s23
            jp=ip                                                         6d26s23
            jq=iq                                                         6d26s23
            jx=ix                                                         6d26s23
            jr=ir                                                         6d26s23
            jz=iz                                                         6d26s23
            xtest=0d0                                                     6d26s23
            do icol=0,nthism                                            6d26s23
             alpha=0d0                                                    6d26s23
             do irow=0,ntot-1                                             6d26s23
              alpha=alpha+bc(jp+irow)*bc(jq+irow)                         6d26s23
             end do                                                       6d26s23
             if(abs(alpha).ne.0d0)then                                  6d26s23
              alpha=bc(idotrz+icol)/alpha                                  6d26s23
             end if                                                     6d26s23
             xdotnew=0d0                                                  6d26s23
             dotrznew=0d0                                                 6d26s23
             do irow=0,ntot-1                                             6d26s23
              bc(jx+irow)=bc(jx+irow)+alpha*bc(jp+irow)                   6d26s23
              bc(jr+irow)=bc(jr+irow)-alpha*bc(jq+irow)                   6d26s23
              bc(jz+irow)=bc(jr+irow)*bc(ihdiag+irow)                     6d26s23
              xdotnew=xdotnew+bc(jr+irow)**2                              6d26s23
              dotrznew=dotrznew+bc(jr+irow)*bc(jz+irow)                   6d26s23
             end do                                                       6d26s23
             bc(idotrznew+icol)=dotrznew                                  6d26s23
             xtest=xtest+dotrznew*bc(iscale+icol)                         6d26s23
             jp=jp+ntot                                                   6d26s23
             jq=jq+ntot                                                   6d26s23
             jx=jx+ntot                                                   6d26s23
             jr=jr+ntot                                                   6d26s23
             jz=jz+ntot                                                   6d26s23
            end do                                                        6d26s23
            xtest=sqrt(xtest/dfloat(nthis))                             6d26s23
            bc(ierror+iter-1)=xtest
            if(xtest.gt.tol)go to 103                                     6d26s23
           if(mynowprog.eq.0)then
            write(6,52)(bc(ierror+izz),izz=0,iter-1)
            write(6,*)('cg iterations converged after '),iter,
     $      (' iterations and '),macit-1,(' restarts')
           end if                                                        6d26s23
           do i=0,nunm                                                   6d26s23
            bc(irhs+i)=bc(isoln+i)+bc(ix+i)                             6d26s23
           end do                                                        6d26s23
          end if
       ibcoff=ix                                                        4d4s23
       ihitc2=0                                                         2d22s23
         ipt=ibcoff
         ibcoff=ipt+nthis0*2
         call enough('paracis.  9',bc,ibc)
         jpt=ipt
         do iat=1,natom
          do ixyz=1,3
           iarg=ixyz+3*(iat-1)+natom4
           if(ipropmat(iarg).ne.-isw)then
            ibc(jpt)=iat
            kpt=jpt+nthis0
            jpt=jpt+1
            ibc(kpt)=ixyz
           end if
          end do
         end do
         jc2=ipt
         ic2=ipt+nthis0
          do i=0,nthis0-1                                               5d9s16
           iadi=irhs+ntot*i                                             5d9s16
           do j=0,nthis0-1                                              5d9s16
            iadj=itmp+ntot*j                                            5d9s16
            sum=0d0                                                     5d9s16
            do k=0,ntot-1                                               5d9s16
             term=bc(iadi+k)*bc(iadj+k)                                 5d4s22
             sum=sum+term                                               5d4s22
            end do                                                      5d9s16
c
c     minus sign is because we solved linear equations with H-E0 rather
c     than E0-H.
c     factor of 2 comes from sr2 I multiplied 1st der matrix in molpro.
c     I think this is because in vtet the correcation is divided by 2,
c     but I need to verify this.
c
            sum2=-sum*2d0                                               6d8s16
            if(nsymb.gt.1)then                                          4d28s22
             iad=ic2mat+ibc(jc2+j)-1+natom*(ibc(ic2+j)-1                5d4s22
     $           +3*(ibc(jc2+i)-1+natom*(ibc(ic2+i)-1)))                5d4s22
c     adding here, because only one symmetry gives nonzero result
             bc(iad)=bc(iad)+sum2                                       8d4s22
            end if                                                      4d28s22
            if(mynowprog.eq.0.and.abs(sum2).gt.1d-10)then                   8d4s22
             if(ihitc2.eq.0)write(6,*)('>C(2) matrix: ')                2d22s23
             ihitc2=1                                                   2d22s23
             write(6,57)cart(ibc(ic2+j)),ibc(jc2+j),                    2d22s23
     $           cart(ibc(ic2+i)),ibc(jc2+i),sum2                       2d1s23
            end if                                                      2d22s23
   57       format('>',x,a1,i2,3x,a1,i2,es15.7,i12)                        5d4s22
           end do                                                       5d9s16
          end do                                                        5d9s16
c
c     if perturbation is not totally symmetric, potential ders will
c     be zero, so cal S is zero as well.
c
         if(isw.ne.-1)then                                               3d1s17
          icals=ibcoff                                                   2d27s17
          ibcoff=icals+ntot*natom                                        2d27s17
          call enough('paracis. 10',bc,ibc)
          do i=0,ntot*natom-1                                            2d27s17
           bc(icals+i)=0d0                                               2d27s17
          end do                                                         2d27s17
          ioff=0                                                         2d27s17
          do iat=1,3                                                     2d27s17
           do ixyz=1,3                                                   2d27s17
            iarg=ixyz+3*(iat-1)+natom4                                   2d27s17
            iarg2=iarg+natom3                                            2d27s17
            if(ipropmat(iarg).ne.-isw)then                                2d27s17
             itrans8=ipropmat(iarg2)                                      2d27s17
             iad1=icals+ntot*(iat-1)                                     2d27s17
             iad2=irhs+ntot*ioff                                         2d27s17
             orig=bc(iad1)
             do i=0,ntot-1                                               2d27s17
              bc(iad1+i)=bc(iad1+i)+trans8*bc(iad2+i)                    2d27s17
             end do                                                      2d27s17
             ioff=ioff+1                                                 2d27s17
            end if                                                       2d27s17
           end do                                                        2d27s17
          end do                                                         2d27s17
         nsymx=0                                                        3d3s23
         do ia=1,natom                                                  3d3s23
          if(iapair(1,ia).gt.0)then                                     3d3s23
           iad1=icals-1+ntot*(ia-1)                                     3d3s23
           iad2=icals-1+ntot*(iapair(1,ia)-1)                           3d3s23
           nsymx=nsymx+1                                                3d3s23
           do i=1,ntot                                                  3d3s23
            avg=0.25d0*(bc(iad1+i)+bc(iad2+i))                          3d3s23
            bc(iad1+i)=avg                                              3d3s23
            bc(iad2+i)=avg                                              3d3s23
           end do                                                       3d3s23
          end if                                                        3d3s23
         end do                                                         3d3s23
c
c     solve linear equations via conjugate gradient
c
         natomm=natom-1                                                  3d3s17
         nun=ntot*natom                                                 3d3s17
         nunm=nun-1                                                     4d4s23
         ix=ibcoff                                                      4d4s23
         ir=ix+nun                                                      4d4s23
         iz=ir+nun                                                      4d4s23
         iq=iz+nun                                                      4d4s23
         ip=iq+nun                                                      4d4s23
         ibeta=ip+nun                                                   4d4s23
         iscale=ibeta+nun                                               4d4s23
         idotrz=iscale+natom                                            4d4s23
         idotrznew=idotrz+natom                                         4d4s23
         ialpha=idotrznew+natom                                         4d4s23
         ibcoff=ialpha+natom                                            4d4s23
         call enough('paracis. 11',bc,ibc)
         jx=ix                                                          4d4s23
         jr=ir                                                          4d4s23
         jcals=icals
         jz=iz                                                          6d26s23
         ss=0d0
         do icol=0,natomm                                                3d3s17
          scale=0d0                                                     4d4s23
          dotrz=0d0                                                     6d26s23
          do irow=0,ntot-1                                               3d3s17
           bc(jx+irow)=0d0                                              4d4s23
           bc(jr+irow)=bc(jcals+irow)                                   4d4s23
           bc(jz+irow)=bc(jr+irow)*bc(ihdiag+irow)                      6d26s23
           dotrz=dotrz+bc(jr+irow)*bc(jz+irow)                          6d26s23
           scale=scale+bc(jcals+irow)**2                                         3d3s17
          end do                                                         3d3s17
          ss=ss+scale                                                   6d26s23
          if(scale.ne.0d0)then                                          4d4s23
           bc(iscale+icol)=1d0/scale                                    4d4s23
          else                                                          4d4s23
           bc(iscale+icol)=1d0                                          4d4s23
          end if                                                        4d4s23
          bc(idotrz+icol)=dotrz                                         6d26s23
          jx=jx+ntot                                                    4d4s23
          jr=jr+ntot                                                    4d4s23
          jz=jz+ntot                                                    6d26s23
          jcals=jcals+ntot                                                 3d3s17
         end do                                                          3d3s17
         ss=sqrt(ss/dfloat(nun))                                        6d26s23
         if(ss.gt.1d-10)then                                            6d26s23
          macit=0                                                        6d26s23
          tol=1d-14                                                         3d3s17
          itmax=10                                                       6d26s23
          ierror=ibcoff                                                  6d26s23
          isoln=ierror+itmax                                             6d26s23
          ibcoff=isoln+ntot*natom                                        6d26s23
          call enough('paracis.err',bc,ibc)                              6d26s23
          do izz=isoln,ibcoff-1                                           6d26s23
           bc(izz)=0d0                                                    6d26s23
          end do                                                         6d26s23
  201     continue                                                       6d26s23
           macit=macit+1                                                 6d26s23
           if(macit.gt.100)then                                          6d26s23
            write(6,*)('too many restarts!!!')
            call dws_synca
            call dws_finalize
            stop 'paracisb'
           end if
           do i=0,nunm                                                   6d26s23
            bc(ix+i)=0d0                                                 6d26s23
            bc(ip+i)=bc(iz+i)                                            6d26s23
           end do                                                        6d26s23
           iter=0                                                        6d26s23
  203      continue                                                      6d26s23
            iter=iter+1                                                  6d26s23
            if(iter.gt.itmax)then                                        6d26s23
             if(mynowprog.eq.0)then                                      6d26s23
              write(6,52)(bc(ierror+izz),izz=0,iter-2)                   6d26s23
   52         format(10es8.1)                                                 5d19s22
             end if                                                      6d26s23
             do i=0,nunm                                                  6d26s23
              bc(isoln+i)=bc(isoln+i)+bc(ix+i)                           6d26s23
             end do                                                       6d26s23
             call mtimesi(bc(iq),bc(isoln),nun,nbk,multh,isw,istoh,icbk,  6d26s23
     $          nvirt,ihcis,noca,natomm,ntot,ibk,bc,ibc)                6d26s23
             jr=ir                                                        6d26s23
             jq=iq                                                        6d26s23
             jz=iz                                                        6d26s23
             jcals=icals                                                  6d26s23
             xtest=0d0                                                    6d26s23
             do icol=0,natomm                                                3d3s17
              xdot=0d0                                                     4d4s23
              dotrz=0d0                                                     6d26s23
              do irow=0,ntot-1                                               3d3s17
               bc(jr+irow)=bc(jcals+irow)-bc(jq+irow)                     6d26s23
               bc(jz+irow)=bc(jr+irow)*bc(ihdiag+irow)                      6d26s23
               xdot=xdot+bc(jr+irow)**2                                   6d26s23
               dotrz=dotrz+bc(jr+irow)*bc(jz+irow)                          6d26s23
              end do                                                         3d3s17
              xtest=xtest+xdot*bc(iscale+icol)                            6d26s23
              bc(idotrz+icol)=dotrz                                       6d26s23
              jr=jr+ntot                                                 6d26s23
              jcals=jcals+ntot                                           6d26s23
              jq=jq+ntot                                                 6d26s23
              jz=jz+ntot                                                 6d26s23
             end do                                                       6d26s23
             xtest=sqrt(xtest/dfloat(natom))                              6d26s23
             go to 201                                                    6d26s23
            end if                                                        6d26s23
            if(iter.ne.1)then                                             6d26s23
             jp=ip                                                        6d26s23
             jz=iz                                                        6d26s23
             do icol=0,natomm                                             6d26s23
              beta=bc(idotrznew+icol)/bc(idotrz+icol)                     6d26s23
              bc(idotrz+icol)=bc(idotrznew+icol)                          6d26s23
              do irow=0,ntot-1                                            6d26s23
               bc(jp+irow)=bc(jz+irow)+beta*bc(jp+irow)                   6d26s23
              end do                                                      6d26s23
              jp=jp+ntot                                                  6d26s23
              jz=jz+ntot                                                  6d26s23
             end do                                                       6d26s23
            end if                                                        6d26s23
            call mtimesi(bc(iq),bc(ip),nun,nbk,multh,isw,istoh,icbk,      6d26s23
     $          nvirt,ihcis,noca,natomm,ntot,ibk,bc,ibc)                6d26s23
            jp=ip                                                         6d26s23
            jq=iq                                                         6d26s23
            jx=ix                                                         6d26s23
            jr=ir                                                         6d26s23
            jz=iz                                                         6d26s23
            xtest=0d0                                                     6d26s23
            do icol=0,natomm                                              6d26s23
             alpha=0d0                                                    6d26s23
             do irow=0,ntot-1                                             6d26s23
              alpha=alpha+bc(jp+irow)*bc(jq+irow)                         6d26s23
             end do                                                       6d26s23
             alpha=bc(idotrz+icol)/alpha                                  6d26s23
             xdotnew=0d0                                                  6d26s23
             dotrznew=0d0                                                 6d26s23
             do irow=0,ntot-1                                             6d26s23
              bc(jx+irow)=bc(jx+irow)+alpha*bc(jp+irow)                   6d26s23
              bc(jr+irow)=bc(jr+irow)-alpha*bc(jq+irow)                   6d26s23
              bc(jz+irow)=bc(jr+irow)*bc(ihdiag+irow)                     6d26s23
              xdotnew=xdotnew+bc(jr+irow)**2                              6d26s23
              dotrznew=dotrznew+bc(jr+irow)*bc(jz+irow)                   6d26s23
             end do                                                       6d26s23
             bc(idotrznew+icol)=dotrznew                                  6d26s23
             xtest=xtest+dotrznew*bc(iscale+icol)                         6d26s23
             jp=jp+ntot                                                   6d26s23
             jq=jq+ntot                                                   6d26s23
             jx=jx+ntot                                                   6d26s23
             jr=jr+ntot                                                   6d26s23
             jz=jz+ntot                                                   6d26s23
            end do                                                        6d26s23
            xtest=sqrt(xtest/dfloat(natom))                               6d26s23
            bc(ierror+iter-1)=xtest
            if(xtest.gt.tol)go to 203                                     6d26s23
           if(mynowprog.eq.0)then
            write(6,52)(bc(ierror+izz),izz=0,iter-1)
            write(6,*)('cg iterations converged after '),iter,
     $      (' iterations and '),macit-1,(' restarts')
           end if                                                        6d26s23
           do i=0,nunm                                                   6d26s23
            bc(isoln+i)=bc(isoln+i)+bc(ix+i)                             6d26s23
           end do                                                        6d26s23
           do i=0,nunm                                                      3d3s17
            bc(icals+i)=bc(isoln+i)                                         6d26s23
           end do                                                            3d3s17
       end if                                                           6d26s23
       ibcoff=ix                                                        4d4s23
       ihitc1=0                                                         2d22s23
          do jat=1,3
           iargs=icals+ntot*(jat-1)
           ioff=0
           do iat=1,3
            do ixyz=1,3
             iarg=ixyz+3*(iat-1)+natom4
             if(ipropmat(iarg).ne.-isw)then
              sum=0d0
              sumb=0d0
              iargbn=itmp+ntot*ioff                                       2d27s17
              if(isw.ne.-1)then
               iargb0=irhs+ntot*(nthis0+jat-1)                           2d27s17
               do i=0,ntot-1                                               2d27s17
c
c     minus sign on argb0 is because we solved linear equations with
c     1/(H-E0) rather than 1/(E0-H). This doesn't effect cal S term
c     because denominator is squared.
c
                term=bc(iargbn+i)*(-bc(iargb0+i)-0.5d0*bc(iargs+i))
                sum=sum-drsign*bc(iargbn+i)                             8d21s24
     $               *(bc(iargb0+i)+0.5d0*bc(iargs+i))                  8d21s24
                sumb=sumb+bc(iargbn+i)*bc(iargb0+i)
  202           format('>',3i3,es15.7,i3,4es15.7,i8)                       3d1s17
               end do                                                      2d27s17
              else
               do i=0,ntot-1                                               2d27s17
                sum=sum-drsign*bc(iargbn+i)*0.5d0*bc(iargs+i)           8d21s24
               end do                                                      2d27s17
              end if
              ioff=ioff+1
c
c     old code had minus sign since dVdx had wrong sign                 8d4s22
c     for factor of 2, see C(2) comment.
c
              if(nsymb.gt.1)then                                        4d28s22
               iad=ic1mat+jat-1+natom*(iat-1+natom*(ixyz-1))            4d28s22
c     adding here, because only one symmetry gives nonzero result
               bc(iad)=bc(iad)+sum*2d0                                  8d4s22
              end if                                                    4d28s22
              if(mynowprog.eq.0.and.abs(sum).gt.1d-10)then              2d22s23
               if(ihitc1.eq.0)write(6,*)('>C(1)' )                      2d22s23
               ihitc1=1                                                 2d22s23
               write(6,58)jat,cart(ixyz),iat,sum*2d0                    2d22s23
              end if                                                    2d22s23
   58         format('>',3x,i3,3x,a1,i2,es15.7,i12)                        5d4s22
             end if
            end do
           end do
          end do
          ihitc0=0                                                      2d22s23
          do iat=1,natom                                                3d3s23
           iargs=icals+ntot*(iat-1)
           iargeb0=irhs+ntot*(nthis0+iat-1)
           jargb0=itmp+ntot*(nthis0+iat-1)                              3d1s17
           do jat=1,natom                                               3d3s23
            iargb0=itmp+ntot*(nthis0+jat-1)
            jargs=icals+ntot*(jat-1)                                    3d1s17
            jargeb0=irhs+ntot*(nthis0+jat-1)                            3d1s17
            sum=0d0
            sumb=0d0
c
c     minus sign on argeb0 is because we solved linear equations with
c     1/(H-E0) rather than 1/(E0-H). This doesn't effect cal S term
c     because denominator is squared.
c     also symmetrize for ease in comparing to molpro
c
            do i=0,ntot-1
             sum=sum+bc(iargb0+i)*(-bc(iargeb0+i)-bc(iargs+i))
     $            +bc(jargb0+i)*(-bc(jargeb0+i)-bc(jargs+i))             3d1s17
            end do
c
c     why no - and factor of 2 here?
c     dont know about sign, but it makes sense there is no factor of 2
c     here, since this is a potential rather than ke term.
c
            if(mynowprog.eq.0.and.abs(sum).gt.1d-10)then                2d22s23
             if(ihitc0.eq.0)write(6,*)('>C(0)')                         2d22s23
             ihitc0=1                                                   2d22s23
             write(6,59)jat,iat,sum                                     2d22s23
             if(nsymb.gt.1)then                                         3d3s23
              iad=ic0mat+jat-1+natom*(iat-1)                             3d3s23
              bc(iad)=sum                                                3d3s23
             end if                                                     3d3s23
            end if                                                      2d22s23
   59       format('>',3x,i3,5x,i3,es15.7)                                  5d9s16
           end do
          end do
         end if
         ibcoff=irhs
        end if
       ibcoff=ihcis(1)
      end do
      ibcoff=ibcoffo
      if(nsymb.gt.1.and.mynowprog.eq.0)then                             4d28s22
       write(6,*)('b4 desymmetrization ')
       write(6,*)('symmetrized C(2) ')
       do ia=1,natom                                                    4d28s22
        iam=ia-1                                                        4d28s22
        do ixyz=1,3                                                      4d28s22
         ixyzm=ixyz-1                                                   4d28s22
         do ja=1,natom                                                  4d28s22
          jam=ja-1                                                      4d28s22
          do jxyz=1,3                                                   4d28s22
           jxyzm=jxyz-1                                                 4d28s22
           iad=ic2mat+jam+natom*(jxyzm+3*(iam+natom*ixyzm))             8d4s22
          end do                                                        4d28s22
         end do                                                         4d28s22
        end do                                                          4d28s22
       end do                                                           4d28s22
       write(6,*)('symmetrized C(1) ')
       do ia=1,natom                                                    4d28s22
        iam=ia-1                                                        4d28s22
        do ja=1,natom                                                   4d28s22
         jam=ja-1                                                       4d28s22
         do ixyz=1,3                                                      4d28s22
          ixyzm=ixyz-1                                                   4d28s22
          iad=ic1mat+jam+natom*(iam+natom*ixyzm)                        4d28s22
         end do                                                         4d28s22
        end do                                                          4d28s22
       end do                                                           4d28s22
       itmp=ibcoff                                                      5d4s22
       ibcoff=itmp+natom3**2                                            5d4s22
       call enough('paracis. 12',bc,ibc)
       do ipass=1,2
        do ixyz=0,2                                                     5d4s22
         do ia=1,natom                                                   5d4s22
          iam=ia-1                                                      5d4s22
          if(iapair(1,ia).eq.0)then                                     5d4s22
           iad1=ic2mat+natom3*(iam+natom*ixyz)                          5d4s22
           iad3=itmp+natom3*(iam+natom*ixyz)                            5d4s22
           do i=0,natom3-1                                              5d4s22
            bc(iad3+i)=bc(iad1+i)                                       5d4s22
           end do                                                       5d4s22
          else if(iapair(1,ia).gt.0)then                                5d4s22
           iad1=ic2mat+natom3*(iam+natom*ixyz)                          5d4s22
           iad2=ic2mat+natom3*(iapair(1,ia)-1+natom*ixyz)               5d4s22
           iad3=itmp+natom3*(iam+natom*ixyz)                            5d4s22
           iad4=itmp+natom3*(iapair(1,ia)-1+natom*ixyz)                 5d4s22
           do i=0,natom3-1                                              5d4s22
            o1=bc(iad1+i)
            o2=bc(iad2+i)
            bc(iad3+i)=0.5d0*(bc(iad1+i)+bc(iad2+i))                    5d4s22
            bc(iad4+i)=0.5d0*(bc(iad1+i)-bc(iad2+i))                    5d4s22
           end do                                                       5d4s22
          end if                                                        5d4s22
         end do                                                         5d4s22
        end do                                                          5d4s22
        do i=0,natom3-1                                                 5d4s22
         do j=0,natom3-1                                                5d4s22
          ji=itmp+j+natom3*i                                            5d4s22
          ij=ic2mat+i+natom3*j                                          5d4s22
          bc(ij)=bc(ji)                                                 5d4s22
         end do                                                         5d4s22
        end do                                                          5d4s22
       end do
       do ixyz=0,2                                                      5d4s22
        do ia=1,natom                                                   5d4s22
         iam=ia-1                                                       5d4s22
         if(iapair(1,ia).eq.0)then                                      5d4s22
          iad1=ic1mat+natom*(iam+natom*ixyz)                            5d4s22
          iad3=itmp+natom*(iam+natom*ixyz)                              5d4s22
          do i=0,natom-1                                                5d4s22
           bc(iad3+i)=bc(iad1+i)                                        5d4s22
          end do                                                        5d4s22
         else if(iapair(1,ia).gt.0)then                                 5d4s22
          iad1=ic1mat+natom*(iam+natom*ixyz)                            5d4s22
          iad2=ic1mat+natom*(iapair(1,ia)-1+natom*ixyz)                 5d4s22
          iad3=itmp+natom*(iam+natom*ixyz)                              5d4s22
          iad4=itmp+natom*(iapair(1,ia)-1+natom*ixyz)                   5d4s22
          do i=0,natom-1                                                5d4s22
           bc(iad3+i)=0.5d0*(bc(iad1+i)+bc(iad2+i))                     5d4s22
           bc(iad4+i)=0.5d0*(bc(iad1+i)-bc(iad2+i))                     5d4s22
          end do                                                        5d4s22
         end if                                                         5d4s22
        end do                                                          5d4s22
       end do                                                           5d4s22
       do i=0,natom3*natom-1                                            5d4s22
        bc(ic1mat+i)=bc(itmp+i)                                         5d4s22
       end do                                                           5d4s22
       ibcoff=itmp                                                      5d4s22
       write(6,*)('>desymmetrized C(2) ')
       do ia=1,natom                                                    4d28s22
        iam=ia-1                                                        4d28s22
        do ixyz=1,3                                                      4d28s22
         ixyzm=ixyz-1                                                   4d28s22
         do ja=1,natom                                                  4d28s22
          jam=ja-1                                                      4d28s22
          do jxyz=1,3                                                   4d28s22
           jxyzm=jxyz-1                                                 4d28s22
           iad=ic2mat+jam+natom*(jxyzm+3*(iam+natom*ixyzm))             5d4s22
           if(abs(bc(iad)).gt.1d-10)then                                3d3s23
            write(6,57)cart(jxyz),ja,cart(ixyz),ia,bc(iad)               4d28s22
           end if                                                       3d3s23
          end do                                                        4d28s22
         end do                                                         4d28s22
        end do                                                          4d28s22
       end do                                                           4d28s22
       write(6,*)('>desymmetrized C(1) ')
       do ja=1,natom                                                    4d28s22
        jam=ja-1                                                        4d28s22
        do ia=1,natom                                                   4d28s22
         iam=ia-1                                                       4d28s22
         do ixyz=1,3                                                    4d28s22
          ixyzm=ixyz-1                                                  4d28s22
          iad=ic1mat+jam+natom*(iam+natom*ixyzm)                        4d28s22
          if(abs(bc(iad)).gt.1d-10)then                                 3d3s23
           write(6,58)ja,cart(ixyz),ia,bc(iad)                           4d28s22
          end if                                                        3d3s23
         end do                                                         4d28s22
        end do                                                          4d28s22
       end do                                                           4d28s22
      end if                                                            4d28s22
      return
      end
