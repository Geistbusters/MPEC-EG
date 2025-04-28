c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine derh01(h0,h0d,idorel,ipsym,multh,ider,iorb,idwsdeb,    5d4s22
     $     itrans,nbasisp,isou,bc,ibc)                                  11d14s22
      implicit real*8 (a-h,o-z)
      include "common.hf"
      include "common.store"
      include "common.spher"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension ieigs(8),iorb(8),h0(*),h0d(*),                          5d4s22
     $    multh(8,8),ipt2(8),iptk(8),itrans(*),                         5d4s22
     $     nbasisp(*),isou(*)                                           4d19s22
      if(idwsdeb.gt.10)write(6,*)('Hi, my name is derh01!')
      ncomp=1
      srh=sqrt(0.5d0)
      if(idorel.ne.0)ncomp=2
      ider=ibcoff
      if(ipsym.eq.1)then
       ioff=1
       ioffmo=1
       jder=ider
       nwds=0                                                           12d5s16
       do isb=1,nsymb                                                   12d5s16
        nwds=nwds+nbasdws(isb)**2                                       5d4s22
       end do                                                           12d5s16
       ibcoff=ider+nwds                                                 5d4s22
       ibcoffo=ibcoff                                                   12d15s16
       do isb=1,nsymb
        nh=nbasdws(isb)                                                 4d7s22
        nhp=nbasisp(isb)*ncomp                                           4d7s22
        if(nh.gt.0)then                                                 5d20s22
         if(idwsdeb.gt.10)then
          write(6,*)('h0 in mo basis '),ioff
          call prntm2(h0(ioffmo),nh,nh,nh)
          if(nsymb.eq.1)then
           call printa(h0(ioffmo),nh,0,1,0,nh,0,1,0,bc(ibcoff))
          end if
         end if
         if(idwsdeb.gt.10)then
          write(6,*)('h0d before squaring: ')
          call mpprnt2(h0d(ioff),nhp)
         end if
         call square(h0d(ioff),nhp)
         iatmp1=ibcoff                                                   4d6s22
         iatmp2=iatmp1+nhp*nhp                                            4d6s22
         iatmp3=iatmp2+nhp*nhp                                            4d6s22
         iatmp4=iatmp3+nhp*nhp                                            4d6s22
         iatmp5=iatmp4+nhp*nhp                                            4d6s22
         iatmp6=iatmp5+nhp*nhp                                            4d6s22
         ibcoff=iatmp6+nhp*nhp                                            4d6s22
         call enough('derh01.  1',bc,ibc)
         do ii=0,nh*nh-1                                                 5d4s22
          bc(iatmp1+ii)=h0(ioffmo+ii)                                    5d4s22
         end do                                                          5d4s22
         do ii=0,nhp*nhp-1                                                4d6s22
          bc(iatmp3+ii)=h0d(ioff+ii)                                     4d6s22
         end do                                                          4d6s22
         nrow=nhp                                                        4d7s22
         if(idwsdeb.gt.10)then
          write(6,*)('orbitals ')
          call prntm2(bc(iorb(isb)),nhp,nh,nhp)
         end if
         do ipass=1,2                                                    4d6s22
          call dgemm('n','n',nrow,nh,nhp,1d0,bc(iatmp3),nrow,            4d7s22
     $        bc(iorb(isb)),nhp,0d0,bc(iatmp4),nrow,
     $         'derh01.1',                                              7d11s22
     d'derh01.  1')
          do i=0,nh-1                                                    4d6s22
           do j=0,nrow-1                                                  4d6s22
            ji=iatmp4+j+nrow*i                                            4d6s22
            ij=iatmp3+i+nh*j                                             4d6s22
            bc(ij)=bc(ji)                                                4d6s22
           end do                                                        4d6s22
          end do                                                         4d6s22
          nrow=nh                                                        4d7s22
         end do                                                          4d6s22
         if(idwsdeb.gt.10)then
          write(6,*)('atmp1 (h0 in mo)')
          call prntm2(bc(iatmp1),nh,nh,nh)
          write(6,*)('atmp3 (dh0 in mo)')
          call prntm2(bc(iatmp3),nh,nh,nh)
          write(6,*)('trans ')
          call prntm2(bc(itrans(isb)),nh,nh,nh)
          if(nsymb.eq.1)then
           itimes2=ibcoff
           ibcoff=itimes2+nh*nh
           do ix2x=0,nh*nh-1
            bc(itimes2+ix2x)=bc(itrans(isb)+ix2x)*2d0
           end do
           call printa(bc(itimes2),nh,0,1,0,nh,0,1,0,bc(ibcoff))
           ibcoff=itimes2
          end if
         end if
         call dgemm('n','n',nh,nh,nh,1d0,bc(iatmp1),nh,bc(itrans(isb)),  4d7s22
     $        nh,0d0,bc(iatmp2),nh,
     $        'derh01.2',
     d'derh01.  2')
         if(idwsdeb.gt.10)then
          write(6,*)('atmp2 =h0*trans')
          call prntm2(bc(iatmp2),nh,nh,nh)
         end if
         do i=0,nh-1                                                     4d6s22
          do j=0,nh-1                                                    4d6s22
           ji=j+nh*i                                                     4d6s22
           ij=i+nh*j                                                     4d6s22
           bc(iatmp1+ji)=bc(iatmp2+ji)+bc(iatmp2+ij)+bc(iatmp3+ji)       4d6s22
          end do                                                         4d6s22
         end do                                                          4d6s22
         if(idwsdeb.gt.10)then
          write(6,*)('atmp1 = der1 ')
          call prntm2(bc(iatmp1),nh,nh,nh)
         end if
         do i=0,nh*nh-1
          bc(jder+i)=bc(iatmp1+i)                                         4d6s22
         end do
         ibcoff=iatmp1                                                   4d6s22
         jder=jder+nh*nh                                                 12d5s16
         ioff=ioff+nhp*nhp                                               4d15s22
         ioffmo=ioffmo+nh*nh                                             5d4s22
        end if                                                          5d20s22
        ibcoff=ibcoffo                                                  3d22s16
       end do
      else
c
c     perturbation is not a1
c
       if(idwsdeb.ne.0)write(6,*)('perturbation is not a1')             7d14s22
       jder=ider
       nx=0                                                             7d5s16
       do isb=1,nsymb                                                    6d21s16
        isk=multh(isb,ipsym)                                            7d19s16
        iptk(isb)=ider+nx                                               6d8s22
        nx=nx+nbasdws(isb)*nbasdws(isk)                                 6d7s22
       end do                                                            6d21s16
       ibcoff=ider+nx                                                   6d8s22
       idera=ider
       call enough('derh01.  2',bc,ibc)
       do iz=ider,ibcoff-1                                              6d8s22
        bc(iz)=0d0                                                      6d8s22
       end do                                                           4d7s22
       ioff=1                                                           4d7s22
       ioffmo=1                                                         5d4s22
       jderab=idera                                                     4d7s22
       ioffd=1                                                          4d7s22
       do isb=1,nsymb                                                   4d7s22
        if(idwsdeb.ne.0)write(6,*)('for isb = '),isb                    7d14s22
        nh=nbasdws(isb)                                                 4d7s22
        nhp=nbasisp(isb)*ncomp                                          4d7s22
        iatmp1=ibcoff                                                   4d7s22
        iatmp2=iatmp1+nhp*nhp                                           4d7s22
        iatmp3=iatmp2+nhp*nhp                                           4d7s22
        iatmp4=iatmp3+nhp*nhp                                           4d7s22
        ibcoff=iatmp4+nhp*nhp                                           4d7s22
        call enough('derh01.  3',bc,ibc)
        do i=0,nh*nh-1                                                  5d4s22
         bc(iatmp1+i)=h0(ioffmo+i)                                      5d4s22
        end do                                                          4d7s22
        if(idwsdeb.ne.0)then                                            7d14s22
         write(6,*)('h0 in mo basis ')
         call prntm2(bc(iatmp1),nh,nh,nh)
        end if                                                          7d14s22
        isk=multh(isb,ipsym)                                            4d7s22
        nk=nbasdws(isk)                                                 4d7s22
        nkp=nbasisp(isk)*ncomp                                          4d7s22
        if(min(nh,nk).gt.0)then                                         7d11s22
         itmpx=ibcoff                                                    4d7s22
         itmpy=itmpx+nhp*nkp                                               4d7s22
         ibcoff=itmpy+nhp*nkp                                              4d7s22
         call enough('derh01.  4',bc,ibc)
         if(idwsdeb.ne.0)then                                           7d14s22
          write(6,*)('multiply atmp1 ')
          call prntm2(bc(iatmp1),nh,nh,nh)
          write(6,*)('times trans ')
          call prntm2(bc(itrans(isb)),nh,nk,nh)
         end if                                                         7d14s22
         call dgemm('n','n',nh,nk,nh,1d0,bc(iatmp1),nh,bc(itrans(isb)),  4d7s22
     $        nh,0d0,bc(itmpx),nh,
     $       'derh01.3',                                                7d11s22
     d'derh01.  3')
         if(idwsdeb.ne.0)then                                           7d14s22
          write(6,*)('h0*t ')
          call prntm2(bc(itmpx),nh,nk,nh)
         end if                                                         7d14s22
         do i=0,nk-1                                                     6d8s22
          do j=0,nh-1                                                    6d8s22
           ji=itmpx+j+nh*i                                               6d8s22
           iadb=iptk(isb)+j+nh*i                                         6d8s22
           iadk=iptk(isk)+i+nk*j                                         6d8s22
           bc(iadb)=bc(iadb)+bc(ji)                                      6d8s22
           bc(iadk)=bc(iadk)+bc(ji)                                      6d8s22
          end do                                                         6d8s22
         end do                                                          6d8s22
         if(idwsdeb.ne.0)then                                           7d14s22
          write(6,*)('isk,isb: '),isk,isb
          write(6,*)('derb so far: ')
          call prntm2(bc(iptk(isb)),nh,nk,nh)
          write(6,*)('derk so far: ')
          call prntm2(bc(iptk(isk)),nk,nh,nk)
         end if                                                         7d14s22
         if(isk.gt.isb)then                                              4d7s22
          ioffd=isou(isb)+1                                              4d19s22
          if(idwsdeb.ne.0)then                                          7d14s22
           write(6,*)('starting h0d '),ioffd
           call prntm2(h0d(ioffd),nhp,nkp,nhp)
          end if                                                        7d14s22
          call dgemm('n','n',nhp,nk,nkp,1d0,h0d(ioffd),nhp,              4d7s22
     $        bc(iorb(isk)),nkp,0d0,bc(itmpx),nhp,
     $        'derh01.4',                                                              4d7s22
     d'derh01.  4')
          if(idwsdeb.ne.0)then                                          7d14s22
           write(6,*)('transformed from the right ')
           call prntm2(bc(itmpx),nhp,nk,nhp)
          end if                                                        7d14s22
          do i=0,nk-1                                                     4d7s22
           do j=0,nhp-1                                                    4d7s22
            ji=itmpx+j+nhp*i                                             4d7s22
            ij=itmpy+i+nk*j                                               4d7s22
            bc(ij)=bc(ji)                                                 4d7s22
           end do                                                         4d7s22
          end do                                                          4d7s22
          call dgemm('n','n',nk,nh,nhp,1d0,bc(itmpy),nk,bc(iorb(isb)),   4d7s22
     $        nhp,0d0,bc(itmpx),nk,'derh01.5',                                     4d7s22
     d'derh01.  5')
          if(idwsdeb.ne.0)then                                          7d14s22
           write(6,*)('dh0 in mo basis ')
           call prntm2(bc(itmpx),nk,nh,nk)
          end if                                                        7d14s22
          do i=0,nk-1                                                    4d7s22
           do j=0,nh-1                                                   4d7s22
            ji=iptk(isb)+j+nh*i                                          6d8s22
            ij2=iptk(isk)+i+nk*j                                         6d8s22
            ij=itmpx+i+nk*j                                              4d7s22
            bc(ji)=bc(ji)+bc(ij)                                         4d7s22
            bc(ij2)=bc(ij2)+bc(ij)                                       6d8s22
           end do                                                        4d7s22
          end do                                                         4d7s22
          if(idwsdeb.ne.0)then                                          7d14s22
           write(6,*)('derb so far '),iptk(isb)
           call prntm2(bc(iptk(isb)),nh,nk,nh)
           write(6,*)('derk so far '),iptk(isk)
           call prntm2(bc(iptk(isk)),nk,nh,nk)
          end if                                                        7d14s22
         end if                                                         7d11s22
        end if                                                          4d7s22
        ibcoff=itmpx                                                    4d7s22
        ioff=ioff+nhp*nhp                                               4d15s22
        ioffmo=ioffmo+nh*nh                                             5d4s22
        ibcoff=iatmp1                                                   4d7s22
       end do                                                           4d7s22
      end if
      return
      end
c mpec2.1 version zeta copyright u.s. government
c mpec2.1 version zeta copyright u.s. government
      subroutine printao(xmat,bc,ibc)
      implicit real*8 (a-h,o-z)
      common/lowersymcm/nsymbgx,iptno(8),ipts(8),nhsz(8),ipao(8)        4d25s18
      dimension ixto(21,4),nx(4),ps(21,4)
      data nx/21,11,6,2/
      data (ixto(j,1),j=1,21)/1,2,3,4,5,6,7,8,9,12,15,18,21,22,23,
     $     27,28,29,30,37,38/
      data (ixto(j,2),j=1,11)/10,13,16,19,25,31,32,33,34,35,40/
      data (ixto(j,3),j=1,6)/11,14,17,20,26,36/
      data (ixto(j,4),j=1,2)/24,39/
      data ps/84*1d0/
      include "common.store"
      include "common.basis"
      include "common.input"
      dimension xmat(*)
      write(6,*)('in printao')
      ps(21,1)=-1d0
      ps(6,2)=-1d0
      ps(7,2)=-1d0
      ps(8,2)=-1d0
      ps(9,2)=-1d0
      ps(11,2)=-1d0
      ps(2,4)=-1d0
      ps(4,4)=-1d0
      srh=sqrt(0.5d0)
      ncomp=1
      if(idorel.ne.0)ncomp=2
      nb=40*ncomp
      itmp1=ibcoff
      itmp2=itmp1+nb*nb
      ibcoff=itmp2+nb*nb
      jtmp1=itmp1-1
      do i=1,nb*nb
       bc(jtmp1+i)=xmat(i)
      end do
      call prntm2(bc(itmp1),nb,nb,nb)
      jgoal=5289449
      idelta=jgoal-itmp1
      i2=idelta/nb
      i1=idelta-nb*i2
      write(6,*)('for '),jgoal,idelta,i1,i2
      jgoal=5289449
      do ipass=1,2
       write(6,*)('for pass '),ipass,srh
       do i=1,26
        ip=i+40
        write(6,*)('move cols '),i
        do j=0,nb-1
         iad=j+nb*(i-1)
         iadp=j+nb*(ip-1)
         bc(itmp2+iad)=bc(itmp1+iad)
         if(idorel.ne.0)then
          iad=j+nb*(i+40-1)
          iadp=j+nb*(ip+40-1)
          bc(itmp2+iad)=bc(itmp1+iad)
         end if
        end do
       end do
       do i=27,30
        ip=i+4
        write(6,*)('combine cols '),i,ip
        do j=0,nb-1
         iad=j+nb*(i-1)
         iadp=j+nb*(ip-1)
         bc(itmp2+iad)=srh*(bc(itmp1+iad)+bc(itmp1+iadp))
         bc(itmp2+iadp)=srh*(bc(itmp1+iad)-bc(itmp1+iadp))
         if(idorel.ne.0)then
          iad=j+nb*(i+40-1)
          iadp=j+nb*(ip+40-1)
          bc(itmp2+iad)=srh*(bc(itmp1+iad)+bc(itmp1+iadp))
          bc(itmp2+iadp)=srh*(bc(itmp1+iad)-bc(itmp1+iadp))
         end if
        end do
       end do
       do i=35,37
        ip=i+3
        write(6,*)('combine cols '),i,ip
        do j=0,nb-1
         iad=j+nb*(i-1)
         iadp=j+nb*(ip-1)
         bc(itmp2+iad)=srh*(bc(itmp1+iad)+bc(itmp1+iadp))
         bc(itmp2+iadp)=srh*(bc(itmp1+iad)-bc(itmp1+iadp))
         if(idorel.ne.0)then
          iad=j+nb*(i+40-1)
          iadp=j+nb*(ip+40-1)
          bc(itmp2+iad)=srh*(bc(itmp1+iad)+bc(itmp1+iadp))
          bc(itmp2+iadp)=srh*(bc(itmp1+iad)-bc(itmp1+iadp))
         end if
        end do
       end do
       write(6,*)('to yield ')
       call prntm2(bc(itmp2),nb,nb,nb)
       do i=0,nb-1
        do j=0,nb-1
         ji=itmp2+j+nb*i
         ij=itmp1+i+nb*j
         bc(ij)=bc(ji)
        end do
       end do
      end do
      write(6,*)('symmetrized ')
      call prntm2(bc(itmp1),nb,nb,nb)
      ibcoff=itmp2
      do isb=1,4
       nxb=nx(isb)*ncomp
       do isk=1,4
        nxk=nx(isk)*ncomp
        itmp=ibcoff
        igoal=itmp+18-1+42*(18-1)
        ibcoff=itmp+nxb*nxk
        rms=0d0
        do ib=0,ncomp-1
         do ik=0,ncomp-1
          do lb=1,nx(isb)
           ioffb=ixto(lb,isb)+40*ib
           ifb=lb+nx(isb)*ib
           do lk=1,nx(isk)
            ioffk=ixto(lk,isk)+40*ik
            ifk=lk+nx(isk)*ik
            iad1=itmp+ifb-1+nxb*(ifk-1)
            iad2=itmp1+ioffb-1+nb*(ioffk-1)
            bc(iad1)=bc(iad2)*ps(lk,isk)*ps(lb,isb)
            rms=rms+bc(iad2)**2
           end do
          end do
         end do
        end do
        rms=sqrt(rms/dfloat(nxb*nxk))
        if(rms.gt.1d-10)then
         write(6,*)('for symmetries '),isb,isk
         write(6,*)('times 2')
         do i=0,nxb*nxk-1
          bc(itmp+i)=bc(itmp+i)*2d0
         end do
         call prntm2(bc(itmp),nxb,nxk,nxb)
        end if
        ibcoff=itmp
       end do
      end do
      ibcoff=itmp1
      return
      end
