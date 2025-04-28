c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine derh0(h0,h0d,ivecs,ixor,morb,idorel,ipsym,multh,ider,  4d7s22
     $     noc,iorb,idwsdeb,dere1,der2e1,h0d2,ider2,itrans,itran2,      4d7s22
     $     nbasisp,isou,bc,ibc)                                         11d14s22
      implicit real*8 (a-h,o-z)
      include "common.hf"
      include "common.store"
      include "common.spher"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension ivecs(8),ieigs(8),iorb(8),morb(8),h0(*),h0d(*),
     $    multh(8,8),noc(*),h0d2(*),ipt2(8),iptk(8),itrans(*),itran2(*),4d7s22
     $     nbasisp(*),isou(*)                                           4d19s22
      if(idwsdeb.gt.10)write(6,*)('in derh0 ')
      igoal=2811778
      ncomp=1
      dere1=0d0                                                         3d21s16
      der2e1=0d0                                                        7d21s16
      srh=sqrt(0.5d0)
      if(idorel.ne.0)ncomp=2
      ider=ibcoff
      if(ipsym.eq.1)then
       jxor=ixor
       ioff=1
       jder=ider
       nwds=0                                                           12d5s16
       do isb=1,nsymb                                                   12d5s16
        nwds=nwds+(nbasdws(isb)*ncomp)**2                               12d5s16
       end do                                                           12d5s16
       ider2=ider+nwds                                                  12d5s16
       ibcoff=ider2+nwds                                                12d5s16
       jder2=ider2                                                      12d5s16
       ibcoffo=ibcoff                                                   12d15s16
       do isb=1,nsymb
        nh=nbasdws(isb)                                                 4d7s22
        nhp=nbasisp(isb)*ncomp                                           4d7s22
        if(idwsdeb.gt.10)then
         write(6,*)('h0 in ao basis '),ioff
         call prntm2(h0(ioff),nhp,nhp,nhp)
        end if
        if(idwsdeb.gt.10)then
         write(6,*)('h0d before squaring: '),bc(igoal)
         call mpprnt2(h0d(ioff),nhp)
         write(6,*)('h0d2 before squaring: ')                            7d1s16
         call mpprnt2(h0d2(ioff),nhp)                                     7d1s16
        end if
        call square(h0d(ioff),nhp)
        call square(h0d2(ioff),nhp)                                      7d1s16
        iatmp1=ibcoff                                                   4d6s22
        iatmp2=iatmp1+nhp*nhp                                            4d6s22
        iatmp3=iatmp2+nhp*nhp                                            4d6s22
        iatmp4=iatmp3+nhp*nhp                                            4d6s22
        iatmp5=iatmp4+nhp*nhp                                            4d6s22
        iatmp6=iatmp5+nhp*nhp                                            4d6s22
        ibcoff=iatmp6+nhp*nhp                                            4d6s22
        call enough('derh0.  1',bc,ibc)
        do ii=0,nhp*nhp-1                                                4d6s22
         bc(iatmp1+ii)=h0(ioff+ii)                                      4d6s22
         bc(iatmp3+ii)=h0d(ioff+ii)                                     4d6s22
         bc(iatmp5+ii)=h0d2(ioff+ii)                                    4d6s22
        end do                                                          4d6s22
        nrow=nhp                                                        4d7s22
        do ipass=1,2                                                    4d6s22
         call dgemm('n','n',nrow,nh,nhp,1d0,bc(iatmp1),nrow,            4d7s22
     $         bc(iorb(isb)),nhp,0d0,bc(iatmp2),nrow,                   4d7s22
     d'derh0.  1')
         call dgemm('n','n',nrow,nh,nhp,1d0,bc(iatmp3),nrow,            4d7s22
     $        bc(iorb(isb)),nhp,0d0,bc(iatmp4),nrow,                    4d7s22
     d'derh0.  2')
         call dgemm('n','n',nrow,nh,nhp,1d0,bc(iatmp5),nrow,            4d7s22
     $        bc(iorb(isb)),nhp,0d0,bc(iatmp6),nrow,                    4d7s22
     d'derh0.  3')
         do i=0,nh-1                                                    4d6s22
          do j=0,nrow-1                                                  4d6s22
           ji=iatmp2+j+nrow*i                                            4d6s22
           ij=iatmp1+i+nh*j                                             4d6s22
           bc(ij)=bc(ji)                                                4d6s22
           ji=iatmp4+j+nrow*i                                            4d6s22
           ij=iatmp3+i+nh*j                                             4d6s22
           bc(ij)=bc(ji)                                                4d6s22
           ji=iatmp6+j+nrow*i                                            4d6s22
           ij=iatmp5+i+nh*j                                             4d6s22
           bc(ij)=bc(ji)                                                4d6s22
          end do                                                        4d6s22
         end do                                                         4d6s22
         nrow=nh                                                        4d7s22
        end do                                                          4d6s22
        call dgemm('n','n',nh,nh,nh,1d0,bc(iatmp1),nh,bc(itrans(isb)),  4d7s22
     $        nh,0d0,bc(iatmp2),nh,                                     4d6s22
     d'derh0.  4')
        do i=0,nh-1                                                     4d6s22
         do j=0,nh-1                                                    4d6s22
          ji=iatmp2+j+nh*i                                              4d6s22
          ij=iatmp4+i+nh*j                                              4d6s22
          bc(ij)=bc(ji)                                                 4d6s22
         end do                                                         4d6s22
        end do                                                          4d6s22
        call dgemm('n','n',nh,nh,nh,2d0,bc(iatmp4),nh,bc(itrans(isb)),  4d7s22
     $        nh,1d0,bc(iatmp5),nh,                                     4d6s22
     d'derh0.  5')
        call dgemm('n','n',nh,nh,nh,2d0,bc(iatmp4),nh,bc(itrans(isb)),  4d7s22
     $        nh,0d0,bc(ibcoff),nh,                                     4d6s22
     d'derh0.  6')
        call dgemm('n','n',nh,nh,nh,2d0,bc(iatmp3),nh,bc(itrans(isb)),  4d6s22
     $        nh,0d0,bc(iatmp4),nh,                                     4d6s22
     d'derh0.  7')
        call dgemm('n','n',nh,nh,nh,1d0,bc(iatmp1),nh,bc(itran2(isb)),  4d6s22
     $        nh,0d0,bc(iatmp6),nh,                                     4d6s22
     d'derh0.  8')
        do i=0,nh-1                                                     4d6s22
         do j=0,nh-1                                                    4d6s22
          ji=j+nh*i                                                     4d6s22
          ij=i+nh*j                                                     4d6s22
          bc(iatmp1+ji)=bc(iatmp2+ji)+bc(iatmp2+ij)+bc(iatmp3+ji)       4d6s22
          bc(ibcoff+ji)=bc(iatmp5+ji)+bc(iatmp6+ji)+bc(iatmp6+ij)
          bc(iatmp5+ji)=bc(iatmp5+ji)+bc(iatmp4+ji)+bc(iatmp4+ij)       4d6s22
     $         +bc(iatmp6+ji)+bc(iatmp6+ij)                              4d6s22
         end do                                                         4d6s22
        end do                                                          4d6s22
        do i=0,nh*nh-1
         bc(jder+i)=bc(iatmp1+i)                                         4d6s22
         bc(jder2+i)=bc(iatmp5+i)                                        4d7s22
        end do
        do i=0,noc(isb)-1                                               4d7s22
         kder=jder+i*(nh+1)                                             4d7s22
         orig=dere1
         dere1=dere1+2d0*bc(kder)                                       4d7s22
        end do                                                          4d7s22
        ibcoff=iatmp1                                                   4d6s22
        jder=jder+nh*nh                                                 12d5s16
        jder2=jder2+nh*nh                                               12d5s16
        ioff=ioff+nhp*nhp                                               4d15s22
        ibcoff=ibcoffo                                                  3d22s16
       end do
      else
c
c     perturbation is not a1
c
       jxor=ixor
       jder=ider
       n2ndsz=0                                                          6d21s16
       n2ndszp=0                                                        4d7s22
       nx=0                                                             7d5s16
       do isb=1,nsymb                                                    6d21s16
        isk=multh(isb,ipsym)                                            7d19s16
        iptk(isb)=nx                                                    7d19s16
        nx=nx+nbasdws(isb)**2                                           4d7s22
        if(isk.gt.isb)then                                              7d5s16
         ipt2(isb)=n2ndszp                                              4d7s22
         ipt2(isk)=n2ndszp                                              4d7s22
         n2ndsz=n2ndsz+nbasdws(isk)*nbasdws(isb)                        4d7s22
         n2ndszp=n2ndszp+nbasisp(isk)*nbasisp(isb)*ncomp*ncomp          4d7s22
        end if                                                          7d5s16
       end do                                                            6d21s16
       ider2=ider+n2ndsz                                                12d5s16
       ibcoff=ider2+nx                                                  12d16s16
       idera=ider
       ider2a=ider2
       call enough('derh0.  2',bc,ibc)
       do iz=idera,ibcoff-1                                             4d7s22
        bc(iz)=0d0                                                      4d7s22
       end do                                                           4d7s22
       ioff=1                                                           4d7s22
       jderab=idera                                                     4d7s22
       ioffd=1                                                          4d7s22
       jder2a=ider2a                                                    4d7s22
       do isb=1,nsymb                                                   4d7s22
        nh=nbasdws(isb)                                                 4d7s22
        nhp=nbasisp(isb)*ncomp                                          4d7s22
        iatmp1=ibcoff                                                   4d7s22
        iatmp2=iatmp1+nhp*nhp                                           4d7s22
        iatmp3=iatmp2+nhp*nhp                                           4d7s22
        iatmp4=iatmp3+nhp*nhp                                           4d7s22
        ibcoff=iatmp4+nhp*nhp                                           4d7s22
        call enough('derh0.  3',bc,ibc)
        do i=0,nhp*nhp-1                                                4d7s22
         bc(iatmp1+i)=h0(ioff+i)                                        4d7s22
         bc(iatmp3+i)=h0d2(ioff+i)                                      4d7s22
        end do                                                          4d7s22
        call square(bc(iatmp3),nhp)                                     4d7s22
        nrow=nhp                                                        4d7s22
        if(min(nh,nhp).gt.0)then                                        11d28s22
         do ipass=1,2                                                    4d7s22
          call dgemm('n','n',nrow,nh,nhp,1d0,bc(iatmp1),nrow,            4d7s22
     $        bc(iorb(isb)),nhp,0d0,bc(iatmp2),nrow,                    4d7s22
     d'derh0.  9')
          do i=0,nh-1                                                    4d7s22
           do j=0,nrow-1                                                 4d7s22
            ji=iatmp2+j+nrow*i                                           4d7s22
            ij=iatmp1+i+nh*j                                             4d7s22
            bc(ij)=bc(ji)                                                4d7s22
           end do                                                        4d7s22
          end do                                                         4d7s22
          call dgemm('n','n',nrow,nh,nhp,1d0,bc(iatmp3),nrow,            4d7s22
     $        bc(iorb(isb)),nhp,0d0,bc(iatmp4),nrow,                    4d7s22
     d'derh0. 10')
          do i=0,nh-1                                                    4d7s22
           do j=0,nrow-1                                                 4d7s22
            ji=iatmp4+j+nrow*i                                           4d7s22
            ij=iatmp3+i+nh*j                                             4d7s22
            bc(ij)=bc(ji)                                                4d7s22
           end do                                                        4d7s22
          end do                                                         4d7s22
          nrow=nh                                                        4d7s22
         end do                                                          4d7s22
        end if                                                          11d28s22
        if(nh.gt.0)then                                                 11d28s22
         call dgemm('n','n',nh,nh,nh,1d0,bc(iatmp1),nh,bc(itran2(isb)),  4d7s22
     $        nh,0d0,bc(iatmp4),nh,                                     4d7s22
     d'derh0. 11')
         do i=0,nh-1                                                     4d7s22
          do j=0,nh-1                                                    4d7s22
           ji=j+nh*i                                                     4d7s22
           ij=i+nh*j                                                     4d7s22
           orig=bc(jder2a+ji)
           bc(jder2a+ji)=bc(jder2a+ji)+bc(iatmp4+ji)+bc(iatmp4+ij)       4d7s22
     $         +bc(iatmp3+ji)                                           4d7s22
          end do                                                         4d7s22
         end do                                                          4d7s22
        end if                                                          11d28s22
        isk=multh(isb,ipsym)                                            4d7s22
        nk=nbasdws(isk)                                                 4d7s22
        nkp=nbasisp(isk)*ncomp                                          4d7s22
        jderak=idera                                                    4d7s22
        kder2a=ider2a                                                   4d7s22
        do i=1,isk-1                                                    4d7s22
         j=multh(i,ipsym)                                               4d7s22
         if(j.gt.i)jderak=jderak+nbasdws(i)*nbasdws(j)                  4d7s22
         kder2a=kder2a+nbasdws(i)*nbasdws(i)                            4d7s22
        end do                                                          4d7s22
        itmpx=ibcoff                                                    4d7s22
        itmpy=itmpx+nhp*nkp                                               4d7s22
        ibcoff=itmpy+nhp*nkp                                              4d7s22
        call enough('derh0.  4',bc,ibc)
        if(min(nh,nk).gt.0)then                                         11d28s22
         call dgemm('n','n',nh,nk,nh,1d0,bc(iatmp1),nh,bc(itrans(isb)),  4d7s22
     $        nh,0d0,bc(itmpx),nh,                                      4d7s22
     d'derh0. 12')
        end if                                                          11d28s22
        if(isk.gt.isb)then                                              4d7s22
         do i=0,nh-1                                                     4d7s22
          do j=0,nk-1                                                    4d7s22
           ij=i+nh*j                                                     4d7s22
           ji=j+nk*i                                                     4d7s22
           bc(itmpy+ji)=bc(itmpx+ij)                                    4d7s22
           bc(jderab+ij)=bc(jderab+ij)+bc(itmpx+ij)                      4d7s22
          end do                                                         4d7s22
         end do                                                          4d7s22
        else                                                            4d7s22
         do i=0,nh-1                                                     4d7s22
          do j=0,nk-1                                                    4d7s22
           ij=i+nh*j                                                     4d7s22
           ji=j+nk*i                                                     4d7s22
           bc(itmpy+ji)=bc(itmpx+ij)                                    4d7s22
           bc(jderak+ji)=bc(jderak+ji)+bc(itmpx+ij)                      4d7s22
          end do                                                         4d7s22
         end do                                                          4d7s22
        end if                                                          4d7s22
        if(min(nk,nh).gt.0)then                                         11d28s22
         call dgemm('n','n',nk,nk,nh,2d0,bc(itmpy),nk,bc(itrans(isb)),   4d7s22
     $       nh,0d0,bc(ibcoff),nk,                                      4d7s22
     d'derh0. 13')
         call dgemm('n','n',nk,nk,nh,2d0,bc(itmpy),nk,bc(itrans(isb)),   4d7s22
     $       nh,1d0,bc(kder2a),nk,                                      4d7s22
     d'derh0. 14')
        end if                                                          11d28s22
        if(isk.gt.isb)then                                              4d7s22
         ioffd=isou(isb)+1                                              4d19s22
         if(min(nhp,nk).gt.0)then                                       11d28s22
          call dgemm('n','n',nhp,nk,nkp,1d0,h0d(ioffd),nhp,              4d7s22
     $        bc(iorb(isk)),nkp,0d0,bc(itmpx),nhp,                      4d7s22
     d'derh0. 15')
          do i=0,nk-1                                                     4d7s22
           do j=0,nhp-1                                                    4d7s22
            ji=itmpx+j+nhp*i                                             4d7s22
            ij=itmpy+i+nk*j                                               4d7s22
            bc(ij)=bc(ji)                                                 4d7s22
           end do                                                         4d7s22
          end do                                                          4d7s22
          call dgemm('n','n',nk,nh,nhp,1d0,bc(itmpy),nk,bc(iorb(isb)),   4d7s22
     $        nhp,0d0,bc(itmpx),nk,                                     4d7s22
     d'derh0. 16')
         end if                                                         11d28s22
         do i=0,nk-1                                                    4d7s22
          do j=0,nh-1                                                   4d7s22
           ji=jderab+j+nh*i                                             4d7s22
           ij=itmpx+i+nk*j                                              4d7s22
           bc(ji)=bc(ji)+bc(ij)                                         4d7s22
          end do                                                        4d7s22
         end do                                                         4d7s22
         jderab=jderab+nh*nk                                            4d7s22
        else                                                            4d7s22
         ioffd=isou(isk)+1                                              4d19s22
         if(min(nkp,nh,nhp).gt.0)then                                   11d28s22
          call dgemm('n','n',nkp,nh,nhp,1d0,h0d(ioffd),nkp,              4d7s22
     $        bc(iorb(isb)),nhp,0d0,bc(itmpy),nkp,                      4d7s22
     d'derh0. 17')
         else                                                           11d28s22
          do iz=itmpy,itmpy+nkp*nh-1                                    11d28s22
           bc(iz)=0d0                                                   11d28s22
          end do                                                        11d28s22
         end if                                                         11d28s22
         do i=0,nkp-1                                                    4d7s22
          do j=0,nh-1                                                   4d7s22
           ji=itmpx+j+nh*i                                              4d7s22
           ij=itmpy+i+nkp*j                                             4d7s22
           bc(ji)=bc(ij)                                                4d7s22
          end do                                                        4d7s22
         end do                                                         4d7s22
         if(min(nh,nk,nkp).gt.0)then                                    11d28s22
          call dgemm('n','n',nh,nk,nkp,1d0,bc(itmpx),nh,bc(iorb(isk)),   4d7s22
     $        nkp,0d0,bc(itmpy),nh,                                     4d7s22
     d'derh0. 18')
         else                                                           11d28s22
          do iz=itmpy,itmpy+nh*nk-1                                     11d28s22
           bc(iz)=0d0                                                   11d28s22
          end do                                                        11d28s22
         end if                                                         11d28s22
         do i=0,nk-1                                                    4d7s22
          do j=0,nh-1                                                   4d7s22
           ji=itmpy+j+nh*i                                              4d7s22
           ij=itmpx+i+nk*j                                              4d7s22
           bc(ij)=bc(ji)                                                4d7s22
          end do                                                        4d7s22
         end do                                                         4d7s22
        end if                                                          4d7s22
        itmpz=ibcoff                                                    4d7s22
        ibcoff=itmpz+nk*nk                                              4d7s22
        call enough('derh0.  5',bc,ibc)
        if(min(nk,nh).gt.0)then                                         11d28s22
         call dgemm('n','n',nk,nk,nh,2d0,bc(itmpx),nk,bc(itrans(isb)),   4d7s22
     $        nh,0d0,bc(itmpz),nk,                                       4d7s22
     d'derh0. 19')
         do i=0,nk-1                                                     4d7s22
          do j=0,nk-1                                                    4d7s22
           ji=j+nk*i                                                     4d7s22
           ij=i+nk*j                                                     4d7s22
           bc(kder2a+ji)=bc(kder2a+ji)+bc(itmpz+ji)+bc(itmpz+ij)         4d7s22
          end do                                                         4d7s22
         end do                                                          4d7s22
        end if                                                          11d28s22
        ibcoff=itmpx                                                    4d7s22
        ioff=ioff+nhp*nhp                                               4d15s22
        jder2a=jder2a+nh*nh                                             4d7s22
        ibcoff=iatmp1                                                   4d7s22
       end do                                                           4d7s22
      end if
      return
      end
