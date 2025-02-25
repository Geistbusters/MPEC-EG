c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine bodcpart(itrans,ovrdk,ovrdd,morb,iorb,idarot,nbasisp,  6d3s22
     $     nbasdws,idorel,bodc,iden1d,idoub,iact,noc,nsymb,ivecs,iden1, 6d20s22
     $     ipuse,multh,isoa1,j2den,nsdlk,isblk,isou,idwsdeb,lprint,bc,  11d14s22
     $     ibc,mden,nstate,isinfo,cdir,cdc,tracer,cdirb,ndentd)         4d10s23
      implicit real*8 (a-h,o-z)
      include "common.store"                                            6d3s22
      logical lprt,lprint                                               7d25s22
      dimension ovrdk(*),ovrdd(*),morb(*),idarot(*),nbasisp(*),
     $     nbasdws(*),iden1d(8,*),idoub(*),iact(*),noc(*),iorb(*),      6d3s22
     $     itrans(*),ivecs(*),iden1(8,*),multh(8,8),isoa1(*),isblk(4,*),3d17s23
     $     idnm(8),j2den(*),isou(*),bodc(*),isinfo(11,*),cdir(*),       4d10s23
     $     cdc(*),cdirb(*),idnmi(8)                                     4d10s23
      lprt=idwsdeb.ne.0                                                 7d14s22
      if(lprt)write(6,*)('Hi, my name is bodcpart ')                    11d28s22
      cipart=0d0                                                        6d20s22
      ovrpart=0d0                                                       6d20s22
      orbpart=0d0                                                       6d20s22
      ciorbpart=0d0
      ibcoffo=ibcoff                                                    6d20s22
      ncomp=1                                                           6d3s22
      if(idorel.ne.0)ncomp=2                                            6d3s22
      do i=1,mden                                                       3d14s23
       bodc(i)=0d0                                                      3d14s23
      end do                                                            3d14s23
      ibc0=ibcoff                                                       11d30s22
      do isb=1,nsymb                                                    11d30s22
       nh=nbasdws(isb)                                                  6d3s22
       jsb=multh(isb,ipuse)                                             6d20s22
       nhj=nbasdws(jsb)                                                 6d20s22
       idnm(jsb)=ibcoff                                                 6d20s22
       idnmi(jsb)=idnm(jsb)+nhj*nh                                      4d10s23
       ibcoff=idnmi(jsb)+nhj*nh                                         4d10s23
      end do                                                            11d30s22
      isuma=ibcoff                                                      3d14s23
      isum1=isuma+mden                                                  3d14s23
      isum2=isum1+mden                                                  3d14s23
      ibcoff=isum2+mden                                                 3d14s23
      call enough('bodcpart.  0.1',bc,ibc)
      do iz=ibc0,ibcoff-1                                               11d30s22
       bc(iz)=0d0                                                       11d30s22
      end do                                                            11d30s22
      do isb=1,nsymb
       nhp=nbasisp(isb)*ncomp                                           6d3s22
       nh=nbasdws(isb)                                                  6d3s22
       jsb=multh(isb,ipuse)                                             6d20s22
       nhj=nbasdws(jsb)                                                 6d20s22
       nhjp=nbasisp(jsb)*ncomp                                          6d20s22
       call enough('bodcpart.  1',bc,ibc)
       if(lprt)write(6,*)('for isb,jsb '),isb,jsb,nh,nhj                11d30s22
       if(nhj.gt.0)then                                                 11d30s22
        idndm=ibcoff                                                    6d6s22
        ibcoff=idndm+nhj*nhj                                             6d20s22
        iorbpart=ibcoff                                                 6d20s22
        ibcoff=iorbpart+nhj*nhj                                         6d20s22
        iovrpart=ibcoff                                                 6d20s22
        ibcoff=iovrpart+nhj*nhj                                         6d20s22
        itmp3=ibcoff                                                    11d30s22
        itmp4=itmp3+nhjp*nhj                                            11d30s22
        ibcoff=itmp4+nhjp*nhj                                           11d30s22
        call enough('bodcpart.  1.1',bc,ibc)
        do iz=idndm,ibcoff-1                                            11d30s22
         bc(iz)=0d0                                                     11d30s22
        end do                                                          11d30s22
        if(lprt)then                                                    11d30s22
         write(6,*)('ovrdd '),isoa1(jsb)+1
         call prntm2(ovrdd(isoa1(jsb)+1),nhjp,nhjp,nhjp)                 6d20s22
        end if                                                          11d30s22
        do i=0,nhjp-1                                                   6d20s22
         do j=0,i-1                                                     6d20s22
          ji=isoa1(jsb)+1+j+nhjp*i                                      6d20s22
          ij=isoa1(jsb)+1+i+nhjp*j                                      6d20s22
          ovrdd(ji)=ovrdd(ij)                                           6d20s22
         end do                                                         6d20s22
        end do                                                          6d20s22
        if(lprt)then                                                    11d30s22
         write(6,*)('transposed ')                                      8d20s24
         call prntm2(ovrdd(isoa1(jsb)+1),nhjp,nhjp,nhjp)                 6d20s22
        end if                                                          11d30s22
        call dgemm('n','n',nhjp,nhj,nhjp,1d0,ovrdd(isoa1(jsb)+1),nhjp,  6d20s22
     $       bc(iorb(jsb)),nhjp,0d0,bc(itmp3),nhjp,                     6d20s22
     d'bodcpart.  6')
        do i=0,nhj-1                                                    6d20s22
         do j=0,nhjp-1                                                  6d20s22
          ji=itmp3+j+nhjp*i                                             6d20s22
          ij=itmp4+i+nhj*j                                              6d20s22
          bc(ij)=bc(ji)                                                 6d6s22
         end do                                                         6d6s22
        end do                                                          6d6s22
        call dgemm('n','n',nhj,nhj,nhjp,1d0,bc(itmp4),nhj,              6d20s22
     $       bc(iorb(jsb)),nhjp,0d0,bc(itmp3),nhj,                      6d20s22
     d'bodcpart.  7')
        if(lprt)then                                                    7d14s22
         write(6,*)('explicit (dv|dv) explicit')
         call prntm2(bc(itmp3),nhj,nhj,nhj)                              6d20s22
        end if                                                          7d14s22
        do i=0,nhj*nhj-1                                                6d20s22
         bc(idndm+i)=bc(itmp3+i)                                        11d30s22
         bc(iovrpart+i)=bc(itmp3+i)                                     6d20s22
        end do                                                          6d6s22
        if(lprt)then                                                    7d14s22
         write(6,*)('dndm so far ')
         call prntm2(bc(idndm),nhj,nhj,nhj)                              6d20s22
        end if                                                          7d14s22
        ibcoff=itmp3                                                    11d30s22
       end if                                                           11d30s22
       if(min(nh,nhj).gt.0)then                                         7d11s22
        if(lprt)then                                                    7d14s22
         write(6,*)('for symmetry block '),isb,jsb
         write(6,*)('idarot ')
         call prntm2(bc(idarot(isb)),nh,nhj,nh)                          6d20s22
         write(6,*)('itrans ')
         call prntm2(bc(itrans(isb)),nh,nhj,nh)                          6d20s22
        end if                                                          7d14s22
        itmp1=ibcoff                                                    6d3s22
        ibcoff=itmp1+nh*nhj                                             6d20s22
        call enough('bodcpart.  2',bc,ibc)
        do i=0,nh*nhj-1                                                 6d20s22
         bc(itmp1+i)=bc(idarot(isb)+i)+bc(itrans(isb)+i)                6d3s22
        end do                                                          6d3s22
        if(lprt)then                                                    7d14s22
         write(6,*)('sum of idarot and trans '),ibcoff
         call prntm2(bc(itmp1),nh,nhj,nh)                                6d20s22
        end if                                                          7d14s22
        jvecs=ivecs(isb)+nhp*nbasdws(isb)                               1d2s23
        if(lprt)then                                                    7d14s22
         write(6,*)('original overlap? '),jvecs,ivecs(isb),isb,nhp,
     $       nbasdws(isb),nbasisp(isb),ncomp
         call prntm2(bc(jvecs),nhp,nhp,nhp)
        end if                                                          7d14s22
        ioff=isou(isb)+1                                                7d1s22
        if(lprt)then                                                    7d14s22
         write(6,*)('ovrdk '),ioff
         call prntm2(ovrdk(ioff),nhp,nhjp,nhp)                           6d20s22
        end if                                                          7d14s22
        if(lprt)then                                                    7d14s22
         write(6,*)('morb*idarot ')
         call dgemm('n','n',nh,nhj,nh,1d0,bc(morb(isb)),nh,
     $       bc(idarot(isb)),nh,0d0,bc(ibcoff),nh,
     d'bodcpart.  3')
         call prntm2(bc(ibcoff),nh,nhj,nh)
         write(6,*)('iorb '),iorb(isb),ibcoff
         write(6,*)('orb*(idarot+trans) ')
        end if                                                          7d14s22
        itmp2=ibcoff                                                     6d3s22
        itmp2i=itmp2+nhp*nhj                                            4d10s23
        ibcoff=itmp2i+nhp*nhj                                           4d10s23
        call enough('bodcpart.  5',bc,ibc)
        call dgemm('n','n',nhp,nhj,nh,1d0,bc(iorb(isb)),nhp,            6d20s22
     $       bc(itmp1),nh,0d0,bc(itmp2),nhp,                            6d3s22
     d'bodcpart.  4')
        call dgemm('n','n',nhp,nhj,nh,1d0,bc(iorb(isb)),nhp,            6d20s22
     $       bc(idarot(isb)),nh,0d0,bc(itmp2i),nhp,                     4d10s23
     d'bodcpart.  4i')                                                  4d10s23
        if(lprt)then                                                    7d14s22
         write(6,*)('this should be der of ao vectors ')
         call prntm2(bc(itmp2),nhp,nhj,nhp)                              6d20s22
         write(6,*)
     $        ('this should be implicit part of der of ao vectors ')
         call prntm2(bc(itmp2i),nhp,nhj,nhp)                              6d20s22
        end if                                                          7d14s22
        itmp3=ibcoff
        itmp3i=itmp3+max(nhp,nhjp)*nhj                                  4d10s23
        itmp4=itmp3i+max(nhp,nhjp)*nhj                                  4d10s23
        ibcoff=itmp4+nhj*max(nhp,nhjp)                                  6d20s22
        itmp5=ibcoff                                                    6d6s22
        itmp5i=itmp5+nhj*nhp                                            4d10s23
        ibcoff=itmp5i+nhj*nhp                                           4d10s23
        call enough('bodcpart.  6',bc,ibc)
        nrow=nhp
        jgrab=jvecs                                                     6d3s22
        do ipass=1,2
         call dgemm('n','n',nrow,nhj,nhp,1d0,bc(jgrab),nrow,bc(itmp2),   6d6s22
     $        nhp,0d0,bc(itmp3),nrow,                                   6d6s22
     d'bodcpart.  5')
         if(lprt)then                                                   7d14s22
          write(6,*)('half'),ibcoff
          call prntm2(bc(itmp3),nrow,nhj,nrow)
         end if                                                         7d14s22
         do i=0,nhj-1                                                    6d3s22
          do j=0,nrow-1                                                 6d3s22
           ji=itmp3+j+nrow*i                                            6d3s22
           ij=itmp4+i+nhj*j                                              6d3s22
           bc(ij)=bc(ji)                                                6d3s22
          end do                                                        6d3s22
         end do                                                         6d3s22
         if(lprt)then                                                   7d14s22
          write(6,*)('transposed '),ibcoff
          call prntm2(bc(itmp4),nhj,nrow,nhj)
         end if                                                         7d14s22
         nrow=nhj                                                        6d3s22
         jgrab=itmp4                                                    6d3s22
        end do
        if(lprt)then                                                    7d14s22
         write(6,*)('implicit (dv|dv) implicit ')                        6d6s22
         call prntm2(bc(itmp4),nhj,nhj,nhj)                                6d20s22
        end if                                                          7d14s22
        do i=0,nhj*nhj-1                                                 6d20s22
         bc(iorbpart+i)=bc(iorbpart+i)+bc(itmp4+i)                      11d30s22
         bc(idndm+i)=bc(idndm+i)+bc(itmp4+i)                            11d30s22
        end do                                                          6d6s22
        if(lprt)then                                                    7d14s22
         write(6,*)('dndm so far ')
         call prntm2(bc(idndm),nhj,nhj,nhj)                              6d20s22
         write(6,*)('transforming ovrdk '),ioff,isb,jsb
         call prntm2(ovrdk(ioff),nhp,nhjp,nhp)
         if(isb.eq.jsb)then
          ssym=0d0
          do i=0,nhjp-1
           do j=0,i-1
            ji=ioff+j+nhp*i
            ij=ioff+i+nhp*j
            ssym=ssym+(ovrdk(ji)+ovrdk(ij))**2
           end do
           ii=ioff+i+nhp*i
           ssym=ssym+ovrdk(ii)**2
          end do
          write(6,*)('skew symmetry test '),ssym
         end if
         write(6,*)('times vectors ')
         call prntm2(bc(iorb(jsb)),nhjp,nhj,nhjp)
        end if                                                          7d14s22
        call dgemm('n','n',nhp,nhj,nhjp,1d0,ovrdk(ioff),nhp,            6d20s22
     $       bc(iorb(jsb)),nhjp,0d0,bc(itmp3),nhp,                      6d20s22
     d'bodcpart.  8')
        if(lprt)then                                                    7d14s22
         write(6,*)('ovrdk times orbs')
         call prntm2(bc(itmp3),nhp,nhj,nhp)
        end if                                                          7d14s22
        call dgemm('n','n',nhp,nhj,nhp,1d0,bc(jvecs),nhp,bc(itmp2),nhp, 6d20s22
     $       0d0,bc(itmp5),nhp,                                         6d6s22
     d'bodcpart.  9')
        call dgemm('n','n',nhp,nhj,nhp,1d0,bc(jvecs),nhp,bc(itmp2i),nhp,4d10s23
     $       0d0,bc(itmp5i),nhp,                                        4d10s23
     d'bodcpart.  9')
        if(lprt)then                                                    7d14s22
         write(6,*)('overlap times ao vect ders ')
         call prntm2(bc(itmp5),nhp,nhj,nhp)
         write(6,*)('overlap times implicit ao vect ders ')
         call prntm2(bc(itmp5i),nhp,nhj,nhp)
        end if                                                          7d14s22
        do i=0,nhj*nhp-1                                                6d20s22
         bc(itmp5+i)=bc(itmp5+i)+bc(itmp3+i)                            6d6s22
        end do                                                          6d6s22
        if(lprt)then                                                    7d14s22
         write(6,*)('summed ')
         call prntm2(bc(itmp5),nhp,nhj,nhp)
        end if                                                          7d14s22
        do i=0,nhj-1                                                    6d20s22
         do j=0,nhp-1                                                   6d6s22
          ji=itmp3+j+nhp*i                                              6d6s22
          ij=itmp4+i+nhj*j                                              6d20s22
          bc(ij)=bc(ji)                                                 6d6s22
         end do                                                         6d6s22
        end do                                                          6d6s22
        if(lprt)then                                                    7d14s22
         write(6,*)('ovrdk*orbs transposed ')
         call prntm2(bc(itmp4),nhj,nhp,nhj)
        end if                                                          7d14s22
        call dgemm('n','n',nhj,nhj,nhp,1d0,bc(itmp4),nhj,bc(itmp2),nhp, 6d20s22
     $       0d0,bc(itmp3),nhj,                                         6d20s22
     d'bodcpart. 10')
        if(lprt)then                                                    7d14s22
         write(6,*)('implicit (dv|dv) explicit')
         call prntm2(bc(itmp3),nhj,nhj,nhj)                              6d20s22
        end if                                                          7d14s22
        do i=0,nhj-1                                                    6d20s22
         do j=0,nhj-1                                                   6d20s22
          ji=itmp3+j+nhj*i                                              6d20s22
          ij=itmp3+i+nhj*j                                              6d20s22
          iad=idndm+j+nhj*i                                             6d20s22
          bc(iad)=bc(iad)+bc(ji)+bc(ij)                                 6d6s22
         end do                                                         6d6s22
        end do                                                          6d6s22
        if(lprt)then                                                    7d14s22
         write(6,*)('total dndm ')
         call prntm2(bc(idndm),nhj,nhj,nhj)                              6d20s22
        end if                                                          7d14s22
        do i=0,nhj-1                                                    6d20s22
         do j=0,nhp-1                                                   6d6s22
          ji=itmp3+i+nhj*j                                              6d20s22
          ij=itmp5+j+nhp*i                                              6d6s22
          bc(ji)=bc(ij)                                                 6d6s22
          ji=itmp3i+i+nhj*j                                             4d10s23
          ij=itmp5i+j+nhp*i                                             4d10s23
          bc(ji)=bc(ij)                                                 6d6s22
         end do                                                         6d6s22
        end do                                                          6d6s22
        if(lprt)then
         write(6,*)('to get full (dn|m), multiply ')
         call prntm2(bc(itmp3),nhj,nhp,nhj)
         write(6,*)('times vectors ')
         call prntm2(bc(iorb(isb)),nhp,nh,nhp)
         write(6,*)('to get implicit (dn|m), multiply ')
         call prntm2(bc(itmp3i),nhj,nhp,nhj)
         write(6,*)('times vectors ')
         sum=0d0
         do i=0,nhp-1
          iad1=itmp3+nhj*i
          iad2=iorb(isb)+i
          term=bc(iad1)*bc(iad2)
          sum=sum+term
          if(abs(term).gt.1d-12)write(6,*)i,bc(iad1),bc(iad2),sum
         end do
        end if
        call dgemm('n','n',nhj,nh,nhp,1d0,bc(itmp3),nhj,bc(iorb(isb)),  6d20s22
     $       nhp,0d0,bc(idnm(jsb)),nhj,                                 6d20s22
     d'bodcpart. 11')
        call dgemm('n','n',nhj,nh,nhp,1d0,bc(itmp3i),nhj,bc(iorb(isb)), 4d10s23
     $       nhp,0d0,bc(idnmi(jsb)),nhj,                                4d10s23
     d'bodcpart. 11i')                                                  4d10s23
        if(lprt)then                                                    7d14s22
         write(6,*)('full (dn|m)')
         call prntm2(bc(idnm(jsb)),nhj,nh,nhj)                           6d20s22
         write(6,*)('full implicit (dn|m)')                             4d10s23
         call prntm2(bc(idnmi(jsb)),nhj,nh,nhj)                         4d10s23
         if(nsymb.eq.1)then
          call printa(bc(idnm(jsb)),nh,0,1,0,nh,0,1,0,bc(ibcoff))
         end if
         write(6,*)('half differentiated densities '),iden1d(isb,2),
     $        iden1d(jsb,3),1861676+2+11*3,bc(1861676+2+11*3)
         call prntm2(bc(iden1d(isb,2)),iact(isb),iact(jsb)*mden,
     $        iact(isb))
         call prntm2(bc(iden1d(jsb,3)),iact(jsb),iact(isb)*mden,
     $        iact(jsb))
         if(nsymb.eq.1)then
          nad=iact(jsb)*iact(jsb)
          iad=iden1d(1,3)
          do k=0,mden-1
           write(6,*)('for density '),k+1
           call printa(bc(iad),iact,idoub,1,0,iact,idoub,1,0,bc(ibcoff))
           iad=iad+nad
          end do
         end if
        end if                                                          7d14s22
        iadd=iden1(jsb,2)                                               3d20s23
        iaddc=1                                                         3d20s23
        do is=1,nstate-1                                                3d20s23
         do isp=is+1,nstate                                             3d20s23
          if(isinfo(2,is).eq.isinfo(2,isp).and.                         4d7s23
     $       isinfo(3,is).eq.isinfo(3,isp))then                         4d7s23
           ksb=multh(jsb,multh(isinfo(1,is),isinfo(1,isp)))              3d20s23
           iheadr=0                                                      3d31s23
           if(iact(ksb).gt.0)then                                        3d20s23
            nad=iact(jsb)*iact(ksb)                                       3d20s23
            do irp=1,isinfo(4,isp)                                        3d20s23
             do ir=1,isinfo(4,is)                                         3d20s23
              if(ksb.eq.isb)then                                          3d20s23
               sumd=0d0
               sumdi=0d0                                                4d10s23
               do ikk=0,iact(ksb)-1                                       3d20s23
                ikkp=ikk+idoub(ksb)                                       3d20s23
                do ijj=0,iact(jsb)-1                                      3d20s23
                 ijjp=ijj+idoub(jsb)                                      3d20s23
                 iad1=iadd+ijj+iact(jsb)*ikk                              3d20s23
                 iad2=idnm(jsb)+ijjp+nhj*ikkp                             3d20s23
                 iad2i=idnmi(jsb)+ijjp+nhj*ikkp                         4d10s23
                 sumd=sumd+bc(iad1)*bc(iad2)                              3d20s23
                 sumdi=sumdi+bc(iad1)*bc(iad2i)                         4d10s23
                end do                                                    3d20s23
               end do                                                     3d20s23
               sumx=sumd-sumdi                                          4d10s23
               cdirb(iaddc)=cdirb(iaddc)+sumdi                          4d10s23
               cdirb(iaddc+ndentd)=cdirb(iaddc+ndentd)+sumx             4d13s23
              end if                                                    4d7s23
              iaddc=iaddc+1                                              3d20s23
              iadd=iadd+nad                                               3d20s23
             end do                                                     4d7s23
            end do                                                       3d20s23
           else                                                         4d7s23
            iaddc=iaddc+isinfo(4,isp)*isinfo(4,is)                      4d7s23
           end if                                                       4d7s23
          end if                                                        4d7s23
         end do                                                         3d20s23
        end do                                                          3d20s23
        do k=0,mden-1                                                   3d14s23
         bc(isuma+k)=0d0                                                3d14s23
         bc(isum1+k)=0d0                                                3d14s23
         bc(isum2+k)=0d0                                                3d14s23
        end do                                                          3d14s23
        if(ipuse.eq.1)then                                              6d20s22
         igoal=isuma+15
         suma=0d0                                                        6d6s22
         sum1=0d0
         sum2=0d0
         sumov=0d0                                                      6d20s22
         sumob=0d0                                                      6d20s22
         do k=0,mden-1                                                  3d14s23
          do i=0,iact(isb)-1                                              6d6s22
           ip=i+idoub(isb)                                                6d6s22
           do j=0,iact(isb)-1                                             6d6s22
            jp=j+idoub(isb)                                               6d6s22
            iad1=idnm(isb)+ip+nh*jp                                      6d20s22
c     vdv
            iad2=iden1d(isb,2)+j+iact(isb)*(i+iact(isb)*k)              3d14s23
            bc(isum1+k)=bc(isum1+k)+bc(iad1)*bc(iad2)                   3d14s23
            iad4=idnm(isb)+jp+nh*ip                                      6d26s22
c     dv v
            iad3=iden1d(isb,3)+j+iact(isb)*(i+iact(isb)*k)              3d14s23
            bc(isum2+k)=bc(isum2+k)+bc(iad3)*bc(iad4)                   3d14s23
            bc(isuma+k)=bc(isuma+k)+bc(iad1)*bc(iad2)+bc(iad3)*bc(iad4) 3d14s23
           end do                                                         6d6s22
          end do                                                          6d6s22
          bc(isuma+k)=bc(isuma+k)*0.5d0                                 3d14s23
         end do
        else                                                            6d20s22
         igoal=isuma+1
         suma=0d0                                                        6d6s22
         sum1=0d0
         sum2=0d0                                                       6d20s22
         do k=0,mden-1                                                  3d14s23
          do i=0,iact(isb)-1                                              6d6s22
           ip=i+idoub(isb)                                                6d6s22
           do j=0,iact(jsb)-1                                             6d6s22
            jp=j+idoub(jsb)                                               6d6s22
            iad1=idnm(jsb)+jp+nhj*ip                                           6d6s22
            iad2=iden1d(isb,2)+i+iact(isb)*(j+iact(jsb)*k)              3d14s23
            bc(isum1+k)=bc(isum1+k)+bc(iad1)*bc(iad2)                   3d14s23
            iad3=iden1d(jsb,3)+j+iact(jsb)*(i+iact(isb)*k)              3d14s23
            bc(isum2+k)=bc(isum2+k)+bc(iad1)*bc(iad3)                   3d14s23
            bc(isuma+k)=bc(isuma+k)+bc(iad1)*bc(iad2)                   3d14s23
     $           +bc(iad3)*bc(iad1)                                     3d27s23
           end do                                                       3d14s23
          end do                                                         6d6s22
          bc(isuma+k)=bc(isuma+k)*0.5d0                                 3d27s23
         end do                                                          6d6s22
        end if                                                          6d20s22
        if(lprt)then                                                    7d14s22
         write(6,*)('suma,sum1,sum2: '),isuma
         call prntm2(bc(isuma),mden,3,mden)                             3d14s23
        end if                                                          7d14s22
        do k=1,mden                                                     3d14s23
         km=k-1                                                         3d14s23
         bodc(k)=bodc(k)+bc(isuma+km)                                   3d14s23
        end do
        ibcoff=itmp1
       end if                                                           7d11s22
       if(nhj.ne.0)then
        dddd=0d0                                                        3d15s23
        do j=0,idoub(jsb)-1                                             6d6s22
         jj=idndm+j*(nhj+1)                                             6d20s22
         dddd=dddd+bc(jj)                                               3d15s23
         jj=iorbpart+j*(nhj+1)                                          6d20s22
         orbpart=orbpart+bc(jj)
         jj=iovrpart+j*(nhj+1)                                          6d20s22
         ovrpart=ovrpart+bc(jj)                                         6d20s22
        end do                                                          6d6s22
        k=1                                                               3d15s23
        do i=1,nstate                                                     3d15s23
         do k1=1,isinfo(4,i)                                              3d15s23
          do k2=1,k1                                                      3d15s23
           if(k1.eq.k2)bodc(k)=bodc(k)+dddd                             3d15s23
           k=k+1                                                          3d15s23
          end do                                                          3d15s23
         end do                                                           3d15s23
        end do                                                            3d15s23
        if(lprt)then                                                    7d14s22
         write(6,*)('undifferentiated density ')                         6d6s22
         call prntm2(bc(iden1(jsb,1)),iact(jsb),iact(jsb)*mden,         3d17s23
     $        iact(jsb))                                                3d17s23
        end if                                                          7d14s22
        if(ipuse.eq.1)then                                              3d16s23
         do k=0,mden-1                                                   3d16s23
          kp=k+1                                                         3d16s23
          kpp=kp+mden                                                   4d10s23
          do i=0,iact(jsb)-1                                              6d20s22
           ip=i+idoub(jsb)                                                6d20s22
           do j=0,iact(jsb)-1                                             6d20s22
            jp=j+idoub(jsb)                                               6d20s22
            iad1=idnm(jsb)+jp+nhj*ip                                    3d16s23
            iad1i=idnmi(jsb)+jp+nhj*ip                                  4d10s23
            iad2=iden1(jsb,1)+j+iact(jsb)*(i+iact(jsb)*k)               3d17s23
            cdir(kp)=cdir(kp)+bc(iad1i)*bc(iad2)                         3d16s23
            cdir(kpp)=cdir(kpp)+(bc(iad1)-bc(iad1i))*bc(iad2)           4d10s23
           end do                                                       3d16s23
          end do                                                        3d16s23
         end do                                                          3d16s23
         if(lprt)then                                                   3d16s23
          write(6,*)('derivative coupling after dnm ')                  3d16s23
          call prntm2(cdir,mden,2,mden)                                 4d10s23
         end if                                                         3d16s23
         do k=0,mden-1                                                  3d16s23
          kp=k+1                                                        3d16s23
          sum3=0d0
          do i=0,iact(isb)-1                                            3d16s23
           iad1=iden1d(isb,3)+i+iact(isb)*(i+iact(isb)*k)               3d16s23
           orig=sum3
           sum3=sum3+bc(iad1)
          end do                                                        3d16s23
          cdir(kp)=cdir(kp)+(sum3/tracer)
         end do                                                         3d16s23
         if(lprt)then                                                   3d16s23
          write(6,*)('derivative coupling trace dv v ')                 3d16s23
          call prntm2(cdir,mden,2,mden)                                 4d10s23
         end if                                                         3d16s23
        else                                                            3d17s23
         do isbtd=1,nsymb                                               3d17s23
          itden=iden1(isbtd,2)                                          3d17s23
          do i=1,nstate                                                  3d17s23
           do ip=1,nstate                                                3d17s23
            jsbtd=multh(isbtd,multh(isinfo(1,i),isinfo(1,ip)))          3d17s23
            if(isbtd.eq.jsb.and.jsbtd.eq.isb)then                       3d17s23
             do irp=1,isinfo(4,ip)                                      3d17s23
              do ir=1,isinfo(4,i)                                       3d17s23
               sum=0d0                                                  3d17s23
               do io=1,iact(isb)                                        3d17s23
                iop=io+idoub(isb)-1                                     3d17s23
                do jo=1,iact(jsb)                                       3d17s23
                 jop=jo+idoub(jsb)-1                                    3d17s23
                 iad=idnm(jsb)+jop+nhj*iop                              3d17s23
                 sum=sum+bc(itden)*bc(iad)                              3d17s23
                 itden=itden+1                                          3d17s23
                end do                                                  3d17s23
               end do                                                   3d17s23
              end do                                                    3d17s23
             end do                                                     3d17s23
            else                                                        3d17s23
             itden=itden+iact(jsb)*iact(isb)*isinfo(4,ip)*isinfo(4,i)   3d17s23
            end if                                                      3d17s23
           end do                                                       3d17s23
          end do                                                        3d17s23
         end do                                                         3d17s23
        end if                                                          3d16s23
        do k=0,mden-1                                                   3d14s23
         kp=k+1                                                         3d14s23
         do i=0,iact(jsb)-1                                              6d20s22
          ip=i+idoub(jsb)                                                6d20s22
          do j=0,iact(jsb)-1                                             6d20s22
           jp=j+idoub(jsb)                                               6d20s22
           iad1=idndm+jp+nhj*ip                                           6d6s22
           iad2=iden1(jsb,1)+j+iact(jsb)*(i+iact(jsb)*k)                3d17s23
           bodc(kp)=bodc(kp)+0.5d0*bc(iad1)*bc(iad2)                    3d14s23
           iad1=iorbpart+jp+nhj*ip                                       6d20s22
           orbpart=orbpart+0.5d0*bc(iad1)*bc(iad2)                       6d20s22
           iad1=iovrpart+jp+nhj*ip                                       6d20s22
           ovrpart=ovrpart+0.5d0*bc(iad1)*bc(iad2)                       6d20s22
          end do                                                         6d6s22
         end do                                                          6d6s22
        end do                                                          3d14s23
        if(iact(jsb).gt.0)then                                           11d30s22
         if(lprt)then                                                    11d30s22
          write(6,*)('double differentiated densities for symmetry '),   11d30s22
     $      jsb,iden1d(jsb,4)                                                         11d30s22
          call prntm2(bc(iden1d(jsb,4)),iact(jsb),iact(jsb)*mden,       3d14s23
     $         iact(jsb))                                               3d14s23
         end if                                                          11d30s22
         do k=0,mden-1                                                  3d14s23
          kp=k+1                                                        3d14s23
          do i=0,iact(jsb)-1                                              11d30s22
           ii=iden1d(jsb,4)+i+iact(jsb)*(i+iact(jsb)*k)                 3d14s23
           cipart=cipart+bc(ii)                                          11d30s22
           bodc(kp)=bodc(kp)+bc(ii)*0.5d0                               3d14s23
          end do                                                         11d30s22
         end do                                                         3d14s23
        end if                                                          11d30s22
       end if                                                           11d30s22
      end do
      if(lprt)then
       write(6,*)('1e contribution to bodc: ')                          3d14s23
       call prntm2(bodc,1,mden,1)                                       3d14s23
       if(nsymb.eq.1)then                                               3d27s23
        bc(ibcoff)=bodc(1)
        bc(ibcoff+1)=bodc(16)
        bc(ibcoff+2)=bodc(21)
        bc(ibcoff+3)=bodc(3)
        bc(ibcoff+4)=bodc(6)
        bc(ibcoff+5)=bodc(13)
        bc(ibcoff+6)=bodc(15)
        bc(ibcoff+7)=bodc(10)
        write(6,*)('in C2v order: ')
        call prntm2(bc(ibcoff),1,8,1)
       end if
      end if
      sumk=0d0                                                          6d20s22
      sumko=0d0                                                         11d29s22
      do kk=0,mden-1                                                    3d15s23
       bc(isum2+kk)=0d0                                                 3d15s23
      end do                                                            3d15s23
      do isb=1,nsymb                                                    6d20s22
       jsb=multh(isb,ipuse)                                             6d20s22
       do i=0,idoub(isb)-1                                              6d20s22
        do j=0,idoub(jsb)-1                                             6d20s22
         ji=idnm(jsb)+j+nbasdws(jsb)*i                                  6d20s22
         sumk=sumk+bc(ji)**2                                            6d20s22
         if(lprt)write(6,*)('dd'),j,i,jsb,isb,bc(ji),sumk
        end do                                                          6d20s22
        do kk=0,mden-1                                                  3d15s23
         do k=0,iact(jsb)-1                                              11d29s22
          ki=idnm(jsb)+k+idoub(jsb)+nbasdws(jsb)*i                       11d29s22
          do l=0,iact(jsb)-1                                             11d29s22
           lk=iden1(jsb,1)+l+iact(jsb)*(k+iact(jsb)*kk)                 3d17s23
           li=idnm(jsb)+l+idoub(jsb)+nbasdws(jsb)*i                      11d29s22
           bc(isum2+kk)=bc(isum2+kk)+bc(ki)*bc(li)*bc(lk)               3d15s23
           if(lprt)write(6,*)('da'),l,k,i,kk,jsb,isb,bc(ki),bc(li),
     $          bc(lk),bc(isum2+kk)
          end do                                                        3d15s23
         end do                                                         11d29s22
        end do                                                          11d29s22
       end do                                                           6d20s22
      end do                                                            6d20s22
      if(lprt)then
       write(6,*)('sumk: '),sumk
       write(6,*)('sumko: ')
       call prntm2(bc(isum2),1,mden,1)
      end if
      do k=1,mden                                                       3d15s23
       km=k-1                                                           3d15s23
       bodc(k)=bodc(k)-bc(isum2+km)                                     3d15s23
      end do                                                            3d15s23
      k=1                                                               3d15s23
      do i=1,nstate                                                     3d15s23
       do k1=1,isinfo(4,i)                                              3d15s23
        do k2=1,k1                                                      3d15s23
         if(k1.eq.k2)bodc(k)=bodc(k)-sumk                               3d15s23
         k=k+1                                                          3d15s23
        end do                                                          3d15s23
       end do                                                           3d15s23
      end do                                                            3d15s23
      if(lprt)then
       write(6,*)('bodc so far')                                        3d15s23
       call prntm2(bodc,1,mden,1)                                       3d15s23
       if(nsymb.eq.1)then                                               3d27s23
        bc(ibcoff)=bodc(1)
        bc(ibcoff+1)=bodc(16)
        bc(ibcoff+2)=bodc(21)
        bc(ibcoff+3)=bodc(3)
        bc(ibcoff+4)=bodc(6)
        bc(ibcoff+5)=bodc(13)
        bc(ibcoff+6)=bodc(15)
        bc(ibcoff+7)=bodc(10)
        write(6,*)('in C2v order: ')
        call prntm2(bc(ibcoff),1,8,1)
       end if
      end if
      if(lprt)write(6,*)('2e densities '),nsdlk                         7d14s22
      sum=0d0                                                           6d20s22
      do k=0,mden-1                                                     3d14s23
       bc(isum2+k)=0d0                                                  3d14s23
      end do                                                            3d14s23
      if(nsymb.eq.1)then
       kkk=15
      else
       kkk=1
      end if
      do is=1,nsdlk                                                     6d20s22
       if(multh(isblk(1,is),isblk(2,is)).eq.ipuse)then                  7d8s22
        nrow=iact(isblk(1,is))*iact(isblk(2,is))                        7d8s22
        ncol=iact(isblk(3,is))*iact(isblk(4,is))                        7d8s22
        if(min(nrow,ncol).gt.0)then                                     7d8s22
         if(lprt)then                                                   7d14s22
          write(6,*)('using density '),(isblk(j,is),j=1,4),j2den(is)
          call prntm2(bc(j2den(is)),nrow,ncol*mden,nrow)
          iad=j2den(is)+nrow*ncol*kkk
          write(6,*)('part we are using '),iad
          call prntm2(bc(iad),nrow,ncol,nrow)
         end if                                                         7d14s22
c     Gijkl*(di|j)*(k|dl)=-Gijkl*(di|j)*(dl|k)
         kden=j2den(is)                                                 7d8s22
         do k=0,mden-1                                                  3d27s23
          do i4=0,iact(isblk(4,is))-1                                    7d8s22
           i4p=i4+idoub(isblk(4,is))                                     7d8s22
           do i3=0,iact(isblk(3,is))-1                                   7d8s22
            i3p=i3+idoub(isblk(3,is))                                    7d8s22
            iad43=idnm(isblk(4,is))+i4p+nbasdws(isblk(4,is))*i3p           7d8s22
            do i2=0,iact(isblk(2,is))-1                                  7d8s22
             i2p=i2+idoub(isblk(2,is))                                   7d8s22
             iad12=idnm(isblk(1,is))+idoub(isblk(1,is))
     $           +nbasdws(isblk(1,is))*i2p                              7d8s22
             do i1=0,iact(isblk(1,is))-1                                 7d8s22
              term=bc(kden+i1)*bc(iad12+i1)*bc(iad43)                    7d8s22
              bc(isum2+k)=bc(isum2+k)-bc(kden+i1)*bc(iad12+i1)*bc(iad43)3d14s23
             end do                                                      7d8s22
             kden=kden+iact(isblk(1,is))                                 7d8s22
            end do                                                       7d8s22
           end do                                                        7d8s22
          end do                                                         7d8s22
         end do                                                         3d14s23
        end if                                                          7d8s22
       end if                                                           7d8s22
      end do                                                            6d20s22
      do k=0,mden-1                                                     3d14s23
       bc(isum2+k)=-bc(isum2+k)*0.5d0                                   3d14s23
      end do                                                            3d14s23
      if(lprt)then                                                      11d28s22
       write(6,*)('active 2e part: '),sum                                6d20s22
       call prntm2(bc(isum2),1,mden,1)                                  3d14s23
       if(nsymb.eq.1)then                                               3d27s23
        bc(ibcoff)=bc(isum2)
        bc(ibcoff+1)=bc(isum2+15)
        bc(ibcoff+2)=bc(isum2+20)
        bc(ibcoff+3)=bc(isum2+2)
        bc(ibcoff+4)=bc(isum2+5)
        bc(ibcoff+5)=bc(isum2+12)
        bc(ibcoff+6)=bc(isum2+14)
        bc(ibcoff+7)=bc(isum2+9)
        write(6,*)('in C2v order: ')
        call prntm2(bc(ibcoff),1,8,1)
       end if
       write(6,*)('ci contribution: '),cipart
       write(6,*)('over contribtion: '),ovrpart
       write(6,*)('orb contribution: '),orbpart
      end if                                                            7d25s22
      ibcoff=ibcoffo                                                    6d20s22
      do k=1,mden                                                       3d14s23
       km=k-1                                                           3d14s23
       bodc(k)=bodc(k)+bc(isum2+km)                                     3d14s23
      end do                                                            3d14s23
      iaddc=1                                                           3d20s23
      icdc=1                                                            3d20s23
      do is=1,nstate                                                    3d20s23
       do ix=is+1,nstate                                                3d21s23
        if(isinfo(2,is).eq.isinfo(2,ix).and.isinfo(3,is).eq.            4d14s23
     $                 isinfo(3,ix))then                                4d14s23
         if(multh(isinfo(1,ix),isinfo(1,is)).eq.ipuse)then               3d20s23
          icc=icdc                                                      3d20s23
          idd=iaddc                                                     3d20s23
          do irp=1,isinfo(4,ix)                                         3d20s23
           do ir=1,isinfo(4,is)                                         3d20s23
            sum=cdirb(idd)+cdc(icc)                                      3d21s23
            cdirb(idd)=sum                                              3d21s23
            idd=idd+1                                                   3d20s23
            icc=icc+1                                                   3d20s23
           end do                                                       3d20s23
          end do                                                        3d20s23
          icdc=icdc+isinfo(4,is)*isinfo(4,ix)                           3d31s23
         end if                                                         3d20s23
          iaddc=iaddc+isinfo(4,is)*isinfo(4,ix)                         3d20s23
        end if                                                          3d20s23
       end do                                                           3d20s23
      end do                                                            3d20s23
      return
      end
