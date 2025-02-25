c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine pshamdnona1(npsb,ipsb,vecb,nsbetab,npsk,ipsk,veck,     6d23s22
     $     nsbetak,numb,iaorb,nalpha,                                   6d23s22
     $                 iborb,nbeta,itmpa,itmpb,iden1,jdenpt,            8d19s14
     $                 itmpc,itmpd,itmpe,itmpf,norb,                    6d23s22
     $                 nsymb,nroot,nowpro,numpro,wgt,                   6d23s22
     $                 l2e,ipxder,bc,ibc,mdenoff)                       3d16s23
      implicit real*8 (a-h,o-z)
      integer*2 ipsb(4,*),ipsk(4,*)                                     6d23s22
      integer*1 iaorb,iborb,itmpc,itmpd,itmpe,itmpf                     11d14s05
      logical lprint,l2e                                                6d22s22
      include "common.cas"
      include "common.store"                                            8d15s06
c
c     compute contribution to 1 and 2 particle densities for off diagonal
c     hamiltonian matrix elements.
c
      dimension vecb(npsb,nroot),veck(npsk,nroot),iaorb(nalpha,1),      6d23s22
     $          iborb(nbeta,1),itmpa(1),itmpb(1),iden1(1),
     $          jdenpt(1),itmp2(2,2),itmpc(1),itmpd(1),nsbetak(*),      6d23s22
     $          itmpe(1),itmpf(1),nsbetab(*),dbg(11),                   6d23s22
     $          wgt(nroot),ipxder(4,8,8,8),igoalx(4)                    6d28s22
      data icall/0/
      save icall
      lprint=.false.                                                     8d19s14
      icall=icall+1
      ism=ibcoff                                                        8d15s06
      irelo=ism+norb                                                    8d15s06
      ibcoff=irelo+norb                                                 8d15s06
      call enough('pshamdnona1.  1',bc,ibc)
      jsm=ism-1                                                         8d15s06
      jrelo=irelo-1                                                     8d15s06
      do i=1,nsymb                                                      8d15s06
       do j=1,iacto(i)                                                  8d15s06
        ibc(jsm+j)=i                                                    8d15s06
        ibc(jrelo+j)=j                                                  8d15s06
  153   format(3i5)                                                     8d15s06
       end do                                                           8d15s06
       jsm=jsm+iacto(i)                                                 8d15s06
       jrelo=jrelo+iacto(i)                                             8d15s06
      end do                                                            8d15s06
      jsm=ism-1                                                         8d15s06
      jrelo=irelo-1                                                     8d15s06
      do ik=1+nowpro,npsk,numpro                                        6d24s22
       ipa=ipsk(1,ik)                                                   6d23s22
       ipb=ipsk(2,ik)                                                   6d23s22
       isbu=ipsk(3,ik)                                                  6d23s22
       nasum2=0                                                         8d15s06
       do isb=1,isbu-1                                                  8d15s06
        nasum2=nasum2+nadet(isb)                                        8d15s06
       end do                                                           8d15s06
       nbsum2=0                                                         8d15s06
       do isb=1,nsbetak(isbu)-1                                         6d23s22
        nbsum2=nbsum2+nbdet(isb)                                        8d15s06
       end do                                                           8d15s06
       ia=ipa+nasum2                                                    8d15s06
       ib=ipb+nbsum2                                                    8d15s06
       do jbb=1,npsb                                                     6d23s22
        vv=0d0                                                          8d18s14
        do k=1,nroot                                                    8d18s14
         vv=vv+wgt(k)*vecb(jbb,k)*veck(ik,k)                             6d23s22
        end do                                                          8d18s14
        jpa=ipsb(1,jbb)                                                 6d23s22
        jpb=ipsb(2,jbb)                                                 6d23s22
        jsbu=ipsb(3,jbb)                                                6d23s22
        nasum1=0                                                        8d15s06
        do isb=1,jsbu-1                                                 8d15s06
         nasum1=nasum1+nadet(isb)                                       8d15s06
        end do                                                          8d15s06
        nbsum1=0                                                        8d15s06
        do isb=1,nsbetab(jsbu)-1                                         8d15s06
         nbsum1=nbsum1+nbdet(isb)                                       8d15s06
        end do                                                          8d15s06
        ja=jpa+nasum1                                                   8d15s06
        jb=jpb+nbsum1                                                   8d15s06
        isuma=0
        do i1=1,norb                                                    8d15s06
         itmpc(i1)=0
         itmpd(i1)=0
         itmpe(i1)=0
         itmpf(i1)=0
        end do
c     bra
        do i1=1,nalpha
         itmpc(iaorb(i1,ja))=1
c     ket
         itmpd(iaorb(i1,ia))=1
        end do
c     bra-ket
        do i1=1,norb                                                    8d15s06
         itmpa(i1)=itmpc(i1)-itmpd(i1)
         isuma=isuma+iabs(itmpa(i1))
        end do
        isumb=0
        do i1=1,nbeta
         itmpe(iborb(i1,jb))=1
         itmpf(iborb(i1,ib))=1
        end do
c     bra-ket
        do i1=1,norb                                                    8d15s06
         itmpb(i1)=itmpe(i1)-itmpf(i1)
         isumb=isumb+iabs(itmpb(i1))
        end do
        isum=isumb+isuma
        if(lprint)write(6,*)('no. difference = '),isum,isuma,isumb
        if(isum.eq.2)then
         if(lprint)write(6,*)('single '),isuma,isumb
         if(isuma.eq.0)then
          if(lprint)write(6,*)('beta single'),i,j
c     bra-ket
          do i1=1,norb                                                  8d15s06
           if(itmpb(i1).eq.1)nbra=i1
           if(itmpb(i1).eq.-1)nket=i1
          end do
          nb=0
          do i1=min(nbra,nket)+1,max(nbra,nket)-1
           if(itmpe(i1).eq.1)nb=nb+1
          end do
          if(mod(nb,2).eq.0)then
           phs=1d0
          else
           phs=-1d0
          end if
          isbket=ibc(jsm+nket)                                             8d15s06
          isbbra=ibc(jsm+nbra)                                           8d15s06
          inket=ibc(jrelo+nket)-1                                       6d23s22
          inbra=ibc(jrelo+nbra)-1                                       6d23s22
          iad1=iden1(isbbra)+inbra+iacto(isbbra)*inket                  6d23s22
          if(lprint)write(6,*)('from to '),inbra,inket,isbket,phs,vv,
     $         iad1,iad1-igoal
          if(l2e)then                                                   6d22s22
           bc(iad1)=bc(iad1)+phs*vv                                     3d16s23
           ixa=inbra+iacto(isbbra)*inket                                6d23s22
           ixb=inket+iacto(isbket)*inbra-ixa                            6d23s22
           nrow=iacto(isbbra)*iacto(isbket)                             6d23s22
           do i1=1,nalpha                                                8d15s06
            isx=ibc(jsm+iaorb(i1,ia))                                    8d15s06
            iisx=ibc(jrelo+iaorb(i1,ia))-1                              6d23s22
            i2eu=ipxder(1,isbbra,isbket,isx)                            6d23s22
            icol=iisx*(iacto(isx)+1)                                    6d23s22
            irow=ixa+ipxder(2,isbbra,isbket,isx)*ixb                    6d28s22
            ix12=irow+nrow*icol                                         6d28s22
            nrr=iacto(isx)**2                                           6d28s22
            ix21=icol+nrr*irow-ix12                                     6d28s22
            iad1=jdenpt(i2eu)+ix12+ipxder(4,isbbra,isbket,isx)*ix21     6d28s22
            bc(iad1)=bc(iad1)+phs*vv                                    6d23s22
           end do                                                        8d15s06
           if(lprint)write(6,*)('beta '),ib
           do i1=1,nbeta                                                 8d15s06
            isx=ibc(jsm+iborb(i1,ib))                                    8d15s06
            iisx=ibc(jrelo+iborb(i1,ib))-1                              6d23s22
            i2eu=ipxder(1,isbbra,isbket,isx)                            6d23s22
            icol=iisx*(iacto(isx)+1)                                    6d23s22
            irow=ixa+ipxder(2,isbbra,isbket,isx)*ixb                    6d28s22
            ix12=irow+nrow*icol                                         6d28s22
            nrr=iacto(isx)**2                                           6d28s22
            ix21=icol+nrr*irow-ix12                                     6d28s22
            iad1=jdenpt(i2eu)+ix12+ipxder(4,isbbra,isbket,isx)*ix21     6d28s22
            bc(iad1)=bc(iad1)+phs*vv                                    6d23s22
            i2eu=ipxder(1,isbbra,isx,isx)                               6d23s22
            ix1=inbra+iacto(isbbra)*iisx                                6d23s22
            ix2=iisx+iacto(isx)*inbra-ix1                               6d23s22
            irow=ix1+ipxder(2,isbbra,isx,isx)*ix2                       6d23s22
            ix3=iisx+iacto(isx)*inket                                   6d23s22
            ix4=inket+iacto(isbket)*iisx-ix3                            6d23s22
            icol=ix3+ipxder(3,isbbra,isx,isx)*ix4                       6d23s22
            nrr=iacto(isbbra)*iacto(isx)                                6d23s22
            ix12=irow+nrr*icol                                          6d28s22
            mrr=iacto(isx)*iacto(isbket)                                6d28s22
            ix21=icol+mrr*irow-ix12                                     6d28s22
            iad1=jdenpt(i2eu)+ix12+ipxder(4,isbbra,isx,isx)*ix21        6d28s22
            bc(iad1)=bc(iad1)-phs*vv                                    6d23s22
           end do                                                        8d15s06
          else                                                          3d16s23
           nad=iacto(isbbra)*iacto(isbket)                              3d16s23
           iad1=iad1+nad*mdenoff                                        3d16s23
           do k1=1,nroot                                                3d16s23
            do k2=1,k1                                                  3d16s23
             bc(iad1)=bc(iad1)+vecb(jbb,k1)*veck(ik,k2)*phs             3d16s23
             iad1=iad1+nad                                              3d16s23
            end do                                                      3d16s23
           end do                                                       3d16s23
          end if                                                        6d22s22
         else
          if(lprint)write(6,*)('alpha single ')
          do i1=1,norb                                                  8d15s06
           if(itmpa(i1).eq.1)nbra=i1
           if(itmpa(i1).eq.-1)nket=i1
          end do
          nb=0
          do i1=min(nbra,nket)+1,max(nbra,nket)-1
           if(itmpc(i1).eq.1)nb=nb+1
          end do
          if(mod(nb,2).eq.0)then
           phs=1d0
          else
           phs=-1d0
          end if
          isbket=ibc(jsm+nket)                                             8d15s06
          isbbra=ibc(jsm+nbra)                                           8d15s06
          inket=ibc(jrelo+nket)-1                                       6d23s22
          inbra=ibc(jrelo+nbra)-1                                       6d23s22
          iad1=iden1(isbbra)+inbra+iacto(isbbra)*inket                  6d23s22
          if(l2e)then                                                   6d22s22
           bc(iad1)=bc(iad1)+phs*vv                                     3d16s23
           if(lprint)write(6,*)('in doubles section '),inbra,isbbra,
     $         inket,isbket
           ixa=inbra+iacto(isbbra)*inket                                6d23s22
           ixb=inket+iacto(isbket)*inbra-ixa                            6d23s22
           nrow=iacto(isbbra)*iacto(isbket)                             6d23s22
           do i1=1,nbeta                                                 8d15s06
            isx=ibc(jsm+iborb(i1,ib))                                    8d15s06
            iisx=ibc(jrelo+iborb(i1,ib))-1                              6d23s22
            i2eu=ipxder(1,isbbra,isbket,isx)                            6d23s22
            icol=iisx*(iacto(isx)+1)                                    6d23s22
            irow=ixa+ipxder(2,isbbra,isbket,isx)*ixb                    6d28s22
            ix12=irow+nrow*icol                                         6d28s22
            nrr=iacto(isx)**2                                           6d28s22
            ix21=icol+nrr*irow-ix12                                     6d28s22
            iad1=jdenpt(i2eu)+ix12+ipxder(4,isbbra,isbket,isx)*ix21     6d28s22
            bc(iad1)=bc(iad1)+phs*vv                                    6d23s22
           end do                                                        8d15s06
           do i1=1,nalpha                                                8d15s06
            isx=ibc(jsm+iaorb(i1,ia))                                    8d15s06
            iisx=ibc(jrelo+iaorb(i1,ia))-1                              6d23s22
            i2eu=ipxder(1,isbbra,isbket,isx)                            6d23s22
            icol=iisx*(iacto(isx)+1)                                    6d23s22
            irow=ixa+ipxder(2,isbbra,isbket,isx)*ixb                    6d28s22
            ix12=irow+nrow*icol                                         6d28s22
            nrr=iacto(isx)**2                                           6d28s22
            ix21=icol+nrr*irow-ix12                                     6d28s22
            iad1=jdenpt(i2eu)+ix12+ipxder(4,isbbra,isbket,isx)*ix21     6d28s22
            if(lprint)write(6,*)('i1 = '),i1,iisx,isx,iad1-igoal
            bc(iad1)=bc(iad1)+phs*vv                                    6d23s22
            i2eu=ipxder(1,isbbra,isx,isx)                               6d23s22
            ix1=inbra+iacto(isbbra)*iisx                                6d23s22
            ix2=iisx+iacto(isx)*inbra-ix1                               6d23s22
            irow=ix1+ipxder(2,isbbra,isx,isx)*ix2                       6d23s22
            ix3=iisx+iacto(isx)*inket                                   6d23s22
            ix4=inket+iacto(isbket)*iisx-ix3                            6d23s22
            icol=ix3+ipxder(3,isbbra,isx,isx)*ix4                       6d23s22
            nrr=iacto(isx)*iacto(isbbra)                                6d24s22
            ix12=irow+nrr*icol                                          6d28s22
            mrr=iacto(isx)*iacto(isbket)                                6d28s22
            ix21=icol+mrr*irow-ix12                                     6d28s22
            iad1=jdenpt(i2eu)+ix12+ipxder(4,isbbra,isx,isx)*ix21        6d28s22
            bc(iad1)=bc(iad1)-phs*vv                                    6d23s22
           end do                                                        8d15s06
          else                                                          3d16s23
           nad=iacto(isbbra)*iacto(isbket)                              3d16s23
           iad1=iad1+nad*mdenoff                                        3d16s23
           do k1=1,nroot                                                3d16s23
            do k2=1,k1                                                  3d16s23
             bc(iad1)=bc(iad1)+vecb(jbb,k1)*veck(ik,k2)*phs             3d16s23
             iad1=iad1+nad                                              3d16s23
            end do                                                      3d16s23
           end do                                                       3d16s23
          end if                                                        6d22s22
         end if
        else if(isum.eq.4.and.l2e)then                                  6d22s22
         if(lprint)write(6,*)('double '),isuma,isumb
         if(isuma.eq.0)then
          if(lprint)write(6,*)('bb')
          ip=1
          im=1
          do i1=1,norb                                                  8d15s06
c     ket
           if(itmpb(i1).eq.1)then
            itmp2(ip,1)=i1
            ip=ip+1
           end if
c     bra
           if(itmpb(i1).eq.-1)then
            itmp2(im,2)=i1
            im=im+1
           end if
          end do
          call phsdet(itmpe,itmpb,itmp2,iphsb,norb)
          call phsdet(itmpf,itmpb,itmp2(1,2),iphsk,norb)
          iphs=iphsb+iphsk
          if(mod(iphs,2).eq.0)then
           phs=1d0
          else
           phs=-1d0
          end if
          isk1=ibc(jsm+itmp2(1,1))                                      8d15s06
          isb1=ibc(jsm+itmp2(1,2))                                      8d15s06
          isk2=ibc(jsm+itmp2(2,1))                                      8d15s06
          isb2=ibc(jsm+itmp2(2,2))                                      8d15s06
          iisxk1=ibc(jrelo+itmp2(1,1))-1                                6d23s22
          iisxb1=ibc(jrelo+itmp2(1,2))-1                                6d23s22
          iisxk2=ibc(jrelo+itmp2(2,1))-1                                6d23s22
          iisxb2=ibc(jrelo+itmp2(2,2))-1                                6d23s22
          i2eu=ipxder(1,isb1,isk1,isb2)                                 6d23s22
          ix1=iisxb1+iacto(isb1)*iisxk1                                 6d23s22
          ix2=iisxk1+iacto(isk1)*iisxb1-ix1                             6d23s22
          irow=ix1+ipxder(2,isb1,isk1,isb2)*ix2                         6d23s22
          ix3=iisxb2+iacto(isb2)*iisxk2                                 6d23s22
          ix4=iisxk2+iacto(isk2)*iisxb2-ix3                             6d23s22
          icol=ix3+ipxder(3,isb1,isk1,isb2)*ix4                         6d23s22
          nrr=iacto(isb1)*iacto(isk1)                                   6d23s22
          ix12=irow+nrr*icol                                            6d28s22
          mrr=iacto(isb2)*iacto(isk2)                                   6d28s22
          ix21=icol+mrr*irow-ix12                                       6d28s22
          iad1=jdenpt(i2eu)+ix12+ipxder(4,isb1,isk1,isb2)*ix21          6d28s22
          bc(iad1)=bc(iad1)+phs*vv                                      6d23s22
          i2eu=ipxder(1,isb1,isk2,isb2)                                 6d23s22
          ix1=iisxb1+iacto(isb1)*iisxk2                                 6d23s22
          ix2=iisxk2+iacto(isk2)*iisxb1-ix1                             6d23s22
          irow=ix1+ipxder(2,isb1,isk2,isb2)*ix2                         6d23s22
          ix3=iisxb2+iacto(isb2)*iisxk1                                 6d23s22
          ix4=iisxk1+iacto(isk1)*iisxb2-ix3                             6d23s22
          icol=ix3+ipxder(3,isb1,isk2,isb2)*ix4                         6d23s22
          nrr=iacto(isb1)*iacto(isk2)                                   6d23s22
          ix12=irow+nrr*icol                                            6d28s22
          mrr=iacto(isb2)*iacto(isk1)                                   6d28s22
          ix21=icol+mrr*irow-ix12                                       6d28s22
          iad1=jdenpt(i2eu)+ix12+ipxder(4,isb1,isk2,isb2)*ix21          6d28s22
          bc(iad1)=bc(iad1)-phs*vv                                      6d23s22
         else if(isumb.eq.0)then
          if(lprint)write(6,*)('aa')
          ip=1
          im=1
          do i1=1,norb                                                  8d15s06
c     ket
           if(itmpa(i1).eq.1)then
            itmp2(ip,1)=i1
            ip=ip+1
           end if
c     bra
           if(itmpa(i1).eq.-1)then
            itmp2(im,2)=i1
            im=im+1
           end if
          end do
          call phsdet(itmpc,itmpa,itmp2,iphsb,norb)                     8d15s06
          call phsdet(itmpd,itmpa,itmp2(1,2),iphsk,norb)                8d15s06
          iphs=iphsb+iphsk
          if(mod(iphs,2).eq.0)then
           phs=1d0
          else
           phs=-1d0
          end if
          if(lprint)write(6,*)('phases '),iphsb,iphsk,iphs,phs
          isk1=ibc(jsm+itmp2(1,1))                                      8d15s06
          isb1=ibc(jsm+itmp2(1,2))                                      8d15s06
          isk2=ibc(jsm+itmp2(2,1))                                      8d15s06
          isb2=ibc(jsm+itmp2(2,2))                                      8d15s06
          iisk1=ibc(jrelo+itmp2(1,1))-1                                 6d23s22
          iisb1=ibc(jrelo+itmp2(1,2))-1                                 6d23s22
          iisk2=ibc(jrelo+itmp2(2,1))-1                                 6d23s22
          iisb2=ibc(jrelo+itmp2(2,2))-1                                 6d23s22
          i2eu=ipxder(1,isb1,isk1,isb2)                                 6d23s22
          ix1=iisb1+iacto(isb1)*iisk1                                   6d23s22
          ix2=iisk1+iacto(isk1)*iisb1-ix1                               6d23s22
          irow=ix1+ipxder(2,isb1,isk1,isb2)*ix2                         6d23s22
          ix3=iisb2+iacto(isb2)*iisk2                                   6d23s22
          ix4=iisk2+iacto(isk2)*iisb2-ix3                               6d23s22
          icol=ix3+ipxder(3,isb1,isk1,isb2)*ix4                         6d23s22
          nrr=iacto(isb1)*iacto(isk1)                                   6d23s22
          ix12=irow+nrr*icol                                            6d28s22
          mrr=iacto(isb2)*iacto(isk2)                                   6d28s22
          ix21=icol+mrr*irow-ix12                                       6d28s22
          iad1=jdenpt(i2eu)+ix12+ipxder(4,isb1,isk1,isb2)*ix21          6d28s22
          bc(iad1)=bc(iad1)+phs*vv                                      6d23s22
          i2eu=ipxder(1,isb1,isk2,isb2)                                 6d23s22
          ix1=iisb1+iacto(isb1)*iisk2                                   6d23s22
          ix2=iisk2+iacto(isk2)*iisb1-ix1                               6d23s22
          irow=ix1+ipxder(2,isb1,isk2,isb2)*ix2                         6d23s22
          ix3=iisb2+iacto(isb2)*iisk1                                   6d23s22
          ix4=iisk1+iacto(isk1)*iisb2-ix3                               6d23s22
          icol=ix3+ipxder(3,isb1,isk2,isb2)*ix4                         6d23s22
          nrr=iacto(isb1)*iacto(isk2)                                   6d23s22
          ix12=irow+nrr*icol                                            6d28s22
          mrr=iacto(isb2)*iacto(isk1)                                   6d28s22
          ix21=icol+mrr*irow-ix12                                       6d28s22
          iad2=jdenpt(i2eu)+ix12+ipxder(4,isb1,isk2,isb2)*ix21          6d28s22
          bc(iad2)=bc(iad2)-phs*vv                                      6d23s22
         else if(isuma.eq.2.and.isumb.eq.2)then
          do i1=1,norb                                                  8d15s06
           if(itmpa(i1).eq.1)nketa=i1
           if(itmpa(i1).eq.-1)nbraa=i1
          end do
          nb=0
          do i1=min(nbraa,nketa)+1,max(nbraa,nketa)-1
           if(itmpc(i1).eq.1)nb=nb+1
          end do
          do i1=1,norb                                                  8d15s06
           if(itmpb(i1).eq.1)nketb=i1
           if(itmpb(i1).eq.-1)nbrab=i1
          end do
          do i1=min(nbrab,nketb)+1,max(nbrab,nketb)-1
           if(itmpe(i1).eq.1)nb=nb+1
          end do
          if(mod(nb,2).eq.0)then
           phs=1d0
          else
           phs=-1d0
          end if
          isba=ibc(jsm+nbraa)                                           6d23s22
          iska=ibc(jsm+nketa)                                           6d23s22
          isbb=ibc(jsm+nbrab)                                           6d23s22
          iskb=ibc(jsm+nketb)                                           6d23s22
          inbraa=ibc(jrelo+nbraa)-1                                     6d23s22
          inketa=ibc(jrelo+nketa)-1                                     6d23s22
          inbrab=ibc(jrelo+nbrab)-1                                     6d23s22
          inketb=ibc(jrelo+nketb)-1                                     6d23s22
          i2eu=ipxder(1,isba,iska,isbb)                                 6d23s22
          ix1=inbraa+iacto(isba)*inketa                                 6d23s22
          ix2=inketa+iacto(iska)*inbraa-ix1                             6d23s22
          irow=ix1+ipxder(2,isba,iska,isbb)*ix2                         6d23s22
          ix3=inbrab+iacto(isbb)*inketb                                 6d23s22
          ix4=inketb+iacto(iskb)*inbrab-ix3                             6d23s22
          icol=ix3+ipxder(3,isba,iska,isbb)*ix4                         6d23s22
          nrr=iacto(isba)*iacto(iska)                                   6d23s22
          ix12=irow+nrr*icol                                            6d28s22
          mrr=iacto(isbb)*iacto(iskb)                                   6d28s22
          ix21=icol+mrr*irow-ix12                                       6d28s22
          iad1=jdenpt(i2eu)+ix12+ipxder(4,isba,iska,isbb)*ix21          6d28s22
          if(lprint)write(6,*)('ab'),nbraa,inbraa,isba,
     $         nketa,inketa,iska,nbrab,inbrab,isbb,nketb,inketb,iskb,
     $         i2eu,iad1-jdenpt(i2eu),iad1-igoal,vv,jdenpt(i2eu),igoal,
     $         ipxder(2,isba,iska,isbb),ipxder(3,isba,iska,isbb)
            if(lprint)write(6,*)iad1-igoal,phs,vv,irow,nrr,icol,
     $           inbrab,iacto(isbb),inketb,iacto(iskb),
     $           ipxder(3,isba,iska,isbb),ix3,ix4
          bc(iad1)=bc(iad1)+phs*vv                                      6d23s22
         end if
        end if
       end do
      end do
      ibcoff=ism                                                        2d16s07
      return
      end
