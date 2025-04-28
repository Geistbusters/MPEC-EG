c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine pshamdbl(nps,ips,vec,numb,iaorb,nalpha,                7d7s22
     $                 iborb,nbeta,itmpa,itmpb,jdenpt,                  7d7s22
     $                 itmpc,itmpd,itmpe,itmpf,debug1,debug2,norb,      8d15s06
     $                 ncsym,nsymb,nsbeta,nroot,dbg,nowpro,numpro,wgt,  11d9s22
     $                 bc,ibc,mdenoff,iden1,igoal)                            3d15s23
      implicit real*8 (a-h,o-z)
      integer*2 ips(4,1)                                                8d22s06
      integer*1 iaorb,iborb,itmpc,itmpd,itmpe,itmpf                     11d14s05
      logical lprint                                                    7d7s22
      include "common.cas"
      include "common.store"                                            8d15s06
c
c     compute full 2e density for bi-linear derivatives.                7d7s22
c     and root separate 1e densities.                                   3d15s23
c
      dimension vec(nps,nroot),iaorb(nalpha,1),                         7d7s22
     $          iborb(nbeta,1),itmpa(1),itmpb(1),                       3d15s23
     $          jdenpt(1),itmp2(2,2),itmpc(1),itmpd(1),igoalx(10),      7d8s22
     $          itmpe(1),itmpf(1),ncsym(8),nsbeta(8),dbg(11),wgt(nroot),3d15s23
     $          iden1(*)                                                3d15s23
      dimension bugf(6)
      data icall/0/
      save icall
      do i=1,6
       bugf(i)=1d0
      end do
      bugf(1)=1d0
      ngx=1
      igoalx(1)=1844364+9+3
      lprint=.false.                                                    8d19s14
      icall=icall+1
      ism=ibcoff                                                        8d15s06
      irelo=ism+norb                                                    8d15s06
      ibcoff=irelo+norb                                                 8d15s06
      ivv=ibcoff                                                        3d14s23
      ntri=(nroot*(nroot+1))/2                                          3d14s23
      ivvp=ivv+ntri                                                     3d14s23
      ibcoff=ivvp+ntri                                                  3d14s23
      call enough('pshamdbl.  1',bc,ibc)
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
      do i=1+nowpro,nps,numpro                                          7d7s22
       ipa=ips(1,i)                                                     8d22s06
       ipb=ips(2,i)                                                     8d22s06
       isbu=ips(3,i)                                                    8d22s06
       nasum2=0                                                         8d15s06
       do isb=1,isbu-1                                                  8d15s06
        nasum2=nasum2+nadet(isb)                                        8d15s06
       end do                                                           8d15s06
       nbsum2=0                                                         8d15s06
       do isb=1,nsbeta(isbu)-1                                          8d15s06
        nbsum2=nbsum2+nbdet(isb)                                        8d15s06
       end do                                                           8d15s06
       ia=ipa+nasum2                                                    8d15s06
       ib=ipb+nbsum2                                                    8d15s06
c
c     for diagonals, we have sum a<b -2*abba
c
       do k=ivv,ivv+ntri-1                                              3d14s23
        bc(k)=0d0                                                       3d14s23
       end do                                                           3d14s23
       jvv=ivv                                                          3d14s23
       do k1=1,nroot                                                    3d14s23
        do k2=1,k1                                                      3d14s23
         bc(jvv)=vec(i,k1)*vec(i,k2)                                    3d14s23
         jvv=jvv+1                                                      3d14s23
        end do                                                          3d14s23
       end do                                                           3d14s23
       do k=ivv,ivv+ntri-1                                              3d14s23
        bc(k)=-2d0*bc(k)                                                3d14s23
       end do                                                           3d14s23
       do i1=1,nalpha                                                   7d7s22
        is14=ibc(jsm+iaorb(i1,ia))
        if14=ibc(jrelo+iaorb(i1,ia))-1                                  7d7s22
        jvv=0                                                           4d9s23
        do k1=1,nroot                                                   4d9s23
         do k2=1,k1                                                     4d9s23
          kp=jvv+mdenoff                                                4d9s23
          iad=iden1(is14)+if14+iacto(is14)*(if14+iacto(is14)*kp)         3d15s23
          bc(iad)=bc(iad)+vec(i,k1)*vec(i,k2)                           4d9s23
          jvv=jvv+1                                                     4d9s23
         end do                                                         4d9s23
        end do                                                          4d9s23
        do i2=i1+1,nalpha                                               7d7s22
         is23=ibc(jsm+iaorb(i2,ia))                                     7d7s22
         if23=ibc(jrelo+iaorb(i2,ia))-1                                 7d7s22
         i2eu=ifind2(is14,is23,is23,is14,icase)
         if(icase.eq.1)then                                             7d7s22
          do k=0,ntri-1                                                 3d14s23
           kp=k+mdenoff                                                 3d14s23
           iad=jdenpt(i2eu)+if14+iacto(is14)*(if23+iacto(is23)*(if23     7d7s22
     $        +iacto(is23)*(if14+iacto(is14)*kp)))                      3d14s23
           bc(iad)=bc(iad)+bc(ivv+k)                                    3d14s23
          end do                                                        3d14s23
         else if(icase.eq.4)then                                        7d7s22
          do k=0,ntri-1                                                 3d14s23
           kp=k+mdenoff                                                 3d14s23
           iad=jdenpt(i2eu)+if23+iacto(is23)*(if14+iacto(is14)*(if23     7d7s22
     $        +iacto(is23)*(if14+iacto(is14)*kp)))                      3d14s23
           bc(iad)=bc(iad)-bc(ivv+k)                                    3d14s23
          end do                                                        3d14s23
         else if(icase.eq.3)then                                        7d7s22
          do k=0,ntri-1                                                 3d14s23
           kp=k+mdenoff                                                 3d14s23
           iad=jdenpt(i2eu)+if14+iacto(is14)*(if23+iacto(is23)*(if23     7d7s22
     $        +iacto(is23)*(if14+iacto(is14)*kp)))                      3d14s23
           bc(iad)=bc(iad)-bc(ivv+k)                                    3d14s23
          end do                                                        3d14s23
         else if(icase.eq.2)then                                        7d7s22
          do k=0,ntri-1                                                 3d14s23
           kp=k+mdenoff                                                 3d14s23
           iad=jdenpt(i2eu)+if23+iacto(is23)*(if14+iacto(is14)*(if14     7d7s22
     $        +iacto(is14)*(if23+iacto(is23)*kp)))                      3d14s23
           bc(iad)=bc(iad)+bc(ivv+k)                                    3d14s23
          end do                                                        3d14s23
         end if                                                         7d7s22
        end do                                                          7d7s22
       end do                                                           7d7s22
       do i1=1,nbeta                                                    7d7s22
        is14=ibc(jsm+iborb(i1,ib))
        if14=ibc(jrelo+iborb(i1,ib))-1                                  7d7s22
        do k=0,ntri-1                                                   3d15s23
         kp=k+mdenoff                                                   3d15s23
         iad=iden1(is14)+if14+iacto(is14)*(if14+iacto(is14)*kp)         3d15s23
         bc(iad)=bc(iad)-0.5d0*bc(ivv+k)                                3d15s23
        end do                                                          3d15s23
        do i2=i1+1,nbeta                                                7d7s22
         is23=ibc(jsm+iborb(i2,ib))                                     7d7s22
         if23=ibc(jrelo+iborb(i2,ib))-1                                 7d7s22
         i2eu=ifind2(is14,is23,is23,is14,icase)
         if(icase.eq.1)then                                             7d7s22
          do k=0,ntri-1                                                 3d14s23
           kp=k+mdenoff                                                 3d14s23
           iad=jdenpt(i2eu)+if14+iacto(is14)*(if23+iacto(is23)*(if23     7d7s22
     $        +iacto(is23)*(if14+iacto(is14)*kp)))                      3d14s23
           bc(iad)=bc(iad)+bc(ivv+k)                                    3d14s23
          end do                                                        3d14s23
         else if(icase.eq.4)then                                        7d7s22
          do k=0,ntri-1                                                 3d14s23
           kp=k+mdenoff                                                 3d14s23
           iad=jdenpt(i2eu)+if23+iacto(is23)*(if14+iacto(is14)*(if23     7d7s22
     $        +iacto(is23)*(if14+iacto(is14)*kp)))                       3d14s23
           bc(iad)=bc(iad)-bc(ivv+k)                                    3d14s23
          end do                                                        3d14s23
         else if(icase.eq.3)then                                        7d7s22
          do k=0,ntri-1                                                 3d14s23
           kp=k+mdenoff                                                 3d14s23
           iad=jdenpt(i2eu)+if14+iacto(is14)*(if23+iacto(is23)*(if23     7d7s22
     $        +iacto(is23)*(if14+iacto(is14)*kp)))                       3d14s23
           bc(iad)=bc(iad)-bc(ivv+k)                                    3d14s23
          end do                                                        3d14s23
         else if(icase.eq.2)then                                        7d7s22
          do k=0,ntri-1                                                 3d14s23
           kp=k+mdenoff                                                 3d14s23
           iad=jdenpt(i2eu)+if23+iacto(is23)*(if14+iacto(is14)*(if14     7d7s22
     $        +iacto(is14)*(if23+iacto(is23)*kp)))                       3d14s23
           bc(iad)=bc(iad)+bc(ivv+k)                                    3d14s23
          end do                                                        3d14s23
         end if                                                         7d7s22
        end do                                                          7d7s22
       end do                                                           7d7s22
       do j=1,nps                                                       3d27s23
        if(j.ne.i)then                                                  3d27s23
         do k=ivv,ivv+ntri-1                                             3d14s23
          bc(k)=0d0                                                      3d14s23
         end do                                                          3d14s23
         jvv=ivv                                                         3d14s23
         do k1=1,nroot                                                   3d14s23
          do k2=1,k1                                                    3d27s23
           bc(jvv)=vec(j,k1)*vec(i,k2)                                  3d27s23
           jvv=jvv+1                                                     3d14s23
          end do                                                         3d14s23
         end do                                                          3d14s23
c     j is bra, i is ket
         jpa=ips(1,j)                                                    8d22s06
         jpb=ips(2,j)                                                    8d22s06
         jsbu=ips(3,j)                                                   8d22s06
         nasum1=0                                                        8d15s06
         do isb=1,jsbu-1                                                 8d15s06
          nasum1=nasum1+nadet(isb)                                       8d15s06
         end do                                                          8d15s06
         nbsum1=0                                                        8d15s06
         do isb=1,nsbeta(jsbu)-1                                         8d15s06
          nbsum1=nbsum1+nbdet(isb)                                       8d15s06
         end do                                                          8d15s06
         ja=jpa+nasum1                                                   8d15s06
         jb=jpb+nbsum1                                                   8d15s06
         if(lprint)write(6,*)('ia,ib,ja,jb '),ia,ib,ja,jb
         isuma=0
         do i1=1,norb                                                    8d15s06
          itmpc(i1)=0
          itmpd(i1)=0
          itmpe(i1)=0
          itmpf(i1)=0
         end do
         do i1=1,nalpha
          itmpc(iaorb(i1,ja))=1
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
         if(lprint)write(6,*)('no. difference = '),isum
         if(isum.eq.2)then
          if(lprint)write(6,*)('single ')
c
c     -2 sum c cbkc
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
           inket=ibc(jrelo+nket)-1                                        7d8s22
           inbra=ibc(jrelo+nbra)-1                                        7d8s22
           do k=0,ntri-1                                                 3d14s23
            kp=k+mdenoff                                                 3d15s23
            iad1=iden1(isbbra)+inbra+iacto(isbbra)                       3d15s23
     $          *(inket+iacto(isbbra)*kp)                               3d15s23
            bc(iad1)=bc(iad1)+bc(ivv+k)*phs
            bc(ivvp+k)=-bc(ivv+k)*2d0*phs                                3d14s23
           end do                                                        3d14s23
           if(isbket.ne.isbbra)then                                          8d15s06
            write(6,*)('single of wrong symmetry '),isbket,isbbra,nket,
     $         nbra
            stop
           end if                                                        8d15s06
           if(lprint)write(6,*)('bra to '),inbra,inket,isbket
           do i1=1,nbeta                                                 7d7s22
            if(iborb(i1,ib).ne.nket)then                                 7d7s22
             is14=ibc(jsm+iborb(i1,ib))                                  7d7s22
             if14=ibc(jrelo+iborb(i1,ib))-1                              7d8s22
             i2eu=ifind2(is14,isbbra,isbket,is14,icase)                  7d7s22
             if(icase.eq.1)then                                          7d7s22
              do k=0,ntri-1                                              3d14s23
               kp=k+mdenoff                                                 3d14s23
               iad=jdenpt(i2eu)+if14+iacto(is14)*(inbra+iacto(isbbra)     7d7s22
     $           *(inket+iacto(isbket)*(if14+iacto(is14)*kp)))           3d14s23
               bc(iad)=bc(iad)+bc(ivvp+k)                                3d14s23
              end do                                                     3d14s23
             else if(icase.eq.4)then                                     7d7s22
              do k=0,ntri-1                                              3d14s23
               kp=k+mdenoff                                                 3d14s23
               iad=jdenpt(i2eu)+inbra+iacto(isbbra)*(if14+iacto(is14)     7d7s22
     $           *(inket+iacto(isbket)*(if14+iacto(is14)*kp)))           3d14s23
               bc(iad)=bc(iad)-bc(ivvp+k)                                3d14s23
              end do                                                     3d14s23
             else if(icase.eq.3)then                                     7d7s22
              do k=0,ntri-1                                              3d14s23
               kp=k+mdenoff                                                 3d14s23
               iad=jdenpt(i2eu)+if14+iacto(is14)*(inbra+iacto(isbbra)     7d7s22
     $           *(if14+iacto(is14)*(inket+iacto(isbket)*kp)))           3d14s23
               bc(iad)=bc(iad)-bc(ivvp+k)                                3d14s23
              end do                                                     3d14s23
             else                                                        7d7s22
              do k=0,ntri-1                                              3d14s23
               kp=k+mdenoff                                                 3d14s23
               iad=jdenpt(i2eu)+inbra+iacto(isbbra)*(if14+iacto(is14)     7d7s22
     $           *(if14+iacto(is14)*(inket+iacto(isbket)*kp)))           3d14s23
               bc(iad)=bc(iad)+bc(ivvp+k)                                3d14s23
              end do                                                     3d14s23
             end if                                                      7d7s22
            end if                                                       7d8s22
           end do                                                        8d15s06
          else
c     bra-ket
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
           isbbra=ibc(jsm+nbra)                                            8d15s06
           if(isbket.ne.isbbra)then                                          8d15s06
            write(6,*)('symmetry of single does not match '),isbket,
     $         isbbra,nket,nbra
            stop
           end if                                                        8d15s06
           inket=ibc(jrelo+nket)-1                                       7d7s22
           inbra=ibc(jrelo+nbra)-1                                       7d7s22
           do k=0,ntri-1                                                 3d14s23
            kp=k+mdenoff                                                 3d15s23
            iad1=iden1(isbbra)+inbra+iacto(isbbra)                       3d15s23
     $           *(inket+iacto(isbbra)*kp)                               3d15s23
            bc(iad1)=bc(iad1)+phs*bc(ivv+k)
            bc(ivvp+k)=-2d0*bc(ivv+k)*phs                                3d14s23
           end do                                                        3d14s23
           do i1=1,nalpha                                                8d15s06
            if(iaorb(i1,ia).ne.nket)then                                 7d8s22
             is14=ibc(jsm+iaorb(i1,ia))                                    8d15s06
             if14=ibc(jrelo+iaorb(i1,ia))-1                              7d8s22
             i2eu=ifind2(is14,isbbra,isbket,is14,icase)                   7d7s22
             if(icase.eq.1)then                                           7d7s22
              do k=0,ntri-1                                              3d14s23
               kp=k+mdenoff                                                 3d14s23
               iad=jdenpt(i2eu)+if14+iacto(is14)*(inbra+iacto(isbbra)      7d7s22
     $           *(inket+iacto(isbket)*(if14+iacto(is14)*kp)))           3d14s23
               bc(iad)=bc(iad)+bc(ivvp+k)                                3d14s23
              end do                                                     3d14s23
             else if(icase.eq.4)then                                      7d7s22
              do k=0,ntri-1                                              3d14s23
               kp=k+mdenoff                                                 3d14s23
               iad=jdenpt(i2eu)+inbra+iacto(isbbra)*(if14+iacto(is14)      7d7s22
     $            *(inket+iacto(isbket)*(if14+iacto(is14)*kp)))          3d14s23
               bc(iad)=bc(iad)-bc(ivvp+k)                                3d14s23
              end do                                                     3d14s23
             else if(icase.eq.3)then                                      7d7s22
              do k=0,ntri-1                                              3d14s23
               kp=k+mdenoff                                                 3d14s23
               iad=jdenpt(i2eu)+if14+iacto(is14)*(inbra+iacto(isbbra)      7d7s22
     $           *(if14+iacto(is14)*(inket+iacto(isbket)*kp)))           3d14s23
               bc(iad)=bc(iad)-bc(ivvp+k)                                3d14s23
              end do                                                     3d14s23
             else                                                         7d7s22
              do k=0,ntri-1                                              3d14s23
               kp=k+mdenoff                                                 3d14s23
               iad=jdenpt(i2eu)+inbra+iacto(isbbra)*(if14+iacto(is14)      7d7s22
     $           *(if14+iacto(is14)*(inket+iacto(isbket)*kp)))           3d14s23
               bc(iad)=bc(iad)+bc(ivvp+k)                                3d14s23
              end do                                                     3d14s23
             end if                                                       7d7s22
            end if                                                       7d8s22
           end do                                                        8d15s06
          end if
         else if(isum.eq.4)then                                          7d7s22
          if(lprint)write(6,*)('double ')
c
c     2 (b1k1b2k2-b1k2b2k1)
          if(isuma.eq.0)then
           if(lprint)write(6,*)('bb')
           ip=1
           im=1
c     bra-ket
           do i1=1,norb                                                  8d15s06
            if(itmpb(i1).eq.1)then
c     bra
             itmp2(ip,1)=i1
             ip=ip+1
            end if
c     ket
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
          else if(isumb.eq.0)then                                        7d8s22
           if(lprint)write(6,*)('aa')
           ip=1
           im=1
           do i1=1,norb                                                  8d15s06
c     bra
            if(itmpa(i1).eq.1)then
             itmp2(ip,1)=i1
             ip=ip+1
            end if
c     ket
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
          else                                                           7d8s22
c     bra-ket
           do i1=1,norb                                                  8d15s06
            if(itmpa(i1).eq.1)nbra1=i1
            if(itmpa(i1).eq.-1)nket1=i1
           end do
           nb=0
           do i1=min(nbra1,nket1)+1,max(nbra1,nket1)-1
            if(itmpc(i1).eq.1)nb=nb+1
           end do
           nbb=nb
           do i1=1,norb                                                  8d15s06
            if(itmpb(i1).eq.1)nbra2=i1
            if(itmpb(i1).eq.-1)nket2=i1
           end do
           do i1=min(nbra2,nket2)+1,max(nbra2,nket2)-1
            if(itmpe(i1).eq.1)nb=nb+1
           end do
           if(mod(nb,2).eq.0)then
            phs=1d0
           else
            phs=-1d0
           end if
           itmp2(1,1)=nbra1                                              7d8s22
           itmp2(1,2)=nket1                                              7d8s22
           itmp2(2,1)=nbra2                                              7d8s22
           itmp2(2,2)=nket2                                              7d8s22
          end if
          do k=0,ntri-1                                                  3d14s23
           bc(ivvp+k)=2d0*bc(ivv+k)*phs                                  3d14s23
          end do                                                         3d14s23
          isb1=ibc(jsm+itmp2(1,1))                                       7d8s22
          isk1=ibc(jsm+itmp2(1,2))                                       7d8s22
          isb2=ibc(jsm+itmp2(2,1))                                       7d8s22
          isk2=ibc(jsm+itmp2(2,2))                                       7d8s22
          iisb1=ibc(jrelo+itmp2(1,1))-1                                  7d8s22
          iisk1=ibc(jrelo+itmp2(1,2))-1                                  7d8s22
          iisb2=ibc(jrelo+itmp2(2,1))-1                                  7d8s22
          iisk2=ibc(jrelo+itmp2(2,2))-1                                  7d8s22
          i2eu=ifind2(isb1,isk1,isb2,isk2,icase)                         8d15s06
          if(icase.eq.1)then                                             7d8s22
           do k=0,ntri-1                                                 3d14s23
            kp=k+mdenoff                                                 3d14s23
            iad=jdenpt(i2eu)+iisb1+iacto(isb1)*(iisk1+iacto(isk1)         7d8s22
     $         *(iisb2+iacto(isb2)*(iisk2+iacto(isk2)*kp)))              3d14s23
            bc(iad)=bc(iad)+bc(ivvp+k)                                   3d14s23
           end do                                                        3d14s23
          else if(icase.eq.4)then                                        7d8s22
           do k=0,ntri-1                                                 3d14s23
            kp=k+mdenoff                                                 3d14s23
            iad=jdenpt(i2eu)+iisk1+iacto(isk1)*(iisb1+iacto(isb1)         7d8s22
     $         *(iisb2+iacto(isb2)*(iisk2+iacto(isk2)*kp)))              3d14s23
            bc(iad)=bc(iad)-bc(ivvp+k)                                   3d14s23
           end do                                                        3d14s23
          else if(icase.eq.3)then                                        7d8s22
           do k=0,ntri-1                                                 3d14s23
            kp=k+mdenoff                                                 3d14s23
            iad=jdenpt(i2eu)+iisb1+iacto(isb1)*(iisk1+iacto(isk1)         7d8s22
     $         *(iisk2+iacto(isk2)*(iisb2+iacto(isb2)*kp)))              3d14s23
            bc(iad)=bc(iad)-bc(ivvp+k)                                   3d14s23
           end do                                                        3d14s23
          else                                                           7d8s22
           do k=0,ntri-1                                                 3d14s23
            kp=k+mdenoff                                                 3d14s23
            iad=jdenpt(i2eu)+iisk1+iacto(isk1)*(iisb1+iacto(isb1)         7d8s22
     $         *(iisk2+iacto(isk2)*(iisb2+iacto(isb2)*kp)))              3d14s23
            bc(iad)=bc(iad)+bc(ivvp+k)                                   3d14s23
           end do                                                        3d14s23
          end if                                                         7d8s22
          if(min(isuma,isumb).eq.0)then
           i2eu=ifind2(isb1,isk2,isb2,isk1,icase)                        8d15s06
           if(icase.eq.1)then                                            7d8s22
            do k=0,ntri-1                                                3d14s23
             kp=k+mdenoff                                                 3d14s23
             iad=jdenpt(i2eu)+iisb1+iacto(isb1)*(iisk2+iacto(isk2)        7d8s22
     $         *(iisb2+iacto(isb2)*(iisk1+iacto(isk1)*kp)))              3d14s23
             bc(iad)=bc(iad)-bc(ivvp+k)                                  3d14s23
            end do                                                       3d14s23
           else if(icase.eq.4)then                                       7d11s22
            do k=0,ntri-1                                                3d14s23
             kp=k+mdenoff                                                 3d14s23
             iad=jdenpt(i2eu)+iisk2+iacto(isk2)*(iisb1+iacto(isb1)        7d8s22
     $         *(iisb2+iacto(isb2)*(iisk1+iacto(isk1)*kp)))              3d14s23
             bc(iad)=bc(iad)+bc(ivvp+k)                                  3d14s23
            end do                                                       3d14s23
           else if(icase.eq.3)then                                       7d11s22
            do k=0,ntri-1                                                3d14s23
             kp=k+mdenoff                                                 3d14s23
             iad=jdenpt(i2eu)+iisb1+iacto(isb1)*(iisk2+iacto(isk2)        7d8s22
     $         *(iisk1+iacto(isk1)*(iisb2+iacto(isb2)*kp)))              3d14s23
             bc(iad)=bc(iad)+bc(ivvp+k)                                  3d14s23
            end do                                                       3d14s23
           else                                                          7d8s22
            do k=0,ntri-1                                                3d14s23
             kp=k+mdenoff                                                 3d14s23
             iad=jdenpt(i2eu)+iisk2+iacto(isk2)*(iisb1+iacto(isb1)        7d8s22
     $         *(iisk1+iacto(isk1)*(iisb2+iacto(isb2)*kp)))              3d14s23
             bc(iad)=bc(iad)-bc(ivvp+k)                                  3d14s23
            end do                                                       3d14s23
           end if                                                        7d8s22
          end if                                                         7d8s22
         end if
        end if                                                          3d27s23
       end do
      end do
      ibcoff=ism                                                        2d16s07
      return
      end
