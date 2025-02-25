c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine pshamd(nps,ips,hdig,vec1,vec2,numb,iaorb,nalpha,       5d16s22
     $                 iborb,nbeta,itmpa,itmpb,iden1,jdenpt,            8d19s14
     $                 itmpc,itmpd,itmpe,itmpf,debug1,debug2,norb,      8d15s06
     $                 ncsym,nsymb,nsbeta,nroot,dbg,nowpro,numpro,wgt,  6d22s22
     $                 l2e,icode,bc,ibc,mdenoff,moreden)                3d27s23
      implicit real*8 (a-h,o-z)
      integer*2 ips(4,1)                                                8d22s06
      integer*1 iaorb,iborb,itmpc,itmpd,itmpe,itmpf                     11d14s05
      logical lprint,l2e                                                6d22s22
      include "common.cas"
      include "common.store"                                            8d15s06
c
c     compute contribution to 1 and 2 particle densities for off diagonal
c     hamiltonian matrix elements.
c
      dimension hdig(1),vec1(nps,nroot),vec2(nps,nroot),iaorb(nalpha,1),5d16s22
     $          iborb(nbeta,1),itmpa(1),itmpb(1),iden1(1),                8d15s06
     $          jdenpt(1),itmp2(2,2),itmpc(1),itmpd(1),                 8d18s14
     $          itmpe(1),itmpf(1),ncsym(8),nsbeta(8),dbg(11),wgt(nroot) 12d30s15
      data icall/0/
      save icall
      lprint=.false.                                                    8d19s14
      icall=icall+1
      ism=ibcoff                                                        8d15s06
      irelo=ism+norb                                                    8d15s06
      ibcoff=irelo+norb                                                 8d15s06
      call enough('pshamd.  1',bc,ibc)
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
      if(icode.eq.2)then                                                6d24s22
       itop=nps-1                                                       6d24s22
      else                                                              6d24s22
       itop=nps                                                         6d24s22
      end if                                                            6d24s22
      ivv=ibcoff                                                        3d27s23
      if(moreden.ne.0)then                                              3d27s23
       ntri=(nroot*(nroot+1))/2                                         3d27s23
       ibcoff=ivv+ntri                                                  3d27s23
       call enough('pshamd.vv',bc,ibc)                                  3d27s23
      end if                                                            3d27s23
      do i=1+nowpro,itop,numpro                                         6d24s22
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
       if(icode.eq.1)then                                               6d23s22
        jbot=1                                                          6d23s22
       else                                                             6d23s22
        jbot=i+1                                                        6d23s22
       end if                                                           6d23s22
       do j=jbot,nps                                                    6d24s22
        if(lprint)write(6,*)('j,i: '),j,i
        vv=0d0                                                          8d18s14
        do k=1,nroot                                                    8d18s14
         vv=vv+wgt(k)*vec1(j,k)*vec2(i,k)                               5d16s22
        end do                                                          8d18s14
        if(icode.eq.2)vv=vv*2d0                                         6d23s22
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
         if(isuma.eq.0)then
          if(lprint)write(6,*)('beta single'),i,j
c     bra-ket
          do i1=1,norb                                                  8d15s06
           if(itmpb(i1).eq.1)nfrom=i1
           if(itmpb(i1).eq.-1)nto=i1
          end do
          nb=0
          do i1=min(nfrom,nto)+1,max(nfrom,nto)-1
           if(itmpe(i1).eq.1)nb=nb+1
          end do
          if(mod(nb,2).eq.0)then
           phs=1d0
          else
           phs=-1d0
          end if
          isbx=ibc(jsm+nto)                                             8d15s06
          isby=ibc(jsm+nfrom)                                           8d15s06
          if(isbx.ne.isby)then                                          8d15s06
           write(6,*)('single of wrong symmetry '),isbx,isby,nto,nfrom  8d15s06
           stop
          end if                                                        8d15s06
          into=ibc(jrelo+nto)                                           8d15s06
          infrom=ibc(jrelo+nfrom)                                       8d15s06
          if(lprint)write(6,*)('from to '),infrom,into,isbx
          iad1=iden1(isby)+infrom-1+iacto(isby)*(into-1)                6d26s22
          if(l2e)then                                                   6d22s22
           bc(iad1)=bc(iad1)+phs*vv*debug1                               12d19s15
           ix=max(into,infrom)                                           8d15s06
           in=min(into,infrom)                                           8d15s06
           ii1=((ix*(ix-1))/2)+in
           na=(iacto(isbx)*(iacto(isbx)+1))/2                            8d15s06
           do i1=1,nalpha                                                8d15s06
            isx=ibc(jsm+iaorb(i1,ia))                                    8d15s06
            i2eu=ifind2(isbx,isbx,isx,isx,icase)                         8d15s06
            iisx=ibc(jrelo+iaorb(i1,ia))                                 8d15s06
            ii2=(iisx*(iisx+1))/2                                        8d15s06
            iad1=jdenpt(i2eu)+ii1-1+na*(ii2-1)                              8d15s06
            bc(iad1)=bc(iad1)+phs*vv*dbg(1)                              8d19s14
           end do                                                        8d15s06
           if(lprint)write(6,*)('beta '),ib
           do i1=1,nbeta                                                 8d15s06
            isx=ibc(jsm+iborb(i1,ib))                                    8d15s06
            i2eu=ifind2(isbx,isbx,isx,isx,icase)                         8d15s06
            iisx=ibc(jrelo+iborb(i1,ib))                                 8d15s06
            ii2=(iisx*(iisx+1))/2                                        8d15s06
            iad1=jdenpt(i2eu)+ii1-1+na*(ii2-1)                           8d18s14
            bc(iad1)=bc(iad1)+phs*vv*dbg(2)                              8d19s14
            i2eu=ifind2(isbx,isx,isbx,isx,icase)                         8d15s06
            if(isbx.eq.isx)then                                          8d15s06
             in=min(iisx,into)                                           8d15s06
             ix=max(iisx,into)                                           8d15s06
             ii3=((ix*(ix-1))/2)+in-1                                    8d15s06
             in=min(iisx,infrom)                                         8d15s06
             ix=max(iisx,infrom)                                         8d15s06
             ii2=((ix*(ix-1))/2)+in-1                                    8d15s06
             iad1=jdenpt(i2eu)+ii3+na*ii2                                8d18s14
            else                                                         8d15s06
             if(icase.eq.1)then                                          8d15s06
              iad1=jdenpt(i2eu)+into-1+iacto(isbx)*(iisx-1+iacto(isx)*   8d18s14
     $           (infrom-1+iacto(isbx)*(iisx-1)))                       8d15s06
             else if(icase.eq.2)then                                     8d15s06
              iad1=jdenpt(i2eu)+iisx-1+iacto(isx)*(into-1+iacto(isbx)*   8d18s14
     $           (iisx-1+iacto(isx)*(infrom-1)))                        8d15s06
             else if(icase.eq.3)then                                     8d15s06
              iad1=jdenpt(i2eu)+into-1+iacto(isbx)*(iisx-1+iacto(isx)*   8d18s14
     $           (iisx-1+iacto(isx)*(infrom-1)))                        8d15s06
             else if(icase.eq.4)then                                     8d15s06
              iad1=jdenpt(i2eu)+iisx-1+iacto(isx)*(into-1+iacto(isbx)*   8d18s14
     $           (infrom-1+iacto(isbx)*(iisx-1)))                       8d15s06
             end if                                                      8d15s06
            end if                                                       8d15s06
            bc(iad1)=bc(iad1)-phs*vv*dbg(3)                              8d19s14
           end do                                                        8d15s06
          else                                                          3d16s23
           nad=iacto(isby)*iacto(isbx)                                  3d16s23
           iad1=iad1+nad*mdenoff                                        3d16s23
           ff=phs                                                       3d16s23
           if(icode.eq.2)ff=ff*2d0                                      3d16s23
           do k1=1,nroot                                                3d16s23
            do k2=1,k1                                                  3d16s23
             bc(iad1)=bc(iad1)+ff*vec1(j,k1)*vec2(i,k2)                 3d16s23
             iad1=iad1+nad                                              3d16s23
            end do                                                      3d16s23
           end do                                                       3d16s23
          end if                                                        6d22s22
         else
c     bra-ket
          do i1=1,norb                                                  8d15s06
           if(itmpa(i1).eq.1)nfrom=i1
           if(itmpa(i1).eq.-1)nto=i1
          end do
          nb=0
          do i1=min(nfrom,nto)+1,max(nfrom,nto)-1
           if(itmpc(i1).eq.1)nb=nb+1
          end do
          if(mod(nb,2).eq.0)then
           phs=1d0
          else
           phs=-1d0
          end if
          isbx=ibc(jsm+nto)                                             8d15s06
          isby=ibc(jsm+nfrom)                                           8d15s06
          if(isbx.ne.isby)then                                          8d15s06
           write(6,*)('symmetry of single does not match '),isbx,isby,
     $         nto,nfrom
           stop
          end if                                                        8d15s06
          into=ibc(jrelo+nto)                                           8d15s06
          infrom=ibc(jrelo+nfrom)                                       8d15s06
           iad1=iden1(isby)+infrom-1+iacto(isby)*(into-1)                6d26s22
          if(l2e)then                                                   6d22s22
           bc(iad1)=bc(iad1)+phs*vv*debug1                               12d19s15
           ix=max(into,infrom)                                           8d15s06
           in=min(into,infrom)                                           8d15s06
           ii1=((ix*(ix-1))/2)+in
           na=(iacto(isbx)*(iacto(isbx)+1))/2                            8d15s06
           do i1=1,nbeta                                                 8d15s06
            isx=ibc(jsm+iborb(i1,ib))                                    8d15s06
            i2eu=ifind2(isbx,isbx,isx,isx,icase)                         8d15s06
            iisx=ibc(jrelo+iborb(i1,ib))                                 8d15s06
            ii2=(iisx*(iisx+1))/2                                        8d15s06
            iad1=jdenpt(i2eu)+ii1-1+na*(ii2-1)                           8d18s14
            bc(iad1)=bc(iad1)+phs*vv*dbg(4)                              8d19s14
           end do                                                        8d15s06
           do i1=1,nalpha                                                8d15s06
            isx=ibc(jsm+iaorb(i1,ia))                                    8d15s06
            i2eu=ifind2(isbx,isbx,isx,isx,icase)                         8d15s06
            iisx=ibc(jrelo+iaorb(i1,ia))                                 8d15s06
            ii2=(iisx*(iisx+1))/2                                        8d15s06
            iad1=jdenpt(i2eu)+ii1-1+na*(ii2-1)                              8d15s06
            bc(iad1)=bc(iad1)+phs*vv*dbg(5)                              8d19s14
            if(isx.eq.isbx)then                                          8d15s06
             i2eu=ifind2(isx,isx,isx,isx,icase)                          8d15s06
             in=min(into,iisx)                                           8d15s06
             ix=max(into,iisx)                                           8d15s06
             ii3=((ix*(ix-1))/2)+in-1                                    8d15s06
             in=min(infrom,iisx)                                         8d15s06
             ix=max(infrom,iisx)                                         8d15s06
             ii2=((ix*(ix-1))/2)+in-1                                    8d15s06
             iad1=jdenpt(i2eu)+ii3+na*ii2                                   8d15s06
            else                                                         8d15s06
             i2eu=ifind2(isbx,isx,isbx,isx,icase)                        8d15s06
             if(icase.eq.1)then                                          8d15s06
              iad1=jdenpt(i2eu)+into-1+iacto(isbx)*(iisx-1+iacto(isx)*      8d15s06
     $           (infrom-1+iacto(isbx)*(iisx-1)))                       8d15s06
             else if(icase.eq.2)then                                     8d15s06
              iad1=jdenpt(i2eu)+iisx-1+iacto(isx)*(into-1+iacto(isbx)*      8d15s06
     $           (iisx-1+iacto(isx)*(infrom-1)))                        8d15s06
             else if(icase.eq.3)then                                     8d15s06
              iad1=jdenpt(i2eu)+into-1+iacto(isbx)*(iisx-1+iacto(isx)*      8d15s06
     $           (iisx-1+iacto(isx)*(infrom-1)))                        8d15s06
             else if(icase.eq.4)then                                     8d15s06
              iad1=jdenpt(i2eu)+iisx-1+iacto(isx)*(into-1+iacto(isbx)*      8d15s06
     $           (infrom-1+iacto(isbx)*(iisx-1)))                       8d15s06
             end if                                                      8d15s06
            end if                                                       8d15s06
            bc(iad1)=bc(iad1)-phs*vv*dbg(6)                              8d19s14
           end do                                                        8d15s06
          else                                                          3d16s23
           nad=iacto(isby)*iacto(isbx)                                  3d16s23
           iad1=iad1+nad*mdenoff                                        3d16s23
           ff=phs                                                       3d16s23
           if(icode.eq.2)ff=2d0                                         3d16s23
           do k1=1,nroot                                                3d16s23
            do k2=1,k1                                                  3d16s23
             bc(iad1)=bc(iad1)+ff*vec1(j,k1)*vec2(i,k2)                 3d16s23
             iad1=iad1+nad                                              3d16s23
            end do                                                      3d16s23
           end do                                                       3d16s23
          end if                                                        6d22s22
         end if
        else if(isum.eq.4.and.l2e)then                                  6d22s22
         if(lprint)write(6,*)('double ')
         if(isuma.eq.0)then
          if(lprint)write(6,*)('bb')
          ip=1
          im=1
          do i1=1,norb                                                  8d15s06
           if(itmpb(i1).eq.1)then
            itmp2(ip,1)=i1
            ip=ip+1
           end if
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
          isx1=ibc(jsm+itmp2(1,1))                                      8d15s06
          isx2=ibc(jsm+itmp2(1,2))                                      8d15s06
          isx3=ibc(jsm+itmp2(2,1))                                      8d15s06
          isx4=ibc(jsm+itmp2(2,2))                                      8d15s06
          iisx1=ibc(jrelo+itmp2(1,1))                                   8d15s06
          iisx2=ibc(jrelo+itmp2(1,2))                                   8d15s06
          iisx3=ibc(jrelo+itmp2(2,1))                                   8d15s06
          iisx4=ibc(jrelo+itmp2(2,2))                                   8d15s06
          i2eu=ifind2(isx1,isx2,isx3,isx4,icase)                        8d15s06
          if(isx1.eq.isx2)then                                          8d15s06
           in=min(iisx1,iisx2)                                          8d15s06
           ix=max(iisx1,iisx2)                                          8d15s06
           ii1=((ix*(ix-1))/2)+in-1                                     8d15s06
           in=min(iisx3,iisx4)                                          8d15s06
           ix=max(iisx3,iisx4)                                          8d15s06
           ii2=((ix*(ix-1))/2)+in-1                                     8d15s06
           na=(iacto(isx1)*(iacto(isx1)+1))/2                           8d15s06
           iad1=jdenpt(i2eu)+ii1+na*ii2                                 8d18s14
          else                                                          8d15s06
           if(icase.eq.1)then                                           8d15s06
            iad1=jdenpt(i2eu)+iisx1-1+iacto(isx1)*(iisx2-1+iacto(isx2)     8d15s06
     $          *(iisx3-1+iacto(isx3)*(iisx4-1)))                       8d15s06
           else if(icase.eq.2)then                                      8d15s06
            iad1=jdenpt(i2eu)+iisx2-1+iacto(isx2)*(iisx1-1+iacto(isx1)     8d15s06
     $          *(iisx4-1+iacto(isx4)*(iisx3-1)))                       8d15s06
           else if(icase.eq.3)then                                      8d15s06
            iad1=jdenpt(i2eu)+iisx1-1+iacto(isx1)*(iisx2-1+iacto(isx2)     8d15s06
     $          *(iisx4-1+iacto(isx4)*(iisx3-1)))                       8d15s06
           else if(icase.eq.4)then                                      8d15s06
            iad1=jdenpt(i2eu)+iisx2-1+iacto(isx2)*(iisx1-1+iacto(isx1)     8d15s06
     $          *(iisx3-1+iacto(isx3)*(iisx4-1)))                       8d15s06
           end if                                                       8d15s06
          end if                                                        8d15s06
          bc(iad1)=bc(iad1)+phs*vv*dbg(7)                               8d19s14
          isx1=ibc(jsm+itmp2(2,1))                                      8d15s06
          isx2=ibc(jsm+itmp2(1,2))                                      8d15s06
          isx3=ibc(jsm+itmp2(1,1))                                      8d15s06
          isx4=ibc(jsm+itmp2(2,2))                                      8d15s06
          iisx1=ibc(jrelo+itmp2(2,1))                                   8d15s06
          iisx2=ibc(jrelo+itmp2(1,2))                                   8d15s06
          iisx3=ibc(jrelo+itmp2(1,1))                                   8d15s06
          iisx4=ibc(jrelo+itmp2(2,2))                                   8d15s06
          i2eu=ifind2(isx1,isx2,isx3,isx4,icase)                        8d15s06
          if(isx1.eq.isx2)then                                          8d15s06
           in=min(iisx1,iisx2)                                          8d15s06
           ix=max(iisx1,iisx2)                                          8d15s06
           ii1=((ix*(ix-1))/2)+in-1                                     8d15s06
           na=(iacto(isx1)*(iacto(isx1)+1))/2                           8d15s06
           in=min(iisx3,iisx4)                                          8d15s06
           ix=max(iisx3,iisx4)                                          8d15s06
           ii2=((ix*(ix-1))/2)+in-1                                     8d15s06
           iad1=jdenpt(i2eu)+ii1+na*ii2                                    8d15s06
          else                                                          8d15s06
           if(icase.eq.1)then                                           8d15s06
            iad1=jdenpt(i2eu)+iisx1-1+iacto(isx1)*(iisx2-1+iacto(isx2)     8d15s06
     $          *(iisx3-1+iacto(isx3)*(iisx4-1)))                       8d15s06
           else if(icase.eq.2)then                                      8d15s06
            iad1=jdenpt(i2eu)+iisx2-1+iacto(isx2)*(iisx1-1+iacto(isx1)     8d15s06
     $           *(iisx4-1+iacto(isx4)*(iisx3-1)))                      8d15s06
           else if(icase.eq.3)then                                      8d15s06
            iad1=jdenpt(i2eu)+iisx1-1+iacto(isx1)*(iisx2-1+iacto(isx2)     8d15s06
     $           *(iisx4-1+iacto(isx4)*(iisx3-1)))                      8d15s06
           else if(icase.eq.4)then                                      8d15s06
            iad1=jdenpt(i2eu)+iisx2-1+iacto(isx2)*(iisx1-1+iacto(isx1)     8d15s06
     $           *(iisx3-1+iacto(isx3)*(iisx4-1)))                      8d15s06
           end if                                                       8d15s06
          end if                                                        8d15s06
          bc(iad1)=bc(iad1)-phs*vv*dbg(8)                               8d19s14
         else if(isumb.eq.0)then
          if(lprint)write(6,*)('aa')
          ip=1
          im=1
          do i1=1,norb                                                  8d15s06
           if(itmpa(i1).eq.1)then
            itmp2(ip,1)=i1
            ip=ip+1
           end if
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
          isx1=ibc(jsm+itmp2(1,1))                                      8d15s06
          isx2=ibc(jsm+itmp2(1,2))                                      8d15s06
          isx3=ibc(jsm+itmp2(2,1))                                      8d15s06
          isx4=ibc(jsm+itmp2(2,2))                                      8d15s06
          iisx1=ibc(jrelo+itmp2(1,1))                                   8d15s06
          iisx2=ibc(jrelo+itmp2(1,2))                                   8d15s06
          iisx3=ibc(jrelo+itmp2(2,1))                                   8d15s06
          iisx4=ibc(jrelo+itmp2(2,2))                                   8d15s06
          i2eu=ifind2(isx1,isx2,isx3,isx4,icase)                        8d15s06
          if(isx1.eq.isx2)then                                          8d15s06
           in=min(iisx1,iisx2)                                          8d15s06
           ix=max(iisx1,iisx2)                                          8d15s06
           ii1=((ix*(ix-1))/2)+in-1                                     8d15s06
           na=(iacto(isx1)*(iacto(isx1)+1))/2                           8d15s06
           in=min(iisx3,iisx4)                                          8d15s06
           ix=max(iisx3,iisx4)                                          8d15s06
           ii2=((ix*(ix-1))/2)+in-1                                     8d15s06
           iad1=jdenpt(i2eu)+ii1+na*ii2                                    8d15s06
          else                                                          8d15s06
           if(icase.eq.1)then                                           8d15s06
            iad1=jdenpt(i2eu)+iisx1-1+iacto(isx1)*(iisx2-1+iacto(isx2)*    8d15s06
     $          (iisx3-1+iacto(isx3)*(iisx4-1)))                        8d15s06
           else if(icase.eq.2)then                                      8d15s06
            iad1=jdenpt(i2eu)+iisx2-1+iacto(isx2)*(iisx1-1+iacto(isx1)*    8d15s06
     $           (iisx4-1+iacto(isx4)*(iisx3-1)))                       8d15s06
           else if(icase.eq.3)then                                      8d15s06
            iad1=jdenpt(i2eu)+iisx1-1+iacto(isx1)*(iisx2-1+iacto(isx2)*    8d15s06
     $           (iisx4-1+iacto(isx4)*(iisx3-1)))                       8d15s06
           else if(icase.eq.4)then                                      8d15s06
            iad1=jdenpt(i2eu)+iisx2-1+iacto(isx2)*(iisx1-1+iacto(isx1)*    8d15s06
     $           (iisx3-1+iacto(isx3)*(iisx4-1)))                       8d15s06
           end if                                                       8d15s06
          end if                                                        8d15s06
          isx1=ibc(jsm+itmp2(2,1))                                      8d15s06
          isx2=ibc(jsm+itmp2(1,2))                                      8d15s06
          isx3=ibc(jsm+itmp2(1,1))                                      8d15s06
          isx4=ibc(jsm+itmp2(2,2))                                      8d15s06
          iisx1=ibc(jrelo+itmp2(2,1))                                   8d15s06
          iisx2=ibc(jrelo+itmp2(1,2))                                   8d15s06
          iisx3=ibc(jrelo+itmp2(1,1))                                   8d15s06
          iisx4=ibc(jrelo+itmp2(2,2))                                   8d15s06
          i2eu=ifind2(isx1,isx2,isx3,isx4,icase)                        8d15s06
          if(isx1.eq.isx2)then                                          8d15s06
           in=min(iisx1,iisx2)                                          8d15s06
           ix=max(iisx1,iisx2)                                          8d15s06
           ii1=((ix*(ix-1))/2)+in-1                                     8d15s06
           na=(iacto(isx1)*(iacto(isx1)+1))/2                           8d15s06
           in=min(iisx3,iisx4)                                          8d15s06
           ix=max(iisx3,iisx4)                                          8d15s06
           ii2=((ix*(ix-1))/2)+in-1                                     8d15s06
           iad2=jdenpt(i2eu)+ii1+na*ii2                                    8d15s06
          else                                                          8d15s06
           if(icase.eq.1)then                                           8d15s06
            iad2=jdenpt(i2eu)+iisx1-1+iacto(isx1)*(iisx2-1+iacto(isx2)*    8d15s06
     $          (iisx3-1+iacto(isx3)*(iisx4-1)))                        8d15s06
           else if(icase.eq.2)then                                      8d15s06
            iad2=jdenpt(i2eu)+iisx2-1+iacto(isx2)*(iisx1-1+iacto(isx1)*    8d15s06
     $           (iisx4-1+iacto(isx4)*(iisx3-1)))                       8d15s06
           else if(icase.eq.3)then                                      8d15s06
            iad2=jdenpt(i2eu)+iisx1-1+iacto(isx1)*(iisx2-1+iacto(isx2)*    8d15s06
     $           (iisx4-1+iacto(isx4)*(iisx3-1)))                       8d15s06
           else if(icase.eq.4)then                                      8d15s06
            iad2=jdenpt(i2eu)+iisx2-1+iacto(isx2)*(iisx1-1+iacto(isx1)*    8d15s06
     $           (iisx3-1+iacto(isx3)*(iisx4-1)))                       8d15s06
           end if                                                       8d15s06
          end if                                                        8d15s06
          bc(iad1)=bc(iad1)+phs*vv*dbg(9)                               8d19s14
          bc(iad2)=bc(iad2)-phs*vv*dbg(10)                              8d19s14
          if(lprint)write(6,*)('integrals '),bc(iad1),bc(iad2)
         else if(isuma.eq.2.and.isumb.eq.2)then
          if(lprint)write(6,*)('ab')
          do i1=1,norb                                                  8d15s06
           if(itmpa(i1).eq.1)nfroma=i1
           if(itmpa(i1).eq.-1)ntoa=i1
          end do
          nb=0
          do i1=min(nfroma,ntoa)+1,max(nfroma,ntoa)-1
           if(itmpc(i1).eq.1)nb=nb+1
          end do
          do i1=1,norb                                                  8d15s06
           if(itmpb(i1).eq.1)nfromb=i1
           if(itmpb(i1).eq.-1)ntob=i1
          end do
          do i1=min(nfromb,ntob)+1,max(nfromb,ntob)-1
           if(itmpe(i1).eq.1)nb=nb+1
          end do
          if(mod(nb,2).eq.0)then
           phs=1d0
          else
           phs=-1d0
          end if
          isfa=ibc(jsm+nfroma)                                          8d15s06
          ista=ibc(jsm+ntoa)                                            8d15s06
          isfb=ibc(jsm+nfromb)                                          8d15s06
          istb=ibc(jsm+ntob)                                            8d15s06
          infroma=ibc(jrelo+nfroma)                                     8d15s06
          intoa=ibc(jrelo+ntoa)                                         8d15s06
          infromb=ibc(jrelo+nfromb)                                     8d15s06
          intob=ibc(jrelo+ntob)                                         8d15s06
          i2eu=ifind2(isfa,ista,isfb,istb,icase)                        8d15s06
          if(isfa.eq.ista)then                                          8d15s06
           in=min(infroma,intoa)                                        8d15s06
           ix=max(infroma,intoa)                                        8d15s06
           ii1=((ix*(ix-1))/2)+in-1                                     8d15s06
           na=(iacto(isfa)*(iacto(isfa)+1))/2                           8d15s06
           in=min(infromb,intob)                                        8d15s06
           ix=max(infromb,intob)                                        8d15s06
           ii2=((ix*(ix-1))/2)+in-1                                     8d15s06
           iad1=jdenpt(i2eu)+ii1+na*ii2                                    8d15s06
          else                                                          8d15s06
           if(icase.eq.1)then                                           8d15s06
            iad1=jdenpt(i2eu)+infroma-1+iacto(isfa)                     8d18s14
     $          *(intoa-1+iacto(ista)*                                  8d18s14
     $          (infromb-1+iacto(isfb)*(intob-1)))                      8d15s06
           else if(icase.eq.2)then                                      8d15s06
            iad1=jdenpt(i2eu)+intoa-1+iacto(ista)*(infroma-1
     $           +iacto(isfa)*                                          8d18s14
     $           (intob-1+iacto(istb)*(infromb-1)))                     8d15s06
           else if(icase.eq.3)then                                      8d15s06
            iad1=jdenpt(i2eu)+infroma-1+iacto(isfa)*(intoa-1            8d18s14
     $           +iacto(ista)*                                          8d18s14
     $           (intob-1+iacto(istb)*(infromb-1)))                     8d15s06
           else if(icase.eq.4)then                                      8d15s06
            iad1=jdenpt(i2eu)+intoa-1+iacto(ista)*(infroma-1            8d18s14
     $           +iacto(isfa)*                                          8d18s14
     $           (infromb-1+iacto(isfb)*(intob-1)))                     8d15s06
           end if                                                       8d15s06
          end if                                                        8d15s06
          bc(iad1)=bc(iad1)+phs*vv*dbg(11)                              8d19s14
         end if
        end if
       end do
      end do
      ibcoff=ism                                                        2d16s07
      return
      end
