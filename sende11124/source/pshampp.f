c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine pshampp(nps,ips,hdig,pham,numb,iaorb,nalpha,
     $                 iborb,nbeta,itmpa,itmpb,ih0e,i2e,                 8d15s06
     $                 itmpc,itmpd,itmpe,itmpf,debug1,debug2,norb,      8d15s06
     $                ncsym,nsymb,nsbeta,dbg,nhere,istart,ldebug,bc,ibc)11d9s22
      implicit real*8 (a-h,o-z)
      integer*2 ips(4,1)                                                8d22s06
      integer*1 iaorb,iborb,itmpc,itmpd,itmpe,itmpf                     11d14s05
      logical lprint,ldebug                                                    12d11s06
      include "common.cas"
      include "common.store"                                            8d15s06
c
c     explicitly form p-space hamiltonian matrix elements               10d23s14
c     do only this processors share of rows                             10d23s14
c
      dimension hdig(1),pham(nhere,nps),iaorb(nalpha,1),                10d23s14
     $          iborb(nbeta,1),itmpa(1),itmpb(1),ih0e(1),                8d15s06
     $          i2e(1),itmp2(2,2),itmpc(1),itmpd(1),                    8d15s06
     $          itmpe(1),itmpf(1),ncsym(8),nsbeta(8),dbg(11)            8d19s14
      data icall/0/
      save icall
      lprint=ldebug
      icall=icall+1
      if(ldebug)then
       write(6,*)('what we have for h0e: ')
       call prntm2(bc(ih0e(1)),iacto,iacto,iacto)
       if(nsymb.eq.1)call printat(bc(i2e(1)),iacto,bc(ibcoff))
      end if
      ism=ibcoff                                                        8d15s06
      irelo=ism+norb                                                    8d15s06
      ibcoff=irelo+norb                                                 8d15s06
      call enough('pshampp.  1',bc,ibc)
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
      do i=istart,istart+nhere-1                                        10d23s14
       isto=i+1-istart                                                  10d23s14
       if(ldebug)write(6,*)('hdig '),i,hdig(i),debug2,isto
       pham(isto,i)=hdig(i)*debug2                                      10d23s14
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
       if(lprint)write(6,*)('for ket'),i,(' alpha:beta '),
     $      (iaorb(i1,ia),i1=1,nalpha),
     $          (':'),(iborb(i1,ib),i1=1,nbeta)
       do j=1,nps                                                       10d23s14
        if(j.ne.i)then                                                  10d23s14
         if(lprint)write(6,*)('j,i '),j,i
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
         do i1=1,norb                                                    8d15s06
          itmpa(i1)=itmpc(i1)-itmpd(i1)
          isuma=isuma+iabs(itmpa(i1))
         end do
         isumb=0
         do i1=1,nbeta
          itmpe(iborb(i1,jb))=1
          itmpf(iborb(i1,ib))=1
         end do
         do i1=1,norb                                                    8d15s06
          itmpb(i1)=itmpe(i1)-itmpf(i1)
          isumb=isumb+iabs(itmpb(i1))
          if(lprint)write(6,*)i1,itmpc(i1),itmpd(i1),itmpa(i1),
     $         itmpe(i1),itmpf(i1),itmpb(i1)
         end do
         isum=isumb+isuma
         ans=0d0
         if(lprint)write(6,*)('no. difference = '),isum
         if(isum.eq.2)then
          if(lprint)write(6,*)('single ')
          if(isuma.eq.0)then
           if(lprint)write(6,*)('beta single'),i,j
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
            stop 'pshampp'                                              9d11s24
           end if                                                        8d15s06
           into=ibc(jrelo+nto)                                           8d15s06
           infrom=ibc(jrelo+nfrom)                                       8d15s06
           if(lprint)write(6,*)('from to '),infrom,into,isbx,ih0e(isbx)
           iad1=ih0e(isbx)+into-1+iacto(isbx)*(infrom-1)                 8d15s06
           if(lprint)write(6,*)('h0 part '),bc(iad1)
           part=bc(iad1)                                                 6d22s22
           ix=max(into,infrom)                                           8d15s06
           in=min(into,infrom)                                           8d15s06
           ii1=((ix*(ix-1))/2)+in
           na=(iacto(isbx)*(iacto(isbx)+1))/2                            8d15s06
           do i1=1,nalpha                                                8d15s06
            isx=ibc(jsm+iaorb(i1,ia))                                    8d15s06
            i2eu=ifind2(isbx,isbx,isx,isx,icase)                         8d15s06
            iisx=ibc(jrelo+iaorb(i1,ia))                                 8d15s06
            ii2=(iisx*(iisx+1))/2                                        8d15s06
            iad1=i2e(i2eu)+ii1-1+na*(ii2-1)                              8d15s06
            if(lprint)write(6,*)('alpha j '),bc(iad1)
            part=part+bc(iad1)                                           6d22s22
           end do                                                        8d15s06
           if(lprint)write(6,*)('beta '),ib
           do i1=1,nbeta                                                 8d15s06
            isx=ibc(jsm+iborb(i1,ib))                                    8d15s06
            i2eu=ifind2(isbx,isbx,isx,isx,icase)                         8d15s06
            iisx=ibc(jrelo+iborb(i1,ib))                                 8d15s06
            ii2=(iisx*(iisx+1))/2                                        8d15s06
            iad1=i2e(i2eu)+ii1-1+na*(ii2-1)                              8d15s06
            if(lprint)write(6,*)('beta j'),bc(iad1)
            part=part+bc(iad1)                                           6d22s22
            i2eu=ifind2(isbx,isx,isbx,isx,icase)                         8d15s06
            if(isbx.eq.isx)then                                          8d15s06
             in=min(iisx,into)                                           8d15s06
             ix=max(iisx,into)                                           8d15s06
             ii3=((ix*(ix-1))/2)+in-1                                    8d15s06
             in=min(iisx,infrom)                                         8d15s06
             ix=max(iisx,infrom)                                         8d15s06
             ii2=((ix*(ix-1))/2)+in-1                                    8d15s06
             iad1=i2e(i2eu)+ii3+na*ii2                                   8d15s06
            else                                                         8d15s06
             if(icase.eq.1)then                                          8d15s06
              iad1=i2e(i2eu)+into-1+iacto(isbx)*(iisx-1+iacto(isx)*      8d15s06
     $           (infrom-1+iacto(isbx)*(iisx-1)))                       8d15s06
             else if(icase.eq.2)then                                     8d15s06
              iad1=i2e(i2eu)+iisx-1+iacto(isx)*(into-1+iacto(isbx)*      8d15s06
     $           (iisx-1+iacto(isx)*(infrom-1)))                        8d15s06
             else if(icase.eq.3)then                                     8d15s06
              iad1=i2e(i2eu)+into-1+iacto(isbx)*(iisx-1+iacto(isx)*      8d15s06
     $           (iisx-1+iacto(isx)*(infrom-1)))                        8d15s06
             else if(icase.eq.4)then                                     8d15s06
              iad1=i2e(i2eu)+iisx-1+iacto(isx)*(into-1+iacto(isbx)*      8d15s06
     $           (infrom-1+iacto(isbx)*(iisx-1)))                       8d15s06
             end if                                                      8d15s06
            end if                                                       8d15s06
            if(lprint)write(6,*)('beta k'),bc(iad1)
            part=part-bc(iad1)                                           6d22s22
           end do                                                        8d15s06
           if(lprint)write(6,*)('phs '),phs,part
           ans=part*phs
          else
           if(lprint)write(6,*)('alpha single '),i,j
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
     $          nto,nfrom
            stop 'pshampp'                                              9d11s24
           end if                                                        8d15s06
           into=ibc(jrelo+nto)                                           8d15s06
           infrom=ibc(jrelo+nfrom)                                       8d15s06
           iad1=ih0e(isbx)+into-1+iacto(isbx)*(infrom-1)                 8d15s06
           part=bc(iad1)
           if(lprint)write(6,*)('h0 part '),part
           ix=max(into,infrom)                                           8d15s06
           in=min(into,infrom)                                           8d15s06
           ii1=((ix*(ix-1))/2)+in
           na=(iacto(isbx)*(iacto(isbx)+1))/2                            8d15s06
           do i1=1,nbeta                                                 8d15s06
            isx=ibc(jsm+iborb(i1,ib))                                    8d15s06
            i2eu=ifind2(isbx,isbx,isx,isx,icase)                         8d15s06
            iisx=ibc(jrelo+iborb(i1,ib))                                 8d15s06
            ii2=(iisx*(iisx+1))/2                                        8d15s06
            iad1=i2e(i2eu)+ii1-1+na*(ii2-1)                              8d15s06
            part=part+bc(iad1)                                           6d23s22
            if(lprint)write(6,*)('beta j'),bc(iad1)
           end do                                                        8d15s06
           do i1=1,nalpha                                                8d15s06
            isx=ibc(jsm+iaorb(i1,ia))                                    8d15s06
            i2eu=ifind2(isbx,isbx,isx,isx,icase)                         8d15s06
            iisx=ibc(jrelo+iaorb(i1,ia))                                 8d15s06
            ii2=(iisx*(iisx+1))/2                                        8d15s06
            iad1=i2e(i2eu)+ii1-1+na*(ii2-1)                              8d15s06
            if(lprint)write(6,*)('alpha j'),bc(iad1)
            part=part+bc(iad1)                                           6d23s22
            if(isx.eq.isbx)then                                          8d15s06
             i2eu=ifind2(isx,isx,isx,isx,icase)                          8d15s06
             in=min(into,iisx)                                           8d15s06
             ix=max(into,iisx)                                           8d15s06
             ii3=((ix*(ix-1))/2)+in-1                                    8d15s06
             in=min(infrom,iisx)                                         8d15s06
             ix=max(infrom,iisx)                                         8d15s06
             ii2=((ix*(ix-1))/2)+in-1                                    8d15s06
             iad1=i2e(i2eu)+ii3+na*ii2                                   8d15s06
            else                                                         8d15s06
             i2eu=ifind2(isbx,isx,isbx,isx,icase)                        8d15s06
             if(icase.eq.1)then                                          8d15s06
              iad1=i2e(i2eu)+into-1+iacto(isbx)*(iisx-1+iacto(isx)*      8d15s06
     $            (infrom-1+iacto(isbx)*(iisx-1)))                       8d15s06
             else if(icase.eq.2)then                                     8d15s06
              iad1=i2e(i2eu)+iisx-1+iacto(isx)*(into-1+iacto(isbx)*      8d15s06
     $           (iisx-1+iacto(isx)*(infrom-1)))                        8d15s06
             else if(icase.eq.3)then                                     8d15s06
              iad1=i2e(i2eu)+into-1+iacto(isbx)*(iisx-1+iacto(isx)*      8d15s06
     $           (iisx-1+iacto(isx)*(infrom-1)))                        8d15s06
             else if(icase.eq.4)then                                     8d15s06
              iad1=i2e(i2eu)+iisx-1+iacto(isx)*(into-1+iacto(isbx)*      8d15s06
     $           (infrom-1+iacto(isbx)*(iisx-1)))                       8d15s06
             end if                                                      8d15s06
            end if                                                       8d15s06
            if(lprint)write(6,*)('alpha k'),bc(iad1)
            part=part-bc(iad1)                                           6d23s22
           end do                                                        8d15s06
           ans=part*phs
          end if
         else if(isum.eq.4)then
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
            iad1=i2e(i2eu)+ii1+na*ii2                                    8d15s06
            val1=bc(iad1)                                                8d15s06
           else                                                          8d15s06
            if(icase.eq.1)then                                           8d15s06
             iad1=i2e(i2eu)+iisx1-1+iacto(isx1)*(iisx2-1+iacto(isx2)     8d15s06
     $          *(iisx3-1+iacto(isx3)*(iisx4-1)))                       8d15s06
             val1=bc(iad1)                                               8d15s06
            else if(icase.eq.2)then                                      8d15s06
             iad1=i2e(i2eu)+iisx2-1+iacto(isx2)*(iisx1-1+iacto(isx1)     8d15s06
     $          *(iisx4-1+iacto(isx4)*(iisx3-1)))                       8d15s06
             val1=bc(iad1)                                               8d15s06
            else if(icase.eq.3)then                                      8d15s06
             iad1=i2e(i2eu)+iisx1-1+iacto(isx1)*(iisx2-1+iacto(isx2)     8d15s06
     $          *(iisx4-1+iacto(isx4)*(iisx3-1)))                       8d15s06
             val1=bc(iad1)                                               8d15s06
            else if(icase.eq.4)then                                      8d15s06
             iad1=i2e(i2eu)+iisx2-1+iacto(isx2)*(iisx1-1+iacto(isx1)     8d15s06
     $          *(iisx3-1+iacto(isx3)*(iisx4-1)))                       8d15s06
             val1=bc(iad1)                                               8d15s06
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
            iad1=i2e(i2eu)+ii1+na*ii2                                    8d15s06
            val2=bc(iad1)                                                8d15s06
           else                                                          8d15s06
            if(icase.eq.1)then                                           8d15s06
             iad1=i2e(i2eu)+iisx1-1+iacto(isx1)*(iisx2-1+iacto(isx2)     8d15s06
     $          *(iisx3-1+iacto(isx3)*(iisx4-1)))                       8d15s06
             val2=bc(iad1)                                               8d15s06
            else if(icase.eq.2)then                                      8d15s06
             iad1=i2e(i2eu)+iisx2-1+iacto(isx2)*(iisx1-1+iacto(isx1)     8d15s06
     $           *(iisx4-1+iacto(isx4)*(iisx3-1)))                      8d15s06
             val2=bc(iad1)                                               8d15s06
            else if(icase.eq.3)then                                      8d15s06
             iad1=i2e(i2eu)+iisx1-1+iacto(isx1)*(iisx2-1+iacto(isx2)     8d15s06
     $           *(iisx4-1+iacto(isx4)*(iisx3-1)))                      8d15s06
             val2=bc(iad1)                                               8d15s06
            else if(icase.eq.4)then                                      8d15s06
             iad1=i2e(i2eu)+iisx2-1+iacto(isx2)*(iisx1-1+iacto(isx1)     8d15s06
     $           *(iisx3-1+iacto(isx3)*(iisx4-1)))                      8d15s06
             val2=bc(iad1)                                               8d15s06
            end if                                                       8d15s06
           end if                                                        8d15s06
           ans=phs*(val1*dbg(7)-val2*dbg(8))                                           8d15s06
           if(lprint)write(6,*)('integrals '),val1,val2,phs
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
            iad1=i2e(i2eu)+ii1+na*ii2                                    8d15s06
           else                                                          8d15s06
            if(icase.eq.1)then                                           8d15s06
             iad1=i2e(i2eu)+iisx1-1+iacto(isx1)*(iisx2-1+iacto(isx2)*    8d15s06
     $          (iisx3-1+iacto(isx3)*(iisx4-1)))                        8d15s06
            else if(icase.eq.2)then                                      8d15s06
             iad1=i2e(i2eu)+iisx2-1+iacto(isx2)*(iisx1-1+iacto(isx1)*    8d15s06
     $           (iisx4-1+iacto(isx4)*(iisx3-1)))                       8d15s06
            else if(icase.eq.3)then                                      8d15s06
             iad1=i2e(i2eu)+iisx1-1+iacto(isx1)*(iisx2-1+iacto(isx2)*    8d15s06
     $           (iisx4-1+iacto(isx4)*(iisx3-1)))                       8d15s06
            else if(icase.eq.4)then                                      8d15s06
             iad1=i2e(i2eu)+iisx2-1+iacto(isx2)*(iisx1-1+iacto(isx1)*    8d15s06
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
            iad2=i2e(i2eu)+ii1+na*ii2                                    8d15s06
           else                                                          8d15s06
            if(icase.eq.1)then                                           8d15s06
             iad2=i2e(i2eu)+iisx1-1+iacto(isx1)*(iisx2-1+iacto(isx2)*    8d15s06
     $          (iisx3-1+iacto(isx3)*(iisx4-1)))                        8d15s06
            else if(icase.eq.2)then                                      8d15s06
             iad2=i2e(i2eu)+iisx2-1+iacto(isx2)*(iisx1-1+iacto(isx1)*    8d15s06
     $           (iisx4-1+iacto(isx4)*(iisx3-1)))                       8d15s06
            else if(icase.eq.3)then                                      8d15s06
             iad2=i2e(i2eu)+iisx1-1+iacto(isx1)*(iisx2-1+iacto(isx2)*    8d15s06
     $           (iisx4-1+iacto(isx4)*(iisx3-1)))                       8d15s06
            else if(icase.eq.4)then                                      8d15s06
             iad2=i2e(i2eu)+iisx2-1+iacto(isx2)*(iisx1-1+iacto(isx1)*    8d15s06
     $           (iisx3-1+iacto(isx3)*(iisx4-1)))                       8d15s06
            end if                                                       8d15s06
           end if                                                        8d15s06
           ans=phs*(bc(iad1)*dbg(9)-bc(iad2)*dbg(10))                                   8d15s06
           if(lprint)write(6,*)('integrals '),bc(iad1),bc(iad2),phs
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
            iad1=i2e(i2eu)+ii1+na*ii2                                    8d15s06
           else                                                          8d15s06
            if(icase.eq.1)then                                           8d15s06
             iad1=i2e(i2eu)+infroma-1+iacto(isfa)*(intoa-1+iacto(ista)*  8d15s06
     $          (infromb-1+iacto(isfb)*(intob-1)))                      8d15s06
            else if(icase.eq.2)then                                      8d15s06
             iad1=i2e(i2eu)+intoa-1+iacto(ista)*(infroma-1+iacto(isfa)*  8d15s06
     $           (intob-1+iacto(istb)*(infromb-1)))                     8d15s06
            else if(icase.eq.3)then                                      8d15s06
             iad1=i2e(i2eu)+infroma-1+iacto(isfa)*(intoa-1+iacto(ista)*  8d15s06
     $           (intob-1+iacto(istb)*(infromb-1)))                     8d15s06
            else if(icase.eq.4)then                                      8d15s06
             iad1=i2e(i2eu)+intoa-1+iacto(ista)*(infroma-1+iacto(isfa)*  8d15s06
     $           (infromb-1+iacto(isfb)*(intob-1)))                     8d15s06
            end if                                                       8d15s06
           end if                                                        8d15s06
           ans=phs*bc(iad1)*dbg(11)                                      8d19s14
           if(lprint)write(6,*)phs,bc(iad1)
          end if
         end if
         if(lprint)write(6,*)('set psham to '),ans,i,j
         pham(isto,j)=ans                                                10d23s14
        end if                                                          10d23s14
       end do
      end do
      ibcoff=ism                                                        2d16s07
      return
      end
