c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine pshamppnona1(nps,ips,pham,numb,iaorb,nalpha,           6d22s22
     $                 iborb,nbeta,itmpa,itmpb,ih0e,i2e,                 8d15s06
     $                 itmpc,itmpd,itmpe,itmpf,debug1,debug2,norb,      8d15s06
     $                 ncsym,nsymb,nsbeta,nsbeta2,dbg,nhere,istart,     6d22s22
     $                 ipxder,ipuse,multh,ips2,hdigo,hdig,bc,ibc)       11d9s22
      implicit real*8 (a-h,o-z)
      integer*2 ips(4,1),ips2(4,1)                                      6d23s22
      integer*1 iaorb,iborb,itmpc,itmpd,itmpe,itmpf                     11d14s05
      logical lprint                                                    12d11s06
      include "common.cas"
      include "common.store"                                            8d15s06
c
c     explicitly form p-space hamiltonian matrix elements               10d23s14
c     do only this processors share of rows                             10d23s14
c
      dimension hdig(*),pham(nhere,nps),iaorb(nalpha,1),ipxder(4,8,8,8),6d27s22
     $          iborb(nbeta,1),itmpa(1),itmpb(1),ih0e(1),multh(8,8),    6d22s22
     $          i2e(1),itmp2(2,2),itmpc(1),itmpd(1),nsbeta2(*),         6d22s22
     $          itmpe(1),itmpf(1),ncsym(8),nsbeta(8),dbg(11),hdigo(*)   6d23s22
      data icall/0/
      save icall
      loop=0
      loopx=797
      lprint=.false.                                                    8d19s14
      icall=icall+1
      ism=ibcoff                                                        8d15s06
      irelo=ism+norb                                                    8d15s06
      ibcoff=irelo+norb                                                 8d15s06
      call enough('pshamppnona1.  1',bc,ibc)
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
      i=0                                                               6d22s22
      do isb=1,nsymb                                                    6d22s22
       jsb=nsbeta2(isb)                                                 6d22s22
       do ipa=1,nadet(isb)                                              6d22s22
        do ipb=1,nbdet(jsb)                                             6d22s22
         i=i+1                                                          6d22s22
         ips2(1,i)=ipa                                                  6d23s22
         ips2(2,i)=ipb                                                  6d23s22
         ips2(3,i)=isb                                                  6d23s22
         if(i.ge.istart.and.i.lt.istart+nhere)then                      6d22s22
          isto=i+1-istart                                               6d22s22
          hdigo(isto)=hdig(i)                                           6d23s22
          nasum2=0                                                         8d15s06
          do isbx=1,isb-1                                                  8d15s06
           nasum2=nasum2+nadet(isbx)                                        8d15s06
          end do                                                           8d15s06
          nbsum2=0                                                         8d15s06
          do isbx=1,jsb-1                                                6d22s22
           nbsum2=nbsum2+nbdet(isbx)                                        8d15s06
          end do                                                           8d15s06
          ia=ipa+nasum2                                                    8d15s06
          ib=ipb+nbsum2                                                    8d15s06
          do j=1,nps                                                       10d23s14
           ans=0d0                                                      6d22s22
           jpa=ips(1,j)                                                    8d22s06
           jpb=ips(2,j)                                                    8d22s06
           jsbu=ips(3,j)                                                   8d22s06
           nasum1=0                                                        8d15s06
           do isbx=1,jsbu-1                                                 8d15s06
            nasum1=nasum1+nadet(isbx)                                       8d15s06
           end do                                                          8d15s06
           nbsum1=0                                                        8d15s06
           do isbx=1,nsbeta(jsbu)-1                                         8d15s06
            nbsum1=nbsum1+nbdet(isbx)                                       8d15s06
           end do                                                          8d15s06
           ja=jpa+nasum1                                                   8d15s06
           jb=jpb+nbsum1                                                   8d15s06
c     j is ket, i is bra
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
c     ket-bra
            itmpa(i1)=itmpc(i1)-itmpd(i1)
            isuma=isuma+iabs(itmpa(i1))
           end do
           isumb=0
           do i1=1,nbeta
            itmpe(iborb(i1,jb))=1
            itmpf(iborb(i1,ib))=1
           end do
           do i1=1,norb                                                    8d15s06
c     ket-bra
            itmpb(i1)=itmpe(i1)-itmpf(i1)
            isumb=isumb+iabs(itmpb(i1))
            if(lprint)write(6,*)i1,itmpc(i1),itmpd(i1),itmpa(i1),
     $           itmpe(i1),itmpf(i1),itmpb(i1)
           end do
           isum=isumb+isuma
           if(lprint)write(6,*)('no. difference = '),isum
           if(isum.eq.2)then
            if(lprint)write(6,*)('single ')
            if(isuma.eq.0)then
             if(lprint)write(6,*)('beta single'),i,j
             do i1=1,norb                                                  8d15s06
              if(itmpb(i1).eq.1)nket=i1
              if(itmpb(i1).eq.-1)nbra=i1
             end do
             nb=0
             do i1=min(nket,nbra)+1,max(nket,nbra)-1
              if(itmpe(i1).eq.1)nb=nb+1
             end do
             if(mod(nb,2).eq.0)then
              phs=1d0
             else
              phs=-1d0
             end if
             isbbra=ibc(jsm+nbra)                                             8d15s06
             isbket=ibc(jsm+nket)                                           8d15s06
             inbra=ibc(jrelo+nbra)-1                                    6d22s22
             inket=ibc(jrelo+nket)-1                                    6d22s22
             iad1=ih0e(isbbra)+inbra+iacto(isbbra)*inket                6d22s22
             if(lprint)write(6,*)('from to '),inbra,isbbra,inket,isbket,
     $             ih0e(isbbra),iad1,iad1-ih0e(isbbra),bc(iad1)
             part=bc(iad1)                                              6d22s22
             nrow=iacto(isbbra)*iacto(isbket)                           6d22s22
             iab=inbra+iacto(isbbra)*inket                              6d22s22
             iba=inket+iacto(isbket)*inbra-iab                          6d22s22
             do i1=1,nalpha                                                8d15s06
              isx=ibc(jsm+iaorb(i1,ia))                                    8d15s06
              i2eu=ipxder(1,isbbra,isbket,isx)                          6d22s22
              iisx=ibc(jrelo+iaorb(i1,ia))-1                            6d22s22
              icol=iisx*(iacto(isx)+1)                                  6d22s22
              irow=iab+ipxder(2,isbbra,isbket,isx)*iba                  6d27s22
              nrr=iacto(isx)**2                                         6d27s22
              ix1=irow+nrow*icol                                        6d27s22
              ix2=icol+nrr*irow-ix1                                     6d27s22
              ix12=ix1+ipxder(4,isbbra,isbket,isx)*ix2                  6d27s22
              iad1=i2e(i2eu)+ix12                                       6d27s22
              if(lprint)write(6,*)('alpha j '),bc(iad1)
              part=part+bc(iad1)                                        6d22s22
             end do                                                        8d15s06
             if(lprint)write(6,*)('beta '),ib,nbeta
             do i1=1,nbeta                                                 8d15s06
              isx=ibc(jsm+iborb(i1,ib))                                    8d15s06
              iisx=ibc(jrelo+iborb(i1,ib))-1                                 8d15s06
              i2eu=ipxder(1,isbbra,isbket,isx)                          6d22s22
              icol=iisx*(iacto(isx)+1)                                  6d22s22
              irow=iab+ipxder(2,isbbra,isbket,isx)*iba                  6d27s22
              ix1=irow+nrow*icol                                        6d27s22
              nrr=iacto(isx)**2                                         6d27s22
              ix2=icol+nrr*irow-ix1                                     6d27s22
              iad1=i2e(i2eu)+ix1+ipxder(4,isbbra,isbket,isx)*ix2        6d27s22
              if(lprint)write(6,*)('beta j '),bc(iad1),i2eu,
     $             i2e(i2eu),iad1
              part=part+bc(iad1)                                        6d22s22
              i2eu=ipxder(1,isbbra,isx,isx)                             6d23s22
              ix1=inbra+iacto(isbbra)*iisx                              6d22s22
              ix2=iisx+iacto(isx)*inbra-ix1                             6d27s22
              ix3=iisx+iacto(isx)*inket                                 6d22s22
              ix4=inket+iacto(isbket)*iisx-ix3                          6d22s22
              nrr=iacto(isbbra)*iacto(isx)                              6d22s22
              irow=ix1+ipxder(2,isbbra,isx,isx)*ix2                     6d27s22
              icol=ix3+ipxder(3,isbbra,isx,isx)*ix4                     6d27s22
              ix12=irow+nrr*icol                                        6d27s22
              mrr=iacto(isx)*iacto(isbket)                              6d27s22
              ix21=icol+mrr*irow-ix12                                   6d27s22
              iad1=i2e(i2eu)+ix12+ipxder(4,isbbra,isx,isx)*ix21         6d27s22
              if(lprint)write(6,*)('beta k '),bc(iad1),i2eu,
     $             iad1-i2e(i2eu),iad1
              part=part-bc(iad1)                                        6d22s22
             end do                                                        8d15s06u
             if(lprint)write(6,*)('phs = '),phs
             ans=part*phs
            else
             if(lprint)write(6,*)('alpha single '),i,j
             do i1=1,norb                                                  8d15s06
              if(itmpa(i1).eq.1)nket=i1
              if(itmpa(i1).eq.-1)nbra=i1
             end do
             nb=0
             do i1=min(nket,nbra)+1,max(nket,nbra)-1
              if(itmpc(i1).eq.1)nb=nb+1
             end do
             if(mod(nb,2).eq.0)then
              phs=1d0
             else
              phs=-1d0
             end if
             isbbra=ibc(jsm+nbra)                                             8d15s06
             isbket=ibc(jsm+nket)                                           8d15s06
             inbra=ibc(jrelo+nbra)-1                                    6d22s22
             inket=ibc(jrelo+nket)-1                                    6d22s22
             iad1=ih0e(isbbra)+inbra+iacto(isbbra)*inket                6d22s22
             part=bc(iad1)
             if(lprint)write(6,*)('h0 part '),bc(iad1)
             iab=inbra+iacto(isbbra)*inket                              6d22s22
             iba=inket+iacto(isbket)*inbra-iab                          6d22s22
             nrow=iacto(isbbra)*iacto(isbket)                           6d22s22
             do i1=1,nbeta                                                  8d15s06
              isx=ibc(jsm+iborb(i1,ib))                                     8d15s06
              iisx=ibc(jrelo+iborb(i1,ib))-1                            6d22s22
              i2eu=ipxder(1,isbbra,isbket,isx)                          6d22s22
              icol=iisx*(iacto(isx)+1)                                  6d23s22
              irow=iab+ipxder(2,isbbra,isbket,isx)*iba                  6d27s22
              ix12=irow+nrow*icol                                       6d27s22
              nrr=iacto(isx)**2                                         6d27s22
              ix21=icol+nrr*irow-ix12                                   6d27s22
              iad1=i2e(i2eu)+ix12+ipxder(4,isbbra,isbket,isx)*ix21      6d27s22
              if(lprint)write(6,*)('beta j '),bc(iad1)
              part=part+bc(iad1)                                        6d22s22
             end do                                                         8d15s06
             do i1=1,nalpha                                                 8d15s06
              isx=ibc(jsm+iaorb(i1,ia))                                     8d15s06
              iisx=ibc(jrelo+iaorb(i1,ia))-1                            6d22s22
              icol=iisx*(iacto(isx)+1)                                  6d22s22
              i2eu=ipxder(1,isbbra,isbket,isx)                          6d22s22
              irow=iab+ipxder(2,isbbra,isbket,isx)*iba                  6d28s22
              ix12=irow+nrow*icol                                       6d28s22
              nrr=iacto(isx)**2                                         6d28s22
              ix21=icol+nrr*irow-ix12                                   6d28s22
              iad1=i2e(i2eu)+ix12+ipxder(4,isbbra,isbket,isx)*ix21      6d28s22
              part=part+bc(iad1)                                        6d22s22
              i2eu=ipxder(1,isbbra,isx,isx)                             6d22s22
              ix1=inbra+iacto(isbbra)*iisx                              6d22s22
              ix2=iisx+iacto(isx)*inbra-ix1                              6d22s22
              ix3=iisx+iacto(isx)*inket                                 6d22s22
              ix4=inket+iacto(isbket)*iisx-ix3                          6d22s22
              nrr=iacto(isx)*iacto(isbbra)                              6d22s22
              irow=ix1+ipxder(2,isbbra,isx,isx)*ix2                     6d28s22
              icol=ix3+ipxder(3,isbbra,isx,isx)*ix4                     6d28s22
              ix12=irow+nrr*icol                                        6d28s22
              mrr=iacto(isx)*iacto(isbket)                              6d28s22
              ix21=icol+mrr*irow-ix12                                   6d28s22
              iad1=i2e(i2eu)+ix12+ipxder(4,isbbra,isx,isx)*ix21         6d28s22
              part=part-bc(iad1)                                        6d22s22
             end do                                                        8d15s06
             ans=part*phs
            end if
           else if(isum.eq.4)then
            if(lprint)write(6,*)('double ')
            if(isuma.eq.0)then
             if(lprint)write(6,*)('bb')
             ip=1
             im=1
c     ket-bra
             do i1=1,norb                                                  8d15s06
              if(itmpb(i1).eq.1)then
c     ket
               itmp2(ip,1)=i1
               ip=ip+1
              end if
              if(itmpb(i1).eq.-1)then
c     bra
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
             isxk1=ibc(jsm+itmp2(1,1))                                      8d15s06
             isxb1=ibc(jsm+itmp2(1,2))                                      8d15s06
             isxk2=ibc(jsm+itmp2(2,1))                                      8d15s06
             isxb2=ibc(jsm+itmp2(2,2))                                      8d15s06
             iisk1=ibc(jrelo+itmp2(1,1))-1                              6d22s22
             iisb1=ibc(jrelo+itmp2(1,2))-1                              6d22s22
             iisk2=ibc(jrelo+itmp2(2,1))-1                              6d22s22
             iisb2=ibc(jrelo+itmp2(2,2))-1                              6d22s22
             i2eu=ipxder(1,isxb1,isxk1,isxb2)                           6d22s22
             ix1=iisb1+iacto(isxb1)*iisk1                               6d22s22
             ix2=iisk1+iacto(isxk1)*iisb1-ix1                           6d22s22
             irow=ix1+ipxder(2,isxb1,isxk1,isxb2)*ix2                   6d22s22
             nrow=iacto(isxk1)*iacto(isxb1)                             6d22s22
             ix3=iisb2+iacto(isxb2)*iisk2                               6d22s22
             ix4=iisk2+iacto(isxk2)*iisb2-ix3                           6d22s22
             icol=ix3+ipxder(3,isxb1,isxk1,isxb2)*ix4                   6d22s22
             ix12=irow+nrow*icol                                        6d28s22
             nrr=iacto(isxb2)*iacto(isxk2)                              6d28s22
             ix21=icol+nrr*irow-ix12                                    6d28s22
             iad1=i2e(i2eu)+ix12+ipxder(4,isxb1,isxk1,isxb2)*ix21       7d1s22
             val1=bc(iad1)                                              6d22s22
             i2eu=ipxder(1,isxb1,isxk2,isxb2)                           6d22s22
             ix1=iisb1+iacto(isxb1)*iisk2                               6d22s22
             ix2=iisk2+iacto(isxk2)*iisb1-ix1                           6d22s22
             irow=ix1+ipxder(2,isxb1,isxk2,isxb2)*ix2                   6d22s22
             nrow=iacto(isxk2)*iacto(isxb1)                             6d22s22
             ix3=iisb2+iacto(isxb2)*iisk1                               6d22s22
             ix4=iisk1+iacto(isxk1)*iisb2-ix3                           6d22s22
             icol=ix3+ipxder(3,isxb1,isxk2,isxb2)*ix4                   6d22s22
             ix12=irow+nrow*icol                                        6d28s22
             nrr=iacto(isxb2)*iacto(isxk1)                              6d28s22
             ix21=icol+nrr*irow-ix12                                    6d28s22
             iad2=i2e(i2eu)+ix12+ipxder(4,isxb1,isxk2,isxb2)*ix21       6d28s22
             val2=bc(iad2)                                              6d22s22
             ans=phs*(val1-val2)                                        6d22s22
             if(lprint)write(6,*)('integrals '),val1,val2,phs
            else if(isumb.eq.0)then
             if(lprint)write(6,*)('aa')
             ip=1
             im=1
c     ket-bra
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
             isxk1=ibc(jsm+itmp2(1,1))                                  6d22s22
             isxb1=ibc(jsm+itmp2(1,2))                                  6d22s22
             isxk2=ibc(jsm+itmp2(2,1))                                  6d22s22
             isxb2=ibc(jsm+itmp2(2,2))                                  6d22s22
             iisxk1=ibc(jrelo+itmp2(1,1))-1                             6d22s22
             iisxb1=ibc(jrelo+itmp2(1,2))-1                             6d22s22
             iisxk2=ibc(jrelo+itmp2(2,1))-1                             6d22s22
             iisxb2=ibc(jrelo+itmp2(2,2))-1                             6d22s22
             i2eu=ipxder(1,isxb1,isxk1,isxb2)                           6d22s22
             ix1=iisxb1+iacto(isxb1)*iisxk1                             6d22s22
             ix2=iisxk1+iacto(isxk1)*iisxb1-ix1                         6d22s22
             irow=ix1+ipxder(2,isxb1,isxk1,isxb2)*ix2                   6d22s22
             nrow=iacto(isxb1)*iacto(isxk1)                             6d22s22
             ix3=iisxb2+iacto(isxb2)*iisxk2                             6d22s22
             ix4=iisxk2+iacto(isxk2)*iisxb2-ix3                         6d22s22
             icol=ix3+ipxder(3,isxb1,isxk1,isxb2)*ix4                   6d22s22
             ix12=irow+nrow*icol                                        6d28s22
             nrr=iacto(isxb2)*iacto(isxk2)                              6d28s22
             ix21=icol+nrr*irow-ix12                                    6d28s22
             iad1=i2e(i2eu)+ix12+ipxder(4,isxb1,isxk1,isxb2)*ix21       6d28s22
             i2eu=ipxder(1,isxb1,isxk2,isxb2)                           6d22s22
             ix1=iisxb1+iacto(isxb1)*iisxk2                             6d22s22
             ix2=iisxk2+iacto(isxk2)*iisxb1-ix1                         6d22s22
             irow=ix1+ipxder(2,isxb1,isxk2,isxb2)*ix2                   6d22s22
             nrow=iacto(isxb1)*iacto(isxk2)                             6d22s22
             ix3=iisxb2+iacto(isxb2)*iisxk1                             6d22s22
             ix4=iisxk1+iacto(isxk1)*iisxb2-ix3                         6d22s22
             icol=ix3+ipxder(3,isxb1,isxk2,isxb2)*ix4                   6d22s22
             ix12=irow+nrow*icol                                        6d28s22
             nrr=iacto(isxb2)*iacto(isxk1)                              6d28s22
             ix21=icol+nrr*irow-ix12                                    6d28s22
             iad2=i2e(i2eu)+ix12+ipxder(4,isxb1,isxk2,isxb2)*ix21       6d28s22
             ans=phs*(bc(iad1)-bc(iad2))                                6d22s22
             if(lprint)write(6,*)('integrals '),bc(iad1),bc(iad2),phs,
     $            isxb1,isxk2,isxb2,isxk1,
     $            i2eu,iad2-i2e(i2eu),irow,nrow,icol,
     $            (ipxder(jz,isxb1,isxk2,isxb2),jz=1,4)
            else if(isuma.eq.2.and.isumb.eq.2)then
c     ket-bra
             do i1=1,norb                                                  8d15s06
              if(itmpa(i1).eq.1)nketa=i1
              if(itmpa(i1).eq.-1)nbraa=i1
             end do
             nb=0
             do i1=min(nketa,nbraa)+1,max(nketa,nbraa)-1
              if(itmpc(i1).eq.1)nb=nb+1
             end do
             do i1=1,norb                                                  8d15s06
              if(itmpb(i1).eq.1)nketb=i1
              if(itmpb(i1).eq.-1)nbrab=i1
             end do
             do i1=min(nketb,nbrab)+1,max(nketb,nbrab)-1
              if(itmpe(i1).eq.1)nb=nb+1
             end do
             if(mod(nb,2).eq.0)then
              phs=1d0
             else
              phs=-1d0
             end if
             iska=ibc(jsm+nketa)                                        6d22s22
             isba=ibc(jsm+nbraa)                                        6d22s22
             iskb=ibc(jsm+nketb)                                        6d22s22
             isbb=ibc(jsm+nbrab)                                        6d22s22
             inketa=ibc(jrelo+nketa)-1                                  6d22s22
             inbraa=ibc(jrelo+nbraa)-1                                  6d22s22
             inketb=ibc(jrelo+nketb)-1                                  6d22s22
             inbrab=ibc(jrelo+nbrab)-1                                  6d22s22
             if(lprint)write(6,*)('ab'),
     $            (ipxder(jq,isba,iska,isbb),jq=1,4)
             i2eu=ipxder(1,isba,iska,isbb)                              6d22s22
             ix1=inbraa+iacto(isba)*inketa                              6d22s22
             ix2=inketa+iacto(iska)*inbraa-ix1                          6d22s22
             irow=ix1+ipxder(2,isba,iska,isbb)*ix2                      6d22s22
             nrow=iacto(isba)*iacto(iska)                               6d22s22
             ix3=inbrab+iacto(isbb)*inketb                              6d22s22
             ix4=inketb+iacto(iskb)*inbrab-ix3                          6d22s22
             icol=ix3+ipxder(3,isba,iska,isbb)*ix4                      6d22s22
             ix12=irow+nrow*icol                                        6d28s22
             nrr=iacto(isbb)*iacto(iskb)                                6d28s22
             ix21=icol+nrr*irow-ix12                                    6d28s22
             iad1=i2e(i2eu)+ix12+ipxder(4,isba,iska,isbb)*ix21          6d28s22
             ans=phs*bc(iad1)                                           6d22s22
             if(lprint)write(6,*)phs,bc(iad1)
            end if
           end if
           if(lprint)write(6,*)('set psham to '),ans,i,j
           pham(isto,j)=ans                                                10d23s14
          end do
         end if                                                         6d22s22
        end do                                                          6d22s22
       end do                                                           6d22s22
      end do
      ibcoff=ism                                                        2d16s07
      return
      end
