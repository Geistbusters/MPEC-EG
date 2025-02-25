c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine psspin(nps,ips,hdig,pham,numb,iaorb,nalpha,
     $                 iborb,nbeta,itmpa,itmpb,itmpc,itmpd,
     $                 itmpe,itmpf,hdigx,nconf,norb,ncsym,nsymb,nsbeta, 11d9s22
     $                 bc,ibc)                                          11d9s22
      implicit real*8 (a-h,o-z)
      integer*8 iarg1,iarg2,iarg3                                       8d22s06
      integer*2 ips(4,1)                                                8d22s06
      integer*1 iaorb,iborb
      logical lprint                                                    5d30s06
      include "common.cas"
      include "common.store"                                            5d30s06
c
c     explicitly form p-space S**2 matrix elements
c
      dimension hdig(1),pham(nps,nps),iaorb(nalpha,1),                  8d22s06
     $          iborb(nbeta,1),itmpa(1),itmpb(1),itmp2(2,2),
     $          itmpc(1),itmpd(1),itmpe(1),itmpf(1),hdigx(1),ncsym(8),  8d14s06
     $          nsbeta(8)                                               12d11s06
      lprint=.false.                                                     12d11s06
      xms=0.5d0*dfloat(nalpha-nbeta)
      nclo=idoub(1)                                                     8d22s06
      do isb=2,nsymb                                                    8d22s06
       nclo=nclo+idoub(isb)                                             8d22s06
      end do                                                            8d22s06
      do i=1,nps
       pham(i,i)=xms*xms+0.5d0*dfloat(2*nclo+nalpha+nbeta)              8d22s06
       ipa=ips(1,i)                                                     8d22s06
       ipb=ips(2,i)                                                     8d22s06
       isbu=ips(3,i)                                                    8d22s06
   87  format('det ',i3,' strings ',3i5)
       nasum=0                                                          8d15s06
       do isb=1,isbu-1                                                  8d15s06
        nasum=nasum+nadet(isb)                                          12d11s06
       end do                                                           8d15s06
       nbsum=0                                                          8d15s06
       do isb=1,nsbeta(isbu)-1                                          8d15s06
        nbsum=nbsum+nbdet(isb)                                          12d11s06
       end do                                                           8d15s06
       ipa=ipa+nasum                                                    8d22s06
       ipb=ipb+nbsum                                                    8d22s06
       isame=0
       do i1=1,nalpha
        do i2=1,nbeta
         if(iborb(i2,ipb).eq.iaorb(i1,ipa))then
          isame=isame+1
          go to 10
         end if
        end do
   10   continue
       end do
   88  format('translated to ',2i5,' using ',2i5,' isame ',i5)
       pham(i,i)=pham(i,i)-dfloat(isame+nclo)                           8d22s06
      end do
      do i=1,nps-1
       ia=ips(1,i)                                                      8d22s06
       ib=ips(2,i)                                                      8d22s06
       isbu=ips(3,i)                                                    8d22s06
       nasum=0                                                          8d15s06
       do isb=1,isbu-1                                                  8d15s06
        nasum=nasum+nadet(isb)                                          12d11s06
       end do                                                           8d15s06
       nbsum=0                                                          8d15s06
       do isb=1,nsbeta(isbu)-1                                          8d15s06
        nbsum=nbsum+nbdet(isb)                                          12d11s06
       end do                                                           8d15s06
       ia0=ia
       ib0=ib
       ia=ia+nasum                                                      8d15s06
       ib=ib+nbsum                                                      8d15s06
       do j=i+1,nps
        ans=0d0
        ja=ips(1,j)                                                     8d22s06
        jb=ips(2,j)                                                     8d22s06
        jsbu=ips(3,j)                                                   8d22s06
        if(lprint)write(6,*)('bra ips: '),(ips(kk,j),kk=1,3)
        if(lprint)write(6,*)('ket ips: '),(ips(kk,i),kk=1,3)
        nasum=0                                                         8d15s06
        do isb=1,jsbu-1                                                 8d15s06
         nasum=nasum+nadet(isb)                                         12d11s06
        end do                                                          8d15s06
        nbsum=0                                                         8d15s06
        do isb=1,nsbeta(jsbu)-1                                         8d15s06
         nbsum=nbsum+nbdet(isb)                                         12d11s06
        end do                                                          8d15s06
        ja0=ja
        jb0=jb
        ja=ja+nasum                                                     8d15s06
        jb=jb+nbsum                                                     8d15s06
        if(lprint)write(6,*)('j,i '),j,i,ia,ib,ja,jb
        isuma=0
        do i1=1,norb                                                    8d14s06
         itmpc(i1)=0
         itmpd(i1)=0
         itmpe(i1)=0
         itmpf(i1)=0
        end do
        do i1=1,nalpha
         itmpc(iaorb(i1,ja))=1
         itmpd(iaorb(i1,ia))=1
        end do
        do i1=1,norb                                                    8d14s06
         itmpa(i1)=itmpc(i1)-itmpd(i1)
         isuma=isuma+iabs(itmpa(i1))
  113    format(4i5)
        end do
        isumb=0
        do i1=1,nbeta
         itmpe(iborb(i1,jb))=1
         itmpf(iborb(i1,ib))=1
        end do
        do i1=1,norb                                                    8d14s06
         itmpb(i1)=itmpe(i1)-itmpf(i1)
         isumb=isumb+iabs(itmpb(i1))
        end do
        isum=isumb+isuma
        if(lprint)write(6,*)('isum '),isum,isuma,isumb
        if(isum.eq.4)then
         if(lprint)write(6,*)('double ')
         if(isuma.eq.2)then
          if(lprint)write(6,*)('ab')
          do i1=1,norb                                                  8d14s06
           if(itmpa(i1).eq.1)then
            nfroma=i1
           end if
           if(itmpa(i1).eq.-1)then
            ntoa=i1
           end if
          end do
          nb=0
          do i1=min(nfroma,ntoa)+1,max(nfroma,ntoa)-1
           if(itmpc(i1).eq.1)nb=nb+1
          end do
          do i1=1,norb                                                  8d14s06
           if(itmpb(i1).eq.1)nfromb=i1
           if(itmpb(i1).eq.-1)ntob=i1
          end do
          do i1=min(nfromb,ntob)+1,max(nfromb,ntob)-1
           if(itmpf(i1).eq.1)nb=nb+1
          end do
          if(mod(nb,2).eq.0)then
           phs=1d0
          else
           phs=-1d0
          end if
          if(lprint)write(6,*)('from '),nfroma,nfromb
          if(lprint)write(6,*)('to '),ntoa,ntob
          if(nfroma.eq.ntob.and.nfromb.eq.ntoa)then
           ans=-phs
          end if
         end if
        end if
        pham(i,j)=ans
        pham(j,i)=ans
       end do
      end do
      iarg1=nps
      return
      end
