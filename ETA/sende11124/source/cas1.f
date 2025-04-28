c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine cas1(nsym,spinm,nelec,nalpha,numa,iaorb,na1p,na1m,iaa1,
     $                iaa1m,na2p,na2m,iaa2,iaa2m,nbeta,numb,iborb,nb1p,
     $                nb1m,ibb1,ibb1m,nb2p,nb2m,ibb2,ibb2m,nconf,bc,ibc)11d9s22
      implicit real*8 (a-h,o-z)
      integer spinm
c
c     generate indicies for computing fci hamiltonian
c
      include "common.store"
      include "common.hf"
      include "common.cas"
      parameter (id=1000,ido=16)
c
c     norb is the number of orbitals
c     nspin is spin multiplicity. (2S+1)
c     nelec is number of electrons
c
      norb=iacto(1)
      nspin=spinm
      spin=0.5d0*(dfloat(nspin)-1d0)
      if(mod(nelec+nspin+1,2).ne.0)then
       write(6,*)('spin is incompatible with number of electrons ')
       stop
      end if
      nx=int(spin+0.01d0)
      write(6,*)('nx = '),nx
      if(mod(nelec,2).eq.0)then
       nalpha=(nelec/2)+nx
       nbeta=nelec-nalpha
      else
       nalpha=((nelec+1)/2)+nx-1
       nbeta=nelec-nalpha
      end if
      write(6,*)('number of alpha electrons '),nalpha
      write(6,*)('number of beta electrons '),nbeta
      if(nbeta.lt.0)then
       write(6,*)('number of electrons incompatible with spin')
       stop
      end if
      check=0.5d0*(nalpha-nbeta)
      write(6,*)('check = '),check
      if(abs(check-spin).gt.1d-6)then
       write(6,*)('check ne spin '),check,spin
       stop
      end if
      write(6,*)('generating alpha strings ')
      write(6,*)('norb = '),norb
      itop=1
      ibot=1
      ix=max(nalpha,norb-nalpha)
      in=min(nalpha,norb-nalpha)
      do i=ix+1,norb
       itop=itop*i
      end do
      do i=1,in
       ibot=ibot*i
      end do
      write(6,*)('numerator '),itop
      write(6,*)('denominator '),ibot
      nconf=itop/ibot
      write(6,*)('number of alpha configurations '),nconf
      iaorb=ibcoff
      needa=nconf*nalpha
      need8=needa/8
      if(need8*8.lt.needa)need8=need8+1
      ibcoff=iaorb+need8
      ialpha=ibcoff
      needa=nconf*norb
      need8=needa/8
      if(need8*8.lt.needa)need8=need8+1
      ibcoff=ialpha+need8
      ix=ibcoff
      ibcoff=ix+norb
      iaa1=ibcoff
      navail=(maxbc-iaa1-1)/4
      ibcoff=iaa1+navail
      write(6,*)('iaorb,... '),iaorb,iaa1,ibcoff,navail
      call enough('cas1.  1',bc,ibc)
      call sgen(ibc(ialpha),nconf,norb,norb,nalpha,nconf2,ibc(iaorb))
      if(nconf.ne.nconf2)then
       write(6,*)('error generating configurations '),nconf,nconf2
       stop
      end if
      call m1gen(ibc(ialpha),nconf,norb,norb,nconf,na1p,
     $     ibc(iaa1),navail,bc(ix),2)
      iaa1m=iaa1+4*na1p
      write(6,*)('iaa1,iaa1m '),iaa1,iaa1m
      navail=navail-na1p
      call m1gen(ibc(ialpha),nconf,norb,norb,nconf,na1m,
     $     ibc(iaa1m),navail,bc(ix),1)
      ibcoff=iaa1m+na1m*4
      write(6,*)('ibcoff '),ibcoff
      write(6,*)('singles '),na1p,na1m
      numa=nconf
      write(6,*)('generating beta strings')
      itop=1
      ibot=1
      ix=max(nbeta,norb-nbeta)
      in=min(nbeta,norb-nbeta)
      do i=ix+1,norb
       itop=itop*i
      end do
      do i=1,in
       ibot=ibot*i
      end do
      write(6,*)('numerator '),itop
      write(6,*)('denominator '),ibot
      nconf=itop/ibot
      write(6,*)('number of beta configurations '),nconf
      ibb1=ibcoff
      need=nconf*(nconf-1)/2
      needs1=need
      ibcoff=ibb1+need*4
      iborb=ibcoff
      needa=nconf*nbeta
      need8=needa/8
      if(need8*8.lt.needa)need8=need8+1
      ibcoff=iborb+need8
      ibeta=ibcoff
      needa=nconf*norb
      need8=needa/8
      if(need8*8.lt.needa)need8=need8+1
      ibcoff=ibeta+need8
      ix=ibcoff
      ibcoff=ix+norb
      call enough('cas1.  2',bc,ibc)
      call sgen(ibc(ibeta),nconf,norb,norb,nbeta,nconf2,ibc(iborb))
      if(nconf.ne.nconf2)then
       write(6,*)('error generating configurations '),nconf,nconf2
       stop
      end if
      call m12gen(ibc(ibeta),nconf,norb,norb,nconf,nb1p,nb1m,
     $     ibc(ibb1),need,bc(ix))
      ibb1m=ibb1+nb1m*4-4
      write(6,*)('m location '),loc(ibc(ibb1m))
      ibb2m=ibb2+nb2m-1
      nb1m=need+1-nb1m
      nb2m=need+1-nb2m
      write(6,*)('singles '),nb1p,nb1m
      numb=nconf
      ibcoff=ibeta
      nconf=numa*numb
      write(6,*)('total number of configurations = '),nconf
      return
      end
