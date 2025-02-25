c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine cas1b(nsym,spinm,nelec,nalpha,nbeta,nproc,nowpro,
     $     ialpha,iaorb,ibeta,iborb,nconfa,nconfb,m1a,m2a,m1b,m2b,
     $     ma1,ma2,mb1,mb2,icall,m12sym,m1c,ma1c,multh,idata1,idata2,   8d12s14
     $     idatac,idatb1,idatb2,idatbc,lprint,mst,mnz,nz,iflag,bc,ibc)  11d9s22
      implicit real*8 (a-h,o-z)
      external second
      integer spinm
      integer*8 ma1(2,8),ma2(2,8),mb1(2,8),mb2(2,8),ma1c(2,36,2),iarg1,
     $     iarg2,iarg3,iarg4,itop                                       4d30s18
      integer*2 idwstest,idwsdcd2(4)                                    12d19s06
      integer*8 idwsdcd8                                                12d19s06
      logical lprint                                                    4d20s18
      equivalence (idwsdcd8,idwsdcd2)                                   12d19s06
      dimension m1(8),m2(8),m12sym(8,2,2),m1c(36,2),idata1(8),idata2(8),8d24s06
     $     idatac(36),il1(8),ih1(8),il2(8),ih2(8),ilc(36),ihc(36),      8d12s14
     $     idatb1(8),idatb2(8),idatbc(36),mst(8,8,2),mnz(3,64,36,2),    4d30s18
     $     nz(36,2)                                                     4d30s18
c
c     storage of integers.
c     sgen: maximum number of orbitals             integer*1=128
c     number of alpha or beta configs in sym block integer*2=32768
c
c     generate indicies for computing fci hamiltonian
c
      include "common.store"
      include "common.hf"
      include "common.cas"
      parameter (id=1000,ido=16)
      data jcall/0/                                                     5d2s18
      save                                                              5d2s18
      jcall=jcall+1                                                     5d2s18
c
c     norb is the number of orbitals
c     nspin is spin multiplicity. (2S+1)
c     nelec is number of electrons
c
      tsgen=0d0
      tsort=0d0
      tcn12=0d0
      tgn12=0d0
      tsrt=0d0
      call second(time2)
      call second(time1)
      tovr=time1-time2
      norb=iacto(1)                                                     8d8s06
      do i=2,nsymb                                                      8d8s06
       norb=norb+iacto(i)                                               8d8s06
      end do                                                            8d8s06
      nspin=spinm
      spin=0.5d0*(dfloat(nspin)-1d0)
      if(mod(nelec+nspin+1,2).ne.0)then
       write(6,*)('norb = '),norb
       write(6,*)('nspin = '),nspin
       write(6,*)('nelec = '),nelec
       write(6,*)('ibcoff = '),ibcoff
       write(6,*)('nsym = '),nsym                                        8d8s06
       write(6,*)('S = '),spin
       write(6,*)('spin is incompatible with number of electrons ')
       stop
      end if
      nos=int(spin*2d0+0.01d0)                                          5d2s07
      nalpha=(nos+nelec)/2                                              8d15s14
      nbeta=nelec-nalpha                                                8d15s14
      if(nbeta.lt.0)then
       write(6,*)('number of alpha electrons '),nalpha
       write(6,*)('number of beta electrons '),nbeta
       write(6,*)('number of electrons incompatible with spin')
       stop
      end if
      check=0.5d0*(nalpha-nbeta)
      if(abs(check).gt.spin+1d-6)then                                   5d17s18
       write(6,*)('ms gt spin '),check,spin                             5d17s18
       call dws_sync                                                    5d17s18
       call dws_finalize                                                5d17s18
       stop
      end if
 1101 format(' symmetry order will be ',8i2)                             8d8s06
      itop=1
      ibot=1
      nalphax=nalpha
      norbx=norb
      ix=max(nalphax,norbx-nalphax)
      in=min(nalphax,norbx-nalphax)
      do i=ix+1,norbx
       itop=itop*i
      end do
      do i=1,in
       ibot=ibot*i
      end do
      nconf=itop/ibot
      idwstest=nconf                                                    5d30s06
      nconfa=nconf
      iaorb=ibcoff
      needa=nconf*nalpha
      need8=needa/8
      if(need8*8.lt.needa)need8=need8+1
      ibcoff=iaorb+need8
      need8a=need8                                                      8d8s06
      ialpha=ibcoff
      needa=nconf*norb
      need8=needa/8
      if(need8*8.lt.needa)need8=need8+1
      ibcoff=ialpha+need8
      call enough('cas1b.  1',bc,ibc)
      call sgen(ibc(ialpha),nconf,norb,norb,nalpha,nconf2,ibc(iaorb))
      call second(time2)
      telap=time2-time1-tovr
      tsgen=tsgen+telap
      if(nconf.ne.nconf2)then
       write(6,*)('for alpha electrons ')
       write(6,*)('error generating configurations '),nconf,nconf2
       write(6,*)('nalphax = '),nalphax
       write(6,*)('norbx = '),norbx
       write(6,*)('spin '),spin
       write(6,*)('nos, nelec '),nos,nelec
       write(6,*)('nalpha,nbeta '),nalpha,nbeta
       stop
      end if
      if(nsymb.ne.1)then
       itmp=ibcoff
       itmp2=itmp+need8
       ibcoff=itmp2+need8a                                              8d8s06
       call enough('cas1b.  2',bc,ibc)
       call sortdet(ibc(ialpha),ibc(iaorb),norb,nconf,ibc(itmp),        8d8s06
     $      ibc(itmp2),nalpha,nadet,multh,icall,lprint,bc,ibc)          11d9s22
       ibcoff=itmp
      else                                                              8d8s06
       nadet(1)=nconf                                                   8d8s06
      end if
      ioffd=0                                                           12d20s06
      do i=1,nsymb                                                      12d20s06
       ioffdeta(i)=ioffd                                                12d20s06
       ioffd=ioffd+nadet(i)                                             12d20s06
      end do                                                            12d20s06
      call dws_sync
      call second(time3)                                                4d30s18
      telap=time3-time2-tovr
      tsort=tsort+telap
      call cnt12(ibc(ialpha),norb,nconf,m1,m2,nadet,nsymb,multh,m1c)    8d9s14
      if(icall.eq.1.and.lprint)then                                     4d20s18
        ioff=0                                                          8d24s06
        do i=1,nsymb                                                      8d23s06
         ioff=ioff+i                                                    8d24s06
        end do                                                            8d23s06
      end if                                                            7d14s06
      call dws_sync
      m1ca=0                                                            8d24s06
      ioff=0                                                            8d24s06
      do i=1,nsymb                                                      8d24s06
       do j=1,i                                                         8d24s06
        m1ca=m1ca+m1c(j+ioff,1)                                         8d24s06
       end do                                                           8d24s06
       ioff=ioff+i                                                      8d29s06
      end do                                                            8d24s06
      do i=1,nsymb                                                      8d8s06
       m12sym(i,1,1)=m1(i)                                              8d8s06
       m12sym(i,2,1)=m2(i)                                              8d8s06
      end do                                                            8d8s06
      m1a=m1(1)                                                         8d8s06
      m2a=m2(1)                                                         8d8s06
      do i=2,nsymb                                                      8d8s06
       m1a=m1a+m1(i)                                                    8d8s06
       m2a=m2a+m2(i)                                                    8d8s06
      end do                                                            8d8s06
      nn=(nsymb*(nsymb+1))/2
      idatas=ibcoff                                                     8d29s06
      iarg1=1                                                           8d24s06
      iarg3=0                                                           8d8s06
      ioff=0                                                            8d24s06
      do isb=1,nsymb                                                    8d24s06
       idata1(isb)=ibcoff                                               1d11s23
       if(m12sym(isb,1,1).gt.0)then                                       8d24s06
        iarg2=m12sym(isb,1,1)                                             8d24s06
        il1(isb)=1                                                      8d12s14
        ih1(isb)=m12sym(isb,1,1)                                        8d12s14
        ibcoff=idata1(isb)+ih1(isb)+1-il1(isb)                          8d24s06
       end if                                                           8d24s06
       idata2(isb)=ibcoff                                               1d11s23
       if(m12sym(isb,2,1).gt.0)then                                       8d24s06
        iarg2=m12sym(isb,2,1)                                             8d24s06
        il2(isb)=1                                                      8d12s14
        ih2(isb)=m12sym(isb,2,1)                                        8d12s14
        ibcoff=idata2(isb)+ih2(isb)+1-il2(isb)                          8d24s06
       end if                                                           8d24s06
       do isk=1,isb                                                     8d24s06
        iad=isk+ioff                                                    8d24s06
        idatac(iad)=ibcoff                                              1d11s23
        if(m1c(iad,1).gt.0)then                                         8d24s06
         iarg2=m1c(iad,1)                                               8d24s06
         ilc(iad)=1                                                     8d12s14
         ihc(iad)=m1c(iad,1)                                            8d12s14
         ibcoff=idatac(iad)+ihc(iad)+1-ilc(iad)                         8d24s06
        end if                                                          8d24s06
       end do                                                           8d24s06
       ioff=ioff+isb                                                    8d29s06
      end do                                                            8d24s06
      call enough('cas1b.  3',bc,ibc)
      call second(time4)
      telap=time4-time3-tovr
      tcn12=tcn12+telap
      nn=(nsymb*(nsymb+1))/2
      iarg4=iarg3                                                       8d29s06
      iarg3=-2                                                          8d8s06
      ix=ibcoff
      ibcoff=ix+norb
      call enough('cas1b.  4',bc,ibc)
      call gen12(ibc(ialpha),norb,nconf,idata1,il1,ih1,                 8d24s06
     $     idata2,il2,ih2,idatac,ilc,ihc,ibc(ix),m1,m2,multh,nadet,     8d9s14
     $     nsymb,bc,ibc)                                                11d14s22
      call second(time5)
      telap=time5-time4-tovr                                            4d30s18
      tgn12=tgn12+telap
      ism=ibcoff                                                        8d29s06
      irelo=ism+norb                                                    8d29s06
      ibcoff=irelo+norb                                                 8d29s06
      call enough('cas1b.  5',bc,ibc)
      jsm=ism-1                                                         8d29s06
      jrelo=irelo-1                                                     8d29s06
      do i=1,nsymb
       do j=1,iacto(i)                                                  8d29s06
        ibc(jsm+j)=i                                                    8d29s06
        ibc(jrelo+j)=j                                                  8d29s06
       end do                                                           8d29s06
       jsm=jsm+iacto(i)                                                 8d29s06
       jrelo=jrelo+iacto(i)                                             8d29s06
      end do                                                            8d29s06
      do isb=1,nsymb                                                    4d26s18
       mst(1,isb,1)=0                                                   5d2s18
       if(m1(isb).gt.0)then                                             4d26s18
        itmpsrt=ibcoff                                                  4d26s18
        ibcoff=itmpsrt+m1(isb)                                          4d26s18
        call enough('cas1b.  6',bc,ibc)
        call srtsym(ibc(idata1(isb)),m1(isb),ibc(ism),nsymb,            4d26s18
     $       ibc(idata1(isb)),ibc(itmpsrt),mst(1,isb,1))                4d26s18
        ibcoff=itmpsrt                                                  4d26s18
       end if                                                           4d26s18
      end do                                                            4d30s18
      do iss=1,iad                                                      4d30s18
       nz(iss,1)=0                                                      5d2s18
       if(m1c(iss,1).gt.0)then                                          4d30s18
        itmpsrt=ibcoff
        ibcoff=itmpsrt+m1c(iss,1)                                        4d30s18
        call enough('cas1b.  7',bc,ibc)
        call srtsymc(ibc(idatac(iss)),m1c(iss,1),ibc(ism),nsymb,           4d30s18
     $      ibc(idatac(iss)),ibc(itmpsrt),mnz(1,1,iss,1),nz(iss,1))     4d30s18
        ibcoff=itmpsrt                                                  4d30s18
       end if                                                           4d30s18
      end do                                                            4d26s18
      ibcoff=ism                                                        4d26s18
      call second(time6)
      telap=time6-time5-tovr
      tsrt=tsrt+telap
      time1=time6
      itop=1
      ibot=1
      nbetax=nbeta
      ix=max(nbetax,norbx-nbetax)
      in=min(nbetax,norbx-nbetax)
      do i=ix+1,norbx
       itop=itop*i
      end do
      do i=1,in
       ibot=ibot*i
      end do
      nconf=itop/ibot
      idwstest=nconf                                                    5d30s06
      nconfb=nconf
      iborb=ibcoff
      needa=nconf*nbeta
      need8=needa/8
      if(need8*8.lt.needa)need8=need8+1
      need8b=need8                                                      8d8s06
      ibcoff=iborb+need8
      ibeta=ibcoff
      needa=nconf*norb
      need8=needa/8
      if(need8*8.lt.needa)need8=need8+1
      ibcoff=ibeta+need8
      call enough('cas1b.  8',bc,ibc)
      if(nbeta.eq.0)then                                                5d2s18
       nconf2=1                                                         5d2s18
      else                                                              5d2s18
      call sgen(ibc(ibeta),nconf,norb,norb,nbeta,nconf2,ibc(iborb))
      end if                                                            5d2s18
      call second(time2)
      telap=time2-time1-tovr
      tsgen=tsgen+telap
      if(nconf.ne.nconf2)then
       write(6,*)('beta electrons ')
       write(6,*)('error generating configurations '),nconf,nconf2
       stop
      end if
      if(nsymb.ne.1)then
       itmp=ibcoff
       itmp2=itmp+need8
       ibcoff=itmp2+need8b                                              8d8s06
       call enough('cas1b.  9',bc,ibc)
       call sortdet(ibc(ibeta),ibc(iborb),norb,nconf,ibc(itmp),         8d8s06
     $      ibc(itmp2),nbeta,nbdet,multh,icall,lprint,bc,ibc)           11d9s22
       ibcoff=itmp
      else                                                              8d8s06
       nbdet(1)=nconf                                                   8d8s06
      end if
      ioffd=0                                                           12d20s06
      do i=1,nsymb                                                      12d20s06
       ioffdetb(i)=ioffd                                                12d20s06
       ioffd=ioffd+nbdet(i)                                             12d20s06
      end do                                                            12d20s06
      call second(time3)                                                4d30s18
      telap=time3-time2-tovr
      tsort=tsort+telap
      call cnt12(ibc(ibeta),norb,nconf,m1,m2,nbdet,nsymb,multh,m1c(1,2))8d9s14
      if(icall.eq.1)then                                                7d14s06
        ioff=0                                                          8d24s06
        do i=1,nsymb                                                      8d23s06
         ioff=ioff+i                                                    8d24s06
        end do                                                            8d23s06
      end if                                                            7d14s06
      do i=1,nsymb                                                      8d8s06
       m12sym(i,1,2)=m1(i)                                              8d8s06
       m12sym(i,2,2)=m2(i)                                              8d8s06
      end do                                                            8d8s06
      m1b=m1(1)                                                         8d8s06
      m2b=m2(1)                                                         8d8s06
      do i=2,nsymb                                                      8d8s06
       m1b=m1b+m1(i)                                                    8d8s06
       m2b=m2b+m2(i)                                                    8d8s06
      end do                                                            8d8s06
      idatas=ibcoff                                                     8d29s06
      iarg1=1                                                           8d24s06
      iarg3=0                                                           8d8s06
      ioff=0                                                            8d24s06
      do isb=1,nsymb                                                    8d24s06
       idatb1(isb)=ibcoff                                               1d11s23
       if(m12sym(isb,1,2).gt.0)then                                     8d29s06
        iarg2=m12sym(isb,1,2)                                             8d24s06
        il1(isb)=1                                                      8d12s14
        ih1(isb)=m12sym(isb,1,2)                                        8d12s14
        ibcoff=idatb1(isb)+ih1(isb)+1-il1(isb)                          8d24s06
       end if                                                           8d24s06
       idatb2(isb)=ibcoff                                               1d11s23
       if(m12sym(isb,2,2).gt.0)then                                     8d29s06
        iarg2=m12sym(isb,2,2)                                           8d29s06
        il2(isb)=1                                                      8d12s14
        ih2(isb)=m12sym(isb,2,2)                                        8d12s14
        ibcoff=idatb2(isb)+ih2(isb)+1-il2(isb)                          8d24s06
       end if                                                           8d24s06
       do isk=1,isb                                                     8d24s06
        iad=isk+ioff                                                    8d24s06
        idatbc(iad)=ibcoff                                              1d11s23
        if(m1c(iad,2).gt.0)then                                         8d29s06
         iarg2=m1c(iad,2)                                               8d29s06
         ilc(iad)=1                                                     8d12s14
         ihc(iad)=m1c(iad,2)                                            8d12s14
         ibcoff=idatbc(iad)+ihc(iad)+1-ilc(iad)                         8d24s06
        end if                                                          8d24s06
       end do                                                           8d24s06
       ioff=ioff+isb                                                    8d29s06
      end do                                                            8d24s06
      call second(time4)
      telap=time4-time3-tovr
      tcn12=tcn12+telap
      call enough('cas1b. 10',bc,ibc)
      iarg4=iarg3                                                       8d29s06
      iarg3=-2                                                          8d8s06
      ix=ibcoff
      ibcoff=ix+norb
      call enough('cas1b. 11',bc,ibc)
      call gen12(ibc(ibeta),norb,nconf,idatb1,il1,ih1,                  8d29s06
     $     idatb2,il2,ih2,idatbc,ilc,ihc,ibc(ix),m1,m2,multh,nbdet,     8d9s14
     $     nsymb,bc,ibc)                                                11d14s22
      call second(time5)
      telap=time5-time4-tovr                                            4d30s18
      tgn12=tgn12+telap
      ism=ibcoff                                                        8d29s06
      irelo=ism+norb                                                    8d29s06
      ibcoff=irelo+norb                                                 8d29s06
      call enough('cas1b. 12',bc,ibc)
      jsm=ism-1                                                         8d29s06
      jrelo=irelo-1                                                     8d29s06
      do i=1,nsymb
       do j=1,iacto(i)                                                  8d29s06
        ibc(jsm+j)=i                                                    8d29s06
        ibc(jrelo+j)=j                                                  8d29s06
       end do                                                           8d29s06
       jsm=jsm+iacto(i)                                                 8d29s06
       jrelo=jrelo+iacto(i)                                             8d29s06
      end do                                                            8d29s06
      do isb=1,nsymb                                                    4d26s18
       mst(1,isb,2)=0                                                   5d2s18
       if(m1(isb).gt.0)then                                             4d26s18
        itmpsrt=ibcoff                                                  4d26s18
        ibcoff=itmpsrt+m1(isb)                                          4d26s18
        call enough('cas1b. 13',bc,ibc)
        call srtsym(ibc(idatb1(isb)),m1(isb),ibc(ism),nsymb,            4d26s18
     $       ibc(idatb1(isb)),ibc(itmpsrt),mst(1,isb,2))                4d26s18
        ibcoff=itmpsrt                                                  4d26s18
       end if                                                           4d26s18
      end do                                                            4d26s18
      do iss=1,iad                                                      4d30s18
       nz(iss,2)=0                                                      5d2s18
       if(m1c(iss,2).gt.0)then                                          4d30s18
        itmpsrt=ibcoff
        ibcoff=itmpsrt+m1c(iss,2)                                        4d30s18
        call enough('cas1b. 14',bc,ibc)
        call srtsymc(ibc(idatbc(iss)),m1c(iss,2),ibc(ism),nsymb,           4d30s18
     $      ibc(idatbc(iss)),ibc(itmpsrt),mnz(1,1,iss,2),nz(iss,2))     4d30s18
        ibcoff=itmpsrt                                                  4d30s18
       end if                                                           4d30s18
      end do                                                            4d26s18
      ibcoff=ism                                                        4d26s18
      call second(time6)
      telap=time6-time5-tovr
      tsrt=tsrt+telap
      return
      end
