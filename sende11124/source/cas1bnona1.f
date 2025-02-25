c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine cas1bnona1(ipuse,nsym,spinm,nelec,nalpha,nbeta,nproc,  6d13s22
     $    nowpro,ialpha,iaorb,ibeta,iborb,nconfa,nconfb,m12sym,multh,   6d13s22
     $     idata1,idata2,idatb1,idatb2,lprint,bc,ibc)                   11d9s22
      implicit real*8 (a-h,o-z)
      external second
      integer spinm
      logical lprint                                                    4d20s18
      dimension m12sym(8,2,2),idata1(8),idata2(8),idatb1(8),idatb2(8)   6d13s22
c
c     storage of integers.
c     sgen: maximum number of orbitals             integer*1=128
c     number of alpha or beta configs in sym block integer*2=32768
c
c     generate indicies for computing fci hamiltonian
c     this version is for bra and ket having different symmetries.      6d10s22
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
       write(6,*)('spin is incompatible with number of electrons ')
       stop
      end if
      nos=int(spin*2d0+0.01d0)                                          5d2s07
      if(nbeta.lt.0)then
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
      call cnt12nona1(ipuse,ibc(ialpha),norb,nconfa,m12sym(1,1,1),      6d13s22
     $     m12sym(1,2,1),nadet,nsymb,multh)                             6d13s22
      do isb=1,nsymb                                                    6d13s22
       idata1(isb)=ibcoff                                               6d13s22
       idata2(isb)=idata1(isb)+m12sym(isb,1,1)                          6d13s22
       ibcoff=idata2(isb)+m12sym(isb,2,1)
      end do                                                            6d13s22
      ix=ibcoff
      ibcoff=ix+norb
      call enough('cas1bnona1.  1',bc,ibc)
      call gen12nona1(ipuse,ibc(ialpha),norb,nconfa,idata1,             6d13s22
     $     m12sym(1,1,1),idata2,m12sym(1,2,1),ibc(ix),multh,nadet,nsymb,11d14s22
     $     bc,ibc)                                                      11d14s22
      ibcoff=ix                                                         6d13s22
      call cnt12nona1(ipuse,ibc(ibeta),norb,nconfb,m12sym(1,1,2),       6d13s22
     $     m12sym(1,2,2),nbdet,nsymb,multh)                             6d13s22
      do isb=1,nsymb                                                    6d13s22
       idatb1(isb)=ibcoff                                               6d13s22
       idatb2(isb)=idatb1(isb)+m12sym(isb,1,2)                          6d13s22
       ibcoff=idatb2(isb)+m12sym(isb,2,2)                               6d13s22
      end do                                                            6d13s22
      ix=ibcoff                                                         6d13s22
      ibcoff=ix+norb                                                    6d13s22
      call enough('cas1bnona1.  2',bc,ibc)
      call gen12nona1(ipuse,ibc(ibeta),norb,nconfb,idatb1,m12sym(1,1,2),6d13s22
     $     idatb2,m12sym(1,2,2),ibc(ix),multh,nbdet,nsymb,bc,ibc)       11d14s22
      ibcoff=ix                                                         6d13s22
      return
      end
