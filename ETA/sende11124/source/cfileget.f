c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine cfileget(idorb,isorb,mdon,mdoo,ibasis,iptrbit,ncsf,    8d10s22
     $     nfcn,multh,ieiginth,ivintinth,isb,nrth,nlzz,lam,lz2,nameci,  8d10s22
     $     nbasisp,nbasdws,nvirt,nsymb,nct,norb,shift,lprint,idorel,bc, 11d14s22
     $     ibc,cfile)                                                   12d19s22
      implicit real*8 (a-h,o-z)                                         8d10s22
      parameter (nspc=22)
      character*(*) cfile                                               12d19s22
      integer*8 ipack8                                                  8d10s22
      integer*4 ipack4(2)                                               8d10s22
      equivalence (ipack8,ipack4)                                       8d10s22
      logical lprint                                                    8d11s22
      dimension nameci(*),multh(8,8),iwavedat(nspc),iorb(8),            8d10s22
     $     nbasisp(*),nbasdws(*),nvirt(*),nbasisp1(8),nbasdws1(8),      8d10s22
     $     nvirt1(8),iptrbit(2,mdoo+1,*),ibasis(3,*),ncsf(*)            8d10s22
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      include "common.store"                                            8d10s22
      if(lprint)then                                                    8d11s22
       write(6,*)('Hi, my name is cfileget ')                            8d10s22
       write(6,*)('input symmetry: '),isb
      end if                                                            8d11s22
      ieiginth=ibcoff                                                   8d10s22
      ivintinth=ieiginth+nrth                                           8d10s22
      ibcoff=ivintinth+nrth*nct                                         8d10s22
      ibcoffo=ibcoff                                                    8d11s22
      ioff=0                                                            8d10s22
      do if=1,nfcn                                                      8d10s22
       nclo=ibasis(1,if)                                                8d10s22
       nclop=nclo+1                                                     8d10s22
       iarg=nclop-mdon                                                  8d10s22
       iic=iptrbit(1,nclop,isb)+ibasis(2,if)-1                          8d10s22
       iio=iptrbit(2,nclop,isb)+ibasis(3,if)-1                          8d10s22
       ioff=ioff+ncsf(iarg)
      end do                                                            8d10s22
      noff=0                                                            8d10s22
      iout=ibcoff                                                       8d10s22
      itargs=10                                                         1d10s23
      itargd=10                                                         1d10s23
      elowest=0d0                                                       2d1s23
      call getwf(cfile,ismult,isymmrci,nroot,iout,iwavedat,mdoox,       12d19s22
     $     mxopm,iorb,0,noff,nspc,2,elowest,idorel,nsymb,               8d10s22
     $     nbasisp,nbasdws,nvirt,multh,maxbx,maxbxb,0,bc,ibc,mddilow,   12d13s22
     $     mddihig,itargs,itargd)                                       1d10s23
      if(lprint)write(6,*)('wf symmetry from file: '),iwavedat(2)       8d11s22
      if(iwavedat(2).ne.isb)then                                        8d10s22
       write(6,*)('but this doesn''t match input symmetry = '),isb
       call dws_synca
       call dws_finalize
       stop
      end if
      nff0=iwavedat(4)+iwavedat(9)                                           7d27s21
      jff0=iwavedat(4)+iwavedat(10)                                          7d27s21
      nec=iwavedat(7)                                                      7d27s21
      ilook=iwavedat(4)+iwavedat(13)                                          7d27s21
      ipack8=ibc(ilook)                                                 5d12s21
      ncsftk=ipack4(1)
      ilook=ilook+1
      ieg=ilook+ncsftk*iwavedat(3)                                      8d11s22
      do i=0,nrth-1                                                     8d11s22
       bc(ieiginth+i)=bc(ieg+i)-shift                                   8d11s22
      end do                                                            8d11s22
      call mapv2v(bc(ilook),ncsftk,ibc(nff0),ibc(jff0),bc(ivintinth),   8d10s22
     $     nct,nrth,nfcn,ibasis,iptrbit(1,1,isb),mdon,mdoo,ncsf,norb,   8d11s22
     $     nec,lprint,bc,ibc)                                           11d14s22
      ibcoff=ibcoffo                                                    8d11s22
      return                                                            8d10s22
      end                                                               8d10s22
