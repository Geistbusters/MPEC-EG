c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine savewf(iwopen,nwavrec,iunc,ethred,nctf,nroot,ntot,     8d16s22
     $     vint,eigint,shift,nlzzci,lambdaci,nsing,ihsdiag,nff1,iff1,   8d16s22
     $     ncsf,mdon,mdoo,nsymb,multh,nvirt,isymmrci,vs,ndoub,nff22,    8d16s22
     $     nfdat,ivdinout,mdoub,ncsf2,ihddiag,nff2,iff2,bc,ibc,nder)    10d7s24
      implicit real*8 (a-h,o-z)                                         8d16s22
      dimension iunc(2),vint(*),eigint(*)                               8d16s22
      integer*4 ipack4                                                  10d7s24
      integer*2 ipack2(2)                                               10d7s24
      equivalence (ipack2,ipack4)                                       10d7s24
      include "common.store"                                            8d16s22
      data icall/0/                                                     8d19s22
      save icall                                                        8d19s22
      icall=icall+1                                                     8d19s22
c
c     save restart data to wavef on unit2.
c     if iwopen=1, then wavef is open and we can simply write to i.
c     if iwopen=0, we need to rewind it, skip over basis set info,
c     and then write to it.
c
      if(iwopen.eq.0)then                                               8d16s22
       open(unit=2,file='wavef',form='unformatted',status='old')        8d16s22
       rewind(2)                                                        8d16s22
       do i=1,nwavrec                                                   8d16s22
        read(2)
       end do                                                           8d16s22
      end if                                                            8d16s22
      write(2)iunc,ethred
      write(2)nctf,ntot
      write(2)(vint(i),i=1,nctf*nroot)                                  8d16s22
      write(2)(eigint(i)+shift,i=1,nroot)                               8d16s22
      ipack2(1)=nlzzci                                                  10d7s24
      ipack2(2)=nder                                                    10d7s24
      write(2)ipack4,lambdaci                                           10d7s24
      if(nsing.ne.0)then                                                8d16s22
       call dumpsing(ihsdiag,nff1,iff1,ncsf,mdon,mdoo,nsymb,multh,nvirt,8d16s22
     $     isymmrci,nroot,vs,0,bc,ibc)                                  11d14s22
      else                                                              8d16s22
       izero=0                                                          8d16s22
       write(2)izero                                                    8d16s22
      end if                                                            8d16s22
      if(ndoub.gt.0)then                                                8d16s22
       if(iunc(2).eq.0)then                                             8d16s22
        call dumpdoubc(ibc(nff22),nfdat,bc(ivdinout),nsymb,ndoub,       8d16s22
     $         mdoub,nroot,mdon,mdoo,ncsf,ncsf2,0,bc,ibc,nder)          10d7s24
       else                                                             8d16s22
        call dumpdoubuc(ibc(ihddiag),ibc(nff2),ibc(iff2),ncsf,          8d16s22
     $      ncsf2,mdon,mdoo,nsymb,multh,nvirt,isymmrci,nroot,0,4,bc,ibc)1d10s23
       end if                                                           8d16s22
      else                                                              8d16s22
       izero=0                                                          8d16s22
       write(2)izero                                                    8d16s22
      end if                                                            8d16s22
      if(nsing.ne.0)then                                                8d16s22
       call dumpsing(ihsdiag,nff1,iff1,ncsf,mdon,mdoo,nsymb,multh,nvirt,8d16s22
     $     isymmrci,nroot,vs,1,bc,ibc)                                  11d14s22
      end if                                                            8d16s22
      if(ndoub.gt.0)then                                                8d16s22
       if(iunc(2).ne.0)then                                             8d16s22
        call dumpdoubuc(ibc(ihddiag),ibc(nff2),ibc(iff2),ncsf,          8d16s22
     $      ncsf2,mdon,mdoo,nsymb,multh,nvirt,isymmrci,nroot,1,4,bc,ibc)1d10s23
       end if                                                           8d16s22
      end if                                                            8d16s22
      close(unit=2)                                                     8d16s22
      iwopen=0                                                          8d16s22
      return                                                            8d16s22
      end                                                               8d16s22
