c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine gencsf1(nec,s0,norb,nos,nod,mdon,mdoo,ncsf,maxopens,   6d12s19
     $     ixsoo,iunit2,ndet,nwavef,nfill,ifill,idox)                   8d28s24
      implicit real*8 (a-h,o-z)
      real*16 fl
c
c     generate csf basis
c
      dimension ncsf(*),ndet(*),ifill(idox,*)                           8d28s24
      COMMON/FACT16/FL(922),NCALL                                       5d27s19
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      is02=nint(s0*2d0)
      nbeta=(nec-is02)/2
      nalpha=nec-nbeta
      xms=0.5d0*(nalpha-nbeta)                                          6d12s19
      if(mynowprog.eq.0)then
       write(6,*)('in gencsf for spin '),s0
       write(6,*)('number of electrons to correlate '),nec
       write(6,*)('number of orbitals '),norb
       write(6,*)('nalpha,nbeta: '),nalpha,nbeta
       write(6,*)                                                       8d28s24
     $      ('  #cs  #os      nae nbe      ncs       noc ndet ncsf')    8d28s24
      end if                                                            11d23s19
      nos=0
      nod=0
      mdon=nbeta
      icsf=1                                                            6d5s19
      ixsoo=-1                                                          6d12s19
      do idoo=0,nbeta
       nae=nalpha-idoo
       nbe=nbeta-idoo
       isoo=nae+nbe
c
c     the +2 allows for double excitations into virtual space           10d18s20
c     and +2 allows for interacting space                               10d18s20
c
       if(isoo+idoo.le.norb+iunit2.and.isoo.le.maxopens+4.and.          10d23s20
     $      min(nae,nbe).ge.0)then                                      10d23s20
        ixsoo=max(ixsoo,isoo)                                           6d12s19
        itop=norb+1
        ibot1=idoo+1
        ibot2=itop-ibot1+1
        ff=fl(itop)-fl(ibot1)-fl(ibot2)
        ndoo=nint(exp(ff))
        itop=norb+1                                                     7d11s19
        ibot1=isoo+1
        ibot2=iabs(itop-ibot1)+1                                        1d17s23
        ff=fl(itop)-fl(ibot1)-fl(ibot2)
        nsoo=nint(exp(ff))
        itop=isoo+1
        ibot1=nae+1
        ibot2=iabs(itop-ibot1)+1                                        1d17s23
        ff=fl(itop)-fl(ibot1)-fl(ibot2)
        ndet(icsf)=nint(exp(ff))                                        11d19s20
        nother=0
        if(ndet(icsf).gt.1)then                                         11d19s20
         do iae=nae+1,nae+1,-1
          ibot1=iae+1
          ibot2=itop-ibot1+1
          ff=fl(itop)-fl(ibot1)-fl(ibot2)
          nhere=nint(exp(ff))
          nother=nother+nhere
         end do
        end if
        ncsf(icsf)=ndet(icsf)-nother                                    11d19s20
        if(mynowprog.eq.0)                                              11d23s19
     $       write(6,2)idoo,isoo,nae,nbe,ndoo,nsoo,ndet(icsf),ncsf(icsf)11d19s20
    2   format(2i5,5x,2i4,5x,i5,i10,2i5)                                8d28s24
        icsf=icsf+1                                                     6d5s19
        nos=nos+nsoo*max(1,isoo)                                        7d5s19
        nod=nod+ndoo*idoo                                               6d12s19
        mdoo=idoo
        mdon=min(mdon,idoo)
       end if
      end do
      return
      end
