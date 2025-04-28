c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine getwf(fname,ismult,isymmrci,nroot,iout,iwavedat,       2d11s22
     $     mdoox,mxopn,iorb,icall,noff,nspc,mwavedat,elowest,idorel,    7d8s21
     $     nsymb,nbasp,nbasdws,nvirt,multh,maxbx,maxbxd,inorm,bc,ibc,   12d13s22
     $     mddilow,mddihig,itargs,itargd)                               12d21s22
      implicit real*8 (a-h,o-z)
      character*(*) fname
      parameter (id=100)                                                11d25s19
      integer*8 ipack8                                                  11d23s19
      integer*4 ipack4(2),jpack4                                        10d7s24
      integer*1 ipack1(4)                                               5d12s21
      integer*2 jpack2(2)
      equivalence (ipack8,ipack4)                                       11d23s19
      equivalence (npack4,ipack1)                                       5d12s21
      equivalence (jpack4,jpack2)                                       10d7s24
      include "common.store"
      include "common.basis"
      dimension multhx(64),nbasdws(8),nbasc(8),nbasp(8),isymx(3,8),
     $     idata(id),iunc(2),iwavedat(nspc,*),iorb(8),nameci(6),        5d21s21
     $     idumx18(18),nvirt(*),multh(8,8),isc2v(4),isd2h(8),isxa(4)    12d6s22
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      data isc2v/4,3,2,1/                                               6d18s22
      data isd2h/4,3,2,1,8,7,6,5/                                       6d18s22
      if(idorel.eq.0)then                                               1d11s20
       ncomp=1                                                          1d11s20
      else                                                              1d11s20
       ncomp=2                                                          1d11s20
      end if                                                            1d11s20
      if(mynowprog.eq.0)then                                            2d3s22
       write(6,*)('in getwf for file '),fname
      end if                                                            2d3s22
      ibcoffo=ibcoff                                                    5d7s21
      if(mynowprog.eq.0)then                                            7d8s21
       if(inorm.ne.0)open(unit=1,file=fname,form='unformatted')         8d10s22
       nwavrec=1                                                         5d4s21
       read(1)ismult,isymmrci,nroot,norb,nsb,nameci                      12d27s19
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)idum
       if(norb.gt.0)then                                                 5d4s21
        nwavrec=nwavrec+1                                                 5d4s21
        read(1)idum                                                      5d4s21
       end if                                                            5d4s21
       if(mynowprog.eq.0)then                                           2d3s22
        write(6,*)('spin multiplicity '),ismult
        write(6,*)('symmetry of wavefunction '),isymmrci
        write(6,*)('number of roots = '),nroot
       end if                                                           2d3s22
       iwavedat(1,1)=ismult                                              5d12s21
       iwavedat(2,1)=isymmrci                                            5d12s21
       iwavedat(3,1)=nroot                                               5d12s21
       iwavedat(4,1)=iout                                                5d12s21
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)nsymbi,idorelx,ngaus,natom,nwcont,numeminus,lmax,nbasallp,8d10s22
     $     multhx,ascale,potdws,idumx18                                 5d25s21
       nextradatatag=idumx18(18)+1                                      2d14s22
       nvec=natom-1                                                      5d26s21
       nexpdata=0                                                        5d26s21
       if(nvec.eq.1)then                                                 5d26s21
        nexpdata=1                                                       5d26s21
       else if(nvec.gt.1)then                                            5d26s21
        nexpdata=3*nvec-3                                                5d26s21
       end if                                                            5d26s21
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)idum                                                      8d10s22
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)(nbasc(isb),isb=1,nsymb)                                   11d20s19
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)idum                                                      8d10s22
       nbasallp=nbasp(1)                                                 12d24s19
       do isb=2,nsymb                                                    12d24s19
        nbasallp=nbasallp+nbasp(isb)                                     12d24s19
       end do                                                            12d24s19
       itmp1=ibcoff
       itmp2=itmp1+natom*3
       iextrad=itmp2+natom*3
       ibcoff=iextrad+nextradatatag                                      5d27s21
       call enough('getwf.  1',bc,ibc)
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)isymx,(bc(itmp1+i),bc(itmp2+i),i=0,natom*3-1),             5d27s21
     $     (bc(iextrad+i),i=0,nextradatatag-1)                          5d27s21
       do i=1,3                                                         1d13s23
        cmx(i,1)=0d0                                                    1d13s23
       end do                                                           1d13s23
       do i=1,natom                                                     1d13s23
        iagrp(i)=1                                                      1d13s23
       end do                                                           1d13s23
       if(nextradatatag-nexpdata.ge.6+natom.and.icall.eq.1)then         5d4s22
        jextrad=iextrad+nexpdata                                         5d27s21
        do i=1,2                                                         5d27s21
         do j=1,3                                                        5d27s21
          cmx(j,i)=bc(jextrad)                                           5d27s21
          jextrad=jextrad+1                                              5d27s21
         end do                                                          5d27s21
        end do                                                           5d27s21
        do i=1,natom                                                     5d27s21
         iagrp(i)=ibc(jextrad)                                           5d27s21
         jextrad=jextrad+1                                               5d27s21
        end do                                                           5d27s21
       end if                                                            5d26s21
       ibcoff=iextrad                                                    5d27s21
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)dum
       rms=0d0                                                           11d22s19
       jtmp1=itmp1-1
       do i=1,natom                                                     1d19s23
        do ixyz=1,3                                                     1d19s23
         diff=xcart(ixyz,i)-bc(jtmp1+ixyz)                              1d19s23
         rms=rms+diff**2                                                1d19s23
        end do                                                          1d19s23
        jtmp1=jtmp1+3                                                   1d19s23
       end do                                                           1d19s23
       rms=sqrt(rms/dfloat(natom*3))                                     11d22s19
       if(rms.lt.1d-8)then                                               11d24s19
        igeod=0                                                          11d24s19
       else                                                              11d24s19
        igeod=1                                                          11d24s19
       end if                                                            11d24s19
       if(inorm.eq.0)write(6,*)('rms geometry difference: '),rms        8d10s22
       ibcoff=itmp1
       nbdat=ngaus*9+nwcont
       ibcoff=itmp1+nbdat
       call enough('getwf.  2',bc,ibc)
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)(bc(itmp1+i),i=0,nbdat-1)
       if(idorel.ne.0)nbasallp=nbasallp/2                                1d13s20
       ibcoff=itmp1+nbasallp*6
       call enough('getwf.  3',bc,ibc)
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)(ibc(itmp1+i),i=0,nbasallp*6-1)
       jout=iout                                                         11d22s19
       do isb=1,nsymb
        need=nbasp(isb)*nbasdws(isb)*ncomp                               1d11s20
        ibcoff=jout+need                                                 11d22s19
        call enough('getwf.  4',bc,ibc)
        if(need.gt.0)then
         nwavrec=nwavrec+1                                                 5d4s21
         read(1)(bc(jout+i),i=0,need-1)                                  5d3s21
        end if                                                           5d3s21
        if(icall.eq.1)then                                               12d16s19
         do i=0,need-1                                                   1d11s20
          bc(iorb(isb)+i)=bc(jout+i)                                     12d16s19
         end do                                                          12d16s19
        end if                                                           12d16s19
        jout=ibcoff                                                      11d22s19
        need=nbasdws(isb)*nbasdws(isb)                                   11d25s19
        ibcoff=jout+need                                                 11d25s19
        if(need.gt.0)then
         nwavrec=nwavrec+1                                                 5d4s21
         read(1)(bc(jout+i),i=0,need-1)                                  5d3s21
        end if                                                           5d3s21
        jout=ibcoff                                                      11d25s19
       end do
       nwavrec=nwavrec+1                                                 5d4s21
       read(1)nec,jff0,mdon,mdoop,nct                                    5d12s21
       iwavedat(21,1)=mdoop                                             8d31s21
       nopenx=nec-2*mdon                                                 5d7s21
       mxopn=max(mxopn,nopenx)                                           5d7s21
       mdoo=mdoop-1                                                      5d3s21
       ipack4(1)=mdoo                                                    11d23s19
       ipack4(2)=mdon                                                    11d23s19
       mdoox=max(mdoo,mdoox)                                             11d25s19
       mdoop=mdoo+1                                                      5d3s21
       ibc(jout)=ipack8                                                  11d23s19
       jout=jout+1                                                       11d23s19
       nhere=0                                                           5d12s21
       if(nhere.gt.id)stop 'id,getwf'
       ipack4(1)=nhere                                                   11d24s19
       ipack4(2)=nec                                                     11d24s19
       ibc(jout)=ipack8                                                  11d24s19
       iwavedat(7,1)=nec                                                 5d12s21
       jout=jout+1                                                       11d24s19
       ipack4(1)=jff0                                                    4d30s21
       ipack4(2)=mdoop                                                   11d23s19
       ibc(jout)=ipack8                                                  11d23s19
       jout=jout+1                                                       11d23s19
       iwavedat(8,1)=jff0                                                5d12s21
       nff0=jout                                                         4d30s21
       iwavedat(9,1)=nff0-iout                                           5d12s21
       iff0=nff0+3*mdoop                                                 5d7s21
       iwavedat(10,1)=iff0-iout                                          5d12s21
       jout=iff0+jff0                                                    4d30s21
       call getbas0(ibc(nff0),mdoop,ibc(iff0),jff0)                       4d30s21
       read(1)iunc,ethred                                                11d23s19
       read(1)nctf,ntot                                                  11d24s19
       iwavedat(13,1)=jout-iout                                          5d12s21
       do i=1,6                                                          12d27s19
        iwavedat(13+i,1)=nameci(i)                                       5d12s21
       end do                                                            12d27s19
       nwavedat=19                                                       12d27s19
       ipack4(1)=nctf                                                    11d24s19
       ipack4(2)=ntot                                                    11d24s19
       ibc(jout)=ipack8                                                  11d24s19
       jout=jout+1                                                       11d24s19
       ivecdws=jout
       read(1)(bc(jout+i),i=0,nctf*nroot-1)                              11d24s19
       jout=jout+nctf*nroot                                              11d24s19
       ibcoff=jout                                                       11d24s19
       ieig=ibcoff                                                       3d20s20
       ibcoff=ieig+nroot                                                 3d20s20
       call enough('getwf.  5',bc,ibc)
       read(1)(bc(ieig+i),i=0,nroot-1)                                   3d20s20
       if(mynowprog.eq.0)then                                           2d15s22
        if(nroot.eq.1)then                                              7d11s23
         write(6,*)('eigenvalue '),bc(ieig)                             7d11s23
        else                                                            7d11s23
         write(6,*)('eigenvalues '),(bc(ieig+i),i=0,nroot-1)               3d20s20
        end if                                                          7d11s23
       end if                                                           2d15s22
       jout=jout+nroot                                                  7d8s21
       elowest=min(elowest,bc(ieig))                                     6d1s21
       ipack1(1)=igeod                                                   5d12s21
       read(1)jpack4,lambdaci                                           10d7s24
       nlzz=jpack2(1)                                                   10d7s24
       nder=jpack2(2)                                                   10d7s24
       ipack1(2)=nlzz                                                    5d12s21
       ipack1(3)=lambdaci                                                5d12s21
       ipack1(4)=lambdaci                                                5d24s21
       if(nlzz.eq.0)then                                                3d28s22
        write(6,*)('general molecule')                                  3d28s22
        iord=0                                                          5d27s21
        if(nsymb.eq.2)then                                               3d28s22
         if(iwavedat(2,1).eq.2)iord=1                                    3d28s22
        else if(nsymb.eq.4)then                                          3d28s22
         if(iwavedat(2,1).gt.2)iord=1                                   3d28s22
        else if(nsymb.eq.8)then                                          3d28s22
         if(iwavedat(2,1).eq.3.or.iwavedat(2,1).eq.4.or.                3d28s22
     $        iwavedat(2,1).eq.7.or.iwavedat(2,1).eq.8)iord=1           3d28s22
        end if                                                           3d28s22
        if(iord.eq.0)then                                               3d28s22
         write(6,*)('this wave function will be taken as pure real')    3d28s22
        else                                                            3d28s22
         write(6,*)('this wave function will be taken as pure imag')    3d28s22
        end if                                                          3d28s22
        ipack1(3)=iord                                                  3d28s22
       end if                                                           3d28s22
       iwavedat(6,1)=npack4                                              5d12s21
       read(1)mff1                                                      7d8s21
       if(mff1.eq.0)then                                                12d21s22
        jtargs=0                                                        12d21s22
       else                                                             12d21s22
        jtargs=1                                                        12d21s22
       end if                                                           12d21s22
       if(itargs.eq.10)then                                             12d21s22
        itargs=jtargs                                                   12d21s22
       else                                                             12d21s22
        if(itargs.ne.jtargs)then                                        12d21s22
         write(6,*)('one wavefunction file has singles turned on! ')     12d21s22
         write(6,*)                                                     12d21s22
     $        ('and one wavefunction file has singles turned off! ')    12d21s22
         write(6,*)                                                     12d21s22
     $       ('wavefuncitons must all have same levels of excitations!')12d21s22
         stop 'getwf'                                                   12d21s22
        end if                                                          12d21s22
       end if                                                           12d21s22
       ibc(jout)=mff1                                                   7d8s21
       jout=jout+1                                                      7d8s21
       if(mff1.gt.0)then                                                7d8s21
        ihsdiag=jout                                                    7d8s21
        nff1=ihsdiag+nsymb*mdoop*2                                      7d9s21
        iff1=nff1+nsymb*mdoop                                           7d8s21
        icsf=iff1+mff1                                                  7d8s21
        jout=icsf+mdoop-mdon                                            7d8s21
        call loadsinga(ibc(nff1),ibc(iff1),mff1,nsymb,mdon,mdoo,        7d8s21
     $       ibc(icsf),1)                                               8d18s22
       end if                                                           7d8s21
       read(1)mff2                                                      7d21s21
       ibc(jout)=mff2                                                   7d21s21
       if(mff2.eq.0)then                                                12d21s22
        jtargd=0                                                        12d21s22
       else if(mff2.gt.0)then                                           12d21s22
        jtargd=1                                                        12d21s22
       else                                                             12d21s22
        jtargd=-1                                                        12d21s22
       end if                                                           12d21s22
       if(itargd.eq.10)then                                             12d21s22
        itargd=jtargd                                                   12d21s22
       else                                                             12d21s22
        if(itargd.ne.jtargd)then                                        12d21s22
         do ipass=1,2                                                    12d21s22
          if(ipass.eq.1)then                                            12d21s22
           ktargd=itargd                                                12d21s22
          else                                                          12d21s22
           ktargd=jtargd                                                12d21s22
          end if                                                        12d21s22
          if(ktargd.gt.0)then
           write(6,*)                                                    12d21s22
     $    ('wavefunction file has uncontracted doubles turned on! ')    12d21s22
          else if(ktargd.lt.0)then                                       12d21s22
           write(6,*)                                                    12d21s22
     $    ('wavefunction file has contracted doubles turned on! ')      12d21s22
          else                                                           12d21s22
           write(6,*)                                                     12d21s22
     $        ('wavefunction file has double turned off! ')             12d21s22
          end if                                                         12d21s22
         end do                                                         12d21s22
         write(6,*)                                                     12d21s22
     $       ('wavefuncitons must all have same levels of excitations!')12d21s22
         stop 'getwf'                                                   12d21s22
        end if                                                          12d21s22
       end if                                                           12d21s22
       jout=jout+1                                                      7d21s21
       if(mff2.gt.0)then                                                7d8s21
        ihddiag=jout                                                    7d8s21
        nff2=ihddiag+nsymb*mdoop*2                                      7d9s21
        iff2=nff2+nsymb*mdoop                                           7d8s21
        icsf=iff2+mff2                                                  7d8s21
        icsf2=icsf+mdoop-mdon                                            7d8s21
        jout=icsf2+mdoop-mdon                                           7d21s21
        call loaddouba(ibc(nff2),ibc(iff2),mff2,nsymb,mdon,mdoo,        7d8s21
     $       ibc(icsf),ibc(icsf2),inorm)                                8d10s22
       else if(mff2.lt.0)then                                           8d3s21
        read(1)ndoub,mdoub                                              8d3s21
        write(6,*)('what we have for mff2 '),mff2
        mff2a=-mff2                                                     8d3s21
        ipack4(1)=ndoub                                                 8d3s21
        ipack4(2)=mdoub                                                 8d3s21
        mdoubstore=jout                                                 8d3s21
        ibc(jout)=ipack8                                                8d3s21
        iff22=jout+1                                                    8d3s21
        nff22=iff22+mff2a                                               8d3s21
        nfdat=nff22+mdoop*nsymb                                         8d3s21
        ivd=nfdat+10*nsymb                                              8d3s21
        icsf=ivd+nroot*(mdoub+ndoub)                                    8d12s21
        icsf2=icsf+mdoop-mdon                                           8d3s21
        jout=icsf2+2*(mdoop-mdon)                                       8d3s21
        nall=ndoub*nroot                                                8d3s21
        call loaddoubc(ibc(iff22),ibc(nff22),ibc(nfdat),                8d3s21
     $       ibc(icsf),ibc(icsf2),mdon,mdoo,nsymb,mff2a,bc(ivd),nall,   8d3s21
     $       isymmrci,nroot,multh,nvirt,inorm)                          8d10s22
        ndoubo=ndoub                                                    11d17s22
        mdoubo=mdoub                                                    11d17s22
        if(nlzz.ne.0.and.lambdaci.gt.0)then                             12d6s22
c
c     need sure to allocate enough space for ndoub,mdoub for other      6d18s22
c     lz component. In spherical case, symmetry with lz=0 should already6d18s22
c     be largest.
c
         if(nsymb.eq.4)then                                             6d18s22
          nsymo=isc2v(isymmrci)                                         6d18s22
         else                                                           6d18s22
          nsymo=isd2h(isymmrci)                                         6d18s22
         end if                                                         6d18s22
         if(nlzz.eq.6)then                                              12d6s22
          if(isymmrci.eq.1.or.isymmrci.eq.4.or.isymmrci.eq.6.or.
     $        isymmrci.eq.7)then                                        12d6s22
           isxa(1)=1                                                    12d6s22
           isxa(2)=4                                                    12d6s22
           isxa(3)=6                                                    12d6s22
           isxa(4)=7                                                    12d6s22
          else                                                          12d6s22
           isxa(1)=2                                                    12d6s22
           isxa(2)=3                                                    12d6s22
           isxa(3)=5                                                    12d6s22
           isxa(4)=8                                                    12d6s22
          end if                                                        12d6s22
          mdxx=0                                                        12d6s22
          ndxx=0                                                        12d6s22
          do isxx=1,4                                                   12d6s22
           isx=isxa(isxx)                                               12d6s22
           call countcc(isx,ibc(nfdat),nsymb,nvirt,mdoubo,ndoubo,        12d6s22
     $        multh)                                                    6d18s22
           mdxx=max(mdxx,mdoubo)                                        12d6s22
           ndxx=max(ndxx,ndoubo)                                        12d6s22
          end do                                                        12d6s22
          mdoubo=mdxx                                                   12d6s22
          ndoubo=ndxx                                                   12d6s22
         else                                                           12d6s22
          call countcc(nsymo,ibc(nfdat),nsymb,nvirt,mdoubo,ndoubo,multh) 6d18s22
         end if
        end if                                                          6d18s22
        ipack4(1)=ndoubo                                                11d17s22
        ipack4(2)=mdoubo                                                11d17s22
        ibc(jout)=ipack8                                                11d17s22
        othermdoub=bc(jout)                                             11d17s22
        ibcoff=jout                                                     8d12s21
       end if                                                           7d8s21
       ibcoff=jout                                                      7d8s21
       iwavedat(5,1)=ibcoff-iout                                         5d12s21
      end if                                                            7d8s21
      nspch=nspc/2                                                      7d8s21
      if(nspch*2.ne.nspc)nspch=nspch+1                                  7d8s21
      call dws_bcast(iwavedat,nspch)                                    7d8s21
      mdoops=iwavedat(21,1)                                             8d31s21
      nsnd=iwavedat(5,1)                                                7d8s21
      nec=iwavedat(7,1)                                                 7d9s21
      ismult=iwavedat(1,1)                                              7d9s21
      isymmrci=iwavedat(2,1)                                            2d3s22
      ibcoff=iout+nsnd                                                  7d8s21
      bc(ibcoff)=dfloat(mxopn)                                          7d9s21
      bc(ibcoff+1)=othermdoub                                           6d18s22
      call dws_bcast(bc(iout),nsnd+2)                                   6d18s22
      mxopn=nint(bc(ibcoff))                                            7d9s21
      ipack8=ibc(ibcoff+1)                                              6d18s22
      ndoubo=ipack4(1)                                                  6d18s22
      mdoubo=ipack4(2)                                                  6d18s22
      ipack8=ibc(iout+iwavedat(13,1))                                   7d8s21
      nctf=ipack4(1)
      nroot=iwavedat(3,1)                                               7d8s21
      npack4=iwavedat(6,1)                                              7d8s21
      nlzz=ipack1(2)                                                    7d8s21
      lambdaci=ipack1(3)
      iwavedat(4,1)=iout                                                7d8s21
      ieg=iwavedat(13,1)+nroot*nctf+iout+1                              7d9s21
      ieig=ieg                                                          7d9s21
      imff1=iwavedat(13,1)+nroot*(nctf+1)+1
      mff1=ibc(imff1+iout)                                              7d8s21
      jout=iout                                                         7d8s21
      do isb=1,nsymb
       need=nbasp(isb)*nbasdws(isb)*ncomp                               7d8s21
       jout=jout+need                                                   7d8s21
       need=nbasdws(isb)*nbasdws(isb)                                   11d25s19
       jout=jout+need                                                   7d8s21
      end do
      ipack8=ibc(jout)                                                  7d8s21
      mdoo=ipack4(1)                                                    7d8s21
      mdon=ipack4(2)                                                    7d8s21
      mdoox=max(mdoo,mdoox)                                             11d25s19
      mdoop=mdoo+1                                                      5d3s21
      inext=imff1+1+iout                                                7d21s21
      if(mff1.gt.0)then                                                 7d8s21
       ihsdiag=imff1+1+iout
       inff1=ihsdiag+nsymb*mdoop*2                                      7d9s21
       iiff1=inff1+nsymb*mdoop
       icsf=iiff1+ibc(imff1+iout)                                        7d8s21
       inext=icsf+mdoop-mdon                                            7d21s21
       call loadsingb(ibc(ihsdiag),ibc(inff1),nsymb,mdon,mdoo,           7d8s21
     $     nvirt,multh,iwavedat(2,1),ibc(icsf),nroot,mynowprog,maxbx,   8d10s22
     $      inorm,bc,ibc,mddilow,mddihig)                               12d13s22
      end if                                                            7d8s21
      mff2=ibc(inext)                                                   8d3s21
      inext=inext+1                                                     7d21s21
      if(mff2.gt.0)then                                                 7d21s21
       ihddiag=inext                                                    7d21s21
       nff2=ihddiag+nsymb*mdoop*2                                       7d21s21
       iff2=nff2+nsymb*mdoop                                            7d21s21
       icsf=iff2+mff2                                                   7d21s21
       icsf2=icsf+mdoop-mdon                                            7d8s21
       call loaddoubb(ibc(ihddiag),ibc(nff2),nsymb,mdon,mdoo,           7d8s21
     $     nvirt,multh,iwavedat(2,1),ibc(icsf),ibc(icsf2),nroot,        7d21s21
     $      mynowprog,maxbxd,inorm,bc,ibc,mddilow,mddihig)              12d13s22
      else if(mff2.lt.0)then                                            8d18s21
       ipack8=ibc(inext)                                                8d27s21
       mdoubstore=inext                                                 8d27s21
       ndoub=ipack4(1)                                                  8d12s21
       mdoub=ipack4(2)                                                  8d3s21
       mff2a=-mff2                                                      8d27s21
       iff22=inext+1                                                    8d27s21
       nff22=iff22+mff2a                                                8d27s21
       nfdat=nff22+mdoop*nsymb                                          8d27s21
       ivd=nfdat+10*nsymb                                               8d27s21
       icsf=ivd+nroot*(mdoub+ndoub)                                     8d27s21
       icsf2=icsf+mdoop-mdon                                            8d27s21
       ivdnon=ivd+nroot*ndoub                                           8d27s21
       if(inorm.ne.0)then                                               8d10s22
        call tofrob(bc(ivd),bc(ivdnon),nroot,ibc(nfdat),nvirt,nsymb,     8d27s21
     $       multh,isymmrci,1,ndoub,mdoub,bc(iff22),bc,ibc)             11d10s22
       end if                                                           8d10s22
      end if                                                            7d21s21
      npass=0                                                           5d12s21
      if(inorm.ne.0)then                                                8d10s22
       if(nlzz.eq.2)then                                                 5d12s21
        npass=1                                                          5d12s21
        if(mynowprog.eq.0)write(6,*)('Lz quantum number '),lambdaci      2d15s22
        if(lambdaci.eq.0)npass=0                                         8d30s21
       else if(nlzz.eq.6)then                                            5d12s21
        npass=2*lambdaci                                                 5d12s21
        if(mynowprog.eq.0)write(6,*)('L quantum number '),lambdaci       2d15s22
       end if                                                            5d12s21
      end if                                                            8d10s22
      if(npass.gt.0)then                                                5d12s21
       do ipass=1,npass                                                 5d12s21
        iout=ibcoff                                                      5d12s21
        jout=iout                                                        5d12s21
        ipp=ipass+1                                                     5d12s21
        if(noff+icall+ipass.gt.mwavedat)then                            5d12s21
         write(6,*)('we are exceeding mwavedat!!! '),noff,icall,ipass,  5d12s21
     $       mwavedat                                                   5d12s21
         stop 'prop'                                                    5d12s21
        end if                                                          5d12s21
        mysym=0                                                         1d13s23
        if(mynowprog.eq.0)then                                          7d8s21
         read(1)mysym                                                    5d12s21
        end if                                                          7d8s21
        ysym=dfloat(mysym)                                              7d8s21
        call dws_bcast(ysym,1)                                          7d8s21
        mysym=nint(ysym)                                                7d8s21
        iwavedat(1,ipp)=ismult                                          5d12s21
        iwavedat(2,ipp)=mysym                                           5d12s21
        iwavedat(3,ipp)=nroot                                           5d12s21
        iwavedat(4,ipp)=iout                                            5d12s21
        iwavedat(21,ipp)=mdoops                                         8d31s21
        ipack1(4)=-1                                                    5d24s21
        iwavedat(6,ipp)=npack4                                          5d12s21
        iwavedat(7,ipp)=nec                                             5d12s21
        jff0z=0                                                         1d13s23
        nctz=0                                                          1d13s23
        if(mynowprog.eq.0)then                                          7d8s21
         read(1)necz,jff0z,mdonz,mdoopz,nctz                             5d12s21
        end if                                                          7d8s21
        bc(iout)=dfloat(jff0z)                                          7d8s21
        bc(iout+1)=dfloat(nctz)                                         7d8s21
        call dws_bcast(bc(iout),2)                                      7d8s21
        jff0z=nint(bc(iout))                                            7d8s21
        nctz=nint(bc(iout+1))                                           7d8s21
        nff0z=jout                                                      5d12s21
        iff0z=nff0z+3*mdoop                                             5d12s21
        jout=iff0z+jff0z                                                5d12s21
        iwavedat(8,ipp)=jff0z                                           5d12s21
        iwavedat(9,ipp)=nff0z-iout
        iwavedat(10,ipp)=iff0z-iout                                     5d12s21
        if(mynowprog.eq.0)then                                          7d8s21
         call getbas0(ibc(nff0z),mdoop,ibc(iff0z),jff0z)                 5d12s21
        end if                                                          7d8s21
        nsnd=jout+1-nff0z                                               7d8s21
        call dws_bcast(bc(nff0z),nsnd)                                  7d8s21
        iwavedat(13,ipp)=jout-iout                                      5d12s21
        do i=1,6                                                        5d12s21
         iwavedat(13+i,ipp)=nameci(i)                                   5d12s21
        end do                                                          5d12s21
        ipack4(1)=nctz                                                  5d12s21
        ibc(jout)=ipack8                                                5d12s21
        jout=jout+1                                                     5d12s21
        do i=0,nctz*nroot-1                                             8d23s21
         bc(jout+i)=0d0                                                 8d23s21
        end do                                                          8d23s21
        jout=jout+nctz*nroot                                            5d12s21
        do i=0,nroot-1                                                  5d12s21
         bc(jout+i)=bc(ieig+i)                                          5d12s21
        end do                                                          5d12s21
        jout=jout+nroot                                                 7d8s21
        ibc(jout)=mff1                                                  7d8s21
        jout=jout+1                                                     7d8s21
        if(mff1.gt.0)then                                               7d8s21
         ihsdiagn=jout                                                    7d8s21
         nff1n=ihsdiagn+nsymb*mdoop*2                                   7d9s21
         iff1n=nff1n+nsymb*mdoop                                           7d8s21
         icsfn=iff1n+mff1                                                  7d8s21
         jout=icsfn+mdoop-mdon                                            7d8s21
         ncpy=jout-nff1n                                                7d8s21
         do i=0,ncpy-1                                                  7d8s21
          ibc(nff1n+i)=ibc(inff1+i)                                     8d19s21
         end do                                                         7d8s21
         ibcoff=jout                                                    11d9s22
         call loadsingb(ibc(ihsdiagn),ibc(inff1),nsymb,mdon,mdoo,       7d8s21
     $     nvirt,multh,mysym,ibc(icsfn),nroot,-1,maxbx,inorm,bc,ibc,    12d13s22
     $        mddilow,mddihig)                                          12d13s22
        end if                                                          7d8s21
        ibc(jout)=mff2                                                  7d21s21
        jout=jout+1                                                     8d18s21
        if(mff2.gt.0)then                                               7d21s21
         ihddiag=jout                                                   7d21s21
         nff2n=ihddiag+nsymb*mdoop*2                                       7d21s21
         iff2=nff2n+nsymb*mdoop                                            7d21s21
         icsf=iff2+mff2                                                   7d21s21
         icsf2=icsf+mdoop-mdon                                            7d8s21
         jout=icsf2+mdoop-mdon                                          7d21s21
         ncpy=jout-nff2                                                 7d21s21
         do i=0,ncpy-1                                                  7d21s21
          bc(nff2n+i)=bc(nff2+i)                                        7d21s21
         end do                                                         7d21s21
         call loaddoubb(ibc(ihddiag),ibc(nff2n),nsymb,mdon,mdoo,        7d21s21
     $     nvirt,multh,mysym,ibc(icsf),ibc(icsf2),nroot,-1,maxbxd,inorm,11d10s22
     $        bc,ibc,mddilow,mddihig)                                   12d13s22
        else if(mff2.lt.0)then                                          8d3s21
         jst=jout                                                       8d3s21
         ipack8=ibc(mdoubstore)                                         8d27s21
         if(ipass.eq.1)then                                             12d6s22
          iff22=mdoubstore+1                                             6d26s22
          nff22=iff22+mff2a                                              8d3s21
          nfdat=nff22+mdoop*nsymb                                        8d3s21
          ivd=nfdat+10*nsymb                                             8d3s21
          icsf=ivd+nroot*(ndoub+mdoub)                                   8d12s21
          icsf2=icsf+mdoop-mdon                                          8d3s21
         end if                                                         12d6s22
         if(nlzz.ne.0)then                                              12d6s22
          ipack4(1)=ndoubo                                              6d18s22
          ipack4(2)=mdoubo                                              6d18s22
          ndoub=ndoubo                                                  6d18s22
          mdoub=mdoubo                                                  6d18s22
         end if                                                           6d18s22
         iff22o=jout+1                                                   8d3s21
         nff22o=iff22o+mff2a                                              8d3s21
         nfdato=nff22o+mdoop*nsymb                                        8d3s21
         ivdo=nfdato+10*nsymb                                             8d3s21
         icsfo=ivdo+nroot*(ndoubo+mdoubo)                                   8d12s21
         icsf2o=icsfo+mdoop-mdon                                          8d3s21
         jout=icsf2o+2*(mdoop-mdon)                                      8d3s21
         ncpy=mff2a+mdoop*nsymb+10*nsymb                                6d24s22
         do i=1,ncpy                                                    6d26s22
          bc(jst+i)=bc(mdoubstore+i)                                    8d3s21
         end do                                                         8d3s21
         ibc(jst)=ipack8                                                6d18s22
         do i=ivdo,icsfo-1                                              6d24s22
          bc(i)=0d0                                                     8d27s21
         end do                                                         8d27s21
         ncpy=3*(mdoop-mdon)                                            6d24s22
         do i=0,ncpy-1                                                  6d24s22
          ibc(icsfo+i)=ibc(icsf+i)                                      6d24s22
         end do                                                         6d24s22
        end if                                                          7d21s21
        ibcoff=jout                                                     7d8s21
        iwavedat(5,ipp)=ibcoff-iout                                     5d14s21
       end do                                                           5d12s21
       noff=noff+npass                                                  5d12s21
      end if                                                            5d12s21
      if(mynowprog.eq.0)then                                            7d8s21
       close(unit=1)
      end if                                                            7d8s21
      return
      end
