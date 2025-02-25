c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine bdens(iden,iwd,ixw1,ixw2,ncsf,nec,nbasdws,nsymb,mdon,  11d2s22
     $     mdoop,ism,irel,irefo,norb,multh,idoubo,maxbx,maxbxd,srh,     11d2s22
     $     sr2,npadddi,lprint,wtt,bc,ibc,nvirt,ident,nh0av,nsdlk,isblk) 4d27s23
c                                                                       11d2s22
c     take wavefunction data stored a la prop and compute one particle  11d2s22
c     density.                                                          11d2s22
c                                                                       11d2s22
      implicit real*8 (a-h,o-z)                                         11d2s22
      logical lprint,lpr,ld2e                                           4d27s23
      integer*8 ipack8                                                  11d2s22
      integer*4 ipack4(2)                                               11d2s22
      equivalence (ipack8,ipack4)                                       11d2s22
      dimension ncsf(*),nbasdws(*),ism(*),irel(*),irefo(*),multh(8,8),  11d2s22
     $     idoubo(*),iwd(*),iden(*),nvirt(*),ident(*),nh0av(*)          4d25s23
      include "common.store"                                            11d2s22
      dimension i2e(512),isblk(4,*)                                     4d27s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension inv(2,8,8,8)                                            4d9s18
      common/fnd2cm/inv                                                 4d9s18
      ibcoffo=ibcoff                                                    8d3s21
      lpr=.true.                                                        11d2s22
      ld2e=.true.
      if(ld2e)then                                                      4d27s23
       nsdlk=0                                                             4d27s23
       i2s=ibcoff                                                       4d27s23
       do is1=1,nsymb                                                   4d27s23
        do is2=1,is1                                                    4d27s23
         itri12=((is1*(is1-1))/2)+is2                                   4d27s23
         is12=multh(is1,is2)                                            4d27s23
         if(is2.eq.is1)then                                             4d27s23
          n12=(nbasdws(is1)*(nbasdws(is1)+1))/2                         4d27s23
         else                                                           4d27s23
          n12=nbasdws(is1)*nbasdws(is2)                                  4d27s23
         end if                                                         4d27s23
         do is3=1,nsymb                                                 4d27s23
          do is4=1,is3                                                  4d27s23
           is34=multh(is3,is4)                                          4d27s23
           if(is12.eq.is34)then                                         4d27s23
            if(is3.eq.is4)then                                          4d27s23
             n34=(nbasdws(is3)*(nbasdws(is3)+1))/2                      4d27s23
            else                                                        4d27s23
             n34=nbasdws(is3)*nbasdws(is4)                              4d27s23
            end if                                                      4d27s23
            if(min(n12,n34).gt.0)then                                   4d27s23
             nsdlk=nsdlk+1                                              4d27s23
             isblk(1,nsdlk)=is3                                         4d27s23
             isblk(2,nsdlk)=is4                                         4d27s23
             isblk(3,nsdlk)=is1                                         4d27s23
             isblk(4,nsdlk)=is2                                         4d27s23
             i2e(nsdlk)=ibcoff                                          4d27s23
             ibcoff=ibcoff+n12*n34                                      4d27s23
            end if                                                      4d27s23
           end if                                                       4d27s23
          end do                                                        4d27s23
         end do                                                         4d27s23
        end do                                                          4d27s23
       end do                                                           4d27s23
       do i1=1,nsymb                                                    4d9s18
        do i2=1,nsymb                                                   4d9s18
         do i3=1,nsymb                                                  4d9s18
          inv(1,i1,i2,i3)=1                                             1d19s23
          do i=1,nsdlk                                                  4d9s18
           if(isblk(1,i).eq.i1.and.                                     4d9s18
     $          isblk(2,i).eq.i2.and.                                   4d9s18
     $          isblk(3,i).eq.i3)then                                   4d9s18
            inv(1,i1,i2,i3)=i                                           4d9s18
            inv(2,i1,i2,i3)=1                                           4d9s18
            go to 1                                                     4d9s18
           else if(isblk(2,i).eq.i1.and.                                4d9s18
     $           isblk(1,i).eq.i2.and.                                  4d9s18
     $           isblk(4,i).eq.i3)then                                  4d9s18
            inv(1,i1,i2,i3)=i                                           4d9s18
            inv(2,i1,i2,i3)=2                                           4d9s18
            go to 1                                                     4d9s18
           else if(isblk(2,i).eq.i2.and.                                4d9s18
     $           isblk(1,i).eq.i1.and.                                  4d9s18
     $           isblk(4,i).eq.i3)then                                  4d9s18
            inv(1,i1,i2,i3)=i                                           4d9s18
            inv(2,i1,i2,i3)=3                                           4d9s18
            go to 1                                                     4d9s18
           else if(isblk(2,i).eq.i1.and.                                4d9s18
     $           isblk(1,i).eq.i2.and.                                  4d9s18
     $           isblk(3,i).eq.i3)then                                  4d9s18
            inv(1,i1,i2,i3)=i                                           4d9s18
            inv(2,i1,i2,i3)=4                                           4d9s18
            go to 1                                                     4d9s18
           end if                                                       4d9s18
          end do                                                        4d9s18
    1     continue                                                      4d9s18
         end do                                                         4d9s18
        end do                                                          4d9s18
       end do                                                           4d9s18
       call enough('bdens.i2e',bc,ibc)                                  4d27s23
       n2espace=ibcoff-i2s                                              4d27s23
       do iz=i2s,ibcoff-1                                               4d27s23
        bc(iz)=0d0                                                      4d27s23
       end do                                                           4d27s23
      end if                                                            4d27s23
      mdoo=mdoop-1                                                      1d27s23
      icsf=1                                                            7d27s21
      icsf2=1                                                           7d27s21
      nff0=iwd(4)+iwd(9)                                                11d2s22
      jff0=iwd(4)+iwd(10)                                               11d2s22
      nroot=iwd(3)                                                      11d2s22
      wtt=wtt+dfloat(nroot)                                             11d2s22
      nec=iwd(7)                                                        11d2s22
      ilook=iwd(4)+iwd(13)                                              11d2s22
      ipack8=ibc(ilook)                                                 5d12s21
      ncsft=ipack4(1)
      ilook=ilook+1
      imff1=ilook+nroot*(ncsft+1)
      if(lpr)write(6,*)('for ket, imff1 = '),imff1
      mff1=ibc(imff1)                                                   7d9s21
      inext=imff1+1                                                     7d21s21
      if(mff1.gt.0)then                                                 7d9s21
       ihsdiag=imff1+1                                                  7d9s21
       nff1=ihsdiag+nsymb*mdoop*2                                       7d9s21
       iff1=nff1+nsymb*mdoop                                            7d9s21
       icsf=iff1+ibc(imff1)                                             7d9s21
       inext=icsf+mdoop-mdon                                            7d21s21
      end if
      mff2=ibc(inext)                                                   7d21s21
      if(lpr)write(6,*)('for ket, mff2 = '),mff2
      if(mff2.gt.0)then                                                 7d27s21
       ihddiag=inext+1                                                  7d21s21
       nff2=ihddiag+nsymb*mdoop*2                                       7d21s21
       iff2=nff2+nsymb*mdoop                                            7d21s21
       icsf=iff2+mff2                                                   7d21s21
       icsf2=icsf+mdoop-mdon                                            7d21s21
       nfdat=1                                                          8d3s21
       ivdk=1                                                           8d3s21
      else if(mff2.lt.0)then                                            8d3s21
       mdoubstore=inext+1
       ipack8=ibc(mdoubstore)                                           8d3s21
       ndoub=ipack4(1)                                                  8d12s21
       mdoub=ipack4(2)                                                  8d3s21
       mff2a=iabs(mff2)                                                 8d3s21
       iff2=mdoubstore+1                                                8d3s21
       nff2=iff2+mff2a                                                  8d3s21
       nfdat=nff2+mdoop*nsymb                                           8d3s21
       ivdk=nfdat+10*nsymb                                               8d3s21
       ivdknon=ivdk+nroot *ndoub                                        8d12s21
       icsf=ivdknon+nroot *mdoub                                        8d12s21
       icsf2=icsf+mdoop-mdon                                            8d3s21
      else                                                              8d3s21
       ihddiag=1                                                        8d3s21
       nff2=1                                                           8d3s21
       iff2=1                                                           8d3s21
       icsf=1                                                           8d3s21
       icsf2=1                                                          8d3s21
       nfdat=1                                                          8d3s21
       ivdk=1                                                           8d3s21
      end if                                                            7d27s21
      call testd(iwd,nsymb,idoubo,irefo,nvirt,nbasdws,multh,mdon,mdoop, 1d27s23
     $     ism,irel,norb,ncsf,ixw1,ixw2,maxbx,bc,ibc,ident,nh0av,ld2e,  4d27s23
     $     i2e,n2espace,nsdlk,isblk)                                    4d27s23
      call hccsfdx(bc(ilook),ncsft,ibc(nff0),ibc(jff0),ncsf,nroot,mdon, 11d2s22
     $     mdoop,ixw1,ixw2,nec,ism,irel,norb,nbasdws,idoubo,nsymb,      4d27s23
     $     iden,multh,bc,ibc,nh0av,ld2e,i2e)                            4d27s23
      if(mff1.gt.0)then                                                 1d27s23
       call hcssd1(ibc(ihsdiag),ibc(nff1),ibc(iff1),ncsf,mdon,mdoo,     1d27s23
     $     nsymb,multh,ixw1,ixw2,iden,nh0av,nroot,ism,irel,irefo,       4d25s23
     $     nvirt,iwd(2),norb,maxbx,bc,ibc,0)                            1d27s23
       call hcisd1(ibc(ihsdiag),ibc(nff1),ibc(iff1),ibc(nff0),ibc(jff0),1d27s23
     $      ncsft,ncsf,mdon,mdoo,nsymb,multh,ixw1,ixw2,iden,nh0av,      4d25s23
     $      nvirt,bc(ilook),bc,ibc,0)                                   1d27s23
      end if                                                            1d27s23
      if(mff2.lt.0)then                                                 1d27s23
       call hcddjkd1(ibc(nff2),ibc(nfdat),bc(ivdknon),nsymb,mdon,mdoo,  1d27s23
     $     nec,multh,iwd(2),nvirt,ncsf,ibc(icsf2),irel,ism,irefo,ixw1,
     $     ixw2,norb,nroot,iden,nh0av,sr2,srh,ibc(iff2),bc(iff2),ibc,   4d25s23
     $     bc,0)                                                        1d30s23
       if(mff1.gt.0)then                                                1d27s23
        call hcdsd1(ibc(ihsdiag),ibc(nff1),ibc(iff1),ibc(nff2),         1d27s23
     $      ibc(nfdat),bc(ivdknon),nsymb,mdon,mdoo,nec,multh,iwd(2),    1d27s23
     $      nvirt,ncsf,ibc(icsf2),irel,ism,irefo,ixw1,ixw2,norb,nroot,  1d27s23
     $      iden,nh0av,maxbx,sr2,ibc(iff2),bc(iff2),ibc,bc,0)           1d30s23
       end if                                                           1d27s23
      end if                                                            1d27s23
      ibcoff=ibcoffo                                                    11d2s22
      return                                                            11d2s22
      end                                                               11d2s22
      subroutine testd(iwket,nsymb,idoubo,irefo,nvirt,nbasdws,multh,    1d27s23
     $     mdon,mdoop,ism,irel,norb,ncsf,ixw1,ixw2,maxbx,bc,ibc,ident,  4d25s23
     $     nh0av,ld2e,i2e,n2espace,nsdlk,isblk)                         4d27s23
      implicit real*8 (a-h,o-z)                                         1d27s23
      integer*8 ipack8
      logical ld2e                                                      4d27s23
      equivalence (ipack8,ipack4)
      dimension iwket(*),ipack4(2),idoubo(*),irefo(*),nvirt(*),         1d27s23
     $     nbasdws(*),multh(8,8),ism(*),irel(*),ncsf(*),                1d27s23
     $     iden(8),ident(8),nh0av(*),i2e(*),isblk(4,*)                  4d27s23
      include "common.store"
      write(6,*)('hi, my name is testd')
      look=1752477
      write(6,*)('what''s under '),look,bc(look)
      ibcoffo=ibcoff
      nec=iwket(7)                                                      7d27s21
      nff0k=iwket(4)+iwket(9)                                           8d25s21
      jff0k=iwket(4)+iwket(10)                                          8d25s21
      nrootk=iwket(3)                                                   8d25s21
      ilook=iwket(4)+iwket(13)                                          8d25s21
      ipack8=ibc(ilook)                                                 5d12s21
      ncsftk=ipack4(1)                                                  5d12s21
      write(6,*)('ncsftk = '),ncsftk
      ilook=ilook+1                                                     5d12s21
      write(6,*)('going for ')
      imff1k=ilook+nrootk*(ncsftk+1)                                    8d19s21
      mff1k=ibc(imff1k)                                                  7d15s21
      write(6,*)('imff1k '),imff1k,mff1k
      inextk=imff1k+1
      ntotk=ncsftk
      write(6,*)('setting ntotk to ncsftk = '),ntotk
      nntot=0
      do isb=1,nsymb
       nntot=nntot+nbasdws(isb)-idoubo(isb)
      end do
      write(6,*)('nntot = '),nntot
      itest=ibcoff
      ibcoff=itest+max(nrootb,nrootk)*nrootk
      call enough('testme.  1',bc,ibc)
      call dgemm('t','n',nrootk,nrootk,ncsftk,1d0,bc(ilook),ncsftk,
     $      bc(ilook),ncsftk,0d0,bc(itest),nrootk,
     d' testme.  1')
      write(6,*)('ortho test for ket after internals: ')
      call prntm2(bc(itest),nrootk,nrootk,nrootk)
      if(mff1k.gt.0)then
       ihsdiagk=imff1k+1                                                  7d9s21
       nff1k=ihsdiagk+nsymb*mdoop*2                                       7d9s21
       iff1k=nff1k+nsymb*mdoop                                            7d9s21
       icsfk=iff1k+ibc(imff1k)                                          10d25s21
       inextk=icsfk+mdoop-mdon                                          10d25s21
       ivx=ibcoff                                                       7d16s21
       ibcoff=ivx+max(nrootb,nrootk)                                                 7d16s21
       call enough('testme.  2',bc,ibc)
       write(6,*)('sotest for ket')
       call sotest(ibc(ihsdiagk),ibc(nff1k),ibc(icsfk),iwket(2),        10d25s21
     $     nvirt,bc(itest),nrootk,mdon,mdoop,nsymb,multh,nsingk,bc(ivx),11d10s22
     $      bc,ibc)                                                     11d10s22
       write(6,*)('nsing from sotest: '),nsingk
       ntotk=ntotk+nsingk
       write(6,*)('adding nsingk = '),nsingk,('to ntotk to get '),
     $      ntotk
       call enough('testme.  3',bc,ibc)
      else                                                              7d15s21
       ihsdiagk=1                                                        7d9s21
       nff1k=1                                                           7d9s21
       iff1k=1                                                           7d9s21
      end if                                                            7d9s21
      mff2k=ibc(inextk)
      write(6,*)('what we have for mff2: '),mff2k
      if(mff2k.gt.0)then
       ihddiagk=inextk+1                                                  7d21s21
       nff2k=ihddiagk+nsymb*mdoop*2                                       7d21s21
       iff2k=nff2k+nsymb*mdoop                                            7d21s21
       icsfk=iff2k+mff2k                                                   7d21s21
       icsf2=icsfk+mdoop-mdon                                            7d21s21
       ivx2=ibcoff
       ibcoff=ivx2+max(nrootb,nrootk)                                                8d3s21
       call enough('testme.  4',bc,ibc)
       write(6,*)('sotest2 for ket')
       call prntm2(bc(itest),nrootk,nrootk,nrootk)
       call sotest2(ibc(ihddiagk),ibc(nff2k),ibc(icsfk),ibc(icsf2),        7d21s21
     $      iwket(2),nvirt,bc(itest),nrootk,mdon,mdoop,nsymb,multh,     7d21s21
     $      ndoubk,bc(ivx2),bc,ibc)                                     11d10s22
       ntotk=ntotk+ndoubk
       write(6,*)('adding ndoubk = '),ndoubk,(' to ntotk to get '),
     $      ntotk
      else if(mff2k.lt.0)then
       mff2a=-mff2k                                                      8d3s21
       mdoubstore=inextk+1                                               8d3s21
       ipack8=ibc(mdoubstore)                                           8d3s21
       ndoubk=ipack4(1)                                                  8d3s21
       mdoubk=ipack4(2)                                                  8d3s21
       write(6,*)('ndoubk,mdoubk: '),ndoubk,mdoubk
       iff22k=mdoubstore+1                                               8d3s21
       nff22k=iff22k+mff2a                                                8d3s21
       nfdatk=nff22k+mdoop*nsymb                                          8d3s21
       ivdk=nfdatk+10*nsymb                                               8d3s21
       icsfk=ivdk+nrootk*(ndoubk+mdoubk)                                     8d12s21
       icsf2=icsfk+mdoop-mdon                                            8d3s21
       ivx2=ibcoff
       ibcoff=ivx2+max(nrootk,nootb)                                                8d3s21
       call enough('testme.  6',bc,ibc)
       ivdknon=ivdk+nrootk*ndoubk                                           8d12s21
       write(6,*)('sotest2c for ket')
       call sotest2c(bc(iff22k),ibc(nff22k),ibc(nfdatk),bc(ivdk),
     $       ibc(icsfk),ibc(icsf2),iwket(2),nvirt,bc(itest),nrootk,mdon,8d25s21
     $       mdoop,nsymb,multh,ndoubxk,bc(ivx2),bc(iff22k),bc,ibc)      11d10s22
       call tofrob(bc(ivdk),bc(ivdknon),nrootk,ibc(nfdatk),nvirt,nsymb, 8d26s21
     $      multh,iwket(2),1,ndoubk,mdoubk,bc(iff22k),bc,ibc)           11d10s22
       ntotk=ntotk+ndoubxk
       write(6,*)('adding ndoubxk = '),ndoubxk,('to ntotk to get '),
     $      ntotk
      else
       ihddiagk=1                                                        7d22s21
       ihddiagb=1
       nff2k=1                                                           7d22s21
       iff2k=1                                                           7d22s21
       nff2b=1                                                           7d22s21
       iff2b=1                                                           7d22s21
       icsf2=1                                                          7d22s21
      end if
      irsum=ibcoff
      ibcoff=irsum+nrootk**2
      do i=irsum,ibcoff-1
       bc(i)=0d0
      end do
      write(6,*)('going after vnew: '),ibcoff,ntotk
      ivnew=ibcoff
      nff10k=ivnew+ntotk*nrootk
      ismv=nff10k+mdoop*3
      irelv=ismv+nntot
      iff10k=irelv+nntot
      write(6,*)('iff10k = '),iff10k,loc(bc(iff10k))
      if(mff1k.ne.0.and.mff2k.eq.0)then                                   7d27s21
       call convert1tor(ibc(ihsdiagk),ibc(nff1k),ibc(iff1k),
     $      ibc(icsfk),iwket(2),nvirt,nrootk,mdon,mdoop,nsymb,multh,
     $      ntotk,bc(ivnew),ibc(nff10k),ism,ibc(ismv),irel,ibc(irelv),
     $      iff10k,norb,irefo,ibc(iff10k),bc(ivx),bc(ilook),ncsftk,
     $      ibc(nff0k),ibc(jff0k),bc,ibc)                               11d14s22
      else if(mff1k.eq.0.and.mff2k.ne.0)then                              7d27s21
       if(mff2k.gt.0)then                                                8d3s21
        call convert2tor(ibc(ihddiagk),ibc(nff2k),ibc(iff2k),           8d25s21
     $       ibc(icsfk),ibc(icsf2),iwket(2),nvirt,nrootk,mdon,mdoop,    8d25s21
     $       nsymb,multh,ntotk,bc(ivnew),ibc(nff10k),ism,ibc(ismv),irel,8d25s21
     $       ibc(irelv),iff10k,norb,irefo,ibc(iff10k),bc(ivx2),         8d25s21
     $       bc(ilook),ncsftk,ibc(nff0k),ibc(jff0k),bc,ibc)             11d14s22
       else                                                             8d3s21
        call convert2torc(bc(iff22k),ibc(nff22k),ibc(nfdatk),
     $        bc(ivdknon),ibc(icsfk),ibc(icsf2),iwket(2),nvirt,nrootk,  8d25s21
     $        mdon,mdoop,nsymb,multh,bc(ivx2),ntotk,bc(ivnew),
     $        ibc(nff10k),ism,ibc(ismv),irel,ibc(irelv),iff10k,norb,
     $        irefo,ibc(iff10k),bc(ilook),ncsftk,ibc(nff0k),ibc(jff0k),
     $        bc(iff22k),bc(irsum),bc,ibc)                              11d14s22
       end if                                                           8d3s21
      else if(mff1k.ne.0.and.mff2k.ne.0)then                              7d27s21
       if(mff2k.gt.0)then                                                8d3s21
        write(6,*)('b4 convert12tor, ntotk = '),ntotk
        call convert12tor(ibc(ihsdiagk),ibc(ihddiagk),ibc(nff1k),       8d25s21
     $      ibc(iff1k),ibc(nff2k),ibc(iff2k),ibc(icsfk),                8d25s21
     $      ibc(icsf2),iwket(2),nvirt,nrootk,mdon,mdoop,nsymb,multh,    8d25s21
     $      ntotk,bc(ivnew),ibc(nff10k),ism,ibc(ismv),irel,ibc(irelv),  8d25s21
     $      iff10k,norb,irefo,ibc(iff10k),bc(ivx2),bc(ilook),ncsftk,       7d22s21
     $      ibc(nff0k),ibc(jff0k),bc,ibc)                               11d14s22
        write(6,*)('after convert12tor, ntotk = '),ntotk
       else
        write(6,*)('calling convert12torc ')
        call convert12torc(bc(iff22k),ibc(nff22k),ibc(nfdatk),
     $        bc(ivdknon),ibc(icsfk),ibc(icsf2),iwket(2),nvirt,nrootk,  8d25s21
     $        mdon,mdoop,nsymb,multh,bc(ivx2),ntotk,bc(ivnew),
     $        ibc(nff10k),ism,ibc(ismv),irel,ibc(irelv),iff10k,norb,
     $        irefo,ibc(iff10k),bc(ilook),ncsftk,ibc(nff0k),ibc(jff0k),
     $        bc(iff22k),bc(irsum),ibc(ihsdiagk),ibc(nff1k),ibc(iff1k), 11d14s22
     $       bc,ibc)                                                    11d14s22
        write(6,*)('back')
       end if
      else                                                              5d18s23
       call convertto(ncsf,nrootk,mdon,mdoop,bc(ivx2),bc(ivnew),        5d18s23
     $      ibc(nff10k),ism,ibc(ismv),irel,ibc(irelv),iff10k,norb,      5d18s23
     $      ibc(iff10k),bc(ilook),ncsftk,ibc(nff0k),ibc(jff0k),bc,ibc)  5d18s23
      end if                                                            7d27s21
      write(6,*)('ibcoff after conversion: '),ibcoff
      write(6,*)('vnew for ortho test: ')
      call prntm2(bc(ivnew),ntotk,nrootk,ntotk)
      call dgemm('t','n',nrootk,nrootk,ntotk,1d0,bc(ivnew),ntotk,
     $      bc(ivnew),ntotk,0d0,bc(ibcoff),nrootk,
     d' testme.  2')
      call prntm2(bc(ibcoff),nrootk,nrootk,nrootk)
      ibd=ibcoff
      do isb=1,nsymb                                                    1d27s23
       write(6,*)('isb, '),isb,nbasdws(isb)
       iden(isb)=ibcoff                                                 1d27s23
       ibcoff=iden(isb)+nh0av(isb)*nh0av(isb)                           4d25s23
      end do                                                            1d27s23
      call enough('testd.1',bc,ibc)                                     1d27s23
      do iz=ibd,ibcoff-1
       bc(iz)=0d0
      end do
      write(6,*)('zeroing '),ibd,ibcoff-1
      nbd=ibcoff-ibd
      if(ld2e)then                                                      4d27s23
       write(6,*)('create copy of irefo ')
       irefocpy=ibcoff                                                  4d27s23
       ibcoff=irefocpy+nsymb                                            4d27s23
       call enough('testd.irefocpy',bc,ibc)                             4d27s23
       jrefocpy=irefocpy-1                                              4d27s23
       do i=1,nsymb                                                     4d27s23
        ibc(jrefocpy+i)=irefo(i)                                        4d27s23
        irefo(i)=nh0av(i)                                               5d18s23
       end do                                                           4d27s23
       write(6,*)('refo now: '),(irefo(i),i=1,nsymb)
      write(6,*)('what''s under '),look,bc(look)
      end if                                                            4d27s23
      write(6,*)('what''s under '),look,bc(look)
      call hccsfdx(bc(ivnew),ntotk,ibc(nff10k),ibc(iff10k),ncsf,nrootk, 1d27s23
     $     mdon,mdoop,ixw1,ixw2,nec,ibc(ismv),ibc(irelv),nntot,         4d27s23
     $     nbasdws,idoubo,nsymb,ident,multh,bc,ibc,nh0av,ld2e,i2e)      4d27s23
      write(6,*)('what''s under '),look,bc(look)
      do i=0,nbd-1
       bc(iden(1)+i)=bc(ident(1)+i)
      end do
      write(6,*)('global summing '),ibd,ibd+nbd-1
      write(6,*)('what''s under '),look,bc(look)
      call dws_gsumf(bc(ibd),nbd)
      write(6,*)('what''s under '),look,bc(look)
      trac=0d0
      do isb=1,nsymb
       if(nbasdws(isb).gt.0)then
        write(6,*)('raw density for symmetry '),isb,iden(isb)
        call prntm2(bc(iden(isb)),nh0av(isb),nh0av(isb),nh0av(isb))     4d25s23
        do i=1,nh0av(isb)-1                                             5d15s23
         do j=0,i-1                                                     1d30s23
          ji=iden(isb)+j+nh0av(isb)*i                                   1d30s23
          ij=iden(isb)+i+nh0av(isb)*j                                   1d30s23
          avg=0.5d0*(bc(ji)+bc(ij))                                     1d30s23
          bc(ji)=avg                                                    1d30s23
          bc(ij)=avg                                                    1d30s23
         end do                                                         1d30s23
        end do                                                          1d30s23
        write(6,*)('symmetrized ')
        call prntm2(bc(iden(isb)),nh0av(isb),nh0av(isb),
     $       nh0av(isb))
        do i=0,nh0av(isb)-1                                             4d25s23
         ii=iden(isb)+i*(nh0av(isb)+1)                                  4d25s23
         trac=trac+bc(ii)
        end do
        write(6,*)('trace so far '),trac
        ieig=ibcoff                                                     1d30s23
        ivec=ieig+nh0av(isb)                                            4d25s23
        isx=ivec+nh0av(isb)*nh0av(isb)                                  4d25s23
        ibcoff=isx+nh0av(isb)                                           4d25s23
        call enough('bdens.d',bc,ibc)                                   1d30s23
        call diagx(nh0av(isb),bc(iden(isb)),bc(ieig),bc(ivec),          4d25s23
     $       ibc(isx),bc,ibc)                                           1d30s23
        write(6,*)('no occupations: ')
        do i=0,nh0av(isb)-1                                             4d25s23
         im=nh0av(isb)-1-i                                              4d25s23
         bc(ivec+im)=bc(ieig+i)                                         4d25s23
        end do                                                          4d25s23
        call prntm2(bc(ivec),1,nh0av(isb),1)                            4d25s23
        ibcoff=ieig                                                     1d30s23
       end if
      end do
      if(ld2e)then                                                      4d27s23
       write(6,*)('reset irefo '),i2e(1),bc(look)
       jrefocpy=irefocpy-1                                              4d27s23
       do i=1,nsymb                                                     4d27s23
        irefo(i)=ibc(jrefocpy+i)                                        4d27s23
       end do                                                           4d27s23
       write(6,*)('refo now: '),(irefo(i),i=1,nsymb)
       ist=i2e(1)
       write(6,*)('now gsumfing '),ist,('to '),ist+n2espace-1
       call dws_gsumf(bc(ist),n2espace)                                 4d27s23
       write(6,*)('look after gsumf '),bc(look)
       do is=1,nsdlk                                                    5d15s23
        if(isblk(1,is).ne.isblk(2,is))then                              5d15s23
         do iss=is+1,nsdlk                                              5d15s23
          if(isblk(1,is).eq.isblk(1,iss).and.isblk(2,is).eq.isblk(2,iss)5d15s23
     $         .and.isblk(3,is).eq.isblk(4,iss))then                    5d15s23
           ncol=nh0av(isblk(3,is))*nh0av(isblk(4,is))                             5d15s23
           nrow=nh0av(isblk(1,is))*nh0av(isblk(2,is))                             5d15s23
           write(6,*)('combining '),(isblk(j,is),j=1,4)                 5d15s23
           call prntm2(bc(i2e(is)),nrow,ncol,nrow)                      5d15s23
           write(6,*)('and '),(isblk(j,iss),j=1,4)                      5d15s23
           call prntm2(bc(i2e(iss)),nrow,ncol,nrow)                     5d15s23
           do i4=0,nh0av(isblk(4,is))-1                                      5d15s23
            do i3=0,nh0av(isblk(3,is))-1                                     5d15s23
             i34=i2e(is)+nrow*(i3+nh0av(isblk(3,is))*i4)                     5d15s23
             i43=i2e(iss)+nrow*(i4+nh0av(isblk(4,is))*i3)                    5d15s23
             do i12=0,nrow-1                                            5d15s23
              sum=bc(i34+i12)+bc(i43+i12)                               5d15s23
              if(i34+12.eq.look.or.i43+i12.eq.look)write(6,*)
     $             ('summing '),bc(i34+i12),bc(i43+i12),sum
              bc(i34+i12)=sum                                           5d15s23
              bc(i43+i12)=sum                                           5d15s23
             end do                                                     5d15s23
            end do                                                      5d15s23
           end do                                                       5d15s23
          end if                                                        5d15s23
         end do                                                         5d15s23
        end if                                                          5d15s23
       end do                                                           5d15s23
       write(6,*)('stamp1 '),bc(look)
       do ib=1,nsdlk                                                     8d19s14
        n1=isblk(1,ib)                                                   8d19s14
        n2=isblk(2,ib)                                                   8d19s14
        n3=isblk(3,ib)
        n4=isblk(4,ib)                                                   8d19s14
        if(n1.eq.n2)then
         nnn=(nh0av(n1)*(nh0av(n1)+1))/2
        else
         nnn=nh0av(n1)*nh0av(n2)
        end if
        if(n3.eq.n4)then
         mmm=(nh0av(n3)*(nh0av(n3)+1))/2
        else
         mmm=nh0av(n3)*nh0av(n4)
        end if
        if(nnn*mmm.gt.0)then
         if(min(n1,n2).eq.min(n3,n4).and.max(n1,n2).eq.max(n3,n4))then   8d20s14
          do j=0,nnn-1                                                    8d20s14
           do k=0,j-1                                                     8d20s14
            iad1=i2e(ib)+k+nnn*j                                        8d20s14
            iad2=i2e(ib)+j+nnn*k                                        8d20s14
            avg=0.5d0*(bc(iad1)+bc(iad2))                                 8d20s14
            if(iad1.eq.look.or.iad2.eq.look)write(6,*)('averaging '),
     $           bc(iad1),bc(iad2),avg,iad1-i2e(ib),iad2-i2e(ib)
            bc(iad1)=avg                                                  8d20s14
            bc(iad2)=avg                                                  8d20s14
           end do                                                         8d20s14
          end do                                                          8d20s14
         end if                                                           8d20s14
         if(.not.(n1.eq.n2.and.n3.eq.n4.and.n1.eq.n3))then                9d18s14
          do io=ib+1,nsdlk                                                 9d18s14
           if(isblk(1,io).eq.n3.and.isblk(2,io).eq.n4.and.                 9d18s14
     $       isblk(3,io).eq.n1.and.isblk(4,io).eq.n2)then               9d18s14
            do j=0,nnn-1                                                    9d18s14
             do k=0,mmm-1                                                   9d18s14
              iad1=i2e(ib)+j+nnn*k                                       9d18s14
              iad2=i2e(io)+k+mmm*j                                       9d18s14
              avg=0.5d0*(bc(iad1)+bc(iad2))                                 9d18s14
              if(iad1.eq.look.or.iad2.eq.look)write(6,*)
     $             ('2nd averaging'),bc(iad1),bc(iad2),avg
              bc(iad1)=avg                                                  9d18s14
              bc(iad2)=avg                                                  9d18s14
             end do                                                         9d18s14
            end do                                                          9d18s14
           end if                                                           9d18s14
          end do                                                            9d18s14
         end if                                                            9d18s14
        end if                                                            9d18s14
       end do                                                            9d18s14
       write(6,*)('stamp2 '),bc(look)
       do ib=1,nsdlk                                                     8d19s14
        n1=isblk(1,ib)                                                   8d19s14
        n2=isblk(2,ib)                                                   8d19s14
        n3=isblk(3,ib)
        n4=isblk(4,ib)                                                   8d19s14
        if(n1.eq.n2)then
         nnn=(nh0av(n1)*(nh0av(n1)+1))/2
        else
         nnn=nh0av(n1)*nh0av(n2)
        end if
        if(n3.eq.n4)then
         mmm=(nh0av(n3)*(nh0av(n3)+1))/2
        else
         mmm=nh0av(n3)*nh0av(n4)
        end if
        if(nnn*mmm.gt.0)then
         if(n1.eq.n2)then                                                6d16s22
          do i2=0,nh0av(n3)-1                                            4d10s17
           do i1=0,i2                                                    4d10s17
            fact=0.5d0                                                   4d10s17
            icol=i2e(ib)+nnn*(((i2*(i2+1))/2)+i1)                       4d10s17
            if(i1.eq.i2)fact=1d0                                         4d10s17
            facth=fact*0.5d0                                             4d10s17
            do i4=0,nh0av(n1)-1                                          4d10s17
             do i3=0,i4-1                                                4d10s17
              if(icol.eq.look)write(6,*)('scaling by facth: '),facth
              bc(icol)=bc(icol)*facth                                    4d10s17
              icol=icol+1                                                4d10s17
             end do                                                      4d10s17
             if(icol.eq.look)write(6,*)('scaling by fact: '),fact
             bc(icol)=bc(icol)*fact                                      4d10s17
             icol=icol+1                                                 4d10s17
            end do                                                       4d10s17
           end do                                                        4d10s17
          end do                                                         4d10s17
         end if                                                          4d10s17
        end if
       end do
       write(6,*)('stamp3 '),bc(look)
       do is=1,nsdlk
        write(6,*)('master: '),is,(isblk(j,is),j=1,4),i2e(is)
        if(isblk(1,is).eq.isblk(2,is))then
         n12=(nh0av(isblk(1,is))*(nh0av(isblk(1,is))+1))/2              4d27s23
         n34=(nh0av(isblk(3,is))*(nh0av(isblk(3,is))+1))/2              4d27s23
         m12=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2
         m34=(irefo(isblk(3,is))*(irefo(isblk(3,is))+1))/2              4d27s23
        else                                                                4d27s23
         n12=nh0av(isblk(1,is))*nh0av(isblk(2,is))                      4d27s23
         n34=nh0av(isblk(3,is))*nh0av(isblk(4,is))                      4d27s23
         m12=irefo(isblk(1,is))*irefo(isblk(2,is))                      4d27s23
         m34=irefo(isblk(3,is))*irefo(isblk(4,is))                      5d18s23
        end if
        nm34=irefo(isblk(3,is))*nvirt(isblk(4,is))
        nm43=irefo(isblk(4,is))*nvirt(isblk(3,is))
        nn34=nvirt(isblk(3,is))*nvirt(isblk(4,is))
        call prntm2(bc(i2e(is)),n12,n34,n12)
        if(min(m12,m34).gt.0)then
         do i4=0,irefo(isblk(4,is))-1
          if(isblk(3,is).eq.isblk(4,is))then
           i3top=i4
          else
           i3top=irefo(isblk(3,is))-1
          end if
          do i3=0,i3top
           if(isblk(3,is).eq.isblk(4,is))then
            i34=((i4*(i4+1))/2)+i3
            j34=i34
           else
            i34=i3+nh0av(isblk(3,is))*i4
            j34=i3+irefo(isblk(3,is))*i4
           end if
           do i2=0,irefo(isblk(2,is))-1
            if(isblk(3,is).eq.isblk(4,is))then
             i1top=i2
            else
             i1top=irefo(isblk(1,is))-1
            end if
            do i1=0,i1top
             if(isblk(3,is).eq.isblk(4,is))then
              i12=((i2*(i2+1))/2)+i1
              j12=i12
             else
              i12=i1+nh0av(isblk(1,is))*i2
              j12=i1+irefo(isblk(1,is))*i2
             end if
             iad1=i2e(is)+i12+n12*i34
             iad2=ibcoff+j12+m12*j34
             bc(iad2)=bc(iad1)
             if(abs(bc(iad2)).gt.1d-10)write(6,*)('getting '),bc(iad2),
     $            ('from '),iad1-i2e(is),bc(iad1),iad1
            end do
           end do
          end do
         end do
         write(6,*)('internal part '),(isblk(j,is),j=1,4)
         call prntm2(bc(ibcoff),m12,m34,m12)
         if(isblk(3,is).ne.isblk(4,is))then
          do i4=0,irefo(isblk(3,is))-1
           do i3=0,irefo(isblk(4,is))-1
            i34=i4+nh0av(isblk(3,is))*i3
            j34=i3+irefo(isblk(4,is))*i4
            do i2=0,irefo(isblk(2,is))-1
             do i1=0,irefo(isblk(1,is))-1
              i12=i1+nh0av(isblk(1,is))*i2
              j12=i1+irefo(isblk(1,is))*i2
              iad1=i2e(is)+i12+n12*i34
              iad2=ibcoff+j12+m12*j34
              bc(iad2)=bc(iad1)
              if(abs(bc(iad2)).gt.1d-10)write(6,*)('getting '),bc(iad2),
     $            ('from '),iad1-i2e(is),bc(iad1),iad1
             end do
            end do
           end do
          end do
          write(6,*)('internal part '),(isblk(j,is),j=1,2),isblk(4,is),
     $         isblk(3,is)
          call prntm2(bc(ibcoff),m12,m34,m12)
         end if
        end if
        if(min(m12,nm34).gt.0)then
         do i4=0,nvirt(isblk(4,is))-1
          i4p=i4+irefo(isblk(4,is))
          do i3=0,irefo(isblk(3,is))-1
           if(isblk(3,is).eq.isblk(4,is))then
            i34=((i4p*(i4p+1))/2)+i3
           else
            i34=i3+nh0av(isblk(3,is))*i4p
           end if
           j34=i3+irefo(isblk(3,is))*i4
           do i2=0,irefo(isblk(2,is))-1
            if(isblk(1,is).eq.isblk(2,is))then
             i1top=i2
            else
             i1top=irefo(isblk(1,is))-1
            end if
            do i1=0,i1top
             if(isblk(1,is).eq.isblk(2,is))then
              i12=((i2*(i2+1))/2)+i1
              j12=i12
              if(i1.eq.i2)then
               ff=2d0
              else
               ff=4d0
              end if
             else
              i12=i1+nh0av(isblk(1,is))*i2
              j12=i1+irefo(isblk(1,is))*i2
              ff=1d0
             end if
             iad1=i2e(is)+i12+n12*i34
             iad2=ibcoff+i4+nvirt(isblk(4,is))*(j12+m12*i3)
             if(abs(bc(iad1)).gt.1d-10)write(6,*)('getting '),bc(iad1),
     $            bc(iad1)*ff,('from '),iad1,i1,i2,ff
             bc(iad2)=bc(iad1)*ff
            end do
           end do
          end do
         end do
         sz=0d0
         do i=0,m12*nm34-1
          sz=sz+bc(ibcoff+i)**2
         end do
         sz=sqrt(sz/dfloat(m12*nm34))
         if(sz.gt.1d-10)then
          write(6,*)('onex part'),(isblk(j,is),j=1,4),('a')
          call prntm2(bc(ibcoff),nvirt(isblk(4,is)),
     $         m12*irefo(isblk(3,is)),nvirt(isblk(4,is)))
         end if
        end if                                                          7d6s23
        write(6,*)('btest '),m12,nm43,isblk(3,is),isblk(4,is)
        if(min(m12,nm43).gt.0.and.isblk(3,is).ne.isblk(4,is))then       7d6s23
         do i4=0,nvirt(isblk(3,is))-1
          i4p=i4+irefo(isblk(3,is))
          do i3=0,irefo(isblk(4,is))-1
           i34=i4p+nh0av(isblk(3,is))*i3
           j34=i3+irefo(isblk(4,is))*i4
           do i2=0,irefo(isblk(2,is))-1
            i1top=irefo(isblk(1,is))-1
            do i1=0,i1top
             i12=i1+nh0av(isblk(1,is))*i2
             j12=i1+irefo(isblk(1,is))*i2
             iad1=i2e(is)+i12+n12*i34
             iad2=ibcoff+i4+nvirt(isblk(3,is))*(j12+m12*i3)
             if(abs(bc(iad1)).gt.1d-10)write(6,*)('getting '),bc(iad1),
     $            ('from '),iad1,i1,i2
             bc(iad2)=bc(iad1)
            end do
           end do
          end do
         end do
         sz=0d0
         do i=0,m12*nm43-1
          sz=sz+bc(ibcoff+i)**2
         end do
         sz=sqrt(sz/dfloat(m12*nm43))
         if(sz.gt.1d-10)then
          write(6,*)('onex part'),(isblk(j,is),j=1,2),isblk(4,is),
     $         isblk(3,is),('b')
          call prntm2(bc(ibcoff),nvirt(isblk(3,is)),
     $         m12*irefo(isblk(4,is)),nvirt(isblk(3,is)))
         end if
        end if
        if(min(m12,nn34).gt.0)then
         fact=1d0
         if(isblk(1,is).eq.isblk(3,is)*50.and.
     $      isblk(2,is).eq.isblk(4,is).or.
     $        (isblk(1,is).eq.isblk(2,is).and.
     $        isblk(3,is).eq.isblk(4,is)))fact=2d0
         do i4=0,nvirt(isblk(4,is))-1
          i4p=i4+irefo(isblk(4,is))
          do i3=0,nvirt(isblk(3,is))-1
           i3p=i3+irefo(isblk(3,is))
           j34=i3+nvirt(isblk(3,is))*i4
           if(isblk(3,is).eq.isblk(4,is))then
            ix=max(i3p,i4p)
            in=min(i3p,i4p)
            i34=((ix*(ix+1))/2)+in
           else
            i34=i3p+nh0av(isblk(3,is))*i4p
           end if
           do i2=0,irefo(isblk(2,is))-1
            if(isblk(1,is).eq.isblk(2,is))then
             i1top=i2
            else
             i1top=irefo(isblk(1,is))-1
            end if
            do i1=0,i1top
             ff=1d0                                                       5d25s23
             if(isblk(1,is).eq.isblk(2,is))then
              i12=((i2*(i2+1))/2)+i1
              j12=i12
              if(i1.ne.i2)then                                          5d25s23
               ff=2d0                                                     5d25s23
              end if                                                      5d25s23
             else
              i12=i1+nh0av(isblk(1,is))*i2
              j12=i1+irefo(isblk(1,is))*i2
             end if
             iad1=i2e(is)+i12+n12*i34
             iad2=ibcoff+j34+nn34*j12
             bc(iad2)=bc(iad1)*fact*ff
             if(abs(bc(iad2)).gt.1d-10)write(6,*)('getting '),bc(iad2),
     $            ('from '),iad1-i2e(is),bc(iad1),iad1,ff
            end do
           end do
          end do
         end do
         write(6,*)('J part '),(isblk(j,is),j=1,4),fact,('a')
         call prntm2(bc(ibcoff),nn34,m12,nn34,m12)
         if(isblk(4,is).ne.isblk(3,is))then
          do i4=0,nvirt(isblk(4,is))-1
           i4p=i4+irefo(isblk(4,is))
           do i3=0,nvirt(isblk(3,is))-1
            i3p=i3+irefo(isblk(3,is))
            j34=i4+nvirt(isblk(4,is))*i3
            if(isblk(3,is).eq.isblk(4,is))then
             ix=max(i3p,i4p)
             in=min(i3p,i4p)
             i34=((ix*(ix+1))/2)+in
            else
             i34=i3p+nh0av(isblk(3,is))*i4p
            end if
            do i2=0,irefo(isblk(2,is))-1
             if(isblk(1,is).eq.isblk(2,is))then
              i1top=i2
             else
              i1top=irefo(isblk(1,is))-1
             end if
             do i1=0,i1top
              if(isblk(1,is).eq.isblk(2,is))then
               i12=((i2*(i2+1))/2)+i1
               j12=i12
              else
               i12=i1+nh0av(isblk(1,is))*i2
               j12=i1+irefo(isblk(1,is))*i2
              end if
              iad1=i2e(is)+i12+n12*i34
              iad2=ibcoff+j34+nn34*j12
              bc(iad2)=bc(iad1)*fact
              if(abs(bc(iad2)).gt.1d-10)write(6,*)('gettingb '),
     $             bc(iad2),('from '),iad1-i2e(is),bc(iad1)
             end do
            end do
           end do
          end do
          write(6,*)('J part '),(isblk(j,is),j=1,2),isblk(4,is),
     $         isblk(3,is),fact,('b')
          call prntm2(bc(ibcoff),nn34,m12,nn34,m12)
         end if
        end if
        m14=irefo(isblk(1,is))*irefo(isblk(4,is))
        m13=irefo(isblk(1,is))*irefo(isblk(3,is))
        nn32=nvirt(isblk(3,is))*nvirt(isblk(2,is))
        nn24=nvirt(isblk(4,is))*nvirt(isblk(2,is))
        m24=irefo(isblk(2,is))*irefo(isblk(4,is))
        m23=irefo(isblk(2,is))*irefo(isblk(3,is))
        nn13=nvirt(isblk(1,is))*nvirt(isblk(3,is))
        nn14=nvirt(isblk(1,is))*nvirt(isblk(4,is))
c     recall K_nm^ab is (nb|ma) with ab distributed across procs.
c     recall K_13^42 is (12|34) with ab distributed across procs.
         ff=4d0
         if(isblk(1,is).eq.isblk(3,is).and.isblk(1,is).ne.isblk(4,is))
     $        ff=1d0
         do i4=0,nvirt(isblk(2,is))-1
          i4p=i4+irefo(isblk(2,is))
          do i3=0,nvirt(isblk(4,is))-1
           i3p=i3+irefo(isblk(4,is))
           j34=i3+nvirt(isblk(4,is))*i4
           do i2=0,irefo(isblk(3,is))-1
            do i1=0,irefo(isblk(1,is))-1
             if(isblk(1,is).eq.isblk(2,is))then
              i12=((i4p*(i4p+1))/2)+i1
              i34=((i3p*(i3p+1))/2)+i2
             else
              i12=i1+nh0av(isblk(1,is))*i4p
              i34=i2+nh0av(isblk(3,is))*i3p
             end if
             j12=i1+irefo(isblk(1,is))*i2
             iad1=i2e(is)+i12+n12*i34
             iad2=ibcoff+j34+nn24*j12
             bc(iad2)=bc(iad1)*ff
             if(abs(bc(iad1)).gt.1d-10)write(6,*)('get '),bc(iad2),
     $            ('from '),iad1,iad1-i2e(is),bc(iad1)
            end do
           end do
          end do
         end do
         write(6,*)('K part '),isblk(1,is),isblk(3,is),
     $        isblk(4,is),isblk(2,is),('a')
         call prntm2(bc(ibcoff),nn24,m13,nn24)
c     recall K_13^42 is (12|34) with ab distributed across procs.
c     but with 1<>2, ie K_23^41
         do i4=0,nvirt(isblk(1,is))-1
          i4p=i4+irefo(isblk(1,is))
          do i3=0,nvirt(isblk(4,is))-1
           i3p=i3+irefo(isblk(4,is))
           j34=i3+nvirt(isblk(4,is))*i4
           do i2=0,irefo(isblk(3,is))-1
            do i1=0,irefo(isblk(2,is))-1
             if(isblk(1,is).eq.isblk(2,is))then
              i12=((i4p*(i4p+1))/2)+i1
              i34=((i3p*(i3p+1))/2)+i2
             else
              i12=i4p+nh0av(isblk(1,is))*i1
              i34=i2+nh0av(isblk(3,is))*i3p
             end if
             j12=i1+irefo(isblk(2,is))*i2
             iad1=i2e(is)+i12+n12*i34
             iad2=ibcoff+j34+nn14*j12
             bc(iad2)=bc(iad1)
             if(abs(bc(iad1)).gt.1d-10)write(6,*)('get '),bc(iad2),
     $            ('from '),iad1,iad1-i2e(is),j34,j12,i12,n12,i34
            end do
           end do
          end do
         end do
         write(6,*)('K part '),isblk(2,is),isblk(3,is),
     $        isblk(4,is),isblk(1,is),('d')
         call prntm2(bc(ibcoff),nn14,m23,nn14)
c     recall K_13^42 is (12|34) with ab distributed across procs.
c     but with 3<>4, ie K_14^32
          do i4=0,nvirt(isblk(2,is))-1
           i4p=i4+irefo(isblk(2,is))
           do i3=0,nvirt(isblk(3,is))-1
            i3p=i3+irefo(isblk(3,is))
            j34=i3+nvirt(isblk(3,is))*i4
            do i2=0,irefo(isblk(4,is))-1
             do i1=0,irefo(isblk(1,is))-1
              if(isblk(1,is).eq.isblk(2,is))then
               i12=((i4p*(i4p+1))/2)+i1
               i34=((i3p*(i3p+1))/2)+i2
              else
               i12=i1+nh0av(isblk(1,is))*i4p
               i34=i3p+nh0av(isblk(3,is))*i2
              end if
              j12=i1+irefo(isblk(1,is))*i2
              iad1=i2e(is)+i12+n12*i34
              iad2=ibcoff+j34+nn32*j12
              bc(iad2)=bc(iad1)
             if(abs(bc(iad1)).gt.1d-10)write(6,*)('get '),bc(iad2),
     $            ('from '),iad1,iad1-i2e(is),i12,n12,i34
             end do
            end do
           end do
          end do
          write(6,*)('K part '),isblk(1,is),isblk(4,is),
     $         isblk(3,is),isblk(2,is),('b')
          call prntm2(bc(ibcoff),nn32,m14,nn32)
c     recall K_13^42 is (12|34) with ab distributed across procs.
c     but with 1<>2 and 3<>4, ie K_24^31
          do i4=0,nvirt(isblk(1,is))-1
           i4p=i4+irefo(isblk(1,is))
           do i3=0,nvirt(isblk(3,is))-1
            i3p=i3+irefo(isblk(3,is))
            j43=i3+nvirt(isblk(3,is))*i4
            do i2=0,irefo(isblk(4,is))-1
             do i1=0,irefo(isblk(2,is))-1
              if(isblk(1,is).eq.isblk(2,is))then
               i12=((i4p*(i4p+1))/2)+i1
               i34=((i3p*(i3p+1))/2)+i2
              else
               i12=i4p+nh0av(isblk(1,is))*i1
               i34=i3p+nh0av(isblk(3,is))*i2
              end if
              j21=i1+irefo(isblk(2,is))*i2
              iad1=i2e(is)+i12+n12*i34
              iad2=ibcoff+j43+nn13*j21
              bc(iad2)=bc(iad1)
             if(abs(bc(iad1)).gt.1d-10)write(6,*)('get '),bc(iad2),
     $            ('from '),iad1,iad1-i2e(is)
             end do
            end do
           end do
          end do
          write(6,*)('K part '),isblk(2,is),isblk(4,is),
     $        isblk(3,is),isblk(1,is),('c')
          call prntm2(bc(ibcoff),nn13,m24,nn13)
       end do
      end if                                                            4d27s23
      ibcoff=ibcoffo
      return
      end
