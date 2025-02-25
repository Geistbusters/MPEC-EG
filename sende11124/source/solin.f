c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine solin(iwfbra,iwfket,nspc,nec,multh,irefo,ih0a,i4so2,   5d27s21
     $     irel,ism,norb,mdon,nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,    9d20s21
     $     irw1,irw2,ih0n,nh0,isopt,nsopt,i4or,ionexr,jmatsr,kmatsr,    10d13s21
     $     kmatsrb,i3xr,iifmx,ntype,npadddi,nbasp,nbaspc,natom,ngaus,   12d20s20
     $     ibdat,iapair,ibstor,isstor,isym,ascale,idorel,iorb,lprint,   11d9s22
     $     bc,ibc,shift,n4vso)                                          2d8s23
      implicit real*8 (a-h,o-z)                                         5d18s21
      logical ldiag,lprint,logb,logk,nlogk,ldebug                       3d29s24
      character*6 labb,labk                                             10d27s20
      integer*1 ipack1(4)                                               5d18s21
      equivalence (ipack1,npack4)                                       5d18s21
      dimension iwfbra(nspc,*),iwfket(nspc,*),multh(8,8),irefo(*),      5d18s21
     $     ih0a(8,2,8),i4so2(8,8,8,2),irel(*),ism(*),nvirt(*),nbasp(*), 12d20s20
     $     nbaspc(*)                                                    12d20s20
      include "common.store"                                            5d18s21
      include "common.print"                                            3d29s24
      common/singcm/iuse,nff
      data icall/0/
      save
      ldebug=iprtr(30).ne.0                                             3d29s24
c
c     if d2h, check on g/u symmetry                                     3d11s22
c                                                                       3d11s22
      if(nsymb.eq.8)then                                                3d11s22
       if(iwfbra(2,1).eq.2.or.iwfbra(2,1).eq.3.or.iwfbra(2,1).eq.5.     3d11s22
     $     or.iwfbra(2,1).eq.8)then                                     3d11s22
        kgub=-1                                                          3d11s22
       else                                                             3d11s22
        kgub=+1                                                          3d11s22
       end if                                                           3d11s22
       if(iwfket(2,1).eq.2.or.iwfket(2,1).eq.3.or.iwfket(2,1).eq.5.     3d11s22
     $     or.iwfket(2,1).eq.8)then                                     3d11s22
        kgu=-1                                                          3d11s22
       else                                                             3d11s22
        kgu=+1                                                          3d11s22
       end if                                                           3d11s22
       if(kgub.ne.kgu)then                                              3d11s22
        return                                                          3d11s22
       end if                                                           3d11s22
      end if                                                            3d11s22
      do i=1,6                                                          10d27s20
       if(iwfbra(13+i,1).ne.0)then                                         10d27s20
        labb(i:i)=char(iwfbra(13+i,1))                                     10d27s20
       else                                                             10d27s20
        labb(i:i)=' '                                                   10d27s20
       end if                                                           10d27s20
       if(iwfket(13+i,1).ne.0)then                                         10d27s20
        labk(i:i)=char(iwfket(13+i,1))                                     10d27s20
       else                                                             10d27s20
        labk(i:i)=' '                                                   10d27s20
       end if                                                           10d27s20
      end do                                                            10d27s20
      ibcoffo=ibcoff                                                    5d18s21
      npack4=iwfbra(6,1)                                                5d18s21
      llb=ipack1(3)                                                     5d18s21
      npack4=iwfket(6,1)                                                5d18s21
      ldiag=loc(iwfbra).eq.loc(iwfket)
      llk=ipack1(3)                                                     5d18s21
      i2sb=iwfbra(1,1)-1                                                5d18s21
      i2sk=iwfket(1,1)-1                                                5d18s21
      if(lprint)then                                                    3d8s22
       write(6,*)('in solin for '),iwfbra(1,1),labb,('and '),           3d11s22
     $     iwfket(1,1),labk                                             3d11s22
      end if                                                            3d8s22
      nrootb=iwfbra(3,1)                                                5d18s21
      nrootk=iwfket(3,1)                                                5d18s21
      llb2=llb*2                                                        5d18s21
c
c     by convention, only consider lambda ge 0 - rest is taken care of
c     by parity adaptation
c     furthermore, if lambda=0, only need sigma ge 0.                   10d20s21
c
       nllb=1                                                           5d27s21
      nsob=nllb*iwfbra(1,1)                                             5d18s21
      llk2=llk*2                                                        5d18s21
       nllk=1                                                           5d27s21
      nsok=nllk*iwfket(1,1)                                             5d18s21
      nroot2=nrootb*nrootk                                              5d18s21
      if(idorel.gt.0)then                                               3d14s22
       iqbot=2                                                          3d14s22
       iqtop=2                                                          3d14s22
      else                                                              3d14s22
       iqbot=iabs(i2sb-i2sk)                                             3d10s22
       iqtop=min(4,i2sb+i2sk)                                           3d15s22
      end if                                                            3d14s22
      nq=((iqtop-iqbot)/2)+1                                            5d17s21
      if(nq.le.0)then                                                   5d20s21
       ibcoff=ibcoffo                                                   5d20s21
       return                                                           5d20s21
      end if                                                            5d20s21
      nn=nsob*nsok                                                        5d17s21
      nf=0                                                              5d27s21
c
c     we are going to fit different msk and msb combinations.
c     the maximum possible number is range of id2: 2qtop+1=iqtop+1
c     because factor of 2 is already in iqtop times the number of
c     msks: (i2sk+1-mskbot)/2
c
       mskbot=-i2sk                                                     10d20s21
      itry2=i2sk+2-mskbot
      itry3=itry2/2
      nf=(iqtop+1)*itry3                                                2d23s24
      ihso=ibcoff                                                       5d14s21
      ibcoff=ihso+nsob*nsok*nroot2                                       5d20s21
      icleb=ibcoff                                                      5d17s21
      ibcoff=icleb+nn*nf                                                5d27s21
      do i=ihso,ibcoff-1                                                5d14s21
       bc(i)=0d0                                                        5d14s21
      end do                                                            5d14s21
      nlind=0                                                           9d5s21
      nreal=0                                                           5d19s21
      nimag=0                                                           5d19s21
      ilind=ibcoff                                                      8d30s21
      idlin=ilind+nf*nf                                                 8d30s21
      itlin=idlin+nf                                                    8d30s21
      nrev=0                                                            3d6s24
      if(llb.eq.1.and.llk.eq.1.and.idorel.lt.0)then                     3d29s24
       if(lprint)                                                       3d29s24
     $     write(6,*)('we have two Pi states, so do reversed as well')  3d29s24
       nrev=1                                                           3d6s24
      end if                                                            3d6s24
      ihsolin=itlin+nf                                                  8d30s21
      ibcoff=ihsolin+nf*nroot2                                          3d6s24
      iaout=ibcoff                                                      5d18s21
      ibcoff=iaout+nroot2                                               3d6s24
      call enough('solin.  1',bc,ibc)
      do irev=0,nrev                                                    3d6s24
       if(ldebug)write(6,*)('for irev = '),irev                         3d29s24
        msbbot=-i2sb                                                     10d20s21
       mlk2=llk2
       numq=(iqtop+2-iqbot)/2                                            2d23s24
       itrial=ibcoff                                                     2d23s24
       ibcoff=itrial+nf*numq                                             2d23s24
       call enough('solin.itrial',bc,ibc)
       do iz=itrial,ibcoff-1                                             2d23s24
        bc(iz)=0d0                                                       2d23s24
       end do                                                            2d23s24
       jtrial=itrial                                                     2d23s24
       do iq=iqbot,iqtop,2                                               2d23s24
        do id2=-iq,iq,2                                                  2d23s24
         do msk2=mskbot,i2sk,2                                           2d23s24
          msb2=msk2-id2                                                  2d23s24
          mlb2=mlk2+id2                                                  2d23s24
          if(irev.eq.1)mlb2=-mlb2                                       3d6s24
          if(mlb2.eq.llb2.and.msb2.ge.msbbot.and.msb2.le.i2sb)then       2d23s24
c                                                                       3d6s24
c     if sigma-sigma, iq=2 ie A1 doesnt exist ... at least for          3d14s24
c     coulomb-spin orbit for there is no totally symmetry so operator.  3d6s24
c                                                                       3d6s24
           if(.not.                                                     3d6s24
     $      (iwfbra(2,1).eq.iwfket(2,1).and.llk2.eq.0.and.iq.eq.2.and.  3d22s24
     $         llb2.eq.0))then                                          3d22s24
            prodc=cleb2(i2sb,msb2,iq,+id2,i2sk,msk2)                      2d23s24
            icol=((iqtop+id2)/2)+(iqtop+1)*((msk2-mskbot)/2)              2d23s24
            if(ldebug)write(6,*)('iq '),iq,('spins '),msb2,id2,msk2,    3d29s24
     $           prodc,icol                                             3d29s24
            bc(jtrial+icol)=prodc                                         2d23s24
           end if
          end if                                                         2d23s24
         end do                                                          2d23s24
        end do                                                           2d23s24
        jtrial=jtrial+nf                                                 2d23s24
       end do
       itrial2=ibcoff                                                    2d23s24
       icol=itrial2+nf*numq                                              2d23s24
       ibcoff=icol+nf                                                    2d23s24
       call enough('solin.trial2',bc,ibc)
       do i=0,nf*numq-1                                                  2d23s24
        bc(itrial2+i)=bc(itrial+i)
       end do
       if(ldebug)then                                                   3d29s24
        write(6,*)('itrial2 '),itrial2                                  3d29s24
        call prntm2(bc(itrial2),nf,numq,nf)                             3d29s24
       end if                                                           3d29s24
       ikeep=0
       do i=0,nf-1                                                       2d23s24
        if(.not.ldebug)then                                             3d29s24
         do j=0,ikeep-1                                                   2d23s24
          dot=0d0                                                         2d23s24
          itrial2i=itrial2+i                                               2d23s24
          itrial2j=itrial2+j                                              2d23s24
          do k=0,numq-1                                                   2d23s24
           dot=dot+bc(itrial2i)*bc(itrial2j)                              2d23s24
           itrial2i=itrial2i+nf                                           2d23s24
           itrial2j=itrial2j+nf                                           2d23s24
          end do                                                          2d23s24
          itrial2i=itrial2+i                                               2d23s24
          itrial2j=itrial2+j                                              2d23s24
          do k=0,numq-1                                                   2d23s24
           bc(itrial2i)=bc(itrial2i)-dot*bc(itrial2j)                     2d23s24
           itrial2i=itrial2i+nf                                           2d23s24
           itrial2j=itrial2j+nf                                           2d23s24
          end do                                                          2d23s24
         end do                                                           2d23s24
        end if                                                          3d29s24
        dot=0d0                                                          2d23s24
        itrial2i=itrial2+i                                               2d23s24
        do k=0,numq-1                                                    2d23s24
         dot=dot+bc(itrial2i)**2                                         2d23s24
         itrial2i=itrial2i+nf                                            2d23s24
        end do                                                           2d23s24
        if(dot.gt.1d-10)then                                             2d23s24
         if(ldebug)write(6,*)('keep fcn '),i                            3d29s24
         dot=1d0/sqrt(dot)                                               2d23s24
         itrial2i=itrial2+i                                               2d23s24
         do k=0,numq-1                                                    2d23s24
          bc(itrial2i)=bc(itrial2i)*dot                                  2d23s24
          itrial2i=itrial2i+nf                                            2d23s24
         end do                                                           2d23s24
         ibc(icol+ikeep)=i                                               2d23s24
         if(ikeep.ne.i)then                                              2d23s24
          itrial2i=itrial2+i                                             2d23s24
          itrial2k=itrial2+ikeep                                         2d23s24
          do k=0,numq-1                                                  2d23s24
           bc(itrial2k)=bc(itrial2i)                                     2d23s24
           bc(itrial2i)=0d0                                              2d23s24
           itrial2k=itrial2k+nf                                          2d23s24
           itrial2i=itrial2i+nf                                          2d23s24
          end do                                                         2d23s24
         end if                                                          2d23s24
         ikeep=ikeep+1                                                   2d23s24
        end if                                                           2d23s24
       end do
       ifkeep=ibcoff                                                     2d23s24
       ibcoff=ifkeep+ikeep                                               2d23s24
       call enough('solin.ifkeep',bc,ibc)
       nfkeep=0                                                          2d23s24
       do i=0,numq-1                                                     2d23s24
        sz=0d0                                                           2d23s24
        imaybsav=itrial2+nf*i                                            2d26s24
        do j=0,ikeep-1                                                   2d23s24
         sz=sz+bc(imaybsav+j)**2                                         2d23s24
        end do                                                           2d23s24
        if(sz.gt.1d-10)then                                              2d23s24
         iq=iqbot+i*2                                                    2d23s24
         ibc(ifkeep+nfkeep)=i                                            2d23s24
         if(nfkeep.ne.i)then                                             2d23s24
          isav=itrial2+nf*nfkeep                                         2d26s24
          do j=0,ikeep-1                                                 2d23s24
           bc(isav+j)=bc(imaybsav+j)                                     2d23s24
          end do                                                         2d23s24
         end if                                                          2d23s24
         nfkeep=nfkeep+1                                                 2d23s24
        end if                                                           2d23s24
       end do                                                            2d23s24
       if(.not.ldebug)then                                              3d29s24
        nfkeep=min(nfkeep,ikeep)                                          2d23s24
        ikeep=min(ikeep,nfkeep)                                           2d26s24
       end if                                                           3d29s24
       do i=0,nfkeep-1                                                   2d26s24
        ii=ibc(ifkeep+i)
        ito=itrial2+ikeep*i                                              2d26s24
        ifrm=itrial+nf*ii                                                2d26s24
        do j=0,ikeep-1                                                   2d26s24
         jj=ibc(icol+j)
         bc(ito+j)=bc(ifrm+jj)                                           2d26s24
        end do                                                           2d26s24
       end do                                                            2d26s24
       do i=0,nfkeep-1                                                   2d26s24
        ii=ibc(icol+i)                                                   2d26s24
        ito=itrial+ikeep*i                                                  2d26s24
        ifrm=itrial2+ikeep*i                                                2d26s24
        do j=0,ikeep-1                                                      2d26s24
         bc(ito+j)=bc(ifrm+j)                                            2d26s24
        end do                                                           2d26s24
       end do                                                            2d26s24
       do jkeep=0,ikeep-1                                                2d23s24
        n2=ibc(icol+jkeep)/(iqtop+1)                                     2d23s24
        n1=ibc(icol+jkeep)-n2*(iqtop+1)                                  2d23s24
        id2=n1*2-iqtop                                                   2d23s24
        msk2=n2*2+mskbot                                                 2d23s24
        msb2=msk2-id2                                                    2d23s24
        if(ldebug)write(6,*)('reconstituted spins: '),msb2,id2,msk2     3d29s24
        mlb2=mlk2+id2                                                    2d23s24
        call genhsoa2(iwfbra,iwfket,nspc,mlb2,mlk2,msb2,msk2,
     $          bc(iaout),nec,multh,irefo,ih0a,i4so2,irel,ism,norb,mdon,5d19s21
     $          iri,idum,nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,    9d20s21
     $           irw2,ih0n,nh0,isopt,nsopt,i4or,ionexr,jmatsr,kmatsr,   10d13s21
     $           kmatsrb,i3xr,iifmx,ntype,npadddi,nbasp,nbaspc,natom,   12d20s20
     $           ngaus,ibdat,iapair,ibstor,isstor,isym,ascale,idorel,   12d23s21
     $           iorb,0,0,0,bc,ibc,n4vso,1,irev)                        3d29s24
         if(ldebug)then                                                 3d29s24
          write(6,*)('back from genhsoa2 '),loop
          call prntm2(bc(iaout),nrootb,nrootk,nrootb)
         end if                                                         3d29s24
        if(iri.ne.0)then
         write(6,*)('on return from genhsoa2, iri = '),iri
        end if
        jhsolin=ihsolin+jkeep                                            2d23s24
        do i=0,nroot2-1                                                 3d6s24
         bc(jhsolin)=bc(iaout+i)                                         2d23s24
         jhsolin=jhsolin+ikeep                                           2d23s24
        end do                                                           2d23s24
        do i=0,nfkeep-1                                                  2d26s24
         bc(itrial2+jkeep+ikeep*i)=bc(itrial+jkeep+ikeep*i)              2d26s24
        end do                                                           2d23s24
       end do                                                            2d23s24
       if(min(nfkeep,ikeep).gt.0)then                                    2d23s24
        iuse=0
        icoef=ibcoff                                                     2d23s24
        iwgt=icoef+nfkeep*nroot2                                        3d6s24
        iscr=iwgt+ikeep                                                  2d23s24
        ibcoff=iscr+2*nfkeep+nfkeep*nfkeep+2*nfkeep*ikeep+ikeep          2d23s24
        call enough('solin.icoeftrial2',bc,ibc)
        do iz=iwgt,iscr-1                                                2d23s24
         bc(iz)=1d0                                                      2d23s24
        end do                                                           2d23s24
        if(ldebug)then                                                  3d29s24
         iuu=6                                                          3d29s24
        else
         iuu=0                                                          3d29s24
        end if                                                          3d29s24
        call lsqfit2(bc(itrial2),ikeep,nfkeep,bc(ihsolin),ikeep,         3d6s24
     $       nroot2,                                                    3d6s24
     $       ikeep,bc(icoef),nfkeep,bc(iscr),bc(iwgt),iuu,rms,bc,ibc)   3d29s24
        if(lprint)then                                                    3d2s22
         if(rms.gt.1d-10)then                                             2d23s24
          write(6,*)('> bad new fit in solin! '),rms                         2d23s24
         end if
         ihit=0                                                          2d23s24
         do icc=0,nfkeep-1                                               2d23s24
          iqq=ibc(ifkeep+icc)*2+iqbot                                     2d23s24
          iqh=iqq/2                                                       3d2s22
          if(iqh.eq.0.and.loc(iwfbra).eq.loc(iwfket))then                 11d21s22
           shiftu=shift                                                   11d21s22
          else                                                            11d21s22
           shiftu=0d0                                                     11d21s22
          end if                                                          11d21s22
          do ik=1,nrootk                                                  3d2s22
           do ib=1,nrootb                                                 3d2s22
            iad=icoef+icc+nfkeep*(ib-1+nrootb*(ik-1))                   3d29s24
            if(ik.eq.ib)bc(iad)=bc(iad)+shiftu                            11d21s22
            if(abs(bc(iad)).gt.1d-14)then                                 2d20s24
             ihit=1                                                      2d23s24
             if(irev.eq.0)then                                          3d29s24
              write(6,28)iqh,ib,iwfbra(1,1),labb,ik,iwfket(1,1),           3d2s22
     $           labk,bc(iad)
   28         format('A'i1,i3,1x,i1,a6,i3,1x,i1,a6,'>=',2es24.16)          3d6s24
             else                                                       3d29s24
              write(6,29)iqh,ib,iwfbra(1,1),labb,ik,iwfket(1,1),           3d2s22
     $           labk,bc(iad)
   29         format('A'i1,'2',i3,1x,i1,a6,i3,1x,i1,a6,'>=',2es24.16)   3d29s24
             end if                                                     3d29s24
            end if                                                        3d2s22
           end do                                                         3d2s22
          end do                                                          3d2s22
         end do                                                            3d2s22
        end if                                                            3d2s22
       end if                                                            2d23s24
       ibcoff=itrial                                                     2d23s24
      end do                                                            3d6s24
      ibcoff=ibcoffo                                                    5d18s21
      return                                                            5d18s21
      end                                                               5d18s21
