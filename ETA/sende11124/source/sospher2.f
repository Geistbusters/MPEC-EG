c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine sospher2(iwfbra,ioffb,iwfket,ioffk,ijjs,ijjm,nspc,nec, 5d18s21
     $     multh,irefo,ih0a,i4so2,irel,ism,norb,mdon,esf,jjmmm,iwww,    8d30s21
     $     nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,irw2,ih0n,nh0,    11d15s21
     $     isopt,nsopt,i4or,ionexr,jmatsr,kmatsr,kmatsrb,i3xr,iifmx,    11d15s21
     $     ntype,npadddi,nbasp,nbaspc,natom,ngaus,ibdat,iapair,ibstor,  12d20s20
     $         isstor,isym,ascale,idorel,iorb,lprint,xmassp,bc,ibc,     11d21s22
     $     shift,n4vso)                                                 2d8s23
      implicit real*8 (a-h,o-z)                                         5d18s21
      logical ldiag,lprint,ldebug                                       4d1s24
      character*6 labb,labk                                             10d27s20
      integer*1 ipack1(4)                                               5d18s21
      integer*2 ipack2(4)                                               8d30s21
      integer*8 ioffb(*),ioffk(*),ijjs(*),ijjm(*),ipack8                8d30s21
      equivalence (ipack1,npack4)                                       5d18s21
      equivalence (ipack8,ipack2)                                       8d30s21
      dimension iwfbra(nspc,*),iwfket(nspc,*),multh(8,8),irefo(*),      5d18s21
     $     ih0a(8,2,8),i4so2(8,8,8,2),irel(*),ism(*),esf(*),itypex(4),  10d8s21
     $     ntype(4),nbasp(*),nbaspc(*),xmassp(*)                        8d16s22
      include "common.store"                                            5d18s21
      include "common.print"                                            3d29s24
      common/singcm/iuse,nff
      data icall/0/
      ldebug=iprtr(30).ne.0                                             3d29s24
      do i=1,6                                                          4d1s24
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
      if(lprint)                                                        3d2s22
     $ write(6,*)('in sospher2 for '),iwfbra(1,1),labb,('and '),        3d2s22
     $     iwfket(1,1),labk                                             6d2s21
      nrootb=iwfbra(3,1)                                                5d18s21
      nrootk=iwfket(3,1)                                                5d18s21
      llb2=llb*2                                                        5d18s21
      nllb=llb2+1                                                       5d18s21
      nsob=nllb*iwfbra(1,1)                                             5d18s21
      llk2=llk*2                                                        5d18s21
      nllk=llk2+1                                                       5d18s21
      nsok=nllk*iwfket(1,1)                                             5d18s21
      nroot2=nrootb*nrootk                                              5d18s21
      iqtop=min(i2sb+i2sk,llb2+llk2)                                    5d18s21
      iqbot=max(iabs(i2sb-i2sk),iabs(llb2-llk2))                        5d19s21
      nq=((iqtop-iqbot)/2)+1                                            5d17s21
      if(nq.le.0)then                                                   5d20s21
       ibcoff=ibcoffo                                                   5d20s21
       return                                                           5d20s21
      end if                                                            5d20s21
      nn=nsob*nsok                                                        5d17s21
      nf=nq                                                             5d18s21
      ihso=ibcoff                                                       5d14s21
      ibcoff=ihso+nsob*nsok*nroot2                                       5d20s21
      icleb=ibcoff                                                      5d17s21
      ibcoff=icleb+nn*nq                                                5d18s21
      call enough('sospher2.  1',bc,ibc)
      do i=ihso,ibcoff-1                                                5d14s21
       bc(i)=0d0                                                        5d14s21
      end do                                                            5d14s21
      iaout=ibcoff                                                      5d18s21
      ibcoff=iaout+nroot2                                               5d18s21
      nreal=0                                                           5d19s21
      nimag=0                                                           5d19s21
      nlind=0                                                           8d30s21
      ilind=ibcoff                                                      8d30s21
      idlin=ilind+nq*nq                                                 8d30s21
      itlin=idlin+nq                                                    8d30s21
      ihsolin=itlin+nq                                                  8d30s21
      ibcoff=ihsolin+nq*nroot2                                          8d30s21
      call enough('sospher2.  2',bc,ibc)
      do mlk2=0,llk2,2                                                  8d30s21
       do msk2=-i2sk,i2sk,2                                              8d30s21
        do id2=-iqtop,iqtop,2                                           8d30s21
         mlb2=mlk2+id2                                                  5d14s21
         msb2=msk2-id2                                                  5d14s21
         if(iabs(mlb2).le.llb2.and.iabs(msb2).le.i2sb)then              8d30s21
          do i=0,nq-1                                                   8d30s21
           bc(itlin+i)=0d0                                              8d30s21
          end do                                                        8d30s21
          do iq=max(iabs(id2),iqbot),iqtop,2                            8d30s21
           prodc=cleb2(llb2,mlb2,iq,-id2,llk2,mlk2)                     5d18s21
     $         *cleb2(i2sb,msb2,iq,+id2,i2sk,msk2)                      5d18s21
           id2h=id2/2                                                   5d18s21
           if(mod(id2h,2).ne.0)prodc=-prodc                             5d18s21
           icol=(iq-iqbot)/2                                            8d30s21
           if(ldebug)write(6,*)('cleb '),mlb2,mlk2,msb2,msk2,id2,       4d1s24
     $          prodc,icol                                              4d1s24
           bc(itlin+icol)=prodc                                         8d30s21
          end do                                                        8d30s21
          sz=0d0                                                        8d30s21
          do i=0,nq-1                                                   8d30s21
           sz=sz+bc(itlin+i)**2                                         8d30s21
          end do                                                        8d30s21
          sz=sqrt(sz/dfloat(nq))                                        8d30s21
          if(sz.gt.1d-10)then                                           8d30s21
           szt=0d0                                                      3d2s22
           do i=0,nq-1                                                  3d2s22
            szt=szt+bc(itlin+i)**2                                      3d2s22
           end do                                                       3d2s22
           szt=sqrt(szt/dfloat(nq))                                     3d2s22
           if(nlind.eq.0.and.szt.gt.1d-5)then                           3d2s22
            do i=0,nq-1                                                 8d30s21
             bc(ilind+i)=bc(itlin+i)                                    8d30s21
            end do                                                      8d30s21
            nlind=1                                                     8d30s21
            call genhsoa2(iwfbra,iwfket,nspc,mlb2,mlk2,msb2,msk2,
     $          bc(iaout),nec,multh,irefo,ih0a,i4so2,irel,ism,norb,mdon,5d19s21
     $          iri,itypex,nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,  9d20s21
     $           irw2,ih0n,nh0,isopt,nsopt,i4or,ionexr,jmatsr,kmatsr,   10d13s21
     $           kmatsrb,i3xr,iifmx,ntype,npadddi,nbasp,nbaspc,natom,   12d20s20
     $           ngaus,ibdat,iapair,ibstor,isstor,isym,ascale,idorel,   12d23s21
     $           iorb,0,0,0,bc,ibc,n4vso,1,0)                           3d27s24
            jhsolin=ihsolin                                             8d30s21
            do i=0,nroot2-1                                             8d30s21
             bc(jhsolin)=bc(iaout+i)                                    8d30s21
             jhsolin=jhsolin+nq                                         8d30s21
            end do                                                      8d30s21
           else                                                         8d30s21
            nlindp=nlind+1
            itrial=ibcoff                                               8d30s21
            ibcoff=itrial+nq*nlindp
            call enough('sospher2.  3',bc,ibc)
            do i=0,nq*nlind-1
             bc(itrial+i)=bc(ilind+i)                                   8d30s21
            end do                                                      8d30s21
            jtrial=itrial+nq*nlind
            do i=0,nq-1                                                 8d30s21
             bc(jtrial+i)=bc(itlin+i)
            end do
            ikeep=0                                                     8d30s21
            do i=0,nlindp-1                                             8d30s21
             iit=itrial+nq*i                                            8d30s21
             if(.not.ldebug)then                                        4d1s24
              do j=0,i-1
               jjt=itrial+nq*j                                           8d30s21
               dot=0d0                                                   8d30s21
               do k=0,nq-1                                               8d30s21
                dot=dot+bc(jjt+k)*bc(iit+k)                              8d30s21
               end do                                                    8d30s21
               do k=0,nq-1                                               8d30s21
                bc(iit+k)=bc(iit+k)-bc(jjt+k)*dot                        8d30s21
               end do                                                    8d30s21
              end do                                                     8d30s21
             end if                                                     4d1s24
             sz=0d0                                                     8d30s21
             do k=0,nq-1                                                8d30s21
              sz=sz+bc(iit+k)**2                                        8d30s21
             end do                                                     8d30s21
             if(sz.gt.1d-10)then                                        8d30s21
              ikeep=ikeep+1                                             8d30s21
              szi=1d0/sqrt(sz)                                          8d30s21
              do k=0,nq-1                                               8d30s21
               bc(iit+k)=bc(iit+k)*szi                                  8d30s21
              end do                                                    8d30s21
             end if                                                     8d30s21
            end do                                                      8d30s21
            if(ikeep.gt.nlind)then                                      8d30s21
             call genhsoa2(iwfbra,iwfket,nspc,mlb2,mlk2,msb2,msk2,
     $          bc(iaout),nec,multh,irefo,ih0a,i4so2,irel,ism,norb,mdon,5d19s21
     $          iri,itypex,nvirt,maxbx,maxbd,srh,sr2,nsymb,irw0,irw1,   9d20s21
     $            irw2,ih0n,nh0,isopt,nsopt,i4or,ionexr,jmatsr,kmatsr,  10d13s21
     $            kmatsrb,i3xr,iifmx,ntype,npadddi,nbasp,nbaspc,natom,  12d20s20
     $            ngaus,ibdat,iapair,ibstor,isstor,isym,ascale,idorel,  12d23s21
     $            iorb,0,0,0,bc,ibc,n4vso,0,0)                          3d27s24
             jhsolin=ihsolin+nlind                                      8d30s21
             do i=0,nroot2-1                                             8d30s21
              bc(jhsolin)=bc(iaout+i)                                    8d30s21
              jhsolin=jhsolin+nq                                         8d30s21
             end do                                                      8d30s21
             jlind=ilind+nq*nlind
             do i=0,nq-1                                                 8d30s21
              bc(jlind+i)=bc(itlin+i)
             end do
             nlind=nlind+1                                              8d30s21
            end if                                                      8d30s21
           end if                                                       8d30s21
           if(nlind.eq.nq.and..not.ldebug)then                          4d1s24
            go to 1066                                                  8d30s21
           end if                                                       8d30s21
          end if                                                        8d30s21
         end if                                                         8d30s21
        end do                                                          8d30s21
       end do                                                           8d30s21
      end do                                                            8d30s21
 1066 continue                                                          8d30s21
      iuse=0                                                            8d30s21
      icoef=ibcoff                                                      8d30s21
      iwgt=icoef+nq*nroot2                                              2d3s22
      iscr=iwgt+nq                                                      8d30s21
      ibcoff=iscr+3*nq*(nq+1)                                           8d30s21
      call enough('sospher2.  4',bc,ibc)
      do i=0,nq-1                                                       8d30s21
       bc(iwgt+i)=1d0                                                   8d30s21
      end do                                                            8d30s21
      do i=0,nq-1                                                       8d30s21
       do j=0,nq-1                                                      8d30s21
        ji=ilind+j+nq*i                                                 8d30s21
        ij=iscr+i+nq*j                                                  8d30s21
        bc(ij)=bc(ji)                                                   8d30s21
       end do                                                           8d30s21
      end do                                                            8d30s21
      do i=0,nq*nq-1                                                    8d30s21
       bc(ilind+i)=bc(iscr+i)                                           8d30s21
      end do                                                            8d30s21
      if(ldebug)then                                                    4d1s24
       iwo=6                                                            4d1s24
      else                                                              4d1s24
       iwo=0                                                            4d1s24
      end if                                                            4d1s24
      call lsqfit2(bc(ilind),nq,nq,bc(ihsolin),nq,nroot2,nq,bc(icoef),
     $     nq,bc(iscr),bc(iwgt),iwo,rms,bc,ibc)                         4d1s24
      if(lprint)then                                                    3d2s22
       if(rms.gt.1d-10)then                                             4d1s24
        write(6,*)('> bad fit in sospher2! '),rms                       4d1s24
       end if                                                           4d1s24
       iqq=iqbot                                                        3d2s22
       do iq=1,nq                                                        3d2s22
        iqh=iqq/2                                                       3d2s22
        if(iqh.eq.0.and.loc(iwfbra).eq.loc(iwfket))then                 11d21s22
         shiftu=shift                                                   11d21s22
        else                                                            11d21s22
         shiftu=0d0                                                     11d21s22
        end if                                                          11d21s22
        do ik=1,nrootk                                                  3d2s22
         do ib=1,nrootb                                                 3d2s22
          iad=icoef+iq-1+nq*(ib-1+nrootb*(ik-1))                        3d2s22
          if(ib.eq.ik)bc(iad)=bc(iad)+shiftu                            11d21s22
          if(abs(bc(iad)).gt.1d-10)then                                 3d2s22
           write(6,28)iqh,ib,iwfbra(1,1),labb,ik,iwfket(1,1),           3d2s22
     $           labk,bc(iad)                                           3d2s22
   28      format('A'i1,i3,1x,i1,a6,i3,1x,i1,a6,'>=',es14.6)            3d2s22
          end if                                                        3d2s22
         end do                                                         3d2s22
        end do                                                          3d2s22
        iqq=iqq+2                                                       3d2s22
       end do                                                            3d2s22
      end if                                                            3d2s22
      do msk2=-i2sk,i2sk,2                                                 5d14s21
       do mlk2=-llk2,llk2,2                                                    5d14s21
        do iq=iqbot,iqtop,2                                             5d18s21
         do id2=-iq,iq,2                                                5d18s21
          mlb2=mlk2+id2                                                  5d14s21
          msb2=msk2-id2                                                  5d14s21
          if(iabs(mlb2).le.llb2.and.iabs(msb2).le.i2sb)then             5d18s21
           prodc=cleb2(llb2,mlb2,iq,-id2,llk2,mlk2)                     5d18s21
     $         *cleb2(i2sb,msb2,iq,+id2,i2sk,msk2)                      5d18s21
           id2h=id2/2                                                   5d18s21
           if(mod(id2h,2).ne.0)prodc=-prodc                             5d18s21
           mlo=mlk2+llk2                                                  5d14s21
           mlo=mlo/2                                                     5d14s21
           mso=msk2+i2sk                                                  5d14s21
           mso=mso/2                                                     5d14s21
           iket=mlo+nllk*mso
           mlo=mlb2+llb2                                                  5d14s21
           mlo=mlo/2                                                     5d14s21
           mso=msb2+i2sb                                                  5d14s21
           mso=mso/2                                                     5d14s21
           ibra=mlo+nllb*mso
           jhso=ihso+ibra+nsob*iket                                       5d14s21
           jcleb=icleb+ibra+nsob*iket+nn*((iq-iqbot)/2)                 5d18s21
           bc(jcleb)=prodc                                              5d18s21
          end if                                                         5d14s21
         end do                                                         5d18s21
        end do                                                          5d14s21
       end do                                                           5d14s21
      end do                                                            5d14s21
      ijjk=ibcoff                                                        5d18s21
      ibcoff=ijjk+nsok*nsok                                             5d18s21
      jjmink=iabs(i2sk-llk2)                                             5d18s21
      jjmaxk=i2sk+llk2                                                     5d18s21
      numjk=((jjmaxk-jjmink)/2)+1                                          5d18s21
      ijjb=ibcoff                                                        5d18s21
      ibcoff=ijjb+nsob*nsob                                             5d18s21
      jjminb=iabs(i2sb-llb2)                                             5d18s21
      jjmaxb=i2sb+llb2                                                     5d18s21
      numjb=((jjmaxb-jjminb)/2)+1                                          5d18s21
      call enough('sospher2.  5',bc,ibc)
      jjmin=max(jjminb,jjmink)                                          5d18s21
      jjmax=min(jjmaxb,jjmaxk)                                          5d18s21
      if(jjmax.ge.jjmin)then                                            5d18s21
       numj=((jjmaxb-jjmin)/2)+1                                        5d18s21
       ieigjj=ibcoff                                                    5d18s21
       ibcoff=ieigjj+numj*nroot2                                        5d18s21
       call enough('sospher2.  6',bc,ibc)
       do i=ieigjj,ibcoff-1                                             5d18s21
        bc(i)=0d0                                                       5d18s21
       end do                                                           5d18s21
       jjj=ijjb
       nrowb=0                                                          5d18s21
       do jj=jjmin,jjmax,2                                              5d18s21
        do mjj=-jj,jj,2                                                  5d18s21
         nrowb=nrowb+1                                                  5d18s21
         do mss=-i2sb,i2sb,2                                               5d18s21
          do mll=-llb2,llb2,2                                               5d18s21
           bc(jjj)=cleb2(llb2,mll,i2sb,mss,jj,mjj)                         5d18s21
           jjj=jjj+1                                                     5d18s21
          end do                                                         5d18s21
         end do                                                          5d18s21
        end do                                                           5d18s21
       end do                                                            5d18s21
       jjj=ijjk
       ncolk=0                                                          5d18s21
       do jj=jjmin,jjmax,2                                               5d18s21
        do mjj=-jj,jj,2                                                  5d18s21
         ncolk=ncolk+1
         do mss=-i2sk,i2sk,2                                               5d18s21
          do mll=-llk2,llk2,2                                               5d18s21
           bc(jjj)=cleb2(llk2,mll,i2sk,mss,jj,mjj)                         5d18s21
           jjj=jjj+1                                                     5d18s21
          end do                                                         5d18s21
         end do                                                          5d18s21
        end do                                                           5d18s21
       end do                                                            5d18s21
       itmp1=ibcoff
       itmp2=itmp1+nn
       ibcoff=itmp2+nn
       call enough('sospher2.  7',bc,ibc)
       do iq=0,nf-1                                                      5d18s21
        jcleb=icleb+iq*nn                                                5d18s21
        call dgemm('n','n',nsob,ncolk,nsok,1d0,bc(jcleb),nsob,          5d18s21
     $       bc(ijjk),nsok,0d0,bc(itmp1),nsob,                           5d18s21
     d' sospher2.  1')
        do i=0,ncolk-1                                                    5d18s21
         do j=0,nsob-1                                                   5d18s21
          ji=itmp1+j+nsob*i                                              5d18s21
          ij=jcleb+i+ncolk*j                                              5d18s21
          bc(ij)=bc(ji)                                                 5d18s21
         end do                                                         5d18s21
        end do                                                          5d18s21
        call dgemm('n','n',ncolk,nrowb,nsob,1d0,bc(jcleb),ncolk,        5d18s21
     $       bc(ijjb),nsob,0d0,bc(itmp1),ncolk,                         5d18s21
     d' sospher2.  2')
        do i=0,nrowb-1                                                    5d18s21
         do j=0,ncolk-1                                                   5d18s21
          ji=itmp1+j+ncolk*i                                              5d18s21
          ij=jcleb+i+nrowb*j                                              5d18s21
          bc(ij)=bc(ji)                                                 5d18s21
         end do                                                         5d18s21
        end do                                                          5d18s21
        kcleb=jcleb                                                      5d18s21
        keigjj=ieigjj                                                    5d18s21
        iq2=iqbot+iq*2                                                  6d2s21
        do jj=jjmin,jjmax,2                                              5d18s21
         try9b=f9j(jj,0,jj,llk2,iq2,llb2,i2sk,iq2,i2sb,1)
         try6b=f6j(llb2,i2sb,jj,i2sk,llk2,iq2,1)
         try6b=try6b/sqrt(dfloat(jj+1)*dfloat(iq2+1))
         isum6=i2sb+llk2+jj+iq2
         if(mod(isum6,2).ne.0)then
          write(6,*)('6j phase error! '),isum6
         end if
         isum6=isum6/2
         if(mod(isum6,2).ne.0)try6b=-try6b
         tryb=try9b*sqrt(dfloat(jj+1)*dfloat(iq2+1)*dfloat(llk2+1)
     $        *dfloat(i2sk+1))
         isumb=i2sk-llk2+3*i2sb+llb2+iq2
         if(mod(isumb,2).ne.0)then
          write(6,*)('9j phase error! '),isumb
         end if
         isumb=isumb/2
         if(mod(isumb,2).ne.0)tryb=-tryb
         trya=f6j(llb2,i2sb,jj,i2sk,llk2,iq2,1)                          6d2s21
         try=trya*sqrt(dfloat(llk2+1)*dfloat(i2sk+1))                    6d2s21
         isum=jj+i2sk+llb2
         if(mod(isum,2).ne.0)then
          write(6,*)('6j phase error! '),jj,i2sk,isum
         end if
         try0=try
         isum=isum/2
         if(mod(isum,2).ne.0)try=-try
         do l=0,nroot2-1                                                5d19s21
          bc(keigjj+numj*l)=bc(keigjj+numj*l)                           5d19s21
     $         +bc(kcleb)*bc(icoef+iq+nf*l)                             5d19s21
         end do                                                         5d19s21
         keigjj=keigjj+1                                                 5d18s21
         mmj=jj+1                                                        5d18s21
         kcleb=kcleb+mmj*(nrowb+1)                                         5d18s21
        end do                                                           5d18s21
       end do                                                            5d18s21
       keigjj=ieigjj                                                     5d18s21
       do jj=jjmin,jjmax,2                                               5d18s21
        if(ldiag)then                                                   5d19s21
         imassp=1                                                       8d16s22
         do l=1,nrootb                                                  5d19s21
          lm=l-1                                                        5d19s21
          leigjj=keigjj+numj*(lm*(nrootb+1))                            5d19s21
          bc(leigjj)=bc(leigjj)+esf(l)                                  5d19s21
          do j=1,l-1                                                    8d16s22
           jm=j-1                                                       8d16s22
           leigjl=keigjj+numj*(jm+nrootb*lm)                            8d16s22
           leiglj=keigjj+numj*(lm+nrootb*jm)                            8d16s22
           bc(leigjl)=bc(leigjl)+xmassp(imassp)                         8d16s22
           bc(leiglj)=bc(leiglj)+xmassp(imassp)                         8d16s22
           imassp=imassp+1                                              8d16s22
          end do                                                        8d16s22
          bc(leigjj)=bc(leigjj)+xmassp(imassp)                          8d16s22
          imassp=imassp+1                                               8d16s22
         end do                                                         5d19s21
        end if                                                          5d19s21
        jju=((jj-jjmmm)/2)+1                                            5d18s21
        ndimj=ijjs(jju)                                                 5d18s21
        ijjl=ijjm(jju)+ndimj*ndimj+ioffk(jju)                           5d19s21
        do lk=0,nrootk-1                                                5d19s21
         ibc(ijjl+lk)=lk+1                                              5d19s21
         ibc(ijjl+lk+ndimj)=iwww                                        5d19s21
         do lb=0,nrootb-1                                               5d19s21
          leigjj=keigjj+numj*(lb+nrootb*lk)                             5d19s21
          ll=ijjm(jju)+ioffb(jju)+lb+ndimj*(ioffk(jju)+lk)               5d19s21
          bc(ll)=bc(leigjj)                                               5d18s21
          ll=ijjm(jju)+ioffk(jju)+lk+ndimj*(ioffb(jju)+lb)               5d19s21
          bc(ll)=bc(leigjj)                                             5d19s21
         end do                                                         5d19s21
        end do                                                          5d19s21
        ioffk(jju)=ioffk(jju)+nrootk                                    5d19s21
        ioffb(jju)=ioffb(jju)+nrootb                                    5d19s21
        keigjj=keigjj+1
       end do                                                            5d18s21
      end if                                                            5d18s21
      ibcoff=ibcoffo                                                    5d18s21
      return                                                            5d18s21
      end                                                               5d18s21
