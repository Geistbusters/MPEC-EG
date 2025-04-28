c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine gtmom(iwavebra,iwaveket,ixw1,ixw2,ncsf,nec,nspc,ixmt,  5d24s21
     $     opdata,iopdata,iosym,nop,opname,nbasdws,nsymb,mdon,mdoop,    5d24s21
     $     ism,irel,irefo,norb,multh,idoubo,itmom,nvirt,maxbx,maxbxd,   8d23s21
     $     srh,sr2,npadddi,lprint,natom,bc,ibc)                         11d10s22
      implicit real*8 (a-h,o-z)                                         5d12s21
      integer*1 ipack1(4)                                               5d12s21
      integer*4 ipack4(2)                                               5d12s21
      integer*8 ipack8,itmom                                            5d25s21
      character*10 opname(*)                                            5d12s21
      character*2 c2                                                    3d8s22
      character*4 ctype                                                 3d30s22
      character*1 crori(2)
      character*6 wnameb,wnamek                                         5d24s21
      logical lri,lprint,lok                                            3d28s22
      equivalence (ipack1,npack4)                                       5d12s21
      equivalence (ipack4,ipack8)                                       5d12s21
      dimension iwavebra(nspc,*),iwaveket(nspc,*),ncsf(*),ixmt(8,*),    8d3s23
     $     opdata(*),iopdata(7,*),iosym(*),nbasdws(*),ism(*),irel(*),   5d27s21
     $     irefo(*),multh(8,8),idoubo(*),itmom(4),nvirt(*),ipmv(2),     1d5s23
     $     idum(2)                                                      1d5s23
      include "common.store"                                            5d12s21
      include "common.basis"                                            3d15s22
      common/singcm/iuse,nff
      data crori/' ','i'/
      data idum/2*1/                                                    12d21s22
      dum=1d0                                                           1d27s23
      ibcoffo=ibcoff                                                    5d26s22
      srh=sqrt(0.5d0)                                                   5d24s21
      npack4=iwavebra(6,1)                                              5d24s21
      wnameb(2:2)=' '                                                   6d3s21
      do i=1,6                                                          5d20s21
       if(iwavebra(13+i,1).ne.0)then                                    5d24s21
        wnameb(i:i)=char(iwavebra(13+i,1))                              5d24s21
        nnnb=i                                                           5d13s21
       else                                                             3d2s22
        wnameb(i:i)=' '                                                 3d2s22
       end if                                                           5d13s21
      end do                                                            5d13s21
      llb=ipack1(3)                                                     5d24s21
      do i=2,2*llb+1                                                      5d24s21
       npack4=iwavebra(6,i)                                             5d24s21
      end do                                                            5d24s21
      npack4=iwaveket(6,1)                                              5d24s21
      wnamek(2:2)=' '                                                   6d3s21
      do i=1,6                                                          5d20s21
       if(iwaveket(13+i,1).ne.0)then                                    5d24s21
        wnamek(i:i)=char(iwaveket(13+i,1))                              5d24s21
        nnnk=i                                                           5d13s21
       else                                                             3d2s22
        wnamek(i:i)=' '                                                 3d2s22
       end if                                                           5d13s21
      end do                                                            5d13s21
      llk=ipack1(3)                                                     5d24s21
      if(lprint)write(6,*)('Hi, my name is gtmom for '),iwavebra(1,1),  3d31s22
     $     wnameb(1:nnnb),(' and'),iwaveket(1,1),wnamek(1:nnnk)         3d31s22
      nlzz=ipack1(2)                                                    5d26s21
      if(nlzz.eq.6)then                                                 5d26s21
       do i=2,2*llk+1                                                      5d24s21
        npack4=iwaveket(6,i)                                             5d24s21
       end do                                                            5d24s21
      else if(nlzz.eq.2.and.ipack1(3).ne.0)then                         5d26s21
       npack4=iwaveket(6,2)
      end if                                                            5d26s21
      if(nlzz.eq.6)then                                                 5d26s21
       if(lprint)write(6,*)('spherical case ')                          3d2s22
       iordb=0                                                           5d20s21
       if(iwavebra(15,1).ne.0)iordb=1                                      5d20s21
       if(mod(llb,2).ne.0)iordb=iordb+1                                     5d19s21
       iordb=mod(iordb,2)                                                  5d19s21
       mdoopb=((nec-iwavebra(1,1)+1)/2)+1                               5d24s21
       iordk=0                                                           5d20s21
       if(iwaveket(15,1).ne.0)iordk=1                                      5d20s21
       if(mod(llk,2).ne.0)iordk=iordk+1                                     5d19s21
       iordk=mod(iordk,2)                                                  5d19s21
       mdoopk=((nec-iwaveket(1,1)+1)/2)+1                               5d24s21
       nllb=2*llb+1                                                     5d24s21
       nllk=2*llk+1                                                     5d24s21
       nrootb=iwavebra(3,1)                                             5d24s21
       nrootk=iwaveket(3,1)                                             5d24s21
       nbk=nrootb*nrootk                                                5d24s21
       if(wnameb(2:2).eq.wnamek(2:2))then                               6d3s21
        ipbot=2                                                         6d3s21
        iptop=3                                                         6d3s21
        itmom(1)=0                                                      6d4s21
        itmom(2)=ibcoff                                                 6d4s21
        itmom(3)=itmom(2)+nbk                                            5d25s21
        irms=itmom(3)+nbk                                                5d25s21
       else                                                             6d3s21
        ipbot=1                                                         6d3s21
        iptop=1                                                         6d3s21
        itmom(1)=ibcoff                                                 6d4s21
        itmom(2)=0                                                      6d4s21
        itmom(3)=0                                                      6d4s21
        irms=itmom(1)+nbk                                               6d4s21
       end if                                                           6d3s21
       ifull=irms+nbk                                                   5d25s21
       icleb=ifull+nrootb*nrootk*nllb*nllk                              5d24s21
       ibcoff=icleb+nllb*nllk                                           5d24s21
       call enough('gtmom.  1',bc,ibc)
       llb2=llb*2                                                       5d24s21
       llk2=llk*2                                                       5d24s21
       do ipass=ipbot,iptop                                             6d3s21
        iheader=0                                                       1d6s23
        do i=0,nbk-1                                                    5d25s21
         bc(itmom(ipass)+i)=0d0                                         5d25s21
         bc(irms+i)=0d0                                                 5d25s21
        end do                                                          5d25s21
        if(ipass.eq.1)then                                              5d24s21
         iq=1                                                           5d24s21
c        +1                 -1
c     (-x-iy)/sqrt(2), (+x-iy)/sqrt(2)
c
c     i.e. real +1=-x/sqrt(2), imag +1=-y/sqrt(2)                       5d24s21
c          real -1=+x/sqrt(2), imag -1= same
         do i=1,nop                                                     5d24s21
          if(opname(i)(1:3).eq.'mux')mur=i
          if(opname(i)(1:3).eq.'muy')mui=i                                 5d24s21
         end do
         ffr=-srh                                                       5d24s21
         ffi=-srh                                                       5d24s21
        else if(ipass.eq.2)then                                         5d24s21
         iq=1                                                           5d24s21
c        +1                    -1
c     (-lx-ily)/sqrt(2),(lx-ily)/sqrt(2)
c     =(+ilx(i)-ly(i))/sqrt(2),(-ilx(i)-ly(i))/sqrt(2)
c
c     i.e. real +1 = -ly(i)/sqrt(2), imag +1=lx(i)/sqrt(2)
c          real -1 = same          , imag -1=-lx(i)/sqrt(2)
         do i=1,nop                                                     5d24s21
          if(opname(i)(1:4).eq.'lx/i')mui=i                                5d24s21
          if(opname(i)(1:4).eq.'ly/i')mur=i                                5d24s21
         end do
         ffr=-srh                                                       5d24s21
         ffi=srh                                                        5d24s21
        else                                                            5d24s21
         iq=2                                                           5d24s21
c        +1     -1
c     -zx-izy,zx-izy
c
c     i.e. real +1 = -zx, imag +1 = -zy
c          real -1 = +zx, imag -1 = same
         do i=1,nop                                                     5d24s21
          if(opname(i)(1:5).eq.'Re Q1')mur=i                               5d24s21
          if(opname(i)(1:5).eq.'Im Q1')mui=i                               5d24s21
         end do
         ffr=1d0                                                        5d26s21
         ffi=1d0                                                        5d26s21
        end if                                                          5d24s21
        iqlow=iabs(llb-llk)                                             5d25s21
        iqhig=llb+llk                                                   5d25s21
        if(iq.ge.iqlow.and.iq.le.iqhig)then                             5d25s21
         iq2=iq*2                                                        5d24s21
         do ipm=-1,1,2                                                   5d25s21
          nimag=0                                                         5d24s21
          do i=0,(nbk+1)*nllb*nllk-1                                       5d24s21
           bc(ifull+i)=0d0                                                 5d24s21
          end do                                                           5d24s21
          do mlb=-llb,llb                                                  5d24s21
           mlba=iabs(mlb)                                                 5d24s21
           do i=1,nllb                                                    5d24s21
            npack4=iwavebra(6,i)                                          5d24s21
            if(ipack1(4).eq.mlba)ibr=i                                     5d24s21
            if(ipack1(4).eq.-mlba)ibi=i                                     5d24s21
           end do                                                         5d24s21
           pbr=1d0                                                           5d14s21
           pbi=1d0                                                           5d14s21
           if(mlb.lt.0)then                                                 5d14s21
            if(mod(mlba+iordb,2).eq.0)then                                   5d20s21
             pfb=1d0                                                          5d14s21
            else                                                              5d14s21
             pfb=-1d0                                                         5d14s21
            end if                                                           5d14s21
            pbr=pfb                                                          5d14s21
            pbi=-pfb                                                         5d14s21
           end if                                                            5d14s21
           ibra=mlb+llb                                                   5d24s21
           mlb2=mlb*2                                                     5d24s21
           mlk=mlb-ipm                                                   5d25s21
           mlka=iabs(mlk)                                                5d24s21
           if(mlka.le.llk)then                                           5d24s21
c                      ipm=-1   ipm=+1
c                       r  i     r  i
c     electric dipole:  -  +     +  +
c     magnetic dipole:  +  -     +  +
c     electric quad:    -  +     +  +
c
            ffru=ffr                                                     5d24s21
            ffiu=ffi                                                     5d24s21
            if(ipm.lt.0)then                                             5d24s21
             if(ipass.eq.2)then                                          5d24s21
              ffiu=-ffiu                                                  5d24s21
             else                                                        5d24s21
              ffru=-ffru                                                 5d24s21
             end if                                                      5d24s21
            end if                                                       5d24s21
            iket=mlk+llk                                                 5d24s21
            jcleb=icleb+ibra+nllb*iket                                   5d24s21
            mlk2=mlk*2                                                   5d24s21
            ipm2=ipm*2                                                   5d24s21
            bc(jcleb)=cleb2(llb2,mlb2,iq2,-ipm2,llk2,mlk2)               5d25s21
            do i=1,nllk                                                     5d24s21
             npack4=iwaveket(6,i)                                          5d24s21
             if(ipack1(4).eq.mlka)ikr=i                                     5d24s21
             if(ipack1(4).eq.-mlka)iki=i                                     5d24s21
            end do                                                         5d24s21
            pkr=1d0                                                           5d14s21
            pki=1d0                                                           5d14s21
            if(mlk.lt.0)then                                                 5d14s21
             if(mod(mlka+iordk,2).eq.0)then                                   5d20s21
              pfk=1d0                                                          5d14s21
             else                                                              5d14s21
              pfk=-1d0                                                         5d14s21
             end if                                                           5d14s21
             pkr=pfk                                                          5d14s21
             pki=-pfk                                                         5d14s21
            end if                                                            5d14s21
             itmp=ibcoff                                                8d23s21
             ibcoff=itmp+nrootb*nrootk                                  8d23s21
             call enough('gtmom.  2',bc,ibc)
c     rr
            if((ibr.ne.ibi.or.iordb.eq.0).and.                           5d24s21
     $          (ikr.ne.iki.or.iordk.eq.0))then                         5d24s21
             isymc=multh(iwavebra(2,ibr),iwaveket(2,ikr))                5d24s21
             if(isymc.eq.iosym(mur))then                                 5d24s21
              ff=ffru*pbr*pkr                                            5d24s21
              if(mlk.ne.0)ff=ff*srh                                      5d24s21
              if(mlb.ne.0)ff=ff*srh                                      5d24s21
              call psioppsi(iwavebra(1,ibr),nrootb,mur,1d0,ixmt,iosym,0,8d23s21
     $             idum,dum,iwaveket(1,ikr),nrootk,bc(itmp),mdon,mdoop, 8d23s21
     $             ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,nvirt,  8d23s21
     $             maxbx,maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,npadddi,11d10s22
     $             bc,ibc)                                              11d10s22
              jfull=ifull+nbk*(ibra+nllb*iket)                           5d24s21
              do i=0,nrootk-1                                           8d23s21
               do j=0,nrootb-1                                          8d23s21
                ji=itmp+j+nrootb*i                                      8d23s21
                ij=jfull+i+nrootk*j                                     8d23s21
                bc(ij)=bc(ij)+bc(ji)*ff                                 8d23s21
               end do                                                   8d23s21
              end do                                                    8d23s21
             else if(isymc.eq.iosym(mui))then                            5d24s21
              write(6,*)('we have rir=i')                                5d24s21
              nimag=nimag+1
             end if                                                      5d24s21
            end if                                                       5d24s21
c     ri
            if((ibr.ne.ibi.or.iordb.eq.0).and.                           5d24s21
     $          (ikr.ne.iki.or.iordk.eq.1))then                         5d24s21
             isymc=multh(iwavebra(2,ibr),iwaveket(2,iki))                5d24s21
             if(isymc.eq.iosym(mur))then                                 5d24s21
              if(lprint)write(6,*)('we have rri=i')                     3d2s22
              nimag=nimag+1
             else if(isymc.eq.iosym(mui))then                            5d24s21
              ff=-pbr*pki*ffiu                                           5d24s21
              if(mlk.ne.0)ff=ff*srh                                      5d24s21
              if(mlb.ne.0)ff=ff*srh                                      5d24s21
              call psioppsi(iwavebra(1,ibr),nrootb,mui,1d0,ixmt,iosym,0,8d23s21
     $             idum,dum,iwaveket(1,iki),nrootk,bc(itmp),mdon,mdoop, 8d23s21
     $             ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,nvirt,  8d23s21
     $             maxbx,maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,npadddi,11d10s22
     $             bc,ibc)                                              11d10s22
              jfull=ifull+nbk*(ibra+nllb*iket)                           5d24s21
              do i=0,nrootk-1                                           8d23s21
               do j=0,nrootb-1                                          8d23s21
                ji=itmp+j+nrootb*i                                      8d23s21
                ij=jfull+i+nrootk*j                                     8d23s21
                bc(ij)=bc(ij)+bc(ji)*ff                                 8d23s21
               end do                                                   8d23s21
              end do                                                    8d23s21
             end if                                                      5d24s21
            end if                                                       5d24s21
c     ir
            if((ibr.ne.ibi.or.iordb.eq.1).and.                           5d24s21
     $          (ikr.ne.iki.or.iordk.eq.0))then                         5d24s21
             isymc=multh(iwavebra(2,ibi),iwaveket(2,ikr))                5d24s21
             if(isymc.eq.iosym(mur))then                                 5d24s21
              if(lprint)write(6,*)('we have irr=i')                     3d2s22
              nimag=nimag+1
             else if(isymc.eq.iosym(mui))then                            5d24s21
              ff=pbi*pkr*ffiu                                            5d24s21
              if(mlk.ne.0)ff=ff*srh                                      5d24s21
              if(mlb.ne.0)ff=ff*srh                                      5d24s21
              call psioppsi(iwavebra(1,ibi),nrootb,mui,1d0,ixmt,iosym,0,8d23s21
     $             idum,dum,iwaveket(1,ikr),nrootk,bc(itmp),mdon,mdoop, 8d23s21
     $             ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,nvirt,  8d23s21
     $             maxbx,maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,npadddi,11d10s22
     $             bc,ibc)                                              11d10s22
              jfull=ifull+nbk*(ibra+nllb*iket)                           5d24s21
              do i=0,nrootk-1                                           8d23s21
               do j=0,nrootb-1                                          8d23s21
                ji=itmp+j+nrootb*i                                      8d23s21
                ij=jfull+i+nrootk*j                                     8d23s21
                bc(ij)=bc(ij)+bc(ji)*ff                                 8d23s21
               end do                                                   8d23s21
              end do                                                    8d23s21
             end if                                                      5d24s21
            end if                                                       5d24s21
c     ii
            if((ibr.ne.ibi.or.iordb.eq.1).and.                           5d24s21
     $          (ikr.ne.iki.or.iordk.eq.1))then                         5d24s21
             isymc=multh(iwavebra(2,ibi),iwaveket(2,iki))                5d24s21
             if(isymc.eq.iosym(mur))then                                 5d24s21
              ff=pbi*pki*ffru                                            5d24s21
              if(mlk.ne.0)ff=ff*srh                                      5d24s21
              if(mlb.ne.0)ff=ff*srh                                      5d24s21
              call psioppsi(iwavebra(1,ibi),nrootb,mur,1d0,ixmt,iosym,0,8d23s21
     $             idum,dum,iwaveket(1,iki),nrootk,bc(itmp),mdon,mdoop, 8d23s21
     $             ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,nvirt,  8d23s21
     $             maxbx,maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,npadddi,11d10s22
     $             bc,ibc)                                              11d10s22
              jfull=ifull+nbk*(ibra+nllb*iket)                           5d24s21
              do i=0,nrootk-1                                           8d23s21
               do j=0,nrootb-1                                          8d23s21
                ji=itmp+j+nrootb*i                                      8d23s21
                ij=jfull+i+nrootk*j                                     8d23s21
                bc(ij)=bc(ij)+bc(ji)*ff                                 8d23s21
               end do                                                   8d23s21
              end do                                                    8d23s21
             else if(isymc.eq.iosym(mui))then                            5d24s21
              if(lprint)write(6,*)('we have iii=i')                     3d2s22
              nimag=nimag+1
             end if                                                      5d24s21
            end if                                                       5d24s21
            ibcoff=itmp                                                 8d23s21
           end if                                                        5d24s21
          end do                                                         5d24s21
          if(nimag.ne.0)then                                              5d24s21
           write(6,*)('we are getting imaginary matrix elements!!! ')     5d24s21
           call dws_synca                                                 5d24s21
           call dws_finalize                                              5d24s21
           stop                                                           5d24s21
          end if                                                          5d24s21
          nn=nllb*nllk                                                    5d24s21
          ifullt=ibcoff                                                   5d24s21
          iscr=ifullt+nbk*nn                                              5d24s21
          iwgt=iscr+2+1+2*nn+nn                                           5d24s21
          icoef=iwgt+nn                                                   5d24s21
          ibcoff=icoef+nbk                                                5d24s21
          call enough('gtmom.  3',bc,ibc)
          do i=0,nn-1                                                     5d24s21
           bc(iwgt+i)=1d0                                                 5d24s21
          end do                                                          5d24s21
          do i=0,nn-1                                                     5d24s21
           do j=0,nbk-1                                                   5d24s21
            ji=ifull+j+nbk*i                                              5d24s21
            ij=ifullt+i+nn*j                                              5d24s21
            bc(ij)=bc(ji)                                                 5d24s21
           end do                                                         5d24s21
          end do                                                          5d24s21
          iunit=0                                                       11d21s22
          iuse=0                                                            5d18s21
          call lsqfit2(bc(icleb),nn,1,bc(ifullt),nn,nbk,nn,bc(icoef),1,     5d24s21
     $     bc(iscr),bc(iwgt),iunit,rms,bc,ibc)                          3d27s23
          nbigb=nllb*iwavebra(1,1)                                      5d25s21
          nbigk=nllk*iwaveket(1,1)                                      5d25s21
          ibig=ibcoff                                                   5d25s21
          jvecrb=ibig+nbigb*nbigk                                       6d3s21
          jvecrk=jvecrb+nbigb*nbigb                                     6d3s21
          ibcoff=jvecrk+nbigk*nbigk                                     6d3s21
          call enough('gtmom.  4',bc,ibc)
          isb2=iwavebra(1,1)-1                                          6d3s21
          isk2=iwaveket(1,1)-1                                          6d3s21
          jjminb=iabs(isb2-llb2)                                         6d3s21
          jjmaxb=isb2+llb2                                               6d3s21
          jjvecrb=jvecrb                                                6d3s21
          do jj=jjminb,jjmaxb,2                                           6d3s21
           do mj=-jj,jj,2                                               6d3s21
            do ml=-llb2,llb2,2                                          6d3s21
             do ms=-isb2,isb2,2                                         6d3s21
              bc(jjvecrb)=cleb2(llb2,ml,isb2,ms,jj,mj)                  6d3s21
              jjvecrb=jjvecrb+1                                         6d3s21
             end do                                                     6d3s21
            end do                                                      6d3s21
           end do                                                       6d3s21
          end do                                                        6d3s21
          jjmink=iabs(isk2-llk2)                                         6d3s21
          jjmaxk=isk2+llk2                                               6d3s21
          jjvecrk=jvecrk                                                6d3s21
          do jj=jjmink,jjmaxk,2                                           6d3s21
           do mj=-jj,jj,2                                               6d3s21
            do ml=-llk2,llk2,2                                          6d3s21
             do ms=-isk2,isk2,2                                         6d3s21
              bc(jjvecrk)=cleb2(llk2,ml,isk2,ms,jj,mj)                  6d3s21
              jjvecrk=jjvecrk+1                                         6d3s21
             end do                                                     6d3s21
            end do                                                      6d3s21
           end do                                                       6d3s21
          end do                                                        6d3s21
          jbig=ibig
          do mlk=-llk2,llk2,2
           iket=(mlk+llk2)/2
           do msk=-isk2,isk2,2
            do mlb=-llb2,llb2,2
             ibra=(mlb+llb2)/2
             do msb=-isb2,isb2,2
              if(msk.eq.msb)then
               jcleb=icleb+ibra+nllb*iket
               value=bc(jcleb)
              else
               value=0d0
              end if
              bc(jbig)=value
              jbig=jbig+1
             end do
            end do
           end do
          end do
          itmpjj=ibcoff                                                 6d3s21
          ibcoff=itmpjj+nbigb*nbigk                                     6d3s21
          call enough('gtmom.  5',bc,ibc)
          call dgemm('n','n',nbigb,nbigk,nbigk,1d0,bc(ibig),nbigb,      6d3s21
     $         bc(jvecrk),nbigk,0d0,bc(itmpjj),nbigb,                   6d3s21
     d' gtmom.  1')
          do i=0,nbigk-1                                                6d3s21
           do j=0,nbigb-1                                               6d3s21
            ji=itmpjj+j+nbigb*i                                         6d3s21
            ij=ibig+i+nbigk*j
            bc(ij)=bc(ji)
           end do
          end do
          call dgemm('n','n',nbigk,nbigb,nbigb,1d0,bc(ibig),nbigk,      6d3s21
     $         bc(jvecrb),nbigb,0d0,bc(itmpjj),nbigk,                   6d3s21
     d' gtmom.  2')
          do i=0,nbigk-1                                                6d3s21
           do j=0,nbigb-1                                               6d3s21
            ji=ibig+j+nbigb*i                                           6d3s21
            ij=itmpjj+i+nbigk*j
            bc(ji)=bc(ij)
           end do
          end do
          jbig=ibig
          rms=0d0
          nrms=0
          do jk=jjmink,jjmaxk,2
           do mjk=-jk,jk,2
            do jb=jjminb,jjmaxb,2
             f1=f6j(jb,jk,iq2,llk2,llb2,isb2,1)                         6d3s21
             f1=f1*sqrt(dfloat(jb+1)*dfloat(jk+1)*dfloat(llk2+1))
             do mjb=-jb,jb,2
              f2=f3j(jb,jk,iq2,-mjb,mjk,ipm2,1)
              trm=f1*f2
              if(abs(trm).gt.1d-10)then
               sum=0d0
               xsum=f3j(jb,jk,iq2,-mjb,mjk,ipm2,1)
     $             *f6j(jb,jk,iq2,llk2,llb2,isb2,1)
                  isum=-mjk+isb2+llk2
               if(mod(isum,2).ne.0)then
                write(6,*)('isum error!!! '),isum
               end if
               isum=isum/2
               if(mod(isum,2).ne.0)xsum=-xsum
               do mlb=-llb2,llb2,2
                do mlk=-llk2,llk2,2
                 do msb=-isb2,isb2,2
                  g1=f3j(jb,llb2,isb2,-mjb,mlb,msb,1)
                  g2=f3j(llk2,jk,isb2,-mlk,mjk,-msb,1)
                  g3=f3j(llk2,llb2,iq2,mlk,-mlb,ipm2,1)
                  isum=-mjk+isb2+iq2-llb2
     $                 +llk2+llb2+iq2
     $                 +llk2+llb2+isb2+mlk+mlb-msb
                  if(mod(isum,2).ne.0)then
                   write(6,*)('isum error!! '),isum
                  end if
                  isum=isum/2
                  ggg=g1*g2*g3
                  if(mod(isum,2).ne.0)ggg=-ggg
                  sum=sum+ggg
                 end do
                end do
               end do
               sum=sum*sqrt(dfloat(jk+1)*dfloat(jb+1)*dfloat(llk2+1))
               xsum=xsum*sqrt(dfloat(jk+1)*dfloat(jb+1)*dfloat(llk2+1))
                  isum=-mjk+isb2+llk2
               if(mod(isum,2).ne.0)then
                write(6,*)('isum error! '),isum
               end if
               isum=isum/2
               if(mod(isum,2).ne.0)trm=-trm
              end if                                                    6d3s21
              if(max(abs(trm),abs(bc(jbig))).gt.1d-10)then
               diff=bc(jbig)-trm
               if(abs(diff).gt.1d-10)
     $              write(6,*)jb,mjb,jk,mjk,trm,bc(jbig),diff,sum,xsum
               rms=rms+diff**2
               nrms=nrms+1
              end if
              jbig=jbig+1
             end do
            end do
           end do
          end do
          ibcoff=ibig
          do i=0,nbk-1                                                  5d25s21
           bc(itmom(ipass)+i)=bc(itmom(ipass)+i)+bc(icoef+i)            5d25s21
           bc(irms+i)=bc(irms+i)+bc(icoef+i)**2                         5d25s21
          end do                                                        5d25s21
          ibcoff=ifullt                                                   5d24s21
         end do                                                          5d25s21
         sz=0d0                                                         5d25s21
         do i=0,nbk-1                                                   5d25s21
          bc(itmom(ipass)+i)=bc(itmom(ipass)+i)*0.5d0                   5d25s21
          bc(irms+i)=bc(irms+i)*0.5d0                                   5d25s21
          bc(irms+i)=sqrt(abs(bc(irms+i)-bc(itmom(ipass)+i)**2))        5d25s21
          sz=sz+bc(irms+i)**2                                           5d25s21
         end do                                                         5d25s21
c
c     itmom has data stored ket,bra. Now transpose it.
c
         itrans=ibcoff                                                   6d4s21
         ibcoff=itrans+nbk                                              6d4s21
         call enough('gtmom.  6',bc,ibc)
         do i=0,nrootb-1                                                6d4s21
          do j=0,nrootk-1                                               6d4s21
           ji=itmom(ipass)+j+nrootk*i                                   6d4s21
           ij=itrans+i+nrootb*j                                         6d4s21
           bc(ij)=bc(ji)                                                6d4s21
          end do                                                        6d4s21
         end do                                                         6d4s21
         do i=0,nbk-1                                                   6d4s21
          bc(itmom(ipass)+i)=bc(itrans+i)                               6d4s21
         end do                                                         6d4s21
         if(lprint)then                                                 3d2s22
          do ik=1,nrootk                                                3d2s22
           do ib=1,nrootb                                               3d2s22
            iad=itmom(ipass)+ib-1+nrootb*(ik-1)                         3d2s22
            if(abs(bc(iad)).gt.1d-10)then                               3d2s22
             if(iheader.eq.0)then                                       1d6s23
              if(ipass.eq.1)then                                        1d6s23
               write(6,*)('>electric dipole')                           1d6s23
              else if(ipass.eq.2)then                                   1d6s23
               write(6,*)('>magnetic dipole')                           1d6s23
              else                                                      1d6s23
               write(6,*)('>electric quadrapole')                       1d6s23
              end if                                                    1d6s23
              iheader=1                                                 1d6s23
             end if                                                     1d6s23
             write(6,28)ib,iwavebra(1,1),wnameb,ik,iwaveket(1,1),       3d2s22
     $           wnamek,bc(iad)                                         3d2s22
   28        format('M',i3,1x,i1,a6,i3,1x,i1,a6,'>=',es14.6)            3d2s22
            end if                                                      3d2s22
           end do                                                       3d2s22
          end do                                                        3d2s22
         end if                                                         3d2s22
         ibcoff=itrans                                                  6d4s21
        end if                                                          5d25s21
       end do                                                           5d24s21
       ibcoff=irms                                                      5d25s21
       ibcoffo=irms                                                     11d21s22
      else if(nlzz.eq.2)then                                            5d26s21
       if(lprint)write(6,*)('linear molecule case ')                    3d2s22
       q0fact=2d0/sqrt(6d0)                                             3d15s22
       if(llb.eq.0)then                                                 5d26s21
        nllb=1                                                          5d26s21
       else                                                             5d26s21
        nllb=2                                                          5d26s21
       end if                                                           5d26s21
       if(llk.eq.0)then                                                 5d26s21
        nllk=1                                                          5d26s21
       else                                                             5d26s21
        nllk=2                                                          5d26s21
       end if                                                           5d26s21
       nrootb=iwavebra(3,1)                                             5d24s21
       nrootk=iwaveket(3,1)                                             5d24s21
       nbk=nrootb*nrootk                                                5d24s21
       ifull=ibcoff                                                     3d16s21
       ibcoff=ifull+nbk*max(nllb*nllk,2)                                3d16s21
       call enough('gtmom.  7',bc,ibc)
       do i=ifull,ibcoff-1                                              5d26s21
        bc(i)=0d0                                                       5d26s21
       end do                                                           5d26s21
       iextob=0                                                         5d26s21
       if(llb.eq.0.and.(iwavebra(2,1).eq.8.or.iwavebra(2,1).eq.4))then  5d26s21
        iextob=1                                                        5d26s21
       end if                                                           5d26s21
       iextok=0                                                         5d26s21
       if(llk.eq.0.and.(iwaveket(2,1).eq.8.or.iwaveket(2,1).eq.4))then  5d26s21
        iextok=1                                                        5d26s21
       end if                                                           5d26s21
       mlba=iabs(llb)                                                   5d26s21
       do i=1,nllb                                                      5d26s21
        npack4=iwavebra(6,i)                                            5d26s21
        if(ipack1(4).eq.mlba)ibr=i                                      5d26s21
        if(ipack1(4).eq.-mlba)ibi=i                                     5d26s21
       end do                                                           5d26s21
       mlka=iabs(llk)                                                   5d26s21
       do i=1,nllk                                                      5d26s21
        npack4=iwaveket(6,i)                                            5d26s21
        if(ipack1(4).eq.mlka)ikr=i                                      5d26s21
        if(ipack1(4).eq.-mlka)iki=i                                     5d26s21
       end do                                                           5d26s21
       do ipass=1,3                                                     5d26s21
        if(ipass.eq.1)then                                              5d24s21
         c2='e1'                                                        3d8s22
         iq=1                                                           5d24s21
c        +1                 -1               0
c     (-x-iy)/sqrt(2), (+x-iy)/sqrt(2)       z
c
c     i.e. real +1=-x/sqrt(2), imag +1=-y/sqrt(2)                       5d24s21
c          real -1=+x/sqrt(2), imag -1= same
         do i=1,nop                                                     5d24s21
          if(opname(i)(1:3).eq.'mux')mur=i
          if(opname(i)(1:3).eq.'muy')mui=i                                 5d24s21
          if(opname(i)(1:3).eq.'mu0')mu0=i                              5d26s21
         end do
         ffr=-srh                                                       5d24s21
         ffi=-srh                                                       5d24s21
        else if(ipass.eq.2)then                                         5d24s21
         c2='m1'                                                        3d8s22
         iq=1                                                           5d24s21
c        +1                    -1                0
c     (-lx-ily)/sqrt(2),(lx-ily)/sqrt(2)         lz=i*lz(i)
c     =(+ilx(i)-ly(i))/sqrt(2),(-ilx(i)-ly(i))/sqrt(2)
c
c     i.e. real +1 = -ly(i)/sqrt(2), imag +1=lx(i)/sqrt(2)
c          real -1 = same          , imag -1=-lx(i)/sqrt(2)
         do i=1,nop                                                     5d24s21
          if(opname(i)(1:4).eq.'lx/i')mui=i                                5d24s21
          if(opname(i)(1:4).eq.'ly/i')mur=i                                5d24s21
          if(opname(i)(1:4).eq.'lz/i')mu0=i                                5d24s21
         end do
         ffr=-srh                                                       5d24s21
         ffi=srh                                                        5d24s21
        else                                                            5d24s21
         iq=2                                                           5d24s21
         c2='e2'                                                        3d8s22
c        +1     -1    0   +/- 2
c     -zx-izy,zx-izy q0  Re Q2 +/- i Im Q2
c
c     i.e. real +1 = -zx, imag +1 = -zy
c          real -1 = +zx, imag -1 = same
         do i=1,nop                                                     5d24s21
          if(opname(i)(1:5).eq.'Re Q1')mur=i                               5d24s21
          if(opname(i)(1:5).eq.'Im Q1')mui=i                               5d24s21
          if(opname(i)(1:2).eq.'Q0')mu0=i                               5d26s21
          if(opname(i)(1:5).eq.'Re Q2')mur2=i                           5d26s21
          if(opname(i)(1:5).eq.'Im Q2')mui2=i                           5d26s21
         end do
         ffr=1d0                                                        5d26s21
         ffi=1d0                                                        5d26s21
        end if                                                          5d24s21
        iq2=iq*2                                                        5d26s21
        mlb=mlba                                                        3d16s21
        pbr=1d0                                                         3d16s21
        pbi=1d0                                                         3d16s21
        if(llb.gt.0)then                                                3d9s22
         if(mod(llb,2).eq.0)then                                        3d9s22
          if(mlb.gt.0)pbi=-1d0                                          3d9s22
         else                                                           3d9s22
          pbi=-1d0                                                      3d9s22
          if(mlb.lt.0)pbr=-1d0                                          3d9s22
         end if                                                         3d9s22
        end if                                                          3d9s22
        do iz=0,nbk*2-1                                                 3d16s21
         bc(ifull+iz)=0d0                                               3d16s21
        end do                                                          3d16s21
        ihit=0                                                          3d16s21
        nimag=0                                                         3d16s21
        ipmax=0                                                         3d16s21
        xnuc=0d0                                                        3d16s21
        do ipm=-iq,iq                                                   5d26s21
c                      ipm=-1   ipm=+1  ipm=-2 ipm=+2
c                       r  i     r  i    r  i   r  i
c     electric dipole:  -  +     +  +
c     magnetic dipole:  +  -     +  +
c     electric quad:    -  +     +  +    +  -   +  +
         fmr=ffr                                                        3d9s22
         fmi=ffi                                                        3d9s22
         if(ipm.ge.0)then                                               5d26s21
         else if(ipm.eq.-1)then                                         5d26s21
          if(ipass.eq.2)then                                            5d26s21
           fmi=-fmi                                                     3d9s22
          else                                                          5d26s21
           fmr=-fmr                                                     3d9s22
          end if                                                        5d26s21
         else if(ipm.eq.-2)then                                         5d26s21
          fmi=-fmi                                                      3d9s22
         end if                                                         5d26s21
         if(loc(iwavebra).eq.loc(iwaveket))then                         3d15s22
          if(ipass.eq.1)then                                             3d15s22
           if(ipm.eq.0)then                                              3d15s22
            do ia=1,natom                                                3d15s22
             xnuc=xnuc+xcart(3,ia)*atnum(1,ia)                           3d15s22
            end do                                                       3d15s22
           end if                                                        3d15s22
          else if(ipass.eq.3)then                                        3d15s22
           if(ipm.eq.0)then                                             3d15s22
            do ia=1,natom                                               3d15s22
             xnuc=xnuc+xcart(3,ia)*xcart(3,ia)*atnum(1,ia)              3d15s22
            end do                                                      3d15s22
            xnuc=xnuc*q0fact                                            3d15s22
           end if                                                       3d15s22
          end if                                                         3d15s22
         end if                                                         3d15s22
         mlk=mlb-ipm                                                    5d25s21
         mlka=iabs(mlk)                                                 5d24s21
         if(mlka.eq.iabs(llk))then                                      5d26s21
          iket=0                                                        5d26s21
          if(mlk.lt.0)iket=1                                            5d26s21
          iketp=iket+1                                                  3d16s21
          ipmv(iketp)=ipm                                               3d16s21
          ipmax=max(ipmax,iketp)                                        3d16s21
          pkr=1d0                                                       5d26s21
          pki=1d0                                                       5d26s21
          if(llk.gt.0)then                                                 3d9s22
           if(mod(llk,2).eq.0)then                                          3d9s22
            if(mlk.gt.0)pki=-1d0                                          3d9s22
           else                                                            3d9s22
            pki=-1d0                                                       3d9s22
            if(mlk.lt.0)pkr=-1d0                                          3d9s22
           end if                                                          3d9s22
          end if                                                           3d9s22
          itmp=ibcoff                                                   3d16s21
          ibcoff=itmp+nbk                                               3d16s21
          call enough('gtmom.  8',bc,ibc)
          if(ipm.eq.0)then                                              5d26s21
           if((ibr.ne.ibi.or.iextob.eq.0).and.                          5d26s21
     $         (ikr.ne.iki.or.iextok.eq.0))then                         5d26s21
            itest=multh(iwavebra(2,ibr),iwaveket(2,ikr))                5d26s21
            if(itest.eq.iosym(mu0))then                                 5d26s21
             ihit=ihit+1                                                5d26s21
             ff=pbr*pkr                                                 5d26s21
             if(mlk.ne.0)ff=ff*srh                                      5d24s21
             if(mlb.ne.0)ff=ff*srh                                      5d24s21
             call psioppsi(iwavebra(1,ibr),nrootb,mu0,1d0,ixmt,iosym,0, 8d23s21
     $         idum,dum,iwaveket(1,ikr),nrootk,bc(itmp),mdon,mdoop,     8d23s21
     $         ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,nvirt,      8d23s21
     $         maxbx,maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,npadddi,bc, 11d10s22
     $            ibc)                                                  11d10s22
             jfull=ifull+nbk*iket                                       3d16s21
             do i=0,nrootk-1                                            8d23s21
              do j=0,nrootb-1                                           8d23s21
               ji=itmp+j+nrootb*i                                       8d23s21
               ij=jfull+i+nrootk*j                                      8d23s21
               bc(ij)=bc(ij)+bc(ji)*ff                                  8d23s21
              end do                                                    8d23s21
             end do                                                     8d23s21
            end if                                                      5d26s21
           end if                                                       5d26s21
           if((ibr.ne.ibi.or.iextob.eq.0).and.                          5d26s21
     $       (ikr.ne.iki.or.iextok.eq.1))then                           5d26s21
            itest=multh(iwavebra(2,ibr),iwaveket(2,iki))                5d26s21
            if(itest.eq.iosym(mu0).and.ipass.eq.2)then                  5d26s21
             ihit=ihit+1                                                5d26s21
             ff=-pbr*pki                                                5d26s21
             if(mlk.ne.0)ff=ff*srh                                      5d24s21
             if(mlb.ne.0)ff=ff*srh                                      5d24s21
             call psioppsi(iwavebra(1,ibr),nrootb,mu0,1d0,ixmt,iosym,0, 8d23s21
     $       idum,dum,iwaveket(1,iki),nrootk,bc(itmp),mdon,mdoop,       8d23s21
     $       ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,nvirt,        8d23s21
     $       maxbx,maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,npadddi,bc,   11d10s22
     $            ibc)                                                  11d10s22
             jfull=ifull+nbk*iket                                       3d16s21
             do i=0,nrootk-1                                            8d23s21
              do j=0,nrootb-1                                           8d23s21
               ji=itmp+j+nrootb*i                                       8d23s21
               ij=jfull+i+nrootk*j                                      8d23s21
               bc(ij)=bc(ij)+bc(ji)*ff                                  8d23s21
              end do                                                    8d23s21
             end do                                                     8d23s21
            end if                                                      5d26s21
           end if                                                       5d26s21
           if((ibr.ne.ibi.or.iextob.eq.1).and.                          5d26s21
     $     (ikr.ne.iki.or.iextok.eq.0))then                             5d26s21
            itest=multh(iwavebra(2,ibi),iwaveket(2,ikr))                5d26s21
            if(itest.eq.iosym(mu0).and.ipass.eq.2)then                                 5d26s21
             ihit=ihit+1                                                5d26s21
             ff=pbi*pkr                                                 5d26s21
             if(mlk.ne.0)ff=ff*srh                                      5d24s21
             if(mlb.ne.0)ff=ff*srh                                      5d24s21
             call psioppsi(iwavebra(1,ibi),nrootb,mu0,1d0,ixmt,iosym,0, 8d23s21
     $     idum,dum,iwaveket(1,ikr),nrootk,bc(itmp),mdon,mdoop,         8d23s21
     $     ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,nvirt,          8d23s21
     $     maxbx,maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,npadddi,bc,ibc) 11d10s22
             jfull=ifull+nbk*iket                                       3d16s21
             do i=0,nrootk-1                                            8d23s21
              do j=0,nrootb-1                                           8d23s21
               ji=itmp+j+nrootb*i                                       8d23s21
               ij=jfull+i+nrootk*j                                      8d23s21
               bc(ij)=bc(ij)+bc(ji)*ff                                  8d23s21
              end do                                                    8d23s21
             end do                                                     8d23s21
            end if                                                      5d26s21
           end if                                                       5d26s21
           if((ibr.ne.ibi.or.iextob.eq.1).and.                          5d26s21
     $     (ikr.ne.iki.or.iextok.eq.1))then                             5d26s21
            itest=multh(iwavebra(2,ibi),iwaveket(2,iki))                5d26s21
            if(itest.eq.iosym(mu0))then                                 5d26s21
             ihit=ihit+1                                                5d26s21
             ff=pbi*pki                                                 5d26s21
             if(mlk.ne.0)ff=ff*srh                                      5d24s21
             if(mlb.ne.0)ff=ff*srh                                      5d24s21
             call psioppsi(iwavebra(1,ibi),nrootb,mu0,1d0,ixmt,iosym,0, 8d23s21
     $     idum,dum,iwaveket(1,iki),nrootk,bc(itmp),mdon,mdoop,         8d23s21
     $     ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,nvirt,          8d23s21
     $     maxbx,maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,npadddi,bc,ibc) 11d10s22
             jfull=ifull+nbk*iket                                       3d16s21
             do i=0,nrootk-1                                            8d23s21
              do j=0,nrootb-1                                           8d23s21
               ji=itmp+j+nrootb*i                                       8d23s21
               ij=jfull+i+nrootk*j                                      8d23s21
               bc(ij)=bc(ij)+bc(ji)*ff                                  8d23s21
              end do                                                    8d23s21
             end do                                                     8d23s21
            end if                                                      5d26s21
           end if                                                       5d26s21
          else if(iabs(ipm).ge.1)then                                   5d26s21
           if(iabs(ipm).eq.1)then                                       5d26s21
            muru=mur                                                    5d26s21
            muiu=mui                                                    5d26s21
           else                                                         5d26s21
            muru=mur2                                                   5d26s21
            muiu=mui2                                                   5d26s21
           end if                                                       5d26s21
           if((ibr.ne.ibi.or.iextob.eq.0).and.                          5d26s21
     $     (ikr.ne.iki.or.iextok.eq.0))then                             5d26s21
            itest=multh(iwavebra(2,ibr),iwaveket(2,ikr))                5d26s21
            if(itest.eq.iosym(muru))then                                 5d26s21
             ihit=ihit+1                                                5d26s21
             ff=pbr*pkr*fmr                                             5d26s21
             if(mlk.ne.0)ff=ff*srh                                      5d24s21
             if(mlb.ne.0)ff=ff*srh                                      5d24s21
             call psioppsi(iwavebra(1,ibr),nrootb,muru,1d0,ixmt,iosym,  8d23s21
     $     0,idum,dum,iwaveket(1,ikr),nrootk,bc(itmp),mdon,             8d23s21
     $     mdoop,ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,          8d23s21
     $     nvirt,maxbx,maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,          11d15s21
     $     npadddi,bc,ibc)                                              11d10s22
             jfull=ifull+nbk*iket                                       3d16s21
             do i=0,nrootk-1                                            8d23s21
              do j=0,nrootb-1                                           8d23s21
               ji=itmp+j+nrootb*i                                       8d23s21
               ij=jfull+i+nrootk*j                                      8d23s21
               bc(ij)=bc(ij)+bc(ji)*ff                                  8d23s21
              end do                                                    8d23s21
             end do                                                     8d23s21
            else if(itest.eq.iosym(muiu))then                            5d26s21
             ihit=ihit+1                                                5d26s21
             nimag=nimag+1
            end if                                                      5d26s21
           end if                                                       5d26s21
           if((ibr.ne.ibi.or.iextob.eq.0).and.                          5d26s21
     $     (ikr.ne.iki.or.iextok.eq.1))then                             5d26s21
            itest=multh(iwavebra(2,ibr),iwaveket(2,iki))                5d26s21
            if(itest.eq.iosym(muru))then                                 5d26s21
             ihit=ihit+1                                                5d26s21
             nimag=nimag+1
            else if(itest.eq.iosym(muiu))then                                 5d26s21
             ihit=ihit+1                                                5d26s21
             ff=-pbr*pki*fmi                                            5d26s21
             if(mlk.ne.0)ff=ff*srh                                      5d24s21
             if(mlb.ne.0)ff=ff*srh                                      5d24s21
             call psioppsi(iwavebra(1,ibr),nrootb,muiu,1d0,ixmt,iosym,  8d23s21
     $      0,idum,dum,iwaveket(1,iki),nrootk,bc(itmp),mdon,            8d23s21
     $      mdoop,ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,         8d23s21
     $      nvirt,maxbx,maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,         11d15s21
     $      npadddi,bc,ibc)                                             11d10s22
             jfull=ifull+nbk*iket                                       3d16s21
             do i=0,nrootk-1                                            8d23s21
              do j=0,nrootb-1                                           8d23s21
               ji=itmp+j+nrootb*i                                       8d23s21
               ij=jfull+i+nrootk*j                                      8d23s21
               bc(ij)=bc(ij)+bc(ji)*ff                                  8d23s21
              end do                                                    8d23s21
             end do                                                     8d23s21
            end if                                                      5d26s21
           end if                                                       5d26s21
           if((ibr.ne.ibi.or.iextob.eq.1).and.                          5d26s21
     $     (ikr.ne.iki.or.iextok.eq.0))then                             5d26s21
            itest=multh(iwavebra(2,ibi),iwaveket(2,ikr))                5d26s21
            if(itest.eq.iosym(muru))then                                 5d26s21
             ihit=ihit+1                                                5d26s21
             nimag=nimag+1
            else if(itest.eq.iosym(muiu))then                                 5d26s21
             ihit=ihit+1                                                5d26s21
             ff=pbi*pkr*fmi                                             5d26s21
             if(mlk.ne.0)ff=ff*srh                                      5d24s21
             if(mlb.ne.0)ff=ff*srh                                      5d24s21
             call psioppsi(iwavebra(1,ibi),nrootb,muiu,1d0,ixmt,iosym,  8d23s21
     $      0,idum,dum,iwaveket(1,ikr),nrootk,bc(itmp),mdon,            8d23s21
     $      mdoop,ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,         8d23s21
     $      nvirt,maxbx,maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,         11d15s21
     $      npadddi,bc,ibc)                                             11d10s22
             jfull=ifull+nbk*iket                                       3d16s21
             do i=0,nrootk-1                                            8d23s21
              do j=0,nrootb-1                                           8d23s21
               ji=itmp+j+nrootb*i                                       8d23s21
               ij=jfull+i+nrootk*j                                      8d23s21
               bc(ij)=bc(ij)+bc(ji)*ff                                  8d23s21
              end do                                                    8d23s21
             end do                                                     8d23s21
            end if                                                      5d26s21
           end if                                                       5d26s21
           if((ibr.ne.ibi.or.iextob.eq.1).and.                          5d26s21
     $     (ikr.ne.iki.or.iextok.eq.1))then                             5d26s21
            itest=multh(iwavebra(2,ibi),iwaveket(2,iki))                5d26s21
            if(itest.eq.iosym(muru))then                                 5d26s21
             ihit=ihit+1                                                5d26s21
             ff=pbi*pki*fmr                                             5d26s21
             if(mlk.ne.0)ff=ff*srh                                      5d24s21
             if(mlb.ne.0)ff=ff*srh                                      5d24s21
             call psioppsi(iwavebra(1,ibi),nrootb,muru,1d0,ixmt,iosym,  8d23s21
     $     0,idum,dum,iwaveket(1,iki),nrootk,bc(itmp),mdon,             8d23s21
     $     mdoop,ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,          8d23s21
     $     nvirt,maxbx,maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,          11d15s21
     $     npadddi,bc,ibc)                                              11d10s22
             jfull=ifull+nbk*iket                                       3d16s21
             do i=0,nrootk-1                                            8d23s21
              do j=0,nrootb-1                                           8d23s21
               ji=itmp+j+nrootb*i                                       8d23s21
               ij=jfull+i+nrootk*j                                      8d23s21
               bc(ij)=bc(ij)+bc(ji)*ff                                  8d23s21
              end do                                                    8d23s21
             end do                                                     8d23s21
            else if(itest.eq.iosym(muiu))then                                 5d26s21
             ihit=ihit+1                                                5d26s21
             nimag=nimag+1
            end if                                                      5d26s21
           end if                                                       5d26s21
          end if                                                        5d26s21
          ibcoff=itmp                                                   3d16s21
         end if                                                         3d16s21
        end do                                                          3d16s21
        if(ihit.ne.0.and.lprint)then                                    3d16s21
         igot=0                                                         3d16s21
         do ik=1,nrootk                                                 3d8s22
          do ib=1,nrootb                                                3d8s22
           if(ib.eq.ik)then                                             3d15s22
            xnucu=xnuc                                                  3d15s22
           else                                                         3d15s22
            xnucu=0d0                                                   3d15s22
           end if                                                       3d15s22
           iad1=ifull+ik-1+nrootk*(ib-1)                                3d16s21
           iad2=iad1+nbk                                                3d16s21
           if(max(abs(bc(iad1)),abs(bc(iad2))).gt.1d-10)then            3d16s21
            if(ipmax.eq.1)then                                          3d16s21
             if(igot.eq.0)write(6,41)(ipmv(ip),ip=1,ipmax)               3d16s21
   41        format(/,'>',29x,'normal, mq=',i3)                         3d16s21
             write(6,40)c2,ib,iwavebra(1,1),wnameb,                      3d16s21
     $             ik,iwaveket(1,1),wnamek,bc(iad1)+xnucu
            else                                                        3d16s21
             if(igot.eq.0)write(6,39)(ipmv(ip),ip=1,ipmax)               3d16s21
   39        format(/,'>',29x,'normal, mq=',i3,x,'exchange, mq=',i3)    3d16s21
             write(6,40)c2,ib,iwavebra(1,1),wnameb,                      3d16s21
     $             ik,iwaveket(1,1),wnamek,bc(iad1)+xnucu,bc(iad2)      3d16s21
   40        format('M(',a2,')',i2,x,i1,a6,x,i2,x,i1,                    3d16s21
     $              a6,'>=',2f14.8,2i8)                                     3d16s21
            end if                                                      3d16s21
            igot=1                                                      3d8s22
           end if                                                       3d16s21
          end do                                                        3d16s21
         end do                                                         3d16s21
        end if                                                          3d16s21
       end do                                                           5d26s21
      else                                                              5d24s21
       npack4=iwavebra(6,1)                                             3d28s22
       iordb=ipack1(3)                                                  3d28s22
       npack4=iwaveket(6,1)                                             3d28s22
       iordk=ipack1(3)                                                  3d28s22
       q0fact=1d0/sqrt(6d0)                                             3d17s22
       nrootb=iwavebra(3,1)                                             5d24s21
       nrootk=iwaveket(3,1)                                             5d24s21
       nbk=nrootb*nrootk                                                5d24s21
       ifull=ibcoff                                                     3d16s21
       ibcoff=ifull+nbk*2                                               3d17s22
       call enough('gtmom.  9',bc,ibc)
       isymneed=multh(iwavebra(2,1),iwaveket(2,1))                      3d17s22
       lzlx=-1                                                          3d22s22
       do i=1,nop
        if(opname(i)(1:4).eq.'LzLx')lzlx=i                              3d22s22
       end do
       if(lzlx.lt.0)then                                                3d22s22
        write(6,*)('oh where is lzlx??? ')
        call dws_synca
        call dws_finalize
        stop                                                            3d22s22
       end if                                                           3d22s22
       if(isymneed.eq.iosym(lzlx))then                                  3d22s22
        itmp=ibcoff                                                     3d22s22
        ibcoff=itmp+nbk                                                 3d22s22
        call enough('gtmom. 10',bc,ibc)
        call psioppsi(iwavebra(1,1),nrootb,lzlx,1d0,ixmt,iosym,0,       3d22s22
     $         idum,dum,iwaveket(1,1),nrootk,bc(itmp),mdon,mdoop,       3d22s22
     $         ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,nvirt,      3d22s22
     $         maxbx,maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,npadddi,bc, 11d10s22
     $       ibc)                                                       11d10s22
        if(lprint)then                                                  3d29s22
         do ik=1,nrootk                                                 3d29s22
          jk=ik-1                                                       3d29s22
          do ib=1,nrootb                                                3d29s22
           jb=ib-1                                                      3d29s22
           iad=itmp+jb+nrootb*jk                                        3d29s22
           if(abs(bc(iad)).gt.1d-10)write(6,401)ib,iwavebra(1,1),       3d29s22
     $          wnameb,ik,iwaveket(1,1),wnamek,bc(iad)                  3d29s22
  401      format('<',i2,x,i1,a6,'|LzLx+LxLz|',i2,x,i1,a6,'>=',es15.7)  3d29s22
          end do                                                        3d29s22
         end do                                                         3d29s22
        end if                                                          3d29s22
        ibcoff=itmp                                                     3d22s22
       end if                                                           3d22s22
       do ipass=1,3                                                     5d26s21
        if(ipass.eq.1)then                                              5d24s21
         c2='e1'                                                        3d8s22
         iq=1                                                           5d24s21
c        +1                 -1               0
c     (-x-iy)/sqrt(2), (+x-iy)/sqrt(2)       z
c
c     i.e. real +1=-x/sqrt(2), imag +1=-y/sqrt(2)                       5d24s21
c          real -1=+x/sqrt(2), imag -1= same
         do i=1,nop                                                     5d24s21
          if(opname(i)(1:3).eq.'mux')mur=i
          if(opname(i)(1:3).eq.'muy')mui=i                                 5d24s21
          if(opname(i)(1:3).eq.'mu0')mu0=i                              5d26s21
         end do
         ffr=-srh                                                       5d24s21
         ffi=-srh                                                       5d24s21
        else if(ipass.eq.2)then                                         5d24s21
         c2='m1'                                                        3d8s22
         iq=1                                                           5d24s21
c        +1                    -1                0
c     (-lx-ily)/sqrt(2),(lx-ily)/sqrt(2)         lz=i*lz(i)
c     =(+ilx(i)-ly(i))/sqrt(2),(-ilx(i)-ly(i))/sqrt(2)
c
c     i.e. real +1 = -ly(i)/sqrt(2), imag +1=lx(i)/sqrt(2)
c          real -1 = same          , imag -1=-lx(i)/sqrt(2)
         do i=1,nop                                                     5d24s21
          if(opname(i)(1:4).eq.'lx/i')mui=i                                5d24s21
          if(opname(i)(1:4).eq.'ly/i')mur=i                                5d24s21
          if(opname(i)(1:4).eq.'lz/i')mu0=i                                5d24s21
         end do
         ffr=-srh                                                       5d24s21
         ffi=srh                                                        5d24s21
        else                                                            5d24s21
         iq=2                                                           5d24s21
         c2='e2'                                                        3d8s22
c        +1     -1    0   +/- 2
c     -zx-izy,zx-izy q0  Re Q2 +/- i Im Q2
c
c     i.e. real +1 = -zx, imag +1 = -zy
c          real -1 = +zx, imag -1 = same
         do i=1,nop                                                     5d24s21
          if(opname(i)(1:5).eq.'Re Q1')mur=i                               5d24s21
          if(opname(i)(1:5).eq.'Im Q1')mui=i                               5d24s21
          if(opname(i)(1:2).eq.'Q0')mu0=i                               5d26s21
          if(opname(i)(1:5).eq.'Re Q2')mur2=i                           5d26s21
          if(opname(i)(1:5).eq.'Im Q2')mui2=i                           5d26s21
         end do
         ffr=1d0                                                        5d26s21
         ffi=1d0                                                        5d26s21
        end if                                                          5d24s21
        iq2=iq*2                                                        5d26s21
        do iz=0,nbk*2-1                                                 3d16s21
         bc(ifull+iz)=0d0                                               3d16s21
        end do                                                          3d16s21
        ihit=0                                                          3d16s21
        nimag=0                                                         3d16s21
        ipmax=0                                                         3d16s21
        xnuc=0d0                                                        3d16s21
        do ipm=-iq,iq                                                   5d26s21
c                      ipm=-1   ipm=+1  ipm=-2 ipm=+2
c                       r  i     r  i    r  i   r  i
c     electric dipole:  -  +     +  +
c     magnetic dipole:  +  -     +  +
c     electric quad:    -  +     +  +    +  -   +  +
         fmr=ffr                                                        3d9s22
         fmi=ffi                                                        3d9s22
         if(ipm.ge.0)then                                               5d26s21
         else if(ipm.eq.-1)then                                         5d26s21
          if(ipass.eq.2)then                                            5d26s21
           fmi=-fmi                                                     3d9s22
          else                                                          5d26s21
           fmr=-fmr                                                     3d9s22
          end if                                                        5d26s21
         else if(ipm.eq.-2)then                                         5d26s21
          fmi=-fmi                                                      3d9s22
         end if                                                         5d26s21
         xnuc=0d0                                                       3d17s22
         if(loc(iwavebra).eq.loc(iwaveket))then                         3d15s22
          if(ipass.eq.1)then                                             3d15s22
           if(ipm.eq.0)then                                              3d15s22
            do ia=1,natom                                                3d15s22
             xnuc=xnuc+xcart(3,ia)*atnum(1,ia)                           3d15s22
            end do                                                       3d15s22
           else if(ipm.gt.0)then                                        3d17s22
            do ia=1,natom                                                3d15s22
             xnuc=xnuc+xcart(1,ia)*atnum(1,ia)                           3d15s22
            end do                                                       3d15s22
           else                                                         3d17s22
            do ia=1,natom                                                3d15s22
             xnuc=xnuc+xcart(2,ia)*atnum(1,ia)                           3d15s22
            end do                                                       3d15s22
           end if                                                        3d15s22
          else if(ipass.eq.3)then                                        3d15s22
           if(ipm.eq.0)then                                             3d15s22
            do ia=1,natom                                               3d15s22
             xnuc=xnuc+(2d0*(xcart(3,ia)**2)-xcart(1,ia)**2             3d17s22
     $            -xcart(2,ia)**2)*atnum(1,ia)                          3d17s22
            end do                                                      3d15s22
            xnuc=xnuc*q0fact                                            3d15s22
           else if(ipm.eq.1)then                                        3d17s22
            do ia=1,natom                                               3d15s22
             xnuc=xnuc-xcart(1,ia)*xcart(3,ia)*atnum(1,ia)              3d17s22
            end do                                                      3d15s22
           else if(ipm.eq.-1)then                                        3d17s22
            do ia=1,natom                                               3d15s22
             xnuc=xnuc-xcart(2,ia)*xcart(3,ia)*atnum(1,ia)              3d17s22
            end do                                                      3d15s22
           else if(ipm.eq.2)then                                        3d17s22
            do ia=1,natom                                               3d15s22
             xnuc=xnuc+0.5d0*(xcart(1,ia)**2-xcart(2,ia)**2)*atnum(1,ia)8d29s22
            end do                                                      3d15s22
           else if(ipm.eq.-2)then                                        3d17s22
            do ia=1,natom                                               3d15s22
             xnuc=xnuc+xcart(1,ia)*xcart(2,ia)*atnum(1,ia)              3d17s22
            end do                                                      3d15s22
           end if                                                       3d15s22
          end if                                                         3d15s22
         end if                                                         3d15s22
         itmp=ibcoff                                                    3d17s22
         ibcoff=itmp+nbk                                                3d17s22
         call enough('gtmom. 11',bc,ibc)
         if(ipm.eq.0)then                                               3d17s22
          if(isymneed.eq.iosym(mu0))then                                3d17s22
           call psioppsi(iwavebra(1,1),nrootb,mu0,1d0,ixmt,iosym,0,     3d17s22
     $         idum,dum,iwaveket(1,1),nrootk,bc(itmp),mdon,mdoop,       3d17s22
     $         ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,nvirt,      8d23s21
     $         maxbx,maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,npadddi,bc, 11d10s22
     $         ibc)                                                     11d10s22
           ihit=ihit+1                                                  3d17s22
           if(lprint)then                                               3d29s22
            lok=(iordb.eq.iordk.and.ipass.ne.2).or.                      3d28s22
     $         (iordb.ne.iordk.and.ipass.eq.2)                          3d28s22
            iok=2
            if(lok)iok=1
            if(ipass.eq.1)then                                           3d28s22
             ctype='mu0 '                                                3d30s22
            else if(ipass.eq.2)then                                      3d28s22
             ctype=' Lz '                                                3d30s22
            else                                                         3d28s22
             ctype=' Q0 '                                                3d30s22
            end if                                                       3d28s22
            if(xnuc.ne.0d0)then                                          3d17s22
             do i=0,nrootb-1                                             3d17s22
              iad=itmp+i*(nrootb+1)                                      3d17s22
              bc(iad)=bc(iad)+xnuc                                       3d17s22
             end do                                                      3d17s22
            end if                                                       3d17s22
            do ik=1,nrootk                                              3d30s22
             jk=ik-1                                                    3d30s22
             do ib=1,nrootb                                             3d30s22
              jb=ib-1                                                   3d30s22
              iad=itmp+jb+nrootb*jk                                     3d30s22
              if(abs(bc(iad)).gt.1d-10)write(6,333)ib,iwavebra(1,1),    3d30s22
     $             wnameb,ctype,ik,iwaveket(1,1),wnamek,bc(iad),        3d30s22
     $             crori(iok)                                           3d30s22
  333         format('<',i2,x,i1,a6,'|',a4,'|',i2,x,i1,a6,'>=',         3d30s22
     $             es22.14,x,a1)                                         3d30s22
             end do                                                     3d30s22
            end do                                                      3d30s22
           end if                                                       3d29s22
          end if                                                        3d17s22
         else if(iabs(ipm).ge.1)then                                    5d26s21
          if(iabs(ipm).eq.1)then                                        5d26s21
           muru=mur                                                     5d26s21
           muiu=mui                                                     5d26s21
          else                                                          5d26s21
           muru=mur2                                                    5d26s21
           muiu=mui2                                                    5d26s21
          end if                                                        5d26s21
          if(isymneed.eq.iosym(muru).and.ipm.gt.0)then                  3d17s22
           ihit=ihit+1                                                  3d17s22
           call psioppsi(iwavebra(1,1),nrootb,muru,1d0,ixmt,iosym,      3d17s22
     $     0,idum,dum,iwaveket(1,1),nrootk,bc(itmp),mdon,               3d17s22
     $     mdoop,ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,          8d23s21
     $     nvirt,maxbx,maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,          11d15s21
     $     npadddi,bc,ibc)                                              11d10s22
           lok=iordb.eq.iordk
           if(xnuc.ne.0d0)then                                          3d17s22
            do i=0,nrootb-1                                             3d17s22
             iad=itmp+i*(nrootb+1)                                      3d17s22
             bc(iad)=bc(iad)+xnuc                                       3d17s22
            end do                                                      3d17s22
           end if                                                       3d17s22
           do i=0,nbk-1                                                 8d29s22
            bc(itmp+i)=bc(itmp+i)*fmr                                   3d17s22
           end do                                                       3d17s22
           if(lprint)then                                               3d29s22
            inz=0
            iok=2                                                       3d30s22
            if(lok)iok=1                                                3d30s22
            if(nsymb.eq.1)then                                          3d30s22
             if(ipass.eq.1)then                                          3d29s22
              ctype='Rm+1'
             else if(ipass.eq.2)then                                     3d29s22
              ctype='RL+1'                                               3d30s22
             else                                                        3d29s22
              write(ctype,303)ipm                                        3d30s22
  303         format('RQ',i2)                                            3d30s22
             end if                                                      3d29s22
            else                                                        3d30s22
             if(ipass.eq.1)then                                          3d29s22
              ctype='mu+1'
             else if(ipass.eq.2)then                                     3d29s22
              ctype='L+1 '                                               3d30s22
             else                                                        3d29s22
              write(ctype,302)ipm                                        3d30s22
  302         format(' Q',i2)                                            3d30s22
             end if                                                      3d29s22
            end if                                                      3d30s22
            do ik=1,nrootk                                              3d30s22
             jk=ik-1                                                    3d30s22
             do ib=1,nrootb                                             3d30s22
              jb=ib-1                                                   3d30s22
              iad=itmp+jb+nrootb*jk                                     3d30s22
              if(abs(bc(iad)).gt.1d-10)then                             3d30s22
               if(inz.eq.0)then                                          3d30s22
                if(ipass.eq.1)then                                      3d30s22
                 write(6,*)('(=-mux/sqrt(2))')
                else if(ipass.eq.2)then                                 3d30s22
                 write(6,*)('(=-Ly/(i*sqrt(2))')
                end if                                                  3d30s22
               end if                                                   3d30s22
               inz=inz+1
               write(6,333)ib,iwavebra(1,1),                            3d30s22
     $             wnameb,ctype,ik,iwaveket(1,1),wnamek,bc(iad),        3d30s22
     $             crori(iok)                                           3d30s22
              end if                                                    3d30s22
             end do                                                     3d30s22
            end do                                                      3d30s22
           end if                                                       3d29s22
          end if                                                        3d17s22
          if(isymneed.eq.iosym(muiu).and.ipm.lt.0)then                  3d17s22
           ihit=ihit+1                                                  3d17s22
           call psioppsi(iwavebra(1,1),nrootb,muiu,1d0,ixmt,iosym,      3d17s22
     $     0,idum,dum,iwaveket(1,1),nrootk,bc(itmp),mdon,               3d17s22
     $     mdoop,ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,          8d23s21
     $     nvirt,maxbx,maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,          11d15s21
     $     npadddi,bc,ibc)                                              11d10s22
           if(xnuc.ne.0d0)then                                          3d17s22
            do i=0,nrootb-1                                             3d17s22
             iad=itmp+i*(nrootb+1)                                      3d17s22
             bc(iad)=bc(iad)+xnuc                                       3d17s22
            end do                                                      3d17s22
           end if                                                       3d17s22
           do i=0,nbk-1                                                 8d29s22
            bc(itmp+i)=bc(itmp+i)*fmi                                   3d17s22
           end do                                                       3d17s22
           if(lprint)then                                               3d29s22
            lok=iordb.ne.iordk
            inz=0
            iok=2                                                       3d30s22
            if(lok)iok=1                                                3d30s22
            if(nsymb.eq.1)then                                          3d30s22
             if(ipass.eq.1)then                                          3d29s22
              ctype='Im-1'                                               3d30s22
             else if(ipass.eq.2)then                                     3d29s22
              ctype='IL-1'                                               3d30s22
             else                                                        3d29s22
              write(ctype,304)ipm                                        3d30s22
  304         format('IQ',i2)                                            3d30s22
             end if                                                      3d29s22
            else                                                        3d30s22
             if(ipass.eq.1)then                                          3d29s22
              ctype='mu-1'                                               3d30s22
             else if(ipass.eq.2)then                                     3d29s22
              ctype='L-1 '                                               3d30s22
             else                                                        3d29s22
              write(ctype,302)ipm                                        3d30s22
             end if                                                      3d29s22
            end if                                                      3d30s22
            do ik=1,nrootk                                              3d30s22
             jk=ik-1                                                    3d30s22
             do ib=1,nrootb                                             3d30s22
              jb=ib-1                                                   3d30s22
              iad=itmp+jb+nrootb*jk                                     3d30s22
              if(abs(bc(iad)).gt.1d-10)then                             3d30s22
               if(inz.eq.0)then                                          3d30s22
                if(ipass.eq.1)then                                      3d30s22
                 write(6,*)('(=-muy/sqrt(2))')
                else if(ipass.eq.2)then                                 3d30s22
                 write(6,*)('(=-Lx/(i*sqrt(2))')
                end if                                                  3d30s22
               end if                                                   3d30s22
               inz=inz+1
               write(6,333)ib,iwavebra(1,1),                            3d30s22
     $             wnameb,ctype,ik,iwaveket(1,1),wnamek,bc(iad),        3d30s22
     $             crori(iok)                                           3d30s22
              end if                                                    3d30s22
             end do                                                     3d30s22
            end do                                                      3d30s22
           end if                                                       3d29s22
          end if                                                        3d17s22
         end if                                                         3d17s22
        end do                                                          3d17s22
       end do                                                           3d17s22
       ibcoff=ifull                                                     3d17s22
      end if                                                            5d24s21
      ibcoff=ibcoffo                                                    5d26s22
      return                                                            5d20s21
      end                                                               5d20s21
