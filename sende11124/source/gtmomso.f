c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine gtmomso(iwfbra,ioffb,iwfket,ioffk,ijjs,ijjm,nspc,nec,  2d17s22
     $     multh,irefo,ixmt,nh0,nopso,iopso,irel,ism,norb,mdon,esf,     2d18s22
     $     nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,irw2,             2d17s22
     $     isopt,nsopt,iifmx,                                           2d17s22
     $     ntype,npadddi,nbasp,nbaspc,natom,ngaus,ibdat,iapair,ibstor,  12d20s20
     $         isstor,isym,ascale,idorel,iorb,opname,ipbra,ipket,       3d1s22
     $     jjlowb,jjlowk,tm,ndim,joffb,joffk,ndimb,tmb,lprint,l2sub,bc, 11d9s22
     $     ibc,n4vso)                                                   2d8s23
c
c     generate matrix elements of spin-orbit transition moments
c
      implicit real*8 (a-h,o-z)                                         5d18s21
      logical ldiag,lprint,logb,logk,nlogk                              3d16s21
      character*6 labb,labk                                             10d27s20
      character*10 opname(*)                                            2d25s22
      character*2 tkind
      character*1 crori(2)                                              3d30s22
      character*19 header                                               1d24s23
      integer*1 ipack1(4)                                               5d18s21
      integer*2 ipack2(4)                                               8d30s21
      integer*8 ioffb(*),ioffk(*),ijjs(*),ijjm(*),ipack8,joffb(*),      3d1s22
     $     joffk(*)                                                     3d1s22
      equivalence (ipack1,npack4)                                       5d18s21
      equivalence (ipack8,ipack2)                                       8d30s21
      dimension iwfbra(nspc,*),iwfket(nspc,*),multh(8,8),irefo(*),      5d18s21
     $     irel(*),ism(*),esf(*),itypex(4),ixmt(8,*),iopso(5,*),        2d18s22
     $     ntype(4),nbasp(*),nbaspc(*),nh0(*),isopt(4,4),isoptm(4),     2d25s22
     $     mur(3),mui(3),fmur(3),fmui(3),trial(3),fmat(3,3),ttmp(3),    3d1s22
     $     ffmat(3,3),tm(ndim,ndimb,2),tmb(ndimb,*),mu0(3),mur2(3),     1d13s23
     $     mui2(3)                                                      3d8s22
      include "common.store"                                            5d18s21
      common/singcm/iuse,nff
      data crori/' ','i'/                                               3d30s22
      data icall/0/
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
      nlzz=ipack1(2)                                                    2d25s22
      llb=ipack1(3)                                                     5d18s21
      npack4=iwfket(6,1)                                                5d18s21
      ldiag=loc(iwfbra).eq.loc(iwfket)
      llk=ipack1(3)                                                     5d18s21
      i2sb=iwfbra(1,1)-1                                                5d18s21
      i2sk=iwfket(1,1)-1                                                5d18s21
      if(lprint)then                                                    3d2s22
       write(6,*)('in gtmomso for '),iwfbra(1,1),labb,('and '),
     $     iwfket(1,1),labk                                             6d2s21
      end if                                                            3d2s22
      nrootb=iwfbra(3,1)                                                5d18s21
      nrootk=iwfket(3,1)                                                5d18s21
      llb2=llb*2                                                        5d18s21
      nllb=llb2+1                                                       5d18s21
      nsob=nllb*iwfbra(1,1)                                             5d18s21
      llk2=llk*2                                                        5d18s21
      nllk=llk2+1                                                       5d18s21
      nsok=nllk*iwfket(1,1)                                             5d18s21
      nroot2=nrootb*nrootk                                              5d18s21
      jlowb=iabs(llb2-i2sb)                                             3d2s22
      jhib=llb2+i2sb                                                    3d2s22
      jlowk=iabs(llk2-i2sk)                                             3d2s22
      jhik=llk2+i2sk                                                    3d2s22
      nn=nsob*nsok                                                        5d17s21
      if(ipbra.eq.ipket)then                                            2d25s22
       ipbot=2                                                          2d25s22
       iptop=3                                                          2d25s22
      else                                                              2d25s22
       ipbot=1                                                          2d25s22
       iptop=1                                                          2d25s22
      end if                                                            2d25s22
      iaout=ibcoff                                                      2d17s22
      ibcoff=iaout+2*nroot2                                             2d25s22
      iaouti=iaout+nroot2                                               2d25s22
      itmpr=ibcoff                                                      2d25s22
      itmpi=itmpr+nroot2                                                2d25s22
      ibcoff=itmpi+nroot2                                               2d25s22
      idata=ibcoff                                                      3d1s22
      ibcoff=idata+3*nroot2                                             3d1s22
      call enough('gtmomso.  1',bc,ibc)
      if(nlzz.eq.6)then                                                 2d25s22
       if(lprint)write(6,*)('spherical system ')                        3d2s22
       ncase1=0                                                         2d28s22
       ncase2=0                                                         2d28s22
       ncase3=0                                                         2d28s22
       do ipass=ipbot,iptop                                             2d25s22
        iheadr=0                                                        1d24s23
        if(ipass.eq.1)then                                              5d24s21
         header='electric dipole   '                                    1d24s23
         tkind='e1'
         iq=1                                                           5d24s21
         irori1=0                                                       2d25s22
         irori2=0                                                       2d25s22
c        +1                 -1
c     (-x-iy)/sqrt(2), (+x-iy)/sqrt(2)
c
c     i.e. real +1=-x/sqrt(2), imag +1=-y/sqrt(2)                       5d24s21
c          real -1=+x/sqrt(2), imag -1= same
         ir=0                                                           2d25s22
         ii=0                                                           2d25s22
         do isop=1,nopso                                                2d25s22
          i=iopso(2,isop)                                               2d25s22
          j=iopso(3,isop)                                               2d25s22
          if(opname(i)(1:3).eq.'mux')then                               2d25s22
           ir=ir+1                                                      2d25s22
           mur(ir)=isop                                                 2d25s22
          end if                                                        2d25s22
          if(opname(i)(1:3).eq.'muy')then                               2d25s22
           ii=ii+1                                                      2d25s22
           mui(ii)=isop                                                 2d25s22
          end if                                                        2d25s22
         end do
         ffr=-srh                                                       5d24s21
         ffi=-srh                                                       5d24s21
        else if(ipass.eq.2)then                                         5d24s21
         header='magnetic dipole   '                                    1d24s23
         tkind='m1'
         irori1=1                                                       2d25s22
         irori2=1                                                       2d25s22
         iq=1                                                           5d24s21
c        +1                    -1
c     (-lx-ily)/sqrt(2),(lx-ily)/sqrt(2)
c     =(+ilx(i)-ly(i))/sqrt(2),(-ilx(i)-ly(i))/sqrt(2)
c
c     i.e. real +1 = -ly(i)/sqrt(2), imag +1=lx(i)/sqrt(2)
c          real -1 = same          , imag -1=-lx(i)/sqrt(2)
         ir=0                                                           2d25s22
         ii=0                                                           2d25s22
         do isop=1,nopso                                                2d25s22
          i=iopso(2,isop)                                               2d25s22
          j=iopso(3,isop)                                               2d25s22
          if(opname(i)(1:4).eq.'ly/i')then                              2d25s22
           ir=ir+1                                                      2d25s22
           mur(ir)=isop                                                 2d25s22
          end if                                                        2d25s22
          if(opname(i)(1:4).eq.'lx/i')then                              2d25s22
           ii=ii+1                                                      2d25s22
           mui(ii)=isop                                                 2d25s22
          end if                                                        2d25s22
         end do
         ffr=-srh                                                       5d24s21
         ffi=srh                                                        5d24s21
        else                                                            5d24s21
         header='electric quadrapole'                                   1d24s23
         tkind='e2'
         iq=2                                                           5d24s21
         irori1=0                                                       2d25s22
         irori2=0                                                       2d26s22
c        +1     -1
c     -zx-izy,zx-izy
c
c     i.e. real +1 = -zx=Re Q1, imag +1 = -zy=Im Q1
c          real -1 = +zx, imag -1 = same
         ir=0                                                           2d25s22
         ii=0                                                           2d25s22
         do isop=1,nopso                                                2d25s22
          i=iopso(2,isop)                                               2d25s22
          j=iopso(3,isop)                                               2d25s22
          if(opname(i)(1:5).eq.'Re Q1')then                             2d25s22
           ir=ir+1                                                      2d25s22
           mur(ir)=isop                                                 2d25s22
          end if                                                        2d25s22
          if(opname(i)(1:5).eq.'Im Q1')then                             2d25s22
           ii=ii+1                                                      2d25s22
           mui(ii)=isop                                                 2d25s22
          end if                                                        2d25s22
         end do
         ffr=1d0                                                        5d26s21
         ffi=1d0                                                        5d26s21
        end if                                                          5d24s21
        iq2=iq*2
        ntrial=0                                                        3d1s22
        jdata=idata                                                     3d1s22
        do mlb=0,llb                                                    3d1s22
         mlb2=mlb*2                                                     3d1s22
         do ipm=-1,1,2
          ffru=ffr                                                       2d26s22
          ffiu=ffi                                                       2d26s22
          if(ipm.lt.0)then                                               2d26s22
           if(ipass.eq.2)then                                            2d26s22
            ffiu=-ffiu                                                   2d26s22
           else                                                          2d26s22
            ffru=-ffru                                                   2d26s22
           end if                                                        2d26s22
          end if                                                         2d26s22
          ipm2=ipm*2                                                    3d1s22
          do mq=-1,1                                                    3d1s22
           mq2=mq*2                                                     3d1s22
           if(mq.eq.0)then                                               2d25s22
            itbot=1
            ittop=1                                                      2d25s22
           else                                                          2d25s22
            itbot=2                                                      2d25s22
            ittop=3                                                      2d25s22
           end if                                                        2d25s22
           itryn=-mq2-ipm2                                                3d1s22
           itrynh=itryn/2
           mlk=mlb+itrynh
           if(iabs(mlk).le.llk)then                                     3d1s22
            mlk2=mlk*2                                                  3d1s22
            do m2sb=-i2sb,i2sb,2                                           2d25s22
             m2sk=m2sb+mq2                                                2d25s22
             if(iabs(m2sk).le.i2sk)then                                   2d25s22
              sc=cleb2(i2sb,m2sb,2,mq2,i2sk,m2sk)                         2d25s22
              if(mq.ne.0)sc=-sc                                           2d25s22
              do isl=iabs(1-iq),1+iq                                     2d25s22
               isl2=isl*2
               jsl=isl-(iq-1)+1                                         3d1s22
               trial(jsl)=0d0                                           3d1s22
               if(isl2.ge.iabs(itryn))then                              3d1s22
                csm=cleb2(2,-mq2,iq2,-ipm2,isl2,itryn)                     2d25s22
                cll=cleb2(llb2,mlb2,isl2,itryn,llk2,mlk2)
                trial(jsl)=csm*cll*sc                                   3d1s22
               end if                                                   3d1s22
              end do                                                    3d1s22
              if(ntrial.eq.0)then                                        3d1s22
               sz=sqrt(trial(1)**2+trial(2)**2+trial(3)**2)              3d1s22
               if(sz.gt.1d-7)then                                        3d1s22
                nkeep=1                                                   3d1s22
                xnorm=1d0/sz                                              3d1s22
                do i=1,3                                                  3d1s22
                 ffmat(1,i)=xnorm*trial(i)                                3d1s22
                end do                                                    3d1s22
               else                                                      3d1s22
                nkeep=0                                                  3d1s22
               end if                                                    3d1s22
              else                                                       3d1s22
               do i=1,3                                                  3d1s22
                ttmp(i)=trial(i)                                         3d1s22
               end do                                                    3d1s22
               do i=1,ntrial                                             3d1s22
                dot=0d0                                                  3d1s22
                do j=1,3                                                 3d1s22
                 dot=dot+ffmat(i,j)*ttmp(j)                              3d1s22
                end do                                                   3d1s22
                do j=1,3                                                 3d1s22
                 ttmp(j)=ttmp(j)-dot*ffmat(i,j)                           3d1s22
                end do                                                   3d1s22
               end do                                                    3d1s22
               sz=sqrt(ttmp(1)**2+ttmp(2)**2+ttmp(3)**2)                 3d1s22
               if(sz.gt.1d-7)then                                        3d1s22
                nkeep=1                                                  3d1s22
                xnorm=1d0/sz                                             3d1s22
                do i=1,3                                                 3d1s22
                 ffmat(ntrial+1,i)=ttmp(i)*xnorm                         3d1s22
                end do                                                   3d1s22
               else                                                      3d1s22
                nkeep=0                                                  3d1s22
               end if                                                    3d1s22
              end if                                                     3d1s22
              if(nkeep.eq.1)then                                         3d1s22
               ntrial=ntrial+1                                           3d1s22
               do i=1,3                                                  3d1s22
                fmat(ntrial,i)=trial(i)                                  3d1s22
               end do                                                    3d1s22
               do iz=0,nroot2-1                                          3d1s22
                bc(jdata+iz*3)=0d0                                       3d1s22
               end do                                                     2d25s22
               do it=itbot,ittop                                          2d25s22
                itp=it+1-l2sub                                          3d2s22
                i=iopso(2,mur(it))                                        2d25s22
                j=iopso(3,mur(it))                                        2d25s22
                isoptm(1)=iopso(1,mur(it))                                2d25s22
                isoptm(2)=isopt(2,j)                                     3d1s22
                isoptm(3)=isopt(3,j)                                      2d25s22
                isoptm(4)=isopt(4,j)                                      2d25s22
                call genhsoa2(iwfbra,iwfket,nspc,mlb2,mlk2,m2sb,m2sk,
     $          bc(iaout),nec,multh,irefo,ih0a,i4so2,irel,ism,norb,mdon,5d19s21
     $          iri,itypex,nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,  9d20s21
     $     irw2,ixmt(1,mur(it)),nh0,isoptm,nsopt,i4or,ionexr,jmatsr,
     $       kmatsr,kmatsrb,i3xr,iifmx,ntype,npadddi,nbasp,nbaspc,natom,   12d20s20
     $           ngaus,ibdat,iapair,ibstor,isstor,isym,ascale,idorel,   12d23s21
     $           iorb,itp,iopso(1,mur(it)),isoptm(2),bc,ibc,n4vso,1,0)  3d27s24
                do iz=0,nroot2-1                                         3d1s22
                 bc(jdata+iz*3)=bc(jdata+iz*3)+bc(iaout+iz)*ffru         3d1s22
                end do                                                    2d25s22
                i=iopso(2,mui(it))                                        2d25s22
                j=iopso(3,mui(it))                                        2d25s22
                isoptm(1)=iopso(1,mui(it))                                2d25s22
                isoptm(2)=mod(1+isopt(2,j),2)                            3d1s22
                isoptm(3)=isopt(3,j)                                      2d25s22
                isoptm(4)=isopt(4,j)                                      2d25s22
                call genhsoa2(iwfbra,iwfket,nspc,mlb2,mlk2,m2sb,m2sk,
     $          bc(iaout),nec,multh,irefo,ih0a,i4so2,irel,ism,norb,mdon,5d19s21
     $          iri,itypex,nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,  9d20s21
     $     irw2,ixmt(1,mui(it)),nh0,isoptm,nsopt,i4or,ionexr,jmatsr,
     $       kmatsr,kmatsrb,i3xr,iifmx,ntype,npadddi,nbasp,nbaspc,natom,   12d20s20
     $           ngaus,ibdat,iapair,ibstor,isstor,isym,ascale,idorel,   12d23s21
     $           iorb,itp,iopso(1,mui(it)),isoptm(2),bc,ibc,n4vso,1,0)  3d27s24
                ffu=ffiu                                                   2d25s22
                if(isopt(2,j).ne.0)ffu=-ffu                              3d1s22
                do iz=0,nroot2-1                                          2d25s22
                 bc(jdata+iz*3)=bc(jdata+iz*3)+bc(iaout+iz)*ffu          3d1s22
                end do                                                    2d25s22
               end do                                                     2d25s22
               jdata=jdata+1                                             3d1s22
               if(ntrial.eq.3)then                                       3d1s22
                iwgt=ibcoff                                              3d1s22
                icoef=iwgt+3                                             3d1s22
                iscr=icoef+3*nroot2                                      3d1s22
                ibcoff=iscr+6+9+18+3                                     3d1s22
                call enough('gtmomso.  2',bc,ibc)
                iuse=0                                                   3d1s22
                do i=0,2                                                 3d1s22
                 bc(iwgt+i)=1d0                                          3d1s22
                end do                                                   3d1s22
                call lsqfit2(fmat,3,3,bc(idata),3,nroot2,3,bc(icoef),    3d1s22
     $              3,bc(iscr),bc(iwgt),0,rms,bc,ibc)                   3d27s23
                go to 27                                                 3d1s22
               end if                                                    3d1s22
              end if
             end if                                                     3d1s22
            end do
           end if                                                       3d1s22
          end do                                                        3d1s22
         end do                                                         3d1s22
        end do                                                          3d1s22
        if(ntrial.gt.0)then                                             3d1s22
         iwgt=ibcoff                                                     3d1s22
         icoef=iwgt+3                                                    3d1s22
         iscr=icoef+3*nroot2                                             3d1s22
         ibcoff=iscr+6+9+18+3                                            3d1s22
         icoef2=ibcoff                                                  3d2s22
         ibcoff=icoef2+ntrial*nroot2                                    3d2s22
         call enough('gtmomso.  3',bc,ibc)
         iuse=0                                                          3d1s22
         do i=0,2                                                        3d1s22
          bc(iwgt+i)=1d0                                                 3d1s22
         end do                                                          3d1s22
         ii=0                                                           3d2s22
         do i=1,3                                                       3d2s22
          sz=0d0                                                        3d2s22
          do j=1,ntrial                                                 3d2s22
           sz=sz+fmat(j,i)**2                                           3d2s22
          end do                                                        3d2s22
          sz=sqrt(sz/dfloat(ntrial))                                    3d2s22
          if(sz.gt.1d-7)then                                            3d2s22
           ii=ii+1                                                      3d2s22
           do j=1,ntrial                                                3d2s22
            ffmat(j,ii)=fmat(j,i)                                       3d2s22
           end do                                                       3d2s22
          end if                                                        3d2s22
         end do                                                         3d2s22
         call lsqfit2(ffmat,3,ntrial,bc(idata),3,nroot2,ntrial,         3d2s22
     $        bc(icoef2),3,bc(iscr),bc(iwgt),0,rms,bc,ibc)              3d27s23
         ii=0                                                           3d2s22
         do i=1,3                                                       3d2s22
          sz=0d0                                                        3d2s22
          do j=1,ntrial                                                 3d2s22
           sz=sz+fmat(j,i)**2                                           3d2s22
          end do                                                        3d2s22
          sz=sqrt(sz/dfloat(ntrial))                                    3d2s22
          if(sz.gt.1d-7)then                                            3d2s22
           iad1=icoef2+ii                                               3d2s22
           iad2=icoef+i-1                                               3d2s22
           do j=0,nroot2-1                                              3d2s22
            bc(iad2+j*3)=bc(iad1+j*3)
           end do
           ii=ii+1                                                      3d2s22
          else                                                          3d2s22
           iad2=icoef+i-1                                               3d2s22
           do j=0,nroot2-1                                              3d2s22
            bc(iad2+j*3)=0d0                                            3d2s22
           end do
          end if                                                        3d2s22
         end do                                                         3d2s22
        else                                                            3d1s22
        end if                                                          3d1s22
   27   continue                                                        3d1s22
        if(ntrial.gt.0)then                                             3d1s22
         if(lprint)then                                                 3d2s22
          do isl=iabs(1-iq),1+iq                                         3d2s22
           jsl=icoef+isl-(iq-1)                                          3d2s22
           do ik=1,nrootk                                                3d2s22
            do ib=1,nrootb                                               3d2s22
             iadt=jsl+3*(ib-1+nrootb*(ik-1))                             3d2s22
             if(abs(bc(iadt)).gt.1d-5)then                               3d2s22
              if(iheadr.eq.0)write(6,*)header                           1d24s23
              iheadr=1                                                  1d24s23
              write(6,28)ib,iwfbra(1,1),labb,tkind,isl,ik,iwfket(1,1),   3d2s22
     $           labk,bc(iadt),bc(iadt)*ascale                          1d24s23
   28         format('<',i3,1x,i1,a6,'|',a2,x,i1,x,'so |',i3,1x,i1,a6,   3d2s22
     $            '>=',f12.6,3x,'times 1/4c^2',2x,es13.6)               1d24s23
             end if                                                      3d2s22
            end do                                                       3d2s22
           end do                                                        3d2s22
          end do                                                         3d2s22
         end if                                                         3d2s22
         if(ipass.le.2)then                                             3d1s22
c     electric dipole or magnetic dipole
          isto=1                                                        3d1s22
         else                                                           3d1s22
          isto=2                                                        3d1s22
         end if                                                         3d1s22
         do jk=jlowk,jhik,2                                               6d4s21
          jjk=1+((jk-jjlowk)/2)                                            6d4s21
          iket=joffk(jjk)+ioffk(jjk)+1                                   6d4s21
          do jb=jlowb,jhib,2                                              6d4s21
           jjb=1+((jb-jjlowb)/2)                                           6d4s21
           do isl=iabs(1-iq),1+iq                                       3d1s22
            isl2=isl*2                                                  3d1s22
            jsl=isl-(iq-1)                                              3d1s22
            fact=f9j(jb,iq2,jk,i2sb,2,i2sk,llb2,isl2,llk2,1)            3d2s22
            fact=fact*sqrt(dfloat(jk+1)*dfloat(jb+1)*dfloat(llk2+1)     3d1s22
     $           *dfloat(i2sk+1)*dfloat(isl2+1))*ascale                 3d2s22
            isum=i2sk-llk2+iq2                                          3d1s22
            if(mod(isum,2).ne.0)isum=isum+1                             3d1s22
            isum=isum/2                                                 3d1s22
            if(mod(isum,2).ne.0)fact=-fact                              3d1s22
            ibra=joffb(jjb)+ioffb(jjb)+1                                  6d4s21
            do jrk=0,nrootk-1                                             6d4s21
             icol=iket+jrk
             do jrb=0,nrootb-1                                            6d4s21
              iadt=icoef+jsl+3*(jrb+nrootb*jrk)                         3d1s22
              irow=ibra+jrb
              tm(irow,icol,isto)=tm(irow,icol,isto)+bc(iadt)*fact       3d1s22
             end do                                                       6d4s21
            end do                                                        6d4s21
            if(loc(iwfbra).ne.loc(iwfket))then                          3d2s22
             fact=f9j(jk,iq2,jb,i2sk,2,i2sb,llk2,isl2,llb2,1)           3d2s22
             fact=fact*sqrt(dfloat(jk+1)*dfloat(jb+1)*dfloat(llb2+1)    3d2s22
     $           *dfloat(i2sb+1)*dfloat(isl2+1))*ascale                 3d2s22
             isum=i2sb-llb2+iq2                                         3d2s22
             if(mod(isum,2).ne.0)isum=isum+1                            3d2s22
             isum=isum/2                                                3d2s22
             if(mod(isum,2).ne.0)fact=-fact                             3d2s22
             ibra=joffb(jjb)+ioffb(jjb)+1                               3d2s22
             if(ipass.gt.1)then                                         3d2s22
c     magnetic dipole or electric quadrapole, i.e. even parity
              do jrk=0,nrootk-1                                          3d2s22
               icol=iket+jrk                                             3d2s22
               do jrb=0,nrootb-1                                         3d2s22
                iadt=icoef+jsl+3*(jrb+nrootb*jrk)                        3d2s22
                irow=ibra+jrb                                            3d2s22
                tm(icol,irow,isto)=tm(icol,irow,isto)+bc(iadt)*fact      3d2s22
               end do                                                    3d2s22
              end do                                                     3d2s22
             else                                                       3d2s22
c     electric dipole, i.e. odd parity
              do jrk=0,nrootk-1                                          3d2s22
               icol=iket+jrk                                             3d2s22
               do jrb=0,nrootb-1                                         3d2s22
                iadt=icoef+jsl+3*(jrb+nrootb*jrk)                        3d2s22
                irow=ibra+jrb                                            3d2s22
                tmb(icol,irow)=tmb(icol,irow)+bc(iadt)*fact             3d2s22
               end do                                                    3d2s22
              end do                                                     3d2s22
             end if                                                     3d2s22
            end if                                                      3d2s22
           end do                                                       3d1s22
          end do                                                        3d1s22
         end do                                                         3d1s22
        end if                                                          3d1s22
       end do                                                           2d25s22
      else if(nlzz.eq.2)then                                            2d25s22
       iaout=ibcoff                                                     3d8s22
       ires=iaout+nroot2*2                                              3d8s22
       ibcoff=ires+nroot2                                               3d8s22
       call enough('gtmomso.  4',bc,ibc)
       mlba=llb                                                         3d8s22
       if(llb.eq.0)then                                                 3d8s22
        nllb=1                                                          3d8s22
       else                                                             3d8s22
        nllb=2                                                          3d8s22
       end if                                                           3d8s22
       mlka=llk                                                         3d8s22
       if(llk.eq.0)then                                                 3d9s22
        nllk=1                                                          3d9s22
       else                                                             3d9s22
        nllk=2                                                          3d9s22
       end if                                                           3d9s22
       do ipass=1,3                                                     5d26s21
        iheadr=0                                                        1d24s23
        if(ipass.eq.1)then                                              5d24s21
         header='electric dipole   '                                    1d24s23
         tkind='e1'                                                        3d8s22
         iq=1                                                           5d24s21
c        +1                 -1               0
c     (-x-iy)/sqrt(2), (+x-iy)/sqrt(2)       z
c
c     i.e. real +1=-x/sqrt(2), imag +1=-y/sqrt(2)                       5d24s21
c          real -1=+x/sqrt(2), imag -1= same
         ir=0                                                           2d25s22
         ii=0                                                           2d25s22
         i0=0                                                           3d8s22
         do isop=1,nopso                                                3d8s22
          i=iopso(2,isop)                                               3d8s22
          if(opname(i)(1:3).eq.'mux')then
           ir=ir+1                                                      3d8s22
           mur(ir)=isop                                                 3d8s22
          end if                                                        3d8s22
          if(opname(i)(1:3).eq.'muy')then
           ii=ii+1                                                      3d8s22
           mui(ii)=isop                                                 3d8s22
          end if                                                        3d8s22
          if(opname(i)(1:3).eq.'mu0')then
           i0=i0+1                                                      3d8s22
           mu0(i0)=isop                                                 3d8s22
          end if                                                        3d8s22
         end do
         ffr=-srh                                                       5d24s21
         ffi=-srh                                                       5d24s21
        else if(ipass.eq.2)then                                         5d24s21
         header='magnetic dipole   '                                    1d24s23
         tkind='m1'                                                        3d8s22
         iq=1                                                           5d24s21
c        +1                    -1                0
c     (-lx-ily)/sqrt(2),(lx-ily)/sqrt(2)         lz=i*lz(i)
c     =(+ilx(i)-ly(i))/sqrt(2),(-ilx(i)-ly(i))/sqrt(2)
c
c     i.e. real +1 = -ly(i)/sqrt(2), imag +1=lx(i)/sqrt(2)
c          real -1 = same          , imag -1=-lx(i)/sqrt(2)
         ir=0                                                           3d8s22
         ii=0                                                           3d8s22
         i0=0                                                           3d8s22
         do isop=1,nopso                                                2d25s22
          i=iopso(2,isop)                                               2d25s22
          if(opname(i)(1:4).eq.'lx/i')then
           ii=ii+1                                                      2d25s22
           mui(ii)=isop                                                 2d25s22
          end if                                                        3d8s22
          if(opname(i)(1:4).eq.'ly/i')then                              3d8s22
           ir=ir+1                                                      2d25s22
           mur(ir)=isop                                                 2d25s22
          end if                                                        3d8s22
          if(opname(i)(1:4).eq.'lz/i')then
           i0=i0+1                                                      3d8s22
           mu0(i0)=isop                                                 3d8s22
          end if                                                        3d8s22
         end do
         ffr=-srh                                                       5d24s21
         ffi=srh                                                        5d24s21
        else                                                            5d24s21
         header='electric quadrapole'                                   1d24s23
         iq=2                                                           5d24s21
         tkind='e2'                                                        3d8s22
c        +1     -1    0   +/- 2
c     -zx-izy,zx-izy q0  Re Q2 +/- i Im Q2
c                          xx-yy      xy
c     i.e. real +1 = -zx, imag +1 = -zy
c          real -1 = +zx, imag -1 = same
         i0=0                                                           3d8s22
         ir=0                                                           3d8s22
         ii=0                                                           3d8s22
         irr=0                                                          3d8s22
         iii=0                                                          3d8s22
         do isop=1,nopso                                                2d25s22
          i=iopso(2,isop)                                               2d25s22
          if(opname(i)(1:5).eq.'Re Q1')then                             3d8s22
           ir=ir+1                                                      3d8s22
           mur(ir)=isop                                                 3d8s22
          end if                                                        3d8s22
          if(opname(i)(1:5).eq.'Im Q1')then                             3d8s22
           ii=ii+1                                                      3d8s22
           mui(ii)=isop                                                 3d8s22
          end if                                                        3d8s22
          if(opname(i)(1:2).eq.'Q0')then                                3d8s22
           i0=i0+1                                                      3d8s22
           mu0(i0)=isop                                                 3d8s22
          end if                                                        3d8s22
          if(opname(i)(1:5).eq.'Re Q2')then                             3d8s22
           irr=irr+1                                                    3d8s22
           mur2(irr)=isop                                               3d8s22
          end if                                                        3d8s22
          if(opname(i)(1:5).eq.'Im Q2')then                             3d8s22
           iii=iii+1                                                    3d8s22
           mui2(iii)=isop                                               3d8s22
          end if                                                        3d8s22
         end do
         ffr=1d0                                                        5d26s21
         ffi=1d0                                                        5d26s21
        end if                                                          5d24s21
        iq2=iq*2                                                        3d17s22
        iqcm=iabs(iq-1)                                                  3d8s22
        iqcx=iq+1                                                        3d8s22
        nf=iqcx+1-iqcm                                                    3d8s22
c
c     first figure out tt vs. ff relation.
c
        phs=1d0                                                         3d17s22
        if(i2sb.ne.i2sk)then
         phs=-phs                                                       3d17s22
        end if
        if(ipass.eq.2)then                                              3d17s22
         phs=-phs                                                       3d17s22
        end if                                                          3d17s22
        if(llb.eq.0.and.labb(4:4).eq.'-')then                           3d17s22
         phs=-phs                                                       3d17s22
        end if                                                          3d17s22
        if(llk.eq.0.and.labk(4:4).eq.'-')then                           3d17s22
         phs=-phs                                                       3d17s22
        end if                                                          3d17s22
        if(llb.eq.0)then                                                3d17s22
         m2sbot=mod(i2sb,2)                                             3d17s22
        else                                                            3d17s22
         m2sbot=-i2sb                                                   3d17s22
        end if                                                          3d17s22
        do lkind=1,2                                                    3d21s22
         if(lkind.eq.1)then                                             3d21s22
         else                                                           3d21s22
         end if                                                         3d21s22
         ndat=0                                                          3d17s22
         mlb=llb                                                         3d17s22
         do m2sbu=m2sbot,i2sb,2                                          3d17s22
          m2sb=m2sbu                                                     3d17s22
          mlb2=mlb*2                                                     3d17s22
          do ipm=-iq,iq                                                  3d17s22
           ipm2=ipm*2                                                     3d8s22
           do mq=-1,1                                                     3d8s22
            mq2=mq*2                                                      3d8s22
            itryn=-mq2-ipm2                                               3d8s22
            itrynh=itryn/2                                                3d8s22
            itrynha=iabs(itrynh)                                          3d8s22
            m2sk=m2sb+mq2                                                 3d17s22
            if(iabs(m2sk).le.i2sk)then                                    3d17s22
             mlk=mlb+itrynh                                               3d8s22
             mlka=iabs(mlk)                                                5d24s21
             if(mlka.eq.llk)then                                         3d17s22
              if((lkind.eq.1.and.(mlk.gt.0.or.(mlk.eq.0.and.m2sk.ge.0)))3d21s22
     $            .or.(lkind.eq.2.and.(mlk.lt.0.or.(mlk.eq.0.and.       3d21s22
     $            m2sk.le.0))))then                                     3d21s22
               ndat=ndat+1                                                3d17s22
              end if                                                      3d17s22
             end if                                                       3d17s22
            end if                                                      3d21s22
           end do                                                        3d17s22
          end do                                                         3d17s22
         end do                                                          3d17s22
         if(ndat.eq.0)go to 1220                                        3d21s22
         nf=iqcx+1-iqcm                                                 3d21s22
         idat=ibcoff                                                     3d17s22
         ifcn=idat+nroot2*ndat*2                                         3d18s22
         nnn=nroot2*ndat                                                 3d18s22
         iwgt=ifcn+nf*ndat                                               3d17s22
         itmp=iwgt+ndat                                                  3d18s22
         ibcoff=itmp+nroot2                                              3d18s22
         call enough('gtmomso.  5',bc,ibc)
         do iz=idat,ibcoff-1                                             3d18s22
          bc(iz)=0d0                                                     3d18s22
         end do                                                          3d18s22
         jdat=idat                                                       3d17s22
         jfcn=ifcn                                                       3d17s22
         jwgt=iwgt                                                       3d17s22
         do m2sbu=m2sbot,i2sb,2                                          3d17s22
          mlb=mlba                                                       3d17s22
          m2sb=m2sbu                                                     3d17s22
          mlb2=mlb*2                                                     3d17s22
          do ipm=-iq,iq                                                  3d17s22
           ipm2=ipm*2                                                     3d8s22
           fmr=ffr                                                        3d8s22
           fmi=ffi                                                        3d8s22
           if(ipm.ge.0)then                                               5d26s21
           else if(ipm.eq.-1)then                                         5d26s21
            if(ipass.eq.2)then                                            5d26s21
             fmi=-fmi                                                     3d8s22
            else                                                          5d26s21
             fmr=-fmr                                                     3d8s22
            end if                                                        5d26s21
           else if(ipm.eq.-2)then                                         5d26s21
            fmi=-fmi                                                      3d8s22
           end if                                                         5d26s21
           do mq=-1,1                                                     3d8s22
            mq2=mq*2                                                      3d8s22
            if(mq.eq.0)then                                               3d8s22
             itbot=1                                                      3d8s22
             ittop=1                                                      3d8s22
            else                                                          3d8s22
             itbot=2                                                      3d8s22
             ittop=3                                                      3d8s22
            end if                                                        3d8s22
            itryn=-mq2-ipm2                                               3d8s22
            itrynh=itryn/2                                                3d8s22
            itrynha=iabs(itrynh)                                          3d8s22
            m2sk=m2sb+mq2                                                 3d17s22
            if(iabs(m2sk).le.i2sk)then                                    3d17s22
             sc=cleb2(i2sb,m2sb,2,mq2,i2sk,m2sk)                          3d17s22
             mlk=mlb+itrynh                                               3d8s22
             mlka=iabs(mlk)                                                5d24s21
             if(mlka.eq.llk)then                                         3d17s22
              if((lkind.eq.1.and.(mlk.gt.0.or.(mlk.eq.0.and.m2sk.ge.0)))3d21s22
     $            .or.(lkind.eq.2.and.(mlk.lt.0.or.
     $            (mlk.eq.0.and.m2sk.le.0))))then                       3d21s22
               mlk2=mlk*2                                                 3d17s22
               do iz=0,nroot2-1                                           3d17s22
                bc(itmp+iz)=0d0                                           3d18s22
               end do                                                     3d17s22
               if(ipm.eq.0)then                                             5d26s21
                do it=itbot,ittop                                          2d25s22
                 itp=it+1-l2sub                                           3d2s22
                 i=iopso(2,mu0(it))                                        2d25s22
                 j=iopso(3,mu0(it))                                        2d25s22
                 isoptm(1)=iopso(1,mu0(it))                                2d25s22
                 if(ipass.eq.2)then                                       3d8s22
                  isoptm(2)=mod(1+isopt(2,j),2)                            3d1s22
                 else                                                     3d8s22
                  isoptm(2)=isopt(2,j)                                     3d1s22
                 end if                                                   3d8s22
                 isoptm(3)=isopt(3,j)                                      2d25s22
                 isoptm(4)=isopt(4,j)                                      2d25s22
                 call genhsoa2(iwfbra,iwfket,nspc,mlb2,mlk2,m2sb,m2sk,
     $          bc(iaout),nec,multh,irefo,ih0a,i4so2,irel,ism,norb,mdon,5d19s21
     $          iri,itypex,nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,  9d20s21
     $     irw2,ixmt(1,mu0(it)),nh0,isoptm,nsopt,i4or,ionexr,jmatsr,
     $       kmatsr,kmatsrb,i3xr,iifmx,ntype,npadddi,nbasp,nbaspc,natom,   12d20s20
     $           ngaus,ibdat,iapair,ibstor,isstor,isym,ascale,idorel,   12d23s21
     $           iorb,itp,iopso(1,mu0(it)),isoptm(2),bc,ibc,n4vso,1,0)  3d27s24
                 szi=0d0                                                  3d9s22
                 do iz=0,nroot2-1                                         3d9s22
                  szi=szi+bc(iaout+iz+nroot2)**2                          3d9s22
                 end do                                                   3d9s22
                 szi=sqrt(szi/dfloat(nroot2))                             3d9s22
                 if(ipass.eq.2.and.isopt(2,j).ne.0)then                   3d8s22
                  do iz=0,nroot2-1                                         3d1s22
                   bc(itmp+iz)=bc(itmp+iz)-bc(iaout+iz)                   3d17s22
                  end do                                                    2d25s22
                 else                                                     3d8s22
                  do iz=0,nroot2-1                                         3d1s22
                   bc(itmp+iz)=bc(itmp+iz)+bc(iaout+iz)                   3d8s22
                  end do                                                    2d25s22
                 end if                                                   3d8s22
                end do                                                    3d8s22
               else if(iabs(ipm).eq.1)then                                3d9s22
                do it=itbot,ittop                                          2d25s22
                 itp=it+1-l2sub                                           3d2s22
                 i=iopso(2,mur(it))                                        2d25s22
                 j=iopso(3,mur(it))                                        2d25s22
                 isoptm(1)=iopso(1,mur(it))                                2d25s22
                 isoptm(2)=isopt(2,j)                                     3d1s22
                 isoptm(3)=isopt(3,j)                                      2d25s22
                 isoptm(4)=isopt(4,j)                                      2d25s22
                 call genhsoa2(iwfbra,iwfket,nspc,mlb2,mlk2,m2sb,m2sk,
     $          bc(iaout),nec,multh,irefo,ih0a,i4so2,irel,ism,norb,mdon,5d19s21
     $          iri,itypex,nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,  9d20s21
     $     irw2,ixmt(1,mur(it)),nh0,isoptm,nsopt,i4or,ionexr,jmatsr,
     $       kmatsr,kmatsrb,i3xr,iifmx,ntype,npadddi,nbasp,nbaspc,natom,   12d20s20
     $           ngaus,ibdat,iapair,ibstor,isstor,isym,ascale,idorel,   12d23s21
     $           iorb,itp,iopso(1,mur(it)),isoptm(2),bc,ibc,n4vso,1,0)  3d27s24
                 szi=0d0                                                  3d9s22
                 do iz=0,nroot2-1                                         3d9s22
                  szi=szi+bc(iaout+iz+nroot2)**2                          3d9s22
                 end do                                                   3d9s22
                 szi=sqrt(szi/dfloat(nroot2))                             3d9s22
                 if(szi.gt.1d-10)write(6,*)
     $               ('imaginary for mtm=+-1 real for it = '),it,szi
                 do iz=0,nroot2-1                                         3d1s22
                  bc(itmp+iz)=bc(itmp+iz)+bc(iaout+iz)*fmr                3d8s22
                 end do                                                    2d25s22
                 i=iopso(2,mui(it))                                        2d25s22
                 j=iopso(3,mui(it))                                        2d25s22
                 isoptm(1)=iopso(1,mui(it))                                2d25s22
                 isoptm(2)=mod(1+isopt(2,j),2)                            3d1s22
                 isoptm(3)=isopt(3,j)                                      2d25s22
                 isoptm(4)=isopt(4,j)                                      2d25s22
                 call genhsoa2(iwfbra,iwfket,nspc,mlb2,mlk2,m2sb,m2sk,
     $          bc(iaout),nec,multh,irefo,ih0a,i4so2,irel,ism,norb,mdon,5d19s21
     $          iri,itypex,nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,  9d20s21
     $     irw2,ixmt(1,mui(it)),nh0,isoptm,nsopt,i4or,ionexr,jmatsr,
     $       kmatsr,kmatsrb,i3xr,iifmx,ntype,npadddi,nbasp,nbaspc,natom,   12d20s20
     $           ngaus,ibdat,iapair,ibstor,isstor,isym,ascale,idorel,   12d23s21
     $           iorb,itp,iopso(1,mui(it)),isoptm(2),bc,ibc,n4vso,1,0)  3d27s24
                 ffu=fmi                                                  3d8s22
                 szi=0d0                                                  3d9s22
                 do iz=0,nroot2-1                                         3d9s22
                  szi=szi+bc(iaout+iz+nroot2)**2                          3d9s22
                 end do                                                   3d9s22
                 szi=sqrt(szi/dfloat(nroot2))                             3d9s22
                 if(szi.gt.1d-10)write(6,*)
     $               ('imaginary for mtm=+-1 imag for it = '),it,szi
                 if(isopt(2,j).ne.0)ffu=-ffu                              3d1s22
                 do iz=0,nroot2-1                                          2d25s22
                  bc(itmp+iz)=bc(itmp+iz)+bc(iaout+iz)*ffu                3d8s22
                 end do                                                    2d25s22
                end do                                                     2d25s22
               else                                                       3d8s22
                do it=itbot,ittop                                          2d25s22
                 itp=it+1-l2sub                                           3d2s22
                 i=iopso(2,mur2(it))                                        2d25s22
                 j=iopso(3,mur2(it))                                        2d25s22
                 isoptm(1)=iopso(1,mur2(it))                                2d25s22
                 isoptm(2)=isopt(2,j)                                     3d1s22
                 isoptm(3)=isopt(3,j)                                      2d25s22
                 isoptm(4)=isopt(4,j)                                      2d25s22
                 call genhsoa2(iwfbra,iwfket,nspc,mlb2,mlk2,m2sb,m2sk,
     $          bc(iaout),nec,multh,irefo,ih0a,i4so2,irel,ism,norb,mdon,5d19s21
     $          iri,itypex,nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,  9d20s21
     $     irw2,ixmt(1,mur2(it)),nh0,isoptm,nsopt,i4or,ionexr,jmatsr,
     $       kmatsr,kmatsrb,i3xr,iifmx,ntype,npadddi,nbasp,nbaspc,natom,   12d20s20
     $           ngaus,ibdat,iapair,ibstor,isstor,isym,ascale,idorel,   12d23s21
     $           iorb,itp,iopso(1,mur2(it)),isoptm(2),bc,ibc,n4vso,1,0) 3d27s24
                 szi=0d0                                                  3d9s22
                 do iz=0,nroot2-1                                         3d9s22
                  szi=szi+bc(iaout+iz+nroot2)**2                          3d9s22
                 end do                                                   3d9s22
                 szi=sqrt(szi/dfloat(nroot2))                             3d9s22
                 if(szi.gt.1d-10)write(6,*)
     $               ('imaginary for mtm=+-2 real for it = '),it,szi
                 do iz=0,nroot2-1                                         3d1s22
                  bc(itmp+iz)=bc(itmp+iz)+bc(iaout+iz)*fmr                3d8s22
                 end do                                                    2d25s22
                 i=iopso(2,mui2(it))                                      3d8s22
                 j=iopso(3,mui2(it))                                      3d8s22
                 isoptm(1)=iopso(1,mui2(it))                              3d8s22
                 isoptm(2)=mod(1+isopt(2,j),2)                            3d1s22
                 isoptm(3)=isopt(3,j)                                      2d25s22
                 isoptm(4)=isopt(4,j)                                      2d25s22
                 call genhsoa2(iwfbra,iwfket,nspc,mlb2,mlk2,m2sb,m2sk,
     $          bc(iaout),nec,multh,irefo,ih0a,i4so2,irel,ism,norb,mdon,5d19s21
     $          iri,itypex,nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,  9d20s21
     $     irw2,ixmt(1,mui2(it)),nh0,isoptm,nsopt,i4or,ionexr,jmatsr,   3d8s22
     $       kmatsr,kmatsrb,i3xr,iifmx,ntype,npadddi,nbasp,nbaspc,natom,   12d20s20
     $           ngaus,ibdat,iapair,ibstor,isstor,isym,ascale,idorel,   12d23s21
     $           iorb,itp,iopso(1,mui2(it)),isoptm(2),bc,ibc,n4vso,1,0) 3d27s24
                 szi=0d0                                                  3d9s22
                 do iz=0,nroot2-1                                         3d9s22
                  szi=szi+bc(iaout+iz+nroot2)**2                          3d9s22
                 end do                                                   3d9s22
                 szi=sqrt(szi/dfloat(nroot2))                             3d9s22
                 if(szi.gt.1d-10)write(6,*)
     $              ('imaginary for mtm=+-1 imag for it = '),it,szi
                 ffu=fmi                                                  3d8s22
                 if(isopt(2,j).ne.0)ffu=-ffu                              3d1s22
                 do iz=0,nroot2-1                                          2d25s22
                  bc(itmp+iz)=bc(itmp+iz)+bc(iaout+iz)*ffu                3d8s22
                 end do                                                    2d25s22
                end do                                                     2d25s22
               end if                                                     3d8s22
               bc(jwgt)=1d0                                               3d17s22
               logb=.false.
               do iqc=iqcm,iqcx                                          3d8s22
                jqc=iqc-iqcm                                             3d8s22
                bc(jfcn+ndat*jqc)=0d0                                    3d17s22
                if(iqc.ge.itrynha)then                                   3d8s22
                 iqc2=iqc*2                                              3d8s22
                 itry=mlk2-mlb2
                 qq=cleb2(2,-mq2,iq2,-ipm2,iqc2,itryn)                   3d8s22
                 bc(jfcn+ndat*jqc)=qq*sc                                 3d17s22
                end if                                                   3d8s22
               end do                                                    3d8s22
               do iz=0,nroot2-1                                          3d18s22
                bc(jdat+ndat*iz)=bc(itmp+iz)                             3d18s22
               end do                                                    3d18s22
               sz=0d0                                                     3d17s22
               do iz=0,nroot2-1                                           3d17s22
                sz=sz+bc(itmp+iz)**2                                      3d17s22
               end do                                                     3d17s22
               sz=sqrt(sz/dfloat(nroot2))
               jdat=jdat+1                                                3d17s22
               jfcn=jfcn+1                                                3d17s22
               jwgt=jwgt+1                                                3d17s22
              end if                                                      3d17s22
             end if                                                       3d17s22
            end if                                                      3d21s22
           end do                                                        3d17s22
          end do                                                         3d17s22
         end do                                                          3d17s22
         icoef=ibcoff                                                    3d17s22
         iscr=icoef+nf*nroot2                                            3d18s22
         ibcoff=iscr+2*nf+nf*nf+2*nf*ndat+ndat                           3d17s22
         call enough('gtmomso.  6',bc,ibc)
         nkeep=-1                                                       3d21s22
         if(nf.gt.0)then                                                3d21s22
          jfcn=ifcn                                                     3d21s22
          kfcn=ifcn                                                     3d21s22
          do i=0,ndat*nf-1                                              3d21s22
           bc(iscr+i)=bc(ifcn+i)                                        3d21s22
          end do                                                        3d21s22
          ikeep=ibcoff                                                  3d21s22
          ibcoff=ikeep+nf                                               3d21s22
          nkeep=0                                                       3d21s22
          do i=0,nf-1                                                   3d21s22
           ix=iscr+ndat*i                                               3d21s22
           do j=0,i-1                                                   3d21s22
            jx=iscr+ndat*j                                              3d21s22
            dot=0d0                                                     3d21s22
            do k=0,ndat-1                                               3d21s22
             dot=dot+bc(jx+k)*bc(ix+k)                                  3d21s22
            end do                                                      3d21s22
            do k=0,ndat-1                                               3d21s22
             bc(ix+k)=bc(ix+k)-dot*bc(jx+k)                             3d21s22
            end do                                                      3d21s22
           end do                                                       3d21s22
           dot=0d0                                                      3d21s22
           do k=0,ndat-1                                                3d21s22
            dot=dot+bc(ix+k)**2                                         3d21s22
           end do                                                       3d21s22
           dot=sqrt(dot)                                                3d21s22
           if(dot.gt.1d-10)then                                         3d21s22
            ibc(ikeep+nkeep)=i                                          3d21s22
            nkeep=nkeep+1                                               3d21s22
            dot=1d0/dot                                                 3d21s22
            do k=0,ndat-1                                               3d21s22
             bc(ix+k)=bc(ix+k)*dot                                      3d21s22
            end do                                                      3d21s22
            do k=0,ndat-1                                               3d21s22
             bc(kfcn+k)=bc(jfcn+k)                                      3d21s22
            end do                                                      3d21s22
            kfcn=kfcn+ndat                                              3d21s22
           end if                                                       3d21s22
           jfcn=jfcn+ndat                                               3d21s22
          end do                                                        3d21s22
          nf=nkeep                                                      3d21s22
         end if                                                         3d21s22
         if(nf.gt.ndat)nf=ndat
         iuse=0                                                          3d17s22
         if(nkeep.gt.0.and.lprint)then                                  2d20s24
          call lsqfit2(bc(ifcn),ndat,nf,bc(idat),ndat,nroot2,ndat,        3d18s22
     $       bc(icoef),nf,bc(iscr),bc(iwgt),0,rms,bc,ibc)               3d27s23
          do ii=0,nf-1                                                  3d23s22
           ibc(ikeep+ii)=ibc(ikeep+ii)+iqcm                             3d23s22
          end do                                                        3d23s22
          do ib=0,nrootb-1                                              3d23s22
           jb=ib+1                                                      3d24s22
           do ik=0,nrootk-1                                             3d23s22
            jk=ik+1                                                     3d24s22
            jcoef=icoef+nf*(ib+nrootb*ik)                               3d23s22
            do ii=0,nf-1                                                3d23s22
             if(abs(bc(jcoef+ii)).gt.1d-10)then                         3d23s22
              if(iheadr.eq.0)write(6,*)header                           1d24s23
              iheadr=1                                                  1d24s23
              write(6,2102)jb,iwfbra(1,1),labb,tkind,ibc(ikeep+ii),     3d23s22
     $            lkind,jk,iwfket(1,1),labk,bc(jcoef+ii),               3d23s22
     $            bc(jcoef+ii)*ascale                                   3d23s22
 2102         format('<',i3,x,i1,a6,'|',a2,x,i1,i2,' so |',i3,x,i1,a6,  3d24s22
     $          '>=',f12.6,es15.7)                                      3d23s22
             end if                                                     3d23s22
            end do                                                      3d23s22
           end do                                                       3d23s22
          end do                                                        3d23s22
         end if                                                         3d21s22
         if(rms.gt.1d-10)then                                            3d17s22
          write(6,*)('bad fit!! ')
         end if                                                          3d17s22
         ibcoff=idat                                                     3d17s22
 1220    continue                                                       3d21s22
        end do                                                          3d21s22
 1059        format(x,2l1,x,'cleb2: 2s''s':3i3,x,'2ms''s',3i3,x,        3d14s22
     $            '2lam''s',2i3,4f10.6)                                 3d14s22
       end do                                                           3d8s22
      else                                                              2d25s22
       iaout=ibcoff                                                     3d28s22
       itmp=iaout+nroot2*2                                              3d28s22
       itmp2=itmp+nroot2                                                3d29s22
       ibcoff=itmp2+nroot2                                              3d29s22
       if(nsymb.eq.1)then                                               3d31s22
        itmpi=ibcoff                                                    3d31s22
        itmp2i=itmpi+nroot2                                             3d31s22
        ibcoff=itmp2i+nroot2                                            3d31s22
       end if                                                           3d31s22
       call enough('gtmomso.  7',bc,ibc)
       do ipass=1,3                                                     5d26s21
        if(ipass.eq.1)then                                              5d24s21
         tkind='e1'                                                        3d8s22
         iq=1                                                           5d24s21
c        +1                 -1               0
c     (-x-iy)/sqrt(2), (+x-iy)/sqrt(2)       z
c
c     i.e. real +1=-x/sqrt(2), imag +1=-y/sqrt(2)                       5d24s21
c          real -1=+x/sqrt(2), imag -1= same
         ir=0                                                           2d25s22
         ii=0                                                           2d25s22
         i0=0                                                           3d8s22
         do isop=1,nopso                                                3d8s22
          i=iopso(2,isop)                                               3d8s22
          if(opname(i)(1:3).eq.'mux')then
           ir=ir+1                                                      3d8s22
           mur(ir)=isop                                                 3d8s22
          end if                                                        3d8s22
          if(opname(i)(1:3).eq.'muy')then
           ii=ii+1                                                      3d8s22
           mui(ii)=isop                                                 3d8s22
          end if                                                        3d8s22
          if(opname(i)(1:3).eq.'mu0')then
           i0=i0+1                                                      3d8s22
           mu0(i0)=isop                                                 3d8s22
          end if                                                        3d8s22
         end do
         ffr=-srh                                                       5d24s21
         ffi=-srh                                                       5d24s21
        else if(ipass.eq.2)then                                         5d24s21
         tkind='m1'                                                        3d8s22
         iq=1                                                           5d24s21
c        +1                    -1                0
c     (-lx-ily)/sqrt(2),(lx-ily)/sqrt(2)         lz=i*lz(i)
c     =(+ilx(i)-ly(i))/sqrt(2),(-ilx(i)-ly(i))/sqrt(2)
c
c     i.e. real +1 = -ly(i)/sqrt(2), imag +1=lx(i)/sqrt(2)
c          real -1 = same          , imag -1=-lx(i)/sqrt(2)
         ir=0                                                           3d8s22
         ii=0                                                           3d8s22
         i0=0                                                           3d8s22
         do isop=1,nopso                                                2d25s22
          i=iopso(2,isop)                                               2d25s22
          if(opname(i)(1:4).eq.'lx/i')then
           ii=ii+1                                                      2d25s22
           mui(ii)=isop                                                 2d25s22
          end if                                                        3d8s22
          if(opname(i)(1:4).eq.'ly/i')then                              3d8s22
           ir=ir+1                                                      2d25s22
           mur(ir)=isop                                                 2d25s22
          end if                                                        3d8s22
          if(opname(i)(1:4).eq.'lz/i')then
           i0=i0+1                                                      3d8s22
           mu0(i0)=isop                                                 3d8s22
          end if                                                        3d8s22
         end do
         ffr=-srh                                                       5d24s21
         ffi=srh                                                        5d24s21
        else                                                            5d24s21
         iq=2                                                           5d24s21
         tkind='e2'                                                        3d8s22
c        +1     -1    0   +/- 2
c     -zx-izy,zx-izy q0  Re Q2 +/- i Im Q2
c                          xx-yy      xy
c     i.e. real +1 = -zx, imag +1 = -zy
c          real -1 = +zx, imag -1 = same
         i0=0                                                           3d8s22
         ir=0                                                           3d8s22
         ii=0                                                           3d8s22
         irr=0                                                          3d8s22
         iii=0                                                          3d8s22
         do isop=1,nopso                                                2d25s22
          i=iopso(2,isop)                                               2d25s22
          if(opname(i)(1:5).eq.'Re Q1')then                             3d8s22
           ir=ir+1                                                      3d8s22
           mur(ir)=isop                                                 3d8s22
          end if                                                        3d8s22
          if(opname(i)(1:5).eq.'Im Q1')then                             3d8s22
           ii=ii+1                                                      3d8s22
           mui(ii)=isop                                                 3d8s22
          end if                                                        3d8s22
          if(opname(i)(1:2).eq.'Q0')then                                3d8s22
           i0=i0+1                                                      3d8s22
           mu0(i0)=isop                                                 3d8s22
          end if                                                        3d8s22
          if(opname(i)(1:5).eq.'Re Q2')then                             3d8s22
           irr=irr+1                                                    3d8s22
           mur2(irr)=isop                                               3d8s22
          end if                                                        3d8s22
          if(opname(i)(1:5).eq.'Im Q2')then                             3d8s22
           iii=iii+1                                                    3d8s22
           mui2(iii)=isop                                               3d8s22
          end if                                                        3d8s22
         end do
         ffr=1d0                                                        5d26s21
         ffi=1d0                                                        5d26s21
        end if                                                          5d24s21
        iq2=iq*2                                                        3d17s22
        phs=1d0                                                         3d17s22
        if(i2sb.ne.i2sk)then
         phs=-phs                                                       3d17s22
        end if
        if(ipass.eq.2)then                                              3d17s22
         phs=-phs                                                       3d17s22
        end if                                                          3d17s22
        if(llb.ne.llk)then                                              3d29s22
         phs=-phs                                                       3d29s22
        end if                                                          3d29s22
        mnz=0                                                           3d30s22
        do m2sb=mod(i2sb,2),i2sb,2                                      3d29s22
         if(m2sb.eq.0)then                                              3d29s22
          msbot=0                                                       3d29s22
         else                                                           3d29s22
          msbot=-i2sk                                                   3d29s22
         end if                                                         3d29s22
         do m2sk=msbot,i2sk,2                                           3d29s22
          if(iabs(m2sk-m2sb).le.2)then                                  3d29s22
           mq=(m2sk-m2sb)/2                                             3d29s22
           if(m2sk.eq.m2sb)then                                               3d8s22
            itbot=1                                                      3d8s22
            ittop=1                                                      3d8s22
           else                                                          3d8s22
            itbot=2                                                      3d8s22
            ittop=3                                                      3d8s22
           end if                                                        3d8s22
           do ipm=-iq,iq                                                3d29s22
            fmr=ffr                                                        3d8s22
            fmi=ffi                                                        3d8s22
            fmro=ffr                                                     3d31s22
            fmio=ffi                                                     3d31s22
            if(ipm.eq.0)then                                            3d31s22
            else if(ipm.eq.+1)then                                       3d31s22
             if(ipass.eq.2)then                                          3d31s22
              fmio=-fmio                                                 3d31s22
             else                                                        3d31s22
              fmro=-fmro                                                 3d31s22
             end if                                                        5d26s21
            else if(ipm.eq.-1)then                                         5d26s21
             if(ipass.eq.2)then                                            5d26s21
              fmi=-fmi                                                     3d8s22
             else                                                          5d26s21
              fmr=-fmr                                                     3d8s22
             end if                                                        5d26s21
            else if(ipm.eq.+2)then                                       3d31s22
             fmio=-fmio                                                  3d31s22
            else if(ipm.eq.-2)then                                         5d26s21
             fmi=-fmi                                                      3d8s22
            end if                                                         5d26s21
            phsu=phs                                                    3d29s22
            if(mod(ipm+mq,2).ne.0)phsu=-phsu                            3d29s22
            if(nsymb.eq.1)then                                          3d31s22
             do iz=0,nroot2*4-1                                         3d31s22
              bc(itmp+iz)=0d0                                            3d18s22
             end do                                                      3d17s22
            else                                                        3d31s22
             do iz=0,nroot2*2-1                                          3d29s22
              bc(itmp+iz)=0d0                                            3d18s22
             end do                                                      3d17s22
            end if                                                      3d31s22
            iok=0                                                       3d30s22
            if(ipm.eq.0)then                                             5d26s21
             do it=itbot,ittop                                          2d25s22
              itp=it+1-l2sub                                            3d2s22
              i=iopso(2,mu0(it))                                        2d25s22
              j=iopso(3,mu0(it))                                        2d25s22
              isoptm(1)=iopso(1,mu0(it))                                2d25s22
              if(ipass.eq.2)then                                        3d8s22
               isoptm(2)=mod(1+isopt(2,j),2)                            3d1s22
              else                                                      3d8s22
               isoptm(2)=isopt(2,j)                                     3d1s22
              end if                                                    3d8s22
              isoptm(3)=isopt(3,j)                                      2d25s22
              isoptm(4)=isopt(4,j)                                      2d25s22
              call genhsoa2(iwfbra,iwfket,nspc,mlb2,mlk2,m2sb,m2sk,
     $          bc(iaout),nec,multh,irefo,ih0a,i4so2,irel,ism,norb,mdon,5d19s21
     $          iri,itypex,nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,  9d20s21
     $     irw2,ixmt(1,mu0(it)),nh0,isoptm,nsopt,i4or,ionexr,jmatsr,
     $       kmatsr,kmatsrb,i3xr,iifmx,ntype,npadddi,nbasp,nbaspc,natom,   12d20s20
     $           ngaus,ibdat,iapair,ibstor,isstor,isym,ascale,idorel,   12d23s21
     $           iorb,itp,iopso(1,mu0(it)),isoptm(2),bc,ibc,n4vso,1,0)  3d27s24
              szi=0d0                                                   3d9s22
              do iz=0,nroot2-1                                          3d9s22
               szi=szi+bc(iaout+iz+nroot2)**2                           3d9s22
              end do                                                    3d9s22
              szi=sqrt(szi/dfloat(nroot2))                              3d9s22
              iause=iaout                                               3d30s22
              itmpu=itmp                                                3d31s22
              if(szi.gt.1d-10)then                                      3d30s22
               iok=iok+1                                                3d30s22
               iause=iaout+nroot2                                       3d30s22
               itmpu=itmpi                                              3d31s22
              end if                                                    3d30s22
              if(ipass.eq.2.and.isopt(2,j).ne.0)then                    3d8s22
               do iz=0,nroot2-1                                         3d1s22
                bc(itmpu+iz)=bc(itmpu+iz)-bc(iause+iz)                    3d30s22
               end do                                                    2d25s22
              else                                                      3d8s22
               do iz=0,nroot2-1                                         3d1s22
                bc(itmpu+iz)=bc(itmpu+iz)+bc(iause+iz)                    3d30s22
               end do                                                    2d25s22
              end if                                                    3d8s22
             end do                                                     3d8s22
            else if(iabs(ipm).eq.1)then                                 3d9s22
             do it=itbot,ittop                                          2d25s22
              itp=it+1-l2sub                                            3d2s22
              i=iopso(2,mur(it))                                        2d25s22
              j=iopso(3,mur(it))                                        2d25s22
              isoptm(1)=iopso(1,mur(it))                                2d25s22
              isoptm(2)=isopt(2,j)                                      3d1s22
              isoptm(3)=isopt(3,j)                                      2d25s22
              isoptm(4)=isopt(4,j)                                      2d25s22
              call genhsoa2(iwfbra,iwfket,nspc,mlb2,mlk2,m2sb,m2sk,
     $          bc(iaout),nec,multh,irefo,ih0a,i4so2,irel,ism,norb,mdon,5d19s21
     $          iri,itypex,nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,  9d20s21
     $     irw2,ixmt(1,mur(it)),nh0,isoptm,nsopt,i4or,ionexr,jmatsr,
     $       kmatsr,kmatsrb,i3xr,iifmx,ntype,npadddi,nbasp,nbaspc,natom,   12d20s20
     $           ngaus,ibdat,iapair,ibstor,isstor,isym,ascale,idorel,   12d23s21
     $           iorb,itp,iopso(1,mur(it)),isoptm(2),bc,ibc,n4vso,1,0)  3d27s24
              szi=0d0                                                   3d9s22
              do iz=0,nroot2-1                                          3d9s22
               szi=szi+bc(iaout+iz+nroot2)**2                           3d9s22
              end do                                                    3d9s22
              szi=sqrt(szi/dfloat(nroot2))                              3d9s22
              iause=iaout                                               3d30s22
              itmpu=itmp                                                3d31s22
              if(szi.gt.1d-10)then                                      3d30s22
               iok=iok+1                                                3d30s22
               iause=iaout+nroot2                                       3d30s22
               itmpu=itmpi                                              3d31s22
              end if                                                    3d30s22
              do iz=0,nroot2-1                                          3d1s22
               bc(itmpu+iz)=bc(itmpu+iz)+bc(iause+iz)*fmr               3d31s22
              end do                                                    2d25s22
              i=iopso(2,mui(it))                                        2d25s22
              j=iopso(3,mui(it))                                        2d25s22
              isoptm(1)=iopso(1,mui(it))                                2d25s22
              isoptm(2)=mod(1+isopt(2,j),2)                             3d1s22
              isoptm(3)=isopt(3,j)                                      2d25s22
              isoptm(4)=isopt(4,j)                                      2d25s22
              call genhsoa2(iwfbra,iwfket,nspc,mlb2,mlk2,m2sb,m2sk,
     $          bc(iaout),nec,multh,irefo,ih0a,i4so2,irel,ism,norb,mdon,5d19s21
     $          iri,itypex,nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,  9d20s21
     $     irw2,ixmt(1,mui(it)),nh0,isoptm,nsopt,i4or,ionexr,jmatsr,
     $       kmatsr,kmatsrb,i3xr,iifmx,ntype,npadddi,nbasp,nbaspc,natom,   12d20s20
     $           ngaus,ibdat,iapair,ibstor,isstor,isym,ascale,idorel,   12d23s21
     $           iorb,itp,iopso(1,mui(it)),isoptm(2),bc,ibc,n4vso,1,0)  3d27s24
              ffu=fmi                                                   3d8s22
              szi=0d0                                                   3d9s22
              do iz=0,nroot2-1                                          3d9s22
               szi=szi+bc(iaout+iz+nroot2)**2                           3d9s22
              end do                                                    3d9s22
              szi=sqrt(szi/dfloat(nroot2))                              3d9s22
              iause=iaout                                               3d30s22
              itmpu=itmp                                                3d31s22
              if(szi.gt.1d-10)then                                      3d30s22
               iok=iok+1                                                3d30s22
               iause=iaout+nroot2                                       3d30s22
               itmpu=itmpi                                              3d31s22
              end if                                                    3d30s22
              if(isopt(2,j).ne.0)ffu=-ffu                               3d1s22
              do iz=0,nroot2-1                                          2d25s22
               bc(itmpu+iz)=bc(itmpu+iz)+bc(iause+iz)*ffu               3d31s22
              end do                                                    2d25s22
             end do                                                     2d25s22
            else                                                        3d8s22
             do it=itbot,ittop                                          2d25s22
              itp=it+1-l2sub                                            3d2s22
              i=iopso(2,mur2(it))                                        2d25s22
              j=iopso(3,mur2(it))                                        2d25s22
              isoptm(1)=iopso(1,mur2(it))                                2d25s22
              isoptm(2)=isopt(2,j)                                      3d1s22
              isoptm(3)=isopt(3,j)                                      2d25s22
              isoptm(4)=isopt(4,j)                                      2d25s22
              call genhsoa2(iwfbra,iwfket,nspc,mlb2,mlk2,m2sb,m2sk,
     $          bc(iaout),nec,multh,irefo,ih0a,i4so2,irel,ism,norb,mdon,5d19s21
     $          iri,itypex,nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,  9d20s21
     $     irw2,ixmt(1,mur2(it)),nh0,isoptm,nsopt,i4or,ionexr,jmatsr,
     $       kmatsr,kmatsrb,i3xr,iifmx,ntype,npadddi,nbasp,nbaspc,natom,   12d20s20
     $           ngaus,ibdat,iapair,ibstor,isstor,isym,ascale,idorel,   12d23s21
     $           iorb,itp,iopso(1,mur2(it)),isoptm(2),bc,ibc,n4vso,1,0) 3d27s24
              szi=0d0                                                   3d9s22
              do iz=0,nroot2-1                                          3d9s22
               szi=szi+bc(iaout+iz+nroot2)**2                           3d9s22
              end do                                                    3d9s22
              szi=sqrt(szi/dfloat(nroot2))                              3d9s22
              iause=iaout                                               3d30s22
              itmpu=itmp                                                3d31s22
              if(szi.gt.1d-10)then                                      3d30s22
               iok=iok+1                                                3d30s22
               iause=iaout+nroot2                                       3d30s22
               itmpu=itmpi                                              3d31s22
              end if                                                    3d30s22
              do iz=0,nroot2-1                                          3d1s22
               bc(itmpu+iz)=bc(itmpu+iz)+bc(iause+iz)*fmr               3d31s22
              end do                                                    2d25s22
              i=iopso(2,mui2(it))                                       3d8s22
              j=iopso(3,mui2(it))                                       3d8s22
              isoptm(1)=iopso(1,mui2(it))                               3d8s22
              isoptm(2)=mod(1+isopt(2,j),2)                             3d1s22
              isoptm(3)=isopt(3,j)                                      2d25s22
              isoptm(4)=isopt(4,j)                                      2d25s22
              call genhsoa2(iwfbra,iwfket,nspc,mlb2,mlk2,m2sb,m2sk,
     $          bc(iaout),nec,multh,irefo,ih0a,i4so2,irel,ism,norb,mdon,5d19s21
     $          iri,itypex,nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,  9d20s21
     $     irw2,ixmt(1,mui2(it)),nh0,isoptm,nsopt,i4or,ionexr,jmatsr,   3d8s22
     $       kmatsr,kmatsrb,i3xr,iifmx,ntype,npadddi,nbasp,nbaspc,natom,   12d20s20
     $           ngaus,ibdat,iapair,ibstor,isstor,isym,ascale,idorel,   12d23s21
     $           iorb,itp,iopso(1,mui2(it)),isoptm(2),bc,ibc,n4vso,1,0) 3d27s24
              szi=0d0                                                   3d9s22
              do iz=0,nroot2-1                                          3d9s22
               szi=szi+bc(iaout+iz+nroot2)**2                           3d9s22
              end do                                                    3d9s22
              szi=sqrt(szi/dfloat(nroot2))                              3d9s22
              iause=iaout                                               3d30s22
              itmpu=itmp                                                3d31s22
              if(szi.gt.1d-10)then                                      3d30s22
               iok=iok+1                                                3d30s22
               iause=iaout+nroot2                                       3d30s22
               itmpu=itmpi                                              3d31s22
              end if                                                    3d30s22
              ffu=fmi                                                   3d8s22
              if(isopt(2,j).ne.0)ffu=-ffu                               3d1s22
              do iz=0,nroot2-1                                          2d25s22
               bc(itmpu+iz)=bc(itmpu+iz)+bc(iause+iz)*ffu               3d31s22
              end do                                                    2d25s22
             end do                                                     2d25s22
            end if                                                      3d8s22
            if(iok.eq.0)then                                            3d30s22
             iok=1                                                      3d30s22
            else                                                        3d30s22
             iok=2                                                      3d30s22
            end if                                                      3d30s22
            rms=0d0                                                     3d31s22
            do ik=1,nrootk                                              3d28s22
             jk=ik-1                                                    3d28s22
             do ib=1,nrootb                                             3d28s22
              jb=ib-1                                                   3d28s22
              iad=itmp+jb+nrootb*jk                                     3d28s22
              iad2=itmp2+jb+nrootb*jk                                     3d28s22
              if(lprint)then                                            3d31s22
               if(nsymb.gt.1)then                                       3d31s22
                if(abs(bc(iad)).gt.1d-12)then                           7d30s22
                 if(mnz.eq.0)write(6,399)                                 3d30s22
  399            format(8x,'2Sigma',4x,'Mtm',x,'Ms',7x,'2Sigma''',4x,   12d23s22
     $              'normal',8x,'reversed phase',5x,'times 1/4c^2')     7d30s22
                 mnz=mnz+1                                                3d30s22
                 write(6,400)ib,iwfbra(1,1),labb,m2sb,tkind,ipm,mq,       3d28s22
     $             ik,iwfket(1,1),labk,m2sk,bc(iad),crori(iok),         3d30s22
     $                phsu,bc(iad)*ascale                               7d30s22
                end if                                                  3d31s22
               else                                                     3d31s22
                iadi=itmpi+jb+nrootb*jk                                 3d31s22
                iad2i=itmp2i+jb+nrootb*jk                               3d31s22
                sn=sqrt(bc(iad)**2+bc(iadi)**2)                         3d31s22
                if(sn.gt.1d-12)then                                     7d30s22
                 if(mnz.eq.0)write(6,399)                                 3d30s22
                 mnz=mnz+1                                              3d31s22
                 if(abs(bc(iad)).gt.1d-12)then                          7d30s22
                  write(6,400)ib,iwfbra(1,1),labb,m2sb,tkind,ipm,mq,       3d28s22
     $             ik,iwfket(1,1),labk,m2sk,bc(iad),crori(1),
     $                phsu,bc(iad)*ascale                               7d30s22
                  if(abs(bc(iadi)).gt.1d-12)then                        7d30s22
                   write(6,404)bc(iadi),-phsu,bc(iadi)*ascale                                 3d31s2
  404              format(37x,'>=',es15.7,x,'i',f4.0,12x,               7d30s22
     $                  es15.7)                                         3d31s22
                  end if                                                3d31s22
                 else                                                   3d31s22
                  write(6,400)ib,iwfbra(1,1),labb,m2sb,tkind,ipm,mq,       3d28s22
     $             ik,iwfket(1,1),labk,m2sk,bc(iadi),crori(2),phsu,     7d30s22
     $                 bc(iadi)*ascale                                  7d30s22
                 end if                                                 3d31s22
                end if                                                  3d31s22a
               end if                                                   3d31s22
              end if                                                    3d29s22
             end do                                                     3d29s22
            end do                                                      3d29s22
           end do                                                       3d29s22
          end if                                                        3d29s22
         end do                                                         3d29s22
        end do                                                          3d29s22
  400          format('<',i2,x,i1,a6,i3,'|',a2,2i3,'|',i2,x,i1,a6,i3,   3d28s22
     $              '>=',es15.7,x,a1,f4.0,12x,es15.7)                   7d30s22
       end do                                                           3d28s22
      end if                                                            2d25s22
      ibcoff=ibcoffo                                                    3d2s22
      return                                                            5d18s21
      end                                                               5d18s21
