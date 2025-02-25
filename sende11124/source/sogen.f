c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine sogen(iwfbra,iwfket,nspc,nec,multh,irefo,ih0a,i4so2,   3d17s22
     $     irel,ism,norb,mdon,nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,    9d20s21
     $     irw1,irw2,ih0n,nh0,isopt,nsopt,i4or,ionexr,jmatsr,kmatsr,    10d13s21
     $     kmatsrb,i3xr,iifmx,ntype,npadddi,nbasp,nbaspc,natom,ngaus,   12d20s20
     $     ibdat,iapair,ibstor,isstor,isym,ascale,idorel,iorb,lprint,   11d9s22
     $     bc,ibc,shift,n4vso)                                          2d8s23
      implicit real*8 (a-h,o-z)                                         5d18s21
      logical ldiag,lprint,logb,logk,nlogk,lok                          3d28s22
      character*6 labb,labk                                             10d27s20
      character*1 crori(2)                                              3d30s22
      integer*1 ipack1(4)                                               5d18s21
      equivalence (ipack1,npack4)                                       5d18s21
      dimension iwfbra(*),iwfket(*),multh(8,8),irefo(*),                3d17s22
     $     ih0a(8,2,8),i4so2(8,8,8,2),irel(*),ism(*),nvirt(*),nbasp(*), 12d20s20
     $     nbaspc(*),isopt(4,4)                                         3d28s22
      include "common.store"                                            5d18s21
      common/singcm/iuse,nff
      data crori/' ','i'/                                               3d30s22
      data icall/0/
      npack4=iwfbra(6)                                                  3d28s22
      iordb=ipack1(3)                                                   3d28s22
      npack4=iwfket(6)                                                  3d28s22
      iordk=ipack1(3)                                                   3d28s22
      do i=1,6                                                          10d27s20
       if(iwfbra(13+i).ne.0)then                                         10d27s20
        labb(i:i)=char(iwfbra(13+i))                                     10d27s20
       else                                                             10d27s20
        labb(i:i)=' '                                                   10d27s20
       end if                                                           10d27s20
       if(iwfket(13+i).ne.0)then                                        3d18s22
        labk(i:i)=char(iwfket(13+i))                                    3d18s22
       else                                                             10d27s20
        labk(i:i)=' '                                                   10d27s20
       end if                                                           10d27s20
      end do                                                            10d27s20
      ibcoffo=ibcoff                                                    5d18s21
      ldiag=loc(iwfbra).eq.loc(iwfket)
      i2sb=iwfbra(1)-1                                                  3d17s22
      i2sk=iwfket(1)-1                                                  3d17s22
      if(lprint)then                                                    3d8s22
       write(6,*)('in sogen for '),iwfbra(1),labb,('and '),             3d17s22
     $     iwfket(1),labk                                               3d17s22
      end if                                                            3d8s22
      inz=0                                                             3d30s22
      nrootb=iwfbra(3)                                                  3d17s22
      nrootk=iwfket(3)                                                  3d17s22
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
      nf=nq                                                             3d24s22
      ndat=iwfbra(1)*iwfket(1)                                          3d17s22
      idat=ibcoff                                                       3d17s22
      ifcn=idat+ndat*nroot2                                             3d17s22
      ibcoff=ifcn+ndat*nq                                               3d17s22
      iaout=ibcoff                                                      5d18s21
      iaouti=iaout+nroot2                                               3d28s22
      ibcoff=iaouti+nroot2                                              3d28s22
      call enough('sogen.  1',bc,ibc)
      isymprod=multh(iwfbra(2),iwfket(2))                               3d28s22
      do msb2=mod(i2sb,2),i2sb,2                                        3d29s22
       if(msb2.eq.0)then                                                3d29s22
        mbot=0                                                          3d29s22
       else                                                             3d29s22
        mbot=-i2sk                                                      3d29s22
       end if                                                           3d29s22
       do msk2=mbot,i2sk,2                                              3d29s22
        if(iabs(msb2-msk2).le.iqtop)then                                3d29s22
         mq2=msk2-msb2                                                  3d29s22
         if(mq2.eq.0)then                                                 3d28s22
          isgoal=0                                                        3d28s22
         else                                                             3d28s22
          isgoal=1                                                        3d28s22
         end if                                                           3d28s22
         itu=-1                                                           3d28s22
c
c     if no symmetry, we need to go for both real and imag parts.
c     fake psisopsi into doing this by temporarily resetting isopt(1
c
         if(nsymb.eq.1)then                                             3d30s22
          if(isgoal.eq.0)then                                            3d30s22
           if(nsopt.eq.4)then                                           3d30s22
            npass=2                                                     3d30s22
            ireset1=1                                                   3d30s22
            ireset2=2                                                   3d30s22
           else                                                         3d30s22
            npass=1                                                     3d30s22
           end if                                                       3d30s22
          else                                                           3d30s22
           npass=2                                                      3d30s22
           if(nsopt.eq.4)then                                           3d30s22
            ireset1=3                                                   3d30s22
            ireset2=4                                                   3d30s22
           else                                                         3d30s22
            ireset1=2                                                   3d30s22
            ireset2=3                                                   3d30s22
           end if                                                       3d30s22
          end if                                                        3d30s22
         else                                                           3d30s22
          npass=1                                                       3d30s22
         end if                                                         3d30s22
         do ipass=1,npass                                               3d30s22
          if(npass.eq.2)then                                            3d30s22
           if(ipass.eq.1)then                                           3d30s22
            isopt(1,ireset2)=2                                          3d30s22
           else                                                         3d30s22
            isopt(1,ireset1)=2                                          3d30s22
           end if                                                       3d30s22
          end if                                                        3d30s22
          do it=1,nsopt                                                    3d28s22
           if(isopt(1,it).eq.isymprod.and.isopt(3,it).eq.isgoal)then       3d28s22
            itu=it                                                         3d28s22
           end if                                                          3d28s22
          end do                                                           3d28s22
          if(itu.gt.0)then                                                 3d28s22
           lri=isopt(2,itu)                                                 3d28s22
c
c     phase for -sig
c
c     (-1)^2J
           if(mod(i2sb,2).eq.0)then                                       3d29s22
            phs=1d0                                                       3d29s22
           else                                                           3d29s22
            phs=-1d0                                                      3d29s22
           end if                                                         3d29s22
c     (-1)^(-S-S')
           iss=i2sb+i2sk                                                 3d29s22
           iss=iss/2                                                     3d29s22
           if(mod(-iss,2).ne.0)phs=-phs                                  3d29s22
c     (-1)^(Sig'-Sig)
           iss=-msb2+msk2                                                 3d29s22
           iss=iss/2                                                     3d29s22
           if(mod(-iss,2).ne.0)phs=-phs                                  3d29s22
c     (-1)^(srho+srho')
           if(iordb.ne.iordk)phs=-phs                                    3d29s22
           write(6,*)('gensoa2 '),msb2,msk2,shift,lri,iordb,iordk
           irori=1                                                      3d19s24
           if(nsymb.eq.1.and.lri.eq.1)irori=2                           3d19s24
           write(6,*)('irori = '),irori                                 3d29s24
           call genhsoa2(iwfbra,iwfket,nspc,mlb2,mlk2,msb2,msk2,           3d17s22
     $          bc(iaout),nec,multh,irefo,ih0a,i4so2,irel,ism,norb,mdon,5d19s21
     $          iri,idum,nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,    9d20s21
     $           irw2,ih0n,nh0,isopt,nsopt,i4or,ionexr,jmatsr,kmatsr,   3d9s22
     $           kmatsrb,i3xr,iifmx,ntype,npadddi,nbasp,nbaspc,natom,   3d9s22
     $           ngaus,ibdat,iapair,ibstor,isstor,isym,ascale,idorel,   3d9s22
     $           iorb,0,0,0,bc,ibc,n4vso,irori,0)                       3d27s24
           lok=(lri.eq.0.and.iordb.eq.iordk).or.
     $       (lri.ne.0.and.iordb.ne.iordk)
           iok=2                                                         3d30s22
           if(lok)iok=1                                                  3d30s22
           if(.not.lok)phs=-phs                                          3d30s22
           if(nsymb.eq.1.and.(msb2.ne.-msb2.or.msk2.ne.-msk2))then      3d30s22
            do ik=1,nrootk                                               3d28s22
             jk=ik-1                                                     3d28s22
             do ib=1,nrootb                                              3d28s22
              jb=ib-1                                                    3d28s22
              iad=iaout+jb+nrootb*jk                                     3d28s22
              if(ib.eq.ik.and.loc(iwfbra).eq.loc(iwfket).and.           3d18s24
     $             msb2.eq.msk2.and.iok.eq.1)then                       3d19s24
               bc(iad)=bc(iad)+shift                                    11d21s22
              end if                                                    11d21s22
              if(abs(bc(iad)).gt.1d-12.and.lprint)then                  3d18s24
               if(inz.eq.0)write(6,399)                                 3d30s22
               inz=inz+1                                                3d30s22
               write(6,401)ib,                                          3d30s22
     $           iwfbra(1),labb,msb2,ik,iwfket(1),labk,msk2,bc(iad),    3d29s22
     $            crori(iok),phs
              end if                                                    3d30s22
             end do
            end do
           else                                                          3d29s22
            do ik=1,nrootk                                               3d28s22
             jk=ik-1                                                     3d28s22
             do ib=1,nrootb                                              3d28s22
              jb=ib-1                                                    3d28s22
              iad=iaout+jb+nrootb*jk                                     3d28s22
              write(6,*)('aoutb '),bc(iad),iad-iaout
              if(ib.eq.ik.and.loc(iwfbra).eq.loc(iwfket)                3d18s24
     $             .and.msb2.eq.msk2)then                               3d18s24
               bc(iad)=bc(iad)+shift                                    11d21s22
              end if                                                    11d21s22
              if(abs(bc(iad)).gt.1d-12.and.lprint)then                  3d30s22
               if(inz.eq.0)write(6,399)                                 3d30s22
  399          format(10x,'2Sig',14x,'2Sig''',6x,'normal',7x,           7d30s22
     $              'reversed Sigmas phase')                            3d19s24
               inz=inz+1                                                3d30s22
               write(6,401)ib,iwfbra(1),                                3d30s22
     $            labb,msb2,ik,iwfket(1),labk,msk2,bc(iad),crori(iok),  3d30s22
     $            phs                                                   3d30s22
  401         format('<',i2,x,i1,a6,i3,'|Hso|',i2,x,i1,a6,i3,'>=',       3d28s22
     $            es22.14,x,a1,f5.0,es22.14,x,a1)                         3d30s22
              end if                                                    3d30s22
             end do                                                      3d28s22
            end do                                                       3d28s22
           end if                                                        3d29s22
          end if                                                        3d30s22
          if(npass.eq.2)then                                            3d30s22
           if(ipass.eq.1)then                                           3d30s22
            isopt(1,ireset2)=1                                          3d30s22
           else                                                         3d30s22
            isopt(1,ireset1)=1                                          3d30s22
           end if                                                       3d30s22
          end if                                                        3d30s22
         end do                                                         3d30s22
        end if                                                          3d29s22
       end do                                                           3d29s22
      end do                                                            3d29s22
      ibcoff=ibcoffo                                                    5d18s21
      return                                                            5d18s21
      end                                                               5d18s21
