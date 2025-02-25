c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine parah0d2(natom,ngaus,ibdat,nbasis,h0,ovr,ovrdk,        6d20s16
     $        isym,iapair,ibstor,isstor,iso,nbb,idwsdeb,idorel,ascale,  4d4s22
     $                multh,ixyz,iatom,ipropsym,idersign,nbasisp,bc,ibc)11d10s22
      implicit real*8 (a-h,o-z)
      external second
      integer*8 nn8,iab,iak                                             4d22s22
      include "common.hf"
      include "common.store"
      include "common.spher"
      integer*8 ibstor,isstor                                           5d6s10
c
c     for d^2 S, need (d^2n|m)+(dn|dm)+(n|d^2m)
c     [if rel, need (d^2n|T|m)+(dn|T|dm)+(n|T|d^2m) as well]
c     for d^2 (n|h|m), need (d^2n|h|m), (dn|dh|m),(dn|h|dm),
c     (n|d^2h|m),(n|dh|dm),(n|h|d^2m)
c     and  (n|d^2m).
c
      include "common.basis"
      logical lpr
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension h0(*),ovr(*),isym(3,*),iapair(3,*),ibstor(*),isstor(*), 2d17s23
     $     iso(*),cartb(3),cartk(3),idxyz(3),multh(8,8),ovrdk(*),       2d17s23
     $     nbasisp(*),igoalx(6),igoaly(2,6)                                                   4d4s22
      data igoalx/1913950,1913954,4*0/
      data igoaly/1913950,1913951,1913954,1913955,8*0/
      ibcoffo=ibcoff                                                    2d19s10
      call second(time1)
      iseek=32*63+1
c
c     build shell pair order
c
      srh=sqrt(0.5d0)                                                   5d3s10
      ascale2=ascale*2d0                                                8d20s15
      nn=(ngaus*(ngaus+1))/2                                            2d19s10
      nn2=nn*2                                                          2d19s10
      if(idorel.eq.0)then                                               8d20s15
c
c     overlap and kinetic energy together
c     and nuclear attraction for each atom
c
       ndo=1+natom                                                       2d19s10
      else                                                              8d20s15
c
c     overlap and kinetic energy together
c     and nuclear attraction
c     and pxVpx, pyVpy, pzVpz for each atom
c
       nbasall2=nbasall*2                                               8d20s15
       ndo=1+natom*4                                                    8d20s15
      end if                                                            8d20s15
      ipair=ibcoff                                                      2d19s10
      ibcoff=ipair+nn*ndo*3                                             2d19s10
      i12=ibcoff                                                        2d19s10
      ibcoff=i12+nn*4                                                   2d19s10
      call enough('parah0d2.  1',bc,ibc)
      j12=i12                                                           2d19s10
      jbdat=ibdat-1                                                     2d19s10
      ngaus7=ngaus*7                                                    5d3s10
      do i1=1,ngaus                                                     2d19s10
       do i2=1,i1                                                       2d19s10
        ibc(j12)=-((ibc(jbdat+i1)+ibc(jbdat+i2)+2)**2)                  2d19s10
        if(iapair(1,ibc(jbdat+i1+ngaus7)).gt.0)ibc(j12)=ibc(j12)*2      5d3s10
        if(iapair(1,ibc(jbdat+i2+ngaus7)).gt.0)ibc(j12)=ibc(j12)*2      5d3s10
        ibc(j12+nn)=i1                                                  2d19s10
        ibc(j12+nn2)=i2                                                 2d19s10
        j12=j12+1                                                       2d19s10
       end do                                                           2d19s10
      end do                                                            2d19s10
      is12=i12+nn*3                                                     2d19s10
      nn8=nn                                                            2d19s10
      call idsortdws(ibc(i12),ibc(is12),nn8)                            1d18s23
      ii1=i12+nn                                                        2d19s10
      ii2=i12+nn2                                                       2d19s10
      jpair=ipair                                                       2d19s10
      do i=1,nn                                                         2d19s10
       j=ibc(is12+i-1)-1                                                2d19s10
    1  format(i5,i8,2i5)                                                2d19s10
       do k=1,ndo                                                       2d19s10
        ibc(jpair)=ibc(ii1+j)                                           2d19s10
        jpair=jpair+1                                                   2d19s10
        ibc(jpair)=ibc(ii2+j)                                           2d19s10
        jpair=jpair+1                                                   2d19s10
        ibc(jpair)=k-1                                                  2d22s10
        jpair=jpair+1                                                   2d19s10
       end do                                                           2d19s10
      end do                                                            2d19s10
      nneed=nn*ndo                                                      2d19s10
       idxyz(1)=0
       idxyz(2)=0
       idxyz(3)=0
       factna=1d0
       idxyz(ixyz)=2                                                    6d20s16
c
c     If we are doing a relativistic calculation, then parahf multiplied
c     nbb by 4, so we have room to store 2 component parts of ham and
c     overlap.
c
      do i=1,nbb                                                        4d18s22
       h0(i)=0d0
       ovr(i)=0d0                                                       2d24s10
       ovrdk(i)=0d0                                                     4d4s16
      end do                                                            2d19s10
      call second(time2)
      telap=time2-time1
      do i=1+mynowprog,nneed,mynprocg                                   2d19s10
       jpair=ipair+3*(i-1)                                              2d19s10
    2  format('I am going to do ',i5,3i3,f18.5)                               2d19s10
       jbra=jbdat+ibc(jpair)                                            2d19s10
       jbra2=jbra+ngaus                                                 2d19s10
       jbra3=jbra2+ngaus                                                2d19s10
       jbra4=jbra3+ngaus                                                2d19s10
       jbra5=jbra4+ngaus                                                2d19s10
       jbra6=jbra5+ngaus                                                2d19s10
       jbra7=jbra6+ngaus                                                2d19s10
       jbra8=jbra7+ngaus                                                5d4s10
       jket=jbdat+ibc(jpair+1)                                          2d19s10
       jket2=jket+ngaus                                                 2d19s10
       jket3=jket2+ngaus                                                2d19s10
       jket4=jket3+ngaus                                                2d19s10
       jket5=jket4+ngaus                                                2d19s10
       jket6=jket5+ngaus                                                2d19s10
       jket7=jket6+ngaus                                                2d19s10
       jket8=jket7+ngaus                                                5d4s10
        nbra=2*ibc(jbra)+1                                              5d3s10
        nbraa=nbra                                                      5d3s10
        if(iapair(1,ibc(jbra8)).gt.0)nbraa=nbra*2                       5d3s10
        nket=2*ibc(jket)+1                                              5d3s10
        nketa=nket                                                      5d3s10
        if(iapair(1,ibc(jket8)).gt.0)nketa=nket*2                       5d3s10
        itmp1=ibcoff                                                    5d3s10
        ibcoff=itmp1+nbraa*nketa                                        4d26s16
        if(ibc(jpair+2).eq.0)ibcoff=ibcoff+nbraa*nketa*3                4d26s16
        nneedz=nbra*nket                                                4d5s16
        nneedzm=nneedz-1                                                2d3s16
        nneedz2=nneedz*2                                                2d3s16
        nneedz2m=nneedz2-1                                              2d3s16
        call enough('parah0d2.  2',bc,ibc)
        idrb=1                                                          2d4s16
        idrk=1                                                          2d4s16
        call oned2ints(ibc(jpair+2),ibc(jket8),ibc(jbra8),
     $     ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),bc(jbra7),
     $     ibc(jket),bc(jket2),bc(jket3),bc(jket5),bc(jket6),bc(jket7),
     $     idxyz,ib1,ib2,ib3,ib4,ib5,ib6,nneedzm,nneedz,nneedz2,        6d20s16
     $     nneedz2m,npass,idrb,idrk,                                    2d27s23
     $     iapair(1,iatom),natom,idorel,xcart,atnum,iatom,              6d20s16
     $     ixyz,mynowprog,i,idersign,bc,ibc)                            11d10s22
        ibu=ib1                                                         5d7s10
        itmpu=itmp1                                                     5d7s10
        do ipass=1,npass                                                5d7s10
         do ik=1,nket                                                   5d7s10
          do ib=1,nbra                                                  5d7s10
           iad1=ik-1+nket*(ib-1)                                        5d7s10
           iad2=ib-1+nbraa*(ik-1)                                       5d7s10
           bc(itmpu+iad2)=bc(ibu+iad1)                                  8d20s15
          end do                                                        5d7s10
         end do                                                         5d7s10
         ibu=ibu+nbra*nket                                              4d5s16
         itmpu=itmpu+nbraa*nketa                                        4d5s16
        end do                                                          5d7s10
        if(nketa.ne.nket)then                                           5d3s10
         cartk(1)=bc(jket5)*dfloat(isym(1,iapair(2,ibc(jket8))))        5d5s10
         cartk(2)=bc(jket6)*dfloat(isym(2,iapair(2,ibc(jket8))))        5d5s10
         cartk(3)=bc(jket7)*dfloat(isym(3,iapair(2,ibc(jket8))))        5d5s10
         idrb=1                                                         2d4s16
         idrk=2                                                         2d27s23
         iak=iapair(1,ibc(jket8))                                       4d22s22
         call oned2ints(ibc(jpair+2),iak,ibc(jbra8),                    4d22s22
     $     ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),bc(jbra7),
     $     ibc(jket),bc(jket2),bc(jket3),cartk,cartk(2),cartk(3),       2d4s16
     $     idxyz,ib1,ib2,ib3,ib4,ib5,ib6,nneedzm,nneedz,nneedz2,        6d20s16
     $        nneedz2m,npass,idrb,idrk,                                 2d27s23
     $     iapair(1,iatom),natom,idorel,xcart,atnum,iatom,              6d20s16
     $        ixyz,mynowprog,i,idersign,bc,ibc)                         11d10s22
         ibu=ib1                                                        5d7s10
         itmpu=itmp1                                                    5d7s10
         do ipass=1,npass                                               5d7s10
          do ik=1,nket                                                    5d3s10
           do ib=1,nbra                                                   5d3s10
            iad1=ik-1+nket*(ib-1)                                         5d3s10
            iad2=ib-1+nbraa*(ik-1+nket)                                  5d3s10
            bc(itmpu+iad2)=bc(ibu+iad1)                                 8d20s15
           end do                                                         5d3s10
          end do                                                          5d3s10
          ibu=ibu+nbra*nket                                             4d5s16
          itmpu=itmpu+nbraa*nketa                                       4d5s16
         end do                                                         5d7s10
        end if                                                          5d3s10
        if(nbraa.ne.nbra)then                                           5d3s10
         cartb(1)=bc(jbra5)*dfloat(isym(1,iapair(2,ibc(jbra8))))        5d5s10
         cartb(2)=bc(jbra6)*dfloat(isym(2,iapair(2,ibc(jbra8))))        5d5s10
         cartb(3)=bc(jbra7)*dfloat(isym(3,iapair(2,ibc(jbra8))))        5d5s10
         idrk=1                                                         2d27s23
         idrb=2                                                         2d27s23
         iab=iapair(1,ibc(jbra8))                                       4d22s22
         call oned2ints(ibc(jpair+2),ibc(jket8),iab,                    4d22s22
     $     ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),cartb(3),       2d4s16
     $     ibc(jket),bc(jket2),bc(jket3),bc(jket5),bc(jket6),bc(jket7),
     $     idxyz,ib1,ib2,ib3,ib4,ib5,ib6,nneedzm,nneedz,nneedz2,        6d20s16
     $        nneedz2m,npass,idrb,idrk,                                 2d27s23
     $     iapair(1,iatom),natom,idorel,xcart,atnum,iatom,              6d20s16
     $        ixyz,mynowprog,i,idersign,bc,ibc)                         11d10s22
         ibu=ib1                                                         5d7s10
         itmpu=itmp1                                                     5d7s10
         do ipass=1,npass                                                5d7s10
          do ik=1,nket                                                    5d3s10
           do ib=1,nbra                                                   5d3s10
            iad1=ik-1+nket*(ib-1)                                         5d3s10
            iad2=ib-1+nbra+nbraa*(ik-1)                                   5d3s10
            bc(itmpu+iad2)=bc(ibu+iad1)                                  8d20s15
           end do                                                        5d7s10
          end do                                                         5d3s10
          ibu=ibu+nbra*nket                                             4d5s16
          itmpu=itmpu+nbraa*nketa                                       4d5s16
         end do                                                          5d3s10
         if(nketa.ne.nket)then                                           5d3s10
          idrb=2                                                        2d27s23
          idrk=2                                                        2d27s23
          iab=iapair(1,ibc(jbra8))                                      4d22s22
          iak=iapair(1,ibc(jket8))                                      4d22s22
          call oned2ints(ibc(jpair+2),iak,iab,                          4d22s22
     $     ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),cartb(3),       2d4s16
     $     ibc(jket),bc(jket2),bc(jket3),cartk,cartk(2),cartk(3),       2d4s16
     $     idxyz,ib1,ib2,ib3,ib4,ib5,ib6,nneedzm,nneedz,nneedz2,        6d20s16
     $         nneedz2m,npass,idrb,idrk,                                2d27s23
     $     iapair(1,iatom),natom,idorel,xcart,atnum,iatom,              6d20s16
     $         ixyz,mynowprog,i,idersign,bc,ibc)                        11d10s22
          ibu=ib1                                                        5d7s10
          itmpu=itmp1                                                    5d7s10
          do ipass=1,npass                                               5d7s10
           do ik=1,nket                                                  5d3s10
            do ib=1,nbra                                                 5d3s10
             iad1=ik-1+nket*(ib-1)                                       5d3s10
             iad2=ib-1+nbra+nbraa*(ik-1+nket)                            5d3s10
             bc(itmpu+iad2)=bc(ibu+iad1)                                 8d20s15
            end do                                                       5d7s10
           end do                                                        5d3s10
           ibu=ibu+nbra*nket                                            4d5s16
           itmpu=itmpu+nbraa*nketa                                      4d5s16
          end do                                                         5d3s10
         end if                                                          5d3s10
        end if                                                           5d3s10
        if(nketa.ne.nket)then
         itmpu=itmp1                                                     5d7s10
         do ipass=1,npass                                                5d7s10
          do ik=1,nket                                                   5d3s10
           do ib=1,nbraa                                                 5d3s10
            iad1=ib-1+nbraa*(ik-1)                                       5d3s10
            iad2=iad1+nbraa*nket                                         5d3s10
            sum=srh*(bc(itmpu+iad1)+bc(itmpu+iad2))                      5d3s10
            dif=srh*(-bc(itmpu+iad1)+bc(itmpu+iad2))                     5d3s10
            bc(itmpu+iad1)=sum                                           5d3s10
            bc(itmpu+iad2)=dif                                           5d3s10
           end do                                                        5d7s10
          end do                                                         5d3s10
          itmpu=itmpu+nbraa*nketa                                       4d5s16
         end do                                                          5d3s10
        end if                                                           5d3s10
        if(nbraa.ne.nbra)then
         itmpu=itmp1                                                     5d7s10
         do ipass=1,npass                                                5d7s10
          do ik=1,nketa                                                   5d3s10
           do ib=1,nbra                                                  5d3s10
            iad1=ib-1+nbraa*(ik-1)                                       5d3s10
            iad2=iad1+nbra                                               5d3s10
            sum=srh*(bc(itmpu+iad1)+bc(itmpu+iad2))                      5d3s10
            dif=srh*(-bc(itmpu+iad1)+bc(itmpu+iad2))                     5d3s10
            bc(itmpu+iad1)=sum                                           5d3s10
            bc(itmpu+iad2)=dif                                           5d3s10
           end do                                                        5d7s10
          end do                                                         5d3s10
          itmpu=itmpu+nbraa*nketa                                       4d5s16
         end do                                                          5d3s10
        end if                                                           5d3s10
        itmpu=itmp1                                                      5d7s10
        if(idorel.eq.0)then                                             8d20s15
         do ipass=1,npass                                                 5d7s10
c
c     npass=1, npass=1: we've got full potential ders
c     npass=6, ipass=1: d2 bra of S, ipass=2: d2 ket of S, ipass=3:
c     d bra d ket of S, ipass=4: d2 bra of T, ipass=5: d2 ket of T                      5d17s16
c     and ipass=6: d bra d ket of T.                                    6d20s16
c     for h0, we need sum of bra and ket ders. For (n|d2/dq2|m), we just6d20s16
c     need ket ders. However, since we only do a triangle, we need
c     to -transpose bra ders to get the full matrix.
c     note: d2 is totally symmetric, but d is not. This means that for
c     ipass=3&6, the contribuion will be diagonal in symmetry, but fall
c     in a block not having the same symmetry of the current block.
c
          do ik=1,nketa                                                   5d3s10
           ikk=ik+ibc(jket4)                                              5d3s10
           isk=isstor(ikk)                                                5d3s10
           do ib=1,nbraa                                                  5d3s10
            ibb=ib+ibc(jbra4)                                             5d3s10
            iad1=ib-1+nbraa*(ik-1)                                        5d3s10
            isb=isstor(ibb)                                             1d11s16
            if(isk.eq.isb)then                                          6d20s16
             iaddk1=iso(isb)+ibstor(ibb)+nbasisp(isb)*(ibstor(ikk)-1)   4d4s22
             iaddk2=iso(isk)+ibstor(ikk)+nbasisp(isk)*(ibstor(ibb)-1)   4d4s22
             in=min(ibstor(ikk),ibstor(ibb))                              5d3s10
             ix=max(ibstor(ikk),ibstor(ibb))                              5d3s10on
             iad=iso(isk)+((ix*(ix-1))/2)+in
             if(ipass.eq.1.and.npass.eq.1.or.ipass.ge.4)then            6d20s16
              if(ibstor(ikk).le.ibstor(ibb).or.                         1d27s16
     $            (isb.ne.isk.and.ibc(jpair).ne.ibc(jpair+1)))then      1d27s16
               h0(iad)=h0(iad)+bc(itmpu+iad1)                               5d3s10
              end if                                                      1d17s12
             else if(ipass.le.3)then                                    5d12s16
              if(ibstor(ikk).le.ibstor(ibb).or.                         1d27s16
     $             (isb.ne.isk.and.ibc(jpair).ne.ibc(jpair+1)))then      4d5s16
               ovr(iad)=ovr(iad)+bc(itmpu+iad1)                         4d5s16
               if(ipass.eq.2)then
                ovrdk(iaddk1)=ovrdk(iaddk1)+bc(itmpu+iad1)
               end if
               if(ipass.eq.1.and.ibstor(ikk).ne.ibstor(ibb))then        12d23s16
                ovrdk(iaddk2)=ovrdk(iaddk2)+bc(itmpu+iad1)
               end if                                                   12d23s16
              end if                                                    4d5s16
             end if                                                       5d7s10
            else                                                          5d3s10
             sz=abs(bc(itmpu+iad1))                                       5d7s10
             if(sz.gt.1d-10.and.npass.eq.2)then                                          5d3s10
              write(6,*)('symmetry transformation failure!!! ')           5d3s10
              write(6,*)ib,ik,ibb,ikk,sz                                  5d7s10
              write(6,*)('ipass,npass: '),ipass,npass                     5d7s10
              call prntm2(bc(itmpu),nbraa,nketa,nbraa)                    5d3s10
              stop                                                        5d3s10
             end if                                                       5d3s10
            end if                                                        5d3s10
           end do                                                         5d7s10
          end do                                                          5d3s10
          itmpu=itmpu+nbraa*nketa                                       4d5s16
         end do                                                           5d3s10
        else                                                            8d20s15
         ia=ibc(jpair+2)                                                 2d22s10
         do ipass=1,npass                                                 5d7s10
          do ik=1,nketa                                                   5d3s10
           ikk=ik+ibc(jket4)                                              5d3s10
           isk=isstor(ikk)                                                5d3s10
           do ib=1,nbraa                                                  5d3s10
            ibb=ib+ibc(jbra4)                                             5d3s10
            iad1=ib-1+nbraa*(ik-1)                                        5d3s10
            isb=isstor(ibb)                                             2d2s16
            iaddk1=iso(isb)+ibstor(ibb)+nbasisp(isb)*2*(ibstor(ikk)-1)  4d4s22
            iaddk2=iso(isk)+ibstor(ikk)+nbasisp(isk)*2*(ibstor(ibb)-1)  4d4s22
            iaddk1ss=iaddk1+nbasisp(isb)*(1+2*nbasisp(isb))             4d4s22
            iaddk2ss=iaddk2+nbasisp(isb)*(1+2*nbasisp(isb))             4d4s22
            if(isk.eq.isb)then                                          2d2s16
             inl=min(ibstor(ikk),ibstor(ibb))                              5d3s10
             ixl=max(ibstor(ikk),ibstor(ibb))                              5d3s10
             iadll=iso(isk)+((ixl*(ixl-1))/2)+inl
             ins=min(ibstor(ikk),ibstor(ibb))+nbasisp(isk)              4d4s22
             ixs=max(ibstor(ikk),ibstor(ibb))+nbasisp(isk)              4d4s22
             iadss=iso(isk)+((ixs*(ixs-1))/2)+ins
             ix=ibstor(ibb)+nbasisp(isk)                                4d4s22
             in=ibstor(ikk)
             iadsl=iso(isk)+((ix*(ix-1))/2)+in
             ix=ibstor(ikk)+nbasisp(isk)                                4d4s22
             in=ibstor(ibb)
             iadls=iso(isk)+((ix*(ix-1))/2)+in
c
c     npass=6: ipass=1 is d2 bra S, ipass=2, d2 ket S,overlap, ipass=3  6d23s16
c           d bra S d ket, ipass=4 d2 bra T, ipass=5 d2 ket T, ipass=6  6d23s16
c           d bra T d ket                                               6d23s16
c     npass=1: Vnuc if ia le natom, dVnucd otherwise.
c     overlap goes into LL part of ovr,
c     t goes into LS part of h0, -t goes into SS part of h0,
c     and ascale*ascale*t goes into SS part of
c     ovr. Vnuc goes into LL part of h0 and dVnucd into SS part of h0.
c
             if(npass.gt.1)then
              if(ipass.le.3)then                                        6d23s16
               if(ibstor(ikk).le.ibstor(ibb).or.                         2d2s16
     $             (isb.ne.isk.and.ibc(jpair).ne.ibc(jpair+1)))then     2d2s16
                ovr(iadll)=ovr(iadll)+bc(itmpu+iad1)                     6d23s16
                if(ipass.eq.2)then                                       6d23s16
                 ovrdk(iaddk1)=ovrdk(iaddk1)+bc(itmpu+iad1)              6d23s16
                else if(ipass.eq.1.and.ibstor(ikk).ne.ibstor(ibb))then  3d7s23
                  ovrdk(iaddk2)=ovrdk(iaddk2)+bc(itmpu+iad1)              6d23s16
                end if                                                   6d23s16
               end if                                                    8d20s15
              else                                                       6d23s16
               h0(iadsl)=h0(iadsl)+bc(itmpu+iad1)                        6d23s16
               if(iadsl.ne.iadls)then                                   3d6s23
                h0(iadls)=h0(iadls)+bc(itmpu+iad1)                        6d23s16
               end if                                                   3d6s23
               if(ibstor(ikk).le.ibstor(ibb).or.
     $             (isb.ne.isk.and.ibc(jpair).ne.ibc(jpair+i)))then     2d2s16
                ovr(iadss)=ovr(iadss)+bc(itmpu+iad1)*ascale2             6d23s16
                h0(iadss)=h0(iadss)-bc(itmpu+iad1)                       8d20s15
                if(ipass.eq.4)then                                       6d23s16
                 ovrdk(iaddk1ss)=ovrdk(iaddk1ss)+bc(itmpu+iad1)*ascale2  6d23s16
                else if(ipass.eq.5)then                                  6d23s16
                 ovrdk(iaddk2ss)=ovrdk(iaddk2ss)+bc(itmpu+iad1)*ascale2  6d23s16
                end if                                                   6d23s16
               end if                                                    8d20s15
              end if                                                     6d23s16
             else if(ia.le.natom)then                                   8d20s15
              if(ibstor(ikk).le.ibstor(ibb).or.
     $             (isb.ne.isk.and.ibc(jpair).ne.ibc(jpair+i)))then     2d2s16
               h0(iadll)=h0(iadll)+bc(itmpu+iad1)
               if(h0(iadll).ne.h0(iadll))then
                write(6,*)('got nan !!! in h0grad ')
               end if
              end if
             else                                                       8d20s15
              if(ibstor(ikk).le.ibstor(ibb).or.
     $             (isb.ne.isk.and.ibc(jpair).ne.ibc(jpair+i)))then     2d2s16
               h0(iadss)=h0(iadss)+ascale*bc(itmpu+iad1)                8d20s15
              end if                                                    8d20s15
             end if                                                     8d20s15
            else                                                          5d3s10
             sz=abs(bc(itmpu+iad1))                                       5d7s10
             if(sz.gt.1d-10.and.npass.eq.2)then                                          5d3s10
              write(6,*)('symmetry transformation failure!!! ')           5d3s10
              write(6,*)ib,ik,ibb,ikk,sz                                  5d7s10
              write(6,*)('ipass,npass: '),ipass,npass                     5d7s10
              call prntm2(bc(itmpu),nbraa,nketa,nbraa)                    5d3s10
              stop                                                        5d3s10
             end if                                                       5d3s10
            end if                                                        5d3s10
           end do                                                         5d7s10
          end do                                                          5d3s10
          itmpu=itmpu+nbraa*nketa                                       4d5s16
         end do                                                           5d3s10
        end if                                                          8d20s15
        ibcoff=itmp1                                                     5d3s10
      end do                                                            2d19s10
      call second(time3)
      call dws_sync                                                     2d22s10
      call dws_gsumf(ovr,nbb)                                           5d3s10
      call dws_gsumf(h0,nbb)                                            5d3s10
      call dws_gsumf(ovrdk,nbb)                                         4d4s16
      call second(time4)
      if(idwsdeb.gt.10)then                                             3d16s12
      do isb=1,nsymb                                                    5d3s10
       write(6,*)('for symmetry block '),isb                            5d3s10
       write(6,*)('overlap: '),iso(isb),ovr(iso(isb)+1),
     $      loc(ovr(iso(isb)+1))
       if(idorel.eq.0)then                                              8d20s15
        nbasz=nbasisp(isb)                                              4d4s22
       else                                                             8d20s15
        nbasz=nbasisp(isb)*2                                            4d4s22
       end if                                                           8d20s15
       call mpprnt2(ovr(iso(isb)+1),nbasz)                              8d20s15
       if(nsymb.eq.1)then                                               4d25s18
        write(6,*)('ovr offset: '),iso(isb),ibcoff
        itmp=ibcoff                                                     4d25s18
        ibcoff=itmp+nbasz*nbasz                                         4d25s18
        ii=iso(isb)+1                                                   4d25s18
        call enough('parah0d2.  3',bc,ibc)
        do i=0,nbasz-1                                                  4d25s18
         do j=0,i                                                       4d25s18
          ij=itmp+i+nbasz*j                                             4d25s18
          ji=itmp+j+nbasz*i                                             4d25s18
          bc(ij)=ovr(ii)                                                4d25s18
          bc(ji)=ovr(ii)                                                4d25s18
          ii=ii+1                                                       4d25s18
         end do                                                         4d25s18
        end do                                                          4d25s18
        write(6,*)('squared ')
        call prntm2(bc(itmp),nbasz,nbasz,nbasz)
       call mpprnt2(ovr(iso(isb)+1),nbasz)                              8d20s15
        call printao(bc(itmp),bc,ibc)                                          4d25s18
        ibcoff=itmp                                                     4d25s18
       end if                                                           4d25s18
       iuse=iso(isb)
       write(6,*)('h0: '),iso(isb)+1                                                2d22s10
       call mpprnt2(h0(iso(isb)+1),nbasz)                               8d20s15
       if(nsymb.eq.1)then                                               4d25s18
        itmp=ibcoff                                                     4d25s18
        ibcoff=itmp+nbasz*nbasz                                         4d25s18
        ii=iso(isb)+1                                                   4d25s18
        call enough('parah0d2.  4',bc,ibc)
        do i=0,nbasz-1                                                  4d25s18
         do j=0,i                                                       4d25s18
          ij=itmp+i+nbasz*j                                             4d25s18
          ji=itmp+j+nbasz*i                                             4d25s18
          bc(ij)=h0(ii)                                                 4d25s18
          bc(ji)=h0(ii)                                                 4d25s18
          ii=ii+1                                                       4d25s18
         end do                                                         4d25s18
        end do                                                          4d25s18
        call printao(bc(itmp),bc,ibc)                                          4d25s18
        ibcoff=itmp                                                     4d25s18
       end if                                                           4d25s18
       write(6,*)('ovrdk '),iso(isb)
       call prntm2(ovrdk(iso(isb)+1),nbasz,nbasz,nbasz)
       if(nsymb.eq.1)then
        call printao(ovrdk(iso(isb)+1),bc,ibc)                                 4d25s18
       end if                                                           4d25s18
       rms=0d0
       end do
      end if                                                            3d16s12
      ibcoff=ibcoffo                                                    2d19s10
      return
      end                                                               2d19s10
