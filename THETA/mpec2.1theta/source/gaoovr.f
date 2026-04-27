c mpec2.1 version theta. For copyright and Disclaimers, see start.f
      subroutine gaoovr(ibdat0,ibdat1,iorbao0,iorbao1,ibstor,isstor,    12d5s24
     $     nsymb,nbasdws,nbasb,ncomp,ngaus,isym,iapair,ascale,bc,ibc)   12d5s24
      implicit real*8 (a-h,o-z)
c
c     compute overlap of ao orbitals at 2 different geometries
c
      integer*8 ibstor(*),isstor(*)                                     12d5s24
      dimension iorbao0(*),iorbao1(*),nbasdws(*),nbasb(*),isym(3,*),    12d5s24
     $     iapair(3,*),iovr(8),cartk(3),cartb(3)                        12d5s24
      data loop,loopx/0,33080000/
      include 'common.store'                                            12d5s24
      integer*8 ibcoffo,iovr,itmp,itmp1
      ibcoffo=ibcoff
      srh=sqrt(0.5d0)                                                   5d3s10
      ascale2=ascale*2d0                                                8d20s15
      do isb=1,nsymb                                                    12d5s24
       nbasbb=nbasb(isb)*ncomp                                          12d5s24
       iovr(isb)=ibcoff                                                 12d5s24
       ibcoff=iovr(isb)+nbasbb*nbasbb                                   12d5s24
      end do                                                            12d5s24
      call enough('iovr.gaoovr',bc,ibc)                                 12d5s24
      do iz=iovr(1),ibcoff-1                                            12d5s24
       bc(iz)=0d0                                                       12d5s24
      end do                                                            12d5s24
      do iket=1,ngaus                                                   12d5s24
       jket=ibdat1+iket-1                                               12d5s24
       jket2=jket+ngaus                                                 12d5s24
       jket3=jket2+ngaus                                                2d19s10
       jket4=jket3+ngaus                                                2d19s10
       jket5=jket4+ngaus                                                2d19s10
       jket6=jket5+ngaus                                                2d19s10
       jket7=jket6+ngaus                                                2d19s10
       jket8=jket7+ngaus                                                5d4s10
       do ibra=1,ngaus                                                   12d5s24
        jbra=ibdat0+ibra-1                                               12d5s24
        jbra2=jbra+ngaus                                                 2d19s10
        jbra3=jbra2+ngaus                                                2d19s10
        jbra4=jbra3+ngaus                                                2d19s10
        jbra5=jbra4+ngaus                                                2d19s10
        jbra6=jbra5+ngaus                                                2d19s10
        jbra7=jbra6+ngaus                                                2d19s10
        jbra8=jbra7+ngaus                                                5d4s10
        nbra=2*ibc(jbra)+1                                              5d3s10
        nbraa=nbra                                                      5d3s10
        if(iapair(1,ibc(jbra8)).gt.0)nbraa=nbra*2                       5d3s10
        nket=2*ibc(jket)+1                                              5d3s10
        nketa=nket                                                      5d3s10
        if(iapair(1,ibc(jket8)).gt.0)nketa=nket*2                       5d3s10
        itmp1=ibcoff                                                    5d3s10
        itmp2=itmp1+nbraa*nketa                                         5d3s10
        ibcoff=itmp2+nbraa*nketa                                        5d3s10
        call enough('gaoovr.  2',bc,ibc)
        call onei(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),    2d19s10
     $      bc(jbra7),ibc(jbra4),ibc(jket),bc(jket2),bc(jket3),         2d19s10
     $      bc(jket5),bc(jket6),bc(jket7),ibc(jket4),dum,dum,nsymb,ib1, 12d5s24
     $      ib2,bc,ibc)                                                 11d9s22
        ibu=ib1                                                         5d7s10
        itmpu=itmp1                                                     5d7s10
        do ipass=1,2                                                    12d5s24
         do ik=1,nket                                                   5d7s10
          do ib=1,nbra                                                  5d7s10
           iad1=ik-1+nket*(ib-1)                                        5d7s10
           iad2=ib-1+nbraa*(ik-1)                                       5d7s10
           bc(itmpu+iad2)=bc(ibu+iad1)                                  8d20s15
          end do                                                        5d7s10
         end do                                                         5d7s10
         ibu=ib2                                                        5d7s10
         itmpu=itmp2                                                    5d7s10
        end do                                                          5d7s10
        if(nketa.ne.nket)then                                           5d3s10
         cartk(1)=bc(jket5)*dfloat(isym(1,iapair(2,ibc(jket8))))        5d5s10
         cartk(2)=bc(jket6)*dfloat(isym(2,iapair(2,ibc(jket8))))        5d5s10
         cartk(3)=bc(jket7)*dfloat(isym(3,iapair(2,ibc(jket8))))        5d5s10
         call onei(ibc(jbra),bc(jbra2),bc(jbra3),bc(jbra5),bc(jbra6),   5d3s10
     $     bc(jbra7),ibc(jbra4),ibc(jket),bc(jket2),bc(jket3),          2d19s10
     $         cartk,cartk(2),cartk(3),ibc(jket4),dum,dum,nsymb,ib1,    12d5s24
     $         ib2,bc,ibc)                                              11d9s22
         ibu=ib1                                                        5d7s10
         itmpu=itmp1                                                    5d7s10
         do ipass=1,2                                                   12d5s24
          do ik=1,nket                                                    5d3s10
           do ib=1,nbra                                                   5d3s10
            iad1=ik-1+nket*(ib-1)                                         5d3s10
            iad2=ib-1+nbraa*(ik-1+nket)                                  5d3s10
            bc(itmpu+iad2)=bc(ibu+iad1)                                 8d20s15
           end do                                                         5d3s10
          end do                                                          5d3s10
          ibu=ib2                                                       5d7s10
          itmpu=itmp2                                                   5d7s10
         end do                                                         5d7s10
        end if                                                          5d3s10
        if(nbraa.ne.nbra)then                                           5d3s10
         cartb(1)=bc(jbra5)*dfloat(isym(1,iapair(2,ibc(jbra8))))        5d5s10
         cartb(2)=bc(jbra6)*dfloat(isym(2,iapair(2,ibc(jbra8))))        5d5s10
         cartb(3)=bc(jbra7)*dfloat(isym(3,iapair(2,ibc(jbra8))))        5d5s10
         call onei(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),        5d3s10
     $        cartb(3),ibc(jbra4),ibc(jket),bc(jket2),bc(jket3),        5d3s10
     $     bc(jket5),bc(jket6),bc(jket7),ibc(jket4),dum,dum,nsymb,ib1,  12d5s24
     $     ib2,bc,ibc)                                                  11d9s22
         ibu=ib1                                                         5d7s10
         itmpu=itmp1                                                     5d7s10
         do ipass=1,2                                                   12d5s24
          do ik=1,nket                                                    5d3s10
           do ib=1,nbra                                                   5d3s10
            iad1=ik-1+nket*(ib-1)                                         5d3s10
            iad2=ib-1+nbra+nbraa*(ik-1)                                   5d3s10
            bc(itmpu+iad2)=bc(ibu+iad1)                                  8d20s15
           end do                                                        5d7s10
          end do                                                         5d3s10
          ibu=ib2                                                        5d7s10
          itmpu=itmp2                                                    5d7s10
         end do                                                          5d3s10
         if(nketa.ne.nket)then                                           5d3s10
         call onei(ibc(jbra),bc(jbra2),bc(jbra3),cartb,cartb(2),        12d5s24
     $        cartb(3),ibc(jbra4),ibc(jket),bc(jket2),bc(jket3),        5d3s10
     $     cartk,cartk(2),cartk(3),ibc(jket4),dum,dum,nsymb,ib1,        12d5s24
     $     ib2,bc,ibc)                                                  11d9s22
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
           ibu=ib2                                                       5d7s10
           itmpu=itmp2                                                   5d7s10
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
          itmpu=itmp2                                                    5d7s10
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
          itmpu=itmp2                                                    5d7s10
         end do                                                          5d3s10
        end if                                                           5d3s10
        if(idorel.eq.0)then                                             8d20s15
         do ik=1,nketa                                                   5d3s10
          ikk=ik+ibc(jket4)                                              5d3s10
          isk=isstor(ikk)                                                5d3s10
          do ib=1,nbraa                                                  5d3s10
           ibb=ib+ibc(jbra4)                                             5d3s10
           iad1=ib-1+nbraa*(ik-1)                                        5d3s10
           if(isk.eq.isstor(ibb))then                                    5d3s10
            iad=iovr(isk)+ibstor(ibb)-1+nbasb(isk)*(ibstor(ikk)-1)      12d5s24
            bc(iad)=bc(itmp1+iad1)                                      12d5s24
           else                                                          5d3s10
            sz=abs(bc(itmp1+iad1))                                       5d7s10
            if(sz.gt.1d-10)then                                          5d3s10
             write(6,*)('symmetry transformation failure!!! ')           5d3s10
             write(6,*)ib,ik,ibb,ikk,sz                                  5d7s10
             call prntm2(bc(itmp1),nbraa,nketa,nbraa)                    5d3s10
             stop 'gaoovr'                                              12d5s24
            end if                                                       5d3s10
           end if                                                        5d3s10
          end do                                                          5d3s10
         end do                                                           5d3s10
        else                                                            8d20s15
         do ipass=1,2                                                   12d5s24
          do ik=1,nketa                                                   5d3s10
           ikk=ik+ibc(jket4)                                              5d3s10
           isk=isstor(ikk)                                                5d3s10
           do ib=1,nbraa                                                  5d3s10
            ibb=ib+ibc(jbra4)                                             5d3s10
            iad1=ib-1+nbraa*(ik-1)                                        5d3s10
            if(isk.eq.isstor(ibb))then                                    5d3s10
             iadll=iovr(isk)+ibstor(ibb)-1+nbasb(isk)*2*(ibstor(ikk)-1) 12d5s24
             iadss=iovr(isk)+ibstor(ibb)-1+nbasb(isk)*(1                12d5s24
     $            +2*(ibstor(ikk)+nbasb(isk)-1))                        12d5s24
             if(ipass.eq.1)then                                         12d5s24
              bc(iadll)=bc(itmpu+iad1)                                  12d5s24
             else
              bc(iadss)=bc(itmpu+iad1)*ascale2                          12d5s24
             end if
            else                                                          5d3s10
             sz=abs(bc(itmpu+iad1))                                       5d7s10
             if(sz.gt.1d-10.and.npass.eq.2)then                                          5d3s10
              write(6,*)('symmetry transformation failure!!! ')           5d3s10
              write(6,*)ib,ik,ibb,ikk,sz                                  5d7s10
              write(6,*)('ipass,npass: '),ipass,npass                     5d7s10
              call prntm2(bc(itmpu),nbraa,nketa,nbraa)                    5d3s10
              stop 'gaoovr'                                                        5d3s10
             end if                                                       5d3s10
            end if                                                        5d3s10
           end do                                                         5d7s10
          end do                                                          5d3s10
          itmpu=itmp2                                                     5d7s10
         end do                                                           5d3s10
        end if                                                          8d20s15
        ibcoff=itmp1                                                     5d3s10
       end do                                                           12d5s24
      end do                                                            12d5s24
      write(6,*)('we have all of integrals ...')
      do isb=1,nsymb                                                    12d5s24
       if(nbasb(isb).gt.0)then                                          12d5s24
        write(6,*)('for symmetry block '),isb                           12d5s24
        nbasbb=nbasb(isb)*ncomp                                         12d5s24
        write(6,*)('ao overlap matrix ')                                12d5s24
        call prntm2(bc(iovr(isb)),nbasbb,nbasbb,nbasbb)                 12d5s24
        itmp=ibcoff                                                     12d5s24
        ibcoff=itmp+nbasbb*nbasdws(isb)                                 12d5s24
        call enough('gaoovr.tmp',bc,ibc)                                12d5s24
        call dgemm('n','n',nbasbb,nbasdws(isb),nbasbb,1d0,bc(iovr(isb)),12d5s24
     $       nbasbb,bc(iorbao1(isb)),nbasbb,0d0,bc(itmp),nbasbb)        12d5s24
        do i=0,nbasbb-1                                                 12d5s24
         do j=0,nbasdws(isb)-1                                          12d5s24
          ji=iovr(isb)+j+nbasdws(isb)*i                                 12d5s24
          ij=itmp+i+nbasbb*j                                            12d5s24
          bc(ji)=bc(ij)                                                 12d5s24
         end do                                                         12d5s24
        end do                                                          12d5s24
        call dgemm('n','n',nbasdws(isb),nbasdws(isb),nbasbb,1d0,        12d5s24
     $       bc(iovr(isb)),nbasdws(isb),bc(iorbao0(isb)),nbasbb,0d0,    12d5s24
     $       bc(itmp),nbasdws(isb))                                     12d5s24
        call prntm2(bc(itmp),nbasdws(isb),nbasdws(isb),nbasdws(isb))    12d5s24
        ibcoff=itmp                                                     12d5s24
       end if                                                           12d5s24
      end do                                                            12d5s24
      ibcoff=ibcoffo
      return
      end
