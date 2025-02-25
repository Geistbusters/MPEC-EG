c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine buildhess(ioooo,ionex,jmats,kmats,h0mo,noc,ihess,
     $     nvirtc,ipsym,multh,idwsdeb,bc,ibc)                           11d14s22
      implicit real*8 (a-h,o-z)
c
c     build orbital rotation hessian matrix
c
      dimension h0mo(1),ioooo(1),ionex(1),jmats(1),kmats(1),noc(1),
     $     ihess(8,8),nvirtc(1),multh(8,8),ipt(8)
      include "common.hf"
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      common/unitcm/iunit                                               11d9s17
      if(idwsdeb.gt.10)then
       write(6,*)('in buildhess '),nvirtc(1),ipsym,loc(ipsym)
       write(6,*)('ibcoff, nsymb '),ibcoff,nsymb
      end if
      ioffh=0
      do isb=1,nsymb
       ipt(isb)=ioffh
       nh=noc(isb)+nvirtc(isb)
       ioffh=ioffh+nh*nh
      end do
      do isb=1,nsymb
       do isa=1,isb
        if(idwsdeb.gt.10)write(6,*)('for isa,isb = '),isa,isb
c
c     the symmetry label isb,isa correspond to the bra and ket
c     occupied orbital symmetries. note that I've coded up things so the
c     order of occupied orbital indicies is bra,ket but the order of
c     virtual orbital indicies is ket,bra.
c     because of the distribution across processors, we need
c     isb,isa to be the virtual labels, which is only different when
c     the perturbation is not totally symmetric.                        3d29s16
c
        ihess(isa,isb)=-1                                               12d5s16
        ihess(isb,isa)=-1                                               12d5s16
        isav=isa                                                        3d29s16
        isbv=isb                                                        3d29s16
        isao=multh(isa,ipsym)                                           3d28s16
        isbo=multh(isb,ipsym)                                           3d28s16
        if(noc(isbo)*noc(isao).gt.0)then                                  3d28s16
         ihess(isa,isb)=ibcoff                                           3d28s16
         call ilimts(nvirtc(isb),nvirtc(isa),mynprocg,mynowprog,ilh,    3d29s16
     $        ihh,                                                      3d28s16
     $        i1sh,i1eh,i2sh,i2eh)                                      3d28s16
         call ilimts(noc(isb),noc(isa),mynprocg,mynowprog,il4o,         1d3s17
     $        ih4o,i1s4o,i1e4o,i2s4o,i2e4o)                             1d3s17
         nhere=ihh+1-ilh                                                3d28s16
         nrowh=noc(isao)*noc(isbo)                                        3d28s16
         ibcoff=ihess(isa,isb)+nhere*noc(isao)*noc(isbo)                  3d28s16
         if(idwsdeb.gt.10)
     $        write(6,*)('space for ihess: '),noc(isao)*noc(isbo),nhere
         ibtop=ibcoff                                                   3d29s16
         call enough('buildhess.  1',bc,ibc)
         do is=1,nsdlk
          if(isblk(1,is).eq.isao.and.isblk(2,is).eq.isbo.and.
     $         isblk(3,is).eq.isbv.and.isblk(4,is).eq.isav)then         3d29s16
           ii=jmats(is)
           if(isao.eq.isbo)then
            nrowj=(noc(isbo)*(noc(isbo)+1))/2
           else
            nrowj=noc(isbo)*noc(isao)
           end if
           kk=ihess(isa,isb)
           i10=i1sh
           i1n=nvirtc(isbv)                                             3d29s16
           do i2=i2sh,i2eh
            if(i2.eq.i2eh)i1n=i1eh
            do i1=i10,i1n
             if(isa.eq.isb)then
              do i4=0,noc(isbo)-1
               do i3=0,noc(isao)-1
                ix=max(i3,i4)
                in=min(i3,i4)
                iadj=ii+((ix*(ix+1))/2)+in
                bc(kk)=-4d0*bc(iadj)
                kk=kk+1
               end do
              end do
              ii=ii+nrowj
             else
              do i4=0,noc(isbo)-1
               do i3=0,noc(isao)-1
                bc(kk)=-4d0*bc(ii)
                kk=kk+1
                ii=ii+1
               end do
              end do
             end if
            end do
            i10=1
           end do
           go to 1
          end if
          if(isblk(1,is).eq.isbo.and.isblk(2,is).eq.isao.and.
     $         isblk(3,is).eq.isbv.and.isblk(4,is).eq.isav)then         3d29s16
           ii=jmats(is)
           kk=ihess(isa,isb)
           i10=i1sh
           i1n=nvirtc(isbv)                                             3d29s16
           do i2=i2sh,i2eh
            if(i2.eq.i2eh)i1n=i1eh
            do i1=i10,i1n
             do i4=0,noc(isbo)-1
              do i3=0,noc(isao)-1
               iadj=ii+i4+noc(isbo)*i3
               bc(kk)=-4d0*bc(iadj)
               kk=kk+1
              end do
             end do
             ii=ii+nrowh
            end do
            i10=1
           end do
           go to 1
          end if
         end do
         if(min(nvirt(isbv),nvirt(isav)).gt.0)then                      11d28s22
          write(6,*)('in buildhess')
          write(6,*)('could not find jmat for '),isao,isbo,isbv,isav       3d29s16
          write(6,*)('nvirt '),nvirt(isbv),nvirt(isav)
          do is=1,nsdlk
           write(6,*)(isblk(j,is),j=1,4)
          end do
          call dws_sync
          call dws_finalize
          stop
         end if                                                         11d28s22
    1    continue
         do is=1,nsdlkk
          if(isblkk(1,is).eq.isao.and.isblkk(2,is).eq.isbo.and.
     $         isblkk(3,is).eq.isbv.and.isblkk(4,is).eq.isav)then       3d29s16
           kk=ihess(isa,isb)
           ii=kmats(is)
           i10=i1sh
           i1n=nvirtc(isbv)                                             3d29s16
           do i2=i2sh,i2eh
            if(i2.eq.i2eh)i1n=i1eh
            do i1=i10,i1n
             do i4=0,noc(isbo)-1
              do i3=0,noc(isao)-1
               bc(kk)=bc(kk)+16d0*bc(ii)
               kk=kk+1
               ii=ii+1
              end do
             end do
            end do
            i10=1
           end do
           go to 2
          end if
         end do
         if(min(nvirt(isbv),nvirt(isav)).gt.0)then                      11d28s22
          write(6,*)('in buildhess')
          write(6,*)('could not find kmat for '),isao,isbo,isbv,isav
          do is=1,nsdlkk
           write(6,*)(isblkk(j,is),j=1,4)
          end do
          call dws_sync
          call dws_finalize
          stop
         end if                                                         11d28s22
    2    continue
         do is=1,nsdlkk
          if(isblkk(1,is).eq.isbo.and.isblkk(2,is).eq.isao.and.
     $         isblkk(3,is).eq.isbv.and.isblkk(4,is).eq.isav)then       3d29s16
           kk=ihess(isa,isb)
           ii=kmats(is)
           i10=i1sh
           i1n=nvirtc(isbv)                                             3d29s16
           do i2=i2sh,i2eh
            if(i2.eq.i2eh)i1n=i1eh
            do i1=i10,i1n
             do i4=0,noc(isbo)-1
              do i3=0,noc(isao)-1
               iadk=ii+i4+noc(isbo)*i3
               bc(kk)=bc(kk)-4d0*bc(iadk)
               kk=kk+1
              end do
             end do
             ii=ii+nrowh
            end do
            i10=1
           end do
           go to 3
          end if
         end do
         if(min(nvirt(isbv),nvirt(isav)).gt.0)then                      11d28s22
          write(6,*)('in buildhess')
          write(6,*)('could not find kmatT for '),isao,isbo,isbv,isav      3d29s16
          do is=1,nsdlkk
           write(6,*)(isblkk(j,is),j=1,4)
          end do
          call dws_sync
          call dws_finalize
          stop
         end if                                                         11d28s22
    3    continue
         if(isa.eq.isb)then
          noo=(noc(isao)*(noc(isao)+1))/2
          nvv=(nvirtc(isav)*(nvirtc(isav)+1))/2                         3d29s16
          ioop=ibcoff
          ivv=ioop+noo
          ibcoff=ivv+nvv
          call enough('buildhess.  2',bc,ibc)
          do i=0,noo+nvv-1
           bc(ioop+i)=0d0
          end do
          do is=1,nsdlk
           if(isblk(1,is).eq.isao.and.isblk(2,is).eq.isao)then
            iscv=isblk(3,is)
            isco=multh(iscv,ipsym)
            call ilimts(noc(iscv),noc(iscv),mynprocg,mynowprog,il,
     $           ih,i1s,i1e,i2s,i2e)
            ii=ioooo(is)
            i10=i1s
            i1n=noc(iscv)
            do i2=i2s,i2e
             if(i2.eq.i2e)i1n=i1e
             do i1=i10,i1n
              if(i1.eq.i2)then
               do i34=0,noo-1
                bc(ioop+i34)=bc(ioop+i34)-8d0*bc(ii+i34)
               end do
              end if
              ii=ii+noo
             end do
             i10=1
            end do
           end if
          end do
          if(idwsdeb.gt.10)then
          write(6,*)('oop after i1=i2 sum ')
          write(6,*)(bc(ioop+i),i=0,noo-1)
          end if
          do isc=1,nsymb
           if(noc(isc).gt.0)then
            do is=1,nsdlk
             if(isblk(1,is).eq.isc.and.isblk(2,is).eq.isao.and.
     $          isblk(3,is).eq.isc)then                                 3d29s16
              call ilimts(noc(isc),noc(isao),mynprocg,mynowprog,il,     3d29s16
     $           ih,i1s,i1e,i2s,i2e)
              ii=ioooo(is)
              if(isc.eq.isao)then
               noo=(noc(isao)*(noc(isao)+1))/2
              else
               noo=noc(isao)*noc(isc)
              end if
              i10=i1s
              i1n=noc(isc)                                              3d29s16
              do i2=i2s,i2e
               if(i2.eq.i2e)i1n=i1e
               i2m=i2-1
               do i1=i10,i1n
                i1m=i1-1
                if(isao.eq.isc)then
                 do il=i2m,noc(isao)-1
                  iad=ioop+((il*(il+1))/2)+i2m
                  ix=max(il,i1m)
                  in=min(il,i1m)
                  iado=ii+((ix*(ix+1))/2)+in
                  bc(iad)=bc(iad)+4d0*bc(iado)
                 end do
                else
                 do il=i2m,noc(isao)-1
                  iad=ioop+((il*(il+1))/2)+i2m
                  iado=ii+i1m+noc(isc)*il
                  bc(iad)=bc(iad)+4d0*bc(iado)
                 end do
                end if
                ii=ii+noo
               end do
               i10=1
              end do
              go to 4
             end if
             if(isblk(1,is).eq.isao.and.isblk(2,is).eq.isc.and.
     $          isblk(3,is).eq.isao)then                                3d29s16
              call ilimts(noc(isao),noc(isc),mynprocg,mynowprog,il,     3d29s16
     $           ih,i1s,i1e,i2s,i2e)
              ii=ioooo(is)
              if(isc.eq.isao)then
               noo=(noc(isao)*(noc(isao)+1))/2
              else
               noo=noc(isao)*noc(isc)
              end if
              i10=i1s
              i1n=noc(isao)                                             3d29s16
              do i2=i2s,i2e
               if(i2.eq.i2e)i1n=i1e
               i2m=i2-1
               do i1=i10,i1n
                i1m=i1-1
                if(isao.eq.isc)then
                 do il=i1m,noc(isao)-1
                  iad=ioop+((il*(il+1))/2)+i1m
                  ix=max(il,i2m)
                  in=min(il,i2m)
                  iado=ii+((ix*(ix+1))/2)+in
                  bc(iad)=bc(iad)+4d0*bc(iado)
                 end do
                else
                 do il=i1m,noc(isao)-1
                  iad=ioop+((il*(il+1))/2)+i1m
                  iado=ii+il+noc(isao)*i2m
                  bc(iad)=bc(iad)+4d0*bc(iado)
                 end do
                end if
                ii=ii+noo
               end do
               i10=1
              end do
              go to 4
             end if
             if(isblk(1,is).eq.isao.and.isblk(2,is).eq.isc.and.
     $          isblk(3,is).eq.isc)then                                 3d29s16
              call ilimts(noc(isc),noc(isao),mynprocg,mynowprog,il,     3d29s16
     $           ih,i1s,i1e,i2s,i2e)
              ii=ioooo(is)
              if(isc.eq.isao)then
               noo=(noc(isao)*(noc(isao)+1))/2
              else
               noo=noc(isao)*noc(isc)
              end if
              i10=i1s
              i1n=noc(isc)                                              3d29s16
              do i2=i2s,i2e
               if(i2.eq.i2e)i1n=i1e
               i2m=i2-1
               do i1=i10,i1n
                i1m=i1-1
                if(isao.eq.isc)then
                 do il=i2m,noc(isao)-1
                  iad=ioop+((il*(il+1))/2)+i2m
                  ix=max(il,i1m)
                  in=min(il,i1m)
                  iado=ii+((ix*(ix+1))/2)+in
                  bc(iad)=bc(iad)+4d0*bc(iado)
                 end do
                else
                 do il=i2m,noc(isao)-1
                  iad=ioop+((il*(il+1))/2)+i2m
                  iado=ii+il+noc(isao)*i2m
                  bc(iad)=bc(iad)+4d0*bc(iado)
                 end do
                end if
                ii=ii+noo
               end do
               i10=1
              end do
              go to 4
             end if
             if(isblk(1,is).eq.isc.and.isblk(2,is).eq.isao.and.          3d29s16
     $          isblk(3,is).eq.isao)then                                3d29s16
              call ilimts(noc(isao),noc(isc),mynprocg,mynowprog,il,     3d29s16
     $           ih,i1s,i1e,i2s,i2e)
              ii=ioooo(is)
              noo=noc(isao)*noc(isc)
              i10=i1s
              i1n=noc(isao)                                             3d29s16
              do i2=i2s,i2e
               if(i2.eq.i2e)i1n=i1e
               i2m=i2-1
               do i1=i10,i1n
                i1m=i1-1
                do il=i1m,noc(isao)-1
                 iad=ioop+((il*(il+1))/2)+i1m
                 iado=ii+i2m+il*noc(isc)
                 bc(iad)=bc(iad)+4d0*bc(iado)
                end do
                ii=ii+noo
               end do
               i10=1
              end do
              go to 4
             end if
            end do
         write(6,*)('in buildhess')
            write(6,*)('missing isc: '),isc,isao,isc,isao               3d29s16
            do is=1,nsdlk
             write(6,*)(isblk(j,is),j=1,4)
            end do
            call dws_sync
            call dws_finalize
            stop
    4       continue
           end if
          end do
          noo=(noc(isao)*(noc(isao)+1))/2
          if(idwsdeb.gt.0)then
          write(6,*)('oop after i1=i3 sum ')
          write(6,*)(bc(ioop+i),i=0,noo-1)
          end if
          call dws_gsumf(bc(ioop),noo)
          if(idwsdeb.gt.10)then
          write(6,*)('2e part of oop: ')
          write(6,*)(bc(ioop+i),i=0,noo-1)
          end if
          joop=ioop
          nh=noc(isao)+nvirtc(isao)
          do i3=0,noc(isao)-1
           do i4=0,i3
            iadh=ipt(isao)+i4+nh*i3+1
            bc(joop)=bc(joop)-4d0*h0mo(iadh)
            joop=joop+1
           end do
          end do
          if(idwsdeb.gt.10)then
          write(6,*)('full ioop ')
          write(6,*)(bc(ioop+i),i=0,noo-1)
          end if
          do is=1,nsdlk
           if(isblk(3,is).eq.isav.and.isblk(4,is).eq.isav)then          3d29s16
            isc=isblk(1,is)
            call ilimts(nvirtc(isav),nvirtc(isav),mynprocg,mynowprog,il,3d29s16
     $           ih,i1s,i1e,i2s,i2e)
            ii=jmats(is)-1
            i10=i1s
            i1n=nvirtc(isav)                                            3d29s16
            nrowj=(noc(isc)*(noc(isc)+1))/2
            do i2=i2s,i2e
             if(i2.eq.i2e)i1n=i1e
             do i1=i10,i1n
              if(i1.ge.i2)then
               jvv=ivv+((i1*(i1-1))/2)+i2-1
               do i34=1,noc(isc)
                iadj=ii+(i34*(i34+1))/2
                bc(jvv)=bc(jvv)+8d0*bc(iadj)
               end do
              end if
              ii=ii+nrowj
             end do
             i10=1
            end do
           end if
          end do
          if(idwsdeb.gt.10)then
          write(6,*)('vv after 34 sum ')
          write(6,*)(bc(ivv+i),i=0,nvv-1)
          end if
          do is=1,nsdlkk
           if(isblkk(3,is).eq.isav.and.isblkk(4,is).eq.isav)then        3d29s16
            isc=isblkk(1,is)
            call ilimts(nvirtc(isav),nvirtc(isav),mynprocg,mynowprog,il,3d29s16
     $           ih,i1s,i1e,i2s,i2e)
            ii=kmats(is)
            i10=i1s
            i1n=nvirtc(isav)                                            3d29s16
            nrowk=noc(isc)*noc(isc)
            do i2=i2s,i2e
             if(i2.eq.i2e)i1n=i1e
             do i1=i10,i1n
              if(i1.ge.i2)then
               jvv=ivv+((i1*(i1-1))/2)+i2-1
               do i34=0,noc(isc)-1
                iadk=ii+i34*(noc(isc)+1)
                bc(jvv)=bc(jvv)-4d0*bc(iadk)
               end do
              end if
              ii=ii+nrowk
             end do
             i10=1
            end do
           end if
          end do
          if(idwsdeb.gt.10)then
          write(6,*)('vv before gsum ')
          write(6,*)(bc(ivv+i),i=0,nvv-1)
          end if
          call dws_gsumf(bc(ivv),nvv)
          if(idwsdeb.gt.10)then
          write(6,*)('2e part of vv: ')
          write(6,*)(bc(ivv+i),i=0,nvv-1)
          write(6,*)('hessian w/o oo,vv&h0 ')
          call prntm2(bc(ihess(isa,isb)),nrowh,nhere,nrowh)
         do icol=0,nhere-1
          iad=ihess(isa,isb)+icol*nrowh
         end do
          end if
          nh=noc(isav)+nvirtc(isav)                                     3d29s16
          i10=i1sh
          i1n=nvirtc(isbv)
          ii=ihess(isa,isb)
          do i2=i2sh,i2eh
           if(i2.eq.i2eh)i1n=i1eh
           do i1=i10,i1n
            ix=max(i1,i2)
            in=min(i1,i2)
            jvv=ivv+((ix*(ix-1))/2)+in-1
            iadh0=ipt(isav)+i2+noc(isav)+nh*(i1+noc(isav)-1)
            h0add=4d0*h0mo(iadh0)
            do i34=0,noc(isao)-1
             iadh=ii+i34*(noc(isao)+1)
             bc(iadh)=bc(iadh)+bc(jvv)+h0add
            end do
            if(i1.eq.i2)then
             do i4=0,noc(isao)-1
              do i3=0,noc(isao)-1
               ix=max(i3,i4)
               in=min(i3,i4)
               joop=ioop+((ix*(ix+1))/2)+in
               bc(ii)=bc(ii)+bc(joop)
               ii=ii+1
              end do
             end do
            else
             ii=ii+nrowh
            end if
           end do
           i10=1
          end do
          ioffh=ioffh+nh*nh
         end if
         if(idwsdeb.gt.10)then
          if(nsymb.ne.1)then                                             4d24s18
           if(nhere.gt.0)then
            iunit=6
            write(iunit,*)('hessian for symmetry block '),isao,isbo,
     $           isbv,isav,isa,isb,ihess(isa,isb)
            call prntm2(bc(ihess(isa,isb)),nrowh,nhere,nrowh)
            iunit=6
           end if                                                        4d24s18
          else                                                           4d24s18
           call prntm2(bc(ihess(isa,isb)),nrowh,nhere,nrowh)
           itmb=ibcoff                                                   4d24s18
           ncol=nvirtc(isb)*nvirtc(isa)                                  4d24s18
           ibcoff=itmb+nrowh*ncol                                        4d24s18
           call enough('buildhess.  3',bc,ibc)
           do i=0,nrowh*ncol-1
            bc(itmb+i)=0d0
           end do                                                        4d24s18
           i10=i1sh                                                      4d24s18
           i1n=nvirtc(isb)                                               4d24s18
           ii=ihess(isa,isb)                                             4d24s18
           do i2=i2sh,i2eh
            i2m=i2-1                                                     4d24s18
            if(i2.eq.i2eh)i1n=i1eh                                       4d24s18
            do i1=i10,i1n                                                4d24s18
             i1m=i1-1                                                    4d24s18
             icol=i1m+nvirtc(isb)*i2m                                    4d24s18
             jj=itmb+nrowh*icol                                          4d24s18
             do i34=0,nrowh-1                                            4d24s18
              bc(jj+i34)=bc(ii+i34)                                      4d24s18
             end do                                                      4d24s18
             ii=ii+nrowh                                                 4d24s18
            end do                                                       4d24s18
            i10=1                                                        4d24s18
           end do
           call dws_gsumf(bc(itmb),nrowh*ncol)                           4d24s18
           call printa(bc(itmb),noc,0,noc,0,nvirtc,noc,nvirtc,noc,
     $         bc(ibcoff))
           ibcoff=itmb
          end if                                                         4d24s18
         end if
         ibcoff=ibtop                                                   3d29s16
        end if
       end do
      end do
      return
      end
