c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine buildhessdx(i4od,jmatd,kmatd,h0mod,noc,ihess,
     $     nvirtc,ipsym,multh,idwsdeb,isblkdd,nsdlkdd,isblkder,nsdlkder,11d14s22
     $     bc,ibc)                                                      11d14s22
      implicit real*8 (a-h,o-z)
c
c     build first derivative of orbital rotation hessian matrix wrt
c     a nuclear center. The orbital rotation hessian symmetries will be
c     sym bra o = sym bra v, sym ket o = sym ket v, and the symmetry of the
c     perturbation is ipsym. Thus the output matrix will have sym ket o *
c     sym ket v = ipsym.
c
      include "common.hf"
      dimension h0mod(1),i4od(1),jmatd(1),kmatd(1),noc(1),              12d5s16
     $     ihess(8,8),nvirtc(1),multh(8,8),ipt(8),isblkdd(4,idbk),      12d5s16
     $     isblkder(4,idbk)                                             12d5s16
      logical log(8)                                                    4d24s18
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      if(mynprocg.eq.1)then
       kgoal=1898561
      else
       kgoal=1783387
      end if
      igoal=1
      if(idwsdeb.gt.10)then
       write(6,*)('in buildhessd '),nvirtc(1),nsdlkdd,nsdlkder,h0mod(1)
       write(6,*)('ibcoff, nsymb '),ibcoff,nsymb
       write(6,*)('isblkdd ')
       do i=1,nsdlkdd
        write(6,*)(isblkdd(j,i),j=1,4)
       end do
       write(6,*)('isblkder ')
       do i=1,nsdlkder
        write(6,*)(isblkder(j,i),j=1,4)
       end do
      end if
      ioffh=0
      do isb=1,nsymb
       isk=multh(isb,ipsym)                                             12d12s16
       if(isk.ge.isb)then                                               12d12s16
        nhb=noc(isb)+nvirtc(isb)                                        12d12s16
        nhk=noc(isk)+nvirtc(isk)                                        12d12s16
        ipt(isb)=ioffh                                                  12d12s16
        ioffh=ioffh+nhb*nhk                                             12d12s16
       end if                                                           12d12s16
      end do
      do isk=1,nsymb
       do isb=1,nsymb
        if(idwsdeb.gt.10)write(6,*)('for isb,isk = '),isb,isk
c
c     the symmetry label isb,isk correspond to the bra and ket
c     occupied orbital symmetries. note that I've coded up things so the
c     order of occupied orbital indicies is bra,ket but the order of
c     virtual orbital indicies is ket,bra for the exchange matrices.
c     because of the distribution across processors, we need
c     isb,isk to be the virtual labels
c
        isbv=isb                                                        3d29s16
        iskv=isk                                                        3d29s16
        ihess(isbv,iskv)=-1                                               12d5s16
        isbo=isb                                                        12d5s16
        isko=multh(isk,ipsym)                                           3d28s16
        if(noc(isko)*noc(isbo).gt.0)then                                  3d28s16
         ihess(isbv,iskv)=ibcoff                                        12d5s16
         call ilimts(nvirtc(iskv),nvirtc(isbv),mynprocg,mynowprog,ilh,  12d5s16
     $        ihh,                                                      3d28s16
     $        i1sh,i1eh,i2sh,i2eh)                                      3d28s16
         nhere=ihh+1-ilh                                                3d28s16
         nrowh=noc(isbo)*noc(isko)                                      12d5s16
         ibcoff=ihess(isbv,iskv)+nhere*noc(isbo)*noc(isko)              12d5s16
         if(idwsdeb.gt.10)then
          write(6,*)('space for ihess: '),noc(isbo)*noc(isko),nhere
         end if
         ibtop=ibcoff                                                   3d29s16
         call enough('buildhessdx.  1',bc,ibc)
         do is=1,nsdlkdd
          if(isblkdd(1,is).eq.isbo.and.isblkdd(2,is).eq.isko.and.       12d5s16
     $         isblkdd(3,is).eq.iskv.and.isblkdd(4,is).eq.isbv)then     12d5s16
           if(idwsdeb.gt.10)write(6,*)('j block 1: '),
     $         (isblk(j,is),j=1,4)
           ii=jmatd(is)
           nrowj=noc(isko)*noc(isbo)
           kk=ihess(isbv,iskv)
           i10=i1sh
           i1n=nvirtc(iskv)                                             3d29s16
           do i2=i2sh,i2eh
            if(i2.eq.i2eh)i1n=i1eh
            do i1=i10,i1n
             do i4=0,noc(isko)-1
              do i3=0,noc(isbo)-1
               bc(kk)=-4d0*bc(ii)
               kk=kk+1
               ii=ii+1
              end do
             end do
            end do
            i10=1
           end do
           go to 1
          end if
         end do
         do is=1,nsdlkdd                                                12d9s16
          if(isblkdd(2,is).eq.isbo.and.isblkdd(1,is).eq.isko.and.       12d9s16
     $         isblkdd(3,is).eq.iskv.and.isblkdd(4,is).eq.isbv)then     12d9s16
           if(idwsdeb.gt.10)write(6,*)('jd block2 '),
     $         (isblkdd(j,is),j=1,4)
           ii=jmatd(is)                                                 12d9s16
           nrowj=noc(isko)*noc(isbo)                                    12d9s16
           kk=ihess(isbv,iskv)                                          12d9s16
           i10=i1sh                                                     12d9s16
           i1n=nvirtc(iskv)                                             12d9s16
           do i2=i2sh,i2eh                                              12d9s16
            if(i2.eq.i2eh)i1n=i1eh                                      12d9s16
            do i1=i10,i1n                                               12d9s16
             do i4=0,noc(isbo)-1                                        12d9s16
              do i3=0,noc(isko)-1                                       12d9s16
               iad=kk+i3*noc(isbo)+i4                                   12d9s16
               bc(iad)=-4d0*bc(ii)                                      12d9s16
               ii=ii+1                                                  12d9s16
              end do                                                    12d9s16
             end do                                                     12d9s16
             kk=kk+nrowj                                                12d9s16
            end do
            i10=1
           end do
           go to 1
          end if
         end do
         write(6,*)('in buildhessdx')
         write(6,*)('could not find jmatd for '),isbo,isko,iskv,isbv       3d29s16
         do is=1,nsdlkdd
          write(6,*)(isblkdd(j,is),j=1,4)
         end do
         call dws_sync
         call dws_finalize
         stop
    1    continue
         if(idwsdeb.gt.10)then
          write(6,*)('dhessian w jmat '),isbo,isko,iskv,isbv
          call prntm2(bc(ihess(isb,isk)),nrowh,nhere,nrowh)
          if(nsymb.eq.1)then
           write(6,*)('assemble full matrix to print as cs ')
           itmp=ibcoff
           nwds=nrowh*nvirtc(1)*nvirtc(1)
           ibcoff=itmp+nwds
           call enough('buildhessdx.  2',bc,ibc)
           do i=0,nwds-1
            bc(itmp+i)=0d0
           end do
           do i=0,nhere-1
            jtmp=itmp+nrowh*(ilh+i-1)
            jhess=ihess(isb,isk)+nrowh*i
            do j=0,nrowh-1
             bc(jtmp+j)=bc(jhess+j)
            end do
           end do
           call dws_gsumf(bc(itmp),nwds)
           ibcoff=itmp
          end if
         end if
         do is=1,nsdlkdd
          if(isblkdd(1,is).eq.isbo.and.isblkdd(2,is).eq.isko.and.        12d5s16
     $         isblkdd(3,is).eq.iskv.and.isblkdd(4,is).eq.isbv)then     12d5s16
           if(idwsdeb.gt.10)write(6,*)('kd block 1'),
     $         (isblkdd(j,is),j=1,4)
           kk=ihess(isbv,iskv)                                          12d5s16
           ii=kmatd(is)
           i10=i1sh
           i1n=nvirtc(iskv)                                             3d29s16
           do i2=i2sh,i2eh
            if(i2.eq.i2eh)i1n=i1eh
            do i1=i10,i1n
             do i4=0,noc(isko)-1                                        12d5s16
              do i3=0,noc(isbo)-1                                       12d5s16
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
         nwds=nrowh*nvirtc(isbv)*nvirtc(iskv)                           12d9s16
         ihessa=ibcoff
         ibcoff=ihessa+nwds                                             12d9s16
         call enough('buildhessdx.  3',bc,ibc)
         do i=0,nwds-1                                                  12d9s16
          bc(ihessa+i)=0d0                                              12d9s16
         end do                                                         12d9s16
         do is=1,nsdlkdd
          if(isblkdd(1,is).eq.isko.and.isblkdd(2,is).eq.isbo.and.       12d9s16
     $         isblkdd(4,is).eq.iskv.and.isblkdd(3,is).eq.isbv)then     12d5s16
          call ilimts(nvirtc(isbv),nvirtc(iskv),mynprocg,mynowprog,ilk,  12d5s16
     $        ihk,                                                      3d28s16
     $        i1sk,i1ek,i2sk,i2ek)                                      3d28s16
           if(idwsdeb.gt.10)write(6,*)('kd block 2'),
     $         (isblkdd(j,is),j=1,4),ilk,ihk
           ii=kmatd(is)
           i10=i1sk
           i1n=nvirtc(isbv)                                             3d29s16
           do i2=i2sk,i2ek
            if(i2.eq.i2ek)i1n=i1ek
            i2m=i2-1                                                    12d9s16
            do i1=i10,i1n
             i1m=i1-1                                                   12d9s16
             do i4=0,noc(isbo)-1                                        12d5s16
              do i3=0,noc(isko)-1                                       12d5s16
               iad=ihessa+i4+noc(isbo)*(i3+noc(isko)*                   12d9s16
     $              (i2m+nvirtc(iskv)*i1m))                             12d9s16
               bc(iad)=bc(iad)+16d0*bc(ii)                              12d9s16
               ii=ii+1
              end do
             end do
            end do
            i10=1
           end do
           call dws_gsumf(bc(ihessa),nwds)                              12d9s16
           i10=i1sh                                                     12d9s16
           i1n=nvirtc(iskv)                                             12d9s16
           ii=ihess(isbv,iskv)
           do i2=i2sh,i2eh                                              12d9s16
            i2m=i2-1                                                    12d9s16
            if(i2.eq.i2eh)i1n=i1eh                                      12d9s16
            do i1=i10,i1n
             i1m=i1-1                                                   12d9s16
             do i34=0,noc(isbo)*noc(isko)-1                             12d9s16
              iad=ihessa+i34+nrowh*(i1m+nvirtc(iskv)*i2m)               12d9s16
              bc(ii+i34)=bc(ii+i34)+bc(iad)                             12d9s16
             end do
             ii=ii+nrowh
            end do
            i10=1                                                       2d21s23
           end do
           ibcoff=ihessa                                                12d9s16
           go to 2
          end if
         end do
         write(6,*)('in buildhessdx')
         write(6,*)('could not find kmatd for '),isbo,isko,iskv,isbv
         do is=1,nsdlkdd
          write(6,*)(isblkdd(j,is),j=1,4)
         end do
         call dws_sync
         call dws_finalize
         stop
    2    continue
         if(idwsdeb.gt.10)then
          write(6,*)('dhessian w jmat and kmat '),isbo,isko,iskv,isbv
          call prntm2(bc(ihess(isb,isk)),nrowh,nhere,nrowh)
          if(nsymb.eq.1)then
           write(6,*)('assemble full matrix to print as cs ')
           itmp=ibcoff
           nwds=nrowh*nvirtc(1)*nvirtc(1)
           ibcoff=itmp+nwds
           call enough('buildhessdx.  4',bc,ibc)
           do i=0,nwds-1
            bc(itmp+i)=0d0
           end do
           do i=0,nhere-1
            jtmp=itmp+nrowh*(ilh+i-1)
            jhess=ihess(isb,isk)+nrowh*i
            do j=0,nrowh-1
             bc(jtmp+j)=bc(jhess+j)
            end do
           end do
           call dws_gsumf(bc(itmp),nwds)
           ibcoff=itmp
          end if
         end if
         do is=1,nsdlkdd
          if(isblkdd(1,is).eq.isko.and.isblkdd(2,is).eq.isbo.and.       12d5s16
     $         isblkdd(3,is).eq.iskv.and.isblkdd(4,is).eq.isbv)then     12d5s16
           if(idwsdeb.gt.10)write(6,*)('kd block 3'),
     $         (isblkdd(j,is),j=1,4)
           kk=ihess(isbv,iskv)
           ii=kmatd(is)
           i10=i1sh
           i1n=nvirtc(iskv)                                             3d29s16
           do i2=i2sh,i2eh
            if(i2.eq.i2eh)i1n=i1eh
            do i1=i10,i1n
             do i4=0,noc(isko)-1
              do i3=0,noc(isbo)-1
               iadk=ii+i4+noc(isko)*i3
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
          if(isblkdd(2,is).eq.isko.and.isblkdd(1,is).eq.isbo.and.       12d9s16
     $         isblkdd(4,is).eq.iskv.and.isblkdd(3,is).eq.isbv)then     12d12s16
           if(idwsdeb.gt.10)write(6,*)('kd block 4'),
     $         (isblkdd(j,is),j=1,4)
           nwds=nrowh*nvirtc(isbv)*nvirtc(iskv)                           12d9s16
           ihessa=ibcoff
           ibcoff=ihessa+nwds                                             12d9s16
           call enough('buildhessdx.  5',bc,ibc)
           do i=0,nwds-1                                                  12d9s16
            bc(ihessa+i)=0d0                                              12d9s16
           end do                                                         12d9s16
           call ilimts(nvirtc(isbv),nvirtc(iskv),mynprocg,mynowprog,ilk,  12d5s16
     $        ihk,                                                      3d28s16
     $        i1sk,i1ek,i2sk,i2ek)                                      3d28s16
           ii=kmatd(is)
           i10=i1sk                                                     12d17s16
           i1n=nvirtc(isbv)                                             12d12s16
           do i2=i2sk,i2ek                                              12d17s16
            i2m=i2-1                                                    12d12s16
            if(i2.eq.i2ek)i1n=i1ek                                      12d17s16
            do i1=i10,i1n
             i1m=i1-1                                                   12d12s16
             do i4=0,noc(isko)-1                                        12d9s16
              do i3=0,noc(isbo)-1                                       12d9s16
               iadk=ii+i3+noc(isbo)*i4                                  12d12s16
               iadh=ihessa+i3+noc(isbo)*                                2d1s17
     $              (i4+noc(isko)*(i2m+nvirtc(iskv)*i1m))               2d1s17
               bc(iadh)=bc(iadh)-4d0*bc(iadk)                           12d12s16
              end do
             end do
             ii=ii+nrowh
            end do
            i10=1
           end do
           call dws_gsumf(bc(ihessa),nwds)                              12d9s16
           i10=i1sh                                                     12d9s16
           i1n=nvirtc(iskv)                                             12d17s16
           ii=ihess(isbv,iskv)
           do i2=i2sh,i2eh                                              12d9s16
            i2m=i2-1                                                    12d9s16
            if(i2.eq.i2eh)i1n=i1eh                                      12d9s16
            do i1=i10,i1n
             i1m=i1-1                                                   12d9s16
             do i34=0,nrowh-1                                           12d17s16
              iad=ihessa+i34+nrowh*(i1m+nvirtc(iskv)*i2m)               2d1s17
              bc(ii+i34)=bc(ii+i34)+bc(iad)                             12d9s16
             end do
             ii=ii+nrowh
            end do
            i10=1                                                       2d1s17
           end do
           ibcoff=ihessa                                                12d9s16
           go to 3
          end if
         end do
         write(6,*)('in buildhessdx')
         write(6,*)('could not find kmatT for '),isbo,isko,iskv,isbv      3d29s16
         do is=1,nsdlkdd
          write(6,*)(isblkdd(j,is),j=1,4)
         end do
         call dws_sync
         call dws_finalize
         stop
    3    continue
         if(idwsdeb.gt.10)then
          write(6,*)('dhessian w/o oo,vv&h0 '),isbo,isko,iskv,isbv,
     $         ihess(isb,isk)
          call prntm2(bc(ihess(isb,isk)),nrowh,nhere,nrowh)
         end if
         if(isbo.eq.isko)then                                           12d9s16
          if(idwsdeb.gt.0)write(6,*)('going for vv ')
          nvv=nvirtc(isbv)*nvirtc(iskv)                                 12d9s16
          ivv=ibcoff                                                    12d9s16
          ibcoff=ivv+nvv                                                12d9s16
          do i=0,nvv-1                                                  12d9s16
           bc(ivv+i)=0d0                                                12d9s16
          end do                                                        12d9s16
          igot=0
          do is=1,nsdlkdd
           if(isblkdd(3,is).eq.isbv.and.isblkdd(4,is).eq.iskv)then          3d29s16
            isc1=isblkdd(1,is)                                          12d9s16
            isc2=isblkdd(2,is)                                          12d9s16
            call ilimts(nvirtc(isbv),nvirtc(iskv),mynprocg,mynowprog,il,3d29s16
     $           ih,i1s,i1e,i2s,i2e)
            ii=jmatd(is)                                                12d8s16
            i10=i1s
            i1n=nvirtc(isbv)                                            3d29s16
            nrowj=noc(isc1)*noc(isc2)                                   12d9s16
            do i2=i2s,i2e
             i2m=i2-1                                                   12d5s16
             if(i2.eq.i2e)i1n=i1e
             do i1=i10,i1n
              i1m=i1-1                                                  12d5s16
              jvv=ivv+i1m+nvirtc(isbv)*i2m                              12d5s16
              do i34=0,noc(isc1)-1                                       12d5s16
               iadj=ii+i34*(noc(isc1)+1)                                 12d5s16
               bc(jvv)=bc(jvv)+8d0*bc(iadj)
              end do
              ii=ii+nrowj
             end do
             i10=1
            end do
            igot=igot+1
           end if
          end do
          if(igot.eq.0)then
           write(6,*)('could not find jmatd for vv: '),isbv,iskv
           do i=1,nsdlkdd
            write(6,*)(isblkdd(j,i),j=1,4)
           end do
           call dws_sync
           call dws_finalize
           stop
          end if
          igot=0
          do is=1,nsdlkdd                                               12d5s16
           if(isblkdd(3,is).eq.isbv.and.isblkdd(4,is).eq.iskv)then      12d5s16
            isc=isblkdd(1,is)                                           12d5s16
            call ilimts(nvirtc(isbv),nvirtc(iskv),mynprocg,mynowprog,il,12d9s16
     $           ih,i1s,i1e,i2s,i2e)
            ii=kmatd(is)
            i10=i1s
            i1n=nvirtc(isbv)                                            3d29s16
            nrowk=noc(isc)*noc(isc)
            do i2=i2s,i2e
             i2m=i2-1                                                   12d5s16
             if(i2.eq.i2e)i1n=i1e
             do i1=i10,i1n
              i1m=i1-1                                                  12d5s16
              jvv=ivv+i1m+nvirtc(isb)*i2m                               12d5s16
              do i34=0,noc(isc)-1
               iadk=ii+i34*(noc(isc)+1)
               bc(jvv)=bc(jvv)-4d0*bc(iadk)
              end do
              ii=ii+nrowk
             end do
             i10=1
            end do
            igot=igot+1
           end if
          end do
          if(igot.eq.0)then
           write(6,*)('could not find kmatd for vv: '),isbv,iskv
           do i=1,nsdlkdd
            write(6,*)(isblkdd(j,i),j=1,4)
           end do
           call dws_sync
           call dws_finalize
           stop
          end if
          call dws_gsumf(bc(ivv),nvv)
          nh=noc(isbv)+nvirtc(isbv)                                     3d29s16
          nhk=noc(iskv)+nvirtc(iskv)                                    12d12s16
c
          i10=i1sh
          i1n=nvirtc(iskv)
          ii=ihess(isb,isk)
          do i2=i2sh,i2eh
           i2m=i2-1                                                     12d5s16
           if(i2.eq.i2eh)i1n=i1eh
           do i1=i10,i1n
            i1m=i1-1                                                    12d5s16
            jvv=ivv+i2m+nvirtc(isbv)*i1m                                12d12s16
            if(iskv.ge.isbv)then                                        12d12s16
             iadh0=ipt(isbv)+i2+noc(isbv)+nh*(i1+noc(iskv)-1)            12d9s16
            else                                                        12d12s16
             iadh0=ipt(iskv)+i1+noc(iskv)+nhk*(i2+noc(isbv)-1)          12d12s16
            end if                                                      12d12s16
            h0add=4d0*h0mod(iadh0)
            do i34=0,noc(isb)-1
             iadh=ii+i34*(noc(isb)+1)
             bc(iadh)=bc(iadh)+bc(jvv)+h0add
            end do
            ii=ii+nrowh
           end do
           i10=1
          end do
         end if
         if(isb.eq.isk)then
          if(idwsdeb.gt.10)write(6,*)('going for oo')
          noo=noc(isbo)*noc(isko)                                       12d5s16
          ioop=ibcoff
          ibcoff=ioop+noo
          call enough('buildhessdx.  7',bc,ibc)
          do i=0,noo-1
           bc(ioop+i)=0d0
          end do
          igot=0                                                        12d9s16
          igotn=ibcoff                                                  2d24s23
          ibcoff=igotn+nsymb                                            2d24s23
          call enough('buildhessdx.gotn',bc,ibc)                        2d24s23
          jgotn=igotn-1                                                 2d24s23
          do is=1,nsymb                                                 2d24s23
           ibc(jgotn+is)=0                                              2d24s23
          end do                                                        2d24s23
          do is=1,nsdlkder
           if(isblkder(1,is).eq.isbo.and.isblkder(2,is).eq.isko)then
            iscv=isblkder(3,is)
            if(ibc(jgotn+iscv).eq.0)then                                2d24s23
             ibc(jgotn+iscv)=1                                          2d24s23
             ii=i4od(is)
             do i2=1,noc(iscv)                                           12d8s16
              do i1=1,noc(iscv)                                          12d8s16
               if(i1.eq.i2)then
                do i34=0,noo-1
                 bc(ioop+i34)=bc(ioop+i34)-8d0*bc(ii+i34)
                end do
               end if
               ii=ii+noo
              end do
             end do
             igot=igot+1                                                 12d9s16
            end if                                                      2d24s23
           else if(isblkder(2,is).eq.isbo.and.                          12d9s16
     $           isblkder(1,is).eq.isko)then                            12d9s16
            iscv=isblkder(3,is)                                         12d9s16
            if(ibc(jgotn+iscv).eq.0)then                                2d24s23
             ibc(jgotn+iscv)=1                                          2d24s23
             ii=i4od(is)                                                 12d9s16
             do i2=1,noc(iscv)                                           12d8s16
              do i1=1,noc(iscv)                                          12d8s16
               if(i1.eq.i2)then                                          12d9s16
                do i4=0,noc(isbo)-1                                      12d9s16
                 do i3=0,noc(isko)-1                                     12d9s16
                  iad=ioop+i4+noc(isbo)*i3                               12d9s16
                  bc(iad)=bc(iad)-8d0*bc(ii)                             12d9s16
                  ii=ii+1                                                12d9s16
                 end do                                                  12d9s16
                end do                                                   12d9s16
               else
                ii=ii+noc(isbo)*noc(isko)                                12d9s16
               end if                                                    12d9s16
              end do                                                     12d9s16
             end do                                                      12d9s16
             igot=igot+1                                                 12d9s16
            end if                                                      2d24s23
           else if(isblkder(4,is).eq.isbo.and.                          12d9s16
     $           isblkder(3,is).eq.isko)then                            12d9s16
            iscv=isblkder(1,is)                                         12d9s16
            if(ibc(jgotn+iscv).eq.0)then                                2d24s23
             ibc(jgotn+iscv)=1                                          2d24s23
             do i2=0,noc(isbo)-1                                         12d9s16
              do i1=0,noc(isko)-1                                        12d9s16
               iad=ioop+i2+noc(isbo)*i1                                  12d9s16
               do i34=0,noc(iscv)-1                                      12d9s16
                iad2=i4od(is)+i34+noc(iscv)*(i34+noc(iscv)               12d9s16
     $              *(i1+noc(isko)*i2))                                 12d9s16
                bc(iad)=bc(iad)-8d0*bc(iad2)                             12d9s16
               end do                                                    12d9s16
              end do                                                     12d9s16
             end do                                                      12d9s16
             igot=igot+1                                                 12d9s16
            end if                                                      2d24s23
           else if(isblkder(3,is).eq.isbo.and.                          12d9s16
     $           isblkder(4,is).eq.isko)then                            12d9s16
            iscv=isblkder(1,is)                                         12d9s16
            if(ibc(jgotn+iscv).eq.0)then                                2d24s23
             ibc(jgotn+iscv)=1                                          2d24s23
             do i2=0,noc(isko)-1                                         12d9s16
              do i1=0,noc(isbo)-1                                        12d9s16
               iad=ioop+i1+noc(isbo)*i2                                  12d9s16
               do i34=0,noc(iscv)-1                                      12d9s16
                iad2=i4od(is)+i34+noc(iscv)*(i34+noc(iscv)               12d9s16
     $              *(i1+noc(isbo)*i2))                                 12d9s16
                bc(iad)=bc(iad)-8d0*bc(iad2)                             12d9s16
               end do                                                    12d9s16
              end do                                                     12d9s16
             end do                                                      12d9s16
             igot=igot+1                                                 12d9s16
            end if                                                      2d24s23
           end if
          end do
          ibcoff=igotn                                                  2d24s23
          if(igot.gt.0)go to 5                                          12d9s16
          write(6,*)('did not find 4od for '),isbo,isko
          do i=1,nsdlkder
           write(6,*)(isblkder(j,i),j=1,4)
          end do
          call dws_sync
          call dws_finalize
          stop
    5     continue                                                      12d9s16
          igot=0                                                        12d9s16
          do isc=1,nsymb
           if(noc(isc).gt.0)then
            jgot=0
            do is=1,nsdlkder
             if(isblkder(1,is).eq.isc.and.isblkder(2,is).eq.isko.and.   12d5s16
     $          isblkder(3,is).eq.isc)then                                 3d29s16
              ii=i4od(is)
              noo=noc(isko)*noc(isc)                                    12d5s16
              do i2=0,noc(isbo)-1                                       12d8s16
               do i1=0,noc(isc)-1                                       12d8s16
                do il=0,noc(isko)-1
                 iad=ioop+i2+noc(isbo)*il
                 iado=ii+i1+noc(isc)*il
c     iado=ii+i1+noc(isc)*(il+noc(isko)*(i1+noc(isc)*i2))
                 bc(iad)=bc(iad)+4d0*bc(iado)
                end do
                ii=ii+noo
               end do
              end do
              igot=igot+1
              jgot=jgot+1                                               2d24s23
              go to 4                                                   2d24s23
             else if(isblkder(3,is).eq.isbo.and.isblkder(4,is).eq.isc.  12d9s16
     $             and.isblkder(1,is).eq.isko)then                      12d9s16
              ii=i4od(is)
              noo=noc(isko)*noc(isc)                                    12d5s16
              do i2=0,noc(isc)-1                                        12d9s16
               do i1=0,noc(isbo)-1                                       12d8s16
                do il=0,noc(isko)-1
                 iad=ioop+i1+noc(isbo)*il                               12d9s16
                 iado=ii+il+noc(isko)*i2                                 12d9s16
c     iado=ii+il+noc(isko)*(i2+noc(isc)*(i1+noc(isbo)*i2))
                 bc(iad)=bc(iad)+4d0*bc(iado)
                end do
                ii=ii+noo
               end do
              end do
              igot=igot+1                                               12d9s16
              jgot=jgot+1                                               2d24s23
              go to 4                                                   2d24s23
             else if(isblkder(2,is).eq.isc.and.isblkder(1,is).eq.isbo.  12d9s16
     $             and.isblkder(3,is).eq.isc)then                       12d9s16
              ii=i4od(is)
              noo=noc(isbo)*noc(isc)                                    12d5s16
              do i2=0,noc(isko)-1                                       12d8s16
               do i1=0,noc(isc)-1                                       12d8s16
                do il=0,noc(isbo)-1
                 iad=ioop+il+noc(isbo)*i2                               12d9s16
c     iado=ii+il+noc(isbo)*(i1+noc(isc)*(i1+noc(isc)*i2)
                 iado=ii+il+noc(isbo)*i1                                12d9s16
                 bc(iad)=bc(iad)+4d0*bc(iado)
                end do
                ii=ii+noo
               end do
              end do
              igot=igot+1
              jgot=jgot+1                                               2d24s23
              go to 4                                                   2d24s23
             else if(isblkder(2,is).eq.isc.and.isblkder(1,is).eq.isbo.  4d24s18
     $             and.isblkder(4,is).eq.isc)then                       4d24s18
              ii=i4od(is)                                               4d24s18
              noo=noc(isbo)*noc(isc)                                    12d5s16
              do i1=0,noc(isc)-1                                        4d24s18
               do i2=0,noc(isko)-1                                      4d24s18
                do il=0,noc(isbo)-1
                 iad=ioop+il+noc(isbo)*i2                               12d9s16
c     iado=ii+il+noc(isbo)*(i1+noc(isc)*(i2+noc(isko)*i1)
                 iado=ii+il+noc(isbo)*i1                                12d9s16
                 bc(iad)=bc(iad)+4d0*bc(iado)
                end do
                ii=ii+noo
               end do
              end do
              igot=igot+1
              jgot=jgot+1                                               2d24s23
              go to 4                                                   2d24s23
             else if(isblkder(2,is).eq.isc.and.isblkder(4,is).eq.isbo.  5d4s22
     $             and.isblkder(3,is).eq.isc)then                       5d4s22
              ii=i4od(is)                                               4d24s18
              noo=noc(isko)*noc(isc)                                    5d4s22
              do il=0,noc(isbo)-1                                       5d4s22
               do i1=0,noc(isc)-1                                        4d24s18
                do i2=0,noc(isko)-1                                      4d24s18
                 iad=ioop+il+noc(isbo)*i2                               12d9s16
c     iado=ii+i2+noc(isko)*(i1+noc(isc)*(i1+noc(isc)*il)
                 iado=ii+i2+noc(isko)*i1                                5d4s22
                 bc(iad)=bc(iad)+4d0*bc(iado)
                end do
                ii=ii+noo
               end do
              end do
              igot=igot+1
              jgot=jgot+1                                               2d24s23
              go to 4                                                   2d24s23
             end if
            end do
            if(jgot.gt.0)go to 4                                        2d24s23
            write(6,*)('in buildhessdx')
            write(6,*)('missing isc: '),isc,isko,isc,isb                12d5s16
            write(6,*)('isbo, isko: '),isbo,isko
            do is=1,nsdlkder
            log(1)=isblkder(1,is).eq.isc.and.isblkder(2,is).eq.isbo
     $            .and.isblkder(3,is).eq.isc
            log(2)=isblkder(2,is).eq.isc.and.isblkder(1,is).eq.isbo
     $            .and.isblkder(3,is).eq.isc
            log(3)=isblkder(1,is).eq.isc.and.isblkder(2,is).eq.isbo
     $            .and.isblkder(3,is).eq.isko
            log(4)=isblkder(2,is).eq.isc.and.isblkder(1,is).eq.isbo
     $            .and.isblkder(3,is).eq.isko
            log(5)=isblkder(3,is).eq.isc.and.isblkder(4,is).eq.isbo
     $            .and.isblkder(1,is).eq.isc
            log(6)=isblkder(4,is).eq.isc.and.isblkder(3,is).eq.isbo
     $            .and.isblkder(1,is).eq.isc
            log(7)=isblkder(3,is).eq.isc.and.isblkder(4,is).eq.isbo
     $            .and.isblkder(2,is).eq.isko
            log(8)=isblkder(4,is).eq.isc.and.isblkder(3,is).eq.isbo
     $            .and.isblkder(2,is).eq.isko
             if(log(1).or.log(2).or.log(3).or.log(4).or.log(5).or.log(6)
     $           .or.log(7).or.log(8))then
              write(6,*)('got it with '),log
             end if
             write(6,*)(isblkder(j,is),j=1,4)
            end do
            call dws_sync
            call dws_finalize
            stop
    4       continue
           end if
          end do
          noo=noc(isbo)*noc(isko)                                       12d5s16
          joop=ioop
          if(isko.ge.isbo)then                                          12d12s16
           nh=noc(isbo)+nvirtc(isbo)
           do i3=0,noc(isko)-1
            do i4=0,noc(isbo)-1
             iadh=ipt(isbo)+i4+nh*i3+1
             bc(joop)=bc(joop)-4d0*h0mod(iadh)
             joop=joop+1
            end do
           end do
          else                                                          12d12s16
           nh=noc(isko)+nvirtc(isko)                                    12d12s16
           do i3=0,noc(isbo)-1                                          12d12s16
            do i4=0,noc(isko)-1                                         12d12s16
             iadh=ipt(isko)+i4+nh*i3+1                                  12d12s16
             bc(joop)=bc(joop)-4d0*h0mod(iadh)                          12d12s16
             joop=joop+1                                                12d12s16
            end do                                                      12d12s16
           end do                                                       12d12s16
          end if                                                        12d12s16
c
          i10=i1sh
          i1n=nvirtc(iskv)
          ii=ihess(isb,isk)
          do i2=i2sh,i2eh
           i2m=i2-1                                                     12d5s16
           if(i2.eq.i2eh)i1n=i1eh
           do i1=i10,i1n
            i1m=i1-1                                                    12d5s16
            if(i1.eq.i2)then
             do i4=0,noc(isko)-1                                        12d12s16
              do i3=0,noc(isbo)-1
               joop=ioop+i3+noc(isbo)*i4                                12d5s16
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
          write(6,*)('dhessian for symmetry block '),isbo,isko,iskv,isbv
          write(6,*)('address in bc: '),ihess(isb,isk)
     $         ,ihess(isb,isk)+nrowh*nhere-ibtop
          write(6,*)('loc '),loc(ihess)
          call prntm2(bc(ihess(isb,isk)),nrowh,nhere,nrowh)
          if(nsymb.eq.1)call printa(bc(ihess(isb,isk)),noc,0,noc,0,
     $         nvirtc,noc,nvirtc,noc,bc(ibcoff))
         end if
         ibcoff=ibtop                                                   3d29s16
        end if
       end do
      end do
      return
      end
