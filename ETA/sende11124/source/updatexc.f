c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine updatexc(vdinout,gdin,hdd,vdold,ndoub,nhsdiag,         1d24s21
     $     vsold,hss,nff1,ncsf,nsymb,mdon,mdoo,nvirt,isymmrci,eig,      1d24s21
     $     nroot,isto,nfdat,vstry,nsing,lprt,ngots,ngotd,               4d21s21
     $     iamconverged,pairs,ldavid,bc,ibc)                            11d10s22
      implicit real*8 (a-h,o-z)                                         6d22s18
c                                                                       6d22s18
c     davidson update for contracted doubles.                           11d10s20
c
      logical lprint,lprt,ldavid                                        12d28s21
      integer*8 nhsdiag(mdoo+1,nsymb,2),irow,icol,il8,ih8,i1,           4d21s21
     $     iamconverged(*)                                              4d21s21
      dimension vdinout(*),gdin(*),hdd(*),vdold(*),vsold(*),hss(*),     11d10s20
     $     nff1(mdoo+1,nsymb),ncsf(*),nvirt(*),eig(*),                  11d10s20
     $     nfdat(5,4,*),vstry(*),pairs(4,*)                             8d19s21
      include "common.basis"                                            7d10s19
      include "common.input"                                            7d10s19
      include "common.store"
      include "common.print"                                            1d13s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data icall/0/
      save
      i1=1                                                              11d10s20
      icall=icall+1
      ibcoffo=ibcoff
      if(iprtr(11).ne.0)then                                            1d13s20
       lprint=.true.                                                     6d26s18
      else                                                              1d13s20
       lprint=.false.                                                    6d26s18
      end if                                                            1d13s20
      if(lprint)then                                                    12d29s21
       write(6,*)('in updatexc ')
       if(ndoub.gt.0)then                                               2d17s21
        write(6,*)('input doubles vector ')
        call prntm2(vdinout,ndoub*nroot,1,ndoub*nroot)                   11d10s20
        write(6,*)('input hcd ')
        call prntm2(gdin,ndoub*nroot,1,ndoub*nroot)                      11d10s20
       end if                                                           2d17s21
       write(6,*)('input eigenvalues ')
       call prntm2(eig,1,nroot,1)
       ivstest=ibcoff                                                   1d24s21
       igstest=ivstest+nsing*nroot                                      1d24s21
       ibcoff=igstest+nsing*nroot                                       1d24s21
       call enough('updatexc.  1',bc,ibc)
       if(nsing.gt.0)then                                               3d19s21
        jvstest=ivstest                                                  1d24s21
        jgstest=igstest                                                  1d24s21
        do isb=1,nsymb                                                   11d10s20
         isbv=multh(isb,isymmrci)                                        11d10s20
         irow=nvirt(isbv)*nroot                                          11d10s20
         irow4=irow                                                      11d10s20
         do ii=mdon+1,mdoo+1                                             11d10s20
          iarg=ii-mdon                                                   11d10s20
          ndblock=ncsf(iarg)*nff1(ii,isb)
          if(ndblock.gt.0)then                                           11d10s20
           icol=ndblock                                                  11d10s20
           icol4=icol                                                     11d10s20
           call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)  11d10s20
           nhere=ih+1-il                                                  11d10s20
           il8=il                                                        11d10s20
           ih8=ih                                                        11d10s20
           ivec=ibcoff                                                   11d10s20
           ibcoff=ivec+irow4*nhere                                       11d10s20
           call ddi_get(bc,ibc,nhsdiag(ii,isb,1),i1,irow,il8,ih8,       11d15s22
     $          bc(ivec))                                               11d15s22
           write(6,*)('sym, nclop '),isb,ii,nff1(ii,isb),iarg,ncsf(iarg)
           write(6,*)('vec ')
           call prntm2(bc(ivec),irow4,nhere,irow4)                       11d10s20
           jvec=ivec                                                     1d24s21
           do if1=1,nff1(ii,isb)                                         1d24s21
            do i=0,ncsf(iarg)-1                                          1d24s21
             do iv=0,nvirt(isbv)-1                                       1d24s21
              do ir=0,nroot-1                                            1d24s21
               iadg=jvstest+nsing*ir                                     1d24s21
               bc(iadg)=bc(jvec+ir)                                      1d24s21
              end do                                                     1d24s21
              jvstest=jvstest+1                                          1d24s21
              jvec=jvec+nroot                                            1d24s21
             end do                                                      1d24s21
            end do                                                       1d24s21
           end do                                                        1d24s21
           write(6,*)('hcs')
           call ddi_get(bc,ibc,nhsdiag(ii,isb,2),i1,irow,il8,ih8,       11d15s22
     $          bc(ivec))                                               11d15s22
           jvec=ivec                                                     1d24s21
           do if1=1,nff1(ii,isb)                                         1d24s21
            do i=0,ncsf(iarg)-1                                          1d24s21
             do iv=0,nvirt(isbv)-1                                       1d24s21
              do ir=0,nroot-1                                            1d24s21
               iadg=jgstest+nsing*ir                                     1d24s21
               bc(iadg)=bc(jvec+ir)                                      1d24s21
              end do                                                     1d24s21
              jgstest=jgstest+1                                          1d24s21
              jvec=jvec+nroot                                            1d24s21
             end do                                                      1d24s21
            end do                                                       1d24s21
           end do                                                        1d24s21
           call prntm2(bc(ivec),irow4,nhere,irow4)                       11d10s20
           ibcoff=ivec                                                   11d10s20
          end if                                                         11d10s20
         end do                                                          11d10s20
        end do                                                           11d10s20
        write(6,*)('vector for mpec2.0')
        call prntm2(bc(ivstest),nsing,nroot,nsing)
        call prntm2(bc(igstest),nsing,nroot,nsing)
       end if                                                           3d19s21
      end if                                                            6d26s18
      do ir=1,nroot                                                     8d19s21
       do l=1,4                                                          8d19s21
        pairs(l,ir)=0d0                                                 8d19s21
       end do                                                           8d19s21
      end do                                                            8d19s21
      if(ndoub.gt.0)then                                                2d17s21
       ioff=1
       ioffd=1                                                           12d11s20
       do jsb=1,nsymb                                                    12d8s20
        jsbv12=multh(jsb,isymmrci)                                       12d8s20
        do jsbv1=1,nsymb                                                 12d8s20
         jsbv2=multh(jsbv1,jsbv12)
         if(jsbv2.ge.jsbv1)then
          if(jsbv12.eq.1.and.nvirt(jsbv1).gt.0)then
           nn=nfdat(3,1,isymmrci)*nroot
           joff=ioff
           do if=1,nfdat(3,1,isymmrci)                                   11d10s20
            do ir=1,nroot
             joffd=ioffd                                                 12d12s20
             if(icall.gt.1.and.ngotd.eq.0)then                          12d29s21
              do iv=0,nvirt(jsbv1)-1                                    12d29s21
               pairs(1,ir)=pairs(1,ir)+gdin(joff)*vdinout(joff)         12d29s21
               gdin(joff)=0d0                                           12d29s21
               joff=joff+1                                              12d29s21
               joffd=joffd+1                                            12d29s21
              end do                                                    12d29s21
             else                                                       12d28s21
              do iv=0,nvirt(jsbv1)-1
               resid=gdin(joff)-eig(ir)*vdinout(joff)                     11d10s20
               pairs(1,ir)=pairs(1,ir)+gdin(joff)*vdinout(joff)          8d19s21
               bot=hdd(joffd)-eig(ir)                                     12d12s20
               bot=bot+1d-14                                              11d10s20
               update=-resid/bot                                               6d22s18
               vdinout(joff)=update                                       11d10s20
               gdin(joff)=0d0                                             11d10s20
               joff=joff+1
               joffd=joffd+1                                              12d12s20
              end do
             end if                                                     12d29s21
            end do
            ioffd=joffd                                                  12d12s20
           end do
           ioff=ioff+nvirt(jsbv1)*nn
          end if
          if(jsbv12.eq.1)then
           nvvs=(nvirt(jsbv1)*(nvirt(jsbv1)+1))/2                         12d8s20
           nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                         12d8s20
           nvvt=nvv
           isw=0
          else                                                           12d8s20
           nvv=nvirt(jsbv1)*nvirt(jsbv2)                                 12d8s20
           nvvs=nvv
           nvvt=nvv
           isw=1
          end if                                                         12d8s20
          if(nvvs.gt.0)then
           do l=1,4
            if(nfdat(3,l,jsb).gt.0)then
             nn=nfdat(3,l,jsb)*nroot
             joff=ioff
             do if=1,nfdat(3,l,jsb)                                      11d10s20
              do ir=1,nroot                                              11d10s20
               joffd=ioffd                                               12d12s20
               if(icall.gt.1.and.ngotd.eq.0)then                        12d29s21
                do iv12=1,nvv                                           12d29s21
                 pairs(l,ir)=pairs(l,ir)+gdin(joff)*vdinout(joff)       12d29s21
                 gdin(joff)=0d0                                         12d29s21
                 joff=joff+1                                            12d29s21
                 joffd=joffd+1                                          12d29s21
                end do                                                  12d29s21
               else                                                     12d29s21
                do iv12=1,nvv                                             11d10s20
                 resid=gdin(joff)-eig(ir)*vdinout(joff)                     11d10s20
                 pairs(l,ir)=pairs(l,ir)+gdin(joff)*vdinout(joff)        8d19s21
                 bot=hdd(joffd)-eig(ir)                                   12d12s20
                 bot=bot+1d-14                                              11d10s20
                 update=-resid/bot                                               6d22s18
                 vdinout(joff)=update                                       11d10s20
                 gdin(joff)=0d0                                             11d10s20
                 joff=joff+1
                 joffd=joffd+1                                            12d12s20
                end do
               end if
              end do
              ioffd=joffd                                                12d12s20
             end do                                                      12d9s20
             ioff=ioff+nn*nvv                                            12d8s20
            end if
           end do
          end if                                                         12d8s20
         end if
        end do                                                           12d8s20
       end do                                                            12d8s20
      end if                                                            2d17s21
      jvsold=1                                                          11d10s20
      jhss=1                                                            11d10s20
      if(nsing.gt.0)then                                                3d19s21
       do isb=1,nsymb                                                    11d10s20
        isbv=multh(isb,isymmrci)                                         11d10s20
        irow=nvirt(isbv)*nroot                                           11d10s20
        irow4=irow                                                       11d10s20
        do ii=mdon+1,mdoo+1                                              11d10s20
         iarg=ii-mdon                                                    11d10s20
         ndblock=ncsf(iarg)*nff1(ii,isb)
         if(ndblock.gt.0)then                                            11d10s20
          icol=ndblock                                                   11d10s20
          icol4=icol                                                      11d10s20
          call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)   11d10s20
          nhere=ih+1-il                                                   11d10s20
          il8=il                                                         11d10s20
          ih8=ih                                                         11d10s20
          ivecs=ibcoff                                                    11d10s20
          igs=ivecs+irow4*nhere                                          11d10s20
          ibcoff=igs+nvirt(isbv)*nhere                                   11d10s20
          call ddi_get(bc,ibc,nhsdiag(ii,isb,1),i1,irow,il8,ih8,        11d15s22
     $         bc(ivecs))                                               11d15s22
          call ddi_get(bc,ibc,nhsdiag(ii,isb,2),i1,irow,il8,ih8,bc(igs))11d15s22
          ic=0                                                           11d10s20
          jvecs=ivecs                                                    11d10s20
          jgs=igs                                                        11d10s20
          jvsold0=jvsold
          do iff=1,nff1(ii,isb)                                          11d10s20
           khss=jhss                                                     11d10s20
           do m=1,ncsf(iarg)                                             11d10s20
            ic=ic+1                                                      11d10s20
            if(ic.ge.il.and.ic.le.ih)then                                11d10s20
             khss=jhss                                                   12d14s20
             if(icall.gt.1.and.ngots.eq.0)then                          12d29s21
              do iv=1,nvirt(isbv)                                         11d10s20
               do ir=1,nroot                                                 11d10s20
                vstry(jvsold)=bc(jvecs)                                 12d29s21
                jgs=jgs+1                                                 11d10s20
                jvecs=jvecs+1                                             11d10s20
                jvsold=jvsold+1                                           11d10s20
               end do                                                     11d10s20
               khss=khss+1                                                12d14s20
              end do                                                       11d10s20
             else                                                       12d29s21
              do iv=1,nvirt(isbv)                                         11d10s20
               do ir=1,nroot                                                 11d10s20
                resid=bc(jgs)-eig(ir)*bc(jvecs)                           11d10s20
                bot=hss(khss)-eig(ir)                                     11d10s20
                bot=bot+1d-14                                              11d10s20
                update=-resid/bot                                         11d10s20
                vstry(jvsold)=update                                      11d10s20
                jgs=jgs+1                                                 11d10s20
                jvecs=jvecs+1                                             11d10s20
                jvsold=jvsold+1                                           11d10s20
               end do                                                     11d10s20
               khss=khss+1                                                12d14s20
              end do                                                       11d10s20
             end if
            end if                                                       11d10s20
           end do                                                        11d10s20
           jhss=khss                                                     11d10s20
          end do                                                         11d10s20
          ibcoff=ivecs                                                   11d10s20
         end if                                                          11d10s20
        end do                                                           11d10s20
       end do                                                            11d10s20
      end if                                                            3d19s21
c
c     print out what we have so far
c
      if(lprt)then                                                      1d25s21
       ntotq=nsing+ndoub
       ivtrialq=ibcoff
       ibcoff=ivtrialq+ntotq*nroot
       call enough('updatexc.  2',bc,ibc)
       do i=ivtrialq,ibcoff-1                                           2d19s21
        bc(i)=0d0                                                       2d19s21
       end do                                                           2d19s21
       jvv=ivtrialq
       jvs=1
       if(nsing.gt.0)then                                               3d19s21
        do isb=1,nsymb
         isbv=multh(isb,isymmrci)
         do nclo=mdon,mdoo
          nclop=nclo+1
          if(nff1(nclop,isb).gt.0)then
           jarg=nclop-mdon
           call ilimts(ncsf(jarg),nff1(nclop,isb),mynprocg,mynowprog,il, 2d19s21
     $         ih,i1s,i1e,i2s,i2e)                                      2d19s21
           nhere=ih+1-il                                                 2d19s21
           do i=0,nhere-1                                                2d19s21
            kvv=jvv+nvirt(isbv)*(il+i-1)                                 2d19s21
            do iv=0,nvirt(isbv)-1
             do ir=0,nroot-1
              bc(kvv+ntotq*ir)=vstry(jvs)
              jvs=jvs+1
             end do
             kvv=kvv+1
            end do
           end do
           jvv=jvv+ncsf(jarg)*nff1(nclop,isb)*nvirt(isbv)                2d19s21
          end if
         end do
        end do
       end if                                                           3d19s21
       if(ndoub.gt.0.and.mynowprog.eq.0)then                            2d19s21
        ioff=1
        do isb=1,nsymb
         isbv12=multh(isb,isymmrci)
         do isbv1=1,nsymb
          isbv2=multh(isbv1,isbv12)
          if(isbv1.le.isbv2)then
           if(isbv1.eq.isbv2)then
            nvvs=(nvirt(isbv1)*(nvirt(isbv1)+1))/2
            nvvt=(nvirt(isbv1)*(nvirt(isbv1)-1))/2
            do if=0,nfdat(3,1,isb)-1
             do ir=0,nroot-1
              do iv=0,nvirt(isbv1)-1
               ivv=((iv*(iv+1))/2)+iv
               ito=jvv+if+nfdat(3,1,isb)*ivv+ntotq*ir
               bc(ito)=vdinout(ioff)
               ioff=ioff+1
              end do
             end do
            end do
            isw=0
           else
            nvvs=nvirt(isbv1)*nvirt(isbv2)
            nvvt=nvvs
            isw=1
           end if
           nvv=nvvs
           ip=+1
           do l=1,4
            if(nfdat(3,l,isb).gt.0)then
             do if=0,nfdat(3,l,isb)-1
              do ir=0,nroot-1
               do iv2=0,nvirt(isbv2)-1
                if(isbv1.eq.isbv2)then
                 itop=iv2-1
                else
                 itop=nvirt(isbv1)-1
                end if
                do iv1=0,itop
                 itri=((iv2*(iv2+ip))/2)+iv1
                 irec=iv1+nvirt(isbv1)*iv2
                 itri=itri+isw*(irec-itri)
                 ito=jvv+if+nfdat(3,l,isb)*itri+ntotq*ir
                 bc(ito)=vdinout(ioff)
                 ioff=ioff+1
                end do
               end do
              end do
             end do
             jvv=jvv+nfdat(3,l,isb)*nvv
            end if
            ip=-1
            nvv=nvvt
           end do
          end if
         end do
        end do
       end if                                                           2d17s21
       write(6,*)('my part of trial vector ')
       call prntm2(bc(ivtrialq),ntotq,nroot,ntotq)
       write(6,*)('trial vector ')
       call dws_gsumf(bc(ivtrialq),ntotq*nroot)                         2d19s21
       call prntm2(bc(ivtrialq),ntotq,nroot,ntotq)
       ibcoff=ivtrialq
      end if                                                            1d25s21
c
c     modified gram-schmidtz orthogonalization of previous (old) vectors
c
      ngot=max(ngots,ngotd)                                             4d14s21
      do ir=1,nroot                                                     11d10s20
       if(lprt)write(6,*)('for new root no.: '),ir
       do iro=1,ngot                                                    1d27s21
        if(lprt)write(6,*)('for old root no.: '),iro
        sum2=0d0                                                        11d10s20
        if(ndoub.gt.0.and.iro.le.ngotd)then                             4d14s21
         ioff=0                                                          11d10s20
         ioffo=0                                                         11d10s20
         do jsb=1,nsymb                                                    12d8s20
          jsbv12=multh(jsb,isymmrci)                                       12d8s20
          do jsbv1=1,nsymb                                                 12d8s20
           jsbv2=multh(jsbv1,jsbv12)
           if(jsbv2.ge.jsbv1)then
            if(jsbv12.eq.1.and.nvirt(jsbv1).gt.0)then
             nn=nfdat(3,1,isymmrci)*nroot
             if(ldavid)then                                             12d28s21
              nno=nfdat(3,1,isymmrci)*ngotd                             12d28s21
              nrooto=ngotd                                              12d28s21
             else                                                       12d28s21
              nno=nn                                                    12d28s21
              nrooto=nroot                                              12d28s21
             end if                                                     12d28s21
             joff=ioff
             joffo=ioffo                                                 11d10s20
             do if=1,nfdat(3,1,isymmrci)                                   11d10s20
              joff=joff+nvirt(jsbv1)*(ir-1)                              11d10s20
              joffo=joffo+nvirt(jsbv1)*(iro-1)                              11d10s20
              do iv=1,nvirt(jsbv1)                                       11d10s20
               sum2=sum2+vdold(joffo+iv)*vdinout(joff+iv)                11d10s20
              end do                                                     11d10s20
              joffo=joffo+nvirt(jsbv1)*(nrooto+1-iro)                   12d28s21
              joff=joff+nvirt(jsbv1)*(nroot+1-ir)                        11d10s20
             end do                                                      11d10s20
             ioff=ioff+nvirt(jsbv1)*nn                                   11d10s20
             ioffo=ioffo+nvirt(jsbv1)*nno                               12d28s21
            end if                                                       11d10s20
            if(jsbv12.eq.1)then
             nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                         12d8s20
            else                                                           12d8s20
             nvv=nvirt(jsbv1)*nvirt(jsbv2)                                 12d8s20
            end if                                                         12d8s20
            if(nvv.gt.0)then
             do l=1,4
              if(nfdat(3,l,jsb).gt.0)then
               nn=nfdat(3,l,jsb)*nroot
               if(ldavid)then                                           12d28s21
                nno=nfdat(3,l,jsb)*ngotd                                12d28s21
                nrooto=ngotd                                            12d28s21
               else                                                     12d28s21
                nno=nn                                                  12d28s21
                nrooto=nroot                                            12d28s21
               end if                                                   12d28s21
               joff=ioff
               joffo=ioffo
               do if=1,nfdat(3,l,jsb)                                      11d10s20
                joff=joff+nvv*(ir-1)                                     11d10s20
                joffo=joffo+nvv*(iro-1)                                  11d10s20
                do iv12=1,nvv                                            11d10s20
                 sum2=sum2+vdold(joffo+iv12)*vdinout(joff+iv12)          11d10s20
                end do                                                   11d10s20
                joff=joff+nvv*(nroot+1-ir)                               11d10s20
                joffo=joffo+nvv*(nrooto+1-iro)                          12d28s21
               end do                                                    11d10s20
               ioff=ioff+nvv*nn                                          11d10s20
               ioffo=ioffo+nvv*nno                                      12d28s21
              end if                                                     11d10s20
             end do                                                      11d10s20
            end if                                                       11d10s20
           end if                                                        11d10s20
          end do                                                         11d10s20
         end do                                                          11d10s20
         if(lprt)write(6,*)('dot for doubles: '),sum2                            11d10s20
        end if                                                          2d17s21
        sum1=0d0                                                        11d10s20
        if(nsing.gt.0)then                                              3d19s21
         jvsold=0                                                        11d10s20
         jvstry=0                                                        11d10s20
         if(ldavid)then                                                 12d28s21
          nrooto=ngots                                                  12d28s21
         else                                                           12d28s21
          nrooto=nroot                                                  12d28s21
         end if                                                         12d28s21
         do isb=1,nsymb                                                  11d10s20
          isbv=multh(isb,isymmrci)                                       11d10s20
          do ii=mdon+1,mdoo+1                                            11d10s20
           iarg=ii-mdon                                                  11d10s20
           ndblock=ncsf(iarg)*nff1(ii,isb)                               11d10s20
           if(ndblock.gt.0)then                                          11d10s20
            icol4=ndblock                                                11d10s20
            call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,   3d19s21
     $           i2e)                                                   3d19s21
            ic=0                                                           11d10s20
            do iff=1,nff1(ii,isb)                                        11d10s20
             do m=1,ncsf(iarg)                                            11d10s20
              ic=ic+1                                                    11d10s20
              if(ic.ge.il.and.ic.le.ih)then                              11d10s20
               jvsoldu=jvsold+iro                                        12d14s20
               jvstryu=jvstry+ir                                         12d14s20
               do iv=1,nvirt(isbv)                                       12d14s20
                sum1=sum1+vsold(jvsoldu)*vstry(jvstryu)                  12d14s20
                jvsoldu=jvsoldu+nrooto                                  12d28s21
                jvstryu=jvstryu+nroot                                    12d14s20
               end do                                                    12d14s20
               jvsold=jvsold+nvirt(isbv)*nrooto                         12d28s21
               jvstry=jvstry+nvirt(isbv)*nroot                           12d14s20
              end if                                                     11d10s20
             end do                                                      11d10s20
            end do                                                       11d10s20
           end if                                                        11d10s20
          end do                                                         11d10s20
         end do                                                          11d10s20
         if(lprt)write(6,*)('my part of dot for singles: '),sum1
         call dws_gsumf(sum1,1)
        end if                                                          3d19s21
        dot=sum1+sum2                                                   11d10s20
        if(lprt)then                                                    1d25s21
         write(6,*)('total dot: '),dot                                   11d10s20
         write(6,*)('subtracting off ...')
        end if                                                          1d25s21
        if(ndoub.gt.0.and.iro.le.ngotd)then                             4d14s21
         ioff=0                                                          11d10s20
         ioffo=0                                                         11d10s20
         do jsb=1,nsymb                                                    12d8s20
          jsbv12=multh(jsb,isymmrci)                                       12d8s20
          do jsbv1=1,nsymb                                                 12d8s20
           jsbv2=multh(jsbv1,jsbv12)
           if(jsbv2.ge.jsbv1)then
            if(jsbv12.eq.1.and.nvirt(jsbv1).gt.0)then
             nn=nfdat(3,1,isymmrci)*nroot
             if(ldavid)then                                             12d28s21
              nno=nfdat(3,1,isymmrci)*ngotd                             12d28s21
              nrooto=ngotd                                              12d28s21
             else                                                       12d28s21
              nno=nn                                                    12d28s21
              nrooto=nroot                                              12d28s21
             end if                                                     12d28s21
             joff=ioff
             joffo=ioffo                                                 11d10s20
             do if=1,nfdat(3,1,isymmrci)                                   11d10s20
              joff=joff+nvirt(jsbv1)*(ir-1)                              11d10s20
              joffo=joffo+nvirt(jsbv1)*(iro-1)                              11d10s20
              do iv=1,nvirt(jsbv1)                                       11d10s20
               vdinout(joff+iv)=vdinout(joff+iv)-dot*vdold(joffo+iv)     11d10s20
              end do                                                     11d10s20
              joffo=joffo+nvirt(jsbv1)*(nrooto+1-iro)                   12d28s21
              joff=joff+nvirt(jsbv1)*(nroot+1-ir)                        11d10s20
             end do                                                      11d10s20
             ioff=ioff+nn*nvirt(jsbv1)                                   11d10s20
             ioffo=ioffo+nno*nvirt(jsbv1)                               12d28s21
            end if                                                       11d10s20
            if(jsbv12.eq.1)then
             nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                         12d8s20
            else                                                           12d8s20
             nvv=nvirt(jsbv1)*nvirt(jsbv2)                                 12d8s20
            end if                                                         12d8s20
            if(nvv.gt.0)then
             do l=1,4
              if(nfdat(3,l,jsb).gt.0)then
               nn=nfdat(3,l,jsb)*nroot
               if(ldavid)then                                           12d28s21
                nno=nfdat(3,l,jsb)*ngotd                                12d28s21
                nrooto=ngotd                                            12d28s21
               else                                                     12d28s21
                nno=nn                                                  12d28s21
                nrooto=nroot                                            12d28s21
               end if                                                   12d28s21
               joff=ioff
               joffo=ioffo
               do if=1,nfdat(3,l,jsb)                                      11d10s20
                joff=joff+nvv*(ir-1)                                     11d10s20
                joffo=joffo+nvv*(iro-1)                                  11d10s20
                do iv12=1,nvv                                            11d10s20
                 vdinout(joff+iv12)=vdinout(joff+iv12)                   11d10s20
     $               -dot*vdold(joffo+iv12)                             11d10s20
                end do                                                   11d10s20
                joff=joff+nvv*(nroot+1-ir)                               11d10s20
                joffo=joffo+nvv*(nrooto+1-iro)                          12d28s21
               end do                                                    11d10s20
               ioff=ioff+nvv*nn                                          11d10s20
               ioffo=ioffo+nvv*nno                                      12d28s21
              end if                                                     11d10s20
             end do                                                      11d10s20
            end if                                                       11d10s20
           end if                                                        11d10s20
          end do                                                         11d10s20
         end do                                                          11d10s20
        end if                                                          2d17s21
        if(nsing.gt.0)then                                              3d19s21
         jvsold=0                                                        11d10s20
         jvstry=0                                                        11d10s20
         do isb=1,nsymb                                                  11d10s20
          isbv=multh(isb,isymmrci)                                       11d10s20
          do ii=mdon+1,mdoo+1                                            11d10s20
           iarg=ii-mdon                                                  11d10s20
           ndblock=ncsf(iarg)*nff1(ii,isb)                               11d10s20
           if(ndblock.gt.0)then                                          11d10s20
            icol4=ndblock                                                11d10s20
            call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,   3d19s21
     $           i2e)                                                   3d19s21
            ic=0                                                           11d10s20
            do iff=1,nff1(ii,isb)                                        11d10s20
             do m=1,ncsf(iarg)                                            11d10s20
              ic=ic+1                                                    11d10s20
              if(ic.ge.il.and.ic.le.ih)then                              11d10s20
               jvsoldu=jvsold+iro                                        12d14s20
               jvstryu=jvstry+ir                                         12d14s20
               do iv=1,nvirt(isbv)                                       12d14s20
                vstry(jvstryu)=vstry(jvstryu)-dot*vsold(jvsoldu)         12d14s20
                jvstryu=jvstryu+nroot                                    12d14s20
                jvsoldu=jvsoldu+nrooto                                  12d28s21
               end do                                                    12d14s20
               jvsold=jvsold+nvirt(isbv)*nrooto                         12d28s21
               jvstry=jvstry+nvirt(isbv)*nroot                           12d14s20
              end if                                                     11d10s20
             end do                                                      11d10s20
            end do                                                       11d10s20
           end if                                                        11d10s20
          end do                                                         11d10s20
         end do                                                          11d10s20
        end if                                                          3d19s21
       end do                                                           11d10s20
      end do                                                            11d10s20
c
c     now modified gram-schmidt among the new vectors
c
      if(lprt)write(6,*)('now among new vectors ')                      1d25s21
      isto=0                                                             3d30s21
      do ir=1,nroot                                                     11d10s20
       if(lprt)write(6,*)('for new root no.: '),ir                      1d25s21
       do iro=1,isto                                                     3d30s21
        if(lprt)write(6,*)('for other new root no.: '),iro,ndoub                    1d25s21
        sum2=0d0                                                        11d10s20
        if(ndoub.gt.0)then                                              2d17s21
         ioff=0                                                          11d10s20
         ioffo=0                                                         11d10s20
         do jsb=1,nsymb                                                    12d8s20
          jsbv12=multh(jsb,isymmrci)                                       12d8s20
          do jsbv1=1,nsymb                                                 12d8s20
           jsbv2=multh(jsbv1,jsbv12)
           if(jsbv2.ge.jsbv1)then
            if(jsbv12.eq.1.and.nvirt(jsbv1).gt.0)then
             nn=nfdat(3,1,isymmrci)*nroot
             joff=ioff
             joffo=ioffo                                                 11d10s20
             do if=1,nfdat(3,1,isymmrci)                                   11d10s20
              joff=joff+nvirt(jsbv1)*(ir-1)                              11d10s20
              joffo=joffo+nvirt(jsbv1)*(iro-1)                              11d10s20
              do iv=1,nvirt(jsbv1)                                       11d10s20
               sum2=sum2+vdinout(joffo+iv)*vdinout(joff+iv)                11d10s20
              end do                                                     11d10s20
              joffo=joffo+nvirt(jsbv1)*(nroot+1-iro)                     11d10s20
              joff=joff+nvirt(jsbv1)*(nroot+1-ir)                        11d10s20
             end do                                                      11d10s20
             ioff=ioff+nvirt(jsbv1)*nn                                   11d10s20
             ioffo=ioffo+nvirt(jsbv1)*nn                                 11d10s20
            end if                                                       11d10s20
            if(jsbv12.eq.1)then
             nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                         12d8s20
            else                                                           12d8s20
             nvv=nvirt(jsbv1)*nvirt(jsbv2)                                 12d8s20
            end if                                                         12d8s20
            if(nvv.gt.0)then
             do l=1,4
              if(nfdat(3,l,jsb).gt.0)then
               nn=nfdat(3,l,jsb)*nroot
               joff=ioff
               joffo=ioffo
               do if=1,nfdat(3,l,jsb)                                      11d10s20
                joff=joff+nvv*(ir-1)                                     11d10s20
                joffo=joffo+nvv*(iro-1)                                  11d10s20
                do iv12=1,nvv                                            11d10s20
                 sum2=sum2+vdinout(joffo+iv12)*vdinout(joff+iv12)          11d10s20
                end do                                                   11d10s20
                joff=joff+nvv*(nroot+1-ir)                               11d10s20
                joffo=joffo+nvv*(nroot+1-iro)                            11d10s20
               end do                                                    11d10s20
               ioff=ioff+nvv*nn                                          11d10s20
               ioffo=ioffo+nvv*nn                                        11d10s20
              end if                                                     11d10s20
             end do                                                      11d10s20
            end if                                                       11d10s20
           end if                                                        11d10s20
          end do                                                         11d10s20
         end do                                                          11d10s20
         if(lprt)write(6,*)('dot for doubles: '),sum2                    1d25s21
        end if                                                          2d17s21
        sum1=0d0                                                        11d10s20
        if(nsing.gt.0)then                                              3d19s21
         jvsold=0                                                        11d10s20
         jvstry=0                                                        11d10s20
         do isb=1,nsymb                                                  11d10s20
          isbv=multh(isb,isymmrci)                                       11d10s20
          do ii=mdon+1,mdoo+1                                            11d10s20
           iarg=ii-mdon                                                  11d10s20
           ndblock=ncsf(iarg)*nff1(ii,isb)                               11d10s20
           if(ndblock.gt.0)then                                          11d10s20
            icol4=ndblock                                                11d10s20
            call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,   3d19s21
     $           i2e)                                                   3d19s21
            ic=0                                                           11d10s20
            do iff=1,nff1(ii,isb)                                        11d10s20
             do m=1,ncsf(iarg)                                            11d10s20
              ic=ic+1                                                    11d10s20
              if(ic.ge.il.and.ic.le.ih)then                              11d10s20
               jvsoldu=jvsold+iro                                         12d14s20
               jvstryu=jvstry+ir                                          12d14s20
               do iv=1,nvirt(isbv)                                        12d14s20
                sum1=sum1+vstry(jvsoldu)*vstry(jvstryu)                   12d14s20
                jvsoldu=jvsoldu+nroot                                     12d14s20
                jvstryu=jvstryu+nroot                                     12d14s20
               end do                                                     12d14s20
               jvsold=jvsold+nvirt(isbv)*nroot                           12d14s20
               jvstry=jvstry+nvirt(isbv)*nroot                           12d14s20
              end if                                                     11d10s20
             end do                                                      11d10s20
            end do                                                       11d10s20
           end if                                                        11d10s20
          end do                                                         11d10s20
         end do                                                          11d10s20
         if(lprt)write(6,*)('my part of dot for singles: '),sum1         1d25s21
         call dws_gsumf(sum1,1)
        end if                                                          3d19s21
        dot=sum1+sum2                                                   11d10s20
        if(lprt)then                                                    1d25s21
         write(6,*)('total dot: '),dot                                   11d10s20
         write(6,*)('project out ...')
        end if                                                          1d25s21
        if(ndoub.gt.0)then                                              2d17s21
         ioff=0                                                          11d10s20
         ioffo=0                                                         11d10s20
         do jsb=1,nsymb                                                    12d8s20
          jsbv12=multh(jsb,isymmrci)                                       12d8s20
          do jsbv1=1,nsymb                                                 12d8s20
           jsbv2=multh(jsbv1,jsbv12)
           if(jsbv2.ge.jsbv1)then
            if(jsbv12.eq.1.and.nvirt(jsbv1).gt.0)then
             nn=nfdat(3,1,isymmrci)*nroot
             joff=ioff
             joffo=ioffo                                                 11d10s20
             do if=1,nfdat(3,1,isymmrci)                                   11d10s20
              joff=joff+nvirt(jsbv1)*(ir-1)                              11d10s20
              joffo=joffo+nvirt(jsbv1)*(iro-1)                              11d10s20
              do iv=1,nvirt(jsbv1)                                       11d10s20
               vdinout(joff+iv)=vdinout(joff+iv)-dot*vdinout(joffo+iv)   11d10s20
              end do                                                     11d10s20
              joffo=joffo+nvirt(jsbv1)*(nroot+1-iro)                     11d10s20
              joff=joff+nvirt(jsbv1)*(nroot+1-ir)                        11d10s20
             end do                                                      11d10s20
             ioff=ioff+nvirt(jsbv1)*nn                                   11d10s20
             ioffo=ioffo+nvirt(jsbv1)*nn                                 11d10s20
            end if                                                       11d10s20
            if(jsbv12.eq.1)then
             nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                         12d8s20
            else                                                           12d8s20
             nvv=nvirt(jsbv1)*nvirt(jsbv2)                                 12d8s20
            end if                                                         12d8s20
            if(nvv.gt.0)then
             do l=1,4
              if(nfdat(3,l,jsb).gt.0)then
               nn=nfdat(3,l,jsb)*nroot
               joff=ioff
               joffo=ioffo
               do if=1,nfdat(3,l,jsb)                                      11d10s20
                joff=joff+nvv*(ir-1)                                     11d10s20
                joffo=joffo+nvv*(iro-1)                                  11d10s20
                do iv12=1,nvv                                            11d10s20
                 vdinout(joff+iv12)=vdinout(joff+iv12)                   11d10s20
     $                -dot*vdinout(joffo+iv12)                           11d10s20
                end do                                                   11d10s20
                joff=joff+nvv*(nroot+1-ir)                               11d10s20
                joffo=joffo+nvv*(nroot+1-iro)                            11d10s20
               end do                                                    11d10s20
               ioff=ioff+nvv*nn                                          11d10s20
               ioffo=ioffo+nvv*nn                                        11d10s20
              end if                                                     11d10s20
             end do                                                      11d10s20
            end if                                                       11d10s20
           end if                                                        11d10s20
          end do                                                         11d10s20
         end do                                                          11d10s20
        end if                                                          2d17s21
        if(nsing.gt.0)then                                              3d19s21
         jvsold=0                                                        11d10s20
         jvstry=0                                                        11d10s20
         do isb=1,nsymb                                                  11d10s20
          isbv=multh(isb,isymmrci)                                       11d10s20
          do ii=mdon+1,mdoo+1                                            11d10s20
           iarg=ii-mdon                                                  11d10s20
           ndblock=ncsf(iarg)*nff1(ii,isb)                               11d10s20
           if(ndblock.gt.0)then                                          11d10s20
            icol4=ndblock                                                11d10s20
            call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,   3d19s21
     $           i2e)                                                   3d19s21
            ic=0                                                           11d10s20
            do iff=1,nff1(ii,isb)                                        11d10s20
             do m=1,ncsf(iarg)                                            11d10s20
              ic=ic+1                                                    11d10s20
              if(ic.ge.il.and.ic.le.ih)then                              11d10s20
               jvsoldu=jvsold+iro                                        12d14s20
               jvstryu=jvstry+ir                                         12d14s20
               do iv=1,nvirt(isbv)                                       11d10s20
                vstry(jvstryu)=vstry(jvstryu)-dot*vstry(jvsoldu)         12d14s20
                jvstryu=jvstryu+nroot                                    12d14s20
                jvsoldu=jvsoldu+nroot                                    12d14s20
               end do                                                    11d10s20
               jvsold=jvsold+nvirt(isbv)*nroot                           12d14s20
               jvstry=jvstry+nvirt(isbv)*nroot                           12d14s20
              end if                                                     11d10s20
             end do                                                      11d10s20
            end do                                                       11d10s20
           end if                                                        11d10s20
          end do                                                         11d10s20
         end do                                                          11d10s20
        end if                                                          3d19s21
       end do                                                           11d10s20
c
c     normalize                                                         11d10s20
c
       if(lprt)write(6,*)('normalize'),isto                                  1d25s21
       sum2=0d0                                                         11d10s20
       if(ndoub.gt.0)then                                               2d17s21
        ioff=0                                                           11d10s20
        do jsb=1,nsymb                                                    12d8s20
         jsbv12=multh(jsb,isymmrci)                                       12d8s20
         do jsbv1=1,nsymb                                                 12d8s20
          jsbv2=multh(jsbv1,jsbv12)
          if(jsbv2.ge.jsbv1)then
           if(jsbv12.eq.1.and.nvirt(jsbv1).gt.0)then
            nn=nfdat(3,1,isymmrci)*nroot
            joff=ioff
            do if=1,nfdat(3,1,isymmrci)                                   11d10s20
             joff=joff+nvirt(jsbv1)*(ir-1)                               11d10s20
             do iv=1,nvirt(jsbv1)                                        11d10s20
              sum2=sum2+vdinout(joff+iv)**2                              11d10s20
             end do                                                      11d10s20
             joff=joff+nvirt(jsbv1)*(nroot+1-ir)                         11d10s20
            end do                                                       11d10s20
            ioff=ioff+nvirt(jsbv1)*nn                                    11d10s20
           end if                                                        11d10s20
           if(jsbv12.eq.1)then
            nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                         12d8s20
           else                                                           12d8s20
            nvv=nvirt(jsbv1)*nvirt(jsbv2)                                 12d8s20
           end if                                                         12d8s20
           if(nvv.gt.0)then
            do l=1,4
             if(nfdat(3,l,jsb).gt.0)then
              nn=nfdat(3,l,jsb)*nroot
              joff=ioff
              do if=1,nfdat(3,l,jsb)                                      11d10s20
               joff=joff+nvv*(ir-1)                                      11d10s20
               do iv12=1,nvv                                             11d10s20
                sum2=sum2+vdinout(joff+iv12)**2
               end do                                                    11d10s20
               joff=joff+nvv*(nroot+1-ir)                                11d10s20
              end do                                                     11d10s20
              ioff=ioff+nvv*nn                                           11d10s20
             end if                                                      11d10s20
            end do                                                       11d10s20
           end if                                                        11d10s20
          end if                                                         11d10s20
         end do                                                          11d10s20
        end do                                                           11d10s20
        if(lprt)write(6,*)('dot for doubles: '),sum2                     1d25s21
       end if                                                           2d17s21
       sum1=0d0                                                         11d10s20
       if(nsing.gt.0)then                                               3d19s21
        jvstry=0                                                         11d10s20
        do isb=1,nsymb                                                   11d10s20
         isbv=multh(isb,isymmrci)                                        11d10s20
         do ii=mdon+1,mdoo+1                                             11d10s20
          iarg=ii-mdon                                                   11d10s20
          ndblock=ncsf(iarg)*nff1(ii,isb)                                11d10s20
          if(ndblock.gt.0)then                                           11d10s20
           icol4=ndblock                                                 11d10s20
           call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)   11d10s20
           ic=0                                                           11d10s20
           do iff=1,nff1(ii,isb)                                         11d10s20
            do m=1,ncsf(iarg)                                            11d10s20
             ic=ic+1                                                     11d10s20
             if(ic.ge.il.and.ic.le.ih)then                               11d10s20
              jvstryu=jvstry+ir                                          12d14s20
              do iv=1,nvirt(isbv)                                        11d10s20
               sum1=sum1+vstry(jvstryu)**2                               12d14s20
               jvstryu=jvstryu+nroot                                     12d14s20
              end do                                                     11d10s20
              jvstry=jvstry+nvirt(isbv)*nroot                            12d14s20
             end if                                                      11d10s20
            end do                                                       11d10s20
           end do                                                        11d10s20
          end if                                                         11d10s20
         end do                                                          11d10s20
        end do                                                           11d10s20
        if(lprt)write(6,*)('my part of dot for singles: '),sum1          1d25s21
        call dws_gsumf(sum1,1)
       end if                                                           3d19s21
       dot=sum1+sum2                                                    11d10s20
       xnorm=1d0/sqrt(dot)                                              11d10s20
       if(lprt)then                                                     1d25s21
        write(6,*)('for root '),ir
        write(6,*)('total dot: '),dot                                    11d10s20
        write(6,*)('xnorm = '),xnorm
       end if                                                           1d25s21
       if(iamconverged(ir).eq.0)then                                    4d21s21
        isto=isto+1                                                       3d30s21
        if(ndoub.gt.0)then                                               2d17s21
         ioff=0                                                           11d10s20
         do jsb=1,nsymb                                                    12d8s20
          jsbv12=multh(jsb,isymmrci)                                       12d8s20
          do jsbv1=1,nsymb                                                 12d8s20
           jsbv2=multh(jsbv1,jsbv12)
           if(jsbv2.ge.jsbv1)then
            if(jsbv12.eq.1.and.nvirt(jsbv1).gt.0)then
             nn=nfdat(3,1,isymmrci)*nroot
             joff=ioff
             joffs=ioff                                                 3d30s21
             do if=1,nfdat(3,1,isymmrci)                                   11d10s20
              joff=joff+nvirt(jsbv1)*(ir-1)                               11d10s20
              joffs=joffs+nvirt(jsbv1)*(isto-1)                          3d30s21
              do iv=1,nvirt(jsbv1)                                        11d10s20
               vdinout(joffs+iv)=vdinout(joff+iv)*xnorm                 3d30s21
              end do                                                      11d10s20
              joff=joff+nvirt(jsbv1)*(nroot+1-ir)                         11d10s20
              joffs=joffs+nvirt(jsbv1)*(nroot+1-isto)                    3d30s21
             end do                                                       11d10s20
             ioff=ioff+nvirt(jsbv1)*nn                                    11d10s20
            end if                                                        11d10s20
            if(jsbv12.eq.1)then
             nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                         12d8s20
            else                                                           12d8s20
             nvv=nvirt(jsbv1)*nvirt(jsbv2)                                 12d8s20
            end if                                                         12d8s20
            if(nvv.gt.0)then
             do l=1,4
              if(nfdat(3,l,jsb).gt.0)then
               nn=nfdat(3,l,jsb)*nroot
               joff=ioff
               joffs=ioff                                                3d30s21
               do if=1,nfdat(3,l,jsb)                                      11d10s20
                joff=joff+nvv*(ir-1)                                     3d30s21
                joffs=joffs+nvv*(isto-1)                                  3d30s21
                do iv12=1,nvv                                             11d10s20
                 vdinout(joffs+iv12)=vdinout(joff+iv12)*xnorm            3d30s21
                end do                                                    11d10s20
                joff=joff+nvv*(nroot+1-ir)                                11d10s20
                joffs=joffs+nvv*(nroot+1-isto)                          3d30s21
               end do                                                     11d10s20
               ioff=ioff+nvv*nn                                           11d10s20
              end if                                                      11d10s20
             end do                                                       11d10s20
            end if                                                        11d10s20
           end if                                                         11d10s20
          end do                                                          11d10s20
         end do                                                           11d10s20
        end if                                                           2d17s21
        if(nsing.gt.0)then                                               3d19s21
         jvstry=0                                                         11d10s20
         do isb=1,nsymb                                                   11d10s20
          isbv=multh(isb,isymmrci)                                        11d10s20
          irow4=nvirt(isbv)*nroot
          do ii=mdon+1,mdoo+1                                             11d10s20
           iarg=ii-mdon                                                   11d10s20
           ndblock=ncsf(iarg)*nff1(ii,isb)                                11d10s20
           if(ndblock.gt.0)then                                           11d10s20
            icol4=ndblock                                                 11d10s20
            call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,   3d30s21
     $           i2e)                                                   3d30s21
            ic=0                                                           11d10s20
            ncol=ih+1-il
            jvstry0=jvstry+1                                              11d10s20
            do iff=1,nff1(ii,isb)                                         11d10s20
             do m=1,ncsf(iarg)                                            11d10s20
              ic=ic+1                                                     11d10s20
              if(ic.ge.il.and.ic.le.ih)then                               11d10s20
               jvstryu=jvstry+ir                                          12d14s20
               jvstrysu=jvstry+isto                                      3d30s21
               do iv=1,nvirt(isbv)                                        11d10s20
                vstry(jvstrysu)=vstry(jvstryu)*xnorm                       12d14s20
                jvstryu=jvstryu+nroot                                     12d14s20
                jvstrysu=jvstrysu+nroot                                 3d30s21
               end do                                                     11d10s20
               jvstry=jvstry+nvirt(isbv)*nroot                            12d14s20
              end if                                                      11d10s20
             end do                                                       11d10s20
            end do                                                        11d10s20
           end if                                                         11d10s20
          end do                                                          11d10s20
         end do                                                           11d10s20
        end if                                                           3d19s21
       end if                                                           3d30s21
      end do                                                            11d10s20
      if(ndoub.gt.0)then                                                2d17s21
       if(isto.ne.nroot)then                                             3d30s21
        ioffo=0                                                         3d30s21
        ioffn=0                                                         3d30s21
        do isb=1,nsymb                                                  3d30s21
         isbv12=multh(isb,isymmrci)                                     3d30s21
         do isbv1=1,nsymb                                               3d30s21
          isbv2=multh(isbv1,isbv12)                                     3d30s21
          if(isbv2.ge.isbv1)then                                        3d30s21
           if(isbv1.eq.isbv2)then                                       3d30s21
            do if=1,nfdat(3,1,isb)                                      3d30s21
             do ir=1,isto                                                3d30s21
              do iv=1,nvirt(isbv1)                                      3d30s21
               vdinout(ioffn+iv)=vdinout(ioffo+iv)                      3d30s21
              end do                                                    3d30s21
              ioffn=ioffn+nvirt(isbv1)                                  3d30s21
              ioffo=ioffo+nvirt(isbv1)                                  3d30s21
             end do                                                     3d30s21
             ioffo=ioffo+nvirt(isbv1)*(nroot-isto)                       3d30s21
            end do                                                      3d30s21
            nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       3d30s21
           else                                                         3d30s21
            nvv=nvirt(isbv1)*nvirt(isbv2)                               3d30s21
           end if                                                       3d30s21
           do l=1,4                                                     3d30s21
            if(nfdat(3,l,isb).gt.0)then                                 3d30s21
             do if=1,nfdat(3,l,isb)                                     3d30s21
              do ir=1,isto                                               3d30s21
               do ivv=1,nvv                                             3d30s21
                vdinout(ioffn+ivv)=vdinout(ioffo+ivv)                   3d30s21
               end do                                                   3d30s21
               ioffn=ioffn+nvv                                          3d30s21
               ioffo=ioffo+nvv                                          3d30s21
              end do                                                    3d30s21
              ioffo=ioffo+nvv*(nroot-isto)                               3d30s21
             end do                                                     3d30s21
            end if                                                      3d30s21
           end do                                                       3d30s21
          end if                                                        3d30s21
         end do                                                         3d30s21
        end do                                                          3d30s21
       end if                                                           3d30s21
       call dws_bcast(vdinout,ndoub*isto)                               3d30s21
      end if                                                            2d17s21
      if(nsing.gt.0)then                                                3d19s21
       jvstry=1                                                          11d10s20
       jvstryn=1                                                        3d30s21
       jhss=1                                                            11d10s20
       do isb=1,nsymb                                                    11d10s20
        isbv=multh(isb,isymmrci)                                         11d10s20
        irow4=nvirt(isbv)*nroot                                         3d30s21
        irown=nvirt(isbv)*isto                                          3d30s21
        irow=irown                                                      3d30s21
        do ii=mdon+1,mdoo+1                                              11d10s20
         iarg=ii-mdon                                                    11d10s20
         ndblock=ncsf(iarg)*nff1(ii,isb)                                 11d10s20
         if(ndblock.gt.0)then                                            11d10s20
          icol4=ndblock                                                  11d10s20
          call ilimts(1,icol4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)   11d10s20
          il8=il                                                         11d10s20
          ih8=ih                                                         11d10s20
          nhere=ih+1-il                                                  11d10s20
          if(isto.ne.nroot)then                                         3d30s21
           jo=jvstry-1                                                  3d30s21
           jn=jvstryn-1                                                 3d30s21
           do if=1,nhere                                                3d30s21
            do iv=1,nvirt(isbv)                                         3d30s21
             do ir=1,isto                                               3d30s21
              vstry(jn+ir)=vstry(jo+ir)                                 3d30s21
             end do                                                     3d30s21
             jn=jn+isto                                                 3d30s21
             jo=jo+nroot                                                3d30s21
            end do                                                      3d30s21
           end do                                                       3d30s21
          end if                                                        3d30s21
          call ddi_put(bc,ibc,nhsdiag(ii,isb,1),i1,irow,il8,ih8,        11d15s22
     $         vstry(jvstryn))                                          11d15s22
          igtmp=ibcoff                                                   12d15s20
          ibcoff=igtmp+irown*nhere                                      3d30s21
          call enough('updatexc.  3',bc,ibc)
          jjvstry=jvstryn                                               3d30s21
          jjgtmp=igtmp                                                   12d15s20
          ic=0                                                           12d15s20
          do iff=1,nff1(ii,isb)                                          11d10s20
           khss=jhss                                                     12d17s20
           do m=1,ncsf(iarg)                                             11d10s20
            ic=ic+1                                                      11d10s20
            if(ic.ge.il.and.ic.le.ih)then                                11d10s20
             khss=jhss                                                   12d14s20
             do iv=1,nvirt(isbv)                                         11d10s20
              do ir=1,isto                                              3d30s21
               bc(jjgtmp)=vstry(jjvstry)*hss(khss)                       12d15s20
               jjgtmp=jjgtmp+1                                           12d15s20
               jjvstry=jjvstry+1                                         12d15s20
              end do                                                     12d15s20
              khss=khss+1                                                12d15s20
             end do                                                      12d15s20
            end if                                                       12d15s20
           end do                                                        12d15s20
           jhss=khss                                                     12d15s20
          end do                                                         12d15s20
          call ddi_put(bc,ibc,nhsdiag(ii,isb,2),i1,irow,il8,ih8,        11d15s22
     $         bc(igtmp))                                               11d15s22
          ibcoff=igtmp                                                   12d15s20
          jvstry=jvstry+irow4*nhere                                      3d30s21
          jvstryn=jvstryn+irown*nhere                                   3d30s21
         end if                                                          11d10s20
        end do                                                           11d10s20
       end do                                                            11d10s20
      end if                                                            3d19s21
      call dws_synca                                                    11d10s20
      return
      end
