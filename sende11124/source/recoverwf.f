c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine recoverwf(iunc,nroot,ntot,nlzzci,lambdaci,nsing,       8d16s22
     $     ihsdiag,nff1,iff1,ncsf,mdon,mdoo,nsymb,multh,nvirt,isymmrci, 8d16s22
     $     vs,vso,ndoub,nff22,nfdat,ivdinout,mdoub,ncsf2,ihddiag,nff2,  8d19s22
     $     iff2,norb,hss,hdd,bc,ibc)                                    11d10s22
      implicit real*8 (a-h,o-z)                                         8d16s22
      integer*8 ihsdiag(mdoo+1,nsymb,2),ihddiag(mdoo+1,nsymb,2),        8d19s22
     $     i18,i28,i38,i48                                              8d19s22
      logical lpr                                                       6d23s23
      dimension iunc(2),iuncr(2),nff1(mdoo+1,nsymb,*),multh(8,8),       8d18s22
     $     nvirt(*),ncsf(*),nfdat(5,4,*),hss(*),vs(*),vso(*),           8d19s22
     $     nff2(mdoo+1,nsymb),hdd(*),ncsf2(4,*)                         8d19s22
      include "common.print"                                            6d23s23
      include "common.store"                                            8d16s22
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      lpr=iprtr(29).ne.0                                                6d23s23
      if(mynowprog.eq.0)then                                            8d16s22
       read(2)iuncr,ethred                                               8d16s22
       if(lpr)then                                                      6d23s23
        write(6,*)('in recoverwf ')
        write(6,*)('iuncr,ethred '),iuncr,ethred
       end if                                                           6d23s23
       ier=0                                                             8d16s22
       do j=1,2                                                          8d16s22
        if(iunc(j).ne.iuncr(j))ier=ier+1                                 8d16s22
       end do                                                            8d16s22
       if(ier.ne.0)then                                                  8d16s22
        write(6,*)('for iunc: ')                                         8d16s22
        write(6,*)('got  :'),iuncr                                       8d16s22
        write(6,*)('want :'),iunc                                        8d16s22
        stop 'restart'                                                   8d16s22
       end if                                                            8d16s22
       read(2)
       read(2)
       read(2)
       read(2)nlzz,lamb
       if(nlzz.ne.nlzzci)ier=ier+1                                      8d16s22
       if(lamb.ne.lambdaci)ier=ier+1                                    8d16s22
       if(ier.ne.0)then                                                  8d16s22
        write(6,*)('for nlzzci: ')                                         8d16s22
        write(6,*)('got  :'),nlzz,lamb                                  8d16s22
        write(6,*)('want :'),nlzzci,lambdaci                            8d16s22
        stop 'restart'                                                   8d16s22
       end if                                                            8d16s22
      end if                                                            8d16s22
      mdoop=mdoo+1                                                      8d16s22
      if(lpr)write(6,*)('nsing '),nsing                                 6d23s23
      if(nsing.ne.0)then                                                8d16s22
       if(mynowprog.eq.0)then                                           8d16s22
        read(2)mff1                                                     8d16s22
        if(lpr)write(6,*)('mff1 ')                                      6d23s23
        nff1r=ibcoff                                                    8d16s22
        iff1r=nff1r+nsymb*mdoop
        icsfr=iff1r+mff1                                                8d16s22
        ibcoff=icsfr+mdoop-mdon                                         8d16s22
        call loadsinga(ibc(nff1r),ibc(iff1r),mff1,nsymb,mdon,mdoo,      8d16s22
     $       ibc(icsfr),2)                                              12d22s22
        if(lpr)write(6,*)('back from loadsinga')                        6d23s23
        call check1(ibc(nff1r),nff1,ibc(iff1r),iff1,mff1,nsymb,mdon,    8d16s22
     $       mdoo,ncsf,ibc(icsfr))                                      12d22s22
        if(lpr)write(6,*)('back from check1')                           6d23s23
        ibcoff=nff1r                                                    8d18s22
       end if                                                           8d16s22
      else                                                              8d18s22
       if(mynowprog.eq.0)read(2)                                        8d18s22
      end if                                                            8d18s22
      if(mynowprog.eq.0)read(2)mff2                                     8d18s22
      icsfr2=ibcoff                                                     4d19s23
      if(ndoub.gt.0)then                                                8d18s22
       if(mynowprog.eq.0)then                                           8d18s22
        if(iunc(2).eq.0.and.mff2.lt.0)then                              8d18s22
         mff2a=-mff2                                                    8d18s22
         read(2)ndoubr,mdoubr                                           8d18s22
         if(ndoubr.ne.ndoub.or.mdoubr.ne.mdoub)then                     8d18s22
          write(6,*)('ndoub,mdoub: ')                                   8d18s22
          write(6,*)('got  :'),ndoubr,mdoubr                            8d18s22
          write(6,*)('want :'),ndoub,mdoub                              8d18s22
          stop 'restart'                                                8d18s22
         end if                                                         8d18s22
         iff22r=ibcoff                                                  8d18s22
         nff22r=iff22r+mff2a                                            8d18s22
         nfdatr=nff22r+mdoop*nsymb                                      8d18s22
         icsf=nfdatr+10*nsymb                                           8d18s22
         icsf2=icsf+mdoop-mdon                                          8d18s22
         ibcoff=icsf2+2*(mdoop-mdon)                                    8d18s22
         call enough('recoverwf.  1',bc,ibc)
         nall=ndoub*nroot                                               8d18s22
         if(lpr)write(6,*)('calling loaddoubc')                         6d23s23
         call loaddoubc(ibc(iff22r),ibc(nff22r),ibc(nfdatr),            8d18s22
     $       ibc(icsf),ibc(icsf2),mdon,mdoo,nsymb,mff2a,bc(ivdinout),   8d18s22
     $       nall,isymmrci,nroot,multh,nvirt,2)                         8d18s22
         if(lpr)write(6,*)('back loaddoubc')                            6d23s23
         call check2(ibc(iff22r),ibc(iff2),ibc(nff22r),ibc(nff22),      2d10s23
     $        ibc(nfdatr),nfdat,ibc(icsf),ncsf,ibc(icsf2),ncsf2,        8d18s22
     $        mff2a,mdon,mdoo,nsymb,norb,bc(iff22r),bc,ibc)             11d14s22
         if(lpr)write(6,*)('back check2')                               6d23s23
         ibcoff=iff22r                                                  8d18s22
        else if(iunc(2).ne.0.and.mff2.gt.0)then                         8d18s22
         icsfr2=ibcoff                                                  8d19s22
         ibcoff=icsfr2+mdoop-mdon                                       8d19s22
         nff2r=ibcoff                                                    8d16s22
         iff2r=nff2r+nsymb*mdoop
         icsfr=iff2r+mff2                                                8d16s22
         ibcoff=icsfr+mdoop-mdon                                        8d19s22
         if(lpr)write(6,*)('calling loaddouba')                         6d23s23
         call loaddouba(ibc(nff2r),ibc(iff2r),mff2,nsymb,mdon,mdoo,      8d16s22
     $       ibc(icsfr),ibc(icsfr2),2)                                  8d19s22
         if(lpr)write(6,*)('back from loaddouba')                       6d23s23
         call check3(ibc(nff2r),nff2,ibc(iff2r),ibc(iff2),mff2,         8d19s22
     $        nsymb,mdon,mdoo,ncsf,ibc(icsfr),ncsf2,ibc(icsfr2))        8d19s22
         if(lpr)write(6,*)('back check3')                               6d23s23
         ibcoff=nff2r                                                    8d18s22
        else                                                            8d18s22
         write(6,*)('iunc(2) = '),iunc(2),('while mff2 = '),mff2,('!!!')8d18s22
         stop 'restart'                                                 8d18s22
        end if                                                          8d18s22
       end if                                                           8d18s22
       if(iunc(2).eq.0)then                                             8d18s22
        nall=ndoub*nroot                                                8d18s22
        call dws_synca
        call dws_bcast(bc(ivdinout),nall)                               8d18s22
       end if                                                           8d18s22
      else                                                              8d18s22
       if(mynowprog.eq.0.and.mff2.ne.0)then                             8d18s22
        write(6,*)('mff2 read in is '),mff2,(', but ndoub = 0!')        8d18s22
        stop 'restart'                                                  8d18s22
       end if                                                           8d18s22
      end if                                                            8d18s22
      if(lpr)write(6,*)('loadsingb? ')                                  6d23s23
      if(nsing.gt.0)then                                                8d18s22
       call loadsingb(ihsdiag,nff1,nsymb,mdon,mdoo,nvirt,multh,         8d18s22
     $      isymmrci,ncsf,nroot,mynowprog,idum,2,bc,ibc,mddilow,mddihig)12d13s22
      end if                                                            8d16s22
      if(lpr)write(6,*)('loaddoubb? ')                                  6d23s23
      if(ndoub.ne.0.and.iunc(2).ne.0)then                               8d19s22
       call loaddoubb(ihddiag,nff2,nsymb,mdon,mdoo,nvirt,               8d19s22
     $      multh,isymmrci,ncsf,ibc(icsfr2),nroot,mynowprog,idum,2,bc,  11d10s22
     $     ibc,mddilow,mddihig)                                         12d13s22
      end if                                                            8d19s22
      close(unit=2)                                                     8d19s22
      if(lpr)write(6,*)('orthogonalize')                                6d23s23
c
c     now we need to orthogonalize the vectors
c
      do ir=1,nroot                                                     8d18s22
       irm=ir-1                                                         8d18s22
       do jr=1,ir-1                                                     8d19s22
        jrm=jr-1                                                        11d2s22
        dot=0d0                                                         8d19s22
        if(nsing.ne.0)then                                               8d18s22
         do isb=1,nsymb                                                  8d18s22
          isbv=multh(isb,isymmrci)                                       8d18s22
          do nclop=mdon+1,mdoo+1                                         8d18s22
           if(min(nff1(nclop,isb,1),nvirt(isbv)).gt.0)then               8d18s22
            iarg=nclop-mdon                                              8d18s22
            ncol=ncsf(iarg)*nff1(nclop,isb,1)                            8d18s22
            nrow=nroot*nvirt(isbv)                                       8d18s22
            call ilimts(ncol,1,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e) 8d18s22
            nhere=ih+1-il                                                8d18s22
            ivtmp=ibcoff                                                 8d18s22
            ibcoff=ivtmp+nhere*nrow                                      8d18s22
            call enough('recoverwf.  2',bc,ibc)
            i18=1                                                        8d18s22
            i28=nrow                                                     8d18s22
            i38=il                                                       8d18s22
            i48=ih                                                       8d18s22
            call ddi_get(bc,ibc,ihsdiag(nclop,isb,1),i18,i28,i38,i48,   12d22s22
     $           bc(ivtmp))                                             12d22s22
            do i=0,nhere-1                                               8d18s22
             do iv=0,nvirt(isbv)-1                                       8d18s22
              iad=ivtmp+irm+nroot*(iv+nvirt(isbv)*i)                     8d18s22
              jad=ivtmp+jrm+nroot*(iv+nvirt(isbv)*i)                     8d18s22
              dot=dot+bc(iad)*bc(jad)                                   8d19s22
             end do                                                      8d18s22
            end do                                                       8d18s22
            ibcoff=ivtmp                                                 8d18s22
           end if
          end do                                                         8d18s22
         end do                                                          8d18s22
        end if                                                           8d18s22
        if(ndoub.gt.0)then                                               8d18s22
         if(iunc(2).eq.0)then                                            8d18s22
          ioff=ivdinout                                                  8d18s22
          do isb=1,nsymb                                                 8d18s22
           isbv12=multh(isb,isymmrci)                                    8d18s22
           do isbv1=1,nsymb                                              8d18s22
            isbv2=multh(isbv1,isbv12)                                    8d18s22
            if(isbv2.ge.isbv1)then                                       8d18s22
             if(isbv1.eq.isbv2)then                                      8d18s22
              do i=0,nfdat(3,1,isb)-1                                    8d18s22
               iad=ioff+nvirt(isbv1)*(irm+nroot*i)                       8d18s22
               jad=ioff+nvirt(isbv1)*(jrm+nroot*i)                       8d18s22
               do iv=mynowprog,nvirt(isbv1)-1,mynprocg                   8d18s22
                dot=dot+bc(iad+iv)*bc(jad+iv)                           8d19s22
               end do                                                    8d18s22
              end do                                                     8d18s22
              ioff=ioff+nvirt(isbv1)*nroot*nfdat(3,1,isb)                8d18s22
              nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                      8d18s22
             else                                                        8d18s22
              nvv=nvirt(isbv1)*nvirt(isbv2)                              8d18s22
             end if                                                      8d18s22
             do l=1,4                                                    8d18s22
              if(nfdat(3,l,isb).gt.0)then
               do i=0,nfdat(3,l,isb)-1                                    8d18s22
                iad=ioff+nvv*(irm+nroot*i)                                8d18s22
                jad=ioff+nvv*(jrm+nroot*i)                                8d18s22
                do ivv=mynowprog,nvv-1,mynprocg                           8d18s22
                 dot=dot+bc(iad+ivv)*bc(jad+ivv)                        8d19s22
                end do                                                    8d18s22
               end do                                                     8d18s22
               ioff=ioff+nfdat(3,l,isb)*nvv*nroot                        8d19s22
              end if                                                     8d18s22
             end do                                                      8d18s22
            end if                                                       8d18s22
           end do                                                        8d18s22
          end do                                                         8d18s22
         else                                                            8d18s22
          call ucdoubop(nsymb,isymmrci,multh,nvirt,mdon,mdoo,nff2,ncsf, 8d24s22
     $        ncsf2,ihddiag,irm,jrm,dot,1,nroot,bc,ibc)                 11d10s22
         end if                                                          8d18s22
        end if                                                           8d18s22
        call dws_gsumf(dot,1)                                            8d18s22
        if(lpr)write(6,*)('dot of '),ir,jr,(' is '),dot
        if(nsing.ne.0)then                                               8d18s22
         do isb=1,nsymb                                                  8d18s22
          isbv=multh(isb,isymmrci)                                       8d18s22
          do nclop=mdon+1,mdoo+1                                         8d18s22
           if(min(nff1(nclop,isb,1),nvirt(isbv)).gt.0)then               8d18s22
            iarg=nclop-mdon                                              8d18s22
            ncol=ncsf(iarg)*nff1(nclop,isb,1)                            8d18s22
            nrow=nroot*nvirt(isbv)                                       8d18s22
            call ilimts(ncol,1,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e) 8d18s22
            nhere=ih+1-il                                                8d18s22
            ivtmp=ibcoff                                                 8d18s22
            ibcoff=ivtmp+nhere*nrow                                      8d18s22
            call enough('recoverwf.  3',bc,ibc)
            i18=1                                                        8d18s22
            i28=nrow                                                     8d18s22
            i38=il                                                       8d18s22
            i48=ih                                                       8d18s22
            call ddi_get(bc,ibc,ihsdiag(nclop,isb,1),i18,i28,i38,i48,   12d22s22
     $           bc(ivtmp))                                             12d22s22
            do i=0,nhere-1                                               8d18s22
             do iv=0,nvirt(isbv)-1                                       8d18s22
              iad=ivtmp+irm+nroot*(iv+nvirt(isbv)*i)                     8d18s22
              jad=ivtmp+jrm+nroot*(iv+nvirt(isbv)*i)                     8d18s22
              bc(iad)=bc(iad)-dot*bc(jad)                               8d19s22
             end do                                                      8d18s22
            end do                                                       8d18s22
            call ddi_put(bc,ibc,ihsdiag(nclop,isb,1),i18,i28,i38,i48,   11d15s22
     $           bc(ivtmp))                                             11d15s22
            ibcoff=ivtmp                                                 8d18s22
           end if
          end do                                                         8d18s22
         end do                                                          8d18s22
        end if                                                           8d18s22
        if(ndoub.gt.0)then                                               8d18s22
         if(iunc(2).eq.0)then                                            8d18s22
          ioff=ivdinout                                                  8d18s22
          do isb=1,nsymb                                                 8d18s22
           isbv12=multh(isb,isymmrci)                                    8d18s22
           do isbv1=1,nsymb                                              8d18s22
            isbv2=multh(isbv1,isbv12)                                    8d18s22
            if(isbv2.ge.isbv1)then                                       8d18s22
             if(isbv1.eq.isbv2)then                                      8d18s22
              do i=0,nfdat(3,1,isb)-1                                    8d18s22
               iad=ioff+nvirt(isbv1)*(irm+nroot*i)                       8d18s22
               jad=ioff+nvirt(isbv1)*(jrm+nroot*i)                       8d18s22
               do iv=0,nvirt(isbv1)-1                                   8d19s22
                bc(iad+iv)=bc(iad+iv)-dot*bc(jad+iv)                    8d19s22
               end do                                                    8d18s22
              end do                                                     8d18s22
              ioff=ioff+nvirt(isbv1)*nroot*nfdat(3,1,isb)                8d18s22
              nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                      8d18s22
             else                                                        8d18s22
              nvv=nvirt(isbv1)*nvirt(isbv2)                              8d18s22
             end if                                                      8d18s22
             do l=1,4                                                    8d18s22
              if(nfdat(3,l,isb).gt.0)then
               do i=0,nfdat(3,l,isb)-1                                    8d18s22
                iad=ioff+nvv*(irm+nroot*i)                                8d18s22
                jad=ioff+nvv*(jrm+nroot*i)                                8d18s22
                do ivv=0,nvv-1                                          8d19s22
                 bc(iad+ivv)=bc(iad+ivv)-dot*bc(jad+ivv)                8d19s22
                end do                                                    8d18s22
               end do                                                     8d18s22
               ioff=ioff+nfdat(3,l,isb)*nvv*nroot                        8d19s22
              end if                                                     8d18s22
             end do                                                      8d18s22
            end if                                                       8d18s22
           end do                                                        8d18s22
          end do                                                         8d18s22
         else                                                            8d18s22
          call ucdoubop(nsymb,isymmrci,multh,nvirt,mdon,mdoo,nff2,ncsf, 8d24s22
     $        ncsf2,ihddiag,irm,jrm,dot,2,nroot,bc,ibc)                 11d10s22
         end if                                                          8d18s22
        end if                                                           8d18s22
       end do                                                           8d19s22
       dot=0d0
       if(nsing.ne.0)then                                               8d18s22
        do isb=1,nsymb                                                  8d18s22
         isbv=multh(isb,isymmrci)                                       8d18s22
         do nclop=mdon+1,mdoo+1                                         8d18s22
          if(min(nff1(nclop,isb,1),nvirt(isbv)).gt.0)then               8d18s22
           iarg=nclop-mdon                                              8d18s22
           ncol=ncsf(iarg)*nff1(nclop,isb,1)                            8d18s22
           nrow=nroot*nvirt(isbv)                                       8d18s22
           call ilimts(ncol,1,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e) 8d18s22
           nhere=ih+1-il                                                8d18s22
           ivtmp=ibcoff                                                 8d18s22
           ibcoff=ivtmp+nhere*nrow                                      8d18s22
           call enough('recoverwf.  4',bc,ibc)
           i18=1                                                        8d18s22
           i28=nrow                                                     8d18s22
           i38=il                                                       8d18s22
           i48=ih                                                       8d18s22
           call ddi_get(bc,ibc,ihsdiag(nclop,isb,1),i18,i28,i38,i48,    11d15s22
     $          bc(ivtmp))                                              11d15s22
           do i=0,nhere-1                                               8d18s22
            do iv=0,nvirt(isbv)-1                                       8d18s22
             iad=ivtmp+irm+nroot*(iv+nvirt(isbv)*i)                     8d18s22
             dot=dot+bc(iad)**2                                         8d18s22
            end do                                                      8d18s22
           end do                                                       8d18s22
           ibcoff=ivtmp                                                 8d18s22
          end if
         end do                                                         8d18s22
        end do                                                          8d18s22
       end if                                                           8d18s22
       if(ndoub.gt.0)then                                               8d18s22
        if(iunc(2).eq.0)then                                            8d18s22
         ioff=ivdinout                                                  8d18s22
         do isb=1,nsymb                                                 8d18s22
          isbv12=multh(isb,isymmrci)                                    8d18s22
          do isbv1=1,nsymb                                              8d18s22
           isbv2=multh(isbv1,isbv12)                                    8d18s22
           if(isbv2.ge.isbv1)then                                       8d18s22
            if(isbv1.eq.isbv2)then                                      8d18s22
             do i=0,nfdat(3,1,isb)-1                                    8d18s22
              iad=ioff+nvirt(isbv1)*(irm+nroot*i)                       8d18s22
              do iv=mynowprog,nvirt(isbv1)-1,mynprocg                   8d18s22
               dot=dot+bc(iad+iv)**2                                    8d18s22
              end do                                                    8d18s22
             end do                                                     8d18s22
             ioff=ioff+nvirt(isbv1)*nroot*nfdat(3,1,isb)                8d18s22
             nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                      8d18s22
            else                                                        8d18s22
             nvv=nvirt(isbv1)*nvirt(isbv2)                              8d18s22
            end if                                                      8d18s22
            do l=1,4                                                    8d18s22
             if(nfdat(3,l,isb).gt.0)then
              do i=0,nfdat(3,l,isb)-1                                    8d18s22
               iad=ioff+nvv*(irm+nroot*i)                                8d18s22
               do ivv=mynowprog,nvv-1,mynprocg                           8d18s22
                dot=dot+bc(iad+ivv)**2                                   8d18s22
               end do                                                    8d18s22
              end do                                                     8d18s22
              ioff=ioff+nfdat(3,l,isb)*nvv*nroot                        8d19s22
             end if                                                     8d18s22
            end do                                                      8d18s22
           end if                                                       8d18s22
          end do                                                        8d18s22
         end do                                                         8d18s22
        else                                                            8d18s22
         call ucdoubop(nsymb,isymmrci,multh,nvirt,mdon,mdoo,nff2,ncsf,  8d24s22
     $        ncsf2,ihddiag,irm,irm,dot,1,nroot,bc,ibc)                 11d10s22
        end if                                                          8d18s22
       end if                                                           8d18s22
       call dws_gsumf(dot,1)                                            8d18s22
       xnorm=1d0/sqrt(dot)                                              8d18s22
       if(lpr)write(6,*)('normalizate '),ir,(' using '),dot,xnorm
       if(nsing.ne.0)then                                               8d18s22
        do isb=1,nsymb                                                  8d18s22
         isbv=multh(isb,isymmrci)                                       8d18s22
         do nclop=mdon+1,mdoo+1                                         8d18s22
          if(min(nff1(nclop,isb,1),nvirt(isbv)).gt.0)then               8d18s22
           iarg=nclop-mdon                                              8d18s22
           ncol=ncsf(iarg)*nff1(nclop,isb,1)                            8d18s22
           nrow=nroot*nvirt(isbv)                                       8d18s22
           call ilimts(ncol,1,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e) 8d18s22
           nhere=ih+1-il                                                8d18s22
           ivtmp=ibcoff                                                 8d18s22
           ibcoff=ivtmp+nhere*nrow                                      8d18s22
           call enough('recoverwf.  5',bc,ibc)
           i18=1                                                        8d18s22
           i28=nrow                                                     8d18s22
           i38=il                                                       8d18s22
           i48=ih                                                       8d18s22
           call ddi_get(bc,ibc,ihsdiag(nclop,isb,1),i18,i28,i38,i48,    11d15s22
     $          bc(ivtmp))                                              11d15s22
           do i=0,nhere-1                                               8d18s22
            do iv=0,nvirt(isbv)-1                                       8d18s22
             iad=ivtmp+irm+nroot*(iv+nvirt(isbv)*i)                     8d18s22
             bc(iad)=bc(iad)*xnorm                                      8d18s22
            end do                                                      8d18s22
           end do                                                       8d18s22
           call ddi_put(bc,ibc,ihsdiag(nclop,isb,1),i18,i28,i38,i48,    11d15s22
     $          bc(ivtmp))                                              11d15s22
           ibcoff=ivtmp                                                 8d18s22
          end if
         end do                                                         8d18s22
        end do                                                          8d18s22
       end if                                                           8d18s22
       if(ndoub.gt.0)then                                               8d18s22
        if(iunc(2).eq.0)then                                            8d18s22
         ioff=ivdinout                                                  8d18s22
         do isb=1,nsymb                                                 8d18s22
          isbv12=multh(isb,isymmrci)                                    8d18s22
          do isbv1=1,nsymb                                              8d18s22
           isbv2=multh(isbv1,isbv12)                                    8d18s22
           if(isbv2.ge.isbv1)then                                       8d18s22
            if(isbv1.eq.isbv2)then                                      8d18s22
             do i=0,nfdat(3,1,isb)-1                                    8d18s22
              iad=ioff+nvirt(isbv1)*(irm+nroot*i)                       8d18s22
              do iv=0,nvirt(isbv1)-1                                    8d18s22
               bc(iad+iv)=bc(iad+iv)*xnorm                              8d18s22
              end do                                                    8d18s22
             end do                                                     8d18s22
             ioff=ioff+nvirt(isbv1)*nroot*nfdat(3,1,isb)                8d18s22
             nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                      8d18s22
            else                                                        8d18s22
             nvv=nvirt(isbv1)*nvirt(isbv2)                              8d18s22
            end if                                                      8d18s22
            do l=1,4                                                    8d18s22
             if(nfdat(3,l,isb).gt.0)then
              do i=0,nfdat(3,l,isb)-1                                    8d18s22
               iad=ioff+nvv*(irm+nroot*i)                                8d18s22
               do ivv=0,nvv-1                                           8d18s22
                bc(iad+ivv)=bc(iad+ivv)*xnorm                           8d18s22
               end do                                                    8d18s22
              end do                                                     8d18s22
              ioff=ioff+nfdat(3,l,isb)*nvv*nroot                        8d19s22
             end if                                                     8d18s22
            end do                                                      8d18s22
           end if                                                       8d18s22
          end do                                                        8d18s22
         end do                                                         8d18s22
         call dws_bcast(bc(ivdinout),nall)                              8d18s22
        else                                                            8d18s22
         call ucdoubop(nsymb,isymmrci,multh,nvirt,mdon,mdoo,nff2,ncsf,  8d24s22
     $        ncsf2,ihddiag,irm,irm,xnorm,3,nroot,bc,ibc)               11d10s22
        end if                                                          8d18s22
       end if                                                           8d18s22
      end do                                                            8d18s22
c
c     now for hss diagonals and vs ...
c
      if(lpr)write(6,*)('hss diagonals as vs ...')
      if(nsing.ne.0)then                                                8d19s22
       ioff=1                                                           8d19s22
       ioffh=1                                                          8d19s22
       do isb=1,nsymb                                                   8d19s22
        isbv=multh(isb,isymmrci)                                        8d19s22
        do nclop=mdon+1,mdoo+1                                          8d19s22
         if(min(nff1(nclop,isb,1),nvirt(isbv)).gt.0)then                8d19s22
          iarg=nclop-mdon                                               8d19s22
          ncol=ncsf(iarg)*nff1(nclop,isb,1)                             8d19s22
          nrow=nroot*nvirt(isbv)                                        8d19s22
          call ilimts(ncsf(iarg),nff1(nclop,isb,1),mynprocg,mynowprog,  8d19s22
     $         il,ih,i1s,i1e,i2s,i2e)                                   8d19s22
          nhere=ih+1-il                                                 8d19s22
          if(nhere.gt.0)then                                            8d19s22
           i18=1                                                        8d19s22
           i28=nrow                                                     8d19s22
           i38=il                                                       8d19s22
           i48=ih                                                       8d19s22
           itmp=ibcoff                                                  8d19s22
           ibcoff=itmp+nhere*nrow                                       8d19s22
           call enough('recoverwf.  6',bc,ibc)
           call ddi_get(bc,ibc,ihsdiag(nclop,isb,1),i18,i28,i38,i48,    11d15s22
     $          vs(ioff))                                               11d15s22
           jtmp=itmp                                                    8d19s22
           i10=i1s                                                      8d19s22
           i1n=ncsf(iarg)                                               8d19s22
           do i2=i2s,i2e                                                8d19s22
            if(i2.eq.i2e)i1n=i1e                                        8d19s22
            joffh=ioffh+nvirt(isbv)*(i2-i2s)                            8d19s22
            do i1=i10,i1n                                               8d19s22
             do i=0,nvirt(isbv)-1                                       8d19s22
              do j=0,nroot-1                                            8d19s22
               bc(jtmp+j)=vs(ioff+j)*hss(joffh+i)                       8d19s22
               vso(ioff+j)=vs(ioff+j)                                   8d19s22
              end do                                                    8d19s22
              jtmp=jtmp+nroot                                           8d19s22
              ioff=ioff+nroot                                           8d19s22
             end do                                                     8d19s22
            end do                                                      8d19s22
            i10=1                                                       8d19s22
           end do                                                       8d19s22
           ioffh=ioffh+nvirt(isbv)*(i2e+1-i2s)                          8d19s22
           call ddi_put(bc,ibc,ihsdiag(nclop,isb,2),i18,i28,i38,i48,    11d15s22
     $          bc(itmp))                                               11d15s22
           ibcoff=itmp                                                  8d19s22
          end if                                                        8d19s22
         end if                                                         8d19s22
        end do                                                          8d19s22
       end do                                                           8d19s22
      end if                                                            8d19s22
c
c     now for hdd diagonals if uc
c
      if(lpr)write(6,*)('hdd diagonals? ')
      if(ndoub.gt.0.and.iunc(2).ne.0)then                               8d19s22
       ioffh=1                                                          8d19s22
       do isb=1,nsymb                                                   8d18s22
        isbv12=multh(isb,isymmrci)                                      8d18s22
        nvisv=0                                                            6d9s21
        nvnotv=0                                                           6d9s21
        do isbv1=1,nsymb                                                8d18s22
         isbv2=multh(isbv1,isbv12)                                      8d18s22
         if(isbv2.ge.isbv1)then                                         8d18s22
          if(isbv1.eq.isbv2)then                                        8d18s22
           nvisv=nvisv+nvirt(isbv1)                                     8d19s22
           nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        8d18s22
          else                                                          8d18s22
           nvv=nvirt(isbv1)*nvirt(isbv2)                                8d18s22
          end if                                                        8d18s22
          nvnotv=nvnotv+nvv                                             8d19s22
         end if                                                         8d19s22
        end do                                                          8d19s22
        do nclop=mdon+1,mdoo+1                                          8d19s22
         if(nff2(nclop,isb).gt.0)then                                   8d19s22
          iarg=nclop-mdon                                               8d19s22
          call ilimts(1,nff2(nclop,isb),mynprocg,mynowprog,il,ih,i1s,   8d19s22
     $         i1e,i2s,i2e)                                             8d19s22
          nhere=ih+1-il                                                 8d19s22
          if(nhere.gt.0)then                                            8d19s22
           nrow=nroot*(nvisv*ncsf2(1,iarg)+nvnotv*ncsf(iarg))           8d19s22
           itmp=ibcoff                                                  8d19s22
           ibcoff=itmp+nrow*nhere                                       8d19s22
           i18=1                                                        8d19s22
           i28=nrow                                                     8d19s22
           i38=il                                                       8d19s22
           i48=ih                                                       8d19s22
           call ddi_get(bc,ibc,ihddiag(nclop,isb,1),i18,i28,i38,i48,    11d15s22
     $          bc(itmp))                                               11d15s22
           jtmp=itmp                                                    8d19s22
           do iff=il,ih                                                 8d19s22
            do isbv1=1,nsymb                                            8d19s22
             isbv2=multh(isbv1,isbv12)                                  8d19s22
             if(isbv2.ge.isbv1)then                                     8d19s22
              if(isbv1.eq.isbv2)then                                    8d19s22
               do iv=0,nvirt(isbv1)-1                                   8d19s22
                do i=0,nroot*ncsf2(1,iarg)-1                              8d19s22
                 bc(jtmp+i)=bc(jtmp+i)*hdd(ioffh+iv)                      8d19s22
                end do                                                  8d19s22
                jtmp=jtmp+nroot*ncsf2(1,iarg)                             8d19s22
               end do                                                     8d19s22
               ioffh=ioffh+nvirt(isbv1)                                   8d19s22
               nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                    8d19s22
              else                                                      8d19s22
               nvv=nvirt(isbv1)*nvirt(isbv2)                            8d19s22
              end if                                                    8d19s22
              do ivv=0,nvv-1                                            8d19s22
               do i=0,nroot*ncsf(iarg)-1                                  8d19s22
                bc(jtmp+i)=bc(jtmp+i)*hdd(ioffh+ivv)                      8d19s22
               end do                                                   8d19s22
               jtmp=jtmp+nroot*ncsf(iarg)                                 8d19s22
              end do                                                    8d19s22
              ioffh=ioffh+nvv                                             8d19s22
             end if                                                     8d19s22
            end do                                                      8d19s22
           end do                                                       8d19s22
           call ddi_put(bc,ibc,ihddiag(nclop,isb,2),i18,i28,i38,i48,    11d15s22
     $          bc(itmp))                                               11d15s22
           ibcoff=itmp                                                  8d19s22
          end if                                                        8d19s22
         end if                                                         8d19s22
        end do                                                          8d19s22
       end do                                                           8d19s22
      end if                                                            8d19s22
      if(lpr)write(6,*)('all done in recoverwf!')                       6d23s23
      call sleep(1)                                                     6d26s23
      call dws_synca                                                    6d26s23
      return
      end                                                               8d16s22
