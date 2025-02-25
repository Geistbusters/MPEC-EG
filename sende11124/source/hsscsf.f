c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hsscsf(mode,nhand,nsing,ih0av,nh0av,nff1,iff1,         11d19s20
     $     ncsf,nct1,mdon,mdoo,jmats,kmats,ntot,                        11d19s20
     $     ioooo,multh,nec,shift,nsymbd,nrootu,tdenss,tovr,ixw1,ixw2,   11d26s20
     $     issdig,ldebug,bc,ibc)                                        11d10s22
      implicit real*8 (a-h,o-z)
      dimension ih0av(*),nh0av(*),iff1(*),nff1(mdoo+1,nsymbd),          11d19s20
     $     ncsf(*),jmats(*),kmats(*),ioooo(*),nct1(*),multh(8,8)
      external second                                                   4d23s20
      integer*8 nhand(mdoo+1,nsymbd,2),i1,i12,in,im                     12d7s20
      logical ldebug                                                    1d25s21
      include "common.store"
      include "common.hf"                                               7d15s19
      include "common.mrci"                                             7d15s19
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
c
c     mode = 1 means compute approximate diagonal of ham matrix
c     mode = ne 1 means compute Cs*Hss
c
      i1=1                                                              11d19s20
       iffoff=1                                                         11d19s20
       do isb=1,nsymb                                                   7d15s19
        isbv=multh(isb,isymmrci)                                        7d15s19
        call hsscsf12(nhand(1,isb,2),ih0av,nh0av,iff1,nff1(1,isb),ncsf, 11d26s20
     $       mdon,mdoo,ioooo,jmats,kmats,nec,nvirt(isbv),shift,isbv,    11d19s20
     $       iffoff,bc,ibc)                                             11d10s22
       end do                                                           7d15s19
       if(ldebug)write(6,*)('what we have for ss diagonals')            1d25s21
       value=0d0                                                        7d8s21
       call dws_synca                                                   12d12s20
       issdig=ibcoff                                                    11d26s20
       jssdig=issdig                                                    11d26s20
       do isb=1,nsymb
        if(ldebug)write(6,*)('for symmetry '),isb,ibcoff
        isbv=multh(isb,isymmrci)                                        7d15s19
        do nclop=1,mdoo+1
         if(nff1(nclop,isb).gt.0)then
          nclo=nclop-1                                                  11d19s20
          iarg=nclop-mdon                                               11d10s20
          nall=nff1(nclop,isb)*ncsf(iarg)                               11d10s20
          if(ldebug)write(6,*)('for nclo = '),nclo,nall                                11d19s20
          call ilimts(1,nall,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)  11d26s20
          ic=0                                                          11d10s20
          idn=nall                                                      11d10s20
          idx=0                                                         11d10s20
          do iff=1,nff1(nclop,isb)                                      11d10s20
           do m=1,ncsf(iarg)                                            11d10s20
            ic=ic+1                                                     11d10s20
            if(ic.ge.il.and.ic.le.ih)then                               11d10s20
             idn=min(idn,iff)                                           11d10s20
             idx=max(idx,iff)                                           11d10s20
            end if                                                      11d10s20
           end do                                                       11d10s20
          end do                                                        11d10s20
          if(ldebug)write(6,*)('idn,idx: '),idn,idx                               11d10s20
          nhere=idx+1-idn                                               11d10s20
          if(ldebug)write(6,*)('nhere '),nhere,il,ih
          if(nhere.gt.0)then                                            2d23s21
           ibcoff=jssdig+nvirt(isbv)*nhere                               12d7s20
           call enough('hsscsf.  1',bc,ibc)
          end if
           i12=nvirt(isbv)                                               12d7s20
           im=idn                                                        11d10s20
           in=idx                                                        11d10s20
           call ddi_get(bc,ibc,nhand(nclop,isb,2),i1,i12,im,in,         11d15s22
     $          bc(jssdig))                                             11d15s22
           if(ldebug)then
            call second(tcur)
            write(6,*)('my part: '),tcur
            call prntm2(bc(jssdig),nvirt(isbv),nhere,nvirt(isbv))         12d7s20
           end if                                                        1d25s21
           nbad=0                                                        2d22s21
           do i=0,nvirt(isbv)*nhere-1                                    2d22s21
            if(bc(jssdig+i).eq.0d0)nbad=nbad+1                           2d22s21
           end do                                                        2d22s21
           if(nbad.ne.0)then                                            2d22s21
            write(6,*)('got '),im,in
            call prntm2(bc(jssdig),nvirt(isbv),nhere,nvirt(isbv))
            im=nff1(nclop,isb)
            call ddi_get(bc,ibc,nhand(nclop,isb,2),i1,i12,i1,im,        11d15s22
     $           bc(ibcoff))                                            11d15s22
            write(6,*)('full block ')
            call prntm2(bc(ibcoff),nvirt(isbv),nff1(nclop,isb),
     $         nvirt(isbv))
            value=1d0                                                   2d19s21
            call dws_gsumf(value,1)                                     2d19s21
            call dws_synca
            call dws_finalize
            stop  'my partss'
           end if                                                       2d19s21
           jssdig=ibcoff                                                 11d26s20
         end if
        end do
       end do
       call dws_gsumf(value,1)                                          7d8s21
       if(value.ne.0d0)then                                             7d8s21
        call dws_synca                                                  7d8s21
        call dws_finalize                                               7d8s21
        stop 'ss diagonal error'                                        7d8s21
       end if                                                           7d8s21
      return                                                            7d15s19
      end                                                               7d15s19
