c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hddcsf(nhand,ndoub,ih0av,nh0av,nff2,iff2,              4d9s21
     $     mdon,mdoo,jmats,kmats,                                       5d28s21
     $     ioooo,multh,nec,shift,nsymbd,idddig,ldebug,bc,ibc)           11d10s22
      implicit real*8 (a-h,o-z)
      dimension ih0av(*),nh0av(*),iff2(*),nff2(mdoo+1,nsymbd),          4d9s21
     $     jmats(*),kmats(*),ioooo(*),multh(8,8)                        5d28s21
      external second                                                   4d23s20
      integer*8 nhand(mdoo+1,nsymbd,2),i1,i12,in,im,i1off               6d8s21
      logical ldebug                                                    1d25s21
      include "common.store"
      include "common.hf"                                               7d15s19
      include "common.mrci"                                             7d15s19
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      call hddcsf12b(nhand(1,1,2),ih0av,nh0av,iff2,nff2,mdon,mdoo,ioooo,6d13s21
     $     jmats,kmats,nec,nvirt,shift,nsymb,multh,bc,ibc)              11d10s22
      if(ldebug)write(6,*)('what we have for dd diagonals')             5d28s21
      call dws_synca                                                    5d28s21
      idddig=ibcoff                                                     5d28s21
      jdddig=idddig                                                     5d28s21
      do isb=1,nsymb
       isbv12=multh(isb,isymmrci)                                       5d28s21
       ioffv=0                                                          6d8s21
       nall=0                                                           6d13s21
       do isbv1=1,nsymb                                                 5d28s21
        isbv2=multh(isbv1,isbv12)                                       5d28s21
        if(isbv2.ge.isbv1)then                                          5d28s21
         if(isbv1.eq.isbv2)then                                         5d28s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)+1))/2                         5d28s21
         else                                                           5d28s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 5d28s21
         end if                                                         5d28s21
         nall=nall+nvv                                                  6d13s21
        end if                                                          6d13s21
       end do                                                           6d13s21
       do nclop=1,mdoo+1
        if(nff2(nclop,isb).gt.0)then                                    5d28s21
         nclo=nclop-1                                                   11d19s20
         iarg=nclop-mdon                                                11d10s20
         call ilimts(1,nff2(nclop,isb),mynprocg,mynowprog,il,ih,i1s,i1e,6d13s21
     $        i2s,i2e)                                                  6d13s21
         idn=il                                                         6d13s21
         idx=ih                                                         6d13s21
         if(ldebug)write(6,*)('idn,idx: '),idn,idx                               11d10s20
         nhere=idx+1-idn                                                11d10s20
         if(ldebug)write(6,*)('nhere '),nhere,il,ih
         if(nhere.gt.0)then                                             2d23s21
          ibcoff=jdddig+nall*nhere                                      5d28s21
          call enough('hddcsf.  1',bc,ibc)
         end if
         i1off=1                                                        6d13s21
         i12=nall                                                       6d13s21
         im=idn                                                         11d10s20
         in=idx                                                         11d10s20
         call ddi_get(bc,ibc,nhand(nclop,isb,2),i1off,i12,im,in,        11d15s22
     $          bc(jdddig))                                             5d28s21
         if(ldebug)then
          call second(tcur)
          write(6,*)('my part: '),i1off,i12,jdddig-idddig
          call prntm2(bc(jdddig),nall,nhere,nall)                       5d28s21
         end if                                                         1d25s21
         nbad=0                                                         2d22s21
         do i=0,nvv*nhere-1                                             5d28s21
          if(bc(jdddig+i).eq.0d0)nbad=nbad+1                            2d22s21
         end do                                                         2d22s21
         if(nbad.ne.0)then                                              2d22s21
          write(6,*)('got '),im,in
          call prntm2(bc(jdddig),nall,nhere,nall)
          im=nff2(nclop,isb)
          call ddi_get(bc,ibc,nhand(nclop,isb,2),i1off,i12,i1,im,       11d15s22
     $           bc(ibcoff))
          write(6,*)('full block ')
          call prntm2(bc(ibcoff),nall,nff2(nclop,isb),nall)             5d28s21
          value=1d0                                                     2d19s21
          call dws_gsumf(value,1)                                       2d19s21
          call dws_synca
          call dws_finalize
          stop  'my partdd'
         end if                                                         2d19s21
         jdddig=ibcoff                                                  11d26s20
        end if
       end do
      end do
      value=0d0                                                         5d28s21
      call dws_gsumf(value,1)                                           5d28s21
      if(value.ne.0d0)then                                              5d28s21
       call dws_synca
       call dws_finalize
       stop
      end if
      return                                                            7d15s19
      end                                                               7d15s19
