c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine loadsingb(ihsdiag,nff1,nsymb,mdon,mdoo,nvirt,multh,    7d8s21
     $     isymmrci,ncsf,nroot,mynowprog,maxbx,inorm,bc,ibc,mddilow,    12d13s22
     $     mddihig)                                                     12d13s22
      implicit real*8 (a-h,o-z)                                         7d8s21
c
c     inorm =
c     0 read from unit 1, but skip over stuff. This for contraction
c       functions.                                                      8d18s22
c     1 read from unit 1, create dms, and store to them. The is for     8d18s22
c       prop.                                                           8d18s22
c     2 read from unit 2 and save to previously created dms. This is    8d18s22
c       for mrci restarts.                                              8d18s22
c                                                                       8d18s22
      integer*8 ihsdiag(mdoo+1,nsymb),i18,i28,i38                       8d21s21
      dimension nff1(mdoo+1,nsymb),nvirt(*),multh(8,8),ncsf(*)          7d8s21
      include "common.store"                                            7d8s21
      data icall/0/                                                     2d14s22
      save icall                                                        2d14s22
      nsing=0
      i18=1                                                             7d8s21
      do isb=1,nsymb                                                    7d8s21
       isbv=multh(isb,isymmrci)                                         7d8s21
       do ii=mdon+1,mdoo+1                                              7d8s21
        if(inorm.eq.1)ihsdiag(ii,isb)=-1                                8d18s22
        if(nff1(ii,isb).gt.0)then                                       7d8s21
         nrow=nroot*nvirt(isbv)                                         7d8s21
         iarg=ii-mdon
         ncol=nff1(ii,isb)*ncsf(iarg)                                   7d8s21
         i28=nrow                                                       7d8s21
         i38=ncol                                                       7d8s21
         maxbx=max(maxbx,nrow*ncol)                                     7d9s21
         if(inorm.eq.1)then                                             8d18s22
          call ddi_create(bc,ibc,i28,i38,ihsdiag(ii,isb))               11d15s22
          if(ihsdiag(ii,isb).lt.mddilow)mddilow=ihsdiag(ii,isb)         12d13s22
          if(ihsdiag(ii,isb).gt.mddihig)mddihig=ihsdiag(ii,isb)         12d13s22
          call ddi_zero(bc,ibc,ihsdiag(ii,isb))                         11d15s22
          call dws_synca                                                 7d9s21
         end if                                                         8d10s22
         if(mynowprog.eq.0)then                                         7d8s21
          itmp=ibcoff                                                   7d8s21
          ibcoff=itmp+nrow*ncol                                         7d8s21
          call enough('loadsingb.  1',bc,ibc)
          if(inorm.eq.0)then                                            8d10s22
           read(1)dum                                                   8d10s22
          else                                                          8d18s22
           if(inorm.eq.1)then                                           8d18s22
            read(1)(bc(itmp+i),i=0,nrow*ncol-1)                           7d8s21
           else                                                         8d18s22
            read(2)(bc(itmp+i),i=0,nrow*ncol-1)                         8d18s22
           end if                                                       8d18s22
           call ddi_put(bc,ibc,ihsdiag(ii,isb),i18,i28,i18,i38,bc(itmp))11d15s22
          end if                                                        8d10s22
          ibcoff=itmp                                                   7d8s21
         end if                                                         7d8s21
         if(inorm.ne.0)then                                             6d26s23
          call sleep(1)                                                 6d26s23
          call dws_synca                                                6d26s23
         end if                                                         6d26s23
         nsing=nsing+ncsf(iarg)*nvirt(isbv)*nff1(ii,isb)                 7d8s21
        end if                                                          7d8s21
       end do                                                           7d8s21
      end do                                                            7d8s21
      return                                                            7d8s21
      end                                                               7d8s21
