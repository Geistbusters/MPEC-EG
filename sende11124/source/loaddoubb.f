c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine loaddoubb(ihddiag,nff2,nsymb,mdon,mdoo,nvirt,multh,    7d8s21
     $     isymmrci,ncsf,ncsf2,nroot,mynowprog,maxbxd,inorm,bc,ibc,     12d13s22
     $     mddilow,mddihig)                                             12d13s22
      implicit real*8 (a-h,o-z)                                         7d8s21
      integer*8 ihddiag(mdoo+1,nsymb),i18,i28,i38                       8d21s21
      dimension nff2(mdoo+1,nsymb),nvirt(*),multh(8,8),ncsf(*),         2d14s22
     $     ncsf2(*)                                                     2d14s22
      include "common.store"                                            7d8s21
      ndoub=0
      i18=1                                                             7d8s21
      do isb=1,nsymb                                                    7d8s21
       isbv12=multh(isb,isymmrci)                                       7d21s21
       nvisv=0                                                          7d21s21
       nvnotv=0                                                         7d21s21
       do isbv1=1,nsymb                                                 7d21s21
        isbv2=multh(isbv1,isbv12)                                       7d21s21
        if(isbv2.ge.isbv1)then                                          7d21s21
         if(isbv1.eq.isbv2)then                                         7d21s21
          nvisv=nvisv+nvirt(isbv1)                                      7d21s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         7d21s21
         else                                                           7d21s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 7d21s21
         end if                                                         7d21s21
         nvnotv=nvnotv+nvv                                              7d21s21
        end if                                                          7d21s21
       end do                                                           7d21s21
       do ii=mdon+1,mdoo+1                                              7d8s21
        if(inorm.eq.1)ihddiag(ii,isb)=-1                                8d19s22
        if(nff2(ii,isb).gt.0)then                                       7d8s21
         iarg=ii-mdon
         nrow=nroot*(ncsf2(iarg)*nvisv+ncsf(iarg)*nvnotv)               2d14s22
         ncol=nff2(ii,isb)                                              7d21s21
         i28=nrow                                                       7d8s21
         i38=ncol                                                       7d8s21
         maxbxd=max(maxbxd,nrow*ncol)                                     7d9s21
         if(inorm.eq.1)then                                             8d19s22
          call ddi_create(bc,ibc,i28,i38,ihddiag(ii,isb))               11d15s22
          if(ihddiag(ii,isb).lt.mddilow)mddilow=ihddiag(ii,isb)         12d13s22
          if(ihddiag(ii,isb).gt.mddihig)mddihig=ihddiag(ii,isb)         12d13s22
          call ddi_zero(bc,ibc,ihddiag(ii,isb))                         11d15s22
         end if                                                         8d19s22
         call dws_synca                                                 7d9s21
         if(mynowprog.eq.0)then                                         7d8s21
          itmp=ibcoff                                                   7d8s21
          ibcoff=itmp+nrow*ncol                                         7d8s21
          call enough('loaddoubb.  1',bc,ibc)
          read(inorm)(bc(itmp+i),i=0,nrow*ncol-1)                       8d19s22
          call ddi_put(bc,ibc,ihddiag(ii,isb),i18,i28,i18,i38,bc(itmp)) 11d15s22
          ibcoff=itmp                                                   7d8s21
         end if                                                         7d8s21
         call sleep(1)
         call dws_synca                                                 7d9s21
         ndoub=ndoub+(ncsf2(iarg)*nvisv                                 2d14s22
     $        +ncsf(iarg)*nvnotv)*nff2(ii,isb)                          2d14s22
        end if                                                          7d8s21
       end do                                                           7d8s21
      end do                                                            7d8s21
      return                                                            7d8s21
      end                                                               7d8s21
