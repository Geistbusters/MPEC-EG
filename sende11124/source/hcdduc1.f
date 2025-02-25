c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcdduc1(xmat,ncsf,ncsf2,iargb,iargk,vd,gd,isbv12,nsymb,6d22s21
     $     multh,nvirt,nrootu,sr2,watch,iflg,name,idncsf2,bc,ibc)       11d10s22
      implicit real*8 (a-h,o-z)                                         6d21s21
      character*(*) name                                                6d24s21
      dimension xmat(*),ncsf(*),ncsf2(idncsf2,*),vd(*),gd(*),multh(8,8),7d26s21
     $     nvirt(*)                                                     6d21s21
      include "common.store"
      watch0=watch
      ixt=ibcoff
      ibcoff=ixt+ncsf(iargb)*ncsf(iargk)
      call enough('hcdduc1.  1',bc,ibc)
      do ib=0,ncsf(iargb)-1
       do ik=0,ncsf(iargk)-1
        ikb=ixt+ik+ncsf(iargk)*ib
        ibk=1+ib+ncsf(iargb)*ik
        bc(ikb)=xmat(ibk)
       end do
      end do
      ik=1                                                              6d21s21
      ib=1                                                              6d22s21
      sz=0d0
      sz2=0d0
      nvisv=0
      nvnotv=0
      nq=0
      do isbv1=1,nsymb                                                  6d21s21
       isbv2=multh(isbv1,isbv12)                                        6d21s21
       if(isbv2.ge.isbv1)then                                           6d21s21
        if(isbv1.eq.isbv2)then                                          6d21s21
         nn=nvirt(isbv1)*nrootu                                         6d21s21
         nvisv=nvisv+nvirt(isbv1)
         if(min(nn,ncsf2(1,iargb),ncsf2(1,iargk)).gt.0)then             7d30s24
          call dgemm('n','n',nn,ncsf2(1,iargb),ncsf2(1,iargk),1d0,       6d25s21
     $        vd(ik),nn,bc(ixt),ncsf(iargk),1d0,gd(ib),nn,              6d25s21
     d' hcdduc1.  1')
         end if                                                         8d16s21
         ik=ik+nn*ncsf2(1,iargk)                                        6d22s21
         ib=ib+nn*ncsf2(1,iargb)                                        6d22s21
         nq=nq+nn*ncsf2(1,iargk)
         nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                          6d21s21
        else                                                            6d21s21
         nvv=nvirt(isbv1)*nvirt(isbv2)                                  6d21s21
        end if                                                          6d21s21
        nn=nvv*nrootu                                                   6d21s21
        nvnotv=nvnotv+nvv
        if(min(nn,ncsf(iargb),ncsf(iargk)).gt.0)then                    7d30s24
         call dgemm('n','n',nn,ncsf(iargb),ncsf(iargk),1d0,              6d22s21
     $       vd(ik),nn,bc(ixt),ncsf(iargk),1d0,gd(ib),nn,               6d25s21
     d' hcdduc1.  2')
        end if                                                          8d19s21
        ik=ik+nn*ncsf(iargk)                                             6d21s21
        nq=nq+nn*ncsf(iargk)
        ib=ib+nn*ncsf(iargb)                                             6d21s21
       end if                                                           6d21s21
      end do                                                            6d21s21
      ibcoff=ixt
      return                                                            6d21s21
      end                                                               6d21s21
