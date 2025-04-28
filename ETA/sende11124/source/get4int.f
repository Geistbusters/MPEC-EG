c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      real*8 function get4int(i4x,isa,isb,isc,isd,iga,igb,igc,igd,      7d6s21
     $      nvirt,iad,bc,ibc)                                           11d10s22
      implicit real*8 (a-h,o-z)                                         7d6s21
      include "common.store"                                            7d6s21
      common/fnd4xcm/inv4x(2,8,8,8)                                     7d4s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      dimension i4x(*),nvirt(*)                                         7d6s21
      i2eu=inv4x(1,isa,isb,isc)                                         7d6s21
      icase=inv4x(2,isa,isb,isc)                                        7d6s21
      if(icase.eq.1)then                                                7d6s21
       is1=isa                                                          7d6s21
       is2=isb                                                          7d6s21
       is3=isc                                                          7d6s21
       is4=isd                                                          7d6s21
       ig1=iga                                                          7d6s21
       ig2=igb                                                          7d6s21
       ig3=igc                                                          7d6s21
       ig4=igd                                                          7d6s21
      else if(icase.eq.2)then                                           7d6s21
       is1=isb                                                          7d6s21
       is2=isa                                                          7d6s21
       is3=isc                                                          7d6s21
       is4=isd                                                          7d6s21
       ig1=igb                                                          7d6s21
       ig2=iga                                                          7d6s21
       ig3=igc                                                          7d6s21
       ig4=igd                                                          7d6s21
      else if(icase.eq.3)then                                           7d6s21
       is1=isa                                                          7d6s21
       is2=isb                                                          7d6s21
       is3=isd                                                          7d6s21
       is4=isc                                                          7d6s21
       ig1=iga                                                          7d6s21
       ig2=igb                                                          7d6s21
       ig3=igd                                                          7d6s21
       ig4=igc                                                          7d6s21
      else if(icase.eq.4)then                                           7d6s21
       is1=isb                                                          7d6s21
       is2=isa                                                          7d6s21
       is3=isd                                                          7d6s21
       is4=isc                                                          7d6s21
       ig1=igb                                                          7d6s21
       ig2=iga                                                          7d6s21
       ig3=igd                                                          7d6s21
       ig4=igc                                                          7d6s21
      else if(icase.eq.5)then                                           7d6s21
       is1=isc                                                          7d6s21
       is2=isd                                                          7d6s21
       is3=isa                                                          7d6s21
       is4=isb                                                          7d6s21
       ig1=igc                                                          7d6s21
       ig2=igd                                                          7d6s21
       ig3=iga                                                          7d6s21
       ig4=igb                                                          7d6s21
      else if(icase.eq.6)then                                           7d6s21
       is4=isa                                                          7d6s21
       is3=isb                                                          7d6s21
       is1=isc                                                          7d6s21
       is2=isd                                                          7d6s21
       ig4=iga                                                          7d6s21
       ig3=igb                                                          7d6s21
       ig1=igc                                                          7d6s21
       ig2=igd                                                          7d6s21
      else if(icase.eq.7)then                                           7d6s21
       is3=isa                                                          7d6s21
       is4=isb                                                          7d6s21
       is2=isc                                                          7d6s21
       is1=isd                                                          7d6s21
       ig3=iga                                                          7d6s21
       ig4=igb                                                          7d6s21
       ig2=igc                                                          7d6s21
       ig1=igd                                                          7d6s21
      else
       is1=isd                                                          7d6s21
       is2=isc                                                          7d6s21
       is3=isb                                                          7d6s21
       is4=isa                                                          7d6s21
       ig1=igd                                                          7d6s21
       ig2=igc                                                          7d6s21
       ig3=igb                                                          7d6s21
       ig4=iga                                                          7d6s21
      end if
      if(is1.eq.is2)then                                                7d6s21
       ix=max(ig1,ig2)                                                  7d6s21
       in=min(ig1,ig2)                                                  7d6s21
       irow=((ix*(ix+1))/2)+in                                          7d6s21
       nrow=(nvirt(is1)*(nvirt(is1)+1))/2                               7d6s21
      else                                                              7d6s21
       irow=ig1+nvirt(is1)*ig2                                          7d6s21
       nrow=nvirt(is1)*nvirt(is2)                                       7d6s21
      end if                                                            7d6s21
      if(is3.eq.is4)then                                                7d6s21
       ix=max(ig3,ig4)                                                  7d6s21
       in=min(ig3,ig4)                                                  7d6s21
       icol=((ix*(ix+1))/2)+in                                          7d6s21
       ncol=(nvirt(is3)*(nvirt(is3)+1))/2                               7d6s21
      else                                                              7d6s21
       icol=ig3+nvirt(is3)*ig4                                          7d6s21
       ncol=nvirt(is3)*nvirt(is4)                                       7d6s21
      end if                                                            7d6s21
      if(is1.eq.is3.and.is2.eq.is4)then                                 7d6s21
       ix=max(irow,icol)                                                7d6s21
       in=min(irow,icol)                                                7d6s21
       iad=((ix*(ix+1))/2)+in                                           7d6s21
       ntot=(nrow*(nrow+1))/2                                           7d6s21
      else                                                              7d6s21
       iad=irow+nrow*icol                                               7d6s21
       ntot=nrow*ncol                                                   7d6s21
       itry=inv4x(2,is1,is2,is3)
       if(itry.ne.1)then
        write(6,*)('bad syms!!! '),itry
        write(6,*)icase,i2eu,is1,is2,is3,is4,ig1,ig2,ig3,ig4,irow,nrow,
     $      icol,isa,isb,isc,isd
       end if
      end if                                                            7d6s21
      call ilimts(1,ntot,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)      7d6s21
      iad=iad+1                                                         7d6s21
      iadu=0
      if(iad.ge.il.and.iad.le.ih)then                                   7d6s21
       get4int=bc(i4x(i2eu)+iad-il)                                     7d6s21
       iadu=i4x(i2eu)+iad-il
      else                                                              7d6s21
       get4int=0d0                                                      7d6s21
       iadu=-1                                                          1d12s24
      end if                                                            7d6s21
      ibc(ibcoff)=iadu
      ibc(ibcoff+1)=i4x(i2eu)
      ibc(ibcoff+2)=i2eu
      ibc(ibcoff+3)=icase
      ibc(ibcoff+4)=iad
      ibc(ibcoff+5)=ntot
      ibc(ibcoff+6)=irow
      ibc(ibcoff+7)=icol
      ibc(ibcoff+8)=nrow
      return                                                            7d6s21
      end                                                               7d6s21
