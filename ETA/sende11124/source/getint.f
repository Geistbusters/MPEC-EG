c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      real*8 function getint(i2e,isa,isb,isc,isd,iga,igb,igc,igd,bc,ibc)11d10s22
      implicit real*8 (a-h,o-z)                                         7d7s19
      dimension i2e(*)                                                  7d7s19
      include "common.mrci"                                             7d7s19
      include "common.store"                                            7d7s19
      i2eu=ifind2(isa,isb,isc,isd,icase)                                7d7s19
      if(isa.eq.isb)then                                                7d7s19
       nrow=(irefo(isa)*(irefo(isa)+1))/2                               7d7s19
       ix=max(iga,igb)                                                  7d7s19
       in=min(iga,igb)                                                  7d7s19
       irow=((ix*(ix+1))/2)+in                                          7d7s19
       ix=max(igc,igd)                                                  7d7s19
       in=min(igc,igd)                                                  7d7s19
       icol=((ix*(ix+1))/2)+in                                          7d7s19
       iad=i2e(i2eu)+irow+nrow*icol                                     7d7s19
      else if(icase.eq.1)then
       iad=i2e(i2eu)+iga+irefo(isa)*(igb+irefo(isb)*                    7d7s19
     $      (igc+irefo(isc)*igd))                                        7d7s19
      else if(icase.eq.2)then
       iad=i2e(i2eu)+igb+irefo(isb)*(iga+irefo(isa)*                    7d7s19
     $      (igd+irefo(isd)*igc))                                        7d7s19
      else if(icase.eq.3)then
       iad=i2e(i2eu)+iga+irefo(isa)*(igb+irefo(isb)*                    7d7s19
     $      (igd+irefo(isd)*igc))                                        7d7s19
      else if(icase.eq.4)then
       iad=i2e(i2eu)+igb+irefo(isb)*(iga+irefo(isa)*                    7d7s19
     $      (igc+irefo(isc)*igd))                                        7d7s19
      end if                                                            7d7s19
      getint=bc(iad)                                                    7d7s19
      return                                                            7d7s19
      end                                                               7d7s19
