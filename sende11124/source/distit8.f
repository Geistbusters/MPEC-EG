c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine distit8(ntot,nhere,istart)
      implicit integer*8 (a-z)
      integer dws_me,dws_np
      common/ddidwscm/dws_me,dws_np                                     8d12s08
c
c     take array of dimension ntot and distribute it across procs.
c     the current proc will have nhere entries, starting at istart
c
      nall=ntot/dws_np
      nleft=ntot-nall*dws_np
      if(dws_me.le.nleft-1)then
       nhere=nall+1
       istart=dws_me*nhere+1
      else
       nhere=nall
       istart=nleft*(nhere+1)+(dws_me-nleft)*nhere+1
      end if
      iend=istart+nhere-1
      return
      end
