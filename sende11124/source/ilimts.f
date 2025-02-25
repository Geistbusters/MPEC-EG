c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine ilimts(nbasl,nbasr,nproc,ip,ilj,ihj,i1s,i1e,i2s,i2e)
      implicit real*8 (a-h,o-z)
c
c     given col dimensions nbasl,nbasr, compute limits for current
c     processor.
c
c     ilimcode=1 means evenly distribute i1 and i2 across procs.
c     ilimcode=2 means i1 is either full dimension or we only have
c                one i2 per proc.
      data icall/0/
      common/ilimcm/ilimcode                                            5d7s19
      save                                                              5d7s19
      nbas2=nbasl*nbasr
      if(nbas2.le.0)then                                                8d2s12
       ilj=1                                                            8d2s12
       ihj=0                                                            8d2s12
       i1s=1                                                            8d2s12
       i1e=0                                                            8d2s12
       i2s=1                                                            8d2s12
       i2e=0                                                            8d2s12
       return                                                           8d2s12
      end if                                                            8d2s12
      nmin=nbas2/nproc
      nminp=nmin+1
      nminm=nmin-1
      nleft=nbas2-nproc*nmin
      if(ip.le.nleft-1)then
       ilj=max(ip*nminp,0)+1
       ihj=ilj+nmin
      else
       ilj=nleft*nminp+(ip-nleft)*nmin+1
       ihj=min(ilj+nminm,nbas2)
      end if
      i1s=1                                                             8d7s7
      i1e=0                                                             8d7s7
      i2s=1                                                             8d7s7
      i2e=0                                                             8d7s7
      if(nbasl.eq.0)return                                              8d7s7
      i2s=ilj/nbasl
      if(i2s*nbasl.eq.ilj)i2s=i2s-1
      i1s=ilj-i2s*nbasl
      i2s=i2s+1
      i2e=ihj/nbasl
      if(i2e*nbasl.eq.ihj)i2e=i2e-1
      i1e=ihj-i2e*nbasl
      i2e=i2e+1
      return
      end
