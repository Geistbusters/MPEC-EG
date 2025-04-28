c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine savetonormal(vec,nps,ips,nroot,nsymb,nadet,nbdet,      7d14s22
     $     nsbeta,vecn)                                                 7d14s22
      implicit real*8 (a-h,o-z)                                         7d14s22
      integer*2 ips                                                     7d14s22
      dimension vec(nps,*),ips(4,*),nadet(*),nbdet(*),nsbeta(*),        7d14s22
     $     ioff(8),vecn(nps,*)                                          7d14s22
      irun=1                                                            7d14s22
      do isb=1,nsymb                                                    7d14s22
       ioff(isb)=irun                                                   7d14s22
       irun=irun+nadet(isb)*nbdet(nsbeta(isb))                          7d14s22
      end do                                                            7d14s22
      do i=1,nps                                                        7d14s22
       ipa=ips(1,i)                                                     7d14s22
       ipb=ips(2,i)                                                     7d14s22
       isa=ips(3,i)                                                     7d14s22
       isb=nsbeta(isa)                                                  7d14s22
       iad=ioff(isa)+ipb-1+nbdet(isb)*(ipa-1)                           7d14s22
       do ir=1,nroot                                                    7d14s22
        vecn(iad,ir)=vec(i,ir)                                               7d14s22
       end do                                                           7d14s22
      end do                                                            7d14s22
      return                                                            7d14s22
      end                                                               7d14s22
