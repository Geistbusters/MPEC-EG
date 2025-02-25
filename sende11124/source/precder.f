c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine precder(idorth,idcont,nsymb,nfdat,nroot,mdon,mdoo,     9d6s23
     $     ncsf2,nff2,nfdatd,bc,ibc,nxs)                                10d16s23
      implicit real*8 (a-h,o-z)
      include "common.store"                                            9d6s23
      dimension idorth(4,*),idcont(*),nfdat(5,4,*),ncsf2(4,*),          9d11s23
     $     nff2(mdoo+1,nsymb,*),nfdatd(8,*),nxs(*)                      10d16s23
      ibcoffo=ibcoff                                                    9d6s23
      do isb=1,nsymb                                                    9d6s23
       do l=1,4                                                         9d6s23
        idorth(l,isb)=ibcoff                                            9d6s23
        ibcoff=idorth(l,isb)+nfdat(2,l,isb)*nfdat(3,l,isb)*nroot        9d6s23
       end do                                                           9d6s23
       lvcv=nfdat(5,1,isb)                                              9d11s23
       nxx=0
       do nclo2=mdon,mdoo                                               9d11s23
        nclo2p=nclo2+1                                                  9d11s23
        iarg2=nclo2p-mdon                                               9d11s23
        nfdatd(isb,nclo2p)=nxx                                          9d18s23
        do if=0,nff2(nclo2p,isb,7)-1
         ipack8=ibc(lvcv)
         nspace=ibc(lvcv+1)
         do l=1,4
          nxx=nxx+ibc(lvcv+1+l)*ncsf2(l,iarg2)
         end do
         lvcv=lvcv+nspace                                               9d11s23
        end do                                                          9d11s23
       end do
       idcont(isb)=ibcoff                                               9d11s23
       ibcoff=idcont(isb)+nxx*nroot                                     9d11s23
       nxs(isb)=nxx*nroot                                               10d16s23
      end do                                                            9d6s23
      call enough('precder',bc,ibc)                                     9d6s23
      do iz=ibcoffo,ibcoff-1                                            9d6s23
       bc(iz)=0d0                                                       9d6s23
      end do                                                            9d6s23
      return                                                            9d6s23
      end                                                               9d6s23
