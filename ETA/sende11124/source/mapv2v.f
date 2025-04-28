c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine mapv2v(vecci,nci,nff0,iff0,vint,ni,nroot,nfcn,ibasis,  8d10s22
     $     iptrbit,mdon,mdoo,ncsf,norb,nec,lprint,bc,ibc)               11d14s22
      implicit real*8 (a-h,o-z)                                         8d10s22
      logical lprint                                                    8d11s22
      integer*8 i0c,i0o                                                 8d11s22
      dimension vecci(nci,*),nff0(mdoo+1,3),iff0(*),vint(ni,*),         8d11s22
     $     ibasis(3,*),iptrbit(2,*),ncsf(*),nab4(2,3)                   8d11s22
      include "common.store"                                            8d11s22
      do i=1,nroot                                                      8d11s22
       do j=1,ni                                                        8d11s22
        vint(j,i)=0d0                                                   8d11s22
       end do                                                           8d11s22
      end do                                                            8d11s22
      ioff=0                                                            8d11s22
      do nclo=mdon,mdoo                                                 8d11s22
       nopen=nec-2*nclo                                                 8d11s22
       nclop=nclo+1                                                     8d11s22
       iarg=nclop-mdon                                                  8d11s22
       iff=nff0(nclop,2)                                                8d11s22
       ipb=nff0(nclop,3)-1                                              8d11s22
       do if=1,nff0(nclop,1)                                            8d11s22
        i0c=iff0(iff)                                                   8d11s22
        iff=iff+1                                                       8d11s22
        i0o=iff0(iff)                                                   8d11s22
        iff=iff+1                                                       8d11s22
        joff=0                                                          8d11s22
        do jf=1,nfcn                                                    8d11s22
         ncloj=ibasis(1,jf)                                             8d11s22
         nopenj=nec-2*ncloj                                             8d11s22
         nclojp=ncloj+1                                                 8d11s22
         jarg=nclojp-mdon                                               8d11s22
         if(ncloj.eq.nclo)then                                          8d11s22
          jjc=iptrbit(1,nclojp)+ibasis(2,jf)-1                          8d11s22
          jjo=iptrbit(2,nclojp)+ibasis(3,jf)-1                          8d11s22
          if(ibc(jjc).eq.i0c.and.ibc(jjo).eq.i0o)then                   8d11s22
           do ir=1,nroot                                                8d11s22
            do j=1,ncsf(jarg)                                           8d11s22
             vint(j+joff,ir)=vecci(j+ipb,ir)                            8d11s22
            end do                                                      8d11s22
           end do                                                       8d11s22
           go to 1                                                      8d11s22
          end if                                                        8d11s22
         end if                                                         8d11s22
         joff=joff+ncsf(jarg)                                           8d11s22
        end do                                                          8d11s22
    1   continue                                                        8d11s22
        ipb=ipb+ncsf(iarg)                                              8d11s22
       end do                                                           8d11s22
      end do                                                            8d11s22
      do ir=1,nroot
       sumci=0d0
       sumi=0d0
       do j=1,nci
        sumci=sumci+vecci(j,ir)**2
       end do
       do j=1,ni
        sumi=sumi+vint(j,ir)**2
       end do
       if(lprint)write(6,*)('for root '),ir,('vecci '),sumci,           8d11s22
     $      ('vint '),sumi                                              8d11s22
c                                                                       8d11s22
c     why do I normalize rather than orthogonalize?                     8d11s22
c     Davidson correction with fixed coefficients.                      8d11s22
c     Otherwise we reorthogonalize after projecting out 2es             8d11s22
c                                                                       8d11s22
       xnorm=1d0/sqrt(sumi)                                             8d11s22
       do j=1,ni                                                        8d11s22
        vint(j,ir)=vint(j,ir)*xnorm                                     8d11s22
       end do                                                           8d11s22
      end do
      return
      end
