c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine maptoold(mdon,mdoo,nff0,iff0,nfcn,ibasis,idorbf,isorbf,1d21s21
     $     iptr,norb,nec,bc,ibc)                                        11d10s22
      implicit real*8 (a-h,o-z)                                         1d21s21
c
c     take new fangled indices nff0 and iff0 and turn them into the
c     old representation using by intcsf
c
      integer*8 itest                                                   1d21s21
      integer*1 iorb(64,2),idorbf(*),isorbf(*)                          1d21s21
      dimension nff0(*),iff0(*),ibasis(3,*),iptr(4,*)                   1d21s21
      include "common.store"
      ioff=1                                                            1d21s21
      noff=1
      do nclo=mdon,mdoo
       nclop=nclo+1
       iptr(1,nclop)=0                                                  1d21s21
       iptr(3,nclop)=0                                                  1d21s21
       if(nff0(nclop).gt.0)then
        do if=1,nff0(nclop)
         io=0                                                           1d21s21
         ic=0
         ioffp=ioff+1                                                   1d21s21
         do i=1,norb
          if(btest(iff0(ioff),i))then                                   1d21s21
           ic=ic+1                                                      1d21s21
           iorb(ic,1)=i                                                 1d21s21
          end if
          if(btest(iff0(ioffp),i))then                                  1d21s21
           io=io+1                                                       1d21s21
           iorb(io,2)=i
          end if                                                        1d21s21
         end do                                                         1d21s21
         noff=noff+1                                                    1d21s21
         ioff=ioff+2                                                    1d21s21
        end do                                                          1d21s21
       end if
      end do
      idorb=ibcoff                                                      1d21s21
      isorb=idorb+nfcn                                                  1d21s21
      ibcoff=isorb+nfcn                                                 1d21s21
      call enough('maptoold.  1',bc,ibc)
      jdorb=idorb                                                       1d21s21
      jsorb=isorb                                                       1d21s21
      kdorb=idorb                                                       1d21s21
      ksorb=isorb                                                       1d21s21
      ks=1                                                              1d21s21
      kd=1                                                              1d21s21
      ioff=1                                                            1d21s21
      noff=0
      do nclo=mdon,mdoo
       nclop=nclo+1
       if(nff0(nclop).gt.0)then
        kdorb0=kdorb                                                      1d21s21
        ksorb0=ksorb                                                    1d21s21
        nopen=nec-2*nclo                                                1d21s21
        jdorb0=jdorb                                                    1d21s21
        jsorb0=jsorb                                                    1d21s21
        iptr(1,nclop)=0                                                 1d21s21
        iptr(2,nclop)=kd                                                1d21s21
        iptr(3,nclop)=0                                                 1d21s21
        iptr(4,nclop)=ks                                                1d21s21
        do if=1,nff0(nclop)
         ibasis(1,if+noff)=nclo                                         1d21s21
         ibc(jdorb)=iff0(ioff)                                          1d21s21
         jdorb=jdorb+1                                                  1d21s21
         ioff=ioff+1                                                    1d21s21
         ibc(jsorb)=iff0(ioff)                                          1d21s21
         jsorb=jsorb+1                                                  1d21s21
         ioff=ioff+1                                                    1d21s21
        end do                                                          1d21s21
        isortd=ibcoff                                                   1d21s21
        isorts=isortd+nff0(nclop)
        ibcoff=isorts+nff0(nclop)
        call enough('maptoold.  2',bc,ibc)
        call idsortdws(ibc(jdorb0),ibc(isortd),nff0(nclop))
        call idsortdws(ibc(jsorb0),ibc(isorts),nff0(nclop))             1d18s23
        is=0                                                            1d21s21
    1   continue                                                        1d21s21
         itest=ibc(jdorb0+is)                                           1d21s21
         ii=ibc(isortd+is)+noff                                         1d21s21
         do i=1,norb
          if(btest(itest,i))then
           idorbf(kd)=i                                                 1d21s21
           kd=kd+1
          end if
         end do
         bc(kdorb)=itest                                                1d21s21
         kdorb=kdorb+1                                                  1d21s21
         nu=kdorb-kdorb0                                                1d21s21
         ic=iptr(2,nclop)+nclo*(nu-1)                                   1d21s21
         ibasis(2,ii)=nu                                                1d21s21
         isp=is+1                                                           1d21s21
         do i=isp,nff0(nclop)-1                                          1d21s21
          if(ibc(jdorb0+i).ne.itest)then                                 1d21s21
           is=i
           go to 1                                                       1d21s21
          else                                                           1d21s21
           ii=ibc(isortd+i)+noff                                        1d21s21
           ibasis(2,ii)=nu                                              1d21s21
          end if                                                         1d21s21
         end do                                                          1d21s21
        iptr(1,nclop)=nu                                                1d21s21
        is=0                                                            1d21s21
    2   continue                                                        1d21s21
         itest=ibc(jsorb0+is)                                           1d21s21
         ii=ibc(isorts+is)+noff                                         1d21s21
         do i=1,norb
          if(btest(itest,i))then
           isorbf(ks)=i                                                 1d21s21
           ks=ks+1                                                      1d21s21
          end if
         end do
         bc(ksorb)=itest                                                1d21s21
         ksorb=ksorb+1                                                  1d21s21
         nu=ksorb-ksorb0                                                1d21s21
         io=iptr(4,nclop)+nopen*(nu-1)                                  1d21s21
         ibasis(3,ii)=nu                                                1d21s21
         isp=is+1                                                           1d21s21
         do i=isp,nff0(nclop)-1                                          1d21s21
          if(ibc(jsorb0+i).ne.itest)then                                 1d21s21
           is=i
           go to 2                                                       1d21s21
          else                                                           1d21s21
           ii=ibc(isorts+i)+noff                                        1d21s21
           ibasis(3,ii)=nu                                              1d21s21
          end if                                                         1d21s21
         end do                                                          1d21s21
        iptr(3,nclop)=nu                                                1d21s21
        ibcoff=isortd
        noff=noff+nff0(nclop)                                           1d21s21
       end if
      end do
      ibcoff=idorb                                                      1d21s21
      return
      end
