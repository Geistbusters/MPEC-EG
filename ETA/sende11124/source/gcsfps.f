c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine gcsfps(nfcn,hdiag,hdps,ecut,ibasis,ipsbase,ncsf,       7d11s19
     $     ipointf,mdon,nvv,hivs,hiv,ncsft,nlzzu,zzdig,zzpsdig,ipointp, 5d14s21
     $     ndigs,npsf)                                                  5d14s21
      implicit real*8 (a-h,o-z)                                         6d11s19
      dimension hdiag(*),hdps(*),ibasis(3,*),ipsbase(3,*),ncsf(*),      7d11s19
     $     ipointf(*),hivs(*),hiv(*),zzdig(nfcn,*),zzpsdig(npsf,*),     5d14s21
     $     ipointp(*)                                                   5d14s21
      ips=1                                                             6d11s19
      ivv=1                                                             9d4s19
      ivvs=1                                                            9d4s19
      irun=0                                                            7d11s19
      iprun=0                                                           7d11s19
      do i=1,nfcn                                                       6d11s19
       edel=hdiag(i)-ecut                                               6d11s19
       nclo=ibasis(1,i)                                                 7d11s19
       nclop=nclo+1
       iarg=nclop-mdon                                                  7d11s19
       if(edel.le.0d0)then                                              6d11s19
        do j=1,ncsf(iarg)                                               7d11s19
         irun=irun+1                                                    7d11s19
         iprun=iprun+1                                                  7d11s19
         ipointf(iprun)=irun                                            7d11s19
        end do                                                          7d11s19
        ipointp(ips)=i                                                  3d24s21
        do j=1,3                                                        6d11s19
         ipsbase(j,ips)=ibasis(j,i)                                     6d11s19
        end do                                                          6d11s19
        hdps(ips)=hdiag(i)                                              6d11s19
        do izz=1,ndigs                                                  5d14s21
         zzpsdig(ips,izz)=zzdig(i,izz)                                  5d14s21
        end do                                                          5d14s21
        do j=0,nvv-1                                                    9d4s19
         do k=0,ncsf(iarg)-1                                            9d4s19
          iadf=ivv+k+ncsft*j                                             9d4s19
          iadt=ivvs+j+nvv*k                                             9d4s19
          hivs(iadt)=hiv(iadf)                                          9d4s19
         end do                                                         9d4s19
        end do                                                          9d4s19
        ivvs=ivvs+nvv*ncsf(iarg)                                        9d4s19
        ips=ips+1                                                       6d11s19
       else                                                             7d11s19
        irun=irun+ncsf(iarg)                                            7d11s19
       end if                                                           6d11s19
       ivv=ivv+ncsf(iarg)                                               9d4s19
      end do                                                            6d11s19
      return                                                            6d11s19
      end                                                               6d11s19
