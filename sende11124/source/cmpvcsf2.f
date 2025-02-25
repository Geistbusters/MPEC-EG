c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine cmpvcsf2(ndet,ncsf,icmp1,icmp2,vcmp,vec,nkeep)         11d8s22
      implicit real*8 (a-h,o-z)                                         11d27s19
c                                                                       11d27s19
c     store only nonzero elements of vec.                               11d27s19
c     compared to cmpvcsf, this uses the transpose of vec.              11d8s22
c                                                                       11d27s19
      dimension vec(ncsf,ndet),icmp1(2,ndet),icmp2(*),vcmp(*)           11d8s22
      ii=1                                                              11d27s19
      do idet=1,ndet                                                    11d27s19
       icmp1(1,idet)=ii
       do icsf=1,ncsf
        if(abs(vec(icsf,idet)).gt.1d-12)then                            11d8s22
         vcmp(ii)=vec(icsf,idet)                                        11d8s22
         icmp2(ii)=icsf                                                 11d27s19
         ii=ii+1                                                        11d27s19
        end if                                                          11d27s19
       end do                                                           11d27s19
       icmp1(2,idet)=ii-1                                               11d27s19
      end do                                                            11d27s19
      nkeep=ii-1
      return                                                            11d27s19
      end                                                               11d27s19
