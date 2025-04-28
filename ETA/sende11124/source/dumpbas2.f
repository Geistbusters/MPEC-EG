c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dumpbas2(nec,mdon,mdoo,nff0,iff0,ncsf)                 7d19s21
      dimension nff0(mdoo+1,3),iff0(*),ncsf(*)                          7d19s21
      mff0=0                                                            7d8s21
      ncc=0                                                             7d19s21
      do ii=mdon+1,mdoo+1                                               7d8s21
       iarg=ii-mdon                                                     7d19s21
       nclo=ii-1
       nopen=nec-nclo*2
       mff0=mff0+nff0(ii,1)*2                                           7d19s21
       ncc=ncc+ncsf(iarg)*nff0(ii,1)                                    7d19s21
      end do                                                            7d8s21
      mdoop=mdoo+1                                                      7d8s21
      write(2)nec,mff0,mdon,mdoop,ncc                                   7d19s21
      write(2)((nff0(j,i),i=1,3),j=1,mdoop)                             7d8s21
      write(2)(iff0(i),i=1,mff0)                                        7d8s21
      return                                                            7d8s21
      end                                                               7d8s21
