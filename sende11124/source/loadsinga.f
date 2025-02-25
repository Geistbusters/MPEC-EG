c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine loadsinga(nff1,iff1,mff1,nsymb,mdon,mdoo,ncsf,iunit)   8d18s22
      dimension nff1(mdoo+1,nsymb,2),iff1(*),ncsf(*)                    7d8s21
      read(iunit)(iff1(i),i=1,mff1)                                     8d18s22
      read(iunit)(((nff1(i,j,k),i=mdon+1,mdoo+1),j=1,nsymb),k=1,2)      9d12s22
      read(iunit)(ncsf(i),i=1,mdoo+1-mdon)                              8d18s22
      return                                                            7d8s21
      end                                                               7d8s21
