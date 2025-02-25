c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine loaddouba(nff2,iff2,mff2,nsymb,mdon,mdoo,ncsf,ncsf2,   8d19s22
     $     inorm)                                                       8d19s22
      dimension nff2(mdoo+1,nsymb,2),iff2(*),ncsf(*),ncsf2(*)           7d21s21
      read(inorm)(iff2(i),i=1,mff2)                                     8d19s22
      read(inorm)(((nff2(i,j,k),i=1,mdoo+1),j=1,nsymb),k=1,2)           8d19s22
      read(inorm)(ncsf(i),i=1,mdoo+1-mdon)                              8d19s22
      read(inorm)(ncsf2(i),i=1,mdoo+1-mdon)                             8d19s22
      return                                                            7d8s21
      end                                                               7d8s21
