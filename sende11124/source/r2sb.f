c mpec2.1 version eta. For copyright and Disclaimers, see start.f
       subroutine r2sb(ibasis,isbasis,nfcn,myfcn)                       1d25s21
c
c     restrict to symmetry basis
c
      dimension ibasis(3,*),isbasis(3,*)                                7d5s19
      myfcn=0                                                           7d5s19
      do if=1,nfcn                                                      7d5s19
        myfcn=myfcn+1                                                   7d5s19
        isbasis(1,myfcn)=ibasis(1,if)                                   7d5s19
        isbasis(2,myfcn)=ibasis(2,if)                                   7d5s19
        isbasis(3,myfcn)=ibasis(3,if)                                   7d5s19
      end do                                                            7d5s19
      return                                                            7d5s19
      end                                                               7d5s19
