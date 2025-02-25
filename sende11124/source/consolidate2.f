c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine consolidate2(ibasis,nfcn,igoal,nhere)                  3d24s21
      dimension ibasis(3,*)                                             3d24s21
      igoalm=igoal-1                                                    3d24s21
      do i=1,nfcn
       if(ibasis(1,i).eq.igoalm)nhere=nhere+1                           3d24s21
      end do                                                            3d24s21
      return                                                            3d24s21
      end                                                               3d24s21
