c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine mov2ps(vecx,vecps,ncsft,ipoint,nps,nroot)              3d25s21
      implicit real*8 (a-h,o-z)                                         3d25s21
      dimension vecx(ncsft,*),vecps(nps,*),ipoint(*)                    3d25s21
      do ir=1,nroot                                                     3d25s21
       do i=1,nps                                                       3d25s21
        vecps(i,ir)=vecx(ipoint(i),ir)                                  3d25s21
       end do                                                           3d25s21
      end do                                                            3d25s21
      return                                                            3d25s21
      end                                                               3d25s21
