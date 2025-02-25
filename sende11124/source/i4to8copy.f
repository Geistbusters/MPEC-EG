c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine i4to8copy(in,iou,nb4,n)                                2d1s21
      dimension in(*),iou(*)                                            2d1s21
      do i=1,n                                                          2d1s21
       iou(nb4+i)=in(i)                                                 2d1s21
      end do                                                            2d1s21
      return                                                            2d1s21
      end                                                               2d1s21
