c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dumpbas0(nff0,mdoop,iff0,jff0)                         2d14s22
      implicit real*8 (a-h,o-z)                                         4d30s21
      dimension  nff0(mdoop,3),iff0(*)                                  5d10s21
      write(2)((nff0(j,i),i=1,3),j=1,mdoop)                             2d14s22
      write(2)(iff0(i),i=1,jff0)                                        2d14s22
      return                                                            4d30s21
      end                                                               4d30s21
