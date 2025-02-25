c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine getbas0(nff0,mdoop,iff0,jff0)                          4d30s21
      implicit real*8 (a-h,o-z)                                         4d30s21
      dimension  nff0(mdoop,3),iff0(*)                                  5d10s21
c
c     it should be noted that iff0 is 32 bit, but if it is copied       8d11s22
c     to a 32 bit integer, the sign bit might be changed, so to avoid   8d11s22
c     this, copy it to a 64 bit integer, then the sign bit is out of    8d11s22
c     harms way.                                                        8d11s22
c
      read(1)((nff0(j,i),i=1,3),j=1,mdoop)                              5d10s21
      read(1)(iff0(i),i=1,jff0)                                         5d3s21
      return                                                            4d30s21
      end                                                               4d30s21
