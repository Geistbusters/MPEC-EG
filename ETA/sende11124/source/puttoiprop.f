c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine puttoiprop(ismult,nec,mdon,mdoop)                      8d16s22
      implicit real*8 (a-h,o-z)                                         8d16s22
      include "common.basis"                                            8d16s22
      include "common.input"                                            8d16s22
      iprop(1)=ismult                                                   8d16s22
      iprop(2)=ismult                                                   8d16s22
      iprop(3)=nec                                                      8d16s22
      iprop(4)=mdon                                                     8d16s22
      iprop(5)=mdoop                                                    8d16s22
      return                                                            8d16s22
      end                                                               8d16s22
