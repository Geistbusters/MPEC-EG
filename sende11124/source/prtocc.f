c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine prtocc(coef,isb,jsb,ia,ib,nadet,nbdet,iaorb,iborb,     4d20s18
     $     nalpha,nbeta,norb,iacto,nsymb,bc,ibc)                        11d9s22
      implicit real*8 (a-h,o-z)                                         4d20s18
      character*80 line
      integer*1 iaorb(nalpha,*),iborb(nbeta,*)                          4d20s18
      dimension nadet(*),nbdet(*),iacto(*)                              4d20s18
      character*1 code(4)                                               4d20s18
      data code/'_','a','b','2'/
      include "common.store"
      iau=ia                                                            4d20s18
      do i=1,isb-1                                                      4d20s18
       iau=iau+nadet(i)                                                 4d20s18
      end do                                                            4d20s18
      ibu=ib                                                            5d7s18
      do i=1,jsb-1
       ibu=ibu+nbdet(i)
      end do                                                            4d20s18
      iorb=ibcoff                                                       4d20s18
      ibcoff=iorb+norb                                                  4d20s18
      jorb=iorb-1                                                       4d20s18
      call enough('prtocc.  1',bc,ibc)
      do i=0,norb-1                                                      4d20s18
       ibc(iorb+i)=1                                                    4d20s18
      end do                                                            4d20s18
      do i=1,nalpha
       ibc(jorb+iaorb(i,iau))=2                                         4d20s18
      end do                                                            4d20s18
      do i=1,nbeta                                                      4d20s18
       ibc(jorb+iborb(i,ibu))=ibc(jorb+iborb(i,ibu))+2                  4d20s18
      end do                                                            4d20s18
      is=1                                                              4d20s18
      ioff=0                                                            4d20s18
      do isym=1,nsymb                                                   4d20s18
       do i=1,iacto(isym)                                               4d20s18
        if(is.gt.80)stop 'prtocc'
        line(is:is)=code(ibc(jorb+ioff+i))                              4d20s18
        is=is+1                                                         4d20s18
       end do                                                           4d20s18
       if(is.gt.80)stop 'prtocc'
       line(is:is)=' '                                                  4d20s18
       is=is+1                                                          4d20s18
       ioff=ioff+iacto(isym)                                            4d20s18
      end do                                                            4d20s18
      ie=is-1                                                           4d20s18
      write(6,1)coef,(line(i:i),i=1,ie)
    1 format(f8.4,1x,50a1)                                              4d20s18
      ibcoff=iorb                                                       4d20s18
      return                                                            4d20s18
      end                                                               4d20s18
