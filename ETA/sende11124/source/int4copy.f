c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine int4copy(i4od,i4odx,n4o,noc,nsblkder,isblkder,bc,ibc)  11d10s22
      implicit real*8 (a-h,o-z)                                         5d16s22
      dimension i4od(*),i4odx(*),noc(*),isblkder(4,*)                   5d16s22
      include "common.store"                                            5d16s22
      n4o=0                                                             5d16s22
      do is=1,nsblkder                                                  5d16s22
       i4odx(is)=ibcoff                                                 5d16s22
       nn=noc(isblkder(1,is))*noc(isblkder(2,is))*noc(isblkder(3,is))   5d16s22
     $      *noc(isblkder(4,is))                                        5d16s22
       ibcoff=i4odx(is)+nn                                              5d16s22
       call enough('int4copy.  1',bc,ibc)
       do i=0,nn-1                                                      5d16s22
        bc(i4odx(is)+i)=bc(i4od(is)+i)                                  5d20s22
       end do                                                           5d16s22
       n4o=n4o+nn                                                       5d16s22
      end do                                                            5d16s22
      return                                                            5d16s22
      end                                                               5d16s22
