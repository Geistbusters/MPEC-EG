c mpec2.1 version theta. For copyright and Disclaimers, see start.f
      subroutine bodcfdb(ibdatf,nbdat,nfname,bc,ibc)                    2d3s25
      implicit real*8 (a-h,o-z)                                         2d3s25
      integer*8 iadd
      include "common.store"                                            2d3s25
      dimension ibdatf(*)                                               2d3s25
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      write(6,*)('Hi, my name is bodcfdb'),nbdat,ibcoff,nfname
      do i=1,nfname                                                     2d3s25
       if(mynowprog.ne.0)then                                           2d3s25
        iadd=ibcoff                                                     2d3s25
        ibcoff=iadd+nbdat                                               2d3s25
        ibdatf(i)=iadd                                                  2d3s25
       end if                                                           2d3s25
       write(6,*)('bcast ibdatf '),ibdatf(i)
       call dws_bcast(bc(ibdatf(i)),nbdat)                              2d3s25
       write(6,*)('back from bcast')
       write(6,*)('bdat ')
       call prntm2(bc(ibdatf(i)),9,28,9)
      end do                                                            2d3s25
      return                                                            2d3s25
      end                                                               2d3s25
