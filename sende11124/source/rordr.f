c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine rordr(i4od,noc,nsblkder,isblkder,i4oduu,bc,ibc)        11d10s22
      implicit real*8 (a-h,o-z)
      dimension i4od(*),noc(*),isblkder(4,*),i4oduu(*)                   5d16s22
      include "common.hf"
      include "common.store"
      do is=1,nsdlk
       if(isblk(1,is).eq.isblk(2,is))then
        nrow=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2
       else
        nrow=noc(isblk(1,is))*noc(isblk(2,is))
       end if
       ncol=noc(isblk(3,is))*noc(isblk(4,is))
       i4oduu(is)=ibcoff                                                5d16s22
       if(min(nrow,ncol).gt.0)then                                      5d16s22
        ibcoff=i4oduu(is)+nrow*ncol                                     5d17s22
        call enough('rordr.  1',bc,ibc)
        do isd=1,nsblkder                                               5d17s22
         if(isblk(1,is).eq.isblkder(1,isd).and.                         5d17s22
     $      isblk(2,is).eq.isblkder(2,isd).and.                         5d17s22
     $      isblk(3,is).eq.isblkder(3,isd))then                         5d17s22
          do i=0,nrow*ncol-1                                            5d17s22
           bc(i4oduu(is)+i)=bc(i4od(isd)+i)                             5d17s22
          end do                                                        5d17s22
          go to 1                                                       5d17s22
         end if                                                         5d17s22
        end do                                                          5d17s22
        write(6,*)('in rordr, looking for '),(isblk(j,is),j=1,4)
        write(6,*)('but could not find it in ')
        do isd=1,nsblkder
         write(6,*)(isblkder(j,isd),j=1,4)
        end do
        call dws_synca
        call dws_finalize
        stop
    1   continue                                                        5d17s22
       end if                                                           5d16s22
      end do                                                            5d16s22
      return                                                            5d17s22
      end                                                               5d17s22
