c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine diaglzz1(tmpa,numa,nalpha,iaorb,ixlzze,debug,nsymb,    12d22s19
     $     ndet,ism,irelo,islz,multh,npass,bc,ibc)                      11d14s22
      implicit real*8 (a-h,o-z)
      integer*1 iaorb(nalpha,numa)
      dimension tmpa(numa)
      include "common.store"
      include "common.cas"
      integer*8 ism(1),irelo(1)                                         5d7s18
      dimension ixlzze(8,*),ndet(1),multh(8,8),islz(*)                  12d31s19
      ione=npass+1                                                      12d31s19
      do i=1,numa
       sum=0d0
       do j=1,nalpha
        jo=irelo(iaorb(j,i))                                            8d8s06
        nj=ism(iaorb(j,i))                                              8d8s06
        iadd=ixlzze(nj,ione)+(jo-1)*(iacto(nj)+1)                       12d31s19
        sum=sum+bc(iadd)                                                7d12s06
        do ipass=1,npass                                                12d31s19
         do jp=1,nalpha
          jpo=irelo(iaorb(jp,i))                                         8d8s06
          njp=ism(iaorb(jp,i))                                           8d8s06
          if(multh(nj,njp).eq.islz(ipass))then                          12d31s19
           iad12=ixlzze(nj,ipass)+(jo-1)+iacto(nj)*(jpo-1)              12d31s19
           term=+2d0*bc(iad12)*bc(iad12)                                 12d22s19
           sum=sum-0.5d0*term                                            12d22s19
          end if                                                         12d22s19
         end do
        end do                                                          12d31s19
       end do
       tmpa(i)=sum
      end do
      return
      end
