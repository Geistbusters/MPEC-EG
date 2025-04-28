c     Copyright © 2023 United States Government as represented by the Administrator of the
c     National Aeronautics and Space Administration.  All Rights Reserved.
c
c Disclaimers
c
c     No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND,
c     EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY
c     THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF
c     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY
c     WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION,
c     IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER,
c     CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS,
c     RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM
c     USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND
c     LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND
c     DISTRIBUTES IT "AS IS."
c
c     Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED
c     STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.
c     IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES,
c     EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON,
c     OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY
c     AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS
c     WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY
c     FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT. 
c 
      implicit real*8 (a-h,o-z)
c
c     mpec2.1 version eta
c
      common/ddidwscm/dws_me,dws_np                                     10d23s14
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      call dws_preinit                                                  1d31s21
      call wrapit(bc,ibc)
      stop                                                              10d31s22
      end
      subroutine wrapit(bc,ibc)                                         11d9s22
      include "common.store"                                            11d9s22
      maxbc4=maxbc                                                      10d31s22
      maxbct=maxbc                                                      11d15s22
      call mkmem(maxbc4)                                                10d31s22
      return                                                            11d9s22
      end                                                               11d9s22
      subroutine dgemm(tl,tr,n1,n2,n3,fab,a,ia,b,ib,fc,c,ic,lab)
      implicit real*8 (a-h,o-z)
      character*(*) lab
      character*1 tl,tr
      logical logl,logr
      dimension a(ia,1),b(ib,1),c(ic,1)
      logl=tl.eq.'n'.or.tl.eq.'N'
      logr=tr.eq.'n'.or.tr.eq.'N'
      if(min(n1,n2,n3,ia,ib,ic).lt.1.or.n1.gt.ic                        11d2s22
     $     .or.(logl.and.n1.gt.ia)                                      11d2s22
     $     .or.(logr.and.n3.gt.ib).or.(.not.logr.and.n2.gt.ib))then     10d7s22
       write(6,*)('bad dgemm args: '),n1,n2,n3,ia,ib,ic
       write(6,*)lab(1:min(20,len(lab)))
       write(6,*)('tl,tr: '),tl,logl,tr,logr
       stop
      end if
      if(fc.eq.0d0)then
       do i=1,n2
        do j=1,n1
         c(j,i)=0d0
        end do
       end do
      else
      do 1 i=1,n2
       do 2 j=1,n1
        c(j,i)=c(j,i)*fc
    2  continue
    1 continue
      end if
      if(logl.and.logr)then
       if(mod(n3,2).eq.0)then                                           5d24s19
        n3top=n3                                                        5d24s19
       else                                                             5d24s19
        n3top=n3-1                                                      5d24s19
       end if                                                           5d24s19
       if(mod(n2,2).eq.0)then                                           5d24s19
        n2top=n2                                                        5d24s19
       else                                                             5d24s19
        n2top=n2-1                                                      5d24s19
       end if                                                           5d24s19
       do 10 i=1,n2top,2                                                5d24s19
        ip=i+1
        do 11 j=1,n3top,2
         jp=j+1
         do 12 k=1,n1
          c(k,i)=c(k,i)+fab*(a(k,j)*b(j,i)+a(k,jp)*b(jp,i))
          c(k,ip)=c(k,ip)+fab*(a(k,j)*b(j,ip)+a(k,jp)*b(jp,ip))
   12    continue
   11   continue
   10  continue
       if(mod(n3,2).ne.0)then
        do i=1,n2top,2
         ip=i+1
         x=b(n3,i)*fab
         y=b(n3,ip)*fab
         do k=1,n1
          c(k,i)=c(k,i)+x*a(k,n3)
          c(k,ip)=c(k,ip)+y*a(k,n3)
         end do
        end do
       end if
       if(mod(n2,2).ne.0)then                                           5d24s19
        do j=1,n3
         x=fab*b(j,n2)
         do k=1,n1                                                      5d24s19
          c(k,n2)=c(k,n2)+a(k,j)*x                                      5d24s19
         end do                                                         5d24s19
        end do                                                          5d24s19
       end if                                                           5d24s19
      else if(logl.and..not.logr)then
       do 20 i=1,n2
        do 21 j=1,n3
         x=b(i,j)*fab
         do 22 k=1,n1
          c(k,i)=c(k,i)+a(k,j)*x
   22    continue
   21   continue
   20  continue
c$$$      else if(tl.eq.'t'.and.tr.eq.'n')then
      else if(.not.logl.and.logr)then
       do 30 i=1,n2
        do 31 j=1,n3
         x=b(j,i)*fab
         do 32 k=1,n1
          c(k,i)=c(k,i)+a(j,k)*x
   32    continue
   31   continue
   30  continue
c$$$      else if(tl.eq.'t'.and.tr.eq.'t')then
      else if(.not.logl.and..not.logr)then
       do 40 i=1,n2
        do 41 j=1,n3
         x=b(i,j)*fab
         do 42 k=1,n1
          c(k,i)=c(k,i)+a(j,k)*x
   42    continue
   41   continue
   40  continue
      else
       write(6,*)('in dgemm, unknown transpose options'),tl,tr
       stop
      end if
      return
      end
