c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine cnt12(istng,no,nconf,m1,m2,ndet,nsymb,itab,m1c)        8d23s06
      implicit real*8 (a-h,o-z)
      integer*1 istng(no,nconf)
      dimension m1(8),m2(8),ndet(8),ioff(8),itab(8,8),m1c(36)           8d23s06
c
c     count no. of singles and double differences
c
      nn=(nsymb*(nsymb+1))/2                                            8d24s06
      do i=1,nn                                                         8d23s06
       m1c(i)=0                                                         8d23s06
      end do                                                            8d23s06
      do i=1,nsymb                                                      8d8s06
       m1(i)=0                                                          8d8s06
       m2(i)=0                                                          8d8s06
       if(i.eq.1)then
        ioff(1)=ndet(i)                                                 8d8s06
       else                                                             8d8s06
        ioff(i)=ioff(i-1)+ndet(i)                                       8d8s06
       end if                                                           8d8s06
      end do                                                            8d8s06
      isymb=1                                                           8d8s06
      do ibra=1,nconf
       ist=isymb                                                        5d7s18
       do ipass=ist,nsymb                                               5d7s18
        if(ibra.gt.ioff(isymb))isymb=isymb+1                             8d8s06
       end do                                                           5d7s18
 1010  format(a3,1x,2i8,5x,10i1)
       isymk=1                                                          8d8s06
       do iket=1,ibra-1
        ist=isymk                                                       5d7s18
        do ipass=ist,nsymb                                              5d7s18
         if(iket.gt.ioff(isymk))isymk=isymk+1                            8d8s06
        end do                                                          5d7s18
        i12=itab(isymb,isymk)                                           8d8s06
        isum=0
        do i=1,no
         if(istng(i,ibra).ne.istng(i,iket))then
          isum=isum+1
         end if
        end do
        if(i12.eq.1)then                                                8d23s06
         if(isum.eq.2)then
          m1(isymb)=m1(isymb)+1                                         8d23s06
         else if(isum.eq.4)then
          m2(isymb)=m2(isymb)+1                                         8d23s06
         end if
        else if(isum.eq.2)then                                          8d23s06
         ii=((isymb*(isymb-1))/2)+isymk                                 8d24s06
         m1c(ii)=m1c(ii)+1                                              8d24s06
        end if                                                          8d23s
       end do
      end do
      return
      end
