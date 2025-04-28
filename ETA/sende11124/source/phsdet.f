c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine phsdet(ia,itmp,nfrom,iphs,no)
      implicit real*8 (a-h,o-z)
      integer*1 ia(1),itmp(1)
      integer nfrom(2)
      io=1
      do i=1,no
       if(ia(i).eq.1)then
        itmp(io)=i
        io=io+1
       end if
      end do
      ne=io-1
      iphs=0
      do ipass=1,2
       if(itmp(ipass).ne.nfrom(ipass))then
        do i=ipass+1,ne
         if(itmp(i).eq.nfrom(ipass))then
          iphs=iphs+1
          itmp(i)=itmp(ipass)
          itmp(ipass)=nfrom(ipass)
          go to 4
         end if
        end do
    4   continue
       end if
      end do
    5 continue
      do i=3,ne-1
       do j=i+1,ne
        if(itmp(i).gt.itmp(j))then
         icpy=itmp(i)
         itmp(i)=itmp(j)
         itmp(j)=icpy
         iphs=iphs+1
         go to 5
        end if
       end do
      end do
      return
      end
