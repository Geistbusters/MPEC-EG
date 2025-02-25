c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine m12gen(ia,id2,id1,no,nconf,n1p,n1m,
     $                  idat1,m1,itmp)
      implicit real*8 (a-h,o-z)
c
c     given det string, determine single and double matrix
c     element differences
c
      integer*1 ia(id1,id2),itmp(no)
      integer*8 idat1
      dimension idat1(4,m1),itmp2(2,2)
      write(6,*)('m1,m2 '),m1,m2
      write(6,*)('in modified m12gen')
      n1p=0
      n1m=m1+1
      n2p=0
      n2m=m2+1
      do iket=1,nconf
       do ibra=1,iket-1
    1   format('string ',10i1)
        isum=0
        do i=1,no
         itmp(i)=ia(i,ibra)-ia(i,iket)
         if(itmp(i).gt.0)then                                           8d8s14
          isum=isum+itmp(i)                                             8d8s14
         else                                                           8d8s14
          isum=isum-itmp(i)                                             8d8s14
         end if                                                         8d8s14
        end do
        if(isum.eq.2)then
         do i=1,no
          if(itmp(i).eq.1)nfrom=i
          if(itmp(i).eq.-1)nto=i
         end do
         nb=0
         do i=min(nfrom,nto)+1,max(nfrom,nto)-1
          if(ia(i,ibra).eq.1)nb=nb+1
         end do
         if(mod(nb,2).eq.0)then
          n1p=n1p+1
          if(n1p.ge.n1m)then
           write(6,*)('singles stack collides!+'),n1p,n1m
           stop
          end if
          idat1(1,n1p)=ibra
          idat1(2,n1p)=iket
          idat1(3,n1p)=nfrom
          idat1(4,n1p)=nto
         else
          n1m=n1m-1
          if(n1m.le.n1p)then
           write(6,*)('singles stack collides!-'),n1p,n1m
           stop
          end if
          idat1(1,n1m)=ibra
          idat1(2,n1m)=iket
          idat1(3,n1m)=nfrom
          idat1(4,n1m)=nto
         end if
        end if
       end do
      end do
      return
      end
