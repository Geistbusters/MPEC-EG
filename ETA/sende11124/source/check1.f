c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine check1(nff1r,nff1,iff1r,iff1,mff1,nsymb,mdon,mdoo,ncsf,8d18s22
     $     ncsfr)                                                       8d18s22
      implicit real*8 (a-h,o-z)                                         8d18s22
      dimension nff1r(mdoo+1,nsymb,*),nff1(mdoo+1,nsymb,*),             8d18s22
     $     iff1r(*),iff1(*),ncsf(*),ncsfr(*)                            8d18s22
      ier=0                                                             8d18s22
      do i=1,mff1                                                       8d18s22
       if(iff1r(i).ne.iff1(i))ier=ier+1                                 8d18s22
      end do                                                            8d18s22
      if(ier.ne.0)then                                                  8d18s22
       write(6,*)('iff1 changes: ')                                     8d18s22
       write(6,*)('got  :'),(iff1r(i),i=1,mff1)                         8d18s22
       write(6,*)('want :'),(iff1(i),i=1,mff1)                          8d18s22
       stop 'restart'                                                   8d18s22
      end if                                                            8d18s22
      do k=1,2                                                          8d18s22
       do j=1,nsymb                                                     8d18s22
        do i=mdon+1,mdoo+1                                              9d29s22
         if(nff1r(i,j,k).ne.nff1(i,j,k))ier=ier+1                       8d18s22
        end do                                                          8d18s22
       end do                                                           8d18s22
      end do                                                            8d18s22
      if(ier.ne.0)then                                                  8d18s22
       write(6,*)('nff1 changes: ')                                     8d18s22
       write(6,*)('got  :'),
     $      (((nff1r(i,j,k),i=mdon+1,mdoo+1),j=1,nsymb),k=1,2)          9d29s22
       write(6,*)('want :'),(((nff1(i,j,k),i=mdon+1,mdoo+1),j=1,nsymb), 9d29s22
     $      k=1,2)                                                      9d29s22
       stop 'restart'                                                   8d18s22
      end if                                                            8d18s22
      do i=1,mdoo+1-mdon                                                8d18s22
       if(ncsf(i).ne.ncsfr(i))ier=ier+1                                 8d18s22
      end do                                                            8d18s22
      if(ier.ne.0)then                                                  8d18s22
       write(6,*)('ncsf changes: ')                                     8d18s22
       write(6,*)('got  :'),(ncsfr(i),i=1,mdoo+1-mdon)                  8d18s22
       write(6,*)('want :'),(ncsf(i),i=1,mdoo+1-mdon)                   8d18s22
       stop 'restart'                                                   8d18s22
      end if                                                            8d18s22
      return                                                            8d18s22
      end                                                               8d18s22
