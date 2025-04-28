c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine check3(nff2r,nff2,iff2r,iff2,mff2,nsymb,mdon,mdoo,ncsf,8d18s22
     $     ncsfr,ncsf2,ncsf2r)                                          8d19s22
      implicit real*8 (a-h,o-z)                                         8d18s22
      dimension nff2r(mdoo+1,nsymb,*),nff2(mdoo+1,nsymb,*),             8d19s22
     $     iff2r(*),iff2(*),ncsf(*),ncsfr(*),ncsf2(4,*),ncsf2r(*)       8d19s22
      ier=0                                                             8d18s22
      do i=1,mff2                                                       8d19s22
       if(iff2r(i).ne.iff2(i))ier=ier+1                                 8d19s22
      end do                                                            8d18s22
      if(ier.ne.0)then                                                  8d18s22
       write(6,*)('iff2 changes: ')                                     8d19s22
       write(6,*)('got  :'),(iff2r(i),i=1,mff2)                         8d19s22
       write(6,*)('want :'),(iff2(i),i=1,mff2)                          8d19s22
       stop 'restart'                                                   8d18s22
      end if                                                            8d18s22
      do k=1,2                                                          8d18s22
       do j=1,nsymb                                                     8d18s22
        do i=1,mdoo+1                                                   8d18s22
         if(nff2r(i,j,k).ne.nff2(i,j,k))ier=ier+1                       8d19s22
        end do                                                          8d18s22
       end do                                                           8d18s22
      end do                                                            8d18s22
      if(ier.ne.0)then                                                  8d18s22
       write(6,*)('nff2 changes: ')                                     8d19s22
       write(6,*)('got  :'),
     $      (((nff2r(i,j,k),i=1,mdoo+1),j=1,nsymb),k=1,2)               8d19s22
       write(6,*)('want :'),(((nff2(i,j,k),i=1,mdoo+1),j=1,nsymb),k=1,2)8d19s22
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
      do i=1,mdoo+1-mdon                                                8d18s22
       if(ncsf2(1,i).ne.ncsf2r(i))ier=ier+1                             8d19s22
      end do                                                            8d18s22
      if(ier.ne.0)then                                                  8d18s22
       write(6,*)('ncsf2 changes: ')                                    8d19s22
       write(6,*)('got  :'),(ncsf2r(i),i=1,mdoo+1-mdon)                 8d19s22
       write(6,*)('want :'),(ncsf2(1,i),i=1,mdoo+1-mdon)                8d19s22
       stop 'restart'                                                   8d18s22
      end if                                                            8d18s22
      return                                                            8d18s22
      end                                                               8d18s22
