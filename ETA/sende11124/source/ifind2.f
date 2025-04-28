c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      integer function ifind2(n1,n2,n3,n4,icase)
      include "common.hf"
      dimension inv(2,8,8,8)                                            4d9s18
      common/fnd2cm/inv                                                 4d9s18
      save                                                              4d9s18
      data icall/0/                                                     4d9s18
      if(nsdlk.eq.1)then                                                4d9s18
       ifind2=1                                                         4d9s18
       icase=1                                                          4d9s18
       return                                                           4d9s18
      end if                                                            4d9s18
      icall=icall+1                                                     4d9s18
      if(icall.eq.1)then                                                4d9s18
       do i1=1,nsymb                                                    4d9s18
        do i2=1,nsymb                                                   4d9s18
         do i3=1,nsymb                                                  4d9s18
          inv(1,i1,i2,i3)=1                                             1d19s23
          do i=1,nsdlk                                                  4d9s18
           if(isblk(1,i).eq.i1.and.                                     4d9s18
     $          isblk(2,i).eq.i2.and.                                   4d9s18
     $          isblk(3,i).eq.i3)then                                   4d9s18
            inv(1,i1,i2,i3)=i                                           4d9s18
            inv(2,i1,i2,i3)=1                                           4d9s18
            go to 1                                                     4d9s18
           else if(isblk(2,i).eq.i1.and.                                4d9s18
     $           isblk(1,i).eq.i2.and.                                  4d9s18
     $           isblk(4,i).eq.i3)then                                  4d9s18
            inv(1,i1,i2,i3)=i                                           4d9s18
            inv(2,i1,i2,i3)=2                                           4d9s18
            go to 1                                                     4d9s18
           else if(isblk(2,i).eq.i2.and.                                4d9s18
     $           isblk(1,i).eq.i1.and.                                  4d9s18
     $           isblk(4,i).eq.i3)then                                  4d9s18
            inv(1,i1,i2,i3)=i                                           4d9s18
            inv(2,i1,i2,i3)=3                                           4d9s18
            go to 1                                                     4d9s18
           else if(isblk(2,i).eq.i1.and.                                4d9s18
     $           isblk(1,i).eq.i2.and.                                  4d9s18
     $           isblk(3,i).eq.i3)then                                  4d9s18
            inv(1,i1,i2,i3)=i                                           4d9s18
            inv(2,i1,i2,i3)=4                                           4d9s18
            go to 1                                                     4d9s18
           end if                                                       4d9s18
          end do                                                        4d9s18
    1     continue                                                      4d9s18
         end do                                                         4d9s18
        end do                                                          4d9s18
       end do                                                           4d9s18
      end if                                                            4d9s18
      ifind2=inv(1,n1,n2,n3)                                            4d9s18
      icase=inv(2,n1,n2,n3)                                             4d9s18
      return
      end
