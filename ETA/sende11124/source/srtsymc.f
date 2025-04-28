c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine srtsymc(idata,m1,ism,nsymb,idata8,jdata8,mnz,nz)       4d30s18
      integer*2 idata(4,m1),itmp2                                       4d26s18
      integer*1 itmp1(2)                                                4d26s18
      integer*8 ism(*),idata8(m1),jdata8(m1)                            4d26s18
      equivalence (itmp2,itmp1)                                         4d26s18
      integer m(8,8),ioff(8,8),mnz(3,64)                                4d30s18
      do i=1,nsymb                                                      4d26s18
       do j=1,nsymb                                                     4d30s18
        m(j,i)=0                                                        4d30s18
       end do                                                            4d26s18
      end do                                                            4d30s18
      do i=1,m1                                                         4d26s18
       itmp2=idata(3,i)                                                 4d26s18
       if(ism(itmp1(1)).lt.0.or.ism(itmp1(1)).gt.nsymb)then
        write(6,*)('sym1 error: '),ism(itmp1(1))
        call dws_sync
        call dws_finalize
        stop
       end if
       if(ism(itmp1(2)).lt.0.or.ism(itmp1(2)).gt.nsymb)then
        write(6,*)('sym2 error: '),ism(itmp1(2))
        call dws_sync
        call dws_finalize
        stop
       end if
       m(ism(itmp1(1)),ism(itmp1(2)))=m(ism(itmp1(1)),ism(itmp1(2)))+1  4d30s18
      end do                                                            4d26s18
      ii=1                                                              4d26s18
      nz=0                                                              4d30s18
      do i=1,nsymb                                                      4d26s18
       do j=1,nsymb                                                     4d30s18
        ioff(j,i)=ii                                                    4d30s18
        ii=ii+m(j,i)                                                    4d30s18
        if(m(j,i).ne.0)then                                             4d30s18
         nz=nz+1                                                        4d30s18
         mnz(1,nz)=j
         mnz(2,nz)=i
         mnz(3,nz)=m(j,i)
        end if
       end do                                                           4d30s18
      end do                                                            4d26s18
      do i=1,m1                                                         4d26s18
       itmp2=idata(3,i)                                                 4d26s18
       isto=ioff(ism(itmp1(1)),ism(itmp1(2)))                           4d30s18
       jdata8(isto)=idata8(i)                                           4d26s18
       ioff(ism(itmp1(1)),ism(itmp1(2)))=                               4d30s18
     $      ioff(ism(itmp1(1)),ism(itmp1(2)))+1                         4d30s18
      end do                                                            4d26s18
      do i=1,m1
       idata8(i)=jdata8(i)
      end do
      return
      end
