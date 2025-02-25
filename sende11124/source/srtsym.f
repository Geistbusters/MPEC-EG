c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine srtsym(idata,m1,ism,nsymb,idata8,jdata8,m)             4d26s18
      integer*2 idata(4,m1),itmp2                                       4d26s18
      integer*1 itmp1(2)                                                4d26s18
      integer*8 ism(*),idata8(m1),jdata8(m1)                            4d26s18
      equivalence (itmp2,itmp1)                                         4d26s18
      integer m(8),ioff(8)                                              4d26s18
      do i=1,nsymb                                                      4d26s18
       m(i)=0                                                           4d26s18
      end do                                                            4d26s18
      do i=1,m1                                                         4d26s18
       itmp2=idata(3,i)                                                 4d26s18
       m(ism(itmp1(1)))=m(ism(itmp1(1)))+1                              4d26s18
      end do                                                            4d26s18
      ii=1                                                              4d26s18
      do i=1,nsymb                                                      4d26s18
       ioff(i)=ii                                                       4d26s18
       ii=ii+m(i)                                                       4d26s18
      end do                                                            4d26s18
      do i=1,m1                                                         4d26s18
       itmp2=idata(3,i)                                                 4d26s18
       isto=ioff(ism(itmp1(1)))                                         4d26s18
       jdata8(isto)=idata8(i)                                           4d26s18
       ioff(ism(itmp1(1)))=ioff(ism(itmp1(1)))+1                        4d26s18
      end do                                                            4d26s18
      do i=1,m1
       idata8(i)=jdata8(i)
      end do
      return
      end
