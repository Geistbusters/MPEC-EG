c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dotdvecs(v1,idorb1,isorb1,ibasis1,iptr1,nfcn1,         10d2s19
     $     v2,idorb2,isorb2,ibasis2,iptr2,nfcn2,mdon,nec,ncsf,dot)      11d15s21
      implicit real*8 (a-h,o-z)                                         10d2s19
      logical ldebug                                                    5d7s21
      integer*1 idorb1(*),isorb1(*),idorb2(*),isorb2(*)                 10d2s19
      dimension v1(*),ibasis1(3,*),iptr1(4,*),v2(*),ibasis2(3,*),       10d2s19
     $     iptr2(4,*),ncsf(*)                                           11d15s21
      include "common.print"                                            5d7s21
      ldebug=iprtr(24).ne.0                                             5d7s21
      if(ldebug)write(6,*)('Hi, my name is dotdvecs')
      dot=0d0                                                           10d2s19
      ioff2=0                                                           10d2s19
      do if2=1,nfcn2
       nc2=ibasis2(1,if2)
       nc2p=nc2+1
       iarg2=nc2p-mdon                                                  10d2s19
       no2=nec-2*nc2                                                    10d2s19
       ic2=iptr2(2,nc2p)+nc2*(ibasis2(2,if2)-1)                         10d2s19
       io2=iptr2(4,nc2p)+no2*(ibasis2(3,if2)-1)                         10d2s19
       ioff1=0                                                          10d2s19
       if(ldebug)then
        write(6,*)('function2 '),if2,('c:o'),
     $       (idorb2(ic2+i),i=0,nc2-1),(':'),
     $       (isorb2(io2+i),i=0,no2-1)                                  5d7s21
       end if
       do if1=1,nfcn1                                                   10d2s19
        nc1=ibasis1(1,if1)                                              10d2s19
        nc1p=nc1+1                                                      10d2s19
        iarg1=nc1p-mdon                                                 10d2s19
        if(nc1.ne.nc2)go to 1                                           10d2s19
         ic1=iptr1(2,nc1p)+nc2*(ibasis1(2,if1)-1)                       10d2s19
         do i=0,nc2-1
          if(idorb1(ic1+i).ne.idorb2(ic2+i))go to 1
         end do
         io1=iptr1(4,nc1p)+no2*(ibasis1(3,if1)-1)                       10d2s19
         do i=0,no2-1                                                   10d2s19
          if(isorb1(io1+i).ne.isorb2(io2+i))go to 1                     10d2s19
         end do                                                         10d2s19
         if(ldebug)then                                                 5d7s21
          write(6,*)('matches function1 '),if1
          write(6,*)('v1 '),ioff1
          call prntm2(v1(ioff1+1),1,ncsf(iarg1),1)
          write(6,*)('v2 '),ioff2
          call prntm2(v2(ioff2+1),1,ncsf(iarg1),1)
         end if                                                         5d7s21
         do i=1,ncsf(iarg1)                                             10d2s19
          dot=dot+v1(ioff1+i)*v2(ioff2+i)                               10d2s19
         end do                                                         10d2s19
         if(ldebug)then                                                 5d7s21
          write(6,*)dot                                                 11d15s21
         end if                                                         5d7s21
    1   continue                                                        10d2s19
        ioff1=ioff1+ncsf(iarg1)                                         10d2s19
       end do                                                           10d2s19
       ioff2=ioff2+ncsf(iarg2)                                          10d2s19
      end do
      return
      end
