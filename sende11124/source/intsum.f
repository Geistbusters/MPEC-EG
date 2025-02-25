c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine intsum(vec,ncsft,nroot,idorb,isorb,ibasis,iptr,        1d15s22
     $     nfcn,mdon,nec,ncsf,norb,irefo,nsymb)                         2d10s22
      implicit real*8 (a-h,o-z)                                         1d15s22
      character*50 oline                                                2d10s22
      character*500 line                                                2d10s22
      integer*1 idorb(*),isorb(*)                                       1d15s22
      dimension vec(ncsft,nroot),ibasis(3,*),iptr(4,*),ncsf(*),irefo(*) 2d11s22
      ioff=0                                                            1d15s22
      do if=1,nfcn                                                      1d15s22
       nc=ibasis(1,if)                                                  1d15s22
       ncp=nc+1                                                         1d15s22
       iarg=ncp-mdon                                                    1d15s22
       no=nec-2*nc                                                      1d15s22
       ic=iptr(2,ncp)+nc*(ibasis(2,if)-1)                               1d15s22
       io=iptr(4,ncp)+no*(ibasis(3,if)-1)                               1d15s22
       do j=1,ncsf(iarg)                                                1d15s22
        xm=0d0                                                          1d15s22
        do ir=1,nroot                                                   1d15s22
         if(abs(vec(ioff+j,ir)).gt.xm)then                              1d15s22
          xm=abs(vec(ioff+j,ir))                                        1d15s22
         end if                                                         1d15s22
        end do                                                          1d15s22
        if(xm.gt.0.1d0)then                                             1d15s22
         do i=1,norb                                                    2d10s22
          oline(i:i)='_'                                                 2d10s22
         end do                                                         2d10s22
         do i=0,nc-1                                                    2d10s22
          oline(idorb(ic+i):idorb(ic+i))='2'                             2d10s22
         end do                                                         2d10s22
         do i=0,no-1                                                    2d10s22
          oline(isorb(io+i):isorb(io+i))='1'                             2d10s22
         end do                                                         2d10s22
         iline=1
         jline=1                                                        2d10s22
         do isb=1,nsymb                                                 2d10s22
          do i=1,irefo(isb)                                             2d10s22
           line(iline:iline)=oline(jline:jline)
           iline=iline+1                                                2d10s22
           jline=jline+1                                                2d11s22
          end do                                                        2d10s22
          line(iline:iline)=' '                                         2d10s22
          iline=iline+1                                                 2d10s22
         end do                                                         2d10s22
         line(iline:iline)=char(ichar('a')+j-1)                         2d11s22
         iline=iline+1                                                  2d11s22
         line(iline:iline)=' '                                          2d11s22
         iline=iline+1                                                  2d11s22
         iend=iline+nroot*12-1                                          2d10s22
         write(line(iline:iend),1)(vec(ioff+j,ir),ir=1,nroot)           2d10s22
    1    format(20f12.6)                                                1d15s22
         write(6,2)(line(i:i),i=1,iend)                                 2d10s22
    2    format(500a1)                                                  2d10s22
        end if                                                          1d15s22
       end do                                                           1d15s22
       ioff=ioff+ncsf(iarg)                                             1d17s22
      end do                                                            1d15s22
      return                                                            1d15s22
      end                                                               1d15s22
