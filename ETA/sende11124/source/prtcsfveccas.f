c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine prtcsfveccas(vec,ncsft,nrootx,nfcn,ibasis,ncsf,mdon,   12d28s19
     $     cprint,xnorm,nec,iunit,norb,irefo,iptrbit,mysym,mdoo,bc,ibc) 11d10s22
      implicit real*8 (a-h,o-z)
      character*80 line,line2                                           7d12s19
      dimension vec(ncsft,nrootx),ibasis(3,*),ncsf(*),irefo(*),         8d2s22
     $     iptrbit(2,mdoo+1,*)                                          8d2s22
      include "common.store"
      iacode=ichar('a')                                                 8d1s19
      do i=1,nrootx
       sum=0d0                                                          7d12s19
       xx=0d0                                                           7d12s19
       do j=1,ncsft
        sum=sum+vec(j,i)**2                                             7d12s19
        xx=max(xx,abs(vec(j,i)))                                        7d12s19
       end do
      end do
      ips=0                                                             7d12s19
      do if=1,nfcn                                                      7d12s19
       nclo=ibasis(1,if)                                                7d12s19
       nclop=nclo+1                                                     7d12s19
       iarg=nclop-mdon                                                  7d12s19
       nopen=nec-nclo*2                                                 7d12s19
       iic=iptrbit(1,nclop,mysym)+ibasis(2,if)-1                        8d2s22
       iio=iptrbit(2,nclop,mysym)+ibasis(3,if)-1                        8d2s22
       do i=1,ncsf(iarg)                                                7d12s19
        ip=i+ips                                                        7d12s19
        xx=abs(vec(ip,1))                                               7d12s19
        do j=2,nrootx                                                   7d12s19
         xx=max(xx,abs(vec(ip,j)))                                      9d6s19
        end do                                                          7d12s19
        if(xx.ge.cprint)then                                            9d6s19
         write(iunit,22)(vec(ip,j),j=1,nrootx)                          7d17s19
   22    format(57es22.14)                                              7d12s19
         do j=1,norb                                                    7d12s19
          line(j:j)='_'                                                 7d12s19
          if(btest(ibc(iic),j))line(j:j)='2'                            8d2s22
          if(btest(ibc(iio),j))line(j:j)='1'                            8d2s22
         end do                                                         7d12s19
   23    format(64a1)                                                   7d12s19
         is=1                                                           7d12s19
         ioff=0                                                         7d12s19
         do isb=1,8                                                     7d12s19
          do j=1,irefo(isb)                                             7d12s19
           jp=j+ioff                                                    7d12s19
           line2(is:is)=line(jp:jp)                                     7d12s19
           is=is+1                                                      7d12s19
          end do                                                        7d12s19
          if(irefo(isb).gt.0)itop=is                                    7d12s19
          line2(is:is)=' '                                              7d12s19
          is=is+1                                                       7d12s19
          ioff=ioff+irefo(isb)                                          7d12s19
         end do                                                         7d12s19
         ndig=1                                                         8d1s19
         itry=i                                                         8d1s19
         nmul=1                                                         8d1s19
  200    continue                                                       8d1s19
         itry=itry/26                                                   8d1s19
         if(itry.gt.0)then                                              8d1s19
          ndig=ndig+1                                                   8d1s19
          nmul=nmul*26                                                  8d1s19
          go to 200                                                     8d1s19
         end if                                                         8d1s19
         itry=i-1
         line(1:5)='     '                                              8d1s19
         do j=1,ndig                                                    6d3s22
          ij=itry/nmul
          nleft=itry-ij*nmul
          itry=nleft
          nmul=nmul/26
          line(j:j)=char(iacode+ij)
         end do
   21    format(i5)                                                     7d12s19
         write(iunit,23)(line2(j:j),j=1,itop),(' '),(line(j:j),j=1,5)   7d17s19
        end if                                                          7d12s19
       end do                                                           7d12s19
       ips=ips+ncsf(iarg)                                               7d12s19
      end do                                                            7d12s19
      return                                                            7d12s19
      end                                                               7d12s19
