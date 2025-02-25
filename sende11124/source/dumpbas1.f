c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dumpbas1(nfcn,iptr,ibasis,idorb,isorb,mdoop,nec,ncsf,  4d28s21
     $     n1,nff0,iff0,mdon,nwavrec,nfcnx,nrestart)                    8d11s22
      implicit real*8 (a-h,o-z)                                         4d28s21
      integer*1 idorb(*),isorb(*)                                       4d28s21
      dimension ibasis(3,*),ncsf(*),n1(*),nff0(mdoop,3),iff0(*),        5d10s21
     $     iptr(4,*)                                                    5d10s21
      do i=1,mdoop                                                      4d28s21
       nff0(i,1)=0                                                      5d10s21
      end do                                                            4d28s21
      jff0=1                                                            4d28s21
      jvv=1                                                             4d28s21
      do i=1,nfcn                                                       4d28s21
       nclo=ibasis(1,i)                                                 4d28s21
       nclop=nclo+1                                                     4d28s21
       if(nff0(nclop,1).eq.0)then                                       5d10s21
        nff0(nclop,2)=jff0                                              5d10s21
        nff0(nclop,3)=jvv                                               5d10s21
       end if                                                           4d28s21
       nff0(nclop,1)=nff0(nclop,1)+1                                    5d10s21
       ic=iptr(2,nclop)+nclo*(ibasis(2,i)-1)                            4d28s21
       nopen=nec-2*nclo                                                 4d28s21
       io=iptr(4,nclop)+nopen*(ibasis(3,i)-1)                           4d28s21
       iff0(jff0)=0                                                     4d28s21
       do j=0,nclo-1                                                    4d28s21
        iff0(jff0)=ibset(iff0(jff0),idorb(ic+j))                        4d28s21
       end do                                                           4d28s21
       jff0=jff0+1                                                      4d28s21
       iff0(jff0)=0                                                     4d28s21
       do j=0,nopen-1                                                    4d28s21
        iff0(jff0)=ibset(iff0(jff0),isorb(io+j))                        4d28s21
       end do                                                           4d28s21
       jff0=jff0+1                                                      4d28s21
       iarg=nclop-mdon
       jvv=jvv+ncsf(iarg)                                               4d28s21
      end do                                                            4d28s21
      jff0=jff0-1                                                       4d30s21
      nwavrec=nwavrec+1                                                 5d4s21
      if(nrestart.lt.0)then                                             8d12s22
       read(2)necr,jff0r,mdonr,mdoopr,nfcnxr                            8d11s22
       ier=0                                                            8d11s22
       if(necr.ne.nec)ier=ier+1                                         8d11s22
       if(jff0r.ne.jff0)ier=ier+1                                       8d11s22
       if(mdonr.ne.mdon)ier=ier+1                                       8d11s22
       if(mdoopr.ne.mdoop)ier=ier+1                                     8d11s22
       if(nfcnxr.ne.nfcnx)ier=ier+1                                     8d11s22
       if(ier.ne.0)then                                                 8d11s22
        write(6,*)('restart data differs for record '),nwavrec          8d16s22
        write(6,*)('gotq  '),necr,jff0r,mdonr,mdoopr,nfcnxr              8d11s22
        write(6,*)('want '),nec,jff0,mdon,mdoop,nfcnx                   8d11s22
        stop 'restart'                                                  8d11s22
       end if                                                           8d11s22
      else                                                              8d11s22
       write(2)nec,jff0,mdon,mdoop,nfcnx                                 5d12s21
      end if                                                            8d11s22
      nwavrec=nwavrec+1                                                 5d4s21
      if(nrestart.lt.0)then                                             8d12s22
       read(2)idum                                                      8d11s22
      else                                                              8d11s22
       write(2)((nff0(j,i),i=1,3),j=1,mdoop)                             5d10s21
      end if                                                            8d11s22
      nwavrec=nwavrec+1                                                 5d4s21
      if(nrestart.lt.0)then                                             8d12s22
       read(2)idum                                                      8d11s22
      else                                                              8d11s22
       write(2)(iff0(i),i=1,jff0)                                        4d30s21
      end if                                                            8d11s22
      return                                                            4d28s21
      end                                                               4d28s21
