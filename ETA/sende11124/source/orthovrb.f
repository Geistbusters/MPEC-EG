c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine orthovrb(ovr,nv,nok,thrs,ivec,bc,ibc,isb,l,nder)       8d31s23
      implicit real*8 (a-h,o-z)                                         11d18s20
      logical ldebug                                                    12d12s19
      dimension ovr(nv,nv)                                              11d18s20
      character*3 scode(4)                                              6d7s23
      include "common.store"
      include "common.print"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data scode/'1ss','3ls','3ms','3hs'/                               6d7s23
      if(iprtr(9).eq.0)then                                             1d6s20
       ldebug=.false.                                                    12d12s19
      else                                                              1d6s20
       ldebug=.true.                                                    1d6s20
      end if                                                            1d6s20
      if(ldebug)write(6,*)('in orthovr ')                               12d12s19
      ibcoffo=ibcoff                                                    8d1s19
      ntri=(nv*(nv+1))/2
      nhere=ntri/mynprocg
      nleft=ntri-mynprocg*nhere
      if(nleft.gt.mynowprog)nhere=nhere+1
      idec=ibcoff
      ibcoff=idec+nhere
      call enough('orthovrb.  1',bc,ibc)
      jdec=idec
      itry=0
      do i=1,nv
       do j=1,i
        jtry=((i*(i-1))/2)+j-1
        if(mod(itry,mynprocg).eq.mynowprog)then
         bc(jdec)=+ovr(j,i)
         jdec=jdec+1
        end if
        itry=itry+1
       end do
      end do
      nok=nv
      if(nder.eq.0)then                                                 8d31s23
       nok=-1
      else                                                              8d31s23
       nderp=nder+1                                                     8d31s23
       ivec=ibcoff                                                      9d1s23
       ieig=ivec+nv*nv*nderp                                            9d1s23
       isym=ieig+nv                                                     9d1s23
       ibcoff=isym+nv                                                   9d1s23
       call enough('orthovrb.  2',bc,ibc)
       call diagxd(nv,ovr,bc(ieig),bc(ivec),ibc(isym),bc,ibc,nder,thrs, 8d31s23
     $      nok,xskip)                                                  9d1s23
       ibcoff=ieig                                                      9d1s23
       if(nok.ne.nv)then
        if(mynowprog.eq.0)write(6,512)isb,scode(l),xskip                9d1s23
       end if                                                           9d1s23
       return                                                           9d1s23
      end if                                                            8d31s23
      if(nok.lt.0)then                                                  1d8s19
       ieig=ibcoff                                                      1d8s19
       ivec=ieig+nv                                                     1d8s19
       isym=ivec+nv*nv                                                  1d8s19
       ibcoff=isym+nv                                                   1d8s19
       call enough('orthovrb.  2',bc,ibc)
       call diagx(nv,ovr,bc(ieig),bc(ivec),ibc(isym),bc,ibc)            11d14s22
       nok=nv                                                           1d8s19
      end if                                                            1d8s19
      if(ldebug)then                                                    12d12s19
       write(6,*)('eigenvalues: ')
       call prntm2(bc(ieig),nok,1,nok)                                     8d1s19
       write(6,*)('vectors ')
       call prntm2(bc(ivec),nv,nv,nok)
      end if                                                            12d12s19
      iskip=0
      do i=0,nok-1
       if(abs(bc(ieig+i)).lt.thrs)iskip=i+1
      end do
      if(iskip.ne.0.and.mynowprog.eq.0)                                 1d8s19
     $     write(6,512)isb,scode(l),bc(ieig+iskip-1)
  512 format('for symmetry ',i1,' coupling ',a3,                        6d7s23
     $     ' largest eigenvalue skipped',es10.2)                        6d7s23
      if(ldebug)call prntm2(bc(ivec),nv,nok,nv)                         12d12s19
      ivec=ivec+nv*iskip
      ieig=ieig+iskip
      nok=nok-iskip
      iveccpy=ivec+nok*nv
      do i=0,nok*nv-1
       bc(iveccpy+i)=bc(ivec+i)
      end do
      do i=0,nok-1
       xx=1d0/sqrt(abs(bc(ieig+i)))
       iad=ivec+nv*i
       do j=0,nv-1
        orig=bc(iad+j)
        bc(iad+j)=bc(iad+j)*xx
       end do
      end do
      ibcoff=iveccpy+nok*nv                                             11d24s20
      return
      end
