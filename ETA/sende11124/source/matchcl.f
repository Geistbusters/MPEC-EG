c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine matchcl(inst,ias,ndet,nopen,isen,ndetj,jas,nopenj,vdet,9d12s19
     $     vtp,v2t,nrootu,nfcn,vout,iextra,lcrenorm,vnorm,fact3jlt,     9d28s20
     $     fact3jst2,fact3jht2)                                         9d28s20
      implicit real*8 (a-h,o-z)                                         9d12s19
      integer*1 ias(nopen,*),jas(nopenj,*),iorb(64)                     9d12s19
      logical ldebug                                                    12d12s19
      dimension inst(2),vdet(max(1,ndet),nrootu),vtp(ndetj,nrootu),     11d4s19
     $     v2t(nfcn,ndetj),vout(nfcn,nrootu),vnorm(*)                   1d22s20
      include "common.print"                                            1d5s20
      if(iprtr(9).eq.0)then                                             1d5s20
       ldebug=.false.                                                    12d12s19
      else                                                              1d5s20
       ldebug=.true.                                                    1d5s20
      end if                                                            1d5s20
      if(ldebug)then                                                    12d12s19
       write(6,*)('in matchcl ')
       write(6,*)('what we have for jas: ')
       do i=1,ndetj
        write(6,*)i,(':'),(jas(j,i),j=1,nopenj)
       end do
      end if                                                            12d12s19
      do i=1,nrootu                                                     9d12s19
       do j=1,ndetj
        vtp(j,i)=0d0                                                    9d12s19
       end do                                                           9d12s19
      end do                                                            9d12s19
      do idet=1,max(1,ndet)                                             11d4s19
       if(ldebug)                                                       12d12s19
     $      write(6,*)('for det '),idet,(':'),(ias(j,idet),j=1,nopen),  12d12s19
     $      (vdet(idet,j),j=1,nrootu)                                   11d18s19
       if(isen.eq.2)then
        nb=inst(2)-inst(1)+1                                            9d13s19
c
c     cases:
c     a: inst(1) le nopen, inst(2) le nopen
c     b: inst(1) le nopen, inst(2) gt nopen
c     c: inst(1) gt nopen, inst(2) gt nopen
        jj=0                                                            9d10s19
        do i=1,min(inst(1)-1,nopen)                                     9d10s19
         jj=jj+1                                                        9d10s19
         iorb(jj)=ias(i,idet)                                           9d10s19
        end do                                                          9d10s19
c     in all cases, inst(1) is next
        jj=jj+1                                                         9d10s19
        iorb(jj)=-1                                                     9d12s19
        do i=inst(1),min(inst(2)-1,nopen)                               9d10s19
         jj=jj+1                                                        9d10s19
         iorb(jj)=ias(i,idet)                                           9d10s19
        end do                                                          9d10s19
        jj=jj+1                                                         9d10s19
        iorb(jj)=-1                                                     9d12s19
        do i=inst(2),nopen                                              9d10s19
         jj=jj+1                                                        9d10s19
         iorb(jj)=ias(i,idet)                                           9d10s19
        end do                                                          9d10s19
       else if(isen.eq.3)then
c
c     inst(1) is open orbital that is excited out of
c     inst(2) is where closed orbital is stuck in.
c
        nb=2*nopen+1-inst(1)-inst(2)+iextra                             11d18s19
        if(inst(2).gt.inst(1))nb=nb+1                                   11d18s19
        if(ias(inst(1),idet).ne.1)go to 2                               9d12s19
        jj=0                                                            9d10s19
        do i=1,min(inst(2)-1,nopen)                                     9d10s19
         if(i.ne.inst(1))then                                           9d10s19
          jj=jj+1                                                       9d10s19
          iorb(jj)=ias(i,idet)                                          9d10s19
         end if                                                         9d10s19
        end do                                                          9d10s19
        jj=jj+1                                                         9d10s19
        iorb(jj)=-1                                                     9d12s19
        ihit=0                                                          9d10s19
        do i=inst(2),nopen                                              9d10s19
         if(i.ne.inst(1))then                                           9d10s19
          ihit=1                                                        9d10s19
          jj=jj+1                                                       9d10s19
          iorb(jj)=ias(i,idet)                                          9d10s19
         end if                                                         9d10s19
        end do                                                          9d10s19
        if(inst(2).gt.nopen.and.ihit.eq.1)then                          9d10s19
         jj=jj+1                                                         9d10s19
         iorb(jj)=-1                                                    9d12s19
        end if                                                          9d10s19
       else                                                             9d12s19
        if(ias(inst(1),idet).lt.0.or.                                    9d12s19
     $     ias(inst(2),idet).lt.0)go to 2                                9d12s19
        jj=0                                                            9d10s19
        do i=1,nopen                                                    9d10s19
         if(i.ne.inst(1).and.i.ne.inst(2))then                          9d10s19
          jj=jj+1                                                       9d10s19
          iorb(jj)=ias(i,idet)                                          9d10s19
         end if                                                         9d10s19
        end do                                                          9d10s19
        nb=nopen-inst(2)+nopen-1-inst(1)                                9d10s19
       end if                                                           9d12s19
       if(ldebug)write(6,*)('promoted det: '),(iorb(ij),ij=1,jj)        12d12s19
       if(jj.ne.nopenj)then                                             9d10s19
        write(6,*)('jj = '),jj,(' ne nopenj = '),nopenj
        call dws_sync
        call dws_finalize
        stop 'matchcl'                                                   9d10s19
       end if                                                           9d10s19
       do i=1,ndetj                                                     9d10s19
        do j=1,nopenj                                                   9d10s19
         if(iorb(j).ne.jas(j,i))go to 1                                 9d10s19
        end do                                                          9d10s19
        if(ldebug)then                                                  12d12s19
         write(6,*)('matches det no. '),i                                9d10s19
         write(6,*)('nb for phase '),nb
        end if                                                          12d12s19
        if(mod(nb,2).eq.0)then                                          9d12s19
         ff=1d0                                                         9d12s19
        else                                                            9d12s19
         ff=-1d0                                                        9d12s19
        end if                                                          9d12s19
        do j=1,nrootu                                                   9d12s19
         prod=ff*vdet(idet,j)                                           9d12s19
         vtp(i,j)=vtp(i,j)+prod                                         9d12s19
        end do                                                          9d12s19
        go to 2                                                         9d10s19
    1   continue                                                        9d10s19
       end do                                                           9d10s19
       write(6,*)('no match found!!! ')                                 9d10s19
       call dws_sync                                                    9d10s19
       call dws_finalize                                                9d10s19
       stop 'matchc'                                                    9d10s19
    2  continue                                                         9d10s19
      end do                                                            9d12s19
      if(ldebug)then                                                    12d12s19
       write(6,*)('new vectors in det basis: ')
       call prntm2(vtp,ndetj,nrootu,ndetj)                               9d12s19
       write(6,*)('mult by ')
       call prntm2(v2t,nfcn,ndetj,nfcn)
      end if                                                            12d12s19
      call dgemm('n','n',nfcn,nrootu,ndetj,1d0,v2t,nfcn,vtp,ndetj,0d0,  9d12s19
     $     vout,nfcn,                                                   9d12s19
     d' matchcl.  1')
      if(ldebug)then                                                    12d12s19
       write(6,*)('in csf basis ')
       call prntm2(vout,nfcn,nrootu,nfcn)                                9d12s19
      end if                                                            12d12s19
      return                                                            9d12s19
      end                                                               9d12s19
