c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine intcsfder(dorb,sorb,mdon,mdoo,ibasis,iptr,icsf,nfcn,   8d23s23
     $     ih0,ioooo,ih0d,i4d,pthrs,nec,multh,isymb,nroot,lprint,vec,   8d23s23
     $     nct,irefo,bc,ibc,iptrbit,ixw1,ixw2,ikeep,dorbf,sorbf,nfcnc,  8d23s23
     $     ibasisc,iptrcb,norb,irel,ism,ih0a,dshift,deig,dvec,tolv,     9d16s23
     $     maxdiis,maxrest)                                             10d16s24
      implicit real*8 (a-h,o-z)                                         8d23s23
      logical lprint,ldeb                                               9d1s23
      include "common.store"
      include "common.print"                                            6d20s24
      integer*1 dorb(*),sorb(*),dorbf(*),sorbf(*)                       8d23s23
      dimension ibasis(3,*),iptr(*),icsf(*),ioooo(*),i4d(*),multh(8,8), 8d23s23
     $     vec(nct,*),irefo(*),ikeep(*),iptrcb(2,*),ibasisc(3,*),       8d23s23
     $     iptrbit(2,mdoo+1,*),deig(*),dvec(nct,*),ih0d(*)              8d17s24
      common/singcm/iuse,nff
      ldeb=iprtr(12).ne.0                                               6d20s24
      if(ldeb)write(6,*)('Hi, my name is intcsfder'),isymb,nfcn,nct,tolv9d1s23
      if(nct.eq.1)then                                                  6d26s24
       dvec(1,1)=0d0                                                    6d26s24
       return                                                           6d26s24
      end if                                                            6d26s24
      do i=1,nfcnc                                                      8d23s23
       ikeep(i)=i                                                       8d23s23
      end do                                                            8d23s23
      ihdig=ibcoff                                                      8d23s23
      ibcoff=ihdig+nfcn                                                 8d23s23
      call enough('intcsfder.hdig',bc,ibc)                              8d23s23
      call hiicsf(bc(ihdig),nfcn,dorb,sorb,mdon,mdoo,ibasis,iptr,       8d23s23
     $     ih0a,ioooo,0d0,nec,bc,ibc)                                    8d23s23
      iprod=ibcoff                                                      8d23s23
      ibcoff=iprod+nct*nroot                                            8d23s23
      call enough('intcsfder.prod',bc,ibc)                              8d23s23
      ieig=ibcoff                                                       8d25s23
      ideig=ieig+nroot                                                  8d25s23
      ibcoff=ideig+nroot                                                8d25s23
      call enough('intcsfder.eig',bc,ibc)                               8d25s23
      call hccsfrilr(bc(iprod),nct,vec,nct,ibasis,ibasis,icsf,nfcn,     8d23s23
     $     nfcn,ih0,ioooo,nroot,mdon,nec,multh,nfcnc,ikeep,ibasisc,     8d23s23
     $     mdoo,iptrbit(1,1,isymb),iptrbit(1,1,isymb),ixw1,ixw2,        8d23s23
     $     iptrcb,bc,ibc)                                               8d23s23
      jeig=ieig-1                                                       8d25s23
      do ir=1,nroot                                                     8d23s23
       dot=0d0                                                          8d23s23
       jprod=iprod-1+nct*(ir-1)                                         8d23s23
       do i=1,nct                                                       8d23s23
        dot=dot+bc(jprod+i)*vec(i,ir)
       end do
       if(ldeb)write(6,*)('dot for root '),ir,(' is '),dot              9d1s23
       bc(jeig+ir)=dot                                                  8d25s23
       sz=0d0
       do i=1,nct
        res=bc(jprod+i)-vec(i,ir)*dot                                   8d23s23
        sz=sz+res**2                                                    8d23s23
       end do
       sz=sqrt(sz/dfloat(nct))                                          8d23s23
      end do                                                            8d23s23
      call hccsfrilr(bc(iprod),nct,vec,nct,ibasis,ibasis,icsf,nfcn,     8d23s23
     $     nfcn,ih0d,i4d,nroot,mdon,nec,multh,nfcnc,ikeep,ibasisc,      8d25s23
     $     mdoo,iptrbit(1,1,isymb),iptrbit(1,1,isymb),ixw1,ixw2,        8d23s23
     $     iptrcb,bc,ibc)                                               8d23s23
      jdeig=ideig-1                                                     8d25s23
      do ir=1,nroot                                                     8d25s23
       dot=0d0                                                          8d25s23
       jprod=iprod-1+nct*(ir-1)                                         8d25s23
       do i=1,nct                                                       8d25s23
        dot=dot+bc(jprod+i)*vec(i,ir)                                   8d25s23
       end do                                                           8d25s23
       if(ldeb)write(6,*)('for root '),ir,('energy dervative '),        9d1s23
     $      dot,dot+dshift                                              9d1s23
       deig(ir)=dot+dshift                                              8d29s23
c
c     replace dH*vec with dE*vec-dH*vec
c
       do i=1,nct                                                       8d25s23
        bc(jprod+i)=dot*vec(i,ir)-bc(jprod+i)                           8d25s23
       end do                                                           8d25s23
       bc(jdeig+ir)=dot                                                 8d25s23
      end do                                                            8d25s23
      nwds=nct*nroot                                                    8d25s23
      irhs=iprod
      idv=ibcoff                                                        8d25s23
      ix=idv+nwds                                                       8d25s23
      ir=ix+nwds                                                        4d10s23
      iz=ir+nwds                                                        4d10s23
      iq=iz+nwds                                                        4d10s23
      iqq=iq+nwds                                                       4d10s23
      ipxr=iqq+nwds                                                     4d10s23
      igdiag=ipxr+nwds                                                  4d10s23
      ibeta=igdiag+nwds                                                 4d10s23
      iscale=ibeta+nroot                                                4d3s23
      idotrz=iscale+nroot                                               4d4s23
      idotrznew=idotrz+nroot                                            4d4s23
      ibest1=idotrznew+nroot                                            4d10s23
      ibest2=ibest1+nroot                                               4d10s23
      ibfit=ibest2+nwds                                                 9d26s24
      ibcoff=ibfit+nwds                                                 9d26s24
      call enough('intcsfder.dv',bc,ibc)
      nwdm=nwds-1                                                       8d25s23
      do iz=0,nwdm                                                      8d25s23
       bc(idv+iz)=0d0                                                   8d25s23
      end do                                                            8d25s23
      jbest1=ibest1-1                                                   4d10s23
      do irt=1,nroot                                                     7d20s22
       scale=0d0                                                         4d3s23
       bc(jbest1+irt)=1d10                                              4d10s23
       jrhs=irhs-1+nct*(irt-1)                                          8d25s23
       do i=1,nct                                                       8d25s23
        scale=scale+bc(jrhs+i)**2                                       8d25s23
       end do
       bc(iscale+irt-1)=1d0/max(1d-4,scale)                             4d10s23
      end do                                                            7d20s22
      jhdig=ihdig-1                                                     8d25s23
      jgdiag=igdiag-1                                                   8d25s23
      do i=1,nfcn
       nclo=ibasis(1,i)
       nclop=nclo+1
       iarg=nclop-mdon                                                  8d25s23
       do j=1,icsf(iarg)                                                8d25s23
        bc(jgdiag+j)=bc(jhdig+i)                                        8d25s23
       end do                                                           8d25s23
       jgdiag=jgdiag+icsf(iarg)                                         8d25s23
      end do
      jgdiag=igdiag+nct                                                 8d25s23
      nrootm=nroot-1                                                    8d25s23
      nctm=nct-1                                                        8d25s23
      do irt=1,nrootm                                                   8d25s23
       do i=0,nctm                                                      8d25s23
        bc(jgdiag+i)=bc(igdiag+i)
       end do                                                           8d25s23
       jgdiag=jgdiag+nct                                                8d25s23
      end do                                                            8d25s23
      jgdiag=igdiag                                                     8d25s23
      do irt=1,nroot                                                    7d20s22
       do i=0,nctm                                                      8d25s23
        ip=i+1                                                          8d25s23
        bc(jgdiag+i)=bc(jgdiag+i)-bc(jeig+irt)                          8d25s23
     $       *(1d0-vec(ip,irt)*vec(ip,irt))                             8d25s23
       end do
       jgdiag=jgdiag+nct                                                8d25s23
      end do                                                            8d25s23
      do i=0,nwdm                                                       8d25s23
       bc(igdiag+i)=1d0/bc(igdiag+i)                                    8d25s23
      end do                                                            8d25s23
      best=1d10                                                         9d18s24
      idiis1=ibcoff                                                     4d12s23
      idiis3=idiis1+maxdiis*nwds                                        4d12s23
      ibcoff=idiis3+maxdiis*nwds                                        4d12s23
      call enough('intcsfder.diis',bc,ibc)
      irest=0                                                           4d14s23
      testlast=1d10                                                     9d26s24
 1010 continue
      jdiis1=idiis1-1                                                   4d13s23
      irest=irest+1
      if(lprint.and.irest.ne.1)write(6,*)('restarting ...'),best
      if(irest.gt.maxrest)then                                          10d16s24
       write(6,*)('too many restarts in intcsfder')
       write(6,*)('best so far is '),best                               9d18s24
       write(6,*)('maxdiis used is '),maxdiis                           10d16s24
       write(6,*)('maxrest used is '),maxrest                           10d16s24
       go to 1822
      end if
      do irt=1,nroot                                                    4d13s23
       sz=0d0
       jdv=idv-1+nct*(irt-1)                                            8d25s23
       do i=1,nct                                                       8d25s23
        sz=sz+bc(jdv+i)**2                                              8d25s23
       end do                                                           4d13s23
       if(sz.ne.0d0)then                                                4d13s23
        sz=1d0/sqrt(sz)                                                 4d13s23
        do i=1,nct                                                      8d25s23
         bc(jdiis1+i)=sz*bc(jdv+i)                                      8d25s23
        end do                                                          4d13s23
       else                                                             4d13s23
        sz=0d0                                                          4d13s23
        jgdiag=igdiag-1+nct*(irt-1)                                     4d13s23
        jrhs=irhs-1+nct*(irt-1)                                         8d25s23
        do i=1,nct                                                      8d25s23
         bc(jdiis1+i)=bc(jgdiag+i)*bc(jrhs+i)                           8d25s23
         sz=sz+bc(jdiis1+i)**2                                          4d13s23
        end do                                                          4d13s23
        sz=1d0/sqrt(sz)                                                 4d13s23
        do i=1,nct                                                      8d25s23
         bc(jdiis1+i)=sz*bc(jdiis1+i)                                   4d13s23
        end do                                                          4d13s23
       end if                                                           4d13s23
       jdiis1=jdiis1+nct                                                8d25s23
      end do                                                            4d13s23
      iterm=0                                                           4d13s23
 1000 continue                                                          4d13s23
       iterm=iterm+1                                                    4d13s23
       if(iterm.gt.maxdiis)then                                         4d13s23
        do i=0,nwdm                                                     9d26s24
         bc(irhs+i)=bc(irhs+i)-bc(ibfit+i)                              9d26s24
         bc(idv+i)=bc(idv+i)+bc(ipxr+i)                                 9d26s24
        end do                                                          9d26s24
        go to 1010                                                      4d14s23
       end if                                                           4d13s23
       iin=idiis1+nwds*(iterm-1)                                        4d13s23
       iout=idiis3+nwds*(iterm-1)                                       4d13s23
c
c     H*trial part
c
       call hccsfrilr(bc(iout),nct,bc(iin),nct,ibasis,ibasis,icsf,nfcn, 8d25s23
     $     nfcn,ih0,ioooo,nroot,mdon,nec,multh,nfcnc,ikeep,ibasisc,      8d25s23
     $     mdoo,iptrbit(1,1,isymb),iptrbit(1,1,isymb),ixw1,ixw2,        8d23s23
     $     iptrcb,bc,ibc)                                               8d23s23
c
c     -e*(I-vvT) part
c
       do ir=0,nrootm                                                   8d25s23
        irp=ir+1                                                        8d25s23
        dot=0d0                                                         8d25s23
        jin=iin-1+nct*ir                                                8d25s23
        do i=1,nct                                                      8d25s23
         dot=dot+vec(i,irp)*bc(jin+i)                                   8d25s23
        end do                                                          8d25s23
        jout=iout-1+nct*ir                                              8d25s23
        do i=1,nct                                                      8d25s23
         bc(jout+i)=bc(jout+i)-bc(ieig+ir)*(bc(jin+i)-dot*vec(i,irp))   8d25s23
        end do                                                          8d25s23
       end do                                                           8d25s23
       nfit=iterm                                                       4d13s23
       ifmat=ibcoff                                                     4d13s23
       icoef=ifmat+nct*nfit                                             4d13s23
       iwgt=icoef+nfit                                                   4d12s23
       iscr=iwgt+nct                                                    8d25s23
       ibcoff=iscr+2*nfit+nfit*nfit+2*nfit*nct+nct                      8d25s23
       call enough('intcsfder.ls',bc,ibc)                                   4d12s23
       nok=0                                                             4d13s23
       do irt=0,nrootm                                                   4d12s23
        irp=irt+1                                                       4d13s23
        do i=0,nfit-1                                                     4d12s23
         iad=ifmat+i*nct                                                8d25s23
         iad3=idiis3+nwds*i+nct*irt                                     4d13s23
         do j=0,nctm
          bc(iad+j)=bc(iad3+j)                                          4d13s23
         end do                                                          4d12s23
        end do                                                           4d12s23
        do i=0,nctm                                                     4d13s23
         bc(iwgt+i)=1d0                                                   4d12s23
        end do                                                            4d12s23
        iuse=0
        jrhs=irhs+nct*irt                                               8d25s23
        call lsqfit2(bc(ifmat),nct,nfit,bc(jrhs),nct,1,nct,             8d25s23
     $      bc(icoef),nfit,bc(iscr),bc(iwgt),0,rms,bc,ibc)              4d12s23
        call dgemm('n','n',nct,1,nfit,1d0,bc(ifmat),nct,bc(icoef),      8d25s23
     $      nfit,0d0,bc(iscr),nct,'intcsfder.dgfit')                        4d13s23
        xx=0d0                                                           4d12s23
        jgdiag=igdiag+nct*irt-1                                         4d13s23
        jscr=iscr-1                                                     4d13s23
        jrhs=irhs-1+nct*irt                                             8d25s23
        jbfit=ibfit-1+nct*irt                                           9d26s24
        do i=1,nct                                                      8d25s23
         bc(jbfit+i)=bc(jscr+i)                                         9d26s24
         diff=bc(jrhs+i)-bc(jscr+i)                                        4d12s23
         bc(jscr+i)=bc(jgdiag+i)*diff                                      4d13s23
         xx=xx+diff**2                                                   4d12s23
        end do                                                           4d12s23
        xtest=sqrt(xx*bc(iscale+irt))                                    4d13s23
        if(ldeb)write(6,*)('xx, scale xtest '),xx,bc(iscale+irt),xtest
        call dws_bcast(xtest,1)                                         6d22s23
        best=min(best,xtest)                                            9d18s24
        if(xtest.lt.tolv)then                                           4d13s23
         nok=nok+1
        end if
        do i=0,nfit-1                                                     4d12s23
         iad=ifmat+i*nct                                                8d25s23
         iad1=idiis1+nwds*i+irt*nct                                     8d25s23
         do j=0,nctm                                                    8d25s23
          bc(iad+j)=bc(iad1+j)                                          8d25s23
         end do                                                          4d12s23
        end do                                                           4d12s23
        if(xtest.lt.tolv.or.iterm.eq.maxdiis)then                        4d14s23
         jp=ipxr+nct*irt                                                 4d13s23
         call dgemm('n','n',nct,1,nfit,1d0,bc(ifmat),nct,bc(icoef),      4d13s23
     $       nfit,0d0,bc(jp),nct,'intcsfder.jp')                            4d13s23
        end if                                                           4d13s23
        do i=0,nfit-1                                                    4d12s23
         dot=0d0                                                         4d12s23
         iad=ifmat+nct*i                                                 4d12s23
         do j=0,nctm                                                     4d12s23
          dot=dot+bc(iad+j)*bc(iscr+j)                                   4d12s23
         end do                                                          4d12s23
         do j=0,nctm
          bc(iscr+j)=bc(iscr+j)-dot*bc(iad+j)                            4d12s23
         end do                                                          4d12s23
        end do                                                           4d12s23
        sz=0d0                                                           4d12s23
        do j=0,nctm                                                      4d12s23
         sz=sz+bc(iscr+j)**2                                             4d12s23
        end do
        sz=1d0/sqrt(sz)                                                  4d12s23
        do j=0,nctm                                                      4d12s23
         bc(iscr+j)=bc(iscr+j)*sz                                        4d12s23
        end do                                                           4d12s23
        jdpt=idiis1+nwds*nfit+nct*irt                                    4d13s23
        do i=0,nctm                                                      4d12s23
         bc(jdpt+i)=bc(iscr+i)                                             4d12s23
        end do                                                           4d12s23
       end do                                                            4d12s23
       ibcoff=ifmat                                                      4d12s23
       if(nok.eq.nroot)then                                              4d13s23
        jp=ipxr-1                                                        4d13s23
        do irt=1,nroot                                                   4d13s23
         jdv=idv-1+nct*(irt-1)                                          8d25s23
         do i=1,nct                                                     8d25s23
          bc(jdv+i)=bc(jdv+i)+bc(jp+i)                                  9d26s24
         end do                                                          4d13s23
         jp=jp+nct                                                      8d25s23
        end do                                                           4d13s23
        if(lprint)
     $       write(6,*)('intcsfder diis iterations converged after '),
     $      iterm,('iterations and '),irest-1,('restarts')
       else                                                             4d13s23
        go to 1000                                                      4d13s23
       end if                                                           4d13s23
       if(ldeb)then                                                     9d1s23
        write(6,*)('derivative of vectors: ')
        call prntm2(bc(idv),nct,nroot,nct)
       end if                                                           9d1s23
 1822 continue                                                          9d26s24
      jdv=idv-1                                                         8d29s23
      do ir=1,nroot                                                     8d29s23
       do i=1,nct                                                       8d29s23
        dvec(i,ir)=bc(jdv+i)                                            8d29s23
       end do                                                           8d29s23
       jdv=jdv+nct                                                      8d31s23
      end do                                                            8d29s23
      ibcoff=ihdig                                                      8d23s23
      return                                                            8d23s23
      end                                                               8d23s23
