c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine lsqfit2(fmat,idf,nf,data,idd,nd,npts,coef,idc,scr,wgt, 10d17s91
     $                   iunit,rmsovr,bc,ibc)                           3d27s23
      implicit real*8 (a-h,o-z)
c
c     perform linear least squares fit.
c     input data:
c     fmat(i,j) - function j at data point i
c     idf - first dimension of fmat                                     10d17s91
c     nf - number of functions
c     data(i,k) - ordinate at data point i for fit k                    10d17s91
c     idd - first dimension of data                                     10d17s91
c     nd - number of data sets to fit                                   10d17s91
c     npts - number of data points
c     idc - first dimension of coef                                     10d17s91
c     wgt(i) - weight for data point i                                  10d13s89
c     iunit - unit to write results                                     10d16s89
c     output data:
c     coef(j,k) - the least squares coefficient for function j for fit k10d17s91
c     rmsovr - overall rms absolute unweighted error                    10d17s91
c     scratch arrays:
c     scr(2*nf+nf*nf+2*nf*npts+npts),
c
      dimension fmat(idf,nf),data(idd,nd),coef(idc,nd),scr(1),wgt(npts) 10d17s91
      dimension isort(3000)
      include "common.store"                                            3d27s23
      common/singcm/iuse,nff
      nff=nf
c
c     perform singular value decomposition of fmat
c
      is=1
      ie=is+nf
      iu=ie+nf
      iv=iu+npts*nf
      iw=iv+nf*nf
      icp=iw+npts
      do 50 i=1,nf
       do 51 j=1,npts
        scr(icp+j-1+(i-1)*npts)=fmat(j,i)*wgt(j)                        10d13s89
   51  continue
   50 continue
      lwork=5*max(npts,nf)                                              3d27s23
      iwork=ibcoff                                                      3d27s23
      ibcoff=iwork+lwork                                                3d27s23
      call enough('lsqfit2.work',bc,ibc)                                3d27s23
      call dgesvd('S','S',npts,nf,scr(icp),npts,scr(is),scr(iu),npts,
     $     scr(iv),nf,bc(iwork),lwork,info)
      itmx=ibcoff                                                       3d27s23
      ibcoff=itmx+nf*nf                                                 3d27s23
      call enough('lsqfit2.tmx',bc,ibc)                                 3d27s23
      do i=0,nf-1                                                       3d27s23
       do j=0,nf-1                                                      3d27s23
        ji=iv+j+nf*i                                                    3d27s23
        ij=itmx+i+nf*j                                                  3d27s23
        bc(ij)=scr(ji)                                                  3d27s23
       end do                                                           3d27s23
      end do                                                            3d27s23
      do i=0,nf*nf-1                                                    3d27s23
       scr(iv+i)=bc(itmx+i)                                             3d27s23
      end do                                                            3d27s23
      ibcoff=iwork                                                      3d27s23
      ibcoff=iwork                                                      3d27s23
c
c     determine number of singular values used
c
      thrs=scr(is)*2.7d-9
      if(iunit.ne.0)call prntm2(scr(is),1,nf,1)
      if(iunit.ne.0)write(iunit,5001)scr(is),scr(is+nf-1)               9d7s93
 5001 format(1x,'range of singular values ',1p2e15.7)
      call flush(6)
      if(iuse.le.0)then                                                 10d15s94
      do 1492 i=1,nf                                                    10d6s00
       scr(iw+i-1)=0d0                                                  10d6s00
 1492 continue                                                          10d6s00
      do 1 i=1,nf
       if(scr(is+i-1).gt.thrs)then                                      11d24s89
        iuse=i                                                          11d24s89
       else if(iunit.ne.0)then                                          9d7s93
        write(iunit,153)                                                11d24s89
  153   format(1x,'column of v cooresponding to small singular value ') 11d24s89
        do 154 j=1,nf                                                   11d24s89
         iadd=iv+(i-1)*nf+j-1                                           11d24s89
         scr(iw+j-1)=scr(iw+j-1)+scr(iadd)**2
         if(abs(scr(iadd)).gt.1d-8)write(iunit,155)j,scr(iadd)          11d24s89
  155    format(1x,i5,1pe15.7)                                          11d24s89
  154   continue                                                        11d24s89
       end if                                                           11d24s89
    1 continue
      if(iuse.ne.nf.and.iunit.ne.0)then                                 9d7s93
       write(iunit,2)nf,iuse                                            10d16s89
    2  format(/1x,'out of ',i5,' functions, only ',i5,' singular',
     $        ' values are used')
       if(nf.gt.3000)then
        write(6,*)('norm of thrown away part ')
        call prntm2(scr(iw),nf,1,nf)
       else
        call dsortdws(scr(iw),isort,nf)                                 1d18s23
        write(6,*)('sorted ')
        do 1493 i=1,nf
         write(6,1494)i,scr(iw+i-1),isort(i)
 1494    format(i5,1pe15.7,2i5)
 1493   continue
       end if
      end if
      end if                                                            10d15s94
c
c     solve for coefficients
c
      rmsovr=0d0
      do 1001 ift=1,nd                                                  10d17s91
       do 3 i=1,iuse
        scr(ie+i-1)=0d0
    3  continue
       do 4 i=1,iuse
        do 5 j=1,npts
         scr(ie+i-1)=scr(ie+i-1)                                        10d17s91
     $        +scr(iu+j-1+(i-1)*npts)*data(j,ift)*wgt(j)                10d17s91
    5   continue
    4  continue
       if(iunit.ne.0)write(iunit,*)
     $      ('coefficients in left singular vector basis ')
       do 6 i=1,iuse
        scr(ie+i-1)=scr(ie+i-1)/scr(is+i-1)
        if(iunit.ne.0)write(iunit,1006)i,scr(ie+i-1)
 1006   format(i5,1pe21.14)
    6  continue
       do 7 i=1,nf
        coef(i,ift)=0d0
    7  continue
       do 8 j=1,iuse
        do 9 i=1,nf
         coef(i,ift)=coef(i,ift)+scr(iv+i-1+(j-1)*nf)*scr(ie+j-1)
    9   continue
    8  continue
       if(iunit.ne.0)then                                               9d7s93
        write(iunit,*)('coefficients for fit no. '),ift
        write(iunit,22)(coef(i,ift),i=1,nf)
   22   format(1x,1p5e15.7)
       end if                                                           9d7s93
c
c     compute various properties of the solution.
c
       rms=0d0
       rmsw=0d0                                                         10d13s89
       rmsx=0d0                                                         10d25s94
       bott=0d0                                                         10d13s89
       if(iunit.ne.0)write(iunit,21)                                    9d7s93
   21  format(1x,' no.       data         fit         diff    ',        10d5s89
     $        '   diff/data         wgt   diff*wgt')                    10d25s94
       do 11 i=1,npts
        bott=bott+wgt(i)                                                10d13s89
        fit=0d0
        do 12 j=1,nf
         fit=fit+fmat(i,j)*coef(j,ift)                                  10d17s91
   12   continue
        diff=data(i,ift)-fit                                            10d17s91
        if(data(i,ift).ne.0d0)then                                      10d17s91
         rat=diff/data(i,ift)                                           10d17s91
        else                                                            10d5s89
         rat=0.d0                                                       10d5s89
        end if                                                          10d5s89
        if(iunit.ne.0)write(iunit,20)i,data(i,ift),fit,diff,rat,wgt(i), 10d25s94
     $       diff*wgt(i)                                                10d25s94
   20   format(1x,i5,1p3e15.7,0p3f10.5)                                 10d13s89
        rms=rms+(data(i,ift)-fit)**2                                    10d17s91
        rmsw=rmsw+(wgt(i)*diff)**2                                      10d13s89
   11  continue
       rmsx=sqrt(rmsw/dfloat(npts))                                     10d25s94
       rmsovr=rmsovr+rmsw
       rms=sqrt(rms/dfloat(npts))
       rmsw=sqrt(rmsw/bott)                                             10d13s89
       if(iunit.ne.0)write(iunit,13)rmsw,rms,rmsx                       10d25s94
   13  format(/1x,'rms error of fit ',1pe15.7,                          10d13s89
     $        /1x,'unweighted rms error ',e15.7,                        10d25s94
     $        /1x,'avg weighted error ',e15.7)                          10d25s94
 1001 continue
      rmsovr=sqrt(rmsovr/dfloat(nd*npts))
      if(iunit.ne.0)write(iunit,*)('overall rms error = '),rmsovr
      scr(1)=rmsw                                                       10d15s94
      return
      end
