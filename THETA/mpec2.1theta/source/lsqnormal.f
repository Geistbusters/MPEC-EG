c mpec2.1 version theta. For copyright and Disclaimers, see start.f
      subroutine lsqnormal(fmat,idf,nf,data,idd,nd,npts,coef,idc,iret,  1d28s25
     $     aratio,bc,ibc)                                               1d28s25
      implicit real*8 (a-h,o-z)
c
c     solve least square problem via normal equations
c     iret is zero for statisfactory completion, nonzero if
c     normal equations are singular.                                    1d28s25
c
      dimension fmat(idf,*),data(idd,*),coef(idc,*)                     1d28s25
      include 'common.store'                                            1d28s25
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      data loop,loopx/0,10/
      save
      iret=0                                                            1d27s25
      ibcoffo=ibcoff                                                    1d27s25
      call ilimts(npts,1,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)      1d27s25
      nhere=ih+1-il                                                     1d27s25
      itrans=ibcoff                                                     1d24s25
      iamat=itrans+nhere*nf                                             1d27s25
      ibmat=iamat+nf*nf                                                 1d24s25
      ipvt=ibmat+nf*nd                                                  1d24s25
      ibcoff=ipvt+nf                                                    1d24s25
      call enough('lsqnormal.trans',bc,ibc)                               1d24s25
      do iz=iamat,ipvt-1                                                1d27s25
       bc(iz)=0d0                                                       1d27s25
      end do                                                            1d27s25
      if(nhere.gt.0)then                                                1d27s25
       ilm=il-1                                                         1d27s25
       do i=1,nf                                                         1d24s25
        do j=1,nhere                                                     1d27s25
         ij=itrans+i-1+nf*(j-1)                                          1d24s25
         bc(ij)=fmat(j+ilm,i)                                           1d27s25
        end do                                                           1d24s25
       end do                                                            1d24s25
       call dgemm('n','n',nf,nf,nhere,1d0,bc(itrans),nf,fmat(il,1),idf, 1d27s25
     $      0d0,bc(iamat),nf,'lsqnormal.amat')                            1d27s25
       call dgemm('n','n',nf,nd,nhere,1d0,bc(itrans),nf,data(il,1),idd, 1d27s25
     $      0d0,bc(ibmat),nf,'lsqnormal.bmat')                            1d27s25
      end if                                                            1d27s25
      ngs=nf*(nf+nd)                                                    1d27s25
      call dws_gsumf(bc(iamat),ngs)                                     1d27s25
      if(bc(iamat).le.1d-20)bc(iamat)=1d0                               2d14s25
      ianorm=ibcoff                                                     1d27s25
      ibcoff=ianaorm+nf                                                 1d27s25
      a1=1d20                                                           1d27s25
      ii=iamat+nf*nf-1
      an=-1d20                                                          1d27s25
      do i=0,nf-1
       ii=iamat+i*(nf+1)                                                1d27s25
       a1=min(a1,bc(ii))                                                1d27s25
       an=max(an,bc(ii))                                                1d27s25
       bc(ianorm+i)=1d0/sqrt(bc(ii))                                    1d27s25
      end do                                                            1d27s25
      aratio=an/a1                                                      1d27s25
      if(aratio.gt.1d10)then                                            1d28s25
       iret=1                                                           1d28s25
       ibcoff=ibcoffo                                                   1d28s25
       return                                                           1d28s25
      end if                                                            1d28s25
      do i=0,nf-1                                                       1d27s25
       jamat=iamat+nf*i                                                 1d27s25
       do j=0,nf-1                                                      1d27s25
        bc(jamat+j)=bc(jamat+j)*bc(ianorm+i)*bc(ianorm+j)               1d27s25
       end do                                                           1d27s25
      end do                                                            1d27s25
      do i=0,nd-1                                                       1d27s25
       jbmat=ibmat+nf*i                                                 1d27s25
       do j=0,nf-1                                                      1d27s25
        bc(jbmat+j)=bc(jbmat+j)*bc(ianorm+j)                            1d27s25
       end do                                                           1d27s25
      end do                                                            1d27s25
      call lusolv(bc(iamat),nf,nf,bc(ibmat),nf,nd,ibc(ipvt),ierr,3)     1d24s25
      do i=0,nd-1                                                       1d27s25
       jbmat=ibmat+nf*i                                                 1d27s25
       do j=0,nf-1                                                      1d27s25
        bc(jbmat+j)=bc(jbmat+j)*bc(ianorm+j)                            1d27s25
       end do                                                           1d27s25
      end do                                                            1d27s25
      do i=1,nd                                                         1d28s25
       jbmat=ibmat-1+nf*(i-1)                                           1d28s25
       do j=1,nf                                                        1d28s25
        coef(j,i)=bc(jbmat+j)                                           1d28s25
       end do                                                           1d28s25
      end do                                                            1d28s25
      ibcoff=ibcoffo                                                    1d28s25
      return                                                            1d28s25
      end                                                               1d28s25
