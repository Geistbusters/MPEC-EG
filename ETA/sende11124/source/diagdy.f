c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine diagdy(hdec,n,nwant,ieig,ivec,bc,ibc)                  11d14s22
      implicit real*8 (a-h,o-z)
c
c     take decimated triangle of hamiltonian matrix and
c     create square matrix to diagonalize separately                    3d19s21
c
      dimension hdec(*)
      include "common.store"
      include "common.print"                                            4d28s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      ieig=ibcoff                                                       3d19s21
      ivec=ieig+nwant                                                   3d19s21
      ibcoff=ivec+n*nwant                                               3d19s21
      itmp=ibcoff
      ibcoff=itmp+n*n
      call enough('diagdy.  1',bc,ibc)
      do i=0,n*n-1
       bc(itmp+i)=0d0
      end do
      ii=0
      jj=1
      do i=0,n-1
       do j=0,i
        if(mod(ii,mynprocg).eq.mynowprog)then
         ji=itmp+j+n*i
         ij=itmp+i+n*j
         bc(ji)=hdec(jj)
         bc(ij)=hdec(jj)
         jj=jj+1
        end if
        ii=ii+1
       end do
      end do
      call dws_gsumf(bc(itmp),n*n)
      if(iprtr(12).ne.0)then                                            4d28s21
       write(6,*)('squared hamiltonian matrix ')
       call prntm2(bc(itmp),n,n,n)
       itcopy=ibcoff
       ibcoff=itcopy+n*n
       call enough('diagdy.  2',bc,ibc)
       do i=0,n*n-1
        bc(itcopy+i)=bc(itmp+i)                                         3d22s22
       end do                                                           3d22s22
      end if                                                            4d28s21
      lwork=100*n                                                       3d19s21
      kwork=ibcoff                                                      3d19s21
      iwork=kwork+lwork                                                 3d19s21
      ifail=iwork+n*5
      ibcoff=ifail+n                                                    3d19s21
      call enough('diagdy.  3',bc,ibc)
      call dsyevx('V','I','L',n,bc(itmp),n,dum,dum,1,nwant,0d0,         3d19s21
     $     nfound,bc(ieig),bc(ivec),n,bc(kwork),lwork,ibc(iwork),       3d19s21
     $     ibc(ifail),info)                                             3d19s21
      if(info.ne.0)then                                                 3d19s21
       write(6,*)('error code from dsyevx in diagdy!! '),info,n
       call prntm2(bc(itcopy),n,n,n)
       stop
      end if
      call phasv(bc(ivec),n,n,nwant)                                    4d16s21
      nsend=nwant*(n+1)                                                 3d19s21
      call dws_bcast(bc(ieig),nsend)                                    3d19s21
      if(iprtr(12).ne.0)then
       write(6,*)('H times all but first element ...')
       nwant2=nwant*2
       ivtmp=ibcoff
       ibcoff=ivtmp+n*nwant2
       call enough('diagdy.  4',bc,ibc)
       do i=0,nwant-1
        iadi=ivtmp+n*i
        do k=0,n-1
         bc(iadi+k)=0d0
        end do
        bc(iadi+i)=1d0
        iadi=iadi+n*nwant
        ivad=ivec+n*i
        do k=0,n-1
         bc(iadi+k)=bc(ivad+k)
        end do
       end do
       write(6,*)('trial vectors ')
       call prntm2(bc(ivtmp),n,nwant2,n)                                3d22s22
       do i=0,nwant2-1
        iadi=ivtmp+n*i
        do j=0,i-1
         iadj=ivtmp+n*j
         dot=0d0
         do k=0,n-1
          dot=dot+bc(iadi+k)*bc(iadj+k)
         end do
         do k=0,n-1
          bc(iadi+k)=bc(iadi+k)-dot*bc(iadj+k)
         end do
        end do
        dot=0d0
        do k=0,n-1
         dot=dot+bc(iadi+k)**2
        end do
        dot=1d0/sqrt(dot)
        do k=0,n-1
         bc(iadi+k)=bc(iadi+k)*dot
        end do
       end do
       write(6,*)('after schmidt ...')
       call prntm2(bc(ivtmp),n,nwant2,n)
       call dgemm('n','n',n,nwant2,n,1d0,bc(itcopy),n,bc(ivtmp),n,
     $      0d0,bc(itmp),n,
     d' diagdy.  1')
       do i=0,nwant2-1
        do j=0,n-1
         ji=itmp+j+n*i
         ij=itcopy+i+nwant2*j
         bc(ij)=bc(ji)
        end do
       end do
       call dgemm('n','n',nwant2,nwant2,n,1d0,bc(itcopy),nwant2,
     $      bc(ivtmp),n,0d0,bc(itmp),nwant2,
     d' diagdy.  2')
       call prntm2(bc(itmp),nwant2,nwant2,nwant2)
       do i=0,n*nwant2-1
        bc(itmp+i)=bc(ivtmp+i)
       end do
      end if
      ibcoff=itmp
      return
      end
