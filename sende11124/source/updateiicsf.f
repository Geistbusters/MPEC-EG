c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine updateiicsf(vec,eig,nroot,ncsft,hdig,nfcn,ibasis,ncsf,
     $     mdon,nps,ipointf,isto,vecn,bc,ibc)                           11d10s22
      implicit real*8 (a-h,o-z)                                         6d22s18
c                                                                       6d22s18
c     davidson update. on input, g contains h*c. on output, it          6d22s18
c     will contain orthogonal component of best guess for best vector   6d22s18
c
      logical lprint                                                    6d26s18
      dimension vec(ncsft,nroot),eig(nroot),hdig(nfcn),vecn(*),         7d12s19
     $     ibasis(3,*),ncsf(*),ipointf(*)                               7d12s19
      include "common.basis"                                            7d10s19
      include "common.input"                                            7d10s19
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data icall/0/
      save
      icall=icall+1
      ibcoffo=ibcoff
      lprint=.false.                                                    6d26s18
      if(lprint)then                                                    6d26s18
       write(6,*)('in updateii ')
       write(6,*)('input vector ')
       call prntm2(vec,ncsft,nroot,ncsft)
       write(6,*)('input hc ')
       call prntm2(vec(1,nroot+1),ncsft,nroot,ncsft)
       write(6,*)('input eigenvalues ')
       call prntm2(eig,1,nroot,1)
      end if                                                            6d26s18
c
c     we are using spin independent diagonals ...
c
      ips=0                                                             7d12s19
      do if=1,nfcn                                                      7d12s19
       nclo=ibasis(1,if)                                                7d12s19
       nclop=nclo+1                                                     7d12s19
       iarg=nclop-mdon                                                  7d12s19
       do i=1,ncsf(iarg)                                                7d12s19
        ip=ips+i                                                        7d12s19
        do j=1,nroot                                                    7d12s19
         resid=vec(ip,j+nroot)-eig(j)*vec(ip,j)                         7d12s19
         bot=hdig(if)-eig(j)                                            7d12s19
         bot=bot+1d-14                                                   6d22s18
         update=-resid/bot                                               6d22s18
         vec(ip,j+nroot)=update                                         7d12s19
        end do                                                           6d22s18
       end do                                                            6d22s18
       ips=ips+ncsf(iarg)
      end do                                                            7d12s19
      if(lprint)then                                                    6d26s18
       write(6,*)('trial update vector ')
       call prntm2(vec(1,nroot+1),ncsft,nroot,ncsft)                    7d12s19
      end if                                                            6d26s18
c
c     zero out pspace components
c
      do i=1,nps
       iad=ipointf(i)                                                   7d12s19
       if(lprint)write(6,*)('p-space configuration '),i,('is csf no. '),7d12s19
     $      iad                                                         7d12s19
       do j=1,nroot
        vec(iad,j)=0d0                                                  7d12s19
        vec(iad,j+nroot)=0d0                                            7d12s19
       end do
      end do
      if(lprint)then                                                    6d26s18
       write(6,*)('new vec: ')
       call prntm2(vec,ncsft,nroot,ncsft)                               7d12s19
       write(6,*)('new g: ')
       call prntm2(vec(1,nroot+1),ncsft,nroot,ncsft)                    7d12s19
      end if                                                            6d26s18
c
c     now orthogonalize
c
c     as one off, modified schmidt
c
      isto=0                                                            3d25s21
      do i=1,nroot*2                                                    3d25s21
       do j=1,i-1                                                         3d25s21
        dot=0d0                                                         3d25s21
        do k=1,ncsft                                                    3d25s21
         dot=dot+vec(k,i)*vec(k,j)                                      3d25s21
        end do                                                          3d25s21
        do k=1,ncsft                                                    3d25s21
         vec(k,i)=vec(k,i)-dot*vec(k,j)                                 3d25s21
        end do                                                          3d25s21
       end do                                                           3d25s21
       dot=0d0                                                          3d25s21
       do k=1,ncsft                                                     3d25s21
        dot=dot+vec(k,i)**2                                             3d25s21
       end do                                                           3d25s21
       if(dot.gt.smallest)then                                          12d9s22
        isto=isto+1                                                     3d25s21
        doti=1d0/sqrt(dot)                                               3d25s21
        do k=1,ncsft                                                     3d25s21
         vec(k,isto)=doti*vec(k,i)                                      3d25s21
        end do                                                           3d25s21
       end if                                                           3d25s21
      end do                                                            3d25s21
      xsto=dfloat(isto)                                                 4d16s21
      call dws_bcast(xsto,1)                                            4d16s21
      isto=nint(xsto)                                                   4d16s21
      call dws_bcast(vec,ncsft*isto)                                    4d16s21
      ibcoff=ibcoffo
      return
      end
