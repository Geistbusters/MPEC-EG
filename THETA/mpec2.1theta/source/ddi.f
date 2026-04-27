cddi1c     implement ddi using mpi commands, call by fortran
cddi1c     half of procs will be compute procs, the other half will
cddi1c     be memory procs, so ddi will work across nodes.
cddi1c     0 and even nodes will be compute nodes.
cddi1c     itag=1 for proc node to memory node directive
cddi1c     other values for memory node to proc node.
cddi1c
cddi1      subroutine dws_preinit                                            1d31s21
cddi1      use mpi
cddi1      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
cddi1     $     mynnode
cddi1      call mpi_init(ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in ddi_int,')
cddi1       write(6,*)('mpi_init returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop
cddi1      end if
cddi1c$$$      write(6,*)('back from mpi_init ')
cddi1      call mpi_comm_size(MPI_COMM_WORLD,isize,ierror)
cddi1      call mpi_comm_rank(MPI_COMM_WORLD,irank,ierror)
cddi1      if(irank.eq.0)write(6,*)('Hi, this is preinit for MPI-1 code')    8d29s22
cddi1      mynowprog=irank                                                   1d31s21
cddi1      mynprocg=isize                                                    1d31s21
cddi1      return                                                            1d31s21
cddi1      end                                                               1d31s21
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine dws_bcast(buf,len)
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer (i-n)
cddi1      integer len
cddi1      integer icont
cddi1      dimension buf(len)
cddi1      include "common.mympi"                                            1d29s21
cddi1      data icall/0/
cddi1      icall=icall+1
cddi1      icont=len
cddi1      isource=0
cddi1      call mpi_bcast(buf,icont,mpi_double_precision,isource,
cddi1     $               my_comm_group,ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in ddi_bcast,')
cddi1       write(6,*)('mpi_bcast returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop
cddi1      end if
cddi1      return
cddi1      end
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine dws_bcasta(buf,len)                                    1d31s21
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer (i-n)
cddi1      integer len
cddi1      integer icont
cddi1      dimension buf(len)
cddi1      icont=len
cddi1      isource=0
cddi1      call mpi_bcast(buf,icont,mpi_double_precision,isource,
cddi1     $               mpi_comm_world,ierror)                             1d31s21
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in ddi_bcasta,')
cddi1       write(6,*)('mpi_bcast returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop
cddi1      end if
cddi1      return
cddi1      end
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine dws_init(ncore,nmn,bc,ibc)                             11d15s22
cddi1      use mpi
cddi1c
cddi1c     mxcproc is the maximum number of compute procs
cddi1c
cddi1      parameter (idm=10)                                                2d21s21
cddi1      implicit integer (i-n)
cddi1      implicit real*8 (a-h,o-z)
cddi1      integer*8 iarg1,iarg2,ddi_np,ddi_me                               11d18s20
cddi1      character*1 file(8)
cddi1      character*8 fname
cddi1      character*255 procname                                            1d29s21
cddi1      character*101 line                                                2d21s21
cddi1      include "common.mympi"                                            1d29s21
c$$$cddi1      dimension npacket(8),istatus(mpi_status_size),
c$$$cddi1     $     mdata(4,id),packet(4)
cddi1      dimension npacket(8),istatus(mpi_status_size),packet(4)           1d26s23
cddi1      dimension minfo(7,idm,2),xinfo(5,idm,2),ninfo(2),nwrt(10,2)       2d21s21
cddi1      equivalence (packet,npacket)
cddi1      data file/'m','o','u','t','.','_','_','_'/
cddi1      data nwrt/7,3,7,7,7,2,3,2,3,7,                                    2d21s21
cddi1     $          0,0,1,3,1,0,0,0,0,3/                                    2d21s21
cddi1      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
cddi1     $     mynnode
cddi1      include "common.store"
cddi1      save
cddi1      data ncall/0/                                                     2d8s21
cddi1      call mpi_comm_size(MPI_COMM_WORLD,isize,ierror)
cddi1      call mpi_comm_rank(MPI_COMM_WORLD,irank,ierror)
cddi1      call mpi_get_processor_name(procname,lenname,ierror)
cddi1c
cddi1c     if nmn = 0, we will use no memory nodes.                          8d12s22
cddi1c     if nmn = 1,                                                       8d12s22
cddi1c     I assume we are hyperthreading with ncore
cddi1c     physical cores per node. Then we have the mapping
cddi1c
cddi1c     node 0
cddi1c     phsyical core 0          irank 0,ncore
cddi1c                   1          irank 1,ncore+1
cddi1c        .
cddi1c        .
cddi1c        .
cddi1c     phsyical core ncore-1    irank ncore-1,ncore*2-1
cddi1c
cddi1c     node 1
cddi1c     physical core 0          irank 2*ncore,3*ncore
cddi1c
cddi1c     physical core ncore-1    irank 3*ncore-1,4*ncore-1
cddi1c
cddi1c
cddi1c     node n
cddi1c     physical core m          irank n*2*ncore+m,"+ncore
cddi1c
cddi1      ncore2=ncore*2
cddi1      if(ncore2.gt.isize)then                                           2d8s21
cddi1       ncore2=isize                                                     2d8s21
cddi1       ncore=ncore2/2                                                   2d8s21
cddi1      end if                                                            2d8s21
cddi1      nu=0                                                              1d29s21
cddi1      ncorem=ncore-1                                                    2d8s21
cddi1      do irankx=0,isize-1                                               1d29s21
cddi1       node=irankx/ncore2
cddi1       n0phys=irankx-node*ncore2
cddi1       if(n0phys.lt.ncore.or.nmn.eq.0)then                              8d12s22
cddi1        nu=nu+1
cddi1        if(nu.gt.mxcproc)stop 'mxcproc'
cddi1        iranks(nu)=irankx
cddi1       end if
cddi1      end do
cddi1      node=irank/ncore2
cddi1      n0phys=irank-node*ncore2
cddi1      if(n0phys.ge.ncore.and.nmn.ne.0)then                              8d12s22
cddi1       immem=1
cddi1      else
cddi1       immem=0
cddi1      end if
cddi1      mnmc=nmn                                                          8d12s22
cddi1      if(nmn.ne.0)then                                                  8d12s22
cddi1       ngroupsize=isize/2                                                1d29s21
cddi1      else                                                              8d12s22
cddi1       ngroupsize=isize                                                 8d12s22
cddi1      end if                                                            8d12s22
cddi1      if(immem.eq.0)then                                                1d29s21
cddi1       mymemp=irank+ncore                                               1d29s21
cddi1       mymast=irank                                                     1d29s21
cddi1       if(irank.ne.0)then
cddi1        write(fname,1066)irank
cddi1        do i=1,8
cddi1         if(fname(i:i).eq.' ')fname(i:i)=file(i)
cddi1        end do
cddi1        fname(1:1)='p'
cddi1        open(unit=6,file=fname)
cddi1       end if
cddi1c
cddi1c     compute node
cddi1c
cddi1       call mpi_comm_group(MPI_COMM_WORLD,igroup_comm_world,ierror)
cddi1       call mpi_group_incl(igroup_comm_world,ngroupsize,iranks,         1d29s21
cddi1     $      igroup_comm_group,ierror)
cddi1       do i=1,ngroupsize                                                1d29s21
cddi1        iranks(i)=iranks(i)+ncore                                       1d29s21
cddi1       end do                                                           1d29s21
cddi1       itag=3
cddi1       call mpi_comm_create_group(MPI_COMM_WORLD,igroup_comm_group,
cddi1     $      itag,my_comm_group,ierror)
cddi1       call mpi_comm_size(my_comm_group,iq,ierror)
cddi1       call mpi_comm_rank(my_comm_group,iq2,ierror)
cddi1       mynowprog=iq2                                                    11d18s20
cddi1       mynprocg=iq                                                      11d18s20
cddi1       return
cddi1      else
cddi1       ninfo(1)=0                                                       2d21s21
cddi1       ninfo(2)=0                                                       2d21s21
cddi1       nuse=1                                                           2d21s21
cddi1       write(fname,1066)irank
cddi1 1066  format(i8)
cddi1       do i=1,8
cddi1        if(fname(i:i).eq.' ')fname(i:i)=file(i)
cddi1       end do
cddi1       open(unit=6,file=fname)
cddi1       do i=1,ngroupsize                                                1d29s21
cddi1        iranks(i)=iranks(i)+ncore                                       1d29s21
cddi1       end do                                                           1d29s21
cddi1c
cddi1c     memory node
cddi1c
cddi1       ibcoff=1
cddi1       ndata=0
cddi1       ipass=0
cddi1    1  continue
cddi1        ipass=ipass+1
cddi1c
cddi1c     wait for incoming message
cddi1c
cddi1       itag=1
cddi1 3329  format(10000i1)
cddi1c
cddi1c     we need to use probe to determine the size of the message.
cddi1c     this is because if we do ddi_putx as two calls, one to
cddi1c     tell me the operation and size, and the second to transfer
cddi1c     the data, another processor might have sent us a message
cddi1c     between the two that will screw things up.
cddi1c
cddi1       call mpi_probe(mpi_any_source,itag,mpi_comm_world,istatus,
cddi1     $      ierror)
cddi1       call mpi_get_count(istatus,mpi_double_precision,nget,ierror)
cddi1       if(nget.le.4)then
c$$$cddi1        write(6,*)('mpi_recva '),nget
cddi1        call mpi_recv(packet,nget,mpi_double_precision,mpi_any_source,
cddi1     $      itag,mpi_comm_world,istatus,ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('mpi_recv returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop 'recva'
cddi1      end if
cddi1       else
cddi1        ipck=ibcoff
cddi1        ibcoff=ipck+nget
cddi1        call enough('ddimem.1',bc,ibc)
c$$$cddi1        write(6,*)('mpi_recvb '),nget
cddi1        call mpi_recv(bc(ipck),nget,mpi_double_precision,mpi_any_source,
cddi1     $      itag,mpi_comm_world,istatus,ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('mpi_recv returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop 'recvb'
cddi1      end if
cddi1        jpck=ipck-1
cddi1        do i=1,4
cddi1         packet(i)=bc(jpck+i)
cddi1        end do
cddi1        jpck=ipck+4
cddi1       end if
cddi1       icode=npacket(1)
cddi1       ninfo(nuse)=ninfo(nuse)+1                                        2d21s21
cddi1       if(ninfo(nuse).gt.idm)then                                       2d21s21
cddi1        ninfo(nuse)=ninfo(nuse)-1                                       2d21s21
cddi1        if(nuse.eq.1)then                                               2d21s21
cddi1         nuse=2                                                         2d21s21
cddi1        else                                                            2d21s21
cddi1         nuse=1                                                         2d21s21
cddi1        end if                                                          2d21s21
cddi1        ninfo(nuse)=1                                                   2d21s21
cddi1       end if                                                           2d21s21
cddi1       ifrom=npacket(6)                                                 2d21s21
cddi1       minfo(1,ninfo(nuse),nuse)=ifrom                                  2d21s21
cddi1       minfo(2,ninfo(nuse),nuse)=icode                                  2d21s21
cddi1c
cddi1c     branch on type of request
cddi1c
cddi1       if(icode.eq.1)then
cddi1c
cddi1c     creating a distributed array
cddi1c
cddi1        ndata=ndata+1
cddi1        minfo(3,ninfo(nuse),nuse)=ndata                                 2d21s21
cddi1        if(ndata.gt.id)then
cddi1         write(6,*)('id exceeded in create call '),ndata,id
cddi1         stop
cddi1        end if
cddi1        do i=1,4
cddi1         ip=i+1
cddi1         mdata(i,ndata)=npacket(ip)
cddi1         minfo(3+i,ninfo(nuse),nuse)=npacket(ip)                        2d21s21
cddi1        end do
cddi1        if(ibcoff.ne.mdata(3,ndata))then
cddi1         write(6,*)('compute and memory node out of sync '),ndata,
cddi1     $        ibcoff,mdata(3,ndata)
cddi1         do i=1,ndata
cddi1          write(6,*)i,(mdata(j,i),j=1,4)
cddi1         end do
cddi1         stop
cddi1        end if
cddi1        ibcoff=ibcoff+mdata(4,ndata)
cddi1        ipck=ibcoff
cddi1        call enough('ddimem.2',bc,ibc)
cddi1       else if(icode.eq.2)then
cddi1c
cddi1c     delete a distributed array
cddi1c
cddi1        nhandle=npacket(2)
cddi1        minfo(3,ninfo(nuse),nuse)=nhandle                               2d21s21
cddi1        if(nhandle.ne.ndata)then
cddi1         write(6,*)('destroying distributed arrays out of order ')
cddi1         write(6,*)('handle: '),nhandle
cddi1         write(6,*)('ndata: '),ndata
cddi1         stop
cddi1        end if
cddi1        ibcoff=mdata(3,ndata)
cddi1        ipck=ibcoff
cddi1        ndata=ndata-1
cddi1       else if(icode.eq.3)then
cddi1c
cddi1c     put part of array
cddi1c
cddi1        nhandle=npacket(2)
cddi1        nrow=npacket(3)
cddi1        istrt=npacket(4)
cddi1        iend=npacket(5)
cddi1        ilj=npacket(8)
cddi1        ifrom=npacket(6)
cddi1        ilow=npacket(7)-1
cddi1        minfo(3,ninfo(nuse),nuse)=nrow                                  2d21s21
cddi1        minfo(4,ninfo(nuse),nuse)=istrt                                 2d21s21
cddi1        minfo(5,ninfo(nuse),nuse)=iend                                  2d21s21
cddi1        minfo(6,ninfo(nuse),nuse)=ilj                                   2d21s21
cddi1        minfo(7,ninfo(nuse),nuse)=ilow                                  2d21s21
cddi1        xinfo(1,ninfo(nuse),nuse)=bc(ipck+4)                            2d21s21
cddi1        jtag=10
cddi1        if(nrow.eq.mdata(1,nhandle))then                                3d6s09
cddi1         iaddx=mdata(1,nhandle)*(istrt-ilj)+mdata(3,nhandle)
cddi1         nmove=nrow*(iend+1-istrt)                                      3d6s09
cddi1         jpck=ipck+4
cddi1         do j=0,nmove-1
cddi1          bc(iaddx+j)=bc(jpck+j)
cddi1         end do
cddi1        else
cddi1         itmp=ipck+4
cddi1         nwds=(iend+1-istrt)*nrow
cddi1         jtmp=itmp
cddi1         do i=istrt,iend
cddi1          iaddx=mdata(1,nhandle)*(i-ilj)+mdata(3,nhandle)+ilow
cddi1          do j=0,nrow-1
cddi1           bc(iaddx+j)=bc(jtmp+j)
cddi1          end do
cddi1          jtmp=jtmp+nrow
cddi1         end do
cddi1         ibcoff=ipck
cddi1        end if
cddi1       else if(icode.eq.4)then
cddi1c
cddi1c     accumulate into an array
cddi1c
cddi1        nhandle=npacket(2)
cddi1        nrow=npacket(3)
cddi1        istrt=npacket(4)
cddi1        iend=npacket(5)
cddi1        ifrom=npacket(6)
cddi1        ilow=npacket(7)-1
cddi1        ilj=npacket(8)
cddi1        minfo(3,ninfo(nuse),nuse)=nrow                                  2d21s21
cddi1        minfo(4,ninfo(nuse),nuse)=istrt                                 2d21s21
cddi1        minfo(5,ninfo(nuse),nuse)=iend                                  2d21s21
cddi1        minfo(6,ninfo(nuse),nuse)=ilj                                   2d21s21
cddi1        minfo(7,ninfo(nuse),nuse)=ilow                                  2d21s21
cddi1        xinfo(1,ninfo(nuse),nuse)=bc(ipck+4)                            2d21s21
cddi1        jtag=11
cddi1        if(nrow.eq.mdata(1,nhandle))then                                3d6s09
cddi1         iaddx=mdata(1,nhandle)*(istrt-ilj)+mdata(3,nhandle)
cddi1         nmove=nrow*(iend+1-istrt)                                      3d6s09
cddi1         itmp=ipck+4
cddi1         xinfo(2,ninfo(nuse),nuse)=bc(iaddx)                            2d21s21
cddi1         do i=0,nmove-1
cddi1          bc(iaddx+i)=bc(iaddx+i)+bc(itmp+i)
cddi1         end do
cddi1         xinfo(3,ninfo(nuse),nuse)=bc(iaddx)                            2d21s21
cddi1         ibcoff=ipck
cddi1        else
cddi1         nwds=(iend+1-istrt)*nrow
cddi1         itmp=ipck+4
cddi1         jtmp=itmp
cddi1         do i=istrt,iend
cddi1          iaddx=mdata(1,nhandle)*(i-ilj)+mdata(3,nhandle)+ilow
cddi1          if(i.eq.istrt)xinfo(2,ninfo(nuse),nuse)=bc(iaddx)                            2d21s21
cddi1          do j=0,nrow-1
cddi1           bc(iaddx+j)=bc(iaddx+j)+bc(jtmp+j)
cddi1          end do
cddi1          if(i.eq.istrt)xinfo(3,ninfo(nuse),nuse)=bc(iaddx)                            2d21s21
cddi1          jtmp=jtmp+nrow
cddi1         end do
cddi1         ibcoff=ipck
cddi1        end if
cddi1       else if(icode.eq.5)then
cddi1c
cddi1c     get part of array
cddi1c
cddi1        nhandle=npacket(2)
cddi1        nrow=npacket(3)
cddi1        istrt=npacket(4)
cddi1        iend=npacket(5)
cddi1        ilj=npacket(8)
cddi1        ifrom=npacket(6)
cddi1        ilow=npacket(7)-1
cddi1        jtag=12
cddi1        minfo(3,ninfo(nuse),nuse)=nrow                                  2d21s21
cddi1        minfo(4,ninfo(nuse),nuse)=istrt                                 2d21s21
cddi1        minfo(5,ninfo(nuse),nuse)=iend                                  2d21s21
cddi1        minfo(6,ninfo(nuse),nuse)=ilj                                   2d21s21
cddi1        minfo(7,ninfo(nuse),nuse)=ilow                                  2d21s21
cddi1        if(nrow.eq.mdata(1,nhandle))then                                3d6s09
cddi1         iaddx=mdata(1,nhandle)*(istrt-ilj)+mdata(3,nhandle)
cddi1         xinfo(1,ninfo(nuse),nuse)=bc(iaddx)
cddi1         nmove=nrow*(iend+1-istrt)                                      3d6s09
cddi1         call mpi_send(bc(iaddx),nmove,mpi_double_precision,ifrom,
cddi1     $        jtag,mpi_comm_world,ierror)
cddi1        else
cddi1         itmp=ibcoff
cddi1         nwds=(iend+1-istrt)*nrow
cddi1         ibcoff=itmp+nwds
cddi1         call enough('ddimem.3',bc,ibc)
cddi1         jtmp=itmp
cddi1         do i=istrt,iend
cddi1          iaddx=mdata(1,nhandle)*(i-ilj)+mdata(3,nhandle)+ilow
cddi1          if(i.eq.istrt)xinfo(1,ninfo(nuse),nuse)=bc(iaddx)             2d21s21
cddi1          do j=0,nrow-1
cddi1           bc(jtmp+j)=bc(iaddx+j)
cddi1          end do
cddi1          jtmp=jtmp+nrow
cddi1         end do
cddi1         call mpi_send(bc(itmp),nwds,mpi_double_precision,ifrom,
cddi1     $         jtag,mpi_comm_world,ierror)
cddi1         ibcoff=itmp
cddi1        end if
cddi1       else if(icode.eq.6)then
cddi1        ipck=ibcoff
cddi1c
cddi1c     initialize load balancing counter
cddi1c
cddi1        idlb=0
cddi1       else if(icode.eq.7)then
cddi1        ipck=ibcoff
cddi1c
cddi1c     give value of balancing counter
cddi1c
cddi1        ito=npacket(2)
cddi1        minfo(3,ninfo(nuse),nuse)=idlb                                  2d21s21
cddi1        itag=10
cddi1        call mpi_send(idlb,1,mpi_integer,ito,itag,mpi_comm_world,
cddi1     $       ierror)
cddi1        idlb=idlb+1
cddi1       else if(icode.eq.8)then
cddi1c
cddi1c     perform global sync
cddi1c
cddi1        call mpi_barrier(mpi_comm_world,ierror)
cddi1        ipck=ibcoff
cddi1       else if(icode.eq.9)then                                          11d21s20
cddi1c
cddi1c     zero out my part of array
cddi1c
cddi1        nhandle=npacket(2)                                              11d21s20
cddi1        minfo(3,ninfo(nuse),nuse)=nhandle                               2d21s21
cddi1        do i=0,mdata(4,nhandle)-1                                       11d21s20
cddi1         bc(mdata(3,nhandle)+i)=0d0                                     11d21s20
cddi1        end do                                                          11d21s20
cddi1       else if(icode.eq.10)then                                         2d19s21
cddi1c
cddi1c     accumulate into an array with message saying we are done.         2d19s21
cddi1c
cddi1        nhandle=npacket(2)
cddi1        nrow=npacket(3)
cddi1        istrt=npacket(4)
cddi1        iend=npacket(5)
cddi1        ifrom=npacket(6)
cddi1        ilow=npacket(7)-1
cddi1        ilj=npacket(8)
cddi1        minfo(3,ninfo(nuse),nuse)=nrow                                  2d21s21
cddi1        minfo(4,ninfo(nuse),nuse)=istrt                                 2d21s21
cddi1        minfo(5,ninfo(nuse),nuse)=iend                                  2d21s21
cddi1        minfo(6,ninfo(nuse),nuse)=ilj                                   2d21s21
cddi1        minfo(7,ninfo(nuse),nuse)=ilow                                  2d21s21
cddi1        xinfo(1,ninfo(nuse),nuse)=bc(ipck+4)                            2d21s21
cddi1        jtag=11                                                         8d11s22
cddi1        if(nrow.eq.mdata(1,nhandle))then                                3d6s09
cddi1         iaddx=mdata(1,nhandle)*(istrt-ilj)+mdata(3,nhandle)
cddi1         nmove=nrow*(iend+1-istrt)                                      3d6s09
cddi1         itmp=ipck+4
cddi1         nwiacc=nwiacc+nmove                                            8d10s22
cddi1         xinfo(2,ninfo(nuse),nuse)=bc(iaddx)                            2d21s21
cddi1         do i=0,nmove-1
cddi1          bc(iaddx+i)=bc(iaddx+i)+bc(itmp+i)
cddi1         end do
cddi1         xinfo(3,ninfo(nuse),nuse)=bc(iaddx)                            2d21s21
cddi1         ibcoff=ipck
cddi1        else
cddi1         nwds=(iend+1-istrt)*nrow
cddi1         nwiacc=nwiacc+nwds                                             8d10s22
cddi1         itmp=ipck+4
cddi1         jtmp=itmp
cddi1         do i=istrt,iend
cddi1          iaddx=mdata(1,nhandle)*(i-ilj)+mdata(3,nhandle)+ilow
cddi1          if(i.eq.istrt)xinfo(2,ninfo(nuse),nuse)=bc(iaddx)                            2d21s21
cddi1          do j=0,nrow-1
cddi1           bc(iaddx+j)=bc(iaddx+j)+bc(jtmp+j)
cddi1          end do
cddi1          if(i.eq.istrt)xinfo(3,ninfo(nuse),nuse)=bc(iaddx)                            2d21s21
cddi1          jtmp=jtmp+nrow
cddi1         end do
cddi1         ibcoff=ipck
cddi1        end if
cddi1       else if(icode.eq.11)then                                         8d10s22
cddi1c
cddi1c     zero iacc word counter
cddi1c
cddi1        nwiacc=0                                                        8d10s22
cddi1       else if(icode.eq.12)then                                         8d10s22
cddi1c
cddi1c     return iacc word counter
cddi1c
cddi1        ito=npacket(2)
cddi1        itag=10
cddi1        call mpi_send(nwiacc,1,mpi_integer,ito,itag,mpi_comm_world,
cddi1     $       ierror)
cddi1       else if(icode.eq.0)then
cddi1c
cddi1c     close down
cddi1c
cddi1        call mpi_barrier(mpi_comm_world,ierror)
cddi1        call mpi_finalize(ierror)
cddi1        stop
cddi1       end if
cddi1       ibcoff=ipck
cddi1       go to 1
cddi1      end if
cddi1      end
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine dws_finalize
cddi1      use mpi
cddi1      implicit integer (i-n)
cddi1      real*8 packet
cddi1      dimension npacket(2)
cddi1      equivalence (packet,npacket)
cddi1      include "common.mympi"                                            1d29s21
cddi1      if(mnmc.ne.0)then                                                 8d12s22
cddi1      itag=1
cddi1      call mpi_comm_rank(my_comm_group,irank,ierror)                     11d18s20
cddi1      ito=mymemp                                                        1d29s21
cddi1      npacket(1)=0
cddi1      call mpi_send(packet,1,mpi_double_precision,ito,itag,
cddi1     $     mpi_comm_world,
cddi1     $     ierror)
cddi1      end if                                                            8d12s22
cddi1      call mpi_barrier(mpi_comm_world,ierror)
cddi1      call mpi_finalize(ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in ddi_finalize,')
cddi1       write(6,*)('mpi_finalize returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop
cddi1      end if
cddi1      return
cddi1      end
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine dws_synca
cddi1      use mpi
cddi1      implicit integer (i-n)
cddi1      real*8 packet
cddi1      dimension npacket(2)
cddi1      equivalence (packet,npacket)
cddi1      include "common.mympi"                                            1d29s21
cddi1      if(mnmc.ne.0)then                                                 8d12s22
cddi1      itag=1
cddi1      call mpi_comm_rank(my_comm_group,irank,ierror)                    11d18s20
cddi1 3329 format(10000i1)
cddi1      ito=mymemp                                                        1d29s21
cddi1      call mpi_barrier(my_comm_group,ierror)
cddi1      npacket(1)=8
cddi1      call mpi_send(packet,1,mpi_double_precision,ito,itag,
cddi1     $     mpi_comm_world,
cddi1     $     ierror)
cddi1      end if                                                            8d12s22
cddi1      call mpi_barrier(mpi_comm_world,ierror)
cddi1      return
cddi1      end
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine dws_sync
cddi1      use mpi
cddi1      implicit integer (i-n)
cddi1      real*8 packet
cddi1      dimension npacket(2)
cddi1      equivalence (packet,npacket)
cddi1      include "common.mympi"                                            1d29s21
cddi1      call mpi_barrier(my_comm_group,ierror)
cddi1      return
cddi1      end
cddi1      subroutine second(time1)
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      data ifirst/0/
cddi1      save
cddi1      time1=mpi_wtime()                                                 5d7s12
cddi1      if(ifirst.eq.0)then
cddi1       time0=time1
cddi1       time1=0d0
cddi1       ifirst=1
cddi1      else
cddi1       time1=time1-time0
cddi1      end if
cddi1      return
cddi1      end
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine dws_gsumf(buff,nwds)
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer (i-n)
cddi1      dimension buff(nwds)                                                 2d22s19
cddi1      include "common.mympi"                                            1d29s21
cddi1      if(nwds.le.0)return                                               10d28s22
cddi1      call mpi_allreduce(mpi_in_place,buff,nwds,mpi_double_precision,   2d22s10
cddi1     $     mpi_sum,my_comm_group,ierror)                                2d22s10
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in dws_gsumf,')
cddi1       write(6,*)('mpi_allreduce returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop
cddi1      end if
cddi1      return
cddi1      end
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine dws_gbor(buff,nwds)                                    1d26s21
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer (i-n)
cddi1      integer buff(nwds)                                                1d26s21
cddi1      include "common.mympi"                                            1d29s21
cddi1      nwds2=nwds*2                                                      1d26s21
cddi1      call mpi_allreduce(mpi_in_place,buff,nwds2,mpi_integer,           1d26s21
cddi1     $     mpi_bor,my_comm_group,ierror)                                1d26s21
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in dws_gbor,')
cddi1       write(6,*)('mpi_allreduce returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop
cddi1      end if
cddi1      return
cddi1      end
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine dws_allgv2(buf,nblock8,ioff8,send)
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer (i-n)
cddi1      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
cddi1     $     mynnode
cddi1      include "common.mympi"                                            1d29s21
cddi1      dimension buf(1),nblock8(1),ioff8(1),send(1)
cccccccddi1      parameter (id=1000)
cddi1      dimension nblock(id),ioff(id)
cddi1      if(mynprocg.gt.id)then
cddi1       write(6,*)('too many procs in dws_allgv!! '),mynprocg,id
cddi1       close(unit=6)
cddi1      end if
cddi1      do i=1,mynprocg
cddi1       nblock(i)=nblock8(i)
cddi1       ioff(i)=ioff8(i)
cddi1      end do
cddi1      idum=nblock(mynowprog+1)
cddi1      call mpi_allgatherv(send,idum,mpi_double_precision,buf,
cddi1     $     nblock,ioff,mpi_double_precision,my_comm_group,ierr)
cddi1      if(ierr.ne.mpi_success)then
cddi1       write(6,*)('in dws_allgv, mpi_allgatherv returned an error'),
cddi1     $      ierror
cddi1       close(unit=6)
cddi1      end if
cddi1      return
cddi1      end
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine dws_all2allvb(bufs,nblock8s,ioff8s,bufr,nblock8r,
cddi1     $     ioff8r)
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer (i-n)
cddi1      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
cddi1     $     mynnode
cddi1      include "common.mympi"                                            1d29s21
cddi1      dimension bufs(1),nblock8s(1),ioff8s(1),bufr(1),nblock8r(1),
cddi1     $     ioff8r(1)
cccccccddi1      parameter (id=1000)
cddi1      dimension nblocks(id),ioffs(id),nblockr(id),ioffr(id)
cddi1      data icall/0/
cddi1      save
cddi1      icall=icall+1
cddi1      if(mynprocg.gt.id)then
cddi1       write(6,*)('too many procs in dws_allgv!! '),mynprocg,id
cddi1       close(unit=6)
cddi1      end if
cddi1      do i=1,mynprocg
cddi1       nblocks(i)=nblock8s(i)
cddi1       ioffs(i)=ioff8s(i)
cddi1       nblockr(i)=nblock8r(i)
cddi1       ioffr(i)=ioff8r(i)
cddi1 3030  format(5i8)
cddi1      end do
cddi1      call mpi_alltoallv(bufs,nblocks,ioffs,mpi_double_precision,
cddi1     $     bufr,nblockr,ioffr,mpi_double_precision,my_comm_group,ierr)
cddi1      if(ierr.ne.mpi_success)then
cddi1       write(6,*)('in dws_all2allv, mpi_alltoallv returned an error'),
cddi1     $      ierror
cddi1       close(unit=6)
cddi1      end if
cddi1      return
cddi1      end
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine dws_all2allvb8(bufs,nblock8s,ioff8s,bufr,nblock8r,
cddi1     $     ioff8r)
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer (i-n)
cddi1      integer*8 nblock8s,ioff8s,nblock8r,ioff8r                         1d10s18
cddi1      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
cddi1     $     mynnode
cddi1      dimension bufs(1),nblock8s(1),ioff8s(1),bufr(1),nblock8r(1),
cddi1     $     ioff8r(1)
cddi1      include "common.mympi"                                            1d29s21
ccccccddi1      parameter (id=1000)
cddi1      dimension nblocks(id),ioffs(id),nblockr(id),ioffr(id)
cddi1      data icall/0/
cddi1      save
cddi1      icall=icall+1
cddi1      if(mynprocg.gt.id)then
cddi1       write(6,*)('too many procs in dws_allgv!! '),mynprocg,id
cddi1       close(unit=6)
cddi1      end if
cddi1      do i=1,mynprocg
cddi1       nblocks(i)=nblock8s(i)
cddi1       ioffs(i)=ioff8s(i)
cddi1       nblockr(i)=nblock8r(i)
cddi1       ioffr(i)=ioff8r(i)
cddi1 3030  format(5i8)
cddi1      end do
cddi1      call mpi_alltoallv(bufs,nblocks,ioffs,mpi_double_precision,
cddi1     $     bufr,nblockr,ioffr,mpi_double_precision,my_comm_group,ierr)
cddi1      if(ierr.ne.mpi_success)then
cddi1       write(6,*)('in dws_all2allv, mpi_alltoallv returned an error'),
cddi1     $      ierror
cddi1       close(unit=6)
cddi1      end if
cddi1      return
cddi1      end
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine dws_12all(buf,len,isource8)
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer (i-n)
cddi1      dimension buf(len)
cddi1      include "common.mympi"                                            1d29s21
cddi1      data icall/0/
cddi1      save icall
cddi1      icall=icall+1
cddi1      icont=len
cddi1      isource=isource8
cddi1      call mpi_bcast(buf,icont,mpi_double_precision,isource,            2d19s10
cddi1     $               my_comm_group,ierror)                              11d18s20
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in dws_12all,')
cddi1       write(6,*)('mpi_bcast returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1      close(unit=6)
cddi1      end if
cddi1      return
cddi1      end
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine ddi_destroy(nhandlex)                                  11d15s22
cddi1      integer*8 nhandle,nhandlex
cddi1      nhandle=nhandlex+1
cddi1      call ddi_destroyx(nhandle)                                        11d15s22
cddi1      return
cddi1      end
cddi1      subroutine ddi_destroyx(nhandle)
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer (i-n)
cddi1      integer*8 irow,icol,nhandle(1),ilow1,ihigh1,ilow2,ihigh2,nhandlex,1d18s21
cddi1     $     ione,iarg3,npno,iarg1,iarg2,idelta(1)                        1d18s21
cddi1      logical no_locks
cddi1      integer (kind=mpi_address_kind) memreq,iadd,iaddx,itrial
cddi1      integer*4 iwin,nneed,ierror,info
c$$$cddi1      include "common.store"
cddi1      include "common.mympi"                                            1d29s21
cccccccccddi1      parameter (id=2000)
c$$$cddi1      dimension mdata(4,id),istatus(mpi_status_size),npacket(8),
c$$$cddi1     $     packet(4)
cddi1      dimension istatus(mpi_status_size),npacket(8),packet(4)
cddi1      equivalence(packet,npacket)
cddi1   51 format('ddi_destroy: ',i5)
cddi1      if(nhandle(1).ne.ndata)then                                       1d18s21
cddi1       write(6,*)('destroying distributed arrays out of order ')
cddi1       write(6,*)('handle: '),nhandle(1)                                1d18s21
cddi1       write(6,*)('ndata: '),ndata
cddi1       stop
cddi1      end if
cddi1      iwinbase=mdata(3,ndata)
cddi1      call mpi_comm_rank(my_comm_group,irank,ierror)
cddi1      ito=mymemp                                                        1d29s21
cddi1      itag=1
cddi1      npacket(1)=2
cddi1      npacket(2)=ndata
cddi1      call mpi_send(packet,1,mpi_double_precision,ito,itag,
cddi1     $     mpi_comm_world,ierror)
cddi1      ndata=ndata-1
cddi1      call dws_sync                                                     3d13s09
cddi1      return
cddi1      end
cddi1      subroutine ddi_zero(bc,ibc,nhandlex)                              7d18s24
cddi1      implicit real*8 (a-h,o-z)                                         7d18s24
cddi1      implicit integer*8 (i-n)                                          7d18s24
cddi1      integer my_comm_group,mymast,mymemp,iranks,mnmc,mdata,ndata,      1d26s23
cddi1     $ iwinbase                                                         1d26s23
cddi1      include "common.mympi"                                            1d29s21
cddi1      nhandle=nhandlex+1                                                11d19s20
cddi1      call ddi_zerox(bc,ibc,nhandle)                                    11d15s22
cddi1      return                                                            7d18s24
cddi1      end                                                               7d18s24
cddi1      subroutine ddi_put(bc,ibc,nhandlex,ilow1,ihigh1,ilow2,ihigh2,     7d18s24
cddi1     $     buff)                                                        7d18s24
cddi1      implicit real*8 (a-h,o-z)                                         7d18s24
cddi1      implicit integer*8 (i-n)                                          7d18s24
cddi1      integer my_comm_group,mymast,mymemp,iranks,mnmc,mdata,ndata,      1d26s23
cddi1     $ iwinbase                                                         1d26s23
cddi1      include "common.mympi"                                            1d29s21
cddi1      nhandle=nhandlex+1
cddi1      call ddi_putx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff)      11d15s22
cddi1      return
cddi1      end                                                               7d18s24
cddi1      subroutine ddi_iacc(bc,ibc,nhandlex,ilow1,ihigh1,ilow2,ihigh2,    7d18s24
cddi1     $   buff,iacc,nacc)                                                11d15s22
cddi1      implicit real*8 (a-h,o-z)                                         7d18s24
cddi1      implicit integer*8 (i-n)                                          7d18s24
cddi1      integer my_comm_group,mymast,mymemp,iranks,mnmc,mdata,ndata,      1d26s23
cddi1     $ iwinbase                                                         1d26s23
cddi1      include "common.mympi"                                            1d29s21
cddi1      nhandle=nhandlex+1                                                2d15s12
cddi1      call ddi_iaccx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff,iacc,11d15s22
cddi1     $   nacc)                                                          11d15s22
cddi1      return                                                            2d15s12
cddi1      end                                                               7d18s24
cddi1      subroutine ddi_acc(bc,ibc,nhandlex,ilow1,ihigh1,ilow2,ihigh2,buff)7d18s24
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer*8 (i-n)
cddi1      integer my_comm_group,mymast,mymemp,iranks,mnmc,mdata,ndata,      1d26s23                  6d17s21
cddi1     $ iwinbase                                                         1d26s23
cddi1      include "common.mympi"                                            1d29s21
cddi1      nhandle=nhandlex+1                                                2d15s12
cddi1      call ddi_accx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff)      11d15s22
cddi1      return                                                            2d15s12
cddi1      end                                                               7d18s24
cddi1      subroutine ddi_create(bc,ibc,irow,icol,nhandlex)                  11d15s22
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer*8 (i-n)
cddi1      integer my_comm_group,mymast,mymemp,iranks,mnmc,mdata,ndata,      1d26s23                  6d17s21
cddi1     $ iwinbase                                                         1d26s23
cddi1      include "common.mympi"                                            1d29s21
cddi1      call ddi_createx(bc,ibc,irow,icol,nhandle)                        11d15s22
cddi1      nhandlex=nhandle-1
cddi1      return
cddi1      end
cddi1      subroutine ddi_getx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff)7d18s24
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer (i-n)
cddi1      integer*8 irow,icol,nhandle(1),ilow1,ihigh1,ilow2,ihigh2,nhandlex,1d18s21
cddi1     $     ione,iarg3,npno,iarg1,iarg2,idelta(1)                        1d18s21
cddi1      logical no_locks
cddi1      integer (kind=mpi_address_kind) memreq,iadd,iaddx,itrial
cddi1      integer*4 iwin,nneed,ierror,info
cddi1      include "common.store"
cddi1      include "common.mympi"                                            1d29s21
cddi1      dimension istatus(mpi_status_size),npacket(8),packet(4)
cddi1      equivalence(packet,npacket)
cddi1      dimension buff(1),ircv(1)                                         1d29s21
cddi1      data loopx/50/                                                    6d17s21
cddi1      data ndata,ncall,loop/3*0/                                        6d17s21
cddi1      save
cddi1      call mpi_comm_rank(my_comm_group,irank,ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in ddi_create,')
cddi1       write(6,*)('mpi_comm_rank returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop
cddi1      end if
cddi1      call mpi_comm_size(my_comm_group,isize,ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in ddi_create,')
cddi1       write(6,*)('mpi_size_rank returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop
cddi1      end if
cddi1      if(nhandle(1).gt.ndata.or.nhandle(1).lt.1)then                    1d18s21
cddi1       write(6,*)('in ddi_get, handle = '),nhandle(1),nhandlex          1d18s21
cddi1       write(6,*)('exceeds ndata = '),ndata
cddi1       stop
cddi1      end if
cddi1      npacket(2)=nhandle(1)                                             1d18s21
cddi1      if(min(ilow1,ilow2).le.0)then
cddi1       write(6,*)('in ddi_get, bad lower limits '),ilow1,ilow2
cddi1       stop
cddi1      end if
cddi1      if(ihigh1.gt.mdata(1,nhandle(1)).or.                              1d18s21
cddi1     $     ihigh2.gt.mdata(2,nhandle(1)))then                           1d18s21
cddi1       write(6,*)('in ddi_get, bad upper limits '),ihigh1,ihigh2
cddi1       write(6,*)('vs '),mdata(1,nhandle(1)),mdata(2,nhandle(1))        1d18s21
cddi1       write(6,*)('nhandle = '),nhandle(1)                              1d18s21
cddi1       stop
cddi1      end if
cddi1      nrow=ihigh1+1-ilow1
cddi1      npacket(3)=nrow
cddi1      npacket(6)=mymast                                                 1d29s21
cddi1      npacket(7)=ilow1
cddi1      npacket(1)=5
cddi1      itag=1
cddi1      jtag=12
cddi1      do iproc=0,isize-1
cddi1       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
cddi1     $      ilj,ihj,i1s,i1e,i2s,i2e)
cddi1       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
cddi1        istrt=max(ilow2,ilj)
cddi1        iend=min(ihigh2,ihj)
cddi1        npacket(4)=istrt
cddi1        npacket(5)=iend
cddi1        npacket(8)=ilj
cddi1        ito=iranks(iproc+1)                                             1d29s21
cddi1        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
cddi1         i1=1+nrow*(istrt-ilow2)                                        3d6s09
cddi1         iaddx=mdata(1,nhandle(1))*(istrt-ilj)                          1d18s21
cddi1     $        +ilow1-1+mdata(3,nhandle(1))                              1d18s21
cddi1         nmove=nrow*(iend+1-istrt)                                      3d6s09
cddi1         if(nmove.gt.0)then
cddi1         call mpi_send(packet,4,mpi_double_precision,ito,itag,
cddi1     $         mpi_comm_world,ierror)
c$$$cddi1        write(6,*)('mpi_recvc '),nmove
cddi1         call mpi_recv(buff(i1),nmove,mpi_double_precision,ito,jtag,
cddi1     $        mpi_comm_world,istatus,ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('mpi_recv returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop 'recvc'
cddi1      end if
cddi1         end if
cddi1        else                                                            3d6s09
cddi1         nwds=(iend+1-istrt)*nrow
cddi1         if(nwds.gt.0)then
cddi1         itmp=ibcoff
cddi1         ibcoff=itmp+nwds
cddi1         call enough('ddi_createx.  5',bc,ibc)
cddi1         call mpi_send(packet,4,mpi_double_precision,ito,itag,
cddi1     $         mpi_comm_world,
cddi1     $       ierror)
c$$$cddi1        write(6,*)('mpi_recvd '),nwds
cddi1         call mpi_recv(bc(itmp),nwds,mpi_double_precision,ito,jtag,
cddi1     $        mpi_comm_world,istatus,ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('mpi_recv returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop 'recvd'
cddi1      end if
cddi1         jtmp=itmp
cddi1         do i=istrt,iend
cddi1          i1=1+nrow*(i-ilow2)
cddi1          do j=0,nrow-1
cddi1           buff(i1+j)=bc(jtmp+j)
cddi1          end do
cddi1          jtmp=jtmp+nrow
cddi1         end do
cddi1         ibcoff=itmp
cddi1         end if
cddi1       end if
cddi1       end if
cddi1      end do
cddi1      return
cddi1      end                                                               7d18s24
cddi1      subroutine ddi_iaccx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,    7d18s24
cddi1     $     buff,ircv,nrcv)                                              7d18s24
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer (i-n)
cddi1      integer*8 irow,icol,nhandle(1),ilow1,ihigh1,ilow2,ihigh2,nhandlex,1d18s21
cddi1     $     ione,iarg3,npno,iarg1,iarg2,idelta(1)                        1d18s21
cddi1      logical no_locks
cddi1      integer (kind=mpi_address_kind) memreq,iadd,iaddx,itrial
cddi1      integer*4 iwin,nneed,ierror,info
cddi1      include "common.store"
cddi1      include "common.mympi"                                            1d29s21
cddi1      dimension istatus(mpi_status_size),npacket(8),packet(4)
cddi1      equivalence(packet,npacket)
cddi1      dimension buff(1),ircv(1)                                         1d29s21
cddi1      data loopx/50/                                                    6d17s21
cddi1      data ndata,ncall,loop/3*0/                                        6d17s21
cddi1      save
cddi1c
cddi1c     important note!!!
cddi1c     buff needs to have 4*mynprocg words of extra storage for this to
cddi1c     work properly, yet there is no way to test for this !!!
cddi1c
cddi1      call mpi_comm_rank(my_comm_group,irank,ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in ddi_create,')
cddi1       write(6,*)('mpi_comm_rank returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop
cddi1      end if
cddi1      call mpi_comm_size(my_comm_group,isize,ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in ddi_create,')
cddi1       write(6,*)('mpi_size_rank returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop
cddi1      end if
cddi1      if(nhand.gt.ndata.or.nhandle(1).lt.1)then
cddi1       write(6,*)('in ddi_accx, handle = '),nhandle(1),nhandlex,nhand
cddi1       write(6,*)('exceeds ndata = '),ndata
cddi1       stop
cddi1      end if
cddi1      npacket(2)=nhandle(1)                                             1d18s21
cddi1      if(min(ilow1,ilow2).le.0)then
cddi1       write(6,*)('in ddi_put, bad lower limits '),ilow1,ilow2
cddi1       stop
cddi1      end if
cddi1      if(ihigh1.gt.mdata(1,nhandle(1)).or.                              1d18s21
cddi1     $     ihigh2.gt.mdata(2,nhandle(1)))then                           1d18s21
cddi1       write(6,*)('in ddi_accx, bad upper limits '),ihigh1,ihigh2,
cddi1     $      mdata(1,nhandle(1)),mdata(2,nhandle(1)),nhandle(1)          1d18s21
cddi1       stop
cddi1      end if
cddi1      nrow=ihigh1+1-ilow1
cddi1      npacket(3)=nrow
cddi1      npacket(6)=mymast                                                 1d29s21
cddi1      npacket(7)=ilow1
cddi1      npacket(1)=10
cddi1      itag=1
cddi1      jtag=11
cddi1      nrcv=0                                                            1d30s21
cddi1      itmp0=ibcoff                                                      1d30s21
cddi1      jtmp0=itmp0                                                       1d30s21
cddi1      do iproc=0,isize-1
cddi1       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
cddi1     $      ilj,ihj,i1s,i1e,i2s,i2e)
cddi1       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
cddi1        istrt=max(ilow2,ilj)
cddi1        iend=min(ihigh2,ihj)
cddi1        npacket(4)=istrt
cddi1        npacket(5)=iend
cddi1        npacket(8)=ilj
cddi1        ito=iranks(iproc+1)                                             1d29s21
cddi1        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
cddi1         i1=1+nrow*(istrt-ilow2)                                        3d6s09
cddi1         iaddx=mdata(1,nhandle(1))*(istrt-ilj)                          1d18s21
cddi1     $        +ilow1-1+mdata(3,nhandle(1))                              1d18s21
cddi1         nmove=nrow*(iend+1-istrt)                                      3d6s09
cddi1         if(nmove.gt.0)then
cddi1          itmp=jtmp0                                                    1d30s21
cddi1          jtmp0=itmp+nmove+4                                            1d30s21
cddi1          ibcoff=jtmp0                                                   1d30s21
cddi1          call enough('ddi_createx.  3',bc,ibc)
cddi1          jtmp=itmp-1
cddi1          do j=1,4
cddi1           bc(jtmp+j)=packet(j)
cddi1          end do
cddi1          jtmp=jtmp+5
cddi1          do j=0,nmove-1
cddi1           bc(jtmp+j)=buff(i1+j)
cddi1          end do
cddi1         end if
cddi1        else                                                            3d6s09
cddi1         nwds=(iend+1-istrt)*nrow
cddi1         if(nwds.gt.0)then
cddi1         itmp=jtmp0                                                     1d30s21
cddi1         jtmp0=itmp+nwds+4                                              1d30s21
cddi1         ibcoff=jtmp0                                                   1d30s21
cddi1         call enough('ddi_createx.  4',bc,ibc)
cddi1         jtmp=itmp-1
cddi1         do j=1,4
cddi1          bc(jtmp+j)=packet(j)
cddi1         end do
cddi1         jtmp=jtmp+5
cddi1         do i=istrt,iend
cddi1          i1=1+nrow*(i-ilow2)
cddi1          do j=0,nrow-1
cddi1           bc(jtmp+j)=buff(i1+j)
cddi1          end do
cddi1          jtmp=jtmp+nrow
cddi1         end do
cddi1         nwdsp=nwds+4
cddi1         end if
cddi1        end if                                                          3d6s09
cddi1       end if
cddi1      end do
cddi1      ntotal=jtmp0-itmp0                                                1d30s21
cddi1      jtmp0=itmp0-1                                                     1d30s21
cddi1      do i=1,ntotal                                                     1d30s21
cddi1       buff(i)=bc(jtmp0+i)                                              1d30s21
cddi1      end do                                                            1d30s21
cddi1      ibcoff=itmp0                                                      1d30s21
cddi1      i1=1
cddi1      myall=0
cddi1      do iproc=0,isize-1
cddi1       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
cddi1     $      ilj,ihj,i1s,i1e,i2s,i2e)
cddi1       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
cddi1        istrt=max(ilow2,ilj)
cddi1        iend=min(ihigh2,ihj)
cddi1        npacket(4)=istrt
cddi1        npacket(5)=iend
cddi1        npacket(8)=ilj
cddi1        ito=iranks(iproc+1)                                             1d29s21
cddi1        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
cddi1         nmove=nrow*(iend+1-istrt)                                      3d6s09
cddi1         if(nmove.gt.0)then
cddi1          nmovep=nmove+4
cddi1          nrcv=nrcv+1                                                   1d30s21
cddi1          myall=myall+nmove
cddi1          call mpi_isend(buff(i1),nmovep,mpi_double_precision,ito,itag, 6d17s21
cddi1     $        mpi_comm_world,ircv(nrcv),ierror)                         1d30s21
cddi1          i1=i1+nmovep                                                  1d30s21
cddi1         end if
cddi1        else                                                            3d6s09
cddi1         nwds=(iend+1-istrt)*nrow
cddi1         if(nwds.gt.0)then
cddi1          nwdsp=nwds+4
cddi1          nrcv=nrcv+1                                                   1d30s21
cddi1          myall=myall+nmove
cddi1          call mpi_isend(buff(i1),nwdsp,mpi_double_precision,ito,itag,  6d17s21
cddi1     $        mpi_comm_world,ircv(nrcv),ierror)                         1d30s21
cddi1          i1=i1+nwdsp                                                   1d30s21
cddi1         end if
cddi1        end if                                                          3d6s09
cddi1       end if
cddi1      end do
cddi1      return
cddi1      end                                                               7d18s24
cddi1      subroutine ddi_accx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff)7d18s24
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer (i-n)
cddi1      integer*8 irow,icol,nhandle(1),ilow1,ihigh1,ilow2,ihigh2,nhandlex,1d18s21
cddi1     $     ione,iarg3,npno,iarg1,iarg2,idelta(1)                        1d18s21
cddi1      logical no_locks
cddi1      integer (kind=mpi_address_kind) memreq,iadd,iaddx,itrial
cddi1      integer*4 iwin,nneed,ierror,info
cddi1      include "common.store"
cddi1      include "common.mympi"                                            1d29s21
cddi1      dimension istatus(mpi_status_size),npacket(8),packet(4)
cddi1      equivalence(packet,npacket)
cddi1      dimension buff(1),ircv(1)                                         1d29s21
cddi1      data loopx/50/                                                    6d17s21
cddi1      data ndata,ncall,loop/3*0/                                        6d17s21
cddi1      save
cddi1      call mpi_comm_rank(my_comm_group,irank,ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in ddi_create,')
cddi1       write(6,*)('mpi_comm_rank returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop
cddi1      end if
cddi1      call mpi_comm_size(my_comm_group,isize,ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in ddi_create,')
cddi1       write(6,*)('mpi_size_rank returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop
cddi1      end if
cddi1      if(nhand.gt.ndata.or.nhandle(1).lt.1)then
cddi1       write(6,*)('in ddi_accx, handle = '),nhandle(1),nhandlex,nhand
cddi1       write(6,*)('exceeds ndata = '),ndata
cddi1       stop
cddi1      end if
cddi1      npacket(2)=nhandle(1)                                             1d18s21
cddi1      if(min(ilow1,ilow2).le.0)then
cddi1       write(6,*)('in ddi_put, bad lower limits '),ilow1,ilow2
cddi1       stop
cddi1      end if
cddi1      if(ihigh1.gt.mdata(1,nhandle(1)).or.                              1d18s21
cddi1     $     ihigh2.gt.mdata(2,nhandle(1)))then                           1d18s21
cddi1       write(6,*)('in ddi_accx, bad upper limits '),ihigh1,ihigh2,
cddi1     $      mdata(1,nhandle(1)),mdata(2,nhandle(1)),nhandle(1)          1d18s21
cddi1       stop
cddi1      end if
cddi1      nrow=ihigh1+1-ilow1
cddi1      npacket(3)=nrow
cddi1      npacket(6)=mymast                                                 1d29s21
cddi1      npacket(7)=ilow1
cddi1      npacket(1)=4                                                      2d21s21
cddi1      itag=1
cddi1      jtag=11
cddi1      ktag=13                                                           2d19s21
cddi1      do iproc=0,isize-1
cddi1       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
cddi1     $      ilj,ihj,i1s,i1e,i2s,i2e)
cddi1       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
cddi1        istrt=max(ilow2,ilj)
cddi1        iend=min(ihigh2,ihj)
cddi1        npacket(4)=istrt
cddi1        npacket(5)=iend
cddi1        npacket(8)=ilj
cddi1        ito=iranks(iproc+1)                                             1d29s21
cddi1        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
cddi1         i1=1+nrow*(istrt-ilow2)                                        3d6s09
cddi1         iaddx=mdata(1,nhandle(1))*(istrt-ilj)                          1d18s21
cddi1     $        +ilow1-1+mdata(3,nhandle(1))                              1d18s21
cddi1         nmove=nrow*(iend+1-istrt)                                      3d6s09
cddi1         if(nmove.gt.0)then
cddi1          itmp=ibcoff
cddi1          ibcoff=itmp+nmove+4
cddi1          jtmp=itmp-1
cddi1          do j=1,4
cddi1           bc(jtmp+j)=packet(j)
cddi1          end do
cddi1          jtmp=jtmp+5
cddi1          do j=0,nmove-1
cddi1           bc(jtmp+j)=buff(i1+j)
cddi1          end do
cddi1          nmovep=nmove+4
cddi1         call mpi_ssend(bc(itmp),nmovep,mpi_double_precision,ito,itag,  2d5s21
cddi1     $        mpi_comm_world,ierror)
cddi1          ibcoff=itmp
cddi1         end if
cddi1        else                                                            3d6s09
cddi1         nwds=(iend+1-istrt)*nrow
cddi1         if(nwds.gt.0)then
cddi1         itmp=ibcoff
cddi1         ibcoff=itmp+nwds+4
cddi1         call enough('ddi_createx.  2',bc,ibc)
cddi1         jtmp=itmp-1
cddi1         do j=1,4
cddi1          bc(jtmp+j)=packet(j)
cddi1         end do
cddi1         jtmp=jtmp+5
cddi1         do i=istrt,iend
cddi1          i1=1+nrow*(i-ilow2)
cddi1          do j=0,nrow-1
cddi1           bc(jtmp+j)=buff(i1+j)
cddi1          end do
cddi1          jtmp=jtmp+nrow
cddi1         end do
cddi1         nwdsp=nwds+4
cddi1         call mpi_ssend(bc(itmp),nwdsp,mpi_double_precision,ito,itag,   2d5s21
cddi1     $        mpi_comm_world,ierror)
cddi1         ibcoff=itmp
cddi1         end if
cddi1        end if                                                          3d6s09
cddi1       end if
cddi1      end do
cddi1      return
cddi1      end                                                               7d18s24
cddi1      subroutine ddi_putx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff)7d18s24
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer (i-n)
cddi1      integer*8 irow,icol,nhandle(1),ilow1,ihigh1,ilow2,ihigh2,nhandlex,1d18s21
cddi1     $     ione,iarg3,npno,iarg1,iarg2,idelta(1)                        1d18s21
cddi1      logical no_locks
cddi1      integer (kind=mpi_address_kind) memreq,iadd,iaddx,itrial
cddi1      integer*4 iwin,nneed,ierror,info
cddi1      include "common.store"
cddi1      include "common.mympi"                                            1d29s21
cddi1      dimension istatus(mpi_status_size),npacket(8),packet(4)
cddi1      equivalence(packet,npacket)
cddi1      dimension buff(1),ircv(1)                                         1d29s21
cddi1      data loopx/50/                                                    6d17s21
cddi1      data ndata,ncall,loop/3*0/                                        6d17s21
cddi1      save
cddi1      call mpi_comm_rank(my_comm_group,irank,ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in ddi_create,')
cddi1       write(6,*)('mpi_comm_rank returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop
cddi1      end if
cddi1      call mpi_comm_size(my_comm_group,isize,ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in ddi_create,')
cddi1       write(6,*)('mpi_size_rank returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop
cddi1      end if
cddi1      if(nhand.gt.ndata.or.nhandle(1).lt.1)then                         1d18s21
cddi1       write(6,*)('in ddi_put, handle = '),nhandle(1),nhandlex,nhand    1d18s21
cddi1       write(6,*)('exceeds ndata = '),ndata
cddi1       stop
cddi1      end if
cddi1      npacket(1)=3
cddi1      npacket(2)=nhandle(1)                                             1d18s21
cddi1      if(min(ilow1,ilow2).le.0)then
cddi1       write(6,*)('in ddi_put, bad lower limits '),ilow1,ilow2
cddi1       stop
cddi1      end if
cddi1      if(ihigh1.gt.mdata(1,nhandle(1)).or                               1d18s21
cddi1     $     .ihigh2.gt.mdata(2,nhandle(1)))then                          1d18s21
cddi1       write(6,*)('in ddi_put, bad upper limits '),ihigh1,ihigh2,
cddi1     $      mdata(1,nhandle(1)),mdata(2,nhandle(1)),nhandle(1)          1d18s21
cddi1       stop
cddi1      end if
cddi1      nrow=ihigh1+1-ilow1
cddi1      npacket(3)=nrow
cddi1      npacket(5)=mymast                                                 1d29s21
cddi1      npacket(7)=ilow1
cddi1      npacket(8)=0
cddi1      itag=1
cddi1      jtag=10
cddi1      do iproc=0,isize-1
cddi1       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
cddi1     $      ilj,ihj,i1s,i1e,i2s,i2e)
cddi1       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
cddi1        istrt=max(ilow2,ilj)
cddi1        iend=min(ihigh2,ihj)
cddi1        npacket(4)=istrt
cddi1        npacket(5)=iend
cddi1        npacket(8)=ilj
cddi1        ito=iranks(iproc+1)                                             1d29s21
cddi1        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
cddi1         i1=1+nrow*(istrt-ilow2)                                        3d6s09
cddi1         iaddx=mdata(1,nhandle(1))*(istrt-ilj)+ilow1-1                  1d18s21
cddi1     $        +mdata(3,nhandle(1))                                      1d18s21
cddi1         nmove=nrow*(iend+1-istrt)                                      3d6s09
cddi138835    format(4i8,1pe15.7)
cddi1         itmp=ibcoff
cddi1         ibcoff=itmp+4+nmove
cddi1         jtmp=itmp-1
cddi1         do j=1,4
cddi1          bc(jtmp+j)=packet(j)
cddi1         end do
cddi1         jtmp=jtmp+5
cddi1         do j=0,nmove-1
cddi1          bc(jtmp+j)=buff(i1+j)
cddi1         end do
cddi1         nmovep=nmove+4
cddi1         call mpi_send(bc(itmp),nmovep,mpi_double_precision,ito,itag,
cddi1     $        mpi_comm_world,ierror)
cddi1         ibcoff=itmp
cddi1        else                                                            3d6s09
cddi1         nwds=(iend+1-istrt)*nrow
cddi1         if(nwds.gt.0)then
cddi1         itmp=ibcoff
cddi1         ibcoff=itmp+nwds+4
cddi1         call enough('ddi_createx.  1',bc,ibc)
cddi1         jtmp=itmp-1
cddi1         do j=1,4
cddi1          bc(jtmp+j)=packet(j)
cddi1         end do
cddi1         jtmp=jtmp+5
cddi1         do i=istrt,iend
cddi1          i1=1+nrow*(i-ilow2)
cddi1          do j=0,nrow-1
cddi1           bc(jtmp+j)=buff(i1+j)
cddi1          end do
cddi1          jtmp=jtmp+nrow
cddi1c          iaddx=mdata(1,nhandle)*(i-ilj)+ilow1-1+mdata(3,nhandle)        3d6s09
cddi1   33    format(5i8)
cddi1         end do
cddi1         nwdsp=nwds+4
cddi1         call mpi_send(bc(itmp),nwdsp,mpi_double_precision,ito,
cddi1     $        itag,mpi_comm_world,ierror)
cddi1         ibcoff=itmp
cddi1         end if
cddi1        end if
cddi1       end if                                                           3d6s09
cddi1      end do
cddi1      return
cddi1      end
cddi1      subroutine ddi_igetx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,    7d18s24
cddi1     $   buff,ircv,nrcv)                                                11d15s22
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer (i-n)
cddi1      integer*8 irow,icol,nhandle(1),ilow1,ihigh1,ilow2,ihigh2,nhandlex,1d18s21
cddi1     $     ione,iarg3,npno,iarg1,iarg2,idelta(1)                        1d18s21
cddi1      logical no_locks
cddi1      integer (kind=mpi_address_kind) memreq,iadd,iaddx,itrial
cddi1      integer*4 iwin,nneed,ierror,info
cddi1      include "common.store"
cddi1      include "common.mympi"                                            1d29s21
cddi1      dimension istatus(mpi_status_size),npacket(8),packet(4)
cddi1      equivalence(packet,npacket)
cddi1      dimension buff(1),ircv(1)                                         1d29s21
cddi1      data loopx/50/                                                    6d17s21
cddi1      data ndata,ncall,loop/3*0/                                        6d17s21
cddi1      save
cddi1      call mpi_comm_rank(my_comm_group,irank,ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in ddi_create,')
cddi1       write(6,*)('mpi_comm_rank returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop
cddi1      end if
cddi1      call mpi_comm_size(my_comm_group,isize,ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in ddi_create,')
cddi1       write(6,*)('mpi_size_rank returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop
cddi1      end if
cddi1      if(nhandle(1).gt.ndata.or.nhandle(1).lt.1)then                    1d18s21
cddi1       write(6,*)('in ddi_igetx, handle = '),nhandle(1),nhandlex          1d18s21
cddi1       write(6,*)('exceeds ndata = '),ndata
cddi1       stop
cddi1      end if
cddi1      npacket(2)=nhandle(1)                                             1d18s21
cddi1      if(min(ilow1,ilow2).le.0)then
cddi1       write(6,*)('in ddi_igetx, bad lower limits '),ilow1,ilow2
cddi1       stop
cddi1      end if
cddi1      if(ihigh1.gt.mdata(1,nhandle(1)).or.                              1d18s21
cddi1     $     ihigh2.gt.mdata(2,nhandle(1)))then                           1d18s21
cddi1       write(6,*)('in ddi_iget, bad upper limits '),ihigh1,ihigh2
cddi1       write(6,*)('vs '),mdata(1,nhandle(1)),mdata(2,nhandle(1))        1d18s21
cddi1       write(6,*)('nhandle = '),nhandle(1)                              1d18s21
cddi1       stop
cddi1      end if
cddi1      nrow=ihigh1+1-ilow1
cddi1      npacket(3)=nrow
cddi1      npacket(6)=mymast                                                 1d29s21
cddi1      npacket(7)=ilow1
cddi1      npacket(1)=5
cddi1      itag=1
cddi1      jtag=12
cddi1      nrcv=0                                                            1d29s21
cddi1      do iproc=0,isize-1
cddi1       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
cddi1     $      ilj,ihj,i1s,i1e,i2s,i2e)
cddi1       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
cddi1        istrt=max(ilow2,ilj)
cddi1        iend=min(ihigh2,ihj)
cddi1        npacket(4)=istrt
cddi1        npacket(5)=iend
cddi1        npacket(8)=ilj
cddi1        ito=iranks(iproc+1)                                             1d29s21
cddi1        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
cddi1         i1=1+nrow*(istrt-ilow2)                                        3d6s09
cddi1         iaddx=mdata(1,nhandle(1))*(istrt-ilj)                          1d18s21
cddi1     $        +ilow1-1+mdata(3,nhandle(1))                              1d18s21
cddi1         nmove=nrow*(iend+1-istrt)                                      3d6s09
cddi1         if(nmove.gt.0)then
cddi1         call mpi_send(packet,4,mpi_double_precision,ito,itag,
cddi1     $         mpi_comm_world,ierror)
cddi1         nrcv=nrcv+1                                                    1d29s21
cddi1         call mpi_irecv(buff(i1),nmove,mpi_double_precision,ito,jtag,   1d29s21
cddi1     $        mpi_comm_world,ircv(nrcv),ierror)
cddi1         end if
cddi1        else                                                            3d6s09
cddi1         nwds=(iend+1-istrt)*nrow
cddi1         if(nwds.gt.0)then
cddi1         itmp=ibcoff
cddi1         ibcoff=itmp+nwds
cddi1         call enough('ddi_createx.  6',bc,ibc)
cddi1         call mpi_send(packet,4,mpi_double_precision,ito,itag,
cddi1     $         mpi_comm_world,
cddi1     $       ierror)
c$$$cddi1        write(6,*)('mpi_recve '),nwds
cddi1         call mpi_recv(bc(itmp),nwds,mpi_double_precision,ito,jtag,     1d29s21
cddi1     $        mpi_comm_world,istatus,ierror)                            1d29s21
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('mpi_recv returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop 'recve'
cddi1      end if
cddi1         jtmp=itmp
cddi1         do i=istrt,iend
cddi1          i1=1+nrow*(i-ilow2)
cddi1          do j=0,nrow-1
cddi1           buff(i1+j)=bc(jtmp+j)
cddi1          end do
cddi1          jtmp=jtmp+nrow
cddi1         end do
cddi1         ibcoff=itmp
cddi1         end if
cddi1       end if
cddi1       end if
cddi1      end do
cddi1      return
cddi1      end                                                               7d18s24
cddi1      subroutine ddi_zerox(bc,ibc,nhandle)                              7d18s24
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer (i-n)
cddi1      integer*8 irow,icol,nhandle(1),ilow1,ihigh1,ilow2,ihigh2,nhandlex,1d18s21
cddi1     $     ione,iarg3,npno,iarg1,iarg2,idelta(1)                        1d18s21
cddi1      logical no_locks
cddi1      integer (kind=mpi_address_kind) memreq,iadd,iaddx,itrial
cddi1      integer*4 iwin,nneed,ierror,info
cddi1      include "common.store"
cddi1      include "common.mympi"                                            1d29s21
cddi1      dimension istatus(mpi_status_size),npacket(8),packet(4)
cddi1      equivalence(packet,npacket)
cddi1      dimension buff(1),ircv(1)                                         1d29s21
cddi1      data loopx/50/                                                    6d17s21
cddi1      data ndata,ncall,loop/3*0/                                        6d17s21
cddi1      save
cddi1      ito=mymemp                                                        1d29s21
cddi1      itag=1
cddi1      npacket(1)=9                                                      11d21s20
cddi1      npacket(2)=nhandle(1)                                             1d18s21
cddi1      call mpi_send(packet,2,mpi_double_precision,ito,itag,             11d21s20
cddi1     $     mpi_comm_world,ierror)
cddi1      return                                                            11d19s20
cddi1      end                                                               7d18s24
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine ddi_createx(bc,ibc,irow,icol,nhandle)                  11d15s22
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer (i-n)
cddi1      integer*8 irow,icol,nhandle(1),ilow1,ihigh1,ilow2,ihigh2,nhandlex,1d18s21
cddi1     $     ione,iarg3,npno,iarg1,iarg2,idelta(1)                        1d18s21
cddi1      logical no_locks
cddi1      integer (kind=mpi_address_kind) memreq,iadd,iaddx,itrial
cddi1      integer*4 iwin,nneed,ierror,info
cddi1      include "common.store"
cddi1      include "common.mympi"                                            1d29s21
cccccccccddi1      parameter (id=2000)
c$$$cddi1      dimension mdata(4,id),istatus(mpi_status_size),npacket(8),
c$$$cddi1     $     packet(4)
cddi1      dimension istatus(mpi_status_size),npacket(8),packet(4)
cddi1      equivalence(packet,npacket)
cddi1      dimension buff(1),ircv(1)                                         1d29s21
cddi1      data loopx/50/                                                    6d17s21
cddi1      data ndata,ncall,loop/3*0/                                        6d17s21
cddi1      save
cddi1      if(ndata.eq.0)then
cddi1       nrun=0
cddi1       iwinbase=1
cddi1      end if
cddi1      ndata=ndata+1
cddi1      if(ndata.gt.id)then
cddi1       write(6,*)('tried to create distributed array '),ndata
cddi1       write(6,*)('but in ddi_create, id = '),id
cddi1       stop
cddi1      end if
cddi1      nhandle(1)=ndata                                                  1d18s21
cddi1      mdata(1,ndata)=irow
cddi1      mdata(2,ndata)=icol
cddi1      call mpi_comm_rank(my_comm_group,irank,ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in ddi_create,')
cddi1       write(6,*)('mpi_comm_rank returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop
cddi1      end if
cddi1      call mpi_comm_size(my_comm_group,isize,ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('in ddi_create,')
cddi1       write(6,*)('mpi_size_rank returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop
cddi1      end if
cddi1      n1=1
cddi1      n2=icol
cddi1      call ilimts(n1,n2,isize,0,ilj,ihj,i1s,i1e,i2s,i2e)
cddi1      nneed=irow*(ihj+1-ilj)                                            3d3s09
cddi1      mdata(3,ndata)=iwinbase
cddi1      iwinbase=iwinbase+nneed                                           3d6s09
cddi1      mdata(4,ndata)=nneed
cddi1      ito=mymemp                                                        1d29s21
cddi1      itag=1
cddi1      npacket(1)=itag
cddi1      do i=1,4
cddi1       ip=i+1
cddi1       npacket(ip)=mdata(i,ndata)
cddi1      end do
cddi1      call mpi_send(packet,4,mpi_double_precision,ito,itag,
cddi1     $     mpi_comm_world,ierror)
cddi1 3356 format('ddi_create: ',i5,5x,3i5,i8)
cddi1      return
c$$$cddi1      entry ddi_zerox(bc,ibc,nhandle)                                   11d15s22
c$$$cddi1      ito=mymemp                                                        1d29s21
c$$$cddi1      itag=1
c$$$cddi1      npacket(1)=9                                                      11d21s20
c$$$cddi1      npacket(2)=nhandle(1)                                             1d18s21
c$$$cddi1      call mpi_send(packet,2,mpi_double_precision,ito,itag,             11d21s20
c$$$cddi1     $     mpi_comm_world,ierror)
c$$$cddi1      return                                                            11d19s20
c$$$cddi1      entry ddi_putx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff)     11d15s22
c$$$cddi1      if(nhand.gt.ndata.or.nhandle(1).lt.1)then                         1d18s21
c$$$cddi1       write(6,*)('in ddi_put, handle = '),nhandle(1),nhandlex,nhand    1d18s21
c$$$cddi1       write(6,*)('exceeds ndata = '),ndata
c$$$cddi1       stop
c$$$cddi1      end if
c$$$cddi1      npacket(1)=3
c$$$cddi1      npacket(2)=nhandle(1)                                             1d18s21
c$$$cddi1      if(min(ilow1,ilow2).le.0)then
c$$$cddi1       write(6,*)('in ddi_put, bad lower limits '),ilow1,ilow2
c$$$cddi1       stop
c$$$cddi1      end if
c$$$cddi1      if(ihigh1.gt.mdata(1,nhandle(1)).or                               1d18s21
c$$$cddi1     $     .ihigh2.gt.mdata(2,nhandle(1)))then                          1d18s21
c$$$cddi1       write(6,*)('in ddi_put, bad upper limits '),ihigh1,ihigh2,
c$$$cddi1     $      mdata(1,nhandle(1)),mdata(2,nhandle(1)),nhandle(1)          1d18s21
c$$$cddi1       stop
c$$$cddi1      end if
c$$$cddi1      nrow=ihigh1+1-ilow1
c$$$cddi1      npacket(3)=nrow
c$$$cddi1      npacket(5)=mymast                                                 1d29s21
c$$$cddi1      npacket(7)=ilow1
c$$$cddi1      npacket(8)=0
c$$$cddi1      itag=1
c$$$cddi1      jtag=10
c$$$cddi1      do iproc=0,isize-1
c$$$cddi1       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
c$$$cddi1     $      ilj,ihj,i1s,i1e,i2s,i2e)
c$$$cddi1       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
c$$$cddi1        istrt=max(ilow2,ilj)
c$$$cddi1        iend=min(ihigh2,ihj)
c$$$cddi1        npacket(4)=istrt
c$$$cddi1        npacket(5)=iend
c$$$cddi1        npacket(8)=ilj
c$$$cddi1        ito=iranks(iproc+1)                                             1d29s21
c$$$cddi1        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
c$$$cddi1         i1=1+nrow*(istrt-ilow2)                                        3d6s09
c$$$cddi1         iaddx=mdata(1,nhandle(1))*(istrt-ilj)+ilow1-1                  1d18s21
c$$$cddi1     $        +mdata(3,nhandle(1))                                      1d18s21
c$$$cddi1         nmove=nrow*(iend+1-istrt)                                      3d6s09
c$$$cddi138835    format(4i8,1pe15.7)
c$$$cddi1         itmp=ibcoff
c$$$cddi1         ibcoff=itmp+4+nmove
c$$$cddi1         jtmp=itmp-1
c$$$cddi1         do j=1,4
c$$$cddi1          bc(jtmp+j)=packet(j)
c$$$cddi1         end do
c$$$cddi1         jtmp=jtmp+5
c$$$cddi1         do j=0,nmove-1
c$$$cddi1          bc(jtmp+j)=buff(i1+j)
c$$$cddi1         end do
c$$$cddi1         nmovep=nmove+4
c$$$cddi1         call mpi_send(bc(itmp),nmovep,mpi_double_precision,ito,itag,
c$$$cddi1     $        mpi_comm_world,ierror)
c$$$cddi1         ibcoff=itmp
c$$$cddi1        else                                                            3d6s09
c$$$cddi1         nwds=(iend+1-istrt)*nrow
c$$$cddi1         if(nwds.gt.0)then
c$$$cddi1         itmp=ibcoff
c$$$cddi1         ibcoff=itmp+nwds+4
c$$$cddi1         call enough('ddi_createx.  1',bc,ibc)
c$$$cddi1         jtmp=itmp-1
c$$$cddi1         do j=1,4
c$$$cddi1          bc(jtmp+j)=packet(j)
c$$$cddi1         end do
c$$$cddi1         jtmp=jtmp+5
c$$$cddi1         do i=istrt,iend
c$$$cddi1          i1=1+nrow*(i-ilow2)
c$$$cddi1          do j=0,nrow-1
c$$$cddi1           bc(jtmp+j)=buff(i1+j)
c$$$cddi1          end do
c$$$cddi1          jtmp=jtmp+nrow
c$$$cddi1c          iaddx=mdata(1,nhandle)*(i-ilj)+ilow1-1+mdata(3,nhandle)        3d6s09
c$$$cddi1   33    format(5i8)
c$$$cddi1         end do
c$$$cddi1         nwdsp=nwds+4
c$$$cddi1         call mpi_send(bc(itmp),nwdsp,mpi_double_precision,ito,
c$$$cddi1     $        itag,mpi_comm_world,ierror)
c$$$cddi1         ibcoff=itmp
c$$$cddi1         end if
c$$$cddi1        end if
c$$$cddi1       end if                                                           3d6s09
c$$$cddi1      end do
c$$$cddi1      return
c$$$cddi1      entry ddi_accx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff)     11d15s22
c$$$cddi1      if(nhand.gt.ndata.or.nhandle(1).lt.1)then
c$$$cddi1       write(6,*)('in ddi_accx, handle = '),nhandle(1),nhandlex,nhand
c$$$cddi1       write(6,*)('exceeds ndata = '),ndata
c$$$cddi1       stop
c$$$cddi1      end if
c$$$cddi1      npacket(2)=nhandle(1)                                             1d18s21
c$$$cddi1      if(min(ilow1,ilow2).le.0)then
c$$$cddi1       write(6,*)('in ddi_put, bad lower limits '),ilow1,ilow2
c$$$cddi1       stop
c$$$cddi1      end if
c$$$cddi1      if(ihigh1.gt.mdata(1,nhandle(1)).or.                              1d18s21
c$$$cddi1     $     ihigh2.gt.mdata(2,nhandle(1)))then                           1d18s21
c$$$cddi1       write(6,*)('in ddi_accx, bad upper limits '),ihigh1,ihigh2,
c$$$cddi1     $      mdata(1,nhandle(1)),mdata(2,nhandle(1)),nhandle(1)          1d18s21
c$$$cddi1       stop
c$$$cddi1      end if
c$$$cddi1      nrow=ihigh1+1-ilow1
c$$$cddi1      npacket(3)=nrow
c$$$cddi1      npacket(6)=mymast                                                 1d29s21
c$$$cddi1      npacket(7)=ilow1
c$$$cddi1      npacket(1)=4                                                      2d21s21
c$$$cddi1      itag=1
c$$$cddi1      jtag=11
c$$$cddi1      ktag=13                                                           2d19s21
c$$$cddi1      do iproc=0,isize-1
c$$$cddi1       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
c$$$cddi1     $      ilj,ihj,i1s,i1e,i2s,i2e)
c$$$cddi1       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
c$$$cddi1        istrt=max(ilow2,ilj)
c$$$cddi1        iend=min(ihigh2,ihj)
c$$$cddi1        npacket(4)=istrt
c$$$cddi1        npacket(5)=iend
c$$$cddi1        npacket(8)=ilj
c$$$cddi1        ito=iranks(iproc+1)                                             1d29s21
c$$$cddi1        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
c$$$cddi1         i1=1+nrow*(istrt-ilow2)                                        3d6s09
c$$$cddi1         iaddx=mdata(1,nhandle(1))*(istrt-ilj)                          1d18s21
c$$$cddi1     $        +ilow1-1+mdata(3,nhandle(1))                              1d18s21
c$$$cddi1         nmove=nrow*(iend+1-istrt)                                      3d6s09
c$$$cddi1         if(nmove.gt.0)then
c$$$cddi1          itmp=ibcoff
c$$$cddi1          ibcoff=itmp+nmove+4
c$$$cddi1          jtmp=itmp-1
c$$$cddi1          do j=1,4
c$$$cddi1           bc(jtmp+j)=packet(j)
c$$$cddi1          end do
c$$$cddi1          jtmp=jtmp+5
c$$$cddi1          do j=0,nmove-1
c$$$cddi1           bc(jtmp+j)=buff(i1+j)
c$$$cddi1          end do
c$$$cddi1          nmovep=nmove+4
c$$$cddi1         call mpi_ssend(bc(itmp),nmovep,mpi_double_precision,ito,itag,  2d5s21
c$$$cddi1     $        mpi_comm_world,ierror)
c$$$cddi1          ibcoff=itmp
c$$$cddi1         end if
c$$$cddi1        else                                                            3d6s09
c$$$cddi1         nwds=(iend+1-istrt)*nrow
c$$$cddi1         if(nwds.gt.0)then
c$$$cddi1         itmp=ibcoff
c$$$cddi1         ibcoff=itmp+nwds+4
c$$$cddi1         call enough('ddi_createx.  2',bc,ibc)
c$$$cddi1         jtmp=itmp-1
c$$$cddi1         do j=1,4
c$$$cddi1          bc(jtmp+j)=packet(j)
c$$$cddi1         end do
c$$$cddi1         jtmp=jtmp+5
c$$$cddi1         do i=istrt,iend
c$$$cddi1          i1=1+nrow*(i-ilow2)
c$$$cddi1          do j=0,nrow-1
c$$$cddi1           bc(jtmp+j)=buff(i1+j)
c$$$cddi1          end do
c$$$cddi1          jtmp=jtmp+nrow
c$$$cddi1         end do
c$$$cddi1         nwdsp=nwds+4
c$$$cddi1         call mpi_ssend(bc(itmp),nwdsp,mpi_double_precision,ito,itag,   2d5s21
c$$$cddi1     $        mpi_comm_world,ierror)
c$$$cddi1         ibcoff=itmp
c$$$cddi1         end if
c$$$cddi1        end if                                                          3d6s09
c$$$cddi1       end if
c$$$cddi1      end do
c$$$cddi1      return
c$$$cddi1      entry ddi_iaccx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff,
c$$$cddi1     $  ircv,nrcv)                                                      11d15s22
c$$$cddi1c
c$$$cddi1c     important note!!!
c$$$cddi1c     buff needs to have 4*mynprocg words of extra storage for this to
c$$$cddi1c     work properly, yet there is no way to test for this !!!
c$$$cddi1c
c$$$cddi1      if(nhand.gt.ndata.or.nhandle(1).lt.1)then
c$$$cddi1       write(6,*)('in ddi_accx, handle = '),nhandle(1),nhandlex,nhand
c$$$cddi1       write(6,*)('exceeds ndata = '),ndata
c$$$cddi1       stop
c$$$cddi1      end if
c$$$cddi1      npacket(2)=nhandle(1)                                             1d18s21
c$$$cddi1      if(min(ilow1,ilow2).le.0)then
c$$$cddi1       write(6,*)('in ddi_put, bad lower limits '),ilow1,ilow2
c$$$cddi1       stop
c$$$cddi1      end if
c$$$cddi1      if(ihigh1.gt.mdata(1,nhandle(1)).or.                              1d18s21
c$$$cddi1     $     ihigh2.gt.mdata(2,nhandle(1)))then                           1d18s21
c$$$cddi1       write(6,*)('in ddi_accx, bad upper limits '),ihigh1,ihigh2,
c$$$cddi1     $      mdata(1,nhandle(1)),mdata(2,nhandle(1)),nhandle(1)          1d18s21
c$$$cddi1       stop
c$$$cddi1      end if
c$$$cddi1      nrow=ihigh1+1-ilow1
c$$$cddi1      npacket(3)=nrow
c$$$cddi1      npacket(6)=mymast                                                 1d29s21
c$$$cddi1      npacket(7)=ilow1
c$$$cddi1      npacket(1)=10
c$$$cddi1      itag=1
c$$$cddi1      jtag=11
c$$$cddi1      nrcv=0                                                            1d30s21
c$$$cddi1      itmp0=ibcoff                                                      1d30s21
c$$$cddi1      jtmp0=itmp0                                                       1d30s21
c$$$cddi1      do iproc=0,isize-1
c$$$cddi1       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
c$$$cddi1     $      ilj,ihj,i1s,i1e,i2s,i2e)
c$$$cddi1       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
c$$$cddi1        istrt=max(ilow2,ilj)
c$$$cddi1        iend=min(ihigh2,ihj)
c$$$cddi1        npacket(4)=istrt
c$$$cddi1        npacket(5)=iend
c$$$cddi1        npacket(8)=ilj
c$$$cddi1        ito=iranks(iproc+1)                                             1d29s21
c$$$cddi1        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
c$$$cddi1         i1=1+nrow*(istrt-ilow2)                                        3d6s09
c$$$cddi1         iaddx=mdata(1,nhandle(1))*(istrt-ilj)                          1d18s21
c$$$cddi1     $        +ilow1-1+mdata(3,nhandle(1))                              1d18s21
c$$$cddi1         nmove=nrow*(iend+1-istrt)                                      3d6s09
c$$$cddi1         if(nmove.gt.0)then
c$$$cddi1          itmp=jtmp0                                                    1d30s21
c$$$cddi1          jtmp0=itmp+nmove+4                                            1d30s21
c$$$cddi1          ibcoff=jtmp0                                                   1d30s21
c$$$cddi1          call enough('ddi_createx.  3',bc,ibc)
c$$$cddi1          jtmp=itmp-1
c$$$cddi1          do j=1,4
c$$$cddi1           bc(jtmp+j)=packet(j)
c$$$cddi1          end do
c$$$cddi1          jtmp=jtmp+5
c$$$cddi1          do j=0,nmove-1
c$$$cddi1           bc(jtmp+j)=buff(i1+j)
c$$$cddi1          end do
c$$$cddi1         end if
c$$$cddi1        else                                                            3d6s09
c$$$cddi1         nwds=(iend+1-istrt)*nrow
c$$$cddi1         if(nwds.gt.0)then
c$$$cddi1         itmp=jtmp0                                                     1d30s21
c$$$cddi1         jtmp0=itmp+nwds+4                                              1d30s21
c$$$cddi1         ibcoff=jtmp0                                                   1d30s21
c$$$cddi1         call enough('ddi_createx.  4',bc,ibc)
c$$$cddi1         jtmp=itmp-1
c$$$cddi1         do j=1,4
c$$$cddi1          bc(jtmp+j)=packet(j)
c$$$cddi1         end do
c$$$cddi1         jtmp=jtmp+5
c$$$cddi1         do i=istrt,iend
c$$$cddi1          i1=1+nrow*(i-ilow2)
c$$$cddi1          do j=0,nrow-1
c$$$cddi1           bc(jtmp+j)=buff(i1+j)
c$$$cddi1          end do
c$$$cddi1          jtmp=jtmp+nrow
c$$$cddi1         end do
c$$$cddi1         nwdsp=nwds+4
c$$$cddi1         end if
c$$$cddi1        end if                                                          3d6s09
c$$$cddi1       end if
c$$$cddi1      end do
c$$$cddi1      ntotal=jtmp0-itmp0                                                1d30s21
c$$$cddi1      jtmp0=itmp0-1                                                     1d30s21
c$$$cddi1      do i=1,ntotal                                                     1d30s21
c$$$cddi1       buff(i)=bc(jtmp0+i)                                              1d30s21
c$$$cddi1      end do                                                            1d30s21
c$$$cddi1      ibcoff=itmp0                                                      1d30s21
c$$$cddi1      i1=1
c$$$cddi1      myall=0
c$$$cddi1      do iproc=0,isize-1
c$$$cddi1       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
c$$$cddi1     $      ilj,ihj,i1s,i1e,i2s,i2e)
c$$$cddi1       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
c$$$cddi1        istrt=max(ilow2,ilj)
c$$$cddi1        iend=min(ihigh2,ihj)
c$$$cddi1        npacket(4)=istrt
c$$$cddi1        npacket(5)=iend
c$$$cddi1        npacket(8)=ilj
c$$$cddi1        ito=iranks(iproc+1)                                             1d29s21
c$$$cddi1        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
c$$$cddi1         nmove=nrow*(iend+1-istrt)                                      3d6s09
c$$$cddi1         if(nmove.gt.0)then
c$$$cddi1          nmovep=nmove+4
c$$$cddi1          nrcv=nrcv+1                                                   1d30s21
c$$$cddi1          myall=myall+nmove
c$$$cddi1          call mpi_isend(buff(i1),nmovep,mpi_double_precision,ito,itag, 6d17s21
c$$$cddi1     $        mpi_comm_world,ircv(nrcv),ierror)                         1d30s21
c$$$cddi1          i1=i1+nmovep                                                  1d30s21
c$$$cddi1         end if
c$$$cddi1        else                                                            3d6s09
c$$$cddi1         nwds=(iend+1-istrt)*nrow
c$$$cddi1         if(nwds.gt.0)then
c$$$cddi1          nwdsp=nwds+4
c$$$cddi1          nrcv=nrcv+1                                                   1d30s21
c$$$cddi1          myall=myall+nmove
c$$$cddi1          call mpi_isend(buff(i1),nwdsp,mpi_double_precision,ito,itag,  6d17s21
c$$$cddi1     $        mpi_comm_world,ircv(nrcv),ierror)                         1d30s21
c$$$cddi1          i1=i1+nwdsp                                                   1d30s21
c$$$cddi1         end if
c$$$cddi1        end if                                                          3d6s09
c$$$cddi1       end if
c$$$cddi1      end do
c$$$cddi1      return
c$$$cddi1      entry ddi_getx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff)     11d15s22
c$$$cddi1      if(nhandle(1).gt.ndata.or.nhandle(1).lt.1)then                    1d18s21
c$$$cddi1       write(6,*)('in ddi_get, handle = '),nhandle(1),nhandlex          1d18s21
c$$$cddi1       write(6,*)('exceeds ndata = '),ndata
c$$$cddi1       stop
c$$$cddi1      end if
c$$$cddi1      npacket(2)=nhandle(1)                                             1d18s21
c$$$cddi1      if(min(ilow1,ilow2).le.0)then
c$$$cddi1       write(6,*)('in ddi_get, bad lower limits '),ilow1,ilow2
c$$$cddi1       stop
c$$$cddi1      end if
c$$$cddi1      if(ihigh1.gt.mdata(1,nhandle(1)).or.                              1d18s21
c$$$cddi1     $     ihigh2.gt.mdata(2,nhandle(1)))then                           1d18s21
c$$$cddi1       write(6,*)('in ddi_get, bad upper limits '),ihigh1,ihigh2
c$$$cddi1       write(6,*)('vs '),mdata(1,nhandle(1)),mdata(2,nhandle(1))        1d18s21
c$$$cddi1       write(6,*)('nhandle = '),nhandle(1)                              1d18s21
c$$$cddi1       stop
c$$$cddi1      end if
c$$$cddi1      nrow=ihigh1+1-ilow1
c$$$cddi1      npacket(3)=nrow
c$$$cddi1      npacket(6)=mymast                                                 1d29s21
c$$$cddi1      npacket(7)=ilow1
c$$$cddi1      npacket(1)=5
c$$$cddi1      itag=1
c$$$cddi1      jtag=12
c$$$cddi1      do iproc=0,isize-1
c$$$cddi1       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
c$$$cddi1     $      ilj,ihj,i1s,i1e,i2s,i2e)
c$$$cddi1       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
c$$$cddi1        istrt=max(ilow2,ilj)
c$$$cddi1        iend=min(ihigh2,ihj)
c$$$cddi1        npacket(4)=istrt
c$$$cddi1        npacket(5)=iend
c$$$cddi1        npacket(8)=ilj
c$$$cddi1        ito=iranks(iproc+1)                                             1d29s21
c$$$cddi1        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
c$$$cddi1         i1=1+nrow*(istrt-ilow2)                                        3d6s09
c$$$cddi1         iaddx=mdata(1,nhandle(1))*(istrt-ilj)                          1d18s21
c$$$cddi1     $        +ilow1-1+mdata(3,nhandle(1))                              1d18s21
c$$$cddi1         nmove=nrow*(iend+1-istrt)                                      3d6s09
c$$$cddi1         if(nmove.gt.0)then
c$$$cddi1         call mpi_send(packet,4,mpi_double_precision,ito,itag,
c$$$cddi1     $         mpi_comm_world,ierror)
c$$$cddi1         call mpi_recv(buff(i1),nmove,mpi_double_precision,ito,jtag,
c$$$cddi1     $        mpi_comm_world,istatus,ierror)
c$$$cddi1         end if
c$$$cddi1        else                                                            3d6s09
c$$$cddi1         nwds=(iend+1-istrt)*nrow
c$$$cddi1         if(nwds.gt.0)then
c$$$cddi1         itmp=ibcoff
c$$$cddi1         ibcoff=itmp+nwds
c$$$cddi1         call enough('ddi_createx.  5',bc,ibc)
c$$$cddi1         call mpi_send(packet,4,mpi_double_precision,ito,itag,
c$$$cddi1     $         mpi_comm_world,
c$$$cddi1     $       ierror)
c$$$cddi1         call mpi_recv(bc(itmp),nwds,mpi_double_precision,ito,jtag,
c$$$cddi1     $        mpi_comm_world,istatus,ierror)
c$$$cddi1         jtmp=itmp
c$$$cddi1         do i=istrt,iend
c$$$cddi1          i1=1+nrow*(i-ilow2)
c$$$cddi1          do j=0,nrow-1
c$$$cddi1           buff(i1+j)=bc(jtmp+j)
c$$$cddi1          end do
c$$$cddi1          jtmp=jtmp+nrow
c$$$cddi1         end do
c$$$cddi1         ibcoff=itmp
c$$$cddi1         end if
c$$$cddi1       end if
c$$$cddi1       end if
c$$$cddi1      end do
c$$$cddi1      return
c$$$cddi1      entry ddi_igetx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,         11d15s22
c$$$cddi1     $   buff,ircv,nrcv)                                                11d15s22
c$$$cddi1      if(nhandle(1).gt.ndata.or.nhandle(1).lt.1)then                    1d18s21
c$$$cddi1       write(6,*)('in ddi_igetx, handle = '),nhandle(1),nhandlex          1d18s21
c$$$cddi1       write(6,*)('exceeds ndata = '),ndata
c$$$cddi1       stop
c$$$cddi1      end if
c$$$cddi1      npacket(2)=nhandle(1)                                             1d18s21
c$$$cddi1      if(min(ilow1,ilow2).le.0)then
c$$$cddi1       write(6,*)('in ddi_igetx, bad lower limits '),ilow1,ilow2
c$$$cddi1       stop
c$$$cddi1      end if
c$$$cddi1      if(ihigh1.gt.mdata(1,nhandle(1)).or.                              1d18s21
c$$$cddi1     $     ihigh2.gt.mdata(2,nhandle(1)))then                           1d18s21
c$$$cddi1       write(6,*)('in ddi_iget, bad upper limits '),ihigh1,ihigh2
c$$$cddi1       write(6,*)('vs '),mdata(1,nhandle(1)),mdata(2,nhandle(1))        1d18s21
c$$$cddi1       write(6,*)('nhandle = '),nhandle(1)                              1d18s21
c$$$cddi1       stop
c$$$cddi1      end if
c$$$cddi1      nrow=ihigh1+1-ilow1
c$$$cddi1      npacket(3)=nrow
c$$$cddi1      npacket(6)=mymast                                                 1d29s21
c$$$cddi1      npacket(7)=ilow1
c$$$cddi1      npacket(1)=5
c$$$cddi1      itag=1
c$$$cddi1      jtag=12
c$$$cddi1      nrcv=0                                                            1d29s21
c$$$cddi1      do iproc=0,isize-1
c$$$cddi1       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
c$$$cddi1     $      ilj,ihj,i1s,i1e,i2s,i2e)
c$$$cddi1       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
c$$$cddi1        istrt=max(ilow2,ilj)
c$$$cddi1        iend=min(ihigh2,ihj)
c$$$cddi1        npacket(4)=istrt
c$$$cddi1        npacket(5)=iend
c$$$cddi1        npacket(8)=ilj
c$$$cddi1        ito=iranks(iproc+1)                                             1d29s21
c$$$cddi1        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
c$$$cddi1         i1=1+nrow*(istrt-ilow2)                                        3d6s09
c$$$cddi1         iaddx=mdata(1,nhandle(1))*(istrt-ilj)                          1d18s21
c$$$cddi1     $        +ilow1-1+mdata(3,nhandle(1))                              1d18s21
c$$$cddi1         nmove=nrow*(iend+1-istrt)                                      3d6s09
c$$$cddi1         if(nmove.gt.0)then
c$$$cddi1         call mpi_send(packet,4,mpi_double_precision,ito,itag,
c$$$cddi1     $         mpi_comm_world,ierror)
c$$$cddi1         nrcv=nrcv+1                                                    1d29s21
c$$$cddi1         call mpi_irecv(buff(i1),nmove,mpi_double_precision,ito,jtag,   1d29s21
c$$$cddi1     $        mpi_comm_world,ircv(nrcv),ierror)
c$$$cddi1         end if
c$$$cddi1        else                                                            3d6s09
c$$$cddi1         nwds=(iend+1-istrt)*nrow
c$$$cddi1         if(nwds.gt.0)then
c$$$cddi1         itmp=ibcoff
c$$$cddi1         ibcoff=itmp+nwds
c$$$cddi1         call enough('ddi_createx.  6',bc,ibc)
c$$$cddi1         call mpi_send(packet,4,mpi_double_precision,ito,itag,
c$$$cddi1     $         mpi_comm_world,
c$$$cddi1     $       ierror)
c$$$cddi1         call mpi_recv(bc(itmp),nwds,mpi_double_precision,ito,jtag,     1d29s21
c$$$cddi1     $        mpi_comm_world,istatus,ierror)                            1d29s21
c$$$cddi1         jtmp=itmp
c$$$cddi1         do i=istrt,iend
c$$$cddi1          i1=1+nrow*(i-ilow2)
c$$$cddi1          do j=0,nrow-1
c$$$cddi1           buff(i1+j)=bc(jtmp+j)
c$$$cddi1          end do
c$$$cddi1          jtmp=jtmp+nrow
c$$$cddi1         end do
c$$$cddi1         ibcoff=itmp
c$$$cddi1         end if
c$$$cddi1       end if
c$$$cddi1       end if
c$$$cddi1      end do
c$$$cddi1      return
cddi1      end
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine ddi_iget(bc,ibc,nhandlex,ilow1,ihigh1,ilow2,ihigh2,    11d15s22
cddi1     $   buff,ircv,nrcv)                                                11d15s22
cddi1c
cddi1c     non-blocking version
cddi1c
cddi1      implicit real*8 (a-h,o-z)                                         1d18s21
cddi1      implicit integer*8 (i-n)                                          1d18s21
cddi1      dimension buff(*),ircv(*)                                         1d29s21
cddi1      nhandle=nhandlex+1                                                1d18s21
cddi1      call ddi_igetx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff,ircv,11d15s22
cddi1     $   nrcv)                                                          11d15s22
cddi1      return                                                            1d18s21
cddi1      end                                                               1d18s21
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine ddi_iaccword1(iarg)                                    7d18s24
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer (i-n)
cddi1      integer*8 iarg,icol,nhandle,istart                                2d4s15
cddi1      integer (kind=mpi_address_kind) memreq,iadd,iaddx,ntrial          3d27s09
cddi1      dimension istatus(mpi_status_size),npacket(2)
cddi1      include "common.mympi"                                            1d29s21
cddi1      data nhandle/-132/
cddi1      equivalence (packet,npacket)
cddi1      save
cddi1      itag=1
cddi1      call mpi_comm_rank(my_comm_group,irank,ierror)
cddi1      ito=mymemp                                                        8d10s22
cddi1      npacket(1)=12                                                     8d10s22
cddi1      npacket(2)=mymast                                                 1d29s21
cddi1      call mpi_send(packet,1,mpi_double_precision,ito,itag,
cddi1     $     mpi_comm_world,ierror)
cddi1      itag=10
c$$$cddi1        write(6,*)('mpi_recvf ')
cddi1      call mpi_recv(icount,1,mpi_integer,ito,itag,mpi_comm_world,
cddi1     $     istatus,ierror)
cddi1      if(ierror.ne.mpi_success)then
cddi1       write(6,*)('mpi_recv returned an error ')
cddi1       write(6,*)('ierror '),ierror
cddi1       write(6,*)('mpi_success '),mpi_success
cddi1       stop 'recvf'
cddi1      end if
cddi1      iarg=icount
cddi1      return
cddi1      end                                                               7d18s24
cddi1      subroutine ddi_iaccword0                                          8d10s22
cddi1      use mpi
cddi1      implicit real*8 (a-h,o-z)
cddi1      implicit integer (i-n)
cddi1      integer*8 iarg,icol,nhandle,istart                                2d4s15
cddi1      integer (kind=mpi_address_kind) memreq,iadd,iaddx,ntrial          3d27s09
cddi1      dimension istatus(mpi_status_size),npacket(2)
cddi1      include "common.mympi"                                            1d29s21
cddi1      data nhandle/-132/
cddi1      equivalence (packet,npacket)
cddi1      save
cddi1      call mpi_comm_rank(my_comm_group,irank,ierror)
cddi1      call mpi_comm_size(my_comm_group,isize,ierror)
cddi1      itag=1
cddi1      ito=mymemp
cddi1      idum=11                                                           8d10s22
cddi1      npacket(1)=idum
cddi1      call mpi_send(packet,1,mpi_double_precision,ito,itag,
cddi1     $      mpi_comm_world,ierror)
cddi1      call dws_sync
cddi1      return
c$$$cddi1      entry ddi_iaccword1(iarg)
c$$$cddi1      itag=1
c$$$cddi1      call mpi_comm_rank(my_comm_group,irank,ierror)
c$$$cddi1      ito=mymemp                                                        8d10s22
c$$$cddi1      npacket(1)=12                                                     8d10s22
c$$$cddi1      npacket(2)=mymast                                                 1d29s21
c$$$cddi1      call mpi_send(packet,1,mpi_double_precision,ito,itag,
c$$$cddi1     $     mpi_comm_world,ierror)
c$$$cddi1      itag=10
c$$$cddi1      call mpi_recv(icount,1,mpi_integer,ito,itag,mpi_comm_world,
c$$$cddi1     $     istatus,ierror)
c$$$cddi1      iarg=icount
c$$$cddi1      return
cddi1      end
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine ddi_iaccword2(nwds)                                    8d10s22
cddi1      implicit real*8 (a-h,o-z)                                         8d10s22
cddi1      integer*8 iwds                                                    8d10s22
cddi1      dimension snd(2)                                                  8d10s22
cddi1      snd(1)=dfloat(nwds)                                               8d10s22
cddi1      call ddi_iaccword1(iwds)                                          8d10s22
cddi1      snd(2)=dfloat(iwds)                                               8d10s22
cddi1      call dws_gsumf(snd,2)                                             8d10s22
cddi1      nwant=nint(snd(1))                                                8d10s22
cddi1      ngot=nint(snd(2))                                                 8d10s22
cddi1      if(nwant.ne.ngot)then                                             8d10s22
cddi1       looper=0                                                         8d10s22
cddi1    1  continue                                                         8d10s22
cddi1       looper=looper+1                                                  8d10s22
cddi1       if(looper.gt.10)then                                             8d10s22
cddi1        write(6,*)('we waited toooo long for iacc to complete! ')       6d22s23
cddi1        write(6,*)('nwant: '),nwant                                     6d22s23
cddi1        write(6,*)('ngot : '),ngot                                      6d22s23
cddi1        call dws_synca
cddi1        call dws_finalize
cddi1        stop 'ddi_iaccword2'                                            8d10s22
cddi1       end if                                                           8d10s22
cddi1       call sleep(1)                                                    8d10s22
cddi1       call ddi_iaccword1(iwds)                                          8d10s22
cddi1       snd(2)=dfloat(iwds)                                               8d10s22
cddi1       call dws_gsumf(snd(2),1)                                         8d10s22
cddi1       ngot=nint(snd(2))                                                8d10s22
cddi1c
cddi1c     keep giving reprieves as long as we are still getting more data
cddi1c       
cddi1       if(looper.eq.1)then                                              6d23s23
cddi1        nlast=ngot                                                      6d23s23
cddi1       else                                                             6d23s23
cddi1        if(nlast.ne.ngot)then                                           6d23s23
cddi1         looper=1                                                       6d23s23
cddi1         nlast=ngot                                                     6d23s23
cddi1        end if                                                          6d23s23
cddi1       end if                                                           6d23s23
cddi1       if(nwant.ne.ngot)go to 1                                         8d10s22
cddi1      end if                                                            8d10s22
cddi1      return                                                            8d10s22
cddi1      end                                                               8d10s22
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine ddi_get(bc,ibc,nhandlex,ilow1,ihigh1,ilow2,ihigh2,buff)11d15s22
cddi1      implicit real*8 (a-h,o-z)                                         1d18s21
cddi1      implicit integer*8 (i-n)                                          1d18s21
cddi1      dimension buff(*)                                                 1d18s21
cddi1      nhandle=nhandlex+1                                                1d18s21
cddi1      call ddi_getx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff)      11d15s22
cddi1      return                                                            1d18s21
cddi1      end                                                               1d18s21
cddi1c mpec2.1 version zeta copyright u.s. government
cddi1      subroutine ddi_done(ircv,nrcv)                                    1d29s21
cddi1      use mpi                                                           1d29s21
cddi1      parameter (id=10)
cddi1      dimension ircv(*),istatus(id)                                      1d29s21
cddi1      do i=1,nrcv                                                       1d29s21
cddi1       call mpi_wait(ircv(i),istatus,ierror)                             1d29s21
cddi1      end do                                                            1d29s21
cddi1      nrcv=0                                                            1d29s21
cddi1      return                                                            1d29s21
cddi1      end                                                               1d29s21
cddi1
c
c     parallelization via mpi.
c     use one-sided memory access via mpi-3 commands.
c     I think this will be limited to a single node.
c
      subroutine dws_init(ncore)                                        3d5s21
      idum=1                                                            3d5s21
      return                                                            3d5s21
      end                                                               3d5s21
      subroutine dws_preinit                                            3d5s21
      use mpi
      implicit real*8 (a-h,o-z)                                         4d10s12
      implicit integer (i-n)
      character*8 file
      integer (kind=mpi_address_kind) isizex                            6d4s10
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      include "common.mympi"                                             3d5s21
      integer*8 itmp
      call mpi_init(ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_int,')
       write(6,*)('mpi_init returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      call mpi_comm_size(MPI_COMM_WORLD,isize,ierror)
      call mpi_comm_rank(MPI_COMM_WORLD,irank,ierror)
      if(irank.eq.0)write(6,*)('Hi, this is preinit for MPI-4 code')    11d1s22
c$$$      if(irank.eq.0)then                                                6d4s10
c$$$       read(5,*)ncpus                                                   6d4s10
c$$$       xcpus=dfloat(ncpus)                                              6d4s10
c$$$      end if                                                            6d4s10
      mynowprog=irank                                                   3d5s21
      mynprocg=isize                                                    3d5s21
      ndata=0                                                           3d5s21
      if(irank.lt.10)then
       write(file,1)irank
    1  format('output.',i1)
      else if(irank.lt.100)then
       write(file,2)irank
    2  format('outpu.',i2)
      else
       write(file,3)irank
    3  format('outp.',i3)
      end if
      if(irank.ne.0)then
      open(unit=6,file=file)
      end if
      return
      end
      subroutine dws_finalize
      use mpi
      implicit integer (i-n)
      call mpi_finalize(ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in dws_finalize,')
       write(6,*)('mpi_finalize returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
      subroutine dws_gbor(buff,nwds)                                    1d26s21
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer buff(nwds)                                                1d26s21
      include "common.mympi"                                            1d29s21
c$$$      common/mpicom/my_comm_group,mymemp
c$$$      write(6,*)('my gsumf input buffer '),buff
      nwds2=nwds*2                                                      1d26s21
      call mpi_allreduce(mpi_in_place,buff,nwds2,mpi_integer,           1d26s21                                                                2d22s10
     $     mpi_bor,MPI_COMM_WORLD,ierror)                                1d26s21
c$$$      write(6,*)('my gsumf output buffer '),buff
      if(ierror.ne.mpi_success)then
       write(6,*)('in dws_gbor,')
       write(6,*)('mpi_allreduce returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
      subroutine dws_bcast(buf,len)
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer icont
      dimension buf(len)
      entry dws_bcasta(buf,len)                                         3d5s21
      icont=len
      isource=0
c$$$      write(6,*)('in dws_bcast for len '),len,loc(buf)
      call mpi_bcast(buf,icont,mpi_double_precision,isource,            2d19s10
     $               MPI_COMM_WORLD,ierror)                             2d19s10
      if(ierror.ne.mpi_success)then
       write(6,*)('in dws_bcast,')
       write(6,*)('mpi_bcast returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
      subroutine dws_sync
      use mpi
      entry dws_synca                                                   3d5s21
      call mpi_barrier(MPI_COMM_WORLD,ierror)                           2d22s10
      if(ierror.ne.mpi_success)then
       write(6,*)('in dws_sync,')
       write(6,*)('mpi_barrier returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
      subroutine dws_gsumf(buff,nwds)
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      dimension buff(nwds)                                                 2d22s19
c$$$      write(6,*)('my gsumf input buffer '),buff
      call mpi_allreduce(mpi_in_place,buff,nwds,mpi_double_precision,   2d22s10
     $     mpi_sum,MPI_COMM_WORLD,ierror)                               2d22s10
c$$$      write(6,*)('my gsumf output buffer '),buff
      if(ierror.ne.mpi_success)then
       write(6,*)('in dws_gsumf,')
       write(6,*)('mpi_allreduce returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
      subroutine dws_win0(mstor,buff,iwin)
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer (kind=mpi_address_kind) size
      common/mycomu/my_comml
      call mpi_info_create(info,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in dws_win0,')
       write(6,*)('mpi_info_create returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
c$$$      call mpi_info_set(info,'no_locks','true',ierror)
      call mpi_info_set(info,'no_locks','false',ierror)                                                                        
      if(ierror.ne.mpi_success)then
       write(6,*)('in dws_win0_create,')
       write(6,*)('mpi_info_set returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      istepsz=8                                                         
      size=mstor*istepsz
      call mpi_win_create(buff,size,istepsz,info,
     $     my_comml,iwin,ierror)                                        3d30s09
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_create,')
       write(6,*)('mpi_win_create returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
      subroutine dws_put0(buff,ioff,nwds,iwin)
      use mpi
      implicit real*8 (a-h,o-z)
      integer (kind=mpi_address_kind) ioffk                             6d9s10
      common/mycomu/my_comml
      itarg=0
      ioffk=ioff                                                        6d9s10
      call mpi_put(buff,nwds,mpi_double_precision,itarg,ioffk,nwds,     6d9s10
     $     mpi_double_precision,iwin,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_put0,')
       write(6,*)('mpi_win_create returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
      subroutine dws_get0(buff,ioff,nwds,iwin)
      use mpi
      implicit real*8 (a-h,o-z)
      integer (kind=mpi_address_kind) ioffk                             6d9s10
      common/mycomu/my_comml
      itarg=0
      ioffk=ioff                                                        6d9s10
      call mpi_get(buff,nwds,mpi_double_precision,itarg,ioffk,nwds,     6d9s10
     $     mpi_double_precision,iwin,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_get0,')
       write(6,*)('mpi_get returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
      subroutine dws_put1(iproc,buff,nwds)                              11d13s12
      use mpi
      implicit real*8 (a-h,o-z)
      integer (kind=mpi_address_kind) ioffk                             6d9s10
      itarg=iproc                                                       11d13s12
      itag=0                                                            11d13s12
      call mpi_send(buff,nwds,mpi_double_precision,itarg,itag,          11d13s12
     $     MPI_COMM_WORLD,ierror)                                       11d13s12
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_put1,')
       write(6,*)('mpi_win_create returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
      subroutine dws_sput1(iproc,buff,nwds)                              11d13s12
      use mpi
      implicit real*8 (a-h,o-z)
      integer (kind=mpi_address_kind) ioffk                             6d9s10
      itarg=iproc                                                       11d13s12
      itag=0                                                            11d13s12
      call mpi_ssend(buff,nwds,mpi_double_precision,itarg,itag,          11d13s12
     $     MPI_COMM_WORLD,ierror)                                       11d13s12
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_sput1,')
       write(6,*)('mpi_win_create returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
      subroutine dws_get1(iproc,buff,nwds)                              11d13s12
      use mpi
      implicit real*8 (a-h,o-z)
      integer (kind=mpi_address_kind) ioffk                             6d9s10
      integer status(mpi_status_size)
      itarg=iproc                                                       11d13s12
      itag=0                                                            11d13s12
      call mpi_recv(buff,nwds,mpi_double_precision,itarg,itag,
     $     MPI_COMM_WORLD,status,ierror)                                       11d13s12                                                                        
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_get1,')
       write(6,*)('mpi_get returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
      subroutine dws_putl(iproc,buff,nwds)                              11d13s12
      use mpi
      implicit real*8 (a-h,o-z)
      integer (kind=mpi_address_kind) ioffk                             6d9s10
      common/mycomu/my_comml
      itarg=iproc                                                       11d13s12
      itag=0                                                            11d13s12
      call mpi_send(buff,nwds,mpi_double_precision,itarg,itag,          11d13s12
     $     my_comml,ierror)                                             11d16s12
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_putl,')
       write(6,*)('mpi_win_create returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
      subroutine dws_getl(iproc,buff,nwds)                              11d13s12
      use mpi
      implicit real*8 (a-h,o-z)
      integer (kind=mpi_address_kind) ioffk                             6d9s10
      integer status(mpi_status_size)
      common/mycomu/my_comml
      itarg=iproc                                                       11d13s12
      itag=0                                                            11d13s12
      call mpi_recv(buff,nwds,mpi_double_precision,itarg,itag,
     $     my_comml,status,ierror)                                      11d16s12
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_getl,')
       write(6,*)('mpi_get returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
      subroutine dws_nowin0(iwin)                                       12d13s11
      use mpi                                                           12d13s11
      implicit real*8 (a-h,o-z)                                         12d13s11
      call mpi_win_free(iwin,ierror)                                    12d13s11
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_nowin0,')
       write(6,*)('mpi_win_free returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
      
      subroutine dws_all2allvl(buf,nblock,ioff)
      use mpi
      implicit real*8 (a-h,o-z)
      dimension buf(1),nblock(1),ioff(1)
      common/mycomu/my_comml
c$$$      write(6,*)('in dws_all2allvl '),buf(1),nblock(1),ioff(1)
c$$$      if(buf(1).ne.-132d0)then
c$$$       call dws_sync
c$$$       call dws_finalize
c$$$       stop
c$$$      end if
      call mpi_alltoallv(mpi_in_place,nblock,ioff,mpi_double_precision,
     $     buf,nblock,ioff,mpi_double_precision,mpi_comm_world,ierr)    6d28s12
      if(ierr.ne.mpi_success)then
       write(6,*)('in dws_all2allvl, mpi_alltoallv returned an error'),
     $      ierror
       close(unit=6)
c$$$       stop
      end if
c$$$      write(6,*)('after words...')
c$$$      if(buf(1).ne.-132d0)then
c$$$       call dws_sync
c$$$       call dws_finalize
c$$$       stop
c$$$      end if
      return
      end
      subroutine dws_all2allvb(bufs,nblock8s,ioff8s,bufr,nblock8r,
     $     ioff8r)
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension bufs(1),nblock8s(1),ioff8s(1),bufr(1),nblock8r(1),
     $     ioff8r(1)
ccccccc      parameter (id=1000)
      include "common.mympi"                                            1d21s21
      dimension nblocks(id),ioffs(id),nblockr(id),ioffr(id)
      data icall/0/
      save
      icall=icall+1
c$$$      write(6,*)('in all2allvb for icall = '),icall
      if(mynprocg.gt.id)then
       write(6,*)('too many procs in dws_allgv!! '),mynprocg,id
       close(unit=6)
c$$$       stop
      end if
c$$$      write(6,*)('in dws_allgv, mynprocg = '),mynprocg
      do i=1,mynprocg
       nblocks(i)=nblock8s(i)
       ioffs(i)=ioff8s(i)
       nblockr(i)=nblock8r(i)
       ioffr(i)=ioff8r(i)
c$$$       write(6,3030)i,nblocks(i),ioffs(i),nblockr(i),ioffr(i)
 3030  format(5i8)
      end do
c$$$      if(icall.eq.1.and.bufs(132).ne.132d0)then
c$$$       call dws_sync
c$$$       call dws_finalize
c$$$       stop
c$$$      end if
c$$$      if(ioffs(3).ne.-132)then
c$$$       call dws_sync
c$$$       call dws_finalize
c$$$       stop
c$$$      end if
      call mpi_alltoallv(bufs,nblocks,ioffs,mpi_double_precision,
     $     bufr,nblockr,ioffr,mpi_double_precision,mpi_comm_world,ierr)
      if(ierr.ne.mpi_success)then
       write(6,*)('in dws_all2allv, mpi_alltoallv returned an error'),
     $      ierror
       close(unit=6)
c$$$       stop
      end if
      return
      end
      subroutine dws_all2allvb8(bufs,nblock8s,ioff8s,bufr,nblock8r,
     $     ioff8r)
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer*8 nblock8s,ioff8s,nblock8r,ioff8r                         1d10s18
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension bufs(1),nblock8s(1),ioff8s(1),bufr(1),nblock8r(1),
     $     ioff8r(1)
ccccccc      parameter (id=1000)
      include "common.mympi"                                            1d21s21
      dimension nblocks(id),ioffs(id),nblockr(id),ioffr(id)
      data icall/0/
      save
      icall=icall+1
c$$$      write(6,*)('in dws_all2allvb8 ')
c$$$      write(6,*)('for call no. '),icall
      if(mynprocg.gt.id)then
       write(6,*)('too many procs in dws_allgv!! '),mynprocg,id
       close(unit=6)
c$$$       stop
      end if
c$$$      write(6,*)('in dws_allgv, mynprocg = '),mynprocg
      do i=1,mynprocg
       nblocks(i)=nblock8s(i)
       ioffs(i)=ioff8s(i)
       nblockr(i)=nblock8r(i)
       ioffr(i)=ioff8r(i)
c$$$       write(6,3030)i,nblocks(i),ioffs(i),nblockr(i),ioffr(i)
 3030  format(5i8)
      end do
c$$$      if(icall.ne.10)then
c$$$       call dws_sync
c$$$       call dws_finalize
c$$$       stop
c$$$      end if
c$$$      if(ioffs(3).ne.-132)then
c$$$       call dws_sync
c$$$       call dws_finalize
c$$$       stop
c$$$      end if
      call mpi_alltoallv(bufs,nblocks,ioffs,mpi_double_precision,
     $     bufr,nblockr,ioffr,mpi_double_precision,mpi_comm_world,ierr)
      if(ierr.ne.mpi_success)then
       write(6,*)('in dws_all2allv, mpi_alltoallv returned an error'),
     $      ierror
       close(unit=6)
c$$$       stop
      end if
      return
      end
      subroutine dws_12all(buf,len,isource8)
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      dimension buf(len)
      data icall/0/
      save icall
      icall=icall+1
      icont=len
      isource=isource8
c$$$      write(6,*)('in dws_12all, send '),icont,(' from '),isource
c$$$      write(6,*)('call no. '),icall
      call mpi_bcast(buf,icont,mpi_double_precision,isource,            2d19s10
     $               MPI_COMM_WORLD,ierror)                             2d19s10
      if(ierror.ne.mpi_success)then
       write(6,*)('in dws_12all,')
       write(6,*)('mpi_bcast returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
      close(unit=6)
c$$$       stop
      end if
      return
      end
      subroutine dws_12all_loc(buf,len,isource8)                        11d16s12
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      common/mycomu/my_comml
      dimension buf(len)
      data icall/0/
      save icall
      icall=icall+1
      icont=len
      isource=isource8
c$$$      write(6,*)('in dws_12all, send '),icont,(' from '),isource
c$$$      write(6,*)('call no. '),icall
      call mpi_bcast(buf,icont,mpi_double_precision,isource,            2d19s10
     $               my_comml,ierror)                                   11d16s12
      if(ierror.ne.mpi_success)then
       write(6,*)('in dws_12all_loc,')                                  11d16s12
       write(6,*)('mpi_bcast returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
      close(unit=6)
c$$$       stop
      end if
      return
      end
      subroutine dws_12all_loc_dum(buf,len,isource8)                        11d16s12
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      common/mycomu/my_comml
      dimension buf(len)
      data icall/0/
      save icall
      icall=icall+1
      return
      end
      subroutine dws_allgv2(buf,nblock8,ioff8,send)
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      include "common.mympi"                                            1d29s21
      dimension buf(1),nblock8(1),ioff8(1),send(1)
cccccc      parameter (id=1000)
      dimension nblock(id),ioff(id)
      if(mynprocg.gt.id)then
       write(6,*)('too many procs in dws_allgv!! '),mynprocg,id
       close(unit=6)
c$$$       stop
      end if
c$$$      write(6,*)('in dws_allgv, mynprocg = '),mynprocg
      do i=1,mynprocg
       nblock(i)=nblock8(i)
       ioff(i)=ioff8(i)
c$$$       write(6,*)i,nblock(i),ioff(i)
      end do
      idum=nblock(mynowprog+1)
c$$$      if(idum.ne.-132)then
c$$$       call dws_sync
c$$$       call dws_finalize
c$$$       stop
c$$$      end if
c$$$      write(6,*)('sending '),idum
      call mpi_allgatherv(send,idum,mpi_double_precision,buf,
     $     nblock,ioff,mpi_double_precision,mpi_comm_world,ierr)
      if(ierr.ne.mpi_success)then
       write(6,*)('in dws_allgv, mpi_allgatherv returned an error'),
     $      ierror
       close(unit=6)
c$$$       stop
      end if
      return
      end
      subroutine second(time1)
      use mpi
      implicit real*8 (a-h,o-z)
      data ifirst/0/
      save
      time1=mpi_wtime()                                                 5d7s12
      if(ifirst.eq.0)then
       time0=time1
       time1=0d0
       ifirst=1
      else
       time1=time1-time0
      end if
      return
      end
      subroutine ddi_zero(bc,ibc,nhandle)                               11d15s22
      use mpi                                                           3d8s21
      implicit real*8 (a-h,o-z)
      integer*8 nhandle
      integer (kind=mpi_address_kind) idisp
      include "common.mympi"                                             3d5s21
      include "common.store"                                            3d8s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      itmp=ibcoff                                                       3d8s21
c$$$      write(6,*)('ddi_zero for '),nhandle
c$$$      write(6,*)('mdata: '),(mdata(j,nhandle),j=1,5)
      ibcoff=itmp+mdata(4,nhandle)                                      3d8s21
      call enough('ddi.zero',bc,ibc)                                    1d18s23
      do i=itmp,ibcoff-1                                                3d8s21
       bc(i)=0d0                                                        3d8s21
      end do                                                            3d8s21
      iassert=0                                                         3d8s21
      call mpi_win_lock(mpi_lock_exclusive,mynowprog,iassert,           3d8s21
     $        mdata(5,nhandle),ierror)                                  3d5s21
      idisp=0                                                           3d8s21
      call mpi_put(bc(itmp),mdata(4,nhandle),mpi_double_precision,      3d8s21
     $     mynowprog,idisp,mdata(4,nhandle),mpi_double_precision,       3d8s21
     $     mdata(5,nhandle),ierror)                                     3d8s21
      call mpi_win_unlock(mynowprog,mdata(5,nhandle),ierror)            3d8s21
      ibcoff=itmp                                                       3d8s21
      return                                                            3d8s21
      end
      subroutine ddi_put(bc,ibc,nhandle,i1,i2,i3,i4,buff)               11d15s22
      use mpi                                                           3d5s21
      implicit real*8 (a-h,o-z)
      integer*8 nhandle,i1,i2,i3,i4
      integer (kind=mpi_address_kind) idisp
      logical lflag
      dimension buff(*)
      include "common.mympi"                                             3d5s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
c$$$      write(6,*)('in ddi_put with args '),nhandle,i1,i2,i3,i4
      if(nhandle.gt.ndata)then
       write(6,*)('you are trying to put to distributed array '),
     $      nhandle
       write(6,*)('but we only have '),ndata,('distributed arrays')
       stop 'ddi_put'
      end if
c$$$      write(6,*)('what we have in mdata: '),(mdata(j,nhandle),j=1,5)
c$$$      if(ndata.ne.132)then
c$$$       call dws_sync
c$$$       call dws_finalize
c$$$       stop
c$$$      end if
      if(i1.eq.1.and.i2.eq.mdata(1,nhandle))then                        3d5s21
       ioff=1                                                           3d5s21
c$$$       write(6,*)('let us query size: ')
c$$$       call mpi_win_get_attr(mdata(5,nhandle),MPI_WIN_SIZE,nsizex,lflag,
c$$$     $      ierror)
c$$$      write(6,*)('lflag, size: '),lflag,nsizex
       do ip=0,mynprocg-1                                               3d5s21
        call ilimts(1,mdata(2,nhandle),mynprocg,ip,il,ih,i1s,i1e,i2s,   3d5s21
     $       i2e)                                                       3d5s21
        istart=il                                                       3d8s21
        if(i3.gt.il)istart=i3                                           3d8s21
        iend=ih                                                         3d8s21
        if(i4.le.iend)iend=i4                                           3d8s21
c$$$        write(6,*)('istart,iend: '),istart,iend
        if(iend.ge.istart)then                                          3d8s21
         iassert=0                                                      3d5s21
         call mpi_win_lock(mpi_lock_exclusive,ip,iassert,               3d5s21
     $        mdata(5,nhandle),ierror)                                  3d5s21
c$$$         write(6,*)('for prod '),ip
c$$$         write(6,*)('limits: '),il,ih
         istart=max(il,i3)                                              3d5s21
         iend=min(ih,i4)                                                3d5s21
         nput=mdata(1,nhandle)*(iend+1-istart)                          3d5s21
         idisp=mdata(1,nhandle)*(istart-il)                             3d5s21
c$$$         write(6,*)('istart,iend '),istart,iend
c$$$         write(6,*)('nput '),nput
c$$$         write(6,*)('idisp '),idisp
c$$$         write(6,*)('window? '),mdata(5,nhandle)
c$$$         write(6,*)('in ddi_put,saving to proc '),ip
c$$$         call prntm2(buff(ioff),mdata(1,nhandle),iend+1-istart,
c$$$     $        mdata(1,nhandle))
         call mpi_put(buff(ioff),nput,mpi_double_precision,ip,idisp,    3d5s21
     $        nput,mpi_double_precision,mdata(5,nhandle),ierror)        3d5s21
         call mpi_win_unlock(ip,mdata(5,nhandle),ierror)                3d5s21
         ioff=ioff+nput                                                 3d5s21
        end if                                                          3d5s21
       end do                                                           3d5s21
      else                                                              3d5s21
       write(6,*)('you asked to put subset of rows ... '),i1,i2
       write(6,*)('out of '),mdata(1,nhandle)
       write(6,*)('I have not coded this yet!'),nhandle
       stop 'ddi_put'
      end if                                                            3d5s21
c$$$      if(ndata.ne.132)then
c$$$       call dws_sync
c$$$       call dws_finalize
c$$$       stop
c$$$      end if
      return                                                            3d5s21
      end
      subroutine ddi_iget(bc,ibc,nhandle,i1,i2,i3,i4,buff,iget,nget)    11d15s22
      implicit real*8 (a-h,o-z)
      integer*8 nhandle,i1,i2,i3,i4                                     3d8s21
      include "common.store"                                            11d15s22
      nget=0                                                            3d8s21
      call ddi_get(bc,ibc,nhandle,i1,i2,i3,i4,buff)                     11d15s22
      return                                                            3d8s21
      end
      subroutine ddi_iacc(bc,ibc,nhandle,i1,i2,i3,i4,buff,iacc,nacc)    11d15s22
      implicit real*8 (a-h,o-z)
      integer*8 nhandle,i1,i2,i3,i4                                     3d8s21
      include "common.store"                                            11d15s22
      nacc=0                                                            3d8s21
      call ddi_acc(bc,ibc,nhandle,i1,i2,i3,i4,buff)                     11d15s22
      return
      end
      subroutine ddi_get(bc,ibc,nhandle,i1,i2,i3,i4,buff)               11d15s22
      use mpi                                                           3d5s21
      implicit real*8 (a-h,o-z)
      integer*8 nhandle,i1,i2,i3,i4
      integer (kind=mpi_address_kind) idisp
      dimension buff(*)
      include "common.mympi"                                             3d5s21
      include "common.store"                                            3d8s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      if(nhandle.gt.ndata)then
       write(6,*)('you are trying to get to distributed array '),
     $      nhandle
       write(6,*)('but we only have '),ndata,('distributed arrays')
       stop 'ddi_get'
      end if
      if(i1.eq.1.and.i2.eq.mdata(1,nhandle))then                        3d5s21
       ioff=1                                                           3d5s21
c$$$       write(6,*)('let us query size: ')
c$$$       call mpi_win_get_attr(mdata(5,nhandle),MPI_WIN_SIZE,nsizex,lflag,
c$$$     $      ierror)
c$$$      write(6,*)('lflag, size: '),lflag,nsizex
       do ip=0,mynprocg-1                                               3d5s21
        call ilimts(1,mdata(2,nhandle),mynprocg,ip,il,ih,i1s,i1e,i2s,   3d5s21
     $       i2e)                                                       3d5s21
c$$$        write(6,*)('from proc '),ip,il,i3,ih,i4
        istart=il                                                       3d8s21
        if(i3.gt.il)istart=i3                                           3d8s21
        iend=ih                                                         3d8s21
        if(i4.le.iend)iend=i4                                           3d8s21
c$$$        write(6,*)('istart,iend: '),istart,iend
        if(iend.ge.istart)then                                          3d8s21
         iassert=0                                                      3d5s21
         call mpi_win_lock(mpi_lock_exclusive,ip,iassert,               3d5s21
     $        mdata(5,nhandle),ierror)                                  3d5s21
c$$$         write(6,*)('for prod '),ip
c$$$         write(6,*)('limits: '),il,ih
         istart=max(il,i3)                                              3d5s21
         iend=min(ih,i4)                                                3d5s21
         nput=mdata(1,nhandle)*(iend+1-istart)                          3d5s21
         idisp=mdata(1,nhandle)*(istart-il)                             3d5s21
c$$$         write(6,*)('istart,iend '),istart,iend
c$$$         write(6,*)('nput '),nput
c$$$         write(6,*)('idisp '),idisp
c$$$         write(6,*)('window? '),mdata(5,nhandle)
         call mpi_get(buff(ioff),nput,mpi_double_precision,ip,idisp,    3d5s21
     $        nput,mpi_double_precision,mdata(5,nhandle),ierror)        3d5s21
c$$$         write(6,*)('I''m getting from proc '),ip
c$$$         call prntm2(buff(ioff),mdata(1,nhandle),iend+1-istart,
c$$$     $        mdata(1,nhandle))
         call mpi_win_unlock(ip,mdata(5,nhandle),ierror)                3d5s21
         ioff=ioff+nput                                                 3d5s21
        end if                                                          3d5s21
       end do                                                           3d5s21
      else                                                              3d5s21
c$$$       write(6,*)('going for limited rows '),i1,i2
       mrow=i2+1-i1                                                     3d8s21
       ioff=1                                                           3d5s21
       do ip=0,mynprocg-1                                               3d5s21
        call ilimts(1,mdata(2,nhandle),mynprocg,ip,il,ih,i1s,i1e,i2s,   3d5s21
     $       i2e)                                                       3d5s21
        istart=il                                                       3d8s21
        if(i3.gt.il)istart=i3                                           3d8s21
        iend=ih                                                         3d8s21
        if(i4.le.iend)iend=i4                                           3d8s21
        if(iend.ge.istart)then                                          3d8s21
         iassert=0                                                      3d5s21
         call mpi_win_lock(mpi_lock_exclusive,ip,iassert,               3d5s21
     $        mdata(5,nhandle),ierror)                                  3d5s21
         mcol=iend+1-istart                                             3d8s21
         nput=mdata(1,nhandle)*(iend+1-istart)                          3d5s21
         idisp=mdata(1,nhandle)*(istart-il)                             3d5s21
         itmp=ibcoff                                                    3d8s21
         ibcoff=itmp+nput                                               3d8s21
         call enough('ddi.get',bc,ibc)                                         11d15s22
         call mpi_get(bc(itmp),nput,mpi_double_precision,ip,idisp,      3d8s21
     $        nput,mpi_double_precision,mdata(5,nhandle),ierror)        3d5s21
c$$$         write(6,*)('what we got from proc '),ip
c$$$         call prntm2(bc(itmp),mdata(1,nhandle),mcol,mdata(1,nhandle))
         call mpi_win_unlock(ip,mdata(5,nhandle),ierror)                3d5s21
         do i=istart,iend                                               3d8s21
          im=i-istart                                                   3d8s21
          jtmp=itmp-1+mdata(1,nhandle)*im                               3d8s21
          joff=ioff-i1+mrow*im                                          3d8s21
          do j=i1,i2                                                    3d8s21
           buff(joff+j)=bc(jtmp+j)                                      3d8s21
          end do                                                        3d8s21
         end do                                                         3d8s21
         ioff=ioff+mrow*mcol                                            3d8s21
         ibcoff=itmp                                                    3d8s21
        end if                                                          3d5s21
       end do                                                           3d5s21
       mtot=i4+1-i3                                                     3d8s21
c$$$       write(6,*)('altogether now ')
c$$$       call prntm2(buff,mrow,mtot,mrow)
c$$$       write(6,*)('you asked to get subset of rows ... '),i1,i2
c$$$       write(6,*)('out of '),mdata(1,nhandle)
c$$$       write(6,*)('I have not coded this yet!')
c$$$       write(6,*)('what about columns? '),i3,i4,mdata(2,nhandle)
c$$$       stop 'ddi_get'
      end if                                                            3d5s21
      return                                                            3d5s21
      end
      subroutine ddi_done(nhandle,i1)
      implicit real*8 (a-h,o-z)
      integer*8 nhandle
c$$$      write(6,*)('you have reached ddi_done')
c$$$      call dws_sync
c$$$      call dws_finalize
c$$$      stop
      idum=1                                                            3d8s21
      return                                                            3d8s21
      end
      subroutine ddi_destroy(nhandle)
      use mpi                                                           3d5s21
      implicit real*8 (a-h,o-z)
      integer*8 nhandle
      include "common.mympi"                                             3d5s21
      if(nhandle.ne.ndata)then
       write(6,*)('trying to destroy distributed arrays out of order')
       write(6,*)nhandle,('vs.'),ndata
       stop 'ddi_destroy'
      end if
      call dws_sync
c$$$      write(6,*)('deleting window '),nhandle,mdata(5,nhandle)
      call mpi_win_free(mdata(5,nhandle),ierror)
      ndata=ndata-1
c$$$      if(ndata.ne.-132)then
c$$$       call dws_sync
c$$$       call dws_finalize
c$$$       stop
c$$$      end if
      return
      end
      subroutine ddi_create(bc,ibc,irow,icol,nhandle)                   11d15s22
      use mpi
      implicit real*8 (a-h,o-z)
      integer*8 nhandle,icol,irow
      include "common.mympi"                                             3d5s21
      character*10 value
      logical lflag                                                                        
      integer(kind=mpi_address_kind) nsizeb,iwinbase,nsizex
      data icall/0/                                                     3d8s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      save icall                                                        3d8s21
      icall=icall+1                                                     3d8s21
      ndata=ndata+1                                                     3d5s21
c$$$      write(6,*)('you have reached ddi_create'),ndata,irow,icol
      if(ndata.gt.id)then
       write(6,*)('you are trying to create distributed array no. '),
     $      ndata
       write(6,*)('but maximum dimensions in common.mympi is '),id
       stop 'ddi_create'
      end if                                                            3d5s21
      nhandle=ndata                                                     3d5s21
      mdata(1,ndata)=irow                                               3d5s21
      mdata(2,ndata)=icol                                               3d5s21
      call ilimts(1,icol,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)      3d5s21
      nhere=ih+1-il                                                     3d5s21
      nneed=nhere*irow                                                  3d5s21
c$$$      write(6,*)('ddi_create '),irow,icol,nhere,nneed
      call mpi_info_create(info,ierror)
c$$$      call mpi_info_set(info,'no_locks','true',ierror)
      call mpi_info_set(info,'no_locks','false',ierror)
c$$$      call mpi_info_get(info,'no_locks',10,value,lflag,ierror)
c$$$      write(6,*)('ierror after info_get '),ierror
c$$$      write(6,*)('lflag from info_get '),lflag
c$$$      if(lflag)then
c$$$       write(6,*)('value = "'),value,('"')
c$$$      end if
      nsizeb=nneed*8                                                    3d5s21
      ndisp=8                                                           3d5s21
c$$$      write(6,*)('creating '),nneed,irow,icol,ndata,icall
      iwinbase=0                                                        8d8s24
      call mpi_win_allocate(nsizeb,ndisp,info,mpi_comm_world,iwinbase,  3d5s21
     $     iwin,ierror)                                                 3d5s21
c$$$      write(6,*)('after win_allocate, ierror = '),ierror
c$$$      write(6,*)('after mpi_win_allocate '),ndata,iwinbase,iwin
c$$$      write(6,*)('let us query size: ')
c$$$      call mpi_win_get_attr(iwin,MPI_WIN_SIZE,nsizex,lflag,ierror)
c$$$      write(6,*)('lflag, size: '),lflag,nsizex
      mdata(3,ndata)=iwinbase                                           3d5s21
      mdata(4,ndata)=nneed                                              3d5s21
      mdata(5,ndata)=iwin                                               3d5s21
c$$$      if(icall.gt.1)then
c$$$       call dws_sync
c$$$       call dws_finalize
c$$$       stop
c$$$      end if
      call dws_sync                                                     3d5s21
c$$$      call dws_sync
c$$$      call dws_finalize
c$$$      stop
      return                                                            3d5s21
      end
      subroutine ddi_acc(bc,ibc,nhandle,i1,i2,i3,i4,buff)               11d15s22
      use mpi                                                           3d5s21
      implicit real*8 (a-h,o-z)
      integer*8 nhandle,i1,i2,i3,i4
      integer (kind=mpi_address_kind) idisp
      logical lflag
      dimension buff(*)
      include "common.mympi"                                             3d5s21
      include "common.store"                                            3d8s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      if(nhandle.gt.ndata)then
       write(6,*)('you are trying to put to distributed array '),
     $      nhandle
       write(6,*)('but we only have '),ndata,('distributed arrays')
       stop 'ddi_acc'
      end if
c$$$      write(6,*)('Hi, my name is ddi_acc!')
c$$$      if(bc(132).ne.-132d0)then
c$$$       call dws_synca
c$$$       call dws_finalize
c$$$       stop
c$$$      end if
      if(i1.eq.1.and.i2.eq.mdata(1,nhandle))then                        3d5s21
c$$$       write(6,*)('ddi_acc block 1')
       ioff=1                                                           3d5s21
c$$$       write(6,*)('let us query size: ')
c$$$       call mpi_win_get_attr(mdata(5,nhandle),MPI_WIN_SIZE,nsizex,lflag,
c$$$     $      ierror)
c$$$      write(6,*)('lflag, size: '),lflag,nsizex
       do ip=0,mynprocg-1                                               3d5s21
c$$$        write(6,*)('for ip = '),ip
        call ilimts(1,mdata(2,nhandle),mynprocg,ip,il,ih,i1s,i1e,i2s,   3d5s21
     $       i2e)                                                       3d5s21
        istart=il                                                       3d8s21
        if(i3.gt.il)istart=i3                                           3d8s21
        iend=ih                                                         3d8s21
        if(i4.le.iend)iend=i4                                           3d8s21
c$$$        write(6,*)('istart,iend: '),istart,iend
        if(iend.ge.istart)then                                          3d8s21
         iassert=0                                                      3d5s21
c$$$         write(6,*)('lock '),mdata(5,nhandle),ip,nhandle
c$$$         if(mdata(5,nhandle).ne.-132.and.bc(132).ne.-132d0)then
c$$$          call dws_synca
c$$$          call dws_finalize
c$$$          stop
c$$$         end if
         call mpi_win_lock(mpi_lock_exclusive,ip,iassert,               3d5s21
     $        mdata(5,nhandle),ierror)                                  3d5s21
c$$$         write(6,*)('after lock, ierror = '),ierror
c$$$         if(ierror.ne.-132.and.bc(132).ne.-132d0)then
c$$$          call dws_synca
c$$$          call dws_finalize
c$$$          stop
c$$$         end if
c$$$         write(6,*)('for prod '),ip
c$$$         write(6,*)('limits: '),il,ih
         istart=max(il,i3)                                              3d5s21
         iend=min(ih,i4)                                                3d5s21
         nput=mdata(1,nhandle)*(iend+1-istart)                          3d5s21
         idisp=mdata(1,nhandle)*(istart-il)                             3d5s21
c$$$         write(6,*)('istart,iend '),istart,iend
c$$$         write(6,*)('nput '),nput
c$$$         write(6,*)('idisp '),idisp
c$$$         write(6,*)('window? '),mdata(5,nhandle)
c$$$         write(6,*)('in ddi_put,saving to proc '),ip
c$$$         call prntm2(buff(ioff),mdata(1,nhandle),iend+1-istart,
c$$$     $        mdata(1,nhandle))
c$$$         write(6,*)('accumulate '),nput,idisp
         call mpi_accumulate(buff(ioff),nput,mpi_double_precision,ip,   3d8s21
     $        idisp,nput,mpi_double_precision,mpi_sum,mdata(5,nhandle), 3d8s21
     $        ierror)                                                   3d8s21
c$$$         write(6,*)('unlock ')
         call mpi_win_unlock(ip,mdata(5,nhandle),ierror)                3d5s21
         ioff=ioff+nput                                                 3d5s21
        end if                                                          3d5s21
       end do                                                           3d5s21
      else                                                              3d5s21
c$$$       write(6,*)('ddi_acc block 2')
       mrow=i2+1-i1                                                     3d8s21
       ioff=1                                                           3d5s21
       do ip=0,mynprocg-1                                               3d5s21
c$$$        write(6,*)('for proc '),ip
        call ilimts(1,mdata(2,nhandle),mynprocg,ip,il,ih,i1s,i1e,i2s,   3d5s21
     $       i2e)                                                       3d5s21
        istart=il                                                       3d8s21
        if(i3.gt.il)istart=i3                                           3d8s21
        iend=ih                                                         3d8s21
        if(i4.le.iend)iend=i4                                           3d8s21
c$$$        write(6,*)('end,start: '),iend,istart
        if(iend.ge.istart)then                                          3d8s21
         iassert=0                                                      3d5s21
c$$$         write(6,*)('lock ')
         call mpi_win_lock(mpi_lock_exclusive,ip,iassert,               3d5s21
     $        mdata(5,nhandle),ierror)                                  3d5s21
         nput=mdata(1,nhandle)*(iend+1-istart)                          3d5s21
         idisp=mdata(1,nhandle)*(istart-il)                             3d5s21
         itmp=ibcoff                                                    3d8s21
         ibcoff=itmp+nput                                               3d8s21
         call enough('ddi.acc',bc,ibc)                                         11d15s22
         do i=itmp,ibcoff-1                                             3d8s21
          bc(i)=0d0                                                     3d8s21
         end do                                                         3d8s21
         do i=istart,iend                                               3d8s21
          im=i-istart                                                   3d8s21
          jtmp=itmp-1+mdata(1,nhandle)*im                               3d8s21
          joff=ioff-i1+mrow*im                                          3d8s21
          do j=i1,i2                                                    3d8s21
           bc(jtmp+j)=buff(joff+j)                                      3d8s21
          end do                                                        3d8s21
         end do                                                         3d8s21
         ioff=ioff+mrow*(iend+1-istart)                                 3d8s21
c$$$         write(6,*)('accumulate '),nput,idisp
         call mpi_accumulate(bc(itmp),nput,mpi_double_precision,ip,     3d8s21
     $        idisp,nput,mpi_double_precision,mpi_sum,mdata(5,nhandle), 3d8s21
     $        ierror)                                                   3d8s21
c$$$         write(6,*)('unlock ')
         call mpi_win_unlock(ip,mdata(5,nhandle),ierror)                3d5s21
         ibcoff=itmp                                                    3d8s21
        end if                                                          3d5s21
       end do                                                           3d5s21
      end if                                                            3d5s21
      return
      end                  
      subroutine ddi_iaccword0
      idum=1
      return
      end
      subroutine ddi_iaccword2
      idum=1
      return
      end

