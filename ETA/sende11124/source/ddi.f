c     implement ddi using mpi commands, call by fortran
c     half of procs will be compute procs, the other half will
c     be memory procs, so ddi will work across nodes.
c     0 and even nodes will be compute nodes.
c     itag=1 for proc node to memory node directive
c     other values for memory node to proc node.
c
      subroutine dws_preinit                                            1d31s21
      use mpi
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      call mpi_init(ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_int,')
       write(6,*)('mpi_init returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
c$$$      write(6,*)('back from mpi_init ')
      call mpi_comm_size(MPI_COMM_WORLD,isize,ierror)
      call mpi_comm_rank(MPI_COMM_WORLD,irank,ierror)
      if(irank.eq.0)write(6,*)('Hi, this is preinit for MPI-1 code')    8d29s22
      mynowprog=irank                                                   1d31s21
      mynprocg=isize                                                    1d31s21
      return                                                            1d31s21
      end                                                               1d31s21
c mpec2.1 version zeta copyright u.s. government
      subroutine dws_bcast(buf,len)
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer len
      integer icont
      dimension buf(len)
      include "common.mympi"                                            1d29s21
      data icall/0/
      icall=icall+1
      icont=len
      isource=0
      call mpi_bcast(buf,icont,mpi_double_precision,isource,
     $               my_comm_group,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_bcast,')
       write(6,*)('mpi_bcast returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
c mpec2.1 version zeta copyright u.s. government
      subroutine dws_bcasta(buf,len)                                    1d31s21
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer len
      integer icont
      dimension buf(len)
      icont=len
      isource=0
      call mpi_bcast(buf,icont,mpi_double_precision,isource,
     $               mpi_comm_world,ierror)                             1d31s21
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_bcasta,')
       write(6,*)('mpi_bcast returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
c mpec2.1 version zeta copyright u.s. government
      subroutine dws_init(ncore,nmn,bc,ibc)                             11d15s22
      use mpi
c
c     mxcproc is the maximum number of compute procs
c
      parameter (idm=10)                                                2d21s21
      implicit integer (i-n)
      implicit real*8 (a-h,o-z)
      integer*8 iarg1,iarg2,ddi_np,ddi_me                               11d18s20
      character*1 file(8)
      character*8 fname
      character*255 procname                                            1d29s21
      character*101 line                                                2d21s21
      include "common.mympi"                                            1d29s21
c$$$      dimension npacket(8),istatus(mpi_status_size),
c$$$     $     mdata(4,id),packet(4)
      dimension npacket(8),istatus(mpi_status_size),packet(4)           1d26s23
      dimension minfo(7,idm,2),xinfo(5,idm,2),ninfo(2),nwrt(10,2)       2d21s21
      equivalence (packet,npacket)
      data file/'m','o','u','t','.','_','_','_'/
      data nwrt/7,3,7,7,7,2,3,2,3,7,                                    2d21s21
     $          0,0,1,3,1,0,0,0,0,3/                                    2d21s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      include "common.store"
      save
      data ncall/0/                                                     2d8s21
      call mpi_comm_size(MPI_COMM_WORLD,isize,ierror)
      call mpi_comm_rank(MPI_COMM_WORLD,irank,ierror)
      call mpi_get_processor_name(procname,lenname,ierror)
c
c     if nmn = 0, we will use no memory nodes.                          8d12s22
c     if nmn = 1,                                                       8d12s22
c     I assume we are hyperthreading with ncore
c     physical cores per node. Then we have the mapping
c
c     node 0
c     phsyical core 0          irank 0,ncore
c                   1          irank 1,ncore+1
c        .
c        .
c        .
c     phsyical core ncore-1    irank ncore-1,ncore*2-1
c
c     node 1
c     physical core 0          irank 2*ncore,3*ncore
c
c     physical core ncore-1    irank 3*ncore-1,4*ncore-1
c
c
c     node n
c     physical core m          irank n*2*ncore+m,"+ncore
c
      ncore2=ncore*2
      if(ncore2.gt.isize)then                                           2d8s21
       ncore2=isize                                                     2d8s21
       ncore=ncore2/2                                                   2d8s21
      end if                                                            2d8s21
      nu=0                                                              1d29s21
      ncorem=ncore-1                                                    2d8s21
      do irankx=0,isize-1                                               1d29s21
       node=irankx/ncore2
       n0phys=irankx-node*ncore2
       if(n0phys.lt.ncore.or.nmn.eq.0)then                              8d12s22
        nu=nu+1
        if(nu.gt.mxcproc)stop 'mxcproc'
        iranks(nu)=irankx
       end if
      end do
      node=irank/ncore2
      n0phys=irank-node*ncore2
      if(n0phys.ge.ncore.and.nmn.ne.0)then                              8d12s22
       immem=1
      else
       immem=0
      end if
      mnmc=nmn                                                          8d12s22
      if(nmn.ne.0)then                                                  8d12s22
       ngroupsize=isize/2                                                1d29s21
      else                                                              8d12s22
       ngroupsize=isize                                                 8d12s22
      end if                                                            8d12s22
      if(immem.eq.0)then                                                1d29s21
       mymemp=irank+ncore                                               1d29s21
       mymast=irank                                                     1d29s21
       if(irank.ne.0)then
        write(fname,1066)irank
        do i=1,8
         if(fname(i:i).eq.' ')fname(i:i)=file(i)
        end do
        fname(1:1)='p'
        open(unit=6,file=fname)
       end if
c
c     compute node
c
       call mpi_comm_group(MPI_COMM_WORLD,igroup_comm_world,ierror)
       call mpi_group_incl(igroup_comm_world,ngroupsize,iranks,         1d29s21
     $      igroup_comm_group,ierror)
       do i=1,ngroupsize                                                1d29s21
        iranks(i)=iranks(i)+ncore                                       1d29s21
       end do                                                           1d29s21
       itag=3
       call mpi_comm_create_group(MPI_COMM_WORLD,igroup_comm_group,
     $      itag,my_comm_group,ierror)
       call mpi_comm_size(my_comm_group,iq,ierror)
       call mpi_comm_rank(my_comm_group,iq2,ierror)
       mynowprog=iq2                                                    11d18s20
       mynprocg=iq                                                      11d18s20
       return
      else
       ninfo(1)=0                                                       2d21s21
       ninfo(2)=0                                                       2d21s21
       nuse=1                                                           2d21s21
       write(fname,1066)irank
 1066  format(i8)
       do i=1,8
        if(fname(i:i).eq.' ')fname(i:i)=file(i)
       end do
       open(unit=6,file=fname)
       do i=1,ngroupsize                                                1d29s21
        iranks(i)=iranks(i)+ncore                                       1d29s21
       end do                                                           1d29s21
c
c     memory node
c
       ibcoff=1
       ndata=0
       ipass=0
    1  continue
        ipass=ipass+1
c
c     wait for incoming message
c
       itag=1
 3329  format(10000i1)
c
c     we need to use probe to determine the size of the message.
c     this is because if we do ddi_putx as two calls, one to
c     tell me the operation and size, and the second to transfer
c     the data, another processor might have sent us a message
c     between the two that will screw things up.
c
       call mpi_probe(mpi_any_source,itag,mpi_comm_world,istatus,
     $      ierror)
       call mpi_get_count(istatus,mpi_double_precision,nget,ierror)
       if(nget.le.4)then
c$$$        write(6,*)('mpi_recva '),nget
        call mpi_recv(packet,nget,mpi_double_precision,mpi_any_source,
     $      itag,mpi_comm_world,istatus,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('mpi_recv returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop 'recva'
      end if
       else
        ipck=ibcoff
        ibcoff=ipck+nget
        call enough('ddimem.1',bc,ibc)
c$$$        write(6,*)('mpi_recvb '),nget
        call mpi_recv(bc(ipck),nget,mpi_double_precision,mpi_any_source,
     $      itag,mpi_comm_world,istatus,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('mpi_recv returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop 'recvb'
      end if
        jpck=ipck-1
        do i=1,4
         packet(i)=bc(jpck+i)
        end do
        jpck=ipck+4
       end if
       icode=npacket(1)
       ninfo(nuse)=ninfo(nuse)+1                                        2d21s21
       if(ninfo(nuse).gt.idm)then                                       2d21s21
        ninfo(nuse)=ninfo(nuse)-1                                       2d21s21
        if(nuse.eq.1)then                                               2d21s21
         nuse=2                                                         2d21s21
        else                                                            2d21s21
         nuse=1                                                         2d21s21
        end if                                                          2d21s21
        ninfo(nuse)=1                                                   2d21s21
       end if                                                           2d21s21
       ifrom=npacket(6)                                                 2d21s21
       minfo(1,ninfo(nuse),nuse)=ifrom                                  2d21s21
       minfo(2,ninfo(nuse),nuse)=icode                                  2d21s21
c
c     branch on type of request
c
       if(icode.eq.1)then
c
c     creating a distributed array
c
        ndata=ndata+1
        minfo(3,ninfo(nuse),nuse)=ndata                                 2d21s21
        if(ndata.gt.id)then
         write(6,*)('id exceeded in create call '),ndata,id
         stop
        end if
        do i=1,4
         ip=i+1
         mdata(i,ndata)=npacket(ip)
         minfo(3+i,ninfo(nuse),nuse)=npacket(ip)                        2d21s21
        end do
        if(ibcoff.ne.mdata(3,ndata))then
         write(6,*)('compute and memory node out of sync '),ndata,
     $        ibcoff,mdata(3,ndata)
         do i=1,ndata
          write(6,*)i,(mdata(j,i),j=1,4)
         end do
         stop
        end if
        ibcoff=ibcoff+mdata(4,ndata)
        ipck=ibcoff
        call enough('ddimem.2',bc,ibc)
       else if(icode.eq.2)then
c
c     delete a distributed array
c
        nhandle=npacket(2)
        minfo(3,ninfo(nuse),nuse)=nhandle                               2d21s21
        if(nhandle.ne.ndata)then
         write(6,*)('destroying distributed arrays out of order ')
         write(6,*)('handle: '),nhandle
         write(6,*)('ndata: '),ndata
         stop
        end if
        ibcoff=mdata(3,ndata)
        ipck=ibcoff
        ndata=ndata-1
       else if(icode.eq.3)then
c
c     put part of array
c
        nhandle=npacket(2)
        nrow=npacket(3)
        istrt=npacket(4)
        iend=npacket(5)
        ilj=npacket(8)
        ifrom=npacket(6)
        ilow=npacket(7)-1
        minfo(3,ninfo(nuse),nuse)=nrow                                  2d21s21
        minfo(4,ninfo(nuse),nuse)=istrt                                 2d21s21
        minfo(5,ninfo(nuse),nuse)=iend                                  2d21s21
        minfo(6,ninfo(nuse),nuse)=ilj                                   2d21s21
        minfo(7,ninfo(nuse),nuse)=ilow                                  2d21s21
        xinfo(1,ninfo(nuse),nuse)=bc(ipck+4)                            2d21s21
        jtag=10
        if(nrow.eq.mdata(1,nhandle))then                                3d6s09
         iaddx=mdata(1,nhandle)*(istrt-ilj)+mdata(3,nhandle)
         nmove=nrow*(iend+1-istrt)                                      3d6s09
         jpck=ipck+4
         do j=0,nmove-1
          bc(iaddx+j)=bc(jpck+j)
         end do
        else
         itmp=ipck+4
         nwds=(iend+1-istrt)*nrow
         jtmp=itmp
         do i=istrt,iend
          iaddx=mdata(1,nhandle)*(i-ilj)+mdata(3,nhandle)+ilow
          do j=0,nrow-1
           bc(iaddx+j)=bc(jtmp+j)
          end do
          jtmp=jtmp+nrow
         end do
         ibcoff=ipck
        end if
       else if(icode.eq.4)then
c
c     accumulate into an array
c
        nhandle=npacket(2)
        nrow=npacket(3)
        istrt=npacket(4)
        iend=npacket(5)
        ifrom=npacket(6)
        ilow=npacket(7)-1
        ilj=npacket(8)
        minfo(3,ninfo(nuse),nuse)=nrow                                  2d21s21
        minfo(4,ninfo(nuse),nuse)=istrt                                 2d21s21
        minfo(5,ninfo(nuse),nuse)=iend                                  2d21s21
        minfo(6,ninfo(nuse),nuse)=ilj                                   2d21s21
        minfo(7,ninfo(nuse),nuse)=ilow                                  2d21s21
        xinfo(1,ninfo(nuse),nuse)=bc(ipck+4)                            2d21s21
        jtag=11
        if(nrow.eq.mdata(1,nhandle))then                                3d6s09
         iaddx=mdata(1,nhandle)*(istrt-ilj)+mdata(3,nhandle)
         nmove=nrow*(iend+1-istrt)                                      3d6s09
         itmp=ipck+4
         xinfo(2,ninfo(nuse),nuse)=bc(iaddx)                            2d21s21
         do i=0,nmove-1
          bc(iaddx+i)=bc(iaddx+i)+bc(itmp+i)
         end do
         xinfo(3,ninfo(nuse),nuse)=bc(iaddx)                            2d21s21
         ibcoff=ipck
        else
         nwds=(iend+1-istrt)*nrow
         itmp=ipck+4
         jtmp=itmp
         do i=istrt,iend
          iaddx=mdata(1,nhandle)*(i-ilj)+mdata(3,nhandle)+ilow
          if(i.eq.istrt)xinfo(2,ninfo(nuse),nuse)=bc(iaddx)                            2d21s21
          do j=0,nrow-1
           bc(iaddx+j)=bc(iaddx+j)+bc(jtmp+j)
          end do
          if(i.eq.istrt)xinfo(3,ninfo(nuse),nuse)=bc(iaddx)                            2d21s21
          jtmp=jtmp+nrow
         end do
         ibcoff=ipck
        end if
       else if(icode.eq.5)then
c
c     get part of array
c
        nhandle=npacket(2)
        nrow=npacket(3)
        istrt=npacket(4)
        iend=npacket(5)
        ilj=npacket(8)
        ifrom=npacket(6)
        ilow=npacket(7)-1
        jtag=12
        minfo(3,ninfo(nuse),nuse)=nrow                                  2d21s21
        minfo(4,ninfo(nuse),nuse)=istrt                                 2d21s21
        minfo(5,ninfo(nuse),nuse)=iend                                  2d21s21
        minfo(6,ninfo(nuse),nuse)=ilj                                   2d21s21
        minfo(7,ninfo(nuse),nuse)=ilow                                  2d21s21
        if(nrow.eq.mdata(1,nhandle))then                                3d6s09
         iaddx=mdata(1,nhandle)*(istrt-ilj)+mdata(3,nhandle)
         xinfo(1,ninfo(nuse),nuse)=bc(iaddx)
         nmove=nrow*(iend+1-istrt)                                      3d6s09
         call mpi_send(bc(iaddx),nmove,mpi_double_precision,ifrom,
     $        jtag,mpi_comm_world,ierror)
        else
         itmp=ibcoff
         nwds=(iend+1-istrt)*nrow
         ibcoff=itmp+nwds
         call enough('ddimem.3',bc,ibc)
         jtmp=itmp
         do i=istrt,iend
          iaddx=mdata(1,nhandle)*(i-ilj)+mdata(3,nhandle)+ilow
          if(i.eq.istrt)xinfo(1,ninfo(nuse),nuse)=bc(iaddx)             2d21s21
          do j=0,nrow-1
           bc(jtmp+j)=bc(iaddx+j)
          end do
          jtmp=jtmp+nrow
         end do
         call mpi_send(bc(itmp),nwds,mpi_double_precision,ifrom,
     $         jtag,mpi_comm_world,ierror)
         ibcoff=itmp
        end if
       else if(icode.eq.6)then
        ipck=ibcoff
c
c     initialize load balancing counter
c
        idlb=0
       else if(icode.eq.7)then
        ipck=ibcoff
c
c     give value of balancing counter
c
        ito=npacket(2)
        minfo(3,ninfo(nuse),nuse)=idlb                                  2d21s21
        itag=10
        call mpi_send(idlb,1,mpi_integer,ito,itag,mpi_comm_world,
     $       ierror)
        idlb=idlb+1
       else if(icode.eq.8)then
c
c     perform global sync
c
        call mpi_barrier(mpi_comm_world,ierror)
        ipck=ibcoff
       else if(icode.eq.9)then                                          11d21s20
c
c     zero out my part of array
c
        nhandle=npacket(2)                                              11d21s20
        minfo(3,ninfo(nuse),nuse)=nhandle                               2d21s21
        do i=0,mdata(4,nhandle)-1                                       11d21s20
         bc(mdata(3,nhandle)+i)=0d0                                     11d21s20
        end do                                                          11d21s20
       else if(icode.eq.10)then                                         2d19s21
c
c     accumulate into an array with message saying we are done.         2d19s21
c
        nhandle=npacket(2)
        nrow=npacket(3)
        istrt=npacket(4)
        iend=npacket(5)
        ifrom=npacket(6)
        ilow=npacket(7)-1
        ilj=npacket(8)
        minfo(3,ninfo(nuse),nuse)=nrow                                  2d21s21
        minfo(4,ninfo(nuse),nuse)=istrt                                 2d21s21
        minfo(5,ninfo(nuse),nuse)=iend                                  2d21s21
        minfo(6,ninfo(nuse),nuse)=ilj                                   2d21s21
        minfo(7,ninfo(nuse),nuse)=ilow                                  2d21s21
        xinfo(1,ninfo(nuse),nuse)=bc(ipck+4)                            2d21s21
        jtag=11                                                         8d11s22
        if(nrow.eq.mdata(1,nhandle))then                                3d6s09
         iaddx=mdata(1,nhandle)*(istrt-ilj)+mdata(3,nhandle)
         nmove=nrow*(iend+1-istrt)                                      3d6s09
         itmp=ipck+4
         nwiacc=nwiacc+nmove                                            8d10s22
         xinfo(2,ninfo(nuse),nuse)=bc(iaddx)                            2d21s21
         do i=0,nmove-1
          bc(iaddx+i)=bc(iaddx+i)+bc(itmp+i)
         end do
         xinfo(3,ninfo(nuse),nuse)=bc(iaddx)                            2d21s21
         ibcoff=ipck
        else
         nwds=(iend+1-istrt)*nrow
         nwiacc=nwiacc+nwds                                             8d10s22
         itmp=ipck+4
         jtmp=itmp
         do i=istrt,iend
          iaddx=mdata(1,nhandle)*(i-ilj)+mdata(3,nhandle)+ilow
          if(i.eq.istrt)xinfo(2,ninfo(nuse),nuse)=bc(iaddx)                            2d21s21
          do j=0,nrow-1
           bc(iaddx+j)=bc(iaddx+j)+bc(jtmp+j)
          end do
          if(i.eq.istrt)xinfo(3,ninfo(nuse),nuse)=bc(iaddx)                            2d21s21
          jtmp=jtmp+nrow
         end do
         ibcoff=ipck
        end if
       else if(icode.eq.11)then                                         8d10s22
c
c     zero iacc word counter
c
        nwiacc=0                                                        8d10s22
       else if(icode.eq.12)then                                         8d10s22
c
c     return iacc word counter
c
        ito=npacket(2)
        itag=10
        call mpi_send(nwiacc,1,mpi_integer,ito,itag,mpi_comm_world,
     $       ierror)
       else if(icode.eq.0)then
c
c     close down
c
        call mpi_barrier(mpi_comm_world,ierror)
        call mpi_finalize(ierror)
        stop
       end if
       ibcoff=ipck
       go to 1
      end if
      end
c mpec2.1 version zeta copyright u.s. government
      subroutine dws_finalize
      use mpi
      implicit integer (i-n)
      real*8 packet
      dimension npacket(2)
      equivalence (packet,npacket)
      include "common.mympi"                                            1d29s21
      if(mnmc.ne.0)then                                                 8d12s22
      itag=1
      call mpi_comm_rank(my_comm_group,irank,ierror)                     11d18s20
      ito=mymemp                                                        1d29s21
      npacket(1)=0
      call mpi_send(packet,1,mpi_double_precision,ito,itag,
     $     mpi_comm_world,
     $     ierror)
      end if                                                            8d12s22
      call mpi_barrier(mpi_comm_world,ierror)
      call mpi_finalize(ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_finalize,')
       write(6,*)('mpi_finalize returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
c mpec2.1 version zeta copyright u.s. government
      subroutine dws_synca
      use mpi
      implicit integer (i-n)
      real*8 packet
      dimension npacket(2)
      equivalence (packet,npacket)
      include "common.mympi"                                            1d29s21
      if(mnmc.ne.0)then                                                 8d12s22
      itag=1
      call mpi_comm_rank(my_comm_group,irank,ierror)                    11d18s20
 3329 format(10000i1)
      ito=mymemp                                                        1d29s21
      call mpi_barrier(my_comm_group,ierror)
      npacket(1)=8
      call mpi_send(packet,1,mpi_double_precision,ito,itag,
     $     mpi_comm_world,
     $     ierror)
      end if                                                            8d12s22
      call mpi_barrier(mpi_comm_world,ierror)
      return
      end
c mpec2.1 version zeta copyright u.s. government
      subroutine dws_sync
      use mpi
      implicit integer (i-n)
      real*8 packet
      dimension npacket(2)
      equivalence (packet,npacket)
      include "common.mympi"                                            1d29s21
      call mpi_barrier(my_comm_group,ierror)
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
c mpec2.1 version zeta copyright u.s. government
      subroutine dws_gsumf(buff,nwds)
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      dimension buff(nwds)                                                 2d22s19
      include "common.mympi"                                            1d29s21
      if(nwds.le.0)return                                               10d28s22
      call mpi_allreduce(mpi_in_place,buff,nwds,mpi_double_precision,   2d22s10
     $     mpi_sum,my_comm_group,ierror)                                2d22s10
      if(ierror.ne.mpi_success)then
       write(6,*)('in dws_gsumf,')
       write(6,*)('mpi_allreduce returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
c mpec2.1 version zeta copyright u.s. government
      subroutine dws_gbor(buff,nwds)                                    1d26s21
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer buff(nwds)                                                1d26s21
      include "common.mympi"                                            1d29s21
      nwds2=nwds*2                                                      1d26s21
      call mpi_allreduce(mpi_in_place,buff,nwds2,mpi_integer,           1d26s21
     $     mpi_bor,my_comm_group,ierror)                                1d26s21
      if(ierror.ne.mpi_success)then
       write(6,*)('in dws_gbor,')
       write(6,*)('mpi_allreduce returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      return
      end
c mpec2.1 version zeta copyright u.s. government
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
      end if
      do i=1,mynprocg
       nblock(i)=nblock8(i)
       ioff(i)=ioff8(i)
      end do
      idum=nblock(mynowprog+1)
      call mpi_allgatherv(send,idum,mpi_double_precision,buf,
     $     nblock,ioff,mpi_double_precision,my_comm_group,ierr)
      if(ierr.ne.mpi_success)then
       write(6,*)('in dws_allgv, mpi_allgatherv returned an error'),
     $      ierror
       close(unit=6)
      end if
      return
      end
c mpec2.1 version zeta copyright u.s. government
      subroutine dws_all2allvb(bufs,nblock8s,ioff8s,bufr,nblock8r,
     $     ioff8r)
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      include "common.mympi"                                            1d29s21
      dimension bufs(1),nblock8s(1),ioff8s(1),bufr(1),nblock8r(1),
     $     ioff8r(1)
cccccc      parameter (id=1000)
      dimension nblocks(id),ioffs(id),nblockr(id),ioffr(id)
      data icall/0/
      save
      icall=icall+1
      if(mynprocg.gt.id)then
       write(6,*)('too many procs in dws_allgv!! '),mynprocg,id
       close(unit=6)
      end if
      do i=1,mynprocg
       nblocks(i)=nblock8s(i)
       ioffs(i)=ioff8s(i)
       nblockr(i)=nblock8r(i)
       ioffr(i)=ioff8r(i)
 3030  format(5i8)
      end do
      call mpi_alltoallv(bufs,nblocks,ioffs,mpi_double_precision,
     $     bufr,nblockr,ioffr,mpi_double_precision,my_comm_group,ierr)
      if(ierr.ne.mpi_success)then
       write(6,*)('in dws_all2allv, mpi_alltoallv returned an error'),
     $      ierror
       close(unit=6)
      end if
      return
      end
c mpec2.1 version zeta copyright u.s. government
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
      include "common.mympi"                                            1d29s21
ccccc      parameter (id=1000)
      dimension nblocks(id),ioffs(id),nblockr(id),ioffr(id)
      data icall/0/
      save
      icall=icall+1
      if(mynprocg.gt.id)then
       write(6,*)('too many procs in dws_allgv!! '),mynprocg,id
       close(unit=6)
      end if
      do i=1,mynprocg
       nblocks(i)=nblock8s(i)
       ioffs(i)=ioff8s(i)
       nblockr(i)=nblock8r(i)
       ioffr(i)=ioff8r(i)
 3030  format(5i8)
      end do
      call mpi_alltoallv(bufs,nblocks,ioffs,mpi_double_precision,
     $     bufr,nblockr,ioffr,mpi_double_precision,my_comm_group,ierr)
      if(ierr.ne.mpi_success)then
       write(6,*)('in dws_all2allv, mpi_alltoallv returned an error'),
     $      ierror
       close(unit=6)
      end if
      return
      end
c mpec2.1 version zeta copyright u.s. government
      subroutine dws_12all(buf,len,isource8)
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      dimension buf(len)
      include "common.mympi"                                            1d29s21
      data icall/0/
      save icall
      icall=icall+1
      icont=len
      isource=isource8
      call mpi_bcast(buf,icont,mpi_double_precision,isource,            2d19s10
     $               my_comm_group,ierror)                              11d18s20
      if(ierror.ne.mpi_success)then
       write(6,*)('in dws_12all,')
       write(6,*)('mpi_bcast returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
      close(unit=6)
      end if
      return
      end
c mpec2.1 version zeta copyright u.s. government
      subroutine ddi_destroy(nhandlex)                                  11d15s22
      integer*8 nhandle,nhandlex
      nhandle=nhandlex+1
      call ddi_destroyx(nhandle)                                        11d15s22
      return
      end
      subroutine ddi_destroyx(nhandle)
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer*8 irow,icol,nhandle(1),ilow1,ihigh1,ilow2,ihigh2,nhandlex,1d18s21
     $     ione,iarg3,npno,iarg1,iarg2,idelta(1)                        1d18s21
      logical no_locks
      integer (kind=mpi_address_kind) memreq,iadd,iaddx,itrial
      integer*4 iwin,nneed,ierror,info
c$$$      include "common.store"
      include "common.mympi"                                            1d29s21
cccccccc      parameter (id=2000)
c$$$      dimension mdata(4,id),istatus(mpi_status_size),npacket(8),
c$$$     $     packet(4)
      dimension istatus(mpi_status_size),npacket(8),packet(4)
      equivalence(packet,npacket)
   51 format('ddi_destroy: ',i5)
      if(nhandle(1).ne.ndata)then                                       1d18s21
       write(6,*)('destroying distributed arrays out of order ')
       write(6,*)('handle: '),nhandle(1)                                1d18s21
       write(6,*)('ndata: '),ndata
       stop
      end if
      iwinbase=mdata(3,ndata)
      call mpi_comm_rank(my_comm_group,irank,ierror)
      ito=mymemp                                                        1d29s21
      itag=1
      npacket(1)=2
      npacket(2)=ndata
      call mpi_send(packet,1,mpi_double_precision,ito,itag,
     $     mpi_comm_world,ierror)
      ndata=ndata-1
      call dws_sync                                                     3d13s09
      return
      end
      subroutine ddi_zero(bc,ibc,nhandlex)                              7d18s24
      implicit real*8 (a-h,o-z)                                         7d18s24
      implicit integer*8 (i-n)                                          7d18s24
      integer my_comm_group,mymast,mymemp,iranks,mnmc,mdata,ndata,      1d26s23
     $ iwinbase                                                         1d26s23
      include "common.mympi"                                            1d29s21
      nhandle=nhandlex+1                                                11d19s20
      call ddi_zerox(bc,ibc,nhandle)                                    11d15s22
      return                                                            7d18s24
      end                                                               7d18s24
      subroutine ddi_put(bc,ibc,nhandlex,ilow1,ihigh1,ilow2,ihigh2,     7d18s24
     $     buff)                                                        7d18s24
      implicit real*8 (a-h,o-z)                                         7d18s24
      implicit integer*8 (i-n)                                          7d18s24
      integer my_comm_group,mymast,mymemp,iranks,mnmc,mdata,ndata,      1d26s23
     $ iwinbase                                                         1d26s23
      include "common.mympi"                                            1d29s21
      nhandle=nhandlex+1
      call ddi_putx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff)      11d15s22
      return
      end                                                               7d18s24
      subroutine ddi_iacc(bc,ibc,nhandlex,ilow1,ihigh1,ilow2,ihigh2,    7d18s24
     $   buff,iacc,nacc)                                                11d15s22
      implicit real*8 (a-h,o-z)                                         7d18s24
      implicit integer*8 (i-n)                                          7d18s24
      integer my_comm_group,mymast,mymemp,iranks,mnmc,mdata,ndata,      1d26s23
     $ iwinbase                                                         1d26s23
      include "common.mympi"                                            1d29s21
      nhandle=nhandlex+1                                                2d15s12
      call ddi_iaccx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff,iacc,11d15s22
     $   nacc)                                                          11d15s22
      return                                                            2d15s12
      end                                                               7d18s24
      subroutine ddi_acc(bc,ibc,nhandlex,ilow1,ihigh1,ilow2,ihigh2,buff)7d18s24
      implicit real*8 (a-h,o-z)
      implicit integer*8 (i-n)
      integer my_comm_group,mymast,mymemp,iranks,mnmc,mdata,ndata,      1d26s23                  6d17s21
     $ iwinbase                                                         1d26s23
      include "common.mympi"                                            1d29s21
      nhandle=nhandlex+1                                                2d15s12
      call ddi_accx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff)      11d15s22
      return                                                            2d15s12
      end                                                               7d18s24
      subroutine ddi_create(bc,ibc,irow,icol,nhandlex)                  11d15s22
      implicit real*8 (a-h,o-z)
      implicit integer*8 (i-n)
      integer my_comm_group,mymast,mymemp,iranks,mnmc,mdata,ndata,      1d26s23                  6d17s21
     $ iwinbase                                                         1d26s23
      include "common.mympi"                                            1d29s21
      call ddi_createx(bc,ibc,irow,icol,nhandle)                        11d15s22
      nhandlex=nhandle-1
      return
      end
      subroutine ddi_getx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff)7d18s24
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer*8 irow,icol,nhandle(1),ilow1,ihigh1,ilow2,ihigh2,nhandlex,1d18s21
     $     ione,iarg3,npno,iarg1,iarg2,idelta(1)                        1d18s21
      logical no_locks
      integer (kind=mpi_address_kind) memreq,iadd,iaddx,itrial
      integer*4 iwin,nneed,ierror,info
      include "common.store"
      include "common.mympi"                                            1d29s21
      dimension istatus(mpi_status_size),npacket(8),packet(4)
      equivalence(packet,npacket)
      dimension buff(1),ircv(1)                                         1d29s21
      data loopx/50/                                                    6d17s21
      data ndata,ncall,loop/3*0/                                        6d17s21
      save
      call mpi_comm_rank(my_comm_group,irank,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_create,')
       write(6,*)('mpi_comm_rank returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      call mpi_comm_size(my_comm_group,isize,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_create,')
       write(6,*)('mpi_size_rank returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      if(nhandle(1).gt.ndata.or.nhandle(1).lt.1)then                    1d18s21
       write(6,*)('in ddi_get, handle = '),nhandle(1),nhandlex          1d18s21
       write(6,*)('exceeds ndata = '),ndata
       stop
      end if
      npacket(2)=nhandle(1)                                             1d18s21
      if(min(ilow1,ilow2).le.0)then
       write(6,*)('in ddi_get, bad lower limits '),ilow1,ilow2
       stop
      end if
      if(ihigh1.gt.mdata(1,nhandle(1)).or.                              1d18s21
     $     ihigh2.gt.mdata(2,nhandle(1)))then                           1d18s21
       write(6,*)('in ddi_get, bad upper limits '),ihigh1,ihigh2
       write(6,*)('vs '),mdata(1,nhandle(1)),mdata(2,nhandle(1))        1d18s21
       write(6,*)('nhandle = '),nhandle(1)                              1d18s21
       stop
      end if
      nrow=ihigh1+1-ilow1
      npacket(3)=nrow
      npacket(6)=mymast                                                 1d29s21
      npacket(7)=ilow1
      npacket(1)=5
      itag=1
      jtag=12
      do iproc=0,isize-1
       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
     $      ilj,ihj,i1s,i1e,i2s,i2e)
       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
        istrt=max(ilow2,ilj)
        iend=min(ihigh2,ihj)
        npacket(4)=istrt
        npacket(5)=iend
        npacket(8)=ilj
        ito=iranks(iproc+1)                                             1d29s21
        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
         i1=1+nrow*(istrt-ilow2)                                        3d6s09
         iaddx=mdata(1,nhandle(1))*(istrt-ilj)                          1d18s21
     $        +ilow1-1+mdata(3,nhandle(1))                              1d18s21
         nmove=nrow*(iend+1-istrt)                                      3d6s09
         if(nmove.gt.0)then
         call mpi_send(packet,4,mpi_double_precision,ito,itag,
     $         mpi_comm_world,ierror)
c$$$        write(6,*)('mpi_recvc '),nmove
         call mpi_recv(buff(i1),nmove,mpi_double_precision,ito,jtag,
     $        mpi_comm_world,istatus,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('mpi_recv returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop 'recvc'
      end if
         end if
        else                                                            3d6s09
         nwds=(iend+1-istrt)*nrow
         if(nwds.gt.0)then
         itmp=ibcoff
         ibcoff=itmp+nwds
         call enough('ddi_createx.  5',bc,ibc)
         call mpi_send(packet,4,mpi_double_precision,ito,itag,
     $         mpi_comm_world,
     $       ierror)
c$$$        write(6,*)('mpi_recvd '),nwds
         call mpi_recv(bc(itmp),nwds,mpi_double_precision,ito,jtag,
     $        mpi_comm_world,istatus,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('mpi_recv returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop 'recvd'
      end if
         jtmp=itmp
         do i=istrt,iend
          i1=1+nrow*(i-ilow2)
          do j=0,nrow-1
           buff(i1+j)=bc(jtmp+j)
          end do
          jtmp=jtmp+nrow
         end do
         ibcoff=itmp
         end if
       end if
       end if
      end do
      return
      end                                                               7d18s24
      subroutine ddi_iaccx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,    7d18s24
     $     buff,ircv,nrcv)                                              7d18s24
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer*8 irow,icol,nhandle(1),ilow1,ihigh1,ilow2,ihigh2,nhandlex,1d18s21
     $     ione,iarg3,npno,iarg1,iarg2,idelta(1)                        1d18s21
      logical no_locks
      integer (kind=mpi_address_kind) memreq,iadd,iaddx,itrial
      integer*4 iwin,nneed,ierror,info
      include "common.store"
      include "common.mympi"                                            1d29s21
      dimension istatus(mpi_status_size),npacket(8),packet(4)
      equivalence(packet,npacket)
      dimension buff(1),ircv(1)                                         1d29s21
      data loopx/50/                                                    6d17s21
      data ndata,ncall,loop/3*0/                                        6d17s21
      save
c
c     important note!!!
c     buff needs to have 4*mynprocg words of extra storage for this to
c     work properly, yet there is no way to test for this !!!
c
      call mpi_comm_rank(my_comm_group,irank,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_create,')
       write(6,*)('mpi_comm_rank returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      call mpi_comm_size(my_comm_group,isize,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_create,')
       write(6,*)('mpi_size_rank returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      if(nhand.gt.ndata.or.nhandle(1).lt.1)then
       write(6,*)('in ddi_accx, handle = '),nhandle(1),nhandlex,nhand
       write(6,*)('exceeds ndata = '),ndata
       stop
      end if
      npacket(2)=nhandle(1)                                             1d18s21
      if(min(ilow1,ilow2).le.0)then
       write(6,*)('in ddi_put, bad lower limits '),ilow1,ilow2
       stop
      end if
      if(ihigh1.gt.mdata(1,nhandle(1)).or.                              1d18s21
     $     ihigh2.gt.mdata(2,nhandle(1)))then                           1d18s21
       write(6,*)('in ddi_accx, bad upper limits '),ihigh1,ihigh2,
     $      mdata(1,nhandle(1)),mdata(2,nhandle(1)),nhandle(1)          1d18s21
       stop
      end if
      nrow=ihigh1+1-ilow1
      npacket(3)=nrow
      npacket(6)=mymast                                                 1d29s21
      npacket(7)=ilow1
      npacket(1)=10
      itag=1
      jtag=11
      nrcv=0                                                            1d30s21
      itmp0=ibcoff                                                      1d30s21
      jtmp0=itmp0                                                       1d30s21
      do iproc=0,isize-1
       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
     $      ilj,ihj,i1s,i1e,i2s,i2e)
       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
        istrt=max(ilow2,ilj)
        iend=min(ihigh2,ihj)
        npacket(4)=istrt
        npacket(5)=iend
        npacket(8)=ilj
        ito=iranks(iproc+1)                                             1d29s21
        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
         i1=1+nrow*(istrt-ilow2)                                        3d6s09
         iaddx=mdata(1,nhandle(1))*(istrt-ilj)                          1d18s21
     $        +ilow1-1+mdata(3,nhandle(1))                              1d18s21
         nmove=nrow*(iend+1-istrt)                                      3d6s09
         if(nmove.gt.0)then
          itmp=jtmp0                                                    1d30s21
          jtmp0=itmp+nmove+4                                            1d30s21
          ibcoff=jtmp0                                                   1d30s21
          call enough('ddi_createx.  3',bc,ibc)
          jtmp=itmp-1
          do j=1,4
           bc(jtmp+j)=packet(j)
          end do
          jtmp=jtmp+5
          do j=0,nmove-1
           bc(jtmp+j)=buff(i1+j)
          end do
         end if
        else                                                            3d6s09
         nwds=(iend+1-istrt)*nrow
         if(nwds.gt.0)then
         itmp=jtmp0                                                     1d30s21
         jtmp0=itmp+nwds+4                                              1d30s21
         ibcoff=jtmp0                                                   1d30s21
         call enough('ddi_createx.  4',bc,ibc)
         jtmp=itmp-1
         do j=1,4
          bc(jtmp+j)=packet(j)
         end do
         jtmp=jtmp+5
         do i=istrt,iend
          i1=1+nrow*(i-ilow2)
          do j=0,nrow-1
           bc(jtmp+j)=buff(i1+j)
          end do
          jtmp=jtmp+nrow
         end do
         nwdsp=nwds+4
         end if
        end if                                                          3d6s09
       end if
      end do
      ntotal=jtmp0-itmp0                                                1d30s21
      jtmp0=itmp0-1                                                     1d30s21
      do i=1,ntotal                                                     1d30s21
       buff(i)=bc(jtmp0+i)                                              1d30s21
      end do                                                            1d30s21
      ibcoff=itmp0                                                      1d30s21
      i1=1
      myall=0
      do iproc=0,isize-1
       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
     $      ilj,ihj,i1s,i1e,i2s,i2e)
       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
        istrt=max(ilow2,ilj)
        iend=min(ihigh2,ihj)
        npacket(4)=istrt
        npacket(5)=iend
        npacket(8)=ilj
        ito=iranks(iproc+1)                                             1d29s21
        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
         nmove=nrow*(iend+1-istrt)                                      3d6s09
         if(nmove.gt.0)then
          nmovep=nmove+4
          nrcv=nrcv+1                                                   1d30s21
          myall=myall+nmove
          call mpi_isend(buff(i1),nmovep,mpi_double_precision,ito,itag, 6d17s21
     $        mpi_comm_world,ircv(nrcv),ierror)                         1d30s21
          i1=i1+nmovep                                                  1d30s21
         end if
        else                                                            3d6s09
         nwds=(iend+1-istrt)*nrow
         if(nwds.gt.0)then
          nwdsp=nwds+4
          nrcv=nrcv+1                                                   1d30s21
          myall=myall+nmove
          call mpi_isend(buff(i1),nwdsp,mpi_double_precision,ito,itag,  6d17s21
     $        mpi_comm_world,ircv(nrcv),ierror)                         1d30s21
          i1=i1+nwdsp                                                   1d30s21
         end if
        end if                                                          3d6s09
       end if
      end do
      return
      end                                                               7d18s24
      subroutine ddi_accx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff)7d18s24
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer*8 irow,icol,nhandle(1),ilow1,ihigh1,ilow2,ihigh2,nhandlex,1d18s21
     $     ione,iarg3,npno,iarg1,iarg2,idelta(1)                        1d18s21
      logical no_locks
      integer (kind=mpi_address_kind) memreq,iadd,iaddx,itrial
      integer*4 iwin,nneed,ierror,info
      include "common.store"
      include "common.mympi"                                            1d29s21
      dimension istatus(mpi_status_size),npacket(8),packet(4)
      equivalence(packet,npacket)
      dimension buff(1),ircv(1)                                         1d29s21
      data loopx/50/                                                    6d17s21
      data ndata,ncall,loop/3*0/                                        6d17s21
      save
      call mpi_comm_rank(my_comm_group,irank,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_create,')
       write(6,*)('mpi_comm_rank returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      call mpi_comm_size(my_comm_group,isize,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_create,')
       write(6,*)('mpi_size_rank returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      if(nhand.gt.ndata.or.nhandle(1).lt.1)then
       write(6,*)('in ddi_accx, handle = '),nhandle(1),nhandlex,nhand
       write(6,*)('exceeds ndata = '),ndata
       stop
      end if
      npacket(2)=nhandle(1)                                             1d18s21
      if(min(ilow1,ilow2).le.0)then
       write(6,*)('in ddi_put, bad lower limits '),ilow1,ilow2
       stop
      end if
      if(ihigh1.gt.mdata(1,nhandle(1)).or.                              1d18s21
     $     ihigh2.gt.mdata(2,nhandle(1)))then                           1d18s21
       write(6,*)('in ddi_accx, bad upper limits '),ihigh1,ihigh2,
     $      mdata(1,nhandle(1)),mdata(2,nhandle(1)),nhandle(1)          1d18s21
       stop
      end if
      nrow=ihigh1+1-ilow1
      npacket(3)=nrow
      npacket(6)=mymast                                                 1d29s21
      npacket(7)=ilow1
      npacket(1)=4                                                      2d21s21
      itag=1
      jtag=11
      ktag=13                                                           2d19s21
      do iproc=0,isize-1
       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
     $      ilj,ihj,i1s,i1e,i2s,i2e)
       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
        istrt=max(ilow2,ilj)
        iend=min(ihigh2,ihj)
        npacket(4)=istrt
        npacket(5)=iend
        npacket(8)=ilj
        ito=iranks(iproc+1)                                             1d29s21
        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
         i1=1+nrow*(istrt-ilow2)                                        3d6s09
         iaddx=mdata(1,nhandle(1))*(istrt-ilj)                          1d18s21
     $        +ilow1-1+mdata(3,nhandle(1))                              1d18s21
         nmove=nrow*(iend+1-istrt)                                      3d6s09
         if(nmove.gt.0)then
          itmp=ibcoff
          ibcoff=itmp+nmove+4
          jtmp=itmp-1
          do j=1,4
           bc(jtmp+j)=packet(j)
          end do
          jtmp=jtmp+5
          do j=0,nmove-1
           bc(jtmp+j)=buff(i1+j)
          end do
          nmovep=nmove+4
         call mpi_ssend(bc(itmp),nmovep,mpi_double_precision,ito,itag,  2d5s21
     $        mpi_comm_world,ierror)
          ibcoff=itmp
         end if
        else                                                            3d6s09
         nwds=(iend+1-istrt)*nrow
         if(nwds.gt.0)then
         itmp=ibcoff
         ibcoff=itmp+nwds+4
         call enough('ddi_createx.  2',bc,ibc)
         jtmp=itmp-1
         do j=1,4
          bc(jtmp+j)=packet(j)
         end do
         jtmp=jtmp+5
         do i=istrt,iend
          i1=1+nrow*(i-ilow2)
          do j=0,nrow-1
           bc(jtmp+j)=buff(i1+j)
          end do
          jtmp=jtmp+nrow
         end do
         nwdsp=nwds+4
         call mpi_ssend(bc(itmp),nwdsp,mpi_double_precision,ito,itag,   2d5s21
     $        mpi_comm_world,ierror)
         ibcoff=itmp
         end if
        end if                                                          3d6s09
       end if
      end do
      return
      end                                                               7d18s24
      subroutine ddi_putx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff)7d18s24
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer*8 irow,icol,nhandle(1),ilow1,ihigh1,ilow2,ihigh2,nhandlex,1d18s21
     $     ione,iarg3,npno,iarg1,iarg2,idelta(1)                        1d18s21
      logical no_locks
      integer (kind=mpi_address_kind) memreq,iadd,iaddx,itrial
      integer*4 iwin,nneed,ierror,info
      include "common.store"
      include "common.mympi"                                            1d29s21
      dimension istatus(mpi_status_size),npacket(8),packet(4)
      equivalence(packet,npacket)
      dimension buff(1),ircv(1)                                         1d29s21
      data loopx/50/                                                    6d17s21
      data ndata,ncall,loop/3*0/                                        6d17s21
      save
      call mpi_comm_rank(my_comm_group,irank,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_create,')
       write(6,*)('mpi_comm_rank returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      call mpi_comm_size(my_comm_group,isize,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_create,')
       write(6,*)('mpi_size_rank returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      if(nhand.gt.ndata.or.nhandle(1).lt.1)then                         1d18s21
       write(6,*)('in ddi_put, handle = '),nhandle(1),nhandlex,nhand    1d18s21
       write(6,*)('exceeds ndata = '),ndata
       stop
      end if
      npacket(1)=3
      npacket(2)=nhandle(1)                                             1d18s21
      if(min(ilow1,ilow2).le.0)then
       write(6,*)('in ddi_put, bad lower limits '),ilow1,ilow2
       stop
      end if
      if(ihigh1.gt.mdata(1,nhandle(1)).or                               1d18s21
     $     .ihigh2.gt.mdata(2,nhandle(1)))then                          1d18s21
       write(6,*)('in ddi_put, bad upper limits '),ihigh1,ihigh2,
     $      mdata(1,nhandle(1)),mdata(2,nhandle(1)),nhandle(1)          1d18s21
       stop
      end if
      nrow=ihigh1+1-ilow1
      npacket(3)=nrow
      npacket(5)=mymast                                                 1d29s21
      npacket(7)=ilow1
      npacket(8)=0
      itag=1
      jtag=10
      do iproc=0,isize-1
       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
     $      ilj,ihj,i1s,i1e,i2s,i2e)
       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
        istrt=max(ilow2,ilj)
        iend=min(ihigh2,ihj)
        npacket(4)=istrt
        npacket(5)=iend
        npacket(8)=ilj
        ito=iranks(iproc+1)                                             1d29s21
        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
         i1=1+nrow*(istrt-ilow2)                                        3d6s09
         iaddx=mdata(1,nhandle(1))*(istrt-ilj)+ilow1-1                  1d18s21
     $        +mdata(3,nhandle(1))                                      1d18s21
         nmove=nrow*(iend+1-istrt)                                      3d6s09
38835    format(4i8,1pe15.7)
         itmp=ibcoff
         ibcoff=itmp+4+nmove
         jtmp=itmp-1
         do j=1,4
          bc(jtmp+j)=packet(j)
         end do
         jtmp=jtmp+5
         do j=0,nmove-1
          bc(jtmp+j)=buff(i1+j)
         end do
         nmovep=nmove+4
         call mpi_send(bc(itmp),nmovep,mpi_double_precision,ito,itag,
     $        mpi_comm_world,ierror)
         ibcoff=itmp
        else                                                            3d6s09
         nwds=(iend+1-istrt)*nrow
         if(nwds.gt.0)then
         itmp=ibcoff
         ibcoff=itmp+nwds+4
         call enough('ddi_createx.  1',bc,ibc)
         jtmp=itmp-1
         do j=1,4
          bc(jtmp+j)=packet(j)
         end do
         jtmp=jtmp+5
         do i=istrt,iend
          i1=1+nrow*(i-ilow2)
          do j=0,nrow-1
           bc(jtmp+j)=buff(i1+j)
          end do
          jtmp=jtmp+nrow
c          iaddx=mdata(1,nhandle)*(i-ilj)+ilow1-1+mdata(3,nhandle)        3d6s09
   33    format(5i8)
         end do
         nwdsp=nwds+4
         call mpi_send(bc(itmp),nwdsp,mpi_double_precision,ito,
     $        itag,mpi_comm_world,ierror)
         ibcoff=itmp
         end if
        end if
       end if                                                           3d6s09
      end do
      return
      end
      subroutine ddi_igetx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,    7d18s24
     $   buff,ircv,nrcv)                                                11d15s22
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer*8 irow,icol,nhandle(1),ilow1,ihigh1,ilow2,ihigh2,nhandlex,1d18s21
     $     ione,iarg3,npno,iarg1,iarg2,idelta(1)                        1d18s21
      logical no_locks
      integer (kind=mpi_address_kind) memreq,iadd,iaddx,itrial
      integer*4 iwin,nneed,ierror,info
      include "common.store"
      include "common.mympi"                                            1d29s21
      dimension istatus(mpi_status_size),npacket(8),packet(4)
      equivalence(packet,npacket)
      dimension buff(1),ircv(1)                                         1d29s21
      data loopx/50/                                                    6d17s21
      data ndata,ncall,loop/3*0/                                        6d17s21
      save
      call mpi_comm_rank(my_comm_group,irank,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_create,')
       write(6,*)('mpi_comm_rank returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      call mpi_comm_size(my_comm_group,isize,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_create,')
       write(6,*)('mpi_size_rank returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      if(nhandle(1).gt.ndata.or.nhandle(1).lt.1)then                    1d18s21
       write(6,*)('in ddi_igetx, handle = '),nhandle(1),nhandlex          1d18s21
       write(6,*)('exceeds ndata = '),ndata
       stop
      end if
      npacket(2)=nhandle(1)                                             1d18s21
      if(min(ilow1,ilow2).le.0)then
       write(6,*)('in ddi_igetx, bad lower limits '),ilow1,ilow2
       stop
      end if
      if(ihigh1.gt.mdata(1,nhandle(1)).or.                              1d18s21
     $     ihigh2.gt.mdata(2,nhandle(1)))then                           1d18s21
       write(6,*)('in ddi_iget, bad upper limits '),ihigh1,ihigh2
       write(6,*)('vs '),mdata(1,nhandle(1)),mdata(2,nhandle(1))        1d18s21
       write(6,*)('nhandle = '),nhandle(1)                              1d18s21
       stop
      end if
      nrow=ihigh1+1-ilow1
      npacket(3)=nrow
      npacket(6)=mymast                                                 1d29s21
      npacket(7)=ilow1
      npacket(1)=5
      itag=1
      jtag=12
      nrcv=0                                                            1d29s21
      do iproc=0,isize-1
       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
     $      ilj,ihj,i1s,i1e,i2s,i2e)
       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
        istrt=max(ilow2,ilj)
        iend=min(ihigh2,ihj)
        npacket(4)=istrt
        npacket(5)=iend
        npacket(8)=ilj
        ito=iranks(iproc+1)                                             1d29s21
        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
         i1=1+nrow*(istrt-ilow2)                                        3d6s09
         iaddx=mdata(1,nhandle(1))*(istrt-ilj)                          1d18s21
     $        +ilow1-1+mdata(3,nhandle(1))                              1d18s21
         nmove=nrow*(iend+1-istrt)                                      3d6s09
         if(nmove.gt.0)then
         call mpi_send(packet,4,mpi_double_precision,ito,itag,
     $         mpi_comm_world,ierror)
         nrcv=nrcv+1                                                    1d29s21
         call mpi_irecv(buff(i1),nmove,mpi_double_precision,ito,jtag,   1d29s21
     $        mpi_comm_world,ircv(nrcv),ierror)
         end if
        else                                                            3d6s09
         nwds=(iend+1-istrt)*nrow
         if(nwds.gt.0)then
         itmp=ibcoff
         ibcoff=itmp+nwds
         call enough('ddi_createx.  6',bc,ibc)
         call mpi_send(packet,4,mpi_double_precision,ito,itag,
     $         mpi_comm_world,
     $       ierror)
c$$$        write(6,*)('mpi_recve '),nwds
         call mpi_recv(bc(itmp),nwds,mpi_double_precision,ito,jtag,     1d29s21
     $        mpi_comm_world,istatus,ierror)                            1d29s21
      if(ierror.ne.mpi_success)then
       write(6,*)('mpi_recv returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop 'recve'
      end if
         jtmp=itmp
         do i=istrt,iend
          i1=1+nrow*(i-ilow2)
          do j=0,nrow-1
           buff(i1+j)=bc(jtmp+j)
          end do
          jtmp=jtmp+nrow
         end do
         ibcoff=itmp
         end if
       end if
       end if
      end do
      return
      end                                                               7d18s24
      subroutine ddi_zerox(bc,ibc,nhandle)                              7d18s24
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer*8 irow,icol,nhandle(1),ilow1,ihigh1,ilow2,ihigh2,nhandlex,1d18s21
     $     ione,iarg3,npno,iarg1,iarg2,idelta(1)                        1d18s21
      logical no_locks
      integer (kind=mpi_address_kind) memreq,iadd,iaddx,itrial
      integer*4 iwin,nneed,ierror,info
      include "common.store"
      include "common.mympi"                                            1d29s21
      dimension istatus(mpi_status_size),npacket(8),packet(4)
      equivalence(packet,npacket)
      dimension buff(1),ircv(1)                                         1d29s21
      data loopx/50/                                                    6d17s21
      data ndata,ncall,loop/3*0/                                        6d17s21
      save
      ito=mymemp                                                        1d29s21
      itag=1
      npacket(1)=9                                                      11d21s20
      npacket(2)=nhandle(1)                                             1d18s21
      call mpi_send(packet,2,mpi_double_precision,ito,itag,             11d21s20
     $     mpi_comm_world,ierror)
      return                                                            11d19s20
      end                                                               7d18s24
c mpec2.1 version zeta copyright u.s. government
      subroutine ddi_createx(bc,ibc,irow,icol,nhandle)                  11d15s22
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer*8 irow,icol,nhandle(1),ilow1,ihigh1,ilow2,ihigh2,nhandlex,1d18s21
     $     ione,iarg3,npno,iarg1,iarg2,idelta(1)                        1d18s21
      logical no_locks
      integer (kind=mpi_address_kind) memreq,iadd,iaddx,itrial
      integer*4 iwin,nneed,ierror,info
      include "common.store"
      include "common.mympi"                                            1d29s21
cccccccc      parameter (id=2000)
c$$$      dimension mdata(4,id),istatus(mpi_status_size),npacket(8),
c$$$     $     packet(4)
      dimension istatus(mpi_status_size),npacket(8),packet(4)
      equivalence(packet,npacket)
      dimension buff(1),ircv(1)                                         1d29s21
      data loopx/50/                                                    6d17s21
      data ndata,ncall,loop/3*0/                                        6d17s21
      save
      if(ndata.eq.0)then
       nrun=0
       iwinbase=1
      end if
      ndata=ndata+1
      if(ndata.gt.id)then
       write(6,*)('tried to create distributed array '),ndata
       write(6,*)('but in ddi_create, id = '),id
       stop
      end if
      nhandle(1)=ndata                                                  1d18s21
      mdata(1,ndata)=irow
      mdata(2,ndata)=icol
      call mpi_comm_rank(my_comm_group,irank,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_create,')
       write(6,*)('mpi_comm_rank returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      call mpi_comm_size(my_comm_group,isize,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('in ddi_create,')
       write(6,*)('mpi_size_rank returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop
      end if
      n1=1
      n2=icol
      call ilimts(n1,n2,isize,0,ilj,ihj,i1s,i1e,i2s,i2e)
      nneed=irow*(ihj+1-ilj)                                            3d3s09
      mdata(3,ndata)=iwinbase
      iwinbase=iwinbase+nneed                                           3d6s09
      mdata(4,ndata)=nneed
      ito=mymemp                                                        1d29s21
      itag=1
      npacket(1)=itag
      do i=1,4
       ip=i+1
       npacket(ip)=mdata(i,ndata)
      end do
      call mpi_send(packet,4,mpi_double_precision,ito,itag,
     $     mpi_comm_world,ierror)
 3356 format('ddi_create: ',i5,5x,3i5,i8)
      return
c$$$      entry ddi_zerox(bc,ibc,nhandle)                                   11d15s22
c$$$      ito=mymemp                                                        1d29s21
c$$$      itag=1
c$$$      npacket(1)=9                                                      11d21s20
c$$$      npacket(2)=nhandle(1)                                             1d18s21
c$$$      call mpi_send(packet,2,mpi_double_precision,ito,itag,             11d21s20
c$$$     $     mpi_comm_world,ierror)
c$$$      return                                                            11d19s20
c$$$      entry ddi_putx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff)     11d15s22
c$$$      if(nhand.gt.ndata.or.nhandle(1).lt.1)then                         1d18s21
c$$$       write(6,*)('in ddi_put, handle = '),nhandle(1),nhandlex,nhand    1d18s21
c$$$       write(6,*)('exceeds ndata = '),ndata
c$$$       stop
c$$$      end if
c$$$      npacket(1)=3
c$$$      npacket(2)=nhandle(1)                                             1d18s21
c$$$      if(min(ilow1,ilow2).le.0)then
c$$$       write(6,*)('in ddi_put, bad lower limits '),ilow1,ilow2
c$$$       stop
c$$$      end if
c$$$      if(ihigh1.gt.mdata(1,nhandle(1)).or                               1d18s21
c$$$     $     .ihigh2.gt.mdata(2,nhandle(1)))then                          1d18s21
c$$$       write(6,*)('in ddi_put, bad upper limits '),ihigh1,ihigh2,
c$$$     $      mdata(1,nhandle(1)),mdata(2,nhandle(1)),nhandle(1)          1d18s21
c$$$       stop
c$$$      end if
c$$$      nrow=ihigh1+1-ilow1
c$$$      npacket(3)=nrow
c$$$      npacket(5)=mymast                                                 1d29s21
c$$$      npacket(7)=ilow1
c$$$      npacket(8)=0
c$$$      itag=1
c$$$      jtag=10
c$$$      do iproc=0,isize-1
c$$$       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
c$$$     $      ilj,ihj,i1s,i1e,i2s,i2e)
c$$$       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
c$$$        istrt=max(ilow2,ilj)
c$$$        iend=min(ihigh2,ihj)
c$$$        npacket(4)=istrt
c$$$        npacket(5)=iend
c$$$        npacket(8)=ilj
c$$$        ito=iranks(iproc+1)                                             1d29s21
c$$$        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
c$$$         i1=1+nrow*(istrt-ilow2)                                        3d6s09
c$$$         iaddx=mdata(1,nhandle(1))*(istrt-ilj)+ilow1-1                  1d18s21
c$$$     $        +mdata(3,nhandle(1))                                      1d18s21
c$$$         nmove=nrow*(iend+1-istrt)                                      3d6s09
c$$$38835    format(4i8,1pe15.7)
c$$$         itmp=ibcoff
c$$$         ibcoff=itmp+4+nmove
c$$$         jtmp=itmp-1
c$$$         do j=1,4
c$$$          bc(jtmp+j)=packet(j)
c$$$         end do
c$$$         jtmp=jtmp+5
c$$$         do j=0,nmove-1
c$$$          bc(jtmp+j)=buff(i1+j)
c$$$         end do
c$$$         nmovep=nmove+4
c$$$         call mpi_send(bc(itmp),nmovep,mpi_double_precision,ito,itag,
c$$$     $        mpi_comm_world,ierror)
c$$$         ibcoff=itmp
c$$$        else                                                            3d6s09
c$$$         nwds=(iend+1-istrt)*nrow
c$$$         if(nwds.gt.0)then
c$$$         itmp=ibcoff
c$$$         ibcoff=itmp+nwds+4
c$$$         call enough('ddi_createx.  1',bc,ibc)
c$$$         jtmp=itmp-1
c$$$         do j=1,4
c$$$          bc(jtmp+j)=packet(j)
c$$$         end do
c$$$         jtmp=jtmp+5
c$$$         do i=istrt,iend
c$$$          i1=1+nrow*(i-ilow2)
c$$$          do j=0,nrow-1
c$$$           bc(jtmp+j)=buff(i1+j)
c$$$          end do
c$$$          jtmp=jtmp+nrow
c$$$c          iaddx=mdata(1,nhandle)*(i-ilj)+ilow1-1+mdata(3,nhandle)        3d6s09
c$$$   33    format(5i8)
c$$$         end do
c$$$         nwdsp=nwds+4
c$$$         call mpi_send(bc(itmp),nwdsp,mpi_double_precision,ito,
c$$$     $        itag,mpi_comm_world,ierror)
c$$$         ibcoff=itmp
c$$$         end if
c$$$        end if
c$$$       end if                                                           3d6s09
c$$$      end do
c$$$      return
c$$$      entry ddi_accx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff)     11d15s22
c$$$      if(nhand.gt.ndata.or.nhandle(1).lt.1)then
c$$$       write(6,*)('in ddi_accx, handle = '),nhandle(1),nhandlex,nhand
c$$$       write(6,*)('exceeds ndata = '),ndata
c$$$       stop
c$$$      end if
c$$$      npacket(2)=nhandle(1)                                             1d18s21
c$$$      if(min(ilow1,ilow2).le.0)then
c$$$       write(6,*)('in ddi_put, bad lower limits '),ilow1,ilow2
c$$$       stop
c$$$      end if
c$$$      if(ihigh1.gt.mdata(1,nhandle(1)).or.                              1d18s21
c$$$     $     ihigh2.gt.mdata(2,nhandle(1)))then                           1d18s21
c$$$       write(6,*)('in ddi_accx, bad upper limits '),ihigh1,ihigh2,
c$$$     $      mdata(1,nhandle(1)),mdata(2,nhandle(1)),nhandle(1)          1d18s21
c$$$       stop
c$$$      end if
c$$$      nrow=ihigh1+1-ilow1
c$$$      npacket(3)=nrow
c$$$      npacket(6)=mymast                                                 1d29s21
c$$$      npacket(7)=ilow1
c$$$      npacket(1)=4                                                      2d21s21
c$$$      itag=1
c$$$      jtag=11
c$$$      ktag=13                                                           2d19s21
c$$$      do iproc=0,isize-1
c$$$       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
c$$$     $      ilj,ihj,i1s,i1e,i2s,i2e)
c$$$       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
c$$$        istrt=max(ilow2,ilj)
c$$$        iend=min(ihigh2,ihj)
c$$$        npacket(4)=istrt
c$$$        npacket(5)=iend
c$$$        npacket(8)=ilj
c$$$        ito=iranks(iproc+1)                                             1d29s21
c$$$        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
c$$$         i1=1+nrow*(istrt-ilow2)                                        3d6s09
c$$$         iaddx=mdata(1,nhandle(1))*(istrt-ilj)                          1d18s21
c$$$     $        +ilow1-1+mdata(3,nhandle(1))                              1d18s21
c$$$         nmove=nrow*(iend+1-istrt)                                      3d6s09
c$$$         if(nmove.gt.0)then
c$$$          itmp=ibcoff
c$$$          ibcoff=itmp+nmove+4
c$$$          jtmp=itmp-1
c$$$          do j=1,4
c$$$           bc(jtmp+j)=packet(j)
c$$$          end do
c$$$          jtmp=jtmp+5
c$$$          do j=0,nmove-1
c$$$           bc(jtmp+j)=buff(i1+j)
c$$$          end do
c$$$          nmovep=nmove+4
c$$$         call mpi_ssend(bc(itmp),nmovep,mpi_double_precision,ito,itag,  2d5s21
c$$$     $        mpi_comm_world,ierror)
c$$$          ibcoff=itmp
c$$$         end if
c$$$        else                                                            3d6s09
c$$$         nwds=(iend+1-istrt)*nrow
c$$$         if(nwds.gt.0)then
c$$$         itmp=ibcoff
c$$$         ibcoff=itmp+nwds+4
c$$$         call enough('ddi_createx.  2',bc,ibc)
c$$$         jtmp=itmp-1
c$$$         do j=1,4
c$$$          bc(jtmp+j)=packet(j)
c$$$         end do
c$$$         jtmp=jtmp+5
c$$$         do i=istrt,iend
c$$$          i1=1+nrow*(i-ilow2)
c$$$          do j=0,nrow-1
c$$$           bc(jtmp+j)=buff(i1+j)
c$$$          end do
c$$$          jtmp=jtmp+nrow
c$$$         end do
c$$$         nwdsp=nwds+4
c$$$         call mpi_ssend(bc(itmp),nwdsp,mpi_double_precision,ito,itag,   2d5s21
c$$$     $        mpi_comm_world,ierror)
c$$$         ibcoff=itmp
c$$$         end if
c$$$        end if                                                          3d6s09
c$$$       end if
c$$$      end do
c$$$      return
c$$$      entry ddi_iaccx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff,
c$$$     $  ircv,nrcv)                                                      11d15s22
c$$$c
c$$$c     important note!!!
c$$$c     buff needs to have 4*mynprocg words of extra storage for this to
c$$$c     work properly, yet there is no way to test for this !!!
c$$$c
c$$$      if(nhand.gt.ndata.or.nhandle(1).lt.1)then
c$$$       write(6,*)('in ddi_accx, handle = '),nhandle(1),nhandlex,nhand
c$$$       write(6,*)('exceeds ndata = '),ndata
c$$$       stop
c$$$      end if
c$$$      npacket(2)=nhandle(1)                                             1d18s21
c$$$      if(min(ilow1,ilow2).le.0)then
c$$$       write(6,*)('in ddi_put, bad lower limits '),ilow1,ilow2
c$$$       stop
c$$$      end if
c$$$      if(ihigh1.gt.mdata(1,nhandle(1)).or.                              1d18s21
c$$$     $     ihigh2.gt.mdata(2,nhandle(1)))then                           1d18s21
c$$$       write(6,*)('in ddi_accx, bad upper limits '),ihigh1,ihigh2,
c$$$     $      mdata(1,nhandle(1)),mdata(2,nhandle(1)),nhandle(1)          1d18s21
c$$$       stop
c$$$      end if
c$$$      nrow=ihigh1+1-ilow1
c$$$      npacket(3)=nrow
c$$$      npacket(6)=mymast                                                 1d29s21
c$$$      npacket(7)=ilow1
c$$$      npacket(1)=10
c$$$      itag=1
c$$$      jtag=11
c$$$      nrcv=0                                                            1d30s21
c$$$      itmp0=ibcoff                                                      1d30s21
c$$$      jtmp0=itmp0                                                       1d30s21
c$$$      do iproc=0,isize-1
c$$$       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
c$$$     $      ilj,ihj,i1s,i1e,i2s,i2e)
c$$$       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
c$$$        istrt=max(ilow2,ilj)
c$$$        iend=min(ihigh2,ihj)
c$$$        npacket(4)=istrt
c$$$        npacket(5)=iend
c$$$        npacket(8)=ilj
c$$$        ito=iranks(iproc+1)                                             1d29s21
c$$$        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
c$$$         i1=1+nrow*(istrt-ilow2)                                        3d6s09
c$$$         iaddx=mdata(1,nhandle(1))*(istrt-ilj)                          1d18s21
c$$$     $        +ilow1-1+mdata(3,nhandle(1))                              1d18s21
c$$$         nmove=nrow*(iend+1-istrt)                                      3d6s09
c$$$         if(nmove.gt.0)then
c$$$          itmp=jtmp0                                                    1d30s21
c$$$          jtmp0=itmp+nmove+4                                            1d30s21
c$$$          ibcoff=jtmp0                                                   1d30s21
c$$$          call enough('ddi_createx.  3',bc,ibc)
c$$$          jtmp=itmp-1
c$$$          do j=1,4
c$$$           bc(jtmp+j)=packet(j)
c$$$          end do
c$$$          jtmp=jtmp+5
c$$$          do j=0,nmove-1
c$$$           bc(jtmp+j)=buff(i1+j)
c$$$          end do
c$$$         end if
c$$$        else                                                            3d6s09
c$$$         nwds=(iend+1-istrt)*nrow
c$$$         if(nwds.gt.0)then
c$$$         itmp=jtmp0                                                     1d30s21
c$$$         jtmp0=itmp+nwds+4                                              1d30s21
c$$$         ibcoff=jtmp0                                                   1d30s21
c$$$         call enough('ddi_createx.  4',bc,ibc)
c$$$         jtmp=itmp-1
c$$$         do j=1,4
c$$$          bc(jtmp+j)=packet(j)
c$$$         end do
c$$$         jtmp=jtmp+5
c$$$         do i=istrt,iend
c$$$          i1=1+nrow*(i-ilow2)
c$$$          do j=0,nrow-1
c$$$           bc(jtmp+j)=buff(i1+j)
c$$$          end do
c$$$          jtmp=jtmp+nrow
c$$$         end do
c$$$         nwdsp=nwds+4
c$$$         end if
c$$$        end if                                                          3d6s09
c$$$       end if
c$$$      end do
c$$$      ntotal=jtmp0-itmp0                                                1d30s21
c$$$      jtmp0=itmp0-1                                                     1d30s21
c$$$      do i=1,ntotal                                                     1d30s21
c$$$       buff(i)=bc(jtmp0+i)                                              1d30s21
c$$$      end do                                                            1d30s21
c$$$      ibcoff=itmp0                                                      1d30s21
c$$$      i1=1
c$$$      myall=0
c$$$      do iproc=0,isize-1
c$$$       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
c$$$     $      ilj,ihj,i1s,i1e,i2s,i2e)
c$$$       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
c$$$        istrt=max(ilow2,ilj)
c$$$        iend=min(ihigh2,ihj)
c$$$        npacket(4)=istrt
c$$$        npacket(5)=iend
c$$$        npacket(8)=ilj
c$$$        ito=iranks(iproc+1)                                             1d29s21
c$$$        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
c$$$         nmove=nrow*(iend+1-istrt)                                      3d6s09
c$$$         if(nmove.gt.0)then
c$$$          nmovep=nmove+4
c$$$          nrcv=nrcv+1                                                   1d30s21
c$$$          myall=myall+nmove
c$$$          call mpi_isend(buff(i1),nmovep,mpi_double_precision,ito,itag, 6d17s21
c$$$     $        mpi_comm_world,ircv(nrcv),ierror)                         1d30s21
c$$$          i1=i1+nmovep                                                  1d30s21
c$$$         end if
c$$$        else                                                            3d6s09
c$$$         nwds=(iend+1-istrt)*nrow
c$$$         if(nwds.gt.0)then
c$$$          nwdsp=nwds+4
c$$$          nrcv=nrcv+1                                                   1d30s21
c$$$          myall=myall+nmove
c$$$          call mpi_isend(buff(i1),nwdsp,mpi_double_precision,ito,itag,  6d17s21
c$$$     $        mpi_comm_world,ircv(nrcv),ierror)                         1d30s21
c$$$          i1=i1+nwdsp                                                   1d30s21
c$$$         end if
c$$$        end if                                                          3d6s09
c$$$       end if
c$$$      end do
c$$$      return
c$$$      entry ddi_getx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff)     11d15s22
c$$$      if(nhandle(1).gt.ndata.or.nhandle(1).lt.1)then                    1d18s21
c$$$       write(6,*)('in ddi_get, handle = '),nhandle(1),nhandlex          1d18s21
c$$$       write(6,*)('exceeds ndata = '),ndata
c$$$       stop
c$$$      end if
c$$$      npacket(2)=nhandle(1)                                             1d18s21
c$$$      if(min(ilow1,ilow2).le.0)then
c$$$       write(6,*)('in ddi_get, bad lower limits '),ilow1,ilow2
c$$$       stop
c$$$      end if
c$$$      if(ihigh1.gt.mdata(1,nhandle(1)).or.                              1d18s21
c$$$     $     ihigh2.gt.mdata(2,nhandle(1)))then                           1d18s21
c$$$       write(6,*)('in ddi_get, bad upper limits '),ihigh1,ihigh2
c$$$       write(6,*)('vs '),mdata(1,nhandle(1)),mdata(2,nhandle(1))        1d18s21
c$$$       write(6,*)('nhandle = '),nhandle(1)                              1d18s21
c$$$       stop
c$$$      end if
c$$$      nrow=ihigh1+1-ilow1
c$$$      npacket(3)=nrow
c$$$      npacket(6)=mymast                                                 1d29s21
c$$$      npacket(7)=ilow1
c$$$      npacket(1)=5
c$$$      itag=1
c$$$      jtag=12
c$$$      do iproc=0,isize-1
c$$$       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
c$$$     $      ilj,ihj,i1s,i1e,i2s,i2e)
c$$$       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
c$$$        istrt=max(ilow2,ilj)
c$$$        iend=min(ihigh2,ihj)
c$$$        npacket(4)=istrt
c$$$        npacket(5)=iend
c$$$        npacket(8)=ilj
c$$$        ito=iranks(iproc+1)                                             1d29s21
c$$$        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
c$$$         i1=1+nrow*(istrt-ilow2)                                        3d6s09
c$$$         iaddx=mdata(1,nhandle(1))*(istrt-ilj)                          1d18s21
c$$$     $        +ilow1-1+mdata(3,nhandle(1))                              1d18s21
c$$$         nmove=nrow*(iend+1-istrt)                                      3d6s09
c$$$         if(nmove.gt.0)then
c$$$         call mpi_send(packet,4,mpi_double_precision,ito,itag,
c$$$     $         mpi_comm_world,ierror)
c$$$         call mpi_recv(buff(i1),nmove,mpi_double_precision,ito,jtag,
c$$$     $        mpi_comm_world,istatus,ierror)
c$$$         end if
c$$$        else                                                            3d6s09
c$$$         nwds=(iend+1-istrt)*nrow
c$$$         if(nwds.gt.0)then
c$$$         itmp=ibcoff
c$$$         ibcoff=itmp+nwds
c$$$         call enough('ddi_createx.  5',bc,ibc)
c$$$         call mpi_send(packet,4,mpi_double_precision,ito,itag,
c$$$     $         mpi_comm_world,
c$$$     $       ierror)
c$$$         call mpi_recv(bc(itmp),nwds,mpi_double_precision,ito,jtag,
c$$$     $        mpi_comm_world,istatus,ierror)
c$$$         jtmp=itmp
c$$$         do i=istrt,iend
c$$$          i1=1+nrow*(i-ilow2)
c$$$          do j=0,nrow-1
c$$$           buff(i1+j)=bc(jtmp+j)
c$$$          end do
c$$$          jtmp=jtmp+nrow
c$$$         end do
c$$$         ibcoff=itmp
c$$$         end if
c$$$       end if
c$$$       end if
c$$$      end do
c$$$      return
c$$$      entry ddi_igetx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,         11d15s22
c$$$     $   buff,ircv,nrcv)                                                11d15s22
c$$$      if(nhandle(1).gt.ndata.or.nhandle(1).lt.1)then                    1d18s21
c$$$       write(6,*)('in ddi_igetx, handle = '),nhandle(1),nhandlex          1d18s21
c$$$       write(6,*)('exceeds ndata = '),ndata
c$$$       stop
c$$$      end if
c$$$      npacket(2)=nhandle(1)                                             1d18s21
c$$$      if(min(ilow1,ilow2).le.0)then
c$$$       write(6,*)('in ddi_igetx, bad lower limits '),ilow1,ilow2
c$$$       stop
c$$$      end if
c$$$      if(ihigh1.gt.mdata(1,nhandle(1)).or.                              1d18s21
c$$$     $     ihigh2.gt.mdata(2,nhandle(1)))then                           1d18s21
c$$$       write(6,*)('in ddi_iget, bad upper limits '),ihigh1,ihigh2
c$$$       write(6,*)('vs '),mdata(1,nhandle(1)),mdata(2,nhandle(1))        1d18s21
c$$$       write(6,*)('nhandle = '),nhandle(1)                              1d18s21
c$$$       stop
c$$$      end if
c$$$      nrow=ihigh1+1-ilow1
c$$$      npacket(3)=nrow
c$$$      npacket(6)=mymast                                                 1d29s21
c$$$      npacket(7)=ilow1
c$$$      npacket(1)=5
c$$$      itag=1
c$$$      jtag=12
c$$$      nrcv=0                                                            1d29s21
c$$$      do iproc=0,isize-1
c$$$       call ilimts(1,mdata(2,nhandle(1)),isize,iproc,                   1d18s21
c$$$     $      ilj,ihj,i1s,i1e,i2s,i2e)
c$$$       if(max(ilow2,ilj).le.min(ihigh2,ihj))then                        3d6s09
c$$$        istrt=max(ilow2,ilj)
c$$$        iend=min(ihigh2,ihj)
c$$$        npacket(4)=istrt
c$$$        npacket(5)=iend
c$$$        npacket(8)=ilj
c$$$        ito=iranks(iproc+1)                                             1d29s21
c$$$        if(nrow.eq.mdata(1,nhandle(1)))then                             1d18s21
c$$$         i1=1+nrow*(istrt-ilow2)                                        3d6s09
c$$$         iaddx=mdata(1,nhandle(1))*(istrt-ilj)                          1d18s21
c$$$     $        +ilow1-1+mdata(3,nhandle(1))                              1d18s21
c$$$         nmove=nrow*(iend+1-istrt)                                      3d6s09
c$$$         if(nmove.gt.0)then
c$$$         call mpi_send(packet,4,mpi_double_precision,ito,itag,
c$$$     $         mpi_comm_world,ierror)
c$$$         nrcv=nrcv+1                                                    1d29s21
c$$$         call mpi_irecv(buff(i1),nmove,mpi_double_precision,ito,jtag,   1d29s21
c$$$     $        mpi_comm_world,ircv(nrcv),ierror)
c$$$         end if
c$$$        else                                                            3d6s09
c$$$         nwds=(iend+1-istrt)*nrow
c$$$         if(nwds.gt.0)then
c$$$         itmp=ibcoff
c$$$         ibcoff=itmp+nwds
c$$$         call enough('ddi_createx.  6',bc,ibc)
c$$$         call mpi_send(packet,4,mpi_double_precision,ito,itag,
c$$$     $         mpi_comm_world,
c$$$     $       ierror)
c$$$         call mpi_recv(bc(itmp),nwds,mpi_double_precision,ito,jtag,     1d29s21
c$$$     $        mpi_comm_world,istatus,ierror)                            1d29s21
c$$$         jtmp=itmp
c$$$         do i=istrt,iend
c$$$          i1=1+nrow*(i-ilow2)
c$$$          do j=0,nrow-1
c$$$           buff(i1+j)=bc(jtmp+j)
c$$$          end do
c$$$          jtmp=jtmp+nrow
c$$$         end do
c$$$         ibcoff=itmp
c$$$         end if
c$$$       end if
c$$$       end if
c$$$      end do
c$$$      return
      end
c mpec2.1 version zeta copyright u.s. government
      subroutine ddi_iget(bc,ibc,nhandlex,ilow1,ihigh1,ilow2,ihigh2,    11d15s22
     $   buff,ircv,nrcv)                                                11d15s22
c
c     non-blocking version
c
      implicit real*8 (a-h,o-z)                                         1d18s21
      implicit integer*8 (i-n)                                          1d18s21
      dimension buff(*),ircv(*)                                         1d29s21
      nhandle=nhandlex+1                                                1d18s21
      call ddi_igetx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff,ircv,11d15s22
     $   nrcv)                                                          11d15s22
      return                                                            1d18s21
      end                                                               1d18s21
c mpec2.1 version zeta copyright u.s. government
      subroutine ddi_iaccword1(iarg)                                    7d18s24
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer*8 iarg,icol,nhandle,istart                                2d4s15
      integer (kind=mpi_address_kind) memreq,iadd,iaddx,ntrial          3d27s09
      dimension istatus(mpi_status_size),npacket(2)
      include "common.mympi"                                            1d29s21
      data nhandle/-132/
      equivalence (packet,npacket)
      save
      itag=1
      call mpi_comm_rank(my_comm_group,irank,ierror)
      ito=mymemp                                                        8d10s22
      npacket(1)=12                                                     8d10s22
      npacket(2)=mymast                                                 1d29s21
      call mpi_send(packet,1,mpi_double_precision,ito,itag,
     $     mpi_comm_world,ierror)
      itag=10
c$$$        write(6,*)('mpi_recvf ')
      call mpi_recv(icount,1,mpi_integer,ito,itag,mpi_comm_world,
     $     istatus,ierror)
      if(ierror.ne.mpi_success)then
       write(6,*)('mpi_recv returned an error ')
       write(6,*)('ierror '),ierror
       write(6,*)('mpi_success '),mpi_success
       stop 'recvf'
      end if
      iarg=icount
      return
      end                                                               7d18s24
      subroutine ddi_iaccword0                                          8d10s22
      use mpi
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      integer*8 iarg,icol,nhandle,istart                                2d4s15
      integer (kind=mpi_address_kind) memreq,iadd,iaddx,ntrial          3d27s09
      dimension istatus(mpi_status_size),npacket(2)
      include "common.mympi"                                            1d29s21
      data nhandle/-132/
      equivalence (packet,npacket)
      save
      call mpi_comm_rank(my_comm_group,irank,ierror)
      call mpi_comm_size(my_comm_group,isize,ierror)
      itag=1
      ito=mymemp
      idum=11                                                           8d10s22
      npacket(1)=idum
      call mpi_send(packet,1,mpi_double_precision,ito,itag,
     $      mpi_comm_world,ierror)
      call dws_sync
      return
c$$$      entry ddi_iaccword1(iarg)
c$$$      itag=1
c$$$      call mpi_comm_rank(my_comm_group,irank,ierror)
c$$$      ito=mymemp                                                        8d10s22
c$$$      npacket(1)=12                                                     8d10s22
c$$$      npacket(2)=mymast                                                 1d29s21
c$$$      call mpi_send(packet,1,mpi_double_precision,ito,itag,
c$$$     $     mpi_comm_world,ierror)
c$$$      itag=10
c$$$      call mpi_recv(icount,1,mpi_integer,ito,itag,mpi_comm_world,
c$$$     $     istatus,ierror)
c$$$      iarg=icount
c$$$      return
      end
c mpec2.1 version zeta copyright u.s. government
      subroutine ddi_iaccword2(nwds)                                    8d10s22
      implicit real*8 (a-h,o-z)                                         8d10s22
      integer*8 iwds                                                    8d10s22
      dimension snd(2)                                                  8d10s22
      snd(1)=dfloat(nwds)                                               8d10s22
      call ddi_iaccword1(iwds)                                          8d10s22
      snd(2)=dfloat(iwds)                                               8d10s22
      call dws_gsumf(snd,2)                                             8d10s22
      nwant=nint(snd(1))                                                8d10s22
      ngot=nint(snd(2))                                                 8d10s22
      if(nwant.ne.ngot)then                                             8d10s22
       looper=0                                                         8d10s22
    1  continue                                                         8d10s22
       looper=looper+1                                                  8d10s22
       if(looper.gt.10)then                                             8d10s22
        write(6,*)('we waited toooo long for iacc to complete! ')       6d22s23
        write(6,*)('nwant: '),nwant                                     6d22s23
        write(6,*)('ngot : '),ngot                                      6d22s23
        call dws_synca
        call dws_finalize
        stop 'ddi_iaccword2'                                            8d10s22
       end if                                                           8d10s22
       call sleep(1)                                                    8d10s22
       call ddi_iaccword1(iwds)                                          8d10s22
       snd(2)=dfloat(iwds)                                               8d10s22
       call dws_gsumf(snd(2),1)                                         8d10s22
       ngot=nint(snd(2))                                                8d10s22
c
c     keep giving reprieves as long as we are still getting more data
c       
       if(looper.eq.1)then                                              6d23s23
        nlast=ngot                                                      6d23s23
       else                                                             6d23s23
        if(nlast.ne.ngot)then                                           6d23s23
         looper=1                                                       6d23s23
         nlast=ngot                                                     6d23s23
        end if                                                          6d23s23
       end if                                                           6d23s23
       if(nwant.ne.ngot)go to 1                                         8d10s22
      end if                                                            8d10s22
      return                                                            8d10s22
      end                                                               8d10s22
c mpec2.1 version zeta copyright u.s. government
      subroutine ddi_get(bc,ibc,nhandlex,ilow1,ihigh1,ilow2,ihigh2,buff)11d15s22
      implicit real*8 (a-h,o-z)                                         1d18s21
      implicit integer*8 (i-n)                                          1d18s21
      dimension buff(*)                                                 1d18s21
      nhandle=nhandlex+1                                                1d18s21
      call ddi_getx(bc,ibc,nhandle,ilow1,ihigh1,ilow2,ihigh2,buff)      11d15s22
      return                                                            1d18s21
      end                                                               1d18s21
c mpec2.1 version zeta copyright u.s. government
      subroutine ddi_done(ircv,nrcv)                                    1d29s21
      use mpi                                                           1d29s21
      parameter (id=10)
      dimension ircv(*),istatus(id)                                      1d29s21
      do i=1,nrcv                                                       1d29s21
       call mpi_wait(ircv(i),istatus,ierror)                             1d29s21
      end do                                                            1d29s21
      nrcv=0                                                            1d29s21
      return                                                            1d29s21
      end                                                               1d29s21

cddi4c
cddi4c     parallelization via mpi.
cddi4c     use one-sided memory access via mpi-3 commands.
cddi4c     I think this will be limited to a single node.
cddi4c
cddi4      subroutine dws_init(ncore)                                        3d5s21
cddi4      idum=1                                                            3d5s21
cddi4      return                                                            3d5s21
cddi4      end                                                               3d5s21
cddi4      subroutine dws_preinit                                            3d5s21
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)                                         4d10s12
cddi4      implicit integer (i-n)
cddi4      character*8 file
cddi4      integer (kind=mpi_address_kind) isizex                            6d4s10
cddi4      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
cddi4     $     mynnode
cddi4      include "common.mympi"                                             3d5s21
cddi4      call mpi_init(ierror)
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in ddi_int,')
cddi4       write(6,*)('mpi_init returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4       stop
cddi4      end if
cddi4      call mpi_comm_size(MPI_COMM_WORLD,isize,ierror)
cddi4      call mpi_comm_rank(MPI_COMM_WORLD,irank,ierror)
cddi4      if(irank.eq.0)write(6,*)('Hi, this is preinit for MPI-4 code')    11d1s22
cddi4c$$$      if(irank.eq.0)then                                                6d4s10
cddi4c$$$       read(5,*)ncpus                                                   6d4s10
cddi4c$$$       xcpus=dfloat(ncpus)                                              6d4s10
cddi4c$$$      end if                                                            6d4s10
cddi4      mynowprog=irank                                                   3d5s21
cddi4      mynprocg=isize                                                    3d5s21
cddi4      ndata=0                                                           3d5s21
cddi4      if(irank.lt.10)then
cddi4       write(file,1)irank
cddi4    1  format('output.',i1)
cddi4      else if(irank.lt.100)then
cddi4       write(file,2)irank
cddi4    2  format('outpu.',i2)
cddi4      else
cddi4       write(file,3)irank
cddi4    3  format('outp.',i3)
cddi4      end if
cddi4      if(irank.ne.0)then
cddi4      open(unit=6,file=file)
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_finalize
cddi4      use mpi
cddi4      implicit integer (i-n)
cddi4      call mpi_finalize(ierror)
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in dws_finalize,')
cddi4       write(6,*)('mpi_finalize returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_gbor(buff,nwds)                                    1d26s21
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      implicit integer (i-n)
cddi4      integer buff(nwds)                                                1d26s21
cddi4      include "common.mympi"                                            1d29s21
cddi4c$$$      common/mpicom/my_comm_group,mymemp
cddi4c$$$      write(6,*)('my gsumf input buffer '),buff
cddi4      nwds2=nwds*2                                                      1d26s21
cddi4      call mpi_allreduce(mpi_in_place,buff,nwds2,mpi_integer,           1d26s21                                                                2d22s10
cddi4     $     mpi_bor,MPI_COMM_WORLD,ierror)                                1d26s21
cddi4c$$$      write(6,*)('my gsumf output buffer '),buff
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in dws_gbor,')
cddi4       write(6,*)('mpi_allreduce returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_bcast(buf,len)
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      implicit integer (i-n)
cddi4      integer icont
cddi4      dimension buf(len)
cddi4      entry dws_bcasta(buf,len)                                         3d5s21
cddi4      icont=len
cddi4      isource=0
cddi4c$$$      write(6,*)('in dws_bcast for len '),len,loc(buf)
cddi4      call mpi_bcast(buf,icont,mpi_double_precision,isource,            2d19s10
cddi4     $               MPI_COMM_WORLD,ierror)                             2d19s10
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in dws_bcast,')
cddi4       write(6,*)('mpi_bcast returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_sync
cddi4      use mpi
cddi4      entry dws_synca                                                   3d5s21
cddi4      call mpi_barrier(MPI_COMM_WORLD,ierror)                           2d22s10
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in dws_sync,')
cddi4       write(6,*)('mpi_barrier returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_gsumf(buff,nwds)
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      implicit integer (i-n)
cddi4      dimension buff(nwds)                                                 2d22s19
cddi4c$$$      write(6,*)('my gsumf input buffer '),buff
cddi4      call mpi_allreduce(mpi_in_place,buff,nwds,mpi_double_precision,   2d22s10
cddi4     $     mpi_sum,MPI_COMM_WORLD,ierror)                               2d22s10
cddi4c$$$      write(6,*)('my gsumf output buffer '),buff
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in dws_gsumf,')
cddi4       write(6,*)('mpi_allreduce returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_win0(mstor,buff,iwin)
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      implicit integer (i-n)
cddi4      integer (kind=mpi_address_kind) size
cddi4      common/mycomu/my_comml
cddi4      call mpi_info_create(info,ierror)
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in dws_win0,')
cddi4       write(6,*)('mpi_info_create returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4       stop
cddi4      end if
cddi4c$$$      call mpi_info_set(info,'no_locks','true',ierror)
cddi4      call mpi_info_set(info,'no_locks','false',ierror)                                                                        
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in dws_win0_create,')
cddi4       write(6,*)('mpi_info_set returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4       stop
cddi4      end if
cddi4      istepsz=8                                                         
cddi4      size=mstor*istepsz
cddi4      call mpi_win_create(buff,size,istepsz,info,
cddi4     $     my_comml,iwin,ierror)                                        3d30s09
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in ddi_create,')
cddi4       write(6,*)('mpi_win_create returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_put0(buff,ioff,nwds,iwin)
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      integer (kind=mpi_address_kind) ioffk                             6d9s10
cddi4      common/mycomu/my_comml
cddi4      itarg=0
cddi4      ioffk=ioff                                                        6d9s10
cddi4      call mpi_put(buff,nwds,mpi_double_precision,itarg,ioffk,nwds,     6d9s10
cddi4     $     mpi_double_precision,iwin,ierror)
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in ddi_put0,')
cddi4       write(6,*)('mpi_win_create returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_get0(buff,ioff,nwds,iwin)
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      integer (kind=mpi_address_kind) ioffk                             6d9s10
cddi4      common/mycomu/my_comml
cddi4      itarg=0
cddi4      ioffk=ioff                                                        6d9s10
cddi4      call mpi_get(buff,nwds,mpi_double_precision,itarg,ioffk,nwds,     6d9s10
cddi4     $     mpi_double_precision,iwin,ierror)
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in ddi_get0,')
cddi4       write(6,*)('mpi_get returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_put1(iproc,buff,nwds)                              11d13s12
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      integer (kind=mpi_address_kind) ioffk                             6d9s10
cddi4      itarg=iproc                                                       11d13s12
cddi4      itag=0                                                            11d13s12
cddi4      call mpi_send(buff,nwds,mpi_double_precision,itarg,itag,          11d13s12
cddi4     $     MPI_COMM_WORLD,ierror)                                       11d13s12
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in ddi_put1,')
cddi4       write(6,*)('mpi_win_create returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_sput1(iproc,buff,nwds)                              11d13s12
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      integer (kind=mpi_address_kind) ioffk                             6d9s10
cddi4      itarg=iproc                                                       11d13s12
cddi4      itag=0                                                            11d13s12
cddi4      call mpi_ssend(buff,nwds,mpi_double_precision,itarg,itag,          11d13s12
cddi4     $     MPI_COMM_WORLD,ierror)                                       11d13s12
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in ddi_sput1,')
cddi4       write(6,*)('mpi_win_create returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_get1(iproc,buff,nwds)                              11d13s12
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      integer (kind=mpi_address_kind) ioffk                             6d9s10
cddi4      integer status(mpi_status_size)
cddi4      itarg=iproc                                                       11d13s12
cddi4      itag=0                                                            11d13s12
cddi4      call mpi_recv(buff,nwds,mpi_double_precision,itarg,itag,
cddi4     $     MPI_COMM_WORLD,status,ierror)                                       11d13s12                                                                        
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in ddi_get1,')
cddi4       write(6,*)('mpi_get returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_putl(iproc,buff,nwds)                              11d13s12
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      integer (kind=mpi_address_kind) ioffk                             6d9s10
cddi4      common/mycomu/my_comml
cddi4      itarg=iproc                                                       11d13s12
cddi4      itag=0                                                            11d13s12
cddi4      call mpi_send(buff,nwds,mpi_double_precision,itarg,itag,          11d13s12
cddi4     $     my_comml,ierror)                                             11d16s12
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in ddi_putl,')
cddi4       write(6,*)('mpi_win_create returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_getl(iproc,buff,nwds)                              11d13s12
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      integer (kind=mpi_address_kind) ioffk                             6d9s10
cddi4      integer status(mpi_status_size)
cddi4      common/mycomu/my_comml
cddi4      itarg=iproc                                                       11d13s12
cddi4      itag=0                                                            11d13s12
cddi4      call mpi_recv(buff,nwds,mpi_double_precision,itarg,itag,
cddi4     $     my_comml,status,ierror)                                      11d16s12
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in ddi_getl,')
cddi4       write(6,*)('mpi_get returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_nowin0(iwin)                                       12d13s11
cddi4      use mpi                                                           12d13s11
cddi4      implicit real*8 (a-h,o-z)                                         12d13s11
cddi4      call mpi_win_free(iwin,ierror)                                    12d13s11
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in ddi_nowin0,')
cddi4       write(6,*)('mpi_win_free returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      
cddi4      subroutine dws_all2allvl(buf,nblock,ioff)
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      dimension buf(1),nblock(1),ioff(1)
cddi4      common/mycomu/my_comml
cddi4c$$$      write(6,*)('in dws_all2allvl '),buf(1),nblock(1),ioff(1)
cddi4c$$$      if(buf(1).ne.-132d0)then
cddi4c$$$       call dws_sync
cddi4c$$$       call dws_finalize
cddi4c$$$       stop
cddi4c$$$      end if
cddi4      call mpi_alltoallv(mpi_in_place,nblock,ioff,mpi_double_precision,
cddi4     $     buf,nblock,ioff,mpi_double_precision,mpi_comm_world,ierr)    6d28s12
cddi4      if(ierr.ne.mpi_success)then
cddi4       write(6,*)('in dws_all2allvl, mpi_alltoallv returned an error'),
cddi4     $      ierror
cddi4       close(unit=6)
cddi4c$$$       stop
cddi4      end if
cddi4c$$$      write(6,*)('after words...')
cddi4c$$$      if(buf(1).ne.-132d0)then
cddi4c$$$       call dws_sync
cddi4c$$$       call dws_finalize
cddi4c$$$       stop
cddi4c$$$      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_all2allvb(bufs,nblock8s,ioff8s,bufr,nblock8r,
cddi4     $     ioff8r)
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      implicit integer (i-n)
cddi4      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
cddi4     $     mynnode
cddi4      dimension bufs(1),nblock8s(1),ioff8s(1),bufr(1),nblock8r(1),
cddi4     $     ioff8r(1)
ccccccccddi4      parameter (id=1000)
cddi4      include "common.mympi"                                            1d21s21
cddi4      dimension nblocks(id),ioffs(id),nblockr(id),ioffr(id)
cddi4      data icall/0/
cddi4      save
cddi4      icall=icall+1
cddi4c$$$      write(6,*)('in all2allvb for icall = '),icall
cddi4      if(mynprocg.gt.id)then
cddi4       write(6,*)('too many procs in dws_allgv!! '),mynprocg,id
cddi4       close(unit=6)
cddi4c$$$       stop
cddi4      end if
cddi4c$$$      write(6,*)('in dws_allgv, mynprocg = '),mynprocg
cddi4      do i=1,mynprocg
cddi4       nblocks(i)=nblock8s(i)
cddi4       ioffs(i)=ioff8s(i)
cddi4       nblockr(i)=nblock8r(i)
cddi4       ioffr(i)=ioff8r(i)
cddi4c$$$       write(6,3030)i,nblocks(i),ioffs(i),nblockr(i),ioffr(i)
cddi4 3030  format(5i8)
cddi4      end do
cddi4c$$$      if(icall.eq.1.and.bufs(132).ne.132d0)then
cddi4c$$$       call dws_sync
cddi4c$$$       call dws_finalize
cddi4c$$$       stop
cddi4c$$$      end if
cddi4c$$$      if(ioffs(3).ne.-132)then
cddi4c$$$       call dws_sync
cddi4c$$$       call dws_finalize
cddi4c$$$       stop
cddi4c$$$      end if
cddi4      call mpi_alltoallv(bufs,nblocks,ioffs,mpi_double_precision,
cddi4     $     bufr,nblockr,ioffr,mpi_double_precision,mpi_comm_world,ierr)
cddi4      if(ierr.ne.mpi_success)then
cddi4       write(6,*)('in dws_all2allv, mpi_alltoallv returned an error'),
cddi4     $      ierror
cddi4       close(unit=6)
cddi4c$$$       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_all2allvb8(bufs,nblock8s,ioff8s,bufr,nblock8r,
cddi4     $     ioff8r)
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      implicit integer (i-n)
cddi4      integer*8 nblock8s,ioff8s,nblock8r,ioff8r                         1d10s18
cddi4      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
cddi4     $     mynnode
cddi4      dimension bufs(1),nblock8s(1),ioff8s(1),bufr(1),nblock8r(1),
cddi4     $     ioff8r(1)
ccccccccddi4      parameter (id=1000)
cddi4      include "common.mympi"                                            1d21s21
cddi4      dimension nblocks(id),ioffs(id),nblockr(id),ioffr(id)
cddi4      data icall/0/
cddi4      save
cddi4      icall=icall+1
cddi4c$$$      write(6,*)('in dws_all2allvb8 ')
cddi4c$$$      write(6,*)('for call no. '),icall
cddi4      if(mynprocg.gt.id)then
cddi4       write(6,*)('too many procs in dws_allgv!! '),mynprocg,id
cddi4       close(unit=6)
cddi4c$$$       stop
cddi4      end if
cddi4c$$$      write(6,*)('in dws_allgv, mynprocg = '),mynprocg
cddi4      do i=1,mynprocg
cddi4       nblocks(i)=nblock8s(i)
cddi4       ioffs(i)=ioff8s(i)
cddi4       nblockr(i)=nblock8r(i)
cddi4       ioffr(i)=ioff8r(i)
cddi4c$$$       write(6,3030)i,nblocks(i),ioffs(i),nblockr(i),ioffr(i)
cddi4 3030  format(5i8)
cddi4      end do
cddi4c$$$      if(icall.ne.10)then
cddi4c$$$       call dws_sync
cddi4c$$$       call dws_finalize
cddi4c$$$       stop
cddi4c$$$      end if
cddi4c$$$      if(ioffs(3).ne.-132)then
cddi4c$$$       call dws_sync
cddi4c$$$       call dws_finalize
cddi4c$$$       stop
cddi4c$$$      end if
cddi4      call mpi_alltoallv(bufs,nblocks,ioffs,mpi_double_precision,
cddi4     $     bufr,nblockr,ioffr,mpi_double_precision,mpi_comm_world,ierr)
cddi4      if(ierr.ne.mpi_success)then
cddi4       write(6,*)('in dws_all2allv, mpi_alltoallv returned an error'),
cddi4     $      ierror
cddi4       close(unit=6)
cddi4c$$$       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_12all(buf,len,isource8)
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      implicit integer (i-n)
cddi4      dimension buf(len)
cddi4      data icall/0/
cddi4      save icall
cddi4      icall=icall+1
cddi4      icont=len
cddi4      isource=isource8
cddi4c$$$      write(6,*)('in dws_12all, send '),icont,(' from '),isource
cddi4c$$$      write(6,*)('call no. '),icall
cddi4      call mpi_bcast(buf,icont,mpi_double_precision,isource,            2d19s10
cddi4     $               MPI_COMM_WORLD,ierror)                             2d19s10
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in dws_12all,')
cddi4       write(6,*)('mpi_bcast returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4      close(unit=6)
cddi4c$$$       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_12all_loc(buf,len,isource8)                        11d16s12
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      implicit integer (i-n)
cddi4      common/mycomu/my_comml
cddi4      dimension buf(len)
cddi4      data icall/0/
cddi4      save icall
cddi4      icall=icall+1
cddi4      icont=len
cddi4      isource=isource8
cddi4c$$$      write(6,*)('in dws_12all, send '),icont,(' from '),isource
cddi4c$$$      write(6,*)('call no. '),icall
cddi4      call mpi_bcast(buf,icont,mpi_double_precision,isource,            2d19s10
cddi4     $               my_comml,ierror)                                   11d16s12
cddi4      if(ierror.ne.mpi_success)then
cddi4       write(6,*)('in dws_12all_loc,')                                  11d16s12
cddi4       write(6,*)('mpi_bcast returned an error ')
cddi4       write(6,*)('ierror '),ierror
cddi4       write(6,*)('mpi_success '),mpi_success
cddi4      close(unit=6)
cddi4c$$$       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine dws_12all_loc_dum(buf,len,isource8)                        11d16s12
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      implicit integer (i-n)
cddi4      common/mycomu/my_comml
cddi4      dimension buf(len)
cddi4      data icall/0/
cddi4      save icall
cddi4      icall=icall+1
cddi4      return
cddi4      end
cddi4      subroutine dws_allgv2(buf,nblock8,ioff8,send)
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      implicit integer (i-n)
cddi4      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
cddi4     $     mynnode
cddi4      include "common.mympi"                                            1d29s21
cddi4      dimension buf(1),nblock8(1),ioff8(1),send(1)
cccccccddi4      parameter (id=1000)
cddi4      dimension nblock(id),ioff(id)
cddi4      if(mynprocg.gt.id)then
cddi4       write(6,*)('too many procs in dws_allgv!! '),mynprocg,id
cddi4       close(unit=6)
cddi4c$$$       stop
cddi4      end if
cddi4c$$$      write(6,*)('in dws_allgv, mynprocg = '),mynprocg
cddi4      do i=1,mynprocg
cddi4       nblock(i)=nblock8(i)
cddi4       ioff(i)=ioff8(i)
cddi4c$$$       write(6,*)i,nblock(i),ioff(i)
cddi4      end do
cddi4      idum=nblock(mynowprog+1)
cddi4c$$$      if(idum.ne.-132)then
cddi4c$$$       call dws_sync
cddi4c$$$       call dws_finalize
cddi4c$$$       stop
cddi4c$$$      end if
cddi4c$$$      write(6,*)('sending '),idum
cddi4      call mpi_allgatherv(send,idum,mpi_double_precision,buf,
cddi4     $     nblock,ioff,mpi_double_precision,mpi_comm_world,ierr)
cddi4      if(ierr.ne.mpi_success)then
cddi4       write(6,*)('in dws_allgv, mpi_allgatherv returned an error'),
cddi4     $      ierror
cddi4       close(unit=6)
cddi4c$$$       stop
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine second(time1)
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      data ifirst/0/
cddi4      save
cddi4      time1=mpi_wtime()                                                 5d7s12
cddi4      if(ifirst.eq.0)then
cddi4       time0=time1
cddi4       time1=0d0
cddi4       ifirst=1
cddi4      else
cddi4       time1=time1-time0
cddi4      end if
cddi4      return
cddi4      end
cddi4      subroutine ddi_zero(bc,ibc,nhandle)                               11d15s22
cddi4      use mpi                                                           3d8s21
cddi4      implicit real*8 (a-h,o-z)
cddi4      integer*8 nhandle
cddi4      integer (kind=mpi_address_kind) idisp
cddi4      include "common.mympi"                                             3d5s21
cddi4      include "common.store"                                            3d8s21
cddi4      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
cddi4     $     mynnode
cddi4      itmp=ibcoff                                                       3d8s21
cddi4c$$$      write(6,*)('ddi_zero for '),nhandle
cddi4c$$$      write(6,*)('mdata: '),(mdata(j,nhandle),j=1,5)
cddi4      ibcoff=itmp+mdata(4,nhandle)                                      3d8s21
cddi4      call enough('ddi.zero',bc,ibc)                                    1d18s23
cddi4      do i=itmp,ibcoff-1                                                3d8s21
cddi4       bc(i)=0d0                                                        3d8s21
cddi4      end do                                                            3d8s21
cddi4      iassert=0                                                         3d8s21
cddi4      call mpi_win_lock(mpi_lock_exclusive,mynowprog,iassert,           3d8s21
cddi4     $        mdata(5,nhandle),ierror)                                  3d5s21
cddi4      idisp=0                                                           3d8s21
cddi4      call mpi_put(bc(itmp),mdata(4,nhandle),mpi_double_precision,      3d8s21
cddi4     $     mynowprog,idisp,mdata(4,nhandle),mpi_double_precision,       3d8s21
cddi4     $     mdata(5,nhandle),ierror)                                     3d8s21
cddi4      call mpi_win_unlock(mynowprog,mdata(5,nhandle),ierror)            3d8s21
cddi4      ibcoff=itmp                                                       3d8s21
cddi4      return                                                            3d8s21
cddi4      end
cddi4      subroutine ddi_put(bc,ibc,nhandle,i1,i2,i3,i4,buff)               11d15s22
cddi4      use mpi                                                           3d5s21
cddi4      implicit real*8 (a-h,o-z)
cddi4      integer*8 nhandle,i1,i2,i3,i4
cddi4      integer (kind=mpi_address_kind) idisp
cddi4      logical lflag
cddi4      dimension buff(*)
cddi4      include "common.mympi"                                             3d5s21
cddi4      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
cddi4     $     mynnode
cddi4c$$$      write(6,*)('in ddi_put with args '),nhandle,i1,i2,i3,i4
cddi4      if(nhandle.gt.ndata)then
cddi4       write(6,*)('you are trying to put to distributed array '),
cddi4     $      nhandle
cddi4       write(6,*)('but we only have '),ndata,('distributed arrays')
cddi4       stop 'ddi_put'
cddi4      end if
cddi4c$$$      write(6,*)('what we have in mdata: '),(mdata(j,nhandle),j=1,5)
cddi4c$$$      if(ndata.ne.132)then
cddi4c$$$       call dws_sync
cddi4c$$$       call dws_finalize
cddi4c$$$       stop
cddi4c$$$      end if
cddi4      if(i1.eq.1.and.i2.eq.mdata(1,nhandle))then                        3d5s21
cddi4       ioff=1                                                           3d5s21
cddi4c$$$       write(6,*)('let us query size: ')
cddi4c$$$       call mpi_win_get_attr(mdata(5,nhandle),MPI_WIN_SIZE,nsizex,lflag,
cddi4c$$$     $      ierror)
cddi4c$$$      write(6,*)('lflag, size: '),lflag,nsizex
cddi4       do ip=0,mynprocg-1                                               3d5s21
cddi4        call ilimts(1,mdata(2,nhandle),mynprocg,ip,il,ih,i1s,i1e,i2s,   3d5s21
cddi4     $       i2e)                                                       3d5s21
cddi4        istart=il                                                       3d8s21
cddi4        if(i3.gt.il)istart=i3                                           3d8s21
cddi4        iend=ih                                                         3d8s21
cddi4        if(i4.le.iend)iend=i4                                           3d8s21
cddi4c$$$        write(6,*)('istart,iend: '),istart,iend
cddi4        if(iend.ge.istart)then                                          3d8s21
cddi4         iassert=0                                                      3d5s21
cddi4         call mpi_win_lock(mpi_lock_exclusive,ip,iassert,               3d5s21
cddi4     $        mdata(5,nhandle),ierror)                                  3d5s21
cddi4c$$$         write(6,*)('for prod '),ip
cddi4c$$$         write(6,*)('limits: '),il,ih
cddi4         istart=max(il,i3)                                              3d5s21
cddi4         iend=min(ih,i4)                                                3d5s21
cddi4         nput=mdata(1,nhandle)*(iend+1-istart)                          3d5s21
cddi4         idisp=mdata(1,nhandle)*(istart-il)                             3d5s21
cddi4c$$$         write(6,*)('istart,iend '),istart,iend
cddi4c$$$         write(6,*)('nput '),nput
cddi4c$$$         write(6,*)('idisp '),idisp
cddi4c$$$         write(6,*)('window? '),mdata(5,nhandle)
cddi4c$$$         write(6,*)('in ddi_put,saving to proc '),ip
cddi4c$$$         call prntm2(buff(ioff),mdata(1,nhandle),iend+1-istart,
cddi4c$$$     $        mdata(1,nhandle))
cddi4         call mpi_put(buff(ioff),nput,mpi_double_precision,ip,idisp,    3d5s21
cddi4     $        nput,mpi_double_precision,mdata(5,nhandle),ierror)        3d5s21
cddi4         call mpi_win_unlock(ip,mdata(5,nhandle),ierror)                3d5s21
cddi4         ioff=ioff+nput                                                 3d5s21
cddi4        end if                                                          3d5s21
cddi4       end do                                                           3d5s21
cddi4      else                                                              3d5s21
cddi4       write(6,*)('you asked to put subset of rows ... '),i1,i2
cddi4       write(6,*)('out of '),mdata(1,nhandle)
cddi4       write(6,*)('I have not coded this yet!'),nhandle
cddi4       stop 'ddi_put'
cddi4      end if                                                            3d5s21
cddi4c$$$      if(ndata.ne.132)then
cddi4c$$$       call dws_sync
cddi4c$$$       call dws_finalize
cddi4c$$$       stop
cddi4c$$$      end if
cddi4      return                                                            3d5s21
cddi4      end
cddi4      subroutine ddi_iget(bc,ibc,nhandle,i1,i2,i3,i4,buff,iget,nget)    11d15s22
cddi4      implicit real*8 (a-h,o-z)
cddi4      integer*8 nhandle,i1,i2,i3,i4                                     3d8s21
cddi4      include "common.store"                                            11d15s22
cddi4      nget=0                                                            3d8s21
cddi4      call ddi_get(bc,ibc,nhandle,i1,i2,i3,i4,buff)                     11d15s22
cddi4      return                                                            3d8s21
cddi4      end
cddi4      subroutine ddi_iacc(bc,ibc,nhandle,i1,i2,i3,i4,buff,iacc,nacc)    11d15s22
cddi4      implicit real*8 (a-h,o-z)
cddi4      integer*8 nhandle,i1,i2,i3,i4                                     3d8s21
cddi4      include "common.store"                                            11d15s22
cddi4      nacc=0                                                            3d8s21
cddi4      call ddi_acc(bc,ibc,nhandle,i1,i2,i3,i4,buff)                     11d15s22
cddi4      return
cddi4      end
cddi4      subroutine ddi_get(bc,ibc,nhandle,i1,i2,i3,i4,buff)               11d15s22
cddi4      use mpi                                                           3d5s21
cddi4      implicit real*8 (a-h,o-z)
cddi4      integer*8 nhandle,i1,i2,i3,i4
cddi4      integer (kind=mpi_address_kind) idisp
cddi4      dimension buff(*)
cddi4      include "common.mympi"                                             3d5s21
cddi4      include "common.store"                                            3d8s21
cddi4      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
cddi4     $     mynnode
cddi4      if(nhandle.gt.ndata)then
cddi4       write(6,*)('you are trying to get to distributed array '),
cddi4     $      nhandle
cddi4       write(6,*)('but we only have '),ndata,('distributed arrays')
cddi4       stop 'ddi_get'
cddi4      end if
cddi4      if(i1.eq.1.and.i2.eq.mdata(1,nhandle))then                        3d5s21
cddi4       ioff=1                                                           3d5s21
cddi4c$$$       write(6,*)('let us query size: ')
cddi4c$$$       call mpi_win_get_attr(mdata(5,nhandle),MPI_WIN_SIZE,nsizex,lflag,
cddi4c$$$     $      ierror)
cddi4c$$$      write(6,*)('lflag, size: '),lflag,nsizex
cddi4       do ip=0,mynprocg-1                                               3d5s21
cddi4        call ilimts(1,mdata(2,nhandle),mynprocg,ip,il,ih,i1s,i1e,i2s,   3d5s21
cddi4     $       i2e)                                                       3d5s21
cddi4c$$$        write(6,*)('from proc '),ip,il,i3,ih,i4
cddi4        istart=il                                                       3d8s21
cddi4        if(i3.gt.il)istart=i3                                           3d8s21
cddi4        iend=ih                                                         3d8s21
cddi4        if(i4.le.iend)iend=i4                                           3d8s21
cddi4c$$$        write(6,*)('istart,iend: '),istart,iend
cddi4        if(iend.ge.istart)then                                          3d8s21
cddi4         iassert=0                                                      3d5s21
cddi4         call mpi_win_lock(mpi_lock_exclusive,ip,iassert,               3d5s21
cddi4     $        mdata(5,nhandle),ierror)                                  3d5s21
cddi4c$$$         write(6,*)('for prod '),ip
cddi4c$$$         write(6,*)('limits: '),il,ih
cddi4         istart=max(il,i3)                                              3d5s21
cddi4         iend=min(ih,i4)                                                3d5s21
cddi4         nput=mdata(1,nhandle)*(iend+1-istart)                          3d5s21
cddi4         idisp=mdata(1,nhandle)*(istart-il)                             3d5s21
cddi4c$$$         write(6,*)('istart,iend '),istart,iend
cddi4c$$$         write(6,*)('nput '),nput
cddi4c$$$         write(6,*)('idisp '),idisp
cddi4c$$$         write(6,*)('window? '),mdata(5,nhandle)
cddi4         call mpi_get(buff(ioff),nput,mpi_double_precision,ip,idisp,    3d5s21
cddi4     $        nput,mpi_double_precision,mdata(5,nhandle),ierror)        3d5s21
cddi4c$$$         write(6,*)('I''m getting from proc '),ip
cddi4c$$$         call prntm2(buff(ioff),mdata(1,nhandle),iend+1-istart,
cddi4c$$$     $        mdata(1,nhandle))
cddi4         call mpi_win_unlock(ip,mdata(5,nhandle),ierror)                3d5s21
cddi4         ioff=ioff+nput                                                 3d5s21
cddi4        end if                                                          3d5s21
cddi4       end do                                                           3d5s21
cddi4      else                                                              3d5s21
cddi4c$$$       write(6,*)('going for limited rows '),i1,i2
cddi4       mrow=i2+1-i1                                                     3d8s21
cddi4       ioff=1                                                           3d5s21
cddi4       do ip=0,mynprocg-1                                               3d5s21
cddi4        call ilimts(1,mdata(2,nhandle),mynprocg,ip,il,ih,i1s,i1e,i2s,   3d5s21
cddi4     $       i2e)                                                       3d5s21
cddi4        istart=il                                                       3d8s21
cddi4        if(i3.gt.il)istart=i3                                           3d8s21
cddi4        iend=ih                                                         3d8s21
cddi4        if(i4.le.iend)iend=i4                                           3d8s21
cddi4        if(iend.ge.istart)then                                          3d8s21
cddi4         iassert=0                                                      3d5s21
cddi4         call mpi_win_lock(mpi_lock_exclusive,ip,iassert,               3d5s21
cddi4     $        mdata(5,nhandle),ierror)                                  3d5s21
cddi4         mcol=iend+1-istart                                             3d8s21
cddi4         nput=mdata(1,nhandle)*(iend+1-istart)                          3d5s21
cddi4         idisp=mdata(1,nhandle)*(istart-il)                             3d5s21
cddi4         itmp=ibcoff                                                    3d8s21
cddi4         ibcoff=itmp+nput                                               3d8s21
cddi4         call enough('ddi.get',bc,ibc)                                         11d15s22
cddi4         call mpi_get(bc(itmp),nput,mpi_double_precision,ip,idisp,      3d8s21
cddi4     $        nput,mpi_double_precision,mdata(5,nhandle),ierror)        3d5s21
cddi4c$$$         write(6,*)('what we got from proc '),ip
cddi4c$$$         call prntm2(bc(itmp),mdata(1,nhandle),mcol,mdata(1,nhandle))
cddi4         call mpi_win_unlock(ip,mdata(5,nhandle),ierror)                3d5s21
cddi4         do i=istart,iend                                               3d8s21
cddi4          im=i-istart                                                   3d8s21
cddi4          jtmp=itmp-1+mdata(1,nhandle)*im                               3d8s21
cddi4          joff=ioff-i1+mrow*im                                          3d8s21
cddi4          do j=i1,i2                                                    3d8s21
cddi4           buff(joff+j)=bc(jtmp+j)                                      3d8s21
cddi4          end do                                                        3d8s21
cddi4         end do                                                         3d8s21
cddi4         ioff=ioff+mrow*mcol                                            3d8s21
cddi4         ibcoff=itmp                                                    3d8s21
cddi4        end if                                                          3d5s21
cddi4       end do                                                           3d5s21
cddi4       mtot=i4+1-i3                                                     3d8s21
cddi4c$$$       write(6,*)('altogether now ')
cddi4c$$$       call prntm2(buff,mrow,mtot,mrow)
cddi4c$$$       write(6,*)('you asked to get subset of rows ... '),i1,i2
cddi4c$$$       write(6,*)('out of '),mdata(1,nhandle)
cddi4c$$$       write(6,*)('I have not coded this yet!')
cddi4c$$$       write(6,*)('what about columns? '),i3,i4,mdata(2,nhandle)
cddi4c$$$       stop 'ddi_get'
cddi4      end if                                                            3d5s21
cddi4      return                                                            3d5s21
cddi4      end
cddi4      subroutine ddi_done(nhandle,i1)
cddi4      implicit real*8 (a-h,o-z)
cddi4      integer*8 nhandle
cddi4c$$$      write(6,*)('you have reached ddi_done')
cddi4c$$$      call dws_sync
cddi4c$$$      call dws_finalize
cddi4c$$$      stop
cddi4      idum=1                                                            3d8s21
cddi4      return                                                            3d8s21
cddi4      end
cddi4      subroutine ddi_destroy(nhandle)
cddi4      use mpi                                                           3d5s21
cddi4      implicit real*8 (a-h,o-z)
cddi4      integer*8 nhandle
cddi4      include "common.mympi"                                             3d5s21
cddi4      if(nhandle.ne.ndata)then
cddi4       write(6,*)('trying to destroy distributed arrays out of order')
cddi4       write(6,*)nhandle,('vs.'),ndata
cddi4       stop 'ddi_destroy'
cddi4      end if
cddi4      call dws_sync
cddi4c$$$      write(6,*)('deleting window '),nhandle,mdata(5,nhandle)
cddi4      call mpi_win_free(mdata(5,nhandle),ierror)
cddi4      ndata=ndata-1
cddi4c$$$      if(ndata.ne.-132)then
cddi4c$$$       call dws_sync
cddi4c$$$       call dws_finalize
cddi4c$$$       stop
cddi4c$$$      end if
cddi4      return
cddi4      end
cddi4      subroutine ddi_create(bc,ibc,irow,icol,nhandle)                   11d15s22
cddi4      use mpi
cddi4      implicit real*8 (a-h,o-z)
cddi4      integer*8 nhandle,icol,irow
cddi4      include "common.mympi"                                             3d5s21
cddi4      character*10 value
cddi4      logical lflag                                                                        
cddi4      integer(kind=mpi_address_kind) nsizeb,iwinbase,nsizex
cddi4      data icall/0/                                                     3d8s21
cddi4      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
cddi4     $     mynnode
cddi4      save icall                                                        3d8s21
cddi4      icall=icall+1                                                     3d8s21
cddi4      ndata=ndata+1                                                     3d5s21
cddi4c$$$      write(6,*)('you have reached ddi_create'),ndata,irow,icol
cddi4      if(ndata.gt.id)then
cddi4       write(6,*)('you are trying to create distributed array no. '),
cddi4     $      ndata
cddi4       write(6,*)('but maximum dimensions in common.mympi is '),id
cddi4       stop 'ddi_create'
cddi4      end if                                                            3d5s21
cddi4      nhandle=ndata                                                     3d5s21
cddi4      mdata(1,ndata)=irow                                               3d5s21
cddi4      mdata(2,ndata)=icol                                               3d5s21
cddi4      call ilimts(1,icol,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)      3d5s21
cddi4      nhere=ih+1-il                                                     3d5s21
cddi4      nneed=nhere*irow                                                  3d5s21
cddi4c$$$      write(6,*)('ddi_create '),irow,icol,nhere,nneed
cddi4      call mpi_info_create(info,ierror)
cddi4c$$$      call mpi_info_set(info,'no_locks','true',ierror)
cddi4      call mpi_info_set(info,'no_locks','false',ierror)
cddi4c$$$      call mpi_info_get(info,'no_locks',10,value,lflag,ierror)
cddi4c$$$      write(6,*)('ierror after info_get '),ierror
cddi4c$$$      write(6,*)('lflag from info_get '),lflag
cddi4c$$$      if(lflag)then
cddi4c$$$       write(6,*)('value = "'),value,('"')
cddi4c$$$      end if
cddi4      nsizeb=nneed*8                                                    3d5s21
cddi4      ndisp=8                                                           3d5s21
cddi4c$$$      write(6,*)('creating '),nneed,irow,icol,ndata,icall
cddi4      iwinbase=0                                                        8d8s24
cddi4      call mpi_win_allocate(nsizeb,ndisp,info,mpi_comm_world,iwinbase,  3d5s21
cddi4     $     iwin,ierror)                                                 3d5s21
cddi4c$$$      write(6,*)('after win_allocate, ierror = '),ierror
cddi4c$$$      write(6,*)('after mpi_win_allocate '),ndata,iwinbase,iwin
cddi4c$$$      write(6,*)('let us query size: ')
cddi4c$$$      call mpi_win_get_attr(iwin,MPI_WIN_SIZE,nsizex,lflag,ierror)
cddi4c$$$      write(6,*)('lflag, size: '),lflag,nsizex
cddi4      mdata(3,ndata)=iwinbase                                           3d5s21
cddi4      mdata(4,ndata)=nneed                                              3d5s21
cddi4      mdata(5,ndata)=iwin                                               3d5s21
cddi4c$$$      if(icall.gt.1)then
cddi4c$$$       call dws_sync
cddi4c$$$       call dws_finalize
cddi4c$$$       stop
cddi4c$$$      end if
cddi4      call dws_sync                                                     3d5s21
cddi4c$$$      call dws_sync
cddi4c$$$      call dws_finalize
cddi4c$$$      stop
cddi4      return                                                            3d5s21
cddi4      end
cddi4      subroutine ddi_acc(bc,ibc,nhandle,i1,i2,i3,i4,buff)               11d15s22
cddi4      use mpi                                                           3d5s21
cddi4      implicit real*8 (a-h,o-z)
cddi4      integer*8 nhandle,i1,i2,i3,i4
cddi4      integer (kind=mpi_address_kind) idisp
cddi4      logical lflag
cddi4      dimension buff(*)
cddi4      include "common.mympi"                                             3d5s21
cddi4      include "common.store"                                            3d8s21
cddi4      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
cddi4     $     mynnode
cddi4      if(nhandle.gt.ndata)then
cddi4       write(6,*)('you are trying to put to distributed array '),
cddi4     $      nhandle
cddi4       write(6,*)('but we only have '),ndata,('distributed arrays')
cddi4       stop 'ddi_acc'
cddi4      end if
cddi4c$$$      write(6,*)('Hi, my name is ddi_acc!')
cddi4c$$$      if(bc(132).ne.-132d0)then
cddi4c$$$       call dws_synca
cddi4c$$$       call dws_finalize
cddi4c$$$       stop
cddi4c$$$      end if
cddi4      if(i1.eq.1.and.i2.eq.mdata(1,nhandle))then                        3d5s21
cddi4c$$$       write(6,*)('ddi_acc block 1')
cddi4       ioff=1                                                           3d5s21
cddi4c$$$       write(6,*)('let us query size: ')
cddi4c$$$       call mpi_win_get_attr(mdata(5,nhandle),MPI_WIN_SIZE,nsizex,lflag,
cddi4c$$$     $      ierror)
cddi4c$$$      write(6,*)('lflag, size: '),lflag,nsizex
cddi4       do ip=0,mynprocg-1                                               3d5s21
cddi4c$$$        write(6,*)('for ip = '),ip
cddi4        call ilimts(1,mdata(2,nhandle),mynprocg,ip,il,ih,i1s,i1e,i2s,   3d5s21
cddi4     $       i2e)                                                       3d5s21
cddi4        istart=il                                                       3d8s21
cddi4        if(i3.gt.il)istart=i3                                           3d8s21
cddi4        iend=ih                                                         3d8s21
cddi4        if(i4.le.iend)iend=i4                                           3d8s21
cddi4c$$$        write(6,*)('istart,iend: '),istart,iend
cddi4        if(iend.ge.istart)then                                          3d8s21
cddi4         iassert=0                                                      3d5s21
cddi4c$$$         write(6,*)('lock '),mdata(5,nhandle),ip,nhandle
cddi4c$$$         if(mdata(5,nhandle).ne.-132.and.bc(132).ne.-132d0)then
cddi4c$$$          call dws_synca
cddi4c$$$          call dws_finalize
cddi4c$$$          stop
cddi4c$$$         end if
cddi4         call mpi_win_lock(mpi_lock_exclusive,ip,iassert,               3d5s21
cddi4     $        mdata(5,nhandle),ierror)                                  3d5s21
cddi4c$$$         write(6,*)('after lock, ierror = '),ierror
cddi4c$$$         if(ierror.ne.-132.and.bc(132).ne.-132d0)then
cddi4c$$$          call dws_synca
cddi4c$$$          call dws_finalize
cddi4c$$$          stop
cddi4c$$$         end if
cddi4c$$$         write(6,*)('for prod '),ip
cddi4c$$$         write(6,*)('limits: '),il,ih
cddi4         istart=max(il,i3)                                              3d5s21
cddi4         iend=min(ih,i4)                                                3d5s21
cddi4         nput=mdata(1,nhandle)*(iend+1-istart)                          3d5s21
cddi4         idisp=mdata(1,nhandle)*(istart-il)                             3d5s21
cddi4c$$$         write(6,*)('istart,iend '),istart,iend
cddi4c$$$         write(6,*)('nput '),nput
cddi4c$$$         write(6,*)('idisp '),idisp
cddi4c$$$         write(6,*)('window? '),mdata(5,nhandle)
cddi4c$$$         write(6,*)('in ddi_put,saving to proc '),ip
cddi4c$$$         call prntm2(buff(ioff),mdata(1,nhandle),iend+1-istart,
cddi4c$$$     $        mdata(1,nhandle))
cddi4c$$$         write(6,*)('accumulate '),nput,idisp
cddi4         call mpi_accumulate(buff(ioff),nput,mpi_double_precision,ip,   3d8s21
cddi4     $        idisp,nput,mpi_double_precision,mpi_sum,mdata(5,nhandle), 3d8s21
cddi4     $        ierror)                                                   3d8s21
cddi4c$$$         write(6,*)('unlock ')
cddi4         call mpi_win_unlock(ip,mdata(5,nhandle),ierror)                3d5s21
cddi4         ioff=ioff+nput                                                 3d5s21
cddi4        end if                                                          3d5s21
cddi4       end do                                                           3d5s21
cddi4      else                                                              3d5s21
cddi4c$$$       write(6,*)('ddi_acc block 2')
cddi4       mrow=i2+1-i1                                                     3d8s21
cddi4       ioff=1                                                           3d5s21
cddi4       do ip=0,mynprocg-1                                               3d5s21
cddi4c$$$        write(6,*)('for proc '),ip
cddi4        call ilimts(1,mdata(2,nhandle),mynprocg,ip,il,ih,i1s,i1e,i2s,   3d5s21
cddi4     $       i2e)                                                       3d5s21
cddi4        istart=il                                                       3d8s21
cddi4        if(i3.gt.il)istart=i3                                           3d8s21
cddi4        iend=ih                                                         3d8s21
cddi4        if(i4.le.iend)iend=i4                                           3d8s21
cddi4c$$$        write(6,*)('end,start: '),iend,istart
cddi4        if(iend.ge.istart)then                                          3d8s21
cddi4         iassert=0                                                      3d5s21
cddi4c$$$         write(6,*)('lock ')
cddi4         call mpi_win_lock(mpi_lock_exclusive,ip,iassert,               3d5s21
cddi4     $        mdata(5,nhandle),ierror)                                  3d5s21
cddi4         nput=mdata(1,nhandle)*(iend+1-istart)                          3d5s21
cddi4         idisp=mdata(1,nhandle)*(istart-il)                             3d5s21
cddi4         itmp=ibcoff                                                    3d8s21
cddi4         ibcoff=itmp+nput                                               3d8s21
cddi4         call enough('ddi.acc',bc,ibc)                                         11d15s22
cddi4         do i=itmp,ibcoff-1                                             3d8s21
cddi4          bc(i)=0d0                                                     3d8s21
cddi4         end do                                                         3d8s21
cddi4         do i=istart,iend                                               3d8s21
cddi4          im=i-istart                                                   3d8s21
cddi4          jtmp=itmp-1+mdata(1,nhandle)*im                               3d8s21
cddi4          joff=ioff-i1+mrow*im                                          3d8s21
cddi4          do j=i1,i2                                                    3d8s21
cddi4           bc(jtmp+j)=buff(joff+j)                                      3d8s21
cddi4          end do                                                        3d8s21
cddi4         end do                                                         3d8s21
cddi4         ioff=ioff+mrow*(iend+1-istart)                                 3d8s21
cddi4c$$$         write(6,*)('accumulate '),nput,idisp
cddi4         call mpi_accumulate(bc(itmp),nput,mpi_double_precision,ip,     3d8s21
cddi4     $        idisp,nput,mpi_double_precision,mpi_sum,mdata(5,nhandle), 3d8s21
cddi4     $        ierror)                                                   3d8s21
cddi4c$$$         write(6,*)('unlock ')
cddi4         call mpi_win_unlock(ip,mdata(5,nhandle),ierror)                3d5s21
cddi4         ibcoff=itmp                                                    3d8s21
cddi4        end if                                                          3d5s21
cddi4       end do                                                           3d5s21
cddi4      end if                                                            3d5s21
cddi4      return
cddi4      end                  
cddi4      subroutine ddi_iaccword0
cddi4      idum=1
cddi4      return
cddi4      end
cddi4      subroutine ddi_iaccword2
cddi4      idum=1
cddi4      return
cddi4      end
cddi4
