c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcbk(gb,ncsftb,veck,ncsftk,nff0b,nff0k,ncsf,ncsf2,     7d22s21
     $     mff0b,                                                       7d22s21
     $     mff0k,nrootz,mdon,mdoop,ixw1,ixw2,nec,ism,irel,irefo,        5d12s21
     $     norb,lxmt,ixmt,multh,phase,nbasdws,idoubo,isymop,n2e,i2eop,  5d20s21
     $     phase2,nsymb,shift,mff1,mff2,iigsb,nff1b,iff1b,ihddiagb,     8d21s21
     $     nff2b,iff2b,isymbra,ihsdiagk,nff1k,iff1k,ihddiagk,nff2k,     7d22s21
     $     iff2k,isymket,nvirt,maxbx,maxbxd,srh,sr2,nfdatb,nfdatk,gdb,  8d25s21
     $     vdk,lpr,bc,ibc)                                              11d10s22
      implicit real*8 (a-h,o-z)                                         7d9s21o
      integer*8 iigsb(mdoop,nsymb),ihsdiagk(mdoop,nsymb),i18,i28,       8d21s21
     $     i38,i48,i58,ihddiagb(mdoop,nsymb),ihddiagk(mdoop,nsymb,2)    8d21s21
      external second                                                   5d4s22
      logical lpr                                                       8d23s21
      dimension nff1b(mdoop,nsymb,*),iff1b(*),nff1k(mdoop,nsymb,*),     7d9s21
     $     iff1k(*),nvirt(*), gb(ncsftb,*),veck(ncsftk,nrootz),mff0k(*),
     $     mff0b(*),nff0b(mdoop,3),nff0k(mdoop,3),ncsf(*),ism(*),       7d9s21
     $     irel(*),irefo(*),multh(8,8),ixmt(8,*),nbasdws(*),            7d9s21
     $     idoubo(*),i2eop(2,3),isymop(*),ixmtf(8),nff2b(mdoop,nsymb,*),7d22s21
     $     nff2k(mdoop,nsymb,*),iff2b(*),iff2k(*),ncsf2(4,*),
     $     nfdatb(5,4,*),nfdatk(5,4,*),vdk(*),gdb(*)                    8d25s21
      include "common.store"                                            7d12s21
      common/cpucom/tovr,top(10),tso(11)                                5d4s22
      data icall/0/
      data loop,loopx/0,1000000/
      save icall
      icall=icall+1
      if(lpr)write(6,*)('in hcbk for call no. '),icall
      call dws_synca                                                    8d11s22
      call ddi_iaccword0                                                8d11s22
      nwiacc=0                                                          8d11s22
      ibcoffo=ibcoff                                                    7d12s21
      if(lpr)write(6,*)('calling hccsfbk')
      call second(time1)                                                5d4s22
      call hccsfbk(gb,ncsftb,veck,ncsftk,nff0b,nff0k,ncsf,mff0b,        7d9s21
     $     mff0k,nrootz,mdon,mdoop,ixw1,ixw2,nec,ism,irel,irefo,        5d12s21
     $     norb,lxmt,ixmt,multh,phase,nbasdws,idoubo,isymop,n2e,i2eop,  5d20s21
     $     phase2,nsymb,shift,ixmtf,bc,ibc)                             11d10s22
      call second(time2)                                                5d4s22
      telap=time2-time1-tovr                                            5d4s22
      top(1)=top(1)+telap                                               5d4s22
      time1=time2                                                       5d4s22
      if(mff1.gt.0)then                                                 7d9s21
      if(lpr)write(6,*)('calling hcsibk')
       call hcsibk(iigsb,nff1b,iff1b,nff0k,mff0k,veck,ncsftk,ncsf,      8d21s21
     $     nec,                                                         7d9s21
     $     mdon,mdoop,nsymb,multh,ixw1,ixw2,lxmt,ixmt,nvirt,1,.false.,
     $     maxbx,isymop,n2e,i2eop,phase,phase2,ism,irel,irefo,norb,     7d9s21
     $     nbasdws,idoubo,ixmtf,nrootz,isymbra,nwiacc,bc,ibc)           1d26s23
      call second(time2)                                                5d4s22
      telap=time2-time1-tovr                                            5d4s22
      top(2)=top(2)+telap                                               5d4s22
      time1=time2                                                       5d4s22
      if(lpr)write(6,*)('calling hcisbk')
       call hcisbk(ihsdiagk,nff1k,iff1k,nff0b,mff0b,gb,ncsftb,ncsf,nec, 7d9s21
     $      mdon,mdoop,nsymb,multh,ixw1,ixw2,lxmt,ixmt,nvirt,nrootz,    7d9s21
     $      .false.,isymop,n2e,i2eop,phase,phase2,ism,irel,irefo,norb,  7d9s21
     $      nbasdws,idoubo,ixmtf,isymket,bc,ibc)                        11d10s22
       call second(time2)                                                5d4s22
       telap=time2-time1-tovr                                            5d4s22
       top(3)=top(3)+telap                                               5d4s22
       time1=time2                                                       5d4s22
      if(lpr)write(6,*)('calling hcssbk')
       call hcssbk(iigsb,nff1b,iff1b,isymbra,ihsdiagk,nff1k,iff1k,      8d21s21
     $      isymket,ncsf,nec,mdon,mdoop,nsymb,multh,ixw1,ixw2,lxmt,ixmt,7d12s21
     $      phase,phase2,nrootz,ism,irel,irefo,idoubo,nvirt,nbasdws,    7d12s21
     $      norb,.false.,maxbx,isymop,n2e,i2eop,ixmtf,bc,ibc,nwiacc)    11d21s22
       call second(time2)                                                5d4s22
       telap=time2-time1-tovr                                            5d4s22
       top(4)=top(4)+telap                                               5d4s22
       time1=time2                                                       5d4s22
      end if                                                            7d9s21
      if(mff2.gt.0)then                                                 8d3s21
       izero=1
       if(n2e.gt.0)then                                                 7d22s21
      if(lpr)write(6,*)('calling hcdibkuc')
        call hcdibkuc(ihddiagb,nff2b,iff2b,nff0k,mff0k,veck,ncsftk,     7d22s21
     $      ncsf,ncsf2,nec,mdon,mdoop,nsymb,multh,ixw1,ixw2,isymop,n2e, 7d22s21
     $      i2eop,phase2,nvirt,1,.false.,maxbxd,srh,isymbra,ism,        7d23s21
     $      irel,idoubo,irefo,nbasdws,ixmt,norb,nrootz,bc,ibc,nwiacc)   11d21s22
        call second(time2)                                                5d4s22
        telap=time2-time1-tovr                                            5d4s22
        top(5)=top(5)+telap                                               5d4s22
        time1=time2                                                       5d4s22
        izero=0
      if(lpr)write(6,*)('calling hcidbkuc')
        call hcidbkuc(ihddiagk,nff2k,iff2k,nff0b,mff0b,gb,ncsftb,       7d22s21
     $      ncsf,ncsf2,nec,mdon,mdoop,nsymb,multh,ixw1,ixw2,isymop,n2e, 7d22s21
     $      i2eop,phase2,nvirt,.false.,maxbxd,nrootz,srh,isymbra,ism,   7d22s21
     $      irel,idoubo,irefo,nbasdws,ixmt,norb,bc,ibc)                 11d10s22
        call second(time2)                                                5d4s22
        telap=time2-time1-tovr                                            5d4s22
        top(6)=top(6)+telap                                               5d4s22
        time1=time2                                                       5d4s22
       end if                                                           7d22s21
      if(lpr)write(6,*)('calling hcdducbk')
       call hcdducbk(ihddiagb,ihddiagk,nff2b,iff2b,nff2k,iff2k,ncsf,    8d24s21
     $     ncsf2,nec,mdon,mdoop,nsymb,multh,ixw1,ixw2,ixmtf,lxmt,phase, 8d24s21
     $     isymop,n2e,ixmt,i2eop,phase2,nvirt,nrootz,ism,irel,irefo,    8d24s21
     $     isymbra,isymket,norb,.false.,maxbxd,sr2,srh,timex,tovr,      8d24s21
     $     idoubo,nbasdws,izero,bc,ibc,nwiacc)                          11d21s22
       call second(time2)                                                5d4s22
       telap=time2-time1-tovr                                            5d4s22
       top(7)=top(7)+telap                                               5d4s22
       time1=time2                                                       5d4s22
      else if(mff2.lt.0)then                                            8d3s21
       izero=1                                                          8d5s21
       if(n2e.gt.0)then                                                 8d3s21
      if(lpr)write(6,*)('calling hcidbk')
        call hcidbk(nff2k,iff2k,iff2k,nsymb,mdon,mdoop,multh,isymket,   8d3s21
     $      nvirt,ncsf,ncsf2,nfdatk,irel,ism,irefo,vdk,ixw1,ixw2,        8d3s21
     $      norb,nff0b,mff0b,gb,ncsftb,nrootz,.false.,sr2,srh,n2e,      8d3s21
     $      ixmt,i2eop,isymop,phase2,idoubo,nbasdws,nec,bc,ibc)         1d13s23
        call second(time2)                                                5d4s22
        telap=time2-time1-tovr                                            5d4s22
        top(6)=top(6)+telap                                             11d17s22
        time1=time2                                                       5d4s22
      if(lpr)write(6,*)('calling hcdibk')
        call hcdibk(nff2b,iff2b,iff2b,nsymb,mdon,mdoop,multh,isymbra,   8d6s21
     $      nvirt,ncsf,ncsf2,nec,nfdatb,irel,ism,irefo,gdb,ixw1,ixw2,    8d6s21
     $      norb,nff0k,mff0k,veck,ncsftk,nrootz,.false.,srh,n2e,        8d6s21
     $      ixmt,i2eop,isymop,phase2,idoubo,nbasdws,izero,bc,ibc)       11d10s22
        call second(time2)                                                5d4s22
        telap=time2-time1-tovr                                            5d4s22
        top(5)=top(5)+telap                                               5d4s22
        time1=time2                                                       5d4s22
        izero=0                                                         8d6s21
       end if                                                           8d3s21
      if(lpr)write(6,*)('calling hcddjkbk')
       call hcddjkbk(nff2b,nfdatb,nff2k,nfdatk,gdb,vdk,nsymb,mdon,mdoop,
     $      nec,multh,isymbra,isymket,nvirt,ncsf,ncsf2,irel,ism,irefo,
     $      ixw1,ixw2,norb,nrootz,ixmtf,lxmt,phase,n2e,ixmt,isymop,
     $      i2eop,phase2,shift,sr2,srh,idoubo,nbasdws,izero,iff2b,
     $      iff2b,iff2k,iff2k,bc,ibc)                                   11d10s22
       call second(time2)                                                5d4s22
       telap=time2-time1-tovr                                            5d4s22
       top(7)=top(7)+telap                                               5d4s22
       time1=time2                                                       5d4s22
      end if                                                            7d22s21
      if(mff1.ne.0.and.mff2.gt.0)then                                   8d3s21
      if(lpr)write(6,*)('calling hcdsucbk')
       call hcdsucbk(ihsdiagk,nff1k,iff1k,ihddiagb,nff2b,iff2b,nsymb,   8d24s21
     $     mdon,mdoop,nec,multh,isymbra,isymket,nvirt,                  8d24s21
     $     ncsf,ncsf2,irel,ism,irefo,ixw1,ixw2,norb,nrootz,ixmtf,phase, 7d29s21
     $     lxmt,isymop,n2e,ixmt,i2eop,phase2,.false.,maxbx,maxbxd,sr2,  7d30s21
     $     idoubo,nbasdws,bc,ibc,nwiacc)                                11d21s22
       call second(time2)                                                5d4s22
       telap=time2-time1-tovr                                            5d4s22
       top(8)=top(8)+telap                                               5d4s22
       time1=time2                                                       5d4s22
      if(lpr)write(6,*)('calling hcsducbk')
       call hcsducbk(iigsb,nff1b,iff1b,ihddiagk,nff2k,iff2k,nsymb,      8d24s21
     $      mdon,mdoop,nec,multh,isymbra,isymket,nvirt,                 8d24s21
     $     ncsf,ncsf2,irel,ism,irefo,ixw1,ixw2,norb,nrootz,ixmtf,phase, 7d29s21
     $     lxmt,isymop,n2e,ixmt,i2eop,phase2,.false.,maxbx,maxbxd,sr2,  7d30s21
     $     idoubo,nbasdws,bc,ibc,nwiacc)                                11d21s22
       call second(time2)                                                5d4s22
       telap=time2-time1-tovr                                            5d4s22
       top(9)=top(9)+telap                                               5d4s22
       time1=time2                                                       5d4s22
      else if(mff1.ne.0.and.mff2.lt.0)then                              8d13s21
      if(lpr)write(6,*)('calling hcdsbk')
       call hcdsbk(ihsdiagk,nff1k,iff1k,nff2b,nfdatb,gdb,               8d26s21
     $      nsymb,mdon,mdoop,nec,multh,isymbra,isymket,nvirt,ncsf,ncsf2,8d16s21
     $      irel,ism,idoubo,irefo,nbasdws,ixw1,ixw2,norb,nrootz,ixmtf,  8d16s21
     $      phase,lxmt,isymop,n2e,ixmt,i2eop,phase2,.false.,maxbx,sr2,  8d16s21
     $      iff2b,iff2b,bc,ibc)                                         11d10s22
       call second(time2)                                                5d4s22
       telap=time2-time1-tovr                                            5d4s22
       top(8)=top(8)+telap                                               5d4s22
       time1=time2                                                       5d4s22
       if(lpr)write(6,*)('calling hcsdbk')
       call hcsdbk(iigsb,nff1b,iff1b,nff2k,nfdatk,vdk,                  8d26s21
     $      nsymb,mdon,mdoop,nec,multh,isymbra,isymket,nvirt,ncsf,ncsf2,8d16s21
     $      irel,ism,idoubo,irefo,nbasdws,ixw1,ixw2,norb,nrootz,ixmtf,  8d16s21
     $      phase,lxmt,isymop,n2e,ixmt,i2eop,phase2,.false.,maxbx,sr2,  8d16s21
     $      iff2k,iff2k,bc,ibc,nwiacc)                                  11d21s22
       call second(time2)                                                5d4s22
       telap=time2-time1-tovr                                            5d4s22
       top(9)=top(9)+telap                                               5d4s22
       time1=time2                                                       5d4s22
      end if                                                            7d27s21
      if(lpr)write(6,*)('all done ')
      call dws_synca
      ibcoff=ibcoffo                                                    7d12s21
      call ddi_iaccword2(nwiacc)                                        8d11s22
      return                                                            7d9s21
      end                                                               7d9s21
      subroutine lookats(strg,igs,nff1,mdon,mdoop,nsymb,multh,isym,
     $     nvirt,ncsf,nroot,bc,ibc)                                     12d12s22
      implicit real*8 (a-h,o-z)
      character*(*) strg
      integer*8 igs(mdoop,nsymb),i18,i28,i48                            12d12s22
      dimension nff1(mdoop,nsymb,2),multh(8,8),nvirt(*),ncsf(*)         12d12s22
      include "common.store"
      write(6,*)strg
      call dws_synca
      i18=1                                                             12d12s22
      do isb=1,nsymb
       isbv=multh(isb,isym)
       do np=mdon+1,mdoop
        nn=np-1
        iarg=np-mdon
        if(min(nvirt(isbv),nff1(np,isb,1)).gt.0)then
         nrow=nroot*nvirt(isbv)                                          12d12s22
         i28=nrow
         ncol=nff1(np,isb,1)*ncsf(iarg)
         i38=ncol
         call ddi_get(bc,ibc,igs(np,isb),i18,i28,i18,i38,bc(ibcoff))
         write(6,*)('for isb, np '),isb,np
         call prntm2(bc(ibcoff),nrow,ncol,nrow)
        end if
       end do
      end do
      return
      end
      subroutine lookatd(strg,gd,nfdat,nsymb,multh,isymbra,nroot,nvirt,
     $     bc,ibc)
      implicit real*8 (a-h,o-z)
      character*(*) strg
      dimension gd(*),nfdat(5,4,*),multh(8,8),nvirt(*)
      include "common.store"
      write(6,*)strg
      call dws_synca
      ioff=1
      do isb=1,nsymb
       isbv12=multh(isb,isymbra)
       do isbv1=1,nsymb
        isbv2=multh(isbv1,isbv12)
        if(isbv2.ge.isbv1)then
         if(isbv1.eq.isbv2)then
          nrow=nvirt(isbv1)*nroot
          ncol=nfdat(2,1,isb)
          if(min(nrow,ncol).gt.0)then
           write(6,*)('visv for '),isb,isbv1,ioff
           itmp=ibcoff
           ibcoff=itmp+nrow*ncol
           call enough('lookatd.1',bc,ibc)
           jtmp=itmp
           do i=0,ncol-1
            do j=0,nrow-1
             bc(jtmp+j)=gd(ioff+j)
            end do
            jtmp=jtmp+nrow
            ioff=ioff+nrow
           end do
           call dws_gsumf(bc(itmp),nrow*ncol)
           call prntm2(bc(itmp),nrow,ncol,nrow)
           ibcoff=itmp
          end if
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2
         else
          nvv=nvirt(isbv1)*nvirt(isbv2)
         end if
         do l=1,4
          if(min(nvv,nfdat(2,l,isb)).gt.0)then
           write(6,*)('vnotv for '),l,isb,isbv1,isbv2,ioff
           nrow=nvv*nroot
           ncol=nfdat(2,l,isb)
           itmp=ibcoff
           ibcoff=itmp+nrow*ncol
           call enough('lookatd.2',bc,ibc)
           jtmp=itmp
           do i=0,ncol-1
            do j=0,nrow-1
             bc(jtmp+j)=gd(ioff+j)
            end do
            jtmp=jtmp+nrow
            ioff=ioff+nrow
           end do
           call dws_gsumf(bc(itmp),nrow*ncol)
           call prntm2(bc(itmp),nrow,ncol,nrow)
           ibcoff=itmp
          end if
         end do
        end if
       end do
      end do
      return
      end
