c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine prop(nfname,ifname,iwavedat,nspc,nsymbx,nbasc,natom,   1d17s20
     $     ngaus,ibdat,isym,iapair,ibstor,isstor,idorel,ascale,multh,   12d16s19
     $     nbasp,nwavedat,lfn,iprop,bc,ibc)                             11d9s22
      implicit real*8 (a-h,o-z)
      character*200 fname                                               11d22s19
      parameter (id=150)                                                5d27s21
      external second                                                   10d8s21
c
c     superintend property calculations of ci wavefunctions
c     for spin orbit matrix elements ...
c     for linear molecules, we need only 1 component of lz pair: the
c     other can be computed from Lz acting on the wf.
c     for spherical systems, will lz acting on the input wf be good enough
c     as well? Certainly for the spin-orbit splitting of the level that
c     would be good enough, but what about coupling to other states?
c
c     diabatic states: do not compute properties, rather generate
c     overlap(s) and rotate to diabatic basis. The transform eigenvectors
c     to diabatic basis and re-write wavef.
c
      parameter (idcsfxv=30)                                            5d14s21
      integer*8 ipack8,icsfxv(idcsfxv)                                  5d14s21
      integer*4 ipack4(2)                                               11d25s19
      integer*1 ipack1(4)                                               5d12s21
      character*10 opname(id)                                           2d17s22
      character*6 esname                                                3d31s20
      character*2 code                                                  11d17s22
      character*3 onoff                                                 2d8s23
      character*(*) lfn                                                 5d25s21
      logical lprint                                                    12d20s19
      equivalence (ipack8,ipack4)                                       11d25s19
      equivalence (ipack1,npack4)                                       5d12s21
      dimension iwavedat(nspc,nwavedat),nbasc(8),nbasp(8),ipt(id),      5d12s21
     $     npt(id),ntype(4),iprop(*),                                   7d30s22
     $     opdata(id),iopdata(7,id),iosym(id),ioprt(id),i2eop(6,id),    5d27s21
     $     ixmt(8,id),isym(3,8),iapair(3,*),ibstor(*),isstor(*),iso(8), 12d16s19
     $     iorb(8),multh(8,8),opnc(id),ih0(8,8),iovr(8,8),irtyp(5),     3d9s20
     $     i4o(8,8,8,4),i4so(8,8,8,2,8),i4so2(8,8,8,2,8),ih0i(8,8),     10d20s21
     $     ih0a(8,2,8),nbasdws2(8),isopt(4,4),negs(3),ih0n(8,4),        10d5s21
     $     i4or(8,8,8,4),ionexr(8,8,8,4),jmatsr(8,8,8,4),nbaspc(8),     12d20s20
     $     kmatsr(8,8,8,4),kmatsrb(8,8,8,4),i3xr(8,8,8,4),nh0(8),       2d15s22
     $     iopso(5,33),ixmtso(8,id),iden(8),ident(8)                    1d30s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      common/cpucom/tovr,top(10),tso(11)                                5d4s22
      common/paddcm/npadddi                                             6d7s22
      include "common.store"
      include "common.mrci"
      include "common.hf"
      include "common.basis"                                            5d27s21
      common/fnd2cm/inv(2,8,8,8)                                        4d9s18
      if(n4vso.eq.0)then                                                2d8s23
       onoff='off'                                                      2d8s23
      else                                                              2d8s23
       onoff='on'                                                       2d8s23
      end if                                                            2d8s23
      ibcoffo=ibcoff                                                    11d24s19
      lprint=mynowprog.eq.0                                             12d20s19
      call second(time1)                                                5d4s22
      call second(time2)                                                5d4s22
      tovr=time2-time1                                                  5d4s22
      elowest=0d0                                                       1d19s23
      do i=1,10                                                         5d4s22
       top(i)=0d0                                                       5d4s22
       tso(i)=0d0                                                       5d4s22
      end do                                                            5d4s22
      itargs=10                                                         12d21s22
      itargd=10                                                         12d21s22
      tso(11)=0d0                                                       5d4s22
      srh=sqrt(0.5d0)                                                   7d22s21
      sr2=sqrt(2d0)                                                     7d23s21
      nfname0=nfname                                                    5d25s21
      xmasstot=0d0                                                      3d31s22
      if(lprint)then                                                    12d20s19
       write(6,*)
       write(6,*)('Hi, my name is prop, and I will be generating'),
     $      (' wave function matrix elements for you.')
       write(6,*)
       write(6,*)('norb : '),norb
       write(6,*)('doub : '),(idoubo(i),i=1,nsymb)                         11d25s19
       write(6,*)('ref  : '),(irefo(i),i=1,nsymb)                          11d25s19
       write(6,*)('nbasc: '),(nbasc(i),i=1,nsymb)                       11d25s19
       write(6,*)('nbasdws: '),(nbasdws(i),i=1,nsymb)                      11d25s19
       write(6,*)('nbasp: '),(nbasp(i),i=1,nsymb)                        12d16s19
       write(6,*)('islow, ishig: '),iprop(1),iprop(2)                   7d26s22
       write(6,*)('nec, mdon, mdoop: '),(iprop(j),j=3,5)                7d26s22
       write(6,*)('current value of ibcoff '),ibcoff
      end if                                                            12d20s19
      nwo=0                                                             12d16s19
      do i=1,nsymb                                                      12d16s19
       iorb(i)=ibcoff                                                   12d16s19
       nhere=nbasp(i)*nbasdws(i)                                        12d16s19
       nbasdws2(i)=nbasdws(i)*2                                         3d23s20
       nbaspc(i)=2*nbasp(i)                                             12d20s20
       if(idorel.ne.0)nhere=nhere*2                                     1d17s20
       ibcoff=iorb(i)+nhere                                             12d16s19
       nwo=nwo+nhere                                                    12d16s19
       nvirt(i)=nbasdws(i)-idoubo(i)-irefo(i)                           3d2s22
      end do                                                            12d16s19
      call enough('prop.  1',bc,ibc)
      noff=0                                                            5d12s21
      if(natom.ne.1.and.iprop(6).eq.0.and.iprop(7).eq.0)then            424s23
       xmasstot=0d0                                                     2d29s24
       maxbx=0                                                          2d29s24
       maxbxd=0                                                         2d29s24
       if(idorel.gt.0)then                                               3d2s22
        l2sub=1                                                          3d2s22
       else                                                              3d2s22
        l2sub=0                                                          3d2s22
       end if                                                            3d2s22
       mxopn=iprop(3)-2*iprop(4)                                         7d26s22
       mxopn=max(mxopn,2)
       if(lprint)write(6,*)('maximum no. of open shells = '),mxopn       12d20s19
       isl=iprop(1)                                                      7d26s22
       ish=iprop(2)                                                      7d26s22
       nec=iprop(3)                                                      7d26s22
       mdon=(nec-mxopn)/2                                                12d17s19
       if(lprint)                                                        2d15s22
     $   write(6,*)('minimum no. of closed shells computed to be '),mdon2d15s22
       if(lprint)write(6,*)('range of spin multiplities: '),isl,ish      12d20s19
       spin=dfloat(ish-1)*0.5d0                                          12d12s19
       xms=dfloat(isl-1)*0.5d0                                           12d12s19
       if(lprint)write(6,*)('i.e. spin = '),spin                         2d15s22
       spinxx=dfloat(mxopn)*0.5d0                                          9d15s21
       if(lprint)
     $   write(6,*)('maximum possible spin given open shells = '),spinxx
       spinx=spin+2d0                                                    9d20s21
       if(lprint)write(6,*)('maximum possible spin given resi: '),spinx  2d15s22
       spin=min(spinx,spinxx)                                            9d20s21
       if(lprint)then                                                    12d20s19
        write(6,*)('high spin: '),spin                                    12d12s19
        write(6,*)('use ms = '),xms                                       12d12s19
       end if                                                            12d20s19
       spinx=spin                                                        9d20s21
       call stepwisecsf(mxopn,spin,spin,dum,idum,idum,0,idum,bc,ibc)    11d9s22
       mdoo=iprop(5)-1                                                   7d26s22
       if(mdoo.eq.mdon)then                                              3d2s22
        mdoo=max(mdoo,nec/2)                                             3d2s22
        if(lprint)write(6,*)('resetting mdoo to '),mdoo                  3d2s22
       end if                                                            3d2s22
       nspin=(iprop(2)+2-iprop(1))/2                                     7d26s22
       icsfid=mdoo+1-mdon                                                7d27s22
       icsfpd=ibcoff                                                     6d11s19
       ibcoff=icsfpd+mdoo+1-mdon                                         6d11s19
       icsf2=ibcoff                                                      8d2s19
       icsf=icsf2+2*icsfid*nspin                                         7d27s22
       ibcoff=icsf+icsfid*nspin                                          7d27s22
       call enough('prop.  2',bc,ibc)
       do i=icsf,ibcoff-1                                                10d26s20
        ibc(i)=0                                                         10d26s20
       end do                                                            10d26s20
       call oplist(id,opname,ipt,npt,opdata,iopdata,iosym,ioprt,i2eop,   12d17s19
     $     nop,opnc,0)                                                  12d20s19
       nbb=0                                                             12d15s19
       if(idorel.eq.0)then                                               1d11s20
        ncomp=1                                                          1d11s20
       else                                                              1d11s20
        ncomp=2                                                          1d11s20
       end if                                                            1d11s20
       do i=1,nop                                                        12d15s19
        do isb=1,nsymb                                                   12d15s19
         isk=multh(isb,iosym(i))                                          12d15s19
         ixmt(isb,i)=ibcoff                                              12d15s19
         ibcoff=ixmt(isb,i)+nbasp(isb)*nbasp(isk)*ncomp*ncomp            1d11s20
         nbb=nbb+nbasp(isb)*nbasp(isk)*ncomp*ncomp                       1d11s20
        end do                                                           12d15s19
       end do                                                            12d15s19
       call enough('prop.  3',bc,ibc)
       do i=ixmt(1,1),ibcoff-1                                           12d15s19
        bc(i)=0d0                                                        12d15s19
       end do                                                            12d15s19
       idwsdeb=0                                                         12d16s19
       if(mynowprog.eq.0)then                                           11d9s22
        jfname=ifname                                                     7d8s21
        nhere=ibc(jfname)                                                 7d26s22
        jfname=jfname+1                                                   7d26s22
        if(nhere.gt.200)then                                              7d26s22
         write(6,*)('name is tooooo long!'),ifname,nhere                               7d26s22
         stop 'prop'                                                      7d26s22
        end if                                                            7d26s22
        do j=1,nhere                                                      7d26s22
         fname(j:j)=char(ibc(jfname))                                     7d26s22
         jfname=jfname+1                                                  7d26s22
        end do                                                            7d26s22
        write(6,*)('grabbing orbitals from file name: '),fname(1:nhere)   7d26s22
        call graborb(iorb,fname(1:nhere),bc,ibc)                        11d10s22
       end if                                                            7d26s22
       needt=0                                                           7d26s22
       do isb=1,nsymb                                                    7d26s22
        need=nbasp(isb)*nbasdws(isb)*ncomp                               1d11s20
        needt=needt+need                                                 7d26s22
        if(mynowprog.ne.0)then                                           7d26s22
         iorb(isb)=ibcoff                                                7d26s22
         ibcoff=iorb(isb)+need                                           7d26s22
        end if                                                           7d26s22
       end do                                                            7d26s22
       call dws_bcast(bc(iorb(1)),needt)                                 7d26s22
       if(mod(natom,2).eq.0)then                                         5d27s21
        nsnd=(natom/2)+6                                                 5d27s21
       else                                                              5d27s21
        nsnd=(natom/2)+7                                                 5d27s21
       end if                                                            5d27s21
       call dws_bcast(cmx,nsnd)                                          7d27s22
       call parap(natom,ngaus,ibdat,ixmt,isym,iapair,ibstor,isstor,iso,  12d16s19
     $     000,idorel,ascale,ipt,npt,opdata,iopdata,iosym,nop,multh,    7d26s22
     $     nbb,nbasp,nbasdws,iorb,opname,1,bc,ibc)                      11d9s22
       idwsdeb=0
       jdenpt=ibcoff                                                     6d12s19
       ibcoff=jdenpt+idbk                                                6d12s19
       call make1dm(nsymb,irefo,nbasdws,.false.,maxddi,jmats,kmats,     7d30s22
     $      jmats,ibc(jdenpt),isblk,isblkk,isblkh,nsdlk,nsdlkk,nsdlkh,  7d30s22
     $      multh,idbk,isblk1,nsdlk1,isblkd,nsdlkd,0)                   7d30s22
       ibcoff=jdenpt                                                     6d12s19
       do i=1,nsdlk                                                      6d30s18
        if(isblk(3,i).lt.isblk(4,i))then                                 7d5s18
         iswap=isblk(4,i)                                                6d30s18
         isblk(4,i)=isblk(3,i)                                           6d30s18
         isblk(3,i)=iswap                                                6d30s18
        end if                                                           6d30s18
       end do                                                            6d30s18
       do is3=1,nsymb                                                    12d16s19
        do is2=1,nsymb                                                   12d16s19
         is23=multh(is2,is3)                                             12d16s19
         do is1=1,nsymb                                                  12d16s19
          is4=multh(is1,is23)                                            12d16s19
          do i=1,nsdlk                                                   12d16s19
           if(isblk(1,i).eq.is1.and.isblk(2,i).eq.is2.and.               12d16s19
     $         isblk(3,i).eq.is3.and.isblk(4,i).eq.is4)then               12d16s19
            inv(1,is1,is2,is3)=i                                         12d16s19
            inv(2,is1,is2,is3)=1                                         12d16s19
            go to 1515                                                   12d16s19
           else if(isblk(2,i).eq.is1.and.isblk(1,i).eq.is2.and.          12d16s19
     $          isblk(4,i).eq.is3.and.isblk(3,i).eq.is4)then            12d16s19
            inv(1,is1,is2,is3)=i                                         12d16s19
            inv(2,is1,is2,is3)=2                                         12d16s19
            go to 1515                                                   12d16s19
           else if(isblk(1,i).eq.is1.and.isblk(2,i).eq.is2.and.          12d16s19
     $          isblk(4,i).eq.is3.and.isblk(3,i).eq.is4)then            12d16s19
            inv(1,is1,is2,is3)=i                                         12d16s19
            inv(2,is1,is2,is3)=3                                         12d16s19
            go to 1515                                                   12d16s19
           else if(isblk(2,i).eq.is1.and.isblk(1,i).eq.is2.and.          12d16s19
     $          isblk(3,i).eq.is3.and.isblk(4,i).eq.is4)then            12d16s19
            inv(1,is1,is2,is3)=i                                         12d16s19
            inv(2,is1,is2,is3)=4                                         12d16s19
            go to 1515                                                   12d16s19
           end if                                                        12d16s19
          end do                                                         12d16s19
 1515     continue                                                       12d16s19
         end do                                                          12d16s19
        end do                                                           12d16s19
       end do                                                            12d16s19
       ixw1=ibcoff                                                       7d26s22
       ixw2=ixw1+nspin                                                   7d26s22
       ibcoff=ixw2+nspin                                                 7d26s22
       call enough('prop.  4',bc,ibc)
       do ispin=iprop(1),iprop(2),2                                      7d26s22
        isad=(ispin-iprop(1))/2                                          7d27s22
        jxw1=ixw1+isad                                                   7d27s22
        jxw2=ixw2+isad                                                   7d27s22
        jcsf=icsf+icsfid*isad                                            7d27s22
        jcsf2=icsf2+icsfid*2*isad                                        7d27s22
        spin=0.5d0*dfloat(ispin-1)                                       2d24s21
        mdoou=(nec-ispin+1)/2                                            2d24s21
        mdooup=mdoou+1                                                   2d24s21
        ibcb4=ibcoff                                                     7d27s22
        call gencsf3(nec,spin,mdon,mdoou,norb,ibc(icsfpd),ispin,         2d24s21
     $     ibc(jcsf2),lprint,kxw1,kxw2,ibc(jcsf),0,ibcb4,               7d27s22
     $      icsfxv(mdooup),bc,ibc)                                      11d11s22
        ibc(jxw1)=kxw1                                                   7d27s22
        ibc(jxw2)=kxw2                                                   7d27s22
       end do                                                            7d26s22
       if(idorel.ne.0.and.idorel.ne.10)then                             4d1s24
        do isb=1,nsymb                                                   9d21s21
         do j=1,3                                                        9d21s21
          if(isym(j,isb).eq.1)then                                       9d21s21
           negs(j)=0                                                     9d21s21
          else                                                           9d21s21
           negs(j)=1                                                     9d21s21
          end if                                                         9d21s21
         end do                                                          9d21s21
         do iopt=1,4                                                      9d21s21
          if(idorel.lt.0.or.iopt.gt.1)then                                9d21s21
           if(iopt.eq.1)then                                              9d21s21
            isopt(1,iopt)=1                                              9d21s21
           else if(iopt.eq.2)then                                         9d21s21
            if(negs(1)*isym(1,isb).eq.-negs(1).and.negs(2)*isym(2,isb).  9d21s21
     $          eq.-negs(2).and.negs(3)*isym(3,isb).eq.negs(3))then     9d21s21
             isopt(1,iopt)=isb                                           9d21s21
            end if                                                       9d21s21
           else if(iopt.eq.3)then                                         9d21s21
            if(negs(1)*isym(1,isb).eq.-negs(1).and.negs(2)*isym(2,isb).  9d21s21
     $          eq.negs(2).and.negs(3)*isym(3,isb).eq.-negs(3))then     9d21s21
             isopt(1,iopt)=isb                                           9d21s21
            end if                                                       9d21s21
           else                                                           9d21s21
            if(negs(1)*isym(1,isb).eq.negs(1).and.negs(2)*isym(2,isb).   9d21s21
     $          eq.-negs(2).and.negs(3)*isym(3,isb).eq.-negs(3))then     9d21s21
             isopt(1,iopt)=isb                                           9d21s21
            end if                                                       9d21s21
           end if                                                         9d21s21
          end if                                                          9d21s21
         end do                                                          9d21s21
        end do                                                           9d21s21
        nsopt=0                                                           9d21s21
        iheadr=0                                                        1d24s23
        do iopt=1,4                                                      9d21s21
         if(idorel.lt.0.or.iopt.gt.1)then                                9d21s21
          if(iopt.eq.1.or.iopt.eq.3)then                                 9d21s21
           fname(1:3)='Re '                                              9d21s21
           irori=0                                                       9d21s21
          else                                                           9d21s21
           fname(1:3)='Im '                                              9d21s21
           irori=1                                                       9d21s21
          end if                                                         9d21s21
          if(iopt.le.2)then                                              9d21s21
           if(idorel.gt.0)then                                          3d29s24
            fname(4:8)='  0  '                                          3d29s24
           else                                                         3d29s24
            fname(4:8)='0,+-2'                                          3d29s24
           end if                                                       3d29s24
           idelms=0                                                      9d21s21
          else                                                           9d21s21
           fname(4:8)='+-1  '                                           3d29s24
           idelms=1                                                      9d21s21
          end if                                                         9d21s21
          if(iopt.eq.1.or.iopt.eq.4)then                                 9d21s21
           fname(9:9)='+'                                               3d29s24
           ispini=0                                                      9d21s21
          else                                                           9d21s21
           fname(9:9)='-'                                               3d29s24
           ispini=1                                                      9d21s21
          end if                                                         9d21s21
          if(idorel.eq.-4)then                                           9d21s21
           fname(12:12)='G'                                             3d29s24
           if(lprint.and.iheadr.eq.0)                                   1d24s23
     $          write(6,*)('>Coulomb-Gaunt Spin-orbit; 4v ints are '),  2d8s23
     $          onoff                                                   2d8s23
           iheadr=1                                                     1d24s23
          else if(idorel.gt.0)then                                       9d21s21
           fname(12:12)='C'                                             3d29s24
           if(lprint.and.iheadr.eq.0)write(6,*)('>Coulomb Spin-orbit;'),2d8s23
     $          ('4v ints are '),onoff                                  2d8s23
           iheadr=1                                                     1d24s23
          else                                                           9d21s21
           fname(12:12)='B'                                             3d29s24
           if(lprint.and.iheadr.eq.0)                                   1d24s23
     $          write(6,*)('>Coulomb-Breit Spin-orbit; 4v ints are '),  2d8s23
     $          onoff                                                   2d8s23
           iheadr=1                                                     1d24s23
          end if                                                         9d21s21
          if(iopt.eq.1)then                                              9d21s21
           if(idorel.eq.-4)then                                          9d21s21
            fname(10:12)=' G '                                          3d29s24
           else                                                          9d21s21
            fname(10:12)=' B '                                          3d29s24
           end if                                                        9d21s21
          else if(iopt.eq.2)then                                         9d21s21
           fname(10:11)='Lz'                                            3d29s24
          else if(iopt.eq.3)then                                         9d21s21
           fname(10:11)='Ly'                                            3d29s24
          else                                                           9d21s21
           fname(10:11)='Lx'                                            3d29s24
          end if                                                         9d21s21
          if(idorel.lt.0)then                                            3d24s22
           ngot=4                                                        3d24s22
          else                                                           3d24s22
           ngot=3                                                        3d24s22
          end if                                                         3d24s22
          if(ngot.gt.0)then                                              9d21s21
           if(nsopt.eq.0.and.lprint)write(6,1849)                        3d2s22
 1849      format(/,'spin-orbit operators',/,4x,'R or I',2x,'DMs',2x,    9d21s21
     $         'MsI',2x,'Op',2x,'sym')                                  9d21s21
           nsopt=nsopt+1                                                   9d21s21
           isopt(1,nsopt)=isopt(1,iopt)                                  9d21s21
           isopt(2,nsopt)=irori                                           9d21s21
           isopt(3,nsopt)=idelms                                          9d21s21
           isopt(4,nsopt)=ispini                                          9d21s21
           if(lprint)write(6,1848)nsopt,fname(1:3),fname(4:8),          3d29s24
     $          fname(9:9),fname(10:12),isopt(1,nsopt)                  3d29s24
          end if                                                         9d21s21
 1848     format(i2,2x,a3,4x,a5,2x,a1,3x,a3,1x,i2)                      3d29s24
         end if                                                          9d21s21
        end do                                                           9d21s21
        call parap42(natom,ngaus,ibdat,ixmtso,isym,iapair,ibstor,isstor, 7d27s22
     $      iso,0,ascale,ipt,npt,opdata,iopdata,iosym,nop,multh,        7d26s22
     $     nbb,nbasp,nbasdws,iorb,opname,isopt,nsopt,id,nopso,          2d17s22
     $      iopso,idoubo,bc,ibc)                                        11d9s22
        call parah042c(natom,ngaus,ibdat,ih0,ih0i,iovr,ibstor,isstor,    3d19s20
     $     0,ascale,nbasp,iorb,iapair,multh,ih0n,isopt,nsopt,srh,       7d26s22
     $      isym,bc,ibc,smsz)                                           11d16s23
        iifmx=ibcoff                                                     10d8s21
        ibcoff=iifmx+32*4                                                10d8s21
        isnd=ibcoff                                                      3d6s20
        nsnd=isnd+mynprocg                                               3d6s20
        ircv=nsnd+mynprocg                                               3d6s20
        nrcv=ircv+mynprocg                                               3d6s20
        ihrow=nrcv+mynprocg                                              3d6s20
        ibcoff=ihrow+mynprocg*64                                         9d24s21
        call enough('prop.  5',bc,ibc)
        call second(timecf0)                                             10d8s21
        call paraeri42cf(natom,ngaus,ibdat,ibstor,isstor,idwsdeb,         3d5s20
     $       ascale,nbasp,iorb,iapair,idoubo,irefo,idorel,ibc(isnd),    9d22s21
     $      ibc(nsnd),ibc(ircv),ibc(nrcv),ibc(ihrow),multh,i4or,ionexr, 10d5s21
     $      jmatsr,kmatsr,kmatsrb,i3xr,                                 11d15s21
     $      shift,ih0n,isopt,nsopt,scals,srh,isym,ibc(iifmx),ntype,nh0, 10d20s21
     $      i4o,bc,ibc)                                                 12d9s22
        call second(timecf1)                                             10d8s21
        telap=timecf1-timecf0
        if(lprint)write(6,*)('time for paraeri42cf'),telap
        call relgammas(nec,mdon,mdoo,irw0,irw1,irw2,spinx,bc,ibc)       11d9s22
       end if                                                            7d26s22
c                                                                       6d17s21
c     npadddi is number of extra words require for ddi_ixxx buffer      6d17s21
c                                                                       6d17s21
       npadddi=4*mynprocg+mynprocg                                               6d17s21
       imasspa=ibcoff                                                   8d16s22
       itmom=imasspa+nfname                                             8d16s22
       ibcoff=itmom+4*nfname*nfname                                      5d25s21
       call enough('prop.  6',bc,ibc)
       do i=imasspa,ibcoff-1                                            8d16s22
        ibc(i)=0                                                         5d25s21
       end do                                                            5d25s21
       jfname=ifname                                                     7d8s21
       ispin=-1                                                          7d26s22
       do jfn=1,nfname                                                   7d26s22
        mddilow=99999999                                                12d13s22
        mddihig=-1                                                      12d13s22
        ibctop=ibcoff                                                    7d26s22
        if(mynowprog.eq.0)then                                          11d9s22
         nhere=ibc(jfname)                                                11d22s19
         jfname=jfname+1                                                  11d22s19
         if(nhere.gt.200)then                                             11d22s19
          write(6,*)('name is tooooo long!')
          stop 'prop'
         end if                                                           11d22s19
         do j=1,nhere
          fname(j:j)=char(ibc(jfname))                                    11d22s19
          jfname=jfname+1
         end do                                                           11d22s19
         write(6,*)('file name: '),fname(1:nhere)                       11d9s22
        end if
        iout=ibcoff                                                      11d24s19
        call getwf(fname(1:nhere),ismult,isymmrci,nroot,iout,            2d11s22
     $       iwavedat,mdoox,mxopn,iorb,0,noff,nspc,nwavedat,            7d26s22
     $       elowest,idorel,nsymb,nbasp,nbasdws,nvirt,multh,maxbx,      7d21s21
     $      maxbxd,1,bc,ibc,mddilow,mddihig,itargs,itargd)              12d21s22
        maxbxi=maxbx+npadddi                                            7d30s22
        maxbxdi=maxbxd+npadddi                                          7d30s22
        call dws_synca                                                    7d9s21
        isad=(iwavedat(1,1)-iprop(1))/2                                  7d27s22
        jxw1=ixw1+isad                                                   7d27s22
        jxw2=ixw2+isad                                                   7d27s22
        kxw1=ibc(jxw1)                                                   7d27s22
        kxw2=ibc(jxw2)                                                   7d27s22
        jcsf=icsf+icsfid*isad                                            7d27s22
        mdoou=(nec-iwavedat(1,1)+1)/2                                    7d27s22
        mdooup=mdoou+1                                                   2d24s21
        jmasspa=imasspa                                                 8d16s22
        call c2eprop(iwavedat,kxw1,kxw2,ibc(jcsf),nec,nspc,ixmt,         7d27s22
     $       opdata,iopdata,iosym,nop,opname,nbasdws,nsymb,mdon,mdooup, 5d13s21
     $       ism,irel,irefo,norb,multh,idoubo,i2eop,ibc(jmasspa),natom, 8d16s22
     $         lfn,nvirt,maxbxi,maxbxdi,srh,sr2,npadddi,xmasstot,bc,ibc)11d10s22
        npack4=iwavedat(6,1)                                             7d27s22
        iwnxt=2                                                          7d27s22
        if(ipack1(2).gt.0)then                                           7d27s22
         if(ipack1(3).ne.0)then
          if(ipack1(2).eq.2)then                                         7d27s22
           iwnxt=3                                                       7d27s22
          else                                                           7d27s22
           iwnxt=2+2*ipack1(3)                                           7d27s22
          end if                                                         7d27s22
          call genwf(iwavedat,kxw1,kxw2,ibc(jcsf),nec,nspc,ixmt,         7d27s22
     $       opdata,iopdata,iosym,nop,opname,nbasdws,nsymb,mdon,mdooup, 5d13s21
     $       ism,irel,irefo,norb,multh,idoubo,isadd,nvirt,maxbxi,       7d30s22
     $        maxbxdi,srh,sr2,npadddi,lprint,genwfeps,bc,ibc)           11d10s22
         end if                                                          7d27s22
        end if                                                           7d27s22
        call gtmom(iwavedat,iwavedat,kxw1,kxw2,ibc(jcsf),nec,nspc,ixmt,  7d27s22
     $       opdata,iopdata,iosym,nop,opname,                           7d27s22
     $        nbasdws,nsymb,mdon,mdooup,ism,irel,irefo,norb,multh,      5d24s21
     $        idoubo,ibc(itmom),nvirt,maxbxi,maxbxdi,srh,sr2,npadddi,   7d30s22
     $        lprint,natom,bc,ibc)                                      11d10s22
        if(idorel.ne.0.and.idorel.ne.10)then                            4d1s24
         if(iwavedat(1,1).ne.1.or.idorel.lt.0)then                       7d27s22
          if(ipack1(2).eq.2)then                                         7d30s22
           call solin(iwavedat,iwavedat,nspc,nec,multh,                    7d27s22
     $         irefo,ih0a,i4o,irel,ism,norb,mdon,nvirt,maxbxi,          7d30s22
     $        maxbxdi,srh,sr2,nsymb,irw0,irw1,irw2,ih0n,nh0,isopt,nsopt,7d30s22
     $         i4or,ionexr,jmatsr,kmatsr,kmatsrb,i3xr,ibc(iifmx),ntype, 11d15s21
     $         npadddi,nbasp,nbaspc,natom,ngaus,ibdat,iapair,ibstor,    12d20s20
     $         isstor,isym,ascale,idorel,iorb,lprint,bc,ibc,shift,n4vso)2d8s23
          else                                                           7d30s22
           call sogen(iwavedat,iwavedat,nspc,nec,multh,                  7d30s22
     $         irefo,ih0a,i4o,irel,ism,norb,mdon,nvirt,maxbxi,           10d22s21
     $        maxbxdi,srh,sr2,nsymb,irw0,irw1,irw2,ih0n,nh0,isopt,nsopt,10d13s21
     $         i4or,ionexr,jmatsr,kmatsr,kmatsrb,i3xr,ibc(iifmx),ntype, 11d15s21
     $         npadddi,nbasp,nbaspc,natom,ngaus,ibdat,iapair,ibstor,    12d20s20
     $         isstor,isym,ascale,idorel,iorb,lprint,bc,ibc,shift,n4vso)2d8s23
          end if                                                         7d30s22
         end if                                                          7d27s22
         if(iwavedat(1,1).ne.1)then                                      7d27s22
          if(ntmso.ne.0)then                                            8d16s22
           call gtmomso(iwavedat,idum,iwavedat,                            7d27s22
     $         idum,idum,idum,nspc,nec,multh,irefo,                     3d8s22
     $     ixmtso,nh0,nopso,iopso,irel,ism,norb,mdon,dum,               12d1s22
     $         nvirt,maxbxi,maxbxdi,srh,sr2,nsymb,irw0,irw1,irw2,       7d30s22
     $         isopt,nsopt,                                             2d17s22
     $         ibc(iifmx),ntype,npadddi,nbasp,nbaspc,natom,ngaus,ibdat, 12d20s20
     $         iapair,ibstor,isstor,isym,ascale,idorel,iorb,opname,0,0, 3d1s22
     $        idum,idum,dum,1,idum,idum,1,dum,lprint,l2sub,bc,ibc,n4vso)2d8s23
          end if                                                        8d16s22
         end if                                                          7d27s22
        end if                                                           7d27s22
        kfname=jfname                                                    7d27s22
        do kfn=jfn+1,nfname                                              7d27s22
         nddilow=99999999                                                12d13s22
         nddihig=-1                                                      12d13s22
         ibctopi=ibcoff                                                    7d26s22
         if(mynowprog.eq.0)then                                         11d9s22
          nhere=ibc(kfname)                                                11d22s19
          kfname=kfname+1                                                  11d22s19
          if(nhere.gt.200)then                                             11d22s19
           write(6,*)('name is tooooo long!')
           stop 'prop'
          end if                                                           11d22s19
          do j=1,nhere
           fname(j:j)=char(ibc(kfname))                                    11d22s19
           kfname=kfname+1
          end do                                                           11d22s19
          write(6,*)('file name: '),fname(1:nhere)                      11d9s22
         end if                                                         11d9s22
         iout=ibcoff                                                      11d24s19
         call getwf(fname(1:nhere),ismult,isymmrci,nroot,iout,            2d11s22
     $       iwavedat(1,iwnxt),mdoox,mxopn,iorb,0,noff,nspc,nwavedat,   7d27s22
     $       elowest,idorel,nsymb,nbasp,nbasdws,nvirt,multh,maxbx,      7d21s21
     $      maxbxd,1,bc,ibc,nddilow,nddihig,itargs,itargd)              12d21s22
         maxbxi=maxbx+npadddi                                            7d30s22
         maxbxdi=maxbxd+npadddi                                          7d30s22
         call dws_synca                                                    7d9s21
         isdelta=iabs(iwavedat(1,iwnxt)-iwavedat(1,1))                   7d27s22
         npack4=iwavedat(6,iwnxt)                                        7d30s22
         if(isdelta.le.2.or.idorel.lt.0)then                             7d27s22
          if(ipack1(2).gt.0)then                                          7d27s22
           if(ipack1(3).ne.0)then                                          7d27s22
            isad=(iwavedat(1,iwnxt)-iprop(1))/2                             7d27s22
            jxw1=ixw1+isad                                                   7d27s22
            jxw2=ixw2+isad                                                   7d27s22
            kxw1=ibc(jxw1)                                                   7d27s22
            kxw2=ibc(jxw2)                                                   7d27s22
            jcsf=icsf+icsfid*isad                                            7d27s22
            mdoou=(nec-iwavedat(1,iwnxt)+1)/2                             7d27s22
            mdooup=mdoou+1                                                   2d24s21
            call genwf(iwavedat(1,iwnxt),kxw1,kxw2,ibc(jcsf),nec,nspc,    7d27s22
     $         ixmt,opdata,iopdata,iosym,nop,opname,nbasdws,nsymb,mdon, 7d27s22
     $         mdooup,ism,irel,irefo,norb,multh,idoubo,isadd,nvirt,     7d27s22
     $         maxbxi,maxbxdi,srh,sr2,npadddi,lprint,genwfeps,bc,ibc)   11d10s22
           end if                                                         7d27s22
          end if                                                          7d27s22
          if(iwavedat(1,iwnxt).eq.iwavedat(1,1))then                      7d27s22
           isad=(iwavedat(1,iwnxt)-iprop(1))/2                             7d27s22
           jxw1=ixw1+isad                                                   7d27s22
           jxw2=ixw2+isad                                                   7d27s22
           kxw1=ibc(jxw1)                                                   7d27s22
           kxw2=ibc(jxw2)                                                   7d27s22
           jcsf=icsf+icsfid*isad                                            7d27s22
           mdoou=(nec-iwavedat(1,iwnxt)+1)/2                              7d27s22
           mdooup=mdoou+1                                                   2d24s21
           call gtmom(iwavedat,iwavedat(1,iwnxt),kxw1,kxw2,ibc(jcsf),   7d30s22
     $          nec,nspc,ixmt,opdata,iopdata,iosym,nop,opname,          7d30s22
     $        nbasdws,nsymb,mdon,mdooup,ism,irel,irefo,norb,multh,      5d24s21
     $        idoubo,ibc(itmom),nvirt,maxbxi,maxbxdi,srh,sr2,npadddi,   7d30s22
     $        lprint,natom,bc,ibc)                                      11d10s22
          end if                                                          7d27s22
          if(idorel.ne.0.and.idorel.ne.10)then                          4d1s24
           if(isdelta.le.2.or.idorel.lt.0)then                           7d27s22
            if(ipack1(2).eq.2)then                                       7d30s22
             call solin(iwavedat,iwavedat(1,iwnxt),nspc,nec,multh,          7d27s22
     $         irefo,ih0a,i4o,irel,ism,norb,mdon,nvirt,maxbxi,          7d30s22
     $        maxbxdi,srh,sr2,nsymb,irw0,irw1,irw2,ih0n,nh0,isopt,nsopt,7d30s22
     $         i4or,ionexr,jmatsr,kmatsr,kmatsrb,i3xr,ibc(iifmx),ntype, 11d15s21
     $         npadddi,nbasp,nbaspc,natom,ngaus,ibdat,iapair,ibstor,    12d20s20
     $         isstor,isym,ascale,idorel,iorb,lprint,bc,ibc,shift,n4vso)2d8s23
            else                                                         7d30s22
             call sogen(iwavedat,iwavedat(1,iwnxt),nspc,nec,multh,       7d30s22
     $         irefo,ih0a,i4o,irel,ism,norb,mdon,nvirt,maxbxi,           10d22s21
     $        maxbxdi,srh,sr2,nsymb,irw0,irw1,irw2,ih0n,nh0,isopt,nsopt,10d13s21
     $         i4or,ionexr,jmatsr,kmatsr,kmatsrb,i3xr,ibc(iifmx),ntype, 11d15s21
     $         npadddi,nbasp,nbaspc,natom,ngaus,ibdat,iapair,ibstor,    12d20s20
     $         isstor,isym,ascale,idorel,iorb,lprint,bc,ibc,shift,n4vso)2d8s23
            end if                                                       7d30s22
           end if                                                        7d27s22
           if(isdelta.le.2)then                                          7d27s22
            if(ntmso.ne.0)then                                          8d16s22
             call gtmomso(iwavedat,idum,iwavedat(1,iwnxt),                  7d27s22
     $         idum,idum,idum,nspc,nec,multh,irefo,                     3d8s22
     $     ixmtso,nh0,nopso,iopso,irel,ism,norb,mdon,dum,               12d1s22
     $         nvirt,maxbxi,maxbxdi,srh,sr2,nsymb,irw0,irw1,irw2,       7d30s22
     $         isopt,nsopt,                                             2d17s22
     $         ibc(iifmx),ntype,npadddi,nbasp,nbaspc,natom,ngaus,ibdat, 12d20s20
     $         iapair,ibstor,isstor,isym,ascale,idorel,iorb,opname,0,0, 3d1s22
     $        idum,idum,dum,1,idum,idum,1,dum,lprint,l2sub,bc,ibc,n4vso)2d8s23
            end if                                                        8d16s22
           end if                                                        7d27s22
          end if                                                         7d27s22
         end if                                                          7d27s22
         ibcoff=ibctopi                                                  7d27s22
         if(nddihig.gt.nddilow)then                                     12d19s22
          do iddi=nddihig,nddilow,-1                                    12d19s22
           ipack8=iddi                                                  12d19s22
           call ddi_destroy(ipack8)                                     12d19s22
          end do                                                        12d19s22
          nddilow=99999999                                                12d13s22
          nddihig=-1                                                      12d13s22
         end if                                                         12d19s22
        end do                                                           7d27s22
        ibcoff=ibctop                                                    7d26s22
        if(mddihig.gt.mddilow)then                                      12d19s22
         do iddi=mddihig,mddilow,-1                                     12d19s22
          ipack8=iddi                                                   12d19s22
          call ddi_destroy(ipack8)                                      12d19s22
         end do                                                         12d19s22
        end if                                                          12d19s22
       end do                                                            7d26s22
       itmp=ibcoff                                                       5d4s22
       ibcoff=itmp+22                                                   3d2s23
       call enough('prop.  7',bc,ibc)
       jtmp=itmp-1                                                       5d4s22
       ktmp=jtmp+10                                                      5d4s22
       do i=1,10                                                         5d4s22
        bc(jtmp+i)=top(i)                                                5d4s22
        bc(ktmp+i)=top(i)**2                                             5d4s22
       end do                                                            5d4s22
       call dws_gsumf(bc(itmp),20)                                       5d4s22
       if(mynowprog.eq.0)then                                            5d4s22
        xxx=1d0/dfloat(mynprocg)                                          5d4s22
        do i=1,10                                                       11d17s22
         bc(jtmp+i)=bc(jtmp+i)*xxx                                        5d4s22
         bc(ktmp+i)=bc(ktmp+i)*xxx                                        5d4s22
         bc(ktmp+i)=sqrt(abs(bc(ktmp+i)-bc(jtmp+i)**2))                   5d4s22
         if(i.eq.1)then                                                 11d17s22
          code='ii'                                                     11d17s22
         else if(i.eq.2)then                                            11d17s22
          code='si'                                                     11d17s22
         else if(i.eq.3)then                                            11d17s22
          code='is'                                                     11d17s22
         else if(i.eq.4)then                                            11d17s22
          code='ss'                                                     11d17s22
         else if(i.eq.5)then                                            11d17s22
          code='di'                                                     11d17s22
         else if(i.eq.6)then                                            11d17s22
          code='id'                                                     11d17s22
         else if(i.eq.7)then                                            11d17s22
          code='dd'                                                     11d17s22
         else if(i.eq.8)then                                            11d17s22
          code='ds'                                                     11d17s22
         else if(i.eq.9)then                                            11d17s22
          code='sd'                                                     11d17s22
         else if(i.eq.10)then                                            11d17s22
          code='dv'
         end if                                                         11d17s22
         write(6,*)('psi*psi timer '),code,bc(jtmp+i),bc(ktmp+i)
        end do                                                            5d4s22
       end if                                                            5d4s22
       ktmp=jtmp+11                                                     3d2s23
       do i=1,11                                                        3d2s23
        bc(jtmp+i)=tso(i)                                               3d2s23
        bc(ktmp+i)=tso(i)**2                                            3d2s23
       end do                                                           3d2s23
       call dws_gsumf(bc(itmp),22)                                      3d2s23
       if(mynowprog.eq.0)then                                           3d2s23
        xxx=1d0/dfloat(mynprocg)                                        3d2s23
        do i=1,11                                                       3d2s23
         bc(jtmp+i)=bc(jtmp+i)*xxx                                      3d2s23
         bc(ktmp+i)=bc(ktmp+i)*xxx                                      3d2s23
         bc(ktmp+i)=sqrt(abs(bc(ktmp+i)-bc(jtmp+i)**2))                 3d2s23
         if(i.eq.1)then                                                 3d2s23
          code='ii'                                                     3d2s23
         else if(i.eq.2)then                                            3d2s23
          code='si'                                                     3d2s23
         else if(i.eq.3)then                                            3d2s23
          code='is'                                                     3d2s23
         else if(i.eq.4)then                                            3d2s23
          code='ss'                                                     3d2s23
         else if(i.eq.5)then                                            3d2s23
          code='di'                                                     3d2s23
         else if(i.eq.6)then                                            3d2s23
          code='id'                                                     3d2s23
         else if(i.eq.7)then                                            3d2s23
          code='dd'                                                     3d2s23
         else if(i.eq.8)then                                            3d2s23
          code='4v'                                                     3d2s23
         else if(i.eq.9)then                                            3d2s23
          code='ds'                                                     3d2s23
         else if(i.eq.10)then                                           3d2s23
          code='sd'                                                     3d2s23
         else if(i.eq.11)then                                           3d2s23
          code='dv'                                                     3d2s23
         end if                                                         3d2s23
         write(6,*)('psisopsi timer '),code,bc(jtmp+i),bc(ktmp+i)       3d2s23
        end do                                                          3d2s23
       end if                                                           3d2s23
       return                                                           7d30s22
      end if                                                            7d30s22
      xmasstot=0d0                                                      3d31s22
      maxbx=0                                                           7d9s21
      maxbxd=0                                                          7d21s21
      mdoox=0                                                           7d8s21
      mxopn=0                                                           7d8s21
      jfname=ifname                                                     7d8s21
      elowest=1d10                                                      7d8s21
      nhere=1                                                           10d2s24
      do i=1,nfname                                                     7d8s21
       if(mynowprog.eq.0)then                                           11d9s22
        nhere=ibc(jfname)                                                11d22s19
        jfname=jfname+1                                                  11d22s19
        if(nhere.gt.200)then                                             11d22s19
         write(6,*)('name is tooooo long!')
         stop 'prop'
        end if                                                           11d22s19
        do j=1,nhere
         fname(j:j)=char(ibc(jfname))                                    11d22s19
         jfname=jfname+1
        end do                                                           11d22s19
        write(6,*)('file name: '),fname(1:nhere)                        11d9s22
       end if                                                           11d9s22
       iout=ibcoff                                                      11d24s19
       call getwf(fname(1:nhere),ismult,isymmrci,nroot,iout,            2d11s22
     $       iwavedat(1,i+noff),mdoox,mxopn,iorb,i,noff,nspc,nwavedat,  6d1s21
     $       elowest,idorel,nsymb,nbasp,nbasdws,nvirt,multh,maxbx,      7d21s21
     $      maxbxd,1,bc,ibc,nddilow,nddihig,itargs,itargd)              12d21s22
      end do                                                            7d8s21
      npadddi=4*mynprocg+mynprocg                                               6d17s21
      maxbx=maxbx+npadddi                                               8d16s21
      maxbxd=maxbxd+npadddi                                             8d16s21
      call dws_synca                                                    7d9s21
      nfname=nfname+noff                                                5d12s21
      call dws_bcast(bc(iorb(1)),nwo)                                   12d16s19
      if(mod(natom,2).eq.0)then                                         5d27s21
       nsnd=(natom/2)+6                                                 5d27s21
      else                                                              5d27s21
       nsnd=(natom/2)+7                                                 5d27s21
      end if                                                            5d27s21
      call dws_bcast(cmx,nsnd)                                          5d27s21
      if(idorel.gt.0)then                                               3d2s22
       l2sub=1                                                          3d2s22
      else                                                              3d2s22
       l2sub=0                                                          3d2s22
      end if                                                            3d2s22
      mxopn=max(mxopn,2)
      if(lprint)write(6,*)('maximum no. of open shells = '),mxopn       12d20s19
      isl=iwavedat(1,1)                                                 12d12s19
      ish=iwavedat(1,1)                                                 12d12s19
      nec=iwavedat(7,1)                                                 12d13s19
      do i=2,nfname                                                     12d12s19
       isl=min(isl,iwavedat(1,i))                                       12d12s19
       ish=max(ish,iwavedat(1,i))                                       12d12s19
       if(iwavedat(7,i).ne.nec)then                                     12d13s19
        write(6,*)('number of electrons differ for filename '),i,
     $      (': want '),nec,('got '),iwavedat(7,i)                      12d13s19
        call dws_synca                                                  5d4s21
        call dws_finalize                                               12d13s19
        stop                                                            12d13s19
       end if                                                           12d13s19
      end do                                                            12d12s19
      mdon=(nec-mxopn)/2                                                12d17s19
      if(lprint)                                                        2d15s22
     $   write(6,*)('minimum no. of closed shells computed to be '),mdon2d15s22
      if(lprint)write(6,*)('range of spin multiplities: '),isl,ish      12d20s19
      spin=dfloat(ish-1)*0.5d0                                          12d12s19
      xms=dfloat(isl-1)*0.5d0                                           12d12s19
      if(lprint)write(6,*)('i.e. spin = '),spin                         2d15s22
      spinxx=dfloat(mxopn)*0.5d0                                          9d15s21
      if(lprint)
     $   write(6,*)('maximum possible spin given open shells = '),spinxx
      spinx=spin+2d0                                                    9d20s21
      if(lprint)write(6,*)('maximum possible spin given resi: '),spinx  2d15s22
      spin=min(spinx,spinxx)                                            9d20s21
      if(lprint)then                                                    12d20s19
       write(6,*)('high spin: '),spin                                    12d12s19
       write(6,*)('use ms = '),xms                                       12d12s19
      end if                                                            12d20s19
      spinx=spin                                                        9d20s21
      call stepwisecsf(mxopn,spin,spin,dum,idum,idum,0,idum,bc,ibc)     11d9s22
      icsfpd=ibcoff                                                     6d11s19
      mdoo=mdoox                                                        12d13s19
      if(mdoo.eq.mdon)then                                              3d2s22
       mdoo=max(mdoo,nec/2)                                             3d2s22
       if(lprint)write(6,*)('resetting mdoo to '),mdoo                  3d2s22
      end if                                                            3d2s22
      ibcoff=icsfpd+mdoo+1-mdon                                         6d11s19
      icsf2=ibcoff                                                      8d2s19
      icsf=icsf2+4*(mdoo+1-mdon)                                        12d14s19
      ibcoff=icsf+mdoo+1-mdon                                           12d14s19
      call enough('prop.  8',bc,ibc)
      do i=icsf,ibcoff-1                                                10d26s20
       ibc(i)=0                                                         10d26s20
      end do                                                            10d26s20
      call oplist(id,opname,ipt,npt,opdata,iopdata,iosym,ioprt,i2eop,   12d17s19
     $     nop,opnc,0)                                                  12d20s19
      nbb=0                                                             12d15s19
      if(idorel.eq.0)then                                               1d11s20
       ncomp=1                                                          1d11s20
      else                                                              1d11s20
       ncomp=2                                                          1d11s20
      end if                                                            1d11s20
      do i=1,nop                                                        12d15s19
       do isb=1,nsymb                                                   12d15s19
        isk=multh(isb,iosym(i))                                          12d15s19
        ixmt(isb,i)=ibcoff                                              12d15s19
        ibcoff=ixmt(isb,i)+nbasp(isb)*nbasp(isk)*ncomp*ncomp            1d11s20
        nbb=nbb+nbasp(isb)*nbasp(isk)*ncomp*ncomp                       1d11s20
       end do                                                           12d15s19
      end do                                                            12d15s19
      call enough('prop.  9',bc,ibc)
      do i=ixmt(1,1),ibcoff-1                                           12d15s19
       bc(i)=0d0                                                        12d15s19
      end do                                                            12d15s19
      idwsdeb=0                                                         12d16s19
      call parap(natom,ngaus,ibdat,ixmt,isym,iapair,ibstor,isstor,iso,  12d16s19
     $     idwsdeb,idorel,ascale,ipt,npt,opdata,iopdata,iosym,nop,multh,12d16s19
     $     nbb,nbasp,nbasdws,iorb,opname,1,bc,ibc)                      11d9s22
      idwsdeb=0
      jdenpt=ibcoff                                                     6d12s19
      ibcoff=jdenpt+idbk                                                6d12s19
      call make1dm(nsymb,irefo,nbasdws,.false.,maxddi,jmats,kmats,jmats,1d17s20
     $   ibc(jdenpt),isblk,isblkk,isblkh,nsdlk,nsdlkk,nsdlkh,multh,idbk,6d12s19
     $     isblk1,nsdlk1,isblkd,nsdlkd,0)                               7d5s21
       call second(timeq)
      ibcoff=jdenpt                                                     6d12s19
      do i=1,nsdlk                                                      6d30s18
       if(isblk(3,i).lt.isblk(4,i))then                                 7d5s18
        iswap=isblk(4,i)                                                6d30s18
        isblk(4,i)=isblk(3,i)                                           6d30s18
        isblk(3,i)=iswap                                                6d30s18
       end if                                                           6d30s18
      end do                                                            6d30s18
      do is3=1,nsymb                                                    12d16s19
       do is2=1,nsymb                                                   12d16s19
        is23=multh(is2,is3)                                             12d16s19
        do is1=1,nsymb                                                  12d16s19
         is4=multh(is1,is23)                                            12d16s19
         do i=1,nsdlk                                                   12d16s19
          if(isblk(1,i).eq.is1.and.isblk(2,i).eq.is2.and.               12d16s19
     $         isblk(3,i).eq.is3.and.isblk(4,i).eq.is4)then               12d16s19
           inv(1,is1,is2,is3)=i                                         12d16s19
           inv(2,is1,is2,is3)=1                                         12d16s19
           go to 2515                                                   12d16s19
          else if(isblk(2,i).eq.is1.and.isblk(1,i).eq.is2.and.          12d16s19
     $          isblk(4,i).eq.is3.and.isblk(3,i).eq.is4)then            12d16s19
           inv(1,is1,is2,is3)=i                                         12d16s19
           inv(2,is1,is2,is3)=2                                         12d16s19
           go to 2515                                                   12d16s19
          else if(isblk(1,i).eq.is1.and.isblk(2,i).eq.is2.and.          12d16s19
     $          isblk(4,i).eq.is3.and.isblk(3,i).eq.is4)then            12d16s19
           inv(1,is1,is2,is3)=i                                         12d16s19
           inv(2,is1,is2,is3)=3                                         12d16s19
           go to 2515                                                   12d16s19
          else if(isblk(2,i).eq.is1.and.isblk(1,i).eq.is2.and.          12d16s19
     $          isblk(3,i).eq.is3.and.isblk(4,i).eq.is4)then            12d16s19
           inv(1,is1,is2,is3)=i                                         12d16s19
           inv(2,is1,is2,is3)=4                                         12d16s19
           go to 2515                                                   12d16s19
          end if                                                        12d16s19
         end do                                                         12d16s19
 2515    continue                                                       12d16s19
        end do                                                          12d16s19
       end do                                                           12d16s19
      end do                                                            12d16s19
c
c     check for wavefunctions we have to generate via angular momentum  5d12s21
c     operators.                                                        5d12s21
c
      mdoop=mdoo+1                                                      5d12s21
      ilastsp=-1                                                        5d21s21
      if(isl.ne.ish)ilastsp=-1                                          5d12s21
      idiab=0                                                           2d11s22
      nno=0                                                             2d11s22
      do i=1,nfname                                                     2d11s22
       npack4=iwavedat(6,i)                                             2d11s22
       if(ipack1(4).ge.0)then                                           2d11s22
        if(i.eq.1)then                                                   2d11s22
         lam=ipack1(3)                                                   2d11s22
         ismult=iwavedat(1,i)                                            2d11s22
         isymmrci=iwavedat(2,i)                                          2d11s22
         nroot=iwavedat(3,i)                                             2d11s22
        else                                                             2d11s22
         if(lam.ne.ipack1(3).or.ismult.ne.iwavedat(1,i).or.              2d11s22
     $     isymmrci.ne.iwavedat(2,i).or.nroot.ne.iwavedat(3,i))nno=nno+12d11s22
        end if                                                           2d11s22
        idiab=idiab+ipack1(1)                                            2d11s22
       end if                                                           2d11s22
      end do                                                            2d11s22
      if(idiab.ne.0)then                                                2d11s22
       if(lprint)write(6,*)('diabatic transformation!'),idiab                           2d11s22
       if(nno.ne.0)then                                                 2d11s22
        write(6,*)('oopsy doopsy: nno = '),nno                          2d11s22
        write(6,*)('all files are not for the same wavefunctions')      2d11s22
        call dws_synca                                                   2d11s22
        call dws_finalize                                                2d11s22
        stop                                                             2d11s22
       end if                                                           2d11s22
       ispin=iwavedat(1,1)                                              2d11s22
       spin=0.5d0*dfloat(ispin-1)                                       2d24s21
       mdoou=(nec-ispin+1)/2                                            2d24s21
       mdooup=mdoou+1                                                   2d24s21
       ibcb4=ibcoff                                                     5d25s21
       call gencsf3(nec,spin,mdon,mdoou,norb,ibc(icsfpd),ispin,         2d24s21
     $     ibc(icsf2),lprint,ixw1,ixw2,ibc(icsf),0,ibcb4,icsfxv(mdooup),11d11s22
     $      bc,ibc)                                                     11d11s22
       nroot=iwavedat(3,1)                                              2d11s22
       idovr=ibcoff                                                      2d11s22
       ibcoff=idovr+nroot*nroot*idiab                                    2d11s22
       ibcmaxbx=ibcoff
       ibcoff=ibcmaxbx+1
       ipack4(1)=maxbx
       ibc(ibcmaxbx)=ipack8
       call enough('prop. 10',bc,ibc)
       jdovr=idovr                                                        2d11s22
       xinv=1d0/dfloat(nec)                                             2d11s22
       do i=2,nfname                                                    2d11s22
        npack4=iwavedat(6,i)                                            2d11s22
        if(ipack1(4).ge.0)then                                          2d11s22
         call dovr(iwavedat(1,1),iwavedat(1,i),bc(jdovr),mdon,mdooup,    2d11s22
     $       ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,nvirt,
     $       ibc(ibcmaxbx),
     $       maxbxd,srh,sr2,multh,nsymb,bc(icsf),npadddi,bc,ibc)        11d10s22
         if(lprint)write(6,*)('overlap with fcn '),i
         do j=0,nroot*nroot-1                                           2d11s22
          bc(jdovr+j)=bc(jdovr+j)*xinv                                  2d11s22
         end do                                                         2d11s22
         if(lprint)call prntm2(bc(jdovr),nroot,nroot,nroot)
         jdovr=jdovr+nroot*nroot                                          2d11s22
        end if                                                          2d11s22
       end do                                                           2d11s22
       iorg=ibcoff
       ibcoff=iorg+nroot*nroot
       ivecon=ibcoff                                                      2d11s22
       ibcoff=ivecon+nroot*nroot                                          2d11s22
       call enough('prop. 11',bc,ibc)
       do i=0,nroot-1                                                   2d11s22
        ii=ivecon+nroot*i                                                 2d11s22
        do j=0,nroot-1                                                  2d11s22
         bc(ii+j)=0d0                                                   2d11s22
        end do                                                          2d11s22
        bc(ii+i)=1d0                                                    2d11s22
       end do                                                           2d11s22
       do i=0,nroot*nroot-1                                             2d11s22
        bc(iorg+i)=bc(idovr+i)                                          2d11s22
       end do                                                           2d11s22
       sz=0d0                                                           2d11s22
       numb=0                                                           2d11s22
       do if=1,idiab                                                    2d11s22
        jdovr=idovr+nroot*nroot*(if-1)                                  2d11s22
        do i=0,nroot-2                                                  2d11s22
         do j=i+1,nroot-1                                               2d11s22
          ji=jdovr+j+nroot*i                                            2d11s22
          ij=jdovr+i+nroot*j                                            2d11s22
          sz=sz+bc(ji)**2+bc(ij)**2                                     2d11s22
          numb=numb+2                                                   2d11s22
         end do                                                         2d11s22
        end do                                                          2d11s22
       end do                                                           2d11s22
       sz=sqrt(sz/dfloat(numb))                                         2d11s22
       call dws_bcast(sz,1)                                             2d15s22
       if(lprint)write(6,*)('starting sz: '),sz                                   2d11s22
       szo=sz                                                           2d11s22
       do isweep=1,300                                                  2d11s22
        do i1=0,nroot-2                                                 2d11s22
         do i2=i1+1,nroot-1                                             2d11s22
          a=0d0                                                         2d11s22
          b=0d0                                                         2d11s22
          c=0d0                                                         2d11s22
          jouse=idovr                                                   2d11s22
          do if=1,idiab                                                 2d11s22
           i11=jouse+i1+nroot*i1                                        2d11s22
           i21=jouse+i2+nroot*i1                                        2d11s22
           i12=jouse+i1+nroot*i2                                        2d11s22
           i22=jouse+i2+nroot*i2                                        2d11s22
           a=a+bc(i11)*bc(i21)-bc(i12)*bc(i22)                          2d11s22
           b=b+bc(i11)**2+bc(i22)**2-bc(i12)**2-bc(i21)**2              2d11s22
           jouse=jouse+nroot*nroot                                      2d11s22
          end do                                                        2d11s22
          c=-a                                                          2d11s22
          dist2=b*b-4d0*a*c                                             2d11s22
          if(dist2.lt.0d0)then                                          2d11s22
           write(6,*)('discriminant is negative!!! '),dist2             2d11s22
           write(6,*)isweep,i1,i2,a,b                                   2d11s22
           stop                                                         2d11s22
          end if                                                        2d11s22
          dist=sqrt(abs(dist2))                                         2d11s22
          q=-0.5d0*(b+sign(dist,b))                                     2d11s22
          t1=q/a                                                        2d11s22
          t2=c/q                                                        2d11s22
          t12=t1**2                                                     2d11s22
          s12=t12/(1d0+t12)                                             2d11s22
          c1=sqrt(1d0-s12)                                              2d11s22
          s1=sign(sqrt(s12),t1)                                         2d11s22
          t22=t2**2                                                     2d11s22
          s22=t22/(1d0+t22)                                             2d11s22
          c2=sqrt(1d0-s22)                                              2d11s22
          s2=sign(sqrt(s22),t2)                                         2d11s22
          phi1=0d0                                                      2d11s22
          phi2=0d0                                                      2d11s22
          jouse=idovr                                                   2d11s22
          do if=1,idiab                                                 2d11s22
           i11=jouse+i1+i1*nroot                                        2d11s22
           i21=jouse+i2+i1*nroot                                        2d11s22
           i12=jouse+i1+i2*nroot                                        2d11s22
           i22=jouse+i2+i2*nroot                                        2d11s22
           phi1=phi1+(-s1*bc(i11)+c1*bc(i21))**2                        2d11s22
     $          +(c1*bc(i12)+s1*bc(i22))**2                             2d11s22
           phi2=phi2+(-s2*bc(i11)+c2*bc(i21))**2                        2d11s22
     $          +(c2*bc(i12)+s2*bc(i22))**2                             2d11s22
           jouse=jouse+nroot*nroot                                      2d11s22
          end do                                                        2d11s22
          if(phi1.lt.phi2)then                                          2d11s22
           cu=c1                                                        2d11s22
           su=s1                                                        2d11s22
          else                                                          2d11s22
           cu=c2                                                        2d11s22
           su=s2                                                        2d11s22
          end if                                                        2d11s22
          jouse=idovr                                                   2d11s22
          do if=1,idiab                                                 2d11s22
           do n=0,nroot-1                                               2d11s22
            i1n=jouse+i1+nroot*n                                        2d11s22
            i2n=jouse+i2+nroot*n                                        2d11s22
            x1=cu*bc(i1n)+su*bc(i2n)                                    2d11s22
            y1=-su*bc(i1n)+cu*bc(i2n)                                   2d11s22
            bc(i1n)=x1                                                  2d11s22
            bc(i2n)=y1                                                  2d11s22
           end do                                                       2d11s22
           jouse=jouse+nroot*nroot                                      2d11s22
          end do                                                        2d11s22
          do n=0,nroot-1                                                2d11s22
           in1=ivecon+n+nroot*i1                                          2d11s22
           in2=ivecon+n+nroot*i2                                          2d11s22
           x1=cu*bc(in1)+su*bc(in2)                                     2d11s22
           y1=-su*bc(in1)+cu*bc(in2)                                    2d11s22
           bc(in1)=x1                                                   2d11s22
           bc(in2)=y1                                                   2d11s22
          end do                                                        2d11s22
         end do                                                         2d11s22
        end do                                                          2d11s22
        sz=0d0                                                          2d11s22
        jouse=idovr                                                     2d11s22
        do if=1,idiab                                                   2d11s22
         do i1=0,nroot-2                                                2d11s22
          do i2=i1+1,nroot-1                                            2d11s22
           i12=jouse+i1+nroot*i2                                        2d11s22
           i21=jouse+i2+nroot*i1                                        2d11s22
           sz=sz+bc(i12)**2+bc(i21)**2                                  2d11s22
          end do                                                        2d11s22
         end do                                                         2d11s22
         jouse=jouse+nroot*nroot                                        2d11s22
        end do                                                          2d11s22
        sz=sqrt(sz/dfloat(numb))                                        2d11s22
        call dws_bcast(sz,1)                                            2d15s22
        if(lprint)write(6,*)('at end of sweep we have '),sz                       2d11s22
        if(sz/szo.gt.0.9999d0)go to 2258                                2d11s22
        szo=sz                                                          2d11s22
       end do                                                           2d11s22
 2258  continue                                                         2d11s22
       det=1d0                                                          2d11s22
       do i1=0,nroot-1                                                  2d11s22
        i11=idovr+i1+nroot*i1                                           2d11s22
        if(bc(i11).lt.0d0)then                                          2d11s22
         do n=0,nroot-1                                                  2d11s22
          in1=ivecon+n+nroot*i1                                            2d11s22
          bc(in1)=-bc(in1)                                               2d11s22
         end do                                                          2d11s22
        end if                                                           2d11s22
        det=det*abs(bc(i11))                                             2d11s22
       end do                                                            2d11s22
       if(lprint)then                                                   2d15s22
        write(6,*)('product of diagonals: '),det                          2d11s22
        write(6,*)('final rotation matrix: ')
        call prntm2(bc(ivecon),nroot,nroot,nroot)
       end if                                                           2d15s22
       call dws_bcast(bc(ivecon),nroot*nroot)                           2d13s22
       ivecno=ibcoff                                                     2d11s22
       ibcoff=ivecno+nroot*nroot                                         2d11s22
       call enough('prop. 12',bc,ibc)
       do i=0,nroot-1                                                    2d11s22
        do j=0,nroot-1                                                   2d11s22
         ji=ivecno+j+nroot*i                                             2d11s22
         ij=ivecon+i+nroot*j                                             2d11s22
         bc(ji)=bc(ij)                                                   2d11s22
        end do                                                           2d11s22
       end do                                                            2d11s22
       if(lprint)then
        write(6,*)('final overlap with first fcn: ')
        call prntm2(bc(idovr),nroot,nroot,nroot)
       end if                                                           2d15s22
       ipack8=ibc(iwavedat(4,1)+iwavedat(13,1))                          2d11s22
       nctf=ipack4(1)                                                    2d11s22
       ieg=iwavedat(13,1)+nroot*nctf+iwavedat(4,1)+1
       do io=0,nroot-1
        ee=bc(ieg+io)                                                     2d11s22
        jdovr=idovr+nroot*io                                             2d11s22
        jvecno=ivecno+nroot*io                                           2d11s22
        do in=0,nroot-1
         bc(jdovr+in)=bc(jvecno+in)*ee                                   2d11s22
        end do                                                           2d11s22
       end do                                                            2d11s22
       ihdd=ibcoff                                                      2d14s22
       ibcoff=ihdd+nroot                                                2d14s22
       call dgemm('n','n',nroot,nroot,nroot,1d0,bc(idovr),nroot,         2d11s22
     $     bc(ivecon),nroot,0d0,bc(ibcoff),nroot,                       2d11s22
     d' prop.  1')
       if(lprint)then                                                   2d15s22
        write(6,*)('eigenvectors')                                        2d11s22
        do i=0,nroot-1                                                    2d11s22
         iv=ivecon+nroot*i                                                2d11s22
         write(6,442)i,(bc(iv+j),j=0,nroot-1)                             2d11s22
  442    format('<d',i2,'|a>',50es22.14)                                  2d11s22
        end do                                                            2d11s22
        write(6,*)('diabatic energies ')                                  2d11s22
       end if                                                           2d15s22
       do i=0,nroot-1                                                    2d11s22
        iham=ibcoff+nroot*i                                              2d11s22
        if(lprint)write(6,441)i,(bc(iham+j),j=0,i)                      2d15s22
  441   format('<d',i2,'|Hd|d''>',50es22.14)                                   2d11s22
        bc(ihdd+i)=bc(iham+i)                                           2d14s22
       end do                                                            2d11s22
       call wrot(iwavedat,mdon,mdooup,nvirt,multh,nsymb,bc(ivecon),      2d13s22
     $     bc(ivecno),bc,ibc)                                           11d9s22
       jdovr=idovr                                                        2d11s22
       do i=2,nfname                                                    2d11s22
        npack4=iwavedat(6,i)                                            2d11s22
        if(ipack1(4).ge.0)then                                          2d11s22
         call dovr(iwavedat(1,1),iwavedat(1,i),bc(jdovr),mdon,mdooup,    2d11s22
     $       ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,nvirt,maxbx,  2d11s22
     $       maxbxd,srh,sr2,multh,nsymb,bc(icsf),npadddi,bc,ibc)        11d10s22
         if(lprint)then                                                 2d15s22
          write(6,*)('final overlap with fcn '),i
          do j=0,nroot*nroot-1                                           2d11s22
           bc(jdovr+j)=bc(jdovr+j)*xinv                                  2d11s22
          end do                                                         2d11s22
          do j=0,nroot-1                                                12d2s22
           iad=jdovr+j
           if(lprint)write(6,447)j,(bc(iad+k*nroot),k=0,nroot-1)        12d2s22
  447      format('<d',i2,'|ref>',50es22.14)                            12d2s22
          end do                                                        12d2s22
         end if                                                         2d15s22
         jdovr=jdovr+nroot*nroot                                          2d11s22
        end if                                                          2d11s22
       end do                                                           2d11s22
       jfname=ifname                                                     7d8s21
       if(mynowprog.eq.0)then                                           11d9s22
        nhere=ibc(jfname)                                                11d22s19
        jfname=jfname+1                                                  11d22s19
        if(nhere.gt.200)then                                             11d22s19
         write(6,*)('name is tooooo long!')
         stop 'prop'
        end if                                                           11d22s19
        fname(1:1)='d'                                                   2d13s22
        do j=1,nhere
         jp=j+1                                                          2d13s22
         fname(jp:jp)=char(ibc(jfname))                                    11d22s19
         jfname=jfname+1
        end do                                                           11d22s19
        nhere=nhere+1                                                    2d13s22
        if(lprint)write(6,*)('file name: '),fname(1:nhere)               2d15s22
       end if                                                           11d9s22
       call putwf(fname(1:nhere),iwavedat,idorel,nsymb,nbasp,idoubo,    2d14s22
     $      irefo,nbasdws,nvirt,irel,ism,multh,maxbx,maxbxd,iapair,     2d14s22
     $      bc(ihdd),mdon,mdooup,bc,ibc)                                11d9s22
       return                                                           1d30s23
      end if                                                            2d11s22
      imasspa=ibcoff                                                    8d16s22
      itmom=imasspa+nfname                                              8d16s22
      ibcoff=itmom+4*nfname*nfname                                      5d25s21
      call enough('prop. 13',bc,ibc)
      do i=imasspa,ibcoff-1                                             8d16s22
       ibc(i)=0                                                         5d25s21
      end do                                                            5d25s21
      do ispin=isl,ish,2                                                2d24s21
       spin=0.5d0*dfloat(ispin-1)                                       2d24s21
       mdoou=(nec-ispin+1)/2                                            2d24s21
       mdooup=mdoou+1                                                   2d24s21
       ibcb4=ibcoff                                                     5d25s21
       call gencsf3(nec,spin,mdon,mdoou,norb,ibc(icsfpd),ispin,         2d24s21
     $     ibc(icsf2),lprint,ixw1,ixw2,ibc(icsf),0,ibcb4,icsfxv(mdooup),11d11s22
     $      bc,ibc)                                                     11d11s22
       nok=0                                                            5d24s21
       iok=ibcoff                                                       5d24s21
       ibcoff=iok+nfname                                                5d25s21
       call enough('prop. 14',bc,ibc)
       is=1                                                             5d24s21
   21  continue                                                         5d24s21
        isadd=1                                                         5d24s21
        if(iwavedat(1,is).eq.ispin)then                                 5d24s21
         npack4=iwavedat(6,is)                                          5d24s21
         if(ipack1(4).lt.0)then                                         5d24s21
          im=is-1                                                          5d12s21
          call genwf(iwavedat(1,im),ixw1,ixw2,ibc(icsf),nec,nspc,ixmt,    5d12s21
     $       opdata,iopdata,iosym,nop,opname,nbasdws,nsymb,mdon,mdooup, 5d13s21
     $       ism,irel,irefo,norb,multh,idoubo,isadd,nvirt,maxbx,maxbxd, 7d27s21
     $         srh,sr2,npadddi,lprint,genwfeps,bc,ibc)                  11d10s22
         else                                                           5d24s21
          jmasspa=imasspa+is-1                                          8d16s22
          if(iprop(7).eq.0)then                                         11d2s22
           call c2eprop(iwavedat(1,is),ixw1,ixw2,ibc(icsf),nec,nspc,    11d16s22
     $         ixmt,opdata,iopdata,iosym,nop,opname,nbasdws,nsymb,mdon, 11d16s22
     $         mdooup,ism,irel,irefo,norb,multh,idoubo,i2eop,           11d16s22
     $         ibc(jmasspa),natom,lfn,nvirt,maxbx,maxbxd,srh,sr2,       11d16s22
     $         npadddi,xmasstot,bc,ibc)                                 11d16s22
          ibc(iok+nok)=is                                               5d24s21
          end if                                                        11d2s22
          nok=nok+1                                                     5d24s21
         end if                                                           5d12s21
        end if                                                          2d24s21
        is=is+isadd                                                     5d24s21
       if(is.le.nfname)go to 21                                         5d24s21
       if(iprop(7).ne.0)then                                            11d2s22
        write(6,*)('computing natural orbitals ...'),nfname                    11d2s22
        ibss=ibcoff                                                     11d2s22
        nbss=0                                                          1d27s23
        do isb=1,nsymb                                                  11d2s22
         nh0(isb)=nbasdws(isb)-idoubo(isb)                              4d25s23
         iden(isb)=ibcoff                                               11d2s22
         ibcoff=iden(isb)+nbasdws(isb)**2                               4d25s23
         write(6,*)('isb '),isb,nbasdws(isb),idoubo(isb),nh0(isb)       4d25s23
         nbss=nbss+nbasdws(isb)**2                                      4d25s23
        end do                                                          11d2s22
        do isb=1,nsymb                                                  1d30s23
         ident(isb)=ibcoff                                              1d30s23
         ibcoff=ident(isb)+nh0(isb)**2                                  4d25s23
        end do
        call enough('prop.noa',bc,ibc)
        do iz=ibss,ibcoff-1                                             11d2s22
         bc(iz)=0d0                                                     11d2s22
        end do                                                          11d2s22
        wtt=0d0                                                         11d2s22
        do iwf=1,nfname                                                 11d2s22
         isymmrci=iwavedat(2,iwf)                                       1d30s23
         call bdens(iden,iwavedat(1,iwf),ixw1,ixw2,ibc(icsf),nec,       11d2s22
     $        nbasdws,nsymb,mdon,mdooup,ism,irel,irefo,norb,multh,      11d2s22
     $        idoubo,maxbx,maxbxd,srh,sr2,npadddi,lprint,wtt,bc,ibc,    1d27s23
     $        nvirt,ident,nh0,nsdlk,isblk)                                          4d25s23
        end do                                                          11d2s22
        call dws_gsumf(bc(ibss),nbss)                                   1d27s23
        write(6,*)('total weight: '),wtt                                11d2s22
        wtti=1d0/wtt                                                    11d2s22
        write(6,*)('scale densities by '),wtti
        trac=0d0                                                        1d27s23
        do isb=1,nsymb                                                  11d2s22
         if(nbasdws(isb).gt.0)then                                      11d2s22
          write(6,*)('starting density for symmetry block '),isb,
     $        iden(isb)
          call prntm2(bc(iden(isb)),nh0(isb),nh0(isb),nh0(isb))         4d25s23
          do i=1,nh0(isb)-1                                             4d25s23
           do j=0,i-1                                                   1d30s23
            ji=iden(isb)+j+nh0(isb)*i                                   4d25s23
            ij=iden(isb)+i+nh0(isb)*j                                   4d25s23
            bc(ji)=bc(ij)                                               1d30s23
           end do                                                       1d30s23
          end do                                                        1d30s23
          write(6,*)('symmetrized: ')
          call prntm2(bc(iden(isb)),nh0(isb),nh0(isb),nh0(isb))         4d25s23
          do i=0,nh0(isb)*nh0(isb)-1                                    4d25s23
           bc(iden(isb)+i)=-bc(iden(isb)+i)*wtti                        11d2s22
          end do                                                        11d2s22
          write(6,*)('-scaled density ')                                 11d2s22
          call prntm2(bc(iden(isb)),nh0(isb),nh0(isb),nh0(isb))         4d25s23
          ieig=ibcoff                                                   11d2s22
          ivec=ieig+nh0(isb)                                            4d25s23
          isymq=ivec+nh0(isb)*nh0(isb)                                  4d25s23
          ibcoff=isymq+nh0(isb)                                         4d25s23
          call enough('prop.nob',bc,ibc)
          call diagx(nh0(isb),bc(iden(isb)),bc(ieig),bc(ivec),          4d25s23
     $         ibc(isymq),bc,ibc)                                       11d14s22
          do i=0,nh0(isb)-1                                             4d25s23
           bc(ieig+i)=-bc(ieig+i)                                       11d2s22
           trac=trac+bc(ieig+i)
          end do                                                        11d2s22
          write(6,*)('no occupancies ')                                 11d2s22
          call prntm2(bc(ieig),1,nh0(isb),1)                            4d25s23
          write(6,*)('trace so far '),trac
          write(6,*)('no in nh0av basis ')                                 11d2s22
          call prntm2(bc(ivec),nh0(isb),nh0(isb),nh0(isb))              4d25s23
          do iz=0,nbasdws(isb)*nbasdws(isb)-1                           4d25s23
           bc(iden(isb)+iz)=0d0                                         4d25s23
          end do                                                        4d25s23
          do i=0,idoubo(isb)-1                                          4d25s23
           ii=iden(isb)+i*(nbasdws(isb)+1)                              4d25s23
           bc(ii)=1d0                                                   4d25s23
          end do                                                        4d25s23
          do i=0,nh0(isb)-1                                             4d25s23
           ip=i+idoubo(isb)                                             4d25s23
           iad1=iden(isb)+idoubo(isb)+nbasdws(isb)*ip                   4d25s23
           iad2=ivec+nh0(isb)*i                                         4d25s23
           do j=0,nh0(isb)-1                                            4d25s23
            bc(iad1+j)=bc(iad2+j)                                       4d25s23
           end do                                                       4d25s23
          end do                                                        4d25s23
          write(6,*)('no in nbasdws basis ')                            4d25s23
          call prntm2(bc(iden(isb)),nbasdws(isb),nbasdws(isb),          4d25s23
     $         nbasdws(isb))                                            4d25s23
          ibcoff=ieig                                                   11d2s22
         end if                                                         11d2s22
        end do                                                          11d2s22
        if(mynowprog.eq.0)then                                           4d25s23
         write(6,*)('writing these orbs to file norbs')                 4d25s23
         jfname=ifname                                                  4d25s23
         nhere=ibc(jfname)                                              4d25s23
         jfname=jfname+1                                                4d25s23
         do j=1,nhere
          fname(j:j)=char(ibc(jfname))                                    11d22s19
          jfname=jfname+1
         end do                                                           11d22s19
         write(6,*)('file name: '),fname(1:nhere)                        11d9s22
         call ntrans(iden,fname(1:nhere),bc,ibc)                        4d25s23
        end if                                                          4d25s23
        call dws_synca                                                  4d25s23
        ibcoff=ibcoffo                                                    11d24s19
        return
       end if                                                           11d2s22
c
       do iket=0,nok-1                                                  5d24s21
        ikwf=ibc(iok+iket)                                              5d24s21
        do ibra=0,iket                                                  5d24s21
         ibwf=ibc(iok+ibra)                                             5d24s21
         jtmom=itmom+4*(ibwf-1+nfname*(ikwf-1))                         5d25s21
         call gtmom(iwavedat(1,ibwf),iwavedat(1,ikwf),ixw1,ixw2,        5d24s21
     $        ibc(icsf),nec,nspc,ixmt,opdata,iopdata,iosym,nop,opname,  5d24s21
     $        nbasdws,nsymb,mdon,mdooup,ism,irel,irefo,norb,multh,      5d24s21
     $        idoubo,ibc(jtmom),nvirt,maxbx,maxbxd,srh,sr2,npadddi,     3d2s22
     $        lprint,natom,bc,ibc)                                      11d10s22
        end do                                                          5d24s21
       end do                                                           5d24s21
      end do                                                            2d24s21
      itmp=ibcoff                                                       5d4s22
      ibcoff=itmp+20                                                    5d4s22
      call enough('prop. 15',bc,ibc)
      jtmp=itmp-1                                                       5d4s22
      ktmp=jtmp+10                                                      5d4s22
      do i=1,10                                                         5d4s22
       bc(jtmp+i)=top(i)                                                5d4s22
       bc(ktmp+i)=top(i)**2                                             5d4s22
      end do                                                            5d4s22
      call dws_gsumf(bc(itmp),20)                                       5d4s22
      if(mynowprog.eq.0)then                                            5d4s22
       xxx=1d0/dfloat(mynprocg)                                          5d4s22
       do i=1,10                                                         5d4s22
        bc(jtmp+i)=bc(jtmp+i)*xxx                                        5d4s22
        bc(ktmp+i)=bc(ktmp+i)*xxx                                        5d4s22
        bc(ktmp+i)=sqrt(abs(bc(ktmp+i)-bc(jtmp+i)**2))                   5d4s22
         if(i.eq.1)then                                                 11d17s22
          code='ii'                                                     11d17s22
         else if(i.eq.2)then                                            11d17s22
          code='si'                                                     11d17s22
         else if(i.eq.3)then                                            11d17s22
          code='is'                                                     11d17s22
         else if(i.eq.4)then                                            11d17s22
          code='ss'                                                     11d17s22
         else if(i.eq.5)then                                            11d17s22
          code='di'                                                     11d17s22
         else if(i.eq.6)then                                            11d17s22
          code='id'                                                     11d17s22
         else if(i.eq.7)then                                            11d17s22
          code='dd'                                                     11d17s22
         else if(i.eq.8)then                                            11d17s22
          code='ds'                                                     11d17s22
         else if(i.eq.9)then                                            11d17s22
          code='sd'                                                     11d17s22
         else if(i.eq.10)then                                            11d17s22
          code='dv'
         end if                                                         11d17s22
        write(6,*)('psioppsi timer '),code,bc(jtmp+i),bc(ktmp+i)
       end do                                                            5d4s22
      end if                                                            5d4s22
      ibcoff=itmp                                                       5d4s22
      if(idorel.ne.0.and.idorel.ne.10)then                              4d1s24
       do isb=1,nsymb                                                   9d21s21
        do j=1,3                                                        9d21s21
         if(isym(j,isb).eq.1)then                                       9d21s21
          negs(j)=0                                                     9d21s21
         else                                                           9d21s21
          negs(j)=1                                                     9d21s21
         end if                                                         9d21s21
        end do                                                          9d21s21
        do iopt=1,4                                                      9d21s21
         if(idorel.lt.0.or.iopt.gt.1)then                                9d21s21
          if(iopt.eq.1)then                                              9d21s21
           isopt(1,iopt)=1                                              9d21s21
          else if(iopt.eq.2)then                                         9d21s21
           if(negs(1)*isym(1,isb).eq.-negs(1).and.negs(2)*isym(2,isb).  9d21s21
     $          eq.-negs(2).and.negs(3)*isym(3,isb).eq.negs(3))then     9d21s21
            isopt(1,iopt)=isb                                           9d21s21
           end if                                                       9d21s21
          else if(iopt.eq.3)then                                         9d21s21
           if(negs(1)*isym(1,isb).eq.-negs(1).and.negs(2)*isym(2,isb).  9d21s21
     $          eq.negs(2).and.negs(3)*isym(3,isb).eq.-negs(3))then     9d21s21
            isopt(1,iopt)=isb                                           9d21s21
           end if                                                       9d21s21
          else                                                           9d21s21
           if(negs(1)*isym(1,isb).eq.negs(1).and.negs(2)*isym(2,isb).   9d21s21
     $          eq.-negs(2).and.negs(3)*isym(3,isb).eq.-negs(3))then     9d21s21
            isopt(1,iopt)=isb                                           9d21s21
           end if                                                       9d21s21
          end if                                                         9d21s21
         end if                                                          9d21s21
        end do                                                          9d21s21
       end do                                                           9d21s21
       nsopt=0                                                           9d21s21
       iheadr=0                                                         1d24s23
       do iopt=1,4                                                      9d21s21
        if(idorel.lt.0.or.iopt.gt.1)then                                9d21s21
         if(iopt.eq.1.or.iopt.eq.3)then                                 9d21s21
          fname(1:3)='Re '                                              9d21s21
          irori=0                                                       9d21s21
         else                                                           9d21s21
          fname(1:3)='Im '                                              9d21s21
          irori=1                                                       9d21s21
         end if                                                         9d21s21
         if(iopt.le.2)then                                              9d21s21
          fname(4:7)='  0 '                                             9d21s21
          idelms=0                                                      9d21s21
         else                                                           9d21s21
          fname(4:7)='+-1 '                                             9d21s21
          idelms=1                                                      9d21s21
         end if                                                         9d21s21
         if(iopt.eq.1.or.iopt.eq.4)then                                 9d21s21
          fname(8:8)='+'                                                9d21s21
          ispini=0                                                      9d21s21
         else                                                           9d21s21
          fname(8:8)='-'                                                9d21s21
          ispini=1                                                      9d21s21
         end if                                                         9d21s21
         if(idorel.eq.-4)then                                           9d21s21
          fname(11:11)='G'                                              9d21s21
          if(lprint.and.iheadr.eq.0)                                    1d24s23
     $         write(6,*)('>Coulomb-Gaunt Spin-orbit; 4v ints are '),   2d8s23
     $         onoff                                                    2d8s23
          iheadr=1                                                      1d24s23
         else if(idorel.gt.0)then                                       9d21s21
          fname(11:11)='C'                                               9d21s21
          if(lprint.and.iheadr.eq.0)write(6,*)('>Coulomb Spin-orbit;'), 2d8s23
     $         (' 4v ints are '),onoff                                  2d8s23
          iheadr=1                                                      1d24s23
         else                                                           9d21s21
          fname(11:11)='B'                                              9d21s21
          if(lprint.and.iheadr.eq.0)                                    1d24s23
     $         write(6,*)('>Coulomb-Breit Spin-orbit; 4v ints are '),   2d8s23
     $         onoff                                                    2d8s23
          iheadr=1                                                      1d24s23
         end if                                                         9d21s21
         if(iopt.eq.1)then                                              9d21s21
          if(idorel.eq.-4)then                                          9d21s21
           fname(9:11)=' G '                                             9d21s21
          else                                                          9d21s21
           fname(9:11)=' B '                                            9d21s21
          end if                                                        9d21s21
         else if(iopt.eq.2)then                                         9d21s21
          fname(9:10)='Lz'                                              9d21s21
         else if(iopt.eq.3)then                                         9d21s21
          fname(9:10)='Ly'                                              9d21s21
         else                                                           9d21s21
          fname(9:10)='Lx'                                              9d21s21
         end if                                                         9d21s21
         if(idorel.lt.0)then                                            3d24s22
          ngot=4                                                        3d24s22
         else                                                           3d24s22
          ngot=3                                                        3d24s22
         end if                                                         3d24s22
         if(ngot.gt.0)then                                              9d21s21
          if(nsopt.eq.0.and.lprint)write(6,1849)                        3d2s22
          nsopt=nsopt+1                                                   9d21s21
          isopt(1,nsopt)=isopt(1,iopt)                                  9d21s21
          isopt(2,nsopt)=irori                                           9d21s21
          isopt(3,nsopt)=idelms                                          9d21s21
          isopt(4,nsopt)=ispini                                          9d21s21
          if(lprint)write(6,1848)nsopt,fname(1:3),fname(4:7),fname(8:8),3d2s22
     $         fname(9:11),isopt(1,nsopt)                               9d21s21
         end if                                                         9d21s21
        end if                                                          9d21s21
       end do                                                           9d21s21
       if(nsopt.eq.0)then                                                9d21s21
        if(lprint)                                                      3d2s22
     $      write(6,*)('there are no spin orbit operators to compute.') 3d2s22
        return                                                          9d21s21
       end if                                                           9d21s21
       call parap42(natom,ngaus,ibdat,ixmt,isym,iapair,ibstor,isstor,   2d15s22
     $      iso,idwsdeb,ascale,ipt,npt,opdata,iopdata,iosym,nop,multh,  2d15s22
     $     nbb,nbasp,nbasdws,iorb,opname,isopt,nsopt,id,nopso,          2d17s22
     $      iopso,idoubo,bc,ibc)                                        11d9s22
       call parah042c(natom,ngaus,ibdat,ih0,ih0i,iovr,ibstor,isstor,    3d19s20
     $     idwsdeb,ascale,nbasp,iorb,iapair,multh,ih0n,isopt,nsopt,srh, 3d11s22
     $      isym,bc,ibc,smsz)                                           11d16s23
       iifmx=ibcoff                                                     10d8s21
       ibcoff=iifmx+32*4                                                10d8s21
       isnd=ibcoff                                                      3d6s20
       nsnd=isnd+mynprocg                                               3d6s20
       ircv=nsnd+mynprocg                                               3d6s20
       nrcv=ircv+mynprocg                                               3d6s20
       ihrow=nrcv+mynprocg                                              3d6s20
       ibcoff=ihrow+mynprocg*64                                         9d24s21
       call enough('prop. 16',bc,ibc)
       call second(timecf0)                                             10d8s21
       call paraeri42cf(natom,ngaus,ibdat,ibstor,isstor,idwsdeb,         3d5s20
     $       ascale,nbasp,iorb,iapair,idoubo,irefo,idorel,ibc(isnd),    9d22s21
     $      ibc(nsnd),ibc(ircv),ibc(nrcv),ibc(ihrow),multh,i4or,ionexr, 10d5s21
     $      jmatsr,kmatsr,kmatsrb,i3xr,                                 11d15s21
     $      shift,ih0n,isopt,nsopt,scals,srh,isym,ibc(iifmx),ntype,nh0, 10d20s21
     $      i4o,bc,ibc)                                                 11d9s22
       call second(timecf1)                                             10d8s21
       telap=timecf1-timecf0
       if(lprint)write(6,*)('time for paraeri42cf'),telap
       call relgammas(nec,mdon,mdoo,irw0,irw1,irw2,spinx,bc,ibc)        11d9s22
       npack4=iwavedat(6,1)                                             5d14s21
       if(ipack1(2).eq.2)then                                           5d14s21
        iuniq=ibcoff                                                     5d27s21
        ibcoff=iuniq+nfname                                             5d27s21
        call enough('prop. 17',bc,ibc)
        nuniq=0                                                         5d27s21
        is=1                                                            5d27s21
 1001   continue                                                        5d27s21
         npack4=iwavedat(6,is)                                          5d27s21
         llam=ipack1(4)                                                 5d27s21
         if(ipack1(3).eq.iabs(llam))then                                5d27s21
          ibc(iuniq+nuniq)=is                                           5d27s21
          nuniq=nuniq+1                                                 5d27s21
          if(ipack1(3).ne.0)is=is+1                                     5d27s21
         end if                                                         5d27s21
         is=is+1                                                        5d27s21
        if(is.le.nfname)go to 1001                                      5d27s21
        do ibra=0,nuniq-1                                               5d27s21
         ibraf=ibc(iuniq+ibra)                                          5d27s21
         do iket=0,ibra                                                 5d27s21
          iketf=ibc(iuniq+iket)                                         5d27s21
          call solin(iwavedat(1,ibraf),iwavedat(1,iketf),nspc,nec,multh,5d27s21
     $         irefo,ih0a,i4o,irel,ism,norb,mdon,nvirt,maxbx,           10d22s21
     $         maxbxd,srh,sr2,nsymb,irw0,irw1,irw2,ih0n,nh0,isopt,nsopt,10d13s21
     $         i4or,ionexr,jmatsr,kmatsr,kmatsrb,i3xr,ibc(iifmx),ntype, 11d15s21
     $         npadddi,nbasp,nbaspc,natom,ngaus,ibdat,iapair,ibstor,    12d20s20
     $         isstor,isym,ascale,idorel,iorb,lprint,bc,ibc,shift,n4vso)2d8s23
          if(ntmso.ne.0)then                                            8d16s22
           call gtmomso(iwavedat(1,ibraf),idum,iwavedat(1,iketf),        3d8s22
     $         idum,idum,idum,nspc,nec,multh,irefo,                     3d8s22
     $     ixmt,nh0,nopso,iopso,irel,ism,norb,mdon,dum,                 12d1s22
     $         nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,irw2,         2d17s22
     $         isopt,nsopt,                                             2d17s22
     $         ibc(iifmx),ntype,npadddi,nbasp,nbaspc,natom,ngaus,ibdat, 12d20s20
     $         iapair,ibstor,isstor,isym,ascale,idorel,iorb,opname,0,0, 3d1s22
     $        idum,idum,dum,1,idum,idum,1,dum,lprint,l2sub,bc,ibc,n4vso)2d8s23
          end if                                                        8d16s22
         end do                                                         5d27s21
        end do                                                          5d27s21
       else if(ipack1(2).eq.6)then                                      3d17s22
        if(lprint)write(6,*)('spherical system ')                       3d2s22
        icharo=ichar('o')                                               5d19s21
        ndime=0                                                         6d4s21
        ndimo=0                                                         6d4s21
        nml0e=0                                                         5d19s21
        nml0o=0                                                         5d19s21
        do iwf=1,nfname                                                  5d18s21
         npack4=iwavedat(6,iwf)                                          5d14s21
         if(ipack1(4).eq.0)then                                         5d24s21
          if(iwavedat(15,iwf).eq.icharo)then                              5d19s21
           nml0o=nml0o+1                                                5d19s21
          else                                                          5d19s21
           nml0e=nml0e+1                                                5d19s21
          end if                                                        5d19s21
          ll=ipack1(3)*2                                                 5d18s21
          isml=iwavedat(1,iwf)-1                                        5d18s21
          jjl=iabs(ll-isml)                                               5d18s21
          jjh=ll+isml                                                     5d18s21
          if(iwavedat(15,iwf).eq.icharo)then                            5d19s21
           if(nml0o.eq.1)then                                               5d18s21
            jjmino=jjl                                                     5d18s21
            jjmaxo=jjh                                                     5d18s21
           else                                                           5d18s21
            jjmino=min(jjmino,jjl)                                          5d18s21
            jjmaxo=max(jjmaxo,jjh)                                          5d18s21
           end if                                                         5d18s21
           ndimo=ndimo+((jjh+2-jjl)/2)*iwavedat(3,iwf)                  6d4s21
          else                                                          5d19s21
           if(nml0e.eq.1)then                                               5d18s21
            jjmine=jjl                                                     5d18s21
            jjmaxe=jjh                                                     5d18s21
           else                                                           5d18s21
            jjmine=min(jjmine,jjl)                                          5d18s21
            jjmaxe=max(jjmaxe,jjh)                                          5d18s21
           end if                                                         5d18s21
           ndime=ndime+((jjh+2-jjl)/2)*iwavedat(3,iwf)                  6d4s21
          end if                                                        5d19s21
         end if
        end do                                                           5d18s21
        if(ndime.gt.0)then                                              1d6s23
         numje=((jjmaxe-jjmine)/2)+1                                        5d18s21
        else                                                            1d6s23
         numje=0                                                        1d6s23
         jjmine=0                                                       1d6s23
         jjmaxe=0                                                       1d6s23
        end if                                                          1d6s23
        if(ndimo.gt.0)then                                              1d6s23
         numjo=((jjmaxo-jjmino)/2)+1                                        5d18s21
        else                                                            1d6s23
         numjo=0                                                        1d6s23
         jjmino=0                                                       1d6s23
         jjmaxo=0                                                       1d6s23
        end if                                                          1d6s23
        iqne=ibcoff                                                     6d4s21
        iqno=iqne+ndime                                                 6d4s21
        ieige=iqno+ndimo                                                6d4s21
        ieigo=ieige+ndime                                               6d4s21
        itmee=ieigo+ndimo                                               6d4s21
        itmoo=itmee+ndime*ndime*2                                       6d4s21
        itmoe=itmoo+ndimo*ndimo*2                                       6d4s21
        itmeo=itmoe+ndime*ndimo                                         6d4s21
        ibcoff=itmeo+ndime*ndimo                                        6d4s21
        iptfe=ibcoff                                                     5d18s21
        ijjse=iptfe+nml0e                                                  5d18s21
        iptfo=ijjse+numje                                                5d18s21
        ijjso=iptfo+nml0o                                                  5d18s21
        ibcoff=ijjso+numjo                                                5d18s21
        call enough('prop. 18',bc,ibc)
        do i=itmee,iptfe-1                                              6d4s21
         bc(i)=0d0                                                      6d4s21
        end do                                                          6d4s21
        do i=0,numje-1                                                   5d18s21
         ibc(ijjse+i)=0                                                  5d18s21
        end do                                                          5d18s21
        do i=0,numjo-1                                                   5d18s21
         ibc(ijjso+i)=0                                                  5d18s21
        end do                                                          5d18s21
        nml0e=0                                                          5d18s21
        nml0o=0                                                          5d18s21
        do iwf=1,nfname                                                  5d18s21
         npack4=iwavedat(6,iwf)                                          5d14s21
         if(ipack1(4).eq.0)then                                         5d24s21
          ll=ipack1(3)*2                                                 5d18s21
          isml=iwavedat(1,iwf)-1                                        5d18s21
          jjl=iabs(ll-isml)                                               5d18s21
          jjh=ll+isml                                                     5d18s21
          nroot=iwavedat(3,iwf)                                          5d18s21
          if(iwavedat(15,iwf).eq.icharo)then                              5d19s21
           ibc(iptfo+nml0o)=iwf                                            5d18s21
           nml0o=nml0o+1                                                   5d18s21
           do jj=jjl,jjh,2                                                5d18s21
            kk=ijjso+((jj-jjmino)/2)                                    5d19s21
            ibc(kk)=ibc(kk)+nroot                                         5d18s21
           end do                                                         5d18s21
          else                                                          5d19s21
           ibc(iptfe+nml0e)=iwf                                            5d18s21
           nml0e=nml0e+1                                                   5d18s21
           do jj=jjl,jjh,2                                                5d18s21
            kk=ijjse+((jj-jjmine)/2)                                    5d19s21
            ibc(kk)=ibc(kk)+nroot                                         5d18s21
           end do                                                         5d18s21
          end if                                                        5d19s21
         end if                                                         5d18s21
        end do                                                           5d18s21
        ijjte=ibcoff                                                    6d4s21
        ijjme=ijjte+numje                                               6d4s21
        ibcoff=ijjme+numje                                                5d18s21
        call enough('prop. 19',bc,ibc)
        njjte=0                                                         6d4s21
        do jj=jjmine,jjmaxe,2                                             5d18s21
         kk=ijjse+((jj-jjmine)/2)                                         5d18s21
         nn=ibc(kk)                                                     5d18s21
         ll=ijjte+((jj-jjmine)/2)                                       6d4s21
         ibc(ll)=njjte                                                  6d4s21
         njjte=njjte+nn                                                 6d4s21
         ll=ijjme+((jj-jjmine)/2)                                         5d18s21
         ibc(ll)=ibcoff                                                 5d18s21
         ibcoff=ibcoff+nn*nn+nn*2                                       5d19s21
        end do                                                          5d18s21
        ijjto=ibcoff                                                    6d4s21
        ijjmo=ijjto+numjo                                               6d4s21
        ibcoff=ijjmo+numjo                                                5d18s21
        njjto=0                                                         6d4s21
        call enough('prop. 20',bc,ibc)
        do jj=jjmino,jjmaxo,2                                             5d18s21
         kk=ijjso+((jj-jjmino)/2)                                         5d18s21
         nn=ibc(kk)                                                     5d18s21
         ll=ijjto+((jj-jjmino)/2)                                       6d4s21
         ibc(ll)=njjto                                                  6d4s21
         njjto=njjto+nn                                                 6d4s21
         ll=ijjmo+((jj-jjmino)/2)                                         5d18s21
         ibc(ll)=ibcoff                                                 5d18s21
         ibcoff=ibcoff+nn*nn+nn*2                                       5d19s21
        end do                                                          5d18s21
        call enough('prop. 21',bc,ibc)
        if(lprint)write(6,*)('for even parity ')                        3d2s22
        ijjobe=ibcoff                                                     5d18s21
        ijjoke=ijjobe+numje                                                5d18s21
        ijjokue=ijjoke+numje                                               5d19s21
        ibcoff=ijjokue+numje                                              5d19s21
        call enough('prop. 22',bc,ibc)
        do i=0,numje-1                                                   5d18s21
         ibc(ijjoke+i)=0                                                 5d18s21
        end do                                                          5d18s21
        do iket=1,nml0e                                                  5d18s21
         iwfk=ibc(iptfe+iket-1)                                          5d18s21
         ieig=iwavedat(4,iwfk)+iwavedat(13,iwfk)                        8d27s21
         ipack8=ibc(ieig)                                               8d27s21
         ncsft=ipack4(1)                                                8d27s21
         ieig=ieig+ncsft*iwavedat(3,iwfk)+1                             8d27s21
         do i=0,numje-1                                                  5d18s21
          ibc(ijjobe+i)=0                                                5d18s21
         end do                                                         5d18s21
         do ibra=1,iket                                                 5d18s21
          do i=0,numje-1                                                 5d19s21
           ibc(ijjokue+i)=ibc(ijjoke+i)                                   5d19s21
          end do                                                        5d19s21
          iwfb=ibc(iptfe+ibra-1)                                         5d18s21
          jtmombk=itmom+4*(iwfb-1+nfname*(iwfk-1))                         5d25s21
          jtmomkb=itmom+4*(iwfk-1+nfname*(iwfb-1))                         5d25s21
          na=0                                                          6d4s21
          nb=0                                                          6d4s21
          do i=0,3                                                      6d4s21
           if(ibc(jtmombk+i).ne.0)na=1                                  6d4s21
           if(ibc(jtmomkb+i).ne.0)nb=1                                  6d4s21
          end do                                                        6d4s21
          if(max(na,nb).gt.0)then                                       6d4s21
           if(na.ne.0)then                                              6d4s21
            call fillrup(iwavedat(1,iwfb),iwavedat(1,iwfk),ibc(jtmombk),6d4s21
     $          jjmine,jjmine,bc(itmee),ndime,ibc(ijjte),ibc(ijjte),    6d4s21
     $          ibc(ijjobe),ibc(ijjokue),ndime,dum,bc,ibc)              1d13s23
           else                                                         6d4s21
            call fillrup(iwavedat(1,iwfk),iwavedat(1,iwfb),ibc(jtmomkb),6d4s21
     $          jjmine,jjmine,bc(itmee),ndime,ibc(ijjte),ibc(ijjte),    6d4s21
     $           ibc(ijjokue),ibc(ijjobe),ndime,dum,bc,ibc)             1d13s23
           end if                                                       6d4s21
          end if                                                        6d4s21
          if(ntmso.ne.0)then                                            8d16s22
c     even parity, tmpb is dum, tm is itmee
           call gtmomso(iwavedat(1,iwfb),ibc(ijjobe),iwavedat(1,iwfk),   5d18s21
     $         ibc(ijjokue),ibc(ijjse),ibc(ijjme),nspc,nec,multh,irefo, 5d19s21
     $     ixmt,nh0,nopso,iopso,irel,ism,norb,mdon,bc(ieig),            12d1s22
     $         nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,irw2,         2d17s22
     $         isopt,nsopt,                                             2d17s22
     $         ibc(iifmx),ntype,npadddi,nbasp,nbaspc,natom,ngaus,ibdat, 12d20s20
     $         iapair,ibstor,isstor,isym,ascale,idorel,iorb,opname,0,0, 3d1s22
     $         jjmine,jjmine,bc(itmee),ndime,ibc(ijjte),ibc(ijjte),     3d1s22
     $         ndime,dum,lprint,l2sub,bc,ibc,n4vso)                     2d8s23
          end if                                                        8d16s22
          jmasspa=imasspa+iwfb-1                                        8d16s22
          call sospher2(iwavedat(1,iwfb),ibc(ijjobe),iwavedat(1,iwfk),   5d18s21
     $         ibc(ijjokue),ibc(ijjse),ibc(ijjme),nspc,nec,multh,irefo, 5d19s21
     $         ih0a,i4o,irel,ism,norb,mdon,bc(ieig),jjmine,iwfk,        10d22s21
     $         nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,irw2,ih0n,nh0,10d13s21
     $         isopt,nsopt,i4or,ionexr,jmatsr,kmatsr,kmatsrb,i3xr,      11d15s21
     $         ibc(iifmx),ntype,npadddi,nbasp,nbaspc,natom,ngaus,ibdat, 12d20s20
     $         iapair,ibstor,isstor,isym,ascale,idorel,iorb,lprint,     8d16s22
     $         bc(ibc(jmasspa)),bc,ibc,shift,n4vso)                     2d8s23
         end do                                                         5d18s21
         do i=0,numje-1                                                  5d19s21
          ibc(ijjoke+i)=ibc(ijjokue+i)                                    5d19s21
         end do                                                         5d19s21
        end do                                                          5d18s21
        if(lprint)write(6,*)('for odd parity ')                         3d2s22
        ijjobo=ibcoff                                                     5d18s21
        ijjoko=ijjobo+numjo                                                5d18s21
        ijjokuo=ijjoko+numjo                                               5d19s21
        ibcoff=ijjokuo+numjo                                              5d19s21
        call enough('prop. 23',bc,ibc)
        do i=0,numjo-1                                                   5d18s21
         ibc(ijjoko+i)=0                                                 5d18s21
        end do                                                          5d18s21
        do iket=1,nml0o                                                  5d18s21
         iwfk=ibc(iptfo+iket-1)                                          5d18s21
         ieig=iwavedat(4,iwfk)+iwavedat(13,iwfk)                        8d27s21
         ipack8=ibc(ieig)                                               8d27s21
         ncsft=ipack4(1)                                                8d27s21
         ieig=ieig+ncsft*iwavedat(3,iwfk)+1                             8d27s21
         do i=0,numje-1                                                  5d18s21
          ibc(ijjobe+i)=0                                                5d18s21
         end do                                                         5d18s21
         do ibra=1,nml0e                                                6d4s21
          iwfb=ibc(iptfe+ibra-1)                                         5d18s21
          jtmombk=itmom+4*(iwfb-1+nfname*(iwfk-1))                         5d25s21
          jtmomkb=itmom+4*(iwfk-1+nfname*(iwfb-1))                         5d25s21
          na=0                                                          6d4s21
          nb=0                                                          6d4s21
          do i=0,3                                                      6d4s21
           if(ibc(jtmombk+i).ne.0)na=1                                  6d4s21
           if(ibc(jtmomkb+i).ne.0)nb=1                                  6d4s21
          end do                                                        6d4s21
          if(max(na,nb).gt.0)then                                       6d4s21
           if(na.ne.0)then                                              6d4s21
            call fillrup(iwavedat(1,iwfb),iwavedat(1,iwfk),ibc(jtmombk),6d4s21
     $          jjmine,jjmino,bc(itmeo),ndime,ibc(ijjte),ibc(ijjto),    6d4s21
     $          ibc(ijjobe),ibc(ijjoko),ndimo,bc(itmoe),bc,ibc)         11d14s22
            if(ntmso.ne.0)then                                          8d16s22
c     bra is even parity, ket is odd parity
c     tm is itmeo, tmb is itmoe
             call gtmomso(iwavedat(1,iwfb),ibc(ijjobe),iwavedat(1,iwfk),12d1s22
     $         ibc(ijjoko),idum,idum,nspc,nec,multh,irefo,              12d1s22
     $     ixmt,nh0,nopso,iopso,irel,ism,norb,mdon,bc(ieig),            12d1s22
     $         nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,irw2,         2d17s22
     $         isopt,nsopt,                                             2d17s22
     $         ibc(iifmx),ntype,npadddi,nbasp,nbaspc,natom,ngaus,ibdat, 12d20s20
     $         iapair,ibstor,isstor,isym,ascale,idorel,iorb,opname,0,1, 12d1s22
     $         jjmine,jjmino,bc(itmeo),ndime,ibc(ijjte),ibc(ijjto),     12d1s22
     $           ndimo,bc(itmoe),lprint,l2sub,bc,ibc,n4vso)             2d8s23
            end if                                                        8d16s22
           else                                                         6d4s21
            call fillrup(iwavedat(1,iwfk),iwavedat(1,iwfb),ibc(jtmomkb),6d4s21
     $          jjmino,jjmine,bc(itmoe),ndimo,ibc(ijjto),ibc(ijjte),    6d4s21
     $           ibc(ijjoko),ibc(ijjobe),ndime,bc(itmeo),bc,ibc)        11d14s22
            if(ntmso.ne.0)then                                          8d16s22
c     bra is odd parity, ket is evend parity
c     tm is itmeo, tmb is itmoe
             call gtmomso(iwavedat(1,iwfk),ibc(ijjoko),                 12d1s22
     $           iwavedat(1,iwfb),ibc(ijjobe),idum,idum,nspc,nec,multh, 12d1s22
     $           irefo,ixmt,nh0,nopso,iopso,irel,ism,norb,mdon,bc(ieig),12d1s22
     $           nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,            12d1s22
     $           irw2,isopt,nsopt,ibc(iifmx),ntype,npadddi,nbasp,nbaspc,12d1s22
     $           natom,ngaus,ibdat,iapair,ibstor,isstor,isym,ascale,    12d1s22
     $           idorel,iorb,opname,1,0,jjmino,jjmine,bc(itmoe),ndimo,  12d1s22
     $           ibc(ijjto),ibc(ijjte),ndime,bc(itmeo),lprint,l2sub,bc, 12d1s22
     $           ibc,n4vso)                                             2d8s23
            end if                                                        8d16s22
           end if                                                       6d4s21
          end if                                                        6d4s21
          ise=iwavedat(1,iwfb)-1                                        6d4s21
          nroote=iwavedat(3,iwfb)                                       6d4s21
          npack4=iwavedat(6,iwfb)                                       6d4s21
          lle=ipack1(3)*2                                               6d4s21
          jjl=iabs(ise-lle)                                             6d4s21
          jjh=ise+lle                                                   6d4s21
          do i=jjl,jjh,2                                                6d4s21
           ii=(i-jjmine)/2                                              6d4s21
           ibc(ijjobe+ii)=ibc(ijjobe+ii)+nroote                         6d4s21
          end do                                                        6d4s21
         end do
         do i=0,numjo-1                                                  5d18s21
          ibc(ijjobo+i)=0                                                5d18s21
         end do                                                         5d18s21
         do ibra=1,iket                                                 5d18s21
          do i=0,numjo-1                                                 5d19s21
           ibc(ijjokuo+i)=ibc(ijjoko+i)                                   5d19s21
          end do                                                        5d19s21
          iwfb=ibc(iptfo+ibra-1)                                         5d18s21
          jtmombk=itmom+4*(iwfb-1+nfname*(iwfk-1))                         5d25s21
          jtmomkb=itmom+4*(iwfk-1+nfname*(iwfb-1))                         5d25s21
          na=0                                                          6d4s21
          nb=0                                                          6d4s21
          do i=0,3                                                      6d4s21
           if(ibc(jtmombk+i).ne.0)na=1                                  6d4s21
           if(ibc(jtmomkb+i).ne.0)nb=1                                  6d4s21
          end do                                                        6d4s21
          if(max(na,nb).gt.0)then                                       6d4s21
           if(na.ne.0)then                                              6d4s21
            call fillrup(iwavedat(1,iwfb),iwavedat(1,iwfk),ibc(jtmombk),6d4s21
     $          jjmino,jjmino,bc(itmoo),ndimo,ibc(ijjto),ibc(ijjto),    6d4s21
     $          ibc(ijjobo),ibc(ijjokuo),ndimo,dum,bc,ibc)              1d13s23
           else                                                         6d4s21
            call fillrup(iwavedat(1,iwfk),iwavedat(1,iwfb),ibc(jtmomkb),6d4s21
     $          jjmino,jjmino,bc(itmoo),ndimo,ibc(ijjto),ibc(ijjto),    6d4s21
     $           ibc(ijjobo),ibc(ijjokuo),ndimo,dum,bc,ibc)             1d13s23
           end if                                                       6d4s21
          end if                                                        6d4s21
          if(ntmso.ne.0)then                                            8d16s22
c     bra is odd parity, ket is odd parity
c     tm is itmoo, tmb is dum
           call gtmomso(iwavedat(1,iwfb),ibc(ijjobo),iwavedat(1,iwfk),   5d18s21
     $         ibc(ijjokuo),idum,idum,nspc,nec,multh,irefo,             2d25s22
     $     ixmt,nh0,nopso,iopso,irel,ism,norb,mdon,bc(ieig),            12d1s22
     $         nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,irw2,         2d17s22
     $         isopt,nsopt,                                             2d17s22
     $         ibc(iifmx),ntype,npadddi,nbasp,nbaspc,natom,ngaus,ibdat, 12d20s20
     $         iapair,ibstor,isstor,isym,ascale,idorel,iorb,opname,1,1, 3d1s22
     $         jjmino,jjmino,bc(itmoo),ndimo,ibc(ijjto),ibc(ijjto),     1d13s23
     $         ndimo,dum,lprint,l2sub,bc,ibc,n4vso)                     2d8s23
          end if                                                        8d16s22
          jmasspa=imasspa+iwfb-1                                        8d16s22
          call sospher2(iwavedat(1,iwfb),ibc(ijjobo),iwavedat(1,iwfk),   5d18s21
     $         ibc(ijjokuo),ibc(ijjso),ibc(ijjmo),nspc,nec,multh,irefo, 5d19s21
     $         ih0a,i4o,irel,ism,norb,mdon,bc(ieig),jjmino,iwfk,        10d22s21
     $         nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,irw2,ih0n,nh0,10d13s21
     $         isopt,nsopt,i4or,ionexr,jmatsr,kmatsr,kmatsrb,i3xr,      11d15s21
     $         ibc(iifmx),ntype,npadddi,nbasp,nbaspc,natom,ngaus,ibdat, 12d20s20
     $         iapair,ibstor,isstor,isym,ascale,idorel,iorb,lprint,     8d16s22
     $         bc(ibc(jmasspa)),bc,ibc,shift,n4vso)                     2d8s23
         end do                                                         5d18s21
         do i=0,numjo-1                                                  5d19s21
          ibc(ijjoko+i)=ibc(ijjokuo+i)                                    5d19s21
         end do                                                         5d19s21
        end do                                                          5d18s21
        if(mynowprog.eq.0)then                                          1d24s23
         call sospher3(numje,ibc(ijjse),ibc(ijjme),jjmine,iwavedat,nspc, 6d1s21
     $        elowest,bc(ieige),ibc(iqne),bc(itmee),ndime,bc(itmeo),    6d4s21
     $       bc(itmoe),ndimo,lprint,bc,ibc)                             11d9s22
         call sospher3(numjo,ibc(ijjso),ibc(ijjmo),jjmino,iwavedat,nspc, 6d1s21
     $       elowest,bc(ieigo),ibc(iqno),bc(itmoo),ndimo,bc(itmoe),     6d4s21
     $       bc(itmeo),ndime,lprint,bc,ibc)                             11d9s22
         call atomll(iwavedat,nspc,bc(ieige),ibc(iqne),bc(itmee),ndime,  6d4s21
     $       bc(ieigo),ibc(iqno),bc(itmoo),ndimo,bc(itmeo),bc(itmoe),   3d2s22
     $       lprint,bc,ibc)                                             11d14s22
        end if                                                          1d24s23
       else                                                             3d17s22
        if(lprint)write(6,*)('general case ')                           3d29s22
        do ibra=1,nfname                                                3d17s22
         do iket=1,ibra                                                 3d17s22
          call sogen(iwavedat(1,ibra),iwavedat(1,iket),nspc,nec,multh,  3d17s22
     $         irefo,ih0a,i4o,irel,ism,norb,mdon,nvirt,maxbx,           10d22s21
     $         maxbxd,srh,sr2,nsymb,irw0,irw1,irw2,ih0n,nh0,isopt,nsopt,10d13s21
     $         i4or,ionexr,jmatsr,kmatsr,kmatsrb,i3xr,ibc(iifmx),ntype, 11d15s21
     $         npadddi,nbasp,nbaspc,natom,ngaus,ibdat,iapair,ibstor,    12d20s20
     $         isstor,isym,ascale,idorel,iorb,lprint,bc,ibc,shift,n4vso)2d8s23
          if(ntmso.ne.0)then                                            8d16s22
           call gtmomso(iwavedat(1,ibra),idum,iwavedat(1,iket),          3d17s22
     $         idum,idum,idum,nspc,nec,multh,irefo,                     3d8s22
     $     ixmt,nh0,nopso,iopso,irel,ism,norb,mdon,dum,                 12d1s22
     $         nvirt,maxbx,maxbxd,srh,sr2,nsymb,irw0,irw1,irw2,         2d17s22
     $         isopt,nsopt,                                             2d17s22
     $         ibc(iifmx),ntype,npadddi,nbasp,nbaspc,natom,ngaus,ibdat, 12d20s20
     $         iapair,ibstor,isstor,isym,ascale,idorel,iorb,opname,0,0, 3d1s22
     $        idum,idum,dum,1,idum,idum,1,dum,lprint,l2sub,bc,ibc,n4vso)2d8s23
          end if                                                        8d16s22
         end do                                                         3d17s22
        end do                                                          3d17s22
       end if                                                           5d14s21
       end if                                                           10d27s20
      itmp=ibcoff                                                       5d4s22
      ibcoff=itmp+22                                                    5d4s22
      call enough('prop. 24',bc,ibc)
      jtmp=itmp-1                                                       5d4s22
      ktmp=jtmp+11                                                      5d4s2
      do i=1,11                                                         5d4s22
       bc(jtmp+i)=tso(i)                                                5d4s22
       bc(ktmp+i)=tso(i)**2                                             5d4s22
      end do                                                            5d4s22
      call dws_gsumf(bc(itmp),22)                                       5d4s22
      if(mynowprog.eq.0)then                                            5d4s22
       xxx=1d0/dfloat(mynprocg)                                          5d4s22
       do i=1,11                                                         5d4s22
        bc(jtmp+i)=bc(jtmp+i)*xxx                                        5d4s22
        bc(ktmp+i)=bc(ktmp+i)*xxx                                        5d4s22
        bc(ktmp+i)=sqrt(abs(bc(ktmp+i)-bc(jtmp+i)**2))                   5d4s22
         if(i.eq.1)then                                                 11d17s22
          code='ii'                                                     11d17s22
         else if(i.eq.2)then                                            11d17s22
          code='si'                                                     11d17s22
         else if(i.eq.3)then                                            11d17s22
          code='is'                                                     11d17s22
         else if(i.eq.4)then                                            11d17s22
          code='ss'                                                     11d17s22
         else if(i.eq.5)then                                            11d17s22
          code='di'                                                     11d17s22
         else if(i.eq.6)then                                            11d17s22
          code='id'                                                     11d17s22
         else if(i.eq.7)then                                            11d17s22
          code='dd'                                                     11d17s22
         else if(i.eq.8)then                                            11d17s22
          code='4v'
         else if(i.eq.9)then                                            11d17s22
          code='ds'                                                     11d17s22
         else if(i.eq.10)then                                            11d17s22
          code='sd'                                                     11d17s22
         else if(i.eq.11)then                                            11d17s22
          code='dv'
         end if                                                         11d17s22
        write(6,*)('psisopsi timer '),code,bc(jtmp+i),bc(ktmp+i)           5d4s22
       end do                                                            5d4s22
      end if                                                            5d4s22
      ibcoff=itmp                                                       5d4s22
      ibcoff=ibcoffo                                                    11d24s19
      return
      end
