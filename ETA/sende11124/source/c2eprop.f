c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine c2eprop(iwavedat,ixw1,ixw2,ncsf,nec,nspc,ixmt,opdata,  5d20s21
     $     iopdata,iosym,nop,opname,nbasdws,nsymb,mdon,mdoop,           5d12s21
     $     ism,irel,irefo,norb,multh,idoubo,i2eop,imassp,natom,lfn,     7d9s21
     $     nvirt,maxbx,maxbxd,srh,sr2,npadddi,xmasstot,bc,ibc)          11d10s22
      implicit real*8 (a-h,o-z)                                         5d12s21
      external second                                                   10d14s22
      integer*1 ipack1(4)                                               5d12s21
      integer*4 ipack4(2)                                               5d12s21
      integer*8 ipack8,imassp                                           5d25s21
      character*(*) lfn                                                 5d26s21
      character*10 opname(*)                                            5d12s21
      character*6 wname                                                 5d13s21
      logical lri,lprint                                                3d2s22
      equivalence (ipack1,npack4)                                       5d12s21
      equivalence (ipack4,ipack8)                                       5d12s21
      dimension iwavedat(nspc),ncsf(*),ixmt(8,*),opdata(*),             5d20s21
     $     iopdata(7,*),iosym(*),nbasdws(*),ism(*),irel(*),irefo(*),    5d27s21
     $     multh(8,8),idoubo(*),i2eop(6,*),nvirt(*),ixmtf(8),ixsum(8,4),2d1s23
     $     i2eops(2,3),iosyms(4)                                        2d1s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      include "common.store"                                            5d12s21
      include "common.basis"                                             5d21s21
      lprint=mynowprog.eq.0                                             3d2s22
      npack4=iwavedat(6)                                                5d20s21
      do i=1,6                                                          5d20s21
       if(iwavedat(13+i).ne.0)then                                      5d20s21
        wname(i:i)=char(iwavedat(13+i))                                 5d20s21
        nnn=i                                                           5d13s21
       end if                                                           5d13s21
      end do                                                            5d13s21
      ll=ipack1(3)                                                      5d21s21
      if(lprint)                                                        3d2s22
     $     write(6,*)('2eprop for '),iwavedat(1),wname(1:nnn)           3d31s22
      if(xmasstot.eq.0d0)then                                           3d31s22
       if(lprint)then                                                    3d2s22
        do ia=1,natom                                                     5d26s21
         write(6,*)('for atom '),ia,atnum(1,ia),atom(ia)                  5d26s21
         if(mynowprog.eq.0)then                                           5d26s21
          call getm(atom(ia),1,xmass,1,lfn//'/atwgts')                   3d29s22
          write(6,*)('xmass = '),xmass                                    5d26s21
          xmasstot=xmasstot+xmass                                         5d26s21
         end if                                                           5d26s21
        end do                                                            5d26s21
       end if
       call dws_bcast(xmasstot,1)                                        5d26s21
      end if                                                            3d31s22
      if(lprint)write(6,*)('xmasstot = '),xmasstot                                5d26s21
      xmassi=219474.635d0/xmasstot                                      5d26s21
      xmassi2=1d0/xmasstot                                              8d16s22
      nlzz=ipack1(2)                                                    5d25s21
      if(nlzz.eq.2)then                                                 5d25s21
       llzz=ll*ll                                                       5d25s21
      else if(nlzz.eq.6)then                                            5d25s21
       llzz=ipack1(4)**2                                                5d25s21
      end if                                                            5d25s21
      nroot=iwavedat(3)                                                 5d20s21
      imassp=ibcoff                                                     5d24s21
      ntri=(nroot*(nroot+1))/2                                          5d24s21
      ibcoff=imassp+ntri                                                5d24s21
      ibcoffo=ibcoff                                                    5d26s22
      isump=ibcoff                                                      5d21s21
      ibcoff=isump+nroot*nroot                                          5d21s21
      call enough('c2eprop.  1',bc,ibc)
      do i=isump,ibcoff-1                                               5d21s21
       bc(i)=0d0                                                        5d21s21
      end do                                                            5d21s21
c
c     Lz^2
c
      if(nlzz.ne.6)then                                                 2d1s23
       n2e=1
       fact2e=-2d0                                                       5d21s21
       do i=1,nop                                                       2d1s23
        if(opname(i)(1:4).eq.'Lz^2')then                                2d1s23
         ikin=i                                                         2d1s23
         go to 112                                                      2d1s23
        end if                                                          2d1s23
       end do                                                           2d1s23
  112  continue                                                         2d1s23
       call psioppsi(iwavedat,nroot,ikin,1d0,ixmt,iosym,n2e,            7d27s21
     $      i2eop(1,ikin),fact2e,iwavedat,nroot,bc(isump),mdon,mdoop,    7d27s21
     $      ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,nvirt,maxbx,   7d27s21
     $      maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,npadddi,bc,ibc)      11d10s22
       if(lprint)then                                                   2d1s23
        call pdump(bc(isump),nroot,iwavedat(1),wname(1:nnn),'Lz^2')     2d1s23
       end if                                                           2d1s23
      end if                                                            2d1s23
c
c     Lx^2+Ly^2(+Lz^2 if atom)
c
      fact2e=-2d0                                                       5d21s21
      n2e=0                                                             2d1s23
      do isb=1,nsymb
       nn=nbasdws(isb)*nbasdws(isb)                                     2d1s23
       ixsum(isb,1)=ibcoff                                              2d1s23
       ibcoff=ixsum(isb,1)+nn                                           2d1s23
       call enough('c2eprop.xsum',bc,ibc)                               2d1s23
       do iz=ixsum(isb,1),ibcoff-1                                      2d1s23
        bc(iz)=0d0                                                      2d1s23
       end do                                                           2d1s23
       iosyms(1)=1                                                      2d1s23
      end do                                                            2d1s23
      do i=0,nroot*nroot-1
       bc(isump+i)=0d0
      end do
      do i=1,nop                                                        2d1s23
       if(opname(i)(1:4).eq.'Lx^2'.or.opname(i)(1:4).eq.'Ly^2'.or.      2d1s23
     $       (nlzz.eq.6.and.opname(i)(1:4).eq.'Lz^2'))then              2d1s23
        n2e=n2e+1                                                       2d1s23
        i2eops(1,n2e)=n2e+1
        i2eops(2,n2e)=n2e+1                                             2d1s23
        iosyms(n2e+1)=iosym(i2eop(1,i))                                 2d1s23
        do isb=1,nsymb                                                  2d1s23
         nn=nbasdws(isb)*nbasdws(isb)                                     2d1s23
         do j=0,nn-1                                                    2d1s23
          bc(ixsum(isb,1)+j)=bc(ixsum(isb,1)+j)+bc(ixmt(isb,i)+j)       2d1s23
         end do                                                         2d1s23
         jsb=multh(isb,iosyms(n2e+1))                                   2d1s23
         nm=nbasdws(isb)*nbasdws(jsb)                                   2d1s23
         ixsum(isb,n2e+1)=ibcoff                                        2d1s23
         ibcoff=ixsum(isb,n2e+1)+nm                                     2d1s23
         call enough('c2eprop.xsnm',bc,ibc)                             2d1s23
         do j=0,nm-1                                                    2d1s23
          bc(ixsum(isb,n2e+1)+j)=bc(ixmt(isb,i2eop(1,i))+j)             2d1s23
         end do                                                         2d1s23
        end do                                                          2d1s23
       end if                                                           2d1s23
      end do
      call psioppsi(iwavedat,nroot,1,1d0,ixsum,iosyms,n2e,              2d1s23
     $      i2eops,fact2e,iwavedat,nroot,bc(isump),mdon,mdoop,          2d1s23
     $      ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,nvirt,maxbx,   7d27s21
     $      maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,npadddi,bc,ibc)      11d10s22
      if(lprint)then                                                    2d1s23
       if(nlzz.eq.6)then                                                 2d1s23
        call pdump(bc(isump),nroot,iwavedat(1),wname(1:nnn),             2d1s23
     $        'L^2')                                                    2d1s23
       else                                                              2d1s23
        call pdump(bc(isump),nroot,iwavedat(1),wname(1:nnn),             2d1s23
     $        'Lx^2+Ly^2')                                              3d30s22
       end if                                                           2d1s23
      end if                                                            2d1s23
c
c     mass polarization
c
      do i=0,nroot*nroot-1
       bc(isump+i)=0d0
      end do
      n2e=3
      fact2e=-1d0                                                       5d21s21
      do i=1,nop                                                        2d1s23
       if(opname(i)(1:3).eq.'kin')then                                  2d1s23
        ikin=i                                                          2d1s23
        go to 113                                                       2d1s23
       end if                                                           2d1s23
      end do                                                            2d1s23
  113 continue                                                          2d1s23
      call psioppsi(iwavedat,nroot,ikin,1d0,ixmt,iosym,n2e,             7d27s21
     $      i2eop(1,ikin),fact2e,iwavedat,nroot,bc(isump),mdon,mdoop,    7d27s21
     $      ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,nvirt,maxbx,   7d27s21
     $      maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,npadddi,bc,ibc)      11d10s22
      if(lprint)then                                                    2d1s23
       call pdump(bc(isump),nroot,iwavedat(1),wname(1:nnn),             2d1s23
     $     'mp w/o mass')                                               2d1s23
      end if                                                            2d1s23
      jmassp=imassp                                                     2d1s23
      do i=0,nroot-1                                                    2d1s23
       jtmp=isump+nroot*i                                               2d1s23
       do j=0,i                                                         2d1s23
        bc(jmassp+j)=bc(jtmp+j)*xmassi2                                 2d1s23
        bc(jtmp+j)=bc(jtmp+j)*xmassi                                    2d1s23
       end do                                                           2d1s23
       jmassp=jmassp+i+1                                                2d1s23
      end do                                                            2d1s23
      if(lprint)then                                                    2d1s23
       call pdump(bc(isump),nroot,iwavedat(1),wname(1:nnn),             2d1s23
     $     'mp in 1/cm')                                                2d1s23
      end if                                                            2d1s23
c
c     Ax^2+Ay^2 if not atom
c
      if(nlzz.ne.6)then                                                 2d1s23
       fact2e=-2d0                                                       5d21s21
       n2e=0                                                             2d1s23
       do isb=1,nsymb
        nn=nbasdws(isb)*nbasdws(isb)                                     2d1s23
        do iz=ixsum(isb,1),ibcoff-1                                      2d1s23
         bc(iz)=0d0                                                      2d1s23
        end do                                                           2d1s23
       end do                                                            2d1s23
       do i=0,nroot*nroot-1
        bc(isump+i)=0d0
       end do
       do i=1,nop                                                        2d1s23
        if(opname(i)(1:4).eq.'Ax^2'.or.opname(i)(1:4).eq.'Ay^2')then    2d1s23
         n2e=n2e+1                                                       2d1s23
         i2eops(1,n2e)=n2e+1
         i2eops(2,n2e)=n2e+1                                             2d1s23
         iosyms(n2e+1)=iosym(i2eop(1,i))                                 2d1s23
         do isb=1,nsymb                                                  2d1s23
          nn=nbasdws(isb)*nbasdws(isb)                                     2d1s23
          do j=0,nn-1                                                    2d1s23
           bc(ixsum(isb,1)+j)=bc(ixsum(isb,1)+j)+bc(ixmt(isb,i)+j)       2d1s23
          end do                                                         2d1s23
          jsb=multh(isb,iosyms(n2e+1))                                   2d1s23
          nm=nbasdws(isb)*nbasdws(jsb)                                   2d1s23
          do j=0,nm-1                                                    2d1s23
           bc(ixsum(isb,n2e+1)+j)=bc(ixmt(isb,i2eop(1,i))+j)             2d1s23
          end do                                                         2d1s23
         end do                                                          2d1s23
        end if                                                           2d1s23
       end do
       call psioppsi(iwavedat,nroot,1,1d0,ixsum,iosyms,n2e,              2d1s23
     $      i2eops,fact2e,iwavedat,nroot,bc(isump),mdon,mdoop,          2d1s23
     $      ixw1,ixw2,ism,irel,irefo,norb,nbasdws,idoubo,nvirt,maxbx,   7d27s21
     $      maxbxd,srh,sr2,multh,nsymb,idum,ncsf,0,npadddi,bc,ibc)      11d10s22
       if(lprint)then                                                    2d1s23
        call pdump(bc(isump),nroot,iwavedat(1),wname(1:nnn),             2d1s23
     $        'Ax^2+Ay^2')                                              3d30s22
       end if                                                           2d1s23
      end if                                                            2d1s23
       ibcoff=isump                                                      5d24s21
      ibcoff=ibcoffo                                                    5d26s22
      return                                                            5d20s21
      end                                                               5d20s21
