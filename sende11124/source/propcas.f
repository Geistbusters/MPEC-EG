c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine propcas(natom,ngaus,ibdat,isym,iapair,ibstor,isstor,   12d19s19
     $     idorel,ascale,multh,nbasp,iorb,ixmtr,noc,islz,nlzz,iflag,    4d3s21
     $     iorbsym,idoorbsym,iorbsymz,bc,ibc,epssym)                    9d19s23
      implicit real*8 (a-h,o-z)                                         12d19s19
      parameter (id=100)                                                12d19s19
      character*10 opname(id)                                           12d19s19
      dimension ipt(id),npt(id),opdata(id),iopdata(7,id),iosym(id),     5d27s21
     $     ioprt(id),i2eop(6,id),ixmt(8,id),isym(3,8),iapair(3,*),      5d20s21
     $     ibstor(*),isstor(*),iso(8),iorb(8),multh(8,8),opnc(id),      12d19s19
     $     nbasp(*),ixmtr(8,6),noc(*),islz(3),iorbsym(*),iorbsymz(*)    4d19s21
      include "common.store"                                            12d19s19
      include "common.mrci"                                             12d19s19
      include "common.hf"                                               12d19s19
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      call oplist(id,opname,ipt,npt,opdata,iopdata,iosym,ioprt,i2eop,   12d17s19
     $     nop,opnc,nlzz)                                               12d31s19
      nbb=0                                                             12d15s19
      if(idoorbsym.eq.2)then                                            4d19s21
       do isb=1,nsymb                                                   4d19s21
        ixmtr(isb,1)=ibcoff                                             4d19s21
        ibcoff=ixmtr(isb,1)+noc(isb)*noc(isb)                           4d19s21
       end do                                                           4d19s21
       if(nlzz.ne.2)then                                                4d19s21
        do isb=1,nsymb                                                   4d19s21
         ixmtr(isb,2)=ibcoff                                             4d19s21
         ibcoff=ixmtr(isb,2)+noc(isb)*noc(isb)                           4d19s21
        end do                                                           4d19s21
       end if                                                           4d19s21
      else                                                              4d19s21
       do i=1,nop                                                        12d31s19
        do isb=1,nsymb                                                   12d22s19
         isk=multh(isb,iosym(i))                                         12d22s19
         ixmtr(isb,i)=ibcoff                                             12d22s19
         ibcoff=ixmtr(isb,i)+noc(isb)*noc(isk)                           12d22s19
        end do                                                           12d22s19
       end do                                                            12d22s19
      end if                                                            4d19s21
      ibcoffo=ibcoff                                                    12d22s19
      if(idorel.eq.0)then                                               1d10s19
       ncomp=1                                                          1d10s19
      else                                                              1d10s19
       ncomp=2                                                          1d10s19
      end if                                                            1d10s19
      do i=1,nop                                                        12d15s19
       do isb=1,nsymb                                                   12d15s19
        isk=multh(isb,iosym(i))                                          12d15s19
        ixmt(isb,i)=ibcoff                                              12d15s19
        ibcoff=ixmt(isb,i)+nbasp(isb)*nbasp(isk)*ncomp*ncomp            1d10s19
        nbb=nbb+nbasp(isb)*nbasp(isk)*ncomp*ncomp                       1d10s19
       end do                                                           12d15s19
      end do                                                            12d15s19
      ipass=0                                                           3d15s21
 2727 continue                                                          3d15s21
      if(ipass.gt.1)then                                                3d15s21
       write(6,*)('failure in propcas ')                                3d15s21
       call dws_synca                                                   3d15s21
       call dws_finalize                                                3d15s21
       stop                                                             3d15s21
      end if                                                            3d15s21
      call enough('propcas.  1',bc,ibc)
      do i=ixmt(1,1),ibcoff-1                                           12d15s19
       bc(i)=0d0                                                        12d15s19
      end do                                                            12d15s19
      idwsdeb=0                                                         12d16s19
      call parap(natom,ngaus,ibdat,ixmt,isym,iapair,ibstor,isstor,iso,  12d16s19
     $     idwsdeb,idorel,ascale,ipt,npt,opdata,iopdata,iosym,nop,multh,12d16s19
     $     nbb,nbasp,nbasdws,iorb,opname,1,bc,ibc)                      11d9s22
      if(idoorbsym.ne.0)then                                            4d3s21
c                                                                       4d3s21
c     set symmetry qn for orbitals.                                     4d3s21
c     idoorbsym=1, just set qn's and proceed.
c     idoorbsym=2, we are in primitive basis, so set qns and return     4d19s21
c     lz^2(and perhaps l^2) to feed into contracg.                      4d19s21
c                                                                       4d19s21
       nbad=0                                                           4d17s21
       if(nlzz.eq.2)then                                                4d19s21
        isav2=1                                                         4d19s21
       else                                                             4d19s21
        isav2=2                                                         4d19s21
       end if                                                           4d19s21
       do isb=1,nsymb                                                   4d3s21
        if(idoorbsym.eq.2)then                                          4d19s21
         nuse=noc(isb)                                                  4d19s21
        else                                                            4d19s21
         nuse=nbasdws(isb)                                              4d19s21
        end if                                                          4d19s21
        idig=ibcoff                                                     4d3s21
        ibcoff=idig+nuse                                                4d19s21
        if(nlzz.gt.2)then                                               4d19s21
         idigz=ibcoff                                                   4d19s21
         ibcoff=idigz+nuse                                              4d19s21
        end if                                                          4d19s21
        call enough('propcas.  2',bc,ibc)
        do i=idig,ibcoff-1                                              4d3s21
         bc(i)=0d0                                                      4d3s21
        end do                                                          4d3s21
        do i=1,nop                                                      4d3s21
         if(opname(i)(3:4).eq.'^2')then                                 4d4s21
          if(idoorbsym.eq.2)then                                        4d19s21
           do j=0,nuse*nuse-1                                           4d19s21
            bc(ixmtr(isb,isav2)+j)=bc(ixmtr(isb,isav2)+j)               4d19s21
     $           +bc(ixmt(isb,i)+j)                                     4d19s21
           end do                                                        4d17s21
          end if                                                        4d19s21
          do j=0,nuse-1                                                 4d19s21
           jj=j*(nuse+1)                                                4d19s21
           bc(idig+j)=bc(idig+j)+bc(ixmt(isb,i)+jj)                     4d4s21
          end do                                                        4d4s21
          if(nlzz.gt.2.and.opname(i)(2:2).eq.'z')then                   4d19s21
           do j=0,nuse-1                                                4d19s21
            jj=j*(nuse+1)                                               4d19s21
            bc(idigz+j)=bc(ixmt(isb,i)+jj)                              4d19s21
           end do                                                        4d4s21
           if(idoorbsym.eq.2)then                                        4d19s21
            do j=0,nuse*nuse-1                                           4d19s21
             bc(ixmtr(isb,1)+j)=bc(ixmtr(isb,1)+j)                      4d19s21
     $           +bc(ixmt(isb,i)+j)                                     4d19s21
            end do                                                        4d17s21
           end if                                                        4d19s21
          end if                                                        4d19s21
         end if                                                         4d4s21
        end do                                                          4d3s21
        idigu=idig                                                      4d17s21
c     nlzz=2, diagonal is lambda squared
c     nlzz=6 diagonal is l*(l+1)=l^2+l=(l+1/2)^2-1/4
        do i=0,nuse-1                                                   4d19s21
         if(nlzz.eq.2)then                                              4d4s21
          xgot=sqrt(abs(bc(idigu+i)))                                    4d4s21
         else                                                           4d4s21
          xgot=sqrt(abs(bc(idigu+i)+0.25d0))-0.5d0                       4d4s21
          zgot=sqrt(abs(bc(idigz+i)))                                   4d19s21
         end if                                                         4d4s21
         lqn=nint(xgot)                                                 4d4s21
         if(nlzz.eq.2)then                                              4d15s21
          test=dfloat(lqn*lqn)                                          4d15s21
         else                                                           4d15s21
          test=dfloat(lqn*(lqn+1))                                      4d15s21
          lzqn=nint(zgot)                                               4d19s21
          testz=dfloat(lzqn*lzqn)                                       4d19s21
         end if                                                         4d15s21
         diff=test-bc(idig+i)                                           4d15s21
         if(nlzz.gt.2)then                                              4d19s21
          diff=max(abs(diff),abs(testz-bc(idigz+i)))                    4d19s21
         end if                                                         4d19s21
         if(abs(diff).gt.epssym)then                                    9d19s23
          write(6,*)('for symmetry '),isb                                 4d3s21
          write(6,*)('bad qn assignment: '),lqn,test,bc(idig+i),epssym  9d19s23
          nbad=nbad+1                                                   4d17s21
         end if                                                         4d4s21
         ibc(iorbsym(isb)+i)=lqn                                        4d4s21
         if(nlzz.eq.2.and.nsymb.eq.8)then                               10d1s21
          if(isb.eq.2.or.isb.eq.3.or.isb.eq.5.or.isb.eq.8)then          10d1s21
           ibc(iorbsym(isb)+i)=ibc(iorbsym(isb)+i)+50                   10d1s21
          end if                                                        10d1s21
         end if                                                         10d1s21
         if(nlzz.gt.2)then                                              4d19s21
          ibc(iorbsymz(isb)+i)=lzqn                                     4d19s21
         end if                                                         4d19s21
        end do                                                          4d4s21
        ibcoff=idig                                                     4d17s21
       end do                                                           4d3s21
       if(nbad.gt.0)then                                                4d17s21
        call dws_synca                                                  4d17s21
        call dws_finalize                                               4d17s21
        stop 'propcas'                                                  10d7s22
       end if                                                           4d17s21
       if(idoorbsym.eq.2)then                                           4d19s21
        ibcoff=ibcoffo                                                  4d19s21
        return                                                          4d19s21
       end if                                                           4d19s21
      end if                                                            4d3s21
      do i=1,nop                                                        12d31s19
       do isb=1,nsymb                                                   12d20s19
        isk=multh(isb,iosym(i))                                         12d22s19
        do icol=0,noc(isk)-1                                            12d22s19
         iad1=ixmt(isb,i)+nbasdws(isb)*icol                             12d22s19
         iad2=ixmtr(isb,i)+noc(isb)*icol                                12d22s19
         do irow=0,noc(isb)-1                                           12d22s19
          bc(iad2+irow)=bc(iad1+irow)                                   12d22s19
         end do                                                         12d22s19
        end do                                                          12d22s19
        if(min(noc(isb),noc(isk)).gt.0)then                             1d11s20
         if(iflag.ne.0)then
          write(6,*)('matrix no. '),i,isb,isk,opname(i)
          call prntm2(bc(ixmtr(isb,i)),noc(isb),noc(isk),noc(isb))
         end if
         thrs=1d-16                                                      4d25s21
         szx=0d0                                                        5d3s21
         nszx=0                                                         5d3s21
         do ii=0,noc(isb)*noc(isk)-1                                    4d25s21
          orig=bc(ixmtr(isb,i)+ii)                                      4d25s21
          if(abs(bc(ixmtr(isb,i)+ii)).lt.thrs)then                      4d25s21
           szx=max(szx,abs(bc(ixmtr(isb,i)+ii)))                        5d3s21
           bc(ixmtr(isb,i)+ii)=0d0                                      4d25s21
           itry=0
          else                                                          4d25s21
           try=orig*orig                                                4d25s21
           itry=nint(try)                                               4d25s21
           diff=try-dfloat(itry*itry)                                   4d25s21
           if(abs(diff).lt.thrs)then                                    4d25s21
            szx=max(szx,abs(diff))                                      5d3s21
            if(orig.gt.0d0)then                                         4d25s21
             bc(ixmtr(isb,i)+ii)=sqrt(dfloat(itry))                     4d25s21
            else                                                        4d25s21
             bc(ixmtr(isb,i)+ii)=-sqrt(dfloat(itry))                    4d25s21
            end if                                                      4d25s21
           else                                                         4d25s21
            itry=-1                                                     4d25s21
           end if                                                       4d25s21
          end if                                                        4d25s21
         end do                                                         4d25s21
         if(szx.gt.thrs)then                                            5d3s21
          write(6,*)('warning: purification in propcas encountered ')   5d3s21
          write(6,*)('errors as large as '),szx                         5d3s21
         end if                                                         5d3s21
         if(noc(isb).eq.noc(isk))then                                   3d15s21
          size=0d0                                                      3d15s21
          do k=1,noc(isk)-1                                             3d15s21
           do j=0,k-1                                                   3d15s21
            jk=ixmtr(isb,i)+j+noc(isb)*k                                3d15s21
            kj=ixmtr(isb,i)+k+noc(isb)*j                                3d15s21
            size=size+bc(jk)**2+bc(kj)**2                               3d15s21
           end do                                                       3d15s21
          end do                                                        3d15s21
         end if                                                         3d15s21
        end if                                                          1d11s20
       end do                                                           12d20s19
      end do                                                            12d20s19
      nhere=nlzz/2                                                      12d31s19
      do i=1,nhere                                                      12d31s19
       islz(i)=iosym(i)                                                 12d31s19
      end do
      ibcoff=ibcoffo                                                    12d19s19
      return                                                            12d19s19
      end                                                               12d19s19
