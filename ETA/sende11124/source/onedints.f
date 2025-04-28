c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine onedints(itype,iak,iab,lb,zetb,xnb,xb,yb,zb,
     $     lk,zetk,xnk,xk,yk,zk,idxyz,jb,jb1,jb2,jb3,nneedzm,           5d17s16
     $     nneedz,                                                      5d12s16
     $     nneedz2,                                                     4d4s16
     $     nneedz2m,npass,idrb,idrk,iapair,natom,idorel,xcart,atnum,
     $     iatom,ixyz,nowpro,ishell,idersign,lbodc,bc,ibc)              11d10s22
      implicit real*8 (a-h,o-z)
      dimension idxyz(3),xcart(3,1),atnum(3,1)
      logical ldeb,lbodc,ldec                                                6d3s22
c
c     compute derivative integrals for a particular shell.
c     idrb,idrk and iapair control signs due to symmetrically
c     equivalent atoms
c     if itype=0, compute overlap and k.e., returning der of bra overlap
c     in jb, der of ket overlap in jb1, der of bra k.e. in jb2, and
c     der of ket k.e. in jb3,
c     otherwise, return full derivative of appropriate potential.
c     if lbodc is true, compute derivative of bra and ket of overlap    6d3s22
c     (and kin, if relativisitc)                                        6d3s22
c
      include "common.store"
      save
      data icall/0/
      icall=icall+1
      ldeb=.false.                                                       4d27s16
      ldec=.false.
      nbbb=2*lb+1
      nkkk=2*lk+1
      if(ldeb)write(6,*)('natom = '),natom,('type '),itype
      if(itype.lt.1)then                                                2d4s16
       jb=ibcoff                                                        4d5s16
       jb1=jb+nneedz                                                    4d5s16
       jb2=jb1+nneedz                                                   4d5s16
       jb3=jb2+nneedz                                                   4d5s16
       ibcoff=jb3+nneedz                                                4d5s16
       npass=4                                                          6d23s16
       if(lbodc)then                                                    6d3s22
        jb4=ibcoff                                                      6d3s22
        ibcoff=jb4+nneedz                                               6d3s22
        if(idorel.eq.0)then                                             6d3s22
         npass=5                                                        6d3s22
        else                                                            6d3s22
         npass=6                                                        6d3s22
         jb5=ibcoff                                                     6d3s22
         ibcoff=jb5+nneedz                                              6d3s22
        end if                                                          6d3s22
       end if                                                           6d3s22
       call enough('onedints.  1',bc,ibc)
       do iz=jb,ibcoff-1                                                4d29s22
        bc(iz)=0d0                                                      4d29s22
       end do                                                           4d29s22
       if(iak.eq.iatom.or.iab.eq.iatom)then                             4d4s16
        if(iak.eq.iab)then                                              4d4s16
         if(lbodc)then                                                  6d3s22
          call oneider(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,        2d4s16
     $       ib1,ib2,idxyz(1),idxyz(2),idxyz(3),idxyz(1),idxyz(2),      6d3s22
     $        idxyz(3),ldec,bc,ibc)                                     11d10s22
          if(idorel.eq.0)then                                           6d3s22
           if(idrk.eq.idrb)then                                          6d3s22
            do iz=0,nneedzm                                             6d3s22
             bc(jb4+iz)=bc(ib1+iz)                                      6d3s22
            end do                                                      6d3s22
           else                                                          6d3s22
            do iz=0,nneedzm                                             6d3s22
             bc(jb4+iz)=-bc(ib1+iz)                                     6d3s22
            end do                                                      6d3s22
           end if                                                        6d3s22
          else                                                          6d3s22
           if(idrk.eq.idrb)then                                          6d3s22
            do iz=0,nneedzm                                             6d3s22
             bc(jb4+iz)=bc(ib1+iz)                                      6d3s22
             bc(jb5+iz)=bc(ib2+iz)                                      6d3s22
            end do                                                      6d3s22
           else                                                          6d3s22
            do iz=0,nneedzm                                             6d3s22
             bc(jb4+iz)=-bc(ib1+iz)                                     6d3s22
             bc(jb5+iz)=-bc(ib2+iz)                                     6d3s22
            end do                                                      6d3s22
           end if                                                        6d3s22
          end if                                                        6d3s22
         end if                                                         6d3s22
         if(ldeb)write(6,*)('oneider the first ')
         call oneider(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,        2d4s16
     $       ib1,ib2,0,0,0,idxyz(1),idxyz(2),idxyz(3),ldec,bc,ibc)      11d10s22
         if(idrk.eq.2)then                                              5d6s16
          do iz=0,nneedzm                                                4d5s16
           bc(jb1+iz)=-bc(ib1+iz)                                         4d5s16
           bc(jb3+iz)=-bc(ib2+iz)                                         4d5s16
          end do                                                         4d5s16
         else
          do iz=0,nneedzm                                                4d5s16
           bc(jb1+iz)=bc(ib1+iz)                                         4d5s16
           bc(jb3+iz)=bc(ib2+iz)                                         4d5s16
          end do                                                         4d5s16
         end if
         ibcoff=ib1                                                     4d5s16
         if(ldeb)write(6,*)('oneider the second ')
         call oneider(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,ib1d,   2d4s16
     $         ib2d,idxyz(1),idxyz(2),idxyz(3),0,0,0,ldec,bc,ibc)       11d10s22
         if(idrb.eq.1)then                                              5d6s16
         do iz=0,nneedzm
          bc(jb+iz)=bc(ib1d+iz)
          bc(jb2+iz)=bc(ib2d+iz)
         end do
         else
          do iz=0,nneedzm                                               2d4s16
           bc(jb+iz)=-bc(ib1d+iz)
           bc(jb2+iz)=-bc(ib2d+iz)
          end do
         end if
         ibcoff=ib1d                                                    3d10s16
        else if(iak.eq.iatom)then                                       2d4s16
         if(ldeb)write(6,*)('oneider the third ')
         call oneider(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,ib1,    2d4s16
     $         ib2,0,0,0,idxyz(1),idxyz(2),idxyz(3),ldec,bc,ibc)        11d10s22
         if(idrk.eq.2)then                                              5d6s16
          do iz=0,nneedzm                                               2d3s16
           bc(jb1+iz)=-bc(ib1+iz)
           bc(jb3+iz)=-bc(ib2+iz)
           bc(jb+iz)=0d0
           bc(jb2+iz)=0d0
          end do
         else
          do iz=0,nneedzm                                               2d3s16
           bc(jb1+iz)=bc(ib1+iz)
           bc(jb3+iz)=bc(ib2+iz)
           bc(jb+iz)=0d0
           bc(jb2+iz)=0d0
          end do
         end if
         ibcoff=ib1
        else if(iab.eq.iatom)then                                       4d4s16
         if(ldeb)write(6,*)('oneider the forth ')
         call oneider(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,ib1,    2d4s16
     $        ib2,idxyz(1),idxyz(2),idxyz(3),0,0,0,ldec,bc,ibc)         11d10s22
         if(idrb.eq.2)then                                              5d6s16
          do iz=0,nneedzm                                               2d3s16
           bc(jb+iz)=-bc(ib1+iz)
           bc(jb2+iz)=-bc(ib2+iz)
           bc(jb1+iz)=0d0
           bc(jb3+iz)=0d0
          end do
         else
          do iz=0,nneedzm                                               2d3s16
           bc(jb+iz)=bc(ib1+iz)
           bc(jb2+iz)=bc(ib2+iz)
           bc(jb1+iz)=0d0
           bc(jb3+iz)=0d0
          end do
         end if
        end if
       else                                                             6d5s14
        do iz=0,nneedz2m
         bc(jb+iz)=0d0
         bc(jb2+iz)=0d0
        end do                                                          4d5s16
       end if                                                           6d4s14
      else                                                              5d7s10
       if(ldeb)write(6,*)('natom = '),natom
       npass=1                                                          5d7s10
       ia=itype                                                         2d4s16
       if(ia.le.natom)then                                              2d2s16
        iau=ia                                                          2d2s16
       else if(ia.le.natom*2)then                                       2d2s16
        iau=ia-natom                                                    2d2s16
       else if(ia.le.natom*3)then                                       2d2s16
        iau=ia-natom*2                                                  2d2s16
       else
        iau=ia-natom*3
       end if                                                           2d2s16
       dersign=1d0                                                      2d3s16
       jau=iau                                                          2d3s16
       if(ldeb)then
        write(6,*)('for ia = '),ia,iau
        write(6,*)('iatom, iapair '),iatom,iapair
       end if
       if(iau.eq.iapair)then                                            2d4s16
        jau=iatom                                                       2d3s16
        if(idersign.eq.2)dersign=-1d0                                   2d3s16
       end if                                                           2d3s16
       if(ldeb)then
        write(6,*)('jau, iatom '),jau,iatom
       end if
       if(jau.eq.iatom.or.iab.eq.iatom.or.iak.eq.iatom)then
        if(ldeb)write(6,*)('inside block ')
        jb=-1
        if(idorel.eq.0)then                                             8d20s15
         if(iab.eq.iak.and.iab.eq.iatom)then
          arg=1d0                                                       5d6s16
          if(idrb.eq.2)arg=-1d0                                         5d6s16
          call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $         jb,idxyz(1),idxyz(2),idxyz(3),0,0,0,0,0,0,0,arg,ldec,    6d30s16
     $         xcart(1,ia),bc,ibc)                                      11d10s22
          if(ldeb)write(6,*)('nattracd the first: '),bc(jb)
          arg=+1d0
          if(idrk.eq.2)arg=-1d0                                         5d6s16
          call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $         ib1d,0,0,0,idxyz(1),idxyz(2),idxyz(3),0,0,0,0,arg,ldec,  6d30s16
     $         xcart(1,ia),bc,ibc)                                      11d10s22
          if(ldeb)write(6,*)('nattracd the 2nd: '),bc(ib1d)
          do iz=0,nneedzm                                               2d3s16
           bc(jb+iz)=bc(jb+iz)+bc(ib1d+iz)
          end do
          ibcoff=ib1d                                                   3d10s16
         else if(iab.eq.iatom)then
          arg=+1d0                                                      2d4s16
          if(idrb.eq.2)arg=-1d0                                         5d6s16
          call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $          jb,idxyz(1),idxyz(2),idxyz(3),0,0,0,0,0,0,0,arg,ldec,   6d30s16
     $         xcart(1,ia),bc,ibc)                                      11d10s22
          if(ldeb)write(6,*)('nattracd the third: '),bc(jb)
         else if(iak.eq.iatom)then
          arg=+1d0
          if(idrk.eq.2)arg=-1d0                                         5d6s16
          call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,       8d20s15
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $          jb,0,0,0,idxyz(1),idxyz(2),idxyz(3),0,0,0,0,arg,ldec,   6d30s16
     $         xcart(1,ia),bc,ibc)                                      11d10s22
          if(ldeb)write(6,*)('nattracd the fourth: '),bc(jb)
         end if
         if(ldeb)write(6,*)('jau, iatom '),jau,iatom
         if(jau.eq.iatom)then
          if(ldeb)write(6,*)('in block '),jb
          ib1sav=-1
          if(jb.gt.0)then
           ib1sav=jb
          end if
c
c     iatom is always a symmetry unique atom, but ia=iau points to all
c     atoms. if iau points to symmetry redundant atom, jau has been
c     set to the symmetry unique atom to pass the iatom=ia test.
c                                                                       4d25s16
          factna=xcart(ixyz,iau)                                        4d25s16
          call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $          jb,0,0,0,0,0,0,ixyz,1,0,0,1d0,ldec,xcart(1,ia),bc,ibc)  11d10s22
          if(ldeb)write(6,*)('nattracd the fifth: '),bc(jb)
          call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $          jb2,0,0,0,0,0,0,ixyz,0,1,0,1d0,ldec,                    6d30s16
     $         xcart(1,ia),bc,ibc)                                      11d10s22
          if(ldeb)write(6,*)('nattracd the fifthb: '),bc(jb2)
          call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $          jb3,0,0,0,0,0,0,ixyz,0,0,1,1d0,ldec,                    6d30s16
     $         xcart(1,ia),bc,ibc)                                      11d10s22
          if(ldeb)write(6,*)('nattracd the fifthc: '),bc(jb3),jb,
     $         jb+nneedzm
          do iz=0,nneedzm                                               6d27s16
           bc(jb+iz)=bc(jb+iz)+bc(jb2+iz)+bc(jb3+iz)                    6d29s16
          end do                                                        6d27s16
          ibcoff=jb2                                                    6d29s16
          if(ldeb)write(6,*)('take dersign to be '),dersign,idersign
          if(dersign.gt.0d0)then                                        2d3s16
           do iz=0,nneedzm                                              2d3s16
            bc(jb+iz)=-bc(jb+iz)                                        2d3s16
           end do                                                       2d3s16
          end if                                                        2d3s16
          if(ib1sav.gt.0)then
           do iz=0,nneedzm                                              2d3s16
            bc(ib1sav+iz)=bc(ib1sav+iz)+bc(jb+iz)
           end do
           ibcoff=jb                                                    3d10s16
           jb=ib1sav                                                    4d5s16
          end if
         end if
        else                                                            8d20s15
         iderx=0                                                        8d20s15
         idery=0                                                        8d20s15
         iderz=0                                                        8d20s15
         if(ia.gt.natom)then                                            6d24s16
          if(ia.le.natom*2)then                                          8d20s15
           iderx=1
          else if(ia.le.natom*3)then                                     8d20s15
           idery=1                                                       8d20s15
          else
           iderz=1
          end if                                                         8d20s15
         end if                                                         6d24s16
         if(ldeb)write(6,*)('iderxyz: '),iderx,idery,iderz
         if(ldeb)write(6,*)('idxyz: '),idxyz
         zm=-atnum(1,iau)                                               8d20s15
         idxp=iderx+idxyz(1)
         idyp=idery+idxyz(2)
         idzp=iderz+idxyz(3)
         if(iab.eq.iak.and.iab.eq.iatom)then                            2d4s16
          arg=+zm                                                       5d6s16
          if(idrb.eq.2)arg=-zm                                          5d6s16
          call derid(lb,zetb,xnb,xb,yb,zb,0,lk,zetk,xnk,xk,yk,zk,0,     2d4s16
     $         0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),
     $         xcart(3,iau),0,0,atnum(2,iau),atnum(3,iau),xcart(1,iau),
     $         xcart(2,iau),xcart(3,iau),0,dum,1,idum,jb,
     $         idxp,idyp,idzp,iderx,idery,iderz,
     $         0,0,0,0,0,0,1,ldeb,arg,bc,ibc)                           11d14s22
          arg=+zm                                                       5d6s16
          if(idrk.eq.2)arg=-zm                                          5d6s16
          call derid(lb,zetb,xnb,xb,yb,zb,0,lk,zetk,xnk,xk,yk,zk,0,     2d4s16
     $          0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),
     $          xcart(3,iau),0,0,atnum(2,iau),atnum(3,iau),xcart(1,iau),
     $          xcart(2,iau),xcart(3,iau),0,dum,1,idum,ib1d,
     $          iderx,idery,iderz,idxp,idyp,idzp,
     $          0,0,0,0,0,0,1,ldeb,arg,bc,ibc)                          11d14s22
          do iz=0,nneedzm                                               2d3s16
           bc(jb+iz)=bc(jb+iz)+bc(ib1d+iz)
          end do
          ibcoff=ib1d                                                   3d10s16
         else if(iab.eq.iatom)then
          arg=+zm                                                       5d6s16
          if(idrb.eq.2)arg=-zm                                          5d6s16
          call derid(lb,zetb,xnb,xb,yb,zb,0,lk,zetk,xnk,xk,yk,zk,0,     2d4s16
     $          0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),
     $          xcart(3,iau),0,0,atnum(2,iau),atnum(3,iau),xcart(1,iau),
     $          xcart(2,iau),xcart(3,iau),0,dum,1,idum,jb,
     $          idxp,idyp,idzp,iderx,idery,iderz,                       2d4s16
     $         0,0,0,0,0,0,1,ldeb,arg,bc,ibc)                           11d14s22
         else if(iak.eq.iatom)then
          arg=+zm                                                       5d6s16
          if(idrk.eq.2)arg=-zm                                          5d6s16
          call derid(lb,zetb,xnb,xb,yb,zb,0,lk,zetk,xnk,xk,yk,zk,0,     2d4s16
     $          0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),
     $          xcart(3,iau),0,0,atnum(2,iau),atnum(3,iau),xcart(1,iau),
     $          xcart(2,iau),xcart(3,iau),0,dum,1,idum,jb,
     $          iderx,idery,iderz,idxp,idyp,idzp,
     $         0,0,0,0,0,0,1,ldeb,arg,bc,ibc)                           11d14s22
         end if
         if(jau.eq.iatom)then                                           2d3s16
          ib1sav=-1
          if(jb.gt.0)then
           ib1sav=jb
          end if
          arg=+2d0*zm*dersign                                           5d6s16
c
c     we use erird instead of derid here because of numerical
c     issues differentiating the extremely tight s function describing
c     the nucleus. erird does the cartesian integrals analtyically via
c     recursion.
c     we do not use erird for the other derivatives, for it can not do
c     mixed ders on a particular center.
c
          call erird(lb,zetb,xnb,xb,yb,zb,                              8d14s24
     $               lk,zetk,xnk,xk,yk,zk,                              8d14s24
     $               0,atnum(2,iau),atnum(3,iau),xcart(1,iau),          8d14s24
     $                              xcart(2,iau),xcart(3,iau),          8d14s24
     $               0,atnum(2,iau),atnum(3,iau),xcart(1,iau),          8d14s24
     $                              xcart(2,iau),xcart(3,iau),          8d14s24
     $         jb,0,arg,iderx,idery,iderz,iderx,idery,iderz,idxyz(1),   8d14s24
     $         idxyz(2),idxyz(3),0,0,0,bc,ibc)                          8d14s24
          if(ib1sav.gt.0)then
           do iz=0,nneedzm                                              2d3s16
            bc(ib1sav+iz)=bc(ib1sav+iz)+bc(jb+iz)
           end do
           ibcoff=jb                                                    3d10s16
           jb=ib1sav
          end if
         end if
        end if
       else                                                             6d5s14
        jb=ibcoff                                                        6d5s14
        ibcoff=jb+nneedz                                                6d5s14
        call enough('onedints.  2',bc,ibc)
        do ineedz=0,nneedzm                                             2d3s16
         bc(jb+ineedz)=0d0                                              6d5s14
        end do                                                          6d5s14
       end if
      end if                                                            5d7s10
      return
      end
