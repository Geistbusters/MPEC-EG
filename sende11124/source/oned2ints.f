c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine oned2ints(itype,iak,iab,lb,zetb,xnb,xb,yb,zb,          6d20s16
     $     lk,zetk,xnk,xk,yk,zk,idxyz,jb,jb1,jb2,jb3,jb4,jb5,nneedzm,   6d20s16
     $     nneedz,                                                      5d12s16
     $     nneedz2,                                                     4d4s16
     $     nneedz2m,npass,idrb,idrk,iapair,natom,idorel,xcart,atnum,    2d27s23
     $     iatom,ixyz,nowpro,ishell,idersign,bc,ibc)                    11d10s22
      implicit real*8 (a-h,o-z)
      dimension idxyz(3),xcart(3,1),atnum(3,1),design(2)                4d22s22
      logical ldeb,lcs,lds                                              2d27s23
c
c     compute derivative integrals for a particular shell.
c     idrb,idrk and iapair control signs due to symmetrically
c     equivalent atoms
c     if itype=0, compute overlap and k.e., returning double der of bra 6d20s16
c     overlap in jb, double der of ket overlap in jb1, der of bra       6d20s16
c     overlap der of ket in jb2, double der of bra k.e. in jb3, double  6d20s16
c     der of ket k.e. in jb4 and der of bra k.e. der of ket in jb5.     6d20s16
c     otherwise, return 2nd derivative of appropriate potential in jb.
c
      include "common.store"
      common/xcom/iadd1,iadd2,iadd3
      save
      data icall/0/
      data design/2*1d0/                                                3d2s23
      icall=icall+1
      ldeb=.false.                                                      4d27s16
      iadd1=0
      iadd2=0
      iadd3=0
       nbbb=2*lb+1
       nkkk=2*lk+1
      if(ldeb)write(6,*)('ldeb is true in oned2ints')
      if(itype.lt.1)then                                                2d4s16
       jb=ibcoff                                                        4d5s16
       jb1=jb+nneedz                                                    4d5s16
       jb2=jb1+nneedz                                                   4d5s16
       jb3=jb2+nneedz                                                   4d5s16
       jb4=jb3+nneedz                                                   6d20s16
       jb5=jb4+nneedz                                                   6d20s16
       ibcoff=jb5+nneedz                                                6d20s16
       call enough('oned2ints.  1',bc,ibc)
       do iz=jb,ibcoff-1                                                4d22s22
        bc(iz)=0d0                                                      4d22s22
       end do                                                           4d22s22
       npass=6                                                          4d4s16
       if(iab.eq.iatom.or.iab.eq.iapair)then                            4d22s22
        call oneider(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,ib1d,    4d22s22
     $         ib2d,idxyz(1),idxyz(2),idxyz(3),0,0,0,ldeb,bc,ibc)       11d10s22
         do iz=0,nneedzm
          bc(jb+iz)=bc(jb+iz)+bc(ib1d+iz)                                4d22s22
          bc(jb3+iz)=bc(jb3+iz)+bc(ib2d+iz)                              4d22s22
         end do
        ibcoff=ib1d                                                     4d22s22
       end if                                                           4d22s22
       if(iak.eq.iatom.or.iak.eq.iapair)then                            4d22s22
        call oneider(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,         4d22s22
     $       ib1,ib2,0,0,0,idxyz(1),idxyz(2),idxyz(3),ldeb,bc,ibc)      11d10s22
         do iz=0,nneedzm                                                 4d22s22
          bc(jb1+iz)=bc(jb1+iz)+bc(ib1+iz)                               4d22s22
          bc(jb4+iz)=bc(jb4+iz)+bc(ib2+iz)                               4d22s22
         end do                                                          4d22s22
        ibcoff=ib1                                                      4d22s22
       end if                                                           4d22s22
       idxh=idxyz(1)/2                                                  4d22s22
       idyh=idxyz(2)/2                                                  4d22s22
       idzh=idxyz(3)/2                                                  4d22s22
       if((iak.eq.iatom.and.iab.eq.iatom).or.                           4d22s22
     $    (iak.eq.iapair.and.iab.eq.iapair))then                        4d22s22
        call oneider(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,         4d22s22
     $       ib1,ib2,idxh,idyh,idzh,idxh,idyh,idzh,ldeb,bc,ibc)         11d10s22
        factme=2d0                                                      4d22s22
        do iz=0,nneedzm                                                 4d22s22
         bc(jb2+iz)=bc(jb2+iz)+bc(ib1+iz)*factme                        4d22s22
         bc(jb5+iz)=bc(jb5+iz)+bc(ib2+iz)*factme                        4d22s22
        end do                                                          4d22s22
        ibcoff=ib1                                                      4d22s22
       else if((iak.eq.iatom.and.iab.eq.iapair).or.                     4d22s22
     $         (iak.eq.iapair.and.iab.eq.iatom))then                    4d22s22
        if(idrb.eq.idrk)then                                            2d27s23
         call oneider(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,         4d22s22
     $       ib1,ib2,idxh,idyh,idzh,idxh,idyh,idzh,ldeb,bc,ibc)         11d10s22
         factme=2d0*design(idersign)                                     4d22s22
         do iz=0,nneedzm                                                 4d22s22
          bc(jb2+iz)=bc(jb2+iz)+bc(ib1+iz)*factme                        4d22s22
          bc(jb5+iz)=bc(jb5+iz)+bc(ib2+iz)*factme                        4d22s22
         end do                                                          4d22s22
         ibcoff=ib1                                                      4d22s22
        end if                                                          2d27s23
       end if                                                           4d22s22
      else                                                              5d7s10
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
        write(6,*)('iab,iak: '),iab,iak
        write(6,*)('for ia = '),ia,iau,atnum(1,ia)
        write(6,*)('iatom, iapair '),iatom,iapair
       end if
       if(iau.eq.iapair)then                                            2d4s16
        jau=iatom                                                       2d3s16
        if(idersign.eq.2)dersign=-1d0                                   2d3s16
       end if                                                           2d3s16
       if(ldeb)then
        write(6,*)('jau, iatom '),jau,iatom
       end if
       if(jau.eq.iatom.or.iab.eq.iatom.or.iak.eq.iatom.or.iab.eq.iapair 4d22s22
     $      .or.iak.eq.iapair)then                                      4d22s22
        if(ldeb)write(6,*)('inside block ')
        jb=ibcoff                                                       4d22s22
        ibcoff=jb+nneedz                                                4d22s22
        call enough('oned2ints.  2',bc,ibc)
        do iz=jb,ibcoff-1                                               4d22s22
         bc(iz)=0d0                                                     4d22s22
        end do                                                          4d22s22
        if(idorel.eq.0)then                                             8d20s15
         tallysum=0d0
         if(iab.eq.iatom.or.iab.eq.iapair)then                          4d22s22
          arg=1d0                                                       5d6s16
          if(idrb.eq.2)arg=arg*design(idersign)                         3d1s23
          call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $         ib1d,idxyz(1),idxyz(2),idxyz(3),0,0,0,0,0,0,0,arg,ldeb,   4d22s22
     $         xcart(1,ia),bc,ibc)                                      11d10s22
          do i=0,nneedzm                                                4d22s22
           bc(jb+i)=bc(ib1d+i)                                          4d25s22
          end do                                                        4d22s22
          ibcoff=ib1d
         end if                                                         4d22s22
         if(iak.eq.iatom.or.iak.eq.iapair)then                          4d22s22
          arg=+1d0
          if(idrk.eq.2)arg=arg*design(idersign)                         3d1s23
          call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $         ib1d,0,0,0,idxyz(1),idxyz(2),idxyz(3),0,0,0,0,arg,ldeb,
     $         xcart(1,ia),bc,ibc)                                      11d10s22
          do i=0,nneedzm                                                4d22s22
           bc(jb+i)=bc(jb+i)+bc(ib1d+i)                                 4d25s22
          end do                                                        4d22s22
          ibcoff=ib1d                                                   4d22s22
         end if                                                         4d22s22
         idxh=idxyz(1)/2                                                  4d22s22
         idyh=idxyz(2)/2                                                  4d22s22
         idzh=idxyz(3)/2                                                  4d22s22
         if((iak.eq.iatom.and.iab.eq.iatom).or.                           4d22s22
     $    (iak.eq.iapair.and.iab.eq.iapair))then                        4d22s22
          arg=2d0
          if(idrk.eq.2)arg=arg*design(idersign)                         3d1s23
          call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $         ib1d,idxh,idyh,idzh,idxh,idyh,idzh,0,0,0,0,arg,ldeb,      4d22s22
     $         xcart(1,ia),bc,ibc)                                      11d10s22
          do i=0,nneedzm                                                4d22s22
           bc(jb+i)=bc(jb+i)+bc(ib1d+i)                                 4d25s22
          end do                                                        4d22s22
          ibcoff=ib1d                                                      4d22s22
         end if                                                         4d22s22
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
c     there are at least 2 and perhaps 4 contributions. The 2 always
c     present are the 1/r**3 and 1/r**5 parts of the 2nd derivative of  6d27s16
c     the potential energy. The perhaps contriubtions are first ders    6d27s16
c     of the potential and first der of bra and/or ket fcn.             6d27s16
c     4d25s16
          factna=xcart(ixyz,iau)                                        4d25s16
          idxh=idxyz(1)/2                                               6d27s16
          idyh=idxyz(2)/2                                               6d27s16
          idzh=idxyz(3)/2                                               6d27s16
          jc1=-1                                                        4d25s22
          jc2=-1                                                        4d25s22
          jc3=-1                                                        4d25s22
          jd1=-1                                                        4d25s22
          jd2=-1                                                        4d25s22
          jd3=-1                                                        4d25s22
          jb1=-1
          jb2=-1
          jb3=-1
          jb4=-1
          jb5=-1
          lcs=.false.                                                   2d27s23
          if(iab.eq.iatom.or.iab.eq.iapair)then                         4d22s22
           lcs=.true.                                                   2d27s23
           if((idrb.eq.1.and.itype.eq.iapair).or.                       3d1s23
     $        (idrb.eq.2.and.itype.eq.iatom))lcs=.false.                3d1s23
           if(lcs)then                                                  2d27s23
            fff=2d0
            call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $         jc1,idxh,idyh,idzh,0,0,0,ixyz,1,0,0,fff,ldeb,xcart(1,ia),11d10s22
     $          bc,ibc)                                                 11d10s22
            call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $         jc2,idxh,idyh,idzh,0,0,0,ixyz,0,1,0,fff,ldeb,xcart(1,ia),11d10s22
     $          bc,ibc)                                                 11d10s22
            call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $         jc3,idxh,idyh,idzh,0,0,0,ixyz,0,0,1,fff,ldeb,xcart(1,ia),11d10s22
     $          bc,ibc)                                                 11d10s22
           end if                                                       2d27s23
          end if                                                        7d1s16
          lds=.false.                                                   2d27s23
          if(iak.eq.iatom.or.iak.eq.iapair)then                         4d22s22
           lds=.true.                                                   2d27s23
           if((idrk.eq.2.and.itype.eq.iatom).or.                        3d1s23
     $         (idrk.eq.1.and.itype.eq.iapair))lds=.false.               2d27s23
           if(lds)then                                                  2d27s23
            fff=2d0
            call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $         jd1,0,0,0,idxh,idyh,idzh,ixyz,1,0,0,fff,ldeb,xcart(1,ia),11d10s22
     $          bc,ibc)                                                 11d10s22
            call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $         jd2,0,0,0,idxh,idyh,idzh,ixyz,0,1,0,fff,ldeb,xcart(1,ia),11d10s22
     $          bc,ibc)                                                 11d10s22
            call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $         jd3,0,0,0,idxh,idyh,idzh,ixyz,0,0,1,fff,ldeb,xcart(1,ia),11d10s22
     $          bc,ibc)                                                 11d10s22
           end if                                                       2d27s23
          end if                                                        7d1s16
          call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $          jb,0,0,0,0,0,0,ixyz,2,0,0,1d0,ldeb,xcart(1,ia),bc,ibc)  11d10s22
          if(lcs)then                                                   2d27s23
           do iz=0,nneedzm                                              7d1s16
            bc(jb+iz)=bc(jb+iz)-(bc(jc1+iz)+bc(jc2+iz)+bc(jc3+iz))      7d1s16
           end do
          end if                                                        4d25s22
          if(lds)then                                                   2d27s23
           do iz=0,nneedzm                                              7d1s16
            bc(jb+iz)=bc(jb+iz)-(bc(jd1+iz)+bc(jd2+iz)+bc(jd3+iz))      7d1s16
           end do
          end if                                                        4d25s22
          if(itype.eq.iapair)then                                       3d1s23
           do iz=0,nneedzm                                              3d1s23
            bc(jb+iz)=bc(jb+iz)*design(idersign)                        3d1s23
           end do                                                       3d1s23
          end if                                                        3d1s23
          call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $          jb1,0,0,0,0,0,0,ixyz,0,2,0,1d0,ldeb,xcart(1,ia),bc,ibc) 11d10s22
          call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $          jb2,0,0,0,0,0,0,ixyz,0,0,2,1d0,ldeb,xcart(1,ia),bc,ibc) 11d10s22
          call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $          jb3,0,0,0,0,0,0,ixyz,1,1,0,2d0,ldeb,xcart(1,ia),bc,ibc) 11d10s22
          call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $          jb4,0,0,0,0,0,0,ixyz,1,0,1,2d0,ldeb,xcart(1,ia),bc,ibc) 11d10s22
          call nattracd(lb,zetb,xnb,xb,yb,zb,lk,zetk,xnk,xk,yk,zk,      2d4s16
     $         atnum(1,ia),xcart(1,ia),xcart(2,ia),xcart(3,ia),         2d4s16
     $          jb5,0,0,0,0,0,0,ixyz,0,1,1,2d0,ldeb,xcart(1,ia),bc,ibc) 11d10s22
          do iz=0,nneedzm                                               6d30s16
           bc(jb+iz)=bc(jb+iz)+bc(jb1+iz)+bc(jb2+iz)+bc(jb3+iz)         6d30s16
     $          +bc(jb4+iz)+bc(jb5+iz)                                  6d30s16
          end do                                                        6d30s16
          if(ldeb)write(6,*)('take dersign to be '),dersign,
     $         idersign
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
         zm=-atnum(1,iau)                                               8d20s15
         idxp=iderx+idxyz(1)
         idyp=idery+idxyz(2)
         idzp=iderz+idxyz(3)
         if(iab.eq.iak.and.(iab.eq.iatom.or.iab.eq.iapair))then         3d6s23
          arg=+zm                                                       5d6s16
          call derid(lb,zetb,xnb,xb,yb,zb,0,lk,zetk,xnk,xk,yk,zk,0,     2d4s16
     $         0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),
     $         xcart(3,iau),0,0,atnum(2,iau),atnum(3,iau),xcart(1,iau),
     $         xcart(2,iau),xcart(3,iau),0,dum,1,idum,jb,
     $         idxp,idyp,idzp,iderx,idery,iderz,
     $         0,0,0,0,0,0,1,ldeb,arg,bc,ibc)                           11d14s22
          if(ldeb)then
           rsum=bc(jb)
           write(6,*)('after derid the first '),bc(jb)
           call prntm2(bc(jb),1,nneedz,1)
          end if
          idxh=iderx+(idxyz(1)/2)                                       7d1s16
          idyh=idery+(idxyz(2)/2)                                       7d1s16
          idzh=iderz+(idxyz(3)/2)                                       7d1s16
          arg=zm*2d0                                                    7d1s16
          call derid(lb,zetb,xnb,xb,yb,zb,0,lk,zetk,xnk,xk,yk,zk,0,     2d4s16
     $         0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),
     $         xcart(3,iau),0,0,atnum(2,iau),atnum(3,iau),xcart(1,iau),
     $         xcart(2,iau),xcart(3,iau),0,dum,1,idum,jb2,
     $         idxh,idyh,idzh,idxh,idyh,idzh,
     $         0,0,0,0,0,0,1,ldeb,arg,bc,ibc)                           11d14s22
          if(ldeb)then
           rsum=rsum+bc(jb2)
           write(6,*)('after derid the firstb '),rsum,bc(jb2)
           call prntm2(bc(jb2),1,nneedz,1)
          end if
          arg=+zm                                                       5d6s16
          call derid(lb,zetb,xnb,xb,yb,zb,0,lk,zetk,xnk,xk,yk,zk,0,     2d4s16
     $          0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),
     $          xcart(3,iau),0,0,atnum(2,iau),atnum(3,iau),xcart(1,iau),
     $          xcart(2,iau),xcart(3,iau),0,dum,1,idum,ib1d,
     $          iderx,idery,iderz,idxp,idyp,idzp,
     $          0,0,0,0,0,0,1,ldeb,arg,bc,ibc)                          11d14s22
          if(ldeb)then
           rsum=rsum+bc(ib1d)
           write(6,*)('after derid the second '),rsum
           call prntm2(bc(ib1d),1,nneedz,1)
          end if
          do iz=0,nneedzm                                               2d3s16
           bc(jb+iz)=bc(jb+iz)+bc(ib1d+iz)+bc(jb2+iz)                   7d1s16
          end do
          if(ldeb)write(6,*)('whats under jb: '),bc(jb)
          ibcoff=ib1d                                                   3d10s16
         else if(iab.eq.iatom)then
          arg=+zm                                                       5d6s16
          call derid(lb,zetb,xnb,xb,yb,zb,0,lk,zetk,xnk,xk,yk,zk,0,     2d4s16
     $          0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),
     $          xcart(3,iau),0,0,atnum(2,iau),atnum(3,iau),xcart(1,iau),
     $          xcart(2,iau),xcart(3,iau),0,dum,1,idum,jb,
     $          idxp,idyp,idzp,iderx,idery,iderz,                       2d4s16
     $         0,0,0,0,0,0,1,ldeb,arg,bc,ibc)                           11d14s22
          if(ldeb)then
           write(6,*)('after derid the third ')
           call prntm2(bc(jb),1,nneedz,1)
          end if
         else if(iak.eq.iatom)then                                      3d6s23
          arg=+zm                                                       5d6s16
          call derid(lb,zetb,xnb,xb,yb,zb,0,lk,zetk,xnk,xk,yk,zk,0,     2d4s16
     $          0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),
     $          xcart(3,iau),0,0,atnum(2,iau),atnum(3,iau),xcart(1,iau),
     $          xcart(2,iau),xcart(3,iau),0,dum,1,idum,jb,
     $          iderx,idery,iderz,idxp,idyp,idzp,
     $         0,0,0,0,0,0,1,ldeb,arg,bc,ibc)                           11d14s22
          if(ldeb)then
           write(6,*)('after derid the fourth '),bc(jb)
           call prntm2(bc(jb),1,nneedz,1)
          end if
         end if
         arg=+zm                                                        3d6s23
         if(iab.ne.iak)then
          if(iak.eq.iapair)then
           call derid(lb,zetb,xnb,xb,yb,zb,0,lk,zetk,xnk,xk,yk,zk,0,     2d4s16
     $          0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),
     $          xcart(3,iau),0,0,atnum(2,iau),atnum(3,iau),xcart(1,iau),
     $          xcart(2,iau),xcart(3,iau),0,dum,1,idum,jbb,
     $          iderx,idery,iderz,idxp,idyp,idzp,
     $         0,0,0,0,0,0,1,ldeb,arg,bc,ibc)                           11d14s22
           if(jb.gt.0)then                                               3d6s23
            do iz=0,nneedzm                                              3d6s23
             bc(jb+iz)=bc(jb+iz)+bc(jbb+iz)                              3d6s23
            end do                                                       3d6s23
           else                                                          3d6s23
            jb=jbb                                                       3d6s23
           end if                                                        3d6s23
          else if(iab.eq.iapair)then                                    3d6s23
           call derid(lb,zetb,xnb,xb,yb,zb,0,lk,zetk,xnk,xk,yk,zk,0,     2d4s16
     $          0,atnum(2,iau),atnum(3,iau),xcart(1,iau),xcart(2,iau),
     $          xcart(3,iau),0,0,atnum(2,iau),atnum(3,iau),xcart(1,iau),
     $          xcart(2,iau),xcart(3,iau),0,dum,1,idum,jbb,             3d6s23
     $          idxp,idyp,idzp,iderx,idery,iderz,                       2d4s16
     $         0,0,0,0,0,0,1,ldeb,arg,bc,ibc)                           11d14s22
           if(jb.gt.0)then                                               3d6s23
            do iz=0,nneedzm                                              3d6s23
             bc(jb+iz)=bc(jb+iz)+bc(jbb+iz)                              3d6s23
            end do                                                       3d6s23
           else                                                          3d6s23
            jb=jbb                                                       3d6s23
           end if                                                        3d6s23
          end if                                                        3d6s23
         end if                                                         3d6s23
         if(jau.eq.iatom)then                                           2d3s16
          ib1sav=-1
          if(jb.gt.0)then
           ib1sav=jb
          end if
          arg=zm                                                        3d6s23
          xuse1=0d0                                                     7d14s15
          xuse2=atnum(2,iau)*2d0                                        7d14s15
          if(ldeb)write(6,*)('before derid, '),bc(ib1sav),ib1sav,jb
          call derid(lb,zetb,xnb,xb,yb,zb,0,lk,zetk,xnk,xk,yk,zk,0,     2d4s16
     $          0,xuse2,atnum(3,iau),xcart(1,iau),xcart(2,iau),
     $          xcart(3,iau),0,0,xuse1,atnum(3,iau),xcart(1,iau),
     $          xcart(2,iau),xcart(3,iau),0,dum,1,idum,jb,
     $          iderx,idery,iderz,iderx,idery,iderz,                    2d4s16
     $          idxyz(1),idxyz(2),idxyz(3),0,0,0,1,ldeb,arg,bc,ibc)     11d14s22
          if(ldeb)then
           write(6,*)('after derid the fifth '),bc(jb),bc(ib1sav),jb
           rsum2=bc(jb)
           call prntm2(bc(jb),1,nneedz,1)
          end if
          if((iab.eq.iatom.and.iau.eq.jau).or.                          3d6s23
     $       (iab.eq.iapair.and.iau.ne.jau))then                        3d6s23
           idxh=idxyz(1)/2                                               6d25s16
           idyh=idxyz(2)/2                                               6d25s16
           idzh=idxyz(3)/2                                               6d25s16
           iderxp=iderx+idxh
           ideryp=idery+idyh
           iderzp=iderz+idzh
           arg=2d0*zm
           twoz=2d0*atnum(2,iau)                                        7d14s15
           call derid(lb,zetb,xnb,xb,yb,zb,0,lk,zetk,xnk,xk,yk,zk,0,     6d25s16
     $          0,0d0,atnum(3,iau),xcart(1,iau),xcart(2,iau),           7d14s15
     $          xcart(3,iau),0,0,twoz,atnum(3,iau),xcart(1,iau),        7d14s15
     $          xcart(2,iau),xcart(3,iau),0,dum,1,idum,lbx,              6d25s16
     $          iderxp,ideryp,iderzp,iderx,idery,iderz,                    6d25s16
     $          0,0,0,idxh,idyh,idzh,1,ldeb,arg,bc,ibc)                 11d14s22
           if(ldeb)then
            rsum3=bc(lbx)
            write(6,*)('after derid the seventh ')
            call prntm2(bc(lbx),1,nneedz,1)
           end if
           do iz=0,nneedzm                                              7d1s16
            bc(jb+iz)=bc(jb+iz)+bc(lbx+iz)                              7d14s15
           end do                                                       7d1s16
          end if
          if((iak.eq.iatom.and.iau.eq.jau).or.                          3d6s23
     $       (iak.eq.iapair.and.iau.ne.jau))then                        3d6s23
           idxh=idxyz(1)/2                                               6d25s16
           idyh=idxyz(2)/2                                               6d25s16
           idzh=idxyz(3)/2                                               6d25s16
           iderxp=iderx+idxh
           ideryp=idery+idyh
           iderzp=iderz+idzh
           arg=2d0*zm
           twoz=atnum(2,iau)*2d0                                        7d14s15
           if(ldeb)write(6,*)('before derid: '),bc(jb),bc(ib1sav),jb,
     $          ib1sav
           call derid(lb,zetb,xnb,xb,yb,zb,0,lk,zetk,xnk,xk,yk,zk,0,     6d25s16
     $          0,0d0,atnum(3,iau),xcart(1,iau),xcart(2,iau),           7d14s15
     $          xcart(3,iau),0,0,twoz,atnum(3,iau),xcart(1,iau),        7d14s15
     $          xcart(2,iau),xcart(3,iau),0,dum,1,idum,lbx,              6d25s16
     $          iderx,idery,iderz,iderxp,ideryp,iderzp,                    6d25s16
     $          0,0,0,idxh,idyh,idzh,1,ldeb,arg,bc,ibc)                 11d14s22
           if(ldeb)then
            rsum30=rsum3
            xlbx0=bc(lbx)
            rsum3=rsum3+bc(lbx)
            write(6,*)('after derid the eighth '),rsum3,rsum+rsum2+rsum3
            call prntm2(bc(lbx),1,nneedz,1)
           end if
           do iz=0,nneedzm                                              7d1s16
            bc(jb+iz)=bc(jb+iz)+bc(lbx+iz)                                 7d1s16
           end do                                                       7d1s16
          end if
          if(ldeb)write(6,*)('ib1sav '),ib1sav,bc(jb),nneedzm
          if(ib1sav.gt.0)then
           do iz=0,nneedzm                                              2d3s16
            bc(ib1sav+iz)=bc(ib1sav+iz)+bc(jb+iz)
           end do
           ibcoff=jb                                                    3d10s16
           jb=ib1sav
          end if                                                        6d25s16
         end if
        end if
       else                                                             6d5s14
        jb=ibcoff                                                        6d5s14
        ibcoff=jb+nneedz                                                6d5s14
        call enough('oned2ints.  3',bc,ibc)
        do ineedz=0,nneedzm                                             2d3s16
         bc(jb+ineedz)=0d0                                              6d5s14
        end do                                                          6d5s14
       end if
      end if                                                            5d7s10
      return
      end
