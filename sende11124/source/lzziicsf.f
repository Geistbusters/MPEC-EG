c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine lzziicsf(hdiag,nfcn,idorb,isorb,mdon,mdoo,ibasis,iptr, 12d24s19
     $     ixlzz,shift,nec,islz,multh,npass,bc,ibc)                     11d10s22
c
c     compute spin independent part of lzlzii or L^2ii                  5d14s21
c     if atom, i.e. npass=3, compute  lzlzii as well.                   5d14s21
c
      implicit real*8 (a-h,o-z)
      integer*1 idorb(*),isorb(*)
      logical ldebug                                                    7d1s19
      dimension hdiag(nfcn,*),iptr(4,*),ixlzz(8,*),ibasis(*),multh(8,8),5d14s21
     $     islz(*)                                                      1d2s20
      include "common.store"
      include "common.mrci"                                             6d10s19
      ione=npass+1                                                      12d31s19
      ionep=ione+1                                                      5d14s21
      ldebug=.false.                                                    7d15s19
      if(ldebug)write(6,*)('hi, i''m lzziicsf'),nfcn,(islz(i),i=1,npass)1d2s20
     $     ,loc(ixlzz)                                                  1d2s20
      if(ldebug)then
       do isb=1,8
        if(irefo(isb).gt.0)then                                         1d2s20
         write(6,*)('ixlzz addresses: '),isb,(ixlzz(isb,j),j=1,npass+1) 1d2s20
         write(6,*)('one electron part: ')
         call prntm2(bc(ixlzz(isb,ione)),irefo(isb),irefo(isb),
     $       irefo(isb))
         if(npass.ne.1)then                                             5d14s21
          write(6,*)('Lz^2 matrix ')                                    5d14s21
          call prntm2(bc(ixlzz(isb,ionep)),irefo(isb),irefo(isb),       5d14s21
     $       irefo(isb))
         end if                                                         5d14s21
         do ipass=1,npass                                                12d31s19
          jsb=multh(islz(ipass),isb)
          if(irefo(jsb).gt.0)then                                       12d31s19
           write(6,*)('for two electron matrix '),ipass,jsb
           call prntm2(bc(ixlzz(isb,ipass)),irefo(isb),irefo(jsb),
     $          irefo(isb))
          end if
         end do                                                          12d31s19
        end if
       end do
       write(6,*)('iptr: ')
       do j=1,4
        write(6,*)j,(iptr(j,k),k=1,3)
       end do
      end if                                                            7d15s19
      if(npass.eq.1)then                                                5d14s21
       ncomp=1                                                          5d14s21
      else                                                              5d14s21
       ncomp=2                                                          5d14s21
      end if                                                            5d14s21
      nclo=-1
      ihdiag=1
      do if=1,nfcn
       jbasis=3*(if-1)
       if(ldebug)                                                       7d1s19
     $    write(6,*)('fcn no. '),if,('basis: '),(ibasis(jbasis+j),j=1,3)7d1s19
       if(ibasis(jbasis+1).ne.nclo)then
        if(nclo.ge.0)then
         ibcoff=ihc
        end if
        nclo=ibasis(jbasis+1)
        nclop=nclo+1
        if(ldebug)write(6,*)('rule sizes: '),iptr(1,nclop),iptr(3,nclop)7d1s19
        ihc=ibcoff
        iho=ihc+iptr(1,nclop)*ncomp                                     5d14s21
        ibcoff=iho+iptr(3,nclop)*ncomp                                  5d14s21
        call enough('lzziicsf.  1',bc,ibc)
        do i=ihc,ibcoff-1                                               5d14s21
         bc(i)=0d0                                                      5d14s21
        end do                                                          5d14s21
        khc=ihc+iptr(1,nclop)                                           5d14s21
        kho=iho+iptr(3,nclop)                                           5d14s21
        if(ldebug)write(6,*)('for closed shells: '),iptr(2,nclop)       7d1s19
        jc=iptr(2,nclop)
        do i=0,iptr(1,nclop)-1
         ip=i+1
         if(ldebug)write(6,*)ip,(idorb(jc+j),j=0,nclo-1)                7d1s19
         do j=0,nclo-1
          iorb=idorb(jc+j)
          isb=ism(iorb)
          ig=irel(iorb)-1+idoubo(isb)                                   12d24s19
          ndm=irefo(isb)+idoubo(isb)                                    12d24s19
          ih=ixlzz(isb,ione)+ig*(ndm+1)                                 1d2s20
          bc(ihc+i)=bc(ihc+i)+2d0*bc(ih)
          if(npass.ne.1)then                                            5d14s21
           kh=ixlzz(isb,ionep)+ig*(ndm+1)                               5d14s21
           bc(khc+i)=bc(khc+i)+2d0*bc(kh)                               5d14s21
          end if                                                        5d14s21
          do l=0,nclo-1
           jorb=idorb(jc+l)
           jsb=ism(jorb)
           jg=irel(jorb)-1+idoubo(jsb)                                  12d24s19
           do ipass=1,npass                                             12d31s19
            if(multh(isb,jsb).eq.islz(ipass))then                       1d2s20
             ndm=idoubo(isb)+irefo(isb)                                  12d24s19
             iad=ixlzz(isb,ipass)+ig+ndm*jg                             1d2s20
             bc(ihc+i)=bc(ihc+i)-2d0*bc(iad)*bc(iad)                     12d24s19
             if(ipass.eq.3)bc(khc+i)=bc(khc+i)-2d0*bc(iad)*bc(iad)      5d14s21
            end if                                                       12d24s19
           end do                                                       1d2s20
          end do
         end do
         jc=jc+nclo
         if(ldebug)write(6,*)bc(ihc+i)                                  7d1s19
        end do
        nopen=nec-2*nclo                                                6d10s19
        if(ldebug)write(6,*)('for open shells: '),iptr(4,nclop)         7d1s19
        jo=iptr(4,nclop)
        do i=0,iptr(3,nclop)-1
         ip=i+1
         if(ldebug)write(6,*)ip,(isorb(jo+j),j=0,nopen-1)               7d1s19
         do j=0,nopen-1
          iorb=isorb(jo+j)
          isb=ism(iorb)
          ig=irel(iorb)-1+idoubo(isb)                                   12d24s19
          ndm=irefo(isb)+idoubo(isb)                                    12d24s19
          ih=ixlzz(isb,ione)+ig*(ndm+1)                                 1d2s20
          bc(iho+i)=bc(iho+i)+bc(ih)
          if(npass.ne.1)then                                            5d14s21
           kh=ixlzz(isb,ionep)+ig*(ndm+1)                               5d14s21
           bc(kho+i)=bc(kho+i)+bc(kh)                                   5d14s21
          end if                                                        5d14s21
         end do
         jo=jo+nopen
         if(ldebug)write(6,*)bc(iho+i)                                  7d1s19
        end do
        jho=iho-1
        jhc=ihc-1
        jkho=kho-1                                                      5d14s21
        jkhc=khc-1                                                      5d14s21
       end if
       hd=bc(jhc+ibasis(jbasis+2))+bc(jho+ibasis(jbasis+3))
       hdzz=0d0                                                         5d14s21
       if(npass.ne.1)hdzz=bc(jkhc+ibasis(jbasis+2))                     5d14s21
     $      +bc(jkho+ibasis(jbasis+3))                                  5d14s21
       if(ldebug)write(6,*)('sum of oo and cc parts: '),hd,hdzz         5d14s21
       jc=iptr(2,nclop)+nclo*(ibasis(jbasis+2)-1)-1                     6d14s19
       jo=iptr(4,nclop)+nopen*(ibasis(jbasis+3)-1)-1                    6d14s19
       if(ldebug)then                                                   7d1s19
       write(6,*)('mixing between '),(idorb(jc+j),j=1,nclo)
       write(6,*)('and '),(isorb(jo+j),j=1,nopen)
       end if                                                           7d1s19
       do i=1,nclo
        iorb=idorb(jc+i)
        isb=ism(iorb)
        ig=irel(iorb)-1+idoubo(isb)                                     12d24s19
        do j=1,nopen                                                    6d11s19
         jorb=isorb(jo+j)                                               6d11s19
         jsb=ism(jorb)                                                  6d11s19
         jg=irel(jorb)-1+idoubo(jsb)                                    12d24s19
         do ipass=1,npass                                               12d31s19
          if(multh(isb,jsb).eq.islz(ipass))then                         1d2s20
           ndm=idoubo(isb)+irefo(isb)                                    12d24s19
           iad=ixlzz(isb,ipass)+ig+ndm*jg                               1d2s20
           hd=hd-2d0*bc(iad)*bc(iad)                                     12d24s19
           if(ipass.eq.3)hdzz=hdzz-2d0*bc(iad)*bc(iad)                  5d14s21
          end if                                                         12d24s19
         end do                                                         1d2s20
        end do                                                          6d11s19
       end do
       if(ldebug)write(6,*)('final spin independent part: '),hd,hdzz    5d14s21
       hdiag(ihdiag,1)=hd                                               5d14s21
       if(npass.ne.1)hdiag(ihdiag,2)=hdzz                               5d14s21
       ihdiag=ihdiag+1
      end do
      return
      end
