c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine makeguess(idorel,idorelg,nsymb,nsymbg,ngaus,ngausg,    2d18s16
     $     bdat,bdatg,nbasdws,nbasdwsg,vecr,vgues,ovr,myguessi,ierr,    3d21s16
     $     idwsdeb,ibstor,isstor,ibstorg,isstorg,ibcode,iptno,ipts,     10d24s17
     $     iapair,isym,iapairg,isymg,nhsz,ipao,morbp,nbasisp,iocode,    3d25s20
     $     idorel4c,ascale,iblstor,ibsstor,nbaslarge,nbassmall,ilsz,    3d27s20
     $     isnorm,nbasisc,nlzz,iorbsym,iorbsymz,morbc,morb,scopy,       4d14s22
     $     nrydb,nvguess,ircode,bc,ibc,ivoffg,icanog)                   5d3s23
      implicit real*8 (a-h,o-z)                                         5d4s23
c
c     generate guess in ao basis from guess orbitals
c     iocode: 0 called by HF and overlap is already squared             2d15s19
c             1 called by cal1int and overlap had not yet been squared  2d15s19
c     if nrydb ne 0, geometry must be the same!                         4d14s22
c
      logical ldebug                                                    5d3s21
      include "common.store"
      include "common.print"                                            5d3s21
      integer*8 iarg8,lw,lg,ipw,ipg,ibstor(*),isstor(*),ibstorg(*),     10d23s17
     $     isstorg(*),iblstor(*),ibsstor(*),lwp                         3d27s20
      dimension bdat(ngaus,9),bdatg(ngausg,9),nbasdws(8),nbasdwsg(8),   4d27s21
     $     vecr(1),vgues(1),ovr(1),iovs(8,8),ibcode(8),iptno(8),ipts(8),10d24s17
     $     iapair(3,*),cartk(3),cartb(3),isym(3,*),iapairg(3,*),        10d27s17
     $     isymg(3,*),nhsz(*),ipao(*),igot(3),morbp(*),nbasisp(*),      5d14s19
     $     shift(3),ider(3,3),nbasisc(*),mrow(8),iorbsym(*),iorbsymz(*),5d5s21
     $     ivoffg(8),ivoffw(8),morbc(*),morb(*),scopy(*),icanog(*)      5d3s23
      equivalence (arg8,iarg8)                                          10d23s17
      data ider/1,0,0, 0,1,0, 0,0,1/                                    3d27s20
      if(iprtr(23).eq.0)then                                            5d3s21
       ldebug=.false.                                                   5d3s21
      else                                                              5d3s21
       ldebug=.true.                                                    5d3s21
      end if                                                            5d3s21
      jpao=0                                                            1d17s23
      igoal=1749477+10
      ircodecx=3                                                        1d4s23
      nvguess=1                                                         10d11s22
      ircode=1                                                          10d12s22
      if(ldebug)write(6,*)('ircode a')                                  1d24s23
      ncomp=1                                                           4d27s21
      if(idorel.ne.0)ncomp=2                                            4d27s21
      srh=sqrt(0.5d0)                                                   10d26s17
      if(ldebug)                                                        5d3s21
     $     write(6,*)('at top of makeguess, ibcoff = '),ibcoff,idwsdeb
      iloopit=0
      ierr=0                                                            2d18s16
      ibdiff=0
      ibcoffo=ibcoff
      if(ldebug)write(6,*)('ngaus vs. ngausg: '),ngaus,ngausg           5d3s21
      if(ngaus.ne.ngausg)then
       if(ldebug)write(6,*)('ngaus ne ngausg block ')                   5d3s21
       ibdiff=1
      else
       ib=ibcoff
       ibg=ib+9                                                         4d27s21
       ibcoff=ibg+9                                                     4d27s21
       jb=ib-1
       jbg=ibg-1
       call enough('makeguess.  1',bc,ibc)
       do i=1,ngaus
        do j=1,9                                                        4d27s21
         bc(jb+j)=bdat(i,j)
         bc(jbg+j)=bdatg(i,j)
        end do
        if(ibc(ib).ne.ibc(ibg))then
         ibdiff=1
         go to 1
        end if                                                          2d18s16
        if(bc(ib+1).ne.bc(ibg+1))then
         ibdiff=1
         go to 1
        end if
        if(ibc(ib+7).ne.ibc(ibg+7))then
         ibdiff=1
         go to 1
        end if
        if(ibdiff.ne.0)stop 'ibdiff'
       end do
      end if
    1 continue
      rms=0d0                                                           1d19s23
      irms=0                                                            1d19s23
      if(ldebug)write(6,*)('iocode,ibdiff,idorel,idorelg '),iocode,     5d3s21
     $     ibdiff,idorel,idorelg                                        5d3s21
      if(iocode.eq.1.or.nrydb.ne.0)then                                 4d14s22
       if(ldebug)write(6,*)('do we have different geometries? ')        5d3s21
       nuw=0                                                            10d12s22
       igeow=ibcoff                                                     10d12s22
       do ig=1,ngaus                                                    10d12s22
        do i=1,nuw                                                      10d12s22
         jgeow=igeow-1+3*(i-1)                                          10d12s22
         dd=0d0                                                         10d12s22
         do ixyz=1,3                                                    10d12s22
          diff=bdat(ig,4+ixyz)-bc(jgeow+ixyz)                           10d12s22
          dd=dd+diff**2                                                 10d12s22
         end do                                                         10d12s22
         dd=sqrt(dd/3d0)                                                10d12s22
         if(dd.lt.1d-8)go to 85                                         10d12s22
        end do                                                          10d12s22
        jgeow=igeow-1+3*nuw                                             10d12s22
        do ixyz=1,3                                                     10d12s22
         bc(jgeow+ixyz)=bdat(ig,4+ixyz)                                 10d12s22
        end do                                                          10d12s22
        nuw=nuw+1                                                       10d12s22
   85   continue                                                        10d12s22
       end do                                                           10d12s22
       if(ldebug)then                                                   1d24s23
        write(6,*)('want atoms: ')
        call prntm2(bc(igeow),3,nuw,3)                                   10d12s22
       end if                                                           1d24s23
       ibcoff=igeow+3*nuw                                               10d12s22
       do igg=1,ngausg                                                  10d12s22
        do i=1,nuw                                                      10d12s22
         jgeow=igeow-1+3*(i-1)                                          10d12s22
         dd=0d0                                                         10d12s22
         do ixyz=1,3                                                    10d12s22
          diff=bdatg(igg,4+ixyz)-bc(jgeow+ixyz)                         10d12s22
          dd=dd+diff**2                                                 10d12s22
         end do                                                         10d12s22
         dd=sqrt(dd/3d0)                                                10d12s22
         if(dd.lt.1d-8)go to 78                                         10d12s22
        end do                                                          10d12s22
        if(ldebug)                                                      1d24s23
     $     write(6,*)('no match for '),(bdatg(igg,4+ixyz),ixyz=1,3),igg 1d24s23
        rms=rms+dd**2                                                   10d12s22
        irms=irms+1                                                     10d12s22
   78   continue                                                        10d12s22
       end do                                                           10d12s22
       ibcoff=igeow                                                     10d12s22
       rms=sqrt(rms/dfloat(max(1,irms)))                                10d12s22
       if(rms.gt.1d-14.and.(iocode.eq.1.or.nrydb.ne.0))then             6d3s22
        write(6,*)('rms difference in geometries = '),rms                10d9s19
        if(iocode.eq.1)then                                             4d14s22
         write(6,*)('you are trying to run ci calculation off of '),
     $      ('orbitals from a different geometry')
         write(6,*)('stopping')                                          10d9s19
        else                                                            4d14s22
         write(6,*)('you are trying to generate Rydberg orbitals off '),4d14s22
     $      ('orbitals from a different geometry!!!')                   4d14s22
        end if
        stop 'makeguess'                                                10d9s19
       end if                                                           10d9s19
      end if                                                            10d9s19
      if(ldebug)then                                                    1d24s23
       write(6,*)('nsymb,nsymbg '),nsymb,nsymbg
       write(6,*)('ibdiff '),ibdiff
       write(6,*)('rms '),rms
      end if                                                            1d24s23
      if(nsymb.ne.nsymbg.or.ibdiff.ne.0.or.rms.gt.1d-14)then            6d3s22
       if(idorel4c.ne.0)then                                            3d27s20
        write(6,*)('we are in 4 component block ')                      3d27s20
        pi=acos(-1d0)                                                   3d27s20
        oneover4pi=0.25d0/pi                                            3d27s20
        if(idorelg.eq.0)then                                            3d27s20
         write(6,*)('using nr orbitals to generate 4c spinors ')        3d27s20
        else                                                            3d27s20
         write(6,*)('using 2c spinors to generate 4c spinors ')         3d27s20
        end if                                                          3d27s20
        nxxt=2*(nbaslarge+nbassmall)
        nxx=2*nbaslarge                                                 3d27s20
        nxx2=nxxt*2*nbaslarge                                           3d27s20
c
c     In 4c space, the overlap matrix is block diagonal in alpha/beta   3d27s20
c     and large/small. Thus the large alpha/beta parts of a spinor      3d27s20
c     are that of the input guess vector.                               3d27s20
c     (but we still need overlaps to figure out symmetry differences)   3d27s20
c     In the nr or 2c cases, the small component is generated via       3d27s20
c     for alpha large, alpha small=(-i/2c)d/dz,                         3d27s20
c     beta small=(1/2c)[d/dy-id/dx] and                                 3d27s20
c     for beta large, alpha small=-(1/2c)[d/dy+id/dx],                  3d27s20
c     beta small=(i/2c)d/dz.                                            3d27s20
c     so we will compute matrix elements with the bra being the small   3d27s20
c     basis fcns and the ket being the ders on the large component      3d27s20
c     basis functions, multiply times the input vectors, then multiply  3d27s20
c     by the inverse of the small component overlap matrix, and Bob's   3d27s20
c     your uncle.
c
        ngt=0                                                           3d27s20
        num1c=0                                                         3d27s20
        num2c=0                                                         3d27s20
        num3c=0                                                         3d27s20
        do isb=1,nsymbg                                                 3d27s20
         iovs(1,isb)=ibcoff                                             3d27s20
         ngt=ngt+morbp(isb)                                             3d27s20
         iovs(2,isb)=iovs(1,isb)+nbaslarge*morbp(isb)                   3d27s20
         iovs(3,isb)=iovs(2,isb)+nbassmall*morbp(isb)                   3d27s20
         iovs(4,isb)=iovs(3,isb)+nbassmall*morbp(isb)                   3d27s20
         ibcoff=iovs(4,isb)+nbassmall*morbp(isb)                        3d27s20
         do i=0,morbp(isb)-1                                            3d27s20
          if(ibc(ibcode(isb)+i).eq.1)then                               3d27s20
           num1c=num1c+1                                                3d27s20
          else if(ibc(ibcode(isb)+i).eq.2)then                          3d27s20
           num2c=num2c+1                                                3d27s20
          else if(ibc(ibcode(isb)+i).eq.3)then                          3d27s20
           num3c=num3c+1                                                3d27s20
          end if                                                        3d27s20
         end do                                                         3d27s20
        end do                                                          3d27s20
        ioff1c=1                                                        3d27s20
        ioff2c=ioff1c+num1c*2*nxxt                                      3d27s20
        ioff3c=ioff2c+num2c*2*nxxt                                      3d27s20
        ioff0c=ioff3c+num3c*2*nxxt                                      3d27s20
        call enough('makeguess.  2',bc,ibc)
        do i=iovs(1,1),ibcoff-1                                         3d27s20
         bc(i)=0d0                                                      3d27s20
        end do                                                          3d27s20
        do iwg=1,ngaus                                                   10d23s17
         arg8=bdat(iwg,1)                                                10d23s17
         lw=iarg8                                                        10d23s17
         lwp=lw+1                                                       3d27s20
         nss=ibc(ilsz+lwp)                                              3d27s20
         arg8=bdat(iwg,4)                                                10d23s17
         ipw=iarg8                                                       10d23s17
         nw=2*lw+1                                                       10d23s17
         arg8=bdat(iwg,8)                                                10d26s17
         japair=iarg8                                                    10d26s17
         cartk(1)=bdat(iwg,5)                                           3d27s20
         cartk(2)=bdat(iwg,6)                                           3d27s20
         cartk(3)=bdat(iwg,7)                                           3d27s20
         iarg8=bdat(iwg,9)                                              4d27s21
         ibdat9=iarg8
         do igg=1,ngausg                                                 10d23s17
          arg8=bdatg(igg,1)                                              10d23s17
          lg=iarg8                                                       10d23s17
          arg8=bdatg(igg,4)                                              10d23s17
          ipg=iarg8                                                      10d23s17
          arg8=bdatg(igg,8)                                              10d27s17
          kapair=iarg8                                                   10d27s17
          ng=2*lg+1                                                      10d23s17
          cartb(1)=bdatg(igg,5)                                         3d27s20
          cartb(2)=bdatg(igg,6)                                         3d27s20
          cartb(3)=bdatg(igg,7)                                         3d27s20
          call onep(lw,bdat(iwg,2),bdat(iwg,3),bdat(iwg,5),bdat(iwg,6), 3d27s20
     $         bdat(iwg,7),lg,bdatg(igg,2),bdatg(igg,3),bdatg(igg,5),   3d27s20
     $         bdatg(igg,6),bdatg(igg,7),idum,ib1,0,0,0, 0,0,0, 0,0,0,  11d9s22
     $         bc,ibc)                                                  11d9s22
c
c     recall we are not making cart to sphere transformation.
c     thus we need to normalize and the indices come out swapped.
c     note this will cause a problem if higher than p functions are
c     used!!!
c
          if(lg.gt.1)then                                               3d27s20
           write(6,*)('oh no! lg = '),lg,(' but I don''t have '),       3d27s20
     $         ('cart to sphere transformation!')                       3d27s20
           stop                                                         3d27s20
          end if                                                        3d27s20
          ilw=2*lw+1                                                        3d27s20
          ilg=2*lg+1                                                    3d27s20
          xnorm=oneover4pi*sqrt(dfloat(ilw*ilg))                        3d27s20
          do ig=1,ng                                                    3d27s20
           isg=isstorg(ipg+ig)                                          3d27s20
           ibb=ibstorg(ipg+ig)-1                                        3d27s20
           do iw=0,nw-1                                                 3d27s20
            ikk=iblstor(iwg)+iw                                         3d27s20
            ito=iovs(1,isg)+ikk+nbaslarge*ibb                            3d27s20
            bc(ito)=bc(ib1)*xnorm                                       3d27s20
            ib1=ib1+1                                                   3d27s20
           end do                                                       3d27s20
          end do                                                        3d27s20
          isnuse=isnorm+iwg-1                                           3d27s20
          ilw=2*lw+3                                                        3d27s20
          xnorm=oneover4pi*sqrt(dfloat(ilw*ilg))                        3d27s20
          do ixyz=1,3                                                   3d27s20
           call onep(lwp,bdat(iwg,2),bc(isnuse),bdat(iwg,5),bdat(iwg,6), 3d27s20
     $         bdat(iwg,7),lg,bdatg(igg,2),bdatg(igg,3),bdatg(igg,5),   3d27s20
     $         bdatg(igg,6),bdatg(igg,7),idum,ib1,0,0,0, 0,0,0,         3d27s20
     $         ider(1,ixyz),ider(2,ixyz),ider(3,ixyz),bc,ibc)           11d9s22
           do ig=1,ng                                                    3d27s20
            isg=isstorg(ipg+ig)                                          3d27s20
            ibb=ibstorg(ipg+ig)-1                                        3d27s20
            do iw=0,nss-1                                               3d27s20
             ikk=ibsstor(iwg)+iw                                        3d27s20
             ito=iovs(1+ixyz,isg)+ikk+nbassmall*ibb                             3d27s20
             bc(ito)=bc(ib1)*xnorm                                       3d27s20
             ib1=ib1+1                                                   3d27s20
            end do                                                      3d27s20
           end do                                                       3d27s20
          end do                                                        3d27s20
         end do                                                         3d27s20
        end do                                                          3d27s20
        ivoff=1                                                         3d27s20
        factr=sqrt(ascale)                                              3d27s20
        ioltmp=ibcoff                                                    3d27s20
        ipvtl=ioltmp+nbaslarge*nbaslarge                                 3d27s20
        ibcoff=ipvtl+nbaslarge                                          3d27s20
        call enough('makeguess.  3',bc,ibc)
        do i=0,nbaslarge-1                                              3d27s20
         iad1=1+nxxt*i                                                  3d27s20
         iad2=ioltmp+nbaslarge*i                                        3d27s20
         do j=0,nbaslarge-1                                             3d27s20
          bc(iad2+j)=ovr(iad1+j)                                        3d27s20
         end do                                                         3d27s20
        end do                                                          3d27s20
        call lusolv(bc(ioltmp),nbaslarge,nbaslarge,dum,1,1,ibc(ipvtl),  3d27s20
     $       ierr,1)                                                    3d27s20
        if(ierr.ne.0)stop                                               3d27s20
        ioff=1+nbaslarge*2*(nxxt+1)                                     3d27s20
        iostmp=ibcoff                                                   3d27s20
        ipvts=iostmp+nbassmall*nbassmall                                3d27s20
        ibcoff=ipvts+nbassmall                                          3d27s20
        call enough('makeguess.  4',bc,ibc)
        do i=0,nbassmall-1                                              3d27s20
         iad1=ioff+nxxt*i                                               3d27s20
         iad2=iostmp+nbassmall*i                                        3d27s20
         do j=0,nbassmall-1                                             3d27s20
          bc(iad2+j)=ovr(iad1+j)                                        3d27s20
         end do                                                         3d27s20
        end do                                                          3d27s20
        call lusolv(bc(iostmp),nbassmall,nbassmall,dum,1,1,ibc(ipvts),  3d27s20
     $       ierr,1)                                                    3d27s20
        if(ierr.ne.0)stop                                               3d27s20
        do isb=1,nsymbg                                                 3d27s20
         if(morbp(isb).gt.0)then                                        3d27s20
          if(idorelg.eq.0)then                                          3d27s20
           nrow=morbp(isb)                                              3d27s20
           isoff=ivoff                                                  3d27s20
          else                                                          3d27s20
           nrow=morbp(isb)*2                                            3d27s20
           isoff=ivoff+morbp(isb)                                       3d27s20
          end if                                                        3d27s20
          itmp=ibcoff                                                   3d27s20
          ibcoff=itmp+nbaslarge*morbp(isb)                               3d27s20
          call enough('makeguess.  5',bc,ibc)
          call dgemm('n','n',nbaslarge,morbp(isb),morbp(isb),1d0,
     $         bc(iovs(1,isb)),nbaslarge,vgues(ivoff),nrow,0d0,         3d27s20
     $         bc(itmp),nbaslarge,
     d' makeguess.  1')
          call lusolv(bc(ioltmp),nbaslarge,nbaslarge,bc(itmp),nbaslarge,3d27s20
     $         morbp(isb),ibc(ipvtl),ierr,2)                            3d27s20
          ivoff=ivoff+nrow*morbp(isb)
          itmps=ibcoff                                                  3d27s20
          ibcoff=itmps+nbassmall*morbp(isb)*3                           3d27s20
          do ixyz=1,3                                                   3d27s20
           jtmps=itmps+nbassmall*morbp(isb)*(ixyz-1)                    3d27s20
           call dgemm('n','n',nbassmall,morbp(isb),morbp(isb),factr,    3d27s20
     $          bc(iovs(1+ixyz,isb)),nbassmall,vgues(isoff),nrow,0d0,   3d27s20
     $          bc(jtmps),nbassmall,                                     3d27s20
     d' makeguess.  2')
           call lusolv(bc(iostmp),nbassmall,nbassmall,bc(jtmps),
     $          nbassmall,morbp(isb),ibc(ipvts),ierr,2)                          3d27s20
          end do                                                        3d27s20
          do i=0,morbp(isb)-1                                           3d27s20
           if(ibc(ibcode(isb)+i).eq.1)then                              3d27s20
            iuse=ioff1c                                                 3d27s20
            ioff1c=ioff1c+2*nxxt                                        3d27s20
           else if(ibc(ibcode(isb)+i).eq.2)then                         3d27s20
            iuse=ioff2c                                                 3d27s20
            ioff2c=ioff2c+2*nxxt                                        3d27s20
           else if(ibc(ibcode(isb)+i).eq.3)then                         3d27s20
            iuse=ioff3c                                                 3d27s20
            ioff3c=ioff3c+2*nxxt                                        3d27s20
           else                                                         3d27s20
            iuse=ioff0c                                                 3d27s20
            ioff0c=ioff0c+2*nxxt                                        3d27s20
           end if                                                       3d27s20
           jtmp=itmp+nbaslarge*i                                        3d27s20
           do ib=0,nbaslarge-1                                          3d27s20
            vecr(iuse+ib)=bc(jtmp+ib)                                   3d27s20
           end do                                                       3d27s20
           juse=iuse+nbaslarge                                          3d27s20
           do ib=0,nbaslarge-1                                          3d27s20
            vecr(juse+ib)=0d0                                           3d27s20
           end do                                                       3d27s20
           juse=juse+nbaslarge                                          3d27s20
           jtmpsx=itmps+nbassmall*i                                     3d27s20
           jtmpsy=itmps+nbassmall*(i+morbp(isb))                        3d27s20
           jtmpsz=itmps+nbassmall*(i+morbp(isb)*2)                       3d27s20
           jusei=juse+nxx2                                              3d27s20
           do ib=0,nbassmall-1                                          3d27s20
            vecr(jusei+ib)=-bc(jtmpsz+ib)                                3d27s20
           end do                                                       3d27s20
           juse=juse+nbassmall                                          3d27s20
           jusei=juse+nxx2                                              3d27s20
           do ib=0,nbassmall-1                                          3d27s20
            vecr(jusei+ib)=-bc(jtmpsx+ib)                               3d27s20
            vecr(juse+ib)=bc(jtmpsy+ib)                                 3d27s20
           end do                                                       3d27s20
           iuse=iuse+nxxt                                               3d27s20
           do ib=0,nbaslarge-1                                          3d27s20
            vecr(iuse+ib)=0d0                                           3d27s20
           end do                                                       3d27s20
           juse=iuse+nbaslarge                                          3d27s20
           do ib=0,nbaslarge-1
            vecr(juse+ib)=bc(jtmp+ib)                                   3d27s20
           end do                                                       3d27s20
           juse=juse+nbaslarge                                          3d27s20
           jusei=juse+nxx2                                              3d27s20
           do ib=0,nbassmall-1                                          3d27s20
            vecr(juse+ib)=-bc(jtmpsy+ib)                                  3d27s20
            vecr(jusei+ib)=-bc(jtmpsx+ib)                                3d27s20
           end do                                                       3d27s20
           jusei=jusei+nbassmall                                        3d27s20
           do ib=0,nbassmall-1                                          3d27s20
            vecr(jusei+ib)=bc(jtmpsz+ib)                                 3d27s20
           end do                                                       3d27s20
          end do                                                        3d27s20
          ibcoff=itmp                                                   3d27s20
         end if                                                         3d27s20
        end do                                                          3d27s20
        itmpr=ibcoff                                                    3d27s20
        itmpi=itmpr+nxxt*nxx                                            3d27s20
        jtmpr=itmpi+nxxt*nxx                                            3d27s20
        jtmpi=jtmpr+nxxt*nxx                                            3d27s20
        ibcoff=jtmpi+nxxt*nxx                                           3d27s20
        call enough('makeguess.  6',bc,ibc)
        call dgemm('n','n',nxxt,nxx,nxxt,1d0,ovr,nxxt,vecr,nxxt,0d0,    3d27s20
     $       bc(itmpr),nxxt,                                            3d27s20
     d' makeguess.  3')
        call dgemm('n','n',nxxt,nxx,nxxt,1d0,ovr,nxxt,vecr(nxx2+1),     3d27s20
     $       nxxt,0d0,bc(itmpi),nxxt,                                   3d27s20
     d' makeguess.  4')
        do i=0,nxx-1                                                    3d27s20
         do j=0,nxxt-1                                                   3d27s20
          ji=itmpr+j+nxxt*i                                              3d27s20
          ij=jtmpr+i+nxx*j                                              3d27s20
          bc(ij)=bc(ji)                                                 3d27s20
          bc(ij+nxx2)=bc(ji+nxx2)                                       3d27s20
         end do                                                         3d27s20
        end do                                                          3d27s20
        call dgemm('n','n',nxx,nxx,nxxt,1d0,bc(jtmpr),nxx,vecr,
     $       nxxt,0d0,bc(itmpr),nxx,
     d' makeguess.  5')
        call dgemm('n','n',nxx,nxx,nxxt,1d0,bc(jtmpi),nxx,
     $       vecr(nxx2+1),nxxt,1d0,bc(itmpr),nxx,
     d' makeguess.  6')
        call dgemm('n','n',nxx,nxx,nxxt,1d0,bc(jtmpi),nxx,vecr,
     $       nxxt,0d0,bc(itmpi),nxx,
     d' makeguess.  7')
        call dgemm('n','n',nxx,nxx,nxxt,-1d0,bc(jtmpr),nxx,
     $       vecr(nxx2+1),nxxt,1d0,bc(itmpi),nxx,
     d' makeguess.  8')
        rms=0d0
        do i=0,nxx-1                                                    3d27s20
         do j=0,i-1                                                     3d27s20
          ji=itmpr+j+nxx*i
          ij=itmpr+i+nxx*j
          rms=rms+bc(ji)**2+bc(ij)**2
         end do
         ii=itmpr+i+nxx*i
         rms=rms+(bc(ii)-1d0)**2
        end do
        rms=sqrt(rms/dfloat(nxx*nxx))
        write(6,*)('real part: '),rms
        call prntm2(bc(itmpr),nxx,nxx,nxx)
        rms=0d0
        do i=0,nxx*nxx-1                                                3d27s20
         rms=rms+bc(itmpi+i)**2                                         3d27s20
        end do                                                          3d27s20
        rms=sqrt(rms/dfloat(nxx*nxx))                                   3d27s20
        write(6,*)('imag part: '),rms                                   3d27s20
        call prntm2(bc(itmpi),nxx,nxx,nxx)
        ibcoff=itmpr                                                    3d27s20
        return                                                          3d27s20
       end if                                                           3d27s20
       if(ldebug)write(6,*)('so basis is different ')                   1d24s23
       if(nsymb.ne.nsymbg)ircodecx=4                                    1d4s23
       if(ldebug)                                                       1d24s23
     $ write(6,*)('compute overlap matrix between different symmetries')10d23s17
c
c     compute "center of charge" for two different geometries in case
c     a simple shift can put them better in alignment.
c
       do i=1,3                                                         5d14s19
        cartb(i)=0d0                                                    5d14s19
        cartk(i)=0d0                                                    5d14s19
       end do                                                           5d14s19
       cofctw=0d0                                                       5d14s19
       cofctg=0d0                                                       5d14s19
       nqgg=0                                                           4d28s21
       nqwg=0                                                           4d28s21
       iqwg=ibcoff                                                      4d28s21
       iqgg=iqwg+3*ngaus                                                4d28s21
       ibcoff=iqgg+3*ngausg                                             4d28s21
       call enough('makeguess.  7',bc,ibc)
       do iwg=1,ngaus                                                   5d14s19
        arg8=bdat(iwg,8)                                                10d26s17
        japair=iarg8                                                    10d26s17
        do i=1,nqwg                                                     4d28s21
         rms=0d0                                                        4d28s21
         jqwg=iqwg+3*(i-1)-1                                            4d28s21
         do ixyz=1,3                                                    4d28s21
          rms=rms+(bdat(iwg,4+ixyz)-bc(jqwg+ixyz))**2                   4d28s21
         end do                                                         4d28s21
         rms=sqrt(rms/3d0)                                              4d28s21
         if(rms.lt.1d-10)go to 88                                       4d28s21
        end do                                                          4d28s21
        jqwg=iqwg+3*nqwg-1                                              4d28s21
        do ixyz=1,3                                                     4d28s21
         bc(jqwg+ixyz)=bdat(iwg,4+ixyz)                                 4d28s21
        end do                                                          4d28s21
        nqwg=nqwg+1                                                     4d28s21
   88   continue                                                        4d28s21
        do ixyz=1,3                                                     5d14s19
         cartb(ixyz)=cartb(ixyz)+bdat(iwg,4+ixyz)*bdat(iwg,2)           5d14s19
        end do                                                          5d14s19
        cofctw=cofctw+bdat(iwg,2)                                       5d14s19
        if(iapair(1,japair).gt.0)then                                   5d14s19
         do i=1,nqwg                                                     4d28s21
          rms=0d0                                                        4d28s21
          jqwg=iqwg+3*(i-1)-1                                           4d28s21
          do ixyz=1,3                                                    4d28s21
           rms=rms+(bdat(iwg,4+ixyz)*dfloat(isym(ixyz,iapair(2,japair)))4d28s21
     $          -bc(jqwg+ixyz))**2                                      4d28s21
          end do                                                         4d28s21
          rms=sqrt(rms/3d0)                                              4d28s21
          if(rms.lt.1d-10)go to 89                                       4d28s21
         end do                                                          4d28s21
         jqwg=iqwg+3*nqwg-1                                             4d28s21
         do ixyz=1,3                                                     4d28s21
          bc(jqwg+ixyz)=bdat(iwg,4+ixyz)                                4d28s21
     $         *dfloat(isym(ixyz,iapair(2,japair)))                     4d28s21
         end do                                                          4d28s21
         nqwg=nqwg+1                                                     4d28s21
   89    continue                                                        4d28s21
         do ixyz=1,3                                                     5d14s19
          cartb(ixyz)=cartb(ixyz)+bdat(iwg,4+ixyz)*bdat(iwg,2)          5d14s19
     $         *dfloat(isym(ixyz,iapair(2,japair)))                     5d14s19
         end do                                                          5d14s19
         cofctw=cofctw+bdat(iwg,2)                                       5d14s19
        end if                                                          5d14s19
       end do                                                           5d14s19
       cofctw=1d0/cofctw                                                5d14s19
       do ixyz=1,3                                                      5d14s19
        cartb(ixyz)=cartb(ixyz)*cofctw                                  5d14s19
       end do                                                           5d14s19
       if(ldebug)then                                                   1d24s23
        write(6,*)('"center of charge" for want '),cofctw
        call prntm2(cartb,1,3,1)
        write(6,*)('number of centers = '),nqwg,iqwg                          4d28s21
       end if                                                           1d24s23
       do igg=1,ngausg                                                   5d14s19
        arg8=bdatg(igg,8)                                                10d26s17
        japairg=iarg8                                                    10d26s17
        do i=1,nqgg                                                     4d28s21
         rms=0d0                                                        4d28s21
         jqgg=iqgg+3*(i-1)-1                                            4d28s21
         do ixyz=1,3                                                    4d28s21
          rms=rms+(bdatg(igg,4+ixyz)-bc(jqgg+ixyz))**2                   4d28s21
         end do                                                         4d28s21
         rms=sqrt(rms/3d0)                                              4d28s21
         if(rms.lt.1d-10)go to 87                                       4d28s21
        end do                                                          4d28s21
        jqgg=iqgg+3*nqgg-1                                              4d28s21
        do ixyz=1,3                                                     4d28s21
         bc(jqgg+ixyz)=bdatg(igg,4+ixyz)                                 4d28s21
        end do                                                          4d28s21
        nqgg=nqgg+1                                                     4d28s21
   87   continue                                                        4d28s21
        do ixyz=1,3                                                     5d14s19
         cartk(ixyz)=cartk(ixyz)+bdatg(igg,4+ixyz)*bdatg(igg,2)           5d14s19
        end do                                                          5d14s19
        cofctg=cofctg+bdat(igg,2)                                       5d14s19
        if(iapairg(1,japairg).gt.0)then                                   5d14s19
         do i=1,nqgg                                                     4d28s21
          rms=0d0                                                        4d28s21
          jqgg=iqgg+3*(i-1)-1                                           4d28s21
          do ixyz=1,3                                                    4d28s21
           rms=rms+(bdatg(igg,4+ixyz)                                   4d28s21
     $         *dfloat(isymg(ixyz,iapairg(2,japairg)))-bc(jqgg+ixyz))**24d28s21
          end do                                                         4d28s21
          rms=sqrt(rms/3d0)                                              4d28s21
          if(rms.lt.1d-10)go to 86                                       4d28s21
         end do                                                          4d28s21
         jqgg=iqgg+3*nqgg-1                                             4d28s21
         do ixyz=1,3                                                     4d28s21
          bc(jqgg+ixyz)=bdatg(igg,4+ixyz)                               4d28s21
     $         *dfloat(isymg(ixyz,iapairg(2,japairg)))                  4d28s21
         end do                                                          4d28s21
         nqgg=nqgg+1                                                     4d28s21
   86    continue                                                        4d28s21
         do ixyz=1,3                                                     5d14s19
          cartk(ixyz)=cartk(ixyz)+bdatg(igg,4+ixyz)*bdatg(igg,2)          5d14s19
     $         *dfloat(isymg(ixyz,iapairg(2,japairg)))                     5d14s19
         end do                                                          5d14s19
         cofctg=cofctg+bdatg(igg,2)                                       5d14s19
        end if                                                          5d14s19
       end do                                                           5d14s19
       cofctg=1d0/cofctg                                                5d14s19
       do ixyz=1,3                                                      5d14s19
        cartk(ixyz)=cartk(ixyz)*cofctg                                  5d14s19
       end do                                                           5d14s19
       if(ldebug)then                                                   1d24s23
        write(6,*)('"center of charge" for got '),cofctg
        call prntm2(cartk,1,3,1)
        write(6,*)('number of centers = '),nqgg,iqgg                          4d28s21
       end if                                                           1d24s23
       if(nqgg.ne.nqwg)then                                             4d28s21
        write(6,*)('number of centers has changed !!!')                 4d28s21
        write(6,*)('I''m sure this means you accidentally used the'),   4d28s21
     $       (' wrong vector file !')                                    4d28s21
        stop 'bad vector file'                                          4d28s21
       else                                                             4d28s21
        rms=0d0                                                         4d28s21
        do i=0,nqgg*3-1                                                 4d28s21
         rms=rms+(bc(iqgg+i)-bc(iqwg+i))**2                             4d28s21
        end do                                                          4d28s21
        rms=sqrt(rms/dfloat(nqgg*3))                                    4d28s21
        if(ldebug)write(6,*)('rms difference in centers = '),rms        1d24s23
       end if                                                           4d28s21
       if(rms.lt.1d-10.and.nsymb.eq.nsymbg)then                         7d1s22
        do ixyz=1,3                                                     4d28s21
         shift(ixyz)=0d0                                                4d28s21
        end do                                                          4d28s21
        if(ldebug)                                                      1d24s23
     $  write(6,*)('centers are at the same geometries: no shifts used')1d24s23
c
c     is the got basis a subset of the want basis?
c
        nsame=0                                                         5d5s21
        imapwg=ibcoff                                                   5d5s21
        if(ngausg.lt.ngaus.and.idorel.eq.idorelg)then                   5d5s21
         iwhit=imapwg+ngausg                                            5d5s21
         ibcoff=iwhit+ngaus                                             5d5s21
         call enough('makeguess.  8',bc,ibc)
         jmapwg=imapwg-1                                                5d5s21
         jwhit=iwhit-1                                                  5d5s21
         do iwg=1,ngaus                                                 5d5s21
          ibc(jwhit+iwg)=0                                              5d5s21
         end do                                                         5d5s21
         do igg=1,ngausg                                                5d5s21
          arg8=bdatg(igg,1)                                             5d5s21
          lg=iarg8                                                      5d5s21
          do iwg=1,ngaus                                                5d5s21
           arg8=bdat(iwg,1)                                             5d5s21
           lw=iarg8                                                     5d5s21
           if(lg.eq.lw)then                                             5d5s21
            rms=sqrt((bdatg(igg,2)-bdat(iwg,2))**2                      5d5s21
     $              +(bdatg(igg,5)-bdat(iwg,5))**2                      5d5s21
     $              +(bdatg(igg,6)-bdat(iwg,6))**2                      5d5s21
     $              +(bdatg(igg,7)-bdat(iwg,7))**2)                     5d5s21
            if(rms.lt.1d-10)then                                        5d5s21
             ibc(jmapwg+igg)=iwg                                        5d5s21
             nsame=1                                                    5d5s21
             ibc(jwhit+iwg)=1                                           5d5s21
             go to 442                                                  5d5s21
            end if                                                      5d5s21
           end if                                                       5d5s21
          end do                                                        5d5s21
          nsame=0                                                       5d5s21
          go to 4442                                                    5d5s21
  442     continue                                                      5d5s21
         end do                                                         5d5s21
        end if                                                          5d5s21
 4442   continue                                                        5d5s21
        if(ldebug)write(6,*)('nsame = '),nsame
        if(nsame.ne.0)then                                              5d5s21
         if(ldebug)then                                                 1d24s23
          write(6,*)('new basis is a subset of the old basis')           5d5s21
          write(6,*)
     $  ('very good. Let us now put old basis in storage of new basis.')
         end if                                                         1d24s23
         ioff=0                                                         5d5s21
         ioffg=0                                                        5d5s21
         do isb=1,nsymb                                                 5d5s21
          ivoffw(isb)=ioff                                              5d5s21
          ivoffg(isb)=ioffg                                             5d5s21
          do i=1,nbasisp(isb)*nbasisc(isb)*ncomp                        5d7s21
           vecr(ioff+i)=0d0                                             5d5s21
          end do                                                        5d5s21
          ioff=ioff+nbasisp(isb)*nbasisc(isb)*ncomp                     5d7s21
          if(ldebug)then                                                5d7s21
           write(6,*)('vgues for sym '),isb,ioffg,ioff,nbasdws(isb),
     $         nbasisp(isb),nbasisc(isb)
           call prntm2(vgues(ioffg+1),morbp(isb)*ncomp,morb(isb),       10d11s22
     $          morbp(isb)*ncomp)                                       10d11s22
           write(6,*)('ivoffg,w: '),ivoffg(isb),ivoffw(isb)
          end if                                                        5d7s21
          ioffg=ioffg+morbp(isb)*morb(isb)*ncomp                        10d11s22
         end do                                                         5d5s21
         do igg=1,ngausg                                                5d5s21
          arg8=bdatg(igg,1)                                             5d5s21
          lg=iarg8                                                      5d5s21
          nl=2*lg+1                                                     5d5s21
          arg8=bdatg(igg,8)                                             5d5s21
          istest=iarg8                                                  5d5s21
          iwg=ibc(jmapwg+igg)                                           5d5s21
          arg8=bdatg(igg,4)                                             5d5s21
          id0g=iarg8                                                     5d5s21
          arg8=bdat(iwg,4)                                              5d5s21
          id0w=iarg8                                                    5d5s21
          do i=1,nl                                                     5d5s21
           idg=id0g+i                                                   5d5s21
           iddg=ibstorg(idg)                                            5d5s21
           idw=id0w+i                                                   5d5s21
           iddw=ibstor(idw)                                             5d5s21
           isdg=isstorg(idg)                                            5d5s21
           isdw=isstor(idw)                                             5d5s21
           if(isdg.ne.isdw)then                                         5d5s21
            write(6,*)('oops, symmetries do not match!!! '),isdg,isdw
            stop 'makeguess'                                            5d5s21
           end if                                                       5d5s21
           idddg=ivoffg(isdg)+iddg                                      5d5s21
           idddw=ivoffw(isdg)+iddw                                      5d5s21
           nrw=nbasisp(isdg)*ncomp                                      5d5s21
           if(ldebug)write(6,*)('vecr no. 1')
           do j=1,morb(isdg)                                            5d5s21
            vecr(idddw)=vgues(idddg)                                    5d5s21
            idddg=idddg+morbp(isdg)*ncomp                               10d11s22
            idddw=idddw+nrw                                             5d5s21
           end do                                                       5d5s21
           if(idorel.ne.0)then                                          5d5s21
            idddg=ivoffg(isdg)+iddg+morbp(isdg)                         5d5s21
            idddw=ivoffw(isdg)+iddw+nbasisp(isdg)                       5d5s21
            do j=1,morb(isdg)                                            5d5s21
             vecr(idddw)=vgues(idddg)                                   5d5s21
             idddg=idddg+morbp(isdg)*ncomp                              10d11s22
             idddw=idddw+nrw                                             5d5s21
            end do                                                       5d5s21
           end if                                                       5d5s21
          end do                                                        5d5s21
         end do                                                         5d5s21
         if(ldebug)write(6,*)('now let us add in guess for remaining')  1d24s23
         nvguess=-ncomp                                                 10d11s22
         ircode=2                                                       10d12s22
         if(ldebug)write(6,*)('ircode b')                               1d24s23
         do isb=1,nsymb                                                 5d5s21
          ivoffg(isb)=morb(isb)                                         5d5s21
         end do                                                         5d5s21
         do iwg=1,ngaus                                                 5d5s21
          if(ibc(jwhit+iwg).eq.0)then                                   5d5s21
           arg8=bdat(iwg,1)                                             5d5s21
           lw=iarg8                                                      5d5s21
           nl=2*lw+1                                                     5d5s21
           arg8=bdat(iwg,8)                                             5d5s21
           istest=iarg8                                                  5d5s21
           if(iapair(1,istest).ne.0)nl=nl*2                             5d5s21
           if(bdat(iwg,2).gt.1d0.and.idorel.ne.0)then                   5d5s21
            write(6,*)
     $   ('this exponential parameter is too large'),
     $          (' for nr limit contraction????'),bdat(iwg,2)           5d5s21
            stop 'makeguess'                                            5d5s21
           end if
           arg8=bdat(iwg,4)
           id0w=iarg8
           do i=1,nl
            idw=id0w+i
            iddw=ibstor(idw)
            isb=isstor(idw)
            idddw=ivoffw(isb)+iddw+ncomp*nbasisp(isb)*ivoffg(isb)       5d7s21
            vecr(idddw)=1d0                                               5d5s21
            if(idorel.ne.0)then                                         5d5s21
             idddw=idddw+nbasisp(isb)                                   5d5s21
             vecr(idddw)=1d0                                            5d5s21
            end if                                                      5d5s21
            ivoffg(isb)=ivoffg(isb)+1                                   5d5s21
           end do                                                       5d5s21
          end if                                                        5d5s21
         end do                                                         5d5s21
         ioffw=0                                                        5d5s21
         do isb=1,nsymb                                                 5d5s21
          nrw=nbasisp(isb)*ncomp                                        5d5s21
          if(nbasisp(isb)*ncomp.eq.nbasisc(isb))then                    5d5s21
           iotmp=ibcoff                                                 5d5s21
           iptmp=iotmp+nrw*nrw                                          5d5s21
           iqtmp=iptmp+nrw*nrw                                          5d5s21
           ibcoff=iqtmp+nrw*nrw                                         5d5s21
           call enough('makeguess.  9',bc,ibc)
           jotmp=iotmp-1                                                5d5s21
           nn=(nrw*(nrw+1))/2                                           5d5s21
           do i=1,nn                                                    5d5s21
            bc(jotmp+i)=ovr(ioffw+i)                                    5d5s21
           end do                                                       5d5s21
           call square(bc(iotmp),nrw)                                   5d5s21
           jptmp=iptmp-1                                                5d5s21
           if(ldebug)write(6,*)('vecr no. 2')                           1d24s23
           do i=morb(isb),ivoffg(isb)-1                                 5d5s21
            ii=ivoffw(isb)+nrw*i                                        5d5s21
            call dgemm('n','n',nrw,i,nrw,1d0,bc(iotmp),nrw,             5d5s21
     $           vecr(ivoffw(isb)+1),nrw,0d0,bc(iptmp),nrw,             5d5s21
     d' makeguess.  9')
            do j=0,i-1                                                  5d5s21
             jj=ivoffw(isb)+nrw*j                                       5d5s21
             jjj=jptmp+nrw*j                                            5d5s21
             dot=0d0                                                    5d5s21
             do k=1,nrw                                                 5d5s21
              dot=dot+vecr(ii+k)*bc(jjj+k)                              5d5s21
             end do                                                     5d5s21
             do k=1,nrw                                                 5d5s21
              vecr(ii+k)=vecr(ii+k)-vecr(jj+k)*dot                      5d5s21
             end do                                                     5d5s21
            end do                                                      5d5s21
            call dgemm('n','n',nrw,1,nrw,1d0,bc(iotmp),nrw,             5d5s21
     $           vecr(ii+1),nrw,0d0,bc(iptmp),nrw,                      5d5s21
     d' makeguess. 10')
            dot=0d0                                                     5d5s21
            do k=1,nrw                                                  5d5s21
             dot=dot+vecr(ii+k)*bc(jptmp+k)                             5d5s21
            end do                                                      5d5s21
            xnorm=1d0/sqrt(dot)                                         5d5s21
            do k=1,nrw                                                  5d5s21
             vecr(ii+k)=vecr(ii+k)*xnorm                                5d5s21
            end do                                                      5d5s21
           end do                                                       5d5s21
           ibcoff=iotmp                                                 5d5s21
           ioffw=ioffw+nrw*nrw                                          5d5s21
          else                                                          5d7s21
           if(ldebug)then                                               5d7s21
            write(6,*)
     $          ('good thing we have original overlap matrix ...'),
     $          ivoffw(isb)
            call prntm2(scopy(ioffw+1),nrw,nrw,nrw)
            call prntm2(vecr(ivoffw(isb)+1),nrw,ivoffg(isb),nrw)
           end if                                                       5d7s21
           iptmp=ibcoff                                                 5d7s21
           iqtmp=iptmp+nrw*nrw                                          5d5s21
           ibcoff=iqtmp+nrw*nrw                                         5d5s21
           call enough('makeguess. 10',bc,ibc)
           if(ldebug)then                                               5d7s21
            call dgemm('n','n',nrw,ivoffg(isb),nrw,1d0,
     $          scopy(ioffw+1),nrw,vecr(ivoffw(isb)+1),nrw,0d0,
     $          bc(iptmp),nrw,
     d' makeguess. 11')
            call dgemm('t','n',ivoffg(isb),ivoffg(isb),nrw,1d0,
     $          vecr(ivoffw(isb)+1),nrw,bc(iptmp),nrw,0d0,
     $          bc(iqtmp),ivoffg(isb),
     d' makeguess. 12')
            write(6,*)('ortho test ...')
            call prntm2(bc(iqtmp),ivoffg(isb),ivoffg(isb),ivoffg(isb))
           end if                                                       5d7s21
           jptmp=iptmp-1                                                5d5s21
           if(ldebug)write(6,*)('vecr no. 3')                           1d24s23
           do i=morb(isb),ivoffg(isb)-1                                 5d5s21
            ii=ivoffw(isb)+nrw*i                                        5d5s21
            call dgemm('n','n',nrw,i,nrw,1d0,scopy(ioffw+1),nrw,        5d7s21
     $           vecr(ivoffw(isb)+1),nrw,0d0,bc(iptmp),nrw,             5d5s21
     d' makeguess. 13')
            do j=0,i-1                                                  5d5s21
             jj=ivoffw(isb)+nrw*j                                       5d5s21
             jjj=jptmp+nrw*j                                            5d5s21
             dot=0d0                                                    5d5s21
             do k=1,nrw                                                 5d5s21
              dot=dot+vecr(ii+k)*bc(jjj+k)                              5d5s21
             end do                                                     5d5s21
             do k=1,nrw                                                 5d5s21
              vecr(ii+k)=vecr(ii+k)-vecr(jj+k)*dot                      5d5s21
             end do                                                     5d5s21
            end do                                                      5d5s21
            call dgemm('n','n',nrw,1,nrw,1d0,scopy(ioffw+1),nrw,        5d7s21
     $           vecr(ii+1),nrw,0d0,bc(iptmp),nrw,                      5d5s21
     d' makeguess. 14')
            dot=0d0                                                     5d5s21
            do k=1,nrw                                                  5d5s21
             dot=dot+vecr(ii+k)*bc(jptmp+k)                             5d5s21
            end do                                                      5d5s21
            if(ldebug)write(6,*)('i: '),i,dot                           5d7s21
            xnorm=1d0/sqrt(dot)                                         5d5s21
            do k=1,nrw                                                  5d5s21
             vecr(ii+k)=vecr(ii+k)*xnorm                                5d5s21
            end do                                                      5d5s21
           end do                                                       5d5s21
           if(ldebug)then                                               5d5s23
            write(6,*)('vec '),ivoffw(isb)
            call prntm2(vecr(ivoffw(isb)+1),nrw,ivoffg(isb),nrw)
            write(6,*)('new code: ortmp')
            iortmp=ibcoff
            iortmp2=iortmp+nrw*ivoffg(isb)
            ibcoff=iortmp2+ivoffg(isb)**2
            call enough('makeguess.ortmp',bc,ibc)
            call dgemm('n','n',nrw,ivoffg(isb),nrw,1d0,scopy(ioffw+1),
     $          nrw,vecr(ivoffw(isb)+1),nrw,0d0,bc(iortmp),nrw,
     $          'makeguess.ortmp')
            call dgemm('t','n',ivoffg(isb),ivoffg(isb),nrw,1d0,
     $          vecr(ivoffw(isb)+1),nrw,bc(iortmp),nrw,0d0,bc(iortmp2),
     $          ivoffg(isb),'makeguess.ortmp2')
            call prntm2(bc(iortmp2),ivoffg(isb),ivoffg(isb),ivoffg(isb))
           end if                                                       5d5s23
           ibcoff=iptmp                                                 5d7s21
           ioffw=ioffw+nrw*nrw                                          5d5s21
          end if                                                        5d5s21
         end do                                                         5d5s21
         ibcoff=ibcoffo                                                 10d7s22
         return                                                         10d7s22
        end if                                                          5d5s21
       end if                                                           10d7s22
        if(ldebug)write(6,*)('different geometry block ...')
        do ixyz=1,3
         diff=cartk(ixyz)-cartb(ixyz)
         shift(ixyz)=0d0                                                 5d14s19
         if(abs(diff).gt.1d-10)then
          nneg=0                                                         5d14s19
          do isb=1,nsymb                                                 5d14s19
           if(isym(ixyz,isb).lt.0)nneg=nneg+1                            5d14s19
          end do                                                         5d14s19
          if(nneg.gt.0)then                                              5d14s19
           write(6,*)('this shift would break symmetry ')                     5d14s19
          else                                                           5d14s19
           write(6,*)('this shift is consistent with symmetry ')         5d14s19
           shift(ixyz)=diff                                              5d14s19
          end if                                                         5d14s19
         end if                                                          5d14s19
        end do
       xnan=-2d0
       do isw=1,nsymb                                                   10d23s17
        nw=nbasisp(isw)                                                 5d9s19
        if(idorel.ne.0)nw=nw*2                                          6d3s22
        do isg=1,nsymbg                                                 10d23s17
         ng=morbp(isg)                                                  5d9s19
         if(idorel.ne.0)ng=ng*2                                         6d3s22
         iovs(isg,isw)=ibcoff                                           10d23s17
         ibcoff=ibcoff+ng*nw                                            10d23s17
         if(ldebug)then
          write(6,*)('for syms '),isg,isw,('ovs has dims '),ng,nw
         end if
         call enough('makeguess. 11',bc,ibc)
         do i=0,ng*nw-1                                                 10d23s17
          bc(iovs(isg,isw)+i)=0d0                                       10d23s17
         end do                                                         10d23s17
        end do                                                          10d23s17
       end do                                                           10d23s17
       do iwg=1,ngaus                                                   10d23s17
        arg8=bdat(iwg,1)                                                10d23s17
        lw=iarg8                                                        10d23s17
        arg8=bdat(iwg,4)                                                10d23s17
        ipw=iarg8                                                       10d23s17
        nw=2*lw+1                                                       10d23s17
        arg8=bdat(iwg,8)                                                10d26s17
        japair=iarg8                                                    10d26s17
        nw0=nw                                                          10d26s17
        if(iapair(1,japair).gt.0)nw=nw*2                                10d26s17
        nwu=nw                                                          6d3s22
        nw0u=nw0                                                        6d3s22
        if(idorel.ne.0)then
         nwu=nw*2                                                       6d3s22
         nw0u=nw0*2                                                     6d3s22
        end if                                                          6d3s22
        do igg=1,ngausg                                                 10d23s17
         arg8=bdatg(igg,1)                                              10d23s17
         lg=iarg8                                                       10d23s17
         arg8=bdatg(igg,4)                                              10d23s17
         ipg=iarg8                                                      10d23s17
         arg8=bdatg(igg,8)                                              10d27s17
         kapair=iarg8                                                   10d27s17
         ng=2*lg+1                                                      10d23s17
         ng0=ng                                                         10d27s17
         if(iapairg(1,kapair).gt.0)ng=ng*2                              10d27s17
         ngu=ng                                                         6d3s22
         ng0u=ng0                                                       6d3s22
         if(idorel.ne.0)then                                            6d3s22
          ngu=ng*2                                                      6d3s22
          ng0u=ng0*2                                                    6d3s22
         end if                                                         6d3s22
         itmp=ibcoff                                                    10d26s17
         ibcoff=itmp+ngu*nwu                                            6d3s22
         call enough('makeguess. 12',bc,ibc)
         do iz=itmp,ibcoff-1                                            6d3s22
          bc(iz)=0d0                                                    6d3s22
         end do                                                         6d3s22
         npass=1                                                        10d27s17
         if(ng.ne.ng0)npass=2                                           10d27s17
         do ipass=1,npass                                               10d27s17
          if(ipass.eq.1)then                                            10d27s17
           cartb(1)=bdatg(igg,5)                                        10d27s17
           cartb(2)=bdatg(igg,6)                                        10d27s17
           cartb(3)=bdatg(igg,7)                                        10d27s17
           ibadd=0                                                      10d27s17
          else                                                          10d27s17
           cartb(1)=bdatg(igg,5)*dfloat(isymg(1,iapairg(2,kapair)))     10d27s17
           cartb(2)=bdatg(igg,6)*dfloat(isymg(2,iapairg(2,kapair)))     10d27s17
           cartb(3)=bdatg(igg,7)*dfloat(isymg(3,iapairg(2,kapair)))     10d27s17
           ibadd=ng0u                                                   6d3s22
          end if                                                        10d27s17
          do ixyz=1,3                                                   5d14s19
           cartk(ixyz)=bdat(iwg,4+ixyz)+shift(ixyz)                     5d14s19
          end do                                                        5d14s19
          call onei(lg,bdatg(igg,2),bdatg(igg,3),cartb,                  10d27s17
     $        cartb(2),cartb(3),ipg,lw,bdat(iwg,2),bdat(iwg,3),         10d27s17
     $        cartk,cartk(2),cartk(3),ipw,dum,dum,idum,ib1,             5d14s19
     $        ib2,bc,ibc)                                               11d9s22
          do ik=0,nw0-1                                                  10d26s17
           do ib=0,ng0-1                                                  10d26s17
            ibp=ib+ibadd                                                 10d27s17
            iadfrm=ib1+ik+nw0*ib                                         10d26s17
            iadto=itmp+ibp+ngu*ik                                       6d3s22
            bc(iadto)=bc(iadfrm)                                         10d26s17
           end do                                                        10d26s17
          end do                                                         10d26s17
          if(idorel.ne.0)then                                           6d3s22
           do ik=0,nw0-1                                                  10d26s17
            ikp=ik+nw0                                                  6d3s22
            do ib=0,ng0-1                                                  10d26s17
             ibp=ib+ibadd+ng0                                           6d3s22
             iadfrm=ib2+ik+nw0*ib                                         10d26s17
             iadto=itmp+ibp+ngu*ikp                                     6d3s22
             bc(iadto)=ascale*bc(iadfrm)                                6d3s22
            end do                                                        10d26s17
           end do                                                         10d26s17
          end if                                                        6d3s22
          if(nw0.ne.nw)then                                              10d26s17
           cartk(1)=bdat(iwg,5)*dfloat(isym(1,iapair(2,japair)))          10d26s17
     $         +shift(1)                                                5d14s19
           cartk(2)=bdat(iwg,6)*dfloat(isym(2,iapair(2,japair)))          10d26s17
     $         +shift(2)                                                5d14s19
           cartk(3)=bdat(iwg,7)*dfloat(isym(3,iapair(2,japair)))          10d26s17
     $         +shift(3)                                                5d14s19
           ikadd=nw0u                                                   6d3s22
           call onei(lg,bdatg(igg,2),bdatg(igg,3),cartb,                10d27s17
     $        cartb(2),cartb(3),ipg,lw,bdat(iwg,2),bdat(iwg,3),         10d27s17
     $        cartk,cartk(2),cartk(3),ipw,dum,dum,idum,ib1,             10d26s17
     $        ib2,bc,ibc)                                               11d9s22
           do ik=0,nw0-1                                                10d27s17
            ikp=ik+ikadd                                                6d3s22
            do ib=0,ng0-1                                               10d27s17
             ibp=ib+ibadd                                               10d27s17
             iadfrm=ib1+ik+nw0*ib                                         10d26s17
             iadto=itmp+ibp+ngu*ikp                                     6d3s22
             bc(iadto)=bc(iadfrm)                                         10d26s17
            end do                                                        10d26s17
           end do                                                         10d26s17
           if(idorel.ne.0)then                                          6d3s22
            do ik=0,nw0-1                                                10d27s17
             ikp=ik+ikadd+nw0                                           6d3s22
             do ib=0,ng0-1                                               10d27s17
              ibp=ib+ibadd+ng0                                          6d3s22
              iadfrm=ib2+ik+nw0*ib                                         10d26s17
              iadto=itmp+ibp+ngu*ikp                                     6d3s22
              bc(iadto)=ascale*bc(iadfrm)                               6d3s22
             end do                                                        10d26s17
            end do                                                         10d26s17
           end if                                                       6d3s22
           do ik=0,nw0u-1                                               6d3s22
            ikp=ik+ikadd                                                6d3s22
            iad1=itmp+ik*ngu+ibadd                                       10d27s17
            iad2=itmp+ikp*ngu+ibadd                                      10d27s17
            do ib=0,ng0u-1                                               6d3s22
             sum=srh*(bc(iad1+ib)+bc(iad2+ib))                           10d26s17
             diff=srh*(-bc(iad1+ib)+bc(iad2+ib))                         10d26s17
             bc(iad1+ib)=sum                                             10d26s17
             bc(iad2+ib)=diff                                            10d26s17
            end do                                                       10d26s17
           end do                                                        10d26s17
          end if                                                        10d27s17
         end do                                                         10d27s17
         if(npass.eq.2)then                                             10d27s17
          do ik=0,nwu-1                                                  10d27s17
           iad1=itmp+ik*ngu                                             6d3s22
           iad2=iad1+ng0u                                               6d3s22
           do ib=0,ng0u-1                                               6d3s22
            sum=srh*(bc(iad1+ib)+bc(iad2+ib))                            10d27s17
            diff=srh*(-bc(iad1+ib)+bc(iad2+ib))                         10d27s17
            bc(iad1+ib)=sum                                             10d27s17
            bc(iad2+ib)=diff                                            10d27s17
           end do                                                       10d27s17
          end do                                                        10d27s17
         end if                                                         10d27s17
c
c     if rel, overlap is now ls 1 ls 2 , while what we want is          6d3s22
c     l1l2 s1s2 ...                                                     6d3s22
c                                                                       6d3s22
         if(idorel.ne.0)then                                            6d3s22
          icpy=ibcoff                                                   6d3s22
          ibcoff=icpy+ngu*nwu                                           6d3s22
          call enough('makeguess. 13',bc,ibc)
          if(nw.ne.nw0)then                                             6d3s22
           do ik1=0,nw0-1                                                 6d3s22
            ik2=ik1+nw0                                                  6d3s22
            ik3=ik2+nw0                                                  6d3s22
            ik4=ik3+nw0                                                  6d3s22
            iad1=itmp+ngu*ik1                                            6d3s22
            iad2=itmp+ngu*ik2                                            6d3s22
            iad3=itmp+ngu*ik3                                            6d3s22
            iad4=itmp+ngu*ik4                                            6d3s22
            jad1=icpy+ngu*ik1                                            6d3s22
            jad2=icpy+ngu*ik2                                            6d3s22
            jad3=icpy+ngu*ik3                                            6d3s22
            jad4=icpy+ngu*ik4                                            6d3s22
            do ib=0,ngu-1                                                6d3s22
             bc(jad1+ib)=bc(iad1+ib)                                     6d3s22
             bc(jad2+ib)=bc(iad3+ib)                                     6d3s22
             bc(jad3+ib)=bc(iad2+ib)                                     6d3s22
             bc(jad4+ib)=bc(iad4+ib)                                     6d3s22
            end do                                                       6d3s22
           end do                                                        6d3s22
          else
           do ibk=0,ngu*nwu-1                                           6d3s22
            bc(icpy+ibk)=bc(itmp+ibk)                                   6d3s22
           end do                                                       6d3s22
          end if                                                        6d3s22
          if(ng0.ne.ng)then
           do ik=0,nwu-1                                                 6d3s22
            iad1=icpy+ngu*ik                                            6d3s22
            iad2=iad1+ng0                                                6d3s22
            iad3=iad2+ng0                                                6d3s22
            iad4=iad3+ng0                                                6d3s22
            jad1=itmp+ngu*ik                                            6d3s22
            jad2=jad1+ng0                                                6d3s22
            jad3=jad2+ng0                                                6d3s22
            jad4=jad3+ng0                                                6d3s22
            do ib=0,ng0-1                                                6d3s22
             bc(jad1+ib)=bc(iad1+ib)                                     6d3s22
             bc(jad2+ib)=bc(iad3+ib)                                     6d3s22
             bc(jad3+ib)=bc(iad2+ib)                                     6d3s22
             bc(jad4+ib)=bc(iad4+ib)                                     6d3s22
            end do                                                       6d3s22
           end do                                                        6d3s22
          else                                                          6d3s22
           do ibk=0,ngu*nwu-1                                           6d3s22
            bc(itmp+ibk)=bc(icpy+ibk)                                   6d3s22
           end do                                                       6d3s22
          end if                                                        6d3s22
          ibcoff=icpy                                                   6d3s22
         end if                                                         6d3s22
         do ik=0,nw-1                                                   10d26s17
          isk=isstor(ipw+ik+1)                                          10d23s17
          ibk=ibstor(ipw+ik+1)                                          10d23s17
          do ib=0,ng-1                                                  10d23s17
           isb=isstorg(ipg+ib+1)                                        10d23s17
           ndd=morbp(isb)                                               6d3s22
           if(idorel.ne.0)ndd=ndd*2                                     6d3s22
           ibb=ibstorg(ipg+ib+1)                                        10d23s17
           iadfrm=itmp+ib+ngu*ik                                         10d26s17
           iadto=iovs(isb,isk)+ibb-1+ndd*(ibk-1)                        6d3s22
           bc(iadto)=bc(iadfrm)                                         10d23s17
          end do                                                        10d23s17
         end do                                                         10d23s17
         if(idorel.ne.0)then                                            6d3s22
          do ik=0,nw-1                                                   10d26s17
           isk=isstor(ipw+ik+1)                                          10d23s17
           ibk=ibstor(ipw+ik+1)+nbasisp(isk)                            6d3s22
           do ib=0,ng-1                                                  10d23s17
            isb=isstorg(ipg+ib+1)                                        10d23s17
            ndd=morbp(isb)*2                                            6d3s22
            ibb=ibstorg(ipg+ib+1)+morbp(isb)                            6d3s22
            iadfrm=itmp+ib+ng+ngu*(ik+nw)                               6d3s22
            iadto=iovs(isb,isk)+ibb-1+ndd*(ibk-1)                       6d3s22
            bc(iadto)=2d0*bc(iadfrm)                                    6d3s22
           end do                                                        10d23s17
          end do                                                         10d23s17
         end if                                                         6d3s22
         ibcoff=itmp                                                    10d26s17
        end do                                                          10d23s17
       end do                                                           10d23s17
       if(ldebug)then                                                   7d5s23
        write(6,*)('raw overlaps ')
        do isg=1,nsymbg
         ng=morbp(isg)*ncomp
         do isw=1,nsymb
          nw=nbasisp(isw)*ncomp
          write(6,*)('for symmetrys g,w '),isg,isw
          call prntm2(bc(iovs(isg,isw)),ng,nw,ng)                        6d3s22
         end do                                                          6d3s22
        end do                                                           6d3s22
       end if                                                           7d5s23
       iarg8=loc(bc(ibcoff))-loc(bdat)                                  4d27s21
       iarg8=iarg8/8                                                    4d27s21
       ipoint2bdat=ibcoff-iarg8                                         4d27s21
       morbpsum=0                                                       4d27s21
       do isg=1,nsymbg                                                  4d27s21
        morbpsum=morbpsum+morbp(isg)                                    4d27s21
       end do                                                           4d27s21
       itmpovr=ibcoff                                                   4d27s21
       jtmpovr=itmpovr                                                  4d27s21
       morbpsum=morbpsum*ncomp                                          10d7s22
       do isw=1,nsymb                                                   4d27s21
        if(nbasisp(isw).gt.0)then                                       4d27s21
         ibcoff=jtmpovr+nbasisp(isw)*ncomp*morbpsum                      4d27s21
         call enough('makeguess. 14',bc,ibc)
         istrt=jtmpovr                                                   4d27s21
         do isg=1,nsymbg                                                 4d27s21
          jovr=iovs(isg,isw)                                             4d27s21
          jstrt=istrt                                                    4d27s21
          do icol=0,nbasisp(isw)*ncomp-1                                 4d27s21
           do irow=0,morbp(isg)*ncomp-1                                 10d7s22
            bc(jstrt+irow)=bc(jovr+irow)                                 4d27s21
           end do                                                        4d27s21
           jstrt=jstrt+morbpsum                                          4d27s21
           jovr=jovr+morbp(isg)*ncomp                                   10d7s22
          end do                                                         4d27s21
          istrt=istrt+morbp(isg)*ncomp                                  10d7s22
         end do                                                          4d27s21
         mrow(isw)=morbpsum                                             4d27s21
         if(ldebug)then
          write(6,*)('all together now ')
          call prntm2(bc(jtmpovr),morbpsum,nbasisp(isw)*ncomp,morbpsum)  10d7s22
         end if
         jtmpovr=jtmpovr+nbasisp(isw)*ncomp*morbpsum                     4d27s21
        end if                                                          4d27s21
       end do                                                           4d27s21
       call contractg(bc(itmpovr),ngaus,ipoint2bdat,ibstor,isstor,      4d27s21
     $       idorel,-1,nbasisp,nbasisc,iapair,0,idum,idum,mrow,0,bc,ibc)11d9s22
       jtmpovr=itmpovr                                                  4d27s21
       do isw=1,nsymb                                                   4d27s21
        if(nbasisc(isw).gt.0)then                                       4d27s21
         if(ldebug)then
          write(6,*)('all together contracted '),isw,nbasisp(isw)
          call prntm2(bc(jtmpovr),morbpsum,nbasisc(isw),morbpsum)        10d7s22
         end if
         istrt=jtmpovr                                                  4d27s21
         do isg=1,nsymbg                                                4d27s21
          jstrt=istrt                                                   4d27s21
          jovr=iovs(isg,isw)                                            4d27s21
          do i=0,nbasisc(isw)-1                                         10d7s22
           do j=0,morbp(isg)*ncomp-1                                    10d7s22
            bc(jovr+j)=bc(jstrt+j)                                      4d27s21
           end do                                                       4d27s21
           jstrt=jstrt+morbpsum                                         4d27s21
           jovr=jovr+morbp(isg)*ncomp                                   10d7s22
          end do                                                        4d27s21
          istrt=istrt+morbp(isg)*ncomp                                  10d7s22
         end do                                                         4d27s21
         jtmpovr=jtmpovr+morbpsum*nbasisp(isw)*ncomp                    10d7s22
        end if                                                          4d27s21
       end do                                                           4d27s21
       ibcoff=itmpovr                                                   10d7s22
       if(ldebug)then                                                   7d5s23
        do isw=1,nsymb
         do isg=1,nsymbg
          if(min(morbp(isg),nbasisc(isw)).gt.0)then                      4d27s21
           write(6,*)('got sym '),isg,(' want sym '),isw
           call prntm2(bc(iovs(isg,isw)),morbp(isg)*ncomp,               10d7s22
     $         nbasisc(isw),morbp(isg)*ncomp)                           10d7s22
          end if                                                         4d27s21
         end do
        end do
       end if                                                           7d5s23
c
c     VTSV=VTS(XSX)V=VT(XSX)SXSXV=VTXS(XSX)SXV=(SXV)T(SXV)
c     we have Sgw. We would like VwTSwgVg=1, Vw and Vg ao vectors.
c     We also require VwTSwwVw=1. Let Uw be ob vectors, i.e. UwTUw=1=
c     UwTXwSwwXwUw=(XwUw)TSww(XwUw), thus Vw=XwwUw. So we want
c     (XwwUw)TSwgVg=1=UwTXwwTSSwgVg, thus Uw=XTSSgwTVg.
c
       ioffw=0                                                          10d23s17
       nbcx=0                                                           10d24s17
       do isg=1,nsymbg                                                  10d24s17
        do i=0,nbasdwsg(isg)-1                                          10d24s17
         nbcx=max(nbcx,ibc(ibcode(isg)+i))                              10d24s17
        end do                                                          10d24s17
       end do
       nbcxp=nbcx+1                                                     10d24s17
       ibcx=ibcoff                                                      10d24s17
       ibcxs=ibcx+nbcxp                                                  10d24s17
       ibcxn=ibcxs+nbcxp                                                 10d24s17
       ibcoff=ibcxn+nbcxp                                               10d24s17
       call enough('makeguess. 15',bc,ibc)
       do i=0,nbcxp-1                                                   10d24s17
        ibc(ibcx+i)=0                                                   10d24s17
       end do                                                           10d24s17
       do isg=1,nsymbg                                                  10d24s17
        do i=0,nbasdwsg(isg)-1                                          10d24s17
         im=ibc(ibcode(isg)+i)                                          10d24s17
         ibc(ibcx+im)=ibc(ibcx+im)+1                                    10d24s17
        end do                                                          10d24s17
       end do                                                           10d24s17
       ivecsav=0                                                        10d24s17
       do isw=1,nsymb                                                   10d23s17
        if(nbasisp(isw).gt.0)then                                       5d14s19
         nbasp=nbasisc(isw)+1                                           4d27s21
         do i=0,nbcxp-1                                                  10d24s17
          ibc(ibcxs+i)=ibcoff                                            10d24s17
          ibcoff=ibcoff+nbasp*ibc(ibcx+i)                                10d24s17
          ibc(ibcxn+i)=0                                                 10d24s17
         end do                                                          10d24s17
         ixi=ibcoff                                                      10d23s17
         ix=ixi+nbasisc(isw)*nbasisc(isw)                               4d27s21
         ieig=ix+nbasisc(isw)*nbasisc(isw)                              4d27s21
         itemp=ieig+nbasisc(isw)                                        4d27s21
         ifv1=itemp+nbasisc(isw)*nbasisc(isw)                           4d27s21
         ibcoff=ifv1+nbasisc(isw)                                       4d27s21
         isi=ibcoff                                                      4d25s18
         ibcoff=isi+nbasisc(isw)*nbasisc(isw)                           4d27s21
         call enough('makeguess. 16',bc,ibc)
         jx=ix                                                           10d23s17
         jcopy=ioffw+1
         do i=1,nbasisc(isw)                                            4d27s21
          do j=1,i                                                       10d23s17
           bc(jx)=ovr(jcopy)                                             10d23s17
           jx=jx+1                                                       10d23s17
           jcopy=jcopy+1                                                 10d23s17
          end do                                                         10d23s17
         end do                                                          10d23s17
         call square(bc(ix),nbasisc(isw))                               4d27s21
         issscpy=ibcoff                                                 4d28s21
         ibcoff=issscpy+nbasisc(isw)*nbasisc(isw)                       4d28s21
         call enough('makeguess. 17',bc,ibc)
         do i=0,nbasisc(isw)*nbasisc(isw)-1                             4d28s21
          bc(issscpy+i)=bc(ix+i)                                        4d28s21
         end do                                                         4d28s21
         ioffw=ioffw+nbasisp(isw)*nbasisp(isw)*ncomp*ncomp              10d11s22
         if(nlzz.ne.0)then                                              4d27s21
          if(nlzz.eq.2)then                                             4d27s21
           do i=0,nbasisc(isw)-1                                         4d27s21
            ibc(ifv1+i)=ibc(iorbsym(isw)+i)                              4d27s21
           end do                                                       4d27s21
          else                                                          4d27s21
           do i=0,nbasisc(isw)-1                                        4d27s21
            ibc(ifv1+i)=ibc(iorbsym(isw)+i)+100*ibc(iorbsymz(isw)+i)    4d27s21
           end do                                                       4d27s21
          end if                                                        4d27s21
          call diagy(nbasisc(isw),bc(ix),bc(ieig),bc(itemp),bc(ifv1),   11d14s22
     $         bc,ibc,0,idum,dum)                                       9d1s23
         else                                                           4d27s21
          call diagx(nbasisc(isw),bc(ix),bc(ieig),bc(itemp),bc(ifv1),   11d14s22
     $         bc,ibc)                                                  11d14s22
         end if                                                         4d27s21
         do i=0,nbasisc(isw)-1                                          4d27s21
          bc(ieig+i)=1d0/sqrt(bc(ieig+i))                                10d23s17
         end do                                                          10d23s17
         do i=0,nbasisc(isw)-1                                          4d27s21
          do j=0,nbasisc(isw)-1                                         4d27s21
           jx=ix+j+nbasisc(isw)*i                                       4d27s21
           jt=itemp+i+nbasisc(isw)*j                                    4d27s21
           bc(jx)=bc(jt)*bc(ieig+j)                                      10d23s17
           kx=isi+j+nbasisc(isw)*i                                      4d27s21
           bc(kx)=bc(jt)*bc(ieig+j)*bc(ieig+j)                           4d25s18
          end do                                                         10d23s17
         end do                                                          10d23s17
         call dgemm('n','n',nbasisc(isw),nbasisc(isw),nbasisc(isw),1d0, 4d27s21
     $       bc(itemp),nbasisc(isw),bc(ix),nbasisc(isw),0d0,            4d27s21
     $       bc(ixi),nbasisc(isw),                                      4d27s21
     d' makeguess. 15')
         call dgemm('n','n',nbasisc(isw),nbasisc(isw),nbasisc(isw),1d0, 4d27s21
     $       bc(itemp),nbasisc(isw),bc(isi),nbasisc(isw),0d0,           4d27s21
     $       bc(ix),nbasisc(isw),                                       4d27s21
     d' makeguess. 16')
        call dgemm('n','n',nbasisc(isw),nbasisc(isw),nbasisc(isw),1d0,
     $       bc(issscpy),nbasisc(isw),bc(ixi),nbasisc(isw),0d0,
     $       bc(ibcoff),nbasisc(isw),
     d' makeguess. 17')
        call dgemm('n','n',nbasisc(isw),nbasisc(isw),nbasisc(isw),1d0,  4d28s21
     $       bc(ixi),nbasisc(isw),bc(ibcoff),nbasisc(isw),0d0,          4d28s21
     $       bc(issscpy),nbasisc(isw),                                  4d28s21
     d' makeguess. 18')
         ioffg=1                                                         10d23s17
         ngot=0
         do isg=1,nsymbg                                                 10d23s17
          if(ldebug)write(6,*)('isg = '),isg,morbp(isg)                 7d5s23
          if(morbp(isg).gt.0)then                                        5d14s19
           do i=1,nbasisc(isw)                                          4d27s21
            do j=1,morbp(isg)*ncomp                                     10d7s22
             iad1=isi+i-1+nbasisc(isw)*(j-1)                            4d27s21
             iad2=iovs(isg,isw)+j-1+morbp(isg)*ncomp*(i-1)              10d7s22
             bc(iad1)=bc(iad2)                                            4d25s18
            end do                                                        4d25s18
           end do                                                         4d25s18
           call dgemm('n','n',nbasisc(isw),morbp(isg)*ncomp,            10d7s22
     $          nbasisc(isw),1d0,                                       10d7s22
     $        bc(ixi),nbasisc(isw),bc(isi),nbasisc(isw),0d0,             4d27s21
     $        bc(itemp),nbasisc(isw),                                   4d27s21
     d' makeguess. 19')
           if(ldebug)then
            write(6,*)('orbs to usea '),ioffg
            call prntm2(vgues(ioffg),morbp(isg)*ncomp,nbasdwsg(isg),     10d7s22
     $        morbp(isg)*ncomp)                                         10d7s22
           end if                                                       7d5s23
           itmp=ibcoff                                                  4d27s21
           ibcoff=itmp+nbasisc(isw)*nbasdwsg(isg)                       4d27s21
           call enough('makeguess. 18',bc,ibc)
           call dgemm('n','n',nbasisc(isw),nbasdwsg(isg),               10d7s22
     $          morbp(isg)*ncomp,1d0,bc(itemp),nbasisc(isw),            10d7s22
     $          vgues(ioffg),morbp(isg)*ncomp,0d0,                      10d7s22
     $          bc(itmp),nbasisc(isw),                                  4d27s21
     d' makeguess. 20')
           if(ldebug)then                                               7d5s23
            call dgemm('t','n',nbasdwsg(isg),nbasdwsg(isg),nbasisc(isw), 4d27s21
     $          1d0,bc(itmp),nbasisc(isw),bc(itmp),nbasisc(isw),0d0,    4d27s21
     $          bc(ibcoff),nbasdwsg(isg),
     d' makeguess. 21')
            write(6,*)('orbs in ob basis: ')
            call prntm2(bc(itmp),nbasisc(isw),nbasdwsg(isg),
     $           nbasisc(isw))
            write(6,*)('ortho test ')
            call prntm2(bc(ibcoff),nbasdwsg(isg),nbasdwsg(isg),
     $          nbasdwsg(isg))
           end if                                                       7d5s23
           itmpv=ibcoff                                                 4d27s21
           ibcoff=itmpv+nbasisc(isw)*nbasdwsg(isg)                      4d27s21
           call enough('makeguess. 19',bc,ibc)
           call dgemm('n','n',nbasisc(isw),nbasdwsg(isg),nbasisc(isw),  4d28s21
     $          1d0,bc(ix),nbasisc(isw),bc(itmp),nbasisc(isw),0d0,      4d28s21
     $          bc(itmpv),nbasisc(isw),                                 4d28s21
     d' makeguess. 22')
           if(ldebug)then
            write(6,*)('Vw')
            call prntm2(bc(itmpv),nbasisc(isw),nbasdwsg(isg),            4d28s21
     $          nbasisc(isw))                                           4d28s21
            write(6,*)('overlap for symmetries '),isg,isw
            call prntm2(bc(iovs(isg,isw)),morbp(isg)*ncomp,nbasisc(isw),   10d7s22
     $        morbp(isg)*ncomp)                                         10d11s22
            write(6,*)('multiply from the left by transpose of '),         10d23s17
     $        ('orbs in ao basis ')
            write(6,*)('orbs to useb '),isg
            call prntm2(vgues(ioffg),morbp(isg)*ncomp,nbasdwsg(isg),       10d7s22
     $        morbp(isg)*ncomp)                                         10d7s22
           end if
           itmp=ibcoff                                                    10d23s17
           ibcoff=itmp+nbasdwsg(isg)*nbasisc(isw)                       4d27s21
           call enough('makeguess. 20',bc,ibc)
           call dgemm('t','n',nbasdwsg(isg),nbasisc(isw),
     $          morbp(isg)*ncomp,1d0,vgues(ioffg),morbp(isg)*ncomp,     10d7s22
     $          bc(iovs(isg,isw)),morbp(isg)*ncomp,0d0,bc(itmp),        10d7s22
     $          nbasdwsg(isg),                                          10d7s22
     d' makeguess. 23')
           if(ldebug)then                                               7d5s23
            write(6,*)('result ')
            call prntm2(bc(itmp),nbasdwsg(isg),nbasisc(isw),            7d5s23
     $           nbasdwsg(isg))                                         7d5s23
            write(6,*)('multiply from the right by orthogonalization '),   10d23s17
     $        ('matrix ')
            call prntm2(bc(ixi),nbasisc(isw),nbasisc(isw),nbasisc(isw))    4d27s21
           end if
           itmp2=ibcoff                                                   10d23s17
           ibcoff=itmp2+nbasdwsg(isg)*nbasisc(isw)                      4d27s21
           call enough('makeguess. 21',bc,ibc)
           call dgemm('n','n',nbasdwsg(isg),nbasisc(isw),nbasisc(isw),  4d27s21
     $         1d0,bc(itmp),nbasdwsg(isg),bc(ixi),nbasisc(isw),0d0,     4d27s21
     $         bc(itmp2),nbasdwsg(isg),                                 5d14s19
     d' makeguess. 24')
           itmp3=ibcoff                                                   10d26s17
           ibcoff=itmp3+nbasdwsg(isg)*nbasdwsg(isg)                       10d26s17
           call enough('makeguess. 22',bc,ibc)
           call dgemm('n','t',nbasdwsg(isg),nbasdwsg(isg),nbasisc(isw),
     $        1d0,bc(itmp2),nbasdwsg(isg),bc(itmp2),nbasdwsg(isg),0d0,
     $        bc(itmp3),nbasdwsg(isg),                                  10d26s17
     d' makeguess. 25')
           if(ldebug)then                                               7d5s23
           write(6,*)('ortho test ')
           call prntm2(bc(itmp3),nbasdwsg(isg),nbasdwsg(isg),             10d26s17
     $        nbasdwsg(isg))                                            10d26s17
           end if                                                       7d5s23
           nkeep=0
           ibct=ibcoff                                                    10d26s17
           ibcoff=ibct+nbasdwsg(isg)                                      10d26s17
           call enough('makeguess. 23',bc,ibc)
           do i=0,nbasdwsg(isg)-1
            iad=itmp3+i*(nbasdwsg(isg)+1)                                 10d26s17
            if(bc(iad).gt.0.5d0)then                                      10d26s17
             ito=itmp+nkeep*nbasisc(isw)                                4d28s21
             ifrm=itmp2+i                                                 10d26s17
             do j=0,nbasisc(isw)-1                                      4d27s21
              bc(ito+j)=bc(ifrm+j*nbasdwsg(isg))                          10d26s17
             end do                                                       10d26s17
             ibc(ibct+nkeep)=ibc(ibcode(isg)+i)                           10d26s17
             nkeep=nkeep+1
            end if                                                        10d26s17
           end do
           if(nkeep.gt.0)then
            nvguess=ncomp                                               10d11s22
            ircode=ircodecx                                             1d4s23
            if(ldebug)then                                              7d5s23
             write(6,*)('ircode c'),ircode                              1d24s23
             write(6,*)('modified schmidt orthogonalization of '),nkeep,
     $          (' vectors ')
             call prntm2(bc(itmp),nbasisc(isw),nkeep,nbasisc(isw))
            end if                                                      7d5s23
            do i=0,nkeep-1                                                10d26s17
             ip=i+1
             iadi=itmp+nbasisc(isw)*i                                   4d27s21
             do j=0,i-1
              iadj=itmp+nbasisc(isw)*j                                  4d27s21
              dot=0d0
              do k=0,nbasisc(isw)-1                                     4d27s21
               dot=dot+bc(iadi+k)*bc(iadj+k)
              end do
              do k=0,nbasisc(isw)-1                                     4d27s21
               bc(iadi+k)=bc(iadi+k)-dot*bc(iadj+k)
              end do
             end do
             sum=0d0
             do k=0,nbasisc(isw)-1                                      4d27s21
              sum=sum+bc(iadi+k)**2
             end do
             ibcc=ibc(ibct+i)                                             10d26s17
             ito=ibc(ibcxs+ibcc)+nbasp*ibc(ibcxn+ibcc)                    10d24s17
             ibc(ibcxn+ibcc)=ibc(ibcxn+ibcc)+1                            10d24s17
             xnorm=1d0/sqrt(sum)                                          10d24s17
             do k=0,nbasisc(isw)-1                                      4d27s21
              bc(iadi+k)=bc(iadi+k)*xnorm                                 10d24s17
              bc(ito+k)=bc(iadi+k)                                        10d24s17
             end do
             if(ldebug)then                                             7d5s23
              write(6,*)('store to address '),ito-itmp
              call prntm2(bc(ito),nbasisc(isw),1,nbasisc(isw))           10d7s22
             end if
             bc(ito+nbasisc(isw))=dfloat(ip*10+isg)                     4d27s21
            end do
           end if
           ngot=ngot+nkeep
           ioffg=ioffg+morbp(isg)*nbasdwsg(isg)*ncomp                   10d11s22
          else                                                          5d14s19
           if(isg.ne.nsymbg)then                                        5d14s19
            ipao(isg+1)=jpao                                            5d14s19
           end if                                                       5d14s19
          end if                                                         5d14s19
         end do                                                          10d23s17
         ivecsav0=ivecsav+1
         jptno=iptno(isw)                                                10d24s17
         jpts=ipts(isw)                                                  10d24s17
 3356    format('      type    #   sym equivalent')                      4d24s18
         irunt=0                                                         4d24s18
         nhitv=0
         if(ldebug)write(6,*)('vecr no. 4')                             1d24s23
         do i=1,nbcx                                                     10d24s17
          iad=ibc(ibcxs+i)-1                                             10d24s17
          if(ldebug)then
           write(6,*)('for nbcx index '),i,iad
           write(6,*)i,iad
           write(6,*)('tmp? ')
           nn=ibc(ibcxn+i)
           call prntm2(bc(iad+1),nbasisc(isw),nn,nbasp)
          end if
          do j=1,ibc(ibcxn+i)                                            10d24s17
           do k=1,nbasisc(isw)                                          7d5s23
            vecr(ivecsav+k)=bc(iad+k)                                    10d24s17
           end do                                                        10d24s17
           iad=iad+nbasp                                                 10d24s17
           val=bc(iad)*0.1d0
           ival=int(val)
           ibc(jptno)=ival
           ibc(jpts)=int(bc(iad)-dfloat(ival*10))
           irunt=irunt+1                                                 4d24s18
           nhitv=nhitv+1
 3357      format(5i5)
           jptno=jptno+1                                                 10d24s17
           jpts=jpts+1                                                   10d24s17
           ivecsav=ivecsav+nbasisc(isw)                                  5d9s19
          end do                                                         10d24s17
         end do                                                          10d24s17
         if(ldebug)then
          write(6,*)('isw = '),isw
          write(6,*)('vecr so far: '),ivecsav0,ivecsav
          call prntm2(vecr(ivecsav0),nbasisc(isw),irunt,nbasisc(isw))    7d5s23
         end if
         iad=ibc(ibcxs)-1                                                10d24s17
         do j=1,ibc(ibcxn)                                               10d24s17
          do k=1,nbasisc(isw)                                           4d28s21
           vecr(ivecsav+k)=bc(iad+k)                                     10d24s17
          end do                                                         10d24s17
          ivecsav=ivecsav+nbasisc(isw)                                  4d28s21
          iad=iad+nbasp                                                  10d24s17
          val=bc(iad)*0.1d0
          ival=int(val)
          ibc(jptno)=ival
          ibc(jpts)=int(bc(iad)-dfloat(ival*10))
          irunt=irunt+1                                                 4d24s18
          jptno=jptno+1                                                 10d24s17
          jpts=jpts+1                                                   10d24s17
         end do                                                          10d24s17
         if(iocode.eq.1)then                                            6d4s19
          nbasdws(isw)=irunt                                              5d14s19
         end if                                                         6d4s19
         if(idorel.ne.0)then                                            10d7s22
          ntop=nbasisc(isw)/2                                           10d7s22
         else                                                           10d7s22
          ntop=nbasisc(isw)                                             10d7s22
         end if                                                         10d7s22
         if(ldebug)then                                                 7d5s23
          write(6,*)('vectors so far '),ivecsav0,ivecsav                        7d5s23
          call prntm2(vecr(ivecsav0),nbasisc(isw),irunt,nbasisc(isw))   7d5s23
          write(6,*)('vecr no. 5')                                      7d5s23
         end if                                                         7d5s23
         do ipad=irunt+1,ntop                                           10d7s22
          do k=1,nbasisc(isw)                                           4d28s21
           vecr(ivecsav+k)=0d0
          end do
          vecr(ivecsav+ipad)=1d0
          do i=0,ipad-2                                                 10d11s22
           ii=ivecsav0-1+nbasisc(isw)*i                                 10d7s22
           dot=0d0                                                      10d7s22
           do k=1,nbasisc(isw)                                          10d7s22
            dot=dot+vecr(ivecsav+k)*vecr(ii+k)                          10d7s22
           end do                                                       10d7s22
           do k=1,nbasisc(isw)                                          10d7s22
            vecr(ivecsav+k)=vecr(ivecsav+k)-dot*vecr(ii+k)              10d7s22
           end do                                                       10d7s22
          end do                                                        10d7s22
          dot=0d0                                                       10d7s22
          do k=1,nbasisc(isw)                                           10d7s22
           dot=dot+vecr(ivecsav+k)**2                                   10d7s22
          end do                                                        10d7s22
          dot=1d0/sqrt(dot)                                             10d7s22
          do k=1,nbasisc(isw)                                           10d7s22
           vecr(ivecsav+k)=vecr(ivecsav+k)*dot                          10d7s22
          end do                                                        10d7s22
          ivecsav=ivecsav+nbasisc(isw)                                  4d28s21
         end do
         itest=ivecsav-ivecsav0
         if(ldebug)then
          write(6,*)('my guess vectors in ob: '),ivecsav0,itest,
     $        loc(vecr(ivecsav0))
          write(6,*)('symmetry block '),isw
          write(6,*)('myguessi = '),myguessi
          call prntm2(vecr(ivecsav0),nbasisc(isw),nbasisc(isw),           4d28s21
     $       nbasisc(isw))                                              4d28s21
         end if                                                         7d5s23
        if(myguessi.ne.2)then                                           4d28s21
c
c     transform to ao basis
c
         itmpx=ibcoff                                                    10d24s17
         ibcoff=itmpx+nbasisc(isw)*nbasisc(isw)                         4d28s21
         call enough('makeguess. 24',bc,ibc)
         if(ldebug)then                                                 3d19s24
          write(6,*)('ixi to use '),ixi
          call prntm2(bc(ixi),nbasisc(isw),nbasisc(isw),nbasisc(isw))
         end if                                                         3d19s24
         call dgemm('n','n',nbasisc(isw),nbasisc(isw),nbasisc(isw),1d0, 4d28s21
     $       bc(ixi),nbasisc(isw),vecr(ivecsav0),nbasisc(isw),0d0,      4d28s21
     $       bc(itmpx),nbasisc(isw),                                    4d28s21
     d' makeguess. 26')
         if(ldebug)write(6,*)('vecr no. 6')                             1d24s23
         do i=0,nbasisc(isw)*nbasisc(isw)-1                             4d28s21
          vecr(ivecsav0+i)=bc(itmpx+i)                                   10d24s17
         end do                                                          10d24s17
         if(ldebug)then                                                 1d24s23
          write(6,*)('my guess vectors in ao basis: '),ivecsav0,
     $        loc(vecr(ivecsav0))
          call prntm2(vecr(ivecsav0),nbasisp(isw)*ncomp,nbasdws(isw),    10d11s22
     $        nbasisp(isw)*ncomp)                                       10d11s22
         end if                                                         1d24s23
         end if                                                         4d28s21
         ibcoff=ixi                                                      10d24s17
        end if                                                          5d14s19
       end do                                                           10d23s17
       ierr=0
       if(ldebug)write(6,*)('return one ')
       return                                                           10d24s17
      end if
      if(ldebug)write(6,*)('idorel, idorelg: '),idorel,idorelg          5d3s21
      if(idorel.eq.idorelg.or.(idorel.ne.0.and.idorelg.ne.0))then       10d25s20
       ioff=0
       if(ldebug)write(6,*)('vecr no. 7'),myguessi
       do isb=1,nsymb
        if(myguessi.eq.2)then
         nh=nbasdwsg(isb)*nbasisc(isb)                                  2d1s23
        else if(myguessi.eq.1)then                                      2d15s19
         nh=nbasdwsg(isb)*nbasisp(isb)                                   2d15s19
         if(idorel.ne.0)nh=nh*2                                         2d15s19
        end if
        do i=1,nh                                                       2d15s19
         vecr(i+ioff)=vgues(i+ioff)
        end do
        ioff=ioff+nh                                                    2d15s19
       end do
      else
       ioffnr=0                                                         2d18s16
       ioffr=0                                                          2d18s16
       do isb=1,nsymb                                                   2d18s16
        if(ldebug)then                                                  5d3s21
        end if
        do j=1,nbasdwsg(isb)                                            2d18s16
         do i=1,nbasisp(isb)                                            3d14s20
          vecr(i+ioffr)=vgues(i+ioffnr)                                 2d18s16
         end do                                                          2d18s16
         ioffr=ioffr+nbasisp(isb)                                       3d14s20
         if(iabs(myguessi).eq.1)then                                    3d14s16
          do i=1,nbasisp(isb)                                           3d14s20
           vecr(i+ioffr)=vgues(i+ioffnr)                                 2d18s16
          end do                                                          2d18s16
         else
          do i=1,nbasdwsg(isb)                                            2d18s16
           vecr(i+ioffr)=0d0
          end do                                                          2d18s16
         end if
         ioffr=ioffr+nbasisp(isb)                                       3d14s20
         ioffnr=ioffnr+nbasisp(isb)                                     3d14s20
        end do                                                          2d18s16
       end do                                                           2d18s16
      end if
c
c     orthogonalize. skip if myguessi=-1. if myguessi=1 need overlap.   3d14s16
c
      if(ldebug)                                                        1d24s23
     $     write(6,*)('myguessi b4 orthogonalization test '),myguessi   1d24s23
      if(myguessi.ge.0)then                                             3d14s16
       ioffs=0                                                           2d18s16
       ioffo=0                                                           2d18s16
       if(ldebug)write(6,*)('orthogonalizing '),myguessi                1d24s23
       ircode=4                                                         10d12s22
       if(ldebug)write(6,*)('ircode d')                                 1d24s23
       if(ldebug)write(6,*)('icanog '),(icanog(isb),isb=1,nsymb)
       do isb=1,nsymb                                                    2d18s16
        nbasdws(isb)=nbasdwsg(isb)                                      2d15s19
        nrow=nbasisc(isb)-icanog(isb)                                   5d3s23
        nrowupdate=nbasisc(isb)                                         5d4s23
        if(myguessi.eq.1)then                                           2d15s19
         nrow=nbasisp(isb)                                              2d15s19
         if(idorel.ne.0)nrow=nrow*2                                     2d15s19
         nrowupdate=nrow                                                5d4s23
        end if                                                          2d15s19
        if(ldebug)then                                                  5d3s21
        write(6,*)('for symmetry block '),isb,ioffo
        call prntm2(vecr(ioffo+1),nrow,nbasdws(isb),nrow)
        call dgemm('t','n',nbasdws(isb),nbasdws(isb),nrow,1d0,
     $       vecr(ioffo+1),nrow,vecr(ioffo+1),nrow,0d0,
     $       bc(ibcoff),nbasdws(isb))
        call prntm2(bc(ibcoff),nbasdws(isb),nbasdws(isb),nbasdws(isb))
        end if
        if(myguessi.eq.1)then                                            2d18s16
         if(ldebug)write(6,*)('vecr no. 8')                             1d24s23
         nn=nrow**2                                                     2d15s19
         iovr=ibcoff                                                     2d18s16
         ibcoff=iovr+nn                                                 2d15s19
         itmp=ibcoff                                                     2d18s16
         ibcoff=itmp+nn                                                 2d15s19
         call enough('makeguess. 25',bc,ibc)
         jovr=iovr-1                                                     2d18s16
          mm=(nrow*(nrow+1))/2                                          2d15s19
          if(iocode.eq.1)then                                           5d7s21
           do i=1,mm                                                     2d15s19
            bc(jovr+i)=ovr(i+ioffs)                                      2d15s19
           end do                                                        2d15s19
           if(ldebug)then                                                5d3s21
            write(6,*)('overlap '),ioffs
            call mpprnt2(bc(iovr),nrow)                                    2d15s19
           end if                                                         5d2s16
           call square(bc(iovr),nrow)                                     2d15s19
          else                                                          5d7s21
           do i=1,nn                                                    5d7s21
            bc(jovr+i)=scopy(i+ioffs)                                   5d7s21
           end do                                                       5d7s21
          end if                                                        5d7s21
         ioffs=ioffs+nn                                                 2d15s19
         if(ldebug)then                                                 5d3s21
         write(6,*)('squared ')
         call prntm2(bc(iovr),nrow,nrow,nrow)                           2d15s19
         end if
         jtmp=itmp-1
         if(nbasdws(isb).gt.0)then                                      5d14s19
         call dgemm('n','n',nrow,nbasdws(isb),nrow,                     2d15s19
     $        1d0,bc(iovr),nrow,vecr(ioffo+1),nrow,                     2d15s19
     $        0d0,bc(itmp),nrow,                                        2d15s19
     d' makeguess. 27')
         if(ldebug)then                                                 5d3s21
         write(6,*)('s*v ')
         call prntm2(bc(itmp),nrow,nbasdws(isb),nrow)
         end if
         itmp2=ibcoff
         ibcoff=itmp2+nbasdws(isb)*nrow                                 2d15s19
         call enough('makeguess. 26',bc,ibc)
         do i=0,nbasdws(isb)-1
          do j=0,nrow-1                                                 2d15s19
           ji=itmp+j+nrow*i                                             2d15s19
           ii=itmp2+i+nbasdws(isb)*j
           bc(ii)=bc(ji)
          end do
         end do
         call dgemm('n','n',nbasdws(isb),nbasdws(isb),nrow,             2d15s19
     $        1d0,bc(itmp2),nbasdws(isb),vecr(ioffo+1),nrow,            2d15s19
     $        0d0,bc(itmp),nbasdws(isb),
     d' makeguess. 28')
         end if                                                         5d14s19
         if(ldebug)then                                                 5d3s21
         write(6,*)('vt*s*v ')
         call prntm2(bc(itmp),nbasdws(isb),nbasdws(isb),nbasdws(isb))
         end if
        end if
        do i=1,nbasdws(isb)                                              2d18s16
         iadi=ioffo+nrow*(i-1)                                          2d15s19
         iadip=iadi+1
         sz=0d0
         do j=1,nrow                                                    2d15s19
          sz=sz+vecr(iadi+j)**2                                          2d18s16
         end do                                                          2d18s16
         if(sz.lt.1d-10)vecr(iadi+i)=1d0                                 2d18s16
         do j=1,i-1                                                      2d18s16
          iad=ioffo+nrow*(j-1)                                          2d15s19
          iadp=iad+1
          if(myguessi.eq.1)then
           call dgemm('n','n',nrow,1,nrow,1d0,bc(iovr),                 2d15s19
     $          nrow,vecr(iadp),nrow,0d0,bc(itmp),nrow,                 2d15s19
     d' makeguess. 29')
           dot=0d0
           do k=1,nrow                                                  2d15s19
            dot=dot+vecr(iadi+k)*bc(jtmp+k)
           end do
          else
           dot=0d0
           do k=1,nrow                                                  2d1s23
            dot=dot+vecr(iadi+k)*vecr(iad+k)
           end do
          end if
          do k=1,nrow                                                   2d15s19
           vecr(iadi+k)=vecr(iadi+k)-dot*vecr(iad+k)
          end do
         end do                                                          2d18s16
         if(myguessi.eq.1)then
          call dgemm('n','n',nrow,1,nrow,1d0,bc(iovr),                  2d15s19
     $         nrow,vecr(iadip),nrow,0d0,bc(itmp),nrow,                 2d15s19
     d' makeguess. 30')
          dot=0d0
          do k=1,nrow                                                   2d15s19
           dot=dot+vecr(iadi+k)*bc(jtmp+k)
          end do
         else
          dot=0d0
          do k=1,nrow                                                   2d1s23
           dot=dot+vecr(iadi+k)*vecr(iadi+k)
          end do
         end if
         if(dot.ne.dot)then
          write(6,*)('got not a number when orthogonalizing ... ')
          write(6,*)('vector '),i,('from symmetry block '),isb
          ierr=1
          return
         end if
         xnorm=1d0/sqrt(dot)                                             2d18s16
         do k=1,nrow                                                    2d15s19
          vecr(iadi+k)=vecr(iadi+k)*xnorm
         end do
        end do                                                           2d18s16
        if(ldebug)then                                                  5d3s21
        write(6,*)('orthogonalized vectors ')
        call prntm2(vecr(ioffo+1),nrow,nbasdws(isb),nrow)               2d15s19
        call dgemm('t','n',nbasdws(isb),nbasdws(isb),nrow,1d0,
     $       vecr(ioffo+1),nrow,vecr(ioffo+1),nrow,0d0,
     $       bc(ibcoff),nbasdws(isb))
        call prntm2(bc(ibcoff),nbasdws(isb),nbasdws(isb),nbasdws(isb))
        end if
        ioffo=ioffo+nrowupdate*nbasdws(isb)                             5d4s23
        if(myguessi.eq.1)ibcoff=iovr                                     2d18s16
       end do
      else                                                              3d14s16
       write(6,*)('not orthogonalizing input vectors in ao basis !!!!') 3d14s16
      end if                                                            3d14s16
      idum=1
      return
      end
