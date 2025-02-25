c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine opto(noc,iyc,ixyc,nv4,ibig2,ibmat10,numpro,nowpro,     3d12s12
     $     igmat,times,lprint,maxd,iorbn,iumat,iamat10,itmat,mgmat,
     $     iorbx,ih0mo,tovr,tchange,i404,ipk,jmats,iqk,kmats,           8d5s14
     $     iden1,igdig,iaseward,thrs,thrsa,xlam,idwsdebx,idavopt,       3d3s20
     $     nbasdwsc,nvirtc,eiglastgood,ncomp,ihessa,ihessb,ihessc,      12d4s17
     $     ihessd,icas,idoub,iacto,iamatx,smallest,nlzz,iorbsym,mlx,    4d5s21
     $     inbr,ipbr,bc,ibc)                                            11d9s22
      implicit real*8 (a-h,o-z)                                         4d2s18
      external second
      integer mgmat,iarg1,iarg2,iarg3,iarg4,noc8                        3d12s12
      logical lprint,lprintx                                            7d13s05
c
c mpec2.1 version zeta copyright u.s. government
c     subroutine opto: perform orbital optimization step.
c
      include "common.store"
      include "common.hf"
      include "common.print"                                            1d14s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension noc(8),iyc(1),ixyc(1),nv4(1),ibig2(1),igmat(8,8),       3d12s12
     $     times(26,3),iorbn(1),iumat(1),itmat(1),mgmat(8,8),iorbx(1),
     $     ipk(1),jmats(1),iqk(1),kmats(1),iden1(1),igdig(1),nbasdwsc(8)8d31s15
     $     ,nvirtc(8),idoub(8),iacto(8),ihessa(8,8),ihessb(8,8),        12d4s17
     $     ihessc(8,8),ihessd(8),iamatx(8),ixycu(8),iycb(8),iorbsym(*)  4d4s21
      save                                                              2d20s14
      data timetot/0d0/                                                 4d23s18
      data iahu,ihess,irelc/0,0,0/                                                      2d20s14
      data icall/0/
      if(iprtr(17).eq.0)then                                            3d3s20
       idwsdeb=0                                                        3d3s20
       lprintx=.false.                                                   3d3s20
      else                                                              3d3s20
       idwsdeb=1000
       lprintx=.true.                                                    3d3s20
       write(6,*)('Hi, My name is opto')
      end if                                                            3d3s20
      ibcoffo=ibcoff                                                    3d3s20
      do isb=1,nsymb                                                    3d3s20
       iycb(isb)=ibcoff                                                 3d3s20
       ibcoff=iycb(isb)+nv4(isb)                                        3d3s20
      end do                                                            3d3s20
      call enough('opto.  1',bc,ibc)
      do i=iycb(1),ibcoff-1                                             3d3s20
       bc(i)=0d0                                                        3d3s20
      end do                                                            3d3s20
      bc(iycb(1))=1d0                                                   3d3s20
      i404=0
      lprintx=.false.
  350 continue
      thc=0d0                                                           4d23s18
      call second(timea)                                                4d23s18
  210 continue
      if(idwsdeb.gt.10)write(6,*)('xlam = '),xlam
      xlmi2=1d0/xlam
      do isb=1,nsymb
       if(noc(isb).gt.0)then                                            9d14s04
        jyc=iyc(isb)-1                                                  9d14s04
        jxyc=ixyc(isb)-1                                                9d14s04
        do 201 i=1,nv4(isb)                                             9d14s04
         bc(jyc+i)=0d0
         bc(jxyc+i)=0d0
  201   continue
        if(isb.eq.1.or.icas.eq.0)then                                   11d29s17
         bc(iyc(isb))=1d0
         bc(ixyc(isb))=1d0                                              11d29s17
        end if                                                          11d29s17
       end if                                                           9d14s04
      end do                                                            9d14s04
      iterd=0
      idav=1
      eigo=1d10
      if(idwsdeb.gt.10)write(6,*)('before davidson, ibcoff = '),ibcoff
  202 continue
      iterd=iterd+1
      call second(time0)                                                3d7s07
c
c      form this processors hess*yc
c
      sumb=0d0                                                          9d16s04
      do isb=1,nsymb
       if(noc(isb).gt.0)then
        ii=ibig2(isb)+nv4(isb)*(idav-1)-1
        do idws=1,nv4(isb)
         bc(ii+idws)=0d0
        end do
       end if
      end do
      ibmat1=ibmat10                                                    9d9s04
      if(icas.ne.0)then                                                 12d28s17
       do isb=1,nsymb                                                   12d28s17
        ixycu(isb)=ibig2(isb)+nv4(isb)*(idav-1)                         12d28s17
       end do                                                           12d28s17
       call second(time1)
       call hessyccas(iyc,ixycu,iamatx,ihessa,ihessb,ihessc,idoub,iacto,12d28s17
     $     noc,nvirtc,nsymb,xlmi2,bc,ibc)                               11d10s22
       do isb=1,nsymb                                                   12d28s17
        call dws_gsumf(bc(ixycu(isb)),nv4(isb))                         12d28s17
       end do                                                           12d28s17
       call second(time2)
       telap=time2-time1-tovr
       thc=thc+telap
      else                                                              12d28s17
      do isb=1,nsymb                                                    9d14s04
       ibmat2=ibmat1+noc(isb)*noc(isb)                                  9d15s04
         noc8=noc(isb)                                                          9d9s04
         if(idwsdeb.gt.10)then
         write(6,*)('bmat1')
         call prntm2(bc(ibmat1),noc8,noc8,noc8)
         write(6,*)('bmat2')                                            7d23s07
         iarg2=nvirtc(isb)                                              8d31s15
         call prntm2(bc(ibmat2),iarg2,noc8,iarg2)                       7d23s07
         end if
       do isa=1,nsymb                                                   9d15s04
        if(noc(isb).gt.0.and.noc(isa).gt.0)then                         9d15s04
         if(idwsdeb.gt.10)then
         writE(6,*)('for symmetry block '),isa,isb                       9d15s04
         end if
         noca=noc(isa)                                                   9d14s04
         nvirta=nvirtc(isa)                                             8d31s15
         nocb=noc(isb)                                                   9d14s04
         nvirtb=nvirtc(isb)                                             8d31s15
         if(isb.ge.isa)then                                             3d12s12
          call ilimts(nvirtc(isb),nvirtc(isa),numpro,nowpro,ilj,ihj,    8d31s15
     $               i1s,i1e,i2s,i2e)                                   9d15s04
         else                                                           3d12s12
          call ilimts(nvirtc(isa),nvirtc(isb),numpro,nowpro,ilj,ihj,    8d31s15
     $               i1s,i1e,i2s,i2e)                                   9d15s04
         end if                                                         3d12s12
         if(idwsdeb.gt.10)then
         writE(6,*)('ibcoff '),ibcoff
         write(6,*)('calling hessycmo '),nowpro
         write(6,*)('whats going into hessycmo ')
         iarg1=1
         iarg2=nv4(isa)
         call prntm2(bc(iyc(isa)),iarg1,iarg2,iarg1)
         end if
         itmp=ibcoff                                                    9d15s04
         ibcoff=itmp+nv4(isb)                                           9d15s04
         call enough('opto.  2',bc,ibc)
         if(isa.eq.isb)then
          nvuse=nv4(isb)                                                9d15s04
         else                                                           9d15s04
          nvuse=0                                                       9d15s04
         end if                                                         9d15s04
         call second(time1)                                             10d19s04
         ihess=ihess+1                                                  2d21s14
         if(isb.ge.isa)then                                             3d12s12
          call hessycmo(nvuse,1,itmp,iyc(isa),noc(isa),noc(isb),        3d12s12
     $                 nvirtc(isa),nvirtc(isb),                         8d31s15
     $                ibmat1,ibmat2,igmat(isa,isb),i2s,i2e,i1s,i1e,     7d26s07
     $                xlmi2,1,bc,ibc)                                   11d10s22
         else                                                           3d12s12
          call hessycmo(nvuse,1,itmp,iyc(isa),noc(isa),noc(isb),        3d12s12
     $                  nvirtc(isa),nvirtc(isb),                        8d31s15
     $                ibmat1,ibmat2,igmat(isb,isa),i2s,i2e,i1s,i1e,     2d5s14
     $                xlmi2,-1,bc,ibc)                                  11d10s22
         end if                                                         3d12s12
         call second(time2)                                             10d19s04
         telap=time2-time1-tovr                                         10d19s04
         times(9,1)=times(9,1)+telap                                    10d19s04
         times(9,2)=times(9,2)+telap**2                                 10d19s04
         times(9,3)=times(9,3)+1d0                                      10d19s04
         if(idwsdeb.gt.10)then
         iarg1=1
         iarg2=nv4(isb)
         end if
         ii=ibig2(isb)+nv4(isb)*(idav-1)-1
         jj=itmp-1
          do idws=1,nv4(isb)
           bc(ii+idws)=bc(ii+idws)+bc(jj+idws)
          end do
         ibcoff=itmp                                                    9d15s04
         if(idwsdeb.gt.10)then
         writE(6,*)('so far ')
         call prntm2(bc(ibig2(isb)+nv4(isb)*(idav-1)),iarg1,iarg2,iarg1)
         end if
         call second(time2)
        end if
       end do                                                           9d15s04
       ibmat1=ibmat1+noc(isb)*nbasdwsc(isb)                             8d31s15
      end do                                                            9d14s04
      end if                                                            12d28s17
c
      if(idav.eq.1)then
       eig=0d0
       sumb=0d0
       do isb=1,nsymb                                                   9d14s04
        if(noc(isb).gt.0)then                                           9d14s04
         jbig=ibig2(isb)-1                                              9d14s04
         jyc=iyc(isb)-1                                                 9d14s04
         do 211 i=1,nv4(isb)                                            9d14s04
          eig=eig+bc(jbig+i)*bc(jyc+i)
  211    continue
         if(icas.eq.0.or.isb.eq.1)then                                  12d28s17
          sumb=sumb+bc(ibig2(isb))
         end if                                                         12d28s17
         if(idwsdeb.gt.10)then
          write(6,*)('after global suma '),isb
         iarg1=nv4(isb)
         iarg2=1                                                        9d14s04
         call prntm2(bc(ibig2(isb)),iarg1,iarg2,iarg1)                  9d14s04
         writE(6,*)('sumb '),sumb
         if(nsymb.eq.1)then
          call printa(bc(ibig2(isb)+1),iacto,idoub,1,0,idoub,0,1,0,
     $         bc(ibcoff))
          call printa(bc(ibig2(isb)+1+iacto(1)*idoub(1)),
     $         nvirtc,noc,1,0,noc,0,1,0,
     $         bc(ibcoff))
         else
          iuse=ibig2(isb)
          if(isb.eq.1)iuse=iuse+1
          call prntm2(bc(iuse),iacto(isb),idoub(isb),iacto(isb))
          call prntm2(bc(iuse+iacto(isb)*idoub(isb)),
     $         nvirtc(isb),noc(isb),nvirtc(isb))
         end if
         end if
        end if                                                          9d14s04
       end do                                                           9d14s04
       if(idwsdeb.gt.10)write(6,*)('eig '),eig
       idavp=idav
      else
       call second(time1)
       itmp1=ibcoff
       ibcoff=itmp1+idav*idav
       do 212 i=1,idav
        jtmp1=itmp1-1+(i-1)*idav
        do 213 j=1,i
         bc(jtmp1+j)=0d0
  213   continue
  212  continue
       sum1=0d0
       do isb=1,nsymb                                                   9d14s04
        if(noc(isb).gt.0)then                                           9d14s04
         jbig=ibig2(isb)+(idav-1)*nv4(isb)                              9d14s04
         if(isb.eq.1.or.icas.eq.0)then                                  12d28s17
          sum1=sum1+bc(jbig)
         end if                                                         12d28s17
         if(idwsdeb.gt.10)then
         write(6,*)('after global sumb '),isb                            9d14s04
         iarg1=nv4(isb)                                                 9d14s04
         iarg2=1                                                        9d14s04
         call prntm2(bc(jbig),iarg1,iarg2,iarg1)                        9d14s04
         if(nsymb.eq.1)then
          call printa(bc(jbig+1),iacto,idoub,1,0,idoub,0,1,0,
     $         bc(ibcoff))
          call printa(bc(jbig+1+iacto(1)*idoub(1)),
     $         nvirtc,noc,1,0,noc,0,1,0,
     $         bc(ibcoff))
         else
          iuse=jbig
          if(isb.eq.1)iuse=iuse+1
          call prntm2(bc(iuse),iacto(isb),idoub(isb),iacto(isb))
          call prntm2(bc(iuse+iacto(isb)*idoub(isb)),
     $         nvirtc(isb),noc(isb),nvirtc(isb))
         end if
         write(6,*)('sum of first element so far '),sum1
         end if
        end if
       end do
       do 214 i=1,idav
        do 215 j=1,i
         sum=0d0
         do isb=1,nsymb
          if(noc(isb).gt.0)then
           jxyc=ixyc(isb)-1+(i-1)*nv4(isb)                               9d14s04
           jbig=ibig2(isb)-1+(j-1)*nv4(isb)                             9d14s04
           do 216 k=1,nv4(isb)
            sum=sum+bc(jbig+k)*bc(jxyc+k)
  216      continue
          end if
         end do
         bc(itmp1+j-1+(i-1)*idav)=sum
         bc(itmp1+i-1+(j-1)*idav)=sum
  215   continue
  214  continue
       call second(time2)
       telap=time2-time1-tovr
       lwork4=idav*8
       itmp2=ibcoff
       iwork=itmp2+idav*idav                                            3d16s12
       iwork4=iwork+lwork4
       ifail4=iwork4+idav*5
       ibcoff=ifail4+idav
       ieig=ibcoff                                                      3d16s12
       ibcoff=ieig+idav                                                 3d16s12
       call second(time1)
       i14=1                                                            9d9s04
       neig4=0
       tol=0d0
       iarg1=idav
       if(idwsdeb.gt.10)then
       end if
       call dws_gsumf(bc(itmp1),idav*idav)                              1d3s14
       xnrm=1d0/dfloat(numpro)                                          1d3s14
       do i=0,idav*idav-1                                               1d3s14
        bc(itmp1+i)=bc(itmp1+i)*xnrm                                    1d3s14
       end do                                                           1d3s14
       if(idwsdeb.gt.20)then
        write(6,*)('diagonalizing ')
        call prntm2(bc(itmp1),idav,idav,idav)
       end if
       call dsyevx('V','A','L',idav,bc(itmp1),idav,dum,dum,i14,idav,tol,
     $             neig4,bc(ieig),bc(itmp2),idav,bc(iwork),lwork4,
     $      bc(iwork4),
     $      bc(ifail4),info4)
       nsend=idav*(idav+1)                                              4d12s22
       call dws_bcast(bc(itmp2),nsend)                                  4d12s22
       if(idwsdeb.gt.20)then
        write(6,*)('eigenvalues: ')
        call prntm2(bc(ieig),1,idav,1)
        write(6,*)('eigenvectors: ')
        call prntm2(bc(itmp2),idav,idav,idav)
       end if
       idelz=0
       if(ncomp.gt.1)then                                               12d11s15
        ieprot=ibcoff                                                    12d11s15
        ibcoff=ieprot+idav                                               12d11s15
        call enough('opto.  3',bc,ibc)
        do i=0,idav-1                                                    12d11s15
         bc(ieprot+i)=0d0                                                12d11s15
        end do                                                           12d11s15
        do isb=1,nsymb                                                   12d11s15
         if(noc(isb).gt.0)then                                           12d11s15
          itmpv=ibcoff                                                   12d11s15
          ibcoff=itmpv+nv4(isb)*idav                                     12d11s15
          call enough('opto.  4',bc,ibc)
          call dgemm('n','n',nv4(isb),idav,idav,1d0,bc(ixyc(isb)),       12d11s15
     $         nv4(isb),bc(itmp2),idav,0d0,bc(itmpv),nv4(isb),            12d11s15
     d' opto.  1')
          nbe=(nbasdwsc(isb)/2)+1                                       6d24s16
          nbe=nbe-noc(isb)                                               12d11s15
          if(icas.eq.0.or.isb.eq.1)then                                 12d28s17
           jtmpv=itmpv+1                                                12d28s17
          else                                                          12d28s17
           jtmpv=itmpv                                                  12d28s17
          end if                                                        12d28s17
          do i=0,idav-1                                                  12d11s15
           do j=0,noc(isb)-1                                               12d11s15
            do k=nbe,nvirtc(isb)                                         12d11s15
             iad=jtmpv+k-1+nvirtc(isb)*j+nv4(isb)*i                     12d28s17
             bc(ieprot+i)=bc(ieprot+i)+bc(iad)**2                        12d11s15
 3352        format(3i5,i8,2es15.7)
            end do                                                       12d11s15
           end do                                                        12d11s15
          end do                                                         12d11s15
         end if                                                          12d11s15
        end do                                                           12d11s15
        if(idav.gt.1)then
         epxm=bc(ieprot)
         do itry=0,idav-2
          itryp=itry+1
          if(bc(ieprot+itry).lt.1d-2)then
           write(6,*)('use eigenvalue no.'),itryp
           idelz=itry
           go to 1066
          end if
         end do
         write(6,*)('at bottom of loop, epxm = '),epxm
         if(epxm.gt.1d-2)then
          write(6,*)('this doesn''t look good - stop optimization now')
          go to 218
         end if
 1066    continue
        end if
       end if                                                           12d11s15
       call dws_sync                                                    3d12s12
       if(info4.ne.0)then
        write(6,*)('on return from dsyevx, info4 = '),info4
        call dws_sync
        call dws_finalize
        stop
       end if
       eig=bc(ieig+idelz)                                               10d19s15
       eigs=eig                                                         1d27s05
       iarg1=1                                                          1d27s05
       call dws_gsumf(eigs,iarg1)                                       3d12s12
       eigs=eigs/dfloat(numpro)                                         1d27s05
       diff=abs(eig-eigs)                                               1d27s05
       rdiff=abs(diff/max(abs(eigs),1d-10))                             3d16s12
       if(rdiff.gt.1d-9)then                                            3d16s12
        write(6,*)('processors have different eigenvalues ')            1d27s05
        write(6,*)('info4 '),info4
        write(6,*)('all eigenvalues ')
        iarg1=1
        iarg2=idav
        call prntm2(bc(ieig),iarg1,iarg2,iarg1)
        val=1d0+diff*0.5d0                                              3d16s12
        call dumb(val)
        dval=val-1d0
        write(6,*)('dval: '),dval
        call dws_sync
        call dws_finalize                                               3d12s12
        stop                                                            1d27s05
       end if                                                           1d27s05
       eig=eigs                                                         1d27s05
       iarg1=idav                                                       1d27s05
       idelz=idelz*idav                                                 10d19s15
       call dws_gsumf(bc(itmp2+idelz),iarg1)                            10d19s15
       do i=1,idav                                                      1d27s05
        bc(itmp2+i-1)=bc(itmp2+idelz+i-1)/dfloat(numpro)                10d19s15
       end do                                                           1d27s05
       if(idwsdeb.gt.10)then
       writE(6,*)('ibcoff '),ibcoff,itmp1
       write(6,*)('eigen vector ')
       iarg1=idav
       iarg2=1
       call prntm2(bc(itmp2),iarg1,iarg2,iarg1)
       end if
       if(info4.ne.0)then
        write(6,*)('on return from dsyevx, info = '),info4
       end if
       xtmp=dfloat(info4)
       iarg1=1
       call dws_gsumf(xtmp,iarg1)                                       3d12s12
       if(xtmp.ne.0d0)then
        write(6,*)('xtmp ne 0'),xtmp                                     9d21s09
        call dws_finalize                                               3d12s12
        stop
       end if
       do isb=1,nsymb
        if(noc(isb).gt.0)then
         do i=1,idav
          iuse=ixyc(isb)+nv4(isb)*(i-1)+1
          if(nsymb.eq.1)then
           iuse=iuse+iacto(1)*idoub(1)
          else
           if(isb.ne.1)iuse=iuse-1
          end if
         end do
         call dgemm('n','n',nv4(isb),1,idav,1d0,bc(ixyc(isb)),nv4(isb),
     $              bc(itmp2),idav,0d0,bc(iyc(isb)),nv4(isb),           9d14s04
     d' opto.  2')
         if(nsymb.eq.1)then
          iuse=iyc(isb)+1
          iuse=iuse+iacto(1)*idoub(1)
         else
          if(isb.eq.1)then
           iuse=iyc(isb)+1
          else
           iuse=iyc(isb)
          end if
         end if
        end if                                                          9d14s04
        if(lprintx)then                                                 3d3s20
         write(6,*)('new iyc: '),isb
         call prntm2(bc(iyc(isb)),1,nv4(isb),1)
        end if                                                          3d3s20
       end do                                                           9d14s04
       ibcoff=itmp1
       if(eigo.lt.1d9)then
        change=eigo-eig
        relc=abs(change/eig)
        tuse=1d-6
        if(idavopt.eq.0)then                                            7d9s13
         if(relc.lt.1d-3)then
          if(abs(eig).gt.1d0)then
           tuse=1d-3
          else if(abs(eig).gt.1d-1)then
           tuse=1d-4
          else if(abs(eig).gt.1d-2)then
           tuse=1d-5
          end if
          if(abs(eig).lt.1d-6)tuse=1d-5
         else if(abs(eig).lt.1d-13)then                                  9d22s09
          relc=0.5d0*tuse                                                9d22s09
         end if
        else
         tuse=0.1d0**idavopt                                            7d9s13
        end if                                                          7d9s13
   66   format(1p4e15.7)                                                1d10s07
        irelc=irelc+1
        if(relc.lt.tuse)go to 218
        xtmp=iterd
        iarg1=1
        call dws_gsumf(xtmp,iarg1)                                      3d12s12
        xtmp=xtmp/dfloat(numpro)
        xtmp2=abs(xtmp-dfloat(iterd))
        if(xtmp2.gt.1d-3)then
         write(6,*)('davidson iterations out of sync '),iterd,xtmp,
     $        xtmp2
         call dws_finalize                                              3d12s12
         stop
        end if
        if(lprintx.or.iprtr(5).ne.0)write(6,219)iterd,idav,eig,change   1d14s20
  219   format(2i5,1p2e15.7)
       end if
       eigo=eig
       if(idav.eq.maxd)then                                             3d12s12
        idav=1
        do isb=1,nsymb                                                  9d14s04
         if(noc(isb).gt.0)then                                          9d14s04
          jxyc=ixyc(isb)-1                                              9d14s04
          jyc=iyc(isb)-1                                                9d14s04
          do 220 i=1,nv4(isb)                                           9d14s04
           bc(jxyc+i)=bc(jyc+i)
  220     continue
         end if                                                         9d14s04
        end do                                                          9d14s04
        go to 202
       end if
       idavp=idav+1
       ibmat1=ibmat10                                                   9d14s04
       if(icas.ne.0)then                                                12d28s17
       do isb=1,nsymb                                                   12d28s17
        ixycu(isb)=ibig2(isb)+nv4(isb)*(idavp-1)                        12d28s17
        do i=0,nv4(isb)-1                                               3d3s20
         bc(iycb(isb)+i)=bc(iyc(isb)+i)                                 3d3s20
        end do                                                          3d3s20
       end do                                                           12d28s17
       call second(time1)
       call hessyccas(iyc,ixycu,iamatx,ihessa,ihessb,ihessc,idoub,iacto,12d28s17
     $     noc,nvirtc,nsymb,xlmi2,bc,ibc)                               11d10s22
       do isb=1,nsymb                                                   12d28s17
        call dws_gsumf(bc(ixycu(isb)),nv4(isb))                         12d28s17
       end do                                                           12d28s17
       call second(time2)                                               4d23s18
       telap=time2-time1-tovr
       thc=thc+telap
       else                                                             12d28s17
       do isb=1,nsymb
        if(noc(isb).gt.0)then
         ii=ibig2(isb)+nv4(isb)*(idavp-1)-1
         do idws=1,nv4(isb)
          bc(ii+idws)=0d0
         end do
        end if
       end do
       ibmat1=ibmat10
       do isb=1,nsymb                                                   9d14s04
        ibmat2=ibmat1+noc(isb)*noc(isb)                                 9d16s04
        do isa=1,nsymb                                                  9d16s04
         if(noc(isb).gt.0.and.noc(isa).gt.0)then                        9d16s04
          call second(time1)
          itmp=ibcoff                                                    9d15s04
          ibcoff=itmp+nv4(isb)                                           9d15s04
          call enough('opto.  5',bc,ibc)
          if(isa.eq.isb)then
           nvuse=nv4(isb)                                                9d15s04
          else                                                           9d15s04
           nvuse=0                                                       9d15s04
          end if                                                         9d15s04
          call second(time1)
         if(idwsdeb.gt.10)then
         writE(6,*)('ibcoff '),ibcoff
         write(6,*)('calling hessycmo '),nowpro
         write(6,*)('whats going into hessycmo ')
         iarg1=1
         iarg2=nv4(isa)
         call prntm2(bc(iyc(isa)),iarg1,iarg2,iarg1)
         write(6,*)('for block '),isa,isb
         end if
          if(isb.ge.isa)then                                            3d12s12
           call ilimts(nvirtc(isb),nvirtc(isa),numpro,nowpro,ilj,ihj,   8d31s15
     $         i1s,i1e,i2s,i2e)                                         9d16s04
           call hessycmo(nvuse,1,itmp,iyc(isa),noc(isa),noc(isb),       3d12s12 C
     $                  nvirtc(isa),nvirtc(isb),                        8d31s15
     $                 ibmat1,ibmat2,igmat(isa,isb),i2s,i2e,i1s,i1e,    2d6s14
     $                 xlmi2,1,bc,ibc)                                  11d10s22
          else                                                          3d12s12
           call ilimts(nvirtc(isa),nvirtc(isb),numpro,nowpro,ilj,ihj,   8d31s15
     $         i1s,i1e,i2s,i2e)                                         9d16s04
           call hessycmo(nvuse,1,itmp,iyc(isa),noc(isa),noc(isb),       3d12s12
     $                   nvirtc(isa),nvirtc(isb),                       8d31s15
     $                 ibmat1,ibmat2,igmat(isb,isa),i2s,i2e,i1s,i1e,    2d6s14
     $                 xlmi2,-1,bc,ibc)                                 11d10s22
          end if                                                        3d12s12
          ii=ibig2(isb)+nv4(isb)*(idavp-1)-1
          jj=itmp-1
          do idws=1,nv4(isb)
           bc(ii+idws)=bc(ii+idws)+bc(jj+idws)
          end do
         if(idwsdeb.gt.10)then
         writE(6,*)('so far ')
         iarg1=1
         iarg2=nv4(isb)
         call prntm2(bc(ibig2(isb)+nv4(isb)*(idavp-1)),iarg1,iarg2,
     $        iarg1)
         end if
          call second(time2)                                             10d19s04
          telap=time2-time1-tovr                                         10d19s04
          times(9,1)=times(9,1)+telap                                    10d19s04
          times(9,2)=times(9,2)+telap**2                                 10d19s04
          times(9,3)=times(9,3)+1d0                                      10d19s04
         end if                                                          9d14s04
        end do                                                          9d16s04
        ibmat1=ibmat1+noc(isb)*nbasdwsc(isb)                            8d31s15
       end do                                                           9d14s04
       end if                                                           12d28s17
       sumb=0d0
       do isb=1,nsymb
        if(noc(isb).gt.0)then
         jbig=ibig2(isb)+(idavp-1)*nv4(isb)
         if(isb.eq.1.or.icas.eq.0)sumb=sumb+bc(jbig)                    12d28s17
        end if
       end do
      end if
      call second(time1)
      ibmat1=ibmat10
      do isb=1,nsymb                                                    9d14s04
       if(noc(isb).gt.0)then                                            9d14s04
        jbig=ibig2(isb)-1+(idavp-1)*nv4(isb)                            9d14s04
        resid=sumb-bc(iyc(isb))*eig                                     9d16s04
        if(idwsdeb.gt.10)then
         write(6,*)('eig '),eig
         write(6,*)('resid '),resid,idavp,sumb,bc(iyc(isb))
        end if
        if(isb.eq.1.or.icas.eq.0)then                                   1d3s18
         if(eig.ne.0d0)then
          bc(iyc(isb))=-resid/(-eig)
         else
          bc(iyc(isb))=0d0
         end if
        end if                                                          1d3s18
        jjbig=jbig+2
        jjyc=iyc(isb)+1
        jjgdig=igdig(isb)
        if(icas.ne.0)then                                               12d28s17
          jbmat=iamatx(isb)                                               12d23s17
          jyc=iyc(isb)                                                    12d4s17
          jxyc=jbig+1                                                   1d3s18
          if(isb.eq.1)then                                                12d4s17
           jyc=jyc+1                                                      12d4s17
           jxyc=jxyc+1                                                    12d4s17
          end if                                                          12d4s17
          nn=idoub(isb)*iacto(isb)+noc(isb)*nvirtc(isb)-1
          if(nsymb.eq.1)then
           ilook=21
          else
           ilook=0
          end if
          do i=0,nn
           resid=bc(jxyc+i)-bc(jyc+i)*eig
           if(resid.eq.0d0)then                                         5d2s18
            bc(jyc+i)=0d0                                               5d2s18
           else                                                         5d2s18
            bot=xlmi2*bc(ihessd(isb)+i)                                  4d2s18
            if(bot.eq.0d0)bot=1d10                                      5d20s22
            bc(jyc+i)=-resid/(bot-eig)                                     12d4s17
           end if                                                       5d2s18
          end do
        else                                                            12d28s17
        do 221 i1=1,noc(isb)                                            3d12s12
         do 230 i2=1,nvirtc(isb)                                        8d31s15
          resid=bc(jjbig)-bc(jjyc)*eig
          jjgdig=igdig(isb)+i1-1+(i2-1)*noc(isb)                        3d12s12
          bot=xlmi2*2d0*(bc(jjgdig)-bc(ibmat1+(i1-1)*(noc(isb)+1)))     3d12s12
          if(idwsdeb.gt.10)then
           write(6,1882)bc(jjbig),bc(jjyc),bc(jjgdig),resid,
     $       bc(ibmat1+(i1-1)*(noc(isb)+1))*2d0,bot
 1882      format(1p6e15.7)
          end if
          bc(jjyc)=-resid/(bot-eig)
          jjbig=jjbig+1
          jjyc=jjyc+1
          jjgdig=jjgdig+1
  230    continue
  221   continue
        end if                                                          12d28s17
        ibmat1=ibmat1+noc(isb)*nbasdwsc(isb)                            8d31s15
       end if
      end do
      do isb=1,nsymb
       if(noc(isb).gt.0.and.idwsdeb.gt.10)then
        write(6,*)('raw trial vector for symmetry '),isb
        iarg1=nv4(isb)
        iarg2=1
        call prntm2(bc(iyc(isb)),iarg1,iarg2,iarg1)
         if(nsymb.eq.1)then
          call printa(bc(iyc(isb)+1),iacto,idoub,1,0,idoub,0,1,0,
     $         bc(ibcoff))
          call printa(bc(iyc(isb)+1+iacto(1)*idoub(1)),
     $         nvirtc,noc,1,0,noc,0,1,0,
     $         bc(ibcoff))
         else
          iuse=iyc(isb)
          if(isb.eq.1)iuse=iuse+1
          call prntm2(bc(iuse),iacto(isb),idoub(isb),iacto(isb))
          call prntm2(bc(iuse+iacto(isb)*idoub(isb)),
     $         nvirtc(isb),noc(isb),nvirtc(isb))
         end if
       end if
      end do
      do 232 j=1,idav
       dot=0d0
       doto=0d0                                                         4d11s22
       ihit=0
       do isb=1,nsymb
        if(noc(isb).gt.0)then
         if(idwsdeb.gt.10)then
         write(6,*)('yc,xyc '),bc(iyc(isb)+1),
     $        bc(ixyc(isb)+1+(j-1)*nv4(isb))
         end if
         jxyc=ixyc(isb)-1+(j-1)*nv4(isb)
         jyc=iyc(isb)-1
         if(idwsdeb.gt.10)then
         write(6,*)('project of symmetry '),isb,jxyc,jyc
         iarg1=nv4(isb)
         iarg2=1
         call prntm2(bc(jyc+1),iarg1,iarg2,iarg1)
         if(nsymb.eq.1)then
          call printa(bc(jyc+2),iacto,idoub,1,0,idoub,0,1,0,
     $         bc(ibcoff))
          call printa(bc(jyc+2+iacto(1)*idoub(1)),
     $         nvirtc,noc,1,0,noc,0,1,0,
     $         bc(ibcoff))
         else
          iuse=jyc+1
          if(isb.eq.1)iuse=iuse+1
          call prntm2(bc(iuse),iacto(isb),idoub(isb),iacto(isb))
          call prntm2(bc(iuse+iacto(isb)*idoub(isb)),
     $         nvirtc(isb),noc(isb),nvirtc(isb))
         end if
         writE(6,*)('against ')
         call prntm2(bc(jxyc+1),iarg1,iarg2,iarg1)
         if(nsymb.eq.1)then
          call printa(bc(jxyc+2),iacto,idoub,1,0,idoub,0,1,0,
     $         bc(ibcoff))
          call printa(bc(jxyc+2+iacto(1)*idoub(1)),
     $         nvirtc,noc,1,0,noc,0,1,0,
     $         bc(ibcoff))
         else
          iuse=jxyc+1
          if(isb.eq.1)iuse=iuse+1
          call prntm2(bc(iuse),iacto(isb),idoub(isb),iacto(isb))
          call prntm2(bc(iuse+iacto(isb)*idoub(isb)),
     $         nvirtc(isb),noc(isb),nvirtc(isb))
         end if
         end if
         if(icas.eq.0.or.isb.eq.1)then                                  12d28s17
          kmin=2                                                        12d28s17
         else                                                           12d28s17
          kmin=1                                                        12d28s17
         end if                                                         12d28s17
         do 233 k=kmin,nv4(isb)                                         12d28s17
          dot=dot+bc(jyc+k)*bc(jxyc+k)
          doto=doto+bc(jxyc+k)**2                                       4d11s22
  233    continue
         if(ihit.eq.0.and.(icas.eq.0.or.isb.eq.1))then                  12d28s17
          dot=dot+bc(jyc+1)*bc(jxyc+1)
          doto=doto+bc(jxyc+1)**2                                       4d11s22
          ihit=1
         end if
        end if
       end do
       if(idwsdeb.gt.10)writE(6,*)('dot '),dot,doto
       do isb=1,nsymb
        if(noc(isb).gt.0)then
         if(idwsdeb.gt.10)then
         write(6,*)('yc,xyc '),bc(iyc(isb)+1),
     $        bc(ixyc(isb)+1+(j-1)*nv4(isb))
         end if
         jyc=iyc(isb)-1
         jxyc=ixyc(isb)-1+(j-1)*nv4(isb)
         do 234 k=1,nv4(isb)
          bc(jyc+k)=bc(jyc+k)-dot*bc(jxyc+k)
  234    continue
         if(idwsdeb.gt.10)then
         write(6,*)('yields '),isb
         call prntm2(bc(jyc+1),iarg1,iarg2,iarg1)
         if(nsymb.eq.1)then
          call printa(bc(jyc+2),iacto,idoub,1,0,idoub,0,1,0,
     $         bc(ibcoff))
          call printa(bc(jyc+2+iacto(1)*idoub(1)),
     $         nvirtc,noc,1,0,noc,0,1,0,
     $         bc(ibcoff))
         else
          iuse=jyc+1
          if(isb.eq.1)iuse=iuse+1
          call prntm2(bc(iuse),iacto(isb),idoub(isb),iacto(isb))
          call prntm2(bc(iuse+iacto(isb)*idoub(isb)),
     $         nvirtc(isb),noc(isb),nvirtc(isb))
         end if
         end if
        end if
       end do
  232 continue
      xnorm=0d0
      idav=idav+1
      ihit=0
      do isb=1,nsymb
       if(noc(isb).gt.0)then
        jyc=iyc(isb)-1
        if(idwsdeb.gt.10)write(6,*)('yc '),bc(iyc(isb)+1)
        if(isb.eq.1.or.icas.eq.0)then                                   12d28s17
         kmin=2                                                         12d28s17
        else                                                            12d28s17
         kmin=1                                                         12d28s17
        end if                                                          12d28s17
        do k=kmin,nv4(isb)                                              12d28s17
         xnorm=xnorm+bc(jyc+k)**2
        end do
        if(ihit.eq.0)then
         xnorm=xnorm+bc(iyc(isb))**2
         ihit=1
        end if
       end if
      end do
      call dws_bcast(xnorm,1)                                           3d3s20
      if(idwsdeb.gt.10)write(6,*)('xnorm '),xnorm,idav
      if(xnorm.lt.smallest*smallest)then                                3d3s20
       if(idwsdeb.gt.10)
     $     write(6,*)('xnorm is too small ... thus converged')          3d3s20
       do isb=1,nsymb                                                   3d3s20
        do i=0,nv4(isb)-1                                               3d3s20
         bc(iyc(isb)+i)=bc(iycb(isb)+i)                                 3d3s20
        end do                                                          3d3s20
       end do                                                           3d3s20
       go to 218                                                        3d3s20
      end if                                                            3d3s20
      if(xnorm.ne.0d0)then                                              3d3s20
       xnorm=1d0/sqrt(xnorm)
      else                                                              3d3s20
       xnorm=1d0                                                        3d3s20
      end if                                                            3d3s20
      do isb=1,nsymb
       if(noc(isb).gt.0)then
        jxyc=ixyc(isb)-1+(idav-1)*nv4(isb)
        jyc=iyc(isb)-1
        do 236 i=1,nv4(isb)
         bc(jyc+i)=bc(jyc+i)*xnorm
         bc(jxyc+i)=bc(jyc+i)
  236   continue
        call dws_bcast(bc(jyc+1),nv4(isb))                              4d12s22
        call dws_bcast(bc(jxyc+1),nv4(isb))                             4d12s22
        if(idwsdeb.gt.10)then
         write(6,*)('yc,xyc '),bc(iyc(isb)+1),
     $        bc(ixyc(isb)+1+(idav-1)*nv4(isb))
        write(6,*)('xyc so far ')
        iarg1=nv4(isb)
        iarg2=idav
        call prntm2(bc(ixyc(isb)),iarg1,iarg2,iarg1)
        end if
       end if                                                           9d14s04
      end do                                                            9d14s04
      call second(time2)
      telap=time2-time1-tovr
      if(iterd.lt.150)go to 202                                          3d22s07
      if(lprintx)write(6,*)('davidson diagonalization not converged ')
      xlam=xlam*2d0                                                     4d20s06
      if(xlam.gt.1000d10)then                                              4d24s06
       write(6,*)('xlam too large: '),xlam                              9d21s09
       write(6,*)('best est of eigenvalue: '),eig,relc                  10d23s17
       call dws_sync                                                    2d20s14
       call dws_finalize                                                3d12s12
       stop
      end if
      go to 210
  218 continue
c
c     davidson iterations converged ...
c
      eps=eig
      call second(time1)
      scal=1d0/(bc(iyc(1))*xlam)
      if(idwsdeb.gt.10)write(6,*)('scal '),scal,bc(iyc(1)),xlam,eig
      dot=0d0
      rms=0d0
      ibmat1=ibmat10                                                    9d14s04
      nv4s=1                                                            9d14s04
      do isb=1,nsymb                                                    9d14s04
       if(noc(isb).gt.0)then                                            9d14s04
        nv4s=nv4s+nv4(isb)-1                                            9d14s04
        ibmat2=ibmat1+noc(isb)*noc(isb)                                 3d12s12
        jtmp1=ibmat2-1                                                    9d9s04
        if(idwsdeb.gt.10)write(6,*)('do 238'),isb
        if(icas.eq.0.or.isb.eq.1)then                                   4d3s18
         iyuse=iyc(isb)                                                 4d3s18
        else                                                            4d3s18
         iyuse=iyc(isb)-1                                               4d3s18
        end if                                                          4d3s18
        nmove=noc(isb)*nvirtc(isb)                                      4d3s18
        if(icas.ne.0)nmove=nmove+iacto(isb)*idoub(isb)                  4d3s18
        do 238 i=1,nmove                                                4d3s18
         ysav=bc(iyuse+i)                                               4d3s18
         bc(iyuse+i)=bc(iyuse+i)*scal                                   4d3s18
         rms=rms+bc(iyuse+i)**2                                         4d3s18
  238   continue
        ibmat1=ibmat1+noc(isb)*nbasdwsc(isb)                            8d31s15
       end if                                                           9d14s04
      end do                                                            9d14s04
      if(lprintx)write(6,*)('raw rms '),rms                             3d3s20
      rms=sqrt(rms/dfloat(nv4s-1))
      if(lprintx)write(6,*)('refined rms '),rms                         3d3s20
      rmslamb=rms                                                       3d22s07
      call second(time2)
      telap=time2-time1-tovr
      times(18,1)=times(18,1)+telap
      times(18,2)=times(18,2)+telap**2
      times(18,3)=times(18,3)+1d0
      if(rms.gt.0.5d0)then
       xlam=xlam*2d0
       go to 210
      end if
      if(rmslamb.lt.1d-2.and.xlam.gt.1d0)then                           3d22s07
       xlamnew=max(1d0,xlam*0.5d0)                                      1d17s20
       if(lprintx)write(6,*)('decrease xlam ')                           3d22s07
      else                                                              3d22s07
       xlamnew=xlam                                                     3d22s07
      end if                                                            3d22s07
      xlam=xlamnew                                                      3d22s07
      call second(timeb)                                                1d23s23
      if(icas.ne.0)then
c
c     for Werner updates, we want to separate occ-virt and act-doub
c     transformation. We will store the former in umat and the latter
c     in tmat.
c
       if(nlzz.ne.0)then                                                4d4s21
        jpbr=ibcoff                                                     4d5s21
        ibcoff=jpbr+mlx                                                 4d5s21
        call enough('opto.  6',bc,ibc)
        do l=0,mlx-1                                                    4d5s21
         do i=0,ibc(inbr+l)-1                                           4d5s21
          bc(ibc(ipbr+l)+i)=0d0                                         4d5s21
         end do                                                         4d5s21
        end do                                                          4d5s21
        rmsko=0d0                                                       4d19s21
        nrmsko=0                                                        4d19s21
        szx=0d0                                                         4d19s21
        do isb=1,nsymb                                                  4d4s21
         do l=0,mlx-1                                                   4d5s21
          ibc(jpbr+l)=ibc(ipbr+l)                                       4d5s21
         end do                                                         4d5s21
         iuse=iyc(isb)                                                  4d4s21
         if(isb.eq.1)iuse=iuse+1                                        4d4s21
         if(min(idoub(isb),iacto(isb)).gt.0)then                        4d4s21
          juse=iuse                                                     4d4s21
          do i=0,idoub(isb)-1                                           4d4s21
           il=ibc(iorbsym(isb)+i)-1                                     4d5s21
           do j=0,iacto(isb)-1                                          4d4s21
            jp=j+idoub(isb)                                             4d4s21
            if(ibc(iorbsym(isb)+i).ne.ibc(iorbsym(isb)+jp))then         4d5s21
             if(abs(bc(juse)).gt.szx)then                               4d19s21
              szx=abs(bc(juse))                                         4d19s21
              iszx=i
              jszx=j
              isbx=isb                                                  5d7s21
             end if                                                     4d19s21
             rmsko=rmsko+bc(juse)**2                                    4d19s21
             nrmsko=nrmsko+1                                            4d19s21
             bc(juse)=0d0                                               4d5s21
            else if(il.ge.0)then                                        4d5s21
             bc(ibc(jpbr+il))=bc(ibc(jpbr+il))+bc(juse)                 4d5s21
             ibc(jpbr+il)=ibc(jpbr+il)+1                                4d5s21
            end if                                                      4d5s21
            juse=juse+1                                                 4d4s21
           end do                                                       4d4s21
          end do                                                        4d4s21
          iuse=iuse+idoub(isb)*iacto(isb)                               4d4s21
         end if                                                         4d4s21
         if(min(nvirtc(isb),noc(isb)).gt.0)then                         4d4s21
          juse=iuse                                                     4d4s21
          do i=0,noc(isb)-1                                             4d4s21
           il=ibc(iorbsym(isb)+i)-1                                     4d5s21
           do j=0,nvirtc(isb)-1                                         4d4s21
            jp=j+noc(isb)                                               4d4s21
            if(ibc(iorbsym(isb)+i).ne.ibc(iorbsym(isb)+jp))then         4d5s21
             if(abs(bc(juse)).gt.szx)then                               4d19s21
              szx=abs(bc(juse))                                         4d19s21
              iszx=i
              jszx=j+1000
              isbx=isb                                                  5d7s21
             end if                                                     4d19s21
             rmsko=rmsko+bc(juse)**2                                    4d19s21
             nrmsko=nrmsko+1                                            4d19s21
             bc(juse)=0d0                                               4d5s21
            else if(il.ge.0)then                                        4d5s21
             bc(ibc(jpbr+il))=bc(ibc(jpbr+il))+bc(juse)                 4d5s21
             ibc(jpbr+il)=ibc(jpbr+il)+1                                4d5s21
            end if                                                      4d5s21
            juse=juse+1                                                 4d4s21
           end do                                                       4d4s21
          end do                                                        4d4s21
         end if                                                         4d4s21
        end do                                                          4d4s21
        if(nrmsko.gt.0)then                                             4d19s21
         rmsko=sqrt(rmsko/dfloat(nrmsko))                               4d19s21
         if(rmsko.gt.1d-1)then                                          5d10s21
          write(6,*)('maximum '),szx,iszx,jszx,isbx                     5d7s21
          write(6,*)('rms sz of kod rotations: '),rmsko                  4d19s21
          call dws_synca                                                5d7s21
          call dws_finalize                                             5d7s21
          stop 'opto nlzz'                                              5d7s21
         end if                                                         5d7s21
        end if                                                          4d19s21
        do l=0,mlx-1                                                    4d5s21
         if(nlzz.eq.2)then                                              4d5s21
          factavg=0.5d0                                                 4d5s21
         else                                                           4d5s21
          factavg=1d0/dfloat(2*l+3)                                     4d5s21
         end if                                                         4d5s21
         do i=0,ibc(inbr+l)-1                                           4d5s21
          bc(ibc(ipbr+l)+i)=bc(ibc(ipbr+l)+i)*factavg                   4d5s21
         end do                                                         4d5s21
         nh=ibc(inbr+l)                                                 4d5s21
        end do                                                          4d5s21
        do isb=1,nsymb                                                  4d4s21
         do l=0,mlx-1                                                   4d5s21
          ibc(jpbr+l)=ibc(ipbr+l)                                       4d5s21
         end do                                                         4d5s21
         iuse=iyc(isb)                                                  4d4s21
         if(isb.eq.1)iuse=iuse+1                                        4d4s21
         if(min(idoub(isb),iacto(isb)).gt.0)then                        4d4s21
          juse=iuse                                                     4d4s21
          do i=0,idoub(isb)-1                                           4d4s21
           il=ibc(iorbsym(isb)+i)-1                                     4d5s21
           do j=0,iacto(isb)-1                                          4d4s21
            jp=j+idoub(isb)                                             4d4s21
            if(ibc(iorbsym(isb)+i).ne.ibc(iorbsym(isb)+jp))then         4d5s21
             bc(juse)=0d0                                               4d5s21
            else if(il.ge.0)then                                        4d5s21
             bc(juse)=bc(ibc(jpbr+il))                                  4d5s21
             ibc(jpbr+il)=ibc(jpbr+il)+1                                4d5s21
            end if                                                      4d5s21
            juse=juse+1                                                 4d4s21
           end do                                                       4d4s21
          end do                                                        4d4s21
          iuse=iuse+idoub(isb)*iacto(isb)                               4d4s21
         end if                                                         4d4s21
         if(min(nvirtc(isb),noc(isb)).gt.0)then                         4d4s21
          juse=iuse                                                     4d4s21
          do i=0,noc(isb)-1                                             4d4s21
           il=ibc(iorbsym(isb)+i)-1                                     4d5s21
           do j=0,nvirtc(isb)-1                                         4d4s21
            jp=j+noc(isb)                                               4d4s21
            if(ibc(iorbsym(isb)+i).ne.ibc(iorbsym(isb)+jp))then         4d5s21
             bc(juse)=0d0                                               4d5s21
            else if(il.ge.0)then                                        4d5s21
             bc(juse)=bc(ibc(jpbr+il))                                  4d5s21
             ibc(jpbr+il)=ibc(jpbr+il)+1                                4d5s21
            end if                                                      4d5s21
            juse=juse+1                                                 4d4s21
           end do                                                       4d4s21
          end do                                                        4d4s21
         end if                                                         4d4s21
        end do                                                          4d4s21
        ibcoff=jpbr                                                     4d5s21
       end if                                                           4d4s21
       do isb=1,nsymb                                                   1d3s18
        iuse=iyc(isb)                                                   1d3s18
        if(isb.eq.1)iuse=iuse+1                                         1d3s18
        if(idoub(isb)*iacto(isb).gt.0)then                              1d3s18
         call taylor(bc(iuse),bc(itmat(isb)),iacto(isb),idoub(isb),iret,11d9s22
     $       bc,ibc)                                                    11d9s22
         if(iret.ne.0)go to 210                                         1d3s18
         iuse=iuse+idoub(isb)*iacto(isb)                                1d3s18
        end if                                                          1d3s18
        if(nvirtc(isb)*noc(isb).gt.0)then
        call taylor(bc(iuse),bc(iumat(isb)),nvirtc(isb),noc(isb),iret,  11d9s22
     $        bc,ibc)                                                   11d9s22
        end if                                                          4d3s18
       end do                                                           1d3s18
      else
      iamat1=iamat10                                                    9d16s04
      do isb=1,nsymb
       if(noc(isb).gt.0)then
        nbas4=nbasdwsc(isb)                                             8d31s15
        idwst1=ibcoff                                                    9d9s04
        ibcoff=idwst1+nbas4*nbas4                                        9d9s04
        idwst2=ibcoff                                                    9d9s04
        ibcoff=idwst2+nbas4*nbas4                                        9d9s04
        call enough('opto.  7',bc,ibc)
        iorb=iorbn(isb)                                                  9d9s04
        call second(time1)
        iunew=ibcoff
        ibcoff=iunew+nbas4*nbas4
        call enough('opto.  8',bc,ibc)
        do 239 i=1,nbas4
         junew=iunew-1+(i-1)*nbas4
         jdwst1=idwst1-1+(i-1)*nbas4
         do 240 j=1,nbas4
          bc(junew+j)=0d0
          bc(jdwst1+j)=0d0
  240    continue
         bc(junew+i)=1d0
  239   continue
        ii=1
        nv=nbas4-noc(isb)
        do 241 i=1,noc(isb)                                             3d12s12
         do 242 j=noc(isb)+1,nbas4                                      3d12s12
          bc(iunew+j-1+(i-1)*nbas4)=bc(iyc(isb)+ii)
          bc(iunew+i-1+(j-1)*nbas4)=-bc(iyc(isb)+ii)
          bc(idwst1+j-1+(i-1)*nbas4)=bc(iyc(isb)+ii)
          bc(idwst1+i-1+(j-1)*nbas4)=-bc(iyc(isb)+ii)
          ii=ii+1
  242    continue
  241   continue
        jdwst1=idwst1-1
        jdwst2=idwst2-1
        do 243 i=1,nbas4*nbas4
         bc(jdwst2+i)=bc(jdwst1+i)
  243   continue
        fact=0.5d0
        fact2=2d0
        iter=0
        idwst3=ibcoff                                                    9d9s04
        ibcoff=idwst3+nbas4*nbas4                                        9d9s04
        call enough('opto.  9',bc,ibc)
c
c     taylor series computation of exponential of matrix of orbital
c     rotations
c
  244   continue
         iter=iter+1
         if(iter.gt.20)then
          if(lprintx)write(6,*)('last rms = '),rms
          xlam=xlam*2d0                                                 11d24s04
          go to 210                                                     11d24s04
         end if
         call dgemm('n','n',nbas4,nbas4,nbas4,1d0,bc(idwst1),nbas4,
     $     bc(idwst2),nbas4,0d0,bc(idwst3),nbas4,
     d' opto.  3')
         jdwst2=idwst2-1
         jdwst3=idwst3-1
         junew=iunew-1
         do 245 i=1,nbas4*nbas4
          bc(jdwst2+i)=bc(jdwst3+i)
          bc(junew+i)=bc(junew+i)+bc(jdwst2+i)*fact
  245    continue
         do 246 i=1,nbas4
          bc(jdwst3+i)=0d0
  246    continue
         do 247 i=1,nbas4
          junew=iunew-1+(i-1)*nbas4
          do 248 j=1,nbas4
           bc(jdwst3+j)=bc(jdwst3+j)+bc(junew+j)**2
  248     continue
  247    continue
         rms=0d0
         do 249 i=1,nbas4
          rms=rms+(bc(jdwst3+i)-1d0)**2
  249    continue
         rms=sqrt(rms/dfloat(nbas4))
         if(rms.gt.1d-9)then
          fact2=fact2+1d0
          fact=fact/fact2
          go to 244
         end if
        call dgemm('n','n',nbas4,nbas4,nbas4,1d0,bc(iumat(isb)),nbas4,    9d9s04
     $      bc(iunew),nbas4,0d0,bc(idwst1),nbas4,
     d' opto.  4')
        if(idwsdeb.gt.10)then
        write(6,*)('new u before schmidt ')
        iarg1=nbas4
        call prntm2(bc(idwst1),iarg1,iarg1,iarg1)
        end if
        do 250 i=1,nbas4
         ji=idwst1-1+(i-1)*nbas4
         do 251 j=1,i-1
          jj=idwst1-1+(j-1)*nbas4
          dot=0d0
          do 252 k=1,nbas4
           dot=dot+bc(ji+k)*bc(jj+k)
  252     continue
          do 253 k=1,nbas4
           bc(ji+k)=bc(ji+k)-dot*bc(jj+k)
  253     continue
  251    continue
         dot=0d0
         do 254 k=1,nbas4
          dot=dot+bc(ji+k)**2
  254    continue
         dot=1d0/sqrt(dot)
         if(dot.gt.1d4)then
          write(6,*)('dot '),i,dot
          call dws_finalize                                             3d12s12
          stop
         end if
         do 255 k=1,nbas4
          bc(ji+k)=bc(ji+k)*dot
  255    continue
  250   continue
        if(idwsdeb.gt.10)then
        write(6,*)('new umat '),idwst1,ibcoff
        iarg1=nbas4
        call prntm2(bc(idwst1),iarg1,iarg1,iarg1)
        end if
        call second(time2)
        telap=time2-time1-tovr
        times(10,1)=times(10,1)+telap
        times(10,2)=times(10,2)+telap**2
        times(10,3)=times(10,3)+1d0
        change=0d0
        idwst3=ibcoff                                                     9d14s04
        ibcoff=idwst3+nbas4*nbas4                                         9d14s04
        do 256 i=1,nbas4
         jdwst1=idwst1-1+(i-1)*nbas4
         jdwst3=idwst3-1+(i-1)*nbas4
         do 257 j=1,nbas4
          bc(jdwst3+j)=bc(jdwst1+j)
  257    continue
         bc(jdwst3+i)=bc(jdwst3+i)-1d0
  256   continue
        if(idwsdeb.gt.10)then
        write(6,*)('whats under dwst1 '),idwst1,ibcoff,iumat(isb)
        iarg1=nbas4
        call prntm2(bc(idwst1),iarg1,iarg1,iarg1)
        end if
        jdwst2=idwst3-1
        do i=1,nbasdwsc(isb)*nbasdwsc(isb)                              8d31s15
         bc(itmat(isb)+i-1)=bc(jdwst2+i)
        end do
       end if
      end do
      end if                                                            1d3s18
      if(icas.eq.0)then                                                 1d3s18
      ntest=0
      do isb=1,nsymb
       if(noc(isb).gt.0)then
        nbas4=nbasdwsc(isb)                                             8d31s15
        ntest=ntest+nbas4*nbas4
        if(idwsdeb.gt.10)then
        writE(6,*)('update umat from ')
        iarg1=nbas4
        call prntm2(bc(iumat(isb)),iarg1,iarg1,iarg1)
        end if
        do i=1,nbas4
         j1=iumat(isb)-1+(i-1)*nbas4
         j2=itmat(isb)-1+(i-1)*nbas4
         do j=1,nbas4
          bc(j1+j)=bc(j2+j)
         end do
         bc(j1+i)=bc(j1+i)+1d0
        end do                                                          9d16s04
        if(idwsdeb.gt.10)then
        write(6,*)('to ')
        call prntm2(bc(iumat(isb)),iarg1,iarg1,iarg1)
        end if
       end if                                                           9d16s04
      end do                                                            9d16s04
      end if                                                            1d3s18
       i404=1
       eiglastgood=eigo                                                 10d19s15
       call second(timec)
       telapc=timec-timeb-tovr
       telapt=timec-timea-tovr
       timetot=timetot+telapt
       telap=timeb-timea-tovr
       ibcoff=ibcoffo                                                   3d3s20
       return
      ibcoff=ibcoffo                                                    3d3s20
      return
      end
