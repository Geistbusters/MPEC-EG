c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine gradcas(natom,ngaus,ibdat,nbasis,isym,iapair,ibstor,   5d4s22
     $     isstor,idorel,ascale,multh,ipropsym,ivecs,ieigs,iorb,        5d4s22
     $     ih0mo,morb,ixinv,ioooo,ionex,jmats,kmats,i3x,idoub,          5d19s22
     $     iact,noc,nbasisp,nbasisc,ivecso,ihessa,ihessb,ihessc,        5d4s22
     $     ihessd,iden1,j2den,lprint,inewl,nlzz,dynw,icasvec,isend,     5d16s22
     $     irecv,iptx,ipf,ixlzz,islz,iptr,dorb,sorb,mdon,mdoo,ibasis,   5d16s22
     $     icsf,nfcn,nct,nec,icsfpd,ismx,irelx,ixw1,ixw2,iptrbit,       5d16s22
     $     iorbsym,iorbsymz,ehf,potdws,inocan,idogrado,nart,denergy,    11d9s22
     $     bc,ibc,nstate,isinfo,nveccnt,icanon)                         5d5s23
c
c     to dos:
c     ok more than 1 doub orbital
c     cas0:
c     csfs!
c     ok dvec when nps ne full space
c     ok initial guess for dvec
c     make it work for multiple roots with dynamic weights
c     > ok coded for nps=everything case, but not yet debugged.
c     > ok next thing to check is if ders with fixed wgts is correct
c     return derivatives
c     ok ders of cannonicalization: need dK and dJ
c     > this method doesn't solve for dV etc, but rather a vector
c     (independent of perturbation) times dV, so I'm not sure it can
c     be used for BODC. And Shepard gets MRCI gradient without any
c     response equations? z matrix method
c     dipole ders, etc
c     ok relativistic
c     ok bodc (sqrt(2) in coordinate definition?)
c     ok for HF C1 bodc
c     ok need higher symmetry Cs
c     ok need bodc for ps is everything.
c     need to compare to molpro (or mpec) f.d.
c     save ders to forbs
c     optimization, gradients only
c     ok mode=4 density for nps ne all case
c     ok initial guess for dvec in nps=all case
c     ok mixed singles for nps ne all nona1 case
c     ok use more than 1 proc
c
      implicit real*8 (a-h,o-z)
      external second                                                   6d14s24
      include "common.hf"
      include "common.store"
      include "common.spher"
      integer*8 ibstor,isstor                                           5d6s10
      integer*1 idogrado(4)                                             7d14s22
      include "common.basis"
      include "common.print"
      logical lprint,lbodc                                              6d2s22
      character*1 cphas                                                 5d4s22
      character*2 stsym(2)                                              4d11s23
      character*4 olab(3),name                                          3d31s23
      character*9 opname                                                3d22s23
      character*6 stname(2)                                             3d22s23
      dimension ipropsym(*),iapair(3,*),multh(8,8),isou(8),nbasisp(*),  5d4s22
     $     isoub(8),itrans(8),noc(*),i4od(idbk),ionexd(idbk),iamatu(8), 5d16s22
     $     isblkxder(4,idbk),iptoh(8,8,8),isblkder(4,idbk),iamatx(8),   5d4s22
     $     iamatb(8),rmssz0(8),i4odx(idbk),i4odu(idbk),idarot(8),       5d16s22
     $    iden1d(8,4),iden1o(8),jdend(idbk),jdeno(idbk),iact(*),pnuc(2),6d2s22
     $     idoub(*),i4oduu(idbk),isblkkder(4,idbk),jmatd(idbk),         5d19s22
     $     kmatd(idbk),i3x(*),ionexdu(idbk),jmatdu(idbk),kmatdu(idbk),  6d13s22
     $     ipxder(4,8,8,8),isoua1(8),isblkderd(4,idbk),ipxderd(4,8,8,8),7d1s22
     $     nz(4),iorb(*),j2denbl(idbk),iden1(*),j2den(*),denergy(*),    8d30s22
     $     idatta(7,3),data(3),iden1i(8,2),isinfo(11,*),iosym(3),       3d31s23
     $     ivel(8,3),ivelp(3),npt(3),nveccnt(*),icanon(*)               5d5s23
      data olab/'d/dx','d/dy','d/dz'/                                   4d15s22
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      common/singcm/iuse,nff
      common/drsigncm/drsign                                            8d20s24
      if(mynowprog.eq.0)then                                            7d25s22
       call second(time0)
       write(6,*)('Hi, my name is gradcas'),time0
       write(6,*)('what we have for idogrado: '),idogrado
       write(6,*)('total number of roots = '),nart                       7d19s22
      end if
      itprt=0                                                           12d19s22
      ides=ibcoff                                                       7d28s22
      ibcoff=ides+nart                                                  7d28s22
      call enough('gradcas.  1',bc,ibc)
      icallcas=10                                                       7d28s22
      iconv=0                                                           7d14s22
      lbodc=idogrado(3).ne.0                                            7d14s22
      if(idorel.eq.0)then                                               5d4s22
       ncomp=1                                                          5d4s22
      else                                                              5d4s22
       ncomp=2                                                          5d4s22
      end if                                                            5d4s22
      if(iprtr(28).eq.0)then
       idwsdeb=0
      else
       idwsdeb=1111
      end if
       idoder=0                                                         4d4s23
      ibodc=ibcoff                                                      1d23s23
      icdc=ibcoff                                                       3d20s23
      if(lbodc)then                                                     6d2s22
       if(lprint)write(6,*)('going for bodc ... ')                      7d14s22
       mden=0                                                            3d14s23
       do i=1,nstate                                                     3d14s23
        mden=mden+((isinfo(4,i)*(isinfo(4,i)+1))/2)                      3d14s23
       end do                                                            3d14s23
       mode=3                                                           6d2s22
       ibcoff=ibodc+natom*mden                                          3d14s23
       icdir=ibcoff                                                     3d16s23
       ibcoff=icdir+natom*mden                                          3d16s23
       call enough('gradcas.  2',bc,ibc)
       do iz=ibodc,ibcoff-1                                             6d7s22
        bc(iz)=0d0                                                      6d7s22
       end do                                                           6d7s22
       ipsymtop=nsymb                                                   6d7s22
c
c     generate 2e densities for bi-linear derivative operator.
c
       ibc0=ibcoff                                                      3d15s23
       nden1i=0                                                         3d15s23
       ndentd=0                                                         3d20s23
       do i=1,nstate-1                                                  3d20s23
        do ip=i+1,nstate                                                3d20s23
         if(isinfo(2,i).eq.isinfo(2,ip).and.isinfo(3,i).eq.isinfo(3,ip))
     $        then                                                      3d27s23
          ndentd=ndentd+isinfo(4,i)*isinfo(4,ip)                         3d27s23
         end if                                                         3d27s23
        end do                                                          3d20s23
       end do                                                           3d20s23
       do isb=1,nsymb                                                   3d15s23
        iden1i(isb,1)=ibcoff                                            3d17s23
        ibcoff=iden1i(isb,1)+iact(isb)*iact(isb)*mden                   3d17s23
        nden1i=nden1i+iact(isb)*iact(isb)*mden                          3d15s23
        iden1i(isb,2)=ibcoff                                            3d17s23
        do i=1,nstate-1                                                 3d17s23
         do ip=i+1,nstate                                               3d17s23
          jsb=multh(isb,multh(isinfo(1,i),isinfo(1,ip)))                3d17s23
          ibcoff=ibcoff+isinfo(4,i)*isinfo(4,ip)*iact(isb)*iact(jsb)    3d17s23
          nden1i=nden1i+isinfo(4,i)*isinfo(4,ip)*iact(isb)*iact(jsb)    3d20s23
         end do                                                         3d17s23
        end do                                                          3d17s23
       end do                                                           3d15s23
       icdc=ibcoff                                                      3d20s23
       ibcoff=icdc+ndentd                                               3d20s23
       do is=1,nsdlk                                                    7d7s22
        nrow=iact(isblk(1,is))*iact(isblk(2,is))
        ncol=iact(isblk(3,is))*iact(isblk(4,is))
        if(min(nrow,ncol).gt.0)then                                     7d7s22
         j2denbl(is)=ibcoff                                             7d7s22
         ibcoff=ibcoff+nrow*ncol*mden                                   3d14s23
        else                                                            7d7s22
         j2denbl(is)=-1                                                 7d7s22
        end if                                                          7d7s22
       end do                                                           7d7s22
       call enough('gradcas.  3',bc,ibc)
       do iz=ibc0,ibcoff-1                                              3d15s23
        bc(iz)=0d0                                                      3d15s23
       end do                                                           3d15s23
       idvec=-1                                                         11d23s22
       call cas0(ih0mo,ioooo,noc,mynprocg,mynowprog,pnuc,               7d7s22
     $      icallcas,iden1i,j2denbl,inewl,nlzz,multh,ehf,ncomp,enobreit,3d15s23
     $      dynw,iconv,icasvec,lprint,eavg2,jmats,kmats,nvirt,          2d22s19
     $      ibc(isend),ibc(irecv),ibc(iptx),ibc(ipf),0,ixlzz,islz,      2d20s20
     $      iptr,dorb,sorb,mdon,mdoo,ibasis,icsf,nfcn,nct,nec,icsfpd,   5d16s22
     $      ibc(ismx),ibc(irelx),iorbn,morb,nbasisp,nryd,ixw1,ixw2,     5d16s22
     $      iptrbit,iorbsym,iorbsymz,0,dum,4,ih0mo,ioooo,idvec,         7d7s22
     $          idenergy,tolv,1,isblk,nsdlk,ivspace,ipxder,             7d7s22
     $            isblk,nsdlk,ipxderd,dwsumi,bc(ides),bc,ibc,bc(icdc),  4d18s23
     $      nveccnt,icanon,ncasvecsize)                                 10d11s24
       igoal=iden1i(1,1)+12*3
       call dws_gsumf(bc(ibc0),nden1i)                                  3d15s23
       do ixyz=1,3                                                      3d31s23
        iosym(ixyz)=ipropsym(ixyz)                                      3d31s23
        do i=1,7                                                        3d31s23
         idatta(i,ixyz)=0                                               3d31s23
        end do                                                          3d31s23
        idatta(3+ixyz,ixyz)=1                                           3d31s23
        data(ixyz)=1d0                                                  3d31s23
        do isb=1,nsymb                                                  3d31s23
         jsb=multh(isb,ipropsym(ixyz))                                  3d31s23
         ivel(isb,ixyz)=ibcoff                                          3d31s23
         ibcoff=ivel(isb,ixyz)+nbasisp(isb)*nbasisp(jsb)*ncomp*ncomp    3d31s23
        end do                                                          3d31s23
        ivelp(ixyz)=ixyz                                                3d31s23
        npt(ixyz)=1                                                     3d31s23
       end do                                                           3d31s23
       nbb=ibcoff-ivel(1,1)                                             3d31s23
       call enough('gradcas.vel',bc,ibc)                                3d31s23
       do iz=ivel(1,1),ibcoff-1                                         4d3s23
        bc(iz)=0d0                                                      4d3s23
       end do                                                           4d3s23
       call parap(natom,ngaus,ibdat,ivel,isym,iapair,ibstor,isstor,idum,3d31s23
     $      idwsdeb,idorel,ascale,ivelp,npt,data,idatta,iosym,3,multh,  3d31s23
     $      nbb,nbasisp,nbasdws,iorb,olab,1,bc,ibc)                     3d31s23
       do ixyz=1,3
        do isb=1,nsymb
         jsb=multh(isb,ipropsym(ixyz))                                  3d31s23
         if(min(iact(isb),iact(jsb)).gt.0)then
          itmp=ibcoff                                                   3d31s23
          ibcoff=itmp+iact(isb)*iact(jsb)                               3d31s23
          call enough('gradcas.tmp',bc,ibc)                             3d31s23
          do i=0,iact(jsb)-1                                            3d31s23
           iad1=itmp+iact(isb)*i                                        3d31s23
           iad2=ivel(isb,ixyz)+idoub(isb)+nbasdws(isb)*(i+idoub(jsb))   3d31s23
           do j=0,iact(isb)-1                                           3d31s23
            bc(iad1+j)=bc(iad2+j)                                       3d31s23
           end do                                                       3d31s23
          end do                                                        3d31s23
          do ij=0,iact(isb)*iact(jsb)-1                                 3d31s23
           bc(ivel(isb,ixyz)+ij)=bc(itmp+ij)                            3d31s23
          end do                                                        3d31s23
          ibcoff=itmp                                                   3d31s23
         end if
        end do
       end do
       tracer=0d0                                                       3d20s23
       ivelvsm=ibcoff                                                     3d31s23
       ivelvdf=ivelvsm+mden*3                                           3d31s23
       ibcoff=ivelvdf+ndentd*3                                          3d31s23
       call enough('gradcas.velv',bc,ibc)                               3d31s23
       do iz=ivelvsm,ibcoff-1                                           3d31s23
        bc(iz)=0d0                                                      3d31s23
       end do                                                           3d31s23
       do isb=1,nsymb
        if(iact(isb).gt.0)then
         if(idwsdeb.ne.0)then                                           3d31s23
          write(6,*)('den1 for symmetry '),isb
          call prntm2(bc(iden1(isb)),iact(isb),iact(isb),iact(isb))
         end if                                                         3d31s23
         do i=0,iact(isb)-1                                             3d20s23
          iad=iden1(isb)+i*(iact(isb)+1)                                3d20s23
          tracer=tracer+bc(iad)                                         3d20s23
         end do                                                         3d20s23
         if(idwsdeb.ne.0)then                                           3d31s23
          write(6,*)('tracer so far '),tracer
          write(6,*)('den1i'),iden1i(isb,1)                              3d17s23
          call prntm2(bc(iden1i(isb,1)),iact(isb),iact(isb)*mden,        3d17s23
     $        iact(isb))                                                3d17s23
         end if                                                         3d31s23
         jvelvsm=ivelvsm                                                3d31s23
         do ixyz=1,3                                                    3d31s23
          if(ipropsym(ixyz).eq.1)then                                   3d31s23
           iad=iden1i(isb,1)                                            3d31s23
           do ixd=0,mden-1                                              3d31s23
            sum=0d0                                                     3d31s23
            do ii=0,iact(isb)*iact(isb)-1                               3d31s23
             sum=sum+bc(ivel(isb,ixyz)+ii)*bc(iad+ii)                   3d31s23
            end do                                                      3d31s23
            bc(jvelvsm+ixd)=bc(jvelvsm+ixd)+sum                         3d31s23
            iad=iad+iact(isb)*iact(isb)                                 3d31s23
           end do                                                       3d31s23
          end if                                                        3d31s23
          jvelvsm=jvelvsm+mden                                          3d31s23
         end do                                                         3d31s23
         jvelvdf=ivelvdf                                                3d31s23
         iad=iden1i(isb,2)                                              3d31s23
         do i=1,nstate-1                                                3d17s23
          do ip=i+1,nstate                                              3d17s23
           if(isinfo(2,i).eq.isinfo(2,ip).and.                          4d7s23
     $          isinfo(3,i).eq.isinfo(3,ip))then                        4d7s23
            jsb=multh(isb,multh(isinfo(1,i),isinfo(1,ip)))               3d17s23
            if(iact(jsb).gt.0)then                                       3d17s23
             nad=iact(isb)*iact(jsb)                                     3d31s23
             do ixyz=1,3                                                 3d31s23
              jad=iad                                                    3d31s23
              kvelvdf=jvelvdf+ndentd*(ixyz-1)                            3d31s23
              if(ipropsym(ixyz).eq.multh(isinfo(1,i),isinfo(1,ip)))then  3d31s23
               do irp=1,isinfo(4,ip)
                do ir=1,isinfo(4,i)
                 sum=0d0                                                 3d31s23
                 do ii=0,iact(isb)*iact(jsb)-1                           3d31s23
                  sum=sum+bc(jad+ii)*bc(ivel(isb,ixyz)+ii)               3d31s23
                 end do                                                  3d31s23
                 bc(kvelvdf)=bc(kvelvdf)+sum                             3d31s23
                 kvelvdf=kvelvdf+1                                       3d31s23
                 jad=jad+nad                                                3d17s23
                end do                                                   3d31s23
               end do                                                    3d31s23
              else                                                      4d7s23
               kvelvdf=kvelvdf+isinfo(4,ip)*isinfo(4,i)                 4d7s23
               jad=jad+isinfo(4,ip)*isinfo(4,i)*nad                     4d7s23
              end if                                                     3d31s23
             end do
             iad=iad+isinfo(4,ip)*isinfo(4,i)*nad                        3d31s23
            end if                                                      4d7s23
            jvelvdf=jvelvdf+isinfo(4,ip)*isinfo(4,i)                    3d31s23
           end if                                                       3d31s23
          end do                                                        3d31s23
         end do                                                         3d31s23
         if(idwsdeb.ne.0)then                                           3d31s23
          iad=iden1i(isb,2)                                              3d17s23
          do i=1,nstate-1                                                3d17s23
           do ip=i+1,nstate                                              3d17s23
            if(isinfo(2,i).eq.isinfo(2,ip).and.                         4d7s23
     $           isinfo(3,i).eq.isinfo(3,ip))then                       4d7s23
             jsb=multh(isb,multh(isinfo(1,i),isinfo(1,ip)))               3d17s23
             if(iact(jsb).gt.0)then                                       3d17s23
              write(6,*)('transitionx density between states '),i,ip,
     $            iad-iden1i(isb,2)
              nad=iact(isb)*iact(jsb)                                      3d17s23
              do irp=1,isinfo(4,ip)
               do ir=1,isinfo(4,i)
                write(6,*)('roots '),ir,irp,iad                                3d17s23
                call prntm2(bc(iad),iact(isb),iact(jsb),iact(isb))         3d17s23
                iad=iad+nad                                                3d17s23
               end do                                                      3d17s23
              end do                                                       3d17s23
             end if                                                       3d17s23
            end if                                                      4d7s23
           end do                                                        3d17s23
          end do                                                         3d17s23
         end if                                                         3d31s23
        end if
       end do
       if(lprint)then                                                   4d7s23
        ihead=0                                                         4d19s23
        jvelvsm=ivelvsm                                                 4d7s23
        do ixyz=1,3                                                     4d7s23
         do i=1,nstate                                                   4d7s23
          nlab=0                                                         4d7s23
          do j=6,11                                                      4d7s23
           if(isinfo(j,i).ne.0)nlab=j                                    4d7s23
          end do                                                         4d7s23
          stname(1)='      '                                             4d7s23
          do j=6,nlab                                                    4d7s23
           jm=j-5                                                        4d7s23
           stname(1)(jm:jm)=char(isinfo(j,i))                            4d7s23
          end do                                                         4d7s23
          stsym(1)='  '                                                 4d11s23
          if(nlzz.eq.2.and.stname(1)(1:1).ne.'S')then                   4d11s23
           if(stname(1)(1:1).eq.'P')then                                4d11s23
            if(isinfo(1,i).eq.2)then                                    4d11s23
             stsym(1)='X '                                              4d11s23
            else if(isinfo(1,i).eq.3)then                               4d11s23
             stsym(1)='Y '                                              4d11s23
            else if(isinfo(1,i).eq.6)then                               4d11s23
             stsym(1)='XZ'                                              4d11s23
            else if(isinfo(1,i).eq.7)then                               4d11s23
             stsym(1)='YZ'                                              4d11s23
            end if                                                      4d11s23
           else                                                         4d11s23
            if(isinfo(1,i).eq.5.or.                                     4d11s23
     $           (isinfo(1,i).eq.1.and.nsymb.eq.4))then                 4d11s23
             stsym(1)='Z '                                              4d11s23
            else if(isinfo(1,i).eq.4)then                               4d11s23
             stsym(1)='XY'                                              4d11s23
            end if                                                      4d11s23
           end if                                                       4d11s23
          end if                                                        4d11s23
          do ir1=1,isinfo(4,i)                                           4d7s23
           do ir2=1,ir1                                                  4d7s23
            if(abs(bc(jvelvsm)).gt.1d-10)then                            4d7s23
             if(ihead.eq.0)then                                         4d19s23
              write(6,*)('> electron velocity matrix elements:')        4d19s23
             end if                                                     4d19s23
             ihead=1                                                    4d19s23
             write(6,537)ir1,isinfo(2,i),stname(1),stsym(1),olab(ixyz), 4d17s23
     $          ir2,isinfo(2,i),stname(1),stsym(1),bc(jvelvsm)          4d17s23
  537        format('<',2i2,a6,a2,'|',a4,'|',2i2,a6,a2,'> = ',          4d11s23
     $            es15.7,2i5)                                           4d11s23
            end if                                                       4d7s23
            jvelvsm=jvelvsm+1                                            4d7s23
           end do                                                        4d7s23
          end do                                                         4d7s23
         end do                                                          4d7s23
        end do                                                          4d7s23
        if(ndentd.gt.0)then
         jvelvdf=ivelvdf                                                4d7s23
         do ixyz=1,3                                                    4d7s23
          do i=1,nstate-1                                               4d7s23
           nlab=0                                                          12d28s19
           do j=6,11                                                       12d31s19
            if(isinfo(j,i).ne.0)nlab=j                                     12d28s19
           end do                                                          12d28s19
           stname(1)='      '                                           4d7s23
           do j=6,nlab                                                  4d7s23
            jm=j-5                                                      4d7s23
            stname(1)(jm:jm)=char(isinfo(j,i))                          4d7s23
           end do                                                       4d7s23
           stsym(1)='  '                                                 4d11s23
           if(nlzz.eq.2.and.stname(1)(1:1).ne.'S')then                   4d11s23
            if(stname(1)(1:1).eq.'P')then                                4d11s23
             if(isinfo(1,i).eq.2)then                                    4d11s23
              stsym(1)='X '                                              4d11s23
             else if(isinfo(1,i).eq.3)then                               4d11s23
              stsym(1)='Y '                                              4d11s23
             else if(isinfo(1,i).eq.6)then                               4d11s23
              stsym(1)='XZ'                                              4d11s23
             else if(isinfo(1,i).eq.7)then                               4d11s23
              stsym(1)='YZ'                                              4d11s23
             end if                                                      4d11s23
            else                                                         4d11s23
             if(isinfo(1,i).eq.5.or.                                     4d11s23
     $           (isinfo(1,i).eq.1.and.nsymb.eq.4))then                 4d11s23
              stsym(1)='Z '                                              4d11s23
             else if(isinfo(1,i).eq.4)then                               4d11s23
              stsym(1)='XY'                                              4d11s23
             end if                                                      4d11s23
            end if                                                       4d11s23
           end if                                                        4d11s23
           do ip=i+1,nstate                                             4d7s23
            if(isinfo(2,ip).eq.isinfo(2,i).and.isinfo(3,ip).eq.         4d7s23
     $          isinfo(3,i))then                                        4d7s23
             if(multh(isinfo(1,i),isinfo(1,ip)).eq.ipropsym(ixyz))then  4d7s23
              nlab=0                                                     4d7s23
              do j=6,11                                                  4d7s23
               if(isinfo(j,ip).ne.0)nlab=j                               4d7s23
              end do                                                     4d7s23
              stname(2)='      '                                         4d7s23
              do j=6,nlab                                                4d7s23
               jm=j-5                                                    4d7s23
               stname(2)(jm:jm)=char(isinfo(j,ip))                       4d7s23
              end do                                                     4d7s23
              stsym(2)='  '                                                 4d11s23
              if(nlzz.eq.2.and.stname(2)(1:1).ne.'S')then                   4d11s23
               if(stname(2)(1:1).eq.'P')then                                4d11s23
                if(isinfo(1,ip).eq.2)then                                    4d11s23
                 stsym(2)='X '                                              4d11s23
                else if(isinfo(1,ip).eq.3)then                               4d11s23
                 stsym(2)='Y '                                              4d11s23
                else if(isinfo(1,ip).eq.6)then                               4d11s23
                 stsym(2)='XZ'                                              4d11s23
                else if(isinfo(1,ip).eq.7)then                               4d11s23
                 stsym(2)='YZ'                                              4d11s23
                end if                                                      4d11s23
               else                                                         4d11s23
                if(isinfo(1,ip).eq.5.or.                                     4d11s23
     $           (isinfo(1,ip).eq.1.and.nsymb.eq.4))then                 4d11s23
                 stsym(2)='Z '                                              4d11s23
                else if(isinfo(1,ip).eq.4)then                               4d11s23
                 stsym(2)='XY'                                              4d11s23
                end if                                                      4d11s23
               end if                                                       4d11s23
              end if                                                        4d11s23
              do irp=1,isinfo(4,ip)                                      4d7s23
               do ir=1,isinfo(4,i)                                       4d7s23
                if(abs(bc(jvelvdf)).gt.1d-10)write(6,537)ir,             4d7s23
     $            isinfo(2,i),stname(1),stsym(1),olab(ixyz),irp,        4d12s23
     $               isinfo(2,ip),stname(2),stsym(2),bc(jvelvdf)        4d12s23
                jvelvdf=jvelvdf+1                                        4d7s23
               end do                                                    4d7s23
              end do                                                     4d7s23
             else                                                        4d7s23
              jvelvdf=jvelvdf+isinfo(4,ip)*isinfo(4,i)                   4d7s23
             end if                                                      4d7s23
            end if                                                      4d7s23
           end do                                                       4d7s23
          end do                                                        4d7s23
         end do                                                         4d7s23
        end if                                                           3d31s23
       end if                                                           4d7s23
       if(idwsdeb.ne.0)then                                             3d31s23
        write(6,*)('bi-linear 2e density ')
        do is=1,nsdlk
         nrow=iact(isblk(1,is))*iact(isblk(2,is))                        7d8s22
         ncol=iact(isblk(3,is))*iact(isblk(4,is))*mden                   3d27s23
         if(min(nrow,ncol).gt.0)then                                     3d27s23
          write(6,*)('address '),j2denbl(is),('symmetry '),(isblk(j,is),
     $       j=1,4)
          call prntm2(bc(j2denbl(is)),nrow,ncol,nrow)
         end if                                                          3d27s23
        end do
       end if                                                           3d31s23
      else                                                              6d2s22
       mode=2                                                           6d2s22
       ipsymtop=1                                                       6d7s22
      end if                                                            6d2s22
      pnuc(1)=potdws                                                    5d16s22
      imsg=ibcoff                                                       6d3s10
      ibcoff=imsg+mynnode                                               6d3s10
      idenergy=1                                                        7d28s22
      if(nlzz.eq.2)ipsymtop=1                                           4d7s23
      do ipsym=1,ipsymtop
       ivspace=0                                                        6d10s22
       do is3=1,nsymb                                                   6d27s22
        do is2=1,nsymb                                                  6d27s22
         do is1=1,nsymb                                                 6d27s22
          do j=1,4                                                      6d27s22
           ipxder(j,is1,is2,is3)=0                                      6d27s22
          end do                                                        6d27s22
         end do                                                         6d27s22
        end do                                                          6d27s22
       end do                                                           6d27s22
       ndvec=0                                                           5d31s22
       nsblkkder=-1                                                      5d4s22
       call gendertype(ipsym,multh,noc,isblkder,isblkxder,nsblkder,      5d4s22
     $        nsblkxder,isblkkder,nsblkkder,ipxder)                     6d13s22
       if(idogrado(2).ne.0)then                                         7d25s22
        call buildhesscas(ioooo,ionex,jmats,kmats,bc(ih0mo),noc,ihessa,   5d4s22
     $      ihessb,ihessc,ihessd,                                       12d4s17
     $      nvirt,ipsym,multh,iden1,j2den,idoub,iact,00,bc,ibc)         11d9s22
       end if                                                           7d25s22
       if(inocan.eq.0.and.idogrado(2).ne.0)then                         7d25s22
        nsblkkder=0                                                     5d26s22
        do isb=1,nsymb                                                  5d26s22
         if(nvirt(isb).gt.0)then                                        5d26s22
          do isa=1,nsymb                                                 5d26s22
           jsa=multh(ipsym,isa)                                         6d21s22
           if(min(noc(jsa),noc(isa)).gt.0)then                          6d7s22
            nsblkkder=nsblkkder+1                                       5d26s22
            isblkkder(1,nsblkkder)=isa                                  5d26s22
            isblkkder(2,nsblkkder)=jsa                                  5d26s22
            isblkkder(3,nsblkkder)=isb                                  5d26s22
            isblkkder(4,nsblkkder)=isb                                  5d26s22
           end if                                                       5d26s22
          end do                                                         5d26s22
         end if                                                         5d26s22
        end do                                                          5d26s22
       end if                                                           5d26s22
       if(ipsym.eq.1.and.idogrado(2).ne.0)then                          12d19s22
c
c     dipole moments
c     quadrapole moments
c     q0, rel(q2) totally symmetric
c     rel(q1):zx,im(q1):zy,im(q2):xy
c
        do ikind=1,8                                                     8d3s23
         if(ikind.le.3)then                                              8d3s23
          ixyz=ikind                                                     8d3s23
          ipkind=ipropsym(ixyz)                                          8d3s23
         else if(ikind.le.5)then                                         8d3s23
          ipkind=1                                                       8d3s23
         else if(ikind.eq.6)then                                         8d3s23
          ipkind=multh(ipropsym(1),ipropsym(3))                          8d3s23
         else if(ikind.eq.7)then                                         8d3s23
          ipkind=multh(ipropsym(2),ipropsym(3))                          8d3s23
         else if(ikind.eq.8)then                                         8d3s23
          ipkind=multh(ipropsym(1),ipropsym(2))                          8d3s23
         end if                                                          8d3s23
         if(ipkind.eq.1)then                                            8d3s23
          iptxx=1                                                        8d30s22
          if(ikind.le.3)then                                            8d3s23
           if(lprint)write(6,575)olab(ixyz)(4:4)                         1d4s23
  575      format('for ',a1,' component of dipole moment')                8d30s22
           npt(1)=1                                                      3d31s23
           iosym(1)=ipropsym(ixyz)                                       3d31s23
           name='mu  '                                                   3d31s23
           name(3:3)=olab(ixyz)(4:4)                                     3d31s23
           do i=1,7                                                       8d30s22
            idatta(i,1)=0                                                3d31s23
           end do                                                         8d30s22
           idatta(ixyz,1)=1                                              3d31s23
           data(1)=-1d0                                                  3d31s23
          else                                                          8d3s23
           if(ikind.eq.4)then                                           8d3s23
            npt(1)=3                                                      3d31s23
            iosym(1)=1                                                  8d3s23
            do j=1,3                                                    8d3s23
             do i=1,7                                                   8d3s23
              idatta(i,j)=0                                             8d3s23
             end do                                                     8d3s23
            end do                                                         8d30s22
            idatta(1,1)=2                                               8d3s23
            data(1)=1d0/sqrt(6d0)                                       8d3s23
            idatta(2,2)=2                                               8d3s23
            data(2)=1d0/sqrt(6d0)                                       8d3s23
            idatta(3,3)=2                                               8d3s23
            data(3)=-2d0/sqrt(6d0)                                       8d3s23
            name='Q0  '                                                 8d3s23
           else if(ikind.eq.5)then                                      8d3s23
            npt(1)=2                                                      3d31s23
            iosym(1)=1                                                  8d3s23
            do j=1,2
             do i=1,7                                                       8d30s22
              idatta(i,j)=0                                                3d31s23
             end do                                                         8d30s22
            end do                                                      8d3s23
            idatta(1,1)=2                                               8d3s23
            data(1)=-0.5d0
            idatta(2,2)=2                                               8d3s23
            data(2)=0.5d0                                                  3d31s23
            name='ReQ2'                                                 8d3s23
           else if(ikind.eq.6)then                                      8d3s23
            name='ImQ2'                                                 8d3s23
            npt(1)=1                                                      3d31s23
            iosym(1)=1                                                  8d3s23
            do i=1,7                                                       8d30s22
             idatta(i,1)=0                                              8d3s23
            end do                                                         8d30s22
            idatta(1,1)=1                                               8d3s23
            idatta(2,1)=1                                               8d3s23
            data(1)=-1d0
           else if(ikind.eq.7)then                                      8d3s23
            name='ReQ1'                                                 8d3s23
            npt(1)=1                                                      3d31s23
            iosym(1)=1                                                  8d3s23
            do i=1,7                                                       8d30s22
             idatta(i,1)=0                                              8d3s23
            end do                                                         8d30s22
            idatta(1,1)=1                                               8d3s23
            idatta(3,1)=1                                               8d3s23
            data(1)=1d0
           else                                                         8d3s23
            name='ImQ1'                                                 8d3s23
            npt(1)=1                                                      3d31s23
            iosym(1)=1                                                  8d3s23
            do i=1,7                                                       8d30s22
             idatta(i,1)=0                                              8d3s23
            end do                                                         8d30s22
            idatta(2,1)=1                                               8d3s23
            idatta(3,1)=1                                               8d3s23
            data(1)=1d0
           end if                                                       8d3s23
          end if                                                        8d3s23
          ifmat=ibcoff                                                   8d30s22
          do isb=1,nsymb                                                 8d30s22
           iamatu(isb)=ibcoff                                            8d30s22
           ibcoff=iamatu(isb)+nbasisp(isb)*nbasisp(isb)*ncomp*ncomp      8d30s22
          end do                                                         8d30s22
          call enough('gradcas.  4',bc,ibc)
          do iz=ifmat,ibcoff-1                                           8d30s22
           bc(iz)=0d0                                                    8d30s22
          end do                                                         8d30s22
          nbb=ibcoff-ifmat                                               8d30s22
          call parap(natom,ngaus,ibdat,iamatu,isym,iapair,ibstor,        8d30s22
     $        isstor,idum,idwsdeb,idorel,ascale,iptxx,npt,data,idatta,  8d30s22
     $        iosym,1,multh,nbb,nbasisp,nbasdws,iorb,name,1,bc,ibc)     11d9s22
          dnuc=0d0
          if(ikind.le.3)then                                            8d3s23
           do ia=1,natom
            dnuc=dnuc+atnum(1,ia)*xcart(ixyz,ia)
           end do
          else if(ikind.eq.4)then                                       8d3s23
           do ia=1,natom
            dnuc=dnuc+atnum(1,ia)*(2d0*xcart(3,ia)*xcart(3,ia)          8d3s23
     $           -xcart(1,ia)*xcart(1,ia)-xcart(2,ia)*xcart(2,ia))      8d3s23
           end do
           dnuc=dnuc/sqrt(6d0)                                          8d3s23
          else if(ikind.eq.5)then                                       8d3s23
           do ia=1,natom
            dnuc=dnuc+atnum(1,ia)*(xcart(1,ia)*xcart(1,ia)              8d3s23
     $           -xcart(2,ia)*xcart(2,ia))                              8d3s23
           end do
           dnuc=dnuc*0.5d0                                              8d3s23
          else if(ikind.eq.6)then                                       8d3s23
           do ia=1,natom
            dnuc=dnuc+atnum(1,ia)*xcart(1,ia)*xcart(2,ia)               8d3s23
           end do
          else if(ikind.eq.7)then                                       8d3s23
           do ia=1,natom
            dnuc=dnuc+atnum(1,ia)*xcart(1,ia)*xcart(3,ia)               8d3s23
           end do
           dnuc=-dnuc                                                   8d3s23
          else                                                          8d3s23
           do ia=1,natom
            dnuc=dnuc+atnum(1,ia)*xcart(2,ia)*xcart(3,ia)               8d3s23
           end do
           dnuc=-dnuc                                                   8d3s23
          end if                                                        8d3s23
          pnuc(2)=dnuc                                                   8d31s22
          ss=0d0                                                         8d30s22
          jfmat=ifmat                                                    8d30s22
          do isb=1,nsymb
           if(nbasdws(isb).gt.0)then                                    11d30s22
            do i=0,nbasdws(isb)*nbasdws(isb)-1                            8d30s22
             bc(jfmat+i)=bc(iamatu(isb)+i)                                8d30s22
            end do                                                        8d30s22
            jfmat=jfmat+nbasdws(isb)*nbasdws(isb)                         8d30s22
           end if                                                       11d30s22
          end do                                                         8d30s22
          call genergy(nsymb,ifmat,idoub,iact,nbasdws,iden1,isblk,0,     8d30s22
     $        j2den,dnuc,mynprocg,mynowprog,idum,dans,ncomp,0,          12d19s22
     $        0d0,lprint,0,bc,ibc)                                      12d19s22
          n4a=0                                                          8d30s22
          do isb=1,nsymb                                                 8d30s22
           idarot(isb)=ibcoff                                            8d30s22
           jsb=multh(isb,ipsym)                                          8d30s22
           ibcoff=idarot(isb)+nbasdws(isb)*nbasdws(jsb)                  8d30s22
           n4a=n4a+nbasdws(isb)*nbasdws(jsb)                             8d30s22
          end do                                                         8d30s22
          ih0da=ibcoff                                                   8d30s22
          ibcoff=ih0da+n4a                                               8d30s22
          do isb=1,nsymb                                                 8d30s22
           jsb=isb                                                       8d30s22
           iamatu(isb)=ibcoff                                            8d30s22
           ibcoff=iamatu(isb)+iact(isb)*idoub(jsb)+noc(jsb)*nvirt(isb)   8d30s22
          end do                                                         8d30s22
          nden=0                                                         8d30s22
          do isb=1,nsymb                                                 8d30s22
           jsb=isb                                                       8d30s22
           ii=iact(isb)*iact(jsb)                                        8d30s22
           iden1d(isb,1)=ibcoff                                          8d30s22
           ibcoff=iden1d(isb,1)+ii                                       8d30s22
           nden=nden+ii                                                  8d30s22
          end do                                                         8d30s22
          do idws=1,nsdlk                                                8d30s22
           n1=isblk(1,idws)                                              8d30s22
           n2=isblk(2,idws)                                              8d30s22
           n3=isblk(3,idws)                                              8d30s22
           n4=isblk(4,idws)                                              8d30s22
           if(n1.eq.n2)then                                              8d30s22
            nnn=(iact(n1)*(iact(n1)+1))/2                                8d30s22
           else                                                          8d30s22
            nnn=iact(n1)*iact(n2)                                        8d30s22
           end if                                                        8d30s22
           if(n3.eq.n4)then                                              8d30s22
            mmm=(iact(n3)*(iact(n3)+1))/2                                8d30s22
           else                                                          8d30s22
            mmm=iact(n3)*iact(n4)                                        8d30s22
           end if                                                        8d30s22
           jdend(idws)=-1                                                8d30s22
           if(nnn*mmm.gt.0)then                                          8d30s22
            jdend(idws)=ibcoff                                           8d30s22
            ibcoff=jdend(idws)+nnn*mmm                                   8d30s22
            nden=nden+nnn*mmm                                            8d30s22
           end if                                                        8d30s22
          end do                                                         8d30s22
          iden1o(1)=ibcoff                                               8d31s22
          ibcoff=iden1o(1)+nden                                          8d31s22
          call buildcasgrad(i4od,ionexd,iamatx,idoub,iact,noc,nvirt,     8d30s22
     $   1,ifmat,multh,idwsdeb,iden1,j2den,0,isblkder,0,isblkxder,      8d30s22
     $        iamatb,namatt,rmssz0,1,1,idum,idum,bc,ibc)                11d9s22
          do ii=0,namatt-1                                               8d30s22
           bc(iamatu(1)+ii)=bc(iamatx(1)+ii)                             8d30s22
          end do                                                         8d30s22
          maxdiis=40                                                     8d30s22
          idiis1=ibcoff                                                  8d30s22
          idiis2=idiis1+maxdiis*nden                                     8d30s22
          idiis3=idiis2+maxdiis*namatt                                   8d30s22
          ibcoff=idiis3+maxdiis*nden                                     8d30s22
          ioder=ibcoff                                                   8d30s22
          ibcoff=ioder+namatt                                            8d30s22
          call enough('gradcas.  5',bc,ibc)
          do iz=ioder,ibcoff-1                                           8d30s22
           bc(iz)=0d0                                                    8d30s22
          end do                                                         8d30s22
          tolo=1d-5                                                      8d30s22
          tolv=1d-5                                                      8d30s22
          bestratio=1d10                                                11d1s24
          nprescf=0                                                     11d1s24
          idiisrest=0                                                   11d1s24
 1122     continue                                                      11d1s24
          iterov=0                                                       8d30s22
          ibctop=ibcoff                                                  8d30s22
          diffo=-1d0                                                    3d23s23
  122     continue                                                       8d30s22
           iterm=iterov                                                  8d30s22
           iterov=iterov+1                                               8d30s22
           if(idwsdeb.ne.0)
     $      write(6,*)('starting orbital-vector iteration no. '),iterov
           if(iterov.gt.maxdiis)then
            write(6,*)('too many iterations!!! ')
            call dws_synca
            call dws_finalize
            stop
           end if
           call orbdercas(iamatu,ihessa,ihessb,ihessc,ihessd,idoub,      8d30s22
     $          iact,noc,ipsym,multh,idwsdeb,bc(ioder),morb,idarot,tolo,11d10s22
     $          bc,ibc)                                                 11d10s22
           call parajkfromhd0(noc,multh,jmats,kmats,ioooo,ionex,itprt,   8d31s22
     $            idwsdeb,nvirt,ipsym,idarot,i4odu,ionexdu,nbasdws,     9d9s22
     $            isblkder,isblkxder,nsblkder,nsblkxder,isblkkder,      5d16s22
     $            -10,jmatd,kmatd,i3x,i4od2b,ionexd2,isblkxder1,        5d16s22
     $            0,isblkder1,0,idorel,igoal,bc,ibc)                    11d10s22
           call sym4o(i4odu,noc,isblkder,isblkxder,nsblkder,nsblkxder,  11d10s22
     $          bc,ibc)                                                 11d10s22
           call derh01b(bc(ih0mo),bc(ifmat),bc(ih0da),idarot,nsymb,      8d31s22
     $          nbasdws,idwsdeb,multh,ipsym,bc,ibc)                     11d14s22
           call trianglez(i4odu,ionexd,noc,nvirt,nsblkder,isblkder,      8d31s22
     $          0,isblkxder,idwsdeb,ipsym.eq.1,bc,ibc)                  11d10s22
           call rordr(i4odu,noc,nsblkder,isblkder,i4oduu,bc,ibc)        11d10s22
           tolgot=tolv                                                  4d27s23
           call cas0(ih0mo,ioooo,noc,mynprocg,mynowprog,pnuc,            8d31s22
     $      icallcas,iden1d,jdend,inewl,nlzz,multh,ehf,ncomp,enobreit,  5d16s22
     $      dynw,iconv,icasvec,lprint,eavg2,jmats,kmats,nvirt,          2d22s19
     $      ibc(isend),ibc(irecv),ibc(iptx),ibc(ipf),0,ixlzz,islz,      2d20s20
     $      iptr,dorb,sorb,mdon,mdoo,ibasis,icsf,nfcn,nct,nec,icsfpd,   5d16s22
     $      ibc(ismx),ibc(irelx),iorbn,morb,nbasisp,nryd,ixw1,ixw2,     5d16s22
     $      iptrbit,iorbsym,iorbsymz,0,dum,2,ih0da,i4oduu,idvec,        11d30s22
     $      idenergy,tolgot,ipsym,isblkder,nsblkder,ivspace,ipxder,     4d27s23
     $            isblkderd,nsblkderd,ipxderd,dwsumi,bc(ides),bc,ibc,   3d20s23
     $          bc(icdc),nveccnt,icanon,ncasvecsize)                    10d11s24
           if(tolgot.ne.tolv.and.lprint)then                            4d27s23
            write(6,*)('desired convergence '),tolv,
     $          ('achieved convergence '),tolgot                        4d27s23
           end if
           if(ndvec.eq.0)then                                            8d31s22
            ndvec=ibcoff+1-idvec                                         8d31s22
            ibctop=ibcoff                                               3d16s23
           end if                                                        8d31s22
           if(abs(dwsumi).gt.1d-14)then                                  8d31s22
            do isb=1,nsymb                                               7d25s22
             if(iact(isb).gt.0)then                                      7d25s22
              if(idwsdeb.gt.0)then                                       7d28s22
               write(6,*)('1e ddensity for symmetry '),isb,
     $          iden1d(isb,1)
               call prntm2(bc(iden1d(isb,1)),iact(isb),iact(isb),         6d2s22
     $          iact(isb))                                              6d2s22
               write(6,*)('undifferentiated density: ')
               call prntm2(bc(iden1(isb)),iact(isb),iact(isb),
     $          iact(isb))
              end if                                                     7d28s22
              do ii=0,iact(isb)*iact(isb)-1                              7d20s22
               bc(iden1d(isb,1)+ii)=bc(iden1d(isb,1)+ii)                 7d20s22
     $          +bc(iden1(isb)+ii)*dwsumi                               7d20s22
              end do                                                     7d20s22
              if(idwsdeb.gt.0)then                                       7d28s22
               write(6,*)('final ddensity ')
               call prntm2(bc(iden1d(isb,1)),iact(isb),iact(isb),         6d2s22
     $         iact(isb))                                               6d2s22
              end if                                                     7d28s22
             end if                                                      7d25s22
            end do                                                       7d25s22
            do idws=1,nsdlk                                                      2d28s07
             n1=isblk(1,idws)                                                 8d19s14
             n2=isblk(2,idws)                                                 8d19s14
             n3=isblk(3,idws)
             n4=isblk(4,idws)                                                 8d19s14
             if(n1.eq.n2)then
              nnn=(iact(n1)*(iact(n1)+1))/2                              5d16s22
             else
              nnn=iact(n1)*iact(n2)                                      5d16s22
             end if
             if(n3.eq.n4)then
              mmm=(iact(n3)*(iact(n3)+1))/2                              5d16s22
             else
              mmm=iact(n3)*iact(n4)                                      5d16s22
             end if
             if(nnn*mmm.gt.0)then
              if(idwsdeb.ne.0)then                                       7d28s22
               write(6,*)('2e ddensity for type '),n1,n2,n3,n4,
     $            jdend(idws)
               call prntm2(bc(jdend(idws)),nnn,mmm,nnn)
               write(6,*)('undifferentiated density: ')
               call prntm2(bc(j2den(idws)),nnn,mmm,nnn)
              end if                                                     7d28s22
              do ii=0,nnn*mmm-1                                          7d20s22
               bc(jdend(idws)+ii)=bc(jdend(idws)+ii)                     7d20s22
     $        +dwsumi*bc(j2den(idws)+ii)                                7d20s22
              end do                                                     7d20s22
              if(idwsdeb.ne.0)then                                       7d28s22
               write(6,*)('final ddensity ')
               call prntm2(bc(jdend(idws)),nnn,mmm,nnn)
              end if                                                     7d28s22
             end if                                                      5d16s22
            end do                                                       5d16s22
           end if                                                        7d25s22
           diff=0d0                                                      5d16s22
           jdiis1=idiis1+nden*iterm                                      5d19s22
           igoal=iden1o(1)+1
           do i=0,nden-1                                                 5d17s22
            term=bc(iden1o(1)+i)-bc(iden1d(1,1)+i)
            diff=diff+(bc(iden1o(1)+i)-bc(iden1d(1,1)+i))**2             6d2s22
            bc(iden1o(1)+i)=bc(iden1d(1,1)+i)                            6d2s22
            bc(jdiis1+i)=bc(iden1d(1,1)+i)                               6d2s22
           end do                                                        5d17s22
           if(iterov.eq.1)then                                           5d19s22
            jdiis3=idiis3+nden                                           5d19s22
            do i=0,nden-1                                                5d19s22
             bc(idiis3+i)=0d0                                            5d19s22
             bc(jdiis3+i)=bc(iden1d(1,1)+i)                              6d2s22
            end do                                                       5d19s22
            diff=1d30                                                   9d9s22
           end if                                                        5d19s22
           if(iterov.ge.2)then                                           5d17s22
            itern=iterm-1                                                5d17s22
            iffmt=ibcoff                                                 5d17s22
            iwgt=iffmt+nden*iterm                                        5d19s22
            icoef=iwgt+nden                                              5d17s22
            iscr=icoef+iterov                                            5d19s22
            idata=iscr+2*iterm+iterm*iterm+2*iterm*nden+nden             5d19s22
            ibcoff=idata+nden                                            5d19s22
            call enough('gradcas.  6',bc,ibc)
            do i=0,iterm-1                                               5d19s22
             jfmat=iffmt+nden*i                                          5d19s22
             iad1=idiis1+i*nden                                          5d19s22
             iad3=idiis3+i*nden                                          5d19s22
             do j=0,nden-1                                               5d19s22
              bc(jfmat+j)=bc(iad1+j)-bc(iad3+j)                          5d17s22
             end do                                                      5d17s22
            end do                                                       5d17s22
            iad1=idiis1+nden*iterm                                       5d19s22
            iad3=idiis3+nden*iterm                                       5d19s22
            do i=0,nden-1                                                5d19s22
             bc(iwgt+i)=1d0                                              5d17s22
             bc(idata+i)=bc(iad3+i)-bc(iad1+i)                           5d17s22
            end do
            do i=0,iterm-1                                               5d19s22
             jfmat=iffmt+nden*i                                          5d19s22
             do j=0,nden-1                                               5d19s22
              bc(jfmat+j)=bc(jfmat+j)+bc(idata+j)                        5d17s22
             end do                                                      5d17s22
            end do                                                       5d17s22
            iuse=0
            if(nden.gt.0)then                                           1d19s23
             if(iterov.gt.nprescf)then                                  11d1s24
              call lsqfit2(bc(iffmt),nden,iterm,bc(idata),nden,1,          5d19s22
     $         nden,bc(icoef),iterov,bc(iscr),bc(iwgt),0,rmsdiis,bc,ibc)3d27s23
             end if                                                     11d1s24
            end if                                                      1d19s23
            if(iterov.gt.nprescf)then                                   11d1s24
             sum=0d0                                                      5d17s22
             do i=0,iterm-1                                               5d19s22
              sum=sum+bc(icoef+i)
             end do
             cn=1d0-sum
             bc(icoef+iterm)=cn                                           5d19s22
             call dgemm('n','n',nden,1,iterov,1d0,bc(idiis1),nden,        5d19s22
     $           bc(icoef),iterov,0d0,bc(iden1d(1,1)),nden,             6d2s22
     d'gradcas.  1')
            end if                                                      11d1s24
            jdiis3=idiis3+iterov*nden                                    5d19s22
            do i=0,nden-1                                                5d19s22
             bc(jdiis3+i)=bc(iden1d(1,1)+i)                              6d2s22
            end do                                                       5d19s22
            ibcoff=iffmt                                                 5d17s22
           end if                                                        5d17s22
           if(nden.eq.0)then                                             5d20s22
            diff=0d0                                                     5d20s22
           else                                                          5d20s22
            diff=sqrt(diff/dfloat(nden))                                 5d17s22
           end if                                                        5d20s22
           ibcoff=ibctop                                                 5d16s22
           call dws_bcast(diff,1)                                        7d14s22
           ratioo=diff/diffo                                            4d12s23
           diffo=diff                                                   3d23s23
           if(iterov.gt.3.and.ratioo.le.bestratio)then                  11d1s24
            bestratio=ratioo                                            11d1s24
            ibestratio=iterov                                           11d1s24
           end if                                                       11d1s24
           if(ratioo.gt.5d0)then                                        11d1s24
            if(lprint)write(6,*)('DIIS iterations diverging')                   10d31s24
            if(idiisrest.eq.0)then                                      11d1s24
             if(lprint)write(6,*)('restart DIIS')                               10d31s24
             idiisrest=idiisrest+1                                      11d1s24
             bestratio=1d10                                             11d1s24
            else if(nprescf.eq.0)then                                   11d1s24
             if(lprint)
     $          write(6,*)('DIIS iterations diverging: switch to scf ')
             nprescf=50
             ibcoff=ibctop
            end if                                                      11d1s24
            if(bestratio.lt.1d0)then                                    11d1s24
             ibestm=ibestratio-1                                        10d31s24
             jdiis3=idiis3+ibestm*nden                                  10d31s24
             do i=0,nden-1                                              10d31s24
              bc(iden1d(1,1)+i)=bc(jdiis3+i)                            10d31s24
             end do                                                     10d31s24
            end if                                                      10d31s24
            go to 1122
           end if                                                       11d1s24
           if(nprescf.ne.0.and.ratioo.gt.0.99d0)then                    11d1s24
            if(lprint)write(6,*)
     $          ('stagnating scf iterations at iteration '),iterov,
     $          ('with last diff '),diff                                11d1s24
            go to 1251                                                  11d1s24
           end if                                                       11d1s24
           if(iterov.gt.35.and.ratioo.gt.0.9d0)then                     4d12s23
            if(lprint)write(6,*)('stagnating diis at iteration '),      4d13s23
     $          iterov,(' with last diff '),diff                        4d13s23
            go to 1251                                                  4d12s23
           end if                                                       4d12s23
           if(diff.gt.1d-9)then                                         4d27s23
            tolo=max(1d-12,min(tolo,diff*0.01d0))                        8d31s22
            tolv=max(1d-12,min(tolv,diff*0.01d0))                        8d31s22
            call buildcasgrad(ioooo,ionex,iamatu,idoub,iact,noc,nvirt,   5d9s22
     $   1,ih0mo,multh,idwsdeb,iden1d,jdend,nsdlk,isblk,                6d16s22
     $          nsdlk1,isblk1,iamatb,namatt,rmssz0,0,ipsym,nsblkderd,   7d1s22
     $          isblkderd,bc,ibc)                                       11d9s22
            if(idwsdeb.gt.10)then                                        5d19s22
             write(6,*)('adding density part ')
             call prntm2(bc(iamatu(1)),namatt,1,namatt)
             write(6,*)('to kinematic part')
             call prntm2(bc(iamatx(1)),namatt,1,namatt)
            end if                                                       5d19s22
            do ii=0,namatt-1                                             5d16s22
             bc(iamatu(1)+ii)=bc(iamatu(1)+ii)+bc(iamatx(1)+ii)          5d16s22
            end do                                                       5d16s22
            go to 122                                                    5d16s22
           end if                                                        5d16s22
           if(lprint)                                                   4d13s23
     $     write(6,*)('diis calculations converged after iteration '),  4d13s23
     $            iterov                                                4d13s23
 1251      continue                                                     4d12s23
          if(lprint)then                                                7d28s22
           write(6,*)('>'),name,(' field gradients: ')                  8d3s23
           if(idogrado(1).eq.0)then                                     7d13s23
            write(1)ikind,0,-1,0                                        8d3s23
           end if                                                       7d13s23
           jdes=ides-1                                                  4d7s23
           do i=1,nstate                                                4d7s23
            nlab=0                                                      4d7s23
            do j=6,11                                                   4d7s23
             if(isinfo(j,i).ne.0)nlab=j                                 4d7s23
            end do                                                      4d7s23
            stname(1)='      '                                           4d7s23
            do j=6,nlab                                                  4d7s23
             jm=j-5                                                      4d7s23
             stname(1)(jm:jm)=char(isinfo(j,i))                          4d7s23
            end do                                                       4d7s23
            stsym(1)='  '                                                 4d11s23
            if(nlzz.eq.2.and.stname(1)(1:1).ne.'S')then                   4d11s23
             if(stname(1)(1:1).eq.'P')then                                4d11s23
              if(isinfo(1,i).eq.2)then                                    4d11s23
               stsym(1)='X '                                              4d11s23
              else if(isinfo(1,i).eq.3)then                               4d11s23
               stsym(1)='Y '                                              4d11s23
              else if(isinfo(1,i).eq.6)then                               4d11s23
               stsym(1)='XZ'                                              4d11s23
              else if(isinfo(1,i).eq.7)then                               4d11s23
               stsym(1)='YZ'                                              4d11s23
              end if                                                      4d11s23
             else                                                         4d11s23
              if(isinfo(1,i).eq.5.or.                                     4d11s23
     $           (isinfo(1,i).eq.1.and.nsymb.eq.4))then                 4d11s23
               stsym(1)='Z '                                              4d11s23
              else if(isinfo(1,i).eq.4)then                               4d11s23
               stsym(1)='XY'                                              4d11s23
              end if                                                      4d11s23
             end if                                                       4d11s23
            end if                                                        4d11s23
            do iir=1,isinfo(4,i)                                        4d7s23
             write(6,517)iir,isinfo(2,i),stname(1),stsym(1),bc(jdes+iir)4d11s23
  517        format(2i2,a6,a2,x,'>',es20.12)                            4d19s23
            end do                                                      4d7s23
            jdes=jdes+isinfo(4,i)                                       4d7s23
           end do                                                       4d7s23
          end if                                                        7d28s22
          if(inocan.eq.0)then                                           9d9s22
           ibcb4=ibcoff                                                 5d20s22
           call parajkfromhd0(noc,multh,jmats,kmats,ioooo,ionex,        7d25s22
     $             itprt,                                               7d25s22
     $            idwsdeb,nvirt,ipsym,idarot,i4odu,ionexdu,nbasdws,     9d9s22
     $            isblkder,isblkxder,nsblkder,nsblkxder,isblkkder,      5d20s22
     $            nsblkkder,jmatdu,kmatdu,i3x,i4od2b,ionexd2,isblkxder1,5d20s22
     $            0,isblkder1,0,idorel,igoalqxz,bc,ibc)                 11d10s22
           call sym4o(i4odu,noc,isblkder,isblkxder,nsblkder,            7d25s22
     $             nsblkxder,bc,ibc)                                    11d10s22
           call trianglez(i4odu,ionexd,noc,nvirt,nsblkder,isblkder,     5d26s22
     $          0,isblkxder,idwsdeb,ipsym.eq.1,bc,ibc)                  11d10s22
           call rordr(i4odu,noc,nsblkder,isblkder,i4oduu,bc,ibc)        11d10s22
           call dcannon(ih0mo,ih0da,mynprocg,mynowprog,ioooo,i4oduu,    5d26s22
     $          kmats,kmatdu,jmats,jmatdu,noc,morb,idarot,idwsdeb,idoub,12d19s22
     $           iact,nvirt,iden1,iden1d,nsblkkder,isblkkder,bc,ibc)    11d14s22
           ibcoff=ibcb4                                                 5d20s22
          end if                                                        5d19s22
          if(mynowprog.eq.0)then                                        7d13s23
           do isb=1,nsymb                                               7d13s23
            if(nbasdws(isb).gt.0)then                                   7d13s23
             nn=nbasdws(isb)**2                                         7d13s23
             write(1)(bc(idarot(isb)+i),i=0,nn-1)                       7d13s23
            end if                                                      7d13s23
           end do                                                       7d13s23
          end if                                                        7d13s23
         end if                                                          8d30s22
        end do                                                           8d30s22
        idvec=-1                                                        11d30s22
       end if                                                           9d9s22
c
c     nuclear derivatives
c
       do ixyz=1,3
        do ia=1,natom
         if(iapair(1,ia).ge.0)then
          ipuse=ipropsym(ixyz)                                           2d3s16
          if(iapair(1,ia).eq.0)then                                      2d3s16
           npass=1                                                       2d3s16
          else                                                           2d3s16
           npass=2
          end if                                                         2d3s16
          do ipass=1,npass
           if(ipuse.eq.ipsym)then                                        3d28s16
            if(lprint)then                                               3d2s17
             write(6,*)('for nucleus no. '),ia
             if(iapair(1,ia).eq.0)then                                      2d3s16
              write(6,*)('symmetrically unique ')
              jbodc=ibodc+ia-1                                           6d7s22
              jbodc2=0                                                  7d7s22
             else                                                           2d3s16
              write(6,*)('symmetrically related to nucleus '),           4d19s23
     $             iapair(1,ia)                                         4d19s23
              jbodc=ibodc+ia-1                                           6d7s22
              jbodc2=ibodc+iapair(1,ia)-1                               7d7s22
             end if
            end if                                                       3d2s17
            ixyzm=ixyz-1
            idersign=1                                                   4d21s22
            if(npass.eq.1)then
             cphas=' '                                                    2d3s16
            else if(ipass.eq.1)then                                       2d3s16
             cphas='+'                                                    2d3s16
            else
             cphas='-'                                                    2d3s16
             idersign=2                                                  4d21s22
            end if                                                        2d3s16
            if(lprint)then                                               4d28s22
             opname(1:4)=olab(ixyz)                                     3d22s23
             write(opname(5:6),202)ia                                   3d22s23
  202        format(i2)                                                 3d22s23
             if(iapair(1,ia).eq.0)then                                   4d28s22
              write(6,1)olab(ixyz),ia,ipuse                              4d28s22
              opname(7:9)='   '                                         3d22s23
              if(idogrado(1).eq.0)then                                  7d28s22
               write(1)ixyz,ia,ipass,0                                  7d28s22
              end if                                                    7d28s22
             else                                                        4d28s22
              write(6,11)olab(ixyz),ia,cphas,iapair(1,ia),ipuse          4d28s22
              opname(7:7)=cphas                                         3d22s23
              write(opname(8:9),202)iapair(1,ia)                        3d22s23
              if(idogrado(1).eq.0)then                                  7d28s22
               write(1)ixyz,ia,ipass,iapair(1,ia)                       7d28s22
              end if                                                    7d28s22
   11         format('>for operator ',a4,i2,a1,i2,' with symmetry ',i1) 7d28s22
             end if                                                      4d28s22
            end if                                                       4d28s22
    1       format('>for operator ',a4,i2,' with symmetry ',i1)           4d28s22
            if(ipuse.eq.1)then
             nnp=0                                                       4d8s22
             nn=0                                                        5d4s22
             do isb=1,nsymb
              nh=nbasdws(isb)                                            5d4s22
              nhp=nbasisp(isb)*ncomp                                     4d8s22
              isou(isb)=nnp                                              4d8s22
              isoub(isb)=nn                                              5d4s22
              isoua1(isb)=nnp                                           6d20s22
              nnp=nnp+nhp*nhp                                            4d8s22
              nn=nn+nh*nh                                                5d4s22
             end do
             nn1=nnp                                                     4d8s22
            else
             nnp=0                                                       4d8s22
             nn=0                                                        5d4s22
             nn1=0                                                      6d20s22
             do isb=1,nsymb
              isoua1(isb)=nn1                                           6d20s22
              nn1=nn1+(nbasisp(isb)*ncomp)**2                           6d20s22
              isk=multh(isb,ipuse)                                          2d3s16
              if(isk.ge.isb)then
               isou(isb)=nnp                                             4d8s22
               isoub(isb)=nn                                             5d4s22
               nnp=nnp+nbasisp(isb)*nbasisp(isk)*ncomp*ncomp             4d8s22
               nn=nn+nbasdws(isb)*nbasdws(isk)                           5d4s22
               isou(isk)=nnp                                             4d8s22
               isoub(isk)=nn                                             5d4s22
               nnp=nnp+nbasisp(isb)*nbasisp(isk)*ncomp*ncomp             4d8s22
               nn=nn+nbasdws(isb)*nbasdws(isk)                           5d4s22
              end if
             end do
             nn1=max(nn1,nnp)                                           6d20s22
            end if
            ih0d=ibcoff
            iovr=ih0d+nnp*2                                              4d8s22
            ibcoff=iovr+nnp*2                                            4d8s22
            iovrdk=ibcoff                                                4d4s16
            npropmat=iovrdk+nnp*2                                        5d4s22
            ibcoff=npropmat+nnp                                          5d4s22
            if(lbodc)ibcoff=ibcoff+nnp                                   6d3s22
            iovrdd=ibcoff                                                6d2s22
            if(lbodc)ibcoff=iovrdd+nn1                                  6d20s22
            call enough('gradcas.  7',bc,ibc)
            do iz=ih0d,ibcoff-1                                          4d14s22
             bc(iz)=0d0                                                  4d14s22
            end do                                                       4d14s22
            derpotn=0d0
            do ixpass=1,npass                                            4d22s16
             sig=1d0                                                     4d22s16
             if(ixpass.eq.1)then                                         4d22s16
              ja=ia                                                      4d22s16
             else                                                        4d22s16
              ja=iapair(1,ia)                                            4d22s16
              if(idersign.eq.2)sig=-1d0                                  4d21s22
             end if                                                      4d22s16
             do iai=1,natom                                              4d22s16
              if(iai.ne.ja)then                                          4d22s16
               dist=0d0                                                  4d22s16
               do jxyz=1,3                                               4d22s16
                dist=dist+(xcart(jxyz,iai)-xcart(jxyz,ja))**2            4d18s16
               end do                                                    4d22s16
               dist=1d0/dist                                             4d22s16
               dists=sqrt(dist)                                          4d22s16
               vij=atnum(1,iai)*atnum(1,ja)*dists                        4d22s16
               dvij=(xcart(ixyz,iai)-xcart(ixyz,ja))*vij*dist            4d22s16
               derpotn=derpotn+drsign*dvij*sig                          8d20s24
              end if                                                     4d22s16
             end do                                                      4d22s16
            end do                                                       4d22s16
            pnuc(2)=derpotn                                              5d16s22
c
c     actually, besides ovrdk, for bodc we need (n'|m') as well
c
            call parah0grad(natom,ngaus,ibdat,nbasis,bc(ih0d),bc(iovr),
     $         bc(iovrdk),isym,iapair,ibstor,isstor,isou,nnp,           5d31s22
     $            idwsdeb,                                              5d12s16
     $            idorel,
     $            ascale,multh,ixyz,ia,ipuse,idersign,nbasisp,lbodc,    6d2s22
     $          bc(iovrdd),isoua1,nn1,bc,ibc)                           11d10s22
            call derofxor1(bc(iovr),isou,ivecs,ipuse,idorel,             5d4s22
     $       ieigs,ixor,multh,idwsdeb,iorb,ivecso,nbasisp,isoub,bc,ibc) 11d14s22
            n4a=0                                                        5d16s22
            do isb=1,nsymb                                               5d16s22
             idarot(isb)=ibcoff                                          5d16s22
             jsb=multh(isb,ipsym)                                       6d7s22
             ibcoff=idarot(isb)+nbasdws(isb)*nbasdws(jsb)               6d7s22
             n4a=n4a+nbasdws(isb)*nbasdws(jsb)                          6d7s22
            end do                                                       5d16s22
            ih0da=ibcoff                                                 5d16s22
            ibcoff=ih0da+n4a                                             5d16s22
            do isb=1,nsymb                                               5d16s22
             jsb=multh(isb,ipsym)                                       6d7s22
             iamatu(isb)=ibcoff                                          5d16s22
             ibcoff=iamatu(isb)+iact(isb)*idoub(jsb)+noc(jsb)*nvirt(isb)6d7s22
            end do                                                       5d16s22
            nden=0                                                       5d17s22
            do isb=1,nsymb                                               5d16s22
             jsb=multh(isb,ipsym)                                       6d7s22
             ii=iact(isb)*iact(jsb)                                     6d7s22
             iden1d(isb,1)=ibcoff                                          5d16s22
             ibcoff=iden1d(isb,1)+ii                                       5d17s22
             nden=nden+ii                                                5d17s22
            end do                                                       5d16s22
            if(ipuse.eq.1)then                                          6d10s22
             do idws=1,nsdlk                                                      2d28s07
              n1=isblk(1,idws)                                                 8d19s14
              n2=isblk(2,idws)                                                 8d19s14
              n3=isblk(3,idws)
              n4=isblk(4,idws)                                                 8d19s14
              if(n1.eq.n2)then
               nnn=(iact(n1)*(iact(n1)+1))/2                              5d16s22
              else
               nnn=iact(n1)*iact(n2)                                      5d16s22
              end if
              if(n3.eq.n4)then
               mmm=(iact(n3)*(iact(n3)+1))/2                              5d16s22
              else
               mmm=iact(n3)*iact(n4)                                      5d16s22
              end if
              jdend(idws)=-1                                              5d17s22
              if(nnn*mmm.gt.0)then
               jdend(idws)=ibcoff                                         5d16s22
               ibcoff=jdend(idws)+nnn*mmm                                 5d17s22
               nden=nden+nnn*mmm                                          5d17s22
              end if                                                      5d16s22
             end do                                                       5d16s22
            else                                                        6d10s22
             nsblkderd=0                                                7d1s22
             do isd=1,nsymb                                             7d1s22
              isdp=multh(isd,ipuse)                                     7d1s22
              do isc=1,nsymb                                            7d1s22
               iscdp=multh(isc,isdp)                                    7d1s22
               do isb=1,nsymb                                           7d1s22
                isa=multh(isb,iscdp)                                    7d1s22
                if(min(iact(isa),iact(isb),iact(isc),iact(isd)).gt.0)   7d1s22
     $               then                                               7d1s22
                 do idws=1,nsblkderd                                    7d1s22
                  if(((isa.eq.isblkderd(1,idws).and.
     $                 isb.eq.isblkderd(2,idws)).or.
     $                (isa.eq.isblkderd(2,idws).and.
     $                 isb.eq.isblkderd(1,idws))).and.
     $               ((isc.eq.isblkderd(3,idws).and.
     $                 isd.eq.isblkderd(4,idws)).or.
     $                (isc.eq.isblkderd(4,idws).and.
     $                 isd.eq.isblkderd(3,idws))))go to 101             7d1s22
                  if(((isa.eq.isblkderd(3,idws).and.
     $                 isb.eq.isblkderd(4,idws)).or.
     $                (isa.eq.isblkderd(4,idws).and.
     $                 isb.eq.isblkderd(3,idws))).and.
     $               ((isc.eq.isblkderd(1,idws).and.
     $                 isd.eq.isblkderd(2,idws)).or.
     $                (isc.eq.isblkderd(2,idws).and.
     $                 isd.eq.isblkderd(1,idws))))go to 101             7d1s22
                 end do                                                 7d1s22
                 nz(1)=isa                                              7d1s22
                 nz(2)=isb                                              7d1s22
                 nz(3)=isc                                              7d1s22
                 nz(4)=isd                                              7d1s22
                 nsblkderd=nsblkderd+1                                    7d1s22
                 do j=1,4                                                 7d1s22
                  isblkderd(j,nsblkderd)=nz(j)                            7d1s22
                 end do                                                   7d1s22
                 ipxderd(1,nz(1),nz(2),nz(3))=nsblkderd
                 ipxderd(2,nz(1),nz(2),nz(3))=0
                 ipxderd(3,nz(1),nz(2),nz(3))=0
                 ipxderd(4,nz(1),nz(2),nz(3))=0
                 ipxderd(1,nz(3),nz(4),nz(1))=nsblkderd
                 ipxderd(2,nz(3),nz(4),nz(1))=0
                 ipxderd(3,nz(3),nz(4),nz(1))=0
                 ipxderd(4,nz(3),nz(4),nz(1))=1
                 if(nz(1).ne.nz(2))then                                   7d1s22
                  ipxderd(1,nz(2),nz(1),nz(3))=nsblkderd
                  ipxderd(2,nz(2),nz(1),nz(3))=1
                  ipxderd(3,nz(2),nz(1),nz(3))=0
                  ipxderd(4,nz(2),nz(1),nz(3))=0
                  ipxderd(1,nz(3),nz(4),nz(2))=nsblkderd
                  ipxderd(2,nz(3),nz(4),nz(2))=0
                  ipxderd(3,nz(3),nz(4),nz(2))=1
                  ipxderd(4,nz(3),nz(4),nz(2))=1
                 end if                                                   7d1s22
                 if(nz(3).ne.nz(4))then                                   7d1s22
                  ipxderd(1,nz(1),nz(2),nz(4))=nsblkderd
                  ipxderd(2,nz(1),nz(2),nz(4))=0
                  ipxderd(3,nz(1),nz(2),nz(4))=1
                  ipxderd(4,nz(1),nz(2),nz(4))=0
                  ipxderd(1,nz(4),nz(3),nz(1))=nsblkderd
                  ipxderd(2,nz(4),nz(3),nz(1))=1
                  ipxderd(3,nz(4),nz(3),nz(1))=0
                  ipxderd(4,nz(4),nz(3),nz(1))=1
                  if(nz(1).ne.nz(2))then
                   ipxderd(1,nz(2),nz(1),nz(4))=nsblkderd
                   ipxderd(2,nz(2),nz(1),nz(4))=1
                   ipxderd(3,nz(2),nz(1),nz(4))=1
                   ipxderd(4,nz(2),nz(1),nz(4))=0
                   ipxderd(1,nz(4),nz(3),nz(2))=nsblkderd
                   ipxderd(2,nz(4),nz(3),nz(2))=1
                   ipxderd(3,nz(4),nz(3),nz(2))=1
                   ipxderd(4,nz(4),nz(3),nz(2))=1
                  end if
                 end if                                                   7d1s22
  101            continue                                               7d1s22
                end if                                                    7d1s22
               end do                                                   7d1s22
              end do                                                    7d1s22
             end do                                                     7d1s22
             do idws=1,nsblkderd                                        7d1s22
              n1=isblkderd(1,idws)                                      7d1s22
              n2=isblkderd(2,idws)                                      7d1s22
              n3=isblkderd(3,idws)                                      7d1s22
              n4=isblkderd(4,idws)                                      7d1s22
              nnn=iact(n1)*iact(n2)                                     6d10s22
              mmm=iact(n3)*iact(n4)                                     6d10s22
              jdend(idws)=-1                                            6d10s22
              if(nnn*mmm.gt.0)then                                      6d10s22
               jdend(idws)=ibcoff                                       6d10s22
               ibcoff=jdend(idws)+nnn*mmm                               6d10s22
               nden=nden+nnn*mmm                                        6d10s22
              end if                                                    6d10s22
             end do                                                     6d10s22
            end if                                                      6d10s22
            ibcoff=ibcoff+nden                                           6d2s22
            do isb=1,nsymb                                               5d17s22
             iden1o(isb)=iden1d(isb,1)+nden                              6d2s22
            end do                                                       5d17s22
            do idws=1,nsdlk                                              5d17s22
             if(jdend(idws).gt.0)then                                    5d17s22
              jdeno(idws)=jdend(idws)+nden                               5d17s22
             else                                                        5d17s22
              jdeno(idws)=-1                                             5d17s22
             end if                                                      5d17s22
            end do
            if(lbodc)then                                                6d2s22
             do im=2,4                                                   6d2s22
              do isb=1,nsymb                                              6d2s22
               if(im.eq.4)then                                          6d15s22
                jsb=isb                                                 6d15s22
               else                                                     6d15s22
                jsb=multh(isb,ipuse)                                     6d7s22
               end if                                                   6d15s22
               nn=iact(isb)*iact(jsb)*mden                              3d14s23
               iden1d(isb,im)=ibcoff                                     6d2s22
               ibcoff=iden1d(isb,im)+nn                                  6d2s22
              end do                                                     6d2s22
             end do                                                      6d2s22
            end if                                                       6d2s22
            icalcas=0                                                   7d14s22
            call transder1(morb,ixor,ixinv,nbasdws,iorb,itrans,ipuse,    5d4s22
     $          multh,idwsdeb,bc,ibc)                                   11d10s22
            call derh01(bc(ih0mo),bc(ih0d),idorel,ipuse,multh,ih0der,    5d4s22
     $          iorb,idwsdeb,itrans,nbasisp,isou,bc,ibc)                11d14s22
c
c     we only need ders of 4o and onex,
c     unless we cannonalized orbitals, then we need ders of j and k as
c     well
c
            call parajkfromhd0(noc,multh,jmats,kmats,ioooo,ionex,itprt,  5d4s22
     $            idwsdeb,nvirt,ipuse,itrans,i4od,ionexd,nbasdws,       5d4s22
     $            isblkder,isblkxder,nsblkder,nsblkxder,isblkkder,      8d4s16
     $            nsblkkder,jmatd,kmatd,i3x,i4od2b,ionexd2,isblkxder1,  5d19s22
     $            0,isblkder1,0,igoal,bc,ibc)                           12d19s22
            do i1=1,nsymb                                                     7d22s14
             do i2=1,nsymb                                                    7d22s14
              do i3=1,nsymb                                                   7d22s14
               iptoh(i3,i2,i1)=0                                              7d22s14
              end do                                                          7d22s14
             end do                                                           7d22s14
            end do                                                            7d22s14
            ii=0
            do isb=1,nsymb                                               9d26s16
             if(nbasdws(isb).gt.0)then                                   9d26s16
              do isc=1,nsymb                                             9d27s16
               if(nbasdws(isc).gt.0)then                                 10d3s16
                isbc=multh(isc,isb)                                      9d27s16
                do isd=1,nsymb                                           9d26s16
                 if(noc(isd).gt.0)then                                   9d26s16
                  isa=multh(ipuse,multh(isbc,isd))                      7d11s22
                  if(nbasdws(isa).gt.0)then                              9d26s16
                   ii=ii+1                                               9d26s16
                   iptoh(isd,isc,isb)=ii                                 9d26s16
                  end if                                                 9d26s16
                 end if                                                  9d26s16
                end do                                                   9d26s16
               end if                                                    9d26s16
              end do                                                     9d26s16
             end if                                                      9d26s16
            end do                                                       9d26s16
c
            call paraerid(natom,ngaus,ibdat,nbasis,ihmat,iorb,noc,       5d4s22
     $            ipair,nhcolt,isym,iapair,ibstor,isstor,multh,iptoh,   4d21s22
     $            itprt,idwsdeb,idorel,ascale,ia,ixyz,ipuse,idersign,   4d21s22
     $            1,nbasisp,bc,ibc)                                     11d10s22
            call parajkfromhd1(ipair,ngaus,ibdat,nhcolt,ihmat,noc,icol,     2d22s16
     $           iorb,iapair,isstor,ibstor,multh,iptoh,imsg,jmats,kmats,      2d23s16
     $         ioooo,ionex,noc,itprt,idwsdeb,ncomp,nvirt,ipuse,itrans,  5d4s22
     $            i4od,ionexd,kmatd,jmatd,i4od2b,ionexd2,der2e,nbasdws, 5d4s22
     $            isblkder,isblkxder,nsblkder,nsblkxder,                5d9s22
     $            isblkxder1,-1,isblkder1,-1,isblkkder,                 5d19s22
     $            nsblkkder,nbasisp,1,i4od,idorel,bc,ibc)               11d10s22
            if(inocan.eq.0.and.idogrado(2).ne.0)then                    7d25s22
             do i1=1,nsymb                                               9d26s16
              do i2=1,nsymb                                              9d26s16
               do i3=1,nsymb                                             9d26s16
                iptoh(i3,i2,i1)=0                                        9d26s16
               end do                                                    9d26s16
              end do                                                     9d26s16
             end do                                                      9d26s16
             ii=0
             do isb=1,nsymb                                              9d26s16
              if(nbasdws(isb).gt.0)then                                  9d26s16
               do isc=1,nsymb                                            9d27s16
                if(noc(isc).gt.0)then                                    9d26s16
                 isbc=multh(isc,isb)                                     9d27s16
                 do isd=1,nsymb                                          9d26s16
                  if(noc(isd).gt.0)then                                  9d26s16
                   isad=multh(isbc,isd)                                  9d27s16
                   isa=multh(isad,ipuse)                                 9d27s16
                   if(nbasdws(isa).gt.0)then                             9d26s16
                    if(isa.le.isb)then
                     ii=ii+1                                             9d26s16
                     iptoh(isd,isc,isb)=ii                               9d26s16
                    end if
                   end if                                                9d26s16
                  end if                                                 9d26s16
                 end do                                                  9d26s16
                end if                                                   9d26s16
               end do                                                    9d26s16
              end if                                                     9d26s16
             end do                                                      9d26s16
             call paraeridj(natom,ngaus,ibdat,nbasis,ihmat,iorb,noc,     5d20s22
     $            ipair,nhcolt,isym,iapair,ibstor,isstor,multh,iptoh,
     $            itprt,idwsdeb,idorel,ascale,ia,ixyz,ipuse,idersign,   11d10s22
     $            nbasisp,bc,ibc)                                       11d10s22
             call parajkfromhd1j(ipair,ngaus,ibdat,nhcolt,ihmat,noc,icol     2d22s16
     $          ,iorb,iapair,isstor,ibstor,multh,iptoh,imsg,jmats,kmats,      2d23s16
     $         ioooo,ionex,noc,itprt,idwsdeb,ncomp,nvirt,ipuse,itrans,  5d20s22
     $            i4od,ionexd,kmatd,jmatd,i4od2b,ionexd2,der2e,nbasdws, 5d20s22
     $            isblkder,isblkxder,nsblkder,nsblkxder,                9d16s16
     $            isblkxder1,0,isblkder1,0,isblkkder,                   6d13s22
     $            nsblkkder,nbasisp,idorel,bc,ibc)                      11d10s22
            end if                                                       5d20s22
            call int4copy(i4od,i4odx,n4o,noc,nsblkder,isblkder,bc,ibc)  11d10s22
c
            tolv=1d-13                                                   6d2s22
            if(ndvec.ne.0)then                                           5d31s22
             idvec=ibcoff                                                5d31s22
             ibcoff=idvec+ndvec                                          5d31s22
             call enough('gradcas.  8',bc,ibc)
             do iz=idvec,ibcoff-1                                       4d6s23
              bc(iz)=0d0                                                4d6s23
             end do                                                     4d6s23
            end if                                                       5d31s22
            call trianglez(i4od,ionexd,noc,nvirt,nsblkder,isblkder,     7d18s22
     $          nsblkxder,isblkxder,idwsdeb,ipuse.eq.1,bc,ibc)          11d10s22
            if(ipuse.eq.1)then                                           5d9s22
             call rordr(i4od,noc,nsblkder,isblkder,i4oduu,bc,ibc)       11d10s22
             enobreit=0d0                                               7d11s22
             if(idogrado(1).ne.0)then                                   7d25s22
              call genergy(nsymb,ih0der,idoub,iact,nbasdws,iden1,isblk,  7d7s22
     $            nsdlk,j2den,derpotn,mynprocg,mynowprog,i4oduu,        7d28s22
     $            denergy(idenergy),ncomp,idwsdeb,enobreit,lprint,0,    12d19s22
     $            bc,ibc)                                               12d19s22
              if(lprint.and.idogrado(4).ne.10)                           7d28s22
     $             write(6,*)('de from genergy> '),denergy(idenergy)    7d28s22
              idenergy=idenergy+1                                       7d28s22
             end if                                                     7d25s22
            end if                                                       5d9s22
            if(idogrado(2).ne.0)then                                    7d25s22
             call buildcasgrad(i4od,ionexd,iamatx,idoub,iact,noc,nvirt,   5d9s22
     $   ipuse,ih0der,multh,idwsdeb,iden1,j2den,nsblkder,isblkder,      5d9s22
     $          nsblkxder,isblkxder,iamatb,namatt,rmssz0,1,1,idum,idum, 11d9s22
     $           bc,ibc)                                                11d9s22
             do ii=0,namatt-1                                             5d16s22
              bc(iamatu(1)+ii)=bc(iamatx(1)+ii)                           5d16s22
             end do                                                       5d16s22
             maxdiis=40                                                   5d19s22
             idiis1=ibcoff                                                5d17s22
             idiis2=idiis1+maxdiis*nden                                   5d17s22
             idiis3=idiis2+maxdiis*namatt                                 5d19s22
             ibcoff=idiis3+maxdiis*nden                                   5d19s22
             ioder=ibcoff                                                 5d9s22
             ibcoff=ioder+namatt                                          5d9s22
             call enough('gradcas.  9',bc,ibc)
             do iz=ioder,ibcoff-1                                         5d17s22
              bc(iz)=0d0                                                  5d17s22
             end do                                                       5d17s22
             tolo=1d-5                                                    6d2s22
             tolv=1d-5
             bestratio=1d10                                             11d1s24
             nprescf=0                                                  11d1s24
             idiisrest=0                                                11d1s24
 2122        continue                                                   11d1s24
             iterov=0                                                     5d16s22
             ibctop=ibcoff                                                5d16s22
             diffo=-1d0                                                 4d10s23
   22        continue                                                     5d16s22
             iterm=iterov                                                 5d17s22
             iterov=iterov+1                                              5d16s22
             if(idwsdeb.ne.0)
     $       write(6,*)('starting orbital-vector iteration no. '),iterov
             if(iterov.gt.maxdiis)then
              write(6,*)('too many iterations!!! ')
              call dws_synca
              call dws_finalize
              stop
             end if
             call orbdercas(iamatu,ihessa,ihessb,ihessc,ihessd,idoub,     5d9s22
     $          iact,noc,ipuse,multh,idwsdeb,bc(ioder),morb,idarot,tolo,11d10s22
     $            bc,ibc)                                               11d10s22
             call parajkfromhd0(noc,multh,jmats,kmats,ioooo,ionex,itprt,  5d4s22
     $            idwsdeb,nvirt,ipuse,idarot,i4odu,ionexdu,nbasdws,     5d20s22
     $            isblkder,isblkxder,nsblkder,nsblkxder,isblkkder,      5d16s22
     $            -10,jmatd,kmatd,i3x,i4od2b,ionexd2,isblkxder1,        5d16s22
     $            0,isblkder1,0,idorel,igoal,bc,ibc)                    11d10s22
             call sym4o(i4odu,noc,isblkder,isblkxder,nsblkder,nsblkxder,11d10s22
     $            bc,ibc)                                               11d10s22
             do i=0,n4o-1                                                 5d16s22
              bc(i4odu(1)+i)=bc(i4odu(1)+i)+bc(i4odx(1)+i)                5d16s22
             end do                                                       5d16s22
             call derh01b(bc(ih0mo),bc(ih0der),bc(ih0da),idarot,nsymb,    5d16s22
     $          nbasdws,idwsdeb,multh,ipuse,bc,ibc)                     11d14s22
             call trianglez(i4odu,ionexd,noc,nvirt,nsblkder,isblkder,    7d18s22
     $          0,isblkxder,idwsdeb,ipuse.eq.1,bc,ibc)                  11d10s22
             tolgot=tolv                                                4d27s23
             if(ipuse.eq.1)then                                           5d9s22
              call rordr(i4odu,noc,nsblkder,isblkder,i4oduu,bc,ibc)     11d10s22
              call cas0(ih0mo,ioooo,noc,mynprocg,mynowprog,pnuc,           5d16s22
     $      icallcas,iden1d,jdend,inewl,nlzz,multh,ehf,ncomp,enobreit,  5d16s22
     $      dynw,iconv,icasvec,lprint,eavg2,jmats,kmats,nvirt,          2d22s19
     $      ibc(isend),ibc(irecv),ibc(iptx),ibc(ipf),0,ixlzz,islz,      2d20s20
     $      iptr,dorb,sorb,mdon,mdoo,ibasis,icsf,nfcn,nct,nec,icsfpd,   5d16s22
     $      ibc(ismx),ibc(irelx),iorbn,morb,nbasisp,nryd,ixw1,ixw2,     5d16s22
     $      iptrbit,iorbsym,iorbsymz,0,dum,mode,ih0da,i4oduu,idvec,     6d2s22
     $      idenergy,tolgot,ipuse,isblkder,nsblkder,ivspace,ipxder,     4d27s23
     $      isblkderd,nsblkderd,ipxderd,dwsumi,bc(ides),bc,ibc,bc(icdc),4d18s23
     $             nveccnt,icanon,ncasvecsize)                          10d11s24
             else
              call cas0(ih0mo,ioooo,noc,mynprocg,mynowprog,pnuc,           5d16s22
     $      icallcas,iden1d,jdend,inewl,nlzz,multh,ehf,ncomp,enobreit,  5d16s22
     $      dynw,iconv,icasvec,lprint,eavg2,jmats,kmats,nvirt,          2d22s19
     $      ibc(isend),ibc(irecv),ibc(iptx),ibc(ipf),0,ixlzz,islz,      2d20s20
     $      iptr,dorb,sorb,mdon,mdoo,ibasis,icsf,nfcn,nct,nec,icsfpd,   5d16s22
     $      ibc(ismx),ibc(irelx),iorbn,morb,nbasisp,nryd,ixw1,ixw2,     5d16s22
     $      iptrbit,iorbsym,iorbsymz,0,dum,mode,ih0da,i4odu,idvec,      6d16s22
     $      idenergy,tolgot,ipuse,isblkder,nsblkder,ivspace,ipxder,     4d27s23
     $      isblkderd,nsblkderd,ipxderd,dwsumi,bc(ides),bc,ibc,bc(icdc),4d18s23
     $             nveccnt,icanon,ncasvecsize)                          10d11s24
             end if                                                       5d9s22
             if(tolgot.ne.tolv.and.lprint)then                            4d27s23
             write(6,*)('desired convergence '),tolv,
     $          ('achieved convergence '),tolgot                        4d27s23
            end if
             if(ndvec.eq.0)then                                         7d25s22
              ndvec=ibcoff+1-idvec                                      7d25s22
              ibctop=ibcoff                                             3d15s23
             end if                                                     7d25s22
             if(abs(dwsumi).gt.1d-14.and.ipuse.eq.1)then                7d25s22
              do isb=1,nsymb                                             7d25s22
               if(iact(isb).gt.0)then                                   7d25s22
                if(idwsdeb.gt.0)then                                    7d28s22
                 write(6,*)('1e density for symmetry '),isb,
     $               iden1d(isb,1)
                 call prntm2(bc(iden1d(isb,1)),iact(isb),iact(isb),         6d2s22
     $            iact(isb))                                            6d2s22
                 write(6,*)('undifferentiated density: ')
                 call prntm2(bc(iden1(isb)),iact(isb),iact(isb),
     $               iact(isb))
                end if                                                  7d28s22
                do ii=0,iact(isb)*iact(isb)-1                           7d20s22
                 bc(iden1d(isb,1)+ii)=bc(iden1d(isb,1)+ii)              7d20s22
     $                +bc(iden1(isb)+ii)*dwsumi                         7d20s22
                end do                                                  7d20s22
                if(idwsdeb.gt.0)then                                    7d28s22
                 write(6,*)('final density ')
                 call prntm2(bc(iden1d(isb,1)),iact(isb),iact(isb),         6d2s22
     $            iact(isb))                                            6d2s22
                end if                                                  7d28s22
               end if                                                   7d25s22
              end do                                                     7d25s22
              do idws=1,nsdlk                                                      2d28s07
               n1=isblk(1,idws)                                                 8d19s14
               n2=isblk(2,idws)                                                 8d19s14
               n3=isblk(3,idws)
               n4=isblk(4,idws)                                                 8d19s14
               if(n1.eq.n2)then
                nnn=(iact(n1)*(iact(n1)+1))/2                              5d16s22
               else
                nnn=iact(n1)*iact(n2)                                      5d16s22
               end if
               if(n3.eq.n4)then
                mmm=(iact(n3)*(iact(n3)+1))/2                              5d16s22
               else
                mmm=iact(n3)*iact(n4)                                      5d16s22
               end if
               if(nnn*mmm.gt.0)then
                if(idwsdeb.ne.0)then                                    7d28s22
                 write(6,*)('2e density for type '),n1,n2,n3,n4,
     $              jdend(idws)
                 call prntm2(bc(jdend(idws)),nnn,mmm,nnn)
                 write(6,*)('undifferentiated density: ')
                 call prntm2(bc(j2den(idws)),nnn,mmm,nnn)
                end if                                                  7d28s22
                do ii=0,nnn*mmm-1                                       7d20s22
                 bc(jdend(idws)+ii)=bc(jdend(idws)+ii)                  7d20s22
     $                +dwsumi*bc(j2den(idws)+ii)                        7d20s22
                end do                                                  7d20s22
                if(idwsdeb.ne.0)then                                    7d28s22
                 write(6,*)('final density ')
                 call prntm2(bc(jdend(idws)),nnn,mmm,nnn)
                end if                                                  7d28s22
               end if                                                      5d16s22
              end do                                                       5d16s22
             end if                                                     7d25s22
             if(idwsdeb.gt.10)then                                       7d14s22
              write(6,*)('densities so far ')
              do isb=1,nsymb                                               5d16s22
               if(iact(isb).gt.0)then
                jsb=multh(isb,ipuse)                                      6d7s22
                write(6,*)('1e density for symmetry '),isb,iden1d(isb,1)
                call prntm2(bc(iden1d(isb,1)),iact(isb),iact(jsb),         6d2s22
     $            iact(isb))                                            6d2s22
                if(lbodc)then
                 do ipx=1,2
                  write(6,*)('bodc den no. '),ipx,iden1d(isb,ipx+1)
                  call prntm2(bc(iden1d(isb,ipx+1)),iact(isb),
     $                 iact(jsb)*mden,iact(isb))
                 end do
                 write(6,*)('bodc den no. 3'),iden1d(isb,4)                             6d16s22
                 call prntm2(bc(iden1d(isb,4)),iact(isb),iact(isb)*mden,3d14s23
     $              iact(isb))
                end if
               end if                                                      5d27s22
              end do                                                       5d16s22
              if(ipuse.eq.1)then                                          6d16s22
               do idws=1,nsdlk                                                      2d28s07
                n1=isblk(1,idws)                                                 8d19s14
                n2=isblk(2,idws)                                                 8d19s14
                n3=isblk(3,idws)
                n4=isblk(4,idws)                                                 8d19s14
                if(n1.eq.n2)then
                 nnn=(iact(n1)*(iact(n1)+1))/2                              5d16s22
                else
                 nnn=iact(n1)*iact(n2)                                      5d16s22
                end if
                if(n3.eq.n4)then
                 mmm=(iact(n3)*(iact(n3)+1))/2                              5d16s22
                else
                 mmm=iact(n3)*iact(n4)                                      5d16s22
                end if
                if(nnn*mmm.gt.0)then
                 write(6,*)('2e density for type '),n1,n2,n3,n4
                 call prntm2(bc(jdend(idws)),nnn,mmm,nnn)
                end if                                                      5d16s22
               end do                                                       5d16s22
              else                                                        6d16s22
               do idws=1,nsblkderd                                        7d1s22
                n1=isblkderd(1,idws)                                      7d1s22
                n2=isblkderd(2,idws)                                      7d1s22
                n3=isblkderd(3,idws)                                      7d1s22
                n4=isblkderd(4,idws)                                      7d1s22
                nnn=iact(n1)*iact(n2)                                     6d10s22
                mmm=iact(n3)*iact(n4)                                     6d10s22
                if(nnn*mmm.gt.0)then                                      6d10s22
                 write(6,*)('2e density for type '),n1,n2,n3,n4,
     $               jdend(idws)
                 call prntm2(bc(jdend(idws)),nnn,mmm,nnn)
                end if                                                    6d10s22
               end do                                                     6d10s22
              end if                                                      6d16s22
             end if                                                      7d14s22
             diff=0d0                                                     5d16s22
             jdiis1=idiis1+nden*iterm                                     5d19s22
             do i=0,nden-1                                                5d17s22
              diff=diff+(bc(iden1o(1)+i)-bc(iden1d(1,1)+i))**2            6d2s22
              bc(iden1o(1)+i)=bc(iden1d(1,1)+i)                           6d2s22
              bc(jdiis1+i)=bc(iden1d(1,1)+i)                              6d2s22
             end do                                                       5d17s22
             if(iterov.eq.1)then                                          5d19s22
              jdiis3=idiis3+nden                                          5d19s22
              do i=0,nden-1                                               5d19s22
               bc(idiis3+i)=0d0                                           5d19s22
               bc(jdiis3+i)=bc(iden1d(1,1)+i)                             6d2s22
              end do                                                      5d19s22
              diff=1d30                                                 8d19s24
             end if                                                       5d19s22
             if(iterov.ge.2)then                                          5d17s22
              itern=iterm-1                                               5d17s22
              ifmat=ibcoff                                                5d17s22
              iwgt=ifmat+nden*iterm                                       5d19s22
              icoef=iwgt+nden                                             5d17s22
              iscr=icoef+iterov                                           5d19s22
              idata=iscr+2*iterm+iterm*iterm+2*iterm*nden+nden            5d19s22
              ibcoff=idata+nden                                           5d19s22
              call enough('gradcas. 10',bc,ibc)
              do i=0,iterm-1                                              5d19s22
               jfmat=ifmat+nden*i                                         5d19s22
               iad1=idiis1+i*nden                                         5d19s22
               iad3=idiis3+i*nden                                         5d19s22
               do j=0,nden-1                                              5d19s22
                bc(jfmat+j)=bc(iad1+j)-bc(iad3+j)                         5d17s22
               end do                                                     5d17s22
              end do                                                      5d17s22
              iad1=idiis1+nden*iterm                                      5d19s22
              iad3=idiis3+nden*iterm                                      5d19s22
              do i=0,nden-1                                               5d19s22
               bc(iwgt+i)=1d0                                             5d17s22
               bc(idata+i)=bc(iad3+i)-bc(iad1+i)                          5d17s22
              end do
              do i=0,iterm-1                                              5d19s22
               jfmat=ifmat+nden*i                                         5d19s22
               do j=0,nden-1                                              5d19s22
                bc(jfmat+j)=bc(jfmat+j)+bc(idata+j)                       5d17s22
               end do                                                     5d17s22
              end do                                                      5d17s22
              iuse=0
              if(nden.gt.0)then                                         1d19s23
               if(iterov.gt.nprescf)then                                11d1s24
                call lsqfit2(bc(ifmat),nden,iterm,bc(idata),nden,1,         5d19s22
     $        nden,bc(icoef),iterov,bc(iscr),bc(iwgt),0,rmsdiis,bc,ibc) 3d27s23
               end if                                                   11d1s24
              end if                                                    1d19s23
              if(iterov.gt.nprescf)then                                 11d1s24
               sum=0d0                                                     5d17s22
               do i=0,iterm-1                                              5d19s22
                sum=sum+bc(icoef+i)
               end do
               cn=1d0-sum
               bc(icoef+iterm)=cn                                          5d19s22
               if(nden.ne.0)then                                         11d29s22
                call dgemm('n','n',nden,1,iterov,1d0,bc(idiis1),nden,    11d29s22
     $           bc(icoef),iterov,0d0,bc(iden1d(1,1)),nden,             6d2s22
     d'gradcas.  2')
               end if                                                    11d29s22
              end if                                                    11d1s24
              jdiis3=idiis3+iterov*nden                                   5d19s22
              do i=0,nden-1                                               5d19s22
               bc(jdiis3+i)=bc(iden1d(1,1)+i)                             6d2s22
              end do                                                      5d19s22
              ibcoff=ifmat                                                5d17s22
             end if                                                       5d17s22
             if(nden.eq.0)then                                            5d20s22
              diff=0d0                                                    5d20s22
             else                                                         5d20s22
              diff=sqrt(diff/dfloat(nden))                                 5d17s22
             end if                                                       5d20s22
             ibcoff=ibctop                                                5d16s22
             call dws_bcast(diff,1)                                      7d14s22
             call dws_bcast(bc(i4odx(1)),n4o)                            7d14s22
             ratioo=diff/diffo                                          4d10s23
             diffo=diff                                                 4d10s23
             if(iterov.gt.3.and.ratioo.le.bestratio)then                11d1s24
              bestratio=ratioo                                          11d1s24
              ibestratio=iterov                                         11d1s24
             end if                                                     11d1s24
             if(ratioo.gt.5d0)then                                      11d1s24
              if(lprint)write(6,*)('DIIS iterations diverging')                   10d31s24
              if(idiisrest.eq.0)then                                    10d31s24
               if(lprint)write(6,*)('restart DIIS')                               10d31s24
               idiisrest=idiisrest+1                                    10d31s24
               bestratio=1d10                                           10d31s24
              else if(nprescf.eq.0)then                                 10d31s24
               if(lprint)
     $           write(6,*)('DIIS iterations diverging: switch to scf ')
               nprescf=50
               ibcoff=ibctop
              end if                                                    10d31s24
              if(bestratio.lt.1d0)then                                    10d31s24
               ibestm=ibestratio-1                                        10d31s24
               jdiis3=idiis3+ibestm*nden                                  10d31s24
               do i=0,nden-1                                              10d31s24
                bc(iden1d(1,1)+i)=bc(jdiis3+i)                            10d31s24
               end do                                                     10d31s24
              end if                                                      10d31s24
              go to 2122
             end if                                                     11d1s24
             if(nprescf.ne.0.and.ratioo.gt.0.99d0)then                    11d1s24
              if(lprint)write(6,*)
     $          ('stagnating scf iterations at iteration '),iterov,
     $          ('with last diff '),diff                                11d1s24
              go to 1215                                                11d1s24
             end if                                                       11d1s24
             if(iterov.gt.35.and.ratioo.gt.0.9d0)then                   4d10s23
              if(lprint)write(6,*)('stagnating diis at iteration '),    4d13s23
     $            iterov,(' with last diff '),diff                      4d13s23
              go to 1215                                                4d10s23
             end if                                                     4d10s23
             if(diff.gt.1d-9)then                                       4d27s23
              tolo=max(1d-12,min(tolo,diff*0.01d0))                     8d31s22
              tolv=max(1d-12,min(tolv,diff*0.01d0))                     8d31s22
              do i=0,n4o-1                                               6d16s22
               bc(i4odu(1)+i)=bc(i4odu(1)+i)+bc(i4odx(1)+i)              6d16s22
              end do                                                     6d16s22
              call buildcasgrad(ioooo,ionex,iamatu,idoub,iact,noc,nvirt,   5d9s22
     $   1,ih0mo,multh,idwsdeb,iden1d,jdend,nsdlk,isblk,                6d16s22
     $          nsdlk1,isblk1,iamatb,namatt,rmssz0,0,ipuse,nsblkderd,   7d1s22
     $          isblkderd,bc,ibc)                                       11d9s22
              if(idwsdeb.gt.10)then                                       5d19s22
               write(6,*)('adding density part ')
               call prntm2(bc(iamatu(1)),namatt,1,namatt)
               write(6,*)('to kinematic part')
               call prntm2(bc(iamatx(1)),namatt,1,namatt)
              end if                                                      5d19s22
              do ii=0,namatt-1                                            5d16s22
               bc(iamatu(1)+ii)=bc(iamatu(1)+ii)+bc(iamatx(1)+ii)         5d16s22
              end do                                                      5d16s22
              go to 22                                                    5d16s22
             end if                                                       5d16s22
             if(lprint)                                                 4d13s23
     $       write(6,*)('diis calculations converged after iteration '),4d13s23
     $            iterov                                                4d13s23
 1215        continue                                                   4d10s23
             if(lprint)then                                             4d7s23
              if(ipsym.eq.1)then                                        4d7s23
               write(6,*)('>energy gradients: ')                          7d28s22
               jdes=ides-1                                                  4d7s23
               do i=1,nstate                                                4d7s23
                nlab=0                                                      4d7s23
                do j=6,11                                                   4d7s23
                 if(isinfo(j,i).ne.0)nlab=j                                 4d7s23
                end do                                                      4d7s23
                stname(1)='      '                                           4d7s23
                do j=6,nlab                                                  4d7s23
                 jm=j-5                                                      4d7s23
                 stname(1)(jm:jm)=char(isinfo(j,i))                          4d7s23
                end do                                                       4d7s23
                stsym(1)='  '                                                 4d11s23
                if(nlzz.eq.2.and.stname(1)(1:1).ne.'S')then                   4d11s23
                 if(stname(1)(1:1).eq.'P')then                                4d11s23
                  if(isinfo(1,i).eq.2)then                                    4d11s23
                   stsym(1)='X '                                              4d11s23
                  else if(isinfo(1,i).eq.3)then                               4d11s23
                   stsym(1)='Y '                                              4d11s23
                  else if(isinfo(1,i).eq.6)then                               4d11s23
                   stsym(1)='XZ'                                              4d11s23
                  else if(isinfo(1,i).eq.7)then                               4d11s23
                   stsym(1)='YZ'                                              4d11s23
                  end if                                                      4d11s23
                 else                                                         4d11s23
                  if(isinfo(1,i).eq.5.or.                                     4d11s23
     $                 (isinfo(1,i).eq.1.and.nsymb.eq.4))then                 4d11s23
                   stsym(1)='Z '                                              4d11s23
                  else if(isinfo(1,i).eq.4)then                               4d11s23
                   stsym(1)='XY'                                              4d11s23
                  end if                                                      4d11s23
                 end if                                                       4d11s23
                end if                                                        4d11s23
                do iir=1,isinfo(4,i)                                        4d7s23
                 write(6,517)iir,isinfo(2,i),stname(1),stsym(1),        4d11s23
     $                bc(jdes+iir)                                      4d11s23
                end do                                                      4d7s23
                jdes=jdes+isinfo(4,i)                                       4d7s23
               end do                                                       4d7s23
              end if                                                    4d7s23
             end if                                                     7d28s22
             if(lbodc)then                                                6d3s22
              ibodc3=ibcoff                                             3d14s23
              icdir3=ibodc3+mden                                        3d14s23
              icdir4=icdir3+mden*2                                      4d10s23
              ibcoff=icdir4+ndentd*2                                    4d10s23
              call enough('gradcas.bodc3',bc,ibc)                       3d14s23
              do iz=ibodc3,ibcoff-1                                     3d14s23
               bc(iz)=0d0                                               3d14s23
              end do                                                    3d14s23
              call dws_gsumf(bc(icdc),ndentd)                           3d20s23
              call bodcpart(itrans,bc(iovrdk),bc(iovrdd),morb,iorb,      6d7s22
     $           idarot,                                                6d7s22
     $          nbasisp,nbasdws,idorel,bc(ibodc3),iden1d,idoub,iact,noc,3d14s23
     $        nsymb,ivecs,iden1i,ipuse,multh,isoua1,j2denbl,nsdlk,isblk,3d15s23
     $            isou,idwsdeb,lprint,bc,ibc,mden,nstate,isinfo,        3d16s23
     $             bc(icdir3),bc(icdc),tracer,bc(icdir4),ndentd)        4d10s23
              if(mynowprog.eq.0)then                                     7d14s22
               if(jbodc2.eq.0)then                                        7d7s22
                jbodcu=jbodc                                            3d15s23
                do k=0,mden-1
                 bc(jbodcu)=bc(jbodcu)+bc(ibodc3+k)                      3d15s23
                 jbodcu=jbodcu+natom                                     3d15s23
                end do                                                  3d15s23
               else                                                       7d7s22
                jbodcu=jbodc                                            3d15s23
                jbodc2u=jbodc2                                          3d15s23
                do k=0,mden-1                                           3d15s23
                 bodc3=bc(ibodc3+k)*0.25d0                              3d15s23
                 bc(jbodcu)=bc(jbodcu)+bodc3                            3d15s23
                 bc(jbodc2u)=bc(jbodc2u)+bodc3                               7d7s22
                 jbodcu=jbodcu+natom                                    3d21s23
                 jbodc2u=jbodc2u+natom                                  3d21s23
                end do                                                  3d15s23
               end if                                                     7d7s22
               jcdir3=icdir3                                            3d22s23
               do i=1,nstate                                            3d22s23
                nlab=0                                                          12d28s19
                do j=6,11                                                       12d31s19
                 if(isinfo(j,i).ne.0)nlab=j                                     12d28s19
                end do                                                          12d28s19
                stname(1)='      '                                      3d22s23
                do j=6,nlab                                             3d22s23
                 jm=j-5                                                 3d22s23
                 stname(1)(jm:jm)=char(isinfo(j,i))                     3d22s23
                end do                                                  3d22s23
                stsym(1)='  '                                                 4d11s23
                if(nlzz.eq.2.and.stname(1)(1:1).ne.'S')then                   4d11s23
                 if(stname(1)(1:1).eq.'P')then                                4d11s23
                  if(isinfo(1,i).eq.2)then                                    4d11s23
                   stsym(1)='X '                                              4d11s23
                  else if(isinfo(1,i).eq.3)then                               4d11s23
                   stsym(1)='Y '                                              4d11s23
                  else if(isinfo(1,i).eq.6)then                               4d11s23
                   stsym(1)='XZ'                                              4d11s23
                  else if(isinfo(1,i).eq.7)then                               4d11s23
                   stsym(1)='YZ'                                              4d11s23
                  end if                                                      4d11s23
                 else                                                         4d11s23
                  if(isinfo(1,i).eq.5.or.                                     4d11s23
     $           (isinfo(1,i).eq.1.and.nsymb.eq.4))then                 4d11s23
                   stsym(1)='Z '                                              4d11s23
                  else if(isinfo(1,i).eq.4)then                               4d11s23
                   stsym(1)='XY'                                              4d11s23
                  end if                                                      4d11s23
                 end if                                                       4d11s23
                end if                                                        4d11s23
                do ir1=1,isinfo(4,i)                                    3d22s23
                 do ir2=1,ir1                                           3d22s23
                  jcdir3p=jcdir3+mden                                   4d10s23
                  sum=bc(jcdir3)+bc(jcdir3p)                            4d10s23
                  if(abs(sum).gt.1d-10)then                             4d10s23
                   write(6,57)ir1,isinfo(2,i),stname(1),stsym(1),opname,4d17s23
     $                 ir2,isinfo(2,i),stname(1),stsym(1),-sum,         4d17s23
     $                 -bc(jcdir3),-bc(jcdir3p)                         4d11s23
   57              format('<',2i2,a6,a2,'|',a9,2i2,a6,a2,'> = ',es15.7, 4d11s23
     $                  ' = ',es15.7,' and ',es15.7)                    4d10s23
                  end if                                                3d22s23
                  jcdir3=jcdir3+1                                       3d22s23
                 end do                                                 3d22s23
                end do                                                  3d22s23
               end do                                                   3d22s23
               if(ndentd.gt.0)then
                jcdir4=icdir4                                           3d22s23
                do i=1,nstate-1                                         3d22s23
                 nlab=0                                                          12d28s19
                 do j=6,11                                                       12d31s19
                  if(isinfo(j,i).ne.0)nlab=j                                     12d28s19
                 end do                                                          12d28s19
                 stname(1)='      '                                      3d22s23
                 do j=6,nlab                                             3d22s23
                  jm=j-5                                                 3d22s23
                  stname(1)(jm:jm)=char(isinfo(j,i))                     3d22s23
                 end do                                                  3d22s23
                 stsym(1)='  '                                                 4d11s23
                 if(nlzz.eq.2.and.stname(1)(1:1).ne.'S')then                   4d11s23
                  if(stname(1)(1:1).eq.'P')then                                4d11s23
                   if(isinfo(1,i).eq.2)then                                    4d11s23
                    stsym(1)='X '                                              4d11s23
                   else if(isinfo(1,i).eq.3)then                               4d11s23
                    stsym(1)='Y '                                              4d11s23
                   else if(isinfo(1,i).eq.6)then                               4d11s23
                    stsym(1)='XZ'                                              4d11s23
                   else if(isinfo(1,i).eq.7)then                               4d11s23
                    stsym(1)='YZ'                                              4d11s23
                   end if                                                      4d11s23
                  else                                                         4d11s23
                   if(isinfo(1,i).eq.5.or.                                     4d11s23
     $           (isinfo(1,i).eq.1.and.nsymb.eq.4))then                 4d11s23
                    stsym(1)='Z '                                              4d11s23
                   else if(isinfo(1,i).eq.4)then                               4d11s23
                    stsym(1)='XY'                                              4d11s23
                   end if                                                      4d11s23
                  end if                                                       4d11s23
                 end if                                                        4d11s23
                 do ip=i+1,nstate                                       3d22s23
                  if(isinfo(2,ip).eq.isinfo(2,i).and.isinfo(3,ip).eq.   3d31s23
     $                 isinfo(3,i))then                                 4d7s23
                   if(multh(isinfo(1,i),isinfo(1,ip)).eq.ipuse)then     4d7s23
                    nlab=0                                                          12d28s19
                    do j=6,11                                                       12d31s19
                     if(isinfo(j,ip).ne.0)nlab=j                                     12d28s19
                    end do                                                          12d28s19
                    stname(2)='      '                                      3d22s23
                    do j=6,nlab                                             3d22s23
                     jm=j-5                                                 3d22s23
                     stname(2)(jm:jm)=char(isinfo(j,ip))                     3d22s23
                    end do                                                  3d22s23
                    stsym(2)='  '                                                 4d11s23
                    if(nlzz.eq.2.and.stname(2)(1:1).ne.'S')then                   4d11s23
                     if(stname(2)(1:1).eq.'P')then                                4d11s23
                      if(isinfo(1,ip).eq.2)then                                    4d11s23
                       stsym(2)='X '                                              4d11s23
                      else if(isinfo(1,ip).eq.3)then                               4d11s23
                       stsym(2)='Y '                                              4d11s23
                      else if(isinfo(1,ip).eq.6)then                               4d11s23
                       stsym(2)='XZ'                                              4d11s23
                      else if(isinfo(1,ip).eq.7)then                               4d11s23
                       stsym(2)='YZ'                                              4d11s23
                      end if                                                      4d11s23
                     else                                                         4d11s23
                      if(isinfo(1,ip).eq.5.or.                                     4d11s23
     $           (isinfo(1,ip).eq.1.and.nsymb.eq.4))then                 4d11s23
                       stsym(2)='Z '                                              4d11s23
                      else if(isinfo(1,ip).eq.4)then                               4d11s23
                       stsym(2)='XY'                                              4d11s23
                      end if                                                      4d11s23
                     end if                                                       4d11s23
                    end if                                                        4d11s23
                    do irp=1,isinfo(4,ip)                                 3d22s23
                     do ir=1,isinfo(4,i)                                  3d22s23
                      jcdir4p=jcdir4+ndentd                             4d10s23
                      sum=bc(jcdir4)+bc(jcdir4p)                        4d10s23
                      if(abs(sum).gt.1d-10)write(6,57)ir,               4d14s23
     $                   isinfo(2,i),stname(1),stsym(1),opname,irp,     4d11s23
     $                     isinfo(2,ip),stname(2),stsym(2),-sum,        4d11s23
     $                     -bc(jcdir4),-bc(jcdir4p)                     4d11s23
                      jcdir4=jcdir4+1                                     3d22s23
                     end do                                               3d22s23
                    end do                                                3d22s23
                   else                                                 4d7s23
                    jcdir4=jcdir4+isinfo(4,ip)*isinfo(4,i)              4d7s23
                   end if                                                3d31s23
                  end if                                                4d7s23
                 end do                                                 3d22s23
                end do                                                  3d22s23
               end if
              end if                                                     7d14s22
              ibcoff=ibodc3                                             3d21s23
             end if                                                       6d3s22
             if(idwsdeb.ne.0)then                                        7d14s22
              do isb=1,nsymb                                               5d16s22
               jsb=multh(isb,ipuse)                                       6d20s22
               if(iact(isb).gt.0)then
                write(6,*)('d of density for sym '),isb
                call prntm2(bc(iden1d(isb,1)),iact(isb),iact(jsb),        6d20s22
     $            iact(isb))
                if(lbodc)then                                              6d2s22
                 write(6,*)('vdv density ')
                 call prntm2(bc(iden1d(isb,2)),iact(isb),iact(jsb),       6d20s22
     $             iact(isb))
                 write(6,*)('dvv density ')
                 call prntm2(bc(iden1d(isb,3)),iact(isb),iact(jsb),       6d20s22
     $             iact(isb))
                 write(6,*)('dvdv density ')
                 call prntm2(bc(iden1d(isb,4)),iact(isb),iact(isb),
     $             iact(isb))
                end if                                                     6d2s22
               end if
              end do                                                       5d16s22
             end if                                                      7d14s22
             if(inocan.eq.0.and.ipuse.eq.1)then                                          5d19s22
              ibcb4=ibcoff                                                5d20s22
              igoalqxz=jmatd(1)
              call parajkfromhd0(noc,multh,jmats,kmats,ioooo,ionex,     7d25s22
     $             itprt,                                               7d25s22
     $            idwsdeb,nvirt,ipuse,idarot,i4odu,ionexdu,nbasdws,     5d20s22
     $            isblkder,isblkxder,nsblkder,nsblkxder,isblkkder,      5d20s22
     $            nsblkkder,jmatdu,kmatdu,i3x,i4od2b,ionexd2,isblkxder1,5d20s22
     $            0,isblkder1,0,igoalqxz,bc,ibc)                        12d19s22
              call sym4o(i4odu,noc,isblkder,isblkxder,nsblkder,         7d25s22
     $             nsblkxder,bc,ibc)                                    11d10s22
              do i=0,n4o-1                                                 5d16s22
               bc(i4odx(1)+i)=bc(i4odu(1)+i)+bc(i4odx(1)+i)               5d26s22
              end do                                                       5d16s22
              call trianglez(i4odx,ionexd,noc,nvirt,nsblkder,isblkder,   5d26s22
     $          0,isblkxder,idwsdeb,ipuse.eq.1,bc,ibc)                  11d10s22
              do isb=1,nsblkkder                                          5d20s22
               nrow=noc(isblkkder(1,isb))*noc(isblkkder(2,isb))           5d20s22
               ncol=nbasdws(isblkkder(3,isb))*nbasdws(isblkkder(4,isb))   5d20s22
               nall=nrow*ncol                                             5d20s22
               if(nall.gt.0)then                                          5d20s22
                call ilimts(nvirt(isblkkder(3,isb)),                      5d20s22
     $            nvirt(isblkkder(4,isb)),mynprocg,mynowprog,il,ih,i1s, 5d20s22
     $            i1e,is,i2e)                                           5d20s22
                nhere=ih+1-il                                             5d20s22
                nn=nhere*nrow-1                                           5d20s22
                do i=0,nn                                                 5d20s22
                 bc(jmatd(isb)+i)=bc(jmatd(isb)+i)+bc(jmatdu(isb)+i)      5d20s22
                 bc(kmatd(isb)+i)=bc(kmatd(isb)+i)+bc(kmatdu(isb)+i)      5d20s22
                end do                                                    5d20s22
               end if                                                     5d20s22
              end do                                                      5d20s22
              ibcoff=ibcb4                                                5d20s22
              call rordr(i4odx,noc,nsblkder,isblkder,i4oduu,bc,ibc)     11d10s22
              call dcannon(ih0mo,ih0da,mynprocg,mynowprog,ioooo,i4oduu,   5d26s22
     $           kmats,kmatd,jmats,jmatd,noc,morb,idarot,idwsdeb,idoub, 5d20s22
     $           iact,nvirt,iden1,iden1d,nsblkkder,isblkkder,bc,ibc)    12d19s22
             end if                                                       5d19s22
             if(mynowprog.eq.0)then                                     7d28s22
              do isb=1,nsymb                                             7d28s22
               if(nbasdws(isb).gt.0)then                                 7d28s22
                nn=nbasdws(isb)**2                                        7d28s22
                write(1)(bc(idarot(isb)+i),i=0,nn-1)                    12d19s22
                write(1)(bc(itrans(isb)+i),i=0,nn-1)                    8d15s23
               end if                                                    7d28s22
              end do                                                     7d28s22
             end if                                                     7d28s22
            end if                                                      7d25s22
c
c
            ibcoff=ih0d                                                  5d4s22
           end if                                                        5d4s22
           if(ipass.lt.npass)then                                       1d17s23
            ipuse=multh(ipuse,iapair(2,ia))                               2d4s16
           end if                                                       1d17s23
          end do                                                         2d3s16
         end if                                                          5d4s22
        end do                                                           5d4s22
       end do
      end do                                                            6d7s22
      if(lbodc.and.lprint)then                                          7d25s22
       if(nlzz.eq.2)then                                                4d7s23
        write(6,*)('>total radial bodc: ')                              4d7s23
       else                                                             4d7s23
        write(6,*)('>total bodc: ')                                       6d7s22
       end if                                                           4d7s23
       moff=0                                                           3d23s23
       do is=1,nstate                                                   3d23s23
        nlab=0                                                          12d28s19
        do j=6,11                                                       12d31s19
         if(isinfo(j,is).ne.0)nlab=j                                     12d28s19
        end do                                                          12d28s19
        stname(1)='      '                                              3d23s23
        do j=6,nlab                                                     3d23s23
         jm=j-5                                                         3d23s23
         stname(1)(jm:jm)=char(isinfo(j,is))                            3d23s23
        end do                                                          3d23s23
        do ir1=1,isinfo(4,is)                                           3d23s23
         do ir2=1,ir1                                                   3d23s23
          write(6,557)ir2,ir1,isinfo(2,is),stname(1)                    3d23s23
  557     format('>states ',2i2,' of ',i2,a6)                           3d23s23
          do ia=1,natom                                                    6d7s22
           jbodc=ibodc+ia-1+natom*moff                                  3d23s23
           write(6,*)('>for nucleus '),ia,bc(jbodc)                      4d19s23
          end do                                                           6d7s22
          moff=moff+1                                                   3d23s23
         end do                                                         3d23s23
        end do
       end do                                                           3d23s23
      end if                                                            6d7s22
      if(mynowprog.eq.0)then                                            10d7s24
       call second(timen)                                               10d7s24
       telap=timen-time0                                                10d7s24
       write(6,*)('total cpu time on this proc for gradcas is '),telap  10d7s24
      end if                                                            10d7s24
      return                                                            5d4s22
      end                                                               5d4s22
