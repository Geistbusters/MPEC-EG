c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine stepwisecsf(nopenx,s,xms,vec,iaorb,iaorb1,iflag,idata, 11d9s22
     $     bc,ibc)                                                      11d9s22
      implicit real*8 (a-h,o-z)
c
c     stepwise determination of csf to det transformation via vector
c     coupling coefficients
c     iflag=0: replicate across procs                                   5d17s21
c     iflag=1: for largest no. of open shells, distribute across procs  5d17s21
c
      parameter (idos=60)                                                 5d28s19
      integer*8 ipower2(idos),itry8,iaorb(*),itmp8                      7s29s19
      integer*2 itmp2(4)                                                7s29s19
      integer*1 iaorb1(*)                                               6d13s19
      logical ldebug                                                    12d12s19
      equivalence (itmp2,itmp8)                                         7s29s19
      dimension vec(*),idata(*)                                         7d30s19
      include "common.print"                                            1d3s20
      save
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      2d1s21
     $     mynnode                                                      2d1s21
      if(iprtr(1).eq.0)then                                             1d3s20
       ldebug=.false.                                                    12d12s19
      else                                                              1d3s20
       ldebug=.true.
      end if                                                            1d3s20
      if(ldebug)write(6,*)('in stepwisecsf for nopenx = '),nopenx,
     $     (' spin '),s,(' ms '),xms
      ipower2(1)=2                                                      5d28s19
      do i=2,idos                                                       5d28s19
       im=i-1                                                           5d28s19
       ipower2(i)=ipower2(im)*2                                         5d28s19
      end do                                                            5d28s19
      isf=nint(s*2d0)
      msf=nint(xms*2d0)
      ivecpt=ibcoff
      ibcoff=ivecpt+nopenx+1
      ivecpp=ibcoff                                                     7s29s19
      ibcoff=ivecpp+nopenx+1                                            7s29s19
      call enough('stepwisecsf.  1',bc,ibc)
c
c     memory map:
c     ibc(ivecpt+io) points to memory info for io open shells
c
c     all subsequent addresses are relative to the base address.
c
c     for each no. of open shells
c     ibc(ims1=ibc(ivecpt+io)),... ibc(ims1+2*io) points to
c     memory info for 2ms= -io,-io+1,...,+io, respectively.
c     for each ms,
c     1st word is no. of dets,
c     then one word per det for occupation code,
c     then pointers for s=iabs(ms),iabs(ms)+1,...,io/2.
c     for each s,
c     first word is no. fcns followed by vectors.
c
c     zero open shells
c
      ibctop=ibcoff
      ibc(ivecpt)=ibcoff
      ibase=ibcoff
      ibc(ibcoff)=1
      ibc(ibcoff+1)=2
      ibc(ibcoff+2)=1
      ibc(ibcoff+3)=0
      ibc(ibcoff+4)=5
      ibc(ibcoff+5)=1
      bc(ibcoff+6)=1d0
      ibc(ibcoff+7)=0                                                   7s29s19
      ibcoff=ibcoff+8                                                   7s29s19
      itmpj=ibcoff
      ibcoff=itmpj+nopenx
      call enough('stepwisecsf.  2',bc,ibc)
      do io=0,2
c
c     add one electron at a time
c
       if(ldebug)write(6,*)('for io = '),io
       itryx=isf+nopenx-io
       itryn=isf+io-nopenx
       ksx=min(io,itryx)
       ksn=mod(io,2)                                                    1d2s20
       if(ldebug)write(6,*)('ksn,ksx '),ksn,ksx
       iop=io+1
       if(ldebug)then
        write(6,*)('for iop= '),iop
        write(6,*)(' ')
        write(6,*)(' ')
       end if
       itryx=isf+nopenx-iop
       itryn=isf+iop-nopenx
       jsx=min(iop,itryx)
       jsn=mod(iop,2)                                                   1d2s20
       if(ldebug)write(6,*)('jsn,jsx '),jsn,jsx
       ibc(ivecpt+iop)=ibcoff
       ibase=ibcoff
       if(ldebug)write(6,*)('base pointer is '),ibase
       msjs=ibase+iop+1
       ibcoff=ibase+2*iop+2
       call enough('stepwisecsf.  3',bc,ibc)
       do msj=-jsx,jsx,2
        if(ldebug)then                                                  12d13s19
         write(6,*)('for msj = '),msj
         write(6,*)(' ')
        end if                                                          12d13s19
        ndetx=0
        do js=max(iabs(msj),jsn),jsx,2
         if(ldebug)write(6,*)('for js = '),js
         njs=0
         isx=min(js+1,ksx)
         isn=max(iabs(js-1),ksn)
         do is=isn,isx,2
          ndetn=0
          nroot=0
          do igamma=-1,1,2
           msi=msj-igamma
           if(iabs(msi).le.is)then
            ibaseo=ibc(ivecpt+io)
            ipms=ibaseo+msi+io+1
            indp=ibc(ipms)                                              6d12s19
            indp=indp+ibaseo
            ndeth=ibc(indp)
            isp=indp+ndeth+is+1
            nroota=ibc(isp)
            nroot=ibc(nroota+ibaseo)
            if(ldebug)write(6,*)('for is = '),is,msi,ndeth,nroot,nroota,
     $           nroota+ibaseo
            ndetn=ndetn+ndeth
           end if
          end do
          njs=njs+nroot
          if(ldebug)write(6,*)('ndetn = '),ndetn
          ndetx=max(ndetx,ndetn)
         end do
         if(ldebug)write(6,*)('number of roots for this js'),njs,js     12d13s19
         if(js.gt.nopenx)then
          write(6,*)('js exceeds nopenx!! '),js,nopenx
          call dws_sync
          call dws_finalize
          stop
         end if
         ibc(itmpj+js)=njs
        end do
        ibc(msjs+msj)=ibcoff-ibase
        ibc(ibcoff)=ndetx
        if(ldebug)then
        write(6,*)('next '),ibase,msjs+msj,msjs+msj-ibase,ibcoff,
     $       ibcoff-ibase
        write(6,*)('storing ndetn under '),msjs+msj,ibcoff
        end if
        jdet=ibcoff
        ibcoff=ibcoff+1
        ibcoff=ibcoff+ndetx
        jptrn=ibcoff
        ibcoff=jptrn+jsx+1                                              1d3s20
        call enough('stepwisecsf.  4',bc,ibc)
        do js=max(iabs(msj),jsn),jsx,2
         kptr=jptrn+js
         ibc(kptr)=ibcoff-ibase
         ibc(ibcoff)=ibc(itmpj+js)
         if(ldebug)write(6,*)('js no. roots stored at '),js,ibcoff-ibase
     $        ,ibc(ibcoff),ibcoff,nopenx,jptrn
         ibcoff=ibcoff+1
         if(ldebug)write(6,*)('vectors stored at '),ibcoff,ibcoff-ibase
         ibcoff=ibcoff+ibc(itmpj+js)*ndetx
         ibcoff=ibcoff+ibc(itmpj+js)                                    7s29s19
        end do
        call enough('stepwisecsf.  5',bc,ibc)
        ndj=0
        do js=max(iabs(msj),jsn),jsx,2
         if(ldebug)write(6,*)('for js = '),js
         kptr=jptrn+js
         nrootja=ibc(kptr)
         nrootj=ibc(nrootja+ibase)
         if(ldebug)write(6,*)('no. roots = '),nrootj,kptr,nrootja+ibase
         ivecn=nrootja+ibase+1
         ilabn=ivecn+nrootj*ndetx                                       7s29s19
         if(ldebug)write(6,*)('new vectors address: '),ivecn,ivecn-ibase
         do i=0,nrootj*ndetx-1
          bc(ivecn+i)=0d0
         end do
         isx=min(js+1,ksx)
         isn=max(iabs(js-1),ksn)
         nj=0
         do is=isn,isx,2
          nroot=0
          do igamma=-1,1,2
           msi=msj-igamma
           if(iabs(msi).le.is)then
            if(ldebug)write(6,*)('for is,msi: '),is,msi
            ibaseo=ibc(ivecpt+io)
            ipms=ibaseo+msi+io+1
            indp=ibc(ipms)                                              6d12s19
            indp=indp+ibaseo
            ndeth=ibc(indp)
            jptr=indp+ndeth
            if(ldebug)write(6,*)('det '),indp,ndeth,jptr
            kptr=jptr+is+1
            nroota=ibc(kptr)
            nroot=ibc(nroota+ibaseo)
            ivec=nroota+ibaseo+1
            ilab=ivec+ndeth*nroot                                       7s29s19
            if(ldebug)then
             write(6,*)('vectors here '),ivec
             call prntm2(bc(ivec),ndeth,nroot,ndeth)
            end if
            iphs=is-1-msj
            iphs=iabs(iphs/2)
            fact=f3j(1,is,js,igamma,msi,-msj,1)
            fact=fact*sqrt(dfloat(js+1))
            if(mod(iphs,2).ne.0)fact=-fact
            if(ldebug)write(6,*)('is, msi '),is,msi,fact
            do j=0,nroot-1                                              7s29s19
             itmp8=ibc(ilab+j)                                          7s29s19
             itmp2(1)=itmp2(2)                                          7s29s19
             itmp2(2)=itmp2(3)                                          7s29s19
             itmp2(3)=itmp2(4)                                          7s29s19
             itmp2(4)=is                                                7s29s19
             iad=ilabn+nj+j                                             7s29s19
             ibc(iad)=itmp8                                             7s29s19
            end do                                                      7s29s19
            do i=1,ndeth
             if(ldebug)write(6,*)('det '),i
             itry8=ibc(indp+i)
             if(igamma.lt.0)itry8=itry8+ipower2(iop)
             do n=1,ndj
              if(itry8.eq.ibc(jdet+n))then
               inrow=n
               go to 1
              end if
             end do
             ndj=ndj+1
             inrow=ndj
             ibc(jdet+ndj)=itry8
    1        continue
             if(ldebug)write(6,*)('matches det '),inrow
             do j=1,nroot
              iadvo=ivec+i-1+ndeth*(j-1)
              iadvn=ivecn+inrow-1+ndetx*(nj+j-1)
              bc(iadvn)=bc(iadvn)+fact*bc(iadvo)
             end do
            end do
            if(ldebug)write(6,*)('for is = '),is,msi,ndeth
           end if
          end do
          nj=nj+nroot
         end do
         if(nj.ne.nrootj)then
          write(6,*)('nj vs nrootj '),nj,nrootj
          call dws_sync
          call dws_finalize
          stop
         end if
        end do
        call prtsigs(ibc(jdet+1),ndj,iop,ibc(ibcoff))
        if(ndj.ne.ndetx)then
         call dws_sync
         call dws_finalize
         stop
        end if
       end do
      end do
c
c     plan b: couple electrons 2 at a time
c
      io=2
      ibaseo=ibc(ivecpt+io)
      itryx=isf+nopenx-io
      jsx=min(io,itryx)
      do is=0,jsx,2                                                       7s29s19
       isp=is+1
       if(ldebug)write(6,*)('for '),isp,('let coupled pair: ')                    7s29s19
       do ms=-is,is,2                                                   7s29s19
        msh=ms/2                                                        7s29s19
        if(ldebug)write(6,*)('for Ms = '),msh                                     7s29s19
        ipms=ibaseo+ms+io+1
        indp=ibc(ipms)+ibaseo
        ndeth=ibc(indp)
        jptr=indp+ndeth
        kptr=jptr+is+1
        nroota=ibc(kptr)
        nroot=ibc(nroota+ibaseo)
        ivec=nroota+ibaseo+1
        do i=1,ndeth
         itry8=ibc(indp+i)                                              7s29s19
         iad=ivec+i-1
         if(ldebug)write(6,*)itry8,(bc(iad+ndeth*j),j=0,nroot-1)
         do j=1,2
          itry8=itry8/2
          itry4=itry8
          if(mod(itry4,2).ne.0)then
           if(ldebug)write(6,*)('-')
          else
           if(ldebug)write(6,*)('+')
          end if
         end do
        end do                                                          7s29s19
       end do                                                           7s29s19
      end do                                                            7s29s19
      iol=2-mod(nopenx,2)                                               7d30s19
c     for nopenx=5, need nopenx rather than nopenx-1
      do io=iol,nopenx-1,2                                              12d12s19
       if(ldebug)write(6,*)('for '),io,('electrons ')
       itryx=isf+nopenx-io
       itryn=isf+io-nopenx
       if(iflag.eq.0)then                                               6d3s21
        itryn=0                                                         6d3s21
       end if                                                           6d3s21
       if(ldebug)write(6,*)('itryn,itryx '),itryn,itryx
       ksx=min(io,itryx)
       ksn=max(mod(io,2),itryn)
       if(ldebug)write(6,*)('ksn,ksx '),ksn,ksx
       iop=io+1
       iopp=iop+1
       itryx=isf+nopenx-iopp
       itryn=isf+iopp-nopenx
       if(ldebug)write(6,*)('itryn,itryx '),itryn,itryx
       jsx=min(iopp,itryx)
       jsn=mod(iopp,2)                                                  3d31s20
       if(ldebug)write(6,*)('jsn,jsx '),jsn,jsx
       ibc(ivecpp+iopp)=ibcoff
       ibase=ibcoff
       if(ldebug)write(6,*)('base pointer is '),ibase,(' for iopp = '),
     $      iopp
       msjs=ibase+iopp+1
       ibcoff=ibase+2*iopp+2
       call enough('stepwisecsf.  6',bc,ibc)
       do msj=-jsx,jsx,2
        if(ldebug)write(6,*)('for msj = '),msj
        ndetx=0
        do js=max(iabs(msj),jsn),jsx,2
         njs=0
         do ipair=0,2,2
          isx=min(js+ipair,ksx)
          isn=max(iabs(js-ipair),ksn)
          do is=isn,isx,2
           ndetn=0
           nroot=0
           do igamma=-ipair,ipair,2
            msi=msj-igamma
            if(iabs(msi).le.is)then
             if(io.le.2)then                                            7s29s19
              ibaseo=ibc(ivecpt+io)                                     7s29s19
             else                                                       7s29s19
              ibaseo=ibc(ivecpp+io)
             end if                                                     7s29s19
             ipms=ibaseo+msi+io+1
             if(ipms-ioffsx.lt.0.or.ipms-ioffsx.gt.maxbc)then           10d31s22
              write(6,*)('ipms = '),io,msi,ibaseo,ipms
              call dws_sync
              call dws_finalize
              stop 'stepwise1'
             end if
             indp=ibc(ipms)                                              6d12s19
             indp=indp+ibaseo
             if(indp-ioffsx.lt.0.or.indp-ioffsx.gt.maxbc)then           10d31s22
              write(6,*)('indp '),indp,ibaseo,ibc(ipms)
              call dws_sync
              call dws_finalize
              stop 'stepwise2'
             end if
             ndeth=ibc(indp)
             isp=indp+ndeth+is+1
             if(isp-ioffsx.lt.0.or.isp-ioffsx.gt.maxbc)then             10d31s22
              write(6,*)('isp '),isp
              call dws_sync
              call dws_finalize
              stop 'stepwise3'
             end if
             nroota=ibc(isp)
             if(nroota.lt.0.or.nroota.gt.maxbc)then                     10d31s22
              write(6,*)('nroota '),nroota                              10d31s22
              call dws_sync
              call dws_finalize
              stop 'stepwise4'
             end if
             nroot=ibc(nroota+ibaseo)
             if(igamma.eq.0)ndeth=ndeth*2                               7s29s19
             ndetn=ndetn+ndeth
            end if
           end do
           njs=njs+nroot
           ndetx=max(ndetx,ndetn)
          end do                                                        7s29s19
         end do
         if(js.gt.nopenx)then
          write(6,*)('js exceeds nopenx!! '),js,nopenx
          call dws_sync
          call dws_finalize
          stop 'stepwise5'
         end if
         ibc(itmpj+js)=njs
        end do
        ibc(msjs+msj)=ibcoff-ibase
        if(io.ge.nopenx-2.and.iflag.ne.0)then                           5d17s21
         call ilimts(1,ndetx,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)  2d1s21
         nhere=ih+1-il                                                  2d1s21
        else                                                            2d1s21
         il=1                                                           2d1s21
         ih=ndetx                                                       2d1s21
         nhere=ndetx                                                    2d1s21
        end if                                                          2d1s21
        ibc(ibcoff)=ndetx
        jdet=ibcoff
        ibcoff=ibcoff+1
        ibcoff=ibcoff+ndetx
        jptrn=ibcoff
        ibcoff=jptrn+jsx+1                                              1d3s20
        do js=max(iabs(msj),jsn),jsx,2
         kptr=jptrn+js
         ibc(kptr)=ibcoff-ibase
         ibc(ibcoff)=ibc(itmpj+js)
         ibcoff=ibcoff+1
         ibcoff=ibcoff+ibc(itmpj+js)*nhere                              2d1s21
         ibcoff=ibcoff+ibc(itmpj+js)                                    7s29s19
        end do
        call enough('stepwisecsf.  7',bc,ibc)
        ndj=0
        do js=max(iabs(msj),jsn),jsx,2
         if(ldebug)write(6,*)('for js = '),js
         kptr=jptrn+js
         nrootja=ibc(kptr)
         nrootj=ibc(nrootja+ibase)
         ivecn=nrootja+ibase+1
         ilabn=ivecn+nrootj*nhere                                       2d1s21
         do i=0,nrootj*nhere-1                                          2d1s21
          bc(ivecn+i)=0d0
         end do
         nj=0
         do ipair=0,2,2                                                 7s29s19
          isx=min(js+ipair,ksx)                                         7s29s19
          isn=max(iabs(js-ipair),ksn)                                   7s29s19
          do is=isn,isx,2
           nroot=0
           do igamma=-ipair,ipair,2                                     7s29s19
            msi=msj-igamma
            if(iabs(msi).le.is)then
             ibasep=ibc(ivecpt+2)                                       7s29s19
             ipmsp=ibasep+igamma+3                                      7s29s19
             indpp=ibc(ipmsp)+ibasep                                    7s29s19
             ndetp=ibc(indpp)                                           7s29s19
             jptrp=indpp+ndetp                                          7s29s19
             kptrp=jptrp+ipair+1                                        7s29s19
             nrootp=ibc(kptrp)                                          7s29s19
             ivecp=nrootp+ibasep+1                                      7s29s19
             if(io.le.2)then                                            7s29s19
              ibaseo=ibc(ivecpt+io)                                     7s29s19
             else                                                       7s29s19
              ibaseo=ibc(ivecpp+io)
             end if                                                     7s29s19
             ipms=ibaseo+msi+io+1
             indp=ibc(ipms)                                              6d12s19
             indp=indp+ibaseo
             ndeth=ibc(indp)
             jptr=indp+ndeth
             kptr=jptr+is+1
             nroota=ibc(kptr)
             nroot=ibc(nroota+ibaseo)
             ivec=nroota+ibaseo+1
             ilab=ivec+ndeth*nroot                                       7s29s19
             iphs=is-ipair-msj
             iphs=iabs(iphs/2)
             fact=f3j(ipair,is,js,igamma,msi,-msj,1)
             fact=fact*sqrt(dfloat(js+1))
             if(mod(iphs,2).ne.0)fact=-fact
             jlabn=ilabn+nj                                             7s29s19
             do j=0,nroot-1                                              7s29s19
              itmp2(1)=is                                               7s29s19
              itmp2(2)=ipair                                            7s29s19
              ibc(jlabn+j)=itmp8                                        7s29s19
             end do                                                      7s29s19
             do j=1,ndetp                                               7s29s19
              if(igamma.eq.0)then                                       7s29s19
               if(j.eq.1)then                                           7s29s19
                itmp8=ipower2(iopp)                                       7s29s19
               else                                                     7s29s19
                itmp8=ipower2(iop)                                        7s29s19
               end if                                                   7s29s19
              else                                                      7s29s19
               if(igamma.gt.0)then                                      7s29s19
                itmp8=0                                                 7s29s19
               else                                                     7s29s19
                itmp8=ipower2(iop)+ipower2(iopp)                            7s29s19
               end if                                                   7s29s19
              end if                                                    7s29s19
              ff=fact*bc(ivecp+j-1)                                     7s29s19
              do i=1,ndeth
               itry8=ibc(indp+i)+itmp8                                  7s29s19
               do n=1,ndj
                if(itry8.eq.ibc(jdet+n))then
                 inrow=n
                 go to 101
                end if
               end do
               ndj=ndj+1
               inrow=ndj
               ibc(jdet+ndj)=itry8
  101          continue
               if(inrow.ge.il.and.inrow.le.ih)then                      2d1s21
                do k=1,nroot
                 iadvo=ivec+i-1+ndeth*(k-1)
                 iadvn=ivecn+inrow-il+nhere*(nj+k-1)                    2d1s21
                 bc(iadvn)=bc(iadvn)+ff*bc(iadvo)
                end do
               end if                                                   2d1s21
              end do
             end do                                                     7s29s19
            end if
           end do
           nj=nj+nroot
          end do
         end do
        end do                                                          7s29s19
        if(nj.ne.nrootj)then
         write(6,*)('nj vs nrootj '),nj,nrootj
         call dws_sync
         call dws_finalize
         stop
        end if
       end do                                                           7s29s19
      end do                                                            7s29s19
      return
      entry fetchcsf2(nopenin,is2,nrooto,vec,iaorb,iaorb1,ms2,icall,    11d15s22
     $     bc,ibc)                                                      11d15s22
      if(ldebug)then                                                    12d12s19
       write(6,*)('in fetchcsf2 ')
       write(6,*)('no. open shells = '),nopenin
       write(6,*)('twice desired spin: '),is2                            9d12s19
       write(6,*)('twice desired ms: '),ms2
      end if                                                            12d12s19
      nrooto=0                                                          9d12s19
      if(ldebug.and.ms2.gt.is2)write(6,*)('bailing for ms2 = '),ms2,
     $     (' is gt is2 = '),is2
      if(ldebug.and.max(is2,ms2).gt.nopenin)write(6,*)('bailing for '),
     $     ('max is2,ms2 = '),max(is2,ms2),(' gt nopenin = '),nopenin
      if(ms2.gt.is2.or.max(is2,ms2).gt.nopenin)return                   9d12s19
      if(nopenin.le.2)then                                                   9d12s19
       ibaseo=ibc(ivecpt+nopenin)                                            9d12s19
      else                                                              9d12s19
       ibaseo=ibc(ivecpp+nopenin)                                            9d12s19
      end if                                                            9d12s19
      if(ldebug)write(6,*)('ibaseo: '),ibaseo,nopenin
      ipms=ibaseo+ms2+nopenin+1                                              9d12s19
      indp=ibc(ipms)                                                    9d12s19
      indp=indp+ibaseo                                                  9d12s19
      ndet=ibc(indp)                                                    9d12s19
      if(ldebug)write(6,*)('no. dets: '),ndet
      itmp=ibcoff                                                       9d12s19
      ibcoff=itmp+ndet*nopenin                                          9d12s19
      if(ibcoff.lt.0)write(6,*)('enoughc '),ibcoff
      call enough('stepwisecsf.  8',bc,ibc)
      call prtsigs(ibc(indp+1),ndet,nopenin,ibc(itmp))                  9d12s19
      do i=0,ndet*nopenin-1                                             9d12s19
       iaorb1(i+1)=ibc(itmp+i)                                          9d12s19
      end do                                                            9d12s19
      ibcoff=itmp                                                       12d1s20
      icall=ndet                                                        9d12s19
      jptr=indp+ndet
      kptr=jptr+is2+1
      nroota=ibc(kptr)
      nrooto=ibc(nroota+ibaseo)
      if(ldebug)write(6,*)('no. roots: '),nrooto
      ivec=nroota+ibaseo+1                                              9d12s19
      if(ldebug)then                                                    12d12s19
       write(6,*)('vectors: '),ivec
       call prntm2(bc(ivec),ndet,nrooto,ndet)                            9d12s19
      end if                                                            12d12s19
      vec(1)=dfloat(ivec)                                               12d1s20
      return                                                            10d2s20
      end
