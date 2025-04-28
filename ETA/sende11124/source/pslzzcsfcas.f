c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine pslzzcsfcas(npsf,hdps,ipsbase,nec,hamps,nps,ncsf,mdon, 12d28s19
     $     ixlzz,multh,islz,nkeep,ism,irel,                             8d1s22
     $     irefo,idoubo,lambdacas,npass,lwrite,mdoo,iptrbit,ixw1,ixw2,  4d12s21
     $     mysym,norb,bc,ibc)                                           11d10s22
      implicit real*8 (a-h,o-z)                                         6d11s19
      external second                                                   7d30s19
      integer*8 ipack,itesta,itestb                                     4d12s21
      integer*2 ipack2(4)                                               11d4s19
      logical lwrite                                                    12d31s19
      integer*1 icode(64),imap(64),nab(2),isame(64),                    8d1s22
     $     nabi(2),nabj(2),nab1(2),nab2(2),itest(64,2)                  4d12s21
      equivalence (ipack,ipack2)                                        11d4s19
      dimension hdps(*),hamps(*),ncsf(*),nother(2),                     8d1s22
     $     ipsbase(3,*),multh(8,8),test(10),ixlzz(8,*),ism(*),irel(*),  12d31s19
     $     irefo(*),idoubo(*),islz(*),nab4(2,3),iptrbit(2,mdoo+1,*)     4d12s21
      include "common.store"                                            6d11s19
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      ione=npass+1                                                      12d31s19
      ihtmp=ibcoff                                                      7d11s19
      jhtmp=ihtmp                                                       7d11s19
      loopit=0                                                          7d10s19
      do i=1,10
       test(i)=0d0
      end do                                                            7d5s19
      ips=0
      do if=1,npsf                                                      6d11s19
       nclo=ipsbase(1,if)                                               6d11s19
       nclop=nclo+1                                                     6d11s19
       iarg=nclo+1-mdon                                                 6d11s19
       ip=ips+1
       nopen=nec-2*nclo
       iic=iptrbit(1,nclop,mysym)+ipsbase(2,if)-1                       4d12s21
       iio=iptrbit(2,nclop,mysym)+ipsbase(3,if)-1                       4d12s21
       jps=ips                                                          12d11s19
       do jf=if,npsf                                                    12d11s19
        ncloj=ipsbase(1,jf)                                              7d11s19
        if(iabs(ncloj-nclo).le.2)then                                   12d9s19
         nclopj=ncloj+1                                                  6d11s19
         nopenj=nec-2*ncloj
         jarg=ncloj+1-mdon                                               6d11s19
         jjo=iptrbit(2,nclopj,mysym)+ipsbase(3,jf)-1                    4d12s21
         jjc=iptrbit(1,nclopj,mysym)+ipsbase(2,jf)-1                    4d12s21
         jp=jps+1
         iflag=0                                                        4d12s21
         call gandc4(ibc(jjc),ibc(jjo),ibc(iic),ibc(iio),nopenj,nopen,  4d12s21
     $         norb,nnot4,nab4)                                         4d12s21
         iok=0                                                          10d2s20
         if(nnot4.gt.1)then                                              10d2s20
          iok=1                                                         10d2s20
         else if(nnot4.eq.1.and.nopen.gt.0)then                          10d2s20
          iok=1                                                         10d2s20
         end if                                                         10d2s20
         if(iok.ne.0)then                                               4d12s21
          itmp=ibcoff                                                   12d9s19
          jhtmp=itmp
          nn=ncsf(jarg)*ncsf(iarg)                                      12d9s19
          nnm=nn-1                                                      4d12s21
          ibcoff=itmp+nn                                                12d9s19
          do i=0,nnm                                                    4d12s21
           bc(jhtmp+i)=0d0                                              12d11s19
          end do                                                        12d11s19
          jhtmp=ibcoff                                                  12d11s19
          if(mod(loopit,mynprocg).eq.mynowprog)then                      12d9s19
           loopit=mynowprog                                             9d10s24
           if(nab4(1,1).eq.0.and.nnot4.ne.1)then                        4d12s21
            write(6,*)('ooops, nab4 is '),nab4                          4d12s21
            stop                                                        4d12s21
           end if                                                       4d12s21
           ibcsav=ibcoff                                                  12d9s19
           call enough('pslzzcsfcas.  1',bc,ibc)
           ihtmpgg=itmp                                                 4d12s21
           if(nnot4.eq.1)then                                             11d13s20
            do i=0,ncsf(jarg)-1                                           11d13s20
             do j=0,i-1                                                   11d13s20
              ji=ihtmpgg+j+ncsf(jarg)*i                                   11d13s20
              bc(ji)=0d0                                                  11d13s20
             end do                                                       11d13s20
             ii=ihtmpgg+i*(ncsf(jarg)+1)                                  11d13s20
             bc(ii)=hdps(if)                                              11d13s20
             do j=i+1,ncsf(jarg)-1                                        11d13s20
              ji=ihtmpgg+j+ncsf(jarg)*i                                   11d13s20
              bc(ji)=0d0                                                  11d13s20
             end do                                                       11d13s20
            end do                                                        11d13s20
            do i1=1,norb-1                                              8d1s22
             if(btest(ibc(iio),i1))then                                 8d1s22
              nab1(2)=i1                                                8d1s22
              do i2=i1+1,norb                                           8d1s22
               if(btest(ibc(iio),i2))then                               8d1s22
                nab1(1)=i2                                              8d1s22
                nopenk=nopen-2                                              11d13s20
                karg=iarg+1                                                 11d13s20
                nab2(1)=nab1(2)                                             11d13s20
                nab2(2)=nab1(1)                                             11d13s20
                jsa=ism(nab1(1))                                              11d13s20
                jga=irel(nab1(1))-1+idoubo(jsa)                           1d21s21
                jsb=ism(nab1(2))                                              11d13s20
                jgb=irel(nab1(2))-1+idoubo(jsb)                           1d21s21
                jsc=ism(nab2(1))                                              11d13s20
                jgc=irel(nab2(1))-1+idoubo(jsc)                           1d21s21
                jsd=ism(nab2(2))                                              11d13s20
                jgd=irel(nab2(2))-1+idoubo(jsd)                           1d21s21
                xint=0d0                                                  1d21s21
                do ipass=1,npass                                            1d2s20
                 if(multh(jsa,jsb).eq.islz(ipass))then                      1d2s20
                  ndma=irefo(jsa)+idoubo(jsa)
                  iada=ixlzz(jsa,ipass)+jga+ndma*jgb                      1d26s21
                  ndmc=irefo(jsc)+idoubo(jsc)
                  iadc=ixlzz(jsc,ipass)+jgc+ndmc*jgd                      1d26s21
                  xint=xint-2d0*bc(iada)*bc(iadc)                         1d22s21
                 end if                                                      12d24s19
                end do                                                      1d2s20
                nqq=karg+mdon-1
                if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                 itesta=ibc(iic)                                             11d13s20
                 itestb=ibc(iio)                                             11d13s20
                 itestb=ibclr(itestb,nab1(1))                             8d1s22
                 itestb=ibclr(itestb,nab1(2))                             8d1s22
                 itesta=ibset(itesta,nab1(2))                             8d1s22
                 call gandc(ibc(jjc),ibc(jjo),itesta,itestb,nopenj,       1d21s21
     $                nopenk,jarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,  1d21s21
     $                iwpb1,iwpk1,ncsfmid1,bc,ibc)                      11d14s22
                 iprod=ibcoff                                               11d13s20
                 itmp1=iprod+ncsf(jarg)*ncsf(jarg)                          11d13s20
                 itmp2=itmp1+ncsf(jarg)*ncsf(karg)                          11d13s20
                 ibcoff=itmp2+ncsf(jarg)*ncsf(karg)                         11d13s20
                 call enough('pslzzcsfcas.  2',bc,ibc)
                 call prodn(iwpb1,iwpk1,ncsf(jarg),ncsf(karg),ncsfmid1,   12d4s20
     $              bc(itmp1),bc,ibc,1d0,0d0)                           2d13s23
                 do k=0,ncsf(karg)-1                                        11d13s20
                  do j=0,ncsf(jarg)-1                                       11d13s20
                   jk=itmp1+j+ncsf(jarg)*k                                  11d13s20
                   kj=itmp2+k+ncsf(karg)*j                                  11d13s20
                   bc(kj)=bc(jk)                                            11d13s20
                  end do                                                    11d13s20
                 end do                                                     11d13s20
                 call dgemm('n','n',ncsf(jarg),ncsf(jarg),ncsf(karg),     1d21s21
     $                1d0,bc(itmp1),ncsf(jarg),bc(itmp2),ncsf(karg),0d0,1d21s21
     $               bc(iprod),ncsf(jarg),                                 11d13s20
     d' pslzzcsfcas.  1')
                 do i=0,ncsf(jarg)-1                                         11d13s20
                  ii=iprod+i*(ncsf(jarg)+1)                                  11d13s20
                  bc(ii)=bc(ii)-1d0                                          11d13s20
                 end do                                                      11d13s20
                 do i=0,ncsf(jarg)*ncsf(jarg)-1                              11d13s20
                  bc(ihtmpgg+i)=bc(ihtmpgg+i)+xint*bc(iprod+i)               11d13s20
                 end do                                                      11d13s20
                 ibcoff=iprod                                                11d13s20
                else                                                        11d13s20
                 do i=0,ncsf(jarg)-1                                        11d13s20
                  ii=ihtmpgg+i*(ncsf(jarg)+1)                               11d13s20
                  bc(ii)=bc(ii)-xint                                        11d13s20
                 end do                                                     11d13s20
                end if                                                      11d13s20
               end if                                                   8d1s22
              end do                                                       11d13s20
             end if                                                     8d1s22
            end do                                                        11d13s20
           else if(nnot4.eq.2)then                                      1d21s21
            do i=0,ncsf(iarg)*ncsf(jarg)-1                              1d21s21
             bc(ihtmpgg+i)=0d0                                          1d21s21
            end do                                                      1d21s21
           else if(nnot4.ge.3)then                                        11d13s20
            if(nnot4.eq.3)then                                          1d26s21
             ipssx=1                                                    1d26s21
            else                                                        1d26s21
             ipssx=3                                                     1d26s21
            end if                                                      1d26s21
            do ipss=1,ipssx                                             12d8s20
             if(ipss.eq.1)then
              iu1=1
              iu2=1
             else if(ipss.eq.2)then
              iu1=1
              iu2=2
             else
              iu1=2
              iu2=1
             end if                                                     12d8s20
             itesta=ibc(jjc)                                               11d13s20
             itestb=ibc(jjo)                                               11d13s20
c     j is bra
             if(btest(itesta,nab4(1,iu1)))then                          1d26s21
              itesta=ibclr(itesta,nab4(1,iu1))                          1d26s21
              itestb=ibset(itestb,nab4(1,iu1))                          1d26s21
              nopenk=nopenj+1                                              11d13s20
              karg=jarg-1                                                  11d13s20
             else if(btest(itestb,nab4(1,iu1)))then                     1d26s21
              itestb=ibclr(itestb,nab4(1,iu1))                          1d26s21
              nopenk=nopenj-1                                              11d13s20
              karg=jarg                                                  11d13s20
             else                                                          11d13s20
              write(6,*)('bit not set for nab4(1,1) = '),nab4(1,1)       11d27s20
              stop 'nab4(1,1)'                                           11d27s20
             end if                                                        11d13s20
             if(btest(itestb,nab4(2,iu2)))then                          1d26s21
              itesta=ibset(itesta,nab4(2,iu2))                          1d26s21
              itestb=ibclr(itestb,nab4(2,iu2))                          1d26s21
              nopenk=nopenk-1                                              11d13s20
              karg=karg+1                                                  11d13s20
             else if(btest(itesta,nab4(2,iu2)))then                     1d26s21
              write(6,*)('already double in nab4(2,1) = '),nab4(2,iu2)     11d27s20
              stop 'nab4(2,1)'                                           11d27s20
             else                                                          11d13s20
              itestb=ibset(itestb,nab4(2,iu2))                          1d26s21
              nopenk=nopenk+1                                              11d13s20
             end if                                                        11d13s20
             nqq=karg+mdon-1                                            1d26s21
             if(nqq.ge.mdon.and.nqq.le.mdoo)then                        1d26s21
              call gandc(ibc(jjc),ibc(jjo),itesta,itestb,nopenj,nopenk,     11d13s20
     $         jarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,    11d13s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
              call gandc(itesta,itestb,ibc(iic),ibc(iio),nopenk,nopen,      11d13s20
     $         karg,iarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,    11d13s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
              if(nnot1.eq.2.and.nnot2.eq.2)then                         1d26s21
               call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),ncsfmid1,
     $             iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod,bc,ibc)       11d10s22
               jsa=ism(nab1(1))                                              11d13s20
               jga=irel(nab1(1))-1+idoubo(jsa)                            1d21s21
               jsb=ism(nab1(2))                                              11d13s20
               jgb=irel(nab1(2))-1+idoubo(jsb)                                           11d13s20
               jsc=ism(nab2(1))                                              11d13s20
               jgc=irel(nab2(1))-1+idoubo(jsc)                                           11d13s20
               jsd=ism(nab2(2))                                              11d13s20
               jgd=irel(nab2(2))-1+idoubo(jsd)                                           11d13s20
               xint=0d0                                                    1d21s21
               do ipass=1,npass                                            1d2s20
                if(multh(jsa,jsb).eq.islz(ipass))then                      1d2s20
                 ndma=irefo(jsa)+idoubo(jsa)
                 iada=ixlzz(jsa,ipass)+jga+ndma*jgb                     1d26s21
                 ndmc=irefo(jsc)+idoubo(jsc)
                 iadc=ixlzz(jsc,ipass)+jgc+ndmc*jgd                     1d26s21
                 xint=xint-2d0*bc(iada)*bc(iadc)                          1d21s21
                end if                                                      12d24s19
               end do                                                      1d2s20
               if(nab1(1).eq.nab2(1).and.nab1(2).eq.nab2(2))              1d21s21
     $              xint=xint*0.5d0                                      1d21s21
               do i=0,ncsf(jarg)*ncsf(iarg)-1                                11d13s20
                bc(ihtmpgg+i)=bc(ihtmpgg+i)+bc(iprod+i)*xint            1d26s21
               end do                                                        11d13s20
               ibcoff=iprod                                                  11d13s20
               if(ipss.eq.2)go to 22                                    1d26s21
              end if                                                    1d26s21
             end if                                                     1d26s21
            end do                                                      1d26s21
   22       continue                                                    1d26s21
           end if                                                         11d13s20
          end if                                                        1d21s21
          loopit=loopit+1                                               4d12s21
         end if                                                         1d21s21
        end if                                                          12d9s19
        jps=jps+ncsf(jarg)
       end do                                                           6d11s19
       ips=ips+ncsf(iarg)
      end do
      nhtmp=jhtmp-ihtmp                                                 7d11s19
      call dws_gsumf(bc(ihtmp),nhtmp)                                   7d11s19
      jhtmp=ihtmp                                                       7d11s19
      npsx=nps                                                          12d24s19
      ntri=(npsx*(npsx+1))/2                                            9d4s19
      do i=1,ntri                                                       7d11s19
       hamps(i)=0d0                                                     7d11s19
      end do                                                            7d11s19
      ih=0                                                              7d11s19
      ips=0
      irowp=10
      icolp=4
      do if=1,npsf                                                      6d11s19
       nclo=ipsbase(1,if)                                               6d11s19
       nclop=nclo+1                                                     6d11s19
       iarg=nclo+1-mdon                                                 6d11s19
       ip=ips                                                           7d11s19
       nopen=nec-2*nclo
       if(nopen.gt.0)then                                               6d11s19
        do ik=1,ncsf(iarg)                                              7d11s19
         ikp=ik+ips                                                     7d11s19
         do ib=1,ncsf(iarg)                                             7d11s19
          ibp=ib+ips                                                    7d11s19
          if(ik.ge.ib)then                                              7d11s19
           itry=((ikp*(ikp-1))/2)+ibp                                   7d11s19
           hamps(itry)=bc(jhtmp)                                         7d11s19
          end if                                                        7d11s19
          jhtmp=jhtmp+1                                                 7d11s19
         end do                                                         7d11s19
        end do                                                          7d11s19
       else                                                             7d11s19
        do i=1,ncsf(iarg)                                                6d11s19
         ipp=ips+i                                                        6d11s19
         itry=((ipp*(ipp-1))/2)+ipp                                     7d11s19
         hamps(itry)=hdps(if)                                            7d11s19
        end do                                                           6d11s19
       end if                                                           7d11s19
       jps=ips+ncsf(iarg)                                               6d19s19
       iic=iptrbit(1,nclop,mysym)+ipsbase(2,if)-1                       4d12s21
       iio=iptrbit(2,nclop,mysym)+ipsbase(3,if)-1                       4d12s21
       do jf=if+1,npsf                                                  6d19s19
        ncloj=ipsbase(1,jf)                                             6d11s19
        nclopj=ncloj+1                                                  12d11s19
        jarg=ncloj+1-mdon                                               6d11s19
        nopenj=nec-2*ncloj
        jjo=iptrbit(2,nclopj,mysym)+ipsbase(3,jf)-1                     8d1s22
        jjc=iptrbit(1,nclopj,mysym)+ipsbase(2,jf)-1                     8d1s22
        call gandc4(ibc(jjc),ibc(jjo),ibc(iic),ibc(iio),nopenj,nopen,   8d1s22
     $         norb,nnot,nab4,bc,ibc)                                   11d14s22
        if(nnot.gt.0)then                                               12d30s19
         do ik=1,ncsf(iarg)                                             12d11s19
          ikp=ik+ips                                                    12d11s19
          do ib=1,ncsf(jarg)                                            12d11s19
           jbp=ib+jps                                                   12d11s19
           itry=((jbp*(jbp-1))/2)+ikp                                   12d11s19
           hamps(itry)=bc(jhtmp)                                        12d11s19
           jhtmp=jhtmp+1                                                12d11s19
          end do                                                        12d11s19
         end do                                                         12d11s19
        end if                                                          12d11s19
        jps=jps+ncsf(jarg)
       end do                                                           7d11s19
       ips=ips+ncsf(iarg)
      end do
      call square(hamps,ips)                                            12d24s19
      ibcoff=ihtmp                                                      12d24s19
      ieig=ibcoff
      ivec=ieig+ips
      isyr=ivec+ips*ips
      ibcoff=isyr+ips
      call enough('pslzzcsfcas.  3',bc,ibc)
      call diagx(ips,hamps,bc(ieig),bc(ivec),ibc(isyr),bc,ibc)          11d14s22
      if(lwrite)then                                                    1d2s20
       if(npass.eq.1)then                                               1d2s20
        write(6,*)('looking for lambda '),lambdacas                     1d2s20
       else                                                             1d2s20
        write(6,*)('looking for L '),lambdacas                          1d2s20
       end if                                                           1d2s20
      end if                                                            1d2s20
      nkeep=0                                                           12d24s19
      nbad=0                                                            12d28s19
      do i=0,nps-1
       if(npass.eq.1)then                                               12d31s19
        try=sqrt(abs(bc(ieig+i)))                                        12d24s19
        itry=nint(try)                                                   12d24s19
        delta=abs(bc(ieig+i)-dfloat(itry*itry))                          12d24s19
       else                                                             12d31s19
        try=sqrt(abs(bc(ieig+i))+0.25d0)-0.5d0                          12d31s19
        itry=nint(try)                                                   12d24s19
        delta=abs(bc(ieig+i)-dfloat(itry*(itry+1)))                     12d31s19
       end if                                                           12d31s19
       if(delta.gt.1d-10)then                                           8d9s22
        nbad=nbad+1                                                     8d9s22
        write(6,*)('bad eigenvalue: '),bc(ieig+i)
        write(6,*)('vector '),ips,nps
        jvec=ivec+ips*(i-1)
        call prntm2(bc(jvec),ips,1,ips)
       end if                                                           8d9s22
       if(delta.gt.1d-10.or.itry.eq.lambdacas)then                       12d28s19
        jvec=ivec+nps*i                                                 12d24s19
        ikeep=ivec+nps*nkeep                                            12d24s19
        do j=0,nps-1                                                    12d24s19
         bc(ikeep+j)=bc(jvec+j)                                         12d24s19
        end do                                                          12d24s19
        nkeep=nkeep+1                                                   12d24s19
       end if                                                           12d24s19
      end do
      if(lwrite)then                                                    12d31s19
       write(6,*)('number of keepers = '),nkeep                          12d24s19
       if(nbad.gt.0)write(6,*)('of which '),nbad,                        12d28s19
     $     (' were not pure eigenvalues')                               2d10s20
      end if                                                            12d31s19
      if(nbad.gt.0)then
       write(6,*)('perhaps your orbitals are not properly paired up?')  8d25s22
       do ipass=1,npass
        write(6,*)('for pass '),ipass
        do isb=1,4
         jsb=multh(isb,islz(ipass))
         if(min(irefo(jsb),irefo(isb)).gt.0)then
          write(6,*)('what we have symmetry block '),isb,jsb
          call prntm2(bc(ixlzz(isb,ipass)),irefo(isb),irefo(jsb),
     $        irefo(isb))
         end if
        end do
       end do
       call dws_synca
       call dws_finalize
       stop
      end if
      jvec=ivec-1                                                       12d24s19
      do j=1,nkeep*nps                                                  12d24s19
       hamps(j)=bc(jvec+j)                                              12d24s19
      end do                                                            12d24s19
      call dws_bcast(hamps,nkeep*nps)                                   12d24s19
      ibcoff=ieig                                                       12d24s19
      return
      end
