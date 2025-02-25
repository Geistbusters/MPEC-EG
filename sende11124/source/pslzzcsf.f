c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine pslzzcsf(npsf,hdps,ipsbase,nec,hamps,nps,ncsf,mdon,    6d11s19
     $     icsfpd,iptr,ixlzz,idorb,isorb,multh,islz,nkeep,npass,lwrite, 1d26s21
     $     mdoo,iptrbit,ixw1,ixw2,mysym,lz2,bc,ibc)                     11d10s22
      implicit real*8 (a-h,o-z)                                         6d11s19
      external second                                                   7d30s19
      integer*8 ipack,itesta,itestb                                     1d21s21
      logical lwrite,ldebug                                             5d21s21
      integer*2 ipack2(4)                                               11d4s19
      character*9 fname                                                 7d17s19
      integer*1 idorb(*),isorb(*),icode(64),imap(64),nab(2),isame(64),  1d21s21
     $     nabi(2),nabj(2),nab1(2),nab2(2),itest(64,2)                  1d21s21
      equivalence (ipack,ipack2)                                        11d4s19
      dimension hdps(npsf,*),hamps(*),ncsf(*),icsfpd(*),iptr(4,*),      5d14s21
     $     nother(2),                                                   5d14s21
     $    ipsbase(3,*),multh(8,8),test(10),ixlzz(8,*),islz(*),nab4(2,3),1d21s21
     $     iptrbit(2,mdoo+1,*)                                          1d21s21
      include "common.store"                                            6d11s19
      include "common.mrci"                                             7d17s19
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      ione=npass+1                                                      12d31s19
      ldebug=.false.
      if(ldebug)write(6,*)('in pslzzcsf '),mdoo,ixw1,ixw2,mysym
      ibcoffo=ibcoff                                                    5d14s21
      ihtmp=ibcoff                                                      7d11s19
      jhtmp=ihtmp                                                       7d11s19
      ibcoff=ihtmp+nps*nps                                              12d30s21
      ihtmpzz=ibcoff                                                    1d19s23
      if(npass.ne.1)then                                                5d14s21
       ibcoff=ihtmpzz+nps*nps                                           5d14s21
       jhtmpzz=ihtmpzz                                                  5d14s21
      end if                                                            5d14s21
      call enough('pslzzcsf.  1',bc,ibc)
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
       ic=iptr(2,nclop)+nclo*(ipsbase(2,if)-1)                           7d11s19
       io=iptr(4,nclop)+nopen*(ipsbase(3,if)-1)                          12d9s19
       iic=iptrbit(1,nclop,mysym)+ipsbase(2,if)-1                       1d21s21
       iio=iptrbit(2,nclop,mysym)+ipsbase(3,if)-1                       1d21s21
       jps=ips                                                          12d11s19
       do jf=if,npsf                                                    12d11s19
        ncloj=ipsbase(1,jf)                                              7d11s19
        if(iabs(ncloj-nclo).le.2)then                                   12d9s19
         nclopj=ncloj+1                                                  6d11s19
         nopenj=nec-2*ncloj
         jarg=ncloj+1-mdon                                               6d11s19
         jo=iptr(4,nclopj)+nopenj*(ipsbase(3,jf)-1)                       7d11s19
         jc=iptr(2,nclopj)+ncloj*(ipsbase(2,jf)-1)                        7d11s19
         jjo=iptrbit(2,nclopj,mysym)+ipsbase(3,jf)-1                    1d21s21
         jjc=iptrbit(1,nclopj,mysym)+ipsbase(2,jf)-1                    1d21s21
         jp=jps+1
         call gcode(idorb(jc),ncloj,isorb(jo),nopenj,
     $       idorb(ic),nclo,isorb(io),nopen,icode,imap,nnot,nx,0,nab)
          iflag=0
          call gandc4(ibc(jjc),ibc(jjo),ibc(iic),ibc(iio),nopenj,nopen,  11d13s20
     $         norb,nnot4,nab4,bc,ibc)                                  11d14s22
          iok=0                                                          10d2s20
          if(nnot4.gt.1)then                                              10d2s20
           iok=1                                                         10d2s20
          else if(nnot4.eq.1.and.nopen.gt.0)then                          10d2s20
           iok=1                                                         10d2s20
          end if                                                         10d2s20
          if(iok.ne.0)then                                              1d21s21
           if(ldebug)then
            write(6,*)('nnot from gandc4 = '),nnot
           call dcbit(ibc(jjc),norb,'jjc')
           call dcbit(ibc(jjo),norb,'jjo')
           call dcbit(ibc(iic),norb,'iic')
           call dcbit(ibc(iio),norb,'iio')
           end if
           itmp=jhtmp                                                   5d14s21
           nn=ncsf(jarg)*ncsf(iarg)                                      12d9s19
           nnm=nn-1                                                      9d20s20
           do i=0,nnm                                                    9d20s20
            bc(jhtmp+i)=0d0                                              12d11s19
           end do                                                        12d11s19
           jhtmp=jhtmp+nn                                               5d14s21
           if(npass.ne.1)then                                           5d14s21
            itmpzz=jhtmpzz                                              5d14s21
            jhtmpzz=jhtmpzz+nn                                          5d14s21
            do i=0,nnm                                                  5d14s21
             bc(itmpzz+i)=0d0                                           5d14s21
            end do                                                      5d14s21
           end if                                                       5d14s21
           if(mod(loopit,mynprocg).eq.mynowprog)then                      12d9s19
            loopit=mynowprog                                            9d10s24
            if(nab4(1,1).eq.0.and.nnot4.ne.1)then                       1d21s21
             write(6,*)('ooops, nab4 is '),nab4
             stop
            end if
            ibcsav=ibcoff                                                9d21s20
            ihtmpgg=itmp                                                11d13s20
            if(nnot4.eq.1)then                                             11d13s20
             if(ldebug)write(6,*)('diagonal contribution: '),hdps(if,1),
     $           hdps(if,2)
             do i=0,ncsf(jarg)-1                                           11d13s20
              do j=0,i-1                                                   11d13s20
               ji=ihtmpgg+j+ncsf(jarg)*i                                   11d13s20
               bc(ji)=0d0                                                  11d13s20
              end do                                                       11d13s20
              ii=ihtmpgg+i*(ncsf(jarg)+1)                                  11d13s20
              bc(ii)=hdps(if,1)                                         5d14s21
              do j=i+1,ncsf(jarg)-1                                        11d13s20
               ji=ihtmpgg+j+ncsf(jarg)*i                                   11d13s20
               bc(ji)=0d0                                                  11d13s20
              end do                                                       11d13s20
             end do                                                        11d13s20
             if(npass.ne.1)then                                         5d14s21
              do i=0,ncsf(jarg)-1                                       5d14s21
               ii=itmpzz+i*(ncsf(jarg)+1)                               5d14s21
               bc(ii)=hdps(if,2)                                        5d14s21
              end do                                                    5d14s21
             end if                                                     5d14s21
             do i1=0,nopen-2                                               11d13s20
              do i2=i1+1,nopen-1                                           11d13s20
               nab1(1)=isorb(io+i2)                                      12d6s20
               nab1(2)=isorb(io+i1)                                      12d6s20
               nopenk=nopen-2                                              11d13s20
               karg=iarg+1                                                 11d13s20
               nab2(1)=nab1(2)                                             11d13s20
               nab2(2)=nab1(1)                                             11d13s20
               jsa=ism(nab1(1))                                              11d13s20
               jga=irel(nab1(1))-1+idoubo(jsa)                          1d21s21
               jsb=ism(nab1(2))                                              11d13s20
               jgb=irel(nab1(2))-1+idoubo(jsb)                          1d21s21
               jsc=ism(nab2(1))                                              11d13s20
               jgc=irel(nab2(1))-1+idoubo(jsc)                          1d21s21
               jsd=ism(nab2(2))                                              11d13s20
               jgd=irel(nab2(2))-1+idoubo(jsd)                          1d21s21
               xint=0d0                                                 1d21s21
               xintzz=0d0                                               5d14s21
               do ipass=1,npass                                            1d2s20
                if(multh(jsa,jsb).eq.islz(ipass))then                      1d2s20
                 ndma=irefo(jsa)+idoubo(jsa)
                 iada=ixlzz(jsa,ipass)+jga+ndma*jgb                     1d26s21
                 ndmc=irefo(jsc)+idoubo(jsc)
                 iadc=ixlzz(jsc,ipass)+jgc+ndmc*jgd                     1d26s21
                 xint=xint-2d0*bc(iada)*bc(iadc)                        1d22s21
                 if(ipass.eq.3)xintzz=xintzz-2d0*bc(iada)*bc(iadc)      5d14s21
                end if                                                      12d24s19
               end do                                                      1d2s20
               nqq=karg+mdon-1
               if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                itesta=ibc(iic)                                             11d13s20
                itestb=ibc(iio)                                             11d13s20
                itestb=ibclr(itestb,isorb(io+i2))                                     11d13s20
                itestb=ibclr(itestb,isorb(io+i1))                                     11d13s20
                itesta=ibset(itesta,isorb(io+i1))                           11d13s20
                call gandc(ibc(jjc),ibc(jjo),itesta,itestb,nopenj,      1d21s21
     $                nopenk,jarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,  1d21s21
     $                iwpb1,iwpk1,ncsfmid1,bc,ibc)                      11d14s22
                iprod=ibcoff                                               11d13s20
                itmp1=iprod+ncsf(jarg)*ncsf(jarg)                          11d13s20
                itmp2=itmp1+ncsf(jarg)*ncsf(karg)                          11d13s20
                ibcoff=itmp2+ncsf(jarg)*ncsf(karg)                         11d13s20
                call enough('pslzzcsf.  2',bc,ibc)
                call prodn(iwpb1,iwpk1,ncsf(jarg),ncsf(karg),ncsfmid1,   12d4s20
     $              bc(itmp1),bc,ibc,1d0,0d0)                           2d13s23
                do k=0,ncsf(karg)-1                                        11d13s20
                 do j=0,ncsf(jarg)-1                                       11d13s20
                  jk=itmp1+j+ncsf(jarg)*k                                  11d13s20
                  kj=itmp2+k+ncsf(karg)*j                                  11d13s20
                  bc(kj)=bc(jk)                                            11d13s20
                 end do                                                    11d13s20
                end do                                                     11d13s20
                call dgemm('n','n',ncsf(jarg),ncsf(jarg),ncsf(karg),    1d21s21
     $                1d0,bc(itmp1),ncsf(jarg),bc(itmp2),ncsf(karg),0d0,1d21s21
     $               bc(iprod),ncsf(jarg),                                 11d13s20
     d' pslzzcsf.  1')
                do i=0,ncsf(jarg)-1                                         11d13s20
                 ii=iprod+i*(ncsf(jarg)+1)                                  11d13s20
                 bc(ii)=bc(ii)-1d0                                          11d13s20
                end do                                                      11d13s20
                 if(ldebug.and.max(abs(xint),abs(xintzz)).gt.1d-10)then
                  write(6,*)('xint '),i1,i2,xint,xintzz,nopenj,nopenk,
     $                 jarg,karg
                  call dcbit(itesta,norb,'itesta')
                  call dcbit(itestb,norb,'itestb')
                  call prntm2(bc(iprod),ncsf(jarg),ncsf(iarg),
     $                 ncsf(jarg))
                 end if
                do i=0,ncsf(jarg)*ncsf(jarg)-1                              11d13s20
                 bc(ihtmpgg+i)=bc(ihtmpgg+i)+xint*bc(iprod+i)               11d13s20
                end do                                                      11d13s20
                if(npass.ne.1)then                                      5d14s21
                 do i=0,ncsf(jarg)*ncsf(jarg)-1                         5d14s21
                  bc(itmpzz+i)=bc(itmpzz+i)+xintzz*bc(iprod+i)          5d14s21
                 end do                                                 5d14s21
                end if                                                  5d14s21
                ibcoff=iprod                                                11d13s20
               else                                                        11d13s20
                if(ldebug.and.max(abs(xint),abs(xintzz)).gt.1d-10)
     $               write(6,*)('karg is out of range '),xint,
     $               xintzz
                do i=0,ncsf(jarg)-1                                        11d13s20
                 ii=ihtmpgg+i*(ncsf(jarg)+1)                               11d13s20
                 bc(ii)=bc(ii)-xint                                        11d13s20
                end do                                                     11d13s20
                if(npass.ne.1)then                                      5d14s21
                 do i=0,ncsf(jarg)-1                                    5d14s21
                  ii=itmpzz+i*(ncsf(jarg)+1)                            5d14s21
                  bc(ii)=bc(ii)-xintzz                                  5d14s21
                 end do                                                 5d14s21
                end if                                                  5d14s21
               end if                                                      11d13s20
              end do                                                       11d13s20
             end do                                                        11d13s20
            else if(nnot4.eq.2)then                                     1d21s21
             do i=0,ncsf(iarg)*ncsf(jarg)-1                             1d21s21
              bc(ihtmpgg+i)=0d0                                         1d21s21
             end do                                                     1d21s21
            else if(nnot4.ge.3)then                                        11d13s20
             if(ldebug)write(6,*)('in nnot4 block '),nnot4
             if(nnot4.eq.3)then                                         1d26s21
              ipssx=1                                                   1d26s21
             else                                                       1d26s21
              ipssx=3                                                   1d26s21
             end if                                                     1d26s21
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
              if(btest(itesta,nab4(1,iu1)))then                         1d26s21
               itesta=ibclr(itesta,nab4(1,iu1))                         1d26s21
               itestb=ibset(itestb,nab4(1,iu1))                         1d26s21
               nopenk=nopenj+1                                              11d13s20
               karg=jarg-1                                                  11d13s20
              else if(btest(itestb,nab4(1,iu1)))then                    1d26s21
               itestb=ibclr(itestb,nab4(1,iu1))                         1d26s21
               nopenk=nopenj-1                                              11d13s20
               karg=jarg                                                  11d13s20
              else                                                          11d13s20
               write(6,*)('bit not set for nab4(1,1) = '),nab4(1,1)       11d27s20
               stop 'nab4(1,1)'                                           11d27s20
              end if                                                        11d13s20
              if(btest(itestb,nab4(2,iu2)))then                         1d26s21
               itesta=ibset(itesta,nab4(2,iu2))                         1d26s21
               itestb=ibclr(itestb,nab4(2,iu2))                         1d26s21
               nopenk=nopenk-1                                              11d13s20
               karg=karg+1                                                  11d13s20
              else if(btest(itesta,nab4(2,iu2)))then                    1d26s21
               write(6,*)('already double in nab4(2,1) = '),nab4(2,iu2)     11d27s20
               stop 'nab4(2,1)'                                           11d27s20
              else                                                          11d13s20
               itestb=ibset(itestb,nab4(2,iu2))                         1d26s21
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
               if(nnot1.eq.2.and.nnot2.eq.2)then                        1d26s21
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
                xintzz=0d0                                              5d14s21
                do ipass=1,npass                                            1d2s20
                 if(multh(jsa,jsb).eq.islz(ipass))then                      1d2s20
                  ndma=irefo(jsa)+idoubo(jsa)
                  iada=ixlzz(jsa,ipass)+jga+ndma*jgb                    1d26s21
                  ndmc=irefo(jsc)+idoubo(jsc)
                  iadc=ixlzz(jsc,ipass)+jgc+ndmc*jgd                    1d26s21
                  xint=xint-2d0*bc(iada)*bc(iadc)                          1d21s21
                  if(ipass.eq.3)xintzz=xintzz-2d0*bc(iada)*bc(iadc)     5d14s21
                 end if                                                      12d24s19
                end do                                                      1d2s20
                if(ldebug.and.max(abs(xint),abs(xintzz)).gt.1d-10)
     $               write(6,*)('initial xint '),xint,xintzz
                if(nab1(1).eq.nab2(1).and.nab1(2).eq.nab2(2))              1d21s21
     $             xint=xint*0.5d0                                      1d21s21
                do i=0,ncsf(jarg)*ncsf(iarg)-1                                11d13s20
                 bc(ihtmpgg+i)=bc(ihtmpgg+i)+bc(iprod+i)*xint           1d26s21
                end do                                                        11d13s20
                if(npass.ne.1)then                                      5d14s21
                 if(nab1(1).eq.nab2(1).and.nab1(2).eq.nab2(2))          5d14s21
     $             xintzz=xintzz*0.5d0                                  5d14s21
                 do i=0,ncsf(jarg)*ncsf(iarg)-1                         5d14s21
                  bc(itmpzz+i)=bc(itmpzz+i)+bc(iprod+i)*xintzz          5d14s21
                 end do                                                 5d14s21
                end if                                                  5d14s21
                if(ldebug.and.max(abs(xint),abs(xintzz)).gt.1d-10)then
                 write(6,*)('final xint '),xint,xintzz
                 write(6,*)('coupling matrix ')
                 call prntm2(bc(iprod),ncsf(jarg),ncsf(iarg),
     $                ncsf(jarg))
                end if
                ibcoff=iprod                                                  11d13s20
                if(ipss.eq.2)go to 22                                   1d26s21
               end if                                                   1d26s21
              end if                                                    1d26s21
             end do                                                     1d26s21
   22        continue                                                   1d26s21
            end if                                                         11d13s20
           end if                                                       1d21s21
          end if                                                        1d21s21
         iok=0                                                          10d2s20
        end if                                                          12d9s19
        jps=jps+ncsf(jarg)
       end do                                                           6d11s19
       ips=ips+ncsf(iarg)
      end do
      nhtmp=jhtmp-ihtmp                                                 7d11s19
      call dws_gsumf(bc(ihtmp),nhtmp)                                   7d11s19
      if(npass.ne.1)then                                                5d14s21
       call dws_gsumf(bc(ihtmpzz),nhtmp)                                5d14s21
      end if                                                            5d14s21
      jhtmp=ihtmp                                                       7d11s19
      jhtmpzz=ihtmpzz                                                   5d14s21
      npsx=nps                                                          12d24s19
      ntri=(npsx*(npsx+1))/2                                            9d4s19
      do i=1,ntri                                                       7d11s19
       hamps(i)=0d0                                                     7d11s19
      end do                                                            7d11s19
      if(npass.ne.1)then                                                5d14s21
       ihampszz=ibcoff                                                  5d14s21
       ibcoff=ihampszz+npsx*npsx                                        5d14s21
       call enough('pslzzcsf.  3',bc,ibc)
       do i=ihampszz,ibcoff-1                                           5d14s21
        bc(i)=0d0                                                       5d14s21
       end do                                                           5d14s21
       jhampszz=ihampszz-1                                              5d14s21
      end if                                                            5d14s21
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
           if(npass.ne.1)bc(jhampszz+itry)=bc(jhtmpzz)                  5d14s21
          end if                                                        7d11s19
          jhtmp=jhtmp+1                                                 7d11s19
          jhtmpzz=jhtmpzz+1                                             5d14s21
         end do                                                         7d11s19
        end do                                                          7d11s19
       else                                                             7d11s19
        if(ldebug)write(6,*)('adding diagonal '),hdps(if,1),hdps(if,2),
     $       if
        do i=1,ncsf(iarg)                                                6d11s19
         ipp=ips+i                                                        6d11s19
         itry=((ipp*(ipp-1))/2)+ipp                                     7d11s19
         hamps(itry)=hdps(if,1)                                         5d14s21
         if(npass.ne.1)bc(jhampszz+itry)=hdps(if,2)                     5d14s21
        end do                                                           6d11s19
       end if                                                           7d11s19
       jps=ips+ncsf(iarg)                                               6d19s19
       ic=iptr(2,nclop)+nclo*(ipsbase(2,if)-1)                           7d11s19
       io=iptr(4,nclop)+nopen*(ipsbase(3,if)-1)                          12d9s19
       do jf=if+1,npsf                                                  6d19s19
        ncloj=ipsbase(1,jf)                                             6d11s19
        nclopj=ncloj+1                                                  12d11s19
        jarg=ncloj+1-mdon                                               6d11s19
        nopenj=nec-2*ncloj
        jo=iptr(4,nclopj)+nopenj*(ipsbase(3,jf)-1)                       7d11s19
        jc=iptr(2,nclopj)+ncloj*(ipsbase(2,jf)-1)                        7d11s19
        call gcode(idorb(jc),ncloj,isorb(jo),nopenj,
     $       idorb(ic),nclo,isorb(io),nopen,icode,imap,nnot,nx,0,nab)
        if(nnot.gt.0)then                                               11d4s19
         do ik=1,ncsf(iarg)                                             12d11s19
          ikp=ik+ips                                                    12d11s19
          do ib=1,ncsf(jarg)                                            12d11s19
           jbp=ib+jps                                                   12d11s19
           jxx=max(jbp,ikp)                                             12d30s21
           jnn=min(jbp,ikp)                                             12d30s21
           itry=((jxx*(jxx-1))/2)+jnn                                   12d30s21
           hamps(itry)=bc(jhtmp)                                        12d11s19
           jhtmp=jhtmp+1                                                12d11s19
          end do                                                        12d11s19
         end do                                                         12d11s19
         if(npass.ne.1)then                                             5d14s21
          do ik=1,ncsf(iarg)                                             12d11s19
           ikp=ik+ips                                                    12d11s19
           do ib=1,ncsf(jarg)                                            12d11s19
            jbp=ib+jps                                                   12d11s19
            jxx=max(jbp,ikp)                                             12d30s21
            jnn=min(jbp,ikp)                                             12d30s21
            itry=((jxx*(jxx-1))/2)+jnn                                  12d30s21
            bc(jhampszz+itry)=bc(jhtmpzz)                               5d14s21
            jhtmpzz=jhtmpzz+1                                           5d14s21
           end do                                                        12d11s19
          end do                                                         12d11s19
         end if                                                         5d14s21
        end if                                                          12d11s19
        jps=jps+ncsf(jarg)
       end do                                                           7d11s19
       ips=ips+ncsf(iarg)
      end do
      call square(hamps,ips)                                            12d24s19
      if(npass.eq.1)then                                                5d14s21
       if(ldebug)then
        write(6,*)('lzz matrix is ps basis ')
        call prntm2(hamps,ips,ips,ips)                                  7d18s22
        icol=0
        do if=1,npsf
         nclo=ipsbase(1,if)                                               6d11s19
         nclop=nclo+1                                                     6d11s19
         iarg=nclo+1-mdon                                                 6d11s19
         nopen=nec-2*nclo
         iic=iptrbit(1,nclop,mysym)+ipsbase(2,if)-1                       1d21s21
         iio=iptrbit(2,nclop,mysym)+ipsbase(3,if)-1                       1d21s21
         irow=0
         do jf=1,if
          ncloj=ipsbase(1,jf)                                               6d11s19
          nclojp=ncloj+1                                                     6d11s19
          jarg=ncloj+1-mdon                                                 6d11s19
          nopenj=nec-2*ncloj
          jjc=iptrbit(1,nclojp,mysym)+ipsbase(2,jf)-1                       1d21s21
          jjo=iptrbit(2,nclojp,mysym)+ipsbase(3,jf)-1                       1d21s21
          do ik=1,ncsf(iarg)
           ikp=ips*(ik+icol-1)
           do ib=1,ncsf(jarg)
            ibp=ib+irow
            iad=ibp+ikp
            if(abs(hamps(iad)).gt.1d-10)then
             call dcbit(ibc(iic),norb,'iic')
             call dcbit(ibc(iio),norb,'iio')
             call dcbit(ibc(jjc),norb,'jjc')
             call dcbit(ibc(jjo),norb,'jjo')
             write(6,*)ib,ik,hamps(iad),ibp,ik+icol
            end if
           end do
          end do
          irow=irow+ncsf(jarg)
         end do
         icol=icol+ncsf(iarg)
        end do
       end if
      else                                                              5d14s21
       if(ldebug)then
        write(6,*)('ll matrix is ps basis ')
        call prntm2(hamps,ips,ips,ips)                                  7d18s22
       end if
      end if                                                            5d14s21
      if(npass.ne.1)then                                                5d14s21
       call square(bc(ihampszz),ips)                                    5d14s21
      else                                                              5d14s21
       ibcoff=ihtmp                                                     5d14s21
      end if                                                            5d14s21
      ieig=ibcoff
      ivec=ieig+ips
      isyr=ivec+ips*ips
      ibcoff=isyr+ips
      call enough('pslzzcsf.  4',bc,ibc)
      call diagx(ips,hamps,bc(ieig),bc(ivec),ibc(isyr),bc,ibc)          11d14s22
      if(ldebug)then
       write(6,*)('eigenvalues ')
       call prntm2(bc(ieig),1,ips,1)
       write(6,*)('eigenvectors ')
       call prntm2(bc(ivec),ips,ips,ips)
      end if
      if(lwrite)then
       if(npass.eq.1)then                                               1d2s20
        write(6,*)('looking for lambda '),lambdaci
       else                                                             1d2s20
        write(6,*)('looking for L '),lambdaci                           1d2s20
       end if                                                           1d2s20
      end if                                                            1d2s20
      nkeep=0                                                           12d24s19
      nbad=0                                                            12d30s19
      do i=0,nps-1
       if(npass.eq.1)then                                               1d2s20
        try=sqrt(abs(bc(ieig+i)))                                        12d24s19
        itry=nint(try)                                                   12d24s19
        delta=abs(bc(ieig+i)-dfloat(itry*itry))                          12d24s19
       else                                                             1d2s20
        try=sqrt(abs(bc(ieig+i))+0.25d0)-0.5d0                          12d31s19
        itry=nint(try)                                                   12d24s19
        delta=abs(bc(ieig+i)-dfloat(itry*(itry+1)))                     12d31s19
       end if                                                           1d2s20
       if(delta.gt.epssymci)nbad=nbad+1                                 9d19s23
       if(delta.gt.epssymci.or.itry.eq.lambdaci)then                    9d19s23
        jvec=ivec+nps*i                                                 12d24s19
        ikeep=ivec+nps*nkeep                                            12d24s19
        do j=0,nps-1                                                    12d24s19
         bc(ikeep+j)=bc(jvec+j)                                         12d24s19
        end do                                                          12d24s19
        nkeep=nkeep+1                                                   12d24s19
       end if                                                           12d24s19
      end do
      if(lwrite)then                                                    1d2s20
       write(6,*)('number of keepers = '),nkeep                          12d24s19
       if(nbad.gt.0)write(6,*)('of which '),nbad,(' are not pure')       12d30s19
      end if                                                            1d2s20
      if(nbad.gt.0)write(6,*)('nbad = '),nbad
      jvec=ivec-1                                                       12d24s19
      do j=1,nkeep*nps                                                  12d24s19
       hamps(j)=bc(jvec+j)                                              12d24s19
      end do                                                            12d24s19
      ibcoff=ieig                                                       5d14s21
      if(npass.ne.1)then                                                5d14s21
       itmp=ibcoff                                                      5d14s21
       ibcoff=itmp+nkeep*nps                                            5d14s21
       call enough('pslzzcsf.  5',bc,ibc)
       nrow=nps                                                         5d14s21
       do ipass=1,2                                                     5d14s21
        call dgemm('n','n',nrow,nkeep,nps,1d0,bc(ihampszz),nrow,        5d14s21
     $       hamps,nps,0d0,bc(itmp),nrow,                               5d14s21
     d' pslzzcsf.  2')
        do i=0,nkeep-1                                                  5d14s21
         do j=0,nrow-1                                                  5d14s21
          ji=itmp+j+nrow*i                                              5d14s21
          ij=ihampszz+i+nkeep*j                                         5d14s21
          bc(ij)=bc(ji)                                                 5d14s21
         end do                                                         5d14s21
        end do                                                          5d14s21
        nrow=nkeep                                                      5d14s21
       end do                                                           5d14s21
       ibcoff=itmp                                                      5d14s21
       ieig=ibcoff                                                      5d14s21
       ivec=ieig+nkeep                                                  5d14s21
       isymz=ivec+nkeep*nkeep                                           5d14s21
       ibcoff=isymz+nkeep                                               5d14s21
       call enough('pslzzcsf.  6',bc,ibc)
       call diagx(nkeep,bc(ihampszz),bc(ieig),bc(ivec),ibc(isymz),bc,   11d14s22
     $      ibc)                                                        11d14s22
        pick=dfloat(lz2*lz2)                                            8d23s21
       jvec=ivec                                                        5d14s21
       kvec=ivec                                                        5d14s21
       nkeepr=0                                                         5d14s21
       do i=0,nkeep-1                                                   5d14s21
        if(abs(bc(ieig+i)-pick).lt.1d-10)then                           8d23s21
         do j=0,nkeep-1                                                 5d14s21
          bc(jvec+j)=bc(kvec+j)                                         5d14s21
         end do                                                         5d14s21
         jvec=jvec+nkeep                                                5d14s21
         nkeepr=nkeepr+1                                                5d14s21
        end if                                                          5d14s21
        kvec=kvec+nkeep                                                 5d14s21
       end do                                                           5d14s21
       if(nkeepr.eq.0)then                                              5d14s21
        write(6,*)('we have no keepers !!! '),pick                           5d14s21
        call prntm2(bc(ieig),1,nkeep,1)
        stop                                                            5d14s21
       end if                                                           5d14s21
       itmp=ibcoff                                                      5d14s21
       ibcoff=itmp+nps*nkeepr                                           5d14s21
       call enough('pslzzcsf.  7',bc,ibc)
       call dgemm('n','n',nps,nkeepr,nkeep,1d0,hamps,nps,bc(ivec),nkeep,5d14s21
     $      0d0,bc(itmp),nps,                                           5d14s21
     d' pslzzcsf.  3')
       jtmp=itmp-1                                                      5d14s21
       do i=1,nps*nkeepr                                                5d14s21
        hamps(i)=bc(jtmp+i)                                             5d14s21
       end do                                                           5d14s21
       nkeep=nkeepr                                                     5d14s21
      end if                                                            5d14s21
      ibcoff=ibcoffo                                                    5d14s21
      xkeep=dfloat(nkeep)                                               5d14s21
      call dws_bcast(xkeep,1)                                           5d14s21
      nkeep=nint(xkeep)                                                 5d14s21
      call dws_bcast(hamps,nkeep*nps)                                   12d24s19
      return
      end
