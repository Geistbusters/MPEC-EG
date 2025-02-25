c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine pshamcsfcas(bc,ibc,npsf,hdps,ipsbase,nec,hamps,nps,    11d10s22
     $     ncsf,mdon,i2e,ih0a,multh,nvv,hiv,hvv,                        11d10s22
     $     nlzzu,npsk,veclzz,ism,irel,irefo,mdoo,mysym,iptrbit,ixw1,    4d12s21
     $     ixw2,norb)                                                   11d14s22
      implicit real*8 (a-h,o-z)                                         6d11s19
      include "common.store"                                            6d11s19
      external second                                                   7d30s19
      integer*8 ipack,itesta,itestb                                     4d12s21
      integer*2 ipack2(4)                                               11d4s19
      character*9 fname                                                 7d17s19
      integer*1 icode(64),imap(64),nab(2),isame(64),                    8d1s22
     $     nabi(2),nabj(2),nab1(2),nab2(2),itest(64,2)                  11d13s20
      equivalence (ipack,ipack2)                                        11d4s19
      dimension hdps(*),hamps(*),ncsf(*),nother(2),                     8d1s22
     $    i2e(*),ipsbase(3,*),ih0a(*),multh(8,8),test(10),hiv(*),hvv(*),12d24s19
     $     veclzz(*),ism(*),irel(*),irefo(*),iptrbit(2,mdoo+1,*),       4d12s21
     $     nab4(2,3)                                                    4d12s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
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
         call gandc4(ibc(jjc),ibc(jjo),ibc(iic),ibc(iio),nopenj,nopen,  11d13s20
     $        norb,nnot4,nab4,bc,ibc)                                   11d14s22
         if(nnot4.gt.0)then                                               11d4s19
          itmp=ibcoff                                                   12d9s19
          jhtmp=itmp
          nn=ncsf(jarg)*ncsf(iarg)                                      12d9s19
          nnm=nn-1                                                      9d20s20
          ibcoff=itmp+nn                                                12d9s19
          call enough('pshamcsfcas.  1',bc,ibc)
          do i=0,nnm                                                    9d20s20
           bc(jhtmp+i)=0d0                                              12d11s19
          end do                                                        12d11s19
          jhtmp=ibcoff                                                  12d11s19
          if(mod(loopit,mynprocg).eq.mynowprog)then                      12d9s19
           loopit=mynowprog                                             9d10s24
           if(nab4(1,1).eq.0.and.nnot4.gt.1)then                        11d13s20
            write(6,*)('ooops, nab4 is '),nab4
            stop
           end if
           ibcsav=ibcoff                                                9d21s20
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
                jga=irel(nab1(1))-1                                           11d13s20
                jsb=ism(nab1(2))                                              11d13s20
                jgb=irel(nab1(2))-1                                           11d13s20
                jsc=ism(nab2(1))                                              11d13s20
                jgc=irel(nab2(1))-1                                           11d13s20
                jsd=ism(nab2(2))                                              11d13s20
                jgd=irel(nab2(2))-1                                           11d13s20
                xint=getintcas(i2e,jsa,jsb,jsc,jsd,jga,jgb,jgc,jgd,     8d1s22
     $               irefo,bc,ibc)                                      11d15s22
                nqq=karg+mdon-1
                if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                 itesta=ibc(iic)                                             11d13s20
                 itestb=ibc(iio)                                             11d13s20
                 itestb=ibclr(itestb,nab1(1))                             8d1s22
                 itestb=ibclr(itestb,nab1(2))                             8d1s22
                 itesta=ibset(itesta,nab1(2))                             8d1s22
                 call gandc(ibc(jjc),ibc(jjo),itesta,itestb,nopenj,     8d1s22
     $                nopenk,jarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,  8d1s22
     $                iwpb1,iwpk1,ncsfmid1,bc,ibc)                      11d14s22
                 iprod=ibcoff                                               11d13s20
                 itmp1=iprod+ncsf(jarg)*ncsf(jarg)                          11d13s20
                 itmp2=itmp1+ncsf(jarg)*ncsf(karg)                          11d13s20
                 ibcoff=itmp2+ncsf(jarg)*ncsf(karg)                         11d13s20
                 call enough('pshamcsfcas.  2',bc,ibc)
                 call prodn(iwpb1,iwpk1,ncsf(jarg),ncsf(karg),ncsfmid1,   12d4s20
     $              bc(itmp1),bc,ibc,1d0,0d0)                           2d13s23
                 do k=0,ncsf(karg)-1                                        11d13s20
                  do j=0,ncsf(jarg)-1                                       11d13s20
                   jk=itmp1+j+ncsf(jarg)*k                                  11d13s20
                   kj=itmp2+k+ncsf(karg)*j                                  11d13s20
                   bc(kj)=bc(jk)                                            11d13s20
                  end do                                                    11d13s20
                 end do                                                     11d13s20
                 call dgemm('n','n',ncsf(jarg),ncsf(jarg),ncsf(karg),   8d1s22
     $                1d0,bc(itmp1),ncsf(jarg),bc(itmp2),ncsf(karg),0d0,8d1s22
     $            bc(iprod),ncsf(jarg),                                 11d13s20
     d' pshamcsfcas.  1')
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
               end if                                                    8d1s22
              end do                                                       11d13s20
             end if                                                     8d1s22
            end do                                                        11d13s20
           else if(nnot4.eq.2)then                                        11d13s20
            call gandc(ibc(jjc),ibc(jjo),ibc(iic),ibc(iio),nopenj,nopen,  11d13s20
     $         jarg,iarg,ncsf,norb,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,    11d13s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
            iprod=ibcoff                                                  11d13s20
            ibcoff=iprod+ncsf(jarg)*ncsf(iarg)                            11d13s20
            call enough('pshamcsfcas.  3',bc,ibc)
            call prodn(iwpb1,iwpk1,ncsf(jarg),ncsf(iarg),ncsfmid1,      12d4s20
     $           bc(iprod),bc,ibc,1d0,0d0)                              2d13s23
            do i=1,norb                                                   11d13s20
             itest(i,1)=0                                                   11d13s20
             itest(i,2)=0                                                    11d13s20
             if(btest(ibc(iic),i))itest(i,1)=2                           8d1s22
             if(btest(ibc(iio),i))itest(i,1)=1                           8d1s22
             if(btest(ibc(jjc),i))itest(i,2)=2                           8d1s22
             if(btest(ibc(jjo),i))itest(i,2)=1                           8d1s22
            end do                                                        11d13s20
            nok=0                                                         11d13s20
            do i=1,norb                                                   11d13s20
             ixn=min(itest(i,1),itest(i,2))
             if(ixn.gt.0)then                                             11d13s20
              nok=nok+1                                                   11d13s20
              itest(nok,1)=ixn                                            11d13s20
              itest(nok,2)=i                                              11d13s20
             end if                                                       11d13s20
            end do                                                        11d13s20
            is=ism(nab4(1,1))                                             11d13s20
            jga=irel(nab4(2,1))-1                                       11d27s20
            jgb=irel(nab4(1,1))-1                                       11d27s20
            ih0u=ih0a(is)+jga+irefo(is)*jgb                               11d13s20
            sum=bc(ih0u)                                                  11d13s20
            do i=1,nok                                                    11d13s20
             js=ism(itest(i,2))                                           11d13s20
             jg=irel(itest(i,2))-1                                        11d13s20
             xint=getintcas(i2e,is,is,js,js,jga,jgb,jg,jg,irefo,bc,ibc) 11d15s22
             if(itest(i,1).eq.2)then                                      11d13s20
              xintk=getintcas(i2e,is,js,is,js,jga,jg,jgb,jg,irefo,bc,   11d15s22
     $            ibc)                                                  11d15s22
              sum=sum+2d0*xint-xintk                                      11d13s20
             else                                                         11d13s20
              sum=sum+xint                                                11d13s20
             end if                                                       11d13s20
            end do                                                        11d13s20
            do i=0,ncsf(jarg)*ncsf(iarg)-1                                11d13s20
             bc(ihtmpgg+i)=sum*bc(iprod+i)                                11d13s20
            end do                                                        11d13s20
            ibcoff=iprod                                                11d13s20
            do i=1,nok                                                    11d13s20
             if(itest(i,1).eq.1)then                                      11d13s20
c     jj is bra, ii is ket.
c
c         J   I    K     thus  JK  KI
c     b   n+1 n    n+1        g   g   (ck|bc)
c     k   m   m+1  m+1         ck  bc
c     c   p   p    p-1
c
              itesta=ibc(jjc)                                              11d13s20
              itestb=ibc(jjo)                                              11d13s20
              nopenk=nopenj                                                11d13s20
c
c     anihilate common
c
              if(btest(itesta,itest(i,2)))then                             11d13s20
               itesta=ibclr(itesta,itest(i,2))                             11d13s20
               itestb=ibset(itestb,itest(i,2))                             11d13s20
               karg=jarg-1                                                11d13s20
               nopenk=nopenk+1                                             11d13s20
              else                                                         11d13s20
               itestb=ibclr(itestb,itest(i,2))                             11d13s20
               karg=jarg                                                  11d13s20
               nopenk=nopenk-1                                             11d13s20
              end if                                                       11d13s20
c
c     create ket
c
              if(btest(itestb,nab4(2,1)))then                           11d27s20
               itesta=ibset(itesta,nab4(2,1))                           11d27s20
               itestb=ibclr(itestb,nab4(2,1))                           11d27s20
               karg=karg+1                                                11d13s20
               nopenk=nopenk-1                                             11d13s20
              else                                                         11d13s20
               itestb=ibset(itestb,nab4(2,1))                           11d27s20
               nopenk=nopenk+1                                             11d13s20
              end if                                                       11d13s20
              nqq=karg+mdon-1
              if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
               call gandc(ibc(jjc),ibc(jjo),itesta,itestb,nopenj,nopenk,     11d13s20
     $         jarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,    11d13s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
               call gandc(itesta,itestb,ibc(iic),ibc(iio),nopenk,nopen,      11d13s20
     $         karg,iarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,    11d13s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
               if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),ncsfmid1,     11d13s20
     $            iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod,bc,ibc)        11d10s22
                jsa=ism(nab1(1))                                              11d13s20
                jga=irel(nab1(1))-1                                           11d13s20
                jsb=ism(nab1(2))                                              11d13s20
                jgb=irel(nab1(2))-1                                           11d13s20
                jsc=ism(nab2(1))                                              11d13s20
                jgc=irel(nab2(1))-1                                           11d13s20
                jsd=ism(nab2(2))                                              11d13s20
                jgd=irel(nab2(2))-1                                           11d13s20
                xint=getintcas(i2e,jsa,jsb,jsc,jsd,jga,jgb,jgc,jgd,     4d13s21
     $               irefo,bc,ibc)                                      11d15s22
                do j=0,ncsf(jarg)*ncsf(iarg)-1                                11d13s20
                 bc(ihtmpgg+j)=bc(ihtmpgg+j)+xint*bc(iprod+j)                11d13s20
                end do                                                       11d13s20
                ibcoff=iprod                                                 11d13s20
               else if(max(nnot1,nnot2).lt.2)then                        11d13s20
                write(6,*)('c:expecting 2s, but got otherwise')
                call dcbit(itesta,norb,'itesta')
                call dcbit(itestb,norb,'itestb')
                stop
               end if
              end if                                                    11d13s20
             end if                                                       11d13s20
            end do                                                        11d13s20
           else if(nnot4.ge.3)then                                        11d13s20
            if(nnot4.eq.3)then                                          1d26s21
             ipssx=1                                                    12d8s20
            else                                                        12d8s20
             ipssx=3                                                    12d18s20
            end if                                                      12d8s20
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
c
c     j is bra
             itesta=ibc(jjc)                                               11d13s20
             itestb=ibc(jjo)                                               11d13s20
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
              write(6,*)('bit not set for nab4(1,1) = '),nab4(1,iu1)       11d27s20
              stop 'nab4(1,1)'                                           11d27s20
             end if                                                        11d13s20
             if(btest(itestb,nab4(2,iu2)))then                          1d26s21
              itesta=ibset(itesta,nab4(2,iu2))                          1d26s21
              itestb=ibclr(itestb,nab4(2,iu2))                          1d26s21
              nopenk=nopenk-1                                              11d13s20
              karg=karg+1                                                  11d13s20
             else if(btest(itesta,nab4(2,iu2)))then                     1d26s21
              write(6,*)('already double in nab4(2,1) = '),nab4(2,iu2)  1d26s21
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
              if(nnot1.eq.2.and.nnot2.eq.2)then                          1d26s21
               call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),ncsfmid1,    1d26s21
     $            iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod,bc,ibc)        11d10s22
               jsa=ism(nab1(1))                                              11d13s20
               jga=irel(nab1(1))-1                                           11d13s20
               jsb=ism(nab1(2))                                              11d13s20
               jgb=irel(nab1(2))-1                                           11d13s20
               jsc=ism(nab2(1))                                              11d13s20
               jgc=irel(nab2(1))-1                                           11d13s20
               jsd=ism(nab2(2))                                              11d13s20
               jgd=irel(nab2(2))-1                                           11d13s20
               xint=getintcas(i2e,jsa,jsb,jsc,jsd,jga,jgb,jgc,jgd,irefo,11d15s22
     $              bc,ibc)                                             11d15s22
               if(nab1(1).eq.nab2(1).and.nab1(2).eq.nab2(2))             1d26s21
     $             xint=xint*0.5d0                                      1d26s21
               do i=0,ncsf(jarg)*ncsf(iarg)-1                                11d13s20
                bc(ihtmpgg+i)=bc(ihtmpgg+i)+bc(iprod+i)*xint            1d26s21
               end do                                                        11d13s20
               ibcoff=iprod                                                  11d13s20
               if(ipss.eq.2)go to 22                                     1d26s21
              end if                                                    1d26s21
             end if                                                     1d26s21
            end do                                                      1d26s21
   22       continue                                                    1d26s21
           end if                                                         11d13s20
           ibcoff=ibcsav                                                 12d9s19
          end if                                                         12d9s19
          loopit=loopit+1                                                12d9s19
         end if                                                          12d9s19
        end if                                                          12d9s19
        jps=jps+ncsf(jarg)
       end do                                                           6d11s19
       ips=ips+ncsf(iarg)
      end do
      nhtmp=jhtmp-ihtmp                                                 7d11s19
      call dws_gsumf(bc(ihtmp),nhtmp)                                   7d11s19
      jhtmp=ihtmp                                                       7d11s19
      npsx=nps+nvv                                                      9d4s19
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
     $        norb,nnot,nab4,bc,ibc)                                    11d14s22
        if(nnot.gt.0)then                                               11d4s19
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
      if(nlzzu.ne.0)then                                                12d24s19
       call square(hamps,nps)                                           12d24s19
       nrow=nps                                                         12d24s19
       do ipass=1,2                                                     12d24s19
        call dgemm('n','n',nrow,npsk,nps,1d0,hamps,nrow,veclzz,nps,     12d24s19
     $       0d0,bc(ihtmp),nrow,                                        12d24s19
     d' pshamcsfcas.  2')
        if(ipass.eq.1)then                                              12d24s19
         do i=0,npsk-1                                                   12d24s19
          do j=0,nrow-1                                                  12d24s19
           ji=ihtmp+j+nrow*i                                             12d24s19
           ij=i+1+npsk*j                                                 12d24s19
           hamps(ij)=bc(ji)                                              12d24s19
          end do                                                         12d24s19
         end do                                                          12d24s19
         nrow=npsk                                                       12d24s19
        else                                                            12d24s19
         ij=1                                                           12d24s19
         do i=0,npsk-1                                                  12d24s19
          do j=0,i                                                      12d24s19
           ji=ihtmp+j+nrow*i                                            12d24s19
           hamps(ij)=bc(ji)                                             12d24s19
           ij=ij+1                                                      12d24s19
          end do                                                        12d24s19
         end do                                                         12d24s19
        end if                                                          12d24s19
       end do                                                           12d24s19
       if(nvv.gt.0)then                                                 12d24s19
        call dgemm('n','n',nvv,npsk,nps,1d0,hiv,nvv,veclzz,nps,0d0,      12d24s19
     $      bc(ihtmp),nvv,                                              12d24s19
     d' pshamcsfcas.  3')
        jhtmp=ihtmp-1                                                    12d24s19
        do i=1,nvv*npsk                                                  12d24s19
         hiv(i)=bc(jhtmp+i)                                              12d24s19
        end do                                                           12d24s19
       end if                                                           12d24s19
       npsu=npsk                                                        12d24s19
       npsx=npsk+nvv                                                    12d24s19
       ntri=(npsx*(npsx+1))/2                                            9d4s19
      else                                                              12d24s19
       npsu=nps                                                         12d24s19
      end if                                                            12d24s19
      do i=1,nvv                                                        9d4s19
       ip=i+npsu                                                        12d24s19
       do j=1,npsu                                                      12d24s19
        itry=((ip*(ip-1))/2)+j                                          9d4s19
        iad=i+nvv*(j-1)                                                 9d4s19
        hamps(itry)=hiv(iad)                                            9d4s19
       end do                                                           9d4s19
       do j=i,nvv                                                       9d4s19
        jp=j+npsu                                                       12d24s19
        itry=((jp*(jp-1))/2)+ip                                         9d4s19
        iad=j+nvv*(i-1)                                                 9d4s19
        hamps(itry)=hvv(iad)                                            9d4s19
       end do                                                           9d4s19
      end do                                                            9d4s19
      jhtmp=ihtmp                                                       7d11s19
      do i=1,ntri                                                       7d11s19
       im=i-1                                                           7d11s19
       if(mod(im,mynprocg).eq.mynowprog)then                            7d11s19
        bc(jhtmp)=hamps(i)                                              7d11s19
        jhtmp=jhtmp+1                                                   7d11s19
       end if                                                           7d11s19
      end do                                                            7d11s19
      nhtmp=jhtmp-ihtmp                                                 7d11s19
      jhtmp=ihtmp-1                                                     7d11s19
      do i=1,nhtmp                                                      7d11s19
       hamps(i)=bc(jhtmp+i)                                             7d11s19
      end do                                                            7d11s19
      ibcoff=ihtmp                                                      7d11s19
      return
      end
