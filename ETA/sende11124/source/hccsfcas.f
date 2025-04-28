c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hccsfcas(vecx,gx,ncsft,ibasis,ncsf,nfcn,ih0a,i2e,      8d1s22
     $     hdig,nrootz,mdon,nec,chc,ism,irel,irefo,                     8d1s22
     $     mysym,iptrbit,ixw1,ixw2,mdoo,norb,bc,ibc)                    11d10s22
      implicit real*8 (a-h,o-z)                                         7d11s19
      external second                                                   8d1s19
      integer*1 icode(64),imap(64),nab(2),isame(64),                    8d1s22
     $     nab1(2),nab2(2)                                              4d12s21
      integer*8 ipack,itesta,itestb                                     4d12s21
      integer*2 ipack2(4)                                               12d1s19
      equivalence (ipack,ipack2)                                        12d1s19
      dimension vecx(ncsft,nrootz),gx(ncsft,nrootz),ih0a(*),            8d1s22
     $    i2e(*),hdig(*),ibasis(3,*),ncsf(*),chc(*),nother(2),          8d1s22
     $     ism(*),irel(*),irefo(*),iptrbit(2,mdoo+1,*),nab4(2,3),       4d12s21
     $     itest(64,2)                                                  4d12s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      include "common.store"                                            7d11s19
      do i=1,nrootz                                                     7d11s19
       do j=1,ncsft                                                     7d11s19
        gx(j,i)=0d0                                                     7d11s19
       end do                                                           7d11s19
      end do                                                            7d11s19
      ips=0                                                             7d11s19
      loopit=0                                                          7d11s19
      do if=1,nfcn                                                      7d11s19
       nclo=ibasis(1,if)                                                7d11s19
       nclop=nclo+1                                                     7d11s19
       iarg=nclop-mdon                                                  7d11s19
       if(mynowprog.eq.0)then                                           7d11s19
        do i=1,nrootz                                                   7d11s19
         do j=1,ncsf(iarg)                                              7d11s19
          gx(j+ips,i)=gx(j+ips,i)+hdig(if)*vecx(j+ips,i)                7d11s19
         end do                                                         7d11s19
        end do                                                          7d11s19
       end if                                                           7d11s19
       nopen=nec-2*nclo                                                 7d11s19
       ip=ips+1                                                         7d11s19
       iic=iptrbit(1,nclop,mysym)+ibasis(2,if)-1                        11d16s20
       iio=iptrbit(2,nclop,mysym)+ibasis(3,if)-1                        11d16s20
       jps=ips                                                          12d11s19
       do jf=if,nfcn                                                    12d9s19
        ncloj=ibasis(1,jf)                                              7d11s19
        if(iabs(ncloj-nclo).le.2)then                                   12d9s19
         nclopj=ncloj+1                                                  6d11s19
         nopenj=nec-2*ncloj
         jarg=ncloj+1-mdon                                               6d11s19
         jjo=iptrbit(2,nclopj,mysym)+ibasis(3,jf)-1                     11d16s20
         jjc=iptrbit(1,nclopj,mysym)+ibasis(2,jf)-1                     11d16s20
         jp=jps+1
         call gandc4(ibc(jjc),ibc(jjo),ibc(iic),ibc(iio),nopenj,nopen,  11d16s20
     $        norb,nnot,nab4,bc,ibc)                                    11d14s22
         if(nnot.gt.0)then                                               11d4s19
          if(mod(loopit,mynprocg).eq.mynowprog)then                      12d9s19
           loopit=mynowprog                                             9d10s24
           if(nnot.gt.1.and.nab4(1,1).eq.0)then                         4d12s21
            write(6,*)('gandc4 failure! '),nnot                         4d12s21
            write(6,*)nab4                                              11d16s20
            stop                                                        11d16s20
           end if                                                       11d16s20
           ibcsav=ibcoff                                                  12d9s19
           jhtmpgg=ibcoff                                               11d16s20
           ibcoff=jhtmpgg+ncsf(jarg)*nrootz                             11d16s20
           if(if.ne.jf)then                                             11d16s20
            ihtmpgg=ibcoff                                              11d16s20
            ibcoff=ihtmpgg+ncsf(iarg)*nrootz                            11d16s20
           end if                                                       11d16s20
           call enough('hccsfcas.  1',bc,ibc)
           if(nnot.eq.1)then                                             11d13s20
            do i=0,nrootz-1                                             11d16s20
             do j=0,ncsf(iarg)-1                                        11d16s20
              ji=jhtmpgg+j+ncsf(jarg)*i                                    11d13s20
              bc(ji)=0d0                                                  11d13s20
             end do                                                       11d13s20
            end do                                                        11d13s20
            do i1=1,norb-1                                              8d1s22
             if(btest(ibc(iio),i1))then                                 8d1s22
              nab1(2)=i1                                                8d1s22
              do i2=i1+1,norb                                           8d1s22
               if(btest(ibc(iio),i2))then                               8d1s22
                nab1(1)=i2                                              8d1s22
                itesta=ibc(iic)                                             11d13s20
                itestb=ibc(iio)                                             11d13s20
                nopenk=nopen-2                                              11d13s20
                karg=iarg+1                                                 11d13s20
                nab2(1)=nab1(2)                                           12d6s20
                nab2(2)=nab1(1)                                           12d6s20
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
                 itestb=ibclr(itestb,nab1(1))                             8d1s22
                 itestb=ibclr(itestb,nab1(2))                             8d1s22
                 itesta=ibset(itesta,nab1(2))                             8d1s22
                 call gandc(ibc(jjc),ibc(jjo),itesta,itestb,nopenj,       12d6s20
     $               nopenk,jarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,   12d6s20
     $               iwpb1,iwpk1,ncsfmid1,bc,ibc)                       11d14s22
                 call gandc(itesta,itestb,ibc(iic),ibc(iio),nopenk,     8d1s22
     $                nopen,karg,iarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,   8d1s22
     $                iwpb2,iwpk2,ncsfmid2,bc,ibc)                      11d14s22
                 iprod=ibcoff                                               11d13s20
                 ibcoff=iprod+ncsf(jarg)*nrootz                           11d16s20
                 call enough('hccsfcas.  2',bc,ibc)
                 call genmatn(ncsf(jarg),ncsf(karg),ncsf(iarg),ncsfmid1,  12d4s20
     $               iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,vecx(ip,1),ncsft, 12d4s20
     $               nrootz,bc(iprod),0d0,bc,ibc)                       11d10s22
                 do ir=0,nrootz-1                                         11d16s20
                  irp=ir+1                                                11d16s20
                  do i=0,ncsf(jarg)-1                                          11d13s20
                   ii=iprod+i+ncsf(jarg)*ir                               11d16s20
                   bc(ii)=bc(ii)-vecx(ip+i,irp)                           11d16s20
                  end do                                                  11d16s20
                 end do                                                      11d13s20
                 do i=0,ncsf(jarg)*nrootz-1                               11d16s20
                  bc(jhtmpgg+i)=bc(jhtmpgg+i)+xint*bc(iprod+i)               11d13s20
                 end do                                                      11d13s20
                 ibcoff=iprod                                                11d13s20
                else                                                        11d13s20
                 do ir=0,nrootz-1                                         11d16s20
                  irp=ir+1                                                11d16s20
                  do i=0,ncsf(jarg)-1                                        11d13s20
                   ii=jhtmpgg+i+ncsf(jarg)*ir                             11d16s20
                   bc(ii)=bc(ii)-xint*vecx(ip+i,irp)                      11d16s20
                  end do                                                  11d16s20
                 end do                                                     11d13s20
                end if                                                      11d13s20
               end if                                                    8d1s22
              end do                                                       11d13s20
             end if                                                     8d1s22
            end do                                                        11d13s20
           else if(nnot.eq.2)then                                        11d13s20
            call gandc(ibc(jjc),ibc(jjo),ibc(iic),ibc(iio),nopenj,      11d16s20
     $            nopen,jarg,iarg,ncsf,norb,ixw1,ixw2,nnot1,nab1,iwpb1, 11d16s20
     $            iwpk1,ncsfmid1,bc,ibc)                                11d14s22
            call gandc(ibc(iic),ibc(iio),ibc(jjc),ibc(jjo),nopen,       11d16s20
     $            nopenj,iarg,jarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,iwpb2, 11d16s20
     $            iwpk2,ncsfmid2,bc,ibc)                                11d14s22
            iprod=ibcoff                                                  11d13s20
            jprod=iprod+ncsf(iarg)*nrootz                               11d16s20
            ibcoff=jprod+ncsf(jarg)*nrootz                              12d4s20
            call enough('hccsfcas.  3',bc,ibc)
            call xtimesn(ncsf(jarg),ncsf(iarg),ncsfmid1,nrootz,iwpb1,   12d4s20
     $            iwpk1,vecx(ip,1),ncsft,bc(jprod),ncsf(jarg),1d0,0d0,  11d10s22
     $           bc,ibc)                                                11d10s22
            call xtimesn(ncsf(iarg),ncsf(jarg),ncsfmid2,nrootz,iwpb2,   12d4s20
     $            iwpk2,vecx(jp,1),ncsft,bc(iprod),ncsf(iarg),1d0,0d0,  11d10s22
     $           bc,ibc)                                                11d10s22
            do i=1,norb                                                   11d13s20
             itest(i,1)=0                                                   11d13s20
             itest(i,2)=0                                                    11d13s20
             if(btest(ibc(iic),i))itest(i,1)=2                          8d1s22
             if(btest(ibc(iio),i))itest(i,1)=1                          8d1s22
             if(btest(ibc(jjc),i))itest(i,2)=2                          8d1s22
             if(btest(ibc(jjo),i))itest(i,2)=1                          8d1s22
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
            do i=0,ncsf(iarg)*nrootz-1                                  11d16s20
             bc(ihtmpgg+i)=sum*bc(iprod+i)                                11d13s20
            end do                                                        11d13s20
            do j=0,ncsf(jarg)*nrootz-1                                  11d16s20
             bc(jhtmpgg+j)=sum*bc(jprod+j)                                11d13s20
            end do                                                        11d13s20
            ibcoff=iprod                                                11d13s20
            do i=1,nok                                                    11d13s20
             if(itest(i,1).eq.1)then                                      11d13s20
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
               call gandc(ibc(jjc),ibc(jjo),itesta,itestb,nopenj,       11d16s20
     $              nopenk,jarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,    11d16s20
     $              iwpb1,iwpk1,ncsfmid1,bc,ibc)                        11d14s22
               call gandc(itesta,itestb,ibc(iic),ibc(iio),nopenk,nopen,      11d13s20
     $         karg,iarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,    11d13s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
               if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                jprod=ibcoff                                            11d16s20
                ibcoff=jprod+ncsf(jarg)*nrootz                          12d4s20
                call enough('hccsfcas.  4',bc,ibc)
                call genmatn(ncsf(jarg),ncsf(karg),ncsf(iarg),ncsfmid1, 12d4s20
     $                iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,vecx(ip,1),ncsft,12d4s20
     $                nrootz,bc(jprod),0d0,bc,ibc)                      11d10s22
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
                do j=0,ncsf(jarg)*nrootz-1                              11d16s20
                 bc(jhtmpgg+j)=bc(jhtmpgg+j)+xint*bc(jprod+j)                11d13s20
                end do                                                       11d13s20
                ibcoff=jprod                                                 11d13s20
                call gandc(ibc(iic),ibc(iio),itesta,itestb,nopen,       11d16s20
     $              nopenk,iarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,    11d16s20
     $              iwpb1,iwpk1,ncsfmid1,bc,ibc)                        11d14s22
                call gandc(itesta,itestb,ibc(jjc),ibc(jjo),nopenk,      11d16s20
     $                nopenj,karg,jarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,  11d16s20
     $                iwpb2,iwpk2,ncsfmid2,bc,ibc)                      11d14s22
                iprod=ibcoff                                            11d16s20
                ibcoff=iprod+ncsf(iarg)*nrootz                          12d4s20
                call enough('hccsfcas.  5',bc,ibc)
                call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),ncsfmid1, 12d4s20
     $                iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,vecx(jp,1),ncsft,12d4s20
     $                nrootz,bc(iprod),0d0,bc,ibc)                      11d10s22
                do j=0,ncsf(iarg)*nrootz-1                              11d16s20
                 bc(ihtmpgg+j)=bc(ihtmpgg+j)+xint*bc(iprod+j)           11d16s20
                end do                                                       11d13s20
                ibcoff=iprod                                            11d16s20
               else if(max(nnot1,nnot2).lt.2)then                        11d13s20
                write(6,*)('c:expecting 2s, but got otherwise')
                call dcbit(itesta,norb,'itesta')
                call dcbit(itestb,norb,'itestb')
                stop
               end if
              end if                                                    11d13s20
             end if                                                       11d13s20
            end do                                                        11d13s20
           else if(nnot.ge.3)then                                        11d13s20
c     j is bra and i is ket
            if(nnot.eq.3)then                                           1d22s21
             ipssx=1                                                    1d22s21
            else                                                        1d22s21
             ipssx=3                                                    1d22s21
            end if                                                      1d22s21
            do j=0,ncsf(jarg)*nrootz-1                                  8d2s22
             bc(jhtmpgg+j)=0d0                                          8d2s22
            end do                                                      8d2s22
            do i=0,ncsf(iarg)*nrootz-1                                  8d2s22
             bc(ihtmpgg+i)=0d0                                          8d2s22
            end do                                                      8d2s22
            do ipss=1,ipssx                                             1d22s21
             if(ipss.eq.1)then                                          1d22s21
              iu1=1                                                     1d22s21
              iu2=1                                                     1d22s21
             else if(ipss.eq.2)then                                     1d22s21
              iu1=1                                                     1d22s21
              iu2=2                                                     1d22s21
             else                                                       1d22s21
              iu1=2                                                     1d22s21
              iu2=1                                                     1d22s21
             end if                                                     12d8s20
             itesta=ibc(jjc)                                               11d13s20
             itestb=ibc(jjo)                                               11d13s20
             if(btest(itesta,nab4(1,iu1)))then                          1d22s21
              itesta=ibclr(itesta,nab4(1,iu1))                          1d22s21
              itestb=ibset(itestb,nab4(1,iu1))                          1d22s21
              nopenk=nopenj+1                                           1d22s21
              karg=jarg-1                                               1d22s21
             else if(btest(itestb,nab4(1,iu1)))then                     1d22s21
              itestb=ibclr(itestb,nab4(1,iu1))                          1d22s21
              nopenk=nopenj-1                                              11d13s20
              karg=jarg                                                  11d13s20
             else                                                          11d13s20
              write(6,*)('bit not set for nab4(1,1) = '),nab4(1,iu1)      11d27s20
              stop 'nab4(1,1)'                                          11d27s20
             end if                                                        11d13s20
             if(btest(itestb,nab4(2,iu2)))then                          1d22s21
              itesta=ibset(itesta,nab4(2,iu2))                          1d22s21
              itestb=ibclr(itestb,nab4(2,iu2))                          1d22s21
              nopenk=nopenk-1                                              11d13s20
              karg=karg+1                                                  11d13s20
             else if(btest(itesta,nab4(2,iu2)))then                     1d22s21
              write(6,*)('already double in nab4(2,1) = '),nab4(2,iu2)  1d22s21
              stop 'nab4(2,1)'                                          11d27s20
             else                                                          11d13s20
              itestb=ibset(itestb,nab4(2,iu2))                          1d22s21
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
              if(nnot1.eq.2.and.nnot2.eq.2)then                         1d22s21
               jprod=ibcoff                                               11d16s20
               ibcoff=jprod+ncsf(jarg)*nrootz                             12d4s20
               call enough('hccsfcas.  6',bc,ibc)
               call genmatn(ncsf(jarg),ncsf(karg),ncsf(iarg),ncsfmid1,    12d4s20
     $            iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,vecx(ip,1),ncsft,    12d4s20
     $            nrootz,bc(jprod),0d0,bc,ibc)                          11d10s22
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
               if(nab1(1).eq.nab2(1).and.nab1(2).eq.nab2(2))              11d16s20
     $              xint=xint*0.5d0                                       11d16s20
               do j=0,ncsf(jarg)*nrootz-1                                 11d16s20
                bc(jhtmpgg+j)=bc(jhtmpgg+j)+bc(jprod+j)*xint            8d2s22
               end do                                                        11d13s20
               ibcoff=jprod                                                  11d13s20
               call gandc(ibc(iic),ibc(iio),itesta,itestb,nopen,nopenk,   11d16s20
     $         iarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,    11d13s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
               call gandc(itesta,itestb,ibc(jjc),ibc(jjo),nopenk,       1d26s21
     $               nopenj,karg,jarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,   1d26s21
     $               iwpb2,iwpk2,ncsfmid2,bc,ibc)                       11d14s22
               iprod=ibcoff                                               11d16s20
               ibcoff=iprod+ncsf(iarg)*nrootz                             12d4s20
               call enough('hccsfcas.  7',bc,ibc)
               call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),ncsfmid1,    12d4s20
     $            iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,vecx(jp,1),ncsft,    12d4s20
     $            nrootz,bc(iprod),0d0,bc,ibc)                          11d10s22
               do i=0,ncsf(iarg)*nrootz-1                                 11d16s20
                bc(ihtmpgg+i)=bc(ihtmpgg+i)+bc(iprod+i)*xint            8d2s22
               end do                                                        11d13s20
               ibcoff=iprod                                                  11d13s20
               if(ipss.eq.2)go to 22                                    1d22s21
              end if                                                    1d26s21
             end if                                                     1d22s21
            end do                                                      1d22s21
   22       continue                                                    1d22s21
           end if                                                         11d13s20
           jht=jhtmpgg                                                   4d12s21
           do ir=0,nrootz-1                                              4d12s21
            irp=ir+1                                                     4d12s21
            do j=0,ncsf(jarg)-1                                          4d12s21
             gx(jp+j,irp)=gx(jp+j,irp)+bc(jht+j)                         4d12s21
            end do                                                       4d12s21
            jht=jht+ncsf(jarg)                                           4d12s21
           end do                                                        4d12s21
           if(jf.ne.if)then                                              12d9s19
            iht=ihtmpgg                                                  11d16s20
            do ir=0,nrootz-1                                             11d16s20
             irp=ir+1                                                    11d16s20
             do i=0,ncsf(iarg)-1                                         11d16s20
              gx(ip+i,irp)=gx(ip+i,irp)+bc(iht+i)                        11d16s20
             end do                                                      11d16s20
             iht=iht+ncsf(iarg)                                          11d16s20
            end do                                                       11d16s20
           end if                                                        12d9s19
           ibcoff=ibcsav                                                 12d9s19
          end if                                                         12d9s19
          loopit=loopit+1                                                12d9s19
         end if                                                          12d9s19
        end if                                                          12d9s19
        jps=jps+ncsf(jarg)
       end do                                                           6d11s19
       ips=ips+ncsf(iarg)
      end do
      nwds=ncsft*nrootz                                                 7d11s19
      call dws_gsumf(gx,nwds)
      do i=1,nrootz                                                     7d12s19
       sum=0d0                                                          7d12s19
       do j=1,ncsft                                                     7d12s19
        sum=sum+gx(j,i)*vecx(j,i)                                       7d12s19
       end do                                                           7d12s19
       chc(i)=sum                                                       7d12s19
      end do                                                            7d12s19
      return
      end                                                               7d11s19
