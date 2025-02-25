c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hccsfd12(vecx,ncsft,ibasis,ncsf,iptrbit,nfcn,ih0a,i2e, 7d28s22
     $     nrootz,mdon,nec,mysym,ixw1,ixw2,mdoo,nh0av,iorbf,bc,ibc,     5d18s23
     $     igoal)
c
c     1 and 2 particle densities
c     and orbital der term
c
      implicit real*8 (a-h,o-z)                                         7d11s19
      external second                                                   8d1s19
      logical ltest,lnew                                                11d16s20
      integer*1 nab(2),nab1(2),nab2(2)                                  7d28s22
      integer*8 ipack,itesta,itestb                                     11d16s20
      integer*2 ipack2(4)                                               12d1s19
      equivalence (ipack,ipack2)                                        12d1s19
      dimension vecx(ncsft,nrootz),ih0a(*),nh0av(*),                    7d28s22
     $    i2e(*),ibasis(3,*),ncsf(*),nother(2),                         7d28s22
     $     iptrbit(2,mdoo+1),nab4(2,3),itest(64,2),iorbf(*)             9d29s22
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      include "common.store"                                            7d11s19
      include "common.mrci"                                             6d19s19
      ltest=.false.
      if(ltest)write(6,*)('hi, my name is hccsf12')
      ngcode=0
      ips=0                                                             7d11s19
      loopit=0                                                          7d11s19
      do if=1,nfcn                                                      7d11s19
       nclo=ibasis(1,if)                                                7d11s19
       nclop=nclo+1                                                     7d11s19
       iarg=nclop-mdon                                                  7d11s19
       nopen=nec-2*nclo                                                 7d11s19
       ip=ips+1                                                         7d11s19
       iic=iptrbit(1,nclop)+ibasis(2,if)-1                              7d29s22
       iio=iptrbit(2,nclop)+ibasis(3,if)-1                              7d29s22
       jps=ips                                                          12d11s19
       do jf=if,nfcn                                                    12d9s19
        ncloj=ibasis(1,jf)                                              7d11s19
        if(iabs(ncloj-nclo).le.2)then                                   12d9s19
         nclopj=ncloj+1                                                  6d11s19
         nopenj=nec-2*ncloj
         jarg=ncloj+1-mdon                                               6d11s19
         jp=jps+1
         jjo=iptrbit(2,nclopj)+ibasis(3,jf)-1                           7d29s22
         jjc=iptrbit(1,nclopj)+ibasis(2,jf)-1                           7d29s22
         call gandc4(ibc(jjc),ibc(jjo),ibc(iic),ibc(iio),nopenj,nopen,  11d16s20
     $        norb,nnot4,nab4,bc,ibc)                                   11d14s22
         nnotu=nnot4                                                    7d28s22
         if(nnotu.gt.0)then                                             11d16s20
          if(mod(loopit,mynprocg).eq.mynowprog)then                      12d9s19
           loopit=mynowprog                                             9d10s24
           if(nnotu.gt.1.and.nab4(1,1).eq.0)then                        11d16s20
            write(6,*)('gandc4 failure! '),nnot4                        11d16s20
            write(6,*)nab4                                              11d16s20
            stop                                                        11d16s20
           end if                                                       11d16s20
           if(ltest)then
            write(6,*)('nnotu = '),nnotu
            call dcbit(ibc(jjc),norb,'jjc')                             7d28s22
            call dcbit(ibc(jjo),norb,'jjo')                             7d28s22
           end if
           ibcsav=ibcoff                                                1d17s23
           if(nnot4.eq.1)then                                             11d13s20
            isum=ibcoff                                                 7d28s22
            ibcoff=isum+nrootz                                          7d28s22
            call enough('hccsfd12.  1',bc,ibc)
            jsum=isum-1                                                 7d28s22
            do ir=1,nrootz                                              7d28s22
             bc(jsum+ir)=0d0                                            7d28s22
             do i=0,ncsf(iarg)-1                                        7d28s22
              bc(jsum+ir)=bc(jsum+ir)+vecx(ip+i,ir)**2                  7d29s22
             end do                                                     7d28s22
            end do                                                      7d28s22
            do i=1,norb                                                 7d28s22
             fact=0d0                                                   3d3s21
             if(btest(ibc(jjc),i))then                                  3d3s21
              fact=2d0                                                  3d3s21
             else if(btest(ibc(jjo),i))then                             3d3s21
              fact=1d0                                                  3d3s21
             end if                                                     3d3s21
             if(fact.ne.0d0)then                                        3d3s21
              jbs=ism(i)                                                3d3s21
              jbg=irel(i)-1                                             3d3s21
              jden=ih0a(jbs)+jbg*(1+nh0av(jbs))                         7d28s22
              do ir=0,nrootz-1                                          3d3s21
               sum=bc(isum+ir)                                          3d3s21
               bc(jden)=bc(jden)+sum*fact                               3d3s21
               jden=jden+nh0av(jbs)*nh0av(jbs)                          3d25s21
              end do                                                    3d3s21
c
c          J  K     j,i
c     cc:  2 -1  unrestricted
c     oo: .5       j ne i
c     co:    -1     ...
c
              do j=1,norb                                               7d28s22
               kbs=ism(j)                                               7d28s22
               kbg=irel(j)-1                                            7d28s22
               if(btest(ibc(jjc),j).or.btest(ibc(jjo),j))then           7d28s22
                ixint=igetint(i2e,kbs,kbs,jbs,jbs,kbg,kbg,jbg,jbg,      7d28s22
     $               nxint,.false.,bc,ibc)                              11d15s22
                ixintk=igetint(i2e,kbs,jbs,jbs,kbs,kbg,jbg,jbg,kbg,      7d28s22
     $               nxintk,.false.,bc,ibc)                             11d15s22
                if(btest(ibc(jjc),j).and.fact.gt.1.5d0)then             7d28s22
                 do ir=1,nrootz                                         7d28s22
                  bc(ixint)=bc(ixint)+2d0*bc(jsum+ir)                   7d28s22
                  bc(ixintk)=bc(ixintk)-bc(jsum+ir)                     7d28s22
                  ixint=ixint+nxint                                     7d28s22
                  ixintk=ixintk+nxintk                                  7d28s22
                 end do                                                 7d28s22
                else if(btest(ibc(jjc),j).and.fact.lt.1.5d0)then        7d29s22
                 do ir=1,nrootz                                         7d28s22
                  bc(ixint)=bc(ixint)+2d0*bc(jsum+ir)                   7d28s22
                  bc(ixintk)=bc(ixintk)-bc(jsum+ir)                     7d28s22
                  ixint=ixint+nxint                                     7d28s22
                  ixintk=ixintk+nxintk                                  7d28s22
                 end do                                                 7d28s22
                else if(btest(ibc(jjo),j).and.fact.lt.1.5d0.and.        7d29s22
     $                j.ne.i)then                                       7d29s22
                 do ir=1,nrootz                                         7d28s22
                  bc(ixint)=bc(ixint)+0.50*bc(jsum+ir)                   7d28s22
                  ixint=ixint+nxint                                     7d28s22
                 end do                                                 7d28s22
                end if                                                  7d28s22
               end if                                                   7d28s22
              end do                                                    7d28s22
             end if                                                     3d3s21
            end do                                                      3d3s21
            ibcoff=isum                                                 7d28s22
            do i1=1,norb-1                                              7d28s22
             if(btest(ibc(iio),i1))then                                 7d28s22
              nab1(2)=i1                                                7d28s22
              do i2=i1+1,norb                                           7d28s22
               if(btest(ibc(iio),i2))then                               7d28s22
                nab1(1)=i2                                              7d28s22
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
                ixint=igetint(i2e,jsa,jsb,jsc,jsd,jga,jgb,jgc,jgd,nxint,7d29s22
     $               .false.,bc,ibc)                                    11d15s22
                nqq=karg+mdon-1
                if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                 itestb=ibclr(itestb,nab1(1))                           7d28s22
                 itestb=ibclr(itestb,nab1(2))                           7d28s22
                 itesta=ibset(itesta,nab1(2))                           7d28s22
                 call gandc(ibc(jjc),ibc(jjo),itesta,itestb,nopenj,       12d6s20
     $               nopenk,jarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,   12d6s20
     $               iwpb1,iwpk1,ncsfmid1,bc,ibc)                       11d14s22
                 call gandc(itesta,itestb,ibc(iic),ibc(iio),nopenk,     7d28s22
     $                nopen,karg,iarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,   7d28s22
     $                iwpb2,iwpk2,ncsfmid2,bc,ibc)                      11d14s22
                 iprod=ibcoff                                               11d13s20
                 ibcoff=iprod+ncsf(jarg)*nrootz                           11d16s20
                 call enough('hccsfd12.  2',bc,ibc)
                 call genmatn(ncsf(jarg),ncsf(karg),ncsf(iarg),ncsfmid1,  12d4s20
     $               iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,vecx(ip,1),ncsft, 12d4s20
     $               nrootz,bc(iprod),0d0,bc,ibc)                       11d10s22
                 do ir=0,nrootz-1                                         11d16s20
                  irp=ir+1                                                11d16s20
                  do i=0,ncsf(jarg)-1                                         11d13s20
                   ii=iprod+i+ncsf(jarg)*ir                               11d16s20
                   bc(ii)=bc(ii)-vecx(ip+i,irp)                           11d16s20
                  end do                                                  11d16s20
                 end do                                                      11d13s20
                 do ir=0,nrootz-1                                         7d28s22
                  irp=ir+1                                                7d28s22
                  jprod=iprod+ncsf(jarg)*ir                               7d28s22
                  do i=0,ncsf(jarg)-1                                     7d28s22
                   bc(ixint)=bc(ixint)+bc(jprod+i)*vecx(jp+i,irp)         7d28s22
                  end do                                                  7d28s22
                  ixint=ixint+nxint                                       7d28s22
                 end do                                                   7d28s22
                 ibcoff=iprod                                                11d13s20
                else                                                        11d13s20
                 do ir=0,nrootz-1                                         11d16s20
                  irp=ir+1                                                11d16s20
                  do i=0,ncsf(jarg)-1                                        11d13s20
                   bc(ixint)=bc(ixint)-vecx(jp+i,irp)*vecx(ip+i,irp)      7d28s22
                  end do                                                  11d16s20
                  ixint=ixint+nxint                                       7d28s22
                 end do                                                     11d13s20
                end if                                                      11d13s20
               end if                                                   7d28s22
              end do                                                       11d13s20
             end if                                                     7d28s22
            end do                                                        11d13s20
           else if(nnot4.eq.2)then                                        11d13s20
            call gandc(ibc(jjc),ibc(jjo),ibc(iic),ibc(iio),nopenj,      11d16s20
     $            nopen,jarg,iarg,ncsf,norb,ixw1,ixw2,nnot1,nab1,iwpb1, 11d16s20
     $            iwpk1,ncsfmid1,bc,ibc)                                11d14s22
            call gandc(ibc(iic),ibc(iio),ibc(jjc),ibc(jjo),nopen,       11d16s20
     $            nopenj,iarg,jarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,iwpb2, 11d16s20
     $            iwpk2,ncsfmid2,bc,ibc)                                11d14s22
            iprod=ibcoff                                                  11d13s20
            jprod=iprod+ncsf(iarg)*nrootz                               11d16s20
            ibcoff=jprod+ncsf(jarg)*nrootz                              12d4s20
            call enough('hccsfd12.  3',bc,ibc)
            call xtimesn(ncsf(iarg),ncsf(jarg),ncsfmid2,nrootz,iwpb2,   12d4s20
     $            iwpk2,vecx(jp,1),ncsft,bc(iprod),ncsf(iarg),1d0,0d0,  11d15s22
     $           bc,ibc)                                                11d15s22
            do i=1,norb                                                   11d13s20
             itest(i,1)=0                                                   11d13s20
             itest(i,2)=0                                                    11d13s20
             if(btest(ibc(iic),i))itest(i,1)=2                           7d29s22
             if(btest(ibc(iio),i))itest(i,1)=1                           7d29s22
             if(btest(ibc(jjc),i))itest(i,2)=2                           7d29s22
             if(btest(ibc(jjo),i))itest(i,2)=1                           7d29s22
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
            ih0u=ih0a(is)+jga+nh0av(is)*jgb                             7d28s22
            nn=nh0av(is)*nh0av(is)                                      7d28s22
            do ir=0,nrootz-1                                            7d28s22
             irp=ir+1                                                   7d28s22
             iiprod=iprod+ncsf(iarg)*ir                                 7d28s22
             do i=0,ncsf(iarg)-1                                        7d28s22
              bc(ih0u)=bc(ih0u)+vecx(ip+i,irp)*bc(iiprod+i)*2d0         7d29s22
             end do                                                     7d28s22
             ih0u=ih0u+nn                                               7d28s22
            end do                                                      7d28s22
            do i=1,nok                                                    11d13s20
             js=ism(itest(i,2))                                           11d13s20
             jg=irel(itest(i,2))-1                                        11d13s20
             ixint=igetint(i2e,is,is,js,js,jga,jgb,jg,jg,nxint,.false., 11d15s22
     $            bc,ibc)                                               11d15s22
             if(itest(i,1).eq.2)then                                      11d13s20
              ixintk=igetint(i2e,is,js,is,js,jga,jg,jgb,jg,nxintk,      7d29s22
     $            .false.,bc,ibc)                                       11d15s22
              do ir=0,nrootz-1                                            7d28s22
               irp=ir+1                                                   7d28s22
               iiprod=iprod+ncsf(iarg)*ir                                 7d28s22
               do ii=0,ncsf(iarg)-1                                        7d28s22
                bc(ixint)=bc(ixint)+4d0*vecx(ip+ii,irp)*bc(iiprod+ii)   7d29s22
                bc(ixintk)=bc(ixintk)-2d0*vecx(ip+ii,irp)*bc(iiprod+ii) 7d29s22
               end do                                                     7d28s22
               ixint=ixint+nxint                                        7d28s22
               ixintk=ixintk+nxintk                                     7d28s22
              end do                                                      7d28s22
             else                                                         11d13s20
              do ir=0,nrootz-1                                            7d28s22
               irp=ir+1                                                   7d28s22
               iiprod=iprod+ncsf(iarg)*ir                                 7d28s22
               do ii=0,ncsf(iarg)-1                                        7d28s22
                bc(ixint)=bc(ixint)+vecx(ip+ii,irp)*bc(iiprod+ii)*2d0   7d29s22
               end do                                                     7d28s22
               ixint=ixint+nxint                                        7d28s22
              end do                                                      7d28s22
             end if                                                       11d13s20
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
                call enough('hccsfd12.  4',bc,ibc)
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
                ixint=igetint(i2e,jsa,jsb,jsc,jsd,jga,jgb,jgc,jgd,nxint,7d29s22
     $               .false.,bc,ibc)                                    11d15s22
                jxint=ixint                                             7d28s22
                do ir=0,nrootz-1                                        7d28s22
                 irp=ir+1                                               7d28s22
                 jjprod=jprod+ncsf(jarg)*ir                             7d28s22
                 do j=0,ncsf(jarg)-1                                    7d28s22
                  bc(jxint)=bc(jxint)+vecx(jp+j,irp)*bc(jjprod+j)*2d0   7d29s22
                 end do                                                 7d28s22
                 jxint=jxint+nxint                                      7d28s22
                end do                                                  7d28s22
                ibcoff=jprod                                                 11d13s20
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
c     j is bra and i is ket
            if(nnot4.eq.3)then                                          1d22s21
             ipssx=1                                                    1d22s21
            else                                                        1d22s21
             ipssx=3                                                    1d22s21
            end if                                                      1d22s21
            fact=0d0                                                    1d22s21
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
             itesta=ibc(jjc)                                                11d13s20
             itestb=ibc(jjo)                                                11d13s20
             if(btest(itesta,nab4(1,iu1)))then                          1d22s21
              itesta=ibclr(itesta,nab4(1,iu1))                          1d22s21
              itestb=ibset(itestb,nab4(1,iu1))                          1d22s21
              nopenk=nopenj+1                                           1d22s21
              karg=jarg-1                                               1d22s21
             else if(btest(itestb,nab4(1,iu1)))then                     1d22s21
              itestb=ibclr(itestb,nab4(1,iu1))                          1d22s21
              nopenk=nopenj-1                                               11d13s20
              karg=jarg                                                   11d13s20
             else                                                           11d13s20
              write(6,*)('bit not set for nab4(1,1) = '),nab4(1,iu1)       11d27s20
              stop 'nab4(1,1)'                                           11d27s20
             end if                                                         11d13s20
             if(btest(itestb,nab4(2,iu2)))then                          1d22s21
              itesta=ibset(itesta,nab4(2,iu2))                          1d22s21
              itestb=ibclr(itestb,nab4(2,iu2))                          1d22s21
              nopenk=nopenk-1                                               11d13s20
              karg=karg+1                                                   11d13s20
             else if(btest(itesta,nab4(2,iu2)))then                     1d22s21
              write(6,*)('already double in nab4(2,1) = '),nab4(2,iu2)  1d22s21
              stop 'nab4(2,1)'                                           11d27s20
             else                                                           11d13s20
              itestb=ibset(itestb,nab4(2,iu2))                          1d22s21
              nopenk=nopenk+1                                               11d13s20
             end if                                                         11d13s20
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
               call enough('hccsfd12.  5',bc,ibc)
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
               ixint=igetint(i2e,jsa,jsb,jsc,jsd,jga,jgb,jgc,jgd,nxint, 7d29s22
     $              .false.,bc,ibc)                                     11d15s22
               xint=2d0                                                 7d28s22
               if(nab1(1).eq.nab2(1).and.nab1(2).eq.nab2(2))              11d16s20
     $              xint=xint*0.5d0                                       11d16s20
               jxint=ixint                                              7d28s22
               do ir=0,nrootz-1                                         7d28s22
                irp=ir+1                                                7d28s22
                jjprod=jprod+ncsf(jarg)*ir                              7d28s22
                do j=0,ncsf(jarg)-1                                     7d28s22
                 bc(jxint)=bc(jxint)+xint*vecx(jp+j,irp)*bc(jjprod+j)   7d29s22
                end do                                                  7d28s22
                jxint=jxint+nxint                                       7d28s22
               end do                                                   7d28s22
               ibcoff=jprod                                                  11d13s20
               if(ipss.eq.2)go to 22                                    1d22s21
               fact=1d0                                                 1d22s21
              end if                                                    1d26s21
             end if                                                     1d22s21
            end do                                                      1d22s21
   22       continue                                                    1d22s21
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
      nwds=ncsft*nrootz                                                 7d11s19
      return
      end                                                               7d11s19
