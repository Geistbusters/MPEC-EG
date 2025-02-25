c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcsid12(nhand,nff1,iff1,nff0,iff0,vec,ncsfv,ncsf,nec,  7d29s22
     $     mdon,mdoo,nsymb,multh,ixw1,ixw2,ih0av,nh0av,ionex,nvirt,     1d24s21
     $     ldebug,maxbx,bc,ibc,igoal)                                         11d10s22
      implicit real*8 (a-h,o-z)
      integer*1 nab1(2),nab2(2)
      integer*8 nhand(mdoo+1,nsymb,2),i1,in,i1c,i1o,i0c,i0o,itesta,     11d26s20
     $     itestb
      logical ldebug                                                    1d25s21
      dimension nff1(mdoo+1,nsymb,2),iff1(*),nff0(mdoo+1,3),            11d25s20
     $     iff0(*),vec(ncsfv,*),ncsf(*),multh(8,8),nab4(2,3),idorb1(32),11d25s20
     $     isorb1(32),idorb(32),isorb(32),itest(32,4),ih0av(*),ionex(*),11d26s20
     $     nvirt(*),iden1x(8,8),nden1x(8,8),nh0av(*)                    7d29s22
      include "common.store"                                            11d25s20
      include "common.mrci"                                             11d25s20
      common/kmfind/invk1(2,8,8,8,2)
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      if(ldebug)write(6,*)('hi, I am the new and improved hcsid12!')       1d25s21
      ibcoffo=ibcoff                                                    12d8s20
      norbx=norb+1                                                      11d25s20
      i1=1                                                              11d25s20
      ifirsttime=0                                                      11d26s20
      do isb=1,nsymb                                                    11d25s20
       isbv=multh(isb,isymmrci)
       do nclo1p=mdon+1,mdoo+1                                           11d25s20
        if(min(nvirt(isbv),nff1(nclo1p,isb,1)).gt.0)then                  11d25s20
         nclo1=nclo1p-1                                                  11d25s20
         nopen1=nec-nclo1*2                                             11d25s20
         iarg=nclo1p-mdon                                                11d25s20
         ncolt=nff1(nclo1p,isb,1)*ncsf(iarg)                            11d10s20
         call ilimts(1,ncolt,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)  12d7s20
         ncolh=ih+1-il                                                  12d7s20
         ighere=ibcoff                                                   11d25s20
         ibcoff=ighere+nvirt(isbv)*ncolh*nroot                          11d10s20
         i1c=nvirt(isbv)*nroot                                          11d10s20
         i1o=il                                                         12d7s20
         in=ih                                                          12d7s20
         call ddi_get(bc,ibc,nhand(nclo1p,isb,1),i1,i1c,i1o,in,         11d15s22
     $        bc(ighere))                                               11d15s22
         if(ldebug)then
          write(6,*)('my S vectors ')
          call prntm2(bc(ighere),nvirt(isbv)*nroot,ncolt,
     $         nvirt(isbv)*nroot)
         end if
         jghere=ighere                                                  12d2s20
         nwds=nvirt(isbv)*ncolh*nroot                                   3d22s21
         nn=ncsf(iarg)*nroot                                            11d26s20
         nnm=nn-1                                                       11d26s20
         idenh=ibcoff                                                   11d25s20
         ndenh=idenh+nn*irefo(isbv)                                     11d28s20
         ibcoff=ndenh+irefo(isbv)                                       11d28s20
         do isa=1,nsymb                                                 11d26s20
          isav=multh(isa,isbv)                                          11d26s20
          do isc=1,isa                                                  11d26s20
           if(isa.eq.isc)then                                           11d26s20
            nnv=(irefo(isa)*(irefo(isa)+1))/2                           11d26s20
           else                                                         11d26s20
            nnv=irefo(isa)*irefo(isc)                                   11d26s20
           end if                                                       11d26s20
           iscav=multh(isc,isav)                                        11d26s20
           iden1x(isc,isa)=ibcoff                                       11d26s20
           nden1x(isc,isa)=iden1x(isc,isa)+nn*nnv*irefo(iscav)          11d28s20
           ibcoff=nden1x(isc,isa)+nnv*irefo(iscav)                      11d28s20
          end do                                                        11d26s20
         end do                                                         11d26s20
         call enough('hcsid12.  1',bc,ibc)
         ibctop=ibcoff-1                                                11d26s20
         jff1=nff1(nclo1p,isb,2)                                        11d25s20
         ifcol=1                                                        12d7s20
         do if1=1,nff1(nclo1p,isb,1)                                    11d26s20
          itcol=ifcol+ncsf(iarg)-1                                      11d10s20
          if(ifcol.le.ih.and.itcol.ge.il)then                           12d7s20
           nnlow=max(0,il-ifcol)                                        12d7s20
           nnhgh=min(ncsf(iarg)-1,ih-ifcol)                             12d12s20
           nnuse=nnhgh+1-nnlow                                          12d7s20
           nnuser=nnuse*nroot                                           11d10s20
           itmpg=ibcoff                                                 12d7s20
           ibcoff=itmpg+nnuser*nvirt(isbv)                              11d10s20
           call enough('hcsid12.  2',bc,ibc)
           do iz=itmpg,ibcoff-1                                         7d29s22
            bc(iz)=0d0                                                  7d29s22
           end do                                                       7d29s22
           do iv=0,nvirt(isbv)-1                                        12d7s20
            do i=nnlow,nnhgh                                            12d7s20
             ip=ifcol+i-il                                              12d7s20
             im=i-nnlow                                                 12d7s20
             do ir=0,nroot-1
              iad1=itmpg+iv+nvirt(isbv)*(ir+nroot*im)                   7d29s22
              iad2=ighere+ir+nroot*(iv+nvirt(isbv)*ip)                  12d14s20
              bc(iad1)=bc(iad2)                                          7d29s22
             end do                                                     11d10s20
            end do                                                      12d7s20
           end do                                                       12d7s20
           do i=idenh,ibctop                                             11d26s20
            bc(i)=0d0                                                    11d26s20
           end do                                                        11d26s20
           i1c=iff1(jff1)                                                11d25s20
           ntest=popcnt(i1c)
           do i=1,norb                                                   11d25s20
            itest(i,1)=0                                                 11d25s20
           end do                                                        11d25s20
           ii=1                                                          11d25s20
           do i=1,norb                                                   11d25s20
            if(btest(i1c,i))then                                         11d25s20
             idorb1(ii)=i                                                11d25s20
             itest(i,1)=2                                                11d25s20
             ii=ii+1                                                     11d25s20
            end if                                                       11d25s20
           end do                                                        11d25s20
           jff1=jff1+1                                                   11d25s20
           i1o=iff1(jff1)                                                11d25s20
           i1o=ibset(i1o,norbx)                                          11d25s20
           ii=1                                                          11d25s20
           do i=1,norb                                                   11d25s20
            if(btest(i1o,i))then                                         11d25s20
             isorb1(ii)=i                                                11d25s20
             itest(i,1)=1                                                11d25s20
             ii=ii+1                                                     11d25s20
            end if                                                       11d25s20
           end do                                                        11d25s20
           jff1=jff1+1                                                   11d25s20
           ngot=0
           do nclop=max(mdon+1,nclo1p-2),min(mdoo+1,nclo1p+2)             11d25s20
            nclo=nclop-1                                                  11d25s20
            nopen=nec-nclo*2                                             11d25s20
            jarg=nclop-mdon                                               11d25s20
            jff0=nff0(nclop,2)                                           11d25s20
            jvs=nff0(nclop,3)                                            11d26s20
            do jf=1,nff0(nclop,1)                                        11d26s20
             i0c=iff0(jff0)                                              11d25s20
             jff0=jff0+1                                                 11d25s20
             i0o=iff0(jff0)                                              11d25s20
             jff0=jff0+1                                                 11d25s20
             call gandc4(i1c,i1o,i0c,i0o,nopen1,nopen,norbx,nnot,nab4,  11d14s22
     $            bc,ibc)                                               11d14s22
             if(nnot.ne.0)then
              if(nab4(1,1).eq.0)then                                     11d27s20
               write(6,*)('gandc4 falure!!! '),nab4
               call dcbit(i1c,norbx,'i1c')
               call dcbit(i0c,norbx,'i0c')
               call dcbit(i1o,norbx,'i1o')
               call dcbit(i0o,norbx,'i0o')
               call dws_synca                                            11d27s20
               call dws_finalize                                         11d27s20
              end if                                                     11d27s20
             end if
             if(nnot.eq.2)then                                           11d25s20
              call gandc(i1c,i1o,i0c,i0o,nopen1,nopen,iarg,jarg,ncsf,    11d25s20
     $            norbx,ixw1,ixw2,nnot1,nab1,iwpb,iwpk,ncsfmid,bc,ibc)  11d14s22
              iprod=ibcoff                                                  11d13s20
              ibcoff=iprod+ncsf(iarg)*nroot                              11d26s20
              call enough('hcsid12.  3',bc,ibc)
              call xtimesn(ncsf(iarg),ncsf(jarg),ncsfmid,                11d26s20
     $            nroot,iwpb,iwpk,vec(jvs,1),ncsfv,bc(iprod),ncsf(iarg),12d7s20
     $             1d0,0d0,bc,ibc)                                      11d10s22
              ii=1                                                          11d25s20
              do i=1,norb                                                   11d25s20
               itest(i,2)=0                                              11d25s20
               if(btest(i0c,i))then                                         11d25s20
                idorb(ii)=i                                                11d25s20
                ii=ii+1                                                     11d25s20
                itest(i,2)=2                                             11d25s20
               end if                                                       11d25s20
              end do                                                     11d25s20
              ii=1                                                          11d25s20
              do i=1,norb                                                   11d25s20
               if(btest(i0o,i))then                                         11d25s20
                isorb(ii)=i                                                11d25s20
                ii=ii+1                                                     11d25s20
                itest(i,2)=1                                             11d25s20
               end if                                                       11d25s20
              end do                                                        11d25s20
              nok=0                                                      11d25s20
              do i=1,norb                                                11d25s20
               ixn=min(itest(i,1),itest(i,2))                            11d25s20
               if(ixn.gt.0)then                                          11d25s20
                nok=nok+1                                                11d25s20
                itest(nok,3)=ixn                                         11d25s20
                itest(nok,4)=i                                           11d25s20
               end if                                                    11d25s20
              end do                                                     11d25s20
              jgb=irel(nab4(2,1))-1                                      11d27s20
              ibc(ndenh+jgb)=1                                           11d28s20
              jdenh=idenh+nn*jgb                                         11d26s20
              do i=0,nnm                                                 11d26s20
               bc(jdenh+i)=bc(jdenh+i)+bc(iprod+i)                       11d26s20
              end do                                                     11d26s20
              do i=1,nok                                                    11d13s20
               js=ism(itest(i,4))                                           11d13s20
               jg=irel(itest(i,4))-1                                        11d13s20
               icolj=((jg*(jg+1)/2))+jg                                  11d26s20
               nnvj=(irefo(js)*(irefo(js)+1))/2                          11d26s20
               iicol=icolj+nnvj*jgb                                      11d28s20
               jden=iden1x(js,js)+nn*iicol                               11d28s20
               ibc(nden1x(js,js)+iicol)=1                                11d28s20
c     is=isbv,jga runs over virts jg,jg,jgb,iv
               if(itest(i,3).eq.2)then                                      11d13s20
                if(isbv.eq.js)then                                       11d26s20
                 ix=max(jg,jgb)                                          11d26s20
                 in=min(jg,jgb)                                          11d26s20
                 icolk=((ix*(ix+1))/2)+in                                11d26s20
                 iicol=icolk+nnvj*jg                                     12d1s20
                 kden=iden1x(js,js)+nn*iicol                             12d2s20
                 ibc(nden1x(js,js)+iicol)=1                                11d28s20
                else                                                     11d26s20
                 nnvk=irefo(js)*irefo(isbv)                              11d26s20
                 if(isbv.gt.js)then                                      11d26s20
                  icolk=jg+irefo(js)*jgb                                  11d26s20
                  iicol=icolk+nnvk*jg                                    11d28s20
                  kden=iden1x(js,isbv)+nn*iicol                          11d28s20
                  ibc(nden1x(js,isbv)+iicol)=1                           11d28s20
                 else                                                    11d26s20
                  icolk=jgb+irefo(isbv)*jg                               11d26s20
                  iicol=icolk+nnvk*jg                                    11d28s20
                  kden=iden1x(isbv,js)+nn*iicol                          12d1s20
                  ibc(nden1x(isbv,js)+iicol)=1                           11d28s20
                 end if                                                  11d26s20
                end if                                                   11d26s20
c     jg,jgb,jg,iv
                do ii=0,nnm                                               11d26s20
                 bc(jden+ii)=bc(jden+ii)+2d0*bc(iprod+ii)                   11d26s20
                 bc(kden+ii)=bc(kden+ii)-bc(iprod+ii)                       11d26s20
                end do                                                   11d26s20
               else                                                         11d13s20
                do ii=0,nnm                                               11d26s20
                 bc(jden+ii)=bc(jden+ii)+bc(iprod+ii)                       11d26s20
                end do                                                   11d26s20
               end if                                                       11d13s20
              end do                                                        11d13s20
              ibcoff=iprod
              do i=1,nok                                                    11d13s20
               if(itest(i,3).eq.1)then                                      11d13s20
                itesta=i1c                                               11d26s20
                itestb=i1o                                               11d26s20
                nopenk=nopen1                                                11d13s20
c
c     anihilate common
c
                if(btest(itesta,itest(i,4)))then                             11d13s20
                 itesta=ibclr(itesta,itest(i,4))                             11d13s20
                 itestb=ibset(itestb,itest(i,4))                             11d13s20
                 karg=iarg-1                                                11d13s20
                 nopenk=nopenk+1                                             11d13s20
                else                                                         11d13s20
                 itestb=ibclr(itestb,itest(i,4))                             11d13s20
                 karg=iarg                                                  11d13s20
                 nopenk=nopenk-1                                             11d13s20
                end if                                                       11d13s20
c
c     create ket
c
                if(btest(itestb,nab4(2,1)))then                          11d27s20
                 itesta=ibset(itesta,nab4(2,1))                          11d27s20
                 itestb=ibclr(itestb,nab4(2,1))                          11d27s20
                 karg=karg+1                                                11d13s20
                 nopenk=nopenk-1                                             11d13s20
                else                                                         11d13s20
                 itestb=ibset(itestb,nab4(2,1))                          11d27s20
                 nopenk=nopenk+1                                             11d13s20
                end if                                                       11d13s20
                nqq=karg+mdon-1
                if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                 call gandc(i1c,i1o,itesta,itestb,nopen1,nopenk,         11d26s20
     $          iarg,karg,ncsf,norbx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,    11d13s20
     $          ncsfmid1,bc,ibc)                                        11d14s22
                 call gandc(itesta,itestb,i0c,i0o,nopenk,nopen,          11d26s20
     $          karg,jarg,ncsf,norbx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,  11d27s20
     $               ncsfmid2,bc,ibc)                                   11d14s22
                 if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                  if(nab1(1).gt.nab1(2))then                             11d27s20
                   icpy=nab1(1)                                          11d27s20
                   nab1(1)=nab1(2)                                       11d27s20
                   nab1(2)=icpy                                          11d27s20
                  end if                                                 11d26s20
                  jsa=ism(nab1(1))                                              11d13s20
                  jga=irel(nab1(1))-1                                           11d13s20
                  jsb=ism(nab1(2))                                              11d13s20
                  jgb=irel(nab1(2))-1                                           11d13s20
                  jsd=ism(nab2(2))                                              11d13s20
                  jgd=irel(nab2(2))-1                                           11d13s20
                  if(jsa.eq.jsb)then                                     11d26s20
                   icolk=((jgb*(jgb+1))/2)+jga                           11d27s20
                   nnvk=(irefo(jsa)*(irefo(jsa)+1))/2                    11d27s20
                  else                                                   11d26s20
                   icolk=jga+irefo(jsa)*jgb                              11d27s20
                   nnvk=irefo(jsa)*irefo(jsb)                            11d26s20
                  end if                                                 11d26s20
                  jdenk=iden1x(jsa,jsb)+nn*(icolk+nnvk*jgd)              11d26s20
                  ibc(nden1x(jsa,jsb)+icolk+nnvk*jgd)=1                  12d1s20
                  call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),        12d7s20
     $                 ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,       12d7s20
     $                 vec(jvs,1),ncsfv,nroot,bc(jdenk),1d0,bc,ibc)     11d10s22
                 end if
                end if                                                    11d13s20
               end if                                                       11d13s20
              end do                                                        11d13s20
             else if(nnot.ge.3)then                                      11d25s20
              if(nnot.eq.3)then                                         1d25s21
               ipssx=1                                                   1d25s21
              else                                                      1d25s21
               ipssx=3                                                   1d25s21
              end if                                                    1d25s21
              do ipss=1,ipssx                                           1d25s21
               if(ipss.eq.1)then                                        1d25s21
                iu1=1                                                   1d25s21
                iu2=1                                                   1d25s21
               else if(ipss.eq.2)then                                   1d25s21
                iu1=1                                                   1d25s21
                iu2=2                                                   1d25s21
               else                                                     1d25s21
                iu1=2                                                   1d25s21
                iu2=1                                                   1d25s21
               end if                                                   1d25s21
c     sing is bra
               itesta=i1c                                                 11d26s20
               itestb=i1o                                                 11d26s20
               if(btest(itesta,nab4(1,iu1)))then                        1d25s21
                itesta=ibclr(itesta,nab4(1,iu1))                        1d25s21
                itestb=ibset(itestb,nab4(1,iu1))                        1d25s21
                nopenk=nopen1+1                                              11d13s20
                karg=iarg-1                                                  11d13s20
               else if(btest(itestb,nab4(1,iu1)))then                   1d25s21
                itestb=ibclr(itestb,nab4(1,iu1))                        1d25s21
                nopenk=nopen1-1                                              11d13s20
                karg=iarg                                                  11d13s20
               end if                                                        11d13s20
               if(btest(itestb,nab4(2,iu2)))then                        1d25s21
                itesta=ibset(itesta,nab4(2,iu2))                        1d25s21
                itestb=ibclr(itestb,nab4(2,iu2))                        1d25s21
                nopenk=nopenk-1                                              11d13s20
                karg=karg+1                                                  11d13s20
               else                                                          11d13s20
                itestb=ibset(itestb,nab4(2,iu2))                        1d25s21
                nopenk=nopenk+1                                              11d13s20
               end if                                                        11d13s20
               call gandc(i1c,i1o,itesta,itestb,nopen1,nopenk,            11d26s20
     $         iarg,karg,ncsf,norbx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,   11d27s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
               call gandc(itesta,itestb,i0c,i0o,nopenk,nopen,             11d26s20
     $         karg,jarg,ncsf,norbx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,   11d27s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
               if(nnot1.eq.2.and.nnot2.eq.2)then                        1d25s21
                if(nab1(1).gt.nab1(2))then                                 11d26s20
                 icpy=nab1(1)                                              11d26s20
                 nab1(1)=nab1(2)                                           11d26s20
                 nab1(2)=icpy                                              11d26s20
                end if                                                     11d26s20
                jsa=ism(nab1(1))                                              11d13s20
                jga=irel(nab1(1))-1                                           11d13s20
                jsb=ism(nab1(2))                                              11d13s20
                jgb=irel(nab1(2))-1                                           11d13s20
                jsd=ism(nab2(2))                                           11d27s20
                jgd=irel(nab2(2))-1                                        11d27s20
                if(jsa.eq.jsb)then                                         11d27s20
                 icolk=((jgb*(jgb+1))/2)+jga                               11d27s20
                 nnvk=(irefo(jsa)*(irefo(jsa)+1))/2                        11d27s20
                else                                                       11d26s20
                 icolk=jga+irefo(jsa)*jgb                                  11d27s20
                 nnvk=irefo(jsa)*irefo(jsb)                                11d27s20
                end if                                                     11d26s20
                iicol=icolk+nnvk*jgd                                       11d28s20
                kden=iden1x(jsa,jsb)+nn*iicol                              11d28s20
                ibc(nden1x(jsa,jsb)+iicol)=1                               11d28s20
                call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),ncsfmid1,    11d26s20
     $            iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,vec(jvs,1),ncsfv,    12d7s20
     $             nroot,bc(kden),1d0,bc,ibc)                           11d10s22
                if(ipss.eq.2)go to 1                                    1d25s21
               end if                                                   1d25s21
              end do                                                    1d25s21
    1         continue                                                  1d25s21
             end if                                                      11d25s20
             jvs=jvs+ncsf(jarg)                                          11d26s20
            end do                                                       11d25s20
           end do                                                        11d25s20
           factg=0d0                                                    12d7s20
           if(irefo(isbv).gt.0)then
            nok=0                                                        12d2s20
            jdenh=idenh                                                  12d2s20
            ikeep=ibcoff                                                 12d2s20
            ibcoff=ikeep+irefo(isbv)                                     12d2s20
            call enough('hcsid12.  4',bc,ibc)
            do i=0,irefo(isbv)-1                                        3d19s21
             if(ibc(ndenh+i).gt.0)then                                  3d19s21
             iad=idenh+nn*i
              ibc(ikeep+nok)=i                                           12d2s20
              do ir=0,nroot-1                                           12d12s20
               do j=nnlow,nnhgh                                          12d12s20
                bc(jdenh)=bc(iad+j)                                      12d12s20
                jdenh=jdenh+1                                            12d12s20
               end do                                                   12d12s20
               iad=iad+ncsf(iarg)                                       12d12s20
              end do                                                     12d2s20
              nok=nok+1                                                  12d2s20
             end if                                                      12d2s20
            end do
            if(nok.gt.0)then                                             12d2s20
c     denh is ncsf(fractional),ir,irefo(fractional)
             nkr=nok*nroot
             itmp=ibcoff                                                 12d2s20
             ibcoff=itmp+nvirt(isbv)*nkr                                7d29s22
             call enough('hcsid12.  5',bc,ibc)
             call dgemm('n','n',nvirt(isbv),nkr,nnuse,2d0,bc(itmpg),    7d11s23
     $            nvirt(isbv),bc(idenh),nnuse,0d0,bc(itmp),nvirt(isbv), 7d29s22
     d' hcsid12.  1')
             do i=0,nok-1                                               7d29s22
              ip=ibc(ikeep+i)                                           12d2s20
              do ir=0,nroot-1                                            7d29s22
               iad1=itmp+nvirt(isbv)*(ir+nroot*i)                       7d29s22
               iad2=ih0av(isbv)+irefo(isbv)+nh0av(isbv)*(ip             7d29s22
     $              +nh0av(isbv)*ir)                                    7d29s22
               do iv=0,nvirt(isbv)-1                                       12d2s20
                bc(iad2+iv)=bc(iad2+iv)+bc(iad1+iv)                     7d29s22
               end do                                                   7d29s22
              end do                                                     12d2s20
             end do                                                      12d2s20
            end if                                                       12d2s20
            ibcoff=ikeep                                                 12d2s20
           end if
           do isa=1,nsymb                                                11d26s20
            isav=multh(isa,isbv)                                         11d26s20
            do isc=1,isa
             isavc=multh(isav,isc)                                       11d26s20
             if(isc.eq.isa)then
              nnv=(irefo(isc)*(irefo(isc)+1))/2
             else
              nnv=irefo(isc)*irefo(isa)
             end if
             ncol=nnv*irefo(isavc)
             if(ncol.gt.0)then
              i2eu=invk1(1,isc,isa,isavc,2)                             12d7s20
              icase=invk1(2,isc,isa,isavc,2)                             12d7s20
              nok=0                                                     12d7s20
              do i=0,ncol-1
               if(ibc(nden1x(isc,isa)+i).ne.0)then                      3d22s21
                ibc(nden1x(isc,isa)+nok)=i                              3d22s21
                nok=nok+1                                               3d22s21
               end if
              end do
              if(nok.gt.0)then                                          12d7s20
               nkr=nok*nroot                                            7d29s22
               itmpi=ibcoff                                             12d7s20
               ibcoff=itmpi+nkr*nvirt(isbv)                             7d29s22
               call enough('hcsid12.  6',bc,ibc)
               jden=iden1x(isc,isa)                                     12d7s20
               jtmpi=itmpi                                              12d7s20
               do i=0,nok-1                                             3d22s21
                icol=ibc(nden1x(isc,isa)+i)                             3d22s21
                ifrm=iden1x(isc,isa)+nn*icol                            12d7s20
                do ir=0,nroot-1                                         12d12s20
                 do j=nnlow,nnhgh                                       12d12s20
                  bc(jden)=bc(ifrm+j)                                   12d12s20
                  jden=jden+1                                           12d12s20
                 end do                                                 12d12s20
                 ifrm=ifrm+ncsf(iarg)                                   12d12s20
                end do                                                  12d12s20
               end do                                                   7d29s22
               call dgemm('n','n',nvirt(isbv),nkr,nnuse,1d0,bc(itmpg),  7d29s22
     $              nvirt(isbv),bc(iden1x(isc,isa)),nnuse,0d0,          7d29s22
     $              bc(itmpi),nvirt(isbv),                              7d29s22
     d' hcsid12.  2')
               nnn=irefo(isa)*irefo(isc)*irefo(isavc)                   7d29s22
               do i=0,nok-1                                             7d29s22
                icol=ibc(nden1x(isc,isa)+i)                             7d29s22
                if(icase.eq.1)then                                      12d7s20
                 icoln=icol                                             12d7s20
                else                                                    12d7s20
                 id=icol/nnv                                            12d7s20
                 ileft=icol-id*nnv                                      12d7s20
                 ia=ileft/irefo(isc)                                    1d25s21
                 ic=ileft-ia*irefo(isc)                                 1d25s21
                 icoln=ia+irefo(isa)*(ic+irefo(isc)*id)                 1d25s21
                 ittry=ic+irefo(isc)*(ia+irefo(isa)*id)                 1d25s21
                end if                                                  12d7s20
                do ir=0,nroot-1                                         7d29s22
                 jint=ionex(i2eu)+nvirt(isbv)*(icoln+nnn*ir)            7d29s22
                 do iv=0,nvirt(isbv)-1                                  7d29s22
                  bc(jint+iv)=bc(jint+iv)+bc(jtmpi+iv)                  7d29s22
                 end do                                                 7d29s22
                 jtmpi=jtmpi+nvirt(isbv)                                7d29s22
                end do                                                  7d29s22
               end do                                                   12d7s20
              end if                                                    12d7s20
             end if
            end do
           end do                                                        11d26s20
          else                                                          12d8s20
           jff1=jff1+2                                                  12d8s20
          end if                                                        12d7s20
          ifcol=itcol+1                                                 12d7s20
         end do                                                         11d25s20
         ibcoff=idenh                                                   3d22s21
        end if                                                          11d25s20
       end do                                                            11d25s20
      end do                                                            11d26s20
      return
      end
