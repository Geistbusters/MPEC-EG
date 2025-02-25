c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcsi(nhand,nff1,iff1,nff0,iff0,vec,ncsfv,ncsf,nec,     11d25s20
     $     mdon,mdoo,nsymb,multh,ixw1,ixw2,ih0av,nh0av,ionex,nvirt,     1d24s21
     $     izero,ldebug,maxbx,ionext,nwiacc,bc,ibc)                     11d10s22
      implicit real*8 (a-h,o-z)
      integer*1 nab1(2),nab2(2)
      integer*8 nhand(mdoo+1,nsymb,2),i1,in,i1c,i1o,i0c,i0o,itesta,     11d26s20
     $     itestb,gandcc,gandco,gandcb                                  10d13s22
      logical ldebug                                                    1d25s21
      dimension nff1(mdoo+1,nsymb,2),iff1(*),nff0(mdoo+1,3),            11d25s20
     $     iff0(*),vec(ncsfv,*),ncsf(*),multh(8,8),nab4(2,3),idorb1(32),11d25s20
     $     isorb1(32),idorb(32),isorb(32),itest(32,4),ih0av(*),ionex(*),11d26s20
     $     nvirt(*),iden1x(8,8),nden1x(8,8),nh0av(*),ionext(*),ioxx(2)  11d1s22
      include "common.store"                                            11d25s20
      include "common.mrci"                                             11d25s20
      common/kmfind/invk1(2,8,8,8,2)
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      if(ldebug)write(6,*)('hi, I am the new and improved hcsi!')       1d25s21
      ibcoffo=ibcoff                                                    12d8s20
      ighere=ibcoff                                                     3d22s21
      iacc=ighere+maxbx                                                 3d22s21
      ibcoff=iacc+mynprocg                                              3d22s21
      call enough('hcsi.  1',bc,ibc)
      nacc=0                                                            3d22s21
      if(ldebug)write(6,*)('izero in hcsi '),izero                      1d25s21
      if(izero.ne.0)then                                                1d24s21
       do isb=1,nsymb                                                    11d26s20
        do ii=mdon+1,mdoo+1                                              11d26s20
         if(nhand(ii,isb,2).ge.0)then                                    11d26s20
          call ddi_zero(bc,ibc,nhand(ii,isb,2))                         11d15s22
         end if                                                          11d26s20
        end do                                                           11d26s20
       end do                                                            11d26s20
      end if                                                            1d24s21
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
         igflag=0                                                       3d22s21
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
         call enough('hcsi.  2',bc,ibc)
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
             gandcc=ieor(i1c,i0c)                                       10d13s22
             gandco=ieor(i1o,i0o)                                       10d13s22
             gandcb=ior(gandcc,gandco)                                  11d1s22
             ndifb=popcnt(gandcb)                                       10d13s22
             if(ndifb.le.4)then                                         11d1s22
              ndifd=popcnt(gandcc)                                       10d13s22
              ndifs=popcnt(gandco)                                       10d13s22
              if(ndifs.eq.2.and.ndifb.eq.2)then
               do i=1,norbx
                if(btest(gandco,i))then
                 if((btest(i0c,i).and.btest(i1o,i)).or.
     $             (btest(i0o,i).and..not.btest(i1c,i)))then
                  nab4(2,1)=i
                 else
                  nab4(1,1)=i
                 end if
                end if
               end do
               call gandc(i1c,i1o,i0c,i0o,nopen1,nopen,iarg,jarg,ncsf,    11d25s20
     $            norbx,ixw1,ixw2,nnot1,nab1,iwpb,iwpk,ncsfmid,bc,ibc)  11d14s22
               iprod=ibcoff                                                  11d13s20
               ibcoff=iprod+ncsf(iarg)*nroot                              11d26s20
               call enough('hcsi.  3',bc,ibc)
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
              else                                                      11d1s22
               nnot=0
               if(ndifs.eq.4.and.ndifb.eq.4)then                                           10d14s22
                nnot=4                                                          10d14s22
                ioxx(1)=1                                                     10d17s22
                ioxx(2)=1                                                     10d17s22
                do i=1,norbx                                                     10d17s22
                 if(btest(gandcb,i))then                                         10d14s22
                  if((btest(i0c,i).and.btest(i1o,i)).or.                         10d17s22
     $              (btest(i0o,i).and..not.btest(i1c,i)))then                   10d14s22
                   nab4(2,ioxx(2))=i                                       10d17s22
                   ioxx(2)=ioxx(2)+1                                       10d17s22
                  else                                                           10d14s22
                   nab4(1,ioxx(1))=i                                       10d17s22
                   ioxx(1)=ioxx(1)+1                                       10d17s22
                  end if                                                         10d14s22
                 end if
                end do
               else if(ndifb.eq.3)then                                           10d14s22
                nnot=3
                ioxx(1)=1                                                        10d14s22
                ioxx(2)=1                                                        10d14s22
                iswap=0                                                          10d17s22
                do i=1,norbx
                 if(btest(gandcb,i))then                                         10d14s22
                  if(btest(gandcc,i).and.                                        10d14s22
     $               ((btest(i1c,i).and..not.btest(i0o,i)).or.                 10d14s22
     $            (btest(i0c,i).and..not.btest(i1o,i))))then            12d19s22
                   if(btest(i0c,i))iswap=1                                        10d17s22
                   nab4(1,1)=i                                                10d17s22
                   nab4(1,2)=i                                                10d17s22
                  else                                                           10d14s22
                   nab4(2,ioxx(2))=i                                             10d14s22
                   ioxx(2)=ioxx(2)+1
                  end if
                 end if                                                          10d14s22
                end do
                if(iswap.ne.0)then                                               10d17s22
                 icpy=nab4(1,1)                                               10d17s22
                 nab4(1,1)=nab4(2,1)                                       10d17s22
                 nab4(2,1)=icpy                                               10d17s22
                 icpy=nab4(1,2)                                               10d17s22
                 nab4(1,2)=nab4(2,2)                                       10d17s22
                 nab4(2,2)=icpy                                               10d17s22
                 nbt=0                                                           10d17s22
                 if(btest(i1c,nab4(1,2)).and..not.btest(i1c,nab4(1,1)))    10d17s22
     $         nbt=1                                                      10d17s22
                else                                                             10d17s22
                 nbt=0                                                           10d17s22
                 if(btest(i0c,nab4(2,2)).and..not.btest(i0c,nab4(2,1)))
     $               nbt=1                                                      10d17s22
                end if                                                           10d17s22
                if(nbt.ne.0)then                                                 10d17s22
                 nab4(1,1)=nab4(1,2)                                             10d17s22
                 nab4(2,1)=nab4(2,2)                                             10d17s22
                end if                                                           10d17s22
               else if(ndifs.eq.0.and.ndifd.eq.2)then                            10d14s22
                nnot=3
                do i=1,norbx
                 if(btest(gandcb,i))then                                         10d14s22
                  if(btest(i1c,i))then
                   nab4(1,1)=i
                   nab4(1,2)=i
                  else                                                           10d14s22
                   nab4(2,1)=i
                   nab4(2,2)=i
                  end if
                 end if                                                          10d14s22
                end do
               end if                                                            10d14s22
               ipssx=0
               if(nnot.eq.3)then                                         1d25s21
                ipssx=1                                                   1d25s21
               else if(nnot.eq.4)then                                   11d1s22
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
    1          continue                                                  1d25s21
              end if                                                      11d25s20
             end if                                                     10d13s22
             jvs=jvs+ncsf(jarg)                                          11d26s20
            end do                                                       11d25s20
           end do                                                        11d25s20
           itmpg=ibcoff                                                 12d7s20
           ibcoff=itmpg+nnuser*nvirt(isbv)                              11d10s20
           itmpgt=ibcoff                                                10d4s22
           ibcoff=itmpgt+nnuser*nvirt(isbv)                             10d4s22
           call enough('hcsi.  4',bc,ibc)
           factg=0d0                                                    12d7s20
           if(irefo(isbv).gt.0)then
            nok=0                                                        12d2s20
            jdenh=idenh                                                  12d2s20
            ikeep=ibcoff                                                 12d2s20
            ibcoff=ikeep+irefo(isbv)                                     12d2s20
            call enough('hcsi.  5',bc,ibc)
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
             itmp=ibcoff                                                 12d2s20
             ibcoff=itmp+nvirt(isbv)*nok                                 12d2s20
             itmpt=ibcoff                                               10d4s22
             ibcoff=itmpt+nvirt(isbv)*nok                               10d4s22
             ident=ibcoff                                               10d4s22
             ibcoff=ident+nnuser*nok                                    10d4s22
             call enough('hcsi.  6',bc,ibc)
             do i=0,nok-1                                               10d4s22
              do j=0,nnuser-1                                           10d4s22
               ji=idenh+j+nnuser*i                                      10d4s22
               ij=ident+i+nok*j                                         10d4s22
               bc(ij)=bc(ji)                                            10d4s22
              end do                                                    10d4s22
             end do                                                     10d4s22
             jtmpt=itmpt                                                10d4s22
             do i=0,nok-1                                               10d4s22
              iad2=ih0av(isbv)+irefo(isbv)+nh0av(isbv)*ibc(ikeep+i)     10d4s22
              do iv=0,nvirt(isbv)-1                                     10d4s22
               bc(jtmpt+iv)=bc(iad2+iv)                                 10d4s22
              end do                                                    10d4s22
              jtmpt=jtmpt+nvirt(isbv)                                   10d4s22
             end do                                                     10d4s22
             call dgemm('n','n',nvirt(isbv),nnuser,nok,1d0,bc(itmpt),   10d4s22
     $            nvirt(isbv),bc(ident),nok,0d0,bc(itmpgt),nvirt(isbv), 10d4s22
     d' hcsi.  1')
             factg=1d0                                                  12d7s20
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
               itmpi=ibcoff                                             12d7s20
               ibcoff=itmpi+nok*nvirt(isbv)                             12d7s20
               itmpit=ibcoff                                            10d4s22
               ident=itmpit+nok*nvirt(isbv)                             10d4s22
               ibcoff=ident+nok*nnuser                                  10d4s22
               call enough('hcsi.  7',bc,ibc)
               jden=iden1x(isc,isa)                                     12d7s20
               jtmpi=itmpi                                              12d7s20
               jtmpit=itmpit                                            10d4s22
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
                jint=ionext(i2eu)+icoln*nvirt(isbv)                     3d23s21
                do iv=0,nvirt(isbv)-1                                   12d7s20
                 bc(jtmpit+iv)=bc(jint+iv)                              10d4s22
                end do                                                  12d7s20
                jtmpit=jtmpit+nvirt(isbv)                               10d4s22
               end do                                                   12d7s20
               do i=0,nnuser-1                                          10d4s22
                do j=0,nok-1                                            10d4s22
                 ji=ident+j+nok*i                                       10d4s22
                 ij=iden1x(isc,isa)+i+nnuser*j                          10d4s22
                 bc(ji)=bc(ij)                                          10d4s22
                end do                                                  10d4s22
               end do                                                   10d4s22
               call dgemm('n','n',nvirt(isbv),nnuser,nok,1d0,           10d4s22
     $              bc(itmpit),nvirt(isbv),bc(ident),nok,factg,         10d4s22
     $              bc(itmpgt),nvirt(isbv),                             10d4s22
     d' hcsi.  2')
               factg=1d0                                                12d7s20
              end if                                                    12d7s20
             end if
            end do
           end do                                                        11d26s20
           if(abs(factg-1d0).lt.1d-3)then                               12d7s20
            if(igflag.eq.0)then                                         3d22s21
             call ddi_done(ibc(iacc),nacc)                                     3d22s21
             do i=0,nwds-1                                                  3d22s21
              bc(ighere+i)=0d0                                          3d22s21
             end do                                                         3d22s21
             igflag=1                                                   3d22s21
            end if                                                      3d22s21
            do iv=0,nvirt(isbv)-1                                       12d7s20
             do i=nnlow,nnhgh                                           12d7s20
              ip=ifcol+i-il                                             12d7s20
              im=i-nnlow                                                12d7s20
              do ir=0,nroot-1
               iad1=itmpgt+iv+nvirt(isbv)*(im+nnuse*ir)                 10d4s22
               iad2=ighere+ir+nroot*(iv+nvirt(isbv)*ip)                 12d14s20
               bc(iad2)=bc(iad2)+bc(iad1)                                12d7s20
              end do                                                    11d10s20
             end do                                                     12d7s20
            end do                                                      12d7s20
           end if                                                       12d7s20
          else                                                          12d8s20
           jff1=jff1+2                                                  12d8s20
          end if                                                        12d7s20
          ifcol=itcol+1                                                 12d7s20
         end do                                                         11d25s20
         if(ifirsttime.eq.0)call dws_synca                              11d26s20
         ifirsttime=1                                                   11d26s20
         i1c=nvirt(isbv)*nroot                                          11d10s20
         i1o=il                                                         12d7s20
         in=ih                                                          12d7s20
         if(igflag.ne.0)then                                            3d22s21
          nwiacc=nwiacc+nroot*nvirt(isbv)*(ih+1-il)                     8d10s22
          call ddi_iacc(bc,ibc,nhand(nclo1p,isb,2),i1,i1c,i1o,in,       11d15s22
     $         bc(ighere),ibc(iacc),nacc)                               11d15s22
         end if                                                         3d22s21
         ibcoff=idenh                                                   3d22s21
        end if                                                          11d25s20
       end do                                                            11d25s20
      end do                                                            11d26s20
      ibcoff=ibcoffo                                                    3d22s21
      call ddi_done(ibc(iacc),nacc)                                     8d31s22
      if(ldebug)then                                                    1d25s21
       write(6,*)('all done in hcsi')
       call dws_synca                                                    1d24s21
       write(6,*)('what we have for hcsi:')
       do isb=1,nsymb                                                    12d7s20
        write(6,*)('for isb = '),isb
        isbv=multh(isb,isymmrci)                                         12d7s20
        i1c=nvirt(isbv)*nroot                                            11d10s20
        nhere=0                                                          12d7s20
        do ii=1,mdoo+1                                                   12d7s20
         iarg=ii-mdon                                                    12d7s20
         ndblock=ncsf(iarg)*nff1(ii,isb,1)                               11d10s20
         nhere=nhere+ndblock
        end do
        mdblock=nhere                                                    11d10s20
        nsing=mdblock*nvirt(isbv)                                        12d7s20
        igs=ibcoff
        ibcoff=igs+nhere*nvirt(isbv)*nroot                               11d10s20
        call enough('hcsi.  8',bc,ibc)
        jgs=igs
        do ii=1,mdoo+1                                                   12d7s20
         iarg=ii-mdon                                                    12d7s20
         ndblock=ncsf(iarg)*nff1(ii,isb,1)                               11d10s20
         if(ndblock.gt.0)then
          itmp=ibcoff                                                    12d7s20
          ibcoff=itmp+ndblock*nvirt(isbv)*nroot                          11d10s20
          call enough('hcsi.  9',bc,ibc)
          in=ndblock                                                     12d7s20
          call ddi_get(bc,ibc,nhand(ii,isb,2),i1,i1c,i1,in,bc(itmp))    11d15s22
          jtmp=itmp                                                      12d7s20
          do iff=1,nff1(ii,isb,1)                                        12d7s20
           do m=0,ncsf(iarg)-1
            do iv=0,nvirt(isbv)-1
             do ir=0,nroot-1
              kgs=jgs+nsing*ir+iv                                        12d14s20
              bc(kgs)=bc(jtmp+ir)                                        12d14s20
             end do
             jtmp=jtmp+nroot                                             12d14s20
            end do
            jgs=jgs+nvirt(isbv)                                          11d10s20
           end do
          end do                                                         12d7s20
          ibcoff=itmp                                                    12d7s20
         end if
        end do                                                           12d7s20
        write(6,*)('all together now')
        call prntm2(bc(igs),nsing,nroot,nsing)                           11d10s20
        ibcoff=igs
       end do                                                            12d7s20
      end if                                                            1d25s21
      return
      end
