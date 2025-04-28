c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcis(vs,nff1,iff1,nff0,iff0,igi,ncsfv,ncsf,nec,mdon,   12d12s20
     $     mdoo,nsymb,multh,ixw1,ixw2,ih0av,nh0av,ionex,nvirt,nrootu,   1d25s21
     $     lprt,ionext,bc,ibc)                                          11d10s22
      implicit real*8 (a-h,o-z)
      integer*1 nab1(2),nab2(2)
      logical lprt                                                      1d25s21
      integer*8 i1,in,i1c,i1o,i0c,i0o,itesta,itestb,igi,gandcc,gandco,  11d1s22
     $     gandcb                                                       11d1s22
      dimension nff1(mdoo+1,nsymb,2),iff1(*),nff0(mdoo+1,3),            11d25s20
     $     iff0(*),ncsf(*),multh(8,8),nab4(2,3),idorb1(32),             12d12s20
     $     isorb1(32),idorb(32),isorb(32),itest(32,4),ih0av(*),ionex(*),11d26s20
     $     nvirt(*),iden1x(8,8),nden1x(8,8),nh0av(*),vs(*),ionext(*),   11d1s22
     $     ioxx(2)                                                      11d1s22
      include "common.store"                                            11d25s20
      include "common.mrci"                                             11d25s20
      common/kmfind/invk1(2,8,8,8,2)
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      if(lprt)write(6,*)('hi, I am the new and improved hcis!')
      nrootm=nroot-1                                                    3d22s21
      i1=1                                                              11d10s20
      ig=ibcoff                                                         12d12s20
      ibcoff=ig+ncsfv*nrootu                                            12d12s20
      call enough('hcis.  1',bc,ibc)
      do i=ig,ibcoff-1                                                  12d12s20
       bc(i)=0d0                                                        12d12s20
      end do                                                            12d12s20
      ibcoffo=ibcoff                                                    12d8s20
      norbx=norb+1                                                      11d25s20
      i1=1                                                              11d25s20
      ivs=1
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
         ivst=ibcoff                                                    12d14s20
         ibcoff=ivst+nvirt(isbv)*nrootu*ncolh                           12d14s20
         call enough('hcis.  2',bc,ibc)
         do i=0,ncolh-1                                                 12d14s20
          do iv=0,nvirt(isbv)-1                                         12d14s20
           do ir=0,nrootu-1                                             12d14s20
            iad1=ivs+ir+nrootu*(iv+nvirt(isbv)*i)                       12d14s20
            iad2=ivst+iv+nvirt(isbv)*(ir+nrootu*i)                      12d14s20
            bc(iad2)=vs(iad1)                                           12d14s20
           end do                                                       12d14s20
          end do                                                        12d14s20
         end do                                                         12d14s20
         do i=0,ncolh*nrootu*nvirt(isbv)-1                              12d14s20
          vs(ivs+i)=bc(ivst+i)                                          12d14s20
         end do                                                         12d14s20
         ibcoff=ivst                                                    12d14s20
         nn=ncsf(iarg)*nrootu                                           12d11s20
         nnm=nn-1                                                       11d26s20
         jff1=nff1(nclo1p,isb,2)                                        11d25s20
         ifcol=1                                                        12d7s20
         do if1=1,nff1(nclo1p,isb,1)                                    11d26s20
          itcol=ifcol+ncsf(iarg)-1                                      11d10s20
          if(ifcol.le.ih.and.itcol.ge.il)then                           12d7s20
           nnlow=max(0,il-ifcol)                                        12d7s20
           nnhgh=min(ncsf(iarg)-1,ih-ifcol)                             12d12s20
           nnuse=nnhgh+1-nnlow                                          12d7s20
           nnuser=nnuse*nrootu                                          12d11s20
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
           do nclop=max(mdon+1,nclo1p-2),min(mdoo+1,nclo1p+2)             11d25s20
            nclo=nclop-1                                                  11d25s20
            nopen=nec-nclo*2                                             11d25s20
            jarg=nclop-mdon                                               11d25s20
            jff0=nff0(nclop,2)                                           11d25s20
            jvs=nff0(nclop,3)                                            11d26s20
            nun=nrootu*ncsf(jarg)                                       10d4s22
            do jf=1,nff0(nclop,1)                                        11d26s20
             i0c=iff0(jff0)                                              11d25s20
             jff0=jff0+1                                                 11d25s20
             i0o=iff0(jff0)                                              11d25s20
             jff0=jff0+1                                                 11d25s20
             gandcc=ieor(i1c,i0c)                                       10d13s22
             gandco=ieor(i1o,i0o)                                       10d13s22
             gandcb=ior(gandcc,gandco)                                  10d20s22
             ndifb=popcnt(gandcb)                                       10d20s22
             if(ndifb.le.4)then                                         10d20s22
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
              iint=ibcoff                                               3d22s21
              ibcoff=iint+nvirt(isbv)                                   3d22s21
              call enough('hcis.  3',bc,ibc)
              ihq=ih0av(isbv)+irefo(isbv)+nh0av(isbv)*jgb               3d22s21
              do iv=0,nvirt(isbv)-1                                     3d22s21
               bc(iint+iv)=bc(ihq+iv)                                   3d22s21
              end do                                                    3d22s21
              do i=1,nok                                                    11d13s20
               js=ism(itest(i,4))                                           11d13s20
               jg=irel(itest(i,4))-1                                        11d13s20
               icolj=((jg*(jg+1)/2))+jg                                  11d26s20
               nnvj=(irefo(js)*(irefo(js)+1))/2                          11d26s20
               iicol=icolj+nnvj*jgb                                      11d28s20
               i2eu=invk1(1,js,js,isbv,2)                               3d22s21
               iad=ionext(i2eu)+iicol*nvirt(isbv)                       3d23s21
               if(itest(i,3).eq.2)then                                      11d13s20
                if(isbv.eq.js)then                                       11d26s20
                 ix=max(jg,jgb)                                          11d26s20
                 in=min(jg,jgb)                                          11d26s20
                 icolk=((ix*(ix+1))/2)+in                                11d26s20
                 iicol=icolk+nnvj*jg                                     12d1s20
                 kad=ionext(i2eu)+iicol*nvirt(isbv)                     3d23s21
                 do iv=0,nvirt(isbv)-1                                    3d22s21
                  bc(iint+iv)=bc(iint+iv)+bc(iad+iv)*2d0-bc(kad+iv)     3d23s21
                 end do                                                   3d22s21
                else                                                     11d26s20
                 nnvk=irefo(js)*irefo(isbv)                              11d26s20
                 i2eu=invk1(1,isbv,js,js,2)                             3d22s21
                 icase=invk1(2,isbv,js,js,2)                            3d22s21
                 if(icase.eq.2)then                                     3d22s21
                  icolk=jg+irefo(js)*jgb                                  11d26s20
                  iicol=icolk+nnvk*jg                                    11d28s20
                 else                                                   3d22s21
                  icolk=jgb+irefo(isbv)*jg                               11d26s20
                  iicol=icolk+nnvk*jg                                    11d28s20
                 end if                                                 3d22s21
                 kad=ionext(i2eu)+iicol*nvirt(isbv)                     3d23s21
                 do iv=0,nvirt(isbv)-1                                    3d22s21
                  bc(iint+iv)=bc(iint+iv)+bc(iad+iv)*2d0-bc(kad+iv)     3d23s21
                 end do                                                 3d22s21
                end if                                                   11d26s20
               else
                do iv=0,nvirt(isbv)-1                                    3d22s21
                 bc(iint+iv)=bc(iint+iv)+bc(iad+iv)                     3d23s21
                end do                                                   3d22s21
               end if                                                   12d11s20
              end do                                                        11d13s20
              iprod=ibcoff                                              3d22s21
              ibcoff=iprod+ncsf(iarg)*ncsf(jarg)                        3d22s21
              call enough('hcis.  4',bc,ibc)
              call prodn(iwpb,iwpk,ncsf(iarg),ncsf(jarg),ncsfmid,       3d22s21
     $             bc(iprod),bc,ibc,1d0,0d0)                            2d13s23
c     Vjr=sum vi iintv*Vsvri*Dij
c     op counts
c     sum v iintv*(sum i Vsvri*Dij)=Nv*nrootu*Nj+Nv*nrootu*Nj*Ni
c                                  =Nv*nrootu*Nj*(1+Ni)
c     sum i Dij*(sumv iintv*Vsvri) =Ni*Nj*nrootu+Nv*nrootu*Nj*Ni)
c                                  =Ni*Nj*nrootu*(1+Nv)
c
              jprod=iprod+nnlow                                         3d22s21
              itmpp=ibcoff                                              3d22s21
              ibcoff=itmpp+nrootu*ncsf(jarg)                            3d22s21
              call enough('hcis.  5',bc,ibc)
              do i=itmpp,ibcoff-1                                       3d22s21
               bc(i)=0d0                                                3d22s21
              end do                                                    3d22s21
              if(nvirt(isbv).gt.nnuse)then                              3d22s21
               nrow=nvirt(isbv)*nrootu                                  3d22s21
               ivdprod=ibcoff                                           3d22s21
               ibcoff=ivdprod+nrow*ncsf(jarg)                           3d22s21
               call enough('hcis.  6',bc,ibc)
               call dgemm('n','n',nrow,ncsf(jarg),nnuse,1d0,            3d22s21
     $              vs(ivs),nrow,bc(jprod),ncsf(iarg),0d0,              3d22s21
     $              bc(ivdprod),nrow,                                   3d22s21
     d' hcis.  1')
               jvdprod=ivdprod                                          3d22s21
               do i=0,nrootu*ncsf(jarg)-1                               3d22s21
                do iv=0,nvirt(isbv)-1                                   3d22s21
                 bc(itmpp+i)=bc(itmpp+i)+bc(jvdprod+iv)*bc(iint+iv)     3d22s21
                end do                                                  3d22s21
                jvdprod=jvdprod+nvirt(isbv)                             3d22s21
               end do                                                   3d22s21
               ibcoff=ivdprod                                           3d22s21
              else                                                      3d22s21
               iivprod=ibcoff                                           3d22s21
               ibcoff=iivprod+nnuser                                    3d22s21
               call enough('hcis.  7',bc,ibc)
               do i=iivprod,ibcoff-1                                    3d22s21
                bc(i)=0d0                                               3d22s21
               end do                                                   3d22s21
               iix=ivs                                                  3d22s21
               do i=0,nnuser-1                                          3d22s21
                do iv=0,nvirt(isbv)-1                                   3d22s21
                 bc(iivprod+i)=bc(iivprod+i)+bc(iint+iv)*vs(iix+iv)     3d22s21
                end do                                                  3d22s21
                iix=iix+nvirt(isbv)                                     3d22s21
               end do                                                   3d22s21
               call dgemm('n','n',nrootu,ncsf(jarg),nnuse,1d0,          3d22s21
     $              bc(iivprod),nrootu,bc(jprod),ncsf(iarg),1d0,        3d22s21
     $              bc(itmpp),nrootu,                                   3d22s21
     d' hcis.  2')
               ibcoff=iivprod
              end if                                                    3d22s21
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
                  call genmat(ncsf(iarg),ncsf(karg),ncsf(jarg),         3d22s21
     $                 ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod, 11d10s22
     $                bc,ibc)                                           11d10s22
                  jsa=ism(nab1(1))                                              11d13s20
                  jga=irel(nab1(1))-1                                           11d13s20
                  jsb=ism(nab1(2))                                              11d13s20
                  jgb=irel(nab1(2))-1                                           11d13s20
                  jsd=ism(nab2(2))                                              11d13s20
                  jgd=irel(nab2(2))-1                                           11d13s20
                  i2eu=invk1(1,jsa,jsb,jsd,2)                           3d22s21
                  icase=invk1(2,jsa,jsb,jsd,2)                          3d22s21
                  if(jsa.eq.jsb)then                                     11d26s20
                   ix=max(jgb,jga)                                      3d22s21
                   in=min(jgb,jga)                                      3d22s21
                   icolk=((ix*(ix+1))/2)+in                             3d22s21
                   nnvk=(irefo(jsa)*(irefo(jsa)+1))/2                    11d27s20
                  else                                                   11d26s20
                   nnvk=irefo(jsa)*irefo(jsb)                            11d26s20
                   if(icase.eq.1)then                                   3d22s21
                    icolk=jga+irefo(jsa)*jgb                            3d22s21
                   else                                                 3d22s21
                    icolk=jgb+irefo(jsb)*jga                            3d22s21
                   end if                                               3d22s21
                  end if                                                 11d26s20
                  kint=ionext(i2eu)+nvirt(isbv)*(icolk+nnvk*jgd)        3d23s21
                  jprod=iprod+nnlow                                     3d22s21
                  if(nvirt(isbv).gt.nnuse)then                              3d22s21
                   nrow=nvirt(isbv)*nrootu                                  3d22s21
                   ivdprod=ibcoff                                           3d22s21
                   ibcoff=ivdprod+nrow*ncsf(jarg)                           3d22s21
                   call enough('hcis.  8',bc,ibc)
                   call dgemm('n','n',nrow,ncsf(jarg),nnuse,1d0,            3d22s21
     $              vs(ivs),nrow,bc(jprod),ncsf(iarg),0d0,              3d22s21
     $              bc(ivdprod),nrow,                                   3d22s21
     d' hcis.  3')
                   jvdprod=ivdprod                                          3d22s21
                   do iz=0,nrootu*ncsf(jarg)-1                               3d22s21
                    do iv=0,nvirt(isbv)-1                                   3d22s21
                     bc(itmpp+iz)=bc(itmpp+iz)+bc(jvdprod+iv)           3d22s21
     $                    *bc(kint+iv)                                  3d23s21
                    end do                                                  3d22s21
                    jvdprod=jvdprod+nvirt(isbv)                             3d22s21
                   end do                                                   3d22s21
                   ibcoff=ivdprod                                           3d22s21
                  else                                                      3d22s21
                   iivprod=ibcoff                                           3d22s21
                   ibcoff=iivprod+nnuser                                    3d22s21
                   call enough('hcis.  9',bc,ibc)
                   do iz=iivprod,ibcoff-1                               3d22s21
                    bc(iz)=0d0                                          3d22s21
                   end do                                                   3d22s21
                   iix=ivs                                                  3d22s21
                   do iz=0,nnuser-1                                          3d22s21
                    do iv=0,nvirt(isbv)-1                                   3d22s21
                     bc(iivprod+iz)=bc(iivprod+iz)+bc(kint+iv)          3d23s21
     $                    *vs(iix+iv)                                   3d22s21
                    end do                                                  3d22s21
                    iix=iix+nvirt(isbv)                                     3d22s21
                   end do                                                   3d22s21
                   call dgemm('n','n',nrootu,ncsf(jarg),nnuse,1d0,          3d22s21
     $              bc(iivprod),nrootu,bc(jprod),ncsf(iarg),1d0,        3d22s21
     $              bc(itmpp),nrootu,                                   3d22s21
     d' hcis.  4')
                   ibcoff=iivprod
                  end if                                                    3d22s21
                  ibcoff=iprod                                          3d22s21
                 end if
                end if                                                    11d13s20
               end if                                                       11d13s20
              end do                                                        11d13s20
              do ir=0,nrootm                                            3d22s21
               iadg=ig+jvs-1+ncsfv*ir                                   3d22s21
               do j=0,ncsf(jarg)-1                                      3d22s21
                jad=itmpp+ir+nrootu*j                                   3d22s21
                bc(iadg+j)=bc(iadg+j)+bc(jad)                           3d22s21
               end do                                                   3d22s21
              end do                                                    3d22s21
              ibcoff=iint                                               3d22s21
             else                                                       10d20s22
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
              ipssx=0                                                   10d20s22
              if(nnot.eq.3)then                                         3d22s21
               ipssx=1                                                  3d22s21
              else if(nnot.eq.4)then                                    11d1s22
               ipssx=2                                                  3d22s21
              end if                                                    3d22s21
              iint=ibcoff                                               3d22s21
              itmpp=iint+nvirt(isbv)                                    3d22s21
              ibcoff=itmpp+nrootu*ncsf(jarg)                            3d22s21
              call enough('hcis. 10',bc,ibc)
              do i=itmpp,ibcoff-1                                       3d22s21
               bc(i)=0d0                                                3d22s21
              end do                                                    3d22s21
              do ipss=1,ipssx                                           3d22s21
               itesta=i1c                                                 11d26s20
               itestb=i1o                                                 11d26s20
               if(btest(itesta,nab4(1,1)))then                            11d27s20
                itesta=ibclr(itesta,nab4(1,1))                            11d27s20
                itestb=ibset(itestb,nab4(1,1))                            11d27s20
                nopenk=nopen1+1                                              11d13s20
                karg=iarg-1                                                  11d13s20
               else if(btest(itestb,nab4(1,1)))then                       11d27s20
                itestb=ibclr(itestb,nab4(1,1))                            11d27s20
                nopenk=nopen1-1                                              11d13s20
                karg=iarg                                                  11d13s20
               end if                                                        11d13s20
               if(btest(itestb,nab4(2,ipss)))then                       3d22s21
                itesta=ibset(itesta,nab4(2,ipss))                        3d22s21
                itestb=ibclr(itestb,nab4(2,ipss))                        3d22s21
                nopenk=nopenk-1                                              11d13s20
                karg=karg+1                                                  11d13s20
               else                                                          11d13s20
                itestb=ibset(itestb,nab4(2,ipss))                       3d22s21
                nopenk=nopenk+1                                              11d13s20
               end if                                                        11d13s20
               call gandc(i1c,i1o,itesta,itestb,nopen1,nopenk,            11d26s20
     $         iarg,karg,ncsf,norbx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,   11d27s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
               call gandc(itesta,itestb,i0c,i0o,nopenk,nopen,             11d26s20
     $         karg,jarg,ncsf,norbx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,   11d27s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
               jsa=ism(nab1(1))                                              11d13s20
               jga=irel(nab1(1))-1                                           11d13s20
               jsb=ism(nab1(2))                                              11d13s20
               jgb=irel(nab1(2))-1                                           11d13s20
               jsd=ism(nab2(2))                                           11d27s20
               jgd=irel(nab2(2))-1                                        11d27s20
               call genmat(ncsf(iarg),ncsf(karg),ncsf(jarg),             3d22s21
     $                 ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod, 11d10s22
     $              bc,ibc)                                             11d10s22
               i2eu=invk1(1,jsa,jsb,jsd,2)                               3d22s21
               icase=invk1(2,jsa,jsb,jsd,2)                              3d22s21
               if(jsa.eq.jsb)then                                         11d27s20
                nnvk=(irefo(jsa)*(irefo(jsa)+1))/2                        11d27s20
                ix=max(jgb,jga)                                         3d22s21
                in=min(jgb,jga)                                         3d22s21
                icolkn=((ix*(ix+1))/2)+in                               3d22s21
               else                                                       11d26s20
                nnvk=irefo(jsa)*irefo(jsb)                                11d27s20
                if(icase.eq.1)then                                       3d22s21
                 icolkn=jga+irefo(jsa)*jgb                               3d22s21
                else                                                     3d22s21
                 icolkn=jgb+irefo(jsb)*jga                               3d22s21
                end if                                                   3d22s21
               end if                                                     11d26s20
               kint=ionext(i2eu)+nvirt(isbv)*(icolkn+nnvk*jgd)          3d23s21
               jprod=iprod+nnlow                                         3d22s21
               if(nvirt(isbv).gt.nnuse)then                              3d22s21
                nrow=nvirt(isbv)*nrootu                                  3d22s21
                ivdprod=ibcoff                                           3d22s21
                ibcoff=ivdprod+nrow*ncsf(jarg)                           3d22s21
                call enough('hcis. 11',bc,ibc)
                call dgemm('n','n',nrow,ncsf(jarg),nnuse,1d0,            3d22s21
     $              vs(ivs),nrow,bc(jprod),ncsf(iarg),0d0,              3d22s21
     $              bc(ivdprod),nrow,                                   3d22s21
     d' hcis.  5')
                jvdprod=ivdprod                                          3d22s21
                do iz=0,nrootu*ncsf(jarg)-1                               3d22s21
                 do iv=0,nvirt(isbv)-1                                   3d22s21
                  bc(itmpp+iz)=bc(itmpp+iz)+bc(jvdprod+iv)               3d22s21
     $                *bc(kint+iv)                                      3d23s21
                 end do                                                  3d22s21
                 jvdprod=jvdprod+nvirt(isbv)                             3d22s21
                end do                                                   3d22s21
                ibcoff=ivdprod                                           3d22s21
               else                                                      3d22s21
                iivprod=ibcoff                                           3d22s21
                ibcoff=iivprod+nnuser                                    3d22s21
                call enough('hcis. 12',bc,ibc)
                do iz=iivprod,ibcoff-1                                   3d22s21
                 bc(iz)=0d0                                              3d22s21
                end do                                                   3d22s21
                iix=ivs                                                  3d22s21
                do iz=0,nnuser-1                                          3d22s21
                 do iv=0,nvirt(isbv)-1                                   3d22s21
                  bc(iivprod+iz)=bc(iivprod+iz)+bc(kint+iv)             3d23s21
     $                    *vs(iix+iv)                                   3d22s21
                 end do                                                  3d22s21
                 iix=iix+nvirt(isbv)                                     3d22s21
                end do                                                   3d22s21
                call dgemm('n','n',nrootu,ncsf(jarg),nnuse,1d0,          3d22s21
     $              bc(iivprod),nrootu,bc(jprod),ncsf(iarg),1d0,        3d22s21
     $              bc(itmpp),nrootu,                                   3d22s21
     d' hcis.  6')
                ibcoff=iivprod
               end if                                                    3d22s21
               ibcoff=iprod                                              3d22s21
              end do
              do ir=0,nrootm                                            3d22s21
               iadg=ig+jvs-1+ncsfv*ir                                   3d22s21
               do j=0,ncsf(jarg)-1                                      3d22s21
                jad=itmpp+ir+nrootu*j                                   3d22s21
                bc(iadg+j)=bc(iadg+j)+bc(jad)                           3d22s21
               end do                                                   3d22s21
              end do                                                    3d22s21
              ibcoff=iint                                               3d22s21
             end if                                                      11d25s20
             end if                                                     10d13s22
             jvs=jvs+ncsf(jarg)                                         3d22s21
            end do                                                       11d25s20
           end do                                                        11d25s20
           ivs=ivs+nvirt(isbv)*nnuser                                   12d11s20
          else                                                          12d8s20
           jff1=jff1+2                                                  12d8s20
          end if                                                        12d7s20
          ifcol=itcol+1                                                 12d7s20
         end do                                                         11d25s20
         ibcoff=ivst                                                    3d22s21
        end if                                                          11d25s20
       end do                                                            11d25s20
      end do                                                            11d26s20
      if(lprt)write(6,*)('all done in hcis')                            1d25s21
      in=ncsfv*nrootu                                                   12d12s20
      call ddi_acc(bc,ibc,igi,i1,i1,i1,in,bc(ig))                       11d15s22
      if(lprt)then                                                      1d25s21
       call dws_gsumf(bc(ig),ncsfv*nrootu)
       write(6,*)('global summed')
       call prntm2(bc(ig),ncsfv,nrootu,ncsfv)
      end if                                                            1d25s21
      ibcoff=ig                                                         12d12s20
      return
      end
