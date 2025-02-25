c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcisbk(ihsdiag,nff1,iff1,nff0,iff0,gi,ncsfv,ncsf,nec,  7d9s21
     $     mdon,mdoop,nsymb,multh,ixw1,ixw2,lxmt,ixmt,nvirt,nrootu,           7d9s21
     $     lprt,isymop,n2e,i2eop,phase,phase2,ism,irel,irefo,norb,      7d9s21
     $     nbasdws,idoubo,ixmtf,isymmrci,bc,ibc)                        11d10s22
      implicit real*8 (a-h,o-z)
      integer*1 nab1(2),nab2(2)
      logical lprt                                                      1d25s21
      integer*8 i1,in,i1c,i1o,i0c,i0o,itesta,itestb,i28,i38,i48,        7d9s21
     $     ihsdiag(mdoop,nsymb,2),gandcc,gandco,gandcb                  2d6s23
      dimension nff1(mdoop,nsymb,2),iff1(*),nff0(mdoop,3),gi(ncsfv,*),  7d9s21
     $     iff0(*),ncsf(*),multh(8,8),nab4(2,3),idorb1(32),             12d12s20
     $     isorb1(32),idorb(32),isorb(32),itest(32,4),ioxx(2),          2d6s23
     $     nvirt(*),iden1x(8,8),nden1x(8,8),ism(*),irel(*),irefo(*),    7d9s21
     $     nbasdws(*),idoubo(*),ixmtf(*),ixmt(8,*),isymop(*),i2eop(2,3) 7d9s21
      include "common.store"                                            11d25s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data icall/0/                                                     12d5s22
      save icall                                                        12d5s22
      icall=icall+1                                                     12d5s22
      if(lprt)write(6,*)('hi, I am the new and improved hcisbk!')
      mdoo=mdoop-1                                                      7d9s21
      nrootm=nrootu-1                                                   7d12s21
      i1=1                                                              11d10s20
      ig=ibcoff                                                         12d12s20
      ibcoff=ig+ncsfv*nrootu                                            12d12s20
      call enough('hcisbk.  1',bc,ibc)
      do i=ig,ibcoff-1                                                  12d12s20
       bc(i)=0d0                                                        12d12s20
      end do                                                            12d12s20
      ibcoffo=ibcoff                                                    12d8s20
      norbx=norb+1                                                      11d25s20
      i1=1                                                              11d25s20
      ivs=1
      do isb=1,nsymb                                                    11d25s20
       isbv=multh(isb,isymmrci)
       isbo=multh(isbv,isymop(lxmt))                                    7d12s21
       do nclo1p=mdon+1,mdoop                                           7d9s21
        if(min(nvirt(isbv),nff1(nclo1p,isb,1)).gt.0)then                  11d25s20
         nclo1=nclo1p-1                                                  11d25s20
         nopen1=nec-nclo1*2                                             11d25s20
         iarg=nclo1p-mdon                                                11d25s20
         ncolt=nff1(nclo1p,isb,1)*ncsf(iarg)                            11d10s20
         call ilimts(1,ncolt,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)  12d7s20
         ncolh=ih+1-il                                                  12d7s20
         ivst=ibcoff                                                    12d14s20
         ibcoff=ivst+nvirt(isbv)*nrootu*ncolh                           12d14s20
         ivtmp=ibcoff                                                   7d9s21
         ibcoff=ivtmp+nvirt(isbv)*nrootu*ncolh                          7d9s21
         call enough('hcisbk.  2',bc,ibc)
         i28=nvirt(isbv)*nrootu                                         7d9s21
         i38=il                                                         7d9s21
         i48=ih                                                         7d9s21
         call ddi_get(bc,ibc,ihsdiag(nclo1p,isb,1),i1,i28,i38,i48,      11d15s22
     $        bc(ivtmp))                                                11d15s22
         do i=0,ncolh-1                                                 12d14s20
          do iv=0,nvirt(isbv)-1                                         12d14s20
           do ir=0,nrootu-1                                             12d14s20
            iad1=ivtmp+ir+nrootu*(iv+nvirt(isbv)*i)                     7d9s21
            iad2=ivst+iv+nvirt(isbv)*(ir+nrootu*i)                      12d14s20
            bc(iad2)=bc(iad1)                                           7d9s21
           end do                                                       12d14s20
          end do                                                        12d14s20
         end do                                                         12d14s20
         ibcoff=ivtmp                                                   7d9s21
         nn=ncsf(iarg)*nrootu                                           12d11s20
         nnm=nn-1                                                       11d26s20
         jff1=nff1(nclo1p,isb,2)                                        11d25s20
         ifcol=1                                                        12d7s20
         jvst=ivst                                                      7d9s21
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
           do nclop=max(mdon+1,nclo1p-2),min(mdoop,nclo1p+2)            7d9s21
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
             gandcc=ieor(i1c,i0c)                                        2d6s23
             gandco=ieor(i1o,i0o)                                        2d6s23
             gandcb=ior(gandcc,gandco)                                     2d6s23
             ndifb=popcnt(gandcb)                                          2d6s23
             if(ndifb.le.4)then                                         2d6s23
              ndifs=popcnt(gandco)                                        2d6s23
              ndifd=popcnt(gandcc)                                        2d6s23
              if(ndifs.eq.2.and.ndifb.eq.2)then                         2d6s23
               do i=1,norbx                                             2d6s23
                if(btest(gandco,i))then                                   2d6s23
                 if((btest(i1o,i).and..not.btest(i0c,i)).or.            2d6s23
     $             (btest(i1c,i).and.btest(i0o,i)))then                 2d6s23
                  nab4(1,1)=i                                             2d6s23
                 else                                                     2d6s23
                  nab4(2,1)=i                                             2d6s23
                 end if                                                   2d6s23
                end if                                                    2d6s23
               end do                                                     2d6s23
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
               jsb=ism(nab4(2,1))                                        7d12s21
               jgb=irel(nab4(2,1))-1+idoubo(jsb)                         7d19s21
               if(jsb.ne.isbo)then                                       7d12s21
                write(6,*)('jsb vs. isbo: '),jsb,isbo
                stop 'hcisbk'
               end if                                                    7d12s21
               iint=ibcoff                                               3d22s21
               ibcoff=iint+nvirt(isbv)                                   3d22s21
               call enough('hcisbk.  3',bc,ibc)
               ihq=ixmtf(jsb)+jgb                                       2d6s23
     $              +nbasdws(jsb)*(irefo(isbv)+idoubo(isbv))            2d6s23
               do iv=0,nvirt(isbv)-1                                     3d22s21
                bc(iint+iv)=bc(ihq+iv*nbasdws(jsb))*phase                7d12s21
               end do                                                    3d22s21
               do i2e=1,n2e                                              7d9s21
                do i=1,nok                                                    11d13s20
                 js=ism(itest(i,4))                                           11d13s20
                 jg=irel(itest(i,4))-1                                        11d13s20
                 if(isymop(i2eop(1,i2e)).eq.1)then                        7d12s21
c     (bv|gg) or (gg|bv)
                  fjj=bc(ixmt(js,i2eop(1,i2e))+idoubo(js)+jg              7d12s21
     $                +nbasdws(js)*(idoubo(js)+jg))                      7d12s21
                  jjj=ixmt(isbo,i2eop(2,i2e))+jgb                        7d19s21
     $               +nbasdws(isbo)*(idoubo(isbv)+irefo(isbv))          7d12s21
                  if(itest(i,3).eq.2)fjj=fjj*2d0                         7d12s21
                  fjj=fjj*phase2                                         7d12s21
                  do iv=0,nvirt(isbv)-1                                  7d12s21
                   bc(iint+iv)=bc(iint+iv)+fjj*bc(jjj+iv*nbasdws(isbo))  7d12s21
                  end do                                                 7d12s21
                 end if                                                  7d12s21
                 if(itest(i,3).eq.2)then                                      11d13s20
                  itry1=multh(js,isbo)
                  itry2=multh(js,isbv)
                  if(itry1.eq.isymop(i2eop(1,i2e)).and.                  7d12s21
     $               itry2.eq.isymop(i2eop(2,i2e)))then                 7d12s21
c     (bg|gv)
                   fkk=bc(ixmt(isbo,i2eop(1,i2e))+jgb                    7d19s21
     $               +nbasdws(isbo)*(idoubo(js)+jg))                    7d12s21
                   kkk=ixmt(js,i2eop(2,i2e))+idoubo(js)+jg                7d12s21
     $                +nbasdws(js)*(idoubo(isbv)+irefo(isbv))           7d12s21
                   fkk=fkk*phase2                                        7d12s21
                   do iv=0,nvirt(isbv)-1                                  7d9s21
                    bc(iint+iv)=bc(iint+iv)-fkk*bc(kkk+iv*nbasdws(js))   7d12s21
                   end do                                                7d12s21
                  end if                                                 7d12s21
                 end if                                                  7d12s21
                end do                                                   7d12s21
               end do                                                    7d9s21
               iprod=ibcoff                                              3d22s21
               ibcoff=iprod+ncsf(iarg)*ncsf(jarg)                        3d22s21
               call enough('hcisbk.  4',bc,ibc)
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
               call enough('hcisbk.  5',bc,ibc)
               do i=itmpp,ibcoff-1                                       3d22s21
                bc(i)=0d0                                                3d22s21
               end do                                                    3d22s21
               if(nvirt(isbv).gt.nnuse)then                              3d22s21
                nrow=nvirt(isbv)*nrootu                                  3d22s21
                ivdprod=ibcoff                                           3d22s21
                ibcoff=ivdprod+nrow*ncsf(jarg)                           3d22s21
                call enough('hcisbk.  6',bc,ibc)
                call dgemm('n','n',nrow,ncsf(jarg),nnuse,1d0,            3d22s21
     $              bc(jvst),nrow,bc(jprod),ncsf(iarg),0d0,             7d9s21
     $              bc(ivdprod),nrow,                                   3d22s21
     d' hcisbk.  1')
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
                call enough('hcisbk.  7',bc,ibc)
                do i=iivprod,ibcoff-1                                    3d22s21
                 bc(i)=0d0                                               3d22s21
                end do                                                   3d22s21
                iix=jvst                                                 7d9s21
                do i=0,nnuser-1                                          3d22s21
                 do iv=0,nvirt(isbv)-1                                   3d22s21
                  bc(iivprod+i)=bc(iivprod+i)+bc(iint+iv)*bc(iix+iv)     7d9s21
                 end do                                                  3d22s21
                 iix=iix+nvirt(isbv)                                     3d22s21
                end do                                                   3d22s21
                call dgemm('n','n',nrootu,ncsf(jarg),nnuse,1d0,          3d22s21
     $              bc(iivprod),nrootu,bc(jprod),ncsf(iarg),1d0,        3d22s21
     $              bc(itmpp),nrootu,                                   3d22s21
     d' hcisbk.  2')
                ibcoff=iivprod
               end if                                                    3d22s21
               if(n2e.gt.0)then                                         2d6s23
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
                   if(nnot1.eq.2.and.nnot2.eq.2.and.n2e.gt.0)then         7d9s21
                    call genmat(ncsf(iarg),ncsf(karg),ncsf(jarg),         3d22s21
     $                  ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod, 11d10s22
     $                bc,ibc)                                           11d10s22
c     (ba|dv)
                    jsa=ism(nab1(1))                                              11d13s20
                    jga=irel(nab1(1))-1                                           11d13s20
                    jsb=ism(nab1(2))                                              11d13s20
                    jgb=irel(nab1(2))-1                                           11d13s20
                    jsd=ism(nab2(2))                                              11d13s20
                    jgd=irel(nab2(2))-1                                           11d13s20
                    kint=ibcoff                                           7d9s21
                    ibcoff=kint+nvirt(isbv)                               7d9s21
                    call enough('hcisbk.  8',bc,ibc)
                    do iz=kint,ibcoff-1                                   7d19s21
                     bc(iz)=0d0                                           7d19s21
                    end do                                                7d19s21
                    do i2e=1,n2e                                          7d19s21
                     if(multh(jsa,jsb).eq.isymop(i2eop(1,i2e)))then          7d14s21
                      fba=phase2                                        2d6s23
     $                    *bc(ixmt(jsb,i2eop(1,i2e))+jgb+idoubo(jsb)    2d6s23
     $                 +nbasdws(jsb)*(jga+idoubo(jsa)))                 7d14s21
                      idv=ixmt(jsd,i2eop(2,i2e))+jgd+idoubo(jsd)           7d12s21
     $                 +nbasdws(jsd)*(idoubo(isbv)+irefo(isbv))         7d12s21
                      do iv=0,nvirt(isbv)-1                                 7d9s21
                       bc(kint+iv)=bc(kint+iv)                          2d6s23
     $                      +fba*bc(idv+iv*nbasdws(jsd))                2d6s23
                      end do                                                7d9s21
                     end if                                                7d19s21
                    end do                                                7d9s21
                    jprod=iprod+nnlow                                     3d22s21
                    if(nvirt(isbv).gt.nnuse)then                              3d22s21
                     nrow=nvirt(isbv)*nrootu                                  3d22s21
                     ivdprod=ibcoff                                           3d22s21
                     ibcoff=ivdprod+nrow*ncsf(jarg)                           3d22s21
                     call enough('hcisbk.  9',bc,ibc)
                     call dgemm('n','n',nrow,ncsf(jarg),nnuse,1d0,            3d22s21
     $                 bc(jvst),nrow,bc(jprod),ncsf(iarg),0d0,             7d9s21
     $                 bc(ivdprod),nrow,                                   3d22s21
     d' hcisbk.  3')
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
                     call enough('hcisbk. 10',bc,ibc)
                     do iz=iivprod,ibcoff-1                               3d22s21
                      bc(iz)=0d0                                          3d22s21
                     end do                                                   3d22s21
                     iix=jvst                                             7d9s21
                     do iz=0,nnuser-1                                          3d22s21
                      do iv=0,nvirt(isbv)-1                                   3d22s21
                       bc(iivprod+iz)=bc(iivprod+iz)+bc(kint+iv)          3d23s21
     $                    *bc(iix+iv)                                   7d9s21
                      end do                                                  3d22s21
                      iix=iix+nvirt(isbv)                                     3d22s21
                     end do                                                   3d22s21
                     call dgemm('n','n',nrootu,ncsf(jarg),nnuse,1d0,          3d22s21
     $                  bc(iivprod),nrootu,bc(jprod),ncsf(iarg),1d0,        3d22s21
     $                  bc(itmpp),nrootu,                                   3d22s21
     d' hcisbk.  4')
                     ibcoff=iivprod
                    end if                                                    3d22s21
                    ibcoff=iprod                                          3d22s21
                   end if
                  end if                                                    11d13s20
                 end if                                                       11d13s20
                end do                                                        11d13s20
               end if                                                   2d6s23
               do ir=0,nrootm                                            3d22s21
                iadg=ig+jvs-1+ncsfv*ir                                   3d22s21
                do j=0,ncsf(jarg)-1                                      3d22s21
                 jad=itmpp+ir+nrootu*j                                   3d22s21
                 bc(iadg+j)=bc(iadg+j)+bc(jad)                           3d22s21
                end do                                                   3d22s21
               end do                                                    3d22s21
               ibcoff=iint                                               3d22s21
              else if(n2e.gt.0)then                                     2d6s23
               nnot=0                                                   2d6s23
               if(ndifs.eq.4.and.ndifb.eq.4)then                        2d6s23
                nnot=4                                                  2d6s23
                ioxx(1)=1                                               2d6s23
                ioxx(2)=1                                               2d6s23
                do i=1,norbx                                            2d6s23
                 if(btest(gandcb,i))then                                2d6s23
                  if((btest(i0c,i).and.btest(i1o,i)).or.                2d6s23
     $                   (btest(i0o,i).and..not.btest(i1c,i)))then      2d6s23
                   nab4(2,ioxx(2))=i                                    2d6s23
                   ioxx(2)=ioxx(2)+1                                    2d6s23
                  else                                                  2d6s23
                   nab4(1,ioxx(1))=i                                    2d6s23
                   ioxx(1)=ioxx(1)+1                                    2d6s23
                  end if                                                2d6s23
                 end if                                                 2d6s23
                end do                                                  2d6s23
               else if(ndifb.eq.3)then                                  2d6s23
                nnot=3                                                  2d6s23
                ioxx(1)=1                                               2d6s23
                ioxx(2)=1                                               2d6s23
                iswap=0                                                 2d6s23
                do i=1,norbx                                            2d6s23
                 if(btest(gandcb,i))then                                2d6s23
                  if(btest(gandcc,i).and.                               2d6s23
     $                   ((btest(i1c,i).and..not.btest(i0o,i)).or.      2d6s23
     $                   (btest(i0c,i).and..not.btest(i1o,i))))then     2d6s23
                   if(btest(i0c,i))iswap=1                              2d6s23
                   nab4(1,1)=i                                          2d6s23
                   nab4(1,2)=i                                          2d6s23
                  else                                                  2d6s23
                   nab4(2,ioxx(2))=i                                    2d6s23
                   ioxx(2)=ioxx(2)+1                                    2d6s23
                  end if                                                2d6s23
                 end if                                                 2d6s23
                end do                                                  2d6s23
                if(iswap.ne.0)then                                      2d6s23
                 icpy=nab4(1,1)                                         2d6s23
                 nab4(1,1)=nab4(2,1)                                    2d6s23
                 nab4(2,1)=icpy                                         2d6s23
                 icpy=nab4(1,2)                                         2d6s23
                 nab4(1,2)=nab4(2,2)                                    2d6s23
                 nab4(2,2)=icpy                                         2d6s23
                 nbt=0                                                  2d6s23
                 if(btest(i1c,nab4(1,2)).and.                           2d6s23
     $                   .not.btest(i1c,nab4(1,1)))nbt=1                2d6s23
                else                                                    2d6s23
                 nbt=0                                                  2d6s23
                 if(btest(i0c,nab4(2,2)).and.                           2d6s23
     $                   .not.btest(i0c,nab4(2,1)))nbt=1                2d6s23
                end if                                                  2d6s23
                if(nbt.ne.0)then                                        2d6s23
                 nab4(1,1)=nab4(1,2)                                    2d6s23
                 nab4(2,1)=nab4(2,2)                                    2d6s23
                end if                                                  2d6s23
               else if(ndifs.eq.0.and.ndifd.eq.2)then                   2d6s23
                nnot=3                                                  2d6s23
                do i=1,norbx                                            2d6s23
                 if(btest(gandcb,i))then                                2d6s23
                  if(btest(i1c,i))then                                  2d6s23
                   nab4(1,1)=i                                          2d6s23
                   nab4(1,2)=i                                          2d6s23
                  else                                                  2d6s23
                   nab4(2,1)=i                                          2d6s23
                   nab4(2,2)=i                                          2d6s23
                  end if                                                2d6s23
                 end if                                                 2d6s23
                end do                                                  2d6s23
               end if                                                   2d6s23
               ipssx=0                                                  2d6s23
               if(nnot.eq.3)then                                         3d22s21
                ipssx=1                                                  3d22s21
               else if(nnot.eq.4)then                                   2d6s23
                ipssx=2                                                  3d22s21
               end if                                                    3d22s21
               iint=ibcoff                                               3d22s21
               itmpp=iint+nvirt(isbv)                                    3d22s21
               ibcoff=itmpp+nrootu*ncsf(jarg)                            3d22s21
               call enough('hcisbk. 11',bc,ibc)
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
c     (ba|dv)
                jsa=ism(nab1(1))                                              11d13s20
                jga=irel(nab1(1))-1                                           11d13s20
                jsb=ism(nab1(2))                                              11d13s20
                jgb=irel(nab1(2))-1                                           11d13s20
                jsd=ism(nab2(2))                                           11d27s20
                jgd=irel(nab2(2))-1                                        11d27s20
                call genmat(ncsf(iarg),ncsf(karg),ncsf(jarg),             3d22s21
     $                 ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod, 11d10s22
     $              bc,ibc)                                             11d10s22
                kint=ibcoff                                              7d9s21
                ibcoff=kint+nvirt(isbv)                                  7d9s21
                call enough('hcisbk. 12',bc,ibc)
                do iz=kint,ibcoff-1                                       7d19s21
                 bc(iz)=0d0                                               7d19s21
                end do                                                   7d19s21
                do i2e=1,n2e                                             7d19s21
                 if(multh(jsa,jsb).eq.isymop(i2eop(1,i2e)))then          7d19s21
                  iba=ixmt(jsb,i2eop(1,i2e))+jgb+idoubo(jsb)             7d14s21
     $              +nbasdws(jsb)*(jga+idoubo(jsa))                     7d14s21
                  fba=phase2*bc(iba)                                     7d14s21
                  idv=ixmt(jsd,i2eop(2,i2e))+jgd+idoubo(jsd)              7d9s21
     $              +nbasdws(jsd)*(irefo(isbv)+idoubo(isbv))            7d9s21
                  do iv=0,nvirt(isbv)-1                                    7d9s21
                   bc(kint+iv)=bc(kint+iv)+fba*bc(idv+iv*nbasdws(jsd))   7d14s21
                  end do                                                   7d9s21
                 end if                                                  7d19s21
                end do                                                   7d9s21
                jprod=iprod+nnlow                                         3d22s21
                if(nvirt(isbv).gt.nnuse)then                              3d22s21
                 nrow=nvirt(isbv)*nrootu                                  3d22s21
                 ivdprod=ibcoff                                           3d22s21
                 ibcoff=ivdprod+nrow*ncsf(jarg)                           3d22s21
                 call enough('hcisbk. 13',bc,ibc)
                 call dgemm('n','n',nrow,ncsf(jarg),nnuse,1d0,            3d22s21
     $              bc(jvst),nrow,bc(jprod),ncsf(iarg),0d0,             7d9s21
     $              bc(ivdprod),nrow,                                   3d22s21
     d' hcisbk.  5')
                 jvdprod=ivdprod                                          3d22s21
                 do iz=0,nrootu*ncsf(jarg)-1                               3d22s21
                  do iv=0,nvirt(isbv)-1                                   3d22s21
                   bc(itmpp+iz)=bc(itmpp+iz)+bc(jvdprod+iv)               3d22s21
     $                  *bc(kint+iv)                                      3d23s21
                  end do                                                  3d22s21
                  jvdprod=jvdprod+nvirt(isbv)                             3d22s21
                 end do                                                   3d22s21
                 ibcoff=ivdprod                                           3d22s21
                else                                                      3d22s21
                 iivprod=ibcoff                                           3d22s21
                 ibcoff=iivprod+nnuser                                    3d22s21
                 call enough('hcisbk. 14',bc,ibc)
                 do iz=iivprod,ibcoff-1                                   3d22s21
                  bc(iz)=0d0                                              3d22s21
                 end do                                                   3d22s21
                 iix=jvst                                                7d9s21
                 do iz=0,nnuser-1                                          3d22s21
                  do iv=0,nvirt(isbv)-1                                   3d22s21
                   bc(iivprod+iz)=bc(iivprod+iz)+bc(kint+iv)             3d23s21
     $                    *bc(iix+iv)                                   7d9s21
                  end do                                                  3d22s21
                  iix=iix+nvirt(isbv)                                     3d22s21
                 end do                                                   3d22s21
                 call dgemm('n','n',nrootu,ncsf(jarg),nnuse,1d0,          3d22s21
     $              bc(iivprod),nrootu,bc(jprod),ncsf(iarg),1d0,        3d22s21
     $              bc(itmpp),nrootu,                                   3d22s21
     d' hcisbk.  6')
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
             end if                                                     2d6s23
             jvs=jvs+ncsf(jarg)                                         3d22s21
            end do                                                       11d25s20
           end do                                                        11d25s20
           jvst=jvst+nvirt(isbv)*nnuser                                 7d9s21
          else                                                          12d8s20
           jff1=jff1+2                                                  12d8s20
          end if                                                        12d7s20
          ifcol=itcol+1                                                 12d7s20
         end do                                                         11d25s20
         ibcoff=ivst                                                    3d22s21
        end if                                                          11d25s20
       end do                                                            11d25s20
      end do                                                            11d26s20
      if(lprt)write(6,*)('all done in hcisbk')                            1d25s21
      jg=ig-1                                                           7d9s21
      do ir=1,nrootu                                                    7d9s21
       do i=1,ncsfv                                                     7d9s21
        gi(i,ir)=gi(i,ir)+bc(jg+i)                                      7d9s21
       end do                                                           7d9s21
       jg=jg+ncsfv                                                      7d9s21
      end do                                                            7d9s21
      ibcoff=ig                                                         12d12s20
      return
      end
