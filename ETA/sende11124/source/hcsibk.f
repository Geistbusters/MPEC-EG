c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcsibk(nhand,nff1,iff1,nff0,iff0,vec,ncsfv,ncsf,nec,     11d25s20
     $     mdon,mdoop,nsymb,multh,ixw1,ixw2,lxmt,ixmt,nvirt,            7d9s21
     $     izero,ldebug,maxbx,isymop,n2e,i2eop,phase,phase2,ism,irel,   7d9s21
     $     irefo,norb,nbasdws,idoubo,ixmtf,nroot,isymmrci,nwiacc,bc,ibc)1d26s23
      implicit real*8 (a-h,o-z)
      integer*1 nab1(2),nab2(2)
      integer*8 nhand(mdoop,nsymb),i1,in,i1c,i1o,i0c,i0o,itesta,        8d21s21
     $     itestb,gandcc,gandco,gandcb                                  2d6s23
      logical ldebug                                                    1d25s21
      dimension nff1(mdoop,nsymb,2),iff1(*),nff0(mdoop,3),              7d9s21
     $     iff0(*),vec(ncsfv,*),ncsf(*),multh(8,8),nab4(2,3),idorb1(32),11d25s20
     $     isorb1(32),idorb(32),isorb(32),itest(32,4),ixmt(8,*),        7d9s21
     $     nvirt(*),iden1x(8,8,3),nden1x(8,8,3),i2eop(2,3),isymop(*),   7d19s21
     $     ism(*),irel(*),irefo(*),nbasdws(*),idoubo(*),ixmtf(*),ioxx(2)2d6s23
      include "common.store"                                            11d25s20
      common/paddcm/npadddi                                             6d7s22
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      if(ldebug)write(6,*)('hi, I am the new and improved hcsibk!')       1d25s21
      mdoo=mdoop-1                                                      7d9s21
c
c     maxbx is for input wavefcns ...
c     but here bra has ket roots, so maxbx may not be large enough.     6d7s22
c
      maxbxb=0                                                          6d7s22
      do isb=1,nsymb                                                    6d7s22
       isbv=multh(isb,isymmrci)                                         6d7s22
       do ii=mdon+1,mdoop                                               6d30s23
        iarg=ii-mdon                                                    6d7s22
        maxbxb=max(maxbxb,ncsf(iarg)*nff1(ii,isb,1)*nvirt(isbv)*nroot)  6d7s22
       end do                                                           6d7s22
      end do                                                            6d7s22
      maxbxb=maxbxb+npadddi                                             6d7s22
      ibcoffo=ibcoff                                                    12d8s20
      ighere=ibcoff                                                     3d22s21
      iacc=ighere+maxbxb                                                6d7s22
      ibcoff=iacc+mynprocg                                              3d22s21
      call enough('hcsibk.  1',bc,ibc)
      nfirst=0                                                          7d14s21
      nlast=0                                                           7d14s21
      nacc=0                                                            3d22s21
      if(ldebug)write(6,*)('izero in hcsi '),izero                      1d25s21
      if(izero.ne.0)then                                                1d24s21
       do isb=1,nsymb                                                    11d26s20
        do ii=mdon+1,mdoop                                              7d9s21
         if(nhand(ii,isb).ge.0)then                                     8d21s21
          call ddi_zero(bc,ibc,nhand(ii,isb))                           11d15s22
         end if                                                          11d26s20
        end do                                                           11d26s20
       end do                                                            11d26s20
      end if                                                            1d24s21
      norbx=norb+1                                                      11d25s20
      i1=1                                                              11d25s20
      ifirsttime=0                                                      11d26s20
      do isb=1,nsymb                                                    11d25s20
       isbv=multh(isb,isymmrci)
       do nclo1p=mdon+1,mdoop                                           7d9s21
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
         isbo=multh(isbv,isymop(lxmt))                                  7d12s21
         idenh=ibcoff                                                   11d25s20
         ndenh=idenh+nn*irefo(isbo)                                     7d12s21
         ibcoff=ndenh+irefo(isbo)                                       7d12s21
         do isa=1,nsymb                                                 11d26s20
          isav=multh(isa,isbo)                                          7d12s21
          do isc=1,nsymb                                                7d14s21
           nnv=irefo(isa)*irefo(isc)                                    7d12s21
           iscav=multh(isc,isav)                                        11d26s20
           do ii2e=1,n2e                                                7d19s21
            iden1x(isc,isa,ii2e)=ibcoff                                 7d19s21
            nden1x(isc,isa,ii2e)=iden1x(isc,isa,ii2e)+nn*               7d19s21
     $           nnv*irefo(iscav)                                       7d19s21
            ibcoff=nden1x(isc,isa,ii2e)+nnv*irefo(iscav)                7d19s21
           end do                                                       7d19s21
          end do                                                        11d26s20
         end do                                                         11d26s20
         call enough('hcsibk.  2',bc,ibc)
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
               iprod=ibcoff                                                  11d13s20
               ibcoff=iprod+ncsf(iarg)*nroot                              11d26s20
               call enough('hcsibk.  3',bc,ibc)
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
               jsb=ism(nab4(2,1))                                        7d12s21
               if(jsb.ne.isbo)then                                       7d12s21
                write(6,*)('jsb vs. isbo: '),jsb,isbo
                stop 'hcsibk'
               end if
               jdenh=idenh+nn*jgb                                         11d26s20
               if(isymop(lxmt).eq.multh(jsb,isbv))then                   7d19s21
                ibc(ndenh+jgb)=1                                         7d19s21
                iad2=ixmtf(isbv)+idoubo(isbv)+irefo(isbv)
     $             +nbasdws(isbv)*(jgb+idoubo(jsb))                     7d19s21
                do i=0,nnm                                                 11d26s20
                 bc(jdenh+i)=bc(jdenh+i)+bc(iprod+i)                       11d26s20
                end do                                                     11d26s20
               end if                                                    7d19s21
               do i=1,nok                                                    11d13s20
                js=ism(itest(i,4))                                           11d13s20
                jg=irel(itest(i,4))-1                                        11d13s20
                iicolj=jgb+irefo(isbo)*(jg+irefo(js)*jg)                  7d12s21
c     (vb|gg)
                ff=1d0                                                   7d12s21
                if(itest(i,3).eq.2)ff=2d0                                7d12s21
                do ii2e=1,n2e                                            7d19s21
                 if(isymop(i2eop(1,ii2e)).eq.1)then                      7d19s21
                  jden=iden1x(isbo,js,ii2e)+nn*iicolj                    7d19s21
                  ibc(nden1x(isbo,js,ii2e)+iicolj)=1                     7d19s21
                  do ii=0,nnm                                               11d26s20
                   bc(jden+ii)=bc(jden+ii)+bc(iprod+ii)*ff                7d12s21
                  end do                                                   11d26s20
                 end if                                                   7d12s21
                end do                                                   7d19s21
c     is=isbv,jga runs over virts jg,jg,jgb,iv
                if(itest(i,3).eq.2)then                                      11d13s20
                 itry1=multh(isbv,js)                                    7d12s21
                 iicol=jg+irefo(js)*(jg+irefo(js)*jgb)                   7d12s21
                 do ii2e=1,n2e                                           7d19s21
                  if(itry1.eq.isymop(i2eop(1,ii2e)))then                 7d19s21
                   kden=iden1x(js,js,ii2e)+nn*iicol                      7d19s21
                   ibc(nden1x(js,js,ii2e)+iicol)=1                       7d19s21
                   do ii=0,nnm                                               11d26s20
                    bc(kden+ii)=bc(kden+ii)-bc(iprod+ii)                       11d26s20
                   end do                                                   11d26s20
                  end if                                                       11d13s20
                 end do                                                  7d19s21
                end if                                                    7d12s21
               end do                                                        11d13s20
               ibcoff=iprod
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
                   if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
c     (ab|vd) = (vd|ab)
                    jsa=ism(nab1(1))                                              11d13s20
                    jga=irel(nab1(1))-1                                           11d13s20
                    jsb=ism(nab1(2))                                              11d13s20
                    jgb=irel(nab1(2))-1                                           11d13s20
                    jsd=ism(nab2(2))                                              11d13s20
                    jgd=irel(nab2(2))-1                                           11d13s20
                    iicol=jgd+irefo(jsd)*(jga+irefo(jsa)*jgb)             7d14s21
                    do ii2e=1,n2e                                         7d19s21
                     if(multh(jsa,jsb).eq.isymop(i2eop(2,ii2e)))then      7d19s21
                      jdenk=iden1x(jsd,jsa,ii2e)+nn*iicol                 7d19s21
                      ibc(nden1x(jsd,jsa,ii2e)+iicol)=1                   7d19s21
                      call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),        12d7s20
     $                 ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,       12d7s20
     $                 vec(jvs,1),ncsfv,nroot,bc(jdenk),1d0,bc,ibc)     11d10s22
                     end if                                                7d14s21
                    end do                                                7d19s21
                   end if
                  end if                                                    11d13s20
                 end if                                                       11d13s20
                end do                                                        11d13s20
               end if                                                   2d6s23
              else if(n2e.gt.0)then                                      2d6s23
               nnot=0                                                    2d6s23
               if(ndifs.eq.4.and.ndifb.eq.4)then                         2d6s23
                nnot=4                                                   2d6s23
                ioxx(1)=1                                                2d6s23
                ioxx(2)=1                                                2d6s23
                do i=1,norbx                                             2d6s23
                 if(btest(gandcb,i))then                                 2d6s23
                  if((btest(i0c,i).and.btest(i1o,i)).or.                 2d6s23
     $               (btest(i0o,i).and..not.btest(i1c,i)))then          2d6s23
                   nab4(2,ioxx(2))=i                                     2d6s23
                   ioxx(2)=ioxx(2)+1                                     2d6s23
                  else                                                   2d6s23
                   nab4(1,ioxx(1))=i                                     2d6s23
                   ioxx(1)=ioxx(1)+1                                     2d6s23
                  end if                                                 2d6s23
                 end if                                                  2d6s23
                end do                                                   2d6s23
               else if(ndifb.eq.3)then                                   2d6s23
                nnot=3                                                   2d6s23
                ioxx(1)=1                                                2d6s23
                ioxx(2)=1                                                2d6s23
                iswap=0                                                  2d6s23
                do i=1,norbx                                             2d6s23
                 if(btest(gandcb,i))then                                 2d6s23
                  if(btest(gandcc,i).and.                                2d6s23
     $               ((btest(i1c,i).and..not.btest(i0o,i)).or.          2d6s23
     $               (btest(i0c,i).and..not.btest(i1o,i))))then         2d6s23
                   if(btest(i0c,i))iswap=1                               2d6s23
                   nab4(1,1)=i                                           2d6s23
                   nab4(1,2)=i                                           2d6s23
                  else                                                   2d6s23
                   nab4(2,ioxx(2))=i                                     2d6s23
                   ioxx(2)=ioxx(2)+1                                     2d6s23
                  end if                                                 2d6s23
                 end if                                                  2d6s23
                end do                                                   2d6s23
                if(iswap.ne.0)then                                       2d6s23
                 icpy=nab4(1,1)                                          2d6s23
                 nab4(1,1)=nab4(2,1)                                     2d6s23
                 nab4(2,1)=icpy                                          2d6s23
                 icpy=nab4(1,2)                                          2d6s23
                 nab4(1,2)=nab4(2,2)                                     2d6s23
                 nab4(2,2)=icpy                                          2d6s23
                 nbt=0                                                   2d6s23
                 if(btest(i1c,nab4(1,2)).and.                            2d6s23
     $                .not.btest(i1c,nab4(1,1)))nbt=1                    2d6s23
                else                                                     2d6s23
                 nbt=0                                                   2d6s23
                 if(btest(i0c,nab4(2,2)).and.                            2d6s23
     $                   .not.btest(i0c,nab4(2,1)))nbt=1                2d6s23
                end if                                                   2d6s23
                if(nbt.ne.0)then                                         2d6s23
                 nab4(1,1)=nab4(1,2)                                     2d6s23
                 nab4(2,1)=nab4(2,2)                                     2d6s23
                end if                                                   2d6s23
               else if(ndifs.eq.0.and.ndifd.eq.2)then                    2d6s23
                nnot=3                                                   2d6s23
                do i=1,norbx                                             2d6s23
                 if(btest(gandcb,i))then                                 2d6s23
                  if(btest(i1c,i))then                                   2d6s23
                   nab4(1,1)=i                                           2d6s23
                   nab4(1,2)=i                                           2d6s23
                  else                                                   2d6s23
                   nab4(2,1)=i                                           2d6s23
                   nab4(2,2)=i                                           2d6s23
                  end if                                                 2d6s23
                 end if                                                  2d6s23
                end do                                                   2d6s23
               end if                                                    2d6s23
               ipssx=0                                                   2d6s23
               if(nnot.eq.3)then                                         1d25s21
                ipssx=1                                                   1d25s21
               else if(nnot.eq.4)then                                    2d6s23
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
c     (ab|vd) = (vd|ab)
                 jsa=ism(nab1(1))                                              11d13s20
                 jga=irel(nab1(1))-1                                           11d13s20
                 jsb=ism(nab1(2))                                              11d13s20
                 jgb=irel(nab1(2))-1                                           11d13s20
                 jsd=ism(nab2(2))                                           11d27s20
                 jgd=irel(nab2(2))-1                                        11d27s20
                 iicol=jgd+irefo(jsd)*(jga+irefo(jsa)*jgb)               7d14s21
                 do ii2e=1,n2e                                           7d19s21
                  if(multh(jsa,jsb).eq.isymop(i2eop(2,ii2e)))then            7d14s21
                   kden=iden1x(jsd,jsa,ii2e)+nn*iicol                    7d19s21
                   ibc(nden1x(jsd,jsa,ii2e)+iicol)=1                     7d19s21
                   call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),        7d19s21
     $                 ncsfmid1,                                        7d19s21
     $            iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,vec(jvs,1),ncsfv,    12d7s20
     $             nroot,bc(kden),1d0,bc,ibc)                           11d10s22
                  end if                                                  7d14s21
                 end do                                                  7d19s21
                 if(ipss.eq.2)go to 1                                    1d25s21
                end if                                                   1d25s21
               end do                                                    1d25s21
    1          continue                                                  1d25s21
              end if                                                    2d6s23
             end if                                                      11d25s20
             jvs=jvs+ncsf(jarg)                                          11d26s20
            end do                                                       11d25s20
           end do                                                        11d25s20
           itmpg=ibcoff                                                 12d7s20
           ibcoff=itmpg+nnuser*nvirt(isbv)                              11d10s20
           call enough('hcsibk.  4',bc,ibc)
           factg=0d0                                                    12d7s20
           if(irefo(isbo).gt.0)then                                     7d12s21
            nok=0                                                        12d2s20
            jdenh=idenh                                                  12d2s20
            ikeep=ibcoff                                                 12d2s20
            ibcoff=ikeep+irefo(isbo)                                    7d12s21
            call enough('hcsibk.  5',bc,ibc)
            do i=0,irefo(isbo)-1                                        7d12s21
             if(ibc(ndenh+i).gt.0)then                                  3d19s21
              iad=idenh+nn*i
              ibc(ikeep+nok)=i+idoubo(isbo)                             7d19s21
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
             call enough('hcsibk.  6',bc,ibc)
             do iv=0,nvirt(isbv)-1                                       12d2s20
              iad1=itmp+nok*iv                                           12d2s20
              iad2=ixmtf(isbv)+idoubo(isbv)+irefo(isbv)+iv              7d12s21
              do i=0,nok-1                                               12d2s20
               ip=ibc(ikeep+i)                                          7d19s21
               bc(iad1+i)=bc(iad2+ip*nbasdws(isbv))                     7d12s21
              end do                                                     12d2s20
             end do                                                      12d2s20
             jdenh=idenh+nnlow                                          12d12s20
             call dgemm('n','n',nnuser,nvirt(isbv),nok,phase,bc(idenh), 7d9s21
     $            nnuser,bc(itmp),nok,factg,bc(itmpg),nnuser,           12d12s20
     d' hcsibk.  1')
             factg=1d0                                                  12d7s20
            end if                                                       12d2s20
            ibcoff=ikeep                                                 12d2s20
           end if
           do isa=1,nsymb                                                11d26s20
            isav=multh(isa,isbo)                                        7d12s21
            do isc=1,nsymb                                              7d14s21
             isavc=multh(isav,isc)                                       11d26s20
             nnv=irefo(isc)*irefo(isa)
             ncol=nnv*irefo(isavc)
             if(ncol.gt.0)then
              do ii2e=1,n2e                                             7d19s21
               nok=0                                                     12d7s20
               do i=0,ncol-1
                if(ibc(nden1x(isc,isa,ii2e)+i).ne.0)then                7d19s21
                 ibc(nden1x(isc,isa,ii2e)+nok)=i                        7d19s21
                 nok=nok+1                                               3d22s21
                end if
               end do
               if(nok.gt.0)then                                         7d19s21
                itmpi=ibcoff                                             12d7s20
                ibcoff=itmpi+nok*nvirt(isbv)                             12d7s20
                call enough('hcsibk.  7',bc,ibc)
                jden=iden1x(isc,isa,ii2e)                               7d19s21
                jtmpi=itmpi                                              12d7s20
                do iz=itmpi,ibcoff-1                                    7d19s21
                 bc(iz)=0d0                                             7d19s21
                end do                                                  7d19s21
                do i=0,nok-1                                             3d22s21
                 icol=ibc(nden1x(isc,isa,ii2e)+i)                       7d19s21
                 ifrm=iden1x(isc,isa,ii2e)+nn*icol                      7d19s21
                 do ir=0,nroot-1                                         12d12s20
                  do j=nnlow,nnhgh                                       12d12s20
                   bc(jden)=bc(ifrm+j)                                   12d12s20
                   jden=jden+1                                           12d12s20
                  end do                                                 12d12s20
                  ifrm=ifrm+ncsf(iarg)                                   12d12s20
                 end do                                                  12d12s20
                 id=icol/nnv                                             12d7s20
                 ileft=icol-id*nnv                                       12d7s20
                 ia=ileft/irefo(isc)                                     1d25s21
                 ic=ileft-ia*irefo(isc)                                  1d25s21
                 ifad=ixmt(isa,i2eop(2,ii2e))+ia+idoubo(isa)            7d19s21
     $                 +nbasdws(isa)*(id+idoubo(isavc))                 7d19s21
                 fad=bc(ifad)                                           7d12s21
                 jint=ixmt(isbv,i2eop(1,ii2e))+idoubo(isbv)+irefo(isbv) 7d19s21
     $                +nbasdws(isbv)*(ic+idoubo(isc))                   7d12s21
                 ktmpi=jtmpi                                             7d9s21
                 do iv=0,nvirt(isbv)-1                                   7d9s21
                  bc(ktmpi)=bc(ktmpi)+fad*bc(jint+iv)                   7d12s21
                  ktmpi=ktmpi+nok                                        7d9s21
                 end do                                                  7d9s21
                 jtmpi=jtmpi+1
                end do                                                   12d7s20
                call dgemm('n','n',nnuser,nvirt(isbv),nok,phase2,        7d9s21
     $              bc(iden1x(isc,isa,ii2e)),nnuser,bc(itmpi),nok,factg,7d19s21
     $              bc(itmpg),nnuser,                                   12d12s20
     d' hcsibk.  2')
                factg=1d0                                                12d7s20
               end if                                                    12d7s20
              end do                                                    7d19s21
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
               iad1=itmpg+im+nnuse*(ir+nroot*iv)                        11d10s20
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
          call ddi_iacc(bc,ibc,nhand(nclo1p,isb),i1,i1c,i1o,in,         11d15s22
     $        bc(ighere),ibc(iacc),nacc)                                11d15s22
          nwiacc=nwiacc+i1c*(ih+1-il)                                   1d26s23
         end if                                                         3d22s21
         ibcoff=idenh                                                   3d22s21
        end if                                                          11d25s20
       end do                                                            11d25s20
      end do                                                            11d26s20
      call ddi_done(ibc(iacc),nacc)                                     3d22s21
      ibcoff=ibcoffo                                                    3d22s21
      if(ldebug)then                                                    1d25s21
       write(6,*)('all done in hcsibk')
       call dws_synca                                                    1d24s21
       write(6,*)('what we have for hcsibk:')
       do isb=1,nsymb                                                    12d7s20
        write(6,*)('for isb = '),isb
        isbv=multh(isb,isymmrci)                                         12d7s20
        i1c=nvirt(isbv)*nroot                                            11d10s20
        nhere=0                                                          12d7s20
        do ii=1,mdoop                                                   7d9s21
         iarg=ii-mdon                                                    12d7s20
         ndblock=ncsf(iarg)*nff1(ii,isb,1)                               11d10s20
         nhere=nhere+ndblock
        end do
        mdblock=nhere                                                    11d10s20
        nsing=mdblock*nvirt(isbv)                                        12d7s20
        igs=ibcoff
        ibcoff=igs+nhere*nvirt(isbv)*nroot                               11d10s20
        call enough('hcsibk.  8',bc,ibc)
        jgs=igs
        do ii=1,mdoop                                                   7d9s21
         iarg=ii-mdon                                                    12d7s20
         ndblock=ncsf(iarg)*nff1(ii,isb,1)                               11d10s20
         if(ndblock.gt.0)then
          itmp=ibcoff                                                    12d7s20
          ibcoff=itmp+ndblock*nvirt(isbv)*nroot                          11d10s20
          call enough('hcsibk.  9',bc,ibc)
          in=ndblock                                                     12d7s20
          call ddi_get(bc,ibc,nhand(ii,isb),i1,i1c,i1,in,bc(itmp))      11d15s22
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
