c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcidbk4(nff0,iff0,gi,ncsfv,i2sb,i2smb,mdoobp,          11d17s21
     $     nff22,iff22,ff22,nfdat,ncsf,ncsf2,vd,i2sk,i2smk,mdookp,      11d17s21
     $     isymmrci,nec,mdon,nsymb,multh,irw1,irw2,                     11d17s21
     $     nvirt,nrootu,ldebug,isymc,irorip,isopt,ism,irel,irefo,       11d17s21
     $     norb,kmatsrb,iifmx,ntype,srh,bc,ibc)                         11d9s22
c
c
      implicit real*8 (a-h,o-z)                                         10d19s20
      external second                                                   3d18s21
      integer*8 itvisvc,itvisvo,itvnotvc,itvnotvo,jrefc,jrefo,ipack8,   12d8s20
     $     itestc,itesto,iff22(*),ipack,gandcc,gandco,gandcb            2d7s23
      integer*2 ipack2(4)                                               11d17s21
      integer*1 idorb(64),isorb(64),nab1(2),iorb(64),imap(64),icode(64),12d8s20
     $     nab2(2)                                                      12d8s20
      logical ldebug,lbail,lchoice                                      3d17s21
      dimension nff22(mdookp,2,nsymb),multh(8,8),nvirt(*),               1d7s21
     $     ncsf(*),ncsf2(4,*),nfdat(5,4,*),irel(*),ism(*),              1d7s21
     $     ipack4(2),itmpt(4),ndenvisv(8),mdenvisv(8),nokf(8),nok(8),   1d6s21
     $     irefo(*),vd(*),ndenvnotv(4,8),mdenvnotv(4,8),ioxx(2),        2d7s23
     $     idenvisv(8),nl(4),npre(4),mpre(4),nokfl(4,8),nokl(4,8),      1d7s21
     $     idenvnotv(4,8),nff0(mdoobp,*),iff0(*),nab4(2,3),ipvint(4,8), 11d17s21
     $     gi(*),ff22(*),kmatsrb(8,8,8),iifmx(*),isopt(*),iwpb1(4),     11d17s21
     $     iwpk1(4),iwpb2(4),iwpk2(4),ipack2a(4),isy(4),igya(4),imy(4)  11d17s21
      equivalence (ipack8,ipack4),(ipack,ipack2)                        11d17s21
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      irori=irorip-1                                                    11d17s21
      if(isopt(3).ne.0)then
       ieoro=1-irori                                                      9d6s21
      else                                                              10d14s21
       ieoro=irori                                                      9d6s21
      end if                                                            9d6s21
      nrootm=nrootu-1                                                   1d6s21
      mdoob=mdoobp-1                                                    11d17s21
      mdook=mdookp-1                                                    11d17s21
      norbx=norb+1                                                      11d19s20
      norbxx=norbx+1                                                    11d19s20
      idoit=0                                                           3d8s21
      ibcoffo=ibcoff                                                    11d22s20
      igg=ibcoff                                                        1d7s21
      ibcoff=igg+ncsfv*nrootu                                           1d7s21
      call enough('hcidbk4.  1',bc,ibc)
      do i=igg,ibcoff-1                                                 1d7s21
       bc(i)=0d0                                                        1d7s21
      end do                                                            1d7s21
      ioffdnon=1                                                        1d7s20
      do isb=1,nsymb                                                    1d7s20
       isbv12=multh(isb,isymmrci)                                       1d7s20
       ibctpp=ibcoff                                                    3d16s21
       do isa=1,nsymb                                                   3d16s21
        ksb=multh(isopt(1),multh(isa,isbv12))                           11d18s21
        nvv=irefo(isa)*irefo(ksb)*ntype                                 11d18s21
        do l=1,4                                                        3d16s21
         ipvint(l,isa)=ibcoff                                           3d16s21
         ibcoff=ipvint(l,isa)+nvv*nrootu*nfdat(2,l,isb)                 3d16s21
        end do                                                          3d16s21
       end do                                                           3d16s21
       call enough('hcidbk4.  2',bc,ibc)
       ibcbtt=ibcoff-1                                                  3d16s21
       nwdtt=ibcoff-ibctpp                                              3d17s21
       do i=ibctpp,ibcbtt                                               3d16s21
        bc(i)=0d0                                                       3d16s21
       end do                                                           3d16s21
       joffdnon=ioffdnon                                                3d16s21
       do isbv1=1,nsymb                                                 3d16s21
        isbv2=multh(isbv1,isbv12)                                       3d16s21
        if(isbv2.ge.isbv1)then                                          3d16s21
         if(isbv12.eq.1)then                                            3d16s21
          call ilimts(nvirt(isbv1),nvirt(isbv1),mynprocg,mynowprog,     3d16s21
     $           il,ih,i1s,i1e,i2s,i2e)                                 1d6s21
          nherev=ih+1-il                                                3d16s21
          if(nherev.gt.0)then                                           3d16s21
           ivl=i2e+1                                                    3d16s21
           ivh=0                                                        3d16s21
           do i2=i2s,i2e                                                3d16s21
            icol=i2+nvirt(isbv1)*(i2-1)                                 3d16s21
            if(icol.ge.il.and.icol.le.ih)then                           3d16s21
             ivl=min(ivl,i2)                                            3d16s21
             ivh=max(ivh,i2)                                            3d16s21
            end if                                                      3d16s21
           end do                                                       3d16s21
           nhere=ivh+1-ivl                                              3d16s21
           if(nhere.gt.0)then                                           3d16s21
            itmp=ibcoff                                                 3d16s21
            ibcoff=itmp+nhere*nrootu*nfdat(2,1,isb)                     3d16s21
            call enough('hcidbk4.  3',bc,ibc)
            jtmp=itmp                                                   3d16s21
            do k=0,nfdat(2,1,isb)-1                                     3d16s21
             do ir=0,nrootm                                             3d16s21
              iadv=joffdnon+ivl-1+nvirt(isbv1)*(ir+nrootu*k)            3d16s21
              do iv=0,nhere-1                                           3d16s21
               bc(jtmp+iv)=vd(iadv+iv)                                  3d16s21
              end do                                                    3d16s21
              jtmp=jtmp+nhere                                           3d16s21
             end do                                                     3d16s21
            end do                                                      3d16s21
            do isa=1,nsymb                                              3d16s21
             ksa=multh(isa,isopt(1))                                    11d18s21
             nrr=irefo(isa)*irefo(ksa)*ntype                            11d18s21
             if(min(nfdat(2,1,isb),nrr).gt.0)then                       11d17s21
c     Gjr = Djkil (il|vv) Vvrk
              itmpv=ibcoff                                               3d16s21
              ibcoff=itmpv+nhere*nrr                                     3d16s21
              call enough('hcidbk4.  4',bc,ibc)
              do iz=itmpv,ibcoff-1                                      8d3s21
               bc(iz)=0d0                                               8d3s21
              end do                                                    8d3s21
              jtmpv=itmpv                                               11d18s21
              do i2=ivl,ivh                                             11d18s21
               irow=i2+nvirt(isbv1)*(i2-1)                              11d18s21
               kint=kmatsrb(isa,isbv1,ksa)+irow-il                      11d18s21
               do iabc=0,nrr-1                                          11d18s21
                bc(jtmpv+iabc)=bc(kint+nherev*iabc)                     11d18s21
               end do                                                   11d18s21
               jtmpv=jtmpv+nrr                                          11d18s21
              end do                                                    11d18s21
              ncol=nrootu*nfdat(2,1,isb)                                 3d16s21
              itmp2=ibcoff                                               3d16s21
              ibcoff=itmp2+nrr*ncol                                      3d16s21
              call enough('hcidbk4.  5',bc,ibc)
              call dgemm('n','n',nrr,ncol,nhere,srh,                     3d16s21
     $               bc(itmpv),nrr,bc(itmp),nhere,0d0,                  3d16s21
     $            bc(itmp2),nrr,                                        3d16s21
     d' hcidbk4.  1')
              jtmp2=itmp2                                                3d16s21
              do k=0,nfdat(2,1,isb)-1                                    3d16s21
               do ir=0,nrootm
                do i=0,nrr-1                                             3d16s21
                 iad3=ipvint(1,isa)+k+nfdat(2,1,isb)*(i+nrr*ir)          3d16s21
                 bc(iad3)=bc(iad3)+bc(jtmp2+i)                           3d16s21
                end do                                                   3d16s21
                jtmp2=jtmp2+nrr                                         3d16s21
               end do                                                    3d16s21
              end do                                                     3d16s21
              ibcoff=itmpv                                              3d16s21
             end if                                                     3d16s21
            end do                                                      3d16s21
            ibcoff=itmp                                                 3d16s21
           end if                                                       3d16s21
          end if                                                        3d16s21
          joffdnon=joffdnon+nvirt(isbv1)*nfdat(2,1,isb)*nrootu           3d16s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                          3d16s21
          isw=0                                                          3d16s21
         else                                                            3d16s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                  3d16s21
          isw=1                                                          3d16s21
         end if                                                          3d16s21
         call ilimts(nvirt(isbv1),nvirt(isbv2),mynprocg,mynowprog,       3d16s21
     $           il,ih,i1s,i1e,i2s,i2e)                                  1d7s21
         nherev=ih+1-il                                                  3d16s21
         do l=1,4                                                        3d16s21
          if(min(nvv,nfdat(2,l,isb)).gt.0)then                           3d16s21
           if(nherev.gt.0)then                                           3d16s21
            itmp=ibcoff                                                  3d16s21
            ibcoff=itmp+nherev*nrootu*nfdat(2,l,isb)                     3d16s21
            call enough('hcidbk4.  6',bc,ibc)
            jtmp=itmp                                                    3d16s21
            do k=0,nfdat(2,l,isb)-1                                      3d16s21
             do ir=0,nrootm                                              3d16s21
              ktmp=jtmp                                                  3d16s21
              i10=i1s                                                    3d16s21
              i1n=nvirt(isbv1)                                           3d16s21
              iadv=joffdnon+nvv*(ir+nrootu*k)                            3d16s21
              do i2=i2s,i2e                                              3d16s21
               i2m=i2-1                                                  3d16s21
               if(i2.eq.i2e)i1n=i1e                                      3d16s21
               ntop=min(i1n,i2m+isw*(nvirt(isbv1)-i2m))                  3d16s21
               do i1=i10,ntop                                            3d16s21
                i1m=i1-1                                                 3d16s21
                itri=((i2m*(i2m-1))/2)+i1m                               3d16s21
                irec=i1m+nvirt(isbv1)*i2m                                3d16s21
                itri=itri+isw*(irec-itri)                                3d16s21
                bc(ktmp)=vd(iadv+itri)                                   3d16s21
                ktmp=ktmp+1                                              3d16s21
               end do                                                    3d16s21
               i10=1                                                     3d16s21
              end do                                                     3d16s21
              mherev=ktmp-jtmp                                           3d16s21
              jtmp=jtmp+mherev                                           3d16s21
             end do                                                      3d16s21
            end do                                                       3d16s21
            if(mherev.gt.0)then                                          3d16s21
             do isa=1,nsymb                                               3d16s21
              ksb=multh(isopt(1),multh(isa,isbv12))                     11d18s21
              nrr=irefo(isa)*irefo(ksb)*ntype                           11d18s21
              if(nrr.gt.0)then                                          3d16s21
               itmpv=ibcoff                                               3d16s21
               ibcoff=itmpv+mherev*nrr                                    3d16s21
               call enough('hcidbk4.  7',bc,ibc)
               do iz=itmpv,ibcoff-1                                     8d3s21
                bc(iz)=0d0                                              8d3s21
               end do                                                   8d3s21
               i10=i1s                                                  11d18s21
               i1n=nvirt(isbv1)                                         11d18s21
               kint=kmatsrb(isa,isbv2,ksb)                              11d18s21
               jtmpv=itmpv                                              11d18s21
               do i2=i2s,i2e                                            11d18s21
                if(i2.eq.i2e)i1n=i1e                                    11d18s21
                i2m=i2-1                                                11d18s21
                ntop=min(i1n,i2m+isw*(nvirt(isbv1)-i2m))                  3d16s21
                kint0=kint                                              2d18s22
                do i1=i10,ntop                                          11d18s21
                 do iabc=0,nrr-1                                        11d18s21
                  bc(jtmpv+iabc)=bc(kint+iabc*nherev)                   11d18s21
                 end do                                                 11d18s21
                 jtmpv=jtmpv+nrr                                        11d18s21
                 kint=kint+1                                            11d18s21
                end do                                                  11d18s21
                kint=kint0+i1n+1-i10                                    2d18s22
                i10=1                                                   11d18s21
               end do                                                   11d18s21
               ncol=nrootu*nfdat(2,l,isb)                                 3d16s21
               itmp2=ibcoff                                               3d16s21
               ibcoff=itmp2+nrr*ncol                                      3d16s21
               call enough('hcidbk4.  8',bc,ibc)
               call dgemm('n','n',nrr,ncol,mherev,1d0,                    3d16s21
     $            bc(itmpv),nrr,bc(itmp),mherev,0d0,                    3d16s21
     $            bc(itmp2),nrr,                                        3d16s21
     d' hcidbk4.  2')
               jtmp2=itmp2                                                3d16s21
               do k=0,nfdat(2,l,isb)-1                                    3d16s21
                do ir=0,nrootm                                            3d16s21
                 do i=0,nrr-1                                             3d16s21
                  iad3=ipvint(l,isa)+k+nfdat(2,l,isb)*(i+nrr*ir)          3d16s21
                  bc(iad3)=bc(iad3)+bc(jtmp2+i)                           3d16s21
                 end do                                                   3d16s21
                 jtmp2=jtmp2+nrr                                          3d16s21
                end do                                                    3d16s21
               end do                                                     3d16s21
               ibcoff=itmpv                                             3d16s21
              end if                                                    3d16s21
             end do                                                      3d16s21
            end if                                                       3d16s21
            ibcoff=itmp                                                  3d16s21
           end if                                                       3d16s21
           joffdnon=joffdnon+nvv*nrootu*nfdat(2,l,isb)                  3d16s21
          end if                                                        3d16s21
         end do                                                         3d16s21
        end if                                                          3d16s21
       end do                                                           3d16s21
       call dws_gsumf(bc(ibctpp),nwdtt)                                 3d17s21
       jvs=igg                                                           1d7s20
       do nclo0=mdon,mdoob                                              11d18s21
        nclo0p=nclo0+1                                                   1d7s20
        if(nff0(nclo0p,1).gt.0)then                                      1d7s20
         iarg0=nclo0p-mdon                                               12d7s20
         nopen=nec-nclo0*2                                               12d7s20
         call wfetch(nopen,mdon,idum,i2sb,i2smb,ndetb,ncsfb,ivecb,      11d17s21
     $        iaorbb,idum,bc,ibc)                                       11d9s22
         ntopx=nclo0+2                                                  3d12s21
         jff0=nff0(nclo0p,2)                                            12d8s20
         do jf=1,nff0(nclo0p,1)                                         12d7s20
          jrefc=iff0(jff0)                                              12d7s20
          jff0=jff0+1                                                   12d7s20
          jrefo=iff0(jff0)                                              12d7s20
          nopen=popcnt(jrefo)
          jff0=jff0+1                                                   12d7s20
          do nclo2=max(mdon,nclo0-2),min(ntopx,mdook)                   11d18s21
           nclo2p=nclo2+1                                                         10d19s20
           if(nff22(nclo2p,1,isb).gt.0)then                                1d7s21
            iarg2=nclo2p-mdon                                                    10d19s20
            jvcv=nfdat(5,1,isb)+nff22(nclo2p,2,isb)                     8d3s21
            ipack8=iff22(jvcv)                                          8d3s21
            nopen2=popcnt(ipack4(2))                                    8d3s21
            nopen2p=nopen2+2                                                 12d7s20
            do if=1,nff22(nclo2p,1,isb)                                   1d7s21
             ipack8=iff22(jvcv)                                         8d3s21
             nclot=popcnt(ipack4(1))                                        11d20s20
             itvnotvc=ipack4(1)                                            12d7s20
             itvisvc=ibset(itvnotvc,norbx)                                 12d7s20
             itvisvo=ipack4(2)                                             12d7s20
             iarg2p=iarg2+1                                                12d8s20
             itvnotvo=ibset(itvisvo,norbx)                                 12d7s20
             itvnotvo=ibset(itvnotvo,norbxx)                               12d7s20
             nspace=iff22(jvcv+1)                                       8d3s21
c
c     v not v
c
             if(mod(idoit,mynprocg).eq.mynowprog)then                   3d12s21
              nhere=0                                                    3d19s21
              do l=1,4                                                   3d17s21
               nl(l)=iff22(jvcv+1+l)                                    8d3s21
               nhere=nhere+nl(l)                                         3d19s21
               if(nl(l).gt.ncsf2(l,iarg2))lchoice=.true.                 3d17s21
              end do                                                     3d17s21
              gandcc=ieor(jrefc,itvnotvc)                               2d6s23
              gandco=ieor(jrefo,itvnotvo)                               2d6s23
              gandcb=ior(gandcc,gandco)                                 2d6s23
              ndifb=popcnt(gandcb)                                      10d20s22
              if(ndifb.le.4)then                                        10d20s22
               ndifs=popcnt(gandco)                                      10d13s22
               ndifd=popcnt(gandcc)                                      10d13s22
               nnot=0                                                   2d6s23
               if(ndifs.eq.4.and.ndifb.eq.4)then                        2d6s23
                nnot=4                                                  2d6s23
                ioxx(1)=1                                               2d6s23
                ioxx(2)=1                                               2d6s23
                do i=1,norbxx                                           2d6s23
                 if(btest(gandcb,i))then                                2d6s23
                  if((btest(itvnotvc,i).and.btest(jrefo,i)).or.         2d6s23
     $                (btest(itvnotvo,i).and..not.btest(jrefc,i)))then  2d6s23
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
                do i=1,norbxx                                           2d6s23
                 if(btest(gandcb,i))then                                2d6s23
                  if(btest(gandcc,i).and.                               2d6s23
     $        ((btest(jrefc,i).and..not.btest(itvnotvo,i)).or.          2d6s23
     $         (btest(itvnotvc,i).and..not.btest(jrefo,i))))then        2d6s23
                   if(btest(itvnotvc,i))iswap=1                         2d6s23
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
                 if(btest(itvnotvc,nab4(2,2)).and.                      2d6s23
     $                .not.btest(itvnotvc,nab4(2,1)))nbt=1              2d6s23
                else                                                    2d6s23
                 nbt=0                                                  2d6s23
                 if(btest(jrefc,nab4(1,2)).and.                         2d6s23
     $                .not.btest(jrefc,nab4(1,1)))nbt=1                 2d6s23
                end if                                                  2d6s23
                if(nbt.ne.0)then                                        2d6s23
                 nab4(1,1)=nab4(1,2)                                    2d6s23
                 nab4(2,1)=nab4(2,2)                                    2d6s23
                end if                                                  2d6s23
               else if(ndifs.eq.0.and.ndifd.eq.2)then                   2d6s23
                nnot=3                                                  2d6s23
                do i=1,norbxx                                           2d6s23
                 if(btest(gandcb,i))then                                2d6s23
                  if(btest(jrefc,i))then                                2d6s23
                   nab4(1,1)=i                                          2d6s23
                   nab4(1,2)=i                                          2d6s23
                  else                                                  2d6s23
                   nab4(2,1)=i                                          2d6s23
                   nab4(2,2)=i                                          2d6s23
                  end if                                                2d6s23
                 end if                                                 2d6s23
                end do                                                  2d6s23
               end if                                                   2d6s23
               if(nnot.ne.0)then                                        2d7s23
                iu1=1
                iu2=1
                itestc=jrefc                                               1d6s21
                itesto=jrefo                                               1d6s21
                if(btest(itestc,nab4(1,iu1)))then                          1d6s21
                 itestc=ibclr(itestc,nab4(1,iu1))                          12d8s20
                 itesto=ibset(itesto,nab4(1,iu1))                          12d8s20
                 nopenk=nopen+1                                            1d6s21
                else if(btest(itesto,nab4(1,iu1)))then                     12d8s20
                 itesto=ibclr(itesto,nab4(1,iu1))                          12d8s20
                 nopenk=nopen-1                                            1d6s21
                else                                                          11d13s20
                 write(6,*)('bit not set for nab4(1,1) = '),nab4(1,iu1)       11d27s20
                 stop 'nab4(1,1)'                                           11d27s20
                end if                                                        11d13s20
                if(btest(itesto,nab4(2,iu2)))then                          12d8s20
                 itestc=ibset(itestc,nab4(2,iu2))                          12d8s20
                 itesto=ibclr(itesto,nab4(2,iu2))                          12d8s20
                 nopenk=nopenk-1                                              11d13s20
                else if(btest(itesto,nab4(2,iu2)))then                     12d8s20
                 write(6,*)('already double in nab4(2,1) = '),          2d7s23
     $                nab4(2,iu2)                                       2d7s23
                 stop 'nab4(2,1)'                                           11d27s20
                else                                                          11d13s20
                 itesto=ibset(itesto,nab4(2,iu2))                          12d8s20
                 nopenk=nopenk+1                                              11d13s20
                end if                                                        11d13s20
                call gandcr(jrefc,jrefo,itestc,itesto,nopen,              11d17s21
     $               nopenk,norbxx,nnot1,nab1,icode,imap,nx,irw1,irw2,  11d15s21
     $               iwpb1,iwpk1,bc,ibc)                                11d14s22
                call gandcr(itestc,itesto,itvnotvc,itvnotvo,nopenk,       11d17s21
     $             nopen2p,norbxx,nnot2,nab2,icode,imap,nx,irw1,irw2,   11d17s21
     $             iwpb2,iwpk2,bc,ibc)                                  11d14s22
                if(nnot1.eq.2.and.nnot2.eq.2)then                              11d13s20
                 iunit=ibcoff                                             11d17s21
                 ibcoff=iunit+ncsf(iarg2)*ncsf(iarg2)                     11d17s21
                 call enough('hcidbk4.  9',bc,ibc)
                 do iz=iunit,ibcoff-1                                     11d17s21
                  bc(iz)=0d0                                              11d17s21
                 end do                                                   11d17s21
                 do iz=0,ncsf(iarg2)-1
                  iad=iunit+iz*(ncsf(iarg2)+1)                            11d17s21
                  bc(iad)=1d0                                             11d17s21
                 end do                                                   11d17s21
                 call spinloop(i2sb,i2smb,i2sk,i2smk,nopen,nopen2p,       11d17s21
     $               nopenk,ncsfb,ncsf(iarg2),itype,imatx,ntypeq,nab1,  11d17s21
     $              iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(iunit),ncsf(iarg2), 11d17s21
     $              ieoro,bc,ibc)                                       11d14s22
                 do iq=0,ntypeq-1                                         11d17s21
                  ipack=ibc(itype+iq)                                     11d17s21
                  itmpxq=imatx+ncsfb*ncsf(iarg2)*iq                      11d18s21
                  do j=1,4                                                11d17s21
                   ipa=ipack2(j)                                          11d17s21
                   ipa=iabs(ipa)                                          11d17s21
                   ipack2a(j)=ipa                                         11d17s21
                   if(ipa.le.norb)then                                    11d17s21
                    if(ipack2(j).gt.0)then                                9d9s21
                     isy(j)=ism(ipack2(j))                                9d9s21
                     igya(j)=irel(ipack2(j))-1                            11d15s21
                     imy(j)=0                                             10d12s21
                    else                                                  9d9s21
                     isy(j)=ism(-ipack2(j))                                9d9s21
                     igya(j)=irel(-ipack2(j))-1                           10d12s21
                     imy(j)=1                                             10d12s21
                    end if                                                9d9s21
                   else                                                   11d17s21
                    igya(j)=-1                                            11d17s21
                    if(ipack2(j).gt.0)then                                11d17s21
                     imy(j)=0                                             11d17s21
                    else                                                  11d17s21
                     imy(j)=1                                             11d17s21
                    end if                                                11d17s21
                   end if                                                 11d17s21
                  end do                                                  11d17s21
                  if(ipack2a(2).gt.norb.and.ipack2a(4).gt.norb)then      11d18s21
c
c     (1v|3v')=+/-(v1|v'3)
c
                   ltest=imy(2)+2*(imy(1)+2*(imy(4)+2*imy(3)))           11d18s21
                   jtest=1-imy(2)+2*(1-imy(1)+2*(1-imy(4)+2*(1-imy(3)))) 11d18s21
                   itestp=ltest+1                                        11d18s21
                   jtestp=jtest+1                                        11d18s21
                   phsj=1d0                                              11d18s21
                   if(isopt(2).ne.0)phsj=-phsj                           11d18s21
                   iusej=-1                                              11d18s21
                   if(iifmx(itestp).ge.0)then                            11d18s21
                    iusej=iifmx(itestp)                                  11d18s21
                   else if(iifmx(jtestp).ge.0)then                       11d18s21
                    iusej=iifmx(jtestp)                                  11d18s21
                    if(isopt(4).ne.0)phsj=-phsj                          11d18s21
                   end if                                                11d18s21
                   iusek=-1                                              11d18s21
                   if(nnot.eq.4)then                                     11d18s21
c
c     (1v'|3v)=(3v|1v')=+/-(v3|v'1)
c
                    ltest=imy(2)+2*(imy(3)+2*(imy(4)+2*imy(1)))           11d18s21
                    jtest=1-imy(2)+2*(1-imy(3)+2*(1-imy(4)              2d7s23
     $                   +2*(1-imy(1))))                                2d7s23
                    itestp=ltest+1                                        11d18s21
                    jtestp=jtest+1                                        11d18s21
                    phsk=-1d0                                              11d18s21
                    if(isopt(2).ne.0)phsk=-phsk                           11d18s21
                    if(iifmx(itestp).ge.0)then                            11d18s21
                     iusek=iifmx(itestp)                                  11d18s21
                    else if(iifmx(jtestp).ge.0)then                       11d18s21
                     iusek=iifmx(jtestp)                                  11d18s21
                     if(isopt(4).ne.0)phsk=-phsk                          11d18s21
                    end if                                                11d18s21
                   end if                                                11d18s21
                   if(max(iusej,iusek).ge.0)then                         11d18s21
                    if(iusej.ge.0)then                                   11d18s21
                     icolj=igya(1)+irefo(isy(1))*(igya(3)+irefo(isy(3))  11d18s21
     $                  *iusej)                                         11d18s21
                     ksb=multh(isopt(1),multh(isy(1),isbv12))            2d3s22
                     nrr1=irefo(isy(1))*irefo(ksb)*ntype                  2d3s22
                    end if                                               11d18s21
                    if(iusek.ge.0)then                                   11d18s21
                     icolk=igya(3)+irefo(isy(3))*(igya(1)+irefo(isy(1))  11d18s21
     $                  *iusek)                                         11d18s21
                     ksb=multh(isopt(1),multh(isy(3),isbv12))              2d3s22
                     nrr3=irefo(isy(3))*irefo(ksb)*ntype                 2d3s22
                    end if                                               11d18s21
                    jprod=itmpxq                                         11d18s21
                    ltmp2=ibcoff                                            3d18s21
                    ltmp3=ltmp2+ncsfb*nhere                                11d17s21
                    ibcoff=ltmp3+nhere*nrootu                               3d18s21
                    call enough('hcidbk4. 10',bc,ibc)
                    do iz=ltmp3,ibcoff-1                                 11d18s21
                     bc(iz)=0d0                                          11d18s21
                    end do                                               11d18s21
                    ktmp2=ltmp2                                          11d18s21
                    ktmp3=ltmp3                                          11d18s21
                    do l=1,4                                                    12d8s20
                     if(nl(l).gt.0)then                                         12d8s20
                      iad1=jvcv+iff22(jvcv+5+l)                            8d3s21
                      iad2=iad1+nl(l)                                       3d19s21
                      itmp2=ktmp2                                           3d18s21
                      call enough('hcidbk4. 11',bc,ibc)
                      call dgemm('n','n',ncsfb,nl(l),ncsf2(l,iarg2),      11d17s21
     $                 1d0,bc(jprod),ncsfb,ff22(iad2),                  11d17s21
     $                 ncsf2(l,iarg2),0d0,bc(itmp2),ncsfb,              11d17s21
     d' hcidbk4.  3')
                      itmp3=ktmp3                                           3d18s21
                      jtmp3=itmp3                                           3d18s21
                      if(iusej.ge.0)then                                 11d18s21
                       nn1=nfdat(2,l,isb)*nrr1                           1d13s23
                       do iii=0,nl(l)-1                                      3d18s21
                        ii=iff22(iad1+iii)-1                                8d3s21
                        jad=ipvint(l,isy(1))+ii+nfdat(2,l,isb)*icolj     11d18s21
                        do ir=0,nrootm                                       3d18s21
                         bc(jtmp3+nhere*ir)=bc(jad+nn1*ir)*phsj          2d3s22
                        end do                                               3d18s21
                        jtmp3=jtmp3+1                                        3d18s21
                       end do                                                3d18s21
                      end if                                                     12d8s20
                      jtmp3=itmp3                                           3d18s21
                      if(iusek.ge.0)then                                 11d18s21
                       nn3=nfdat(2,l,isb)*nrr3                           1d13s23
                       do iii=0,nl(l)-1                                      3d18s21
                        ii=iff22(iad1+iii)-1                                8d3s21
                        jad=ipvint(l,isy(3))+ii+nfdat(2,l,isb)*icolk     11d18s21
                        do ir=0,nrootm                                       3d18s21
                         bc(jtmp3+nhere*ir)=bc(jtmp3+nhere*ir)+          11d18s21
     $                       bc(jad+nn3*ir)*phsk                        2d3s22
                        end do                                               3d18s21
                        jtmp3=jtmp3+1                                        3d18s21
                       end do                                                3d18s21
                      end if                                                     12d8s20
                     end if                                              11d18s21
                     jprod=jprod+ncsfb*ncsf2(l,iarg2)                      11d17s21
                     ktmp2=ktmp2+ncsfb*nl(l)                               11d17s21
                     ktmp3=ktmp3+nl(l)                                      3d18s21
                    end do                                                      12d8s20
                    call dgemm('n','n',ncsfb,nrootu,nhere,1d0,             11d17s21
     $                 bc(ltmp2),ncsfb,bc(ltmp3),nhere,1d0,             11d17s21
     $                 bc(jvs),ncsfv,                                   3d18s21
     d' hcidbk4.  4')
                    ibcoff=ltmp2                                          11d18s21
                   end if                                                11d18s21
                  else                                                   11d18s21
                   write(6,*)('don''t know how to deal with '),ipack2a,  11d18s21
     $                 ('!!!')                                          11d18s21
                   stop 'hcidbk4'                                        11d18s21
                  end if                                                 11d18s21
                 end do                                                   11d17s21
                 ibcoff=iunit                                             11d17s21
                else
                 write(6,*)('nnots not both 2!!! '),nnot1,nnot2
                 call dws_synca
                 call dws_finalize
                 stop 'hcidbk4'                                         2d7s23
                end if                                                     12d8s20
               end if                                                       12d8s20
              end if                                                     3d8s21
             end if                                                     2d7s23
             idoit=idoit+1                                              3d8s21
             jvcv=jvcv+nspace                                           3d19s21
            end do
           end if                                                       1d7s20
          end do                                                        1d7s20
          jvs=jvs+ncsfb                                                 11d17s21
         end do                                                         1d7s20
        end if                                                          1d7s20
       end do                                                           1d7s20
       ibcoff=ibctpp                                                    3d16s21
       do isbv1=1,nsymb                                                 1d7s21
        isbv2=multh(isbv1,isbv12)                                       1d7s21
        if(isbv2.ge.isbv1)then                                          1d7s21
         if(isbv12.eq.1)then                                            1d7s21
          ioffdnon=ioffdnon+nrootu*nvirt(isbv1)*nfdat(2,1,isb)          1d7s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         1d7s21
         else                                                           1d7s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 1d7s21
         end if                                                         1d7s21
         do l=1,4                                                       1d7s21
          ioffdnon=ioffdnon+nvv*nrootu*nfdat(2,l,isb)                   1d7s21
         end do                                                         1d7s21
        end if                                                          1d7s21
       end do                                                           1d7s21
      end do
      jrefc=1                                                           1d8s21
      jrefo=ncsfv*nrootu                                                1d8s21
      jgg=igg-1                                                         8d3s21
      do i=1,ncsfv*nrootu                                               8d3s21
       gi(i)=gi(i)+bc(jgg+i)                                            8d3s21
      end do                                                            8d3s21
      ibcoff=ibcoffo                                                    11d23s20
      return                                                            10d19s20
      end                                                               10d19s20
