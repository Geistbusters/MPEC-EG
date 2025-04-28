c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcidbk(nff22,iff22,ff22,nsymb,mdon,mdoop,multh,
     $     isymmrci,nvirt,ncsf,ncsf2,nfdat,irel,ism,irefo,vd,ixw1,      8d3s21
     $     ixw2,norb,nff0,iff0,gi,ncsfv,nrootu,ldebug,sr2,srh,n2e,      8d3s21
     $     ixmt,i2eop,isymop,phase2,idoubo,nbasdws,nec,bc,ibc)          1d13s23
c
c
      implicit real*8 (a-h,o-z)                                         10d19s20
      external second                                                   3d18s21
      integer*8 itvisvc,itvisvo,itvnotvc,itvnotvo,jrefc,jrefo,ipack8,   12d8s20
     $     itestc,itesto,iff22(*),gandcc,gandco,gandcb                  2d6s23
      integer*1 idorb(64),isorb(64),nab1(2),iorb(64),imap(64),icode(64),12d8s20
     $     nab2(2)                                                      12d8s20
      logical ldebug,lbail,lchoice                                      3d17s21
      dimension nff22(mdoop,2,nsymb),multh(8,8),nvirt(*),               1d7s21
     $     ncsf(*),ncsf2(4,*),nfdat(5,4,*),irel(*),ism(*),              1d7s21
     $     ipack4(2),itmpt(4),ndenvisv(8),mdenvisv(8),nokf(8),nok(8),   1d6s21
     $     irefo(*),vd(*),ndenvnotv(4,8),mdenvnotv(4,8),ixmt(8,*),      8d3s21
     $     idenvisv(8),nl(4),npre(4),mpre(4),nokfl(4,8),nokl(4,8),      1d7s21
     $     idenvnotv(4,8),nff0(mdoop,*),iff0(*),i2eop(2,3),isymop(*),   8d3s21
     $     nab4(2,3),ipvint(4,8),idoubo(*),nbasdws(*),gi(*),            8d3s21
     $     ff22(*),ioxx(2)                                              2d6s23
      equivalence (ipack8,ipack4)                                       11d20s20
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data loop,loopx/0,8121/
      mdoo=mdoop-1                                                      8d3s21
      nrootm=nrootu-1                                                   1d6s21
      norbx=norb+1                                                      11d19s20
      norbxx=norbx+1                                                    11d19s20
      idoit=0                                                           3d8s21
      ibcoffo=ibcoff                                                    11d22s20
      igg=ibcoff                                                        1d7s21
      ibcoff=igg+ncsfv*nrootu                                           1d7s21
      call enough('hcidbk.  1',bc,ibc)
      do i=igg,ibcoff-1                                                 1d7s21
       bc(i)=0d0                                                        1d7s21
      end do                                                            1d7s21
      ioffdnon=1                                                        1d7s20
      do isb=1,nsymb                                                    1d7s20
       isbv12=multh(isb,isymmrci)                                       1d7s20
       ibctpp=ibcoff                                                    3d16s21
       do isa=1,nsymb                                                   3d16s21
        ksb=multh(isa,isbv12)                                           3d16s21
        nvv=irefo(isa)*irefo(ksb)                                       3d16s21
        do l=1,4                                                        3d16s21
         ipvint(l,isa)=ibcoff                                           3d16s21
         ibcoff=ipvint(l,isa)+nvv*nrootu*nfdat(2,l,isb)                 3d16s21
        end do                                                          3d16s21
       end do                                                           3d16s21
       call enough('hcidbk.  2',bc,ibc)
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
            call enough('hcidbk.  3',bc,ibc)
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
             nrr=irefo(isa)*irefo(isa)                                  3d16s21
             if(min(nfdat(2,1,isb),nrr).gt.0)then                       11d17s21
c     Gjr = Djkil (il|vv) Vvrk
              itmpv=ibcoff                                               3d16s21
              ibcoff=itmpv+nhere*nrr                                     3d16s21
              call enough('hcidbk.  4',bc,ibc)
              do iz=itmpv,ibcoff-1                                      8d3s21
               bc(iz)=0d0                                               8d3s21
              end do                                                    8d3s21
              itmpra=ibcoff                                             8d3s21
              ibcoff=itmpra+irefo(isa)                                  8d3s21
              call enough('hcidbk.  5',bc,ibc)
              do ii2e=1,n2e                                             8d3s21
               if(multh(isa,isbv1).eq.isymop(i2eop(1,ii2e)))then        8d3s21
                jtmpv=itmpv                                             8d3s21
                do i2=ivl,ivh                                           8d3s21
                 iv=i2-1+idoubo(isbv1)+irefo(isbv1)                     8d3s21
                 do ia=0,irefo(isa)-1                                   8d3s21
                  iaa=ia+idoubo(isa)                                    8d5s21
                  iav=ixmt(isa,i2eop(1,ii2e))+iaa+nbasdws(isa)*iv       8d3s21
                  bc(itmpra+ia)=bc(iav)                                 8d5s21
                 end do                                                 8d3s21
                 do ia=0,irefo(isa)-1                                   8d3s21
                  fact=phase2*bc(itmpra+ia)                             8d3s21
                  do ib=0,irefo(isa)-1                                  8d3s21
                   bc(jtmpv+ib)=bc(jtmpv+ib)+fact*bc(itmpra+ib)         8d3s21
                  end do                                                8d3s21
                  jtmpv=jtmpv+irefo(isa)                                8d3s21
                 end do                                                 8d3s21
                end do                                                  8d3s21
               end if                                                   8d3s21
              end do                                                    8d3s21
              ibcoff=itmpra                                             8d3s21
              ncol=nrootu*nfdat(2,1,isb)                                 3d16s21
              itmp2=ibcoff                                               3d16s21
              ibcoff=itmp2+nrr*ncol                                      3d16s21
              call enough('hcidbk.  6',bc,ibc)
              call dgemm('n','n',nrr,ncol,nhere,srh,                     3d16s21
     $               bc(itmpv),nrr,bc(itmp),nhere,0d0,                  3d16s21
     $            bc(itmp2),nrr,                                        3d16s21
     d' hcidbk.  1')
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
            call enough('hcidbk.  7',bc,ibc)
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
              ksb=multh(isa,isbv12)                                       3d16s21
              nrr=irefo(isa)*irefo(ksb)                                  3d16s21
              if(nrr.gt.0)then                                          3d16s21
               itmpv=ibcoff                                               3d16s21
               ibcoff=itmpv+mherev*nrr                                    3d16s21
               call enough('hcidbk.  8',bc,ibc)
               do iz=itmpv,ibcoff-1                                     8d3s21
                bc(iz)=0d0                                              8d3s21
               end do                                                   8d3s21
               itmpi2=ibcoff                                            8d3s21
               itmpk1=itmpi2+irefo(isa)*(i2e+1-i2s)                     8d3s21
               ibcoff=itmpk1+irefo(ksb)*nvirt(isbv1)                    8d3s21
               call enough('hcidbk.  9',bc,ibc)
               do ii2e=1,n2e                                            8d3s21
                if(multh(isa,isbv2).eq.isymop(i2eop(1,ii2e)).and.       8d3s21
     $             multh(ksb,isbv1).eq.isymop(i2eop(2,ii2e)))then       8d3s21
                 do i2=i2s,i2e                                          8d3s21
                  iv2=i2-1+idoubo(isbv2)+irefo(isbv2)                   8d3s21
                  jtmpi2=itmpi2+irefo(isa)*(i2-i2s)                     8d3s21
                  do ia=0,irefo(isa)-1                                  8d3s21
                   iaa=ia+idoubo(isa)                                   8d3s21
                   iav=ixmt(isa,i2eop(1,ii2e))+iaa+nbasdws(isa)*iv2     8d3s21
                   bc(jtmpi2+ia)=phase2*bc(iav)                         8d3s21
                  end do                                                8d3s21
                 end do                                                 8d3s21
                 do iv1=0,nvirt(isbv1)-1                                8d3s21
                  jtmpk1=itmpk1+irefo(ksb)*iv1                          8d3s21
                  ivv=iv1+idoubo(isbv1)+irefo(isbv1)                    8d3s21
                  do ib=0,irefo(ksb)-1                                  8d3s21
                   ibb=ib+idoubo(ksb)                                   8d3s21
                   ibv=ixmt(ksb,i2eop(2,ii2e))+ibb+nbasdws(ksb)*ivv     8d3s21
                   bc(jtmpk1+ib)=bc(ibv)                                8d3s21
                  end do                                                8d3s21
                 end do                                                 8d3s21
                 i10=i1s                                                    3d16s21
                 i1n=nvirt(isbv1)                                           3d16s21
                 jtmpv=itmpv                                                3d16s21
                 do i2=i2s,i2e                                              3d16s21
                  jtmpi2=itmpi2+irefo(isa)*(i2-i2s)                     8d3s21
                  i2m=i2-1                                                  3d16s21
                  if(i2.eq.i2e)i1n=i1e                                      3d16s21
                  ntop=min(i1n,i2m+isw*(nvirt(isbv1)-i2m))                  3d16s21
                  do i1=i10,ntop                                            3d16s21
                   jtmpk1=itmpk1+irefo(ksb)*(i1-1)                      8d3s21
                   i1m=i1-i10                                               3d16s21
                   do ib=0,irefo(ksb)-1                                 8d3s21
                    do ia=0,irefo(isa)-1                                8d3s21
                     bc(jtmpv+ia)=bc(jtmpv+ia)                          8d3s21
     $                    +bc(jtmpi2+ia)*bc(jtmpk1+ib)                  8d3s21
                    end do                                              8d3s21
                    jtmpv=jtmpv+irefo(isa)                              8d3s21
                   end do                                               8d3s21
                  end do                                                8d3s21
                  i10=1                                                 8d3s21
                 end do                                                 8d3s21
                end if                                                  8d3s21
               end do                                                   8d3s21
               ibcoff=itmpi2                                            8d3s21
               ncol=nrootu*nfdat(2,l,isb)                                 3d16s21
               itmp2=ibcoff                                               3d16s21
               ibcoff=itmp2+nrr*ncol                                      3d16s21
               call enough('hcidbk. 10',bc,ibc)
               call dgemm('n','n',nrr,ncol,mherev,1d0,                    3d16s21
     $            bc(itmpv),nrr,bc(itmp),mherev,0d0,                    3d16s21
     $            bc(itmp2),nrr,                                        3d16s21
     d' hcidbk.  2')
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
       do nclo0=mdon,mdoo                                                1d7s20
        nclo0p=nclo0+1                                                   1d7s20
        if(nff0(nclo0p,1).gt.0)then                                      1d7s20
         iarg0=nclo0p-mdon                                               12d7s20
         nopen=nec-nclo0*2                                               12d7s20
         ntopx=nclo0+2                                                  3d12s21
         jff0=nff0(nclo0p,2)                                            12d8s20
         do jf=1,nff0(nclo0p,1)                                         12d7s20
          jrefc=iff0(jff0)                                              12d7s20
          jff0=jff0+1                                                   12d7s20
          jrefo=iff0(jff0)                                              12d7s20
          nopen=popcnt(jrefo)
          jff0=jff0+1                                                   12d7s20
          do nclo2=max(mdon,nclo0-2),min(ntopx,mdoo)                    1d7s20
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
              idoit=mynowprog                                           9d18s24
              lchoice=.false.                                            3d17s21
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
               ipssx=0                                                  2d6s23
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
                 if(btest(jrefc,nab4(1,2)).and.                         2d6s23
     $                .not.btest(jrefc,nab4(1,1)))nbt=1                 2d6s23
                else                                                    2d6s23
                 nbt=0                                                  2d6s23
                 if(btest(itvnotvc,nab4(2,2)).and.                      2d6s23
     $                .not.btest(itvnotvc,nab4(2,1)))nbt=1              2d6s23
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
               if(nnot.eq.3)then                                           12d8s20
                ipssx=1                                                    12d8s20
               else if(nnot.eq.4)then                                    2d6s23
                ipssx=2                                                    12d8s20
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
                itestc=jrefc                                               1d6s21
                itesto=jrefo                                               1d6s21
                if(btest(itestc,nab4(1,iu1)))then                          1d6s21
                 itestc=ibclr(itestc,nab4(1,iu1))                          12d8s20
                 itesto=ibset(itesto,nab4(1,iu1))                          12d8s20
                 nopenk=nopen+1                                            1d6s21
                 karg=iarg0-1                                               12d8s20
                else if(btest(itesto,nab4(1,iu1)))then                     12d8s20
                 itesto=ibclr(itesto,nab4(1,iu1))                          12d8s20
                 nopenk=nopen-1                                            1d6s21
                 karg=iarg0                                                1d6s21
                else                                                          11d13s20
                 write(6,*)('bit not set for nab4(1,1) = '),nab4(1,iu1)       11d27s20
                 stop 'nab4(1,1)'                                           11d27s20
                end if                                                        11d13s20
                if(btest(itesto,nab4(2,iu2)))then                          12d8s20
                 itestc=ibset(itestc,nab4(2,iu2))                          12d8s20
                 itesto=ibclr(itesto,nab4(2,iu2))                          12d8s20
                 nopenk=nopenk-1                                              11d13s20
                 karg=karg+1                                                  11d13s20
                else if(btest(itesto,nab4(2,iu2)))then                     12d8s20
                 write(6,*)('already double in nab4(2,1) = '),          2d6s23
     $                nab4(2,iu2)                                       2d6s23
                 stop 'nab4(2,1)'                                           11d27s20
                else                                                          11d13s20
                 itesto=ibset(itesto,nab4(2,iu2))                          12d8s20
                 nopenk=nopenk+1                                              11d13s20
                end if                                                        11d13s20
                nqq=karg+mdon-1                                          4d20s21
                if(nqq.ge.mdon.and.nqq.le.mdoo)then                      4d20s21
                 call gandc(jrefc,jrefo,itestc,itesto,nopen,nopenk,          12d8s20
     $         iarg0,karg,ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1, 12d8s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
                 call gandc(itestc,itesto,itvnotvc,itvnotvo,nopenk,
     $              nopen2p,karg,iarg2,ncsf,norbxx,ixw1,ixw2,nnot2,nab2,
     $              iwpb2,iwpk2,ncsfmid2,bc,ibc)                        11d14s22
                 if(nnot1.eq.2.and.nnot2.eq.2)then                              11d13s20
                  if(nab1(2).lt.nab2(2))then                               1d7s21
                   icpy=nab1(1)                                            1d7s21
                   nab1(1)=nab2(1)                                         1d7s21
                   nab2(1)=icpy                                            1d7s21
                  end if                                                      12d8s20
                  jsa=ism(nab1(1))                                          1d7s21
                  jga=irel(nab1(1))-1                                       1d7s21
                  jsb=ism(nab2(1))                                          1d7s21
                  jgb=irel(nab2(1))-1                                       1d7s21
                  icol=jga+irefo(jsa)*jgb                                 1d6s21
                  nrr=irefo(jsa)*irefo(jsb)                                   1d6s21
                  iprod=ibcoff                                            3d17s21
                  if(lchoice)then                                         3d17s21
                   call genmat(ncsf(iarg0),ncsf(karg),ncsf(iarg2),            12d8s20
     $            ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod,bc,   11d10s22
     $                ibc)                                              11d10s22
                  end if                                                  3d17s21
                  ltmp2=ibcoff                                            3d18s21
                  ltmp3=ltmp2+ncsf(iarg0)*nhere                           3d18s21
                  ibcoff=ltmp3+nhere*nrootu                               3d18s21
                  call enough('hcidbk. 11',bc,ibc)
                  jprod=iprod                                               1d6s21
                  iargo=1                                                 3d12s21
                  ktmp2=ltmp2                                             3d18s21
                  ktmp3=ltmp3                                             3d18s21
                  do l=1,4                                                    12d8s20
                   if(nl(l).gt.0)then                                         12d8s20
                    nn=nfdat(2,l,isb)*nrr                                 3d18s21
                    iad1=jvcv+iff22(jvcv+5+l)                            8d3s21
                    iad2=iad1+nl(l)                                       3d19s21
                    itmp2=ktmp2                                           3d18s21
                    call enough('hcidbk. 12',bc,ibc)
                    if(lchoice)then                                       3d17s21
                     call dgemm('n','n',ncsf(iarg0),nl(l),
     $                   ncsf2(l,iarg2),1d0,bc(jprod),ncsf(iarg0),
     $                   ff22(iad2),ncsf2(l,iarg2),0d0,bc(itmp2),
     $                   ncsf(iarg0),' hcidbk.  3')
                    else                                                  3d17s21
                     call genmatn2(ncsf(iarg0),ncsf(karg),ncsf(iarg2),     3d12s21
     $                 ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,       3d12s21
     $                 ff22(iad2),ncsf2(l,iarg2),nl(l),bc(itmp2),0d0,   8d3s21
     $                 iargo,ncsf2(l,iarg2),1d0,bc,ibc)                 12d9s22
                    end if                                                3d17s21
                    itmp3=ktmp3                                           3d18s21
                    jtmp3=itmp3                                           3d18s21
                    do iii=0,nl(l)-1                                      3d18s21
                     ii=iff22(iad1+iii)-1                                8d3s21
                     jad=ipvint(l,jsa)+ii+nfdat(2,l,isb)*icol             3d18s21
                     do ir=0,nrootm                                       3d18s21
                      bc(jtmp3+nhere*ir)=bc(jad+nn*ir)                    3d18s21
                     end do                                               3d18s21
                     jtmp3=jtmp3+1                                        3d18s21
                    end do                                                3d18s21
                   end if                                                     12d8s20
                   jprod=jprod+ncsf(iarg0)*ncsf2(l,iarg2)                   1d6s21
                   iargo=iargo+ncsf2(l,iarg2)                             3d12s21
                   ktmp2=ktmp2+ncsf(iarg0)*nl(l)                          3d18s21
                   ktmp3=ktmp3+nl(l)                                      3d18s21
                  end do                                                      12d8s20
                  call dgemm('n','n',ncsf(iarg0),nrootu,nhere,1d0,        3d18s21
     $                 bc(ltmp2),ncsf(iarg0),bc(ltmp3),nhere,1d0,       3d18s21
     $                 bc(jvs),ncsfv,                                   3d18s21
     d' hcidbk.  4')
                  ibcoff=iprod                                              1d6s21
                 else
                  write(6,*)('nnots not both 2!!! '),nnot1,nnot2
                  call dws_synca
                  call dws_finalize
                  stop
                 end if                                                     12d8s20
                end if                                                   4d20s21
               end do                                                      12d8s20
              end if                                                       12d8s20
             end if                                                     3d8s21
             idoit=idoit+1                                              3d8s21
             jvcv=jvcv+nspace                                           3d19s21
            end do
           end if                                                       1d7s20
          end do                                                        1d7s20
          jvs=jvs+ncsf(iarg0)                                           1d7s20
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
