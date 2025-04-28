c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcdibkuc(nhand,nff2,iff2,nff0,iff0,vec,ncsfv,ncsf,     7d22s21
     $     ncsf2,nec,mdon,mdoop,nsymb,multh,ixw1,ixw2,isymop,n2e,i2eop,  7d22s21
     $     phase2,nvirt,izero,ldebug,maxbxd,srh,isymmrci,ism,irel,      7d22s21
     $     idoubo,irefo,nbasdws,ixmt,norb,nroot,bc,ibc)                 11d10s22
      implicit real*8 (a-h,o-z)
      integer*1 nab1(2),nab2(2)
      integer*8 nhand(mdoop,nsymb),i1,in,i2c,i2o,i0c,i0o,               8d21s21
     $     itesta,itestb,last8(2)                                       6d17s21
      logical ldebug                                                    1d25s21
      dimension nff2(mdoop,nsymb,2),iff2(*),nff0(mdoop,3),              7d22s21
     $     iff0(*),vec(ncsfv,*),ncsf(*),multh(8,8),nab4(2,3),idorb1(32),11d25s20
     $     isorb1(32),idorb(32),isorb(32),itest(32,4),isymop(*),        7d22s21
     $     nvirt(*),iden(8,8),nden(8,8),ncsf2(*),mden(8,8),itmpx(4,8),  7d22s21
     $     jtmpx(4,8),ntmpx(4,8),mtmpx(4),i2eop(2,*),ism(*),irel(*),    7d22s21
     $     idoubo(*),irefo(*),nbasdws(*),ixmt(8,*)                      7d22s21
      include "common.store"                                            11d25s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data loopx/2600000/
      mdoo=mdoop-1                                                      7d22s21
      if(ldebug)write(6,*)('hi, I am the new and improved hcdibkuc!')       1d25s21
      xnan=-2d0
      nrootm=nroot-1                                                    6d8s21
      last8(1)=-1                                                       6d17s21
      loop=0
      ibcoffo=ibcoff                                                    12d8s20
      iacc=ibcoff                                                       6d17s21
      ighere=iacc+mynprocg                                              3d22s21
      ibcoff=ighere+maxbxd                                              6d17s21
      call enough('hcdibkuc.  1',bc,ibc)
      nacc=0                                                            3d22s21
      if(ldebug)write(6,*)('izero in hcdiuc '),izero                      1d25s21
      i1=1                                                              11d25s20
      if(izero.ne.0)then                                                1d24s21
       do isb=1,nsymb                                                    11d26s20
        do ii=mdon+1,mdoop                                              7d22s21
         if(nhand(ii,isb).ge.0)then                                     8d21s21
          call ddi_zero(bc,ibc,nhand(ii,isb))                           11d15s22
         end if                                                          11d26s20
        end do                                                          6d8s21
       end do                                                            11d26s20
       call dws_synca                                                    6d11s21
      end if                                                            1d24s21
      norbx=norb+1                                                      11d25s20
      norbxx=norbx+1                                                    6d7s21
      ifirsttime=0                                                      11d26s20
      jff2=1                                                            6d8s21
      do isb=1,nsymb                                                    11d25s20
       isbv12=multh(isb,isymmrci)                                       6d7s21
       do nclo2p=mdon+1,mdoop                                           7d22s21
        if(nff2(nclo2p,isb,1).gt.0)then                                 7d22s21
         nclo2=nclo2p-1                                                  11d25s20
         if(isb.eq.isymmrci)then                                         12d7s20
          ntop=nclo2p+3                                                  12d7s20
         else                                                            12d7s20
          ntop=nclo2p+2                                                  12d7s20
         end if                                                          12d7s20
         nopen2=nec-nclo2*2                                             11d25s20
         iarg=nclo2p-mdon                                                11d25s20
         nvall=0                                                        6d8s21
         do isbv1=1,nsymb                                               6d8s21
          isbv2=multh(isbv12,isbv1)                                     6d8s21
          if(isbv2.ge.isbv1)then                                        6d8s21
           if(isbv1.eq.isbv2)then
            nvall=nvall+ncsf2(iarg)*nvirt(isbv1)                        7d22s21
            nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       6d8s21
           else                                                         6d8s21
            nvv=nvirt(isbv1)*nvirt(isbv2)                               6d8s21
           end if                                                       6d8s21
           nvall=nvall+nvv*ncsf(iarg)                                   6d8s21
          end if                                                        6d8s21
         end do                                                         6d8s21
         nwds=nvall*nroot*nff2(nclo2p,isb,1)                            7d22s21
         igflag=0                                                       3d22s21
         nn=ncsf(iarg)*nroot                                            11d26s20
         ncsfm=ncsf(iarg)-1                                             6d8s21
         nnm=nn-1                                                       11d26s20
         idenh=ibcoff                                                   11d25s20
         do isa=1,nsymb                                                 11d26s20
          isc=multh(isa,isbv12)                                          11d26s20
          iden(isc,isa)=ibcoff                                          6d7s21
          nden(isc,isa)=iden(isc,isa)+nn*irefo(isa)*irefo(isc)          6d7s21
          ibcoff=nden(isc,isa)+irefo(isa)*irefo(isc)                    6d7s21
         end do                                                         11d26s20
         call enough('hcdibkuc.  2',bc,ibc)
         ibctop=ibcoff-1                                                11d26s20
         jghere=ighere                                                  6d8s21
         do if2=1,nff2(nclo2p,isb,1)                                    7d22s21
          do i=idenh,ibctop                                             11d26s20
           bc(i)=0d0                                                    11d26s20
          end do                                                        11d26s20
          i2c=iff2(jff2)                                                11d25s20
          jff2=jff2+1                                                   11d25s20
          i2o=iff2(jff2)                                                11d25s20
          i2o=ibset(i2o,norbx)                                          6d7s21
          i2o=ibset(i2o,norbxx)                                         6d7s21
          jff2=jff2+1                                                   11d25s20
          do nclop=max(mdon+1,nclo2p-2),min(mdoop,ntop)                 7d22s21
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
            call gandc4(i2c,i2o,i0c,i0o,nopen2,nopen,norbxx,nnot,nab4,  11d14s22
     $           bc,ibc)                                                11d14s22
            if(nnot.ne.0)then
             if(nab4(1,1).eq.0)then                                     11d27s20
              write(6,*)('gandc4 falure!!! '),nab4
              call dcbit(i2c,norbxx,'i2c')
              call dcbit(i0c,norbxx,'i0c')
              call dcbit(i2o,norbxx,'i2o')
              call dcbit(i0o,norbxx,'i0o')
              call dws_synca                                            11d27s20
              call dws_finalize                                         11d27s20
             end if                                                     11d27s20
             if(nnot.eq.3)then                                           1d25s21
              ipssx=1                                                    1d25s21
             else                                                        1d25s21
              ipssx=3                                                    1d25s21
             end if                                                      1d25s21
             do ipss=1,ipssx                                             1d25s21
              if(ipss.eq.1)then                                          1d25s21
               iu1=1                                                     1d25s21
               iu2=1                                                     1d25s21
              else if(ipss.eq.2)then                                     1d25s21
               iu1=1                                                     1d25s21
               iu2=2                                                     1d25s21
              else                                                       1d25s21
               iu1=2                                                     1d25s21
               iu2=1                                                     1d25s21
              end if                                                     1d25s21
c     doub is bra
              itesta=i2c                                                 11d26s20
              itestb=i2o                                                 11d26s20
              if(btest(itesta,nab4(1,iu1)))then                          1d25s21
               itesta=ibclr(itesta,nab4(1,iu1))                          1d25s21
               itestb=ibset(itestb,nab4(1,iu1))                          1d25s21
               nopenk=nopen2+1                                              11d13s20
               karg=iarg-1                                                  11d13s20
              else if(btest(itestb,nab4(1,iu1)))then                     1d25s21
               itestb=ibclr(itestb,nab4(1,iu1))                          1d25s21
               nopenk=nopen2-1                                              11d13s20
               karg=iarg                                                  11d13s20
              end if                                                        11d13s20
              if(btest(itestb,nab4(2,iu2)))then                          1d25s21
               itesta=ibset(itesta,nab4(2,iu2))                          1d25s21
               itestb=ibclr(itestb,nab4(2,iu2))                          1d25s21
               nopenk=nopenk-1                                              11d13s20
               karg=karg+1                                                  11d13s20
              else                                                          11d13s20
               itestb=ibset(itestb,nab4(2,iu2))                          1d25s21
               nopenk=nopenk+1                                              11d13s20
              end if                                                        11d13s20
              call gandc(i2c,i2o,itesta,itestb,nopen2,nopenk,            11d26s20
     $         iarg,karg,ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,   11d27s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
              call gandc(itesta,itestb,i0c,i0o,nopenk,nopen,             11d26s20
     $         karg,jarg,ncsf,norbxx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,   11d27s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
              if(nnot1.eq.2.and.nnot2.eq.2)then                          1d25s21
               if(nab1(1).ne.norbx.or.nab2(1).ne.norbxx)then
               end if
               jsa=ism(nab2(2))                                         6d8s21
               jga=irel(nab2(2))-1                                      6d8s21
               jsb=ism(nab1(2))                                         6d8s21
               jgb=irel(nab1(2))-1                                      6d8s21
               jden=iden(jsa,jsb)+nn*(jga+irefo(jsa)*jgb)               6d8s21
               mdenx=nden(jsa,jsb)+jga+irefo(jsa)*jgb                    6d8s21
               ibc(mdenx)=1                                              6d8s21
               call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),ncsfmid1,  6d8s21
     $              iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,vec(jvs,1),ncsfv,  6d8s21
     $              nroot,bc(jden),1d0,bc,ibc)                          11d10s22
               if(ipss.eq.2)go to 1                                      1d25s21
              end if                                                     1d25s21
             end do                                                      1d25s21
    1        continue                                                    1d25s21
            end if                                                       11d25s20
            jvs=jvs+ncsf(jarg)                                           11d26s20
           end do                                                        11d25s20
          end do                                                         11d25s20
          nokt=0                                                        6d8s21
          do isa=1,nsymb                                                6d8s21
           isc=multh(isa,isbv12)                                        6d8s21
           nok=0                                                        6d8s21
           jden=iden(isa,isc)                                           6d8s21
           itmp=ibcoff                                                  6d8s21
           ibcoff=itmp+nn*irefo(isa)*irefo(isc)                         6d8s21
           call enough('hcdibkuc.  3',bc,ibc)
           kden=itmp                                                    6d8s21
           do i=0,irefo(isa)*irefo(isc)-1                               6d8s21
            if(ibc(nden(isa,isc)+i).ne.0)then                           6d8s21
             ibc(nden(isa,isc)+nok)=i                                   6d8s21
             nok=nok+1                                                  6d8s21
             do j=0,nnm                                                 6d8s21
              bc(kden+j)=bc(jden+j)                                     6d8s21
             end do                                                     6d8s21
             kden=kden+nn                                               6d8s21
            end if                                                      6d8s21
            jden=jden+nn                                                6d8s21
           end do                                                       6d8s21
           mden(isa,isc)=nok                                            6d8s21
           nokt=nokt+nok                                                6d8s21
           if(nok.gt.0)then                                             6d8s21
c     dimensions are ncsf,nroot,iab. transpose to iab,nroot,ncsf.
            do i=0,nok-1                                                6d8s21
             do j=0,nrootm                                              6d8s21
              do k=0,ncsfm                                              6d8s21
               kji=itmp+k+ncsf(iarg)*(j+nroot*i)                        6d8s21
               ijk=iden(isa,isc)+i+nok*(j+nroot*k)                      6d8s21
               bc(ijk)=bc(kji)                                          6d8s21
              end do                                                    6d8s21
             end do                                                     6d8s21
            end do                                                      6d8s21
           end if                                                       6d8s21
           ibcoff=itmp                                                  6d8s21
          end do                                                        6d8s21
          if(nokt.gt.0)then                                             6d8s21
           if(igflag.eq.0)then                                          6d8s21
            call ddi_done(ibc(iacc),nacc)                               6d8s21
            do i=0,nwds-1                                               6d8s21
             bc(ighere+i)=0d0                                           6d8s21
            end do                                                      6d8s21
            igflag=1                                                    6d8s21
           end if                                                       6d8s21
           do isbv1=1,nsymb                                              6d8s21
            isbv2=multh(isbv1,isbv12)                                    6d8s21
            if(isbv2.ge.isbv1)then                                       6d8s21
             call ilimts(nvirt(isbv1),nvirt(isbv2),mynprocg,mynowprog,   6d8s21
     $          il,ih,i1s,i1e,i2s,i2e)                                  6d8s21
             if(isbv1.eq.isbv2)then                                      6d8s21
              nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                      6d8s21
              nvx=nvirt(isbv1)*ncsf2(iarg)                              7d22s21
             else                                                        6d8s21
              nvv=nvirt(isbv1)*nvirt(isbv2)                              6d8s21
              nvx=0                                                     6d8s21
             end if                                                      6d8s21
             nvx=nvx+nvv*ncsf(iarg)                                     6d8s21
             nhere=ih+1-il                                               6d8s21
             if(nhere.gt.0)then                                          6d8s21
              nall=nhere*nn                                              6d8s21
              itmp=ibcoff                                                6d8s21
              ibcoff=itmp+nall                                           6d8s21
              call enough('hcdibkuc.  4',bc,ibc)
              do i=itmp,ibcoff-1                                         6d8s21
               bc(i)=0d0                                                 6d8s21
              end do                                                     6d8s21
              nhit=0                                                     6d8s21
              do isa=1,nsymb                                              6d8s21
               isc=multh(isa,isbv12)                                      6d8s21
               if(mden(isa,isc).gt.0)then                                6d8s21
                nhit=1                                                   6d8s21
                itmpi=ibcoff                                             6d8s21
                ibcoff=itmpi+nhere*mden(isa,isc)                         6d8s21
                call enough('hcdibkuc.  5',bc,ibc)
                do iz=itmpi,ibcoff-1                                    7d22s21
                 bc(iz)=0d0                                             7d22s21
                end do                                                  7d22s21
                do i=0,mden(isa,isc)-1                                  6d8s21
                 ic=ibc(nden(isa,isc)+i)/irefo(isa)                     7d22s21
                 ia=ibc(nden(isa,isc)+i)-ic*irefo(isa)                  7d22s21
                 ic=ic+idoubo(isc)                                      7d22s21
                 ia=ia+idoubo(isa)                                      7d22s21
                 i10=i1s                                                7d22s21
                 i1n=nvirt(isbv1)                                       7d22s21
                 icolt=itmpi+nhere*i                                    6d8s21
                 do i2=i2s,i2e                                          7d22s21
                  i2p=i2+idoubo(isbv2)+irefo(isbv2)-1                   7d22s21
                  if(i2.eq.i2e)i1n=i1e                                  7d22s21
                  do i1=i10,i1n                                         7d22s21
                   i1p=i1+idoubo(isbv1)+irefo(isbv1)-1                  7d22s21
                   do ii2e=1,n2e                                        7d22s21
                    if(multh(isbv1,isc).eq.isymop(i2eop(1,ii2e)))then   7d22s21
                     ivc=ixmt(isbv1,i2eop(1,ii2e))+i1p                  7d22s21
     $                    +nbasdws(isbv1)*ic                            7d22s21
                     iva=ixmt(isbv2,i2eop(2,ii2e))+i2p+nbasdws(isbv2)*ia7d23s21
                     bc(icolt)=bc(icolt)+phase2*bc(ivc)*bc(iva)         7d23s21
                    end if                                              7d22s21
                   end do                                               7d22s21
                   icolt=icolt+1                                        7d22s21
                  end do                                                7d22s21
                  i10=1                                                 7d22s21
                 end do                                                 7d22s21
                end do                                                  6d8s21
                call dgemm('n','n',nhere,nn,mden(isa,isc),1d0,           6d8s21
     $              bc(itmpi),nhere,bc(iden(isa,isc)),mden(isa,isc),1d0,6d8s21
     $              bc(itmp),nhere,                                     6d8s21
     d' hcdibkuc.  1')
                ibcoff=itmpi                                             6d8s21
               end if                                                    6d8s21
              end do                                                      6d8s21
              if(nhit.ne.0)then                                          6d8s21
               kghere=jghere                                            6d8s21
               isw=1                                                    6d8s21
               if(isbv1.eq.isbv2)then                                   6d8s21
                isw=0                                                   6d8s21
                do ivp=i2s,i2e                                          6d8s21
                 iv=ivp-1                                               6d8s21
                 icol=ivp+nvirt(isbv1)*iv                               6d8s21
                 if(icol.ge.il.and.icol.le.ih)then                      6d8s21
                  lghere=kghere+nroot*ncsf2(iarg)*iv                    7d22s21
                  icol=icol-il                                          6d8s21
                  jtmp=itmp+icol                                        6d11s21
                  do i=0,nroot*ncsf2(iarg)-1                            7d22s21
                   bc(lghere+i)=bc(jtmp+nhere*i)*srh                    6d9s21
                  end do                                                6d8s21
                 end if                                                 6d8s21
                end do                                                  6d8s21
                kghere=kghere+nvirt(isbv1)*nroot*ncsf2(iarg)            7d22s21
               end if                                                   6d8s21
               i10=i1s                                                  6d8s21
               i1n=nvirt(isbv1)                                         6d8s21
               jtmp=itmp                                                6d8s21
               do i2=i2s,i2e                                            6d8s21
                iv2=i2-1                                                6d8s21
                if(i2.eq.i2e)i1n=i1e                                    6d8s21
                itop=min(i1n,iv2+isw*(nvirt(isbv1)-iv2))                6d8s21
                do i1=i10,itop                                          6d8s21
                 iv1=i1-1                                               6d8s21
                 itri=((iv2*(iv2-1))/2)+iv1                             6d8s21
                 irec=iv1+nvirt(isbv1)*iv2                              6d8s21
                 ivv=itri+isw*(irec-itri)                               6d8s21
                 lghere=kghere+nn*ivv                                   6d8s21
                 ltmp=jtmp+i1-i10                                       6d11s21
                 do i=0,nnm                                             6d8s21
                  bc(lghere+i)=bc(ltmp+i*nhere)                         6d11s21
                 end do                                                 6d8s21
                end do                                                  6d8s21
                jtmp=jtmp+i1n+1-i10                                     6d11s21
                i10=1                                                   6d8s21
               end do                                                   6d8s21
              end if                                                     6d8s21
              ibcoff=itmp
             end if                                                      6d8s21
             jghere=jghere+nroot*nvx                                    6d8s21
            end if                                                       6d8s21
           end do                                                        6d8s21
          else                                                          6d8s21
           jghere=jghere+nvall*nroot                                    6d8s21
          end if                                                        6d8s21
         end do                                                         6d7s21
         if(ifirsttime.eq.0)call dws_synca                              6d8s21
         ifirsttime=1                                                   6d8s21
         i2o=1                                                          6d8s21
         i2c=nvall*nroot                                                6d8s21
         i0o=1                                                          6d8s21
         i0c=nff2(nclo2p,isb,1)                                         7d22s21
         if(igflag.ne.0)then                                            3d22s21
          nrow=i2c
          nacc=0                                                        6d17s21
          call ddi_ acc(bc,ibc,nhand(nclo2p,isb),i2o,i2c,i0o,i0c,       12d13s22
     $         bc(ighere),ibc(iacc),nacc)                               12d13s22
          last8(1)=nhand(nclo2p,isb)                                    8d21s21
          last8(2)=i0c                                                  6d17s21
         else                                                           6d17s21
          nacc=0                                                        6d17s21
         end if                                                         3d22s21
         ibcoff=idenh                                                   3d22s21
        end if                                                          12d7s20
       end do                                                           11d25s20
      end do                                                            11d25s20
      call ddi_done(ibc(iacc),nacc)                                     3d22s21
      ibcoff=ibcoffo                                                    6d11s21
      return
      end
