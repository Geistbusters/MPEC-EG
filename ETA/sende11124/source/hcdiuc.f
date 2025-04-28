c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcdiuc(nhand,nff2,iff2,nff0,iff0,vec,ncsfv,ncsf,ncsf2, 6d7s21
     $     nec,mdon,mdoo,nsymb,multh,ixw1,ixw2,kmats,nvirt,             6d7s21
     $     izero,ldebug,maxbxd,srh,nwiacc,bc,ibc)                       11d10s22
      implicit real*8 (a-h,o-z)
      integer*1 nab1(2),nab2(2)
      integer*8 nhand(mdoo+1,nsymb,2),i1,in,i2c,i2o,i0c,i0o,            6d8s21
     $     itesta,itestb,last8(2)                                       6d17s21
      logical ldebug                                                    1d25s21
      dimension nff2(mdoo+1,nsymb),iff2(*),nff0(mdoo+1,3),              6d8s21
     $     iff0(*),vec(ncsfv,*),ncsf(*),multh(8,8),nab4(2,3),idorb1(32),11d25s20
     $     isorb1(32),idorb(32),isorb(32),itest(32,4),kmats(*),         6d7s21
     $     nvirt(*),iden(8,8),nden(8,8),ncsf2(4,*),mden(8,8),itmpx(4,8),6d11s21
     $     jtmpx(4,8),ntmpx(4,8),mtmpx(4)                               6d11s21
      include "common.store"                                            11d25s20
      include "common.mrci"                                             11d25s20
      common/kmfind/invk1(2,8,8,8,2)
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data loopx/2600000/
      if(ldebug)write(6,*)('hi, I am the new and improved hcdiuc!')       1d25s21
      nrootm=nroot-1                                                    6d8s21
      last8(1)=-1                                                       6d17s21
      loop=0
      ibcoffo=ibcoff                                                    12d8s20
      iacc=ibcoff                                                       6d17s21
      ighere=iacc+mynprocg                                              3d22s21
      ibcoff=ighere+maxbxd                                              6d17s21
      call enough('hcdiuc.  1',bc,ibc)
      nacc=0                                                            3d22s21
      if(ldebug)write(6,*)('izero in hcdiuc '),izero                      1d25s21
      i1=1                                                              11d25s20
      if(izero.ne.0)then                                                1d24s21
       do isb=1,nsymb                                                    11d26s20
        do ii=mdon+1,mdoo+1                                              11d26s20
         if(nhand(ii,isb,2).ge.0)then                                   6d8s21
          call ddi_zero(bc,ibc,nhand(ii,isb,2))                         11d15s22
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
       do nclo2p=mdon+1,mdoo+1                                           11d25s20
        if(nff2(nclo2p,isb).gt.0)then                                    6d8s21
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
            nvall=nvall+ncsf2(1,iarg)*nvirt(isbv1)                      6d8s21
            nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       6d8s21
           else                                                         6d8s21
            nvv=nvirt(isbv1)*nvirt(isbv2)                               6d8s21
           end if                                                       6d8s21
           nvall=nvall+nvv*ncsf(iarg)                                   6d8s21
          end if                                                        6d8s21
         end do                                                         6d8s21
         nwds=nvall*nroot*nff2(nclo2p,isb)                              6d8s21
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
         call enough('hcdiuc.  2',bc,ibc)
         ibctop=ibcoff-1                                                11d26s20
         jghere=ighere                                                  6d8s21
         do if2=1,nff2(nclo2p,isb)                                      6d8s21
          do i=idenh,ibctop                                             11d26s20
           bc(i)=0d0                                                    11d26s20
          end do                                                        11d26s20
          i2c=iff2(jff2)                                                11d25s20
          jff2=jff2+1                                                   11d25s20
          i2o=iff2(jff2)                                                11d25s20
          i2o=ibset(i2o,norbx)                                          6d7s21
          i2o=ibset(i2o,norbxx)                                         6d7s21
          jff2=jff2+1                                                   11d25s20
          do nclop=max(mdon+1,nclo2p-2),min(mdoo+1,ntop)                6d7s21
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
           call enough('hcdiuc.  3',bc,ibc)
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
              nvx=nvirt(isbv1)*ncsf2(1,iarg)                            6d8s21
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
              call enough('hcdiuc.  4',bc,ibc)
              do i=itmp,ibcoff-1                                         6d8s21
               bc(i)=0d0                                                 6d8s21
              end do                                                     6d8s21
              nhit=0                                                     6d8s21
              do isa=1,nsymb                                              6d8s21
               isc=multh(isa,isbv12)                                      6d8s21
               if(mden(isa,isc).gt.0)then                                6d8s21
                nhit=1                                                   6d8s21
                i2eu=invk1(1,isa,isc,isbv1,1)                            6d8s21
                icase=invk1(2,isa,isc,isbv1,1)                            6d8s21
                itmpi=ibcoff                                             6d8s21
                ibcoff=itmpi+nhere*mden(isa,isc)                         6d8s21
                call enough('hcdiuc.  5',bc,ibc)
                if(icase.eq.1)then                                       6d8s21
                 do i=0,mden(isa,isc)-1                                  6d8s21
                  icoli=kmats(i2eu)+nhere*ibc(nden(isa,isc)+i)           6d8s21
                  icolt=itmpi+nhere*i                                    6d8s21
                  do j=0,nhere-1                                         6d8s21
                   bc(icolt+j)=bc(icoli+j)                               6d8s21
                  end do                                                 6d8s21
                 end do                                                  6d8s21
                else                                                     6d8s21
                 write(6,*)('icase = 2 !!!')
                end if                                                   6d8s21
                call dgemm('n','n',nhere,nn,mden(isa,isc),1d0,           6d8s21
     $              bc(itmpi),nhere,bc(iden(isa,isc)),mden(isa,isc),1d0,6d8s21
     $              bc(itmp),nhere,                                     6d8s21
     d' hcdiuc.  1')
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
                  lghere=kghere+nroot*ncsf2(1,iarg)*iv                  6d8s21
                  icol=icol-il                                          6d8s21
                  jtmp=itmp+icol                                        6d11s21
                  do i=0,nroot*ncsf2(1,iarg)-1                          6d8s21
                   bc(lghere+i)=bc(jtmp+nhere*i)*srh                    6d9s21
                  end do                                                6d8s21
                 end if                                                 6d8s21
                end do                                                  6d8s21
                kghere=kghere+nvirt(isbv1)*nroot*ncsf2(1,iarg)          6d8s21
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
         i0c=nff2(nclo2p,isb)                                           6d8s21
         if(igflag.ne.0)then                                            3d22s21
          nrow=i2c
          nacc=0                                                        6d17s21
          call ddi_iacc(bc,ibc,nhand(nclo2p,isb,2),i2o,i2c,i0o,i0c,     11d15s22
     $         bc(ighere),ibc(iacc),nacc)                               11d15s22
          nwiacc=nwiacc+i0c*i2c                                         8d11s22
          last8(1)=nhand(nclo2p,isb,2)                                  6d17s21
          last8(2)=i0c                                                  6d17s21
         else                                                           6d17s21
          nacc=0                                                        6d17s21
         end if                                                         3d22s21
         ibcoff=idenh                                                   3d22s21
        end if                                                          12d7s20
       end do                                                           11d25s20
      end do                                                            11d25s20
      call ddi_done(ibc(iacc),nacc)                                     3d22s21
      ibcoff=ibcoffo                                                    3d22s21
      if(ldebug)then                                                    1d25s21
       write(6,*)('all done in hcdiuc')
       call dws_synca                                                    1d24s21
       write(6,*)('what we have for hcdi:')
       do isb=1,nsymb                                                    12d7s20
        write(6,*)('for isb = '),isb
        isbv12=multh(isb,isymmrci)                                         12d7s20
        do l=1,4                                                        6d11s21
         mtmpx(l)=0                                                     6d11s21
        end do                                                          6d11s21
        do ii=mdon+1,mdoo+1                                             6d11s21
         if(nff2(ii,isb).gt.0)then                                      6d11s21
          iarg=ii-mdon                                                  6d11s21
          write(6,*)ii,(ncsf2(l,iarg),l=1,4),nff2(ii,isb)
          do l=1,4                                                      6d11s21
           mtmpx(l)=mtmpx(l)+ncsf2(l,iarg)*nff2(ii,isb)                 6d11s21
          end do                                                        6d11s21
         end if                                                         6d11s21
        end do                                                          6d11s21
        write(6,*)('mtmpx: '),mtmpx
        nvisv=0                                                         6d8s21
        nvnotv=0                                                        6d8s21
        ibcsrt=ibcoff                                                   6d11s21
        nvisv=0                                                         6d12s21
        do isbv1=1,nsymb                                                6d8s21
         isbv2=multh(isbv1,isbv12)                                      6d8s21
         if(isbv2.ge.isbv1)then                                         6d8s21
          nvisv0=0                                                      6d11s21
          if(isbv1.eq.isbv2)then                                        6d8s21
           nvisv=nvisv+nvirt(isbv1)                                     6d8s21
           nvisv0=nvirt(isbv1)                                          6d11s21
           write(6,*)('nvisv '),nvisv,nvirt(isbv1),isbv1
           nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        6d9s21
          else                                                          6d8s21
           nvv=nvirt(isbv1)*nvirt(isbv2)                                6d8s21
          end if                                                        6d8s21
          do l=1,4                                                      6d11s21
           if(l.eq.1)then                                               6d11s21
            nuse=nvv+nvisv0                                             6d11s21
           else                                                         6d11s21
            nuse=nvv                                                    6d11s21
           end if                                                       6d11s21
           ntmpx(l,isbv1)=nuse*mtmpx(l)                                 6d11s21
           itmpx(l,isbv1)=ibcoff                                        6d11s21
           jtmpx(l,isbv1)=ibcoff                                        6d11s21
           ibcoff=itmpx(l,isbv1)+nroot*ntmpx(l,isbv1)                   6d11s21
          end do                                                        6d11s21
          nvnotv=nvnotv+nvv                                             6d8s21
         end if                                                         6d8s21
        end do                                                          6d8s21
        call enough('hcdiuc.  6',bc,ibc)
        xnan=-2d0                                                       6d11s21
        do i=ibcsrt,ibcoff-1                                            6d11s21
         bc(i)=xnan                                                     6d11s21
        end do                                                          6d11s21
        write(6,*)('nvisv,nvnotv '),nvisv,nvnotv,xnan
        do ii=1,mdoo+1                                                   12d7s20
         if(nff2(ii,isb).gt.0)then                                      6d8s21
          iarg=ii-mdon                                                    12d7s20
          nrow=nroot*(nvisv*ncsf2(1,iarg)+nvnotv*ncsf(iarg))
          itmp=ibcoff                                                   6d8s21
          ibcoff=itmp+nrow*nff2(ii,isb)                                 6d8s21
          write(6,*)('for nclop = '),ii,nrow,nff2(ii,isb)                                 6d8s21
          call enough('hcdiuc.  7',bc,ibc)
          do i=itmp,ibcoff-1
           bc(i)=xnan
          end do
          i2o=1                                                         6d8s21
          i2c=nrow                                                      6d8s21
          i0o=1                                                         6d8s21
          i0c=nff2(ii,isb)                                              6d8s21
          call ddi_get(bc,ibc,nhand(ii,isb,2),i2o,i2c,i0o,i0c,bc(itmp)) 11d15s22
          do i=0,nff2(ii,isb)
           do j=0,nrow-1
            ji=itmp+j+nrow*i
            if(abs(bc(ji)).gt.1d20)write(6,*)('getting '),bc(ji),
     $           j,i,nhand(ii,isb,2),ii,isb
           end do
          end do
          jtmp=itmp                                                     6d11s21
          do if2=0,nff2(ii,isb)-1                                       6d11s21
           do isbv1=1,nsymb                                              6d11s21
            isbv2=multh(isbv1,isbv12)                                    6d11s21
            if(isbv2.ge.isbv1)then                                       6d11s21
             if(isbv1.eq.isbv2)then                                      6d11s21
              nvvs=(nvirt(isbv1)*(nvirt(isbv1)+1))/2                    6d11s21
              nvvt=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                    6d11s21
              do iv=0,nvirt(isbv1)-1                                    6d11s21
               ivv=((iv*(iv+1))/2)+iv                                   6d11s21
               do i=0,ncsf2(1,iarg)-1                                    6d11s21
                do ir=0,nroot-1                                         6d11s21
                 iad=jtmpx(1,isbv1)+i+mtmpx(1)*(ivv+nvvs*ir)            6d11s21
                 bc(iad)=bc(jtmp+ir)                                    6d11s21
                 if(isb.eq.4.and.isbv1.eq.2.and.
     $                iad.eq.itmpx(1,isbv1)+735)then
                  write(6,*)('getting it '),bc(iad),('from '),
     $                 jtmp+ir-itmp,ii
                 end if
                end do                                                  6d11s21
                jtmp=jtmp+nroot                                         6d11s21
               end do                                                   6d11s21
              end do                                                    6d11s21
              isw=0                                                     6d11s21
             else                                                       6d11s21
              nvvs=nvirt(isbv1)*nvirt(isbv2)                            6d11s21
              nvvt=nvvs                                                 6d11s21
              isw=1                                                     6d11s21
             end if                                                     6d11s21
             do iv2=0,nvirt(isbv2)-1                                    6d11s21
              itop=iv2-1                                                6d11s21
              itop=itop+isw*(nvirt(isbv1)-1-itop)                       6d11s21
              do iv1=0,itop                                             6d11s21
               nvv=nvvs                                                   6d11s21
               ip=+1                                                      6d11s21
               do l=1,4                                                   6d11s21
                if(ntmpx(l,isbv1).gt.0)then                               6d11s21
                itri=((iv2*(iv2+ip))/2)+iv1                              6d11s21
                irec=iv1+nvirt(isbv1)*iv2                                6d11s21
                ivv=itri+isw*(irec-itri)                                 6d11s21
                 do i=0,ncsf2(l,iarg)-1                                   6d11s21
                  do ir=0,nroot-1                                       6d11s21
                   iad=jtmpx(l,isbv1)+i+mtmpx(l)*(ivv+nvv*ir)           6d11s21
                   bc(iad)=bc(jtmp+ir)                                  6d11s21
                 if(isb.eq.4.and.isbv1.eq.2.and.
     $                iad.eq.itmpx(1,isbv1)+735)then
                  write(6,*)('getting it '),bc(iad),('from '),
     $                 jtmp+ir-itmp,ii
                 end if
                  end do                                                6d11s21
                  jtmp=jtmp+nroot                                       6d11s21
                 end do                                                 6d11s21
                end if                                                   6d11s21
                nvv=nvvt                                                  6d11s21
                ip=-1                                                     6d11s21
               end do                                                   6d11s21
              end do                                                    6d11s21
             end do                                                     6d11s21
             do l=1,4
              jtmpx(l,isbv1)=jtmpx(l,isbv1)+ncsf2(l,iarg)               6d11s21
             end do
            end if                                                       6d11s21
           end do                                                       6d11s21
          end do                                                        6d11s21
          ibcoff=itmp                                                   6d8s21
         end if                                                         6d8s21
        end do
          do isbv1=1,nsymb                                              6d11s21
           isbv2=multh(isbv1,isbv12)                                    6d11s21
           if(isbv2.ge.isbv1)then                                       6d11s21
            if(isbv1.eq.isbv2)then                                      6d11s21
             nvvs=(nvirt(isbv1)*(nvirt(isbv1)+1))/2                     6d11s21
             nvvt=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                     6d11s21
            else                                                        6d11s21
             nvvs=nvirt(isbv1)*nvirt(isbv2)                             6d11s21
             nvvt=nvvs                                                  6d11s21
            end if                                                      6d11s21
            nvv=nvvs                                                    6d11s21
            do l=1,4                                                      6d11s21
             if(ntmpx(l,isbv1).gt.0)then                                6d11s21
              if(l.eq.1)then                                               6d11s21
               write(6,*)('for singlet coupled pairs of symmetry '),
     $            isb,isbv1,isbv2
              else
               write(6,*)('for triplet couple pairs of spin '),l-1,
     $              ('symmetry '),isb,isbv1,isbv2
              end if
             end if                                                     6d11s21
             nvv=nvvt                                                   6d11s21
            end do                                                      6d11s21
           end if                                                       6d11s21
          end do                                                        6d11s21
       end do                                                            12d7s20
      end if                                                            1d25s21
      ibcoff=ibcoffo                                                    6d11s21
      return
      end
