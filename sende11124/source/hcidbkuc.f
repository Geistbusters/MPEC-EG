c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcidbkuc(nhand,nff2,iff2,nff0,iff0,gb,ncsfv,ncsf,ncsf2,7d22s21
     $     nec,mdon,mdoop,nsymb,multh,ixw1,ixw2,isymop,n2e,i2eop,phase2, 7d22s21
     $     nvirt,ldebug,maxbxd,nrootu,srh,isymmrci,ism,irel,idoubo,     7d22s21
     $     irefo,nbasdws,ixmt,norb,bc,ibc)                              11d10s22
      implicit real*8 (a-h,o-z)
      integer*1 nab1(2),nab2(2)
      integer*8 nhand(mdoop,nsymb,2),i1,in,i2c,i2o,i0c,i0o,             7d22s21
     $     itesta,itestb                                                7d22s21
      logical ldebug                                                    1d25s21
      dimension nff2(mdoop,nsymb),iff2(*),nff0(mdoop,3),                7d22s21
     $     iff0(*),ncsf(*),multh(8,8),nab4(2,3),idorb1(32),             6d19s21
     $     isorb1(32),idorb(32),isorb(32),itest(32,4),isymop(*),         6d7s21
     $     nvirt(*),iden(8,8),nden(8,8),ncsf2(*),mden(8,8),itmpx(4,8),  7d22s21
     $     jtmpx(4,8),ntmpx(4,8),mtmpx(4),i2eop(2,*),ism(*),irel(*),    7d22s21
     $     idoubo(*),irefo(*),nbasdws(*),ixmt(8,*),gb(ncsfv,*)          7d22s21
      include "common.store"                                            11d25s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data loopx/2600000/
      if(ldebug)write(6,*)('hi, I am the new and improved hciducbk!')       1d25s21
      nrootm=nrootu-1                                                    6d8s21
      loop=0
      ibcoffo=ibcoff                                                    12d8s20
      igi=ibcoff                                                        6d19s21
      ivhere=igi+nrootu*ncsfv                                           6d19s21
      ibcoff=ivhere+maxbxd                                              6d18s21
      call enough('hcidbkuc.  1',bc,ibc)
      do i=igi,ivhere-1                                                 6d19s21
       bc(i)=0d0                                                        6d19s21
      end do                                                            6d19s21
      nacc=0                                                            3d22s21
      i1=1                                                              11d25s20
      norbx=norb+1                                                      11d25s20
      norbxx=norbx+1                                                    6d7s21
      jff2=1                                                            6d8s21
      do isb=1,nsymb                                                    11d25s20
       isbv12=multh(isb,isymmrci)                                       6d7s21
       do nclo2p=mdon+1,mdoop                                           7d22s21
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
            nvall=nvall+ncsf2(iarg)*nvirt(isbv1)                        7d22s21
            nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       6d8s21
           else                                                         6d8s21
            nvv=nvirt(isbv1)*nvirt(isbv2)                               6d8s21
           end if                                                       6d8s21
           nvall=nvall+nvv*ncsf(iarg)                                   6d8s21
          end if                                                        6d8s21
         end do                                                         6d8s21
         nwds=nvall*nrootu*nff2(nclo2p,isb)                              6d8s21
         i2o=1                                                          6d18s21
         i2c=nvall*nrootu                                                6d18s21
         i0o=1                                                          6d18s21
         i0c=nff2(nclo2p,isb)                                           6d18s21
         call ddi_ get(bc,ibc,nhand(nclo2p,isb,1),i2o,i2c,i0o,i0c,      12d13s22
     $        bc(ivhere))                                               12d13s22
         sz=0d0
         do iz=0,nwds-1
          sz=sz+bc(ivhere+iz)**2
         end do
         sz=sqrt(sz/dfloat(nwds))
c
c     ivd is root, ncsf,nvv,iff2.
c
         igflag=0                                                       3d22s21
         nn=ncsf(iarg)*nrootu                                            11d26s20
         ncsfm=ncsf(iarg)-1                                             6d8s21
         nnm=nn-1                                                       11d26s20
         jvhere=ivhere                                                  6d18s21
         do if2=1,nff2(nclo2p,isb)                                      6d8s21
          i2c=iff2(jff2)                                                11d25s20
          ntest=popcnt(i2c)
          if(ntest.ne.nclo2)then
           call dws_synca
           call dcbit(i2c,norb,'i2c')
           write(6,*)('vs. nclo2 '),nclo2
           call dws_finalize
           stop
          end if
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
           jgs=nff0(nclop,3)+igi-1                                      6d19s21
           do jf=1,nff0(nclop,1)                                        11d26s20
            i0c=iff0(jff0)                                              11d25s20
            ntest=popcnt(i0c)
            if(ntest.ne.nclo)then
             call dws_synca
             call dcbit(i0c,norb,'i0c')
             write(6,*)('vs. nclo '),nclo
             call dws_finalize
             stop
            end if
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
               mcolx=jga+irefo(jsa)*jgb                                 6d18s21
               call genmat(ncsf(iarg),ncsf(karg),ncsf(jarg),ncsfmid1,   6d18s21
     $              iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod,bc,ibc)      11d10s22
               ittt=ibcoff                                              6d18s21
               ibcoff=ittt+ncsf(iarg)*nrootu                             6d18s21
               call enough('hcidbkuc.  2',bc,ibc)
               do i=ittt,ibcoff-1                                       6d18s21
                bc(i)=0d0                                               6d18s21
               end do                                                   6d18s21
               sv=0d0
               sx=0d0
               kvhere=jvhere                                            6d18s21
               do isbv1=1,nsymb                                              6d8s21
                isbv2=multh(isbv1,isbv12)                                    6d8s21
                if(isbv2.ge.isbv1)then                                       6d8s21
                 call ilimts(nvirt(isbv1),nvirt(isbv2),mynprocg,
     $                mynowprog,il,ih,i1s,i1e,i2s,i2e)                  6d18s21
                 if(isbv1.eq.isbv2)then                                      6d8s21
                  nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                      6d8s21
                  nvx=nvirt(isbv1)*ncsf2(iarg)                          7d22s21
                 else                                                        6d8s21
                  nvv=nvirt(isbv1)*nvirt(isbv2)                              6d8s21
                  nvx=0                                                     6d8s21
                 end if                                                      6d8s21
                 nvx0=nvx                                               6d21s21
                 nvx=nvx+nvv*ncsf(iarg)                                     6d8s21
                 nhere=ih+1-il                                               6d8s21
                 if(nhere.gt.0)then                                          6d8s21
                  nall=nhere*nn                                              6d8s21
                  itmp=ibcoff                                                6d8s21
                  ixint=itmp+nall                                           6d8s21
                  ibcoff=ixint+nhere                                    7d22s21
                  call enough('hcidbkuc.  3',bc,ibc)
                  do i=itmp,ibcoff-1                                         6d8s21
                   bc(i)=0d0                                                 6d8s21
                  end do                                                     6d8s21
                  nhit=0                                                     6d8s21
                  jgap=jga+idoubo(jsa)                                  7d22s21
                  jgbp=jgb+idoubo(jsb)                                  7d22s21
                  do ii2e=1,n2e                                         7d22s21
                   if(multh(jsb,isbv1).eq.isymop(i2eop(1,ii2e)))then    7d22s21
                    jxint=ixint                                         7d22s21
                    i10=i1s                                             7d22s21
                    i1n=nvirt(isbv1)                                    7d22s21
                    do i2=i2s,i2e                                       7d22s21
                     i2p=i2+idoubo(isbv2)+irefo(isbv2)-1                7d22s21
                     iav=ixmt(jsa,i2eop(2,ii2e))+jgap+nbasdws(jsa)*i2p  7d22s21
                     fav=phase2*bc(iav)                                 7d22s21
                     if(i2.eq.i2e)i1n=i1e                               7d22s21
                     do i1=i10,i1n                                      7d22s21
                      i1p=i1+idoubo(isbv1)+irefo(isbv1)-1               7d22s21
                      ibv=ixmt(jsb,i2eop(1,ii2e))+jgbp+nbasdws(jsb)*i1p 7d23s21
                      bc(jxint)=bc(jxint)+fav*bc(ibv)                   7d22s21
                      jxint=jxint+1                                     7d22s21
                     end do                                             7d22s21
                     i10=1                                              7d22s21
                    end do                                              7d22s21
                   end if                                               7d22s21
                  end do                                                7d22s21
                  isw=1                                                    6d8s21
                  if(isbv1.eq.isbv2)then                                   6d8s21
                   isw=0                                                   6d8s21
                   do ivp=i2s,i2e                                          6d8s21
                    iv=ivp-1                                               6d8s21
                    icol=ivp+nvirt(isbv1)*iv                               6d8s21
                    if(icol.ge.il.and.icol.le.ih)then                      6d8s21
                     jtmpi=ixint+icol-il                                7d22s21
                     xint=bc(jtmpi)*srh                                 6d18s21
                     lvhere=kvhere+nrootu*ncsf2(iarg)*iv                7d22s21
                     do i=0,nrootu*ncsf2(iarg)-1                        7d22s21
                      bc(ittt+i)=bc(ittt+i)+bc(lvhere+i)*xint           6d18s21
                     end do                                                6d8s21
                    end if                                                 6d8s21
                   end do                                                  6d8s21
                   kvhere=kvhere+nvirt(isbv1)*nrootu*ncsf2(iarg)        7d22s21
                  end if                                                   6d8s21
                  i10=i1s                                                  6d8s21
                  i1n=nvirt(isbv1)                                         6d8s21
                  jtmpi=ixint                                           7d22s21
                  do i2=i2s,i2e                                            6d8s21
                   iv2=i2-1                                                6d8s21
                   if(i2.eq.i2e)i1n=i1e                                    6d8s21
                   itop=min(i1n,iv2+isw*(nvirt(isbv1)-iv2))                6d8s21
                   do i1=i10,itop                                          6d8s21
                    iv1=i1-1                                               6d8s21
                    itri=((iv2*(iv2-1))/2)+iv1                             6d8s21
                    irec=iv1+nvirt(isbv1)*iv2                              6d8s21
                    ivv=itri+isw*(irec-itri)                               6d8s21
                    lvhere=kvhere+nn*ivv                                6d18s21
                    do i=0,nnm                                             6d8s21
                     bc(ittt+i)=bc(ittt+i)+bc(lvhere+i)*bc(jtmpi)       6d18s21
                    end do                                              6d18s21
                    jtmpi=jtmpi+1                                       6d18s21
                   end do                                                 6d8s21
                   if(itop.ge.i10)then                                  6d21s21
                    jtmpi=jtmpi+i1n-itop                                 6d21s21
                   else                                                 6d21s21
                    jtmpi=jtmpi+i1n+1-i10                                6d21s21
                   end if                                               6d21s21
                   i10=1                                                   6d8s21
                  end do                                                   6d8s21
                  kvhere=kvhere+nvv*nrootu*ncsf(iarg)                    6d18s21
                  ibcoff=itmp                                           7d22s21
                 else                                                   6d18s21
                  kvhere=kvhere+nrootu*nvx                               6d18s21
                 end if                                                     6d8s21
                end if                                                  6d18s21
               end do                                                   6d18s21
               itmp=ibcoff                                              6d18s21
               ibcoff=itmp+ncsf(jarg)*nrootu                             6d18s21
               call dgemm('n','n',nrootu,ncsf(jarg),ncsf(iarg),1d0,      6d18s21
     $              bc(ittt),nrootu,bc(iprod),ncsf(iarg),0d0,            6d18s21
     $              bc(itmp),nrootu,                                     6d18s21
     d' hcidbkuc.  1')
               jtmp=itmp-1                                              6d18s21
               do j=0,ncsf(jarg)-1                                      6d18s21
                do ir=1,nrootu                                           6d18s21
                 iad=jgs+j+ncsfv*(ir-1)                                 6d19s21
                 bc(iad)=bc(iad)+bc(jtmp+ir)                            6d19s21
                end do                                                  6d18s21
                jtmp=jtmp+nrootu                                         6d18s21
               end do                                                   6d18s21
               ibcoff=iprod                                             6d18s21
               if(ipss.eq.2)go to 1                                      1d25s21
              end if                                                     1d25s21
             end do                                                      1d25s21
    1        continue                                                    1d25s21
            end if                                                       11d25s20
            jgs=jgs+ncsf(jarg)                                           11d26s20
           end do                                                        11d25s20
          end do                                                         11d25s20
          jvhere=jvhere+nrootu*nvall                                     6d18s21
         end do                                                         6d7s21
        end if                                                          12d7s20
       end do                                                           11d25s20
      end do                                                            11d25s20
      i1=1                                                              6d19s21
      in=ncsfv*nrootu                                                   6d19s21
      jgi=igi-1                                                         7d22s21
      sz=0d0
      do ir=1,nrootu                                                    7d22s21
       do k=1,ncsfv                                                     7d22s21
        gb(k,ir)=gb(k,ir)+bc(jgi+k)                                     7d22s21
       end do                                                           7d22s21
       jgi=jgi+ncsfv                                                    7d22s21
      end do                                                            7d22s21
      ibcoff=ibcoffo                                                    6d11s21
      return
      end
