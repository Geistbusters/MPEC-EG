c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcdduc(ihddiag,nff2,iff2,ncsf,ncsf2,nec,mdon,mdoo,     6d22s21
     $     nsymb,multh,ixw1,ixw2,                                       7d6s21
     $     ih0av,nh0av,ioooo,jmats,kmats,i4x,i4xb,nvirt,nrootu,ism,     7d6s21
     $    irel,irefo,isymmrci,norb,lprt,maxbxd,sr2,srh,timex,tovr,      8d11s22
     $     nwiacc,bc,ibc)                                               11d10s22
      implicit real*8 (a-h,o-z)
      external second                                                   2d18s21
      integer*8 ihddiag(mdoo+1,nsymb,2),i18,i28,i38,i48,i2c,i2o,j2c,j2o,12d14s20
     $     itestc,itesto,last8(2),k2o,l2o                               6d28s21
      integer*1 nab1(2),nab2(2),nab1b(2),nab2b(2)                       2d19s21
      logical lkeep,lprt,lchoice                                        3d17s21
      dimension nff2(mdoo+1,nsymb,2),ncsf(*),multh(8,8),ih0av(*),       12d12s20
     $     nh0av(*),jmats(*),kmats(*),nvirt(*),ism(*),irel(*),
     $     itest(32,3),isorb2(32),idorb2(32),jsorb2(32),jdorb2(32),     6d22s21
     $     nab4(2,3),idenj(8),idenk(8,2),irefo(*),iff2(*),ndenj(8),     7d1s21
     $     ndenk(8),ioooo(*),idenjk(8),ndenjk(8),ncsf2(4,*),            6d23s21
     $     ivec(7,8),ivecp(7,8),nvecp(8),ircv(7,8),nrcv(7,8),           6d22s21
     $     mdenjk(8),idenjkc(8),mdenj(8),idenjc(8),mdenk(8),            6d23s21
     $     idenkc(8),ww(5),i4xb(*),                                     7d6s21
     $     itransv(7,8),i4x(*),mtmpx(4),ntmpx(4,8),itmpx(4,8),jtmpx(4,8)
      include "common.store"                                            12d12s20
      common/fnd2cm/inv(2,8,8,8)                                        4d9s18
      common/kmfind/invk1(2,8,8,8,2)                                    6d30s18
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data loopx/660000/
      do isb=1,nsymb                                                    11d15s21
       do i=1,7                                                         11d15s21
        nrcv(i,isb)=0                                                   11d15s21
       end do                                                           11d15s21
      end do                                                            11d15s21
      if(lprt)                                                          1d25s21
     $write(6,*)('hello    ! my name is hcdduc'),isymmrci,nec,mdon,mdoo,
     $ nsymb,ixw1,ixw2,norb
      iacc=ibcoff                                                       2d12s21
      ibcoff=iacc+mynprocg                                              2d12s21
      do isb=1,nsymb                                                    2d12s21
       nvecp(isb)=0                                                     2d12s21
       do i=1,7                                                         6d25s21
        ivec(i,isb)=ibcoff                                              2d12s21
        ircv(i,isb)=ivec(i,isb)+maxbxd                                  6d21s21
        ibcoff=ircv(i,isb)+mynprocg                                     2d12s21
       end do                                                           2d12s21
      end do                                                            2d12s21
      last8(1)=-1                                                       2d8s21
      igg=ibcoff                                                        2d12s21
      ibcoff=igg+maxbxd                                                 6d21s21
      iygoal=1
      nacc=0                                                            1d30s21
      itransgg=0                                                        2d12s21
      call enough('hcdduc.  1',bc,ibc)
      nrootm=nrootu-1                                                   12d15s20
      nsing=0                                                           12d15s20
      loop=0
      norbx=norb+1                                                      12d14s20
      norbxx=norbx+1                                                    12d14s20
      norbxxx=norbxx+1                                                  6d21s21
      idoit=0                                                           6d21s21
      do isb=1,nsymb                                                    12d12s20
       isbv12=multh(isb,isymmrci)                                       6d21s21
       mvisv=0                                                          6d21s21
       mvnotv=0                                                         6d21s21
       do isbv1=1,nsymb                                                 6d21s21
        isbv2=multh(isbv1,isbv12)                                       6d21s21
        if(isbv2.ge.isbv1)then                                          6d21s21
         if(isbv1.eq.isbv2)then                                         6d21s21
          mvisv=mvisv+nvirt(isbv1)                                      6d21s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         6d21s21
         else                                                           6d21s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 6d21s21
         end if                                                         6d21s21
         mvnotv=mvnotv+nvv                                              6d21s21
        end if                                                          6d21s21
       end do                                                           6d21s21
       do jsb=1,nsymb                                                   2d12s21
        nvecp(jsb)=0                                                    2d12s21
       end do                                                           2d12s21
       do ncloip=mdon+1,mdoo+1                                           12d12s20
        if(min(nff2(ncloip,isb,1),mvisv+mvnotv).gt.0)then               6d21s21
         ncloi=ncloip-1                                                  12d12s20
         iarg=ncloip-mdon                                                12d12s20
         mrow=(ncsf2(1,iarg)*mvisv+ncsf(iarg)*mvnotv)*nrootu            6d21s21
         nopeni=nec-2*ncloi                                              12d12s20
         nopenim=nopeni-1                                               6d21s21
         do jsb=1,nsymb                                                   12d12s20
          ijsb=multh(jsb,isb)                                           6d22s21
          jsbv12=multh(jsb,isymmrci)                                    6d21s21
          nvisv=0                                                       6d21s21
          nvnotv=0                                                      6d21s21
          do jsbv1=1,nsymb                                              6d21s21
           jsbv2=multh(jsbv1,jsbv12)                                    6d21s21
           if(jsbv2.ge.jsbv1)then                                       6d21s21
            if(jsbv1.eq.jsbv2)then                                      6d21s21
             nvisv=nvisv+nvirt(jsbv1)                                   6d21s21
             nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                      6d21s21
            else                                                        6d21s21
             nvv=nvirt(jsbv1)*nvirt(jsbv2)                              6d21s21
            end if                                                      6d21s21
            nvnotv=nvnotv+nvv                                           6d21s21
           end if                                                       6d21s21
          end do                                                        6d21s21
          if(isbv12.eq.1)then                                           6d21s21
           nctop=ncloip+3                                               6d21s21
          else                                                          6d21s21
           nctop=ncloip+2                                               6d21s21
          end if                                                        6d21s21
          nmrow=nrootu*(ncsf2(1,iarg)*nvisv+ncsf(iarg)*nvnotv)          6d30s21
          nmrowk=nrootu*ncsf(iarg)*(nvisv+nvnotv)                       7d1s21
          nmrowk2=nmrowk*2                                              7d2s21
          do nclojp=max(mdon+1,ncloip-2),min(mdoo+1,nctop)              6d21s21
           if(min(nff2(nclojp,jsb,1),nvisv+nvnotv).gt.0)then            6d21s21
            ncloj=nclojp-1                                               12d12s20
            jarg=nclojp-mdon                                             12d12s20
            nopenj=nec-2*ncloj                                           12d14s20
            nrow=(ncsf2(1,jarg)*nvisv+ncsf(jarg)*nvnotv)*nrootu         6d21s21
            ivecbc=ibcoff                                               2d12s21
            call ilimts(1,nff2(nclojp,jsb,1),mynprocg,mynowprog,ilf,    6d22s21
     $           ihf,i1s,i1e,i2s,i2e)                                   6d22s21
            nfhere=ihf+1-ilf                                            6d22s21
c
c     has vector alread been loaded?
c
            nminc=mdoo+1                                                2d12s21
            do i=1,nvecp(jsb)                                           2d12s21
             if(nclojp.eq.ivecp(i,jsb))then                             2d12s21
              nuse=i                                                    2d12s21
              go to 22                                                  2d12s21
             end if                                                     2d12s21
             if(ivecp(i,jsb).lt.nminc)then                              2d12s21
              nminc=ivecp(i,jsb)                                        2d12s21
              iminc=i                                                   2d12s21
             end if                                                     2d12s21
            end do                                                      2d12s21
c
c     no ... are we building up library?
c
            if(nvecp(jsb).lt.7)then                                     6d22s21
c     build
             nuse=nvecp(jsb)+1                                          2d12s21
             nvecp(jsb)=nvecp(jsb)+1                                    2d12s21
             ivecp(nvecp(jsb),jsb)=nclojp                               2d12s21
            else
c
c     smallest one out
c
             nuse=iminc                                                 2d12s21
             ivecp(iminc,jsb)=nclojp                                    2d12s21
            end if
            i18=1                                                       6d21s21
            i28=nrow                                                    6d21s21
            i38=ilf                                                     6d22s21
            i48=ihf                                                     6d22s21
            itransv(nuse,jsb)=0                                         2d12s21
            call ddi_iget(bc,ibc,ihddiag(nclojp,jsb,1),i18,i28,i38,i48, 11d15s22
     $           bc(ivec(nuse,jsb)),ibc(ircv(nuse,jsb)),nrcv(nuse,jsb)) 2d12s21
   22       continue                                                    2d12s21
            if2o=nff2(ncloip,isb,2)                                      12d14s20
            idenstrt=ibcoff                                             6d22s21
            do isa=1,nsymb                                              6d22s21
             ncold=(irefo(isa)*(irefo(isa)+1))/2                         6d22s21
             isc=multh(isa,ijsb)                                        6d24s21
             ncol=irefo(isa)*irefo(isc)                                 6d24s21
             ndenk(isa)=ibcoff                                          6d24s21
             ibcoff=ndenk(isa)+ncol                                     6d24s21b*
             if(isc.ge.isa)then                                         6d24s21
              if(isc.eq.isa)then                                        6d24s21
               ncolu=ncold                                              6d24s21
              else                                                      6d24s21
               ncolu=ncol                                               6d24s21
              end if                                                    6d24s21
              ndenj(isa)=ibcoff                                         6d29s21
              ibcoff=ndenj(isa)+ncolu                                   6d29s21
             end if                                                     6d24s21
            end do                                                      6d22s21
            mdenh=ibcoff                                                6d23s21
            ibcoff=ibcoff+1                                             6d23s21
            nncol=ibcoff-idenstrt                                       6d23s21
            idenh=ibcoff                                                6d22s21
            ibcoff=idenh+nrow                                           6d22s21
            do isa=1,nsymb                                              6d22s21
             ncold=(irefo(isa)*(irefo(isa)+1))/2                         6d22s21
             isc=multh(isa,ijsb)                                        6d23s21
             ncol=irefo(isa)*irefo(isc)                                 6d23s21
             idenk(isa,1)=ibcoff                                          6d23s21
             idenk(isa,2)=idenk(isa,1)+nmrowk*ncol                      7d2s21
             ibcoff=idenk(isa,2)+nmrowk*ncol                            7d2s21
             if(isc.ge.isa)then                                         6d22s21
              if(isc.eq.isa)then                                        6d22s21
               ncolu=ncold                                              6d22s21
              else                                                      6d22s21
               ncolu=ncol                                               6d22s21
              end if                                                    6d22s21
              idenj(isa)=ibcoff                                         6d22s21
              ibcoff=idenj(isa)+nmrow*ncolu                             6d30s21
             end if                                                     6d22s21
            end do                                                      6d22s21
            idenend=ibcoff-1                                            6d22s21
            ndenend=ibcoff-idenstrt                                     6d22s21
            call enough('hcdduc.  2',bc,ibc)
            jgg=igg                                                     12d14s20
            jvecer=0                                                    2d12s21
            do if2=1,nff2(ncloip,isb,1)                                  12d14s20
             call second(timex)
             do i=idenstrt,idenend                                      6d22s21
              bc(i)=0d0                                                 6d22s21
             end do                                                     6d22s21
             i2c=iff2(if2o)                                                11d25s20
             ntest=popcnt(i2c)
             if(ntest.ne.ncloi)then
              write(6,*)('ntest = '),ntest,if2,ncloi,if2o
              call dcbit(i2c,32,'dorb')
              call dws_synca
              call dws_finalize
              stop
             end if
             do i=1,norb                                                   11d25s20
              itest(i,1)=0                                                 11d25s20
             end do                                                        11d25s20
             ii=1                                                          11d25s20
             do i=1,norb                                                   11d25s20
              if(btest(i2c,i))then                                         11d25s20
               idorb2(ii)=i                                                11d25s20
               itest(i,1)=2                                                11d25s20
               ii=ii+1                                                     11d25s20
              end if                                                       11d25s20
             end do                                                        11d25s20
             if2o=if2o+1                                                   11d25s20
             i2o=iff2(if2o)                                                11d25s20
             i2o=ibset(i2o,norbx)                                          11d25s20
             l2o=ibset(i2o,norbxxx)                                     6d28s21
             i2o=ibset(i2o,norbxx)                                      6d21s21
             ii=1                                                          11d25s20
             do i=1,norb                                                6d22s21
              if(btest(i2o,i))then                                         11d25s20
               isorb2(ii)=i                                                11d25s20
               itest(i,1)=1                                                11d25s20
               ii=ii+1                                                     11d25s20
              end if                                                       11d25s20
             end do                                                        11d25s20
             isorb2(ii)=norbx                                           6d21s21
             ii=ii+1                                                    6d21s21
             isorb2(ii)=norbxx                                          6d21s21
             if2o=if2o+1                                                   11d25s20
             jf2o=nff2(nclojp,jsb,2)                                      12d14s20
             jvec=ivec(nuse,jsb)                                        2d12s21
             do jf2=1,nff2(nclojp,jsb,1)                                  12d14s20
              if(jf2.ge.ilf.and.jf2.le.ihf)then                          6d22s21
               jvd=ivec(nuse,jsb)+nrow*(jf2-ilf)                        6d22s21
               j2c=iff2(jf2o)                                                11d25s20
               jf2o=jf2o+1                                                   11d25s20
               j2o=iff2(jf2o)                                                11d25s20
               j2o=ibset(j2o,norbx)                                      6d21s21
               j2o=ibset(j2o,norbxx)                                          11d25s20
               jf2o=jf2o+1                                                   11d25s20
               if(isb.eq.jsb)then                                       6d22s21
                call gandc4(i2c,i2o,j2c,j2o,nopeni,nopenj,norbxx,nnot,    12d14s20
     $             nab4,bc,ibc)                                         11d14s22
c     jvd is root, ncsf,nvv,jff2.
c     re-order it to nvv,nroot,ncsf,iff2                                6d21s21
                if(itransv(nuse,jsb).eq.0.and.nnot.gt.0)then              2d12s21
                 call ddi_done(ibc(ircv(nuse,jsb)),nrcv(nuse,jsb))        2d12s21
                 itmp=ibcoff                                              6d21s21
                 ibcoff=itmp+nrow                                         6d21s21
                 call enough('hcdduc.  3',bc,ibc)
                 szz=0d0
                 do iff=0,nfhere-1                                       6d22s21
                  jjvd=ivec(nuse,jsb)+nrow*iff                          6d26s21
                  jtmp=itmp                                               6d21s21
                  do jsbv1=1,nsymb                                               6d9s21
                   jsbv2=multh(jsbv1,jsbv12)                                     6d10s21
                   if(jsbv2.ge.jsbv1)then                                        6d9s21
                    if(jsbv1.eq.jsbv2)then                                       6d9s21
                     do iv=0,nvirt(jsbv1)-1                                      6d9s21
                      do i=0,ncsf2(1,jarg)-1                                    6d9s21
                       do ir=0,nrootm                                           6d9s21
                        kvd=jtmp+iv+nvirt(jsbv1)*(ir+nrootu*i)                   6d16s21
                        bc(kvd)=bc(jjvd+ir)                             6d26s21
                        szz=szz+bc(kvd)**2
                       end do                                                   6d9s21
                       jjvd=jjvd+nrootu                                 6d26s21
                      end do                                                    6d9s21
                     end do                                                     6d9s21
                     nhere=nvirt(jsbv1)*ncsf2(1,jarg)*nrootu                     6d9s21
                     jtmp=jtmp+nhere                                               6d16s21
                     nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                       6d9s21
                    else                                                         6d9s21
                     nvv=nvirt(jsbv1)*nvirt(jsbv2)                               6d9s21
                    end if                                                       6d9s21
                    do ivv=0,nvv-1                                              6d17s21
                     do i=0,ncsf(jarg)-1                                        6d17s21
                      do ir=0,nrootm                                            6d17s21
                       kvd=jtmp+ivv+nvv*(ir+nrootu*i)                            6d17s21
                        szz=szz+bc(kvd)**2
                       bc(kvd)=bc(jjvd+ir)                              6d26s21
                      end do                                                    6d9s21
                      jjvd=jjvd+nrootu                                  6d26s21
                     end do                                                     6d17s21
                    end do                                                      6d17s21
                    jtmp=jtmp+nvv*nrootu*ncsf(jarg)                               6d17s21
                   end if                                                        6d9s21
                  end do                                                         6d9s21
                  jjvd=ivec(nuse,jsb)+nrow*iff                          6d26s21
                  do i=0,nrow-1                                           6d21s21
                   bc(jjvd+i)=bc(itmp+i)                                6d26s21
                  end do                                                  6d21s21
                 end do                                                         6d16s21
                 ibcoff=itmp                                              6d21s21
                 itransv(nuse,jsb)=1                                      2d12s21
                end if                                                     2d12s21
                if(itransgg.eq.0.and.nnot.gt.0)then                       6d21s21
                 call ddi_done(ibc(iacc),nacc)                             2d12s21
                 itransgg=1                                                2d12s21
                 do i=igg,igg+mrow*nff2(ncloip,isb,1)-1                   6d21s21
                  bc(i)=0d0                                                2d12s21
                  if(i.eq.iygoal)write(6,*)('zeroing '),iygoal
                 end do                                                    2d12s21
                end if                                                     2d12s21
                if(nnot.gt.1.and.nab4(1,1).eq.0)then
                 write(6,*)('after gandc4, nnot = '),nnot
                 write(6,*)('nab4: '),nab4
                 stop
                end if
                if(nnot.eq.1)then                                         6d21s21
c     we already have everything, independent of csf, except for 4v.
                 nn=ncsf(iarg)*ncsf(iarg)                                 6d21s21
                 nnm=nn-1                                                 6d21s21
                 imat4o=ibcoff                                            6d21s21
                 ibcoff=imat4o+nn                                         6d21s21
                 call enough('hcdduc.  4',bc,ibc)
                 do i=imat4o,ibcoff-1                                     6d21s21
                  bc(i)=0d0                                               6d21s21
                 end do                                                   6d21s21
                 do i1=1,nopeni-1                                         6d21s21
                  if(i1.eq.nopenim)then                                   6d21s21
                   nab1(2)=norbx                                          6d21s21
                  else                                                    6d21s21
                   nab1(2)=isorb2(i1)                                     6d21s21
                   lsb=ism(nab1(2))                                              11d13s20
                   lgb=irel(nab1(2))-1                                           11d13s20
                   lsc=lsb                                                6d21s21
                   lgc=lgb                                                6d21s21
                  end if                                                  6d21s21
                  do i2=i1+1,nopeni                                       6d21s21
                   if(i2.ge.nopenim)then                                   6d21s21
                    nab1(1)=norbx+i2-nopenim                              6d21s21
                   else
                    nab1(1)=isorb2(i2)                                    6d21s21
                    lsa=ism(nab1(1))                                              11d13s20
                    lga=irel(nab1(1))-1                                           11d13s20
                    lsd=lsa                                               6d21s21
                    lgd=lga                                               6d21s21
                   end if                                                 6d21s21
                   nopenk=nopeni-2                                              11d13s20
                   karg=iarg+1                                                 11d13s20
                   nab2(1)=nab1(2)                                             11d13s20
                   nab2(2)=nab1(1)                                             11d13s20
                   nqq=karg+mdon-1                                      6d23s21
                   if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                    itestc=i2c                                            12d14s20
                    itesto=ibclr(i2o,nab1(1))                             6d21s21
                    itesto=ibclr(itesto,nab1(2))                          6d21s21
                    itestc=ibset(itestc,nab1(2))                          6d21s21
                    call gandc(itestc,itesto,j2c,j2o,nopenk,nopenj,       6d21s21
     $         karg,iarg,ncsf,norbxx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,   12d14s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
                    call gandc(i2c,i2o,itestc,itesto,nopeni,nopenk,       12d14s20
     $         iarg,karg,ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,   12d14s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
                    call genmat(ncsf(iarg),ncsf(karg),ncsf(jarg),         6d21s21
     $                ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod,  11d10s22
     $                   bc,ibc)                                        11d10s22
                   else                                                   6d21s21
                    iprod=ibcoff                                          6d21s21
                    ibcoff=iprod+ncsf(iarg)*ncsf(iarg)                    6d21s21
                    call enough('hcdduc.  5',bc,ibc)
                    do i=iprod,ibcoff-1                                   6d21s21
                     bc(i)=0d0                                            6d21s21
                    end do                                                6d21s21
                   end if                                                 6d21s21
                   do i=0,ncsf(iarg)-1                                    6d21s21
                    ii=iprod+i*(ncsf(iarg)+1)                             6d21s21
                    bc(ii)=bc(ii)-1d0                                     6d21s21
                   end do                                                 6d21s21
                   if(i1.lt.nopenim.and.i2.lt.nopenim)then                6d21s21
                    xint=getint(ioooo,lsa,lsb,lsc,lsd,lga,lgb,lgc,lgd,  11d15s22
     $                  bc,ibc)                                         11d15s22
                    do i=0,nnm                                           6d21s21
                     bc(imat4o+i)=bc(imat4o+i)+xint*bc(iprod+i)         6d24s21
                    end do                                               6d21s21
                   else if(i1.lt.nopenim.and.i2.eq.nopenim)then         6d23s21
                   end if                                                 6d21s21
                   ibcoff=iprod                                          6d22s21
                  end do                                                  6d21s21
                 end do                                                   6d21s21
                 jgg=igg+mrow*(jf2-1)                                    6d21s21
                 call hcdduc1(bc(imat4o),ncsf,ncsf2,iarg,iarg,bc(jvd),   6d22s21
     $               bc(jgg),jsbv12,nsymb,multh,nvirt,nrootu,sr2,       6d24s21
     $                bc(iygoal),jgg,'b',4,bc,ibc)                      11d10s22
                 ibcoff=imat4o                                           6d22s21
                else if(nnot.eq.2)then                                          12d14s20
                 call gandc(i2c,i2o,j2c,j2o,nopeni,nopenj,iarg,jarg,    6d22s21
     $                ncsf,norbxx,ixw1,ixw2,nnot,nab1,iwpb1,iwpk1,      6d22s21
     $                ncsfmid1,bc,ibc)                                  11d14s22
                 iprod=ibcoff                                            6d22s21
                 ibcoff=iprod+ncsf(iarg)*ncsf(jarg)                      6d22s21
                 call enough('hcdduc.  6',bc,ibc)
                 call prodn(iwpb1,iwpk1,ncsf(iarg),ncsf(jarg),ncsfmid1,  6d22s21
     $               bc(iprod),bc,ibc,1d0,0d0)                          2d13s23
                 if(nab1(1).gt.nab1(2))then                              6d22s21
                  icopy=nab1(1)                                          6d22s21
                  nab1(1)=nab1(2)                                        6d22s21
                  nab1(2)=icopy                                          6d22s21
                 end if                                                  6d22s21
                 isa=ism(nab1(1))                                        6d22s21
                 iga=irel(nab1(1))-1                                     6d22s21
                 igb=irel(nab1(2))-1                                     6d22s21
                 ncolj=((igb*(igb+1))/2)+iga                             6d22s21
                 ih0=ih0av(isa)+iga+nh0av(isa)*igb                       6d22s21
                 sum=bc(ih0)
                 do i=1,norb                                                   11d25s20
                  itest(i,2)=0                                                 11d25s20
                 end do                                                        11d25s20
                 ii=1                                                          11d25s20
                 do i=1,norb                                                   11d25s20
                  if(btest(j2c,i))then                                         11d25s20
                   itest(i,2)=2                                                11d25s20
                   ii=ii+1                                                     11d25s20
                  end if                                                       11d25s20
                 end do                                                        11d25s20
                 ii=1                                                          11d25s20
                 do i=1,norb                                                   11d25s20
                  if(btest(j2o,i))then                                         11d25s20
                   itest(i,2)=1                                                11d25s20
                   ii=ii+1                                                     11d25s20
                  end if                                                       11d25s20
                 end do                                                        11d25s20
                 nok=0                                                         11d13s20
                 do i=1,norb                                                   11d13s20
                  ixn=min(itest(i,1),itest(i,2))
                  if(ixn.gt.0)then                                             11d13s20
                   nok=nok+1                                                   11d13s20
                   itest(nok,3)=ixn                                            11d13s20
                   itest(nok,2)=i                                              11d13s20
                  end if                                                       11d13s20
                 end do                                                        11d13s20
                 do i=1,nok                                                    11d13s20
                  js=ism(itest(i,2))                                           11d13s20
                  jg=irel(itest(i,2))-1                                        11d13s20
                  if(itest(i,3).eq.2)then                                  12d14s20
                   orig=sum                                             6d25s21
                   xj=getint(ioooo,isa,isa,js,js,iga,igb,jg,jg,bc,ibc)  11d15s22
                   xk=getint(ioooo,isa,js,isa,js,iga,jg,igb,jg,bc,ibc)  11d15s22
                   sum=sum+xj*2d0
     $                   -xk
                  else                                                   6d22s21
                   orig=sum                                             6d25s21
                   xj=getint(ioooo,isa,isa,js,js,iga,igb,jg,jg,bc,ibc)  11d15s22
                   sum=sum+xj
                  end if                                                 6d22s21
                 end do                                                    12d14s20
                 do i=0,ncsf(iarg)*ncsf(jarg)-1                          6d22s21
                  bc(iprod+i)=bc(iprod+i)*sum                            6d22s21
                 end do                                                  6d22s21
                 jgg=igg+mrow*(if2-1)                                    6d21s21
                 call hcdduc1(bc(iprod),ncsf,ncsf2,iarg,jarg,bc(jvd),    6d22s21
     $               bc(jgg),jsbv12,nsymb,multh,nvirt,nrootu,sr2,       6d24s21
     $                bc(iygoal),jgg,'d',4,bc,ibc)                      11d10s22
                 ibcoff=iprod                                            6d22s21
                 do i=1,nok                                                    11d13s20
                  if(itest(i,3).eq.1)then                                  12d14s20
                   itestc=i2c                                              12d14s20
                   itesto=i2o                                              12d14s20
                   nopenk=nopeni                                                11d13s20
c
c     anihilate common
c
                   if(btest(itestc,itest(i,2)))then                             11d13s20
                    itestc=ibclr(itestc,itest(i,2))                             11d13s20
                    itesto=ibset(itesto,itest(i,2))                             11d13s20
                    karg=iarg-1                                                11d13s20
                    nopenk=nopenk+1                                             11d13s20
                   else                                                         11d13s20
                    itesto=ibclr(itesto,itest(i,2))                             11d13s20
                    karg=iarg                                                  11d13s20
                    nopenk=nopenk-1                                             11d13s20
                   end if                                                       11d13s20
c
c     create ket
c
                   if(btest(itesto,nab4(2,1)))then                        12d14s20
                    itestc=ibset(itestc,nab4(2,1))                        12d14s20
                    itesto=ibclr(itesto,nab4(2,1))                        12d14s20
                    karg=karg+1                                                11d13s20
                    nopenk=nopenk-1                                             11d13s20
                   else                                                         11d13s20
                    itesto=ibset(itesto,nab4(2,1))                        12d14s20
                    nopenk=nopenk+1                                             11d13s20
                   end if                                                       11d13s20
                   nqq=karg+mdon-1
                   if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                    call gandc(i2c,i2o,itestc,itesto,nopeni,nopenk,       2d19s21
     $            iarg,karg,ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,     6d23s21
     $                iwpk1,ncsfmid1,bc,ibc)                            11d14s22
                    call gandc(itestc,itesto,j2c,j2o,nopenk,nopenj,       2d19s21
     $            karg,jarg,ncsf,norbxx,ixw1,ixw2,nnot2,nab2,iwpb2,     6d23s21
     $               iwpk2,ncsfmid2,bc,ibc)                             11d14s22
                    if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                     call genmat(ncsf(iarg),ncsf(karg),ncsf(jarg),       6d22s21
     $                ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,        6d23s21
     $                iprod,bc,ibc)                                     11d10s22
                     lsa=ism(nab1(1))                                      12d14s20
                     lga=irel(nab1(1))-1                                   12d14s20
                     lsb=ism(nab1(2))                                      12d14s20
                     lgb=irel(nab1(2))-1                                   12d14s20
                     lsc=ism(nab2(1))                                      12d14s20
                     lgc=irel(nab2(1))-1                                   12d14s20
                     lsd=ism(nab2(2))                                      12d14s20
                     lgd=irel(nab2(2))-1                                   12d14s20
                     xint=getint(ioooo,lsa,lsb,lsc,lsd,lga,lgb,lgc,lgd, 11d15s22
     $                    bc,ibc)                                       11d15s22
                     do ii=0,ncsf(iarg)*ncsf(jarg)-1                      6d22s21
                      bc(iprod+ii)=bc(iprod+ii)*xint                    6d22s21
                     end do                                              6d22s21
                     call hcdduc1(bc(iprod),ncsf,ncsf2,iarg,jarg,       6d22s21
     $                    bc(jvd),bc(jgg),jsbv12,nsymb,multh,nvirt,     6d24s21
     $                    nrootu,sr2,bc(iygoal),jgg,'e',4,bc,ibc)       11d10s22
                     ibcoff=iprod                                        6d22s21
                    end if                                                  12d14s20
                   end if                                                 4d14s21
                  end if                                                   12d14s20
                 end do                                                    12d14s20
                else if(nnot.ge.3)then                                     12d14s20
                 if(nnot.eq.3)then                                          12d8s20
                  ipssx=1                                                   12d8s20
                 else                                                       12d8s20
                  ipssx=2                                                   12d8s20
                 end if                                                     12d8s20
                 jgg=igg+mrow*(if2-1)                                    6d21s21
                 do ipss=1,ipssx                                            12d8s20
                  if(ipss.eq.1)then
                   iu1=1
                   iu2=1
                  else if(ipss.eq.2)then
                   iu1=1
                   iu2=2
                  else
                   iu1=2
                   iu2=1
                  end if                                                    12d8s20
                  itestc=i2c                                              12d14s20
                  itesto=i2o                                              12d14s20
                  if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                   itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                   itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopeni+1                                      6d23s21
                   karg=iarg-1                                             12d8s20
                  else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                   itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopeni-1                                      6d23s21
                   karg=iarg                                              12d14s20
                  end if                                                        11d13s20
                  if(btest(itesto,nab4(2,iu2)))then                         12d8s20
                   itestc=ibset(itestc,nab4(2,iu2))                         12d8s20
                   itesto=ibclr(itesto,nab4(2,iu2))                             11d27s20
                   nopenk=nopenk-1                                              11d13s20
                   karg=karg+1                                                  11d13s20
                  else                                                          11d13s20
                   itesto=ibset(itesto,nab4(2,iu2))                         12d8s20
                   nopenk=nopenk+1                                              11d13s20
                  end if                                                        11d13s20
                  nqq=karg+mdon-1                                         4d14s21
                  if(nqq.ge.mdon.and.nqq.le.mdoo)then                     4d14s21
                   call gandc(i2c,i2o,itestc,itesto,nopeni,nopenk,         2d19s21
     $         iarg,karg,ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,        6d23s21
     $               iwpk1,ncsfmid1,bc,ibc)                             11d14s22
                   call gandc(itestc,itesto,j2c,j2o,nopenk,nopenj,         2d19s21
     $           karg,jarg,ncsf,norbxx,ixw1,ixw2,nnot2,nab2,iwpb2,      6d23s21
     $               iwpk2,ncsfmid2,bc,ibc)                             11d14s22
                   if(nnot1.eq.2.and.nnot2.eq.2)then                    6d22s21
                    lsa=ism(nab1(1))                                       12d14s20
                    lga=irel(nab1(1))-1                                    12d14s20
                    lsb=ism(nab1(2))                                       12d14s20
                    lgb=irel(nab1(2))-1                                    12d14s20
                    lsc=ism(nab2(1))                                     1d8s21
                    lgc=irel(nab2(1))-1                                  1d8s21
                    lsd=ism(nab2(2))                                     1d8s21
                    lgd=irel(nab2(2))-1                                  1d8s21
                    xint=getint(ioooo,lsa,lsb,lsc,lsd,lga,lgb,lgc,lgd,  11d15s22
     $                   bc,ibc)                                        11d15s22
                    if(min(nab1(1),nab1(2)).eq.min(nab2(1),nab2(2))     4d20s21
     $                    .and.                                         4d20s21
     $                max(nab1(1),nab1(2)).eq.max(nab2(1),nab2(2)))     1d12s21
     $                  xint=xint*0.5d0                                 1d12s21
                    call genmat(ncsf(iarg),ncsf(karg),ncsf(jarg),       6d22s21
     $                ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,        6d23s21
     $                iprod,bc,ibc)                                     11d10s22
                    do i=0,ncsf(iarg)*ncsf(jarg)-1                      6d22s21
                     bc(iprod+i)=bc(iprod+i)*xint                       6d22s21
                    end do                                              6d22s21
                    call hcdduc1(bc(iprod),ncsf,ncsf2,iarg,jarg,        6d22s21
     $                    bc(jvd),bc(jgg),jsbv12,nsymb,multh,nvirt,     6d24s21
     $                    nrootu,sr2,bc(iygoal),jgg,'f',4,bc,ibc)       11d10s22
                    ibcoff=iprod                                        6d22s21
                    if(ipss.eq.2)go to 3                                   12d18s20
                   end if                                                  12d14s20
                  end if                                                  12d14s20
    3             continue                                              6d22s21
                 end do                                                   12d14s20
                end if                                                     12d14s20
               end if                                                   6d22s21
               j2o=ibclr(j2o,norbx)                                      1d8s21
               j2o=ibset(j2o,norbxxx)                                    1d8s21
c                                                                       1d8s21
               call gandc4(i2c,l2o,j2c,j2o,nopeni,nopenj,norbxxx,       6d28s21
     $               nnot,nab4,bc,ibc)                                  11d14s22
               if(nnot.gt.0)then                                         1d8s21
c     jvd is root, ncsf,nvv,jff2.
c     re-order it to nvv,nroot,ncsf,iff2                                6d21s21
                if(itransv(nuse,jsb).eq.0)then                          6d22s21
                 call ddi_done(ibc(ircv(nuse,jsb)),nrcv(nuse,jsb))        2d12s21
                 itmp=ibcoff                                              6d21s21
                 ibcoff=itmp+nrow                                         6d21s21
                 call enough('hcdduc.  7',bc,ibc)
                 do iff=0,nfhere-1                                       6d22s21
                  jjvd=ivec(nuse,jsb)+nrow*iff                          6d26s21
                  jtmp=itmp                                               6d21s21
                  do jsbv1=1,nsymb                                               6d9s21
                   jsbv2=multh(jsbv1,jsbv12)                            6d24s21
                   if(jsbv2.ge.jsbv1)then                               6d24s21
                    if(jsbv1.eq.jsbv2)then                              6d24s21
                     do iv=0,nvirt(jsbv1)-1                             6d24s21
                      do i=0,ncsf2(1,jarg)-1                                    6d9s21
                       do ir=0,nrootm                                           6d9s21
                        kvd=jtmp+iv+nvirt(jsbv1)*(ir+nrootu*i)          6d24s21
                        bc(kvd)=bc(jjvd+ir)                             6d26s21
                       end do                                                   6d9s21
                       jjvd=jjvd+nrootu                                 6d26s21
                      end do                                                    6d9s21
                     end do                                                     6d9s21
                     nhere=nvirt(jsbv1)*ncsf2(1,jarg)*nrootu            6d24s21
                     jtmp=jtmp+nhere                                               6d16s21
                     nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2              6d24s21
                    else                                                         6d9s21
                     nvv=nvirt(jsbv1)*nvirt(jsbv2)                      6d24s21
                    end if                                                       6d9s21
                    do ivv=0,nvv-1                                              6d17s21
                     do i=0,ncsf(jarg)-1                                        6d17s21
                      do ir=0,nrootm                                            6d17s21
                       kvd=jtmp+ivv+nvv*(ir+nrootu*i)                            6d17s21
                       bc(kvd)=bc(jjvd+ir)                              6d26s21
                      end do                                                    6d9s21
                      jjvd=jjvd+nrootu                                  6d26s21
                     end do                                                     6d17s21
                    end do                                                      6d17s21
                    jtmp=jtmp+nvv*nrootu*ncsf(jarg)                               6d17s21
                   end if                                                        6d9s21
                  end do                                                         6d9s21
                  jjvd=ivec(nuse,jsb)+nrow*iff                          6d26s21
                  do i=0,nrow-1                                           6d21s21
                   bc(jjvd+i)=bc(itmp+i)                                6d26s21
                  end do                                                  6d21s21
                 end do                                                         6d16s21
                 ibcoff=itmp                                              6d21s21
                 itransv(nuse,jsb)=1                                      2d12s21
                end if                                                     2d12s21
                if(nab4(1,1).eq.0)then                                   1d8s21
                 write(6,*)('gandc4 falure!!! '),nab4
                 call dcbit(j2c,norb,'j2c')
                 call dcbit(j2o,norbxxx,'j2o')
                 call dcbit(i2c,norb,'i2c')
                 call dcbit(l2o,norbxxx,'l2o')
                 call dws_synca                                            11d27s20
                 call dws_finalize                                         11d27s20
                 stop 'hcdduc:f'
                end if                                                     11d27s20
                if(nnot.eq.2)then                                        1d8s21
                 if(isb.ne.jsb)then                                      1d9s21
                  write(6,*)('hey, what''s up with isb ne jsb ? ')       1d9s21
                  call dcbit(i2c,norbxxx,'i2c')
                  call dcbit(l2o,norbxxx,'l2o')
                  call dcbit(j2c,norbxxx,'j2c')
                  call dcbit(j2o,norbxxx,'j2o')
                  stop 'hcdduc:g'                                                   1d9s21
                 end if                                                  1d9s21
                 if(itransgg.eq.0.and.nnot.gt.0)then                       6d21s21
                  call ddi_done(ibc(iacc),nacc)                             2d12s21
                  itransgg=1                                                2d12s21
                  do i=igg,igg+mrow*nff2(ncloip,isb,1)-1                   6d21s21
                   bc(i)=0d0                                                2d12s21
                  end do                                                    2d12s21
                 end if                                                     2d12s21
                 call gandc(i2c,l2o,j2c,j2o,nopeni,nopenj,iarg,jarg,    6d28s21
     $              ncsf,norbxxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,      1d8s21
     $              ncsfmid1,bc,ibc)                                    11d14s22
                 iprod=ibcoff                                           6d22s21
                 ibcoff=iprod+ncsf(iarg)*ncsf(jarg)                     6d22s21
                 call enough('hcdduc.  8',bc,ibc)
                 call prodn(iwpb1,iwpk1,ncsf(iarg),ncsf(jarg),ncsfmid1,  6d22s21
     $               bc(iprod),bc,ibc,1d0,0d0)                          2d13s23
                 rms=0d0
                 do i=0,ncsf(iarg)-1
                  do j=0,i-1
                   ji=iprod+j+ncsf(iarg)*i
                   ij=iprod+i+ncsf(iarg)*j
                   rms=rms+bc(ji)**2+bc(ij)**2
                  end do
                  ii=iprod+i*(ncsf(iarg)+1)
                  rms=rms+(bc(ii)-1d0)**2
                 end do
                 rms=sqrt(rms/dfloat(ncsf(iarg)*ncsf(iarg)))
                 if(iarg.ne.jarg.or.rms.gt.1d-10)then
                  write(6,*)('hden prod is not unit matrix!!!'),iarg,
     $                 jarg,rms
                  call prntm2(bc(iprod),ncsf(iarg),ncsf(jarg),
     $                 ncsf(iarg))
                  call dws_synca
                  call dws_finalize
                  stop
                 end if
                 do kkk=0,nrow-1                                        6d28s21
                  bc(idenh+kkk)=bc(jvd+kkk)                             6d28s21
                 end do                                                 6d28s21
                 bc(mdenh)=1d0
                 do i=1,norb                                                   11d25s20
                  itest(i,2)=0                                                 11d25s20
                 end do                                                        11d25s20
                 ii=1                                                          11d25s20
                 do i=1,norb                                                   11d25s20
                  if(btest(j2c,i))then                                         11d25s20
                   itest(i,2)=2                                                11d25s20
                   ii=ii+1                                                     11d25s20
                  end if                                                       11d25s20
                 end do                                                        11d25s20
                 ii=1                                                          11d25s20
                 do i=1,norb                                                   11d25s20
                  if(btest(j2o,i))then                                         11d25s20
                   itest(i,2)=1                                                11d25s20
                   ii=ii+1                                                     11d25s20
                  end if                                                       11d25s20
                 end do                                                        11d25s20
                 nok=0                                                    1d8s21
                 do i=1,norb                                             1d8s21
                  if(itest(i,1).gt.0)then                                1d8s21
                   nok=nok+1                                             1d8s21
                   itest(nok,3)=itest(i,1)                               1d8s21
                   itest(nok,2)=i                                        1d8s21
                  end if                                                 1d8s21
                 end do                                                  1d8s21
                 itmpj=ibcoff                                           6d22s21
                 itmpk=itmpj+ncsf(iarg)*ncsf(jarg)                      6d22s21
                 itmpk2=itmpk+ncsf(iarg)*ncsf(jarg)                     7d2s21
                 ibcoff=itmpk2+ncsf(iarg)*ncsf(jarg)                    7d2s21
                 call enough('hcdduc.  9',bc,ibc)
                 do ii=0,ncsf(iarg)*ncsf(jarg)-1                        6d22s21
                  bc(itmpj+ii)=bc(iprod+ii)*2d0                         6d22s21
                  bc(itmpk+ii)=-bc(iprod+ii)                            6d22s21
                  bc(itmpk2+ii)=-bc(iprod+ii)                           7d2s21
                 end do                                                 6d22s21
                 do jj=ncsf2(1,jarg),ncsf(jarg)-1                       7d2s21
                  jtmpk2=itmpk2+ncsf(iarg)*jj                           7d2s21
                  do ii=0,ncsf(iarg)-1                                  7d2s21
                   bc(jtmpk2+ii)=-bc(jtmpk2+ii)                         7d2s21
                  end do                                                7d2s21
                 end do                                                 7d2s21
                 sz=0d0
                 do kkk=0,mrow-1
                  sz=sz+bc(jvd+kkk)**2
                 end do
                 sz=sqrt(sz/dfloat(mrow))
                 kdump=0
                 do i=1,nok                                              1d8s21
                  lsa=ism(itest(i,2))                                    1d8s21
                  lga=irel(itest(i,2))-1                                 1d8s21
                  icolj=((lga*(lga+1))/2)+lga                            1d8s21
                  icolk=lga+irefo(lsa)*lga                               1d8s21
                  bc(ndenj(lsa)+icolj)=1d0                              6d24s21
                  jden=idenj(lsa)+nmrow*icolj                           6d30s21
                  if(itest(i,3).eq.2)then                                1d8s21
                   bc(ndenk(lsa)+icolk)=1d0                             6d23s21
                   call hcdduc1(bc(itmpj),ncsf,ncsf2,iarg,jarg,           6d22s21
     $                    bc(jvd),bc(jden),jsbv12,nsymb,multh,nvirt,    6d24s21
     $                    nrootu,sr2,bc(iygoal),0,'h',4,bc,ibc)         11d10s22
                   kden=idenk(lsa,1)+nmrowk*icolk                       7d2s21
                   kden2=idenk(lsa,2)+nmrowk*icolk                      7d2s21
                   call hcdduc2(bc(itmpk),ncsf,ncsf2,iarg,jarg,bc(jvd), 7d2s21
     $                  bc(kden),jsbv12,nsymb,multh,nvirt,nrootu,sr2,   7d2s21
     $                  bc(iygoal),0,'hk1',4,bc,ibc)                    11d10s22
                   call hcdduc2(bc(itmpk2),ncsf,ncsf2,iarg,jarg,bc(jvd),7d2s21
     $                  bc(kden2),jsbv12,nsymb,multh,nvirt,nrootu,sr2,  7d2s21
     $                  bc(iygoal),0,'hk2',4,bc,ibc)                    11d10s22
                  else                                                  6d22s21
                   call hcdduc1(bc(iprod),ncsf,ncsf2,iarg,jarg,         6d22s21
     $                    bc(jvd),bc(jden),jsbv12,nsymb,multh,nvirt,    6d24s21
     $                    nrootu,sr2,bc(iygoal),0,'j',4,bc,ibc)         11d10s22
                  end if                                                6d22s21
                 end do                                                 6d22s21
                 ibcoff=iprod                                           6d22s21
                 do i=1,nok                                              1d8s21
                  if(itest(i,3).eq.1)then                                  12d14s20
                   itestc=i2c                                           6d22s21
                   itesto=l2o                                           6d28s21
                   nopenk=nopeni                                        6d22s21
c
c     i is bra, j is ket.
c
c         I   J    K     thus  IK  KJ
c     b   n+1 n    n+1        g   g   (ck|bc)
c     k   m   m+1  m+1         ck  bc
c     c   p   p    p-1
c
c
c     anihilate common
c
                   if(btest(itestc,itest(i,2)))then                             11d13s20
                    itestc=ibclr(itestc,itest(i,2))                             11d13s20
                    itesto=ibset(itesto,itest(i,2))                             11d13s20
                    karg=iarg-1                                                11d13s20
                    nopenk=nopenk+1                                             11d13s20
                   else                                                         11d13s20
                    itesto=ibclr(itesto,itest(i,2))                             11d13s20
                    karg=iarg                                                  11d13s20
                    nopenk=nopenk-1                                             11d13s20
                   end if                                                       11d13s20
c
c     create ket
c
                   if(btest(itesto,nab4(2,1)))then                        12d14s20
                    itestc=ibset(itestc,nab4(2,1))                        12d14s20
                    itesto=ibclr(itesto,nab4(2,1))                        12d14s20
                    karg=karg+1                                                11d13s20
                    nopenk=nopenk-1                                             11d13s20
                   else                                                         11d13s20
                    itesto=ibset(itesto,nab4(2,1))                        12d14s20
                    nopenk=nopenk+1                                             11d13s20
                   end if                                                       11d13s20
                   nqq=karg+mdon-1
                   if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                    call gandc(i2c,l2o,itestc,itesto,nopeni,nopenk,     6d28s21
     $            iarg,karg,ncsf,norbxxx,ixw1,ixw2,nnot1,nab1,iwpb1,    6d22s21
     $               iwpk1,ncsfmid1,bc,ibc)                             11d14s22
                    call gandc(itestc,itesto,j2c,j2o,nopenk,nopenj,     6d22s21
     $            karg,jarg,ncsf,norbxxx,ixw1,ixw2,nnot2,nab2,iwpb2,     12d14s20
     $                iwpk2,ncsfmid2,bc,ibc)                            11d14s22
                    if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                      ia=nab2(2)
                      ib=nab1(1)
                     lsa=ism(ia)                                          1d8s21
                     lga=irel(ia)-1                                       1d8s21
                     lsb=ism(ib)                                          1d8s21
                     lgb=irel(ib)-1                                       1d8s21
                     lsab=multh(lsa,lsb)                                  1d9s21
                     icol=lga+irefo(lsa)*lgb                              1d8s21
                     if(lsab.ne.ijsb)then                                 1d9s21
                      write(6,*)('hey, lsab = '),lsab,(' ne ijsb = '),
     $                  ijsb
                      stop 'hcdduc:i'
                     end if                                               1d9s21
                     call genmat(ncsf(iarg),ncsf(karg),ncsf(jarg),       6d22s21
     $                ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,        6d23s21
     $                iprod,bc,ibc)                                     11d10s22
                     bc(ndenk(lsa)+icol)=1d0                            6d23s21
                     kden=idenk(lsa,1)+nmrowk*icol                      7d2s21
                     call hcdduc2(bc(iprod),ncsf,ncsf2,iarg,jarg,         6d22s21
     $                    bc(jvd),bc(kden),jsbv12,nsymb,multh,nvirt,    6d24s21
     $                    nrootu,sr2,bc(iygoal),0,'k1',4,bc,ibc)        11d10s22
                     do jjj=ncsf2(1,jarg),ncsf(jarg)-1                  7d1s21
                      jprod=iprod+ncsf(iarg)*jjj                        7d1s21
                      do iii=0,ncsf(iarg)-1                             7d1s21
                       bc(jprod+iii)=-bc(jprod+iii)                     7d1s21
                      end do                                            7d1s21
                     end do                                             7d1s21
                     kden=idenk(lsa,2)+nmrowk*icol                      7d2s21
                     call hcdduc2(bc(iprod),ncsf,ncsf2,iarg,jarg,       7d2s21
     $                    bc(jvd),bc(kden),jsbv12,nsymb,multh,nvirt,    6d24s21
     $                    nrootu,sr2,bc(iygoal),0,'k2',4,bc,ibc)        11d10s22
                     ibcoff=iprod                                       7d1s21
                    end if                                               6d22s21
                   end if
                  end if                                                6d22s21
                 end do                                                 6d22s21
                else if(nnot.gt.3)then                                  6d22s21
                 if(nnot.eq.3)then                                          12d8s20
                  ipssx=1                                                   12d8s20
                 else                                                       12d8s20
                  ipssx=3                                                 12d18s20
                 end if                                                     12d8s20
                 do ipss=1,ipssx                                            12d8s20
                  if(ipss.eq.1)then
                   iu1=1
                   iu2=1
                  else if(ipss.eq.2)then
                   iu1=1
                   iu2=2
                  else
                   iu1=2
                   iu2=1
                  end if                                                    12d8s20
                  itestc=i2c                                              12d8s20
                  itesto=l2o                                            6d28s21
                  if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                   itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                   itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopeni+1                                              11d13s20
                   karg=iarg-1                                             12d8s20
                  else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                   itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopeni-1                                              11d13s20
                   karg=iarg
                  end if                                                        11d13s20
                  if(btest(itesto,nab4(2,iu2)))then                         12d8s20
                   itestc=ibset(itestc,nab4(2,iu2))                         12d8s20
                   itesto=ibclr(itesto,nab4(2,iu2))                             11d27s20
                   nopenk=nopenk-1                                              11d13s20
                   karg=karg+1                                                  11d13s20
                  else                                                          11d13s20
                   itesto=ibset(itesto,nab4(2,iu2))                         12d8s20
                   nopenk=nopenk+1                                              11d13s20
                  end if                                                        11d13s20
                  nqq=karg+mdon-1                                       4d20s21
                  if(nqq.ge.mdon.and.nqq.le.mdoo)then                   4d20s21
                   call gandc(i2c,l2o,itestc,itesto,nopeni,nopenk,      6d28s21
     $         iarg,karg,ncsf,norbxxx,ixw1,ixw2,nnot1,nab1,iwpb1,       6d22s21
     $                iwpk1,ncsfmid1,bc,ibc)                            11d14s22
                   call gandc(itestc,itesto,j2c,j2o,nopenk,nopenj,      6d22s21
     $         karg,jarg,ncsf,norbxxx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,   12d18s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
                   if(nnot1.eq.2.and.nnot2.eq.2)then                       12d19s20
                    call genmat(ncsf(iarg),ncsf(karg),ncsf(jarg),       6d22s21
     $                ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,        6d23s21
     $                iprod,bc,ibc)                                     11d10s22
                    if(nab1(2).gt.norb)then                               1d8s21
                      ia=nab1(1)                                          1d8s21
                      ib=nab2(2)                                          1d8s21
                     lsa=ism(ia)                                           1d8s21
                     lsb=ism(ib)                                           1d8s21
                     lga=irel(ia)-1                                        1d8s21
                     lgb=irel(ib)-1                                        1d8s21
                     icol=lga+irefo(lsa)*lgb                              1d8s21
                     bc(ndenk(lsa)+icol)=1d0                            6d23s21
                     kden=idenk(lsa,1)+nmrowk*icol                      7d2s21
                     call hcdduc2(bc(iprod),ncsf,ncsf2,iarg,jarg,       7d2s21
     $                    bc(jvd),bc(kden),jsbv12,nsymb,multh,nvirt,    6d24s21
     $                    nrootu,sr2,bc(iygoal),0,'l',4,bc,ibc)         11d10s22
                     do jjj=ncsf2(1,jarg),ncsf(jarg)-1                  7d1s21
                      jprod=iprod+ncsf(iarg)*jjj                        7d1s21
                      do iii=0,ncsf(iarg)-1                             7d1s21
                       bc(jprod+iii)=-bc(jprod+iii)                     7d1s21
                      end do                                            7d1s21
                     end do                                             7d1s21
                     kden=idenk(lsa,2)+nmrowk*icol                      7d2s21
                     call hcdduc2(bc(iprod),ncsf,ncsf2,iarg,jarg,       7d2s21
     $                    bc(jvd),bc(kden),jsbv12,nsymb,multh,nvirt,    6d24s21
     $                    nrootu,sr2,bc(iygoal),0,'l2',4,bc,ibc)        11d10s22
                    else                                                  1d8s21
                     ia=min(nab1(1),nab1(2))                              1d8s21
                     ib=max(nab1(1),nab1(2))                              1d8s21
                     lsa=ism(ia)                                           1d8s21
                     lsb=ism(ib)                                           1d8s21
                     lga=irel(ia)-1                                        1d8s21
                     lgb=irel(ib)-1                                        1d8s21
                     if(lsa.eq.lsb)then                                   1d8s21
                      icol=((lgb*(lgb+1))/2)+lga                          1d10s21
                     else                                                 1d8s21
                      icol=lgb+irefo(lsb)*lga                             1d9s21
                     end if                                               1d8s21
                     bc(ndenj(lsa)+icol)=1d0                            6d23s21
                     jden=idenj(lsa)+nmrow*icol                         6d30s21
                     call hcdduc1(bc(iprod),ncsf,ncsf2,iarg,jarg,         6d22s21
     $                    bc(jvd),bc(jden),jsbv12,nsymb,multh,nvirt,    6d24s21
     $                    nrootu,sr2,bc(iygoal),0,'m',4,bc,ibc)         11d10s22
                    end if                                                1d8s21
                    ibcoff=iprod                                        6d22s21
                    if(ipss.eq.2)go to 4                                6d22s21
                   end if                                               6d22s21
                  end if                                                6d22s21
                 end do                                                 6d22s21
    4            continue                                               6d22s21
                end if                                                  6d22s21
               end if                                                   6d22s21
              else                                                      6d22s21
               jf2o=jf2o+2                                              6d22s21
              end if                                                    6d22s21
              if(isb.eq.jsb.and.ncloi.eq.ncloj.and.if2.eq.jf2)then      7d1s21
               if(itransgg.eq.0)then                                    7d1s21
                call ddi_done(ibc(iacc),nacc)                             2d12s21
                itransgg=1                                                2d12s21
                do i=igg,igg+mrow*nff2(ncloip,isb,1)-1                   6d21s21
                 bc(i)=0d0                                                2d12s21
                end do                                                    2d12s21
               end if                                                     2d12s21
               ivdtmp=ibcoff                                            7d1s21
               itmp=ivdtmp+nrow                                         7d1s21
               ibcoff=itmp+nrow                                         7d1s21
               call enough('hcdduc. 10',bc,ibc)
               i18=1                                                    7d1s21
               i28=nrow                                                 7d1s21
               i38=if2                                                  7d1s21
               call ddi_get(bc,ibc,ihddiag(ncloip,isb,1),i18,i28,i38,   11d15s22
     $              i38,bc(itmp))                                       11d15s22
c     4x contributions: Gvv'ir+/-=[(vu|v'u')+/-(vu'|v'u)]Vuu'ir+/-
c     /sqrt((1+delta vv')(1+delta uu'))
               jjvd=itmp                                                7d1s21
               do i=0,nrow-1                                            7d2s21
                bc(ivdtmp+i)=0d0                                        7d2s21
               end do                                                   7d2s21
               nna=nrootu*ncsf(iarg)                                    7d3s21
               nnam=nna-1                                               7d3s21
               nn=nrootu*ncsf2(1,iarg)                                  7d3s21
               nnm=nn-1                                                 7d3s21
               do jsbv1=1,nsymb                                         7d1s21
                jsbv2=multh(jsbv1,jsbv12)                               7d1s21
                if(jsbv2.ge.jsbv1)then                                  7d1s21
                 if(jsbv1.eq.jsbv2)then                                 7d2s21
                  jjvdp=jjvd+nrootu*ncsf2(1,iarg)*nvirt(jsbv1)          7d2s21
                  nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                 7d1s21
                  isw=0                                                 7d2s21
                 else                                                   7d2s21
                  jjvdp=jjvd                                            7d2s21
                  isw=1                                                 7d2s21
                  nvv=nvirt(jsbv1)*nvirt(jsbv2)                         7d1s21
                 end if                                                 7d2s21
                 ktmp=ivdtmp                                              7d1s21
                 do ksbv1=1,nsymb                                       7d2s21
                  ksbv2=multh(ksbv1,jsbv12)                             7d2s21
                  if(ksbv2.ge.ksbv1)then                                7d2s21
                   call ilimts(nvirt(ksbv2),nvirt(jsbv2),mynprocg,      7d2s21
     $                 mynowprog,ilj,ihj,i1s,i1e,i2s,i2e)               7d2s21
                   nherej=ihj+1-ilj                                     7d2s21
                   i2euj=inv(1,ksbv1,jsbv1,ksbv2)                       7d2s21
                   icasej=inv(2,ksbv1,jsbv1,ksbv2)                      7d2s21
                   if(icasej.eq.1)then                                  7d3s21
                    iscj=0                                              7d3s21
                   else                                                 7d3s21
                    iscj=1                                              7d3s21
                   end if                                               7d3s21
                   if(ksbv1.eq.jsbv1)then                               7d3s21
                    nnrowj=(nvirt(ksbv1)*(nvirt(ksbv1)+1))/2             7d3s21
                    jksj=0                                              7d3s21
                   else                                                 7d3s21
                    nnrowj=nvirt(ksbv1)*nvirt(jsbv1)                     7d3s21
                    jksj=1                                              7d3s21
                   end if                                               7d3s21
                   call ilimts(nvirt(ksbv2),nvirt(jsbv1),mynprocg,      7d2s21
     $                 mynowprog,ilk,ihk,i1s,i1e,i2s,i2e)               7d2s21
                   nherek=ihk+1-ilk                                     7d2s21
                   i2euk=inv(1,ksbv1,jsbv2,ksbv2)                       7d2s21
                   icasek=inv(2,ksbv1,jsbv2,ksbv2)                      7d2s21
                   if(icasek.eq.1)then                                  7d4s21
                    isck=0                                              7d3s21
                   else                                                 7d3s21
                    isck=1                                              7d3s21
                   end if                                               7d3s21
                   if(ksbv1.eq.jsbv2)then                               7d3s21
                    nnrowk=(nvirt(ksbv1)*(nvirt(ksbv1)+1))/2             7d3s21
                    jksk=0                                              7d3s21
                   else                                                 7d3s21
                    nnrowk=nvirt(ksbv1)*nvirt(jsbv2)                     7d3s21
                    jksk=1                                              7d3s21
                   end if                                               7d3s21
                   if(jsbv12.eq.1)then                                  7d2s21
                    mvv=(nvirt(ksbv1)*(nvirt(ksbv1)-1))/2               7d3s21
                    ktmpp=ktmp+nrootu*nvirt(ksbv1)*ncsf2(1,iarg)        7d3s21
                    do jv=0,nvirt(jsbv1)-1                                7d1s21
                     jjvdl=jjvd+nn*jv                                   7d3s21
                     do kv2=0,nvirt(ksbv2)-1                             7d2s21
                      ktmpl=ktmp+nn*kv2                                 7d3s21
                      ntop=(kv2+isw*(nvirt(ksbv1)-kv2))-1               7d3s21
                      xinta=get4int(i4xb,ksbv1,jsbv1,ksbv2,jsbv2,       7d6s21
     $                      kv2,jv,kv2,jv,nvirt,idum,bc,ibc)            11d10s22
                      if(xinta.ne.0d0)then                              7d6s21
                       do ii=0,nnm                                      7d3s21
                        bc(ktmpl+ii)=bc(ktmpl+ii)+xinta*bc(jjvdl+ii)    7d6s21
                       end do                                           7d3s21
                      end if                                            7d6s21
                      do kv1=0,ntop                                     7d2s21
                       xinta=get4int(i4xb,ksbv1,jsbv1,ksbv2,jsbv2,      7d6s21
     $                      kv1,jv,kv2,jv,nvirt,idum,bc,ibc)*sr2        11d10s22
                       if(xinta.ne.0d0)then                             7d6s21
                        ktri=((kv2*(kv2-1))/2)+kv1                      7d3s21
                        krec=kv1+nvirt(ksbv1)*kv2                       7d3s21
                        kcol=ktri+isw*(krec-ktri)                       7d3s21
                        ktmpl=ktmpp+nna*kcol                            7d3s21
                        do ii=0,nnm                                      7d3s21
                         bc(ktmpl+ii)=bc(ktmpl+ii)+xinta*bc(jjvdl+ii)   7d6s21
                        end do                                           7d3s21
                       end if                                           7d6s21
                      end do                                            7d6s21
                     end do                                             7d3s21
                    end do                                              7d3s21
                    do kv=0,nvirt(ksbv1)-1                              7d3s21
                     ktmpl=ktmp+nn*kv                                   7d3s21
                     do jv2=0,nvirt(jsbv2)-1                            7d3s21
                      ntop=(jv2+isw*(nvirt(jsbv1)-jv2))-1               7d3s21
                       icol=kv+1+nvirt(ksbv2)*jv2                        7d4s21
                       ixint=i4x(i2euj)+nnrowj*(icol-ilj)                7d3s21
                      do jv1=0,ntop                                     7d6s21
                       xinta=get4int(i4xb,ksbv1,jsbv1,ksbv2,jsbv2,       7d6s21
     $                      kv,jv1,kv,jv2,nvirt,idum,bc,ibc)*sr2        11d10s22
                       if(xinta.ne.0d0)then                              7d6s21
                        jtri=((jv2*(jv2-1))/2)+jv1                      7d3s21
                        jrec=jv1+nvirt(jsbv1)*jv2                       7d3s21
                        jcol=jtri+isw*(jrec-jtri)                       7d3s21
                        jjvdl=jjvdp+nna*jcol                            7d3s21
                        do ii=0,nnm                                     7d3s21
                         bc(ktmpl+ii)=bc(ktmpl+ii)+xinta*bc(jjvdl+ii)    7d3s21
                        end do                                          7d3s21
                       end if                                           7d6s21
                      end do                                            7d6s21
                     end do                                             7d3s21
                    end do                                              7d3s21
                   else                                                 7d3s21
                    ktmpp=ktmp                                          7d3s21
                    mvv=nvirt(ksbv1)*nvirt(ksbv2)                       7d3s21
                   end if                                               7d3s21
                   do kv2=0,nvirt(ksbv2)-1                              7d3s21
                    ktop=(kv2+isw*(nvirt(ksbv1)-kv2))-1                 7d3s21
                    ktri=((kv2*(kv2-1))/2)                              7d3s21
                    krec=nvirt(ksbv1)*kv2                               7d3s21
                    krow=ktri+isw*(krec-ktri)                           7d3s21
                    do jv2=0,nvirt(jsbv2)-1                             7d3s21
                     jtop=(jv2+isw*(nvirt(jsbv1)-jv2))-1                7d3s21
                     jtri=((jv2*(jv2-1))/2)                             7d3s21
                     jrec=nvirt(jsbv1)*jv2                              7d3s21
                     jrow=jtri+isw*(jrec-jtri)                          7d3s21
                     do kv1=0,ktop                                      7d6s21
                      ktmpl=ktmpp+nna*(kv1+krow)                        7d6s21
                      do jv1=0,jtop                                     7d6s21
                       xinta=get4int(i4xb,ksbv1,jsbv1,ksbv2,jsbv2,      7d6s21
     $                      kv1,jv1,kv2,jv2,nvirt,idum,bc,ibc)          11d10s22
                       if(xinta.ne.0d0)then                             7d6s21
                        jjvdl=jjvdp+nna*(jv1+jrow)                      7d3s21
                        do ii=0,nnam                                    7d3s21
                         bc(ktmpl+ii)=bc(ktmpl+ii)+xinta*bc(jjvdl+ii)    7d3s21
                        end do                                          7d3s21
                       end if                                           7d6s21
                      end do                                            7d3s21
                     end do                                             7d3s21
                    end do                                              7d3s21
                    do jv1=0,nvirt(jsbv1)-1                             7d3s21
                     jbot=(jv1+1)*(1-isw)                               7d3s21
                     do jv2=jbot,nvirt(jsbv2)-1                         7d6s21
                      jtri=((jv2*(jv2-1))/2)+jv1                        7d6s21
                      jrec=jv1+nvirt(jsbv1)*jv2                         7d6s21
                      jcol=jtri+isw*(jrec-jtri)                         7d6s21
                      jjvdl=jjvdp+nna*jcol                              7d6s21
                      do kv1=0,ktop
                       xinta=get4int(i4xb,ksbv1,jsbv2,ksbv2,jsbv1,      7d6s21
     $                      kv1,jv2,kv2,jv1,nvirt,idum,bc,ibc)          11d10s22
                       if(xinta.ne.0d0)then                             7d6s21
                        ktmpl=ktmpp+nna*(kv1+krow)                      7d3s21
                        do ii=0,nnm                                     7d4s21
                         bc(ktmpl+ii)=bc(ktmpl+ii)+xinta*bc(jjvdl+ii)   7d6s21
                        end do                                          7d3s21
                        do ii=nn,nnam                                   7d4s21
                         bc(ktmpl+ii)=bc(ktmpl+ii)-xinta*bc(jjvdl+ii)   7d6s21
                        end do                                          7d3s21
                       end if                                           7d6s21
                      end do                                            7d3s21
                     end do                                             7d3s21
                    end do                                              7d3s21
                   end do                                               7d3s21
                   ktmp=ktmpp+mvv*nna                                   7d3s21
                  end if                                                7d3s21
                 end do                                                 7d3s21
                 jjvd=jjvdp+nvv*nna                                     7d3s21
                end if                                                  7d1s21
               end do                                                   7d1s21
               jtmp=ivdtmp                                              7d3s21
               kgg=igg+mrow*(if2-1)                                     7d3s21
               do jsbv1=1,nsymb                                         7d3s21
                jsbv2=multh(jsbv1,jsbv12)                               7d3s21
                if(jsbv2.ge.jsbv1)then                                  7d3s21
                 if(jsbv1.eq.jsbv2)then                                 7d3s21
                  jtmpp=jtmp+nvirt(jsbv1)*nn                            7d3s21
                  kggp=kgg+nvirt(jsbv1)*nn                              7d3s21
                  nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                 7d3s21
                  do iv=0,nvirt(jsbv1)-1                                7d3s21
                   do ii=0,nnm                                          7d3s21
                    ij=kgg+iv+nvirt(jsbv1)*ii                           7d3s21
                    ji=jtmp+ii+nn*iv                                    7d3s21
                    bc(ij)=bc(ij)+bc(ji)                                7d3s21
                   end do                                               7d3s21
                  end do                                                7d3s21
                 else                                                   7d3s21
                  jtmpp=jtmp                                            7d3s21
                  kggp=kgg                                              7d3s21
                  nvv=nvirt(jsbv1)*nvirt(jsbv2)                         7d3s21
                 end if                                                 7d3s21
                 nvvm=nvv-1                                             7d3s21
                 do ivv=0,nvvm
                  do ii=0,nnam                                          7d3s21
                   ij=kggp+ivv+nvv*ii                                   7d3s21
                   ji=jtmpp+ii+nna*ivv                                  7d3s21
                   bc(ij)=bc(ij)+bc(ji)                                 7d3s21
                  end do                                                7d3s21
                 end do                                                 7d3s21
                 jtmp=jtmpp+nvv*nna                                     7d3s21
                 kgg=kggp+nvv*nna                                       7d3s21
                end if                                                  7d3s21
               end do                                                   7d3s21
c
               jjvd=itmp                                                7d1s21
               jtmp=ivdtmp                                              7d1s21
               do jsbv1=1,nsymb                                         7d1s21
                jsbv2=multh(jsbv1,jsbv12)                               7d1s21
                if(jsbv2.ge.jsbv1)then                                  7d1s21
                 if(jsbv1.eq.jsbv2)then                                 7d1s21
                  do iv=0,nvirt(jsbv1)-1                                7d1s21
                   do i=0,ncsf2(1,iarg)-1                               7d1s21
                    do ir=0,nrootm                                      7d1s21
                     kvd=jtmp+iv+nvirt(jsbv1)*(ir+nrootu*i)             7d1s21
                     bc(kvd)=bc(jjvd+ir)                                7d1s21
                    end do                                              7d1s21
                    jjvd=jjvd+nrootu                                    7d1s21
                   end do                                               7d1s21
                  end do                                                7d1s21
                  nhere=nvirt(jsbv1)*ncsf2(1,iarg)*nrootu               7d1s21
                  jtmp=jtmp+nhere                                       7d1s21
                  nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                 7d1s21
                 else                                                   7d1s21
                  nvv=nvirt(jsbv1)*nvirt(jsbv2)                         7d1s21
                 end if                                                 7d1s21
                 do ivv=0,nvv-1                                         7d1s21
                  do i=0,ncsf(iarg)-1                                   7d1s21
                   do ir=0,nrootm                                       7d1s21
                    kvd=jtmp+ivv+nvv*(ir+nrootu*i)                      7d1s21
                    bc(kvd)=bc(jjvd+ir)                                 7d1s21
                   end do                                               7d1s21
                   jjvd=jjvd+nrootu                                     7d1s21
                  end do                                                7d1s21
                 end do                                                 7d1s21
                 jtmp=jtmp+nvv*nrootu*ncsf(iarg)                        7d1s21
                end if                                                  7d1s21
               end do                                                   7d1s21
               ibcoff=itmp                                              7d1s21
               do i=1,norb                                               7d1s21
                if(btest(i2c,i).or.btest(l2o,i))then                    7d1s21
                 lsa=ism(i)                                             7d1s21
                 lga=irel(i)-1                                          7d1s21
                 icolj=((lga*(lga+1))/2)+lga                            1d8s21
                 if(btest(i2c,i))then
                  icolk=lga+irefo(lsa)*lga                              7d2s21
                  kgg=igg+mrow*(if2-1)                                    6d21s21
                  kvd=ivdtmp                                              7d1s21
                  do ksbv1=1,nsymb                                      7d1s21
                   ksbv2=multh(isbv12,ksbv1)                            7d1s21
                   if(ksbv2.ge.ksbv1)then                               7d1s21
                    i2eu1=inv(1,lsa,lsa,ksbv1)                          7d1s21
                    k2eu1=invk1(1,lsa,lsa,ksbv1,1)                      7d2s21
                    call ilimts(nvirt(ksbv1),nvirt(ksbv1),mynprocg,     7d1s21
     $                    mynowprog,ilk,ihk,i1sk,i1ek,i2sk1,i2ek1)      7d1s21
                    n1n=nvirt(ksbv1)                                    7d1s21
                    n1x=-1                                              8d8s24
                    do iv1=i2sk1,i2ek1                                  7d1s21
                     irow=iv1+nvirt(ksbv1)*(iv1-1)                      7d1s21
                     if(irow.ge.ilk.and.irow.le.ihk)then                7d1s21
                      n1n=min(n1n,iv1)                                  7d1s21
                      n1x=max(n1x,iv1)                                  7d1s21
                     end if                                             7d1s21
                    end do                                              7d1s21
                    nhere1=ihk+1-ilk                                    7d1s21
                    jint1=jmats(i2eu1)+nhere1*icolj-ilk                 7d1s21
                    kint1=kmats(k2eu1)+nhere1*icolk-ilk                 7d2s21
                    i2eu2=inv(1,lsa,lsa,ksbv2)                          7d1s21
                    k2eu2=invk1(1,lsa,lsa,ksbv2,1)                      7d2s21
                    call ilimts(nvirt(ksbv2),nvirt(ksbv2),mynprocg,     7d1s21
     $                    mynowprog,ilk,ihk,i1sk,i1ek,i2sk2,i2ek2)      7d1s21
                    nhere2=ihk+1-ilk                                    7d1s21
                     jint2=jmats(i2eu2)+nhere2*icolj-ilk                 7d1s21
                    kint2=kmats(k2eu2)+nhere2*icolk-ilk                 7d2s21
                    n2n=nvirt(ksbv2)                                    7d1s21
                    n2x=-1                                              8d8s24
                    do iv2=i2sk2,i2ek2                                  7d1s21
                     irow=iv2+nvirt(ksbv2)*(iv2-1)                      7d1s21
                     if(irow.ge.ilk.and.irow.le.ihk)then                7d1s21
                      n2n=min(n2n,iv2)                                  7d1s21
                      n2x=max(n2x,iv2)                                  7d1s21
                     end if                                             7d1s21
                    end do                                              7d1s21
                    if(ksbv1.eq.ksbv2)then                              7d1s21
                     do j=0,ncsf2(1,iarg)-1                             7d1s21
                      do ir=0,nrootm                                    7d1s21
                       do ivp=n1n,n1x                                   7d1s21
                        iv=ivp-1                                        7d1s21
                        ivv=ivp+nvirt(ksbv1)*iv                         7d1s21
                        bc(kgg+iv)=bc(kgg+iv)-2d0*(2d0*bc(jint1+ivv)    7d2s21
     $                       -bc(kint1+ivv)                             7d2s21
     $                       )*bc(kvd+iv)                               7d2s21
                       end do                                           7d1s21
                       kgg=kgg+nvirt(ksbv1)                             7d1s21
                       kvd=kvd+nvirt(ksbv1)                             7d1s21
                      end do                                            7d1s21
                     end do                                             7d1s21
                     isw=0                                              7d1s21
                     kvv=(nvirt(ksbv1)*(nvirt(ksbv1)-1))/2              7d1s21
                    else                                                7d1s21
                     isw=1                                              7d1s21
                     kvv=nvirt(ksbv1)*nvirt(ksbv2)                      7d1s21
                    end if                                              7d1s21
                    do j=0,ncsf(iarg)-1                                 7d1s21
                     do ir=0,nrootm                                     7d1s21
                      do iv2=0,nvirt(ksbv2)-1                           7d1s21
                       ntop=iv2+isw*(nvirt(ksbv1)-iv2)                  7d1s21
                       do iv1p=n1n,min(ntop,n1x)                        7d1s21
                        iv1=iv1p-1                                      7d1s21
                        ivv=iv1p+nvirt(ksbv1)*iv1                       7d1s21
                        itri=((iv2*(iv2-1))/2)+iv1                      7d1s21
                        irec=iv1+nvirt(ksbv1)*iv2                       7d1s21
                        irow=itri+isw*(irec-itri)                       7d1s21
                        bc(kgg+irow)=bc(kgg+irow)-(2d0*bc(jint1+ivv)     7d1s21
     $                       -bc(kint1+ivv)                             7d2s21
     $                       )*bc(kvd+irow)                             7d2s21
                       end do                                           7d1s21
                      end do                                            7d1s21
                      do iv2p=n2n,n2x                                   7d1s21
                       iv2=iv2p-1                                       7d1s21
                       ntop=(iv2+isw*(nvirt(ksbv1)-iv2))-1               7d1s21
                       ivv=iv2p+nvirt(ksbv2)*iv2                        7d1s21
                       xint=bc(jint2+ivv)*2d0                           7d2s21
     $                      -bc(kint2+ivv)                              7d2s21
                       do iv1=0,ntop                                    7d1s21
                        itri=((iv2*(iv2-1))/2)+iv1                      7d1s21
                        irec=iv1+nvirt(ksbv1)*iv2                       7d1s21
                        irow=itri+isw*(irec-itri)                       7d1s21
                        bc(kgg+irow)=bc(kgg+irow)-xint*bc(kvd+irow)     7d1s21
                       end do                                           7d1s21
                      end do                                            7d1s21
                      kgg=kgg+kvv                                       7d1s21
                      kvd=kvd+kvv                                       7d1s21
                     end do                                             7d1s21
                    end do                                              7d1s21
                   end if                                               7d1s21
                  end do                                                7d1s21
                 else
                  kgg=igg+mrow*(if2-1)                                    6d21s21
                  kvd=ivdtmp                                            7d1s21
                  do ksbv1=1,nsymb                                      7d1s21
                   ksbv2=multh(isbv12,ksbv1)                            7d1s21
                   if(ksbv2.ge.ksbv1)then                               7d1s21
                    i2eu1=inv(1,lsa,lsa,ksbv1)                          7d1s21
                    k2eu1=invk1(1,lsa,lsa,ksbv1,1)                      7d2s21
                    call ilimts(nvirt(ksbv1),nvirt(ksbv1),mynprocg,     7d1s21
     $                    mynowprog,ilk,ihk,i1sk,i1ek,i2sk1,i2ek1)      7d1s21
                    n1n=nvirt(ksbv1)                                    7d1s21
                    n1x=-1                                              8d8s24
                    do iv1=i2sk1,i2ek1                                  7d1s21
                     irow=iv1+nvirt(ksbv1)*(iv1-1)                      7d1s21
                     if(irow.ge.ilk.and.irow.le.ihk)then                7d1s21
                      n1n=min(n1n,iv1)                                  7d1s21
                      n1x=max(n1x,iv1)                                  7d1s21
                     end if                                             7d1s21
                    end do                                              7d1s21
                    nhere1=ihk+1-ilk                                    7d1s21
                    jint1=jmats(i2eu1)+nhere1*icolj-ilk                 7d1s21
                    icolk=lga+irefo(lsa)*lga                            7d2s21
                    kint1=kmats(k2eu1)+nhere1*icolk-ilk                 7d2s21
                    i2eu2=inv(1,lsa,lsa,ksbv2)                          7d1s21
                    call ilimts(nvirt(ksbv2),nvirt(ksbv2),mynprocg,     7d1s21
     $                    mynowprog,ilk,ihk,i1sk,i1ek,i2sk2,i2ek2)       7d1s21
                    nhere2=ihk+1-ilk                                    7d1s21
                    n2n=nvirt(ksbv2)                                    7d1s21
                    n2x=-1                                              8d8s24
                    do iv2=i2sk2,i2ek2                                  7d1s21
                     irow=iv2+nvirt(ksbv2)*(iv2-1)                      7d1s21
                     if(irow.ge.ilk.and.irow.le.ihk)then                7d1s21
                      n2n=min(n2n,iv2)                                  7d1s21
                      n2x=max(n2x,iv2)                                  7d1s21
                     end if                                             7d1s21
                    end do                                              7d1s21
                    jint2=jmats(i2eu2)+nhere2*icolj-ilk                 7d1s21
                    if(ksbv1.eq.ksbv2)then                              7d1s21
                     do j=0,ncsf2(1,iarg)-1                             7d1s21
                      do ir=0,nrootm                                    7d1s21
                       do ivp=n2n,n2x                                   7d1s21
                        iv=ivp-1                                        7d1s21
                        ivv=ivp+nvirt(ksbv1)*iv                         7d1s21
                        bc(kgg+iv)=bc(kgg+iv)-(2d0*bc(jint1+ivv)        7d2s21
     $                       -bc(kint1+ivv)                             7d2s21
     $                       )*bc(kvd+iv)                               7d2s21
                       end do                                           7d1s21
                       kgg=kgg+nvirt(ksbv1)                             7d1s21
                       kvd=kvd+nvirt(ksbv1)                             7d1s21
                      end do                                            7d1s21
                     end do                                             7d1s21
                     isw=0                                              7d1s21
                     kvv=(nvirt(ksbv1)*(nvirt(ksbv1)-1))/2              7d1s21
                    else                                                7d1s21
                     isw=1                                              7d1s21
                     kvv=nvirt(ksbv1)*nvirt(ksbv2)                      7d1s21
                    end if                                              7d1s21
                    do j=0,ncsf(iarg)-1                                 7d1s21
                     do ir=0,nrootm                                     7d1s21
                      do iv2=0,nvirt(ksbv2)-1                           7d1s21
                       ntop=iv2+isw*(nvirt(ksbv1)-iv2)                  7d1s21
                       do iv1p=n1n,min(ntop,n1x)                        7d1s21
                        iv1=iv1p-1                                      7d1s21
                        ivv=iv1p+nvirt(ksbv1)*iv1                       7d1s21
                        itri=((iv2*(iv2-1))/2)+iv1                      7d1s21
                        irec=iv1+nvirt(ksbv1)*iv2                       7d1s21
                        irow=itri+isw*(irec-itri)                       7d1s21
                        bc(kgg+irow)=bc(kgg+irow)-bc(jint1+ivv)         7d1s21
     $                        *bc(kvd+irow)
                       end do                                           7d1s21
                      end do                                            7d1s21
                      do iv2p=n2n,n2x                                   7d1s21
                       iv2=iv2p-1                                       7d1s21
                       ntop=(iv2+isw*(nvirt(ksbv1)-iv2))-1              7d1s21
                       ivv=iv2p+nvirt(ksbv2)*iv2                        7d1s21
                       xint=bc(jint2+ivv)
                       do iv1=0,ntop                                    7d1s21
                        itri=((iv2*(iv2-1))/2)+iv1                      7d1s21
                        irec=iv1+nvirt(ksbv1)*iv2                       7d1s21
                        irow=itri+isw*(irec-itri)                       7d1s21
                        bc(kgg+irow)=bc(kgg+irow)-xint*bc(kvd+irow)     7d1s21
                       end do                                           7d1s21
                      end do                                            7d1s21
                      kgg=kgg+kvv                                       7d1s21
                      kvd=kvd+kvv                                       7d1s21
                     end do                                             7d1s21
                    end do                                              7d1s21
                   end if                                               7d1s21
                  end do                                                7d1s21
                 end if
                end if                                                  7d1s21
               end do                                                   7d1s21
               ibcoff=ivdtmp                                            7d1s21
              end if                                                    7d1s21
c     !jf2
             end do                                                      12d14s20
             call second(timex)
             call dws_gsumf(bc(idenstrt),nncol)                         6d23s21
             ibctop=ibcoff                                              6d23s21
             mdenhc=0                                                   6d23s21
             xxx=bc(mdenh)                                              6d23s21
             ibc(mdenh)=nint(xxx)                                       6d23s21
             if(ibc(mdenh).ne.0)then                                    6d23s21
              mdenhc=ibcoff                                             6d23s21
              ibcoff=mdenhc+mrow                                        6d29s21
              call enough('hcdduc. 11',bc,ibc)
              do kkk=0,mrow-1                                           6d29s21
               bc(mdenhc+kkk)=bc(idenh+kkk)                             6d23s21
              end do                                                    6d23s21
             end if                                                     6d23s21
             do isa=1,nsymb                                             6d22s21
              ncold=(irefo(isa)*(irefo(isa)+1))/2                         6d22s21
              if(ncold.gt.0)then
              end if
              isc=multh(isa,ijsb)                                       6d23s21
              ncol=irefo(isa)*irefo(isc)                                6d22s21
              mdenk(isa)=0                                              7d1s21
              if(ncol.gt.0)then
               idenkc(isa)=ibcoff                                       7d1s21
               jdenkc=ibcoff                                            6d23s21
               kkkk=0                                                   6d23s21
               do kkk=0,ncol-1                                          6d22s21
                xxx=bc(ndenk(isa)+kkk)                                  6d22s21
                ibc(ndenk(isa)+kkk)=nint(xxx)                           6d22s21
                if(ibc(ndenk(isa)+kkk).ne.0)then                        6d23s21
                 ibc(ndenk(isa)+kkkk)=kkk                               6d23s21
                 iad=idenk(isa,1)+nmrowk*kkk                            7d2s21
                 ibcoff=ibcoff+nmrowk                                   7d2s21
                 call enough('hcdduc. 12',bc,ibc)
                 do i=0,nmrowk-1                                        7d2s21
                  bc(jdenkc+i)=bc(iad+i)                                6d23s21
                 end do                                                 6d23s21
                 jdenkc=jdenkc+nmrowk                                   7d2s21
                 iad=idenk(isa,2)+nmrowk*kkk                            7d2s21
                 ibcoff=ibcoff+nmrowk                                   7d2s21
                 call enough('hcdduc. 13',bc,ibc)
                 do i=0,nmrowk-1                                        7d2s21
                  bc(jdenkc+i)=bc(iad+i)                                6d23s21
                 end do                                                 6d23s21
                 jdenkc=jdenkc+nmrowk                                   7d2s21
                 ibcoff=jdenkc                                          7d1s21
                 kkkk=kkkk+1                                            6d23s21
                end if                                                  6d23s21
               end do                                                    6d22s21
               ibcoff=jdenkc                                            6d23s21
               mdenk(isa)=kkkk                                          6d23s21
              end if
              if(isa.le.isc)then                                        6d29s21
               mdenj(isa)=0                                             6d29s21
               if(isc.eq.isa)then                                       6d22s21
                ncolu=ncold                                             6d22s21
               else                                                     6d22s21
                ncolu=ncol                                              6d22s21
               end if                                                   6d22s21
               if(ncolu.gt.0)then
                idenjc(isa)=ibcoff                                      6d23s21
                jdenjc=ibcoff                                           6d23s21
                kkkk=0                                                  6d23s21
                do kkk=0,ncolu-1                                         6d22s21
                 xxx=bc(ndenj(isa)+kkk)                                 6d23s21
                 ibc(ndenj(isa)+kkk)=nint(xxx)                          6d23s21
                 if(ibc(ndenj(isa)+kkk).ne.0)then                       6d23s21
                  ibc(ndenj(isa)+kkkk)=kkk                              6d23s21
                  iad=idenj(isa)+nmrow*kkk                              6d30s21
                  ibcoff=ibcoff+nmrow                                   6d30s21
                  call enough('hcdduc. 14',bc,ibc)
                  do i=0,nmrow-1                                        6d30s21
                   bc(jdenjc+i)=bc(iad+i)                               6d23s21
                  end do                                                6d23s21
                  jdenjc=jdenjc+nmrow                                   6d30s21
                  kkkk=kkkk+1                                           6d23s21
                 end if                                                 6d23s21
                end do                                                   6d22s21
                ibcoff=jdenjc                                           6d23s21
                mdenj(isa)=kkkk                                         6d23s21
               end if
              end if                                                    6d22s21
             end do                                                     6d22s21
             nwds=ibcoff-ibctop                                         6d23s21
             if(nwds.gt.0)then                                          6d23s21
             call second(timex)
              call dws_gsumf(bc(ibctop),nwds)                           6d23s21
              if(itransgg.eq.0)then                                      2d12s21
               call ddi_done(ibc(iacc),nacc)                             2d12s21
               itransgg=1                                                2d12s21
               do i=igg,igg+mrow*nff2(ncloip,isb,1)-1                    6d22s21
                bc(i)=0d0                                                2d12s21
               end do                                                    2d12s21
              end if                                                     2d12s21
              ioff=igg+mrow*(if2-1)                                     6d23s21
             call second(timex)
              do isbv1=1,nsymb                                                1d12s21
               isbv2=multh(isbv1,isbv12)                                      1d12s21
               if(isbv2.ge.isbv1)then                                         1d12s21
                if(isbv12.eq.1)then                                           1d12s21
                 isw=0                                                        1d12s21
                 mvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        1d12s21
                 ioffp=ioff+nvirt(isbv1)*nrootu*ncsf2(1,iarg)           6d23s21
                else                                                          1d12s21
                 isw=1                                                        1d12s21
                 mvv=nvirt(isbv1)*nvirt(isbv2)                                1d12s21
                 ioffp=ioff
                end if                                                        1d12s21
                mmvv=mvv*nrootu                                               1d12s21
                joff=0                                                  6d23s21
                joffk=0                                                 7d2s21
                nqvisv=0
                nqvnotv=0
                do jsbv1=1,nsymb                                              1d12s21
                 jsbv2=multh(jsbv1,jsbv12)                                    1d12s21
                 if(jsbv2.ge.jsbv1)then                                       1d12s21
                  if(jsbv12.eq.1)then                                         1d12s21
                   jsw=0                                                      1d12s21
                   nvv=(nvirt(jsbv1)*(nvirt(jsbv2)-1))/2                      1d12s21
                   joffp=joff+nvirt(jsbv1)*nrootu*ncsf2(1,iarg)         6d24s21
                   joffkp=joffk+nvirt(jsbv1)*nrootu*ncsf(iarg)          7d2s21
                   nqvisv=nqvisv+nvirt(jsbv1)
                  else                                                        1d12s21
                   jsw=1                                                      1d12s21
                   nvv=nvirt(jsbv1)*nvirt(jsbv2)                              1d12s21
                   joffp=joff                                                 1d13s21
                   joffkp=joffk                                         7d2s21
                  end if                                                      1d12s21
                  nqvnotv=nqvnotv+nvv
                  nnvv=nvv*nrootu                                             1d12s21
c     isbv             lr    jsbv     match  br   bl     common=v',l=v,r=v"
c     Gvvrj =djiab (ab|vv") Vvv"ri     3&4
c     Gvv'rj=djiab (ab|vv") Vv'v"ri     4     2    n vv'
c     Gvv'rj=djiab (ab|vv") Vv"v'ri     2     1    n vv'
c     Gv'vrj=djiab (ab|vv") Vv'v"ri     3     2    t vv'
c     Gv'vrj=djiab (ab|vv") Vv"v'ri     1     1    t vv'
c
c     for diagonal: Gvv'rj=djiab (ab|vv) Vvv'ri match=2
c                   Gv'vrj=djiab (ab|vv) Vv'vri match=3
c                   Gvv rj=djiab (ab|vv ) Vvv ri match = 1 to 4.
                  if(jsbv2.eq.isbv1.or.jsbv2.eq.isbv2.or.jsbv1.eq.isbv1
     $                 .or.jsbv1.eq.isbv2)then                          6d23s21
                   do imatch=1,4                                              1d13s21
                    kase=0                                                    1d14s21
                    if(jsbv2.eq.isbv1.and.imatch.eq.1)then                    1d13s21
                     isbvc=jsbv2                                               1d13s21
                     isbr=jsbv1                                                1d13s21
                     isbl=isbv2                                                1d13s21
                     itf=0                                                     1d13s21
                     jtf=0                                                     1d13s21
                     ibl=2                                              6d30s21
                     ibr=1                                                     1d13s21
                     kase=1
                    end if                                                    1d13s21
                    if(jsbv2.eq.isbv2.and.imatch.eq.2)then                    1d13s21
                     isbvc=jsbv2                                               1d13s21
                     isbr=jsbv1                                                1d13s21
                     isbl=isbv1                                                1d13s21
                     itf=1                                                     1d13s21
                     jtf=0                                                    1d15s21
                     ibl=1                                              6d30s21
                     ibr=1                                                     1d13s21
                     kase=2
                    end if                                                    1d13s21
                    if(jsbv1.eq.isbv1.and.imatch.eq.3)then                    1d14s21
                     isbvc=jsbv1                                               1d13s21
                     isbr=jsbv2                                                1d13s21
                     isbl=isbv2                                                1d13s21
                     jtf=1                                                     1d13s21
                     itf=0                                                    1d15s21
                     ibl=2                                              6d30s21
                     ibr=2                                                     1d13s21
                     kase=3
                    end if                                                    1d13s21
                    if(jsbv1.eq.isbv2.and.imatch.eq.4)then                    1d13s21
                     isbvc=jsbv1                                               1d13s21
                     isbr=jsbv2                                                1d13s21
                     isbl=isbv1                                                1d13s21
                     itf=1                                                     1d13s21
                     jtf=1                                                     1d13s21
                     ibl=1                                                     1d13s21
                     ibr=2                                                     1d13s21
                     kase=4                                                    1d13s21
                    end if                                                     1d13s21
                    iioffp=ioffp                                               1d12s21
                    jjoffp=joffp                                            1d14s21
                    if(kase.ne.0)then                                   6d23s21
                     call ilimts(nvirt(isbl),nvirt(isbr),mynprocg,      6d30s21
     $                       mynowprog,il,ih,i1s,i1e,i2s,i2e)                  1d12s21
                     nhere=ih+1-il                                            1d12s21
                     nhere2=i2e+1-i2s                                         1d12s21
                     if(min(nhere,nvirt(isbvc)).gt.0)then               6d29s21
                      nr=nvirt(isbvc)*nrootu*ncsf(iarg)                 6d23s21
                      itmpp=ibcoff                                      6d23s21
                      ibcoff=itmpp+nr*nvirt(isbl)                       6d30s21
                      call enough('hcdduc. 15',bc,ibc)
                      fact=0d0                                          6d23s21
                      do lsb=1,nsymb                                    6d23s21
                       lsc=multh(lsb,ijsb)                              6d29s21
                       if(mdenj(lsb).ne.0.and.lsb.le.lsc)then           6d29s21
                        nm=nhere2*mdenj(lsb)                            6d30s21
                        itmpi=ibcoff                                    6d23s21
                        ivtmp=itmpi+nvirt(isbl)*nm                      6d30s21
                        ibcoff=ivtmp+nr*nm                              6d30s21
                        call enough('hcdduc. 16',bc,ibc)
                        do kkk=itmpi,ibcoff-1                           6d23s21
                         bc(kkk)=0d0                                    6d23s21
                        end do                                          6d23s21
                        i2eu=inv(1,lsc,lsb,isbl)                        6d30s21
                        icase=inv(2,lsc,lsb,isbl)                       6d30s21
                        do j=0,mdenj(lsb)-1                             6d23s21
                         jj=ibc(ndenj(lsb)+j)                           6d23s21
                         jint=jmats(i2eu)+nhere*jj                      6d24s21
                         ivd=idenjc(lsb)+nmrow*j+joff                   6d30s21
                         ivdp=idenjc(lsb)+nmrow*j+joffp                 6d30s21
                         i10=i1s                                         6d23s21
                         i1n=nvirt(isbl)                                6d30s21
                         do i2=i2s,i2e                                   6d23s21
                          if(i2.eq.i2e)i1n=i1e                           6d23s21
                          do i1=i10,i1n                                 6d23s21
                           iv1=i1-1                                      6d23s21
                           iad=itmpi+iv1+nvirt(isbl)*(i2-i2s+nhere2*j)  6d30s21
                           bc(iad)=bc(jint)
                           jint=jint+1                                  6d23s21
                          end do                                         6d23s21
                          i10=1                                          6d23s21
                         end do                                          6d23s21
                         if(ibr.eq.1)then                               6d23s21
                          do i=0,ncsf(iarg)-1                           6d23s21
                           do ir=0,nrootm                               6d23s21
                            do iv2=0,nvirt(isbvc)-1                     6d30s21
                             ntop=iv2+jsw*(nvirt(isbr)-iv2)             6d30s21
                             do iv1p=i2s,min(i2e,ntop)                  6d30s21
                              iv1=iv1p-1                                6d30s21
                              itri=((iv2*(iv2-1))/2)+iv1                6d23s21
                              irec=iv1+nvirt(isbr)*iv2                  6d30s21
                              iad=ivdp+itri+jsw*(irec-itri)+nvv*(ir     6d23s21
     $                             +nrootu*i)                           6d23s21
                              jad=ivtmp+iv1p-i2s+nhere2*(j+mdenj(lsb)*  6d30s21
     $                             (ir+nrootu*(i+ncsf(iarg)*iv2)))      6d30s21
                              bc(jad)=bc(iad)                           6d23s21
                             end do                                     6d23s21
                            end do                                      6d23s21
                           end do                                       6d23s21
                          end do                                        6d23s21
                         else                                           6d23s21
                          tfr=1d0                                       6d30s21
                          do i=0,ncsf(iarg)-1                           6d23s21
                           if(i.eq.ncsf2(1,iarg))tfr=-1d0               6d30s21
                           do ir=0,nrootm                               6d23s21
                            do iv2p=i2s,i2e                             6d30s21
                             iv2=iv2p-1                                 6d30s21
                             ntop=iv2+jsw*(nvirt(isbvc)-iv2)            6d30s21
                             do iv1=0,ntop-1                            6d23s21
                              itri=((iv2*(iv2-1))/2)+iv1                6d23s21
                              irec=iv1+nvirt(isbvc)*iv2                 6d30s21
                              iad=ivdp+itri+jsw*(irec-itri)+nvv*(ir      6d23s21
     $                             +nrootu*i)                           6d23s21
                              jad=ivtmp+iv2p-i2s+nhere2*(j+mdenj(lsb)*  6d30s21
     $                             (ir+nrootu*(i+ncsf(iarg)*iv1)))      6d30s21
                              bc(jad)=bc(iad)*tfr                       6d30s21
                             end do                                     6d23s21
                            end do                                      6d23s21
                           end do                                       6d23s21
                          end do                                        6d23s21
                         end if                                         6d23s21
                         if(jsbv12.eq.1)then                            6d23s21
                          do i=0,ncsf2(1,iarg)-1                        6d23s21
                           do ir=0,nrootm                               6d23s21
                            do iv2p=i2s,i2e                             6d30s21
                             iv2=iv2p-1                                 6d30s21
                             iad=ivd+iv2+nvirt(jsbv2)*(ir+nrootu*i)     6d23s21
                             jad=ivtmp+iv2p-i2s+nhere2*(j+mdenj(lsb)*(  6d30s21
     $                            ir+nrootu*(i+ncsf(iarg)*iv2)))        6d30s21
                             bc(jad)=bc(iad)*srh                        6d30s21
                            end do                                      6d23s21
                           end do                                       6d23s21
                          end do                                        6d23s21
                         end if                                         6d23s21
                        end do                                          6d23s21
                        if(min(nvirt(isbl),nr,nm).gt.0)then             8d8s24
                         call dgemm('n','n',nvirt(isbl),nr,nm,1d0,       6d30s21
     $                       bc(itmpi),nvirt(isbl),bc(ivtmp),nm,fact,   6d30s21
     $                       bc(itmpp),nvirt(isbl),                     6d30s21
     d' hcdduc.  1')
                         fact=1d0                                        6d23s21
                        end if                                          8d8s24
                        ibcoff=itmpi                                    6d29s21
                       end if                                           6d23s21
                       if(mdenk(lsb).ne.0)then                          7d1s21
                        nm=nhere2*mdenk(lsb)                            7d1s21
                        itmpi=ibcoff                                    6d23s21
                        ivtmp=itmpi+nvirt(isbl)*nm                      6d30s21
                        ibcoff=ivtmp+nr*nm                              6d30s21
                        call enough('hcdduc. 17',bc,ibc)
                        do kkk=itmpi,ibcoff-1                           6d23s21
                         bc(kkk)=0d0                                    6d23s21
                        end do                                          6d23s21
                        i2eu=invk1(1,lsb,lsc,isbl,1)                        6d30s21
                        icase=invk1(2,lsb,lsc,isbl,1)                       6d30s21
                        do j=0,mdenk(lsb)-1                             6d23s21
                         jj=ibc(ndenk(lsb)+j)                           6d23s21
                         kint=kmats(i2eu)+nhere*jj                      7d1s21
                         i10=i1s                                         6d23s21
                         i1n=nvirt(isbl)                                6d30s21
                         do i2=i2s,i2e                                   6d23s21
                          if(i2.eq.i2e)i1n=i1e                           6d23s21
                          do i1=i10,i1n                                 6d23s21
                           iv1=i1-1                                      6d23s21
                           iad=itmpi+iv1+nvirt(isbl)*(i2-i2s+nhere2*j)  6d30s21
                           bc(iad)=bc(kint)
                           kint=kint+1                                  6d23s21
                          end do                                         6d23s21
                          i10=1                                          6d23s21
                         end do                                          6d23s21
                         if(ibr.eq.1)then                               6d23s21
                          ivd=idenkc(lsb)+nmrowk2*j+joffk               7d2s21
                          ivdp=idenkc(lsb)+nmrowk2*j+joffkp             7d2s21
                          do i=0,ncsf(iarg)-1                           6d23s21
                           do ir=0,nrootm                               6d23s21
                            do iv2=0,nvirt(isbvc)-1                     6d30s21
                             ntop=iv2+jsw*(nvirt(isbr)-iv2)             6d30s21
                             do iv1p=i2s,min(i2e,ntop)                  6d30s21
                              iv1=iv1p-1                                6d30s21
                              itri=((iv2*(iv2-1))/2)+iv1                6d23s21
                              irec=iv1+nvirt(isbr)*iv2                  6d30s21
                              iad=ivdp+itri+jsw*(irec-itri)+nvv*(ir     6d23s21
     $                             +nrootu*i)                           6d23s21
                              jad=ivtmp+iv1p-i2s+nhere2*(j+mdenk(lsb)*  7d1s21
     $                             (ir+nrootu*(i+ncsf(iarg)*iv2)))      6d30s21
                              bc(jad)=bc(iad)                           6d23s21
                             end do                                     6d23s21
                            end do                                      6d23s21
                           end do                                       6d23s21
                          end do                                        6d23s21
                         else                                           6d23s21
                          ivd=idenkc(lsb)+nmrowk2*j+joffk+nmrowk        7d2s21
                          ivdp=idenkc(lsb)+nmrowk2*j+joffkp+nmrowk      7d2s21
                          tf=1d0                                        7d2s21
                          do i=0,ncsf(iarg)-1                           6d23s21
                           if(i.eq.ncsf2(1,iarg))tf=-1d0                7d2s21
                           do ir=0,nrootm                               6d23s21
                            do iv2p=i2s,i2e                             6d30s21
                             iv2=iv2p-1                                 6d30s21
                             ntop=iv2+jsw*(nvirt(isbvc)-iv2)            6d30s21
                             do iv1=0,ntop-1                            6d23s21
                              itri=((iv2*(iv2-1))/2)+iv1                6d23s21
                              irec=iv1+nvirt(isbvc)*iv2                 6d30s21
                              iad=ivdp+itri+jsw*(irec-itri)+nvv*(ir      6d23s21
     $                             +nrootu*i)                           6d23s21
                              jad=ivtmp+iv2p-i2s+nhere2*(j+mdenk(lsb)*  7d1s21
     $                             (ir+nrootu*(i+ncsf(iarg)*iv1)))      6d30s21
                              bc(jad)=bc(iad)
                             end do                                     6d23s21
                            end do                                      6d23s21
                           end do                                       6d23s21
                          end do                                        6d23s21
                         end if                                         6d23s21
                         if(jsbv12.eq.1)then                            6d23s21
                          do i=0,ncsf(iarg)-1                           7d2s21
                           do ir=0,nrootm                               6d23s21
                            do iv2p=i2s,i2e                             6d30s21
                             iv2=iv2p-1                                 6d30s21
                             iad=ivd+iv2+nvirt(jsbv2)*(ir+nrootu*i)     6d23s21
                             jad=ivtmp+iv2p-i2s+nhere2*(j+mdenk(lsb)*(  7d1s21
     $                            ir+nrootu*(i+ncsf(iarg)*iv2)))        6d30s21
                             bc(jad)=bc(iad)*srh                        6d30s21
                            end do                                      6d23s21
                           end do                                       6d23s21
                          end do                                        6d23s21
                         end if                                         6d23s21
                        end do                                          6d23s21
                        if(min(nvirt(isbl),nr,nm).gt.0)then             8d8s24
                         call dgemm('n','n',nvirt(isbl),nr,nm,1d0,       6d30s21
     $                       bc(itmpi),nvirt(isbl),bc(ivtmp),nm,fact,   6d30s21
     $                       bc(itmpp),nvirt(isbl),                     6d30s21
     d' hcdduc.  2')
                         fact=1d0                                        6d23s21
                        end if                                          8d8s24
                        ibcoff=itmpi                                    6d29s21
                       end if                                           6d23s21
                      end do                                            6d23s21
                      if(fact.ne.0d0)then                               6d29s21
                       if(ibl.eq.1)then                                  6d23s21
                        do iv2=0,nvirt(isbv2)-1                         6d30s21
                         ntop=(iv2+isw*(nvirt(isbv1)-iv2))-1            6d30s21
                         do i=0,ncsf(iarg)-1                             6d23s21
                          do ir=0,nrootm                                 6d23s21
                           do iv1=0,ntop                                6d30s21
                            itri=((iv2*(iv2-1))/2)+iv1                   6d23s21
                            irec=iv1+nvirt(isbv1)*iv2                    6d23s21
                            iad1=ioffp+itri+isw*(irec-itri)             6d29s21
     $                           +mvv*(ir+nrootu*i)                     6d29s21
                            iad2=itmpp+iv1+nvirt(isbv1)*(ir+nrootu*     6d30s21
     $                          (i+ncsf(iarg)*iv2))                     6d30s21
                            bc(iad1)=bc(iad1)+bc(iad2)                   6d24s21
                           end do                                        6d23s21
                          end do                                         6d23s21
                         end do                                          6d23s21
                        end do
                       else                                              6d23s21
                        do iv1=0,nvirt(isbv1)-1                         6d30s21
                         nbot=(iv1+1)*(1-isw)                           6d30s21
                         tfl=1d0                                        6d30s21
                         do i=0,ncsf(iarg)-1                             6d23s21
                          if(i.eq.ncsf2(1,iarg))tfl=-1d0                6d30s21
                          do ir=0,nrootm                                 6d23s21
                           do iv2=nbot,nvirt(isbv2)-1                   6d30s21
                            itri=((iv2*(iv2-1))/2)+iv1                   6d23s21
                            irec=iv1+nvirt(isbv1)*iv2                    6d23s21
                            iad1=ioffp+itri+isw*(irec-itri)             6d29s21
     $                           +mvv*(ir+nrootu*i)                     6d29s21
                            iad2=itmpp+iv2+nvirt(isbv2)*(ir+nrootu*     6d29s21
     $                          (i+ncsf(iarg)*iv1))                     6d30s21
                            bc(iad1)=bc(iad1)+bc(iad2)*tfl              6d30s21
                           end do                                        6d23s21
                          end do                                         6d23s21
                         end do                                          6d23s21
                        end do
                       end if                                            6d23s21
                       if(isbv12.eq.1)then                               6d23s21
                        do iv2=0,nvirt(isbv2)-1                         6d30s21
                         do i=0,ncsf2(1,iarg)-1                         6d30s21
                          do ir=0,nrootm                                 6d23s21
                           iad1=ioff+iv2+nvirt(isbv2)*(ir+nrootu*i)      6d23s21
                           iad2=itmpp+iv2+nvirt(isbv2)*(ir+nrootu*       6d23s21
     $                          (i+ncsf(iarg)*iv2))                     6d30s21
                           bc(iad1)=bc(iad1)+bc(iad2)*srh                6d24s21
                          end do                                         6d23s21
                         end do                                          6d23s21
                        end do
                       end if                                            6d23s21
                      end if                                            6d29s21
                      ibcoff=itmpp                                      6d29s21
                     end if                                             6d23s21
                    end if                                              6d23s21
                   end do                                               6d23s21
                  end if                                                6d23s21
                  joff=joffp+nvv*nrootu*ncsf(iarg)                      6d24s21
                  joffk=joffkp+nvv*nrootu*ncsf(iarg)                    7d2s21
                 end if                                                 6d23s21
                end do                                                  6d23s21
                ioff=ioffp+mvv*nrootu*ncsf(iarg)                        6d24s21
               end if                                                   6d23s21
              end do                                                    6d23s21
              call second(timex)
              if(ijsb.eq.1.and.mdenhc.gt.0)then                         6d26s21
               ioff=0                                                   6d24s21
               jgg=igg+mrow*(if2-1)                                     6d24s21
               do isbv1=1,nsymb                                         6d24s21
                isbv2=multh(isbv1,isbv12)                               6d24s21
                if(isbv2.ge.isbv1)then                                  6d24s21
                 if(isbv1.eq.isbv2)then                                 6d24s21
                  ioffp=ioff+nvirt(isbv1)*nrootu*ncsf2(1,iarg)          6d24s21
c     nvv,nroot,ncsf,iff2                                               6d24s21
                  nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                 6d24s21
                  isw=0                                                 6d24s21
                  do i=0,ncsf2(1,iarg)-1                                6d28s21
                   do ir=0,nrootm                                       6d28s21
                    iad1=mdenhc+ioff+nvirt(isbv1)*(ir+nrootu*i)         6d28s21
                    lgg=jgg+ioff+nvirt(isbv1)*(ir+nrootu*i)             6d29s21
                    do iv2=1,nvirt(isbv1)-1                             6d28s21
                     if(mod(idoit,mynprocg).eq.mynowprog)then           6d28s21
                      ih0=ih0av(isbv1)+irefo(isbv1)+nh0av(isbv1)        6d28s21
     $                   *(iv2+irefo(isbv1))                            6d28s21
                      kgg=jgg+ioffp+((iv2*(iv2-1))/2)+nvv*(ir+nrootu*i) 6d28s21
                      kad1=mdenhc+ioffp+((iv2*(iv2-1))/2)               6d29s21
     $                     +nvv*(ir+nrootu*i)                           6d29s21
                      do iv1=0,iv2-1                                    6d28s21
                       bc(kgg+iv1)=bc(kgg+iv1)+bc(ih0+iv1)*(bc(iad1+iv1)6d28s21
     $                      +bc(iad1+iv2))*sr2
                       bc(lgg+iv1)=bc(lgg+iv1)+bc(ih0+iv1)*bc(kad1+iv1) 6d29s21
     $                      *sr2
                       bc(lgg+iv2)=bc(lgg+iv2)+bc(ih0+iv1)*bc(kad1+iv1) 6d29s21
     $                      *sr2
                      end do                                            6d28s21
                     end if                                             6d28s21
                     idoit=idoit+1                                      6d28s21
                    end do                                              6d28s21
                   end do                                               6d28s21
                  end do                                                6d28s21
                  tf=1d0                                                6d28s21
                  do i=0,ncsf(iarg)-1                                    6d28s21
                   if(i.eq.ncsf2(1,iarg))tf=-1d0                        6d28s21
                   do ir=0,nrootm                                        6d28s21
                    iad1=mdenhc+ioffp+nvv*(ir+nrootu*i)                  6d28s21
                    lgg=jgg+ioffp+nvv*(ir+nrootu*i)                      6d28s21
                    do iv2=1,nvirt(isbv2)-1                             6d29s21
                     if(mod(idoit,mynprocg).eq.mynowprog)then           6d28s21
                      iv2m=iv2-1                                        6d28s21
                      iv2p=iv2+1                                        6d28s21
                      lad1=iad1+iv2                                     6d28s21
                      kgg=lgg+((iv2*(iv2-1))/2)                         6d28s21
                      lh0=ih0av(isbv1)+irefo(isbv1)+nh0av(isbv1)        6d28s21
     $                      *(iv2+irefo(isbv1))                         6d28s21
                      do iv1=0,iv2m                                     6d28s21
                       iv1m=iv1-1                                       6d28s21
                       ih0=ih0av(isbv1)+irefo(isbv1)+nh0av(isbv1)       6d28s21
     $                      *(iv1+irefo(isbv1))                         6d28s21
c     Gvv'ri=dij Hvv" Vv'v"rj = Hvv"  Dv'v"ri
                       do iv=iv2p,nvirt(isbv1)-1                        6d28s21
                        kad1=lad1+((iv*(iv-1))/2)                       6d28s21
                        bc(kgg+iv1)=bc(kgg+iv1)+bc(ih0+iv)*bc(kad1)*tf  6d28s21
                       end do                                           6d28s21
                       kad1=iad1+((iv1*(iv1-1))/2)                      6d28s21
c     Gvv'ri=dij Hv'v" Vv"vrj = Hv'v" Dv"vri                            6d28s21
                       do iv=0,iv1m                                     6d28s21
                        bc(kgg+iv1)=bc(kgg+iv1)                         6d28s21
     $                       +bc(lh0+iv)*bc(kad1+iv)*tf
                       end do                                           6d28s21
                      end do                                            6d28s21
                     end if                                             6d28s21
                     idoit=idoit+1                                      6d28s21
                    end do                                              6d28s21
                   end do                                               6d28s21
                  end do                                                6d28s21
                 else                                                   6d24s21
                  ioffp=ioff                                            6d24s21
                  nvv=nvirt(isbv1)*nvirt(isbv2)                         6d24s21
                  isw=1                                                 6d24s21
                 end if                                                 6d24s21
                 do i=0,ncsf(iarg)-1                                    6d28s21
                  do ir=0,nrootm                                        6d28s21
                   iad1=mdenhc+ioffp+nvv*(ir+nrootu*i)                  6d28s21
                   lgg=jgg+ioffp+nvv*(ir+nrootu*i)                      6d28s21
                   do iv2=0,nvirt(isbv2)-1                              6d28s21
                    ntop=(iv2+isw*(nvirt(isbv1)-iv2))-1                 6d28s21
                    if(mod(idoit,mynprocg).eq.mynowprog)then            6d28s21
                     ivv2a=((iv2*(iv2-1))/2)                            6d28s21
                     ivv2b=nvirt(isbv1)*iv2                             6d28s21
                     ivv2=ivv2a+isw*(ivv2b-ivv2a)                       6d28s21
                     kgg=lgg+ivv2                                       6d28s21
                     kad1=iad1+ivv2                                     6d28s21
                     lh0=ih0av(isbv2)+irefo(isbv2)+nh0av(isbv2)*        6d28s21
     $                   (irefo(isbv2)+iv2)                             6d28s21
                     do iv1=0,ntop                                      6d28s21
                      ih0=ih0av(isbv1)+irefo(isbv1)+nh0av(isbv1)*        6d28s21
     $                   (irefo(isbv1)+iv1)                             6d28s21
c     Gvv'ri=dij Hvv" Vv"v'rj = Hvv"  Dv"v'ri
c      xx              x x
c       x              x x
c                        x
                      do iv=0,ntop                                      6d28s21
                       bc(kgg+iv1)=bc(kgg+iv1)+bc(ih0+iv)*bc(kad1+iv)   6d28s21
                      end do                                            6d28s21
                      bc(kgg+iv1)=bc(kgg+iv1)-bc(ih0+iv1)*bc(kad1+iv1)  6d28s21
                      nbot=(iv1+1)*(1-isw)                              6d28s21
c     Gvv'ri=dij Hv'v" Vvv"rj = Hv'v" Dvv"ri
c      xx               xx
c       x                x
c                        x
                      do iv=nbot,nvirt(isbv2)-1                         6d28s21
                       itri=((iv*(iv-1))/2)+iv1                         6d28s21
                       irec=iv1+nvirt(isbv1)*iv                         6d28s21
                       lad1=iad1+itri+isw*(irec-itri)                   6d28s21
                       bc(kgg+iv1)=bc(kgg+iv1)+bc(lh0+iv)*bc(lad1)      6d28s21
                      end do                                            6d28s21
                      itri=((iv2*(iv2-1))/2)+iv1                        6d28s21
                      irec=iv1+nvirt(isbv1)*iv2                         6d28s21
                      lad1=iad1+itri+isw*(irec-itri)                    6d28s21
                      bc(kgg+iv1)=bc(kgg+iv1)-bc(lh0+iv2)*bc(lad1)      6d28s21
                     end do                                             6d28s21
                    end if                                              6d28s21
                    idoit=idoit+1                                       6d28s21
                   end do                                               6d28s21
                  end do                                                6d28s21
                 end do                                                 6d28s21
                 ioff=ioffp+nvv*nrootu*ncsf(iarg)                       6d24s21
c     !isbv2.ge.isbv1
                end if                                                  6d24s21
c     !isbv1
               end do                                                   6d24s21
c     !ijsb
              end if                                                    6d24s21
c     !nwds
             end if                                                     6d23s21
             ibcoff=ibctop                                              6d25s21
c     !if2
            end do                                                      6d22s21
            ibcoff=idenstrt                                             6d22s21
           end if                                                       6d22s21
c     !nclojp
          end do                                                        6d22s21
c     !jsb
         end do                                                         6d22s21
         i18=1                                                          12d15s20
         i28=mrow                                                       6d22s21
         i38=1                                                          12d15s20
         i48=nff2(ncloip,isb,1)                                         6d24s21
         if(itransgg.ne.0)then                                          2d12s21
c
c     igg is nvv,nroot,ncsf,iff2                                6d21s21
c     reorder it to is root, ncsf,nvv,jff2.
c
          itmp=ibcoff                                                    1d27s21
          ibcoff=itmp+mrow                                              6d22s21
          call enough('hcdduc. 18',bc,ibc)
          do if=0,nff2(ncloip,isb,1)-1                                   6d22s21
           jgg=igg+mrow*if                                              6d22s21
           jtmp=itmp                                                    6d22s21
           do isbv1=1,nsymb                                             6d22s21
            isbv2=multh(isbv1,isbv12)                                   6d22s21
            if(isbv2.ge.isbv1)then                                      6d22s21
             if(isbv1.eq.isbv2)then                                     6d22s21
              do i=0,ncsf2(1,iarg)-1                                    6d22s21
               do ir=0,nrootm
                do iv=0,nvirt(isbv1)-1                                  6d22s21
                 ktmp=jtmp+ir+nrootu*(i+ncsf2(1,iarg)*iv)                 6d22s21
                 bc(ktmp)=bc(jgg+iv)                                    6d22s21
                end do                                                  6d22s21
                jgg=jgg+nvirt(isbv1)                                    6d22s21
               end do                                                   6d22s21
              end do                                                    6d22s21
              jtmp=jtmp+nvirt(isbv1)*nrootu*ncsf2(1,iarg)               6d22s21
              nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                     6d22s21
             else                                                       6d22s21
              nvv=nvirt(isbv1)*nvirt(isbv2)                             6d22s21
             end if                                                     6d22s21
             do i=0,ncsf(iarg)-1                                        6d22s21
              do ir=0,nrootm                                            6d22s21
               do ivv=0,nvv-1                                           6d22s21
                ktmp=jtmp+ir+nrootu*(i+ncsf(iarg)*ivv)                  6d22s21
                bc(ktmp)=bc(jgg+ivv)                                    6d22s21
               end do                                                   6d22s21
               jgg=jgg+nvv                                              6d22s21
              end do                                                    6d22s21
             end do                                                     6d22s21
             jtmp=jtmp+nvv*nrootu*ncsf(iarg)                            6d25s21
            end if                                                      6d22s21
           end do                                                       6d22s21
           jgg=igg+mrow*if                                              6d22s21
           do i=0,mrow-1                                                6d22s21
            bc(jgg+i)=bc(itmp+i)                                        6d22s21
           end do                                                       6d22s21
          end do                                                        6d22s21
          ibcoff=itmp                                                   6d22s21
          nacc=0
          call ddi_iacc(bc,ibc,ihddiag(ncloip,isb,2),i18,i28,i38,i48,   11d15s22
     $         bc(igg),ibc(iacc),nacc)                                  11d15s22
          nwiacc=nwiacc+i48*i28                                         8d10s22
c
c     this shouldn't be neccessary, but it seems to be needed on my
c     laptop under gcc7.
c
          last8(1)=ihddiag(ncloip,isb,2)                                2d17s21
          last8(2)=i48                                                  2d17s21
          itransgg=0                                                     2d12s21
         else                                                           2d12s21
          nacc=0                                                        2d12s21
         end if                                                         2d12s21
        end if                                                          6d22s21
c     !ncloip
       end do                                                           6d22s21
c     !isb
      end do                                                            6d22s21
      call ddi_done(ibc(iacc),nacc)                                     7d1s21
      ibcoff=iacc                                                       2d12s21
      if(lprt)then                                                      1d25s21
       call dws_synca                                                    1d24s21
       ndoub=0                                                          6d14s21
       do isb=1,nsymb                                                   6d14s21
        isbv12=multh(isb,isymmrci)                                      6d14s21
        nvisv=0                                                         6d14s21
        nvnotv=0                                                        6d14s21
        do isbv1=1,nsymb                                                6d14s21
         isbv2=multh(isbv1,isbv12)                                      6d14s21
         if(isbv2.ge.isbv1)then                                         6d14s21
          if(isbv2.eq.isbv1)then                                        6d14s21
           nvisv=nvisv+nvirt(isbv1)                                     6d14s21
           nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        6d14s21
          else                                                          6d14s21
           nvv=nvirt(isbv1)*nvirt(isbv2)                                6d14s21
          end if                                                        6d14s21
          nvnotv=nvnotv+nvv                                             6d14s21
         end if                                                         6d14s21
        end do                                                          6d14s21
        do ii=mdon+1,mdoo+1                                             6d14s21
         if(nff2(ii,isb,1).gt.0)then                                      6d14s21
          iarg=ii-mdon                                                  6d14s21
          ndoub=ndoub+(ncsf2(1,iarg)*nvisv+ncsf(iarg)*nvnotv)           6d14s21
     $         *nff2(ii,isb,1)                                          6d22s21
         end if                                                         6d14s21
        end do                                                          6d14s21
       end do                                                           6d14s21
       igdmaster=ibcoff                                                 6d15s21
       ibcoff=igdmaster+ndoub*nrootu
       jgdmaster=igdmaster
       call enough('hcdduc. 19',bc,ibc)
       nsp=0
       ntp=0
       do isb=1,nsymb                                                    12d7s20
        isbv12=multh(isb,isymmrci)                                         12d7s20
        do l=1,4                                                        6d11s21
         mtmpx(l)=0                                                     6d11s21
        end do                                                          6d11s21
        do ii=mdon+1,mdoo+1                                             6d11s21
         if(nff2(ii,isb,1).gt.0)then                                    6d22s21
          iarg=ii-mdon                                                  6d11s21
          do l=1,4                                                      6d11s21
           mtmpx(l)=mtmpx(l)+ncsf2(l,iarg)*nff2(ii,isb,1)               6d22s21
          end do                                                        6d11s21
         end if                                                         6d11s21
        end do                                                          6d11s21
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
           itmpx(l,isbv1)=jgdmaster                                     6d15s21
           jtmpx(l,isbv1)=jgdmaster                                     6d15s21
           jgdmaster=itmpx(l,isbv1)+nrootu*ntmpx(l,isbv1)                   6d11s21
          end do                                                        6d11s21
          nvnotv=nvnotv+nvv                                             6d8s21
         end if                                                         6d8s21
        end do                                                          6d8s21
        call enough('hcdduc. 20',bc,ibc)
        xnan=-2d0                                                       6d11s21
        do i=ibcsrt,ibcoff-1                                            6d11s21
         bc(i)=xnan                                                     6d11s21
        end do                                                          6d11s21
        do ii=1,mdoo+1                                                   12d7s20
         if(nff2(ii,isb,1).gt.0)then                                      6d8s21
          iarg=ii-mdon                                                    12d7s20
          nrow=nrootu*(nvisv*ncsf2(1,iarg)+nvnotv*ncsf(iarg))
          itmp=ibcoff                                                   6d8s21
          ibcoff=itmp+nrow*nff2(ii,isb,1)                                 6d8s21
          call enough('hcdduc. 21',bc,ibc)
          do i=itmp,ibcoff-1
           bc(i)=xnan
          end do
          i2o=1                                                         6d8s21
          i2c=nrow                                                      6d8s21
          j2o=1                                                         6d8s21
          k2o=nff2(ii,isb,1)                                              6d8s21
          call ddi_get(bc,ibc,ihddiag(ii,isb,2),i2o,i2c,j2o,k2o,        11d15s22
     $         bc(itmp))                                                11d15s22
          jtmp=itmp                                                     6d11s21
          do if2=0,nff2(ii,isb,1)-1                                       6d11s21
           do isbv1=1,nsymb                                              6d11s21
            isbv2=multh(isbv1,isbv12)                                    6d11s21
            if(isbv2.ge.isbv1)then                                       6d11s21
             if(isbv1.eq.isbv2)then                                      6d11s21
              nvvs=(nvirt(isbv1)*(nvirt(isbv1)+1))/2                    6d11s21
              nvvt=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                    6d11s21
              do iv=0,nvirt(isbv1)-1                                    6d11s21
               ivv=((iv*(iv+1))/2)+iv                                   6d11s21
               do i=0,ncsf2(1,iarg)-1                                    6d11s21
                do ir=0,nrootu-1                                         6d11s21
                 iad=jtmpx(1,isbv1)+i+mtmpx(1)*(ivv+nvvs*ir)            6d11s21
                 bc(iad)=bc(jtmp+ir)                                    6d11s21
                 if(abs(bc(iad)).gt.1d-10)then
                  nsp=nsp+1
                 end if
                end do                                                  6d11s21
                jtmp=jtmp+nrootu                                         6d11s21
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
                  do ir=0,nrootu-1                                       6d11s21
                   iad=jtmpx(l,isbv1)+i+mtmpx(l)*(ivv+nvv*ir)           6d11s21
                   bc(iad)=bc(jtmp+ir)                                  6d11s21
                   if(abs(bc(iad)).gt.1d-10)then
                    if(l.eq.1)then
                     nsp=nsp+1
                    else
                     ntp=ntp+1
                    end if
                   end if
                  end do                                                6d11s21
                  jtmp=jtmp+nrootu                                       6d11s21
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
       end do                                                            12d7s20
       write(6,*)('the whole kit and kaboodle'),nsp,ntp
       call prntm2(bc(igdmaster),ndoub,nrootu,ndoub)
       ibcoff=igdmaster                                                 6d22s21
      end if                                                            1d25s21
      return                                                            12d12s20
      end                                                               12d12s20
