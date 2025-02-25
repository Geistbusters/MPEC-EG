c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcdibk(nff2,iff22,ff22,nsymb,mdon,mdoop,               8d6s21
     $     multh,isymmrci,nvirt,ncsf,ncsf2,nec,nfdat,                   8d6s21
     $     irel,ism,irefo,gb,ixw1,ixw2,norb,                            8d6s21
     $     nff0,iff0,vec,ncsfv,nroot,lprt,srh,n2e,ixmt,i2eop,isymop,    8d6s21
     $     phase2,idoubo,nbasdws,izero,bc,ibc)                          11d10s22
      implicit real*8 (a-h,o-z)                                         10d19s20
      integer*8 itvisvc,itvisvo,itvnotvc,itvnotvo,jrefc,jrefo,ipack8,   12d8s20
     $     itestc,itesto,last8(2),iff22(*),gandcc,gandco,gandcb         2d6s23
      integer*1 idorb(64),isorb(64),nab1(2),iorb(64),imap(64),icode(64),12d8s20
     $     nab2(2)                                                      12d8s20
      logical ldebug,lbail,lprt                                         1d25s21
      dimension nff2(mdoop,2,nsymb),multh(8,8),nvirt(*),                8d6s21
     $     ncsf(*),ncsf2(4,*),nfdat(5,4,*),irel(*),ism(*),jvcv(8),      12d8s20
     $     ipack4(2),itmpt(4),ff22(*),ixmt(8,*),i2eop(2,3),isymop(*),   8d6s21
     $     irefo(*),gb(*),idoubo(*),nbasdws(*),                         8d6s21
     $     idenvisv(8),nl(4),npre(4),mpre(4),                           12d8s20
     $     idenvnotvk(4,8,8),nff0(mdoop,*),iff0(*),                     8d6s21
     $     vec(ncsfv,*),nab4(2,3),ioxx(2)                               2d6s23
      equivalence (ipack8,ipack4)                                       11d20s20
      include "common.store"
      common/kmfind/invk1(2,8,8,8,2)                                    6d30s18
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      if(lprt)write(6,*)('Hi, my name is hcdibk'),ibcoff
      if(izero.ne.0)then                                                8d6s21
       ioff=0                                                           8d6s21
       do isb=1,nsymb                                                   8d6s21
        isbv12=multh(isb,isymmrci)                                      8d6s21
        nft=nfdat(2,1,isb)+nfdat(2,2,isb)+nfdat(2,3,isb)+nfdat(2,4,isb) 8d6s21
        do isbv1=1,nsymb                                                8d6s21
         isbv2=multh(isbv1,isbv12)                                      8d6s21
         if(isbv2.ge.isbv1)then                                         8d6s21
          if(isbv1.eq.isbv2)then                                        8d6s21
           nn=nvirt(isbv1)*nfdat(2,1,isb)                               8d6s21
           nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        8d6s21
          else                                                          8d6s21
           nn=0                                                         8d6s21
           nvv=nvirt(isbv1)*nvirt(isbv2)                                8d6s21
          end if                                                        8d6s21
          nn=nn+nvv*nft                                                 8d6s21
          nn=nn*nroot                                                   8d6s21
          do iz=1,nn                                                    8d6s21
           gb(ioff+iz)=0d0                                              8d6s21
          end do                                                        8d6s21
          ioff=ioff+nn                                                  8d6s21
         end if                                                         8d6s21
        end do                                                          8d6s21
       end do                                                           8d6s21
      end if                                                            8d6s21
      last8(1)=-1                                                       6d17s21
      ldebug=.false.                                                    10d19s20
      mdoo=mdoop-1                                                      8d6s21
      norbx=norb+1                                                      11d19s20
      norbxx=norbx+1                                                    11d19s20
      loop=0
      idoit=0                                                           3d8s21
      ibcoffo=ibcoff                                                    11d22s20
      nn=nroot*nfdat(2,1,isymmrci)
      do jsb=1,nsymb                                                    10d19s20
       jvcv(jsb)=nfdat(5,1,jsb)                                         12d8s20
       jsbv12=multh(jsb,isymmrci)                                       12d8s20
       do l=1,4
        nn=nroot*nfdat(2,l,jsb)                                         12d8s20
        do isa=1,nsymb                                                  11d23s20
         isb=multh(isa,jsbv12)                                          12d8s20
         idenvnotvk(l,jsb,isa)=ibcoff                                   12d8s20
         ibcoff=idenvnotvk(l,jsb,isa)+irefo(isa)*irefo(isb)*nn          12d8s20
        end do                                                          11d23s20
       end do                                                           11d23s20
      end do                                                            10d19s20
      call enough('hcdibk.  1',bc,ibc)
      nden=ibcoff-ibcoffo                                               3d8s21
      do i=ibcoffo,ibcoff-1                                             11d22s20
       bc(i)=0d0                                                        11d22s20
      end do                                                            11d22s20
      ibclast=ibcoff                                                    10d17s20
      do nclo2=mdon,mdoo                                                   10d19s20
       nclo2p=nclo2+1                                                         10d19s20
       iarg2=nclo2p-mdon                                                    10d19s20
       nopen2=nec-2*nclo2p
       nopen2p=nopen2+2                                                 12d7s20
       do isb=1,nsymb                                                   10d19s20
        if(nff2(nclo2p,1,isb).gt.0)then                                  8d6s21
         if(isb.eq.isymmrci)then                                         12d7s20
          ntop=nclo2p+3                                                  12d7s20
         else                                                            12d7s20
          ntop=nclo2p+2                                                  12d7s20
         end if                                                          12d7s20
         do nclo0p=max(mdon+1,nclo2p-2),min(mdoop,ntop)                  8d6s21
          if(nff0(nclo0p,1).gt.0)then                                   8d6s21
           nclo0=nclo0p-1                                                 12d7s20
           iarg0=nclo0p-mdon                                               12d7s20
           nopen=nec-nclo0*2                                              12d7s20
           jvcvu=nfdat(5,1,isb)+nff2(nclo2p,2,isb)                       8d6s21
           do if=1,nff2(nclo2p,1,isb)                                      11d20s20
            ipack8=iff22(jvcvu)                                         8d6s21
            nclot=popcnt(ipack4(1))                                        11d20s20
            itvnotvc=ipack4(1)                                            12d7s20
            itvisvc=ibset(itvnotvc,norbx)                                 12d7s20
            itvisvo=ipack4(2)                                             12d7s20
            iarg2p=iarg2+1                                                12d8s20
            itvnotvo=ibset(itvisvo,norbx)                                 12d7s20
            itvnotvo=ibset(itvnotvo,norbxx)                               12d7s20
            nspace=iff22(jvcvu+1)                                       8d6s21
            if(mod(idoit,mynprocg).eq.mynowprog)then                      3d8s21
             ibctmpt=ibcoff                                                12d8s20
             do l=1,4                                                      12d8s20
              nl(l)=iff22(jvcvu+1+l)                                    8d6s21
              if(nl(l).gt.0)then                                           12d8s20
               iad1=jvcvu+iff22(jvcvu+5+l)                              8d6s21
               iad2=iad1+nl(l)                                             3d19s21
               itmpt(l)=ibcoff                                             12d8s20
               ibcoff=itmpt(l)+ncsf2(l,iarg2)*nl(l)                        12d8s20
               call enough('hcdibk.  2',bc,ibc)
               do i=0,nl(l)-1                                                    11d20s20
                do j=0,ncsf2(l,iarg2)-1                                        11d20s20
                 ji=iad2+j+ncsf2(l,iarg2)*i                                    11d20s20
                 ij=itmpt(l)+i+nl(l)*j                                     11d23s20
                 bc(ij)=ff22(ji)                                        8d6s21
                end do                                                        11d20s20
               end do                                                         11d20s20
              end if                                                       12d8s20
             end do                                                        12d8s20
             jff0=nff0(nclo0p,2)                                           12d8s20
             jvs=nff0(nclo0p,3)                                            12d7s20
             do jf=1,nff0(nclo0p,1)                                        12d7s20
              jrefc=iff0(jff0)                                             12d7s20
              ntestc=popcnt(jrefc)
              jff0=jff0+1                                                  12d7s20
              jrefo=iff0(jff0)                                             12d7s20
              ntesto=popcnt(jrefo)
              if(ntestc.ne.nclo0.or.ntesto.ne.nopen)then
               write(6,*)('for internal fcn '),jf
               write(6,*)('ntestc vs. nclo0 '),ntestc,nclo0
               write(6,*)('ntesto vs. nopen '),ntesto,nopen
               call dws_synca
               call dws_finalize
               stop
              end if
              jff0=jff0+1                                                  12d7s20
c
c     v not v
c
              gandcc=ieor(itvnotvc,jrefc)                               2d6s23
              gandco=ieor(itvnotvo,jrefo)                               2d6s23
              gandcb=ior(gandcc,gandco)                                 2d6s23
              ndifb=popcnt(gandcb)                                      2d6s23
              if(ndifb.le.4)then                                        2d6s23
               ndifs=popcnt(gandco)                                     2d6s23
               ndifd=popcnt(gandcc)                                     2d6s23
               nnot=0                                                   2d6s23
               ipssx=0                                                  2d6s23
               if(ndifs.eq.4.and.ndifb.eq.4)then                        2d6s23
                nnot=4                                                  2d6s23
                ioxx(1)=1                                               2d6s23
                ioxx(2)=1                                               2d6s23
                do i=1,norbxx                                           2d6s23
                 if(btest(gandcb,i))then                                2d6s23
                  if((btest(jrefc,i).and.btest(itvnotvo,i)).or.         2d6s23
     $                (btest(jrefo,i).and..not.btest(itvnotvc,i)))then  2d6s23
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
     $        ((btest(itvnotvc,i).and..not.btest(jrefo,i)).or.          2d6s23
     $         (btest(jrefc,i).and..not.btest(itvnotvo,i))))then        2d6s23
                   if(btest(jrefc,i))iswap=1                            2d6s23
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
                 if(btest(itvnotvc,nab4(1,2)).and.                      2d6s23
     $                .not.btest(itvnotvc,nab4(1,1)))nbt=1              2d6s23
                else                                                    2d6s23
                 nbt=0                                                  2d6s23
                 if(btest(jrefc,nab4(2,2)).and.                         2d6s23
     $                .not.btest(jrefc,nab4(2,1)))nbt=1                 2d6s23
                end if                                                  2d6s23
                if(nbt.ne.0)then                                        2d6s23
                 nab4(1,1)=nab4(1,2)                                    2d6s23
                 nab4(2,1)=nab4(2,2)                                    2d6s23
                end if                                                  2d6s23
               else if(ndifs.eq.0.and.ndifd.eq.2)then                   2d6s23
                nnot=3                                                  2d6s23
                do i=1,norbxx                                           2d6s23
                 if(btest(gandcb,i))then                                2d6s23
                  if(btest(itvnotvc,i))then                             2d6s23
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
               else if(nnot.eq.4)then                                   2d6s23
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
                itestc=itvnotvc                                              12d8s20
                itesto=itvnotvo                                              12d8s20
                if(btest(itestc,nab4(1,iu1)))then                          12d8s20
                 itestc=ibclr(itestc,nab4(1,iu1))                          12d8s20
                 itesto=ibset(itesto,nab4(1,iu1))                          12d8s20
                 nopenk=nopen2p+1                                              11d13s20
                 karg=iarg2-1                                               12d8s20
                else if(btest(itesto,nab4(1,iu1)))then                     12d8s20
                 itesto=ibclr(itesto,nab4(1,iu1))                          12d8s20
                 nopenk=nopen2p-1                                              11d13s20
                 karg=iarg2                                                 12d8s20
                end if                                                        11d13s20
                if(btest(itesto,nab4(2,iu2)))then                          12d8s20
                 itestc=ibset(itestc,nab4(2,iu2))                          12d8s20
                 itesto=ibclr(itesto,nab4(2,iu2))                          12d8s20
                 nopenk=nopenk-1                                              11d13s20
                 karg=karg+1                                                  11d13s20
                else                                                          11d13s20
                 itesto=ibset(itesto,nab4(2,iu2))                          12d8s20
                 nopenk=nopenk+1                                              11d13s20
                end if                                                        11d13s20
                call gandc(itvnotvc,itvnotvo,itestc,itesto,nopen2p,     8d6s21
     $               nopenk,iarg2,karg,ncsf,norbxx,ixw1,ixw2,nnot1,nab1,8d6s21
     $               iwpb1,iwpk1,ncsfmid1,bc,ibc)                       11d14s22
                call gandc(itestc,itesto,jrefc,jrefo,nopenk,nopen,          12d8s20
     $         karg,iarg0,ncsf,norbxx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2, 12d8s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
                if(nnot1.eq.2.and.nnot2.eq.2)then                              11d13s20
                 jsa=ism(nab1(2))                                            12d8s20
                 jga=irel(nab1(2))-1                                         12d8s20
                 jsb=ism(nab2(2))                                            12d8s20
                 jgb=irel(nab2(2))-1
                 itmp=ibcoff                                                 12d8s20
                 ibcoff=itmp+ncsf(iarg2)*nroot                               12d8s20
                 call enough('hcdibk.  3',bc,ibc)
                 call genmatn(ncsf(iarg2),ncsf(karg),ncsf(iarg0),            12d8s20
     $            ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,vec(jvs,1), 12d8s20
     $            ncsfv,nroot,bc(itmp),0d0,bc,ibc)                      11d10s22
                 jtmp=itmp                                                 12d9s20
                 do l=1,4                                                    12d8s20
                  if(nl(l).gt.0)then                                         12d8s20
                   jden=idenvnotvk(l,isb,jsb)+nroot*nfdat(2,l,isb)*        12d8s20
     $            (jgb+irefo(jsb)*jga)                                  12d8s20
                   itmp2=ibcoff                                               12d8s20
                   ibcoff=itmp2+nl(l)*nroot                                   12d8s20
                   call enough('hcdibk.  4',bc,ibc)
                   call dgemm('n','n',nl(l),nroot,ncsf2(l,iarg2),1d0,         12d8s20
     $            bc(itmpt(l)),nl(l),bc(jtmp),ncsf(iarg2),0d0,          12d9s20
     $            bc(itmp2),nl(l),                                      12d8s20
     d' hcdibk.  1')
                   iadu=jvcvu+iff22(jvcvu+5+l)                          8d6s21
                   jtmp2=itmp2                                                12d8s20
                   do ir=0,nroot-1                                            12d8s20
                    kden=jden-1+nfdat(2,l,isb)*ir                          12d9s20
                    do i=0,nl(l)-1                                            12d8s20
                     ip=iff22(iadu+i)                                   8d6s21
                     bc(kden+ip)=bc(kden+ip)+bc(jtmp2+i)                      12d8s20
                    end do                                                    12d8s20
                    jtmp2=jtmp2+nl(l)                                         12d8s20
                   end do                                                     12d8s20
                   ibcoff=itmp2                                              12d8s20
                  end if                                                     12d8s20
                  jtmp=jtmp+ncsf2(l,iarg2)                                 12d9s20
                 end do                                                      12d8s20
                 ibcoff=itmp                                                 12d8s20
                else
                 write(6,*)('nnots not both 2!!! '),nnot1,nnot2
                 call dws_synca
                 call dws_finalize
                 stop
                end if                                                     12d8s20
               end do                                                      12d8s20
              end if                                                       12d8s20
              jvs=jvs+ncsf(iarg0)                                          12d8s20
             end do                                                        12d8s20
            end if                                                        3d8s21
            idoit=idoit+1                                                 3d8s21
            jvcvu=jvcvu+nspace                                            3d19s21
           end do
          end if                                                        8d6s21
         end do
        end if                                                          8d6s21
       end do
      end do
      nn=nroot*nfdat(2,1,isymmrci)
      call dws_gsumf(bc(ibcoffo),nden)                                  3d8s21
      do jsb=1,nsymb                                                    10d19s20
       jsbv12=multh(jsb,isymmrci)                                       12d8s20
       do l=1,4
        nn=nroot*nfdat(2,l,jsb)                                         12d8s20
        do isa=1,nsymb                                                  11d23s20
         isb=multh(isa,jsbv12)                                          12d8s20
         if(min(irefo(isa),irefo(isb),nn).gt.0)then
          nrr=irefo(isa)*irefo(isb)                                     12d8s20
          sz=0d0
          do i=0,nrr*nn-1
           sz=sz+bc(idenvnotvk(l,jsb,isa)+i)**2                         12d8s20
          end do
          sz=sqrt(sz/dfloat(nrr*nn))                                    12d8s20
          itmp=ibcoff                                                   12d8s20
          ibcoff=itmp+nrr*nn                                            12d8s20
          call enough('hcdibk.  5',bc,ibc)
          do irr=0,nrr-1                                                 12d8s20
           do ir=0,nroot-1
            do i=0,nfdat(2,l,jsb)-1
             iad1=idenvnotvk(l,jsb,isa)+i+nfdat(2,l,jsb)                12d8s20
     $             *(ir+nroot*irr)                                      12d8s20
             iad2=itmp+irr+nrr*(ir+nroot*i)                              12d8s20
             bc(iad2)=bc(iad1)                                           12d8s20
            end do                                                       12d8s20
           end do                                                        12d8s20
          end do                                                         12d8s20
          nrow=nrr*nroot                                                12d8s20
          do ic=0,nrow*nfdat(2,l,jsb)-1                                 8d6s21
           bc(idenvnotvk(l,jsb,isa)+ic)=bc(itmp+ic)                     8d6s21
           if(bc(itmp+ic).ne.bc(itmp+ic))write(6,*)('itmp+ic '),
     $          bc(itmp+ic),itmp+ic
          end do                                                        8d6s21
          ibcoff=itmp
         end if
        end do                                                          11d23s20
       end do                                                           11d23s20
      end do                                                            10d19s20
      ioff=0                                                            12d8s20
      do jsb=1,nsymb                                                    12d8s20
       jsbv12=multh(jsb,isymmrci)                                       12d8s20
       do jsbv1=1,nsymb                                                 12d8s20
        jsbv2=multh(jsbv1,jsbv12)                                       12d8s20
        if(jsbv2.ge.jsbv1)then                                          12d8s20
         call ilimts(nvirt(jsbv1),nvirt(jsbv2),mynprocg,mynowprog,il,   12d8s20
     $          ih,i1s,i1e,i2s,i2e)                                     10d19s20
         nhere=ih+1-il                                                  12d8s20
         if(jsbv1.eq.jsbv2.and.nvirt(jsbv1).gt.0)then                   12d8s20
          i10=i1s                                                       12d8s20
          i1n=nvirt(jsbv1)                                              12d8s20
          nvv=0                                                         12d8s20
          ivn=i2e+1                                                     12d8s20
          do i2=i2s,i2e                                                 12d8s20
           if(i2.eq.i2e)i1n=i1e                                         12d8s20
           if(i10.le.i2.and.i1n.ge.i2)then                              12d8s20
            nvv=nvv+1                                                   12d8s20
            ivn=min(ivn,i2)                                             12d8s20
           end if                                                       12d8s20
           i10=1                                                        12d8s20
          end do                                                        12d8s20
          nn=nroot*nfdat(2,1,isymmrci)                                  8d6s21
          if(min(nn,nvv).gt.0)then                                      11d17s21
           iputto=ioff+ivn                                               12d8s20
           do isa=1,nsymb                                                12d8s20
            if(irefo(isa).gt.0)then                                     12d8s20
             nrr=irefo(isa)*irefo(isa)                                  12d8s20
             itmpi=ibcoff                                               12d8s20
             ibcoff=itmpi+nvv*nrr                                       12d8s20
             call enough('hcdibk.  6',bc,ibc)
             do iz=itmpi,ibcoff-1                                       8d6s21
              bc(iz)=0d0                                                8d6s21
             end do                                                     8d6s21
             do ii2e=1,n2e                                              8d6s21
              if(multh(isa,jsbv1).eq.isymop(i2eop(1,ii2e)))then         8d6s21
               itmpva=ibcoff                                            8d6s21
               ibcoff=itmpva+irefo(isa)                                 8d6s21
               call enough('hcdibk.  7',bc,ibc)
               i10=i1s                                                    12d8s20
               i1n=nvirt(jsbv1)
               jtmpi=itmpi                                                12d8s20
               do i2=i2s,i2e
                if(i2.eq.i2e)i1n=i1e
                if(i2.ge.i10.and.i2.le.i1n)then                           1d28s21
                 iv=i2-1+idoubo(jsbv1)+irefo(jsbv1)                     8d6s21
                 do ia=0,irefo(isa)                                     8d6s21
                  iap=ia+idoubo(isa)                                    8d6s21
                  iva=ixmt(jsbv1,i2eop(1,ii2e))+iv+nbasdws(jsbv1)*iap   8d6s21
                  bc(itmpva+ia)=bc(iva)                                 8d6s21
                 end do                                                 8d6s21
                 jtmpi=itmpi+i2-ivn                                        12d8s20
                 do i4=0,irefo(isa)-1                                   8d6s21
                  fact=phase2*bc(itmpva+i4)                             8d6s21
                  do i3=0,irefo(isa)-1                                  8d6s21
                   bc(jtmpi+i3*nvv)=bc(jtmpi+i3*nvv)+fact*bc(itmpva+i3) 8d9s21
                  end do                                                8d6s21
                  jtmpi=jtmpi+irefo(isa)*nvv                            8d9s21
                 end do                                                    12d8s20
                end if                                                    1d28s21
                i10=1                                                     12d8s20
               end do                                                     12d8s20
               ibcoff=itmpva                                            8d6s21
              end if                                                    8d6s21
             end do                                                     8d6s21
             call dgemm('n','n',nvv,nn,nrr,srh,bc(itmpi),nvv,           3d8s21
     $             bc(idenvnotvk(1,isymmrci,isa)),nrr,1d0,gb(iputto),   3d8s21
     $             nvirt(jsbv1),                                        12d8s20
     d' hcdibk.  2')
             ibcoff=itmpi                                               12d8s20
            end if                                                      12d8s20
           end do                                                       12d8s20
          end if                                                        12d8s20
          ioff=ioff+nvirt(jsbv1)*nn                                     12d8s20
         end if                                                         12d8s20
         if(jsbv1.eq.jsbv2)then
          nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                         12d8s20
          isw=0                                                         12d8s20
         else                                                           12d8s20
          nvv=nvirt(jsbv1)*nvirt(jsbv2)                                 12d8s20
          isw=1                                                         12d8s20
         end if                                                         12d8s20
         if(nvv.gt.0)then                                               12d8s20
          if(nhere.gt.0)then                                            12d8s20
           do l=1,4                                                     12d8s20
            itmpt(l)=ibcoff                                             12d8s20
            ibcoff=itmpt(l)+nhere*nroot*nfdat(2,l,jsb)                  8d6s21
           end do                                                       12d8s20
           call enough('hcdibk.  8',bc,ibc)
           do i=itmpt(1),ibcoff-1                                       12d8s20
            bc(i)=0d0                                                   12d8s20
           end do                                                       12d8s20
           do isa=1,nsymb                                                12d8s20
            isb=multh(isa,jsbv12)                                       12d8s20
            if(min(irefo(isa),irefo(isb)).gt.0)then                     12d8s20
             nrr=irefo(isa)*irefo(isb)                                  12d8s20
             if(jsbv12.eq.1)then                                        12d8s20
              jsw=0                                                     8d6s21
             else                                                       8d6s21
              jsw=1                                                     8d6s21
             end if                                                     8d6s21
             itmpi=ibcoff                                               12d8s20
             ibcoff=itmpi+nhere*nrr                                     12d8s20
             call enough('hcdibk.  9',bc,ibc)
             do iz=itmpi,ibcoff-1                                       8d6s21
              bc(iz)=0d0                                                8d6s21
             end do                                                     8d6s21
             jtmpi=itmpi                                                8d6s21
             do ii2e=1,n2e                                              8d6s21
              if(multh(isa,jsbv2).eq.isymop(i2eop(1,ii2e)).and.         8d6s21
     $            multh(isb,jsbv1).eq.isymop(i2eop(2,ii2e)))then        8d6s21
               itmp2a=ibcoff                                            8d6s21
               itmp1b=itmp2a+irefo(isa)*(i2e+1-i2s)                     8d6s21
               ibcoff=itmp1b+irefo(isb)*nvirt(jsbv1)                    8d6s21
               call enough('hcdibk. 10',bc,ibc)
               do i2=i2s,i2e                                            8d6s21
                jtmp2a=itmp2a+irefo(isa)*(i2-i2s)                       8d6s21
                iv=i2-1+idoubo(jsbv2)+irefo(jsbv2)                      8d6s21
                do ia=0,irefo(isa)-1                                    8d6s21
                 iaa=ia+idoubo(isa)                                     8d6s21
                 iva=ixmt(jsbv2,i2eop(1,ii2e))+iv+nbasdws(jsbv2)*iaa    8d6s21
                 bc(jtmp2a+ia)=bc(iva)*phase2                           8d6s21
                end do                                                  8d6s21
               end do                                                   8d6s21
               do i1=1,nvirt(jsbv1)                                     8d6s21
                iv=i1-1                                                 8d6s21
                jtmp1b=itmp1b+irefo(isb)*iv                             8d6s21
                iv=iv+idoubo(jsbv1)+irefo(jsbv1)                        8d6s21
                do ib=0,irefo(isb)-1                                    8d6s21
                 ibb=ib+idoubo(isb)                                     8d6s21
                 ivb=ixmt(jsbv1,i2eop(2,ii2e))+iv+nbasdws(jsbv1)*ibb    8d6s21
                 bc(jtmp1b+ib)=bc(ivb)                                  8d6s21
                end do                                                  8d6s21
               end do                                                   8d6s21
               i10=i1s                                                    12d8s20
               i1n=nvirt(jsbv1)                                           12d8s20
               jtmpi=itmpi                                               12d8s20
               do i2=i2s,i2e                                             12d8s20
                i2m=i2-1                                                 12d8s20
                if(i2.eq.i2e)i1n=i1e                                     12d8s20
                jtmp2a=itmp2a+irefo(isa)*(i2-i2s)                       8d6s21
                itop=min(i1n,i2m+jsw*(i1n-i2m))                         8d6s21
                do i1=i10,itop                                          8d6s21
                 jtmp1b=itmp1b+irefo(isb)*(i1-1)                        8d6s21
                 ktmpi=jtmpi                                             12d8s20
                 do i4=0,irefo(isb)-1                                   8d6s21
                  do i3=0,irefo(isa)-1                                  8d6s21
                   bc(ktmpi)=bc(ktmpi)+bc(jtmp1b+i4)*bc(jtmp2a+i3)      8d6s21
                   ktmpi=ktmpi+nhere                                      12d8s20
                  end do                                                8d6s21
                 end do                                                  12d8s20
                 jtmpi=jtmpi+1                                           12d8s20
                end do                                                   12d8s20
                i10=1                                                    12d8s20
               end do                                                    12d8s20
               ibcoff=itmp2a
              end if                                                    8d6s21
             end do                                                     8d6s21
             nvvh=jtmpi-itmpi                                           12d8s20
             kint=itmpi                                                 12d8s20
             if(nvvh.gt.0)then                                          12d8s20
              do l=1,4                                                       12d8s20
               if(nfdat(2,l,jsb).gt.0)then                              8d6s21
                nn=nfdat(2,l,jsb)*nroot                                 8d6s21
                call dgemm('n','n',nvvh,nn,nrr,1d0,bc(kint),nhere,       12d8s20
     $              bc(idenvnotvk(l,jsb,isa)),nrr,1d0,                  12d8s20
     $              bc(itmpt(l)),nhere,                                 12d8s20
     d' hcdibk.  3')
               end if                                                    12d8s20
              end do                                                     12d8s20
             end if                                                      12d8s20
             ibcoff=itmpi                                               8d6s21
            end if                                                      12d8s20
           end do                                                       12d8s20
           joff=ioff                                                    12d8s20
           do l=1,4                                                     12d8s20
            if(nfdat(2,l,jsb).gt.0)then                                 8d6s21
             nn=nfdat(2,l,jsb)*nroot                                    8d6s21
             i10=i1s                                                      12d8s20
             i1n=nvirt(jsbv1)                                             12d8s20
             jtmp=itmpt(l)                                              12d8s20
             do i2=i2s,i2e                                                12d8s20
              i2m=i2-1                                                  12d8s20
              if(i2.eq.i2e)i1n=i1e                                        12d8s20
              if(jsbv12.eq.1)then                                         12d8s20
               itop=min(i2-1,i1n)                                       12d9s20
              else                                                        12d8s20
               itop=i1n                                                   12d8s20
              end if                                                      12d8s20
              do i1=i10,itop                                            12d8s20
               i1m=i1-1                                                 12d8s20
               itri=((i2m*(i2m-1))/2)+i1                                12d8s20
               irec=i1+nvirt(jsbv1)*i2m                                 12d8s20
               itri=itri+isw*(irec-itri)                                12d8s20
               ito=joff+itri                                            12d8s20
               ktmp=jtmp                                                12d8s20
               do i34=1,nn                                              12d8s20
                gb(ito)=gb(ito)+bc(ktmp)                                12d8s20
                ito=ito+nvv                                             12d8s20
                ktmp=ktmp+nhere                                         12d9s20
               end do                                                   12d8s20
               jtmp=jtmp+1                                              12d8s20
              end do                                                    12d8s20
              i10=1                                                     12d8s20
             end do                                                     12d8s20
             joff=joff+nvv*nn                                           12d8s20
            end if                                                      12d8s20
           end do                                                       12d8s20
           ibcoff=itmpt(1)                                              12d8s20
          end if                                                        12d8s20
          do l=1,4                                                      12d8s20
           ioff=ioff+nvv*nfdat(2,l,jsb)*nroot                           8d6s21
          end do                                                        12d8s20
         end if                                                         12d8s20
        end if                                                          12d8s20
       end do                                                           12d8s20
      end do                                                            12d8s20
      call dws_synca
      ibcoff=ibcoffo                                                    11d23s20
      return                                                            10d19s20
      end                                                               10d19s20
