c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcdi(nff2,nsymb,mdon,mdoo,                             12d7s20
     $     multh,isymmrci,nvirt,ncsf,nec,ncsf2,nfdat,                   11d19s20
     $     irel,ism,kmats,irefo,ismult,gb,ixw1,ixw2,norb,               12d7s20
     $     nff0,iff0,vec,ncsfv,nroot,lprt,srh,bc,ibc)                   11d10s22
      implicit real*8 (a-h,o-z)                                         10d19s20
      integer*8 itvisvc,itvisvo,itvnotvc,itvnotvo,jrefc,jrefo,ipack8,   12d8s20
     $     itestc,itesto,last8(2),gandcc,gandco,gandcb                  11d1s22
      integer*1 idorb(64),isorb(64),nab1(2),iorb(64),imap(64),icode(64),12d8s20
     $     nab2(2)                                                      12d8s20
      logical ldebug,lbail,lprt                                         1d25s21
      dimension nff2(mdoo+1,nsymb,*),multh(8,8),nvirt(*),               11d25s20
     $     ncsf(*),ncsf2(4,*),nfdat(5,4,*),irel(*),ism(*),jvcv(8),      12d8s20
     $     ipack4(2),itmpt(4),                                          12d8s20
     $     kmats(*),irefo(*),gb(*),
     $     idenvisv(8),nl(4),npre(4),mpre(4),                           12d8s20
     $     idenvnotvk(4,8,8),nff0(mdoo+1,*),iff0(*),                    12d8s20
     $     vec(ncsfv,*),nab4(2,3),ioxx(2)                               11d1s22
      equivalence (ipack8,ipack4)                                       11d20s20
      include "common.store"
      common/kmfind/invk1(2,8,8,8,2)                                    6d30s18
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      if(lprt)then
       write(6,*)('Hi, my name is hcdi'),ibcoff
       write(6,*)('input vec: ')
       call prntm2(vec,ncsfv,nroot,ncsfv)
      end if
      last8(1)=-1                                                       6d17s21
      ldebug=.false.                                                    10d19s20
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
      call enough('hcdi.  1',bc,ibc)
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
        if(isb.eq.isymmrci)then                                         12d7s20
         ntop=nclo2p+3                                                  12d7s20
        else                                                            12d7s20
         ntop=nclo2p+2                                                  12d7s20
        end if                                                          12d7s20
        do nclo0p=max(mdon+1,nclo2p-2),min(mdoo+1,ntop)                 12d7s20
         nclo0=nclo0p-1                                                 12d7s20
         iarg0=nclo0p-mdon                                               12d7s20
         nopen=nec-nclo0*2                                              12d7s20
         jvcvu=jvcv(isb)
         do if=1,nff2(nclo2p,isb,7)                                      11d20s20
          ipack8=ibc(jvcvu)                                             12d8s20
          nclot=popcnt(ipack4(1))                                        11d20s20
          itvnotvc=ipack4(1)                                            12d7s20
          itvisvc=ibset(itvnotvc,norbx)                                 12d7s20
          itvisvo=ipack4(2)                                             12d7s20
          iarg2p=iarg2+1                                                12d8s20
          itvnotvo=ibset(itvisvo,norbx)                                 12d7s20
          itvnotvo=ibset(itvnotvo,norbxx)                               12d7s20
          if(nclot.ne.nclo2)then                                         11d20s20
           write(6,*)('for isb = '),isb,(' if = '),if,(' nclot = '),
     $          nclot,(' is ne nclo2 = '),nclo2
           write(6,*)ibc(jvcv(isb)),jvcv(isb)
           stop
          end if
          nspace=ibc(jvcvu+1)                                           3d19s21
          if(mod(idoit,mynprocg).eq.mynowprog)then                      3d8s21
          ibctmpt=ibcoff                                                12d8s20
          do l=1,4                                                      12d8s20
           nl(l)=ibc(jvcvu+1+l)                                         3d19s21
           if(nl(l).gt.0)then                                           12d8s20
            iad1=jvcvu+ibc(jvcvu+5+l)                                   3d19s21
            iad2=iad1+nl(l)                                             3d19s21
            itmpt(l)=ibcoff                                             12d8s20
            ibcoff=itmpt(l)+ncsf2(l,iarg2)*nl(l)                        12d8s20
            call enough('hcdi.  2',bc,ibc)
            do i=0,nl(l)-1                                                    11d20s20
             do j=0,ncsf2(l,iarg2)-1                                        11d20s20
              ji=iad2+j+ncsf2(l,iarg2)*i                                    11d20s20
              ij=itmpt(l)+i+nl(l)*j                                     11d23s20
              bc(ij)=bc(ji)                                               11d20s20
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
            stop
           end if
           jff0=jff0+1                                                  12d7s20
c
c     v not v
c
           gandcc=ieor(itvnotvc,jrefc)                                  10d13s22
           gandco=ieor(itvnotvo,jrefo)                                  10d13s22
           gandcb=ior(gandcc,gandco)                                    10d21s22
           ndifb=popcnt(gandcb)                                         10d21s22
           if(ndifb.le.4)then                                           10d21s22
            ndifd=popcnt(gandcc)                                         10d13s22
            ndifs=popcnt(gandco)                                         10d13s22
            nnot=0                                                      10d21s22
            ipssx=0
            if(ndifs.eq.4.and.ndifb.eq.4)then                                           10d14s22
             nnot=4                                                          10d14s22
             ioxx(1)=1                                                     10d17s22
             ioxx(2)=1                                                     10d17s22
             do i=1,norbxx                                                     10d17s22
              if(btest(gandcb,i))then                                         10d14s22
               if((btest(jrefc,i).and.btest(itvnotvo,i)).or.                         10d17s22
     $             (btest(jrefo,i).and..not.btest(itvnotvc,i)))then                   10d14s22
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
             do i=1,norbxx
              if(btest(gandcb,i))then                                         10d14s22
               if(btest(gandcc,i).and.                                        10d14s22
     $             ((btest(itvnotvc,i).and..not.btest(jrefo,i)).or.                 10d14s22
     $             (btest(jrefc,i).and..not.btest(itvnotvo,i))))then                     10d14s22
               if(btest(jrefc,i))iswap=1                                        10d17s22
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
              if(btest(itvnotvc,nab4(1,2)).and.                         10d21s22
     $             .not.btest(itvnotvc,nab4(1,1)))                      10d21s22
     $             nbt=1                                                      10d17s22
             else                                                             10d17s22
              nbt=0                                                           10d17s22
              if(btest(jrefc,nab4(2,2)).and..not.btest(jrefc,nab4(2,1)))
     $             nbt=1                                                      10d17s22
             end if                                                           10d17s22
             if(nbt.ne.0)then                                                 10d17s22
              nab4(1,1)=nab4(1,2)                                             10d17s22
              nab4(2,1)=nab4(2,2)                                             10d17s22
             end if                                                           10d17s22
            else if(ndifs.eq.0.and.ndifd.eq.2)then                            10d14s22
             nnot=3
             do i=1,norbxx
              if(btest(gandcb,i))then                                         10d14s22
               if(btest(itvnotvc,i))then
                nab4(1,1)=i
                nab4(1,2)=i
               else                                                           10d14s22
                nab4(2,1)=i
                nab4(2,2)=i
               end if
              end if                                                          10d14s22
             end do
            end if                                                            10d14s22
            if(nnot.eq.3)then                                           11d1s22
             ipssx=1                                                    12d8s20
            else if(nnot.eq.4)then                                      11d1s22
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
              write(6,*)('already double in nab4(2,1) = '),nab4(2,iu2)  12d8s20
              stop 'nab4(2,1)'                                           11d27s20
             else                                                          11d13s20
              itesto=ibset(itesto,nab4(2,iu2))                          12d8s20
              nopenk=nopenk+1                                              11d13s20
             end if                                                        11d13s20
             call gandc(itvnotvc,itvnotvo,itestc,itesto,nopen2p,nopenk,      12d8s20
     $         iarg2,karg,ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1, 12d8s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
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
              call enough('hcdi.  3',bc,ibc)
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
                call enough('hcdi.  4',bc,ibc)
                call dgemm('n','n',nl(l),nroot,ncsf2(l,iarg2),1d0,         12d8s20
     $            bc(itmpt(l)),nl(l),bc(jtmp),ncsf(iarg2),0d0,          12d9s20
     $            bc(itmp2),nl(l),                                      12d8s20
     d' hcdi.  1')
                iadu=jvcvu+ibc(jvcvu+5+l)                               3d19s21
                jtmp2=itmp2                                                12d8s20
                do ir=0,nroot-1                                            12d8s20
                 kden=jden-1+nfdat(2,l,isb)*ir                          12d9s20
                 do i=0,nl(l)-1                                            12d8s20
                  ip=ibc(iadu+i)                                           12d8s20
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
              stop
             end if                                                     12d8s20
            end do                                                      12d8s20
           end if                                                       10d13s22
           jvs=jvs+ncsf(iarg0)                                          12d8s20
          end do                                                        12d8s20
          end if                                                        3d8s21
          idoit=idoit+1                                                 3d8s21
          jvcvu=jvcvu+nspace                                            3d19s21
         end do
        end do
        jvcv(isb)=jvcvu
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
         if(min(irefo(isa),irefo(isb),nn,nfdat(3,l,jsb)).gt.0)then      3d16s21
          nrr=irefo(isa)*irefo(isb)                                     12d8s20
          sz=0d0
          do i=0,nrr*nn-1
           sz=sz+bc(idenvnotvk(l,jsb,isa)+i)**2                         12d8s20
          end do
          sz=sqrt(sz/dfloat(nrr*nn))                                    12d8s20
          itmp=ibcoff                                                   12d8s20
          ibcoff=itmp+nrr*nn                                            12d8s20
          call enough('hcdi.  5',bc,ibc)
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
          call dgemm('n','n',nrow,nfdat(3,l,jsb),nfdat(2,l,jsb),1d0,    12d8s20
     $         bc(itmp),nrow,bc(nfdat(4,l,jsb)),nfdat(2,l,jsb),0d0,     12d8s20
     $         bc(idenvnotvk(l,jsb,isa)),nrow,                          12d8s20
     d' hcdi.  2')
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
          nn=nroot*nfdat(3,1,isymmrci)                                  12d9s20
          if(min(nn,nvv).gt.0)then                                      11d17s21
           iputto=ioff+ivn                                               12d8s20
           do isa=1,nsymb                                                12d8s20
            if(irefo(isa).gt.0)then                                     12d8s20
             nrr=irefo(isa)*irefo(isa)                                  12d8s20
             itmpi=ibcoff                                               12d8s20
             ibcoff=itmpi+nvv*nrr                                       12d8s20
             call enough('hcdi.  6',bc,ibc)
             i10=i1s                                                    12d8s20
             i1n=nvirt(jsbv1)
             i2eu=invk1(1,isa,isa,jsbv1,1)                              12d8s20
             kint=kmats(i2eu)                                           12d8s20
             jtmpi=itmpi                                                12d8s20
             do i2=i2s,i2e
              if(i2.eq.i2e)i1n=i1e
              if(i2.ge.i10.and.i2.le.i1n)then                           1d28s21
               irow=i2+nvirt(jsbv1)*(i2-1)                               12d8s20
               kintu=kint+irow-il                                        12d8s20
               jtmpi=itmpi+i2-ivn                                        12d8s20
               do i34=0,nrr-1                                            12d8s20
                bc(jtmpi)=bc(kintu)                                      12d8s20
                jtmpi=jtmpi+nvv                                          12d8s20
                kintu=kintu+nhere                                        12d8s20
               end do                                                    12d8s20
              end if                                                    1d28s21
              i10=1                                                     12d8s20
             end do                                                     12d8s20
             call dgemm('n','n',nvv,nn,nrr,srh,bc(itmpi),nvv,           3d8s21
     $             bc(idenvnotvk(1,isymmrci,isa)),nrr,1d0,gb(iputto),   3d8s21
     $             nvirt(jsbv1),                                        12d8s20
     d' hcdi.  3')
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
            ibcoff=itmpt(l)+nhere*nroot*nfdat(3,l,jsb)                  12d8s20
           end do                                                       12d8s20
           call enough('hcdi.  7',bc,ibc)
           do i=itmpt(1),ibcoff-1                                       12d8s20
            bc(i)=0d0                                                   12d8s20
           end do                                                       12d8s20
           do isa=1,nsymb                                                12d8s20
            isb=multh(isa,jsbv12)                                       12d8s20
            if(min(irefo(isa),irefo(isb)).gt.0)then                     12d8s20
             nrr=irefo(isa)*irefo(isb)                                  12d8s20
             i2eu=invk1(1,isa,isb,jsbv1,1)                              12d8s20
             icase=invk1(2,isa,isb,jsbv1,1)                              12d8s20
             if(jsbv12.eq.1)then                                        12d8s20
              itmpi=ibcoff                                               12d8s20
              ibcoff=itmpi+nhere*nrr                                     12d8s20
              call enough('hcdi.  8',bc,ibc)
              i10=i1s                                                    12d8s20
              i1n=nvirt(jsbv1)                                           12d8s20
              jtmpi=itmpi                                               12d8s20
              do i2=i2s,i2e                                             12d8s20
               i2m=i2-1                                                 12d8s20
               if(i2.eq.i2e)i1n=i1e                                     12d8s20
               do i1=i10,min(i1n,i2-1)                                  1d28s21
                kint=kmats(i2eu)+i1+nvirt(jsbv1)*i2m-il                 12d8s20
                ktmpi=jtmpi                                             12d8s20
                do i34=1,nrr                                            12d8s20
                 bc(ktmpi)=bc(kint)                                     12d8s20
                 ktmpi=ktmpi+nhere                                      12d8s20
                 kint=kint+nhere                                        12d8s20
                end do                                                  12d8s20
                jtmpi=jtmpi+1                                           12d8s20
               end do                                                   12d8s20
               i10=1                                                    12d8s20
              end do                                                    12d8s20
              nvvh=jtmpi-itmpi                                          12d8s20
              kint=itmpi                                                12d8s20
             else                                                       12d8s20
              kint=kmats(i2eu)                                          12d8s20
              nvvh=nhere                                                12d8s20
             end if                                                     12d8s20
             if(nvvh.gt.0)then                                          12d8s20
              do l=1,4                                                       12d8s20
               if(nfdat(3,l,jsb).gt.0)then                                  12d8s20
                nn=nfdat(3,l,jsb)*nroot                                  12d8s20
                call dgemm('n','n',nvvh,nn,nrr,1d0,bc(kint),nhere,       12d8s20
     $              bc(idenvnotvk(l,jsb,isa)),nrr,1d0,                  12d8s20
     $              bc(itmpt(l)),nhere,                                 12d8s20
     d' hcdi.  4')
               end if                                                    12d8s20
              end do                                                     12d8s20
             end if                                                      12d8s20
             if(jsbv12.eq.1)ibcoff=itmpi                                12d9s20
            end if                                                      12d8s20
           end do                                                       12d8s20
           joff=ioff                                                    12d8s20
           do l=1,4                                                     12d8s20
            if(nfdat(3,l,jsb).gt.0)then                                 12d8s20
             nn=nfdat(3,l,jsb)*nroot                                    12d8s20
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
           ioff=ioff+nvv*nfdat(3,l,jsb)*nroot                           12d8s20
          end do                                                        12d8s20
         end if                                                         12d8s20
        end if                                                          12d8s20
       end do                                                           12d8s20
      end do                                                            12d8s20
      call dws_gsumf(gb,ioff)                                           8d28s22
      if(lprt)then                                                      1d25s21
       igbtest=ibcoff
       ibcoff=igbtest+ioff                                               12d9s20
       call enough('hcdi.  9',bc,ibc)
       ioff=1
       do jsb=1,nsymb                                                    12d8s20
        jsbv12=multh(jsb,isymmrci)                                       12d8s20
        do jsbv1=1,nsymb                                                 12d8s20
         jsbv2=multh(jsbv1,jsbv12)
         if(jsbv2.ge.jsbv1)then
          write(6,*)('for symmetries '),jsb,jsbv1,jsbv2,ioff
          igoff=igbtest
          do i=igbtest,ibcoff-1
           bc(i)=0d0                                                     12d9s20
          end do                                                         12d9s20
          if(jsbv12.eq.1.and.nvirt(jsbv1).gt.0)then
           ndoub=nfdat(3,1,isymmrci)*((nvirt(jsbv1)*(nvirt(jsbv1)+1))/2)
           write(6,*)('visv'),ioff
           nn=nfdat(3,1,isymmrci)*nroot
           call prntm2(gb(ioff),nvirt(jsbv1),nn,nvirt(jsbv1))
           do if=0,nfdat(3,1,isymmrci)-1
            do ir=0,nroot-1
             do iv=0,nvirt(jsbv1)-1
              ivv=((iv*(iv+1))/2)+iv                                     12d9s20
              iad=igoff+if+nfdat(3,1,isymmrci)*ivv+ndoub*ir              12d9s20
              bc(iad)=gb(ioff)
              ioff=ioff+1
             end do
            end do
           end do
          end if
          if(jsbv12.eq.1)then
           nvvs=(nvirt(jsbv1)*(nvirt(jsbv1)+1))/2                         12d8s20
           nvvt=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                         12d8s20
           isw=0
          else                                                           12d8s20
           nvvs=nvirt(jsbv1)*nvirt(jsbv2)                                 12d8s20
           nvvt=nvvs
           isw=1
          end if                                                         12d8s20
          if(nvvt.gt.0)then
           write(6,*)('vnotv')
           ip=+1
           nvv=nvvs
           do l=1,4
            if(nfdat(3,l,jsb).gt.0)then
             write(6,*)('for l = '),l,ioff
             nn=nfdat(3,l,jsb)*nroot
             call prntm2(gb(ioff),nvvt,nn,nvvt)
             ndoub=nfdat(3,l,jsb)*nvv
             do if=0,nfdat(3,l,jsb)-1
              do ir=0,nroot-1
               do iv2=0,nvirt(jsbv2)-1
                if(jsbv12.eq.1)then
                 itop=iv2-1
                else
                 itop=nvirt(jsbv1)-1
                end if
                do iv1=0,itop
                 itri=((iv2*(iv2+ip))/2)+iv1
                 irec=iv1+nvirt(jsbv1)*iv2
                 itri=itri+isw*(irec-itri)
                 iad=igoff+if+nfdat(3,l,jsb)*itri+ndoub*ir
                 bc(iad)=gb(ioff)
                 ioff=ioff+1
                end do
               end do
              end do
             end do                                                      12d9s20
             write(6,*)('global summed result for mpec2.0')
             call prntm2(bc(igbtest),ndoub,nroot,ndoub)
             igoff=igbtest
            end if
            ip=-1
            nvv=nvvt
           end do
          end if                                                         12d8s20
         end if
        end do                                                           12d8s20
       end do                                                            12d8s20
      end if                                                            1d25s21
      ibcoff=ibcoffo                                                    11d23s20
      return                                                            10d19s20
      end                                                               10d19s20
