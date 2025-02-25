c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcdibk4(nff2,iff22,ff22,nfdat,ncsfd,ncsfd2,gb,i2sb,    11d15s21
     $     i2smb,mdoobp,isymmrci,nff0,iff0,vec,ncsfv,i2sk,i2smk,mdookp, 11d15s21
     $     nec,mdon,nsymb,multh,irw1,irw2,nvirt,nroot,lprt,             11d15s21
     $     isymc,irorip,isopt,ism,irel,irefo,norb,kmatsrb,iifmx,ntype,  11d16s21
     $     srh,bc,ibc)                                                  11d9s22
      implicit real*8 (a-h,o-z)                                         10d19s20
      integer*8 itvisvc,itvisvo,itvnotvc,itvnotvo,jrefc,jrefo,ipack8,   12d8s20
     $     itestc,itesto,last8(2),iff22(*),ipack,gandcc,gandco,gandcb   2d7s23
      integer*2 ipack2(4)                                               11d15s21
      integer*1 idorb(64),isorb(64),nab1(2),iorb(64),imap(64),icode(64),12d8s20
     $     nab2(2)                                                      12d8s20
      logical ldebug,lbail,lprt                                         1d25s21
      dimension nff2(mdoobp,2,nsymb),multh(8,8),nvirt(*),                8d6s21
     $     ncsfd(*),ncsfd2(4,*),nfdat(5,4,*),irel(*),ism(*),jvcv(8),      12d8s20
     $     ipack4(2),itmpt(4),ff22(*),isopt(*),kmatsrb(8,8,*),iifmx(*), 11d15s21
     $     irefo(*),gb(*),idenvisv(8),nl(4),npre(4),mpre(4),imy(4),     11d15s21
     $     idenvnotvk(4,8,8),nff0(mdookp,*),iff0(*),isy(4),igya(4),     11d15s21
     $     vec(ncsfv,*),nab4(2,3),iwpb1(4),iwpk1(4),iwpb2(4),iwpk2(4),  11d15s21
     $     ipack2a(4),ndenvnotvk(8,8),noks(8,8),ioxx(2)                 2d7s23
      equivalence (ipack8,ipack4),(ipack,ipack2)                        11d15s21
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      if(lprt)write(6,*)('Hi, my name is hcdibk4'),ibcoff
      xnan=-2d0
      mdoo=mdoobp-1                                                     11d15s21
      irori=irorip-1                                                    11d15s21
      if(isopt(3).ne.0)then
       ieoro=1-irori                                                      9d6s21
      else                                                              10d14s21
       ieoro=irori                                                      9d6s21
      end if                                                            9d6s21
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
        nn=nroot*nfdat(2,l,jsb)*ntype                                   11d15s21
        do isa=1,nsymb                                                  11d23s20
         isb=multh(isopt(1),multh(isa,jsbv12))                          11d16s21
         idenvnotvk(l,isa,jsb)=ibcoff                                   11d15s21
         ibcoff=idenvnotvk(l,isa,jsb)+irefo(isa)*irefo(isb)*nn          11d15s21
        end do                                                          11d23s20
       end do                                                           11d23s20
       do isa=1,nsymb                                                   11d15s21
        isb=multh(isopt(1),multh(isa,jsbv12))                           11d16s21
        ndenvnotvk(isa,jsb)=ibcoff                                      11d15s21
        ibcoff=ndenvnotvk(isa,jsb)+irefo(isa)*irefo(isb)*ntype          11d15s21
       end do                                                           11d15s21
      end do                                                            10d19s20
      call enough('hcdibk4.  1',bc,ibc)
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
         do nclo0p=max(mdon+1,nclo2p-2),min(mdookp,ntop)                11d16s21
          if(nff0(nclo0p,1).gt.0)then                                   8d6s21
           nclo0=nclo0p-1                                                 12d7s20
           nopen=nec-nclo0*2                                              12d7s20
           call wfetch(nopen,mdon,idum,i2sk,i2smk,ndetk,ncsfk,iveck,    11d15s21
     $          iaorbk,idum,bc,ibc)                                     11d9s22
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
               ibcoff=itmpt(l)+ncsfd2(l,iarg2)*nl(l)                    11d15s21
               call enough('hcdibk4.  2',bc,ibc)
               do i=0,nl(l)-1                                                    11d20s20
                do j=0,ncsfd2(l,iarg2)-1                                11d15s21
                 ji=iad2+j+ncsfd2(l,iarg2)*i                            11d15s21
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
                 if(btest(jrefc,nab4(2,2)).and.                         2d6s23
     $                .not.btest(jrefc,nab4(2,1)))nbt=1                 2d6s23
                else                                                    2d6s23
                 nbt=0                                                  2d6s23
                 if(btest(itvnotvc,nab4(1,2)).and.                      2d6s23
     $                .not.btest(itvnotvc,nab4(1,1)))nbt=1              2d6s23
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
               if(nnot.ne.0)then                                        2d7s23
                iu1=1
                iu2=1
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
                call gandcr(itvnotvc,itvnotvo,itestc,itesto,nopen2p,     11d15s21
     $               nopenk,norbxx,nnot1,nab1,icode,imap,nx,irw1,irw2,  11d15s21
     $               iwpb1,iwpk1,bc,ibc)                                11d14s22
                call gandcr(itestc,itesto,jrefc,jrefo,nopenk,nopen,       11d15s21
     $         norbxx,nnot2,nab2,icode,imap,nx,irw1,irw2,iwpb2,iwpk2,   11d14s22
     $              bc,ibc)                                             11d14s22
                if(nnot1.eq.2.and.nnot2.eq.2)then                              11d13s20
                 call spinloop(i2sb,i2smb,i2sk,i2smk,nopen2p,nopen,      11d15s21
     $               nopenk,ncsfd(iarg2),nroot,itype,imatx,ntypeq,nab1, 11d16s21
     $              iwpb1,iwpk1,nab2,iwpb2,iwpk2,vec(jvs,1),ncsfv,ieoro,11d14s22
     $              bc,ibc)                                             11d14s22
                 do iq=0,ntypeq-1                                        11d16s21
                  ipack=ibc(itype+iq)                                    11d15s21
                  itmpxq=imatx+ncsfd(iarg2)*nroot*iq                     11d15s21
                  do j=1,4                                               9d9s21
                   ipa=ipack2(j)                                         11d15s21
                   ipa=iabs(ipa)                                         11d15s21
                   ipack2a(j)=ipa                                        11d15s21
                   if(ipa.le.norb)then                                   11d15s21
                    if(ipack2(j).gt.0)then                                9d9s21
                     isy(j)=ism(ipack2(j))                                9d9s21
                     igya(j)=irel(ipack2(j))-1                            11d15s21
                     imy(j)=0                                             10d12s21
                    else                                                  9d9s21
                     isy(j)=ism(-ipack2(j))                                9d9s21
                     igya(j)=irel(-ipack2(j))-1                           10d12s21
                     imy(j)=1                                             10d12s21
                    end if                                                9d9s21
                   else                                                  11d15s21
                    igya(j)=-1                                           11d15s21
                    if(ipack2(j).gt.0)then                               11d15s21
                     imy(j)=0                                            11d15s21
                    else                                                 11d15s21
                     imy(j)=1                                            11d15s21
                    end if                                               11d15s21
                   end if                                                11d15s21
                  end do                                                 9d9s21
                  if(ipack2a(1).gt.norb.and.ipack2a(3).gt.norb)then      11d16s21
                   if(nnot.eq.3)then                                     11d16s21
                    ngothru=1                                            11d16s21
                   else                                                  11d16s21
                    ngothru=2                                            11d16s21
                   end if                                                11d16s21
                   do igothru=1,ngothru                                  11d16s21
                    if(igothru.eq.1)then                                 11d16s21
c
c     type is (v2|v'4)
c
                     ltest=imy(1)+2*(imy(2)+2*(imy(3)+2*imy(4)))           11d15s21
                     jtest=1-imy(1)+2*(1-imy(2)+2*(1-imy(3)              11d16s21
     $                   +2*(1-imy(4))))                                11d16s21
                     phs=1d0                                               11d15s21
                    else                                                 11d16s21
c
c     type is (v4|v'2)
c
                     ltest=imy(1)+2*(imy(4)+2*(imy(3)+2*imy(2)))           11d15s21
                     jtest=1-imy(1)+2*(1-imy(4)+2*(1-imy(3)              11d16s21
     $                   +2*(1-imy(2))))                                11d16s21
                     phs=-1d0                                               11d15s21
                    end if                                               11d16s21
                    itestp=ltest+1                                        11d15s21
                    jtestp=jtest+1                                        11d15s21
                    iuse=-1                                               11d15s21
                    if(iifmx(itestp).ge.0)then                             11d15s21
                     iuse=iifmx(itestp)                                    11d15s21
                    else if(iifmx(jtestp).ge.0)then                       11d15s21
                     iuse=iifmx(jtestp)                                   11d15s21
                     if(isopt(4).ne.0)phs=-phs                            11d15s21
                    end if                                                11d15s21
                    if(iuse.ge.0)then                                     11d15s21
                     jtmp=itmpxq                                          11d15s21
                     if(igothru.eq.1)then                                11d16s21
                      icol=igya(2)+irefo(isy(2))*(igya(4)+irefo(isy(4))    11d15s21
     $                  *iuse)                                          11d15s21
                      isuse=isy(2)                                       11d16s21
                     else                                                11d16s21
                      icol=igya(4)+irefo(isy(4))*(igya(2)+irefo(isy(2))  11d16s21
     $                  *iuse)                                          11d15s21
                      isuse=isy(4)                                       11d16s21
                     end if                                              11d16s21
                     bc(ndenvnotvk(isuse,isb)+icol)=1d0                  2d18s22
                     do l=1,4                                                    12d8s20
                      if(nl(l).gt.0)then                                         12d8s20
                       jden=idenvnotvk(l,isuse,isb)+                     11d16s21
     $                   nroot*nfdat(2,l,isb)*icol                      11d15s21
                       itmp2=ibcoff                                               12d8s20
                       ibcoff=itmp2+nl(l)*nroot                                   12d8s20
                       call enough('hcdibk4.  3',bc,ibc)
                       call dgemm('n','n',nl(l),nroot,ncsfd2(l,iarg2),   11d16s21
     $                     phs,bc(itmpt(l)),nl(l),bc(jtmp),ncsfd(iarg2),11d16s21
     $                     0d0,bc(itmp2),nl(l),                         11d16s21
     d' hcdibk4.  1')
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
                      jtmp=jtmp+ncsfd2(l,iarg2)                           11d15s21
                     end do                                                      12d8s20
                    end if                                                11d15s21
                   end do                                                11d16s21
                  else                                                   11d15s21
                   write(6,*)('confused about where virts are!!!'),
     $                 ipack2a                                          11d15s21
                   stop 'hcdibk4'                                        11d15s21
                  end if                                                 11d15s21
                 end do                                                  11d15s21
                 ibcoff=itype                                            11d15s21
                else
                 write(6,*)('nnots not both 2!!! '),nnot1,nnot2
                 call dws_synca
                 call dws_finalize
                 stop 'hcdibk4'
                end if                                                     12d8s20
               end if                                                   2d7s23
              end if                                                       12d8s20
              jvs=jvs+ncsfk                                             11d15s21
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
      call dws_gsumf(bc(ibcoffo),nden)                                  3d8s21
      do jsb=1,nsymb                                                    10d19s20
       jsbv12=multh(jsb,isymmrci)                                       12d8s20
       do l=1,4
        nn=nroot*nfdat(2,l,jsb)                                         12d8s20
        do isa=1,nsymb                                                  11d23s20
         isb=multh(isopt(1),multh(isa,jsbv12))                          11d16s21
         if(min(irefo(isa),irefo(isb),nn).gt.0)then
          nrr=irefo(isa)*irefo(isb)*ntype                               11d15s21
          nok=0                                                         11d15s21
          do i=0,nrr-1                                                  11d15s21
           if(bc(ndenvnotvk(isa,jsb)+i).gt.0.1d0)nok=nok+1              2d18s22
          end do                                                        11d15s21
          itmp=ibcoff                                                   11d15s21
          ibcoff=itmp+nok*nn                                            11d15s21
          call enough('hcdibk4.  4',bc,ibc)
          jtmp=itmp                                                     11d15s21
c
c     dims are nfdat,ir,ic,id,type.
c     transpose to (compressed ic,id,type),ir,nfdat
c
          do i=0,nrr-1                                                  11d15s21
           if(bc(ndenvnotvk(isa,jsb)+i).gt.0.1d0)then                   2d18s22
            do ir=0,nroot-1                                              11d15s21
             do j=0,nfdat(2,l,jsb)-1                                     11d15s21
              iad1=idenvnotvk(l,isa,jsb)+j+nfdat(2,l,jsb)*(ir+nroot*i)  11d15s21
              iad2=jtmp+nok*(ir+nroot*j)                                 11d15s21
              bc(iad2)=bc(iad1)                                          11d15s21
             end do                                                     11d15s21
            end do                                                      11d15s21
            jtmp=jtmp+1                                                 11d15s21
           end if                                                       11d15s21
          end do                                                        11d15s21
          if(nok.gt.0)then
           sz=0d0                                                        11d16s21
           do i=0,nok*nn-1                                               11d15s21
            bc(idenvnotvk(l,isa,jsb)+i)=bc(itmp+i)                       11d15s21
            sz=sz+bc(itmp+i)**2
           end do                                                        11d15s21
           sz=sqrt(sz/dfloat(nok*nn))
          end if
          ibcoff=itmp                                                   11d15s21
         end if
        end do                                                          11d23s20
       end do                                                           11d23s20
       do isa=1,nsymb                                                    11d15s21
        isb=multh(isopt(1),multh(isa,jsbv12))                           11d16s21
        if(min(irefo(isa),irefo(isb)).gt.0)then                          11d15s21
         nrr=irefo(isa)*irefo(isb)*ntype                                 11d15s21
         nrrm=nrr-1                                                      11d15s21
         nok=0                                                           11d15s21
         do i=0,nrrm                                                    11d15s21
          if(ibc(ndenvnotvk(isa,jsb)+i).gt.0)then                       11d15s21
           ibc(ndenvnotvk(isa,jsb)+nok)=i                               11d15s21
           nok=nok+1                                                    11d15s21
          end if                                                        11d15s21
         end do                                                         11d15s21
         noks(isa,jsb)=nok                                              11d15s21
        else                                                             11d15s21
         noks(isa,jsb)=0                                                11d15s21
        end if                                                           11d15s21
       end do                                                            11d15s21
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
          if(nvv.gt.0)then                                              12d8s20
           iputto=ioff+ivn                                               12d8s20
           do isa=1,nsymb                                                12d8s20
            if(min(nn,noks(isa,jsb),irefo(isa)).gt.0)then               11d17s22
             isb=multh(isa,isopt(1))                                    11d16s21
             nrr=noks(isa,jsb)                                          11d15s21
             itmpi=ibcoff                                               12d8s20
             ibcoff=itmpi+nvv*nrr                                       12d8s20
             nrrm=nrr-1                                                 11d15s21
             call enough('hcdibk4.  5',bc,ibc)
             i10=i1s                                                    12d8s20
             i1n=nvirt(jsbv1)
             jtmpi=itmpi                                                12d8s20
             do i2=i2s,i2e
              if(i2.eq.i2e)i1n=i1e
              if(i2.ge.i10.and.i2.le.i1n)then                           1d28s21
               irow=i2+nvirt(jsbv1)*(i2-1)-il                           11d15s21
               kint=kmatsrb(isa,jsbv2,isb)+irow                         11d16s21
               jtmpi=itmpi+i2-ivn                                        12d8s20
               do i34=0,nrrm
                icol=nhere*ibc(ndenvnotvk(isa,jsb)+i34)                 2d5s22
                bc(jtmpi+i34*nvv)=bc(kint+icol)                         11d15s21
               end do                                                   11d15s21
              end if                                                    1d28s21
              i10=1                                                     12d8s20
             end do                                                     12d8s20
             call dgemm('n','n',nvv,nn,noks(isa,jsb),srh,bc(itmpi),     11d15s21
     $            nvv,bc(idenvnotvk(1,isa,jsb)),noks(isa,jsb),          11d15s21
     $            1d0,gb(iputto),nvirt(jsbv1),                          11d15s21
     d' hcdibk4.  2')
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
           call enough('hcdibk4.  6',bc,ibc)
           do i=itmpt(1),ibcoff-1                                       12d8s20
            bc(i)=0d0                                                   12d8s20
           end do                                                       12d8s20
           do isa=1,nsymb                                                12d8s20
            isb=multh(isopt(1),multh(isa,jsbv12))                       11d16s21
            if(min(noks(isa,jsb),irefo(isa),irefo(isb)).gt.0)then       11d16s21
             nrr=noks(isa,jsb)                                          11d15s21
             nrrm=nrr-1                                                 11d15s21
             if(jsbv12.eq.1)then                                        12d8s20
              jsw=0                                                     8d6s21
             else                                                       8d6s21
              jsw=1                                                     8d6s21
             end if                                                     8d6s21
             itmpi=ibcoff                                               12d8s20
             ibcoff=itmpi+nhere*nrr                                     12d8s20
             call enough('hcdibk4.  7',bc,ibc)
             do iz=itmpi,ibcoff-1
              bc(iz)=xnan
             end do
             jtmpi=itmpi                                                8d6s21
             i10=i1s                                                    11d15s21
             i1n=nvirt(jsbv1)                                           11d15s21
             do i2=i2s,i2e                                              11d15s21
              i2m=i2-1                                                  11d15s21
              if(i2.eq.i2e)i1n=i1e                                      11d15s21
              itop=i2m+jsw*(nvirt(jsbv1)-i2m)                           11d16s21
              do i1=i10,min(itop,i1n)                                   11d15s21
               kint=kmatsrb(isa,jsbv2,isb)+i1+nvirt(jsbv1)*i2m-il       11d15s21
               do i34=0,nrrm                                            11d15s21
                icol=ibc(ndenvnotvk(isa,jsb)+i34)                       11d15s21
                bc(jtmpi+i34*nhere)=bc(kint+nhere*icol)                 11d15s21
               end do                                                   11d15s21
               jtmpi=jtmpi+1                                            11d15s21
              end do                                                    11d15s21
              i10=1                                                     11d15s21
             end do                                                     11d15s21
             nvvh=jtmpi-itmpi                                           12d8s20
             kint=itmpi                                                 12d8s20
             if(nvvh.gt.0)then                                          12d8s20
              do l=1,4                                                       12d8s20
               if(nfdat(2,l,jsb).gt.0)then                              8d6s21
                nn=nfdat(2,l,jsb)*nroot                                 8d6s21
                call dgemm('n','n',nvvh,nn,nrr,1d0,bc(kint),nhere,       12d8s20
     $              bc(idenvnotvk(l,isa,jsb)),nrr,1d0,                  11d15s21
     $              bc(itmpt(l)),nhere,                                 12d8s20
     d' hcdibk4.  3')
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
