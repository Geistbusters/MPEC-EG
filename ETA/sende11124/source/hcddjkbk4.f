c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcddjkbk4(nff22b,iff22b,ff22b,nfdatb,ncsfb,ncsfb2,gd,  11d19s21
     $     i2sb,i2smb,mdoobp,isymbra,                                   11d19s21
     $     nff22k,iff22k,ff22k,nfdatk,ncsfk,ncsfk2,vd,i2sk,i2smk,       11d19s21
     $     mdookp,isymket,                                              11d19s21
     $     nec,mdon,nsymb,multh,irw0,irw1,irw2,nvirt,nrootu,ldebug,     11d19s21
     $     isymc,irorip,isopt,ism,irel,irefo,norb,ih0n,nh0,i4or,jmatsr, 11d19s21
     $     kmatsr,iifmx,ntype,sr2,srh,iall,idv4,ndv4,iuall,l2e,bc,ibc,  2d8s23
     $     n4vso)                                                       2d8s23
      implicit real*8 (a-h,o-z)                                         1d8s21
c                                                                       1d8s21
c     0 and 2 virt contribution to gd=hdd*vd                            1d8s21
c                                                                       1d8s21
      logical lprt,lchoice,ldebug,lkeep                                 2d10s22
      integer*1 nab1(2),nab2(2),nab1b(2),nab2b(2),imap(64),icode(64),   11d19s21
     $     ipackc1(8)                                                   11d19s21
      integer*2 ipack2(4)                                               11d19s21
      integer*8 ipack8,itc,ito,jtc,jto,itestc,itesto,itesta,itestb,     1d17s21
     $     jtesta,jtestb,iff22b(*),iff22k(*),ipack,ipackc,gandcc,gandco,2d7s23
     $     gandcb,gandcc2,gandco2,gandcb2                               2d8s23
      external second                                                   2d18s21
      equivalence (ipack8,ipack4),(ipack,ipack2),(ipackc,ipackc1)       11d19s21
      dimension nff22b(mdoobp,2,nsymb),nfdatb(5,4,*),vd(*),gd(*),          8d9s21
     $     multh(8,8),nvirt(*),ncsfb(*),ncsfb2(4,*),irel(*),ism(*),       1d8s21
     $     irefo(*),ncsfk(*),ncsfk2(4,*),isy(4),imy(4),igya(4),         11d19s21
     $     ipack4(2),nprei(4),mprei(4),nli(4),nprej(4),ipack2a(4),      11d22s21
     $     mprej(4),nlj(4),nab4(2,3),idenj(4,8),ndenj(4,4,8),           11d23s21
     $     mdenj(4,8),idenk(4,8),ndenk(4,8),mdenk(4,8),loope(6),                 1d9s21
     $     idenhvv(4,4),mdenhvv(4),itest(32,3),iden1e(4),mden1e(4),     11d26s21
     $     ldenj(4,8),ldenk(4,8),iden1ef(4),idenhvvf(4),nh0(*),         11d22s21
     $     idenjf(4,4,8),idenkf(4,4,8),ndenjf(4,8),ndenkf(4,4,8),
     $     nokj(4,4,8),nokk(4,4,8),isorb(32),idorb(32),idenjn(4,4,8),   11d23s21
     $     idenkn(4,4,8),iden1en(4,4),idenhvvn(4,4),ff22b(*),itmpl(4),  11d26s21
     $     nff22k(mdookp,2,nsymb),nfdatk(5,4,*),ff22k(*),mcsf(2),       11d19s21
     $     ih0n(*),i4or(8,8,8),jmatsr(8,8,8),kmatsr(8,8,8),iifmx(*),    11d19s21
     $     isopt(*),iwpb1(4),iwpk1(4),iwpb2(4),iwpk2(4),iall(8,8,8),    11d30s21
     $     iden4(4,4),idv4(2,4,8),idv4l(2,4),idenhvvx(4),idv4y(4,8),    2d8s22
     $     idv4x(4,8),idenhvvy(4),itermy(2,2),termy(2,2),itermx(2,2),   2d8s22
     $     termx(2,2),idv4yl(4),ioxx(2)                                 2d7s23
      data loopx/1900000000/
      include "common.store"                                            1d8s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data icall/0/
      save icall
      icall=icall+1
      lprt=.false.
      loop=0
      irori=irorip-1                                                    11d17s21
      if(isopt(3).ne.0)then
       ieoro=1-irori                                                      9d6s21
      else                                                              10d14s21
       ieoro=irori                                                      9d6s21
      end if                                                            9d6s21
      mdoob=mdoobp-1                                                    11d19s21
      mdook=mdookp-1                                                    11d19s21
      ibcoffo=ibcoff
      norbx=norb+1                                                      1d8s21
      norbxx=norbx+1                                                    1d8s21
      norbxxx=norbxx+1                                                  1d8s21
      norbxxxx=norbxxx+1                                                11d30s21
      ioffdnon=1                                                        1d8s21
      igoal=1+5*2
      igoalj=1
      igoaljj=1
      igoalp=2
      igoalpp=2
      igoal1=1
      igoal4=1
      igoali=2
      igoalii=2
      igoalk=2
      igoal2=2
      igoal1j=2
      igoal2j=2
      nrootm=nrootu-1                                                   1d11s21
c
c     space for partially transformed vectors for 4v integrals
c
      ib4=ibcoff                                                        12d3s21
      do isb=1,nsymb                                                    12d3s21
       isbv12=multh(isb,isymket)                                        8d9s21
       nvisv=0                                                          12d3s21
       nvnotv=0                                                         12d3s21
       do isbv1=1,nsymb                                                 12d3s21
        isbv2=multh(isbv1,isbv12)                                       12d3s21
        if(isbv2.ge.isbv1)then                                          12d3s21
         if(isbv1.eq.isbv2)then                                         12d3s21
          nvisv=nvisv+nvirt(isbv1)                                      12d3s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         12d3s21
         else                                                           12d3s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 12d3s21
         end if                                                         12d3s21
         nvnotv=nvnotv+nvv                                              12d3s21
        end if                                                          12d3s21
       end do                                                           12d3s21
       nvisv0=nvisv*nrootu*2                                            2d8s22
       nvisv=nvisv*nrootu*ntype                                         12d3s21
       nvnotv0=nvnotv*nrootu*2                                          2d8s22
       nvnotv=nvnotv*nrootu*ntype                                        12d3s21
       do lb=1,4                                                        12d3s21
        idv4(1,lb,isb)=ibcoff                                           12d3s21
        ibcoff=idv4(1,lb,isb)+nvnotv*nfdatb(2,lb,isb)                   12d3s21
        idv4x(lb,isb)=ibcoff                                            2d9s22
        ibcoff=idv4x(lb,isb)+nvnotv0*nfdatb(2,lb,isb)                   2d9s22
        if(nvisv.ne.0)then                                              12d3s21
         idv4(2,lb,isb)=ibcoff                                          12d3s21
         ibcoff=idv4(2,lb,isb)+nvisv*nfdatb(2,lb,isb)                   12d3s21
         idv4y(lb,isb)=ibcoff                                           2d8s22
         ibcoff=idv4y(lb,isb)+nvisv0*nfdatb(2,lb,isb)                   2d8s22
        end if                                                          12d3s21
       end do                                                           12d3s21
      end do                                                            12d3s21
      call enough('hcddjkbk4.  1',bc,ibc)
      do iz=ib4,ibcoff-1                                                12d3s21
       bc(iz)=0d0                                                       12d3s21
      end do                                                            12d3s21
      ndv4=ibcoff-ib4                                                   12d3s21
      nhvij=0                                                           2d7s22
      nhvii=0                                                           2d7s22
      do isb=1,nsymb                                                    1d8s21
       isbv12=multh(isb,isymket)                                        8d9s21
       joffdnon=1                                                       1d8s21
       nhvabp=0                                                         2d7s22
       nhvabm=0                                                         2d7s22
       do jsb=1,nsymb                                                   8d10s21
        ijsb=multh(jsb,isb)                                             1d9s21
        jsbv12=multh(jsb,isymbra)                                       8d9s21
        ibc0=ibcoff                                                     1d8s21
        nden4=0                                                         11d30s21
        do l=1,4                                                        1d8s21
         if(nfdatk(2,l,isb).gt.0)then                                   8d26s21
          do ll=1,4                                                     1d8s21
           if(nfdatb(2,ll,jsb).gt.0)then                                8d26s21
            nll=nfdatk(2,l,isb)*nfdatb(2,ll,jsb)                        8d26s21
            iden1en(l,ll)=ibcoff                                          2d24s21
            ibcoff=iden1en(l,ll)+nll                                    11d23s21
             idenhvvn(l,ll)=ibcoff                                       11d26s21
             ibcoff=idenhvvn(l,ll)+nll                                  2d8s22
             iden4(l,ll)=ibcoff                                         2d8s22
             ibcoff=iden4(l,ll)+nll*ntype                                11d30s21
            do lsa=1,nsymb                                              1d8s21
c
c     lsa*lsb*v*v'=isymc=symbra*symket
c     isb*symket=v1*v2
c     jsb*symbra=v2*v3.
c     isb*isymket*jsb*isymbra=v1*v3=ijsb*isymc
             lsb=multh(lsa,ijsb)                                        11d22s21
             nn=irefo(lsa)*irefo(lsb)*ntype                             11d19s21
             ndenkf(l,ll,lsa)=ibcoff                                    2d25s21
             ibcoff=ndenkf(l,ll,lsa)+nn                                     1d11s21
             idenkn(l,ll,lsa)=ibcoff                                    1d9s21
             ibcoff=idenkn(l,ll,lsa)+nn*nll                             2d24s21
             idenjn(l,ll,lsa)=ibcoff                                      1d11s21
             ibcoff=idenjn(l,ll,lsa)+nn*nll                             11d23s21
             ndenj(l,ll,lsa)=ibcoff                                     11d23s21
             ibcoff=ndenj(l,ll,lsa)+nn                                  11d23s21
            end do                                                      1d8s21
           end if                                                       1d8s21
          end do                                                        1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
        call enough('hcddjkbk4.  2',bc,ibc)
        nwds=ibcoff-ibc0                                                1d8s21
        do i=ibc0,ibcoff-1                                              1d8s21
         bc(i)=0d0                                                      1d8s21
        end do                                                          1d8s21
        ibclast=ibcoff                                                  1d17s21
        idoit=0                                                         2d23s21
        nlprt=0
        do ncloi=mdon,mdook                                             11d19s21
         ncloip=ncloi+1                                                 1d8s21
         if(nff22k(ncloip,1,isb).gt.0)then                              8d26s21
          iarg=ncloip-mdon                                              1d8s21
          nopeni=nec-2*ncloip                                           1d8s21
          nopenip=nopeni+2                                              1d8s21
          ibcst=ibcoff                                                  1d8s21
          ibcnd=ibcoff-1                                                1d8s21
          ivcv=nfdatk(5,1,isb)+nff22k(ncloip,2,isb)                     8d26s21
          do if=1,nff22k(ncloip,1,isb)                                  8d26s21
           ipack8=iff22k(ivcv)                                          8d26s21
           itc=ipack4(1)                                                1d8s21
           ito=ipack4(2)                                                1d8s21
           ito=ibset(ito,norbx)                                         1d8s21
           ito=ibset(ito,norbxx)                                        1d8s21
           nspacei=iff22k(ivcv+1)                                       8d26s21
           do l=1,4                                                     3d19s21
            nli(l)=iff22k(ivcv+1+l)                                     8d26s21
           end do                                                       3d19s21
           do i=1,norb                                                  1d8s21
            itest(i,1)=0                                                1d8s21
           end do                                                       1d8s21
           do i=1,norb                                                  1d8s21
            if(btest(itc,i))itest(i,1)=2                                1d8s21
            if(btest(ito,i))itest(i,1)=1                                1d8s21
           end do                                                       1d8s21
           itest(norbx,1)=1                                             1d8s21
           itest(norbxx,1)=1                                            1d8s21
           do ncloj=max(ncloi-2,mdon),min(ncloi+2,mdoob)                11d19s21
            nclojp=ncloj+1                                                1d8s21
            if(nff22b(nclojp,1,jsb).gt.0)then                           8d26s21
             if(mod(idoit,mynprocg).eq.mynowprog)then                   3d1s21
              jarg=nclojp-mdon                                             1d8s21
              nopenj=nec-2*nclojp                                          1d8s21
              nopenjp=nopenj+2                                           1d8s21
              jvcv=nfdatb(5,1,jsb)+nff22b(nclojp,2,jsb)                  8d26s21
              do jf=1,nff22b(nclojp,1,jsb)                               8d26s21
               ipack8=iff22b(jvcv)                                       8d26s21
               jtc=ipack4(1)                                               1d8s21
               jto=ipack4(2)                                               1d8s21
               jto=ibset(jto,norbx)                                        1d8s21
               jto=ibset(jto,norbxx)                                       1d8s21
               nspacej=iff22b(jvcv+1)                                    8d26s21
               do l=1,4                                                  3d19s21
                nlj(l)=iff22b(jvcv+1+l)                                  8d26s21
               end do                                                    3d19s21
c     j is bra, i is ket
               if(lprt)nlprt=nlprt+1
               if(lprt)write(6,*)('yes, we have lprt!'),isbv12,jsbv12,
     $             nlprt
               if(isbv12.eq.jsbv12)then                                  8d9s21
c                                                                       1d8s21
                gandcc=ieor(itc,jtc)                                      10d13s22
                gandco=ieor(ito,jto)                                      10d13s22
                gandcb=ior(gandcc,gandco)                                 10d20s22
                ndifb=popcnt(gandcb)                                      10d20s22
                if(ndifb.le.4)then                                        10d20s22
                 ndifs=popcnt(gandco)                                      10d13s22
                 ndifd=popcnt(gandcc)                                      10d13s22
                 iunit=ibcoff                                           11d19s21
                 ibcoff=iunit+ncsfk(iarg)*ncsfk(iarg)                   11d19s21
                 call enough('hcddjkbk4.  3',bc,ibc)
                 do iz=iunit,ibcoff-1                                   11d19s21
                  bc(iz)=0d0                                            11d19s21
                 end do                                                 11d19s21
                 do iz=0,ncsfk(iarg)-1                                  11d19s21
                  iad=iunit+iz*(ncsfk(iarg)+1)                          11d19s21
                  bc(iad)=1d0                                           11d19s21
                 end do                                                 11d19s21
                 if(ndifs.eq.0.and.ndifd.eq.0)then                      2d7s23
                  ii=1                                                           11d21s20
                  do i=1,norb                                                    11d21s20
                   if(btest(ipack4(1),i))then                                    11d21s20
                    idorb(ii)=i                                                  11d21s20
                    ii=ii+1                                                      11d21s20
                   end if                                                        11d21s20
                  end do                                                         11d21s20
                  ii=1                                                           11d21s20
                  do i=1,norb                                                    11d21s20
                   if(btest(ipack4(2),i))then                                    11d21s20
                    isorb(ii)=i                                                  11d21s20
                    ii=ii+1                                                      11d21s20
                   end if                                                        11d21s20
                  end do                                                         11d21s20
                  idelta=i2smk-i2smb                                      11d5s21
                  idelta=iabs(idelta)                                     11d5s21
                  if(idelta.eq.0)then                                        9d13s21
                   if(i2sb.eq.i2sk)then                                      9d13s21
                    sum=0d0                                                  9d13s21
                    do i=1,norb                                               9d1s21
                     if(btest(itc,i))then                                    9d1s21
                      is=ism(i)                                               9d1s21
                      ig=irel(i)-1                                            9d1s21
                      ig0=ig                                                 10d8s21
                      do ipass=1,2                                            9d1s21
                       h0=-2d0                                               10d14s21
                       if(ih0n(is).gt.0)then                                 10d14s21
                        iadn=ih0n(is)+(nh0(is)+1)*ig0                        10d14s21
                        h0=bc(iadn)                                           10d14s21
                        if(ipass.eq.2.and.isopt(4).ne.0)h0=-h0          2d24s22
                        orig=sum
                        sum=sum+h0                                           10d14s21
                       end if                                                10d14s21
                       isp=ig/irefo(is)                                      10d8s21
                       ispm=1-isp                                            10d8s21
                       do j=1,norb                                            9d1s21
                        if(btest(itc,j).and.l2e.eq.0)then               2d17s22
                         js=ism(j)                                            9d1s21
                         jg=irel(j)-1                                         9d1s21
                         jg0=jg                                              10d8s21
                         do jpass=1,2                                         9d1s21
                          jsp=jg/irefo(js)                                   10d8s21
                          jspm=1-jsp                                         10d8s21
                          ltest=jsp+2*(jsp+2*(isp+2*isp))                    10d8s21
                          itestp=ltest+1                                     10d8s21
                          jtest=jspm+2*(jspm+2*(ispm+2*ispm))                    10d8s21
                          jtestp=jtest+1                                     10d8s21
                          if(iifmx(itestp).ge.0)then                    11d19s21
                           iad=i4or(js,is,is)+jg0+irefo(js)*(jg0        11d19s21
     $                         +irefo(js)*(ig0+irefo(is)*(ig0           11d19s21
     $                         +irefo(is)*iifmx(itestp))))              11d19s21
                           sum=sum+bc(iad)*0.5d0                             10d8s21
                          else if(iifmx(jtestp).ge.0)then                    10d8s21
                           iad=i4or(js,is,is)+jg0+irefo(js)*(jg0        11d19s21
     $                          +irefo(js)*(ig0+irefo(is)*(ig0          11d19s21
     $                          +irefo(is)*iifmx(jtestp))))             11d19s21
                           xint=bc(iad)                                      10d8s21
                           if(isopt(4).ne.0)then                             10d8s21
                            xint=-xint                                       10d8s21
                           end if
                           sum=sum+xint*0.5d0                                10d8s21
                          end if                                             10d8s21
                          ltest=jsp+2*(isp+2*(isp+2*jsp))                    10d8s21
                          itestp=ltest+1                                     10d8s21
                          jtest=jspm+2*(ispm+2*(ispm+2*jspm))                    10d8s21
                          jtestp=jtest+1                                     10d8s21
                          if(iifmx(itestp).ge.0)then                         10d8s21
                           iad=i4or(is,is,js)+jg0+irefo(js)*(ig0        11d19s21
     $                         +irefo(is)*(ig0+irefo(is)*(jg0           11d19s21
     $                         +irefo(js)*iifmx(itestp))))              11d19s21
                           sum=sum-bc(iad)*0.5d0                             10d8s21
                          else if(iifmx(jtestp).ge.0)then                    10d8s21
                           iad=i4or(is,is,js)+jg0+irefo(js)*(ig0        11d19s21
     $                          +irefo(is)*(ig0+irefo(is)*(jg0          11d19s21
     $                          +irefo(js)*iifmx(jtestp))))             11d19s21
                           xint=bc(iad)                                      10d8s21
                           if(isopt(4).ne.0)then                             10d8s21
                            xint=-xint                                       10d8s21
                           end if
                           sum=sum-xint*0.5d0                               10d8s21
                          end if                                             10d8s21
                          jg=jg+irefo(js)                                     9d1s21
                         end do                                               9d1s21
                        end if                                                9d1s21
                       end do                                                 9d1s21
                       ig=ig+irefo(is)                                        9d1s21
                      end do                                                  9d1s21
                     end if                                                   9d1s21
                    end do                                                    9d1s21
                    do l=1,4                                              2d24s21
                     if(nli(l).gt.0)then                                  2d24s21
                      itmpt=ibcoff                                        2d24s21
                      itmpm=itmpt+ncsfk2(l,iarg)*nli(l)                   11d19s21
                      ibcoff=itmpm+nli(l)*nli(l)                          2d24s21
                      call enough('hcddjkbk4.  4',bc,ibc)
                      iad1=ivcv+iff22k(ivcv+5+l)                          8d26s21
                      iad2=iad1+nli(l)                                    3d19s21
                      do iii=0,nli(l)-1                                   2d24s21
                       do j=0,ncsfk2(l,iarg)-1                             2d24s21
                        ji=iad2+j+ncsfk2(l,iarg)*iii                       2d24s21
                        ij=itmpt+iii+nli(l)*j                             2d24s21
                        bc(ij)=ff22k(ji)                                  8d26s21
                       end do                                             2d24s21
                      end do                                              2d24s21
                      call dgemm('n','n',nli(l),nli(l),ncsfk2(l,iarg),  11d19s21
     $                     sum,                                         11d19s21
     $                  bc(itmpt),nli(l),ff22k(iad2),ncsfk2(l,iarg),0d0,11d19s21
     $                   bc(itmpm),nli(l),                              2d24s21
     d' hcddjkbk4.  1')
                      jtmpm=itmpm                                         2d24s21
                      do j=0,nli(l)-1                                     2d24s21
                       jj=iff22k(iad1+j)-1                                8d26s21
                       jden=iden1en(l,l)+nfdatk(2,l,isb)*jj-1           11d23s21
                       do k=0,nli(l)-1                                    2d24s21
                        kk=iff22k(iad1+k)                                 8d26s21
                        bc(jden+kk)=bc(jden+kk)+bc(jtmpm+k)               2d24s21
                       end do                                             2d24s21
                       jtmpm=jtmpm+nli(l)                                 2d24s21
                      end do                                              2d24s21
                      ibcoff=itmpt                                        2d24s21
                     end if                                               2d24s21
                    end do                                                2d24s21
                   end if                                                 11d19s21
                  end if                                                  11d19s21
                  nx1=0                                                     9d15s21
                  do i=1,norb                                               9d15s21
                   if(btest(ito,i))then                                  11d19s21
                    nx1=nx1+1                                               9d15s21
                    imap(nx1)=i                                            9d15s21
                   end if                                                   9d15s21
                  end do                                                    9d15s21
                  nx1=nx1+1                                              11d5s21
                  imap(nx1)=norbx                                        11d19s21
                  nx1=nx1+1                                              11d5s21
                  imap(nx1)=norbxx                                       11d19s21
                  if(idelta.le.2)then                                        9d13s21
                   sigc=1d0                                                  9d21s21
                   ntype1a=0                                                 10d14s21
                   ibcb4=ibcoff                                           11d5s21
                   call getcup(irw0,nopenip,i2sb,i2smb,i2sk,i2smk,      11d26s21
     $                  ntypeg,ioutg,imatg,mcsf,imap,bc,ibc)            11d14s22
                   if(ntypeg.gt.0)then
                    ncsfa=mcsf(1)                                            10d18s21
                    ncsfc=mcsf(2)                                            10d19s21
                    do iq=0,ntypeg-1
                     ipackc=ibc(ioutg+iq)                                    10d19s21
                     if(ipackc1(1).gt.0)then                                 10d19s21
                      mmm=ipackc1(1)                                         10d19s21
                      if(mmm.le.norb)then                                 11d5s21
                       isy(1)=ism(mmm)                                  11d29s21
                       igya(1)=irel(mmm)-1                                    10d14s21
                      else
                       igya(1)=-1                                         11d5s21
                      end if                                              11d5s21
                      imy(1)=0                                               10d14s21
                     else                                                    10d14s21
                      mmm=-ipackc1(1)                                        10d19s21
                      if(mmm.le.norb)then                                 11d5s21
                       isy(1)=ism(mmm)                                  11d29s21
                       igya(1)=irel(mmm)-1                                    10d14s21
                      else                                                11d5s21
                       igya(1)=-1                                         11d5s21
                      end if                                              11d5s21
                      imy(1)=1                                               10d14s21
                     end if                                                  10d14s21
                     if(ipackc1(2).gt.0)then                                 10d19s21
                      mmm=ipackc1(2)                                         10d19s21
                      if(mmm.le.norb)then                                 11d5s21
                       isy(2)=ism(mmm)                                  11d29s21
                       igya(2)=irel(mmm)-1                                    10d14s21
                      else                                                11d5s21
                       igya(2)=-1                                         11d5s21
                      end if                                              11d5s21
                      imy(2)=0                                               10d14s21
                     else                                                    10d14s21
                      mmm=-ipackc1(2)                                        10d19s21
                      if(mmm.le.norb)then                                 11d5s21
                       isy(2)=ism(mmm)                                  11d29s21
                       igya(2)=irel(mmm)-1                                    10d14s21
                      else                                                11d5s21
                       igya(2)=-1                                         11d5s21
                      end if                                              11d5s21
                      imy(2)=1                                               10d14s21
                     end if                                                  10d14s21
                     imatq=imatg+mcsf(1)*mcsf(2)*iq
                     if(min(igya(1),igya(2)).ge.0)then                  11d29s21
                      if(ih0n(isy(1)).gt.0)then                         11d29s21
                       iad=ih0n(isy(1))+igya(1)+nh0(isy(1))*igya(2)     11d29s21
                       sum=bc(iad)                                      11d29s21
                       if(imy(1).ne.0.and.isopt(4).ne.0)sum=-sum        11d29s21
                      else                                              11d29s21
                       iad=ih0n(isy(2))+igya(2)+nh0(isy(2))*igya(1)     11d29s21
                       sum=bc(iad)                                      11d29s21
                       if(isopt(2).ne.0)sum=-sum                        11d29s21
                       if(imy(2).ne.0.and.isopt(4).ne.0)sum=-sum        11d29s21
                      end if                                            11d29s21
                      do i=1,norb                                       11d29s21
                       if((btest(itc,i).or.btest(ito,i))                2d17s22
     $                      .and.l2e.eq.0)then                          2d17s22
                        js=ism(i)                                       11d29s21
                        jg=irel(i)-1                                    11d29s21
                        do jpass=0,1                                    11d29s21
                         phs=-1d0                                        11d29s21
                         if(btest(ito,i))phs=-0.5d0                      11d29s21
                         ltest=jpass+2*(imy(1)+2*(imy(2)+2*jpass))      11d29s21
                         jtest=1-jpass+2*(1-imy(1)+2*(1-imy(2)          11d29s21
     $                        +2*(1-jpass)))                            11d29s21
                         itestp=ltest+1                                 11d29s21
                         jtestp=jtest+1                                 11d29s21
                         juse=-1                                        11d29s21
                         if(iifmx(itestp).ge.0)then                     11d29s21
                          juse=iifmx(itestp)                            11d29s21
                         else if(iifmx(jtestp).ge.0)then                11d29s21
                          juse=iifmx(jtestp)                            11d29s21
                          if(isopt(4).ne.0)phs=-phs                     11d29s21
                         end if                                         11d29s21
                         if(juse.ge.0)then                              11d29s21
                          iad=i4or(isy(1),isy(2),js)+jg+irefo(js)*(     11d29s21
     $                        igya(1)+irefo(isy(1))*(igya(2)            11d29s21
     $                        +irefo(isy(2))*(jg+irefo(js)*juse)))      11d29s21
                          orig=sum
                          sum=sum+phs*bc(iad)                           11d29s21
                         end if                                         11d29s21
                        end do                                          11d29s21
                       end if                                           11d29s21
                      end do                                            11d29s21
                      do lk=1,4                                          11d26s21
                       if(nli(lk).gt.0)then                              11d26s21
                        iad1=ivcv+iff22k(ivcv+5+lk)                        11d26s21
                        iad2=iad1+nli(lk)                                  11d26s21
                        itmpx=ibcoff                                       11d26s21
                        itmpy=itmpx+ncsfb(jarg)*nli(lk)                    11d26s21
                        ibcoff=itmpy+ncsfb(jarg)*nli(lk)                   11d26s21
                        call enough('hcddjkbk4.  5',bc,ibc)
                        call dgemm('n','n',ncsfb(jarg),nli(lk),
     $                      ncsfk2(lk,iarg),sum,bc(imatq),ncsfb(jarg),  11d29s21
     $                    ff22k(iad2),ncsfk2(lk,iarg),0d0,bc(itmpx),    11d26s21
     $                    ncsfb(jarg),                                  11d26s21
     d' hcddjkbk4.  2')
                        do j=0,ncsfb(jarg)-1                               11d26s21
                         do k=0,nli(lk)-1                                  11d26s21
                          jk=itmpx+j+ncsfb(jarg)*k                         11d26s21
                          kj=itmpy+k+nli(lk)*j                             11d26s21
                          bc(kj)=bc(jk)                                    11d26s21
                         end do                                            11d26s21
                        end do                                             11d26s21
                        jtmpy=itmpy                                        11d26s21
                        do lb=1,4                                          11d26s21
                         if(nlj(lb).gt.0)then                              11d26s21
                          jad1=jvcv+iff22b(jvcv+5+lb)                          8d26s21
                          jad2=jad1+nlj(lb)                                    3d19s21
                          itmpz=ibcoff                                     11d26s21
                          ibcoff=itmpz+nli(lk)*nlj(lb)                     11d26s21
                          call enough('hcddjkbk4.  6',bc,ibc)
                          call dgemm('n','n',nli(lk),nlj(lb),
     $                            ncsfb2(lb,jarg),1d0,bc(jtmpy),nli(lk),      11d26s21
     $                      ff22b(jad2),ncsfb2(lb,jarg),0d0,bc(itmpz),  11d26s21
     $                      nli(lk),                                    11d26s21
     d' hcddjkbk4.  3')
                          jtmp=itmpz                                       11d26s21
                          do iii=0,nlj(lb)-1                             11d29s21
                           jj=iff22b(jad1+iii)-1                         11d29s21
                           jjden=iden1en(lk,lb)                          11d29s21
     $                          +nfdatk(2,lk,isb)*jj-1                   11d29s21
                           do j=0,nli(lk)-1                              11d29s21
                            kk=iff22k(iad1+j)                            11d29s21
                            bc(jjden+kk)=bc(jjden+kk)+bc(jtmp+j)         11d29s21
                           end do                                             2d24s21
                           jtmp=jtmp+nli(lk)                                  11d19s21
                          end do                                              2d24s21
                          ibcoff=itmpz                                     11d26s21
                         end if                                            11d26s21
                         jtmpy=jtmpy+nli(lk)*ncsfb2(lb,jarg)               11d26s21
                        end do                                             11d26s21
                        ibcoff=itmpx                                       11d26s21
                       end if                                            11d26s21
                       imatq=imatq+ncsfb(jarg)*ncsfk2(lk,iarg)           11d26s21
                      end do                                             11d26s21
                     end if                                             11d29s21
                    end do                                                   9d20s21
                   end if                                                    9d20s21
                   ibcoff=ibcb4                                           11d5s21
                  end if                                                     9d13s21
                  if(idelta.le.4.and.l2e.eq.0)then                      2d17s22
c
c     intermediate state = bra or ket
c
                   call spinloop1(i2sb,i2smb,i2sk,i2smk,nopenip,          11d19s21
     $                ncsfb(jarg),imap,ncsfk(iarg),itypeg,imatg,ntypeg, 11d19s21
     $                irw0,bc(iunit),ncsfk(iarg),ieoro,bc,ibc)          11d14s22
                   do ii=0,ntypeg-1
                    ipack=ibc(itypeg+ii)
                    jmat=imatg+ncsfb(jarg)*ncsfk(iarg)*ii               11d29s21
                    do j=1,4                                                 10d20s21
                     if(ipack2(j).gt.0)then                             11d29s21
                      imy(j)=0                                          11d29s21
                      ipack2a(j)=ipack2(j)                              11d29s21
                     else                                               11d29s21
                      imy(j)=1                                          11d29s21
                      ipack2a(j)=-ipack2(j)                             11d29s21
                     end if                                             11d29s21
                     if(ipack2a(j).le.norb)then                           11d5s21
                      isy(j)=ism(ipack2a(j))                            11d29s21
                      igya(j)=irel(ipack2a(j))-1                        11d29s21
                     else                                               11d29s21
                      igya(j)=-1                                        11d29s21
                     end if                                                  10d20s21
                    end do                                                   10d20s21
                    if(min(igya(1),igya(2),igya(3),igya(4)).ge.0)then   11d29s21
                     ltest=imy(1)+2*(imy(2)+2*(imy(3)+2*imy(4)))        11d19s21
                     jtest=1-imy(1)+2*(1-imy(2)+2*(1-imy(3)             11d19s21
     $                    +2*(1-imy(4))))                               11d19s21
                     itestp=ltest+1                                           10d20s21
                     jtestp=jtest+1                                           10d20s21
                     iuse=-1                                               11d5s21
                     phs=0.5d0                                          11d29s21
                     if(iifmx(itestp).ge.0)then                               10d20s21
                      iuse=iifmx(itestp)                                   11d5s21
                     else if(iifmx(jtestp).ge.0)then                       11d5s21
                      iuse=iifmx(jtestp)                                   11d5s21
                      if(isopt(4).ne.0)phs=-phs                            11d5s21
                     end if                                                11d5s21
                     if(iuse.ge.0)then                                  11d29s21
                      iad=i4or(isy(2),isy(3),isy(4))+igya(1)             11d29s21
     $                    +irefo(isy(1))*(igya(2)+irefo(isy(2))*(igya(3)11d29s21
     $                    +irefo(isy(3))*(igya(4)+irefo(isy(4))*iuse))) 11d29s21
                      xint=bc(iad)*phs                                    11d5s21
                      do lk=1,4                                          11d29s21
                       if(nli(lk).gt.0)then                              11d29s21
                        iad1=ivcv+iff22k(ivcv+5+lk)                        11d26s21
                        iad2=iad1+nli(lk)                                  11d26s21
                        itmpx=ibcoff                                       11d26s21
                        itmpy=itmpx+ncsfb(jarg)*nli(lk)                    11d26s21
                        ibcoff=itmpy+ncsfb(jarg)*nli(lk)                   11d26s21
                        call enough('hcddjkbk4.  7',bc,ibc)
                        call dgemm('n','n',ncsfb(jarg),nli(lk),
     $                       ncsfk2(lk,iarg),xint,bc(jmat),ncsfb(jarg),    11d29s21
     $                    ff22k(iad2),ncsfk2(lk,iarg),0d0,bc(itmpx),    11d26s21
     $                    ncsfb(jarg),                                  11d26s21
     d' hcddjkbk4.  4')
                        do j=0,ncsfb(jarg)-1                               11d26s21
                         do k=0,nli(lk)-1                                  11d26s21
                          jk=itmpx+j+ncsfb(jarg)*k                         11d26s21
                          kj=itmpy+k+nli(lk)*j                             11d26s21
                          bc(kj)=bc(jk)                                    11d26s21
                         end do                                            11d26s21
                        end do                                             11d26s21
                        jtmpy=itmpy                                        11d26s21
                        do lb=1,4                                          11d26s21
                         if(nlj(lb).gt.0)then                              11d26s21
                          jad1=jvcv+iff22b(jvcv+5+lb)                          8d26s21
                          jad2=jad1+nlj(lb)                                    3d19s21
                          itmpz=ibcoff                                     11d26s21
                          ibcoff=itmpz+nli(lk)*nlj(lb)                     11d26s21
                          call enough('hcddjkbk4.  8',bc,ibc)
                          call dgemm('n','n',nli(lk),nlj(lb),
     $                            ncsfb2(lb,jarg),1d0,bc(jtmpy),nli(lk),      11d26s21
     $                      ff22b(jad2),ncsfb2(lb,jarg),0d0,bc(itmpz),  11d26s21
     $                      nli(lk),                                    11d26s21
     d' hcddjkbk4.  5')
                          jtmp=itmpz                                       11d26s21
                          do iii=0,nlj(lb)-1                                   2d24s21
                           jj=iff22b(jad1+iii)-1                              8d26s21
                           jjden=iden1en(lk,lb)                          11d29s21
     $                              +nfdatk(2,lk,isb)*jj-1              11d26s21
                           do j=0,nli(lk)-1                                    2d24s21
                            kk=iff22k(iad1+j)                                 8d26s21
                            bc(jjden+kk)=bc(jjden+kk)+bc(jtmp+j)             11d19s21
                           end do                                             2d24s21
                           jtmp=jtmp+nli(lk)                                  11d19s21
                          end do                                              2d24s21
                          ibcoff=itmpz                                     11d26s21
                         end if                                            11d26s21
                         jtmpy=jtmpy+nli(lk)*ncsfb2(lb,jarg)               11d26s21
                        end do                                             11d26s21
                        ibcoff=itmpx                                       11d26s21
                       end if                                            11d29s21
                       jmat=jmat+ncsfb(jarg)*ncsfk2(lk,iarg)             11d29s21
                      end do                                             11d29s21
                     end if                                             11d29s21
                    end if                                                11d5s21
 2021               continue                                              11d9s21
                   end do                                                    10d20s21
                   ibcoff=itypeg                                          11d19s21
                   do i1=1,norbx                                          11d19s21
                    if(btest(ito,i1))then                                   9d15s21
                     do i2=i1+1,norbxx                                    11d19s21
                      if(btest(ito,i2))then                                 9d15s21
                       itestc=itc                                              5d7s21
                       itesto=ito                                              5d7s21
                       if(btest(itestc,i1))then                               9d13s21
                        itestc=ibclr(itestc,i1)                               9d13s21
                        itesto=ibset(itesto,i1)                               9d13s21
                        nopenkk=nopenip+1                                           1d22s21
                       else if(btest(itesto,i1))then                          9d13s21
                        itesto=ibclr(itesto,i1)                               9d13s21
                        nopenkk=nopenip-1                                              11d13s20
                       else                                                          11d13s20
                        write(6,*)('bita not set for i1 = '),i1
                        stop 'nab4a(1,1)'                                          11d27s20
                       end if                                                        11d13s20
                       if(btest(itesto,i2))then                               9d13s21
                        itestc=ibset(itestc,i2)                               9d13s21
                        itesto=ibclr(itesto,i2)                               9d13s21
                        nopenkk=nopenkk-1                                              11d13s20
                       else if(btest(itestc,i2))then                         9d13s21
                        go to 1848
                       else                                                          11d13s20
                        itesto=ibset(itesto,i2)                               9d13s21
                        nopenkk=nopenkk+1                                              11d13s20
                       end if                                                        11d13s20
                       call gandcr(jtc,jto,itestc,itesto,nopenip,         11d19s21
     $               nopenkk,norbxx,nnot1,nab1,icode,imap,nx1,irw1,     9d16s21
     $                 irw2,iwpb1,iwpk1,bc,ibc)                         11d14s22
                       call gandcr(itestc,itesto,itc,ito,nopenkk,       11d23s21
     $                      nopenip,norbxx,nnot2,nab2,icode,imap,nx2,   11d23s21
     $                      irw1,irw2,iwpb2,iwpk2,bc,ibc)               11d14s22
                       call spinloop(i2sb,i2smb,i2sk,i2smk,nopenip,       11d19s21
     $                    nopenjp,nopenkk,ncsfb(jarg),ncsfk(iarg),itype,11d19s21
     $                    imatx,ntypeg,nab1,iwpb1,iwpk1,nab2,iwpb2,     11d19s21
     $                    iwpk2,bc(iunit),ncsfk(iarg),ieoro,bc,ibc)     11d14s22
                       do ii=0,ntypeg-1                                    11d5s21
                        ipack=ibc(itype+ii)
c     we have 4o or (1v|v'4), or 4v
                        nx=0                                            11d26s21
                        do j=1,4                                        11d26s21
                         if(ipack2(j).gt.0)then                         11d26s21
                          imy(j)=0                                      11d26s21
                          ipack2a(j)=ipack2(j)                          11d26s21
                         else                                           11d26s21
                          imy(j)=1                                      11d26s21
                          ipack2a(j)=-ipack2(j)                         11d26s21
                         end if                                         11d26s21
                         if(ipack2a(j).le.norb)then                     11d26s21
                          isy(j)=ism(ipack2a(j))                        11d26s21
                          igya(j)=irel(ipack2a(j))-1                    11d26s21
                         else                                           11d26s21
                          igya(j)=-1                                    11d26s21
                          nx=nx+1                                       11d26s21
                         end if                                         11d26s21
                        end do                                          11d26s21
                        ltest=imy(1)+2*(imy(2)+2*(imy(3)+2*imy(4)))     11d26s21
                        jtest=1-imy(1)+2*(1-imy(2)+2*(1-imy(3)          11d26s21
     $                       +2*(1-imy(4))))                            11d26s21
                        itestp=ltest+1                                  11d26s21
                        jtestp=jtest+1                                  11d26s21
                        phs=1d0                                         11d26s21
                        juse=-1                                         11d26s21
                        if(iifmx(itestp).ge.0)then                      11d26s21
                         juse=iifmx(itestp)                             11d26s21
                        else if(iifmx(jtestp).ge.0)then                 11d26s21
                         juse=iifmx(jtestp)                             11d26s21
                         if(isopt(4).ne.0)phs=-phs                      11d26s21
                        end if                                          11d26s21
                        if(juse.ge.0.and.nx.le.2)then                   11d26s21
                         if(nx.eq.0)then                                 11d26s21
                          iad=i4or(isy(2),isy(3),isy(4))+igya(1)        11d26s21
     $                        +irefo(isy(1))*(igya(2)+irefo(isy(2))*(   11d26s21
     $                        igya(3)+irefo(isy(3))*(igya(4)            11d26s21
     $                        +irefo(isy(4))*juse)))                    11d26s21
                          phs=phs*bc(iad)                               11d26s21
                         else if(nx.eq.2)then                            11d26s21
                          jcol=igya(1)+irefo(isy(1))*(igya(4)           11d26s21
     $                         +irefo(isy(4))*juse)                     11d26s21
                         end if                                          11d26s21
                         jmat=imatx+ii*ncsfb(jarg)*ncsfk(iarg)           11d26s21
                         do lk=1,4                                       11d26s21
                          if(nli(lk).gt.0)then                           11d26s21
                           iad1=ivcv+iff22k(ivcv+5+lk)                        11d26s21
                           iad2=iad1+nli(lk)                                  11d26s21
                           itmpx=ibcoff                                       11d26s21
                           itmpy=itmpx+ncsfb(jarg)*nli(lk)                    11d26s21
                           ibcoff=itmpy+ncsfb(jarg)*nli(lk)                   11d26s21
                           call enough('hcddjkbk4.  9',bc,ibc)
                           call dgemm('n','n',ncsfb(jarg),nli(lk),
     $                    ncsfk2(lk,iarg),phs,bc(jmat),ncsfb(jarg),     11d29s21
     $                    ff22k(iad2),ncsfk2(lk,iarg),0d0,bc(itmpx),    11d26s21
     $                    ncsfb(jarg),                                  11d26s21
     d' hcddjkbk4.  6')
                           do j=0,ncsfb(jarg)-1                               11d26s21
                            do k=0,nli(lk)-1                                  11d26s21
                             jk=itmpx+j+ncsfb(jarg)*k                         11d26s21
                             kj=itmpy+k+nli(lk)*j                             11d26s21
                             bc(kj)=bc(jk)                                    11d26s21
                            end do                                            11d26s21
                           end do                                             11d26s21
                           jtmpy=itmpy                                        11d26s21
                           do lb=1,4                                          11d26s21
                            if(nlj(lb).gt.0)then                              11d26s21
                             jad1=jvcv+iff22b(jvcv+5+lb)                          8d26s21
                             jad2=jad1+nlj(lb)                                    3d19s21
                             itmpz=ibcoff                                     11d26s21
                             ibcoff=itmpz+nli(lk)*nlj(lb)                     11d26s21
                             call enough('hcddjkbk4. 10',bc,ibc)
                             call dgemm('n','n',nli(lk),nlj(lb),
     $                            ncsfb2(lb,jarg),1d0,bc(jtmpy),nli(lk),      11d26s21
     $                      ff22b(jad2),ncsfb2(lb,jarg),0d0,bc(itmpz),  11d26s21
     $                      nli(lk),                                    11d26s21
     d' hcddjkbk4.  7')
                             jtmp=itmpz                                       11d26s21
                             if(nx.eq.0)then                            11d26s21
                              do iii=0,nlj(lb)-1                                   2d24s21
                               jj=iff22b(jad1+iii)-1                              8d26s21
                               jjden=iden1en(lk,lb)                     11d26s21
     $                              +nfdatk(2,lk,isb)*jj-1              11d26s21
                               do j=0,nli(lk)-1                                    2d24s21
                                kk=iff22k(iad1+j)                                 8d26s21
                                bc(jjden+kk)=bc(jjden+kk)+bc(jtmp+j)             11d19s21
                               end do                                             2d24s21
                               jtmp=jtmp+nli(lk)                                  11d19s21
                              end do                                              2d24s21
                             end if                                     11d26s21
                             ibcoff=itmpz                                     11d26s21
                            end if                                            11d26s21
                            jtmpy=jtmpy+nli(lk)*ncsfb2(lb,jarg)               11d26s21
                           end do                                             11d26s21
                           ibcoff=itmpx                                       11d26s21
                          end if                                         11d26s21
                          jmat=jmat+ncsfb(jarg)*ncsfk2(lk,iarg)          11d26s21
                         end do                                          11d26s21
                        end if                                          11d26s21
                       end do                                             11d19s21
                       ibcoff=itype                                       11d19s21
 1848                  continue                                         11d22s21
                      end if                                              11d19s21
                     end do                                               11d19s21
                    end if                                                11d19s21
                   end do                                                 11d19s21
                  end if                                                  11d19s21
                 else if(ndifs.eq.2.and.ndifb.eq.2)then                  2d6s23
                  do i=1,norbxx
                   if(btest(gandco,i))then
                    if((btest(jtc,i).and.btest(ito,i)).or.                 10d21s22
     $                   (btest(jto,i).and..not.btest(itc,i)))then            10d21s22
                     nab4(1,1)=i
                    else
                     nab4(2,1)=i
                    end if
                   end if
                  end do
                  call gandcr(jtc,jto,itc,ito,nopenjp,nopenip,          11d19s21
     $               norbxx,nnot1b,nab1b,icode,imap,nx,irw1,irw2,       11d19s21
     $                 iwpb1,iwpk1,bc,ibc)                              11d14s22
                  call gencup(i2sb,i2smb,i2sk,i2smk,nopenjp,nopenip,    11d19s21
     $                 nab1b,iwpb1,iwpk1,ioutg,imatg,ntypeg,mcsf,       11d19s21
     $                 bc(iunit),ncsfk(iarg),ncsfk(iarg),bc,ibc)        11d14s22
                  do ii=0,ntypeg-1                                      11d19s21
                   ipackc=ibc(ioutg+ii)                                 11d19s21
                   imatu=imatg+ncsfb(jarg)*ncsfk(iarg)*ii               11d23s21
                   if(ipackc1(1).gt.0)then                                   9d10s21
                    isy(3)=ism(ipackc1(1))                              11d19s21
                    igya(3)=irel(ipackc1(1))-1                          11d19s21
                    imy(3)=0                                                 10d12s21
                   else
                    isy(3)=ism(-ipackc1(1))                             11d19s21
                    igya(3)=irel(-ipackc1(1))-1                              10d12s21
                    imy(3)=1                                                 10d12s21
                   end if
                   if(ipackc1(2).gt.0)then                                   9d10s21
                    isy(4)=ism(ipackc1(2))                              11d19s21
                    igya(4)=irel(ipackc1(2))-1                          11d19s21
                    imy(4)=0                                                 10d12s21
                   else
                    isy(4)=ism(-ipackc1(2))                             11d19s21
                    igya(4)=irel(-ipackc1(2))-1                              10d12s21
                    imy(4)=1                                                 10d12s21
                   end if
                   h0=-2d0                                                   10d13s21
                   sum=0d0                                                   10d14s21
                   if(ih0n(isy(3)).gt.0)then                            11d19s21
                    iad=ih0n(isy(3))+igya(3)+nh0(isy(3))*igya(4)        11d19s21
                    h0=bc(iad)                                                10d13s21
                    if(imy(3).ne.0.and.isopt(4).ne.0)h0=-h0                   10d13s21
                    sum=h0                                                   10d14s21
                   else if(ih0n(isy(4)).gt.0)then                       11d19s21
                    iad=ih0n(isy(4))+igya(4)+nh0(isy(4))*igya(3)        11d19s21
                    h0=bc(iad)                                                10d13s21
                    if(isopt(2).ne.0)h0=-h0                                  10d13s21
                    if(imy(4).ne.0.and.isopt(4).ne.0)h0=-h0                  10d13s21
                    sum=h0                                                   10d14s21
                   end if                                                    10d13s21
                   do i=1,norb                                               9d10s21
                    if(btest(jtc,i).and.btest(itc,i).and.l2e.eq.0)then  2d17s22
                     js=ism(i)                                               9d10s21
                     jg=irel(i)-1                                            9d10s21
                     jg0=jg                                                  10d12s21
                     do jpass=1,2                                            9d10s21
                      ltest=jpass-1+2*(jpass-1+2*(imy(3)+2*imy(4)))          10d12s21
                      itestp=ltest+1                                         10d12s21
                      jtest=2-jpass+2*(2-jpass+2*(1-imy(3)+2*           11d19s21
     $                     (1-imy(4))))                                 11d19s21
                      jtestp=jtest+1                                         10d12s21
                      xintj=-2d0                                             10d12s21
                      if(iifmx(itestp).ge.0)then                             10d12s21
                       iad=i4or(js,isy(3),isy(4))+jg0+irefo(js)*(jg0    11d19s21
     $                     +irefo(js)*(igya(3)+irefo(isy(3))*(igya(4)   11d19s21
     $                     +irefo(isy(4))*iifmx(itestp))))              11d19s21
                       xintj=bc(iad)                                         10d12s21
                       orig=sum
                       sum=sum+xintj                                         10d12s21
                      else if(iifmx(jtestp).ge.0)then                        10d12s21
                       iad=i4or(js,isy(3),isy(4))+jg0+irefo(js)*(jg0    11d19s21
     $                      +irefo(js)*(igya(3)+irefo(isy(3))*(igya(4)  11d19s21
     $                      +irefo(isy(4))*iifmx(jtestp))))             11d19s21
                       xintj=bc(iad)                                         10d12s21
                       if(isopt(4).ne.0)xintj=-xintj                         10d12s21
                       orig=sum
                       sum=sum+xintj                                         10d12s21
                      end if                                                 10d12s21
                      ltest=jpass-1+2*(imy(4)+2*(imy(3)+2*(jpass-1)))        10d12s21
                      itestp=ltest+1                                         10d12s21
                      jtest=2-jpass+2*(1-imy(4)+2*(1-imy(3)             11d19s21
     $                     +2*(2-jpass)))                               11d19s21
                      jtestp=jtest+1                                         10d12s21
                      xintj=-2d0                                             10d12s21
                      if(iifmx(itestp).ge.0)then                             10d12s21
                       iad=i4or(isy(4),isy(3),js)+jg0+irefo(js)*(igya(4)11d19s21
     $                +irefo(isy(4))*(igya(3)+irefo(isy(3))*(jg0        11d19s21
     $                     +irefo(js)*iifmx(itestp))))                  11d19s21
                       xintk=bc(iad)                                         10d12s21
                       orig=sum
                       sum=sum-xintk                                         10d12s21
                      else if(iifmx(jtestp).ge.0)then                        10d12s21
                       iad=i4or(isy(4),isy(3),js)+jg0+irefo(js)*(igya(4)11d19s21
     $                 +irefo(isy(4))*(igya(3)+irefo(isy(3))*(jg0       11d19s21
     $                      +irefo(js)*iifmx(jtestp))))                 11d19s21
                       xintk=bc(iad)                                         10d12s21
                       if(isopt(4).ne.0)xintk=-xintk                         10d12s21
                       orig=sum
                       sum=sum-xintk                                         10d12s21
                      end if                                                 10d12s21
                      jg=jg+irefo(js)                                        9d10s21
                     end do                                                  9d10s21
                    end if                                                   9d10s21
                   end do                                                    9d10s21
                   if(ipackc1(1).gt.0)then                                   10d12s21
                    in=ipackc1(1)                                            10d12s21
                    imy(1)=1                                                 10d12s21
                   else                                                      9d10s21
                    in=-ipackc1(1)                                           10d12s21
                    imy(1)=0                                                 10d12s21
                   end if                                                    9d13s21
                   if(btest(jtc,in).and.l2e.eq.0)then                   2d17s22
                    ltest=imy(1)+2*(imy(1)+2*(imy(3)+2*imy(4)))              10d12s21
                    itestp=ltest+1                                           10d12s21
                    jtest=1-imy(1)+2*(1-imy(1)+2*(1-imy(3)              11d19s21
     $                   +2*(1-imy(4))))                                11d19s21
                    jtestp=jtest+1                                           10d12s21
                    xintj=-2d0                                               10d12s21
                    if(iifmx(itestp).ge.0)then                               10d12s21
                     iad=i4or(isy(3),isy(3),isy(4))+igya(3)             11d19s21
     $                   +irefo(isy(3))*(igya(3)+irefo(isy(3))*(igya(3) 11d19s21
     $                   +irefo(isy(3))*(igya(4)+irefo(isy(4))          11d19s21
     $                  *iifmx(itestp))))                                   10d12s21
                     xintj=bc(iad)                                           10d12s21
                     orig=sum
                     sum=sum+xintj                                           10d13s21
                    else if(iifmx(jtestp).ge.0)then                          10d12s21
                     iad=i4or(isy(3),isy(3),isy(4))+igya(3)             11d19s21
     $                    +irefo(isy(3))*(igya(3)+irefo(isy(3))*(igya(3)11d19s21
     $                    +irefo(isy(3))*(igya(4)+irefo(isy(4))         11d19s21
     $                    *iifmx(jtestp))))                             11d19s21
                     xintj=bc(iad)                                           10d12s21
                     if(isopt(4).ne.0)xintj=-xintj                           10d12s21
                     orig=sum
                     sum=sum+xintj                                           10d13s21
                    end if                                                   10d12s21
                    ltest=imy(1)+2*(imy(4)+2*(imy(3)+2*imy(1)))              10d12s21
                    itestp=ltest+1                                           10d12s21
                    jtest=1-imy(1)+2*(1-imy(4)+2*(1-imy(3)+2*           11d19s21
     $                   (1-imy(1))))                                   11d19s21
                    jtestp=jtest+1                                           10d12s21
                    xintk=-2d0                                               10d12s21
                    if(iifmx(itestp).ge.0)then                               10d12s21
                     iad=i4or(isy(4),isy(3),isy(3))+igya(3)             11d19s21
     $                   +irefo(isy(3))*(igya(4)+irefo(isy(4))*(igya(3) 11d19s21
     $                   +irefo(isy(3))*(igya(3)+irefo(isy(3))          11d19s21
     $                 *iifmx(itestp))))                                   10d12s21
                     xintk=bc(iad)                                           10d12s21
                     orig=sum
                     sum=sum-xintk                                           10d13s21
                    else if(iifmx(jtestp).ge.0)then                          10d12s21
                     iad=i4or(isy(4),isy(3),isy(3))+igya(3)             11d19s21
     $                    +irefo(isy(3))*(igya(4)+irefo(isy(4))*(igya(3)11d19s21
     $                    +irefo(isy(3))*(igya(3)+irefo(isy(3))         11d19s21
     $                    *iifmx(jtestp))))                             11d19s21
                     xintk=bc(iad)                                           10d12s21
                     if(isopt(4).ne.0)xintk=-xintk                           10d12s21
                     orig=sum
                     sum=sum-xintk                                           10d13s21
                    end if                                                   10d12s21
                   end if                                                    9d13s21
                   if(ipackc1(2).gt.0)then                                   10d12s21
                    im=ipackc1(2)                                            10d12s21
                    imy(2)=1                                                 10d13s21
                   else                                                      10d12s21
                    im=-ipackc1(2)                                           10d12s21
                    imy(2)=0                                                 10d13s21
                   end if                                                    10d12s21
                   if(btest(itc,im).and.l2e.eq.0)then                   2d17s22
                    ltest=imy(2)+2*(imy(2)+2*(imy(3)+2*imy(4)))              10d12s21
                    itestp=ltest+1                                           10d12s21
                    jtest=1-imy(2)+2*(1-imy(2)+2*(1-imy(3)              11d22s21
     $                   +2*(1-imy(4))))                                11d22s21
                    jtestp=jtest+1                                           10d12s21
                    xintj=-2d0                                               10d12s21
                    if(iifmx(itestp).ge.0)then                          11d19s21
                     iad=i4or(isy(4),isy(3),isy(4))+igya(4)             11d19s21
     $                   +irefo(isy(4))*(igya(4)+irefo(isy(4))*(igya(3) 11d19s21
     $                   +irefo(isy(3))*(igya(4)+irefo(isy(4))          11d19s21
     $                   *iifmx(itestp))))                              11d19s21
                     xintj=bc(iad)                                           10d12s21
                     orig=sum
                     sum=sum+xintj                                           10d13s21
                    else if(iifmx(jtestp).ge.0)then                          10d12s21
                     iad=i4or(isy(4),isy(3),isy(4))+igya(4)             11d19s21
     $                    +irefo(isy(4))*(igya(4)+irefo(isy(4))*(igya(3)11d19s21
     $                    +irefo(isy(3))*(igya(4)+irefo(isy(4))         11d19s21
     $                    *iifmx(jtestp))))                             11d19s21
                     xintj=bc(iad)                                           10d12s21
                     if(isopt(4).ne.0)xintj=-xintj                           10d12s21
                     orig=sum
                     sum=sum+xintj                                           10d13s21
                    end if                                                   10d12s21
                    ltest=imy(2)+2*(imy(4)+2*(imy(3)+2*imy(2)))              10d12s21
                    itestp=ltest+1                                           10d12s21
                    jtest=1-imy(2)+2*(1-imy(4)+2*(1-imy(3)              11d19s21
     $                   +2*(1-imy(2))))                                11d19s21
                    jtestp=jtest+1                                           10d12s21
                    xintk=-2d0                                               10d12s21
                    if(iifmx(itestp).ge.0)then                               10d12s21
                     iad=i4or(isy(4),isy(3),isy(4))+igya(4)             11d19s21
     $                   +irefo(isy(4))*(igya(4)+irefo(isy(4))*(igya(3) 11d19s21
     $                   +irefo(isy(3))*(igya(4)+irefo(isy(4))          11d19s21
     $              *iifmx(itestp))))                                   10d12s21
                     xintk=bc(iad)                                           10d12s21
                     orig=sum
                     sum=sum-xintk                                           10d13s21
                    else if(iifmx(jtestp).ge.0)then                          10d12s21
                     iad=i4or(isy(4),isy(3),isy(4))+igya(4)             11d19s21
     $                    +irefo(isy(4))*(igya(4)+irefo(isy(4))*(igya(3)11d19s21
     $                    +irefo(isy(3))*(igya(4)                       11d19s21
     $               +irefo(isy(4))*iifmx(jtestp))))                       10d12s21
                     xintk=bc(iad)                                           10d12s21
                     if(isopt(4).ne.0)xintk=-xintk                           10d12s21
                     orig=sum
                     sum=sum-xintk                                           10d13s21
                    end if                                                   10d12s21
                   end if                                                    9d13s21
                   jmatu=imatu                                          11d19s21
                   do lk=1,4                                            11d26s21
                    if(nli(lk).gt.0)then                                11d26s21
                     iad1=ivcv+iff22k(ivcv+5+lk)                        11d26s21
                     iad2=iad1+nli(lk)                                  11d26s21
                     itmpx=ibcoff                                       11d26s21
                     itmpy=itmpx+ncsfb(jarg)*nli(lk)                    11d26s21
                     ibcoff=itmpy+ncsfb(jarg)*nli(lk)                   11d26s21
                     call enough('hcddjkbk4. 11',bc,ibc)
                     call dgemm('n','n',ncsfb(jarg),nli(lk),
     $                    ncsfk2(lk,iarg),sum,bc(jmatu),ncsfb(jarg),    11d26s21
     $                    ff22k(iad2),ncsfk2(lk,iarg),0d0,bc(itmpx),    11d26s21
     $                    ncsfb(jarg),                                  11d26s21
     d' hcddjkbk4.  8')
                     do j=0,ncsfb(jarg)-1                               11d26s21
                      do k=0,nli(lk)-1                                  11d26s21
                       jk=itmpx+j+ncsfb(jarg)*k                         11d26s21
                       kj=itmpy+k+nli(lk)*j                             11d26s21
                       bc(kj)=bc(jk)                                    11d26s21
                      end do                                            11d26s21
                     end do                                             11d26s21
                     jtmpy=itmpy                                        11d26s21
                     do lb=1,4                                          11d26s21
                      if(nlj(lb).gt.0)then                              11d26s21
                       jad1=jvcv+iff22b(jvcv+5+lb)                          8d26s21
                       jad2=jad1+nlj(lb)                                    3d19s21
                       itmpz=ibcoff                                     11d26s21
                       ibcoff=itmpz+nli(lk)*nlj(lb)                     11d26s21
                       call enough('hcddjkbk4. 12',bc,ibc)
                       call dgemm('n','n',nli(lk),nlj(lb),
     $                      ncsfb2(lb,jarg),1d0,bc(jtmpy),nli(lk),      11d26s21
     $                      ff22b(jad2),ncsfb2(lb,jarg),0d0,bc(itmpz),  11d26s21
     $                      nli(lk),                                    11d26s21
     d' hcddjkbk4.  9')
                       jtmp=itmpz                                       11d26s21
                       sz=0d0
                       do iii=0,nlj(lb)-1                                   2d24s21
                        jj=iff22b(jad1+iii)-1                              8d26s21
                        jjden=iden1en(lk,lb)+nfdatk(2,lk,isb)*jj-1           11d23s21
                        do j=0,nli(lk)-1                                    2d24s21
                         kk=iff22k(iad1+j)                                 8d26s21
                         bc(jjden+kk)=bc(jjden+kk)+bc(jtmp+j)             11d19s21
                         sz=sz+bc(jtmp+j)**2
                        end do                                             2d24s21
                        jtmp=jtmp+nli(lk)                                  11d19s21
                       end do                                              2d24s21
                       sz=sqrt(sz/dfloat(nlj(lb)*nli(lk)))
                       ibcoff=itmpz                                     11d26s21
                      end if                                            11d26s21
                      jtmpy=jtmpy+nli(lk)*ncsfb2(lb,jarg)               11d26s21
                     end do                                             11d26s21
                     ibcoff=itmpx                                       11d26s21
                    end if                                              11d26s21
                    jmatu=jmatu+ncsfb(jarg)*ncsfk2(lk,iarg)             11d26s21
                   end do                                               11d26s21
                  end do                                                11d19s21
                  ibcoff=ioutg                                          11d19s21
                  nok=0                                                         11d13s20
                  do i=1,norb                                                   11d13s20
                   itest(i,2)=0                                                    11d13s20
                   if(btest(jtc,i))itest(i,2)=2                              5d7s21
                   if(btest(jto,i))itest(i,2)=1                              5d7s21
                   ixn=min(itest(i,1),itest(i,2))
                   if(ixn.gt.0)then                                             11d13s20
                    nok=nok+1                                                   11d13s20
                    itest(nok,3)=ixn                                            11d13s20
                    itest(nok,2)=i                                              11d13s20
                   end if                                                       11d13s20
                  end do                                                        11d13s20
                  do i=1,nok                                             1d8s21
                   if(itest(i,3).eq.1.and.l2e.eq.0)then                 2d17s22
                    itestc=jtc                                              12d14s20
                    itesto=jto                                              12d14s20
                    nopenk=nopenjp                                                11d13s20
c
c     anihilate common
c
                    if(btest(itestc,itest(i,2)))then                             11d13s20
                     itestc=ibclr(itestc,itest(i,2))                             11d13s20
                     itesto=ibset(itesto,itest(i,2))                             11d13s20
                     nopenk=nopenk+1                                             11d13s20
                    else                                                         11d13s20
                     itesto=ibclr(itesto,itest(i,2))                             11d13s20
                     nopenk=nopenk-1                                             11d13s20
                    end if                                                       11d13s20
c
c     create ket
c
                    if(btest(itesto,nab4(2,1)))then                        12d14s20
                     itestc=ibset(itestc,nab4(2,1))                        12d14s20
                     itesto=ibclr(itesto,nab4(2,1))                        12d14s20
                     nopenk=nopenk-1                                             11d13s20
                    else                                                         11d13s20
                     itesto=ibset(itesto,nab4(2,1))                        12d14s20
                     nopenk=nopenk+1                                             11d13s20
                    end if                                                       11d13s20
                    call gandcr(jtc,jto,itestc,itesto,nopenjp,          11d19s21
     $              nopenk,norbxx,nnot1,nab1,icode,imap,nx1,irw1,irw2,  11d19s21
     $              iwpb1,iwpk1,bc,ibc)                                 11d14s22
                    call gandcr(itestc,itesto,itc,ito,nopenk,nopenip,   11d19s21
     $         norbxx,nnot2,nab2,icode,imap,nx2,irw1,irw2,iwpb2,iwpk2,  11d14s22
     $                   bc,ibc)                                        11d14s22
                    if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                     call spinloop(i2sb,i2smb,i2sk,i2smk,nopenjp,       11d19s21
     $                   nopenip,nopenk,ncsfb(jarg),ncsfk(iarg),itype,  11d19s21
     $                   imatx,ntypeg,nab1,iwpb1,iwpk1,nab2,iwpb2,iwpk2,11d19s21
     $                   bc(iunit),ncsfk(iarg),ieoro,bc,ibc)            11d14s22
                     do ii=0,ntypeg-1                                   11d22s21
                      ipack=ibc(itype+ii)                               11d22s21
                      do j=1,4                                               9d9s21
                       if(ipack2(j).gt.0)then                                9d9s21
                        isy(j)=ism(ipack2(j))                                9d9s21
                        igya(j)=irel(ipack2(j))-1                             9d9s21
                        imy(j)=0                                             10d13s21
                       else                                                  9d9s21
                        isy(j)=ism(-ipack2(j))                                9d9s21
                        igya(j)=irel(-ipack2(j))-1                           10d13s21
                        imy(j)=1                                             10d13s21
                       end if                                                9d9s21
                      end do                                                 9d9s21
                      ltest=imy(1)+2*(imy(2)+2*(imy(3)+2*imy(4)))            10d13s21
                      itestp=ltest+1                                         10d13s21
                      jtest=1-imy(1)+2*(1-imy(2)+2*(1-imy(3)+2*         11d19s21
     $                     (1-imy(4))))                                 11d19s21
                      jtestp=jtest+1                                         10d13s21
                      xint=0d0                                               9d10s21
                      ihit=0
                      xintj=-2d0                                             10d13s21
                      if(iifmx(itestp).ge.0)then                             10d13s21
                       iad=i4or(isy(2),isy(3),isy(4))+igya(1)           11d19s21
     $                     +irefo(isy(1))*(igya(2)+irefo(isy(2))        11d19s21
     $                     *(igya(3)+irefo(isy(3))*(igya(4)             11d19s21
     $                     +irefo(isy(4))*iifmx(itestp))))              11d19s21
                       xintj=bc(iad)                                         10d13s21
                       ihit=1                                                10d13s21
                       xint=xintj                                            10d13s21
                      else if(iifmx(jtestp).ge.0)then                        10d13s21
                       iad=i4or(isy(2),isy(3),isy(4))+igya(1)           11d19s21
     $                      +irefo(isy(1))*(igya(2)+irefo(isy(2))       11d19s21
     $                      *(igya(3)+irefo(isy(3))*                    11d19s21
     $                (igya(4)+irefo(isy(4))*iifmx(jtestp))))           10d13s21
                       xintj=bc(iad)                                         10d13s21
                       if(isopt(4).ne.0)xintj=-xintj                         10d13s21
                       ihit=1                                                10d13s21
                       xint=xintj                                            10d13s21
                      end if                                                 10d13s21
                      ltest=imy(1)+2*(imy(4)+2*(imy(3)+2*imy(2)))            10d13s21
                      itestp=ltest+1                                         10d13s21
                      jtest=1-imy(1)+2*(1-imy(4)+2*(1-imy(3)+2*         11d19s21
     $                     (1-imy(2))))                                 11d19s21
                      jtestp=jtest+1                                         10d13s21
                      xintk=-2d0                                             10d13s21
                      if(iifmx(itestp).ge.0)then                             10d13s21
                       iad=i4or(isy(4),isy(3),isy(2))+igya(1)           11d19s21
     $                     +irefo(isy(1))*(igya(4)+irefo(isy(4))        11d19s21
     $                     *(igya(3)+irefo(isy(3))*                     11d19s21
     $                (igya(2)+irefo(isy(2))*iifmx(itestp))))           10d13s21
                       xintk=bc(iad)                                         10d13s21
                       ihit=1                                                10d13s21
                       xint=xint-xintk                                       10d13s21
                      else if(iifmx(jtestp).ge.0)then                        10d13s21
                       iad=i4or(isy(4),isy(3),isy(2))+igya(1)           11d22s21
     $                      +irefo(isy(1))*(igya(4)+irefo(isy(4))*      11d22s21
     $                      (igya(3)+irefo(isy(3))*                     11d22s21
     $                (igya(2)+irefo(isy(2))*iifmx(jtestp))))           10d13s21
                       xintk=bc(iad)                                         10d13s21
                       if(isopt(4).ne.0)xintk=-xintk                         10d13s21
                       ihit=1                                                10d13s21
                       xint=xint-xintk                                       10d13s21
                      end if                                                 10d13s21
                      if(ihit.ne.0)then                                 11d22s21
                       jmat=imatx+ii*ncsfb(jarg)*ncsfk(iarg)            11d22s21
                       do lk=1,4                                        11d26s21
                        if(nli(lk).gt.0)then                            11d26s21
                         iad1=ivcv+iff22k(ivcv+5+lk)                       8d26s21
                         iad2=iad1+nli(lk)                                 3d19s21
                         itmpx=ibcoff                                   11d26s21
                         itmpy=itmpx+ncsfb(jarg)*nli(lk)                11d26s21
                         ibcoff=itmpy+ncsfb(jarg)*nli(lk)               11d26s21
                         call enough('hcddjkbk4. 13',bc,ibc)
                         call dgemm('n','n',ncsfb(jarg),nli(lk),        11d26s21
     $                        ncsfk2(lk,iarg),xint,bc(jmat),ncsfb(jarg),11d26s21
     $                        ff22k(iad2),ncsfk2(lk,iarg),0d0,bc(itmpx),11d26s21
     $                        ncsfb(jarg),                              11d26s21
     d' hcddjkbk4. 10')
                         do j=0,ncsfb(jarg)-1                           11d26s21
                          do k=0,nli(lk)-1                              11d26s21
                           jk=itmpx+j+ncsfb(jarg)*k                     11d26s21
                           kj=itmpy+k+nli(lk)*j                         11d26s21
                           bc(kj)=bc(jk)                                11d26s21
                          end do                                        11d26s21
                         end do                                         11d26s21
                         jtmpy=itmpy                                    11d26s21
                         do lb=1,4                                      11d26s21
                          if(nlj(lb).gt.0)then                          11d26s21
                           jad1=jvcv+iff22b(jvcv+5+lb)                  11d26s21
                           jad2=jad1+nlj(lb)                            11d26s21
                           itmpz=ibcoff                                 11d26s21
                           ibcoff=itmpz+nli(lk)*nlj(lb)                 11d26s21
                           call enough('hcddjkbk4. 14',bc,ibc)
                           call dgemm('n','n',nli(lk),nlj(lb),
     $                          ncsfb2(lb,jarg),1d0,bc(jtmpy),nli(lk),  11d26s21
     $                          ff22b(jad2),ncsfb2(lb,jarg),0d0,        11d26s21
     $                          bc(itmpz),nli(lk),                      11d26s21
     d' hcddjkbk4. 11')
                           jtmpm=itmpz                                  11d26s21
                           do iii=0,nlj(lb)-1                           11d26s21
                            jj=iff22b(jad1+iii)-1                           8d26s21
                            jjden=iden1en(lk,lb)+nfdatk(2,lk,isb)*jj-1  11d26s21
                            do j=0,nli(lk)-1                            11d26s21
                             kk=iff22k(iad1+j)                              8d26s21
                             bc(jjden+kk)=bc(jjden+kk)+bc(jtmpm+j)             2d24s21
                            end do                                             2d24s21
                            jtmpm=jtmpm+nli(lk)                         11d26s21
                           end do                                              2d24s21
                           ibcoff=itmpz                                 11d26s21
                          end if                                        11d26s21
                          jtmpy=jtmpy+nli(lk)*ncsfb2(lb,jarg)           11d26s21
                         end do                                         11d26s21
                         ibcoff=itmpx                                   11d26s21
                        end if                                          11d26s21
                        jmat=jmat+ncsfb(jarg)*ncsfk2(lk,iarg)           11d26s21
                       end do                                           11d26s21
                      end if                                                  12d14s20
                     end do                                             11d22s21
                     ibcoff=itype                                       11d22s21
                    end if                                                   12d14s20
                   end if                                               11d22s21
                  end do                                                 1d8s21
                 else if(l2e.eq.0)then                                  2d23s22
                  nnot=0                                                  2d6s23
                  if(ndifs.eq.4.and.ndifb.eq.4)then                                           10d14s
                   nnot=4                                                          10d14s22
                   ioxx(1)=1                                                     10d17s22
                   ioxx(2)=1                                                     10d17s22
                   do i=1,norbxx                                                     10d17s22
                    if(btest(gandcb,i))then                                         10d14s22
                     if((btest(itc,i).and.btest(jto,i)).or.                         10d17s22
     $                (btest(ito,i).and..not.btest(jtc,i)))then                   10d14s22
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
     $                 ((btest(jtc,i).and..not.btest(ito,i)).or.                 10d14s22
     $         (btest(itc,i).and..not.btest(jto,i))))then                     10d14s22
                      if(btest(itc,i))iswap=1                                        10d17s22
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
                    if(btest(itc,nab4(2,2)).and.                        2d6s23
     $                   .not.btest(itc,nab4(2,1)))nbt=1                2d6s23
                   else                                                             10d17s22
                    nbt=0                                                           10d17s22
                    if(btest(jtc,nab4(1,2)).and.                        2d6s23
     $                   .not.btest(jtc,nab4(1,1)))nbt=1                2d6s23
                   end if                                                           10d17s22
                   if(nbt.ne.0)then                                                 10d17s22
                    nab4(1,1)=nab4(1,2)                                             10d17s22
                    nab4(2,1)=nab4(2,2)                                             10d17s22
                   end if                                                           10d17s22
                  else if(ndifs.eq.0.and.ndifd.eq.2)then                            10d14s22
                   nnot=3
                   do i=1,norbxx
                    if(btest(gandcb,i))then                                         10d14s22
                     if(btest(jtc,i))then
                      nab4(1,1)=i
                      nab4(1,2)=i
                     else                                                           10d14s22
                      nab4(2,1)=i
                      nab4(2,2)=i
                     end if
                    end if                                                          10d14s22
                   end do
                  end if                                                            10d14s22
                  if(nnot.ne.0)then                                     2d7s23
                   iu1=1
                   iu2=1
                   itestc=jtc                                            11d23s21
                   itesto=jto                                            11d23s21
                   if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                    itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                    itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                    nopenk=nopenjp+1                                     11d23s21
                   else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                    itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                    nopenk=nopenjp-1                                     11d23s21
                   else                                                          11d13s20
                    write(6,*)('bitb not set for nab4(1,1) ='),           3d1s21
     $                   nab4(1,iu1)                                    3d1s21
                    stop 'nab4b(1,1)'                                           11d27s20
                   end if                                                        11d13s20
                   if(btest(itesto,nab4(2,iu2)))then                         12d8s20
                    itestc=ibset(itestc,nab4(2,iu2))                         12d8s20
                    itesto=ibclr(itesto,nab4(2,iu2))                             11d27s20
                    nopenk=nopenk-1                                              11d13s20
                   else if(btest(itesto,nab4(2,iu2)))then                    12d8s20
                    write(6,*)('already double in nab4(2,1) = '),          12d18s20
     $                 nab4(2,iu2)                                       12d18s20
                    stop 'nab4c(2,1)'                                           11d27s20
                   else                                                          11d13s20
                    itesto=ibset(itesto,nab4(2,iu2))                         12d8s20
                    nopenk=nopenk+1                                              11d13s20
                   end if                                                        11d13s20
                   call gandcr(jtc,jto,itestc,itesto,nopenjp,nopenk,     11d22s21
     $         norbxx,nnot1,nab1,icode,imap,nx1,irw1,irw2,iwpb1,iwpk1,  1d18s23
     $                 bc,ibc)                                          11d14s22
                   call gandcr(itestc,itesto,itc,ito,nopenk,nopenip,      2d25s21
     $         norbxx,nnot2,nab2,icode,imap,nx1,irw1,irw2,iwpb2,iwpk2,  1d18s23
     $                 bc,ibc)                                          11d14s22
                   if(nnot1.eq.2.and.nnot2.eq.2)then                       12d19s20
                    call spinloop(i2sb,i2smb,i2sk,i2smk,nopenjp,nopenip, 11d22s21
     $               nopenk,ncsfb(jarg),ncsfk(iarg),itype,imatx,ntypeq, 11d22s21
     $                 nab1,iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(iunit),     11d22s21
     $                 ncsfk(iarg),ieoro,bc,ibc)                        11d14s22
                    do ii=0,ntypeq-1                                      11d22s21
                     ipack=ibc(itype+ii)                                 11d22s21
                     ltest=0                                                10d12s21
                     jtest=0                                                10d12s21
                     mtest=1                                                10d12s21
                     do j=1,4                                               9d9s21
                      if(ipack2(j).gt.0)then                                9d9s21
                       isy(j)=ism(ipack2(j))                                9d9s21
                       igya(j)=irel(ipack2(j))-1                             9d9s21
                       jtest=jtest+mtest                                    10d12s21
                       imy(j)=0                                             10d12s21
                      else                                                  9d9s21
                       isy(j)=ism(-ipack2(j))                                9d9s21
                       igya(j)=irel(-ipack2(j))-1                           10d12s21
                       ltest=ltest+mtest                                    10d12s21
                       imy(j)=1                                             10d12s21
                      end if                                                9d9s21
                      mtest=mtest*2                                         10d12s21
                     end do                                                 9d9s21
                     itestp=ltest+1                                         10d12s21
                     jtestp=jtest+1                                         10d12s21
                     ihit=0
                     xintj=-2d0                                             10d12s21
                     xint=0d0                                               11d9s21
                     jchoice=-1
                     kchoice=-1
                     if(iifmx(itestp).ge.0)then                             10d12s21
                      iad=i4or(isy(2),isy(3),isy(4))+igya(1)             11d22s21
     $                   +irefo(isy(1))*                                11d22s21
     $          (igya(2)+irefo(isy(2))*(igya(3)+irefo(isy(3))*(         10d12s21
     $                igya(4)+irefo(isy(4))*iifmx(itestp))))            10d12s21
                      xintj=bc(iad)                                          10d12s21
                      ihit=1                                                10d12s21
                      xint=xintj                                            10d12s21
                      jchoice=iifmx(itestp)
                     else if(iifmx(jtestp).ge.0)then                        10d12s21
                      iad=i4or(isy(2),isy(3),isy(4))+igya(1)             11d22s21
     $                    +irefo(isy(1))*                               11d22s21
     $          (igya(2)+irefo(isy(2))*(igya(3)+irefo(isy(3))*(         10d12s21
     $                igya(4)+irefo(isy(4))*iifmx(jtestp))))            10d12s21
                      xintj=bc(iad)                                          10d12s21
                      jchoice=100+iifmx(jtestp)
                      if(isopt(4).ne.0)then                                 10d12s21
                       xintj=-xintj                                           10d12s21
                      end if
                      ihit=1                                                10d12s21
                      xint=xintj                                            10d12s21
                     end if                                                 10d12s21
                     xintk=-2d0                                             10d12s21
                     if(nnot.eq.4)then                                      10d12s21
                      ltest=imy(1)+2*(imy(4)+2*(imy(3)+2*imy(2)))            10d12s21
                      itestp=ltest+1                                         10d12s21
                      jtest=1-imy(1)+2*(1-imy(4)+2*(1-imy(3)+2           11d22s21
     $                    *(1-imy(2))))                                 11d22s21
                      jtestp=jtest+1                                         10d12s21
                      if(iifmx(itestp).ge.0)then                             10d12s21
                       iad=i4or(isy(4),isy(3),isy(2))+igya(1)            11d22s21
     $                    +irefo(isy(1))*                               11d22s21
     $          (igya(4)+irefo(isy(4))*(igya(3)+irefo(isy(3))*(         10d12s21
     $                igya(2)+irefo(isy(2))*iifmx(itestp))))            10d12s21
                       xintk=bc(iad)                                          10d12s21
                       xint=xint-xintk                                       10d12s21
                       kchoice=iifmx(itestp)
                       ihit=1                                               11d9s21
                      else if(iifmx(jtestp).ge.0)then                        10d12s21
                       iad=i4or(isy(4),isy(3),isy(2))+igya(1)            11d22s21
     $                     +irefo(isy(1))*                              11d22s21
     $          (igya(4)+irefo(isy(4))*(igya(3)+irefo(isy(3))*(         10d12s21
     $                igya(2)+irefo(isy(2))*iifmx(jtestp))))            10d12s21
                       xintk=bc(iad)                                          10d12s21
                       if(isopt(4).ne.0)then                                 10d12s21
                        xintk=-xintk                                           10d12s21
                       end if
                       xint=xint-xintk                                       10d12s21
                       kchoice=100+iifmx(jtestp)
                       ihit=1                                               11d9s21
                      end if                                                 10d12s21
                     end if                                                 10d12s21
                     if(ihit.ne.0)then                                      9d13s21
                      xint0=xint
                      if(nab4(1,1).eq.nab4(1,2).and.nab4(2,1).eq.        11d22s21
     $                    nab4(2,2))                                    11d22s21
     $                 xint=xint*0.5d0                                  9d9s21
                      jmat=imatx+ii*ncsfb(jarg)*ncsfk(iarg)              11d22s21
                      do lk=1,4                                          11d23s21
                       if(nli(lk).gt.0)then                              11d23s21
                        iad1=ivcv+iff22k(ivcv+5+lk)                       8d26s21
                        iad2=iad1+nli(lk)                                  3d19s21
                        itmpx=ibcoff                                     11d22s21
                        itmpy=itmpx+ncsfb(jarg)*nli(lk)                  11d23s21
                        ibcoff=itmpy+ncsfb(jarg)*nli(lk)                 11d23s21
                        call enough('hcddjkbk4. 15',bc,ibc)
                        call dgemm('n','n',ncsfb(jarg),nli(lk),          11d23s21
     $                      ncsfk2(lk,iarg),xint,bc(jmat),ncsfb(jarg),  11d23s21
     $                      ff22k(iad2),ncsfk2(lk,iarg),0d0,bc(itmpx),  11d23s21
     $                      ncsfb(jarg),                                11d23s21
     d' hcddjkbk4. 12')
                        do j=0,ncsfb(jarg)-1                             11d23s21
                         do k=0,nli(lk)-1                                 11d22s21
                          kj=itmpy+k+nli(lk)*j                            11d22s21
                          jk=itmpx+j+ncsfb(jarg)*k                       11d23s21
                          bc(kj)=bc(jk)                                  11d22s21
                         end do                                          11d22s21
                        end do                                           11d22s21
                        jtmpy=itmpy                                      11d23s21
                        do lb=1,4                                        11d23s21
                         if(nlj(lb).gt.0)then                            11d23s21
                          jad1=jvcv+iff22b(jvcv+5+lb)                    11d23s21
                          jad2=jad1+nlj(lb)                              11d23s21
                          itmpz=ibcoff                                   11d26s21
                          ibcoff=itmpz+nli(lk)*nlj(lb)                   11d26s21
                          call enough('hcddjkbk4. 16',bc,ibc)
                          call dgemm('n','n',nli(lk),nlj(lb),            11d23s21
     $                        ncsfb2(lb,jarg),1d0,bc(jtmpy),nli(lk),    11d23s21
     $                        ff22b(jad2),ncsfb2(lb,jarg),0d0,bc(itmpz),11d26s21
     $                        nli(lk),                                  11d23s21
     d' hcddjkbk4. 13')
                          jtmpm=itmpz                                    11d26s21
                          do iii=0,nlj(lb)-1                                   2d24s21
                           jj=iff22b(jad1+iii)-1                           8d26s21
                           jjden=iden1en(lk,lb)+nfdatk(2,lk,isb)*jj-1    11d23s21
                           do j=0,nli(lk)-1                              11d23s21
                            kk=iff22k(iad1+j)                              8d26s21
                            bc(jjden+kk)=bc(jjden+kk)+bc(jtmpm+j)             2d24s21
                           end do                                             2d24s21
                           jtmpm=jtmpm+nli(lk)                                 2d24s21
                          end do                                              2d24s21
                          ibcoff=itmpz                                   11d26s21
                         end if                                          11d23s21
                         jtmpy=jtmpy+nli(lk)*ncsfb2(lb,jarg)             11d23s21
                        end do                                           11d23s21
                        ibcoff=itmpx                                     11d23s21
                       end if                                            11d23s21
                       jmat=jmat+ncsfb(jarg)*ncsfk2(lk,iarg)             11d23s21
                      end do                                             11d23s21
                     end if                                                  12d18s20
                    end do                                               11d22s21
                   end if                                               2d7s23
                   ibcoff=itype                                         11d22s21
                  end if                                                4d20s21
                 end if                                                  1d8s21
                 ibcoff=iunit                                           11d22s21
                end if                                                   1d8s21
               end if                                                    1d8s21
               jto=ibclr(jto,norbx)                                      1d8s21
               jto=ibset(jto,norbxxx)                                    1d8s21
c                                                                       1d8s21
               gandcc=ieor(itc,jtc)                                      10d13s22
               gandco=ieor(ito,jto)                                      10d13s22
               gandcb=ior(gandcc,gandco)                                 10d20s22
               ndifb=popcnt(gandcb)                                      10d20s22
               if(ndifb.le.4)then                                        10d20s22
                ndifs=popcnt(gandco)                                      10d13s22
                ndifd=popcnt(gandcc)                                      10d13s22
                iunit=ibcoff                                            11d22s21
                ibcoff=iunit+ncsfk(iarg)*ncsfk(iarg)                    11d22s21
                call enough('hcddjkbk4. 17',bc,ibc)
                do iz=iunit,ibcoff-1                                    11d22s21
                 bc(iz)=0d0                                             11d22s21
                end do                                                  11d22s21
                do iz=0,ncsfk(iarg)-1                                   11d22s21
                 iad=iunit+iz*(ncsfk(iarg)+1)                           11d22s21
                 bc(iad)=1d0                                            11d22s21
                end do                                                  11d22s21
                if(ndifs.eq.2.and.ndifb.eq.2)then                       2d7s23
                 do i=1,norbxxx                                         2d7s23
                  if(btest(gandco,i))then                               2d7s23
                   if((btest(jtc,i).and.btest(ito,i)).or.               2d7s23
     $             (btest(jto,i).and..not.btest(itc,i)))then             2d7s23
                    nab4(1,1)=i                                         2d7s23
                   else                                                 2d7s23
                    nab4(2,1)=i                                         2d7s23
                   end if                                               2d7s23
                  end if                                                2d7s23
                 end do                                                 2d7s23
                 if(isb.ne.jsb)then                                      1d9s21
                  write(6,*)('hey, what''s up with isb ne jsb ? ')       1d9s21
                  stop 'hcddjk:g'                                                   1d9s21
                 end if                                                  1d9s21
                 nx1=0                                                  11d22s21
                 do l=1,4                                               2d23s21
                  nsz=max(nlj(l),ncsfb2(l,jarg))*nli(l)                 11d23s21
                  nx1=max(nx1,nsz)                                      11d22s21
                 end do                                                 2d23s21
                 itmp1=ibcoff                                           2d23s21
                 itmp2=itmp1+nx1                                        11d22s21
                 ibcoff=itmp2+nx1                                       11d22s21
                 call enough('hcddjkbk4. 18',bc,ibc)
                 nok=0                                                         11d13s20
                 do i=1,norb                                                   11d13s20
                  itest(i,2)=0                                                    11d13s20
                  if(btest(jtc,i))itest(i,2)=2                              5d7s21
                  if(btest(jto,i))itest(i,2)=1                              5d7s21
                  ixn=min(itest(i,1),itest(i,2))
                  if(ixn.gt.0)then                                             11d13s20
                   nok=nok+1                                                   11d13s20
                   itest(nok,3)=ixn                                            11d13s20
                   itest(nok,2)=i                                              11d13s20
                  end if                                                       11d13s20
                 end do                                                        11d13s20
                 call gandcr(jtc,jto,itc,ito,nopenjp,nopenip,norbxxx,   11d22s21
     $                nnot1,nab1,icode,imap,nx1,irw1,irw2,iwpb1,iwpk1,  11d14s22
     $                bc,ibc)                                           11d14s22
                 call gencup(i2sb,i2smb,i2sk,i2smk,nopenjp,nopenip,
     $               nab1,iwpb1,iwpk1,ioutg,imatg,ntypeg,mcsf,bc(iunit),11d22s21
     $               ncsfk(iarg),ncsfk(iarg),bc,ibc)                    11d14s22
                 do ii=0,ntypeg-1                                       11d22s21
                  ipackc=ibc(ioutg+ii)                                  11d22s21
                  phs=1d0                                               11d22s21
                  if(ipackc1(1).gt.0)then                               11d22s21
                   imy(3)=0                                             11d22s21
                   nhvabp=nhvabp+1                                      2d7s22
                  else                                                  11d22s21
                   imy(3)=1                                             11d22s21
                   nhvabm=nhvabm+1                                      2d7s22
                   if(isopt(4).ne.0)phs=-phs                            11d22s21
                  end if                                                11d22s21
                  if(ipackc1(2).gt.0)then                               11d22s21
                   imy(4)=0                                             11d22s21
                  else                                                  11d22s21
                   imy(4)=1                                             11d22s21
                  end if                                                11d22s21
                  jmat=imatg+ncsfb(jarg)*ncsfk(iarg)*ii                 11d22s21
                  jmat0=jmat
                  do lk=1,4                                             11d26s21
                   if(nli(lk).gt.0)then                                 11d26s21
                    iad1=ivcv+iff22k(ivcv+5+lk)                         11d26s21
                    iad2=iad1+nli(lk)                                   11d26s21
                    itmpx=ibcoff                                        11d26s21
                    itmpy=itmpx+nli(lk)*ncsfb(jarg)                     11d26s21
                    ibcoff=itmpy+nli(lk)*ncsfb(jarg)                    11d26s21
                    call enough('hcddjkbk4. 19',bc,ibc)
                    call dgemm('n','n',ncsfb(jarg),nli(lk),             11d26s21
     $                   ncsfk2(lk,iarg),1d0,bc(jmat),ncsfb(jarg),      11d26s21
     $                   ff22k(iad2),ncsfk2(lk,iarg),0d0,bc(itmpx),     11d26s21
     $                   ncsfb(jarg),                                   11d26s21
     d' hcddjkbk4. 14')
                    do j=0,ncsfb(jarg)-1                                11d26s21
                     do k=0,nli(lk)-1                                   11d26s21
                      kj=itmpy+k+nli(lk)*j                              11d26s21
                      jk=itmpx+j+ncsfb(jarg)*k                          11d26s21
                      bc(kj)=bc(jk)                                     11d26s21
                     end do                                             11d26s21
                    end do                                              11d26s21
                    jtmpy=itmpy                                         11d26s21
                    do lb=1,4                                           11d26s21
                     if(nlj(lb).gt.0)then                               11d26s21
                      jad1=jvcv+iff22b(jvcv+5+lb)                       11d26s21
                      jad2=jad1+nlj(lb)                                 11d26s21
                      itmpz=ibcoff                                      11d26s21
                      ibcoff=itmpz+nli(lk)*nlj(lb)                      11d26s21
                      call enough('hcddjkbk4. 20',bc,ibc)
                      call dgemm('n','n',nli(lk),nlj(lb),
     $                     ncsfb2(lb,jarg),1d0,bc(jtmpy),nli(lk),       11d26s21
     $                     ff22b(jad2),ncsfb2(lb,jarg),0d0,bc(itmpz),   11d26s21
     $                     nli(lk),                                     11d26s21
     d' hcddjkbk4. 15')
                      ktmpn=itmpz                                       11d26s21
                      do iii=0,nlj(lb)-1                                11d26s21
                       jj=iff22b(jad1+iii)-1                               8d26s21
                       jjden=idenhvvn(lk,lb)+nfdatk(2,lk,isb)*jj-1      11d26s21
                       do j=0,nli(lk)-1                                 11d26s21
                        kk=iff22k(iad1+j)                                  8d26s21
                        bc(jjden+kk)=bc(jjden+kk)+bc(ktmpn+j)*phs         11d23s21
                       end do                                              2d24s21
                       ktmpn=ktmpn+nli(lk)                              11d26s21
                      end do                                               2d24s21
c     imy(3) is v and imy(4) is v' ...
c     we need order to be v,v'
                      do i=1,norb                                         11d23s21
                       if(btest(itc,i).and.btest(jtc,i))then              11d23s21
                        js=ism(i)                                               9d10s21
                        jg=irel(i)-1                                            9d10s21
                        do jpass=0,1                                       11d22s21
c
c     (ii|vv')=(vv'|ii)
c
                         ltest=imy(3)+2*(imy(4)+2*(jpass+2*jpass))        11d23s21
                         itestp=ltest+1                                         10d12s21
                         jtest=1-imy(3)+2*(1-imy(4)+2*(1-jpass            11d23s21
     $                      +2*(1-jpass)))                              11d23s21
                         jtestp=jtest+1                                    11d22s21
                         iuse=-1                                           11d22s21
                         phsx=1d0                                          11d22s21
                         if(iifmx(itestp).ge.0)then                        11d22s21
                          iuse=iifmx(itestp)                               11d22s21
                         else if(iifmx(jtestp).ge.0)then                   11d22s21
                          iuse=iifmx(jtestp)                               11d22s21
                          if(isopt(4).ne.0)phsx=-phsx                      11d22s21
                         end if                                            11d22s21
                         if(iuse.ge.0)then                                 11d22s21
                          icol=jg+irefo(js)*(jg+irefo(js)*iuse)            11d22s21
                          jden=idenjn(lk,lb,js)+nfdatk(2,lk,isb)        11d26s21
     $                       *nfdatb(2,lb,jsb)*icol                      11d23s21
                          ktmpn=itmpz                                   11d26s21
                          do iii=0,nlj(lb)-1                                    2d24s21
                           jj=iff22b(jad1+iii)-1                               8d26s21
                           jjden=jden+jj*nfdatk(2,lk,isb)-1                11d23s21
                           do j=0,nli(lk)-1                                     2d24s21
                            kk=iff22k(iad1+j)                                  8d26s21
                            bc(jjden+kk)=bc(jjden+kk)+bc(ktmpn+j)*phsx     11d22s21
                           end do                                              2d24s21
                           ktmpn=ktmpn+nli(lk)                                  2d24s21
                          end do                                               2d24s21
                         end if                                            11d22s21
c
c     (iv'|vi)=(vi|iv')=+/-(iv|v'i)
c
                         ltest=jpass+2*(imy(3)+2*(imy(4)+2*jpass))         11d22s21
                         itestp=ltest+1                                         10d12s21
                         jtest=1-jpass+2*(1-imy(3)+2*(1-imy(4)+2           11d22s21
     $                        *(1-jpass)))                                 11d22s21
                         jtestp=jtest+1                                    11d22s21
                         phsx=-1d0                                         11d22s21
                         if(isopt(2).ne.0)phsx=-phsx                      11d23s21
                         iuse=-1                                           11d22s21
                         if(iifmx(itestp).ge.0)then                        11d22s21
                          iuse=iifmx(itestp)                               11d22s21
                         else if(iifmx(jtestp).ge.0)then                   11d22s21
                          iuse=iifmx(jtestp)                               11d22s21
                          if(isopt(4).ne.0)phsx=-phsx                      11d22s21
                         end if                                            11d22s21
                         if(iuse.ge.0)then                                 11d22s21
                          icol=jg+irefo(js)*(jg+irefo(js)*iuse)            11d22s21
                          jden=idenkn(lk,lb,js)+nfdatk(2,lk,isb)        11d26s21
     $                      *nfdatb(2,lb,jsb)*icol                       11d22s21
                          ktmpn=itmpz                                   11d26s21
                          do iii=0,nlj(lb)-1                                    2d24s21
                           jj=iff22b(jad1+iii)-1                               8d26s21
                           jjden=jden+jj*nfdatk(2,lk,isb)-1                11d23s21
                           do j=0,nli(lk)-1                                     2d24s21
                            kk=iff22k(iad1+j)                                  8d26s21
                            bc(jjden+kk)=bc(jjden+kk)+bc(ktmpn+j)*phsx     11d22s21
                           end do                                              2d24s21
                           ktmpn=ktmpn+nli(lk)                                  2d24s21
                          end do                                               2d24s21
                         end if                                            11d22s21
                        end do                                             11d22s21
                       end if                                             11d23s21
                      end do                                              11d22s21
                      ibcoff=itmpz                                      11d26s21
                     end if                                             11d26s21
                     jtmpy=jtmpy+nli(lk)*ncsfb2(lb,jarg)                11d26s21
                    end do                                              11d26s21
                    ibcoff=itmpx                                        11d26s21
                   end if                                               11d26s21
                   jmat=jmat+ncsfb(jarg)*ncsfk2(lk,iarg)                11d26s21
                  end do                                                11d26s21
                 end do                                                 11d22s21
                 ibcoff=ioutg                                           11d22s21
                 do i=1,nok                                              1d8s21
                  if(itest(i,3).eq.1.and.l2e.eq.0)then                  2d17s22
                   itestc=jtc                                           11d23s21
                   itesto=jto                                           11d23s21
                   nopenk=nopenjp                                       11d23s21
c
c     j is bra, i is ket.
c
c         J   I    K     thus  JK  KI
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
                    nopenk=nopenk+1                                             11d13s20
                   else                                                         11d13s20
                    itesto=ibclr(itesto,itest(i,2))                             11d13s20
                    nopenk=nopenk-1                                             11d13s20
                   end if                                                       11d13s20
c
c     create ket
c
                   if(btest(itesto,nab4(2,1)))then                        12d14s20
                    itestc=ibset(itestc,nab4(2,1))                        12d14s20
                    itesto=ibclr(itesto,nab4(2,1))                        12d14s20
                    nopenk=nopenk-1                                             11d13s20
                   else                                                         11d13s20
                    itesto=ibset(itesto,nab4(2,1))                        12d14s20
                    nopenk=nopenk+1                                             11d13s20
                   end if                                                       11d13s20
                   call gandcr(jtc,jto,itestc,itesto,nopenjp,nopenk,    11d22s21
     $            norbxxx,nnot1,nab1,icode,imap,nx1,irw1,irw2,iwpb1,    11d22s21
     $                  iwpk1,bc,ibc)                                   11d14s22
                   call gandcr(itestc,itesto,itc,ito,nopenk,nopenip,    11d22s21
     $            norbxxx,nnot2,nab2,icode,imap,nx1,irw1,irw2,iwpb2,    11d22s21
     $               iwpk2,bc,ibc)                                      11d14s22
                   if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                    call spinloop(i2sb,i2smb,i2sk,i2smk,nopenjp,nopenip,11d22s21
     $               nopenk,ncsfb(jarg),ncsfk(iarg),itype,imatx,ntypeq, 11d22s21
     $                  nab1,iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(iunit),    11d22s21
     $                  ncsfk(iarg),ieoro,bc,ibc)                       11d14s22
                    do ii=0,ntypeq-1                                    11d22s21
                     ipack=ibc(itype+ii)                                11d22s21
                     do j=1,4                                               9d9s21
                      if(ipack2(j).gt.0)then                                9d9s21
                       ipack2a(j)=ipack2(j)                             11d22s21
                       imy(j)=0                                             10d13s21
                      else                                                  9d9s21
                       ipack2a(j)=-ipack2(j)                            11d22s21
                       imy(j)=1                                             10d13s21
                      end if                                                9d9s21
                      if(ipack2a(j).le.norb)then                        11d22s21
                       isy(j)=ism(ipack2a(j))                                9d9s21
                       igya(j)=irel(ipack2a(j))-1                             9d9s21
                      end if                                            11d22s21
                     end do                                                 9d9s21
                     if(ipack2a(2).eq.norbx.and.ipack2a(3).eq.norbxxx)  11d23s21
     $                    then                                          11d23s21
c
c     (1v'|v4)=(v4|1v')=+/-(4v|v'1)
c
                      phsk=1d0                                          11d22s21
                      if(isopt(2).ne.0)phsk=-phsk                       11d23s21
                      kuse=-1                                           11d22s21
                      ltest=imy(4)+2*(imy(3)+2*(imy(2)+2*imy(1)))       11d23s21
                      jtest=1-imy(4)+2*(1-imy(3)+2*(1-imy(2)            11d22s21
     $                     +2*(1-imy(1))))                              11d22s21
                      itestp=ltest+1                                    11d22s21
                      jtestp=jtest+1                                    11d22s21
                      if(iifmx(itestp).ge.0)then                        11d22s21
                       kuse=iifmx(itestp)                               11d22s21
                      else if(iifmx(jtestp).ge.0)then                   11d22s21
                       kuse=iifmx(jtestp)                               11d22s21
                       if(isopt(4).ne.0)phsk=-phsk                      11d22s21
                      end if                                            11d22s21
c
c     (14|vv')=(vv'|14)
c
                      phsj=-1d0                                         11d22s21
                      juse=-1                                           11d22s21
                      ltest=imy(3)+2*(imy(2)+2*(imy(1)+2*imy(4)))       11d23s21
                      jtest=1-imy(3)+2*(1-imy(2)+2*(1-imy(1)            11d23s21
     $                     +2*(1-imy(4))))                              11d23s21
                      itestp=ltest+1                                    11d22s21
                      jtestp=jtest+1                                    11d22s21
                      if(iifmx(itestp).ge.0)then                        11d22s21
                       juse=iifmx(itestp)                               11d22s21
                      else if(iifmx(jtestp).ge.0)then                   11d22s21
                       juse=iifmx(jtestp)                               11d22s21
                       if(isopt(4).ne.0)phsj=-phsj                      11d22s21
                      end if                                            11d22s21
                     else                                               11d22s21
                      write(6,*)('don''t know how to handle typec '),
     $                     ipack2
                      call dws_synca
                      call dws_finalize
                      stop
                     end if                                             11d22s21
                     jmat=imatx+ncsfb(jarg)*ncsfk(iarg)*ii              11d22s21
                     if(min(juse,kuse).ge.0)then                        11d22s21
                      jcol=igya(1)+irefo(isy(1))*(igya(4)                11d22s21
     $                   +irefo(isy(4))*juse)                           11d22s21
                      kcol=igya(4)+irefo(isy(4))*(igya(1)                11d22s21
     $                   +irefo(isy(1))*kuse)                           11d22s21
                      do lk=1,4                                          11d22s21
                       if(nli(lk).gt.0)then                              11d22s21
                        iad1=ivcv+iff22k(ivcv+5+lk)                      11d22s21
                        iad2=iad1+nli(lk)                                11d22s21
                        itmpx=ibcoff                                     11d22s21
                        itmpy=itmpx+nli(lk)*ncsfb(jarg)                  11d22s21
                        ibcoff=itmpy+nli(lk)*ncsfb(jarg)                11d26s21
                        call dgemm('n','n',ncsfb(jarg),nli(lk),
     $                      ncsfk2(lk,iarg),1d0,bc(jmat),ncsfb(jarg),      11d22s21
     $                      ff22k(iad2),ncsfk2(lk,iarg),0d0,bc(itmpx),  11d23s21
     $                      ncsfb(jarg),                                11d22s21
     d' hcddjkbk4. 16')
                        do j=0,ncsfb(jarg)-1                             11d22s21
                         do k=0,nli(lk)-1                                11d22s21
                          kj=itmpy+k+nli(lk)*j                           11d22s21
                          jk=itmpx+j+ncsfb(jarg)*k                       11d22s21
                          bc(kj)=bc(jk)                                  11d22s21
                         end do                                          11d22s21
                        end do                                           11d22s21
                        jtmpy=itmpy                                      11d22s21
                        do lb=1,4                                        11d22s21
                         if(nlj(lb).gt.0)then                            11d22s21
                          jad1=jvcv+iff22b(jvcv+5+lb)                    11d22s21
                          jad2=jad1+nlj(lb)                              11d22s21
                          itmpz=ibcoff                                  11d26s21
                          ibcoff=itmpz+nli(lk)*nlj(lb)                  11d26s21
                          call enough('hcddjkbk4. 21',bc,ibc)
                          call dgemm('n','n',nli(lk),nlj(lb),
     $                        ncsfb2(lb,jarg),1d0,bc(jtmpy),nli(lk),    11d22s21
     $                       ff22b(jad2),ncsfb2(lb,jarg),               11d23s21
     $                      0d0,bc(itmpz),nli(lk),                      11d26s21
     d' hcddjkbk4. 17')
                          if(juse.ge.0)then                             11d23s21
                           ktmpn=itmpz                                  11d26s21
                           jden=idenjn(lk,lb,isy(4))                    11d23s21
     $                        +nfdatk(2,lk,isb)*nfdatb(2,lb,jsb)*jcol   11d23s21
                           do iii=0,nlj(lb)-1                                    2d24s21
                            jj=iff22b(jad1+iii)-1                               8d26s21
                            jjden=jden+nfdatk(2,lk,isb)*jj-1                 11d22s21
                            do j=0,nli(lk)-1                                     2d24s21
                             kk=iff22k(iad1+j)                                  8d26s21
                             bc(jjden+kk)=bc(jjden+kk)+bc(ktmpn+j)*phsj     11d22s21
                            end do                                              2d24s21
                            ktmpn=ktmpn+nli(lk)                                  2d24s21
                           end do                                               2d24s21
                          end if                                         11d22s21
                          if(kuse.ge.0)then                             11d26s21
                           ktmpn=itmpz                                  11d26s21
                           jden=idenkn(lk,lb,isy(4))                            11d22s21
     $                         +nfdatk(2,lk,isb)*nfdatb(2,lb,jsb)*kcol      11d22s21
                           do iii=0,nlj(lb)-1                                    2d24s21
                            jj=iff22b(jad1+iii)-1                               8d26s21
                            jjden=jden+nfdatk(2,lk,isb)*jj-1            11d26s21
                            do j=0,nli(lk)-1                                     2d24s21
                             kk=iff22k(iad1+j)                                  8d26s21
                             bc(jjden+kk)=bc(jjden+kk)+bc(ktmpn+j)*phsk   11d22s21
                            end do                                              2d24s21
                            ktmpn=ktmpn+nli(lk)                                  2d24s21
                           end do                                               2d24s21
                          end if
                          ibcoff=itmpz                                  11d26s21
                         end if                                          11d22s21
                         jtmpy=jtmpy+nli(lk)*ncsfb2(lb,jarg)             11d22s21
                        end do                                           11d22s21
                        ibcoff=itmpx                                     11d22s21
                       end if                                            11d22s21
                       jmat=jmat+ncsfb(jarg)*ncsfk2(lk,iarg)             11d22s21
                      end do                                             11d22s21
                     end if                                             11d22s21
                    end do                                              11d22s21
                    ibcoff=itype                                        11d22s21
                   end if                                                  12d14s20
                  end if                                                   12d14s20
                 end do                                                  12d18s20
c
c     coupling coefficients for 4 v ints ...
c
                 if(l2e.eq.0.and.n4vso.ne.0)then                        2d8s23
                  jto=ibclr(jto,norbxx)                                  11d30s21
                  jto=ibset(jto,norbxxxx)                                11d30s21
                  gandcc2=ieor(itc,jtc)                                 2d8s23
                  gandco2=ieor(ito,jto)                                 2d8s23
                  gandcb2=ior(gandcc2,gandco2)                          2d8s23
                  ndifb2=popcnt(gandcb2)                                2d8s23
                  ndifs2=popcnt(gandco2)                                2d8s23
                  ndifd2=popcnt(gandcc2)                                2d8s23
                  nnot=0                                                2d7s23
                  if(ndifs2.eq.4.and.ndifb2.eq.4)then                   2d8s23
                   nnot=4                                                          10d14s22
                   ioxx(1)=1                                                     10d17s22
                   ioxx(2)=1                                                     10d17s22
                   do i=1,norbxxxx                                      2d7s23
                    if(btest(gandcb2,i))then                            2d8s23
                     if((btest(itc,i).and.btest(jto,i)).or.                         10d17s22
     $                (btest(ito,i).and..not.btest(jtc,i)))then                   10d14s22
                      nab4(2,ioxx(2))=i                                       10d17s22
                      ioxx(2)=ioxx(2)+1                                       10d17s22
                     else                                                           10d14s22
                      nab4(1,ioxx(1))=i                                       10d17s22
                      ioxx(1)=ioxx(1)+1                                       10d17s22
                     end if                                                         10d14s22
                    end if
                   end do
                  end if                                                2d7s23
                  if(nnot.ne.4)then                                      11d30s21
                   write(6,*)('expecting 4, but got '),nnot,('???')
                   call dws_synca                                        11d30s21
                   call dws_finalize                                     11d30s21
                   stop                                                  11d30s21
                  end if                                                 11d30s21
                  nden4=1                                                11d30s21
                  iu1=1
                  iu2=1
                  itestc=jtc                                             11d23s21
                  itesto=jto                                             11d23s21
                  if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                   itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                   itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopenjp+1                                      11d23s21
                  else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                   itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopenjp-1                                      11d23s21
                  else                                                          11d13s20
                   write(6,*)('bitc not set for nab4(1,1) ='),          2d17s22
     $                  nab4(1,iu1)                                     2d17s22
                   stop 'nab4d(1,1)'                                           11d27s20
                  end if                                                        11d13s20
                  if(btest(itesto,nab4(2,iu2)))then                         12d8s20
                   itestc=ibset(itestc,nab4(2,iu2))                         12d8s20
                   itesto=ibclr(itesto,nab4(2,iu2))                             11d27s20
                   nopenk=nopenk-1                                              11d13s20
                  else if(btest(itesto,nab4(2,iu2)))then                    12d8s20
                   write(6,*)('already double in nab4(2,1) = '),          12d18s20
     $                 nab4(2,iu2)                                       12d18s20
                   stop 'nab4e(2,1)'                                           11d27s20
                  else                                                          11d13s20
                   itesto=ibset(itesto,nab4(2,iu2))                         12d8s20
                   nopenk=nopenk+1                                              11d13s20
                  end if                                                        11d13s20
                  call gandcr(jtc,jto,itestc,itesto,nopenjp,nopenk,       11d22s21
     $         norbxxxx,nnot1,nab1,icode,imap,nx1,irw1,irw2,iwpb1,iwpk1,11d14s22
     $                 bc,ibc)                                          11d14s22
                  call gandcr(itestc,itesto,itc,ito,nopenk,nopenip,       2d23s21
     $         norbxxxx,nnot2,nab2,icode,imap,nx1,irw1,irw2,iwpb2,iwpk2,11d14s22
     $                 bc,ibc)                                          11d14s22
                  if(nnot1.eq.2.and.nnot2.eq.2)then                       12d19s20
                   call spinloop(i2sb,i2smb,i2sk,i2smk,nopenjp,nopenip,  11d22s21
     $                nopenk,ncsfb(jarg),ncsfk(iarg),itype,imatx,ntypeq,11d22s21
     $                nab1,iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(iunit),      11d22s21
     $                ncsfk(iarg),ieoro,bc,ibc)                         11d14s22
                   do ii=0,ntypeq-1                                      11d22s21
                    ipack=ibc(itype+ii)                                  11d22s21
                    jmat=imatx+ii*ncsfb(jarg)*ncsfk(iarg)                11d30s21
                    do j=1,4                                             11d30s21
                     if(ipack2(j).gt.0)then                              11d30s21
                      imy(j)=0                                           11d30s21
                     else                                                11d30s21
                      imy(j)=1                                           11d30s21
                     end if                                              11d30s21
                    end do                                               11d30s21
                    ltest=imy(1)+2*(imy(2)+2*(imy(3)+2*imy(4)))          11d30s21
                    jtest=1-imy(1)+2*(1-imy(2)+2*(1-imy(3)              2d17s22
     $                   +2*(1-imy(4))))                                2d17s22
                    itestp=ltest+1                                       11d30s21
                    jtestp=jtest+1                                       11d30s21
                    phs=1d0                                              11d30s21
                    iuse=-1                                              11d30s21
                    if(iifmx(itestp).ge.0)then                           11d30s21
                     iuse=iifmx(itestp)                                  11d30s21
                    else if(iifmx(jtestp).ge.0)then                      11d30s21
                     iuse=iifmx(jtestp)                                  11d30s21
                     if(isopt(4).ne.0)phs=-phs                           11d30s21
                    end if                                               11d30s21
                    if(iuse.ge.0)then                                    11d30s21
                     do lk=1,4                                           11d22s21
                      if(nli(lk).gt.0)then                               11d22s21
                       iad1=ivcv+iff22k(ivcv+5+lk)                       11d22s21
                       iad2=iad1+nli(lk)                                 11d22s21
                       itmpx=ibcoff                                      11d22s21
                       itmpy=itmpx+nli(lk)*ncsfb(jarg)                   11d22s21
                       ibcoff=itmpy+nli(lk)*ncsfb(jarg)                  11d30s21
                       call dgemm('n','n',ncsfb(jarg),nli(lk),
     $                   ncsfk2(lk,iarg),phs,bc(jmat),ncsfb(jarg),      11d22s21
     $                   ff22k(iad2),ncsfk2(lk,iarg),0d0,bc(itmpx),     11d30s21
     $                   ncsfb(jarg),                                   11d30s21
     d' hcddjkbk4. 18')
                       do j=0,ncsfb(jarg)-1                              11d30s21
                        do k=0,nli(lk)-1                                 11d30s21
                         kj=itmpy+k+nli(lk)*j                            11d30s21
                         jk=itmpx+j+ncsfb(jarg)*k                        11d30s21
                         bc(kj)=bc(jk)                                   11d30s21
                        end do                                           11d30s21
                       end do                                            11d30s21
                       jtmpy=itmpy                                       11d30s21
                       do lb=1,4                                         11d30s21
                        if(nlj(lb).gt.0)then                             11d30s21
                         jad1=jvcv+iff22b(jvcv+5+lb)                     11d30s21
                         jad2=jad1+nlj(lb)                               11d30s21
                         itmpz=ibcoff                                    11d30s21
                         ibcoff=itmpz+nli(lk)*nlj(lb)                    11d30s21
                         call enough('hcddjkbk4. 22',bc,ibc)
                         call dgemm('n','n',nli(lk),nlj(lb),             11d30s21
     $                    ncsfb2(lb,jarg),1d0,bc(jtmpy),nli(lk),        11d30s21
     $                    ff22b(jad2),ncsfb2(lb,jarg),                  11d30s21
     $                    0d0,bc(itmpz),nli(lk),                        11d30s21
     d' hcddjkbk4. 19')
                         ktmpn=itmpz                                     11d30s21
                         jden=iden4(lk,lb)+nfdatk(2,lk,isb)              11d30s21
     $                        *nfdatb(2,lb,jsb)*iuse                     11d30s21
                         do iii=0,nlj(lb)-1                              11d30s21
                          jj=iff22b(jad1+iii)-1                          11d30s21
                          jjden=jden+nfdatk(2,lk,isb)*jj-1               11d30s21
                          do j=0,nli(lk)-1                               11d30s21
                           kk=iff22k(iad1+j)                             11d30s21
                           bc(jjden+kk)=bc(jjden+kk)+bc(ktmpn+j)         11d30s21
                          end do                                         11d30s21
                          ktmpn=ktmpn+nli(lk)                            11d30s21
                         end do                                          11d30s21
                         ibcoff=itmpz                                    11d30s21
                        end if                                           11d30s21
                        jtmpy=jtmpy+nli(lk)*ncsfb2(lb,jarg)              11d30s21
                       end do                                            11d30s21
                       ibcoff=itmpx                                      11d30s21
                      end if                                             11d30s21
                      jmat=jmat+ncsfb(jarg)*ncsfk2(lk,iarg)              11d30s21
                     end do                                              11d30s21
                    end if                                               11d30s21
                   end do                                                11d30s21
                   ibcoff=itype                                          11d30s21
                  end if                                                2d17s22
                 end if                                                 11d30s21
                else if(l2e.eq.0)then                                   2d17s22
                 nnot=0                                                  2d6s23
                 if(ndifs.eq.4.and.ndifb.eq.4)then                                           10d14s2
                  nnot=4                                                          10d14s22
                  ioxx(1)=1                                                     10d17s22
                  ioxx(2)=1                                                     10d17s22
                  do i=1,norbxxx                                                     10d17s22
                   if(btest(gandcb,i))then                                         10d14s22
                    if((btest(itc,i).and.btest(jto,i)).or.                         10d17s22
     $               (btest(ito,i).and..not.btest(jtc,i)))then                   10d14s22
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
                  do i=1,norbxxx
                   if(btest(gandcb,i))then                                         10d14s22
                    if(btest(gandcc,i).and.                                        10d14s22
     $                ((btest(jtc,i).and..not.btest(ito,i)).or.                 10d14s22
     $        (btest(itc,i).and..not.btest(jto,i))))then                     10d14s22
                     if(btest(itc,i))iswap=1                                        10d17s22
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
                   if(btest(itc,nab4(2,2)).and.                         2d6s23
     $                  .not.btest(itc,nab4(2,1)))nbt=1                 2d6s23
                  else                                                             10d17s22
                   nbt=0                                                           10d17s22
                   if(btest(jtc,nab4(1,2)).and.                         2d6s23
     $                  .not.btest(jtc,nab4(1,1)))nbt=1                 2d6s23
                  end if                                                           10d17s22
                  if(nbt.ne.0)then                                                 10d17s22
                   nab4(1,1)=nab4(1,2)                                             10d17s22
                   nab4(2,1)=nab4(2,2)                                             10d17s22
                  end if                                                           10d17s22
                 else if(ndifs.eq.0.and.ndifd.eq.2)then                            10d14s22
                  nnot=3
                  do i=1,norbxxx
                   if(btest(gandcb,i))then                                         10d14s22
                    if(btest(jtc,i))then
                     nab4(1,1)=i
                     nab4(1,2)=i
                    else                                                           10d14s22
                     nab4(2,1)=i
                     nab4(2,2)=i
                    end if
                   end if                                                          10d14s22
                  end do
                 end if                                                            10d14s22
                 if(nnot.ne.0)then                                      2d8s23
                  iu1=1
                  iu2=1
                  itestc=jtc                                             11d23s21
                  itesto=jto                                             11d23s21
                  if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                   itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                   itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopenjp+1                                      11d23s21
                  else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                   itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopenjp-1                                      11d23s21
                  else                                                          11d13s20
                   write(6,*)('bitc not set for nab4(1,1) ='),          2d8s23
     $                  nab4(1,iu1)                                     2d8s23
                   stop 'nab4f(1,1)'                                           11d27s20
                  end if                                                        11d13s20
                  if(btest(itesto,nab4(2,iu2)))then                         12d8s20
                   itestc=ibset(itestc,nab4(2,iu2))                         12d8s20
                   itesto=ibclr(itesto,nab4(2,iu2))                             11d27s20
                   nopenk=nopenk-1                                              11d13s20
                  else if(btest(itesto,nab4(2,iu2)))then                    12d8s20
                   write(6,*)('already double in nab4(2,1) = '),          12d18s20
     $                 nab4(2,iu2)                                       12d18s20
                   stop 'nab4g(2,1)'                                           11d27s20
                  else                                                          11d13s20
                   itesto=ibset(itesto,nab4(2,iu2))                         12d8s20
                   nopenk=nopenk+1                                              11d13s20
                  end if                                                        11d13s20
                  call gandcr(jtc,jto,itestc,itesto,nopenjp,nopenk,       11d22s21
     $         norbxxx,nnot1,nab1,icode,imap,nx1,irw1,irw2,iwpb1,iwpk1, 11d14s22
     $                bc,ibc)                                           11d14s22
                  call gandcr(itestc,itesto,itc,ito,nopenk,nopenip,       2d23s21
     $         norbxxx,nnot2,nab2,icode,imap,nx1,irw1,irw2,iwpb2,iwpk2, 11d14s22
     $                bc,ibc)                                           11d14s22
                  if(nnot1.eq.2.and.nnot2.eq.2)then                       12d19s20
                   call spinloop(i2sb,i2smb,i2sk,i2smk,nopenjp,nopenip,  11d22s21
     $                nopenk,ncsfb(jarg),ncsfk(iarg),itype,imatx,ntypeq,11d22s21
     $                nab1,iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(iunit),      11d22s21
     $                ncsfk(iarg),ieoro,bc,ibc)                         11d14s22
                   do ii=0,ntypeq-1                                      11d22s21
                    ipack=ibc(itype+ii)                                  11d22s21
                    do j=1,4                                             11d22s21
                     if(ipack2(j).gt.0)then                              11d22s21
                      imy(j)=0                                           11d22s21
                      ipack2a(j)=ipack2(j)                               11d22s21
                     else                                                11d22s21
                      imy(j)=1                                           11d22s21
                      ipack2a(j)=-ipack2(j)                              11d22s21
                     end if                                              11d22s21
                     if(ipack2a(j).le.norb)then                          11d22s21
                      igya(j)=irel(ipack2a(j))-1                         11d22s21
                      isy(j)=ism(ipack2a(j))                             11d22s21
                     end if                                              11d22s21
                    end do                                               11d22s21
                    if(ipack2a(3).eq.norbxxx.and.ipack2a(4).eq.norbx)    11d23s21
     $                   then                                           11d23s21
c
c     (12|vv')=(vv'|12) case
c
                     juse=-1                                             11d22s21
                     ltest=imy(3)+2*(imy(4)+2*(imy(1)+2*imy(2)))         11d23s21
                     jtest=1-imy(3)+2*(1-imy(4)+2*(1-imy(1)              2d17s22
     $                    +2*(1-imy(2)))                                2d17s22
     $                    )                                              11d22s21
                     itestp=ltest+1                                      11d22s21
                     jtestp=jtest+1                                      11d22s21
                     isja=1                                              11d23s21
                     isjb=2                                              11d23s21
                     phsj=1d0                                             11d22s21
                     if(iifmx(itestp).ge.0)then                          11d22s21
                      juse=iifmx(itestp)                                 11d22s21
                     else if(iifmx(jtestp).ge.0)then                     11d22s21
                      juse=iifmx(jtestp)                                 11d22s21
                      if(isopt(4).ne.0)phsj=-phsj                        11d22s21
                     end if                                              11d22s21
c
c     (1v'|v2)=(v2|1v')=+/-(2v|v'1) case
c
                     phsk=-1d0                                           11d22s21
                     if(isopt(2).ne.0)phsk=-phsk                         11d23s21
                     ltest=imy(2)+2*(imy(3)+2*(imy(4)+2*imy(1)))         11d22s21
                     jtest=1-imy(2)+2*(1-imy(3)+2*(1-imy(4)              11d22s21
     $                   +2*(1-imy(1))))                                11d22s21
                     itestp=ltest+1                                      11d22s21
                     jtestp=jtest+1                                      11d22s21
                     iska=2                                               11d23s21
                     iskb=1                                              11d23s21
                     kuse=-1                                             11d22s21
                     if(iifmx(itestp).ge.0)then                          11d22s21
                      kuse=iifmx(itestp)                                 11d22s21
                     else if(iifmx(jtestp).ge.0)then                     11d22s21
                      kuse=iifmx(jtestp)                                 11d22s21
                      if(isopt(4).ne.0)phsk=-phsk                        11d22s21
                     end if                                              11d22s21
                    else                                                 11d22s21
                     write(6,*)('don''t know how to handle typed '),
     $                   ipack2
                     call dws_synca
                     call dws_finalize
                     stop
                    end if                                               11d22s21
                    if(min(juse,kuse).ge.0)then                          11d22s21
                     jmat=imatx+ii*ncsfb(jarg)*ncsfk(iarg)               11d22s21
                     jcol=igya(isja)+irefo(isy(isja))*(igya(isjb)        11d23s21
     $                   +irefo(isy(isjb))*juse)                        11d23s21
                     kcol=igya(iska)+irefo(isy(iska))*(igya(iskb)        11d23s21
     $                   +irefo(isy(iskb))*kuse)                        11d23s21
                     do lk=1,4                                           11d22s21
                      if(nli(lk).gt.0)then                               11d22s21
                       iad1=ivcv+iff22k(ivcv+5+lk)                       11d22s21
                       iad2=iad1+nli(lk)                                 11d22s21
                       itmpx=ibcoff                                      11d22s21
                       itmpy=itmpx+nli(lk)*ncsfb(jarg)                   11d22s21
                       ibcoff=itmpy+nli(lk)*ncsfb(jarg)                  11d22s21
                       call dgemm('n','n',ncsfb(jarg),nli(lk),
     $                      ncsfk2(lk,iarg),1d0,bc(jmat),ncsfb(jarg),      11d22s21
     $                      ff22k(iad2),ncsfk2(lk,iarg),0d0,bc(itmpx),  11d23s21
     $                      ncsfb(jarg),                                11d22s21
     d' hcddjkbk4. 20')
                       do j=0,ncsfb(jarg)-1                              11d22s21
                        do k=0,nli(lk)-1                                 11d22s21
                         kj=itmpy+k+nli(lk)*j                            11d22s21
                         jk=itmpx+j+ncsfb(jarg)*k                        11d22s21
                         bc(kj)=bc(jk)                                   11d22s21
                        end do                                           11d22s21
                       end do                                            11d22s21
                       jtmpy=itmpy                                       11d22s21
                       do lb=1,4                                         11d22s21
                        if(nlj(lb).gt.0)then                             11d22s21
                         jad1=jvcv+iff22b(jvcv+5+lb)                     11d22s21
                         jad2=jad1+nlj(lb)                               11d22s21
                         itmpz=ibcoff                                    11d26s21
                         ibcoff=itmpz+nli(lk)*nlj(lb)                    11d26s21
                         call enough('hcddjkbk4. 23',bc,ibc)
                         call dgemm('n','n',nli(lk),nlj(lb),
     $                        ncsfb2(lb,jarg),1d0,bc(jtmpy),nli(lk),    11d22s21
     $                       ff22b(jad2),ncsfb2(lb,jarg),               11d23s21
     $                      0d0,bc(itmpz),nli(lk),                      11d22s21
     d' hcddjkbk4. 21')
                         if(juse.ge.0)then                               11d23s21
                          ktmpn=itmpz                                    11d26s21
                          jden=idenjn(lk,lb,isy(isja))                   11d23s21
     $                         +nfdatk(2,lk,isb)*nfdatb(2,lb,jsb)*jcol      11d22s21
                          do iii=0,nlj(lb)-1                                    2d24s21
                           jj=iff22b(jad1+iii)-1                               8d26s21
                           jjden=jden+nfdatk(2,lk,isb)*jj-1              11d23s21
                           do j=0,nli(lk)-1                                     2d24s21
                            kk=iff22k(iad1+j)                                  8d26s21
                            bc(jjden+kk)=bc(jjden+kk)+bc(ktmpn+j)*phsj     11d22s21
                           end do                                              2d24s21
                           ktmpn=ktmpn+nli(lk)                                  2d24s21
                          end do                                               2d24s21
                         end if                                          11d22s21
                         if(kuse.ge.0)then                               11d26s21
                          ktmpn=itmpz                                    11d26s21
                          jden=idenkn(lk,lb,isy(iska))                   11d23s21
     $                         +nfdatk(2,lk,isb)*nfdatb(2,lb,jsb)*kcol      11d22s21
                          do iii=0,nlj(lb)-1                                    2d24s21
                           jj=iff22b(jad1+iii)-1                               8d26s21
                           jjden=jden+nfdatk(2,lk,isb)*jj-1              11d23s21
                           do j=0,nli(lk)-1                                     2d24s21
                            kk=iff22k(iad1+j)                                  8d26s21
                            bc(jjden+kk)=bc(jjden+kk)+bc(ktmpn+j)*phsk   11d22s21
                           end do                                              2d24s21
                           ktmpn=ktmpn+nli(lk)                                  2d24s21
                          end do                                               2d24s21
                         end if                                          11d22s21
                         ibcoff=itmpz                                    11d26s21
                        end if                                           11d26s21
                        jtmpy=jtmpy+nli(lk)*ncsfb2(lb,jarg)              11d22s21
                       end do                                            11d22s21
                       ibcoff=itmpx                                      11d22s21
                      end if                                             11d26s21
                      jmat=jmat+ncsfb(jarg)*ncsfk2(lk,iarg)              11d22s21
                     end do                                              11d22s21
                    end if                                               11d22s21
                   end do                                                11d22s21
                   ibcoff=itype                                          11d22s21
                  end if                                                 11d22s21
                 end if                                                 2d8s23
                end if                                                  11d22s21
                ibcoff=iunit                                            11d22s21
               end if                                                    1d8s21
               jvcv=jvcv+nspacej                                         3d19s21
              end do                                                       1d8s21
             end if                                                     3d1s21
             idoit=idoit+1                                               3d1s21
            end if                                                      1d8s21
           end do                                                       1d8s21
           ivcv=ivcv+nspacei                                            3d19s21
          end do                                                        1d8s21
          ibcoff=ibcst                                                  1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
        call dws_gsumf(bc(ibc0),nwds)                                   2d23s21
        do l=1,4                                                        1d8s21
         if(nfdatk(2,l,isb).gt.0)then                                   8d26s21
          do ll=1,4                                                     1d8s21
           if(nfdatb(2,ll,jsb).gt.0)then                                 1d8s21
            nll=nfdatk(2,l,isb)*nfdatb(2,ll,jsb)                        8d26s21
            do lsa=1,nsymb                                              1d8s21
             lsb=multh(lsa,ijsb)                                        11d22s21
             nn=irefo(lsa)*irefo(lsb)                                   1d9s21
             if(min(nn,nll).gt.0)then                                   2d24s21
              idenkf(l,ll,lsa)=idenkn(l,ll,lsa)                         2d25s21
             end if                                                     2d24s21
            end do                                                      1d8s21
           end if                                                       1d8s21
          end do                                                        1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
        nall=0                                                          5d10s21
        do ll=1,4                                                       1d8s21
         if(nfdatb(2,ll,jsb).gt.0)then                                  8d26s21
          do l=1,4                                                      1d8s21
           if(nfdatk(2,l,isb).gt.0)then                                 8d26s21
            nnn=nfdatk(2,l,isb)*nfdatb(2,ll,jsb)                        11d23s21
            sz=0d0                                                      11d23s21
            do i=0,nnn-1                                                11d23s21
             sz=sz+bc(iden1en(l,ll)+i)**2                               11d23s21
            end do                                                      11d23s21
            sz=sqrt(sz/dfloat(nfdatk(2,l,isb)*nfdatb(2,ll,jsb)))        11d23s21
            if(sz.gt.1d-10)then                                         11d23s21
             nall=nall+1                                                11d23s21
            end if                                                      11d23s21
            sz=0d0                                                      8d21s21
            do i=0,nnn-1                                                11d22s21
             sz=sz+bc(idenhvvn(l,ll)+i)**2                                 8d21s21
            end do                                                      8d21s21
            sz=sqrt(sz/dfloat(nfdatk(2,l,isb)*nfdatb(2,ll,jsb)))        8d26s21
            if(sz.gt.1d-10)then                                         11d23s21
             nall=nall+1                                                11d23s21
             if(isb.eq.jsb)then                                         2d7s22
              nhvii=nhvii+1                                             2d7s22
             else                                                       2d7s22
              nhvij=nhvij+1                                             2d7s22
             end if                                                     2d7s22
            end if                                                      11d23s21
            do lsa=1,nsymb                                              1d8s21
             lsb=multh(lsa,ijsb)                                        11d22s21
             nokk(l,ll,lsa)=0                                           1d15s21
             if(min(irefo(lsa),irefo(lsb)).gt.0)then                    1d13s21
              nnn=nfdatk(2,l,isb)*nfdatb(2,ll,jsb)                      8d26s21
              nn=irefo(lsa)*irefo(lsb)*ntype                            11d22s21
              jden=idenjn(l,ll,lsa)                                     11d23s21
              nokj(l,ll,lsa)=0                                          11d23s21
              kden=idenjn(l,ll,lsa)                                     11d23s21
              do icol=0,nn-1                                            1d9s21
               sz=0d0                                                    1d8s21
               do i=0,nnn-1                                             1d9s21
                sz=sz+bc(jden+i)**2                                     1d13s21
               end do                                                    1d8s21
               sz=sqrt(sz/dfloat(nnn))                                  1d9s21
               if(sz.gt.1d-10.and.l2e.eq.0)then                         2d17s22
                ibc(ndenj(l,ll,lsa)+nokj(l,ll,lsa))=icol                11d23s21
                nokj(l,ll,lsa)=nokj(l,ll,lsa)+1                         11d23s21
                do i=0,nnn-1                                            1d12s21
                 bc(kden+i)=bc(jden+i)                                  1d12s21
                end do                                                  1d12s21
                kden=kden+nnn                                           1d12s21
                nall=nall+1                                             5d10s21
               end if                                                   1d13s21
               jden=jden+nnn                                            1d9s21
              end do                                                    1d9s21
              nn=irefo(lsa)*irefo(lsb)*ntype                            11d22s21
              jden=idenkf(l,ll,lsa)                                     1d9s21
              kden=idenkf(l,ll,lsa)                                     1d12s21
              do icol=0,nn-1                                            1d9s21
               sz=0d0                                                    1d8s21
               do i=0,nnn-1                                             1d9s21
                sz=sz+bc(jden+i)**2                                     1d9s21
               end do                                                    1d8s21
               sz=sqrt(sz/dfloat(nnn))                                  1d9s21
               if(sz.gt.1d-10.and.l2e.eq.0)then                         2d17s22
                ibc(ndenkf(l,ll,lsa)+nokk(l,ll,lsa))=icol               1d12s21
                nokk(l,ll,lsa)=nokk(l,ll,lsa)+1                         1d12s21
                do i=0,nnn-1                                            1d12s21
                 bc(kden+i)=bc(jden+i)                                  1d12s21
                end do                                                  1d12s21
                kden=kden+nnn                                           1d12s21
                nall=nall+1                                             5d10s21
               end if                                                    1d8s21
               jden=jden+nnn                                            1d9s21
              end do                                                    1d9s21
             end if                                                     1d8s21
            end do                                                      1d8s21
           end if                                                       1d8s21
          end do                                                        1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
        if(jsb.eq.isb)then                                              12d3s21
         nden4=1                                                        12d3s21
        else                                                            12d3s21
         nden4=0                                                        12d3s21
        end if                                                          12d3s21
        if(max(nden4,nall).gt.0)then                                    2d23s22
         ioff=ioffdnon                                                   1d12s21
         if(jsb.eq.isb)then                                             12d3s21
          do lb=1,4                                                      12d3s21
           idv4l(1,lb)=idv4(1,lb,jsb)                                   12d3s21
           idv4l(2,lb)=idv4(2,lb,jsb)                                   12d3s21
          end do                                                         12d3s21
         end if
         do isbv1=1,nsymb                                                1d12s21
          isbv2=multh(isbv1,isbv12)                                      1d12s21
          if(isbv2.ge.isbv1)then                                         1d12s21
           if(isbv12.eq.1)then                                           1d12s21
            isw=0                                                        1d12s21
            mvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        1d12s21
            ioffp=ioff+nvirt(isbv1)*nrootu*nfdatk(2,1,isb)              8d26s21
           else                                                          1d12s21
            isw=1                                                        1d12s21
            mvv=nvirt(isbv1)*nvirt(isbv2)                                1d12s21
            ioffp=ioff
           end if                                                        1d12s21
           mmvv=mvv*nrootu                                               1d12s21
           if(nden4.ne.0.and.l2e.eq.0)then                              2d23s22
            igoal4va=1
            igoal4vb=1
c
c     vectors are v,r,nfdatk
c     density is nfdatk,nfdatb,type, so
c     v*d is v,r,nfdatb,type
            do lb=1,4                                                   11d30s21
             nrow=mmvv                                                  11d30s21
             ncol=nfdatb(2,lb,jsb)*ntype                                11d30s21
             if(isw.eq.0)then                                           11d30s21
              nrowx=nvirt(isbv1)*nrootu                                 11d30s21
             end if                                                     11d30s21
             call enough('hcddjkbk4. 24',bc,ibc)
             if(nfdatb(2,lb,jsb).gt.0)then                              11d30s21
              fact=0d0                                                  11d30s21
              iioffp=ioffp                                              12d2s21
              do lk=1,4                                                 11d30s21
               if(nfdatk(2,lk,isb).gt.0)then                            11d30s21
                if(lk.eq.1.and.isw.eq.0.and.nrowx.gt.0)then             11d30s21
c     vd is nvv,nroot,nfdat
                 call dgemm('n','n',nrowx,ncol,nfdatk(2,lk,isb),srh,    11d30s21
     $                vd(ioff),nrowx,bc(iden4(lk,lb)),nfdatk(2,lk,isb), 11d30s21
     $                0d0,bc(idv4l(2,lb)),nrowx,                        12d3s21
     d' hcddjkbk4. 22')
c     dv4 nvv,nroot,ntype
                end if                                                  11d30s21
                if(nrow.gt.0)then                                       11d30s21
                 call dgemm('n','n',nrow,ncol,nfdatk(2,lk,isb),1d0,      11d30s21
     $                vd(iioffp),nrow,bc(iden4(lk,lb)),nfdatk(2,lk,isb),12d2s21
     $                fact,bc(idv4l(1,lb)),nrow,                        12d3s21
     d' hcddjkbk4. 23')
                 fact=1d0                                                11d30s21
                end if                                                  11d30s21
                iioffp=iioffp+mmvv*nfdatk(2,lk,isb)                     12d2s21
               end if                                                   11d30s21
              end do                                                    11d30s21
             end if                                                     11d30s21
            end do                                                      11d30s21
           end if                                                       11d30s21
           joff=joffdnon                                                 1d12s21
           do jsbv1=1,nsymb                                              1d12s21
            jsbv2=multh(jsbv1,jsbv12)                                    1d12s21
            if(jsbv2.ge.jsbv1)then                                       1d12s21
             if(jsbv12.eq.1)then                                         1d12s21
              jsw=0                                                      1d12s21
              nvv=(nvirt(jsbv1)*(nvirt(jsbv2)-1))/2                      1d12s21
              joffp=joff+nvirt(jsbv1)*nrootu*nfdatb(2,1,jsb)            8d26s21
             else                                                        1d12s21
              jsw=1                                                      1d12s21
              nvv=nvirt(jsbv1)*nvirt(jsbv2)                              1d12s21
              joffp=joff                                                 1d13s21
             end if                                                      1d12s21
             nnvv=nvv*nrootu                                             1d12s21
c
c     my stab at 4v ...                                                 11d30s21
c     Gv"v'''rj=djit [(v"v|v'''v')(t)-(v"v'|v'''v)(t)]Vvv'ri          11d30s21
c                      j i j   i       j i  j   i
c                      1 1 2   2       1 2  2   1
c                  = [(v"v|v'''v')(t)-(v"v'|v'''v)(t)]Dvv'rjt          11d30s21
c                      j i j   i       j i  j   i
c                      1 1 2   2       1 2  2   1
c     with Dvv'rjt=Vvv'ri*djit
c
             if(nden4.ne.0.and.l2e.eq.0)then                            2d23s22
              call ilimts(nh0(jsbv2),nh0(isbv2),mynprocg,mynowprog,     11d30s21
     $            jl,jh,j1s,j1e,j2s,j2e)                                11d30s21
              nherejx=nh0(jsbv2)*nh0(isbv2)
              do it=0,ntype-1                                           11d30s21
               i10=j1s                                                  11d30s21
               i1n=nh0(jsbv2)                                           2d10s22
               nrowj=nh0(jsbv1)*nh0(isbv1)                              11d30s21
cc
               if(iuall.ne.0)then                                       1d24s22
                j4=iall(isbv1,jsbv2,isbv2)+nrowj*nherejx*it               11d30s21
                do i2=j2s,j2e                                            11d30s21
                 if(i2.eq.j2e)i1n=j1e                                    11d30s21
                 iv2=i2-1-irefo(isbv2)                                   12d3s21
                 if(iv2.ge.0)then                                        12d3s21
                  iv1t=(iv2+isw*(nvirt(isbv1)-iv2))-1                     11d30s21
                  irec=nvirt(isbv1)*iv2                                   11d30s21
                  itri=((iv2*(iv2-1))/2)                                  11d30s21
                  itri=itri+isw*(irec-itri)                               11d30s21
                  do i1=i10,i1n                                           11d30s21
                   jv2=i1-1-irefo(jsbv2)                                 12d3s21
                   if(jv2.ge.0)then                                      12d3s21
                    icol=jv2+irefo(jsbv2)+nh0(jsbv2)*(iv2+irefo(isbv2)) 2d8s22
     $                  +1                                              2d8s22
                    j4u=j4+nrowj*(icol-1)                                 12d3s21
                    jv1t=(jv2+jsw*(nvirt(jsbv1)-jv2))-1                    11d30s21
                    jrec=nvirt(jsbv1)*jv2                                  11d30s21
                    jtri=((jv2*(jv2-1))/2)                                 11d30s21
                    jtri=jtri+jsw*(jrec-jtri)                              11d30s21
                    if(jsw.eq.0)then                                      11d30s21
                     idens=idv4l(1,1)+mmvv*nfdatb(2,1,jsb)*it             12d3s21
                     do jr=0,nrootu*nfdatb(2,1,jsb)-1                     11d30s21
                      ifrm=idens+itri+mvv*jr                               11d30s21
                      ito=joff+jv2+nvirt(jsbv1)*jr                        11d30s21
                      do iv1=0,iv1t                                          11d30s21
                       jj4=j4u+jv2+irefo(jsbv1)                           11d30s21
     $                     +nh0(jsbv1)*(iv1+irefo(isbv1))               11d30s21
                       gd(ito)=gd(ito)+bc(ifrm+iv1)*bc(jj4)*srh           11d30s21
                      end do                                              11d30s21
                     end do                                               11d30s21
                     if(isw.eq.0)then                                     11d30s21
                      idenx=idv4l(2,1)+nvirt(isbv1)*nfdatb(2,1,jsb)*it    12d3s21
     $                   *nrootu                                        2d3s22
                      do jr=0,nrootu*nfdatb(2,1,jsb)-1                     11d30s21
                       ifrm=idenx+iv2+nvirt(isbv1)*jr                     11d30s21
                       ito=joff+jv2+nvirt(jsbv1)*jr                        11d30s21
                       jj4=j4u+jv2+irefo(jsbv1)                           11d30s21
     $                     +nh0(jsbv1)*(iv2+irefo(isbv1))               11d30s21
                       gd(ito)=gd(ito)+bc(ifrm)*bc(jj4)*srh               11d30s21
                      end do                                              11d30s21
                     end if                                               11d30s21
                    end if                                                11d30s21
                    jjoffp=joffp                                          12d2s21
                    do lb=1,4                                              11d30s21
                     idens=idv4l(1,lb)+mmvv*nfdatb(2,lb,jsb)*it           12d3s21
                     do jr=0,nrootu*nfdatb(2,lb,jsb)-1                     11d30s21
                      ifrm=idens+itri+mvv*jr                               11d30s21
                      ito=jjoffp+jtri+nvv*jr                              12d2s21
                      do iv1=0,iv1t                                          11d30s21
                       jj4=j4u+irefo(jsbv1)+nh0(jsbv1)*(iv1             2d8s22
     $                      +irefo(isbv1))                              2d8s22
                       do jv1=0,jv1t                                         11d30s21
                        gd(ito+jv1)=gd(ito+jv1)+bc(ifrm+iv1)*bc(jj4+jv1)   11d30s21
                       end do                                             11d30s21
                      end do                                              11d30s21
                     end do                                               11d30s21
                     if(isw.eq.0)then                                     11d30s21
                      idenx=idv4l(2,lb)+nvirt(isbv1)*nfdatb(2,lb,jsb)*it  12d3s21
     $                   *nrootu                                        2d3s22
                      do jr=0,nrootu*nfdatb(2,lb,jsb)-1                     11d30s21
                       ifrm=idenx+iv2+nvirt(isbv1)*jr                     11d30s21
                       ito=jjoffp+jtri+nvv*jr                             12d2s21
                       jj4=j4u+irefo(jsbv1)+nh0(jsbv1)*(iv2             2d8s22
     $                      +irefo(isbv1))                              2d8s22
                       do jv1=0,jv1t                                         11d30s21
                        gd(ito+jv1)=gd(ito+jv1)+bc(ifrm)*bc(jj4+jv1)      11d30s21
                       end do                                             11d30s21
                      end do                                              11d30s21
                     end if                                               11d30s21
                     jjoffp=jjoffp+nnvv*nfdatb(2,lb,jsb)                  12d2s21
                    end do                                                11d30s21
                   end if                                                 11d30s21
                  end do                                                  11d30s21
                 end if                                                  12d3s21
                 i10=1                                                   11d30s21
                end do                                                   11d30s21
               end if
cc
              end do                                                    11d30s21
              call ilimts(nh0(jsbv2),nh0(isbv1),mynprocg,mynowprog,     11d30s21
     $            kl,kh,k1s,k1e,k2s,k2e)                                11d30s21
              nherekx=nh0(jsbv2)*nh0(isbv1)                             12d3s21
              if(iuall.ne.0)then                                        2d8s22
              end if                                                    2d8s22
              do it=0,ntype-1                                           11d30s21
               itu=-1                                                   12d2s21
               phs=1d0                                                  12d2s21
               do jj=1,16                                               12d2s21
                if(iifmx(jj).eq.it)then                                 12d2s21
                 imy(4)=(jj-1)/8                                        12d2s21
                 ltest=jj-1-8*imy(4)                                    12d2s21
                 imy(3)=ltest/4                                         12d2s21
                 ltest=ltest-4*imy(3)                                   12d2s21
                 imy(2)=ltest/2                                         12d2s21
                 imy(1)=ltest-2*imy(2)                                  12d2s21
                 ltest=imy(1)+2*(imy(4)+2*(imy(3)+2*imy(2)))            12d2s21
                 jtest=1-imy(1)+2*(1-imy(4)+2*(1-imy(3)+2*(1-imy(2))))  12d2s21
                 itestp=ltest+1                                         12d2s21
                 jtestp=jtest+1                                         12d2s21
                 if(iifmx(itestp).ge.0)then                             12d2s21
                  itu=iifmx(itestp)                                     12d2s21
                 else if(iifmx(jtestp).ge.0)then                        12d2s21
                  itu=iifmx(jtestp)                                     12d2s21
                  if(isopt(4).ne.0)phs=-phs                             12d2s21
                 end if                                                 12d2s21
                 go to 1887                                             12d2s21
                end if                                                  12d2s21
               end do                                                   12d2s21
 1887          continue                                                 12d2s21
               if(itu.ge.0)then                                         12d2s21
                i10=k1s                                                  11d30s21
                i1n=nh0(jsbv2)                                          2d10s22
                nrowk=nh0(jsbv1)*nh0(isbv2)                              11d30s21
cc
                if(iuall.ne.0)then
                 k4=iall(isbv2,jsbv2,isbv1)+nrowk*nherekx*itu            12d3s21
                 do i2=k2s,k2e                                            11d30s21
                  if(i2.eq.k2e)i1n=k1e                                    11d30s21
                  iv1=i2-1-irefo(isbv1)                                  12d3s21
                  if(iv1.ge.0)then                                       12d3s21
                   iv2b=(iv1+1)*(1-isw)                                    11d30s21
                   do i1=i10,i1n                                           11d30s21
                    jv2=i1-1-irefo(jsbv2)                                12d3s21
                    if(jv2.ge.0)then
                     icol=jv2+irefo(jsbv2)                               12d3s21
     $                   +nh0(jsbv2)*(iv1+irefo(isbv1))+1               12d3s21
                     k4u=k4+nrowk*(icol-1)                                12d3s21
                     jv1t=(jv2+jsw*(nvirt(jsbv1)-jv2))-1                    11d30s21
                     jrec=nvirt(jsbv1)*jv2                                  11d30s21
                     jtri=((jv2*(jv2-1))/2)                                 11d30s21
                     jtri=jtri+jsw*(jrec-jtri)                              11d30s21
                     if(jsw.eq.0)then                                      11d30s21
                      idens=idv4l(1,1)+mmvv*nfdatb(2,1,jsb)*it            12d3s21
                      do jr=0,nrootu*nfdatb(2,1,jsb)-1                     11d30s21
                       ito=joff+jv2+nvirt(jsbv1)*jr                        11d30s21
                       do iv2=iv2b,nvirt(isbv2)-1                          11d30s21
                        irec=iv1+nvirt(isbv1)*iv2                          11d30s21
                        itri=((iv2*(iv2-1))/2)+iv1                         11d30s21
                        itri=itri+isw*(irec-itri)                          11d30s21
                        ifrm=idens+itri+mvv*jr                             11d30s21
                        kk4=k4u+jv2+irefo(jsbv1)                           11d30s21
     $                     +nh0(jsbv1)*(iv2+irefo(isbv2))               11d30s21
                        gd(ito)=gd(ito)-bc(ifrm)*bc(kk4)*srh*phs           12d2s21
                       end do                                              11d30s21
                      end do                                               11d30s21
                      if(isw.eq.0)then                                     11d30s21
                       idenx=idv4l(2,1)+nvirt(isbv1)*nfdatb(2,1,jsb)*it   12d3s21
     $                    *nrootu                                       2d3s22
                       do jr=0,nrootu*nfdatb(2,1,jsb)-1                     11d30s21
                        ifrm=idenx+iv1+nvirt(isbv1)*jr                     11d30s21
                        ito=joff+jv2+nvirt(jsbv1)*jr                        11d30s21
                        kk4=k4u+jv2+irefo(jsbv1)                           11d30s21
     $                     +nh0(jsbv1)*(iv1+irefo(isbv2))               11d30s21
                        gd(ito)=gd(ito)-bc(ifrm)*bc(kk4)*srh*phs           12d2s21
                       end do                                              11d30s21
                      end if                                               11d30s21
                     end if                                                11d30s21
                     jjoffp=joffp                                         12d2s21
                     do lb=1,4                                              11d30s21
                      idens=idv4l(1,lb)+mmvv*nfdatb(2,lb,jsb)*it          12d3s21
                      do jr=0,nrootu*nfdatb(2,lb,jsb)-1                     11d30s21
                       ito=jjoffp+jtri+nvv*jr                             12d2s21
                       do iv2=iv2b,nvirt(isbv2)-1                          11d30s21
                        irec=iv1+nvirt(isbv1)*iv2                          11d30s21
                        itri=((iv2*(iv2-1))/2)+iv1                         11d30s21
                        itri=itri+isw*(irec-itri)                          11d30s21
                        ifrm=idens+itri+mvv*jr                               11d30s21
                        kk4=k4u+irefo(jsbv1)                             12d3s21
     $                      +nh0(jsbv1)*(iv2+irefo(isbv2))              12d3s21
                        do jv1=0,jv1t                                         11d30s21
                         gd(ito+jv1)=gd(ito+jv1)
     $                        -bc(ifrm)*bc(kk4+jv1)*phs                 2d8s22
                        end do                                             11d30s21
                       end do                                              11d30s21
                      end do                                               11d30s21
                      if(isw.eq.0)then                                     11d30s21
                       idenx=idv4l(2,lb)+nvirt(isbv1)*nfdatb(2,lb,jsb)  2d8s22
     $                     *it*nrootu                                   2d8s22
                       do jr=0,nrootu*nfdatb(2,lb,jsb)-1                     11d30s21
                        ifrm=idenx+iv1+nvirt(isbv1)*jr                     11d30s21
                        ito=jjoffp+jtri+nvv*jr                            12d2s21
                        kk4=k4u+irefo(jsbv1)                             12d3s21
     $                      +nh0(jsbv1)*(iv1+irefo(isbv2))              12d3s21
                        do jv1=0,jv1t                                         11d30s21
                         gd(ito+jv1)=gd(ito+jv1)                        2d8s22
     $                        -bc(ifrm)*bc(kk4+jv1)*phs                 2d8s22
                        end do                                             11d30s21
                       end do                                              11d30s21
                      end if                                               11d30s21
                      jjoffp=jjoffp+nnvv*nfdatb(2,lb,jsb)                 12d2s21
                     end do                                                11d30s21
                    end if                                                 11d30s21
                   end do                                                  11d30s21
                  end if                                                 12d3s21
                  i10=1                                                   11d30s21
                 end do                                                   11d30s21
                end if
cc
               end if                                                   12d2s21
              end do                                                    11d30s21
             end if                                                     11d30s21
c
c
c     Gvvrj =djiab (ab|vv") Vvv"ri
c     Gvv'rj=djiab (ab|vv") Vv'v"ri
c     Gvv'rj=djiab (ab|vv") Vv"v'ri
c     Gv'vrj=djiab (ab|vv") Vv'v"ri
c     Gv'vrj=djiab (ab|vv") Vv"v'ri
             if(jsbv2.eq.isbv1.or.jsbv2.eq.isbv2.or.jsbv1.eq.isbv1.or.  2d23s22
     $         jsbv1.eq.isbv2)then                                      2d23s22
              do imatch=1,4                                              1d13s21
               kase=0                                                    1d14s21
               if(jsbv2.eq.isbv1.and.imatch.eq.1)then                    1d13s21
                isbvc=jsbv2                                               1d13s21
                isbl=jsbv1                                                1d13s21
                isbr=isbv2                                                1d13s21
                itf=0                                                     1d13s21
                jtf=0                                                     1d13s21
                ibl=1                                                     1d13s21
                ibr=1                                                     1d13s21
                kase=1
               end if                                                    1d13s21
               if(jsbv2.eq.isbv2.and.imatch.eq.2)then                    1d13s21
                isbvc=jsbv2                                               1d13s21
                isbl=jsbv1                                                1d13s21
                isbr=isbv1                                                1d13s21
                itf=1                                                     1d13s21
                jtf=0                                                    1d15s21
                ibl=1                                                     1d13s21
                ibr=2                                                     1d13s21
                kase=2
               end if                                                    1d13s21
               if(jsbv1.eq.isbv1.and.imatch.eq.3)then                    1d14s21
                isbvc=jsbv1                                               1d13s21
                isbl=jsbv2                                                1d13s21
                isbr=isbv2                                                1d13s21
                jtf=1                                                     1d13s21
                itf=0                                                    1d15s21
                ibl=2                                                     1d13s21
                ibr=1                                                     1d13s21
                kase=3
               end if                                                    1d13s21
               if(jsbv1.eq.isbv2.and.imatch.eq.4)then                    1d13s21
                isbvc=jsbv1                                               1d13s21
                isbl=jsbv2                                                1d13s21
                isbr=isbv1                                                1d13s21
                itf=1                                                     1d13s21
                jtf=1                                                     1d13s21
                ibl=2                                                     1d13s21
                ibr=2                                                     1d13s21
                kase=4                                                    1d13s21
               end if                                                     1d13s21
               iioffp=ioffp                                               1d12s21
               tfi=1d0                                                   1d14s21
               do li=1,4                                                 1d14s21
                if(min(kase,nfdatk(2,li,isb)).gt.0)then                 11d22s21
                 jjoffp=joffp                                            1d14s21
                 tfj=1d0                                                 1d14s21
                 do lj=1,4                                               1d14s21
                  if(nfdatb(2,lj,jsb).gt.0)then                         8d26s21
                   tf=tfi*tfj
                   call ilimts(nvirt(isbl),nvirt(isbr),mynprocg,            1d13s21
     $                mynowprog,il,ih,i1s,i1e,i2s,i2e)                  1d12s21
                   nhere=ih+1-il                                            1d12s21
                   nhere2=i2e+1-i2s                                         1d12s21
                   if(min(nhere,nvirt(isbvc)).gt.0)then                 8d11s21
                    if(min(lj,li).eq.1.and.max(lj,li).gt.1)then                   1d11s21
                     phs=-1d0                                                   1d11s21
                    else                                                        1d11s21
                     phs=+1d0                                                   1d11s21
                    end if                                                      1d11s21
                    mn=nfdatk(2,li,isb)*nfdatb(2,lj,jsb)                8d26s21
                    intden=ibcoff                                            1d12s21
                    ibcoff=intden+mn*nhere                                   1d12s21
                    call enough('hcddjkbk4. 25',bc,ibc)
                    fact=0d0                                                1d13s21
                    if(multh(isbl,isbr).eq.isopt(1))then                11d26s21
                     fact=1d0                                           11d26s21
                     if(ih0n(isbl).gt.0)then                            11d26s21
                      i10=i1s                                           8d10s21
                      i1n=nvirt(isbl)                                   8d10s21
                      jntden=intden                                     11d26s21
                      do i2=i2s,i2e                                     8d10s21
                       i2p=i2+irefo(isbr)-1                             11d26s21
                       ih=ih0n(isbl)+irefo(isbl)-1+nh0(isbl)*i2p        11d26s21
                       if(i2.eq.i2e)i1n=i1e                             8d10s21
                       do i1=i10,i1n                                    8d10s21
                        do j=0,mn-1                                     11d26s21
                         bc(jntden+j)=bc(ih+i1)*bc(idenhvvn(li,lj)+j)    11d26s21
                        end do                                          11d26s21
                        jntden=jntden+mn                                11d26s21
                       end do                                           11d26s21
                       i10=1                                            11d26s21
                      end do                                            11d26s21
                     else                                               11d26s21
                      ffh=1d0                                            11d26s21
                      if(isopt(2).ne.0)ffh=-1d0                         11d26s21
                      if(isopt(4).ne.0.and.isopt(3).ne.0)ffh=-ffh       11d30s21
                      i10=i1s                                           8d10s21
                      i1n=nvirt(isbl)                                   8d10s21
                      jntden=intden                                     11d26s21
                      do i2=i2s,i2e                                     8d10s21
                       i2p=i2+irefo(isbr)-1                             11d26s21
                       ih=ih0n(isbr)+i2p+nh0(isbr)*(irefo(isbl)-1)      11d26s21
                       if(i2.eq.i2e)i1n=i1e                             8d10s21
                       do i1=i10,i1n                                    8d10s21
                        hh=bc(ih+i1*nh0(isbr))*ffh                      11d26s21
                        do j=0,mn-1                                     11d26s21
                         bc(jntden+j)=hh*bc(idenhvvn(li,lj)+j)          11d26s21
                        end do                                          11d26s21
                        jntden=jntden+mn                                11d26s21
                       end do                                           11d26s21
                       i10=1                                            11d26s21
                      end do                                            11d26s21
                     end if                                             11d26s21
                    end if                                              11d26s21
                    do lsa=1,nsymb                                              1d8s21
                     lsb=multh(lsa,ijsb)                                11d22s21
                     if(min(irefo(lsa),irefo(lsb)).gt.0                 2d23s22
     $                    .and.l2e.eq.0)then                            2d23s22
                      if(nokj(li,lj,lsa).gt.0)then                      11d23s21
                       itmpi=ibcoff                                         1d12s21
                       ibcoff=itmpi+nokj(li,lj,lsa)*nhere               11d23s21
                       call enough('hcddjkbk4. 26',bc,ibc)
                       do iz=itmpi,ibcoff-1                               8d10s21
                        bc(iz)=0d0                                        8d10s21
                       end do                                             8d10s21
                       do j=0,nokj(li,lj,lsa)-1                         11d23s21
                        jj=ibc(ndenj(li,lj,lsa)+j)                      11d23s21
                        jtmpi=itmpi+j                                     11d22s21
                        jint=jmatsr(isbr,lsa,lsb)+nhere*jj                11d22s21
                        i10=i1s                                           8d10s21
                        i1n=nvirt(isbl)                                   8d10s21
                        do i2=i2s,i2e                                     8d10s21
                         if(i2.eq.i2e)i1n=i1e                             8d10s21
                         do i1=i10,i1n                                    8d10s21
                          bc(jtmpi)=bc(jtmpi)+bc(jint)                    11d22s21
                          jint=jint+1                                       11d22s21
                          jtmpi=jtmpi+nokj(li,lj,lsa)                   11d23s21
                         end do                                           8d10s21
                         i10=1                                            8d10s21
                        end do                                            8d10s21
                       end do                                             8d10s21
                       call dgemm('n','n',mn,nhere,nokj(li,lj,lsa),1d0, 11d23s21
     $                  bc(idenjn(li,lj,lsa)),mn,bc(itmpi),             11d23s21
     $                      nokj(li,lj,lsa),fact,bc(intden),mn,         11d23s21
     d' hcddjkbk4. 24')
                       fact=1d0                                             1d12s21
                       ibcoff=itmpi                                         1d13s21
                      end if                                                1d12s21
                      if(nokk(li,lj,lsa).gt.0)then                              1d12s21
                       itmpi=ibcoff                                         1d12s21
                       ibcoff=itmpi+nokk(li,lj,lsa)*nhere                       1d12s21
                       call enough('hcddjkbk4. 27',bc,ibc)
                       do iz=itmpi,ibcoff-1                             8d10s21
                        bc(iz)=0d0                                      8d10s21
                       end do                                           8d10s21
                       nn=irefo(lsa)*irefo(lsb)                         11d22s21
                       do j=0,nokk(li,lj,lsa)-1                                1d12s21
                        jj=ibc(ndenkf(li,lj,lsa)+j)                            1d12s21
                        i10=i1s                                         8d10s21
                        i1n=nvirt(isbl)                                 8d10s21
                        jtmpi=itmpi+j                                   8d10s21
                        kint=kmatsr(isbl,isbr,lsb)+nhere*jj             11d22s21
                        do i2=i2s,i2e                                   8d10s21
                         if(i2.eq.i2e)i1n=i1e                           8d10s21
                         do i1=i10,i1n                                  8d10s21
                          bc(jtmpi)=bc(kint)                              11d22s21
                          jtmpi=jtmpi+nokk(li,lj,lsa)                   8d10s21
                          kint=kint+1                                     11d22s21
                         end do                                         8d10s21
                         i10=1                                          8d10s21
                        end do                                          8d10s21
                       end do                                           8d10s21
                       call dgemm('n','n',mn,nhere,nokk(li,lj,lsa),1d0,   1d14s21
     $               bc(idenkf(li,lj,lsa)),mn,bc(itmpi),nokk(li,lj,lsa),1d14s21
     $                    fact,bc(intden),mn,                           1d14s21
     d' hcddjkbk4. 25')
                       fact=1d0                                             1d12s21
                       ibcoff=itmpi                                         1d13s21
                      end if                                                1d12s21
                     end if                                                 1d13s21
                    end do                                                  1d13s21
                    if(fact.gt.0.5d0)then                                    1d12s21
                     nrow=nfdatb(2,lj,jsb)*nvirt(isbl)                  8d26s21
                     nmul=nfdatk(2,li,isb)*nhere2                       8d26s21
                     itmp1=ibcoff                                            1d12s21
                     ibcoff=itmp1+nrow*nmul                                 1d12s21
                     call enough('hcddjkbk4. 28',bc,ibc)
                     do i=itmp1,ibcoff-1                                    1d12s21
                      bc(i)=0d0                                             1d12s21
                     end do                                                 1d12s21
                     i10=i1s                                                1d12s21
                     i1n=nvirt(isbl)                                        1d13s21
                     jntden=intden                                          1d12s21
                     do i2=i2s,i2e                                          1d12s21
                      i2n=i2-i2s                                            1d12s21
                      if(i2.eq.i2e)i1n=i1e                                  1d12s21
                      do i1=i10,i1n                                         1d12s21
                       i1m=i1-1                                             1d12s21
                       do j=0,nfdatb(2,lj,jsb)-1                        8d26s21
                        do i=0,nfdatk(2,li,isb)-1                       8d26s21
                         iad=itmp1+j+nfdatb(2,lj,jsb)*(i1m+nvirt(isbl)*(8d26s21
     $                       i2n+nhere2*i))                             1d14s21
                         bc(iad)=bc(jntden+i)                               1d13s21
                        end do                                              1d12s21
                        jntden=jntden+nfdatk(2,li,isb)                  8d26s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                      i10=1                                                 1d12s21
                     end do                                                 1d12s21
                     ncol=nvirt(isbvc)*nrootu                               1d12s21
                     itmp2=ibcoff                                           1d12s21
                     ibcoff=itmp2+nmul*ncol                                 1d12s21
                     call enough('hcddjkbk4. 29',bc,ibc)
                     nx=nhere2*nfdatk(2,li,isb)*nrootu                  8d26s21
                     do i=itmp2,ibcoff-1                                    1d12s21
                      bc(i)=0d0                                             1d12s21
                     end do                                                 1d12s21
                     if(ibr.eq.1)then                                       1d13s21
                      do i2=i2s,i2e                                          1d12s21
                       i2n=i2-i2s
                       iv2=i2-1                                              1d12s21
                       ntop=iv2-1                                            1d12s21
                       ntop=ntop+isw*(nvirt(isbvc)-1-ntop)                   1d12s21
                       do i=0,nfdatk(2,li,isb)-1                        8d26s21
                        do ir=0,nrootm                                       1d12s21
                         iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                         jtmp2=itmp2+i2n+nhere2*(i+nfdatk(2,li,isb)*ir) 8d26s21
                         do iv1=0,ntop                                       1d12s21
                          irec=iv1+nvirt(isbv1)*iv2                          1d13s21
                          itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                          itri=itri+isw*(irec-itri)                          1d13s21
                          bc(jtmp2+iv1*nx)=vd(iuse+itri)                     1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     else                                                   1d13s21
                      do i2=i2s,i2e                                          1d12s21
                       i2n=i2-i2s
                       iv1=i2-1                                              1d12s21
                       nbot=iv1+1                                            1d12s21
                       nbot=nbot+isw*(0-nbot)                               1d13s21
                       do i=0,nfdatk(2,li,isb)-1                        8d26s21
                        do ir=0,nrootm                                       1d12s21
                         iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                         jtmp2=itmp2+i2n+nhere2*(i+nfdatk(2,li,isb)*ir) 8d26s21
                         do iv2=nbot,nvirt(isbv2)-1                         1d13s21
                          irec=iv1+nvirt(isbv1)*iv2                          1d13s21
                          itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                          itri=itri+isw*(irec-itri)                          1d13s21
                          bc(jtmp2+iv2*nx)=vd(iuse+itri)                     1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     end if                                                 1d13s21
                     if(li.eq.1.and.isbv12.eq.1)then                     1d15s21
                      iioff=iioffp-nvirt(isbv1)*nrootu*nfdatk(2,1,isb)  8d26s21
                      do i2=i2s,i2e                                      1d15s21
                       i2n=i2-i2s                                        1d15s21
                       iv1=i2-1                                          1d15s21
                       do i=0,nfdatk(2,1,isb)-1                         8d26s21
                        do ir=0,nrootm                                   1d15s21
                         iuse=iioff+iv1+nvirt(isbv1)*(ir+nrootu*i)       1d15s21
                         jtmp2=itmp2+i2n+nhere2*(i+nfdatk(2,1,isb)*(ir  8d26s21
     $                       +nrootu*iv1))                              1d15s21
                         bc(jtmp2)=vd(iuse)*srh                          1d15s21
                        end do                                           1d15s21
                       end do                                            1d15s21
                      end do                                             1d15s21
                     end if                                              1d15s21
                     iprod=ibcoff                                           1d12s21
                     ibcoff=iprod+nrow*ncol                                 1d12s21
                     call enough('hcddjkbk4. 30',bc,ibc)
                     tff=tf*phs                                          1d14s21
                     call dgemm('n','n',nrow,ncol,nmul,tff,              1d14s21
     $                 bc(itmp1),nrow,bc(itmp2),nmul,0d0,                1d12s21
     $                bc(iprod),nrow,                                   1d12s21
     d' hcddjkbk4. 26')
                     if(ibl.eq.1)then                                       1d13s21
                      do iv2=0,nvirt(isbvc)-1                                1d12s21
                       ntop=iv2-1                                            1d12s21
                       ntop=ntop+jsw*(nvirt(jsbv1)-1-ntop)                   1d13s21
                       do ir=0,nrootm                                        1d12s21
                        do iv1=0,ntop                                        1d12s21
                         irec=iv1+nvirt(jsbv1)*iv2                           1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                         itri=itri+jsw*(irec-itri)                          1d13s21
                         iad=jjoffp+itri+nvv*ir                            1d13s21
                         jprod=iprod+nfdatb(2,lj,jsb)*(iv1              8d26s21
     $                        +nvirt(jsbv1)*(ir+nrootu*iv2))            8d26s21
                         do j=0,nfdatb(2,lj,jsb)-1                      8d26s21
                          gd(iad+j*nnvv)=gd(iad+j*nnvv)+bc(jprod+j)          1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     else                                                   1d13s21
                      do iv1=0,nvirt(isbvc)-1                                1d12s21
                       nbot=iv1+1                                           1d13s21
                       nbot=nbot+jsw*(0-nbot)                               1d13s21
                       do ir=0,nrootm                                        1d12s21
                        do iv2=nbot,nvirt(jsbv2)-1                          1d13s21
                         irec=iv1+nvirt(jsbv1)*iv2                           1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                         itri=itri+jsw*(irec-itri)                          1d13s21
                         iad=jjoffp+itri+nvv*ir                            1d13s21
                         jprod=iprod+nfdatb(2,lj,jsb)*(iv2              8d26s21
     $                        +nvirt(jsbv2)*(ir+nrootu*iv1))            8d26s21
                         do j=0,nfdatb(2,lj,jsb)-1                      8d26s21
                          gd(iad+j*nnvv)=gd(iad+j*nnvv)+bc(jprod+j)          1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     end if                                                 1d13s21
                     if(lj.eq.1.and.jsbv12.eq.1)then                     1d14s21
                      nv=nvirt(jsbv1)*nrootu                             1d14s21
                      do iv2=0,nvirt(isbvc)-1                            1d14s21
                       do ir=0,nrootm                                    1d14s21
                        iad=joff+iv2+nvirt(jsbv1)*ir                     1d14s21
                        jprod=iprod+nfdatb(2,lj,jsb)*(iv2+nvirt(jsbv1)  8d26s21
     $                      *(ir+nrootu*iv2))                           1d14s21
                        do j=0,nfdatb(2,lj,jsb)-1                       8d26s21
                         gd(iad+j*nv)=gd(iad+j*nv)+bc(jprod+j)*srh       1d14s21
                        end do                                           1d14s21
                       end do                                            1d14s21
                      end do                                             1d14s21
                     end if                                              1d14s21
                    end if                                                   1d12s21
                    ibcoff=intden                                            1d12s21
                   end if                                                1d14s21
                   jjoffp=jjoffp+nnvv*nfdatb(2,lj,jsb)                  8d26s21
                  end if                                                  1d14s21
                  if(jtf.ne.0)tfj=-1d0                                   1d15s21
                 end do                                                  1d14s21
                end if                                                   1d14s21
                iioffp=iioffp+mmvv*nfdatk(2,li,isb)                     8d26s21
                if(kase.ne.0)then                                       1d13s23
                 if(itf.ne.0)tfi=-1d0                                     1d14s21
                end if                                                  1d13s23
               end do                                                    1d14s21
              end do                                                     1d13s21
             end if                                                      1d13s21
             if(jsbv12.eq.1)joff=joff                                    1d12s21
     $           +nvirt(jsbv1)*nrootu*nfdatb(2,1,jsb)                   8d26s21
             do l=1,4                                                    1d12s21
              joff=joff+nnvv*nfdatb(2,l,jsb)                            8d26s21
             end do                                                      1d12s21
            end if                                                       1d12s21
           end do                                                        1d12s21
           if(isb.eq.jsb.and.l2e.eq.0)then                              2d23s22
            do lb=1,4                                                    12d3s21
             idv4l(1,lb)=idv4l(1,lb)+mmvv*nfdatb(2,lb,jsb)*ntype         12d3s21
             idv4yl(lb)=idv4yl(lb)                                      2d8s22
     $            +nvirt(isbv1)*nrootu*nfdatb(2,lb,jsb)*2               2d8s22
            end do                                                       12d3s21
           end if                                                       12d3s21
           if(isbv12.eq.1)then
            nvnr=nvirt(isbv1)*nrootu                                    12d3s21
            ioff=ioff+nvnr*nfdatk(2,1,isb)                              12d3s21
            if(isb.eq.jsb.and.l2e.eq.0)then                             2d23s22
             do lb=1,4                                                    12d3s21
              iorig=idv4l(2,lb)
              idv4l(2,lb)=idv4l(2,lb)+nvnr*nfdatb(2,lb,jsb)*ntype        12d3s21
             end do                                                       12d3s21
            end if                                                      12d3s21
           end if                                                       12d3s21
           do l=1,4                                                      1d12s21
            ioff=ioff+mmvv*nfdatk(2,l,isb)                              8d26s21
           end do                                                        1d12s21
          end if                                                         1d12s21
         end do                                                          1d12s21
c     we are here if isbv1=jsbv1 and isbv2=jsbv2. This occurs when
c     isbv12=jsbv12=isb*isymket=jsb*isymbra.
c     ijsb=isb*jsb
         lprt=nlprt.ne.0
         if(lprt)write(6,*)('test for den1 '),isbv12,jsbv12
         if(isbv12.eq.jsbv12)then                                       8d11s21
          if(lprt)write(6,*)('den1 block '),loop
          ioff=ioffdnon                                                  1d11s21
          joff=joffdnon                                                 8d11s21
          do isbv1=1,nsymb                                               1d11s21
           isbv2=multh(isbv1,isbv12)                                     1d11s21
           if(isbv2.ge.isbv1)then                                        1d11s21
            nvv=0                                                        4d28s21
            if(isbv1.eq.isbv2)then                                       1d11s21
             ioffs=ioff                                                  1d12s21
             joffs=joff                                                 8d11s21
             nnn=nvirt(isbv1)*nrootu                                     1d11s21
             call ilimts(nvirt(isbv1),nrootu,mynprocg,mynowprog,il,ih,   1d11s21
     $           i1s,i1e,i2s,i2e)                                       1d11s21
             nhere=ih+1-il                                               1d11s21
             if(min(nfdatb(2,1,jsb),nfdatk(2,1,isb),nhere).gt.0)then    8d26s21
              itmp=ibcoff                                                1d11s21
              ibcoff=itmp+nhere*nfdatb(2,1,jsb)                         8d26s21
              call enough('hcddjkbk4. 31',bc,ibc)
              call dgemm('n','n',nhere,nfdatb(2,1,jsb),nfdatk(2,1,isb), 8d26s21
     $           1d0,vd(ioff+il-1),nnn,bc(iden1en(1,1)),nfdatk(2,1,isb),11d23s21
     $            0d0,bc(itmp),nhere,                                   1d11s21
     d' hcddjkbk4. 27')
              do i=0,nfdatb(2,1,jsb)-1                                   8d11s21
               i10=i1s                                                   1d11s21
               i1n=nvirt(isbv1)                                          1d11s21
               jtmp=itmp+nhere*i                                         1d11s21
               do i2=i2s,i2e                                             1d11s21
                if(i2.eq.i2e)i1n=i1e                                     1d11s21
                iad=joff-1+nvirt(isbv1)*(i2-1+nrootu*i)                  1d11s21
                do i1=i10,i1n                                            1d11s21
                 gd(iad+i1)=gd(iad+i1)+bc(jtmp)                          1d11s21
                 jtmp=jtmp+1                                             1d11s21
                end do                                                   1d11s21
                i10=1                                                    1d11s21
               end do                                                    1d11s21
              end do                                                     1d11s21
              ibcoff=itmp                                                1d11s21
             end if                                                      1d11s21
             ioff=ioff+nnn*nfdatk(2,1,isb)                              8d26s21
             joff=joff+nnn*nfdatb(2,1,jsb)                              8d26s21
             nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       1d11s21
            else                                                         1d11s21
             nvv=nvirt(isbv1)*nvirt(isbv2)                               1d11s21
            end if                                                       1d11s21
            nnn=nvv*nrootu                                               1d11s21
            call ilimts(nvv,nrootu,mynprocg,mynowprog,il,ih,            1d11s21
     $            i1s,i1e,i2s,i2e)                                       1d11s21
            nhere=ih+1-il                                               1d11s21
            if(nhere.gt.0)then                                          1d11s21
             jjoff=joff                                                 11d23s21
             do lb=1,4                                                   11d23s21
              if(nfdatb(2,lb,jsb).gt.0)then                             11d23s21
               itmp=ibcoff                                                1d11s21
               ibcoff=itmp+nhere*nfdatb(2,lb,jsb)                       11d23s21
               call enough('hcddjkbk4. 32',bc,ibc)
               do iz=itmp,ibcoff-1                                      11d23s21
                bc(iz)=0d0                                              11d23s21
               end do                                                   11d23s21
               iioff=ioff                                               11d23s21
               do lk=1,4                                                11d23s21
                if(nfdatk(2,lk,isb).gt.0)then                           11d23s21
                 call dgemm('n','n',nhere,nfdatb(2,lb,jsb),              11d23s21
     $               nfdatk(2,lk,isb),1d0,vd(iioff+il-1),nnn,           11d23s21
     $               bc(iden1en(lk,lb)),nfdatk(2,lk,isb),1d0,bc(itmp),  11d23s21
     $               nhere,                                             11d23s21
     d' hcddjkbk4. 28')
                 iioff=iioff+nnn*nfdatk(2,lk,isb)                       11d23s21
                end if                                                  11d23s21
               end do                                                   11d23s21
               do i=0,nfdatb(2,lb,jsb)-1                                 8d26s21
                i10=i1s                                                   1d11s21
                i1n=nvv                                                  1d11s21
                jtmp=itmp+nhere*i                                         1d11s21
                do i2=i2s,i2e                                             1d11s21
                 if(i2.eq.i2e)i1n=i1e                                     1d11s21
                 iad=jjoff-1+nvv*(i2-1+nrootu*i)                        11d23s21
                 do i1=i10,i1n                                            1d11s21
                  gd(iad+i1)=gd(iad+i1)+bc(jtmp)                          1d11s21
                  jtmp=jtmp+1                                             1d11s21
                 end do                                                   1d11s21
                 i10=1                                                    1d11s21
                end do                                                    1d11s21
               end do                                                     1d11s21
               ibcoff=itmp                                               1d12s21
               jjoff=jjoff+nnn*nfdatb(2,lb,jsb)                         11d23s21
              end if                                                    11d23s21
             end do                                                      11d23s21
            end if                                                      11d23s21
            do l=1,4                                                     1d11s21
             ioff=ioff+nnn*nfdatk(2,l,isb)                              8d26s21
             joff=joff+nnn*nfdatb(2,l,jsb)                              8d26s21
            end do                                                       1d11s21
           end if                                                        1d11s21
          end do                                                         1d11s21
         end if                                                          1d11s21
        end if                                                          5d10s21
        ibcoff=ibc0                                                     1d8s21
        do jsbv1=1,nsymb                                                1d8s21
         jsbv2=multh(jsbv1,jsbv12)                                      1d8s21
         if(jsbv2.ge.jsbv1)then                                         1d8s21
          if(jsbv12.eq.1)then                                           1d8s21
           joffdnon=joffdnon+nrootu*nvirt(jsbv1)*nfdatb(2,1,jsb)        8d26s21
           nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                        1d8s21
          else                                                          1d8s21
           nvv=nvirt(jsbv1)*nvirt(jsbv2)                                1d8s21
          end if                                                        1d8s21
          nvv=nvv*nrootu                                                1d8s21
          do l=1,4                                                      1d8s21
           joffdnon=joffdnon+nvv*nfdatb(2,l,jsb)                        8d26s21
          end do                                                        1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
       end do                                                           1d8s21
       do isbv1=1,nsymb                                                 1d8s21
        isbv2=multh(isbv1,isbv12)                                       1d8s21
        if(isbv2.ge.isbv1)then                                          1d8s21
         if(isbv12.eq.1)then                                            1d8s21
          ioffdnon=ioffdnon+nvirt(isbv1)*nrootu*nfdatk(2,1,isb)         8d26s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         1d8s21
         else                                                           1d8s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 1d8s21
         end if                                                         1d8s21
         nvv=nvv*nrootu                                                 1d8s21
         do l=1,4                                                       1d8s21
          ioffdnon=ioffdnon+nvv*nfdatk(2,l,isb)                         8d26s21
         end do                                                         1d8s21
        end if                                                          1d8s21
       end do                                                           1d8s21
      end do                                                            1d8s21
      return                                                            1d8s21
      end                                                               1d8s21
