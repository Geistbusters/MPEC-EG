c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hccsfbk4(gb,ncsftb,nff0b,mff0b,mdoobp,i2sb,i2smb,      8d31s21
     $     veck,ncsftk,nff0k,mff0k,mdookp,i2sk,i2smk,                   8d31s21
     $     nrootz,mdon,nec,ism,irel,irefo,norb,multh,nsymb,ih0,         8d31s21
     $     i4o,isymc,idorbb,isorbb,idorbk,isorbk,icode,imap,irorip,     9d16s21
     $     irw0,irw1,irw2,ih0n,nh0,i4or,iifmx,ntypeso,isopt,l2e,bc,ibc) 11d9s22
      implicit real*8 (a-h,o-z)                                         7d11s19
c
      external second                                                   8d1s19
      integer*8 ipack,itesta,itestb,i0kc,i0ko,j0bc,j0bo,ipackc,ipackm,  9d14s21
     $     ipacks,ipacksm,ipack8,gandcc,gandco,gandcb                   2d7s23
      integer*4 ipack48(2)                                              9d20s21
      integer*2 ipack2(4),ipack2m(4),ipack2s(4),ipack2sm(4)             9d14s21
      integer*1 nab1(2),nab2(2),idorbb(*),isorbb(*),idorbk(*),isorbk(*),8d30s21
     $     icode(*),imap(*),ipackc1(8),icode1(64),icode2(64),isorb(64), 9d6s21
     $     idorb(64),nabx(2),imap1(64),imap2(64),icode1a(64),           9d9s21
     $     icode1b(64),icode2a(64),icode2b(64),imap1a(64),imap1b(64),   9d9s21
     $     imap2a(64),imap2b(64)                                        9d9s21
      equivalence (ipack,ipack2),(ipackc,ipackc1),(ipackm,ipack2m),     9d14s21
     $     (ipacks,ipack2s),(ipacksm,ipack2sm),(ipack8,ipack48)         9d20s21
      logical ldebug,lpr,lkeep                                                    1d3s20
      dimension gb(ncsftb,*),veck(ncsftk,nrootz),mff0k(*),mff0b(*),     5d12s21
     $     nff0b(mdoobp,3),nff0k(mdookp,3),nother(2),ism(*),            8d31s21
     $     irel(*),irefo(*),nab4(2,3),multh(8,8),itest(64,2),ixmtf(8),  8d30s21
     $     ih0(*),i4o(8,8,*),isy(4),igy(4),irefo2(8),ips(2),iwpb(4),    9d16s21
     $     iwpk(4),iwpb1(4),iwpb2(4),iwpk1(4),iwpk2(4),ih0n(*),nh0(*),  10d13s21
     $     i4or(8,8,*),iifmx(32),iqptr(32),isopt(4),igya(4),imy(4),     10d12s21
     $     mcsf(2),ioxx(2),nab4xx(2,2)                                              2d7s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      data ips/1,-1/                                                    9d14s21
      data icall/0/
      include "common.store"                                            7d11s19
      include "common.print"                                            1d3s20
      common/n0virtcm/n0virt                                            10d25s21
      save
      if(nsymb.eq.1)then
       igoalr=3
       igoal=1
      else
       igoalr=1
       igoal=1
      end if
      loopit=0                                                          2d24s22
      ibcoffo=ibcoff                                                    10d22s21
      nh0same=0                                                         10d12s21
      nh0sd=0
      nh0ss=0
      nh0dd=0
      nh0ds=0
      nh0dif=0                                                          10d12s21
      nh0sign=0                                                         10d12s21
      do i=1,nsymb                                                      9d1s21
       irefo2(i)=irefo(i)*2                                             8d31s21
      end do                                                            8d31s21
      ldebug=.false.                                                    5d12s21
       icall=icall+1
      lpr=.false.
      if(icall.eq.-17)ldebug=.true.
      irori=irorip-1                                                    9d6s21
      if(isopt(3).ne.0)then
       ieoro=1-irori                                                      9d6s21
      else                                                              10d14s21
       ieoro=irori                                                      9d6s21
      end if                                                            9d6s21
      if(ldebug.or.lpr)then                                                    1d3s20
       write(6,*)('Hi, my name is hccsfbk4 '),ncsftb,ncsftk,nrootz,
     $     mdon,nec,norb
       write(6,*)('ieoro '),ieoro,irori,isopt(3),loc(irori)
       write(6,*)('call '),icall
       write(6,*)('ism: '),(ism(i),i=1,norb)
       write(6,*)('rel: '),(irel(i),i=1,norb)
       write(6,*)('n0virt= '),n0virt,mdon,mdookp
       do nclop=mdon+1,mdookp
        iarg=nclop-mdon
        iff=nff0k(nclop,2)
        ipk=nff0k(nclop,3)
        nclo=nclop-1                                                    7d22s21
        nopen=nec-nclo*2                                                10d22s21
        write(6,*)('calling wfetch '),nopen
        call wfetch(nopen,mdon,idum,i2sk,i2smk,ndetk,ncsfk,iveck,iaorbk,8d30s21
     $      idum,bc,ibc)                                                11d9s22
        write(6,*)('for nclo '),nclo,nff0k(nclop,1),nff0b(nclop,1),
     $       ipk,ncsfk
        do if=1,nff0k(nclop,1)
         iff0x=iff
         i0kc=mff0k(iff)
         if(popcnt(i0kc).ne.nclo)then
          write(6,*)('we want nclo = '),nclo
          write(6,*)('we''ve got ')
          call dcbit(i0kc,norb,'kc')
          write(6,*)('this is for if = '),if
          call dws_synca
          call dws_finalize
          stop
         end if
         iff=iff+1
         i0ko=mff0k(iff)
         iff=iff+1
         sz=0d0                                                         7d16s21
         do ir=1,nrootz                                                 7d16s21
          do i=0,ncsfk-1                                                10d22s21
           sz=sz+veck(ipk+i,ir)**2                                      7d16s21
          end do                                                        7d16s21
         end do                                                         7d16s21
         sz=sqrt(sz/dfloat(nrootz*ncsfk))                               10d22s21
         if(sz.gt.1d-10)then
          write(6,*)('ket '),nclop,if,ipk
          call dcbit(i0kc,norb,'c')
          do i=n0virt,norb
           if(btest(i0kc,i))then
            write(6,*)irel(i),('s'),ism(i)
           end if
          end do
          call dcbit(i0ko,norb,'o')
          do i=n0virt,norb
           if(btest(i0ko,i))then
            write(6,*)irel(i),('s'),ism(i)
           end if
          end do
 1493     continue
          call prntm2(veck(ipk,1),ncsfk,nrootz,ncsftk)
         end if                                                         7d16s21
         ipk=ipk+ncsfk                                                  10d22s21
        end do
       end do
       write(6,*)('at end of block')
       ldebug=.true.
      end if                                                            1d3s20
      if(lpr)then
      write(6,*)('what do we have for ih0n? '),isymc
      do lsb=1,nsymb
       lsk=multh(lsb,isymc)
       if(ih0n(lsb).gt.0)then
        write(6,*)('for lsb, lsk '),lsb,lsk
        call prntm2(bc(ih0n(lsb)),nh0(lsb),nh0(lsk),nh0(lsb))
       end if
      end do
      end if
      mdoo=mdoop-1                                                      5d10s21
      idelta=i2smb-i2smk                                                8d31s21
      do nclokp=mdon+1,mdookp                                           8d31s21
       nclok=nclokp-1                                                     5d7s21
       iargk=nclokp-mdon                                                  7d11s19
       nopenk=nec-2*nclok                                                 7d11s19
       call wfetch(nopenk,mdon,idum,i2sk,i2smk,ndetk,ncsfk,iveck,iaorbk,8d30s21
     $      idum,bc,ibc)                                                11d9s22
       do nclobp=max(mdon+1,nclokp-2),min(mdoobp,nclokp+2)              8d31s21
        nclob=nclobp-1                                                  6d11s19
        nopenb=nec-2*nclob
        jargb=nclobp-mdon                                               6d11s19
        call wfetch(nopenb,mdon,idum,i2sb,i2smb,ndetb,ncsfb,ivecb,      8d30s21
     $       iaorbb,idum,bc,ibc)                                        11d9s22
        iff0=nff0k(nclokp,2)                                               5d7s21
        ipk=nff0k(nclokp,3)                                                 5d7s21
        do if=1,nff0k(nclokp,1)                                            5d7s21
         i0kc=mff0k(iff0)                                                  5d7s21
         ii=1                                                           8d30s21
         do i=1,norb                                                    8d30s21
          if(btest(i0kc,i))then                                         8d30s21
           idorbk(ii)=i                                                 8d30s21
           ii=ii+1                                                      8d30s21
          end if                                                        8d30s21
         end do                                                         8d30s21
         iff0=iff0+1                                                     5d7s21
         i0ko=mff0k(iff0)                                                  5d7s21
         ii=1                                                           8d30s21
         do i=1,norb                                                    8d30s21
          if(btest(i0ko,i))then                                         8d30s21
           isorbk(ii)=i                                                 8d30s21
           ii=ii+1                                                      8d30s21
          end if                                                        8d30s21
         end do                                                         8d30s21
         iff0=iff0+1                                                     5d7s21
         if(mod(loopit,mynprocg).eq.mynowprog)then                      3d13s23
          loopit=mynowprog                                              3d13s23
          jff0=nff0b(nclobp,2)                                               5d7s21
          jpb=nff0b(nclobp,3)                                                 5d7s21
          do jf=1,nff0b(nclobp,1)                                         5d7s21
           if(lpr)write(6,*)('for '),jf,nclob,if,nclok,loop
           j0bc=mff0b(jff0)                                                  5d7s21
           ii=1                                                           8d30s21
           do i=1,norb                                                    8d30s21
            if(btest(j0bc,i))then                                         8d30s21
             idorbb(ii)=i                                                 8d30s21
             ii=ii+1                                                      8d30s21
            end if                                                        8d30s21
           end do                                                         8d30s21
           jff0=jff0+1                                                     5d7s21
           j0bo=mff0b(jff0)                                                  5d7s21
           if(lpr)then
            call dcbit(i0kc,norb,'kc')
            call dcbit(i0ko,norb,'ko')
            call dcbit(j0bc,norb,'bc')
            call dcbit(j0bo,norb,'bo')
           end if
           ii=1                                                           8d30s21
           do i=1,norb                                                    8d30s21
            if(btest(j0bo,i))then                                         8d30s21
             isorbb(ii)=i                                                 8d30s21
             ii=ii+1                                                      8d30s21
            end if                                                        8d30s21
           end do                                                         8d30s21
           jff0=jff0+1                                                     5d7s21
           gandcc=ieor(j0bc,i0kc)                                        2d6s23
           gandco=ieor(j0bo,i0ko)                                        2d6s23
           gandcb=ior(gandcc,gandco)                                     2d6s23
           ndifb=popcnt(gandcb)                                          2d6s23
           idelta=i2smk-i2smb                                            9d13s21
           idelta=iabs(idelta)                                           9d13s21
           if(lpr)write(6,*)('ndifb '),ndifb,idelta,loopit,
     $         mynprocg
           if(ndifb.le.4)then                                            2d7s23
            ndifs=popcnt(gandco)                                        2d6s23
            ndifd=popcnt(gandcc)                                        2d6s23
            if(lpr)write(6,*)('ndifs,ndifd '),ndifs,ndifd
            if(ndifs.eq.0.and.ndifd.eq.0)then                           2d7s23
             idelta=i2smk-i2smb                                            9d13s21
             idelta=iabs(idelta)                                           9d13s21
             if(idelta.eq.0)then                                        9d13s21
              if(i2sb.eq.i2sk)then                                      9d13s21
               sum=0d0                                                  9d13s21
               do i=1,norb                                               9d1s21
                if(btest(j0bc,i))then                                    9d1s21
                 is=ism(i)                                               9d1s21
                 ig=irel(i)-1                                            9d1s21
                 ig0=ig                                                 10d8s21
                 do ipass=1,2                                            9d1s21
                  h0=-2d0                                               10d14s21
                  if(ih0n(is).gt.0)then                                 10d14s21
                   iadn=ih0n(is)+(nh0(is)+1)*ig0                        10d14s21
                   h0=bc(iadn)                                           10d14s21
                   if(ipass.eq.2.and.isopt(4).ne.0)h0=-h0               2d24s22
                   sum=sum+h0                                           10d14s21
                  end if                                                10d14s21
                  isp=ig/irefo(is)                                      10d8s21
                  ispm=1-isp                                            10d8s21
                  do j=1,norb                                            9d1s21
                   if(btest(j0bc,j).and.l2e.eq.0)then                   2d17s22
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
                     if(iifmx(itestp).ge.0)then                         10d8s21
                       iad=i4or(js,is,is)+jg0+irefo(js)*(jg0+irefo(js)   10d8s21
     $                  *(ig0+irefo(is)*(ig0+irefo(is)*iifmx(itestp)))) 10d8s21
                       sum=sum+bc(iad)*0.5d0                             10d8s21
                     else if(iifmx(jtestp).ge.0)then                    10d8s21
                       iad=i4or(js,is,is)+jg0+irefo(js)*(jg0+irefo(js)   10d8s21
     $                  *(ig0+irefo(is)*(ig0+irefo(is)*iifmx(jtestp)))) 10d8s21
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
                       iad=i4or(is,is,js)+jg0+irefo(js)*(ig0+irefo(is)   10d8s21
     $                  *(ig0+irefo(is)*(jg0+irefo(js)*iifmx(itestp)))) 10d8s21
                       sum=sum-bc(iad)*0.5d0                             10d8s21
                      else if(iifmx(jtestp).ge.0)then                    10d8s21
                       iad=i4or(is,is,js)+jg0+irefo(js)*(ig0+irefo(is)   10d8s21
     $                   *(ig0+irefo(is)*(jg0+irefo(js)*iifmx(jtestp)))) 10d8s21
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
               do ir=1,nrootz                                           10d20s21
                do i=0,ncsfb-1                                           9d13s21
                 gb(jpb+i,ir)=gb(jpb+i,ir)+sum*veck(ipk+i,ir)           10d20s21
                end do
               end do                                                   9d13s21
              end if                                                    9d13s21
             end if                                                     9d13s21
             ii=0                                                       9d13s21
             do i=1,norb                                                9d13s21
              if(btest(j0bo,i))then                                     9d13s21
               ii=ii+1                                                  9d13s21
               icode1(ii)=3                                             9d13s21
               imap1(ii)=i
              end if                                                    9d13s21
             end do                                                     9d13s21
             nx1=ii                                                     9d13s21
              mxmat=norb*32                                             9d15s21
              itype=ibcoff
              imatx=itype+mxmat
              ibcoff=imatx+mxmat*ncsfb*ncsfk
              if(ibcoff.lt.0)write(6,*)('enougha')
              call enough('hccsfbk4.  1',bc,ibc)
              do iz=imatx,ibcoff-1                                      9d15s21
               bc(iz)=0d0                                               9d15s21
              end do                                                    9d15s21
              ntype=0
             if(idelta.le.2)then                                        9d13s21
              sigc=1d0                                                  9d21s21
              ntype1a=0                                                 10d14s21
              call getcup(irw0,nopenb,i2sb,i2smb,i2sk,i2smk,ntypeg,     10d18s21
     $             ioutg,imatg,mcsf,imap1,bc,ibc)                       11d14s22
              if(ntypeg.gt.0)then
               ncsfa=mcsf(1)                                            10d18s21
               ncsfc=mcsf(2)                                            10d19s21
               do iq=0,ntypeg-1
                ipackc=ibc(ioutg+iq)                                    10d19s21
                if(ipackc1(1).gt.0)then                                 10d19s21
                 mmm=ipackc1(1)                                         10d19s21
                 igy(1)=mmm
                 lsb=ism(mmm)                                           10d14s21
                 igya(1)=irel(mmm)-1                                    10d14s21
                 imy(1)=0                                               10d14s21
                else                                                    10d14s21
                 mmm=-ipackc1(1)                                        10d19s21
                 igy(1)=mmm
                 lsb=ism(mmm)                                           10d14s21
                 igya(1)=irel(mmm)-1                                    10d14s21
                 imy(1)=1                                               10d14s21
                end if                                                  10d14s21
                if(ipackc1(2).gt.0)then                                 10d19s21
                 mmm=ipackc1(2)                                         10d19s21
                 igy(2)=mmm
                 lsk=ism(mmm)                                           10d14s21
                 igya(2)=irel(mmm)-1                                    10d14s21
                 imy(2)=0                                               10d14s21
                else                                                    10d14s21
                 mmm=-ipackc1(2)                                        10d19s21
                 igy(2)=mmm
                 lsk=ism(mmm)                                           10d14s21
                 igya(2)=irel(mmm)-1                                    10d14s21
                 imy(2)=1                                               10d14s21
                end if                                                  10d14s21
                imatq=imatg+mcsf(1)*mcsf(2)*iq
                sum=0d0                                                  10d14s21
                h0=-2d0                                                   10d13s21
                if(ih0n(lsb).gt.0)then                                    10d13s21
                 iad=ih0n(lsb)+igya(1)+nh0(lsb)*igya(2)                  10d14s21
                 h0=bc(iad)                                                10d13s21
                 if(imy(1).ne.0.and.isopt(4).ne.0)h0=-h0                 10d14s21
                 sum=h0                                                   10d14s21
                else if(ih0n(lsk).gt.0)then                               10d13s21
                 iad=ih0n(lsk)+igya(2)+nh0(lsk)*igya(1)                  10d14s21
                 h0=bc(iad)                                                10d13s21
                 if(isopt(2).ne.0)h0=-h0                                  10d13s21
                 if(imy(2).ne.0.and.isopt(4).ne.0)h0=-h0                  10d13s21
                 sum=h0                                                   10d14s21
                end if                                                    10d13s21
                do i=1,norb                                             10d14s21
                 if((btest(j0bc,i).or.btest(j0bo,i)).and.l2e.eq.0)then  2d17s22
                  is=ism(i)                                               9d13s21
                  ig=irel(i)-1                                            9d13s21
                  if(btest(j0bc,i))then                                 10d14s21
                   ff=-1d0                                              10d14s21
                  else                                                  10d14s21
                   ff=-0.5d0                                            10d14s21
                  end if                                                10d14s21
                  do isp=0,1                                             10d14s21
                   ltest=isp+2*(imy(1)+2*(imy(2)+2*isp))                10d14s21
                   itestp=ltest+1                                        10d14s21
                   jtest=1-isp+2*(1-imy(1)+2*(1-imy(2)+2*(1-isp)))      10d14s21
                   jtestp=jtest+1                                        10d14s21
                    if(iifmx(itestp).ge.0)then                           10d14s21
                     iad=i4or(lsb,lsb,is)+ig+irefo(is)*(igya(1)          10d14s21
     $                  +irefo(lsb)*(igya(2)+irefo(lsb)*                10d14s21
     $                  (ig+irefo(is)*iifmx(itestp))))                  10d14s21
                     xint=bc(iad)                                        10d14s21
                     sum=sum+xint*ff                                     10d14s21
                    else if(iifmx(jtestp).ge.0)then                           10d14s21
                     iad=i4or(lsb,lsb,is)+ig+irefo(is)*(igya(1)          10d14s21
     $                  +irefo(lsb)*(igya(2)+irefo(lsb)*                10d14s21
     $                  (ig+irefo(is)*iifmx(jtestp))))                  10d14s21
                     xint=bc(iad)                                        10d14s21
                     if(isopt(4).ne.0)xint=-xint                         10d14s21
                     sum=sum+xint*ff                                     10d14s21
                    end if                                               10d14s21
                  end do                                                10d14s21
                 end if                                                 10d14s21
                end do                                                  10d14s21
                call dgemm('n','n',ncsfb,nrootz,ncsfk,sum,bc(imatq),    10d20s21
     $                ncsfb,veck(ipk,1),ncsftk,1d0,gb(jpb,1),ncsftb,    10d20s21
     d' hccsfbk4.  1')
               end do                                                   9d20s21
              end if                                                    9d20s21
             end if                                                     9d13s21
             if(idelta.le.4.and.l2e.eq.0)then                           2d17s22
c
c     intermediate state = bra or ket
c
              nx1=0                                                     9d15s21
              do i=1,norb                                               9d15s21
               if(btest(j0bo,i))then                                    9d15s21
                nx1=nx1+1                                               9d15s21
                imap1(nx1)=i                                            9d15s21
                icode1(nx1)=3                                           9d15s21
               end if                                                   9d15s21
              end do                                                    9d15s21
              call spinloop1(i2sb,i2smb,i2sk,i2smk,nopenb,ncsfb,imap1,  10d20s21
     $             nrootz,itypeg,imatg,ntypeg,irw0,veck(ipk,1),ncsftk,  10d20s21
     $             ieoro,bc,ibc)                                        11d14s22
              do i=0,ntypeg-1
               ipack=ibc(itypeg+i)
               do j=1,4                                                 10d20s21
                if(ipack2(j).gt.0)then                                  10d20s21
                 isy(j)=ism(ipack2(j))                                  10d20s21
                 igya(j)=irel(ipack2(j))-1                               10d20s21
                 igy(j)=ipack2(j)
                 imy(j)=0                                               10d20s21
                else                                                    10d20s21
                 isy(j)=ism(-ipack2(j))                                 10d20s21
                 igy(j)=-ipack2(j)                                      11d23s21
                 igya(j)=irel(-ipack2(j))-1                             10d20s21
                 imy(j)=1                                               10d20s21
                end if                                                  10d20s21
               end do                                                   10d20s21
               ltest=imy(1)+2*(imy(2)+2*(imy(3)+2*imy(4)))              10d20s21
               itestp=ltest+1                                           10d20s21
               jtest=1-imy(1)+2*(1-imy(2)+2*(1-imy(3)+2*(1-imy(4))))    10d20s21
               jtestp=jtest+1                                           10d20s21
               xint=0d0                                                 10d20s21
               ihit=0                                                   10d20s21
                if(iifmx(itestp).ge.0)then                               10d20s21
                 iad=i4or(isy(2),isy(3),isy(4))+igya(1)+irefo(isy(1))    10d20s21
     $                *(igya(2)+irefo(isy(2))*(igya(3)+irefo(isy(3))*   10d20s21
     $                (igya(4)+irefo(isy(4))*iifmx(itestp))))           10d20s21
                 xint=bc(iad)                                            10d20s21
                 ihit=1                                                  10d20s21
                else if(iifmx(jtestp).ge.0)then                          10d20s21
                 iad=i4or(isy(2),isy(3),isy(4))+igya(1)+irefo(isy(1))    10d20s21
     $                *(igya(2)+irefo(isy(2))*(igya(3)+irefo(isy(3))*   10d20s21
     $                (igya(4)+irefo(isy(4))*iifmx(jtestp))))           10d20s21
                 xint=bc(iad)                                            10d20s21
                 if(isopt(4).ne.0)xint=-xint                             10d20s21
                 ihit=1                                                  10d20s21
                end if                                                   10d20s21
               if(ihit.ne.0)then
                xint0=xint                                              11d29s21
                xint=xint*0.5d0
                jmat=imatg+ncsfb*nrootz*i
                do ir=1,nrootz                                          10d20s21
                 jtmpxq=jmat+ncsfb*(ir-1)                               10d20s21
                 do j=0,ncsfb-1                                         10d20s21
                  gb(jpb+j,ir)=gb(jpb+j,ir)+xint*bc(jtmpxq+j)           10d20s21
                 end do                                                 10d20s21
                end do                                                  10d20s21
               end if
              end do                                                    10d20s21
              ibcoff=itypeg
              do i1=1,norb                                              9d13s21
               if(btest(i0ko,i1))then                                   9d15s21
                do i2=i1+1,norb                                             9d13s21
                 if(btest(i0ko,i2))then                                 9d15s21
                  itesta=i0kc                                              5d7s21
                  itestb=i0ko                                              5d7s21
                  if(btest(itesta,i1))then                               9d13s21
                   itesta=ibclr(itesta,i1)                               9d13s21
                   itestb=ibset(itestb,i1)                               9d13s21
                   nopenkk=nopenb+1                                           1d22s21
                  else if(btest(itestb,i1))then                          9d13s21
                   itestb=ibclr(itestb,i1)                               9d13s21
                   nopenkk=nopenb-1                                              11d13s20
                  else                                                          11d13s20
                   write(6,*)('bit not set for i1 = '),i1
                   stop 'nab4(1,1)'                                          11d27s20
                  end if                                                        11d13s20
                  if(btest(itestb,i2))then                               9d13s21
                   itesta=ibset(itesta,i2)                               9d13s21
                   itestb=ibclr(itestb,i2)                               9d13s21
                   nopenkk=nopenkk-1                                              11d13s20
                  else if(btest(itesta,i2))then                         9d13s21
                   go to 1848
                  else                                                          11d13s20
                   itestb=ibset(itestb,i2)                               9d13s21
                   nopenkk=nopenkk+1                                              11d13s20
                  end if                                                        11d13s20
                  call gandcr(j0bc,j0bo,itesta,itestb,nopenb,            9d13s21
     $               nopenkk,norb,nnot1,nab1,icode1,imap1,nx1,irw1,     9d16s21
     $                 irw2,iwpb1,iwpk1,bc,ibc)                         11d14s22
                  call gandcr(itesta,itestb,i0kc,i0ko,nopenkk,nopenk,        5d7s21
     $                 norb,nnot2,nab2,icode2,imap2,nx2,irw1,irw2,      9d16s21
     $                 iwpb2,iwpk2,bc,ibc)                              11d14s22
                  call spinloop(i2sb,i2smb,i2sk,i2smk,nopenb,nopenk,    10d19s21
     $                 nopenkk,ncsfb,nrootz,itype,imatx,ntype,nab1,     10d19s21
     $                 iwpb1,iwpk1,nab2,iwpb2,iwpk2,veck(ipk,1),ncsftk, 10d19s21
     $                 ieoro,bc,ibc)                                    11d14s22
                  do i=0,ntype-1
                   ipack=ibc(itype+i)
                   do j=1,4                                               9d9s21
                    if(ipack2(j).gt.0)then                                9d9s21
                     isy(j)=ism(ipack2(j))                                9d9s21
                     igya(j)=irel(ipack2(j))-1                             9d9s21
                     igy(j)=ipack2(j)
                     imy(j)=0                                             10d13s21
                    else                                                  9d9s21
                     isy(j)=ism(-ipack2(j))                                9d9s21
                     igy(j)=-ipack2(j)
                     igya(j)=irel(-ipack2(j))-1                           10d13s21
                     imy(j)=1                                             10d13s21
                    end if                                                9d9s21
                   end do                                                 9d9s21
                   ltest=imy(1)+2*(imy(2)+2*(imy(3)+2*imy(4)))            10d13s21
                   itestp=ltest+1                                         10d13s21
                   jtest=1-imy(1)+2*(1-imy(2)+2*(1-imy(3)+2*(1-imy(4))))  10d13s21
                   jtestp=jtest+1                                         10d13s21
                   xint=0d0                                               9d10s21
                   ihit=0
                   if(iifmx(itestp).ge.0)then                             10d13s21
                    iad=i4or(isy(2),isy(3),isy(4))+igya(1)+irefo(isy(1))  10d13s21
     $                *(igya(2)+irefo(isy(2))*(igya(3)+irefo(isy(3))*   10d13s21
     $                (igya(4)+irefo(isy(4))*iifmx(itestp))))           10d13s21
                    xint=bc(iad)                                         10d13s21
                    ihit=1                                                10d13s21
                   else if(iifmx(jtestp).ge.0)then                        10d13s21
                    iad=i4or(isy(2),isy(3),isy(4))+igya(1)+irefo(isy(1))  10d13s21
     $                *(igya(2)+irefo(isy(2))*(igya(3)+irefo(isy(3))*   10d13s21
     $                (igya(4)+irefo(isy(4))*iifmx(jtestp))))           10d13s21
                    xint=bc(iad)                                         10d13s21
                    if(isopt(4).ne.0)xint=-xint                         10d19s21
                    ihit=1                                                10d13s21
                   end if                                                 10d13s21
                   if(ihit.ne.0)then                                      9d10s21
                    itmpxq=imatx+ncsfb*nrootz*i                           10d19s21
                    do ir=1,nrootz                                        10d19s21
                     jtmpxq=itmpxq+ncsfb*(ir-1)                           10d19s21
                     do j=0,ncsfb-1                                       10d19s21
                      gb(jpb+j,ir)=gb(jpb+j,ir)+xint*bc(jtmpxq+j)         10d19s21
                     end do                                               10d19s21
                    end do                                                10d19s21
                   end if                                                 9d9s21
                  end do
                  ibcoff=itype                                            9d6s21
 1848             continue
                 end if
                end do                                                   9d13s21
               end if                                                   9d13s21
              end do                                                    9d13s21
             end if                                                     9d13s21
            else if(ndifs.eq.2.and.ndifb.eq.2)then                      2d6s23
             do i=1,norb                                                2d6s23
              if(btest(gandco,i))then                                   2d6s23
               if((btest(j0bo,i).and..not.btest(i0kc,i)).or.            2d6s23
     $             (btest(j0bc,i).and.btest(i0ko,i)))then               2d6s23
                nab4(1,1)=i                                             2d6s23
               else                                                     2d6s23
                nab4(2,1)=i                                             2d6s23
               end if                                                   2d6s23
              end if                                                    2d6s23
             end do                                                     2d6s23
             call gandcr(j0bc,j0bo,i0kc,i0ko,nopenb,                     9d10s21
     $          nopenk,norb,nnot1,nab1,icode1,imap1,nx1,irw1,irw2,      9d16s21
     $            iwpb,iwpk,bc,ibc)                                     11d14s22
             call gencup(i2sb,i2smb,i2sk,i2smk,nopenb,nopenk,nab1,iwpb, 9d17s21
     $            iwpk,ioutg,imatg,ntypeg,mcsf,veck(ipk,1),ncsftk,      10d19s21
     $            nrootz,bc,ibc)                                        11d14s22
             ntype1a=ntypeg                                             10d12s21
             iouta=ioutg                                                10d12s21
             imata=imatg                                                10d12s21
             do ii=0,ntype1a-1                                          9d10s21
              ipackc=ibc(iouta+ii)                                      9d10s21
              if(ipackc1(1).gt.0)then                                   9d10s21
               lsb=ism(ipackc1(1))                                      10d12s21
               lgb=irel(ipackc1(1))-1                                   10d12s21
               igya(3)=lgb                                              10d12s21
               igy(3)=ipackc1(1)
               imy(3)=0                                                 10d12s21
              else
               lsb=ism(-ipackc1(1))                                     10d12s21
               igya(3)=irel(-ipackc1(1))-1                              10d12s21
               igy(3)=-ipackc1(1)
               lgb=irel(-ipackc1(1))-1+irefo(lsb)                       10d12s21
               imy(3)=1                                                 10d12s21
              end if
              if(ipackc1(2).gt.0)then                                   9d10s21
               lsk=ism(ipackc1(2))                                      10d12s21
               lgk=irel(ipackc1(2))-1                                   10d12s21
               igy(4)=ipackc1(2)
               igya(4)=lgk                                              10d12s21
               imy(4)=0                                                 10d12s21
              else
               lsk=ism(-ipackc1(2))                                     10d12s21
               igya(4)=irel(-ipackc1(2))-1                              10d12s21
               lgk=irel(-ipackc1(2))-1+irefo(lsk)                       10d12s21
               igy(4)=-ipackc1(2)
               imy(4)=1                                                 10d12s21
              end if
              h0=-2d0                                                   10d13s21
              sum=0d0                                                   10d14s21
              if(ih0n(lsb).gt.0)then                                    10d13s21
               iad=ih0n(lsb)+igya(3)+nh0(lsb)*igya(4)                   10d14s21
               h0=bc(iad)                                                10d13s21
               if(imy(3).ne.0.and.isopt(4).ne.0)h0=-h0                   10d13s21
               sum=h0                                                   10d14s21
              else if(ih0n(lsk).gt.0)then                               10d13s21
               iad=ih0n(lsk)+igya(4)+nh0(lsk)*igya(3)                   10d14s21
               h0=bc(iad)                                                10d13s21
               if(isopt(2).ne.0)h0=-h0                                  10d13s21
               if(imy(4).ne.0.and.isopt(4).ne.0)h0=-h0                  10d13s21
               sum=h0                                                   10d14s21
              end if                                                    10d13s21
              if(lpr.and.abs(sum).gt.1d-10)write(6,*)('h0 '),h0
              do i=1,norb                                               9d10s21
               if(btest(j0bc,i).and.btest(i0kc,i).and.l2e.eq.0)then     2d17s22
                js=ism(i)                                               9d10s21
                jg=irel(i)-1                                            9d10s21
                jg0=jg                                                  10d12s21
                do jpass=1,2                                            9d10s21
                 ltest=jpass-1+2*(jpass-1+2*(imy(3)+2*imy(4)))          10d12s21
                 itestp=ltest+1                                         10d12s21
                 jtest=2-jpass+2*(2-jpass+2*(1-imy(3)+2*(1-imy(4))))    10d12s21
                 jtestp=jtest+1                                         10d12s21
                 xintj=-2d0                                             10d12s21
                 if(iifmx(itestp).ge.0)then                             10d12s21
                  iad=i4or(js,lsb,lsk)+jg0+irefo(js)*(jg0+irefo(js)*    10d12s21
     $                (igya(3)+irefo(lsb)*(igya(4)+irefo(lsk)*          10d12s21
     $                iifmx(itestp))))                                  10d12s21
                  xintj=bc(iad)                                         10d12s21
                  sum=sum+xintj                                         10d12s21
                 else if(iifmx(jtestp).ge.0)then                        10d12s21
                  iad=i4or(js,lsb,lsk)+jg0+irefo(js)*(jg0+irefo(js)*    10d12s21
     $                (igya(3)+irefo(lsb)*(igya(4)+irefo(lsk)*          10d12s21
     $                iifmx(jtestp))))                                  10d12s21
                  xintj=bc(iad)                                         10d12s21
                  if(isopt(4).ne.0)xintj=-xintj                         10d12s21
                  sum=sum+xintj                                         10d12s21
                 end if                                                 10d12s21
                 ltest=jpass-1+2*(imy(4)+2*(imy(3)+2*(jpass-1)))        10d12s21
                 itestp=ltest+1                                         10d12s21
                 jtest=2-jpass+2*(1-imy(4)+2*(1-imy(3)+2*(2-jpass)))    10d12s21
                 jtestp=jtest+1                                         10d12s21
                 xintj=-2d0                                             10d12s21
                 if(iifmx(itestp).ge.0)then                             10d12s21
                  iad=i4or(lsk,lsb,js)+jg0+irefo(js)*(igya(4)           10d12s21
     $                +irefo(lsk)*(igya(3)+irefo(lsb)*(jg0+irefo(js)*   10d12s21
     $                iifmx(itestp))))                                  10d12s21
                  xintk=bc(iad)                                         10d12s21
                  sum=sum-xintk                                         10d12s21
                 else if(iifmx(jtestp).ge.0)then                        10d12s21
                  iad=i4or(lsk,lsb,js)+jg0+irefo(js)*(igya(4)           10d12s21
     $                 +irefo(lsk)*(igya(3)+irefo(lsb)*(jg0+irefo(js)*  10d12s21
     $                iifmx(jtestp))))                                  10d12s21
                  xintk=bc(iad)                                         10d12s21
                  if(isopt(4).ne.0)xintk=-xintk                         10d12s21
                  sum=sum-xintk                                         10d12s21
                 end if                                                 10d12s21
                 jg=jg+irefo(js)                                        9d10s21
                end do                                                  9d10s21
               end if                                                   9d10s21
              end do                                                    9d10s21
              if(ipackc1(1).gt.0)then                                   10d12s21
               in=ipackc1(1)                                            10d12s21
               lgbm=lgb+irefo(lsb)                                      9d13s21
               imy(1)=1                                                 10d12s21
              else                                                      9d10s21
               in=-ipackc1(1)                                           10d12s21
               lgbm=lgb-irefo(lsb)                                      9d13s21
               imy(1)=0                                                 10d12s21
              end if                                                    9d13s21
              if(btest(j0bc,in).and.l2e.eq.0)then                       2d17s22
               ltest=imy(1)+2*(imy(1)+2*(imy(3)+2*imy(4)))              10d12s21
               itestp=ltest+1                                           10d12s21
               jtest=1-imy(1)+2*(1-imy(1)+2*(1-imy(3)+2*(1-imy(4))))    10d12s21
               jtestp=jtest+1                                           10d12s21
               xintj=-2d0                                               10d12s21
               if(iifmx(itestp).ge.0)then                               10d12s21
                iad=i4or(lsb,lsb,lsk)+igya(3)+irefo(lsb)*(igya(3)       10d12s21
     $              +irefo(lsb)*(igya(3)+irefo(lsb)*(igya(4)+irefo(lsk) 10d12s21
     $              *iifmx(itestp))))                                   10d12s21
                xintj=bc(iad)                                           10d12s21
                sum=sum+xintj                                           10d13s21
               else if(iifmx(jtestp).ge.0)then                          10d12s21
                iad=i4or(lsb,lsb,lsk)+igya(3)+irefo(lsb)*(igya(3)       10d12s21
     $               +irefo(lsb)*(igya(3)+irefo(lsb)*(igya(4)           10d12s21
     $               +irefo(lsk)*iifmx(jtestp))))                       10d12s21
                xintj=bc(iad)                                           10d12s21
                if(isopt(4).ne.0)xintj=-xintj                           10d12s21
                sum=sum+xintj                                           10d13s21
               end if                                                   10d12s21
               ltest=imy(1)+2*(imy(4)+2*(imy(3)+2*imy(1)))              10d12s21
               itestp=ltest+1                                           10d12s21
               jtest=1-imy(1)+2*(1-imy(4)+2*(1-imy(3)+2*(1-imy(1))))    10d12s21
               jtestp=jtest+1                                           10d12s21
               xintk=-2d0                                               10d12s21
               if(iifmx(itestp).ge.0)then                               10d12s21
                iad=i4or(lsk,lsb,lsb)+igya(3)+irefo(lsb)*(igya(4)       10d12s21
     $              +irefo(lsk)*(igya(3)+irefo(lsb)*(igya(3)+irefo(lsb) 10d12s21
     $              *iifmx(itestp))))                                   10d12s21
                xintk=bc(iad)                                           10d12s21
                sum=sum-xintk                                           10d13s21
               else if(iifmx(jtestp).ge.0)then                          10d12s21
                iad=i4or(lsk,lsb,lsb)+igya(3)+irefo(lsb)*(igya(4)       10d12s21
     $               +irefo(lsk)*(igya(3)+irefo(lsb)*(igya(3)           10d12s21
     $               +irefo(lsb)*iifmx(jtestp))))                       10d12s21
                xintk=bc(iad)                                           10d12s21
                if(isopt(4).ne.0)xintk=-xintk                           10d12s21
                sum=sum-xintk                                           10d13s21
               end if                                                   10d12s21
              end if                                                    9d13s21
              if(ipackc1(2).gt.0)then                                   10d12s21
               im=ipackc1(2)                                            10d12s21
               lgkm=lgk+irefo(lsk)                                      10d12s21
               imy(2)=1                                                 10d13s21
              else                                                      10d12s21
               im=-ipackc1(2)                                           10d12s21
               lgkm=lgk-irefo(lsk)                                      10d12s21
               imy(2)=0                                                 10d13s21
              end if                                                    10d12s21
              if(btest(i0kc,im).and.l2e.eq.0)then                       2d17s22
               ltest=imy(2)+2*(imy(2)+2*(imy(3)+2*imy(4)))              10d12s21
               itestp=ltest+1                                           10d12s21
               jtest=1-imy(2)+2*(1-imy(2)+2*(1-imy(3)+2*(1-imy(4))))    10d12s21
               jtestp=jtest+1                                           10d12s21
               xintj=-2d0                                               10d12s21
               if(iifmx(itestp).ge.0)then                               10d12s21
                iad=i4or(lsk,lsb,lsk)+igya(4)+irefo(lsk)*(igya(4)       10d12s21
     $              +irefo(lsk)*(igya(3)+irefo(lsb)*(igya(4)+irefo(lsk) 10d12s21
     $              *iifmx(itestp))))                                   10d12s21
                xintj=bc(iad)                                           10d12s21
                sum=sum+xintj                                           10d13s21
               else if(iifmx(jtestp).ge.0)then                          10d12s21
                iad=i4or(lsk,lsb,lsk)+igya(4)+irefo(lsk)*(igya(4)       10d12s21
     $               +irefo(lsk)*(igya(3)+irefo(lsb)*(igya(4)           10d12s21
     $               +irefo(lsk)*iifmx(jtestp))))                       10d12s21
                xintj=bc(iad)                                           10d12s21
                if(isopt(4).ne.0)xintj=-xintj                           10d12s21
                sum=sum+xintj                                           10d13s21
               end if                                                   10d12s21
               ltest=imy(2)+2*(imy(4)+2*(imy(3)+2*imy(2)))              10d12s21
               itestp=ltest+1                                           10d12s21
               jtest=1-imy(2)+2*(1-imy(4)+2*(1-imy(3)+2*(1-imy(2))))    10d12s21
               jtestp=jtest+1                                           10d12s21
               xintk=-2d0                                               10d12s21
               if(iifmx(itestp).ge.0)then                               10d12s21
                iad=i4or(lsk,lsb,lsk)+igya(4)+irefo(lsk)*(igya(4)       10d12s21
     $              +irefo(lsk)*(igya(3)+irefo(lsb)*(igya(4)+irefo(lsk) 10d12s21
     $              *iifmx(itestp))))                                   10d12s21
                xintk=bc(iad)                                           10d12s21
                sum=sum-xintk                                           10d13s21
               else if(iifmx(jtestp).ge.0)then                          10d12s21
                iad=i4or(lsk,lsb,lsk)+igya(4)+irefo(lsk)*(igya(4)       10d12s21
     $               +irefo(lsk)*(igya(3)+irefo(lsb)*(igya(4)           10d12s21
     $               +irefo(lsk)*iifmx(jtestp))))                       10d12s21
                xintk=bc(iad)                                           10d12s21
                if(isopt(4).ne.0)xintk=-xintk                           10d12s21
                sum=sum-xintk                                           10d13s21
               end if                                                   10d12s21
              end if                                                    9d13s21
              imat=imata+ncsfb*nrootz*ii                                10d19s21
              do ir=1,nrootz                                            10d19s21
               jtmp=imat+ncsfb*(ir-1)                                   10d19s21
               do i=0,ncsfb-1                                           10d19s21
                gb(jpb+i,ir)=gb(jpb+i,ir)+sum*bc(jtmp+i)                10d19s21
               end do                                                   10d19s21
              end do                                                    10d19s21
             end do                                                     9d10s21
             nok=0                                                         11d13s20
             do i=1,norb                                                   11d13s20
              itest(i,1)=0                                                   11d13s20
              itest(i,2)=0                                                    11d13s20
              if(btest(i0kc,i))itest(i,1)=2                              5d7s21
              if(btest(i0ko,i))itest(i,1)=1                              5d7s21
              if(btest(j0bc,i))itest(i,2)=2                              5d7s21
              if(btest(j0bo,i))itest(i,2)=1                              5d7s21
              ixn=min(itest(i,1),itest(i,2))
              if(ixn.gt.0)then                                             11d13s20
               nok=nok+1                                                   11d13s20
               itest(nok,1)=ixn                                            11d13s20
               itest(nok,2)=i                                              11d13s20
              end if                                                       11d13s20
             end do                                                        11d13s20
             if(lpr)write(6,*)('nok '),nok
             do ii=1,nok                                                    11d13s20
              if(itest(ii,1).eq.1.and.l2e.eq.0)then                     2d17s22
               if(lpr)write(6,*)('try orb '),itest(ii,2)
               itesta=j0bc                                               5d7s21
               itestb=j0bo                                               5d7s21
               nopenkk=nopenb                                                11d13s20
c
c     anihilate common
c
               if(btest(itesta,itest(ii,2)))then                             11d13s20
                itesta=ibclr(itesta,itest(ii,2))                             11d13s20
                itestb=ibset(itestb,itest(ii,2))                             11d13s20
                karg=jargb-1                                                11d13s20
                nopenkk=nopenkk+1                                             11d13s20
               else                                                         11d13s20
                itestb=ibclr(itestb,itest(ii,2))                             11d13s20
                karg=jargb                                                  11d13s20
                nopenkk=nopenkk-1                                             11d13s20
               end if                                                       11d13s20
c
c     create ket
c
               if(btest(itestb,nab4(2,1)))then                           11d27s20
                itesta=ibset(itesta,nab4(2,1))                           11d27s20
                itestb=ibclr(itestb,nab4(2,1))                           11d27s20
                karg=karg+1                                                11d13s20
                nopenkk=nopenkk-1                                             11d13s20
               else                                                         11d13s20
                itestb=ibset(itestb,nab4(2,1))                           11d27s20
                nopenkk=nopenkk+1                                             11d13s20
               end if                                                       11d13s20
               call gandcr(j0bc,j0bo,itesta,itestb,nopenb,                 5d10s21
     $              nopenkk,norb,nnot1,nab1,icode1,imap1,nx1,irw1,irw2, 9d16s21
     $              iwpb1,iwpk1,bc,ibc)                                 11d14s22
               call gandcr(itesta,itestb,i0kc,i0ko,nopenkk,nopenk,           5d10s21
     $         norb,nnot2,nab2,icode2,imap2,nx2,irw1,irw2,iwpb2,iwpk2,  11d14s22
     $              bc,ibc)                                             11d14s22
               if(lpr)write(6,*)('nnot1,nnot2 '),nnot1,nnot2
               if(nnot1.eq.2.and.nnot2.eq.2)then                        9d10s21
                call spinloop(i2sb,i2smb,i2sk,i2smk,nopenb,nopenk,      10d19s21
     $               nopenkk,ncsfb,nrootz,itype,imatx,ntype,nab1,iwpb1, 10d19s21
     $               iwpk1,nab2,iwpb2,iwpk2,veck(ipk,1),ncsftk,ieoro,bc,11d14s22
     $              ibc)                                                11d14s22
                do i=0,ntype-1
                 ipack=ibc(itype+i)
                 do j=1,4                                               9d9s21
                  if(ipack2(j).gt.0)then                                9d9s21
                   isy(j)=ism(ipack2(j))                                9d9s21
                   igya(j)=irel(ipack2(j))-1                             9d9s21
                   igy(j)=ipack2(j)
                   imy(j)=0                                             10d13s21
                  else                                                  9d9s21
                   isy(j)=ism(-ipack2(j))                                9d9s21
                   igy(j)=-ipack2(j)
                   igya(j)=irel(-ipack2(j))-1                           10d13s21
                   imy(j)=1                                             10d13s21
                  end if                                                9d9s21
                 end do                                                 9d9s21
                 ltest=imy(1)+2*(imy(2)+2*(imy(3)+2*imy(4)))            10d13s21
                 itestp=ltest+1                                         10d13s21
                 jtest=1-imy(1)+2*(1-imy(2)+2*(1-imy(3)+2*(1-imy(4))))  10d13s21
                 jtestp=jtest+1                                         10d13s21
                 xint=0d0                                               9d10s21
                 ihit=0
                 xintj=-2d0                                             10d13s21
                 if(iifmx(itestp).ge.0)then                             10d13s21
                  iad=i4or(isy(2),isy(3),isy(4))+igya(1)+irefo(isy(1))  10d13s21
     $                *(igya(2)+irefo(isy(2))*(igya(3)+irefo(isy(3))*   10d13s21
     $                (igya(4)+irefo(isy(4))*iifmx(itestp))))           10d13s21
                  xintj=bc(iad)                                         10d13s21
                  ihit=1                                                10d13s21
                  xint=xintj                                            10d13s21
                 else if(iifmx(jtestp).ge.0)then                        10d13s21
                  iad=i4or(isy(2),isy(3),isy(4))+igya(1)+irefo(isy(1))  10d13s21
     $                *(igya(2)+irefo(isy(2))*(igya(3)+irefo(isy(3))*   10d13s21
     $                (igya(4)+irefo(isy(4))*iifmx(jtestp))))           10d13s21
                  xintj=bc(iad)                                         10d13s21
                  if(isopt(4).ne.0)xintj=-xintj                         10d13s21
                  ihit=1                                                10d13s21
                  xint=xintj                                            10d13s21
                 end if                                                 10d13s21
                 ltest=imy(1)+2*(imy(4)+2*(imy(3)+2*imy(2)))            10d13s21
                 itestp=ltest+1                                         10d13s21
                 jtest=1-imy(1)+2*(1-imy(4)+2*(1-imy(3)+2*(1-imy(2))))  10d13s21
                 jtestp=jtest+1                                         10d13s21
                 xintk=-2d0                                             10d13s21
                 if(iifmx(itestp).ge.0)then                             10d13s21
                  iad=i4or(isy(4),isy(3),isy(2))+igya(1)+irefo(isy(1))  10d13s21
     $                *(igya(4)+irefo(isy(4))*(igya(3)+irefo(isy(3))*   10d13s21
     $                (igya(2)+irefo(isy(2))*iifmx(itestp))))           10d13s21
                  xintk=bc(iad)                                         10d13s21
                  ihit=1                                                10d13s21
                  xint=xint-xintk                                       10d13s21
                 else if(iifmx(jtestp).ge.0)then                        10d13s21
                  iad=i4or(isy(4),isy(3),isy(2))+igya(1)+irefo(isy(1))  10d13s21
     $                *(igya(4)+irefo(isy(4))*(igya(3)+irefo(isy(3))*   10d13s21
     $                (igya(2)+irefo(isy(2))*iifmx(jtestp))))           10d13s21
                  xintk=bc(iad)                                         10d13s21
                  if(isopt(4).ne.0)xintk=-xintk                         10d13s21
                  ihit=1                                                10d13s21
                  xint=xint-xintk                                       10d13s21
                 end if                                                 10d13s21
                 if(ihit.ne.0)then                                      9d10s21
                  itmpxq=imatx+ncsfb*nrootz*i                           10d19s21
                  do ir=1,nrootz                                        10d19s21
                   jtmpxq=itmpxq+ncsfb*(ir-1)                           10d19s21
                   do j=0,ncsfb-1                                       10d19s21
                    gb(jpb+j,ir)=gb(jpb+j,ir)+xint*bc(jtmpxq+j)         10d19s21
                   end do                                               10d19s21
                  end do                                                10d19s21
                 end if                                                 9d9s21
                end do
                ibcoff=itype                                            9d6s21
               end if                                                       11d13s20
              end if                                                    9d10s21
             end do                                                        11d13s20
            else if(l2e.eq.0)then                                       2d7s23
             nnotxx=nnot
             do i=1,2
              do j=1,2
               nab4xx(j,i)=nab4(j,i)
              end do
             end do
             nnot=0                                                     2d6s23
             if(ndifs.eq.4.and.ndifb.eq.4)then                          2d6s23
              nnot=4                                                    2d6s23
              ioxx(1)=1                                                 2d6s23
              ioxx(2)=1                                                 2d6s23
              do i=1,norb                                               2d6s23
               if(btest(gandcb,i))then                                  2d6s23
                if((btest(i0kc,i).and.btest(j0bo,i)).or.                2d8s23
     $               (btest(i0ko,i).and..not.btest(j0bc,i)))then        2d8s23
                 nab4(2,ioxx(2))=i                                      2d6s23
                 ioxx(2)=ioxx(2)+1                                      2d6s23
                else                                                    2d6s23
                 nab4(1,ioxx(1))=i                                      2d6s23
                 ioxx(1)=ioxx(1)+1                                      2d6s23
                end if                                                  2d6s23
               end if                                                   2d6s23
              end do                                                    2d6s23
             else if(ndifb.eq.3)then                                    2d6s23
              nnot=3                                                    2d6s23
              ioxx(1)=1                                                 2d6s23
              ioxx(2)=1                                                 2d6s23
              iswap=0                                                   2d6s23
              do i=1,norb                                               2d6s23
               if(btest(gandcb,i))then                                  2d6s23
                if(btest(gandcc,i).and.                                 2d6s23
     $        ((btest(j0bc,i).and..not.btest(i0ko,i)).or.               2d6s23
     $         (btest(i0kc,i).and..not.btest(j0bo,i))))then             2d6s23
                 if(btest(i0kc,i))iswap=1                               2d6s23
                 nab4(1,1)=i                                            2d6s23
                 nab4(1,2)=i                                            2d6s23
                else                                                    2d6s23
                 nab4(2,ioxx(2))=i                                      2d6s23
                 ioxx(2)=ioxx(2)+1                                      2d6s23
                end if                                                  2d6s23
               end if                                                   2d6s23
              end do                                                    2d6s23
              if(iswap.ne.0)then                                        2d6s23
               icpy=nab4(1,1)                                           2d6s23
               nab4(1,1)=nab4(2,1)                                      2d6s23
               nab4(2,1)=icpy                                           2d6s23
               icpy=nab4(1,2)                                           2d6s23
               nab4(1,2)=nab4(2,2)                                      2d6s23
               nab4(2,2)=icpy                                           2d6s23
               nbt=0                                                    2d6s23
               if(btest(i0kc,nab4(2,2)).and..not.btest(i0kc,nab4(2,1))) 2d6s23
     $              nbt=1                                               2d6s23
              else                                                      2d6s23
               nbt=0                                                    2d6s23
               if(btest(j0bc,nab4(1,2)).and..not.btest(j0bc,nab4(1,1))) 2d6s23
     $       nbt=1                                                      2d6s23
              end if                                                    2d6s23
              if(nbt.ne.0)then                                          2d6s23
               nab4(1,1)=nab4(1,2)                                      2d6s23
               nab4(2,1)=nab4(2,2)                                      2d6s23
              end if                                                    2d6s23
             else if(ndifs.eq.0.and.ndifd.eq.2)then                     2d6s23
              nnot=3                                                    2d6s23
              do i=1,norb                                               2d6s23
               if(btest(gandcb,i))then                                  2d6s23
                if(btest(j0bc,i))then                                   2d8s23
                 nab4(1,1)=i                                            2d6s23
                 nab4(1,2)=i                                            2d6s23
                else                                                    2d6s23
                 nab4(2,1)=i                                            2d6s23
                 nab4(2,2)=i                                            2d6s23
                end if                                                  2d6s23
               end if                                                   2d6s23
              end do                                                    2d6s23
             end if                                                     2d6s23
             if(nnot.ne.0)then                                          2d7s23
              iu1=1                                                     1d22s21
              iu2=1                                                     1d22s21
              itesta=j0bc                                                5d7s21
              itestb=j0bo                                                5d7s21
              if(btest(itesta,nab4(1,iu1)))then                          1d22s21
               itesta=ibclr(itesta,nab4(1,iu1))                          1d22s21
               itestb=ibset(itestb,nab4(1,iu1))                          1d22s21
               nopenkk=nopenb+1                                           1d22s21
               karg=jargb-1                                               1d22s21
              else if(btest(itestb,nab4(1,iu1)))then                     1d22s21
               itestb=ibclr(itestb,nab4(1,iu1))                          1d22s21
               nopenkk=nopenb-1                                              11d13s20
               karg=jargb                                                  11d13s20
              else                                                          11d13s20
               write(6,*)('bit not set for nab4(1,1) = '),nab4(1,iu1)      11d27s20
               stop 'nab4(1,1)'                                          11d27s20
              end if                                                        11d13s20
              if(btest(itestb,nab4(2,iu2)))then                          1d22s21
               itesta=ibset(itesta,nab4(2,iu2))                          1d22s21
               itestb=ibclr(itestb,nab4(2,iu2))                          1d22s21
               nopenkk=nopenkk-1                                              11d13s20
               karg=karg+1                                                  11d13s20
              else if(btest(itesta,nab4(2,iu2)))then                     1d22s21
               write(6,*)('already double in nab4(2,1) = '),nab4(2,iu2)  1d22s21
               stop 'nab4(2,1)'                                          11d27s20
              else                                                          11d13s20
               itestb=ibset(itestb,nab4(2,iu2))                          1d22s21
               nopenkk=nopenkk+1                                              11d13s20
              end if                                                        11d13s20
c     iarg=nclop-mdon. thus iarg+mdon-1=nclo,
c     so nqq=nec-nopenkk
              nqq=nec-nopenkk                                           9d6s21
              nqqp=nqq+1
               call gandcr(j0bc,j0bo,itesta,itestb,nopenb,                 5d10s21
     $              nopenkk,norb,nnot1,nab1,icode1,imap1,nx1,irw1,irw2, 9d16s21
     $             iwpb1,iwpk1,bc,ibc)                                  11d14s22
               call gandcr(itesta,itestb,i0kc,i0ko,nopenkk,nopenk,           5d10s21
     $         norb,nnot2,nab2,icode2,imap2,nx2,irw1,irw2,iwpb2,iwpk2,  11d14s22
     $              bc,ibc)                                             11d14s22
               if(nnot1.eq.2.and.nnot2.eq.2)then                         1d22s21
                call spinloop(i2sb,i2smb,i2sk,i2smk,nopenb,nopenk,      10d19s21
     $               nopenkk,ncsfb,nrootz,itype,imatx,ntype,nab1,iwpb1, 10d19s21
     $               iwpk1,nab2,iwpb2,iwpk2,veck(ipk,1),ncsftk,ieoro,bc,11d14s22
     $              ibc)                                                11d14s22
                do i=0,ntype-1
                 ipack=ibc(itype+i)
                 itmpxq=imatx+ncsfb*nrootz*i                            10d19s21
                 ltest=0                                                10d12s21
                 jtest=0                                                10d12s21
                 mtest=1                                                10d12s21
                 do j=1,4                                               9d9s21
                  if(ipack2(j).gt.0)then                                9d9s21
                   isy(j)=ism(ipack2(j))                                9d9s21
                   igya(j)=irel(ipack2(j))-1                             9d9s21
                   igy(j)=ipack2(j)                                       10d12s21
                   jtest=jtest+mtest                                    10d12s21
                   imy(j)=0                                             10d12s21
                  else                                                  9d9s21
                   isy(j)=ism(-ipack2(j))                                9d9s21
                   igy(j)=-ipack2(j)
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
                  iad=i4or(isy(2),isy(3),isy(4))+igya(1)+irefo(isy(1))* 10d12s21
     $          (igya(2)+irefo(isy(2))*(igya(3)+irefo(isy(3))*(         10d12s21
     $                igya(4)+irefo(isy(4))*iifmx(itestp))))            10d12s21
                  xintj=bc(iad)                                          10d12s21
                  ihit=1                                                10d12s21
                  xint=xintj                                            10d12s21
                  jchoice=iifmx(itestp)
                 else if(iifmx(jtestp).ge.0)then                        10d12s21
                  iad=i4or(isy(2),isy(3),isy(4))+igya(1)+irefo(isy(1))* 10d12s21
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
                  jtest=1-imy(1)+2*(1-imy(4)+2*(1-imy(3)+2*(1-imy(2))))  10d12s21
                  jtestp=jtest+1                                         10d12s21
                  if(iifmx(itestp).ge.0)then                             10d12s21
                   iad=i4or(isy(4),isy(3),isy(2))+igya(1)+irefo(isy(1))* 10d12s21
     $          (igya(4)+irefo(isy(4))*(igya(3)+irefo(isy(3))*(         10d12s21
     $                igya(2)+irefo(isy(2))*iifmx(itestp))))            10d12s21
                   xintk=bc(iad)                                          10d12s21
                   xint=xint-xintk                                       10d12s21
                   kchoice=iifmx(itestp)
                   ihit=1                                               11d9s21
                  else if(iifmx(jtestp).ge.0)then                        10d12s21
                   iad=i4or(isy(4),isy(3),isy(2))+igya(1)+irefo(isy(1))* 10d12s21
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
                  if(nab4(1,1).eq.nab4(1,2).and.nab4(2,1).eq.nab4(2,2)) 9d9s21
     $                 xint=xint*0.5d0                                  9d9s21
                  do ir=1,nrootz                                        10d19s21
                   jtmpxq=itmpxq+ncsfb*(ir-1)                           10d19s21
                   do j=0,ncsfb-1                                       10d19s21
                    gb(jpb+j,ir)=gb(jpb+j,ir)+xint*bc(jtmpxq+j)         10d19s21
                   end do                                               10d19s21
                  end do                                                10d19s21
                 end if                                                 9d9s21
                end do
                ibcoff=itype                                            9d6s21
               end if                                                    1d26s21
              end if                                                    2d7s23
            end if                                                      9d6s21
           end if                                                          12d9s19
           jpb=jpb+ncsfb                                                 8d31s21
          end do                                                         5d7s21
         end if                                                         3d13s23
         loopit=loopit+1                                                3d13s23
         ipk=ipk+ncsfk                                                  8d31s21
        end do                                                          5d7s21
       end do                                                           6d11s19
      end do
      ibcoff=ibcoffo                                                    10d22s21
      if(lpr)write(6,*)('goal at end: '),gb(igoal,igoalr)
      if(ldebug)then
       write(6,*)('gb at end: ')
       call prntm2(gb,ncsftb,nrootz,ncsftb)
       itmp=ibcoff
       ibcoff=itmp+ncsftb*nrootz
       call enough('hccsfbk4.  2',bc,ibc)
       jtmp=itmp-1
       do ir=1,nrootz
        do k=1,ncsftb
         bc(jtmp+k)=gb(k,ir)
        end do
        jtmp=jtmp+ncsftb
       end do
       call dws_gsumf(bc(itmp),ncsftb*nrootz)
       write(6,*)('global summed ')
       call prntm2(bc(itmp),ncsftb,nrootz,ncsftb)
      end if
      if(lpr)then
       ipb=1                                                            7d16s21
       rsum=0d0
       do nclop=mdon+1,mdoobp
        iarg=nclop-mdon
        iff=nff0b(nclop,2)
        ipb=nff0b(nclop,3)
        nopen=nec-(nclop-1)*2
        call wfetch(nopen,mdon,idum,i2sb,i2smb,ndetb,ncsfb,ivecb,iaorbb,8d30s21
     $      idum,bc,ibc)                                                11d9s22
        do if=1,nff0b(nclop,1)
         i0bc=mff0b(iff)
         iff=iff+1
         i0bo=mff0b(iff)
         iff=iff+1
         sz=0d0
         do ir=1,nrootz
          do i=0,ncsfb-1
           iad=itmp+ipb+i-1+ncsftb*(ir-1)
           sz=sz+bc(iad)**2
          end do
         end do
         sz=sqrt(sz/dfloat(nrootz*ncsfb))
         if(sz.gt.1d-14)then
          call dcbit(i0bc,norb,'c')
          mvirt=0
          do i=n0virt,norb                                                   7d22s21
           if(btest(i0bc,i))then
            write(6,*)irel(i),('s'),ism(i)
            mvirt=mvirt+2
           end if
          end do
          call dcbit(i0bo,norb,'o')
          do i=n0virt,norb                                                   7d22s21
           if(btest(i0bo,i))then
            write(6,*)irel(i),('s'),ism(i)
            mvirt=mvirt+1
           end if
          end do
          if(mvirt.eq.2)then                                             7d16s21
           write(6,*)('bra '),nclop,if,ipb,ncsfb,ipb-igoal,igoalr
           iad=itmp+ipb-1
           call prntm2(bc(iad),ncsfb,nrootz,ncsftb)
           do i=0,ncsfb-1
            iad=itmp+ipb+i-1
            rsum=rsum+bc(iad)**2
           end do
           write(6,*)('rsum so far: '),rsum
          end if
         end if
         ipb=ipb+ncsfb
        end do
       end do
      end if
      return
      end                                                               7d11s19
