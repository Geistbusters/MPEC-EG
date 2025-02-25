c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcdsbk4(nff22,iff22,ff22,nfdat,ncsfb,ncsfb2,gd,i2sb,   12d6s21
     $     i2smb,mdoobp,isymbra,                                        12d6s21
     $     ihsdiagk,nff1,iff1,i2sk,i2smk,mdookp,ncsfk,isymket,          12d6s21
     $     nec,mdon,nsymb,multh,irw1,irw2,nvirt,nrootu,irorip,isopt,ism,12d7s21
     $     irel,irefo,norb,ih0n,nh0,ionex,i3x,iifmx,ntype,maxbx,sr2,l2e,11d14s22
     $     bc,ibc)                                                      11d14s22
      implicit real*8 (a-h,o-z)                                         12d18s20
c
c     to do ...
c     consolidate densities, and perhaps njhere in density calculation
c     did I swap iarg and jarg in 4th gandc4?
c
      integer*8 ihsdiagk(mdookp,nsymb,2),i18,i28,                       12d7s21
     $   i38,i48,i1c,i1o,j2c,j2o,itestc,itesto,ipack8,last8(2),iff22(*),8d16s21
     $     ipackc,ipack,gandcc,gandco,gandcb                            2d7s23
      integer*1 nab1(2),nab2(2),imap(64),icode(64),ipackc1(8)           12d6s21
      integer*2 ipack2(4)                                               12d6s21
      external second                                                   2d18s21
      logical l3x,lprt,lnew,ldebug,lchoice                              3d17s21
      equivalence (ipack8,ipack4),(ipackc,ipackc1),(ipack,ipack2)       12d6s21
      dimension nff1(mdookp,nsymb,2),iff1(*),nff22(mdoobp,2,nsymb),     12d6s21
     $     nfdat(5,4,*),gd(*),multh(8,8),nvirt(*),ncsfb(*),ncsfb2(4,*), 12d7s21
     $     ncsfk(*),irel(*),ism(*),irefo(*),itest(32,3),                12d6s21
     $     nab4(2,3),ipack4(2),nl(4),npre(4),mpre(4),mcsf(2),           12d6s21
     $     id1visv(8,8),nd1visv(8,8),nokdc(8,8,4),nok3v(4),             12d22s20
     $     idhvnotv(4),ndhvnotv(4),id1vnotv(4,8,8),nd1vnotv(4,8,8),     12d21s20
     $     id3vnotv3(4),nd3vnotv3(4),loff(4),idkeep(2),ndkeep(2),       12d7s21
     $     mdkeep(2),keep(2),idhvnotvf(4),ibmat(8),ivmat(8),            2d4s21
     $     mdhvnotv(4),md3vnotv3(4),nok4f(4),nok4(4),nok3vf(4),         12d7s21
     $     nok33f(4,2),nok33(4,2),ff22(*),isopt(*),ih0n(*),nh0(*),      12d6s21
     $     ionex(8,8,8),i3x(8,8,8),iifmx(*),iwpb1(4),iwpk1(4),iwpb2(4), 12d6s21
     $     iwpk2(4),imy(4),igya(4),isy(4),ipack2a(4),i3xk(2,16),        12d9s21
     $     phss(2),ioxx(2)                                              2d7s23
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data loopx/43140000/
      data phss/1d0,-1d0/                                               12d9s21
      data icall/0/
      save
      icall=icall+1
      irori=irorip-1                                                    11d17s21
      if(isopt(3).ne.0)then
       ieoro=1-irori                                                      9d6s21
      else                                                              10d14s21
       ieoro=irori                                                      9d6s21
      end if                                                            9d6s21
      loop=0
      nrootm=nrootu-1                                                   1d4s21
      idoit=0                                                           3d1s21
      mdoob=mdoobp-1                                                      8d13s21
      mdook=mdookp-1                                                      8d13s21
      igoal=1
      igoal1=5
      igoal2=5
      igoal3=5
      igoal4=2
      igoal5=2
      igoal6=5
      igoal7=2
      ivgoal=5
      igoalv=5
      igoalvb=5
      igoalvc=5
      igoalt=1
      last8(1)=-1                                                       2d8s21
      ircv=ibcoff                                                       1d29s21
      iacc=ircv+mynprocg                                                1d30s21
      ivs=iacc+mynprocg                                                 1d30s21
      igg=ivs+maxbx                                                     1d30s21
      ibcoff=igg+maxbx                                                  1d30s21
      nacc=0                                                            1d30s21
      call enough('hcdsbk4.  1',bc,ibc)
      itransgg=0                                                        1d30s21
      loop=0
      nsing=0                                                           12d23s20
      norbx=norb+1                                                      12d18s20
      norbxx=norbx+1                                                    12d18s20
      norbxxx=norbxx+1                                                  12d18s20
      isymbk=multh(isymbra,isymket)                                     8d13s21
      do jsb=1,nsymb                                                    12d18s20
       ksbv=multh(jsb,isymket)                                          8d13s21
c
c     let us form bvrkn=[(vv"|nv')+p(k)(vv'|nv")]Vv'v"rk,               2d3s21
c     v ne v' and v"                                                    2d3s21
c     and space for vvrkn=Vvrj*Djkn                                     2d4s21
c     gandc4(Vsv,Vdv'v"), idx=1 (iv'|vv"), idx=2 (iv"|vv')
c
       ioffvd=1                                                         2d3s21
       ibcvmat=ibcoff                                                   12s15s21
       nbmat=0                                                          12d6s21
       do isb=1,nsymb                                                   3d2s21
        isbv12=multh(isb,isymbra)                                       8d13s21
        isn=multh(isopt(1),multh(isbv12,ksbv))                          12d7s21
        ivmat(isb)=ibcoff                                               2d3s21
        if(min(irefo(isn),nvirt(ksbv)).gt.0)then                        12d6s21
         nftrip=nfdat(2,2,isb)+nfdat(2,3,isb)+nfdat(2,4,isb)            12d6s21
         nfh=nfdat(2,1,isb)+nftrip                                      12d6s21
         nn=nvirt(ksbv)*nrootu*nfh*ntype                                12d6s21
         nbmat=nbmat+nn*irefo(isn)                                      12d6s21
         ibcoff=ivmat(isb)+nn*irefo(isn)                                3d2s21
        end if                                                          12d6s21
       end do                                                           3d2s21
       call enough('hcdsbk4.  2',bc,ibc)
       do i=ibcvmat,ibcoff-1                                            12s15s21
        bc(i)=0d0                                                       3d2s21
       end do                                                           3d2s21
       do it=1,16                                                       1d18s23
        i3xk(1,it)=0                                                    1d18s23
        i3xk(2,it)=0                                                    1d18s23
       end do                                                           1d18s23
       do nclo1p=mdon+1,mdookp                                           8d13s21
        if(min(nff1(nclo1p,jsb,1),nvirt(ksbv)).gt.0)then                8d26s21
         nclo1=nclo1p-1                                                  12d12s20
         jarg=nclo1p-mdon                                                12d12s20
         nopen1=nec-2*nclo1                                              12d12s20
         nggk=nvirt(ksbv)*nrootu                                        8d13s21
         nggg=nff1(nclo1p,jsb,1)*ncsfk(jarg)                            12d6s21
         ncolt=nggg*nggb                                                8d13s21
         x2=dfloat(nff1(nclo1p,jsb,1))/dfloat(mynprocg)                 3d2s21
         call ilimts(1,nff1(nclo1p,jsb,1),mynprocg,mynowprog,           3d2s21
     $        i1l,i1h,i11s,i11e,i12s,i12e)                              12d19s20
         i1l=1+ncsfk(jarg)*(i1l-1)                                      3d2s21
         i1h=ncsfk(jarg)*i1h                                            3d2s21
         ibcgg=ibcoff                                                   1d30s21
         i18=1
         i28=nggk                                                       8d13s21
         i38=1                                                          12d18s20
         i48=nggg                                                       12d29s20
         call ddi_iget(bc,ibc,ihsdiagk(nclo1p,jsb,1),i18,i28,i38,i48,   11d15s22
     $        bc(ivs),ibc(ircv),nrcv)                                   11d15s22
         itransvs=0                                                     1d29s21
         ioffdnon=1                                                     12d21s20
         koffdnon=1                                                     8d16s21
         do isb=1,nsymb                                                 12d18s20
          nfh=nfdat(2,1,isb)+nfdat(2,2,isb)+nfdat(2,3,isb)              2d3s21
     $         +nfdat(2,4,isb)                                          2d3s21
          isbv12=multh(isb,isymbra)                                     8d13s21
          loff(1)=0                                                     1d4s21
          do l=2,4                                                      1d4s21
           lm=l-1                                                       1d4s21
           loff(l)=loff(lm)+nfdat(2,lm,isb)                             1d4s21
          end do                                                        1d4s21
          nnt=loff(4)+nfdat(2,4,isb)                                    1d4s21
          nnt=nnt*ntype                                                 12d6s21
          do nclo2p=max(mdon+1,nclo1p-2),min(mdoobp,nclo1p+3)           12d6s21
           if(nff22(nclo2p,1,isb).gt.0)then                             12d18s20
            nclo2=nclo2p-1                                              12d18s20
            iarg=nclo2p-mdon                                            12d18s20
            iargp=iarg+1                                                12d18s20
            nopen2=nec-2*nclo2p                                         12d18s20
            nopen2p=nopen2+2                                            12d18s20
            idhvisv=ibcoff                                              12d18s20
            nnj=ncsfk(jarg)*nfdat(2,1,isb)                              12d6s21
            lgoal=multh(jsb,isb)                                        12d22s20
            lgoal3=multh(isopt(1),multh(ksbv,isbv12))                   12d8s21
            do isc=1,nsymb                                              12d18s20
             do l=1,4                                                   12d18s20
              nnl=ncsfk(jarg)*nfdat(2,l,isb)*irefo(isc)                 12d8s21
              if(isc.eq.lgoal)then                                      12d21s20
               idhvnotvf(l)=ibcoff                                      1d6s21
               idhvnotv(l)=idhvnotvf(l)+nnl                             1d6s21
               ndhvnotv(l)=idhvnotv(l)+nnl                              12d21s20
               mdhvnotv(l)=ndhvnotv(l)+irefo(isc)                       12d31s20
               ibcoff=mdhvnotv(l)+nfdat(2,l,isb)                        12d21s20
              end if                                                    12d21s20
             end do                                                     12d18s20
             iscv=multh(isc,ksbv)                                       8d26s21
             jscv=multh(isc,lgoal)                                      12d21s20
             do isd=1,nsymb                                             8d16s21
              iscdv=multh(iscv,isd)                                     12d18s20
              jscdv=multh(jscv,isd)                                     12d21s20
              nn=irefo(isd)*irefo(isc)                                  8d13s21
              nnn=nn*irefo(jscdv)*ntype                                 12d6s21
              do l=1,4                                                  12d19s20
               id1vnotv(l,isd,isc)=ibcoff                               12d19s20
               nd1vnotv(l,isd,isc)=id1vnotv(l,isd,isc)+nnn*ncsfk(jarg)  12d6s21
     $              *nfdat(2,l,isb)                                     12d19s20
               ibcoff=nd1vnotv(l,isd,isc)+nnn                           12d21s20
              end do                                                    12d19s20
             end do                                                     12d18s20
            end do                                                      12d18s20
            do ipass=1,1                                                3d1s21
             idkeep(ipass)=ibcoff                                       1d4s21
             ndkeep(ipass)=idkeep(ipass)+irefo(lgoal3)*nnt*ncsfk(jarg)  12d6s21
             mdkeep(ipass)=ndkeep(ipass)+irefo(lgoal3)*nnt               8d13s21
             ibcoff=mdkeep(ipass)+irefo(lgoal3)*nnt                      8d13s21
             call enough('hcdsbk4.  3',bc,ibc)
            end do                                                      1d4s21
            do l=1,4
             nnl=ntype*ncsfk(jarg)*nfdat(2,l,isb)*irefo(lgoal3)         12d6s21
             id3vnotv3(l)=ibcoff                                        12d7s21
             nd3vnotv3(l)=id3vnotv3(l)+nnl                              12d7s21
             md3vnotv3(l)=nd3vnotv3(l)+ntype*irefo(lgoal3)              12d8s21
             ibcoff=md3vnotv3(l)+nfdat(2,l,isb)*ntype                   12d8s21
            end do                                                      12d22s20
            ibcb4=ibcoff-1                                              12d18s20
            if1o=nff1(nclo1p,jsb,2)                                      12d14s20
            jvs=ivs                                                     12d18s20
            do if1=1,nff1(nclo1p,jsb,1)                                  12d14s20
             ist=1+ncsfk(jarg)*(if1-1)                                  12d7s21
             ien=ncsfk(jarg)*if1                                        12d7s21
             istu=max(ist,i1l)                                          1d5s21
             ienu=min(ien,i1h)                                          1d5s21
             if(ienu.ge.istu)then                                       1d5s21
              i11s=istu-ist+1                                           1d5s21
              njhere=ienu+1-istu                                        1d5s21
             else                                                       1d5s21
              njhere=0                                                  1d5s21
             end if                                                     1d5s21
             idoit=idoit+1                                              3d1s21
             joffdnon=ioffdnon                                          12d21s20
             loffdnon=koffdnon                                          8d16s21
             do i=idhvisv,ibcb4                                         12d18s20
              bc(i)=0d0                                                 12d18s20
             end do                                                     12d18s20
             i1c=iff1(if1o)                                                11d25s20
             ntest=popcnt(i1c)
             do i=1,norbxx                                              12d19s20
              itest(i,1)=0                                                 11d25s20
             end do                                                        11d25s20
             do i=1,norb                                                   11d25s20
              if(btest(i1c,i))then                                         11d25s20
               itest(i,1)=2                                                11d25s20
              end if                                                       11d25s20
             end do                                                        11d25s20
             if1o=if1o+1                                                   11d25s20
             i1o=iff1(if1o)                                                11d25s20
             if1o=if1o+1                                                   11d25s20
             i1o=ibset(i1o,norbx)                                          11d25s20
             do i=1,norb                                                   11d25s20
              if(btest(i1o,i))then                                         11d25s20
               itest(i,1)=1                                                11d25s20
              end if                                                       11d25s20
             end do                                                        11d25s20
             itest(norbx,1)=1                                           12d19s20
             ivcv=nfdat(5,1,isb)                                        12d18s20
             jvcv=ivcv+nff22(nclo2p,2,isb)                              12d18s20
             do if2=1,nff22(nclo2p,1,isb)                               12d18s20
              ipack8=iff22(jvcv)                                        8d16s21
              j2c=ipack4(1)                                             12d18s20
              j2o=ipack4(2)                                             12d18s20
              nclo=popcnt(ipack4(1))                                    12d18s20
              nspace=iff22(jvcv+1)                                      8d16s21
              lchoice=.false.                                           3d18s21
              do l=1,4                                                   3d17s21
               nl(l)=iff22(jvcv+1+l)                                    8d16s21
               if(nl(l).gt.ncsfb2(l,iarg))lchoice=.true.                1d13s23
              end do                                                     3d17s21
              j2o=ibset(j2o,norbx)                                      12d18s20
              j2o=ibset(j2o,norbxx)                                     12d18s20
              if(njhere.gt.0)then                                       2d26s21
               gandcc=ieor(i1c,j2c)                                     12d21s22
               gandco=ieor(i1o,j2o)                                     12d21s22
               gandcb=ior(gandcc,gandco)                                12d21s22
               ndifb=popcnt(gandcb)                                     12d21s22
               if(ndifb.le.4)then                                       12d21s22
                ndifd=popcnt(gandcc)                                     10d13s22
                ndifs=popcnt(gandco)                                     10d13s22
                if(ndifs.eq.2.and.ndifb.eq.2)then
                 do i=1,norbxx
                  if(btest(gandco,i))then                                          10d14s22
                   if((btest(j2o,i).and..not.btest(i1c,i)).or.
     $                (btest(j2c,i).and.btest(i1o,i)))then                           10d14s22
                    nab4(1,1)=i
                   else                                                           10d14s22
                    nab4(2,1)=i
                   end if
                  end if                                                          10d14s22
                 end do                                                           10d14s22
c
c     i1o is rest,a,x
c     j2o is rest,x,xx
c     so nab1(1) is xx,nab1(2) is a.                                    12d6s21
c     a*xx should have sym isymbra*isymket
c     rest*a*ksbv=isymket, rest*a*jsb*isymket=isymket
c     rest*isbv12=isymbra, rest*isb*isymbra=isymbra,
c     or rest*a*jsb=1=rest*isb or a*jsb=isb or a=jsb*isb
c
                 do i=1,norb                                             12d18s20
                  itest(i,2)=0                                           12d18s20
                 end do                                                  12d18s20
                 do i=1,norb                                             12d18s20
                  if(btest(j2c,i))then                                   12d18s20
                   itest(i,2)=2                                          12d18s20
                  end if                                                 12d18s20
                  if(btest(j2o,i))then                                   12d18s20
                   itest(i,2)=1                                          12d18s20
                  end if                                                 12d18s20
                 end do                                                  12d18s20
                 itest(norbx,2)=1                                        12d19s20
                 itest(norbxx,2)=1                                       12d19s20
                 nok=0                                                         11d13s20
                 do i=1,norbxx                                           12d19s20
                  ixn=min(itest(i,1),itest(i,2))
                  if(ixn.gt.0)then                                             11d13s20
                   nok=nok+1                                                   11d13s20
                   itest(nok,3)=ixn                                            11d13s20
                   itest(nok,2)=i                                              11d13s20
                  end if                                                       11d13s20
                 end do                                                        11d13s20
                 iunit=ibcoff                                            12d6s21
                 ibcoff=iunit+ncsfk(jarg)*ncsfk(jarg)                    12d6s21
                 call enough('hcdsbk4.  4',bc,ibc)
                 do iz=iunit,ibcoff-1                                    12d6s21
                  bc(iz)=0d0                                             12d6s21
                 end do                                                  12d6s21
                 do iz=0,ncsfk(jarg)-1                                   12d6s21
                  iad=iunit+iz*(ncsfk(jarg)+1)                           12d6s21
                  bc(iad)=1d0                                            12d6s21
                 end do                                                  12d6s21
                 itrans=ibcoff                                           12d7s21
                 ibcoff=itrans+ncsfb(iarg)*ncsfk(jarg)                   12d7s21
                 call enough('hcdsbk4.  5',bc,ibc)
                 call gandcr(j2c,j2o,i1c,i1o,nopen2p,nopen1,norbxx,     2d7s23
     $            nnot1,nab1,icode,imap,nx,irw1,irw2,iwpb1,iwpk1,bc,ibc)2d7s23
                 call gencup(i2sb,i2smb,i2sk,i2smk,nopen2p,nopen1,nab1,  12d6s21
     $               iwpb1,iwpk1,ioutg,imatg,ntypeg,mcsf,bc(iunit),     12d6s21
     $               ncsfk(jarg),ncsfk(jarg),bc,ibc)                    11d14s22
                 do ixi=0,ntypeg-1                                        12d6s21
                  ipackc=ibc(ioutg+ixi)                                   12d6s21
                  imatu=imatg+ncsfb(iarg)*ncsfk(jarg)*ixi                12d7s21
                  do i=0,ncsfb(iarg)-1                                   12d6s21
                   do k=0,ncsfk(jarg)-1                                  12d6s21
                    ki=itrans+k+ncsfk(jarg)*i                            12d6s21
                    ik=imatu+i+ncsfb(iarg)*k                             12d6s21
                    bc(ki)=bc(ik)                                        12d6s21
                   end do                                                12d6s21
                  end do                                                 12d6s21
                  if(ipackc1(1).gt.0)then                                   9d10s21
                   imy(3)=0                                              12d6s21
                  else                                                   12d6s21
                   imy(3)=1                                              12d6s21
                  end if                                                 12d6s21
                  if(ipackc1(2).gt.0)then                                12d6s21
                   imy(4)=0                                              12d6s21
                   ipack2a(4)=ipackc1(2)
                  else                                                   12d6s21
                   imy(4)=1                                              12d6s21
                   ipack2a(4)=-ipackc1(2)                                12d6s21
                  end if                                                 12d6s21
                  isy(4)=ism(ipack2a(4))                                 12d6s21
                  igya(4)=irel(ipack2a(4))-1                             12d6s21
                  phsh=1d0                                               12d6s21
                  isy(3)=multh(isy(4),isopt(1))                          12d7s21
                  if(ih0n(isy(3)).gt.0)then                              12d6s21
                   if(imy(3).ne.0.and.isopt(4).ne.0)phsh=-phsh           12d6s21
                  else                                                   12d6s21
                   if(isopt(2).ne.0)phsh=-phsh                           12d6s21
                   if(imy(4).ne.0.and.isopt(4).ne.0)phsh=-phsh           12d6s21
                  end if                                                 12d6s21
                  jmat=itrans                                            12d6s21
                  do l=1,4                                                12d18s20
                   if(nl(l).gt.0)then                                     12d18s20
                    nnl=ncsfk(jarg)*nfdat(2,l,isb)                       12d6s21
                    iad1=jvcv+iff22(jvcv+5+l)                             8d16s21
                    iad2=iad1+nl(l)                                       3d19s21
                    itmp=ibcoff                                           12d18s20
                    ibcoff=itmp+ncsfk(jarg)*nl(l)                        12d6s21
                    call enough('hcdsbk4.  6',bc,ibc)
                    call dgemm('n','n',ncsfk(jarg),nl(l),ncsfb2(l,iarg), 12d6s21
     $               1d0,bc(jmat),ncsfk(jarg),ff22(iad2),ncsfb2(l,iarg),8d16s21
     $                 0d0,bc(itmp),ncsfk(jarg),                        12d6s21
     d' hcdsbk4.  1')
c     jsbv=jsb*isymbra, ksbv=jsb*isymket                                8d13s21
c     isbv12=isb*isymbra, ksbv12=jsb*isymket                            8d13s21
c     isbv(1or2)=ksbv12*jsbv=jsb*isymbra*isb*isymket
c     isbv(1or2)=isbv12*ksbv=jsb*isymket*isb*isymbra                    8d13s21
                    jden=idhvnotv(l)+nnl*igya(4)                         12d6s21
                    ibc(ndhvnotv(l)+igya(4))=1                           12d6s21
                    jtmp=itmp                                             12d18s20
                    do iii=0,nl(l)-1                                      12d18s20
                     ip=iff22(iad1+iii)-1                                 8d16s21
                     jjden=jden+ncsfk(jarg)*ip                             12d18s20
                     ibc(mdhvnotv(l)+ip)=1                                12d21s20
                     do j=0,ncsfk(jarg)-1                                  12d18s20
                      bc(jjden+j)=bc(jjden+j)+bc(jtmp+j)*phsh            12d6s21
                     end do                                               12d18s20
                     jtmp=jtmp+ncsfk(jarg)                               12d6s21
                    end do                                                12d18s20
                    do i=1,norb                                          12d6s21
                     if(btest(j2c,i).and.btest(i1c,i).and.l2e.eq.0)then  2d17s22
                      js=ism(i)                                               9d10s21
                      jg=irel(i)-1                                            9d10s21
                      nn=irefo(js)*irefo(js)                             12d6s21
                      do jpass=0,1                                       11d22s21
c
c     (v4|ii)=(ii|v4)                                                           12d6s21
c
                       ltest=jpass+2*(jpass+2*(imy(3)+2*imy(4)))         12d6s21
                       jtest=1-jpass+2*(1-jpass+2*(1-imy(3)              12d7s21
     $                     +2*(1-imy(4))))                              12d7s21
                       itestp=ltest+1                                         10d12s21
                       jtestp=jtest+1                                    11d22s21
                       iuse=-1                                           11d22s21
                       phsj=1d0                                          11d22s21
                       if(iifmx(itestp).ge.0)then                        11d22s21
                        iuse=iifmx(itestp)                               11d22s21
                       else if(iifmx(jtestp).ge.0)then                   11d22s21
                        iuse=iifmx(jtestp)                               11d22s21
                        if(isopt(4).ne.0)phsj=-phsj                      11d22s21
                       end if                                            11d22s21
                       if(iuse.ge.0)then                                 11d22s21
                        icolj=jg+irefo(js)*(jg+irefo(js)*(igya(4)+       12d6s21
     $                     irefo(isy(4))*iuse))                         12d6s21
                        jdenj=id1vnotv(l,js,js)+nnl*icolj                    12d19s20
                        ibc(nd1vnotv(l,js,js)+icolj)=1                       12d19s20
                        jtmp=itmp                                              12d18s20
                        do iii=0,nl(l)-1                                         12d18s20
                         ii=iff22(iad1+iii)-1                              8d16s21
                         jdj=jdenj+ncsfk(jarg)*ii                             12d18s20
                         do j=0,ncsfk(jarg)-1                                   12d18s20
                          bc(jdj+j)=bc(jdj+j)+bc(jtmp+j)*phsj            12d6s21
                         end do                                                12d18s20
                         jtmp=jtmp+ncsfk(jarg)                           12d7s21
                        end do                                                 12d18s20
                       end if                                            12d6s21
c
c why not (vi|i4)=(i4|vi)
c     (iv|4i)=(4i|iv)=+/-(i4|vi)
c
                       ltest=jpass+2*(imy(4)+2*(imy(3)+2*jpass))         12d6s21
                       jtest=1-jpass+2*(1-imy(4)+2*(1-imy(3)             12d6s21
     $                     +2*(1-jpass)))                               12d6s21
                       itestp=ltest+1                                    12d6s21
                       jtestp=jtest+1                                    12d6s21
                       iuse=-1                                           12d6s21
                       phsk=-1d0                                         12d6s21
                       if(iifmx(itestp).ge.0)then                        12d6s21
                        iuse=iifmx(itestp)                               12d6s21
                       else if(iifmx(jtestp).ge.0)then                   12d6s21
                        iuse=iifmx(jtestp)                               12d6s21
                        if(isopt(4).ne.0)phsk=-phsk                      12d6s21
                       end if                                            12d6s21
                       if(iuse.ge.0)then                                 12d6s21
                        icolk=jg+irefo(js)*(igya(4)+irefo(isy(4))*(jg+    12d6s21
     $                     irefo(js)*iuse))                             12d6s21
                        jdenj=id1vnotv(l,js,isy(4))+nnl*icolk                    12d19s20
                        ibc(nd1vnotv(l,js,isy(4))+icolk)=1               12d6s21
                        jtmp=itmp                                              12d18s20
                        do iii=0,nl(l)-1                                         12d18s20
                         ii=iff22(iad1+iii)-1                              8d16s21
                         jdj=jdenj+ncsfk(jarg)*ii                             12d18s20
                         do j=0,ncsfk(jarg)-1                                   12d18s20
                          bc(jdj+j)=bc(jdj+j)+bc(jtmp+j)*phsk            12d6s21
                         end do                                                12d18s20
                         jtmp=jtmp+ncsfk(jarg)                           12d7s21
                        end do                                                 12d18s20
                       end if                                            12d6s21
                      end do                                             12d6s21
                     end if                                              12d6s21
                    end do                                               12d6s21
                    if(btest(i1c,ipack2a(4)).and.l2e.eq.0)then           2d17s22
c
c     y term
c         --   --
c     (v4|44)=(44|v4)
c
                     ltest=1-imy(4)+2*(1-imy(4)+2*(imy(3)+2*imy(4)))     12d6s21
                     jtest=imy(4)+2*(imy(4)+2*(1-imy(3)+2*(1-imy(4))))   12d6s21
                     itestp=ltest+1                                      12d6s21
                     jtestp=jtest+1                                      12d6s21
                     phsj=1d0                                            12d6s21
                     iuse=-1                                             12d6s21
                     if(iifmx(itestp).ge.0)then                          12d6s21
                      iuse=iifmx(itestp)                                 12d6s21
                     else if(iifmx(jtestp).ge.0)then                     12d6s21
                      iuse=iifmx(jtestp)                                 12d6s21
                      if(isopt(4).ne.0)phsj=-phsj                        12d7s21
                     end if                                              12d6s21
                     if(iuse.ge.0)then                                   12d6s21
                      icolj=igya(4)+irefo(isy(4))*(igya(4)+irefo(isy(4)) 12d6s21
     $                   *(igya(4)+irefo(isy(4))*iuse))                 12d6s21
                      jdenj=id1vnotv(l,isy(4),isy(4))+nnl*icolj          12d6s21
                      ibc(nd1vnotv(l,isy(4),isy(4))+icolj)=1             12d6s21
                      jtmp=itmp                                              12d18s20
                      do iii=0,nl(l)-1                                         12d18s20
                       ii=iff22(iad1+iii)-1                              8d16s21
                       jdj=jdenj+ncsfk(jarg)*ii                             12d18s20
                       do j=0,ncsfk(jarg)-1                                   12d18s20
                        bc(jdj+j)=bc(jdj+j)+bc(jtmp+j)*phsj              12d6s21
                       end do                                                12d18s20
                       jtmp=jtmp+ncsfk(jarg)                             12d7s21
                      end do                                                 12d18s20
                     end if                                              12d6s21
c       - -    -   -
c     (v4|44)=(44|v4)
c
                     ltest=1-imy(4)+2*(imy(4)+2*(imy(3)+2*(1-imy(4))))   12d6s21
                     jtest=imy(4)+2*(1-imy(4)+2*(1-imy(3)+2*imy(4)))     12d6s21
                     itestp=ltest+1                                      12d6s21
                     jtestp=jtest+1                                      12d6s21
                     phsk=-1d0                                            12d6s21
                     iuse=-1                                             12d6s21
                     if(iifmx(itestp).ge.0)then                          12d6s21
                      iuse=iifmx(itestp)                                 12d6s21
                     else if(iifmx(jtestp).ge.0)then                     12d6s21
                      iuse=iifmx(jtestp)                                 12d6s21
                      if(isopt(4).ne.0)phsk=-phsk                        12d6s21
                     end if                                              12d6s21
                     if(iuse.ge.0)then                                   12d6s21
                      icolk=igya(4)+irefo(isy(4))*(igya(4)+irefo(isy(4)) 12d6s21
     $                   *(igya(4)+irefo(isy(4))*iuse))                 12d6s21
                      jdenj=id1vnotv(l,isy(4),isy(4))+nnl*icolk          12d6s21
                      ibc(nd1vnotv(l,isy(4),isy(4))+icolk)=1             12d6s21
                      jtmp=itmp                                              12d18s20
                      do iii=0,nl(l)-1                                         12d18s20
                       ii=iff22(iad1+iii)-1                              8d16s21
                       jdj=jdenj+ncsfk(jarg)*ii                             12d18s20
                       do j=0,ncsfk(jarg)-1                                   12d18s20
                        bc(jdj+j)=bc(jdj+j)+bc(jtmp+j)*phsk              12d6s21
                       end do                                                12d18s20
                       jtmp=jtmp+ncsfk(jarg)                             12d7s21
                      end do                                                 12d18s20
                     end if                                              12d6s21
                    end if                                               12d6s21
                    ibcoff=itmp                                          12d6s21
                   end if                                                12d6s21
                   jmat=jmat+ncsfk(jarg)*ncsfb2(l,iarg)                  12d6s21
                  end do                                                 12d7s21
                 end do
                 ibcoff=ioutg                                            12d6s21
                 iadd=ipack2a(4)                                         12d9s21
                 do i=1,norb                                             12d9s21
                  if(btest(i1o,i).and.btest(j2o,i).and.l2e.eq.0)then     2d17s22
                   itestc=j2c                                            12d6s21
                   itesto=j2o                                            12d6s21
                   nopenk=nopen2p                                        12d6s21
c
c     anihilate common
c
                   if(btest(itestc,i))then                               12d9s21
                    itestc=ibclr(itestc,i)                               12d9s21
                    itesto=ibset(itesto,i)                               12d9s21
                    nopenk=nopenk+1                                             11d13s20
                   else                                                         11d13s20
                    itesto=ibclr(itesto,i)                               12d9s21
                    nopenk=nopenk-1                                             11d13s20
                   end if                                                       11d13s20
c
c     create ket
c
                   if(btest(itesto,nab4(2,1)))then                            12d9s21
                    itestc=ibset(itestc,nab4(2,1))                            12d9s21
                    itesto=ibclr(itesto,nab4(2,1))                            12d9s21
                    nopenk=nopenk-1                                             11d13s20
                   else                                                         11d13s20
                    itesto=ibset(itesto,nab4(2,1))                            12d9s21
                    nopenk=nopenk+1                                             11d13s20
                   end if                                                       11d13s20
                   call gandcr(j2c,j2o,itestc,itesto,nopen2p,nopenk,     12d6s21
     $            norbxx,nnot1,nab1,icode,imap,nx1,irw1,irw2,iwpb1,     12d6s21
     $                  iwpk1,bc,ibc)                                   11d14s22
                   call gandcr(itestc,itesto,i1c,i1o,nopenk,nopen1,      12d6s21
     $            norbxx,nnot2,nab2,icode,imap,nx1,irw1,irw2,iwpb2,     12d6s21
     $               iwpk2,bc,ibc)                                      11d14s22
                   if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                    call spinloop(i2sb,i2smb,i2sk,i2smk,nopen2p,nopen1,  12d6s21
     $               nopenk,ncsfb(iarg),ncsfk(jarg),itype,imatx,ntypeq, 11d22s21
     $                  nab1,iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(iunit),    11d22s21
     $                  ncsfk(jarg),ieoro,bc,ibc)                       11d14s22
                    do ixi=0,ntypeq-1                                     12d6s21
                     ipack=ibc(itype+ixi)                                12d7s21
                     imatu=imatx+ncsfb(iarg)*ncsfk(jarg)*ixi             12d7s21
                     do k=0,ncsfk(jarg)-1                                12d7s21
                      do l=0,ncsfb(iarg)-1                               12d7s21
                       lk=imatu+l+ncsfb(iarg)*k                          12d7s21
                       kl=itrans+k+ncsfk(jarg)*l                         12d7s21
                       bc(kl)=bc(lk)                                     12d7s21
                      end do                                             12d7s21
                     end do                                              12d7s21
                     do j=1,4                                               9d9s21
                      if(ipack2(j).gt.0)then                                9d9s21
                       ipack2a(j)=ipack2(j)                              12d6s21
                       imy(j)=0                                             10d13s21
                      else                                                  9d9s21
                       ipack2a(j)=-ipack2(j)                             12d6s21
                       imy(j)=1                                             10d13s21
                      end if                                                9d9s21
                      if(ipack2a(j).le.norb)then                         12d6s21
                       isy(j)=ism(ipack2a(j))                                9d9s21
                       igya(j)=irel(ipack2a(j))-1                             9d9s21
                      end if                                             12d6s21
                     end do                                                 9d9s21
                     if(ipack2a(3).gt.norb)then                          12d6s21
c
c     (12|v4)
c
                      ltest=imy(1)+2*(imy(2)+2*(imy(3)+2*imy(4)))        12d7s21
                      jtest=1-imy(1)+2*(1-imy(2)+2*(1-imy(3)             12d7s21
     $                    +2*(1-imy(4))))                               12d7s21
                      itestp=ltest+1                                     12d7s21
                      jtestp=jtest+1                                     12d7s21
                      juse=-1                                            12d7s21
                      phsj=1d0                                           12d7s21
                      if(iifmx(itestp).ge.0)then                         12d7s21
                       juse=iifmx(itestp)                                12d7s21
                      else if(iifmx(jtestp).ge.0)then                    12d7s21
                       juse=iifmx(jtestp)                                12d7s21
                       if(isopt(4).ne.0)phsj=-phsj                       12d7s21
                      end if                                             12d7s21
c
c     (14|v2)
c
                      ltest=imy(1)+2*(imy(4)+2*(imy(3)+2*imy(2)))        12d7s21
                      jtest=1-imy(1)+2*(1-imy(4)+2*(1-imy(3)             12d7s21
     $                    +2*(1-imy(2))))                               12d7s21
                      itestp=ltest+1                                     12d7s21
                      jtestp=jtest+1                                     12d7s21
                      kuse=-1                                            12d7s21
                      phsk=-1d0                                           12d7s21
                      if(iifmx(itestp).ge.0)then                         12d7s21
                       kuse=iifmx(itestp)                                12d7s21
                      else if(iifmx(jtestp).ge.0)then                    12d7s21
                       kuse=iifmx(jtestp)                                12d7s21
                       if(isopt(4).ne.0)phsk=-phsk                       12d7s21
                      end if                                             12d7s21
                      if(max(juse,kuse).ge.0)then                        12d7s21
                       jcol=igya(1)+irefo(isy(1))*(igya(2)+irefo(isy(2)) 12d7s21
     $                     *(igya(4)+irefo(isy(4))*juse))               12d7s21
                       kcol=igya(1)+irefo(isy(1))*(igya(4)+irefo(isy(4)) 12d7s21
     $                     *(igya(2)+irefo(isy(2))*kuse))               12d7s21
                       jmat=itrans                                        12d7s21
                       do l=1,4
                        if(nl(l).gt.0)then                                     12d18s20
                         nnl=ncsfk(jarg)*nfdat(2,l,isb)                         12d19s20
                         iad1=jvcv+iff22(jvcv+5+l)                        12d7s21
                         iad2=iad1+nl(l)                                  3d19s21
                         itmp=ibcoff                                           12d18s20
                         ibcoff=itmp+ncsfk(jarg)*nl(l)                          12d18s20
                         call enough('hcdsbk4.  7',bc,ibc)
                         call dgemm('n','n',ncsfk(jarg),nl(l),            3d19s21
     $                      ncsfb2(l,iarg),1d0,bc(jmat),ncsfk(jarg),     3d19s21
     $                      ff22(iad2),ncsfb2(l,iarg),0d0,bc(itmp),      8d16s21
     $                      ncsfk(jarg),                                 3d19s21
     d' hcdsbk4.  2')
                         if(juse.ge.0)then                                12d7s21
                          jden=id1vnotv(l,isy(1),isy(2))+nnl*jcol         12d7s21
                          ibc(nd1vnotv(l,isy(1),isy(2))+jcol)=1           12d7s21
                          jtmp=itmp                                              12d18s20
                          mden=mdhvnotv(l)                                 12d21s20
                          do iii=0,nl(l)-1                                         12d18s20
                           ii=iff22(iad1+iii)-1                           8d16s21
                           jdh=jden+ncsfk(jarg)*ii                             12d18s20
                           ibc(mden+ii)=1                                  12d19s20
                           do j=0,ncsfk(jarg)-1                                   12d18s20
                            bc(jdh+j)=bc(jdh+j)+bc(jtmp+j)*phsj           12d7s21
                           end do                                                12d18s20
                           jtmp=jtmp+ncsfk(jarg)                          12d7s21
                          end do                                                 12d18s20
                         end if                                           12d7s21
                         if(kuse.ge.0)then                                12d7s21
                          jden=id1vnotv(l,isy(1),isy(4))+nnl*kcol         12d7s21
                          ibc(nd1vnotv(l,isy(1),isy(4))+kcol)=1           12d7s21
                          jtmp=itmp                                              12d18s20
                          mden=mdhvnotv(l)                                 12d21s20
                          do iii=0,nl(l)-1                                         12d18s20
                           ii=iff22(iad1+iii)-1                           8d16s21
                           jdh=jden+ncsfk(jarg)*ii                             12d18s20
                           ibc(mden+ii)=1                                  12d19s20
                           do j=0,ncsfk(jarg)-1                                   12d18s20
                            bc(jdh+j)=bc(jdh+j)+bc(jtmp+j)*phsk           12d7s21
                           end do                                                12d18s20
                           jtmp=jtmp+ncsfk(jarg)                          12d7s21
                          end do                                                 12d18s20
                         end if                                           12d7s21
                         ibcoff=itmp                                      12d19s20
                        end if                                              12d19s20
                        jmat=jmat+ncsfk(jarg)*ncsfb2(l,iarg)              3d19s21
                       end do                                             12d19s20
                      end if                                             12d7s21
                     else                                                12d6s21
                      write(6,*)('don''t know how to handle ipack2! '),  12d6s21
     $                    ipack2                                        12d6s21
                      stop 'hcdsbk4'                                     12d6s21
                     end if                                              12d6s21
                    end do                                               12d6s21
                    ibcoff=itype                                         12d6s21
                   end if                                                12d6s21
                  end if                                                 12d7s21
                 end do                                                  12d18s20
                 ibcoff=iunit                                            12d7s21
                else if(l2e.eq.0)then                                   2d7s23
                 nnot=0                                                 2d7s23
                 if(ndifs.eq.4.and.ndifb.eq.4)then                                           10d14s2
                  nnot=4                                                          10d14s22
                  ioxx(1)=1                                                     10d17s22
                  ioxx(2)=1                                                     10d17s22
                  do i=1,norbxx                                                     10d17s22
                   if(btest(gandcb,i))then                                         10d14s22
                    if((btest(i1c,i).and.btest(j2o,i)).or.              2d8s23
     $                 (btest(i1o,i).and..not.btest(j2c,i)))then        2d8s23
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
     $                 ((btest(j2c,i).and..not.btest(i1o,i)).or.                 10d14s22
     $                 (btest(i1c,i).and..not.btest(j2o,i))))then                     10d14s22
                     if(btest(i1c,i))iswap=1                                        10d17s22
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
                   if(btest(i1c,nab4(2,2)).and.                         10d26s22
     $                  .not.btest(i1c,nab4(2,1)))nbt=1                 10d26s22
                  else                                                             10d17s22
                   nbt=0                                                           10d17s22
                   if(btest(j2c,nab4(1,2)).and.                         10d26s22
     $                  .not.btest(j2c,nab4(1,1)))nbt=1                 10d26s22
                  end if                                                           10d17s22
                  if(nbt.ne.0)then                                                 10d17s22
                   nab4(1,1)=nab4(1,2)                                             10d17s22
                   nab4(2,1)=nab4(2,2)                                             10d17s22
                  end if                                                           10d17s22
                 else if(ndifs.eq.0.and.ndifd.eq.2)then                            10d14s22
                  nnot=3
                  do i=1,norbxx
                   if(btest(gandcb,i))then                                         10d14s22
                    if(btest(j2c,i))then
                     nab4(1,1)=i
                     nab4(1,2)=i
                    else                                                           10d14s22
                     nab4(2,1)=i
                     nab4(2,2)=i
                    end if
                   end if                                                          10d14s22
                  end do
                 end if                                                            10d14s22
                 if(nnot.ne.0)then                                      2d7s23
                  iunit=ibcoff                                            12d6s21
                  ibcoff=iunit+ncsfk(jarg)*ncsfk(jarg)                    12d6s21
                  call enough('hcdsbk4.  8',bc,ibc)
                  do iz=iunit,ibcoff-1                                    12d6s21
                   bc(iz)=0d0                                             12d6s21
                  end do                                                  12d6s21
                  do iz=0,ncsfk(jarg)-1                                   12d6s21
                   iad=iunit+iz*(ncsfk(jarg)+1)                           12d6s21
                   bc(iad)=1d0                                            12d6s21
                  end do                                                  12d6s21
                  itrans=ibcoff                                           12d7s21
                  ibcoff=itrans+ncsfb(iarg)*ncsfk(jarg)                   12d7s21
                  call enough('hcdsbk4.  9',bc,ibc)
                  iu1=1
                  iu2=1
                  itestc=j2c                                              12d7s21
                  itesto=j2o                                              12d7s21
                  if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                   itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                   itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopen2p+1                                              11d13s20
                  else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                   itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopen2p-1                                              11d13s20
                  end if                                                        11d13s20
                  if(btest(itesto,nab4(2,iu2)))then                         12d8s20
                   itestc=ibset(itestc,nab4(2,iu2))                         12d8s20
                   itesto=ibclr(itesto,nab4(2,iu2))                             11d27s20
                   nopenk=nopenk-1                                              11d13s20
                  else                                                          11d13s20
                   itesto=ibset(itesto,nab4(2,iu2))                         12d8s20
                   nopenk=nopenk+1                                              11d13s20
                  end if                                                        11d13s20
                  call gandcr(j2c,j2o,itestc,itesto,nopen2p,nopenk,       12d6s21
     $            norbxx,nnot1,nab1,icode,imap,nx1,irw1,irw2,iwpb1,     12d6s21
     $                  iwpk1,bc,ibc)                                   11d14s22
                  call gandcr(itestc,itesto,i1c,i1o,nopenk,nopen1,        12d6s21
     $            norbxx,nnot2,nab2,icode,imap,nx1,irw1,irw2,iwpb2,     12d6s21
     $               iwpk2,bc,ibc)                                      11d14s22
                  if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                   call spinloop(i2sb,i2smb,i2sk,i2smk,nopen2p,nopen1,    12d6s21
     $               nopenk,ncsfb(iarg),ncsfk(jarg),itype,imatx,ntypeq, 11d22s21
     $                  nab1,iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(iunit),    11d22s21
     $                  ncsfk(jarg),ieoro,bc,ibc)                       11d14s22
                   do ixi=0,ntypeq-1                                      12d6s21
                    ipack=ibc(itype+ixi)                                  12d6s21
                    imatu=imatx+ncsfb(iarg)*ncsfk(jarg)*ixi               12d7s21
                    do k=0,ncsfk(jarg)-1                                  12d7s21
                     do l=0,ncsfb(iarg)-1                                 12d7s21
                      lk=imatu+l+ncsfb(iarg)*k                            12d7s21
                      kl=itrans+k+ncsfk(jarg)*l                           12d7s21
                      bc(kl)=bc(lk)                                       12d7s21
                     end do                                               12d7s21
                    end do                                                12d7s21
                    do j=1,4                                               9d9s21
                     if(ipack2(j).gt.0)then                                9d9s21
                      ipack2a(j)=ipack2(j)                                12d6s21
                      imy(j)=0                                             10d13s21
                     else                                                  9d9s21
                      ipack2a(j)=-ipack2(j)                               12d6s21
                      imy(j)=1                                             10d13s21
                     end if                                                9d9s21
                     if(ipack2a(j).le.norb)then                           12d6s21
                      isy(j)=ism(ipack2a(j))                                9d9s21
                      igya(j)=irel(ipack2a(j))-1                             9d9s21
                     end if                                               12d6s21
                    end do                                                 9d9s21
                    if(ipack2a(3).gt.norb)then                            12d6s21
c
c     (12|v4)
c
                     ltest=imy(1)+2*(imy(2)+2*(imy(3)+2*imy(4)))          12d7s21
                     jtest=1-imy(1)+2*(1-imy(2)+2*(1-imy(3)             2d7s23
     $                    +2*(1-imy(4))))                               2d7s23
                     itestp=ltest+1                                       12d7s21
                     jtestp=jtest+1                                       12d7s21
                     juse=-1                                              12d7s21
                     phsj=1d0                                             12d7s21
                     if(iifmx(itestp).ge.0)then                           12d7s21
                      juse=iifmx(itestp)                                  12d7s21
                     else if(iifmx(jtestp).ge.0)then                      12d7s21
                      juse=iifmx(jtestp)                                  12d7s21
                      if(isopt(4).ne.0)phsj=-phsj                         12d7s21
                     end if                                               12d7s21
                     kuse=-1                                              12d7s21
                     if(nnot.eq.4)then                                    12d7s21
c
c     (14|v2)
c
                      ltest=imy(1)+2*(imy(4)+2*(imy(3)+2*imy(2)))         12d7s21
                      jtest=1-imy(1)+2*(1-imy(4)+2*(1-imy(3)              12d7s21
     $                    +2*(1-imy(2))))                                12d7s21
                      itestp=ltest+1                                      12d7s21
                      jtestp=jtest+1                                      12d7s21
                      phsk=-1d0                                           12d7s21
                      if(iifmx(itestp).ge.0)then                          12d7s21
                       kuse=iifmx(itestp)                                 12d7s21
                      else if(iifmx(jtestp).ge.0)then                     12d7s21
                       kuse=iifmx(jtestp)                                 12d7s21
                       if(isopt(4).ne.0)phsk=-phsk                        12d7s21
                      end if                                              12d7s21
                     end if                                               12d7s21
                     if(max(juse,kuse).ge.0)then                          12d7s21
                      jcol=igya(1)+irefo(isy(1))*(igya(2)+irefo(isy(2))   12d7s21
     $                     *(igya(4)+irefo(isy(4))*juse))               12d7s21
                      kcol=igya(1)+irefo(isy(1))*(igya(4)+irefo(isy(4))   12d7s21
     $                     *(igya(2)+irefo(isy(2))*kuse))               12d7s21
                      jmat=itrans                                         12d7s21
                      do l=1,4
                       if(nl(l).gt.0)then                                     12d18s20
                        nnl=ncsfk(jarg)*nfdat(2,l,isb)                         12d19s20
                        iad1=jvcv+iff22(jvcv+5+l)                         12d7s21
                        iad2=iad1+nl(l)                                   12d7s21
                        itmp=ibcoff                                           12d18s20
                        ibcoff=itmp+ncsfk(jarg)*nl(l)                          12d18s20
                        call enough('hcdsbk4. 10',bc,ibc)
                        call dgemm('n','n',ncsfk(jarg),nl(l),             12d7s21
     $                      ncsfb2(l,iarg),1d0,bc(jmat),ncsfk(jarg),     3d19s21
     $                      ff22(iad2),ncsfb2(l,iarg),0d0,bc(itmp),      8d16s21
     $                      ncsfk(jarg),                                 3d19s21
     d' hcdsbk4.  3')
                        if(juse.ge.0)then                                 12d7s21
                         jden=id1vnotv(l,isy(1),isy(2))+nnl*jcol          12d7s21
                         ibc(nd1vnotv(l,isy(1),isy(2))+jcol)=1            12d7s21
                         jtmp=itmp                                              12d18s20
                         mden=mdhvnotv(l)                                 12d21s20
                         do iii=0,nl(l)-1                                         12d18s20
                          ii=iff22(iad1+iii)-1                            12d7s21
                          jdh=jden+ncsfk(jarg)*ii                             12d18s20
                          ibc(mden+ii)=1                                  12d19s20
                          do j=0,ncsfk(jarg)-1                                   12d18s20
                           bc(jdh+j)=bc(jdh+j)+bc(jtmp+j)*phsj            12d7s21
                          end do                                                12d18s20
                          jtmp=jtmp+ncsfk(jarg)                           12d7s21
                         end do                                                 12d18s20
                        end if                                            12d7s21
                        if(kuse.ge.0)then                                 12d7s21
                         jden=id1vnotv(l,isy(1),isy(4))+nnl*kcol          12d7s21
                         ibc(nd1vnotv(l,isy(1),isy(4))+kcol)=1            12d7s21
                         jtmp=itmp                                              12d18s20
                         mden=mdhvnotv(l)                                 12d21s20
                         do iii=0,nl(l)-1                                         12d18s20
                          ii=iff22(iad1+iii)-1                            12d7s21
                          jdh=jden+ncsfk(jarg)*ii                             12d18s20
                          ibc(mden+ii)=1                                  12d19s20
                          do j=0,ncsfk(jarg)-1                                   12d18s20
                           bc(jdh+j)=bc(jdh+j)+bc(jtmp+j)*phsk            12d7s21
                          end do                                                12d18s20
                          jtmp=jtmp+ncsfk(jarg)                           12d7s21
                         end do                                                 12d18s20
                        end if                                            12d7s21
                        ibcoff=itmp                                       12d7s21
                       end if                                              12d19s20
                       jmat=jmat+ncsfk(jarg)*ncsfb2(l,iarg)               12d7s21
                      end do                                              12d7s21
                     end if                                               12d7s21
                    else                                                  12d7s21
                     write(6,*)('don''t know how to handle ipack2! '),    12d7s21
     $                    ipack2                                        12d6s21
                     stop 'hcdsbk4'                                       12d7s21
                    end if                                                12d7s21
                   end do                                                 12d7s21
                   ibcoff=itype                                           12d7s21
                  end if                                                  12d7s21
                 end if                                                 2d7s23
                end if                                                   12d18s20
               end if                                                   2d7s23
               j2o=ibclr(j2o,norbx)                                      12d18s20
               j2o=ibset(j2o,norbxxx)                                     12d18s20
               gandcc=ieor(i1c,j2c)                                     12d21s22
               gandco=ieor(i1o,j2o)                                     12d21s22
               gandcb=ior(gandcc,gandco)                                12d21s22
               ndifb=popcnt(gandcb)                                     12d21s22
               if(ndifb.le.4.and.l2e.eq.0)then                          2d7s23
                ndifd=popcnt(gandcc)                                     10d13s22
                ndifs=popcnt(gandco)                                     10d13s22
                nnot=0                                                  2d7s23
                if(ndifs.eq.4.and.ndifb.eq.4)then                                           10d14s2
                 nnot=4                                                          10d14s22
                 ioxx(1)=1                                                     10d17s22
                 ioxx(2)=1                                                     10d17s22
                 do i=1,norbxxx                                                     10d17s22
                  if(btest(gandcb,i))then                                         10d14s22
                   if((btest(i1c,i).and.btest(j2o,i)).or.               2d7s23
     $                 (btest(i1o,i).and..not.btest(j2c,i)))then        2d7s23
                    nab4(2,ioxx(2))=i                                       10d17s22
                    ioxx(2)=ioxx(2)+1                                       10d17s22
                   else                                                           10d14s22
                    nab4(1,ioxx(1))=i                                       10d17s22
                    ioxx(1)=ioxx(1)+1                                       10d17s22
                   end if                                                         10d14s22
                  end if
                 end do
                 if(nnot.eq.4)then                                        2d7s23
                  iunit=ibcoff                                            12d6s21
                  ibcoff=iunit+ncsfk(jarg)*ncsfk(jarg)                    12d6s21
                  call enough('hcdsbk4. 11',bc,ibc)
                  do iz=iunit,ibcoff-1                                    12d6s21
                   bc(iz)=0d0                                             12d6s21
                  end do                                                  12d6s21
                  do iz=0,ncsfk(jarg)-1                                   12d6s21
                   iad=iunit+iz*(ncsfk(jarg)+1)                           12d6s21
                   bc(iad)=1d0                                            12d6s21
                  end do                                                  12d6s21
                  itrans=ibcoff                                           12d7s21
                  ibcoff=itrans+ncsfb(iarg)*ncsfk(jarg)                   12d7s21
                  call enough('hcdsbk4. 12',bc,ibc)
                  iu1=1
                  iu2=1
                  itestc=j2c                                              12d8s20
                  itesto=j2o                                              12d8s20
                  if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                   itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                   itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopen2p+1                                              11d13s20
                  else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                   itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopen2p-1                                              11d13s20
                  end if                                                        11d13s20
                  if(btest(itesto,nab4(2,iu2)))then                         12d8s20
                   itestc=ibset(itestc,nab4(2,iu2))                         12d8s20
                   itesto=ibclr(itesto,nab4(2,iu2))                             11d27s20
                   nopenk=nopenk-1                                              11d13s20
                  else                                                          11d13s20
                   itesto=ibset(itesto,nab4(2,iu2))                         12d8s20
                   nopenk=nopenk+1                                              11d13s20
                  end if                                                        11d13s20
                  call gandcr(j2c,j2o,itestc,itesto,nopen2p,nopenk,       12d7s21
     $              norbxxx,nnot1,nab1,icode,imap,nx1,irw1,irw2,iwpb1,     12d6s21
     $                  iwpk1,bc,ibc)                                   11d14s22
                  call gandcr(itestc,itesto,i1c,i1o,nopenk,nopen1,        12d6s21
     $            norbxxx,nnot2,nab2,icode,imap,nx1,irw1,irw2,iwpb2,     12d6s21
     $               iwpk2,bc,ibc)                                      11d14s22
                  if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                   call spinloop(i2sb,i2smb,i2sk,i2smk,nopen2p,nopen1,    12d6s21
     $               nopenk,ncsfb(iarg),ncsfk(jarg),itype,imatx,ntypeq, 11d22s21
     $                  nab1,iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(iunit),    11d22s21
     $                  ncsfk(jarg),ieoro,bc,ibc)                       11d14s22
                   do ixi=0,ntypeq-1                                       12d6s21
                    ipack=ibc(itype+ixi)                                   12d6s21
                    imatu=imatx+ncsfb(iarg)*ncsfk(jarg)*ixi                12d7s21
                    do k=0,ncsfk(jarg)-1                                  12d7s21
                     do l=0,ncsfb(iarg)-1                                 12d7s21
                      lk=imatu+l+ncsfb(iarg)*k                            12d7s21
                      kl=itrans+k+ncsfk(jarg)*l                           12d7s21
                      bc(kl)=bc(lk)                                       12d7s21
                     end do                                               12d7s21
                    end do                                                12d7s21
                    do j=1,4                                               9d9s21
                     if(ipack2(j).gt.0)then                                9d9s21
                      ipack2a(j)=ipack2(j)                                12d6s21
                      imy(j)=0                                              10d13s21
                     else                                                  9d9s21
                      ipack2a(j)=-ipack2(j)                               12d6s21
                      imy(j)=1                                             10d13s21
                     end if                                                9d9s21
                     if(ipack2a(j).le.norb)then                           12d6s21
                      isy(j)=ism(ipack2a(j))                                9d9s21
                      igya(j)=irel(ipack2a(j))-1                             9d9s21
                     end if                                               12d6s21
                    end do                                                 9d9s21
                    if(ipack2a(2).le.norb)then                            12d7s21
                     if(isy(2).ne.lgoal3)then
                      write(6,*)('symmetries do not match!!! '),isy(2),
     $                   lgoal3
                      stop 'hcdsbk4'
                     end if
c
c but 3x is stored badc, indexed b,c,d
c                a b c d     c b a d = bcda
c          _     _               _
c      d   d s   d s d       d s d
c     (v'2|v"v)=(v"v|v'2) & (v'v|v"2)
c
                     ltest=imy(3)+2*(imy(4)+2*(imy(1)+2*imy(2)))          12d7s21
                     jtest=1-imy(3)+2*(1-imy(4)+2*(1-imy(1)             2d7s23
     $                    +2*(1-imy(2))))                               2d7s23
                     itestp=ltest+1                                       12d7s21
                     jtestp=jtest+1                                       12d7s21
                     iuse=-1                                              12d7s21
                     phs=1d0                                              12d7s21
                     if(iifmx(itestp).ge.0)then                           12d7s21
                      iuse=iifmx(itestp)                                  12d7s21
                     else if(iifmx(jtestp).ge.0)then                      12d7s21
                      iuse=iifmx(jtestp)                                  12d7s21
                      if(isopt(4).ne.0)phs=-phs                           12d7s21
                     end if                                               12d7s21
                     kuse=-1                                              12d9s21
                     phsk=-1d0                                            12d9s21
                     ltest=imy(1)+2*(imy(4)+2*(imy(3)+2*imy(2)))          12d7s21
                     jtest=1-imy(1)+2*(1-imy(4)+2*(1-imy(3)             2d7s23
     $                    +2*(1-imy(2))))                               2d7s23
                     itestp=ltest+1                                       12d7s21
                     jtestp=jtest+1                                       12d7s21
                     if(iifmx(itestp).ge.0)then                           12d7s21
                      kuse=iifmx(itestp)                                  12d7s21
                     else if(iifmx(jtestp).ge.0)then                      12d7s21
                      kuse=iifmx(jtestp)                                  12d7s21
                      if(isopt(4).ne.0)phsk=-phsk                         12d9s21
                     end if                                               12d7s21
                     if(iuse.ge.0)then                                    12d7s21
                      icol=igya(2)+irefo(isy(2))*iuse                     12d7s21
                      i3xk(1,iuse+1)=kuse                                 12d9s21
                      if(phsk*phs.gt.0d0)then                             12d9s21
                       i3xk(2,iuse+1)=1                                   12d9s21
                      else                                                12d9s21
                       i3xk(2,iuse+1)=2                                   12d9s21
                      end if                                              12d9s21
                      jmat=itrans                                         12d7s21
                      do l=1,4                                               12d19s20
                       if(nl(l).gt.0)then                                     12d18s20
                        nnl=ncsfk(jarg)*nfdat(2,l,isb)                         12d19s20
                        iad1=jvcv+iff22(jvcv+5+l)                         8d16s21
                        iad2=iad1+nl(l)                                    3d19s21
                        itmp=ibcoff                                           12d18s20
                        ibcoff=itmp+ncsfk(jarg)*nl(l)                          12d18s20
                        call enough('hcdsbk4. 13',bc,ibc)
                        call dgemm('n','n',ncsfk(jarg),nl(l),             12d7s21
     $                     ncsfb2(l,iarg),                              12d7s21
     $            1d0,bc(jmat),ncsfk(jarg),ff22(iad2),ncsfb2(l,iarg),   12d7s21
     $                     0d0,bc(itmp),ncsfk(jarg),                    12d7s21
     d' hcdsbk4.  4')
                        jden=id3vnotv3(l)+nnl*icol                        12d7s21
                        ibc(nd3vnotv3(l)+icol)=1                          12d22s20
                        mden=md3vnotv3(l)                                 12d7s21
                        jtmp=itmp                                              12d18s20
                        do iii=0,nl(l)-1                                         12d18s20
                         ii=iff22(iad1+iii)-1                             12d7s21
                         jdh=jden+ncsfk(jarg)*ii                              12d18s20
                         ibc(mden+ii)=1                                      12d19s20
                         do j=0,ncsfk(jarg)-1                                   12d18s20
                          bc(jdh+j)=bc(jdh+j)+bc(jtmp+j)*phs                       12d18s20
                         end do                                                12d18s20
                         jtmp=jtmp+ncsfk(jarg)                                  12d18s20
                        end do                                                 12d18s20
                        ibcoff=itmp                                          12d19s20
                       end if                                                12d19s20
                       jmat=jmat+ncsfk(jarg)*ncsfb2(l,iarg)               12d7s21
                      end do                                               12d7s21
                     end if                                               12d7s21
                    else                                                  12d7s21
                     write(6,*)('I don''t know how to deal with '),
     $                    ipack2
                     stop 'hcdsbk4'
                    end if                                                12d7s21
                   end do                                                 12d7s21
                   ibcoff=iunit                                           12d7s21
                  end if                                                  12d7s21
                 end if                                                 2d7s23
                end if                                                  2d7s23
               end if                                                   8d15s21
              end if                                                    3d2s21
              jvcv=jvcv+nspace                                          3d19s21
             end do                                                     12d18s20
c
             if(njhere.gt.0)then                                        3d2s21
              do l=1,4
               nok4f(l)=0                                               12d21s20
               nok4(l)=0                                                11d17s22
               do isc=1,nsymb                                           11d17s22
                do isd=1,nsymb                                          11d17s22
                 nokdc(isd,isc,l)=0                                     12d21s20
                end do                                                  11d17s22
               end do                                                   11d17s22
               if(nfdat(2,l,isb).gt.0)then                               12d19s20
                do i=0,nfdat(2,l,isb)-1                                  12d21s20
                 if(ibc(mdhvnotv(l)+i).ne.0)then                         12d21s20
                  ibc(mdhvnotv(l)+nok4f(l))=i                            12d21s20
                  nok4f(l)=nok4f(l)+1                                    12d21s20
                 end if                                                  12d21s20
                end do                                                   12d21s20
                if(nok4f(l).gt.0)then                                    12d21s20
                 nnl=ncsfk(jarg)*nfdat(2,l,isb)                         12d7s21
                 nok4(l)=0                                               12d21s20
                 jdhvnotv=idhvnotv(l)                                    12d21s20
                 jdhvnotvf=idhvnotvf(l)                                  1d6s21
                 do i=0,irefo(lgoal)-1                                  12d8s21
                  if(ibc(ndhvnotv(l)+i).ne.0)then                        12d21s20
                   ibc(ndhvnotv(l)+nok4(l))=i                            12d21s20
                   nok4(l)=nok4(l)+1                                     12d21s20
                   iad=idhvnotv(l)+nnl*i
                   do k=0,nok4f(l)-1                                     12d21s20
                    iad=idhvnotv(l)+ncsfk(jarg)*ibc(mdhvnotv(l)+k)+nnl*i12d7s21
                    do j=0,ncsfk(jarg)-1                                12d7s21
                     bc(jdhvnotvf+j)=bc(iad+j)                           1d6s21
                    end do                                               12d21s20
                    jdhvnotvf=jdhvnotvf+ncsfk(jarg)                     12d7s21
                    iad=iad+i11s-1                                       1d6s21
                    do j=0,njhere-1                                      12d21s20
                     bc(jdhvnotv+j)=bc(iad+j)                            12d21s20
                    end do                                               12d21s20
                    jdhvnotv=jdhvnotv+njhere                             12d21s20
                   end do                                                12d21s20
   10              format(i5,5x,20i2)
                  end if
                 end do
                 do isc=1,nsymb                                          12d21s20
                  iscv=multh(isc,lgoal)                                   12d21s20
                  do isd=1,nsymb                                        8d16s21
                   iscdv=multh(iscv,isd)                                     12d18s20
                   nn=irefo(isd)*irefo(isc)                                 12d18s20
                   nnn=nn*irefo(iscdv)*ntype                            12d7s21
                   nokdc(isd,isc,l)=0                                     12d21s20
                   jd1vnotv=id1vnotv(l,isd,isc)                           12d21s20
                   do i=0,nnn-1
                    if(ibc(nd1vnotv(l,isd,isc)+i).ne.0)then
                     icol=i                                             12d29s20
                     ibc(nd1vnotv(l,isd,isc)+nokdc(isd,isc,l))=icol      12d29s20
                     nokdc(isd,isc,l)=nokdc(isd,isc,l)+1                 12d21s20
                     sz=0d0
                     do k=0,nok4f(l)-1                                    12d21s20
                      iad=id1vnotv(l,isd,isc)+i11s-1                      12d21s20
     $                   +ncsfk(jarg)*ibc(mdhvnotv(l)+k)+nnl*i          12d7s21
                      do j=0,njhere-1                                     12d21s20
                       bc(jd1vnotv+j)=bc(iad+j)                           12d21s20
                      end do                                              12d21s20
                      jd1vnotv=jd1vnotv+njhere                            12d21s20
                     end do                                               12d21s20
                    end if
                   end do
                  end do
                 end do                                                  12d21s20
                end if
               end if                                                    12d21s20
              end do                                                     12d21s20
             end if                                                     3d2s21
             if(itransvs.eq.0)then                                      1d29s21
              call ddi_done(ibc(ircv),nrcv)                                  1d29s21
              itransvs=1                                                1d29s21
              itmp=ibcoff                                                    1d27s21
              ibcoff=itmp+nggk                                          8d13s21
              call enough('hcdsbk4. 14',bc,ibc)
              jjvs=ivs                                                        1d27s21
              nggm=nggk-1                                               8d13s21
              do i=0,nggg-1                                                  1d27s21
               jvs0=jjvs                                                      1d27s21
               do iv=0,nvirt(ksbv)-1                                    8d13s21
                iad=itmp+iv                                                  1d27s21
                do ir=0,nrootm                                               1d27s21
                 bc(iad+ir*nvirt(ksbv))=bc(jjvs+ir)                     8d13s21
                end do                                                       1d27s21
                jjvs=jjvs+nrootu                                               1d27s21
               end do                                                        1d27s21
               do j=0,nggm                                                   1d27s21
                bc(jvs0+j)=bc(itmp+j)                                        1d27s21
               end do                                                        1d27s21
              end do                                                         1d27s21
              ibcoff=itmp                                                    1d27s21
             end if                                                     1d29s21
             if(njhere.gt.0)then                                        8d15s21
              ipass=1                                                   3d2s21
              jdkeep=idkeep(ipass)                                      1d4s21
              nok=0                                                     1d4s21
              do i=0,irefo(lgoal3)*ntype-1                              12d7s21
               do l=1,4                                                 1d4s21
                if(nfdat(2,l,isb).gt.0)then                             1d4s21
                 if(ibc(nd3vnotv3(l)+i).ne.0)then                       12d7s21
                  do k=0,nfdat(2,l,isb)-1                               1d4s21
                   if(ibc(md3vnotv3(l)+k).ne.0)then                     12d7s21
                    ibc(ndkeep(ipass)+nok)=i                            1d5s21
                    ibc(mdkeep(ipass)+nok)=k+loff(l)                    1d5s21
                    iad=id3vnotv3(l)+ncsfk(jarg)*(k                     12d7s21
     $                   +nfdat(2,l,isb)*i)                             1d4s21
     $                   +i11s-1                                        3d2s21
                    ibc(ibcoff+nok)=l
                    do j=0,njhere-1                                     3d2s21
                     bc(jdkeep+j)=bc(iad+j)                             1d4s21
                    end do                                              1d4s21
                    jdkeep=jdkeep+njhere                                3d2s21
                    nok=nok+1                                           1d4s21
                   end if                                               1d4s21
                  end do                                                1d4s21
                 end if                                                 1d4s21
                end if                                                  1d4s21
               end do                                                   1d4s21
              end do                                                     1d4s21
              keep(ipass)=nok                                           1d4s21
              if(nok.gt.0)then                                          3d2s21
               nnk=nvirt(ksbv)*nrootu                                    8d13s21
               itmpk=ibcoff                                             2d3s21
               itmpd=itmpk+nnk*nok                                       2d3s21
               ibcoff=itmpd+njhere*nok                                  3d2s21
               call enough('hcdsbk4. 15',bc,ibc)
               jvss=jvs+nnk*(i11s-1)                                     3d2s21
               call dgemm('n','n',nnk,nok,njhere,1d0,                   8d13s21
     $               bc(jvss),nnk,bc(idkeep(ipass)),njhere,0d0,         8d13s21
     $               bc(itmpk),nnk,                                     8d13s21
     d' hcdsbk4.  5')
               jtmpk=itmpk                                              2d3s21
               do i=0,nok-1                                             2d3s21
                do j=0,njhere-1                                         3d2s21
                 ji=idkeep(ipass)+j+njhere*i                            3d2s21
                 ij=itmpd+i+nok*j                                       2d3s21
                 bc(ij)=bc(ji)                                          2d3s21
                end do                                                  2d3s21
                in=ibc(ndkeep(ipass)+i)                                 2d3s21
                k=ibc(mdkeep(ipass)+i)                                  2d3s21
                jvmat=ivmat(isb)+nnk*(k+nfh*in)                         8d13s21
                do j=0,nnk-1                                            8d13s21
                 bc(jvmat+j)=bc(jvmat+j)+bc(jtmpk+j)                    2d4s21
                end do                                                  8d13s21
                jtmpk=jtmpk+nnk                                         8d13s21
               end do                                                   2d3s21
               ibcoff=itmpk                                             8d13s21
              end if                                                    2d3s21
              do isbv1=1,nsymb                                           12d21s20
               isbv2=multh(isbv1,isbv12)                                 12d21s20
               if(isbv1.le.isbv2)then                                    12d21s20
c     S is ket, D is bra
                if(isbv12.eq.1.and.ksbv.eq.isbv1)then                   8d13s21
                 nokf=nok4f(1)                                           2d25s21
                 nrow=njhere*nokf                                          12d21s20
                 intden=ibcoff                                             12d21s20
                 ibcoff=intden+nrow*nvirt(ksbv)                         8d13s21
                 call enough('hcdsbk4. 16',bc,ibc)
                 fact=0d0                                                  12d21s20
                 if(min(nok4(1),nok4f(1)).gt.0)then                      2d25s21
                  nok=nok4(1)                                            2d25s21
                  nokf=nok4f(1)                                          2d25s21
                  itmp=ibcoff                                              12d21s20
                  ibcoff=itmp+nok*nvirt(ksbv)                           8d13s21
                  call enough('hcdsbk4. 17',bc,ibc)
                  jtmp=itmp                                                12d21s20
                  iosym=multh(ksbv,isopt(1))                            12d7s21
                  if(ih0n(ksbv).gt.0)then                               12d7s21
                   do iv=0,nvirt(ksbv)-1                                 8d13s21
                    ivp=iv+irefo(ksbv)                                   12d7s21
                    do i=0,nok-1                                            12d21s20
                     iadh=ih0n(ksbv)+ivp+nh0(ksbv)*ibc(ndhvnotv(1)+i)   12d7s21
                     bc(jtmp+i)=bc(iadh)
                    end do                                                  12d21s20
                    jtmp=jtmp+nok                                           12d21s20
                   end do                                                   12d21s20
                  else if(ih0n(iosym).gt.0)then                         12d7s21
                   phsh=1d0                                             12d7s21
                   do iv=0,nvirt(ksbv)-1                                 8d13s21
                    ivp=iv+irefo(ksbv)                                   12d7s21
                    do i=0,nok-1                                            12d21s20
                     iadh=ih0n(iosym)+ibc(ndhvnotv(1)+i)                12d7s21
     $                    +nh0(iosym)*ivp                               12d7s21
                     bc(jtmp+i)=bc(iadh)*phsh                           12d7s21
                    end do                                                  12d21s20
                    jtmp=jtmp+nok                                           12d21s20
                   end do                                                   12d21s20
                  end if                                                12d7s21
                  call dgemm('n','n',nrow,nvirt(ksbv),nok,sr2,          8d13s21
     $              bc(idhvnotv(1)),nrow,bc(itmp),nok,fact,bc(intden),  2d25s21
     $                nrow,                                             2d25s21
     d' hcdsbk4.  6')
                  fact=1d0                                                 12d21s20
                  ibcoff=itmp                                              12d21s20
                 end if                                                    12d21s20
                 nok=0                                                     12d21s20
                 do isc=1,nsymb
                  iscv=multh(isc,lgoal)                                 12d7s21
                  do isd=1,nsymb                                        8d16s21
                   iscdv=multh(iscv,isd)                                    12d19s20
                   nn=irefo(isd)*irefo(isc)                                 12d18s20
                   nnn=nn*irefo(iscdv)                                       12d18s20
                   if(min(nokdc(isd,isc,1),nok4f(1)).gt.0)then          12d7s21
                    nok=nokdc(isd,isc,1)                                 2d25s21
                    iprod=multh(multh(isd,isc),multh(iscdv,ksbv))
                    nokf=nok4f(1)                                        2d25s21
                    ncol=nok*nokf                                            12d21s20
c     (d|c)(scdv|ksbv)
                    itmp=ibcoff                                            12d21s20
                    ibcoff=itmp+nok*nvirt(ksbv)                            12d21s20
                    call enough('hcdsbk4. 18',bc,ibc)
                    do iz=itmp,ibcoff-1                                 8d15s21
                     bc(iz)=0d0                                         8d15s21
                    end do                                              8d15s21
                    irdirc=irefo(isd)*irefo(isc)                        8d15s21
                    do i=0,nok-1                                        12d7s21
                     jtmp=itmp+i                                        12d7s21
                     idt=ibc(nd1vnotv(1,isd,isc)+i)/nnn                 12d7s21
                     left=ibc(nd1vnotv(1,isd,isc)+i)-idt*nnn            12d7s21
                     idv=left/irdirc                                    12d7s21
                     idc=left-irdirc*idv                                12d7s21
                     ic=idc/irefo(isd)                                  12d7s21
                     id=idc-irefo(isd)*ic                               12d7s21
                     idd=ionex(isc,isd,iscdv)+ic+irefo(isc)*(id+        12d8s21
     $                    irefo(isd)*(idv+irefo(iscdv)*nvirt(ksbv)*idt))12d8s21
                     do kv=0,nvirt(ksbv)-1                              8d15s21
                      bc(jtmp+kv*nok)=bc(jtmp+kv*nok)+bc(idd+nnn*kv)    12d7s21
                     end do                                             8d15s21
                    end do                                              8d15s21
                    call dgemm('n','n',nrow,nvirt(ksbv),nok,sr2,          2d25s21
     $           bc(id1vnotv(1,isd,isc)),nrow,bc(itmp),nok,fact,        2d25s21
     $                bc(intden),nrow,                                  12d21s20
     d' hcdsbk4.  7')
                    fact=1d0                                                 12d21s20
                    ibcoff=itmp                                              12d21s20
                   end if                                               8d15s21
                  end do                                                   12d21s20
                 end do                                                    12d21s20
                 if(fact.gt.0.5d0)then                                     12d21s20
                  itrans=ibcoff                                            12d21s20
                  ibcoff=itrans+nvirt(ksbv)*nrow                           12d21s20
                  call enough('hcdsbk4. 19',bc,ibc)
                  do i=0,nvirt(ksbv)-1                                     12d21s20
                   do j=0,nrow-1                                           12d21s20
                    ji=intden+j+nrow*i                                   12d29s20
                    ij=itrans+i+nvirt(ksbv)*j                              12d21s20
                    bc(ij)=bc(ji)                                          12d21s20
                   end do                                                  12d21s20
                  end do                                                   12d21s20
                  do if=0,nokf-1                                           12d21s20
                   do ir=0,nrootu-1                                        12d21s20
                    iadvd=joffdnon+nvirt(ksbv)*(ir                       12d21s20
     $                   +nrootu*ibc(mdhvnotv(1)+if))                    2d25s21
                    do j=0,njhere-1                                        12d21s20
                     iad=itrans+nvirt(ksbv)*(j+njhere*if)                  12d21s20
                     iads=jvs+nvirt(ksbv)*(ir+nrootu*(i11s+j-1))         1d28s21
                     do iv=0,nvirt(ksbv)-1                                  12d21s20
                      gd(iadvd+iv)=gd(iadvd+iv)+bc(iad+iv)*bc(iads+iv)   1d27s21
                     end do                                                12d21s20
                    end do                                                 12d21s20
                   end do                                                  12d21s20
                  end do                                                   12d21s20
                 end if                                                    12d21s20
                 ibcoff=intden                                             12d21s20
                end if                                                     12d21s20
                if(isbv12.eq.1)then                                      12d22s20
                 joffdnon=joffdnon+nfdat(2,1,isb)*nvirt(isbv1)*nrootu     12d21s20
                 nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                    12d21s20
                 isw=0                                                   12d21s20
                else                                                      12d21s20
                 nvv=nvirt(isbv1)*nvirt(isbv2)                            12d21s20
                 isw=1                                                   12d21s20
                end if                                                    12d21s20
                do l=1,4                                                 12d21s20
                 if(nfdat(2,l,isb).gt.0)then                             12d21s20
                  if(ksbv.eq.isbv1.or.ksbv.eq.isbv2)then                   12d22s20
                   if(isbv1.eq.ksbv)then                                    12d21s20
                    isbvu=isbv2                                             12d21s20
                    tf=1d0                                                  12d21s20
                   else                                                     12d21s20
                    isbvu=isbv1                                             12d21s20
                    tf=-1d0                                                 12d21s20
                   end if                                                   12d21s20
                   if(l.eq.1)then                                         12d21s20
                    factt=1d0                                             12d21s20
                    tf2=1d0
                   else                                                   12d21s20
                    tf2=-1d0                                             12d29s20
                    factt=tf                                              12d21s20
                   end if                                                 12d21s20
                   if(nok4f(l).gt.0.and.nvirt(isbvu).gt.0)then          3d2s21
                    nrow=njhere*nok4f(l)                                  12d21s20
                    intden=ibcoff                                         12d21s20
                    ibcoff=intden+nrow*nvirt(isbvu)                       12d21s20
                    call enough('hcdsbk4. 20',bc,ibc)
                    fact=0d0                                              12d21s20
                    if(nok4(l).gt.0)then                                 1d1s21
c     Gr jv j = D j k i H0iv2 V jv iv2 r k
                     iosym=multh(isbvu,isopt(1))                        12d7s21
                     itmp=ibcoff                                          12d21s20
                     ibcoff=itmp+nok4(l)*nvirt(isbvu)                     12d21s20
                     call enough('hcdsbk4. 21',bc,ibc)
                     jtmp=itmp                                            12d21s20
                     if(ih0n(isbvu).gt.0)then                            12d7s21
                      do iv=0,nvirt(isbvu)-1                               12d21s20
                       ivp=iv+irefo(isbvu)                              12d7s21
                       ih0=ih0n(isbvu)+ivp
                       do i=0,nok4(l)-1                                 12d7s21
                        ii=ibc(ndhvnotv(l)+i)                              12d22s20
                        bc(jtmp+i)=bc(ih0+nh0(isbvu)*ii)
                       end do                                           12d7s21
                       jtmp=jtmp+nok4(l)                                   12d21s20
                      end do                                            12d7s21
                     else if(ih0n(iosym).gt.0)then                      12d7s21
                      phsh=1d0                                          12d7s21
                      do iv=0,nvirt(isbvu)-1                               12d21s20
                       ivp=iv+irefo(isbvu)                              12d7s21
                       ih0=ih0n(iosym)+nh0(iosym)*ivp                   12d7s21
                       do i=0,nok4(l)-1                                 12d7s21
                        ii=ibc(ndhvnotv(l)+i)                              12d22s20
                        bc(jtmp+i)=bc(ih0+ii)*phsh                      12d7s21
                       end do                                           12d7s21
                       jtmp=jtmp+nok4(l)                                   12d21s20
                      end do                                            12d7s21
                     end if                                             12d7s21
                     call dgemm('n','n',nrow,nvirt(isbvu),nok4(l),factt,  12d21s20
     $                  bc(idhvnotv(l)),nrow,bc(itmp),nok4(l),fact,     12d23s20
     $                  bc(intden),nrow,                                12d21s20
     d' hcdsbk4.  8')
                     ibcoff=itmp                                          12d21s20
                     fact=1d0                                             12d21s20
                    end if                                                12d21s20
                    do isc=1,nsymb                                        12d22s20
                     iscv=multh(isc,lgoal)                                12d22s20
                     do isd=1,nsymb                                     12d7s21
                      irdrc=irefo(isd)*irefo(isc)                       12d7s21
                      iscdv=multh(iscv,isd)                               12d22s20
                      if(nokdc(isd,isc,l).gt.0)then                       12d22s20
                       nn=irefo(isc)*irefo(isd)                          12d22s20
                       nnn=nn*irefo(iscdv)                              12d7s21
                       itmp=ibcoff                                        12d22s20
                       ibcoff=itmp+nokdc(isd,isc,l)*nvirt(isbvu)          12d22s20
                       call enough('hcdsbk4. 22',bc,ibc)
                       do i=0,nokdc(isd,isc,l)-1                        12d7s21
                        idt=ibc(nd1vnotv(l,isd,isc)+i)/nnn              12d7s21
                        left=ibc(nd1vnotv(l,isd,isc)+i)-nnn*idt         12d7s21
                        idv=left/irdrc                                  12d7s21
                        idc=left-irdrc*idv                              12d7s21
                        ic=idc/irefo(isd)                               12d7s21
                        id=idc-irefo(isd)*ic                            12d7s21
                        idvk=ionex(isc,isd,iscdv)+ic+irefo(isc)*(id     12d8s21
     $            +irefo(isd)*(idv+irefo(iscdv)*nvirt(isbvu)*idt))      12d8s21
                        jtmp=itmp+i                                     12d8s21
                        do iv=0,nvirt(isbvu)-1                          12d7s21
                         bc(jtmp+iv*nokdc(isd,isc,l))=bc(idvk+iv*nnn)   12d7s21
                        end do                                          8d16s21
                       end do                                           8d16s21
                       call dgemm('n','n',nrow,nvirt(isbvu),
     $                   nokdc(isd,isc,l),factt,bc(id1vnotv(l,isd,isc)),12d22s20
     $                    nrow,bc(itmp),nokdc(isd,isc,l),fact,          12d22s20
     $                    bc(intden),nrow,                              12d22s20
     d' hcdsbk4.  9')
                       fact=1d0                                           12d22s20
                       ibcoff=itmp                                        12d22s20
                      end if                                              12d22s20
                     end do                                               12d22s20
                    end do                                                12d22s20
                    if(fact.gt.0.5d0)then                                 12d22s20
                     itrans=ibcoff                                        12d22s20
                     ibcoff=itrans+nrow*nvirt(isbvu)                      12d22s20
                     call enough('hcdsbk4. 23',bc,ibc)
                     if(isbv1.eq.isbv2)then                              12d29s20
                      do i=0,nvirt(isbvu)-1                                12d22s20
                       do j=0,nrow-1                                       12d22s20
                        ji=intden+j+nrow*i                                 12d22s20
c     vjk
                        ij=itrans+i+nvirt(isbvu)*j                         12d22s20
                        bc(ij)=bc(ji)                                      12d22s20
                       end do                                              12d22s20
                      end do                                               12d22s20
                     else                                                12d29s20
                      do iv=0,nvirt(isbvu)-1                             12d29s20
                       do k=0,nok4f(l)-1                                  12d29s20
                        do j=0,njhere-1                                   12d29s20
                         jki=intden+j+njhere*(k+nok4f(l)*iv)             12d29s20
                         kij=itrans+k+nok4f(l)*(iv+nvirt(isbvu)*j)       12d29s20
c     kvj
                         bc(kij)=bc(jki)                                 12d29s20
                        end do                                           12d29s20
                       end do                                            12d29s20
                      end do                                              12d29s20
                     end if                                              12d29s20
                     ncol=nrootu*nvirt(ksbv)                              12d22s20
                     itmpsv=ibcoff                                         12d22s20
                     ibcoff=itmpsv+njhere*ncol                          8d26s21
                     call enough('hcdsbk4. 24',bc,ibc)
                     if(isbv1.ne.isbv2)then                               12d22s20
                      do j=0,njhere-1                                      12d22s20
                       do i=0,ncol-1                                       12d22s20
                        ij=jvs+i+ncol*(j+i11s-1)                           12d22s20
                        ji=itmpsv+j+njhere*i                               12d22s20
                        bc(ji)=bc(ij)                                      12d22s20
                       end do                                              12d22s20
                      end do                                               12d22s20
                     end if                                               12d22s20
                     if(isbv1.eq.isbv2)then                               12d22s20
                      do ir=0,nrootm                                     1d27s21
                       do j=0,njhere-1                                     12d22s20
                        iad=jvs+nvirt(ksbv)*(ir+nrootu*(j+i11s-1))      8d26s21
                        jad=itmpsv+nvirt(ksbv)*(ir+nrootu*j)            8d26s21
                        do jv=0,nvirt(ksbv)-1                           8d26s21
                         bc(jad+jv)=bc(iad+jv)                           1d28s21
                        end do                                            12d22s20
                       end do                                             12d22s20
                      end do                                              12d22s20
                      do k=0,nok4f(l)-1                                   12d22s20
                       kk=ibc(mdhvnotv(l)+k)                              12d22s20
                       do j=0,njhere-1                                    12d22s20
                        do ir=0,nrootu-1                                  12d22s20
                         iadv=joffdnon+nvv*(ir+nrootu*kk)                 12d22s20
                         iadd=itrans+nvirt(isbvu)*(j+njhere*k)            12d22s20
                         iadsv=itmpsv+nvirt(ksbv)*(ir+nrootu*j)            12d22s20
                         do iv2=0,nvirt(isbv2)-1                          12d22s20
                          do iv1=0,iv2-1                                  12d22s20
                           gd(iadv+iv1)=gd(iadv+iv1)                      12d22s20
     $                        +bc(iadd+iv1)*bc(iadsv+iv2)*tf2           12d29s20
                           gd(iadv+iv1)=gd(iadv+iv1)                     12d29s20
     $                        +bc(iadd+iv2)*bc(iadsv+iv1)               12d31s20
                          end do                                          12d22s20
                          iadv=iadv+iv2                                   12d22s20
                         end do                                           12d22s20
                        end do                                            12d22s20
                       end do                                             12d22s20
                      end do                                              12d22s20
                     else if(isbvu.eq.isbv1)then                          12d22s20
c     r v2 j, j f v1, v1v2 r f
c     f v1 r v2, f v1 j, j r v2
                      mrow=nvirt(isbv1)*nok4f(l)                          12d22s20
                      mcol=nvirt(isbv2)*nrootu                            12d22s20
                      itmpdv=ibcoff                                        12d22s20
                      itmpdg=itmpdv+mrow*mcol                             12d22s20
                      ibcoff=itmpdg+mrow*mcol                             12d22s20
                      call enough('hcdsbk4. 25',bc,ibc)
                      call dgemm('n','n',mrow,mcol,njhere,1d0,          3d2s21
     $                     bc(itrans),mrow,bc(itmpsv),njhere,0d0,       3d2s21
     $                     bc(itmpdg),mrow,                             3d2s21
     d' hcdsbk4. 10')
                      do iv2=0,nvirt(isbv2)-1                             12d22s20
                       do ir=0,nrootu-1                                   12d22s20
                        do iv1=0,nvirt(isbv1)-1                           12d22s20
                         do k=0,nok4f(l)-1                                12d22s20
                          kk=ibc(mdhvnotv(l)+k)                           12d22s20
                          iad1=itmpdg+k+nok4f(l)*(iv1+nvirt(isbv1)*(iv2  1d27s21
     $                        +nvirt(isbv2)*ir))                        1d27s21
                          iad2=joffdnon+iv1+nvirt(isbv1)*(iv2            12d22s20
     $                        +nvirt(isbv2)*(ir+nrootu*kk))             12d22s20
                          gd(iad2)=gd(iad2)+bc(iad1)                      12d22s20
                         end do                                           12d22s20
                        end do                                            12d22s20
                       end do                                             12d22s20
                      end do                                              12d22s20
                     else                                                 12d22s20
c     r v1 j, j f v2, v1v2 r f
c     f v1 r v2, f v2 j, j r v1
                      mrow=nrootu*nvirt(isbv1)                            12d22s20
                      mcol=nok4f(l)*nvirt(isbv2)                          12d22s20
                      itmpdv=ibcoff                                       12d22s20
                      itmpgv=itmpdv+mrow*mcol                             12d22s20
                      ibcoff=itmpgv+mrow*mcol                             12d22s20
                      call enough('hcdsbk4. 26',bc,ibc)
                      call dgemm('n','n',mcol,mrow,njhere,1d0,          3d2s21
     $                     bc(itrans),mcol,bc(itmpsv),njhere,0d0,       3d2s21
     $                     bc(itmpgv),mcol,                             3d2s21
     d' hcdsbk4. 11')
                      do ir=0,nrootm                                     1d28s21
                       do iv1=0,nvirt(isbv1)-1                             12d22s20
                        do iv2=0,nvirt(isbv2)-1                           12d22s20
                         do k=0,nok4f(l)-1                                12d22s20
                          kk=ibc(mdhvnotv(l)+k)                           12d22s20
                          iad1=itmpgv+k+nok4f(l)*(iv2+nvirt(isbv2)        12d22s20
     $                       *(iv1+nvirt(isbv1)*ir))                    1d28s21
                          iad2=joffdnon+iv1+nvirt(isbv1)*(iv2             12d22s20
     $                       +nvirt(isbv2)*(ir+nrootu*kk))              12d22s20
                          gd(iad2)=gd(iad2)+bc(iad1)                      12d22s20
                         end do                                           12d22s20
                        end do                                            12d22s20
                       end do                                             12d22s20
                      end do                                              12d22s20
                     end if                                               12d22s20
                    end if                                                12d22s20
                    ibcoff=intden                                         12d22s20
                   end if                                                 12d21s20
                  end if                                                 12d22s20
                  joffdnon=joffdnon+nvv*nfdat(2,l,isb)*nrootu            12d22s20
                 end if                                                  12d21s20
                end do                                                   12d21s20
               end if                                                     12d21s20
              end do                                                     12d21s20
             end if                                                     3d2s21
             jvs=jvs+ncsfk(jarg)*nvirt(ksbv)*nrootu                      8d13s21
            end do                                                      12d18s20
            ibcoff=idhvisv                                              12d24s20
           end if                                                       12d18s20
          end do                                                        12d18s20
          do isbv1=1,nsymb                                              12d21s20
           isbv2=multh(isbv1,isbv12)                                    12d21s20
           if(isbv1.le.isbv2)then                                       12d21s20
            if(isbv1.eq.isbv2)then                                      12d21s20
             nvvs=(nvirt(isbv1)*(nvirt(isbv1)+1))/2                     12d21s20
             nvvt=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                     12d21s20
            else                                                        12d21s20
             nvvs=nvirt(isbv1)*nvirt(isbv2)                             12d21s20
             nvvt=nvvs                                                  12d21s20
            end if                                                      12d21s20
            ioffdnon=ioffdnon+(nvvs*nfdat(2,1,isb)                      12d21s20
     $           +nvvt*(nfdat(2,2,isb)+nfdat(2,3,isb)+nfdat(2,4,isb)))  12d21s20
     $             *nrootu                                              12d21s20
           end if                                                       12d21s20
          end do                                                        12d21s20
         end do                                                         12d18s20
         itransgg=0                                                     1d30s21
         ibcoff=ibcgg                                                   1d30s21
        end if                                                          12d18s20
       end do                                                           12d18s20
       ioffdnon=1                                                       2d4s21
       if(l2e.eq.0)then                                                 2d17s22
        call dws_gsumf(bc(ibcvmat),nbmat)                                3d2s21
        call dws_gbor(i3xk,16)                                          12d21s22
        do i=1,ntype                                                    1d18s23
         if(i3xk(2,i).eq.0)i3xk(2,i)=1                                  1d26s23
         if(i3xk(2,i).lt.0.or.i3xk(2,i).gt.2)then                       1d18s23
          write(6,*)('bad i3xk!!! ')
          write(6,*)('for type '),i,(' of '),ntype,
     $         ('get '),i3xk(2,i)
          call dws_synca
          call dws_finalize
          stop 'hcdsbk4'
         end if                                                         1d18s23
        end do                                                          1d18s23
        do isb=1,nsymb                                                   2d4s21
         nfh=nfdat(2,1,isb)+nfdat(2,2,isb)+nfdat(2,3,isb)+nfdat(2,4,isb) 2d4s21
         isbv12=multh(isb,isymbra)                                       8d16s21
         isn=multh(isopt(1),multh(isbv12,ksbv))                          12d7s21
         if(min(irefo(isn),nvirt(ksbv)).gt.0)then                        2d25s21
          nn=nvirt(ksbv)*nrootu*nfh                                       2d4s21
c
c     form Gv'v"rk=[(vv"|nv')-(vv'|nv")]vvrkn
c
          do isbv1=1,nsymb                                               2d4s21
           isbv2=multh(isbv1,isbv12)                                     2d4s21
           if(isbv2.ge.isbv1)then                                        2d4s21
            call ilimts(irefo(isn),nvirt(isbv1),mynprocg,mynowprog,il,
     $          ih,i1s,i1e,i2s,i2e)                                         2d4s21
            nhere=ih+1-il                                                12d7s21
            ilm=il-1                                                     12d7s21
            if(isbv1.eq.isbv2)then                                       2d4s21
             i10=i1s                                                     2d4s21
             i1n=irefo(isn)                                              2d4s21
             do i2=i2s,i2e                                               12d7s21
              iv1=i2-1                                                   12d7s21
              if(i2.eq.i2e)i1n=i1e                                       2d4s21
              do i1=i10,i1n                                              2d4s21
               i1m=i1-1
               do k=0,nfdat(2,1,isb)-1                                   2d4s21
                do ir=0,nrootm                                           2d4s21
                 sum=0d0                                                 2d4s21
                 iadv=ioffdnon+iv1+nvirt(isbv1)*(ir+nrootu*k)            2d4s21
                 do it=0,ntype-1                                         12d7s21
                  jvmat=ivmat(isb)+nvirt(ksbv)*(ir+nrootu*(k+nfh*        12d7s21
     $                (i1m+irefo(isn)*it)))                             12d7s21
                  ixv=i3x(ksbv,isbv1,isn)+nvirt(ksbv)*(iv1               12d8s21
     $                +nvirt(isbv1)*(i1m+irefo(isn)*iv1-ilm             12d8s21
     $                 +nhere*it))                                      12d8s21
                  do kv=0,nvirt(ksbv)-1                                    2d4s21
                   sum=sum+bc(jvmat+kv)*bc(ixv+kv)                       12d8s21
                  end do                                                 2d4s21
                 end do                                                  12d7s21
                 sum=sum*sr2                                             12d9s21
                 gd(iadv)=gd(iadv)+sum                                   2d4s21
                end do                                                   2d4s21
               end do                                                    2d4s21
              end do                                                     2d4s21
              i10=1                                                      2d4s21
             end do                                                      2d4s21
             ioffdnon=ioffdnon+nvirt(isbv1)*nrootu*nfdat(2,1,isb)        2d4s21
             nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       2d4s21
             isw=0                                                       2d4s21
            else                                                         2d4s21
             nvv=nvirt(isbv1)*nvirt(isbv2)                               2d4s21
             isw=1                                                       2d4s21
            end if                                                       2d4s21
            i10=i1s                                                      2d4s21
            i1n=irefo(isn)                                               2d4s21
            do i2=i2s,i2e                                                12d7s21
             iv1=i2-1                                                    12d7s21
             ibots=i2                                                    2d4s21
             ibotn=0                                                     2d4s21
             ibot=ibots+isw*(ibotn-ibots)                                2d4s21
             if(i2.eq.i2e)i1n=i1e                                        2d4s21
             do i1=i10,i1n                                               2d4s21
              i1m=i1-1                                                   8d16s21
              do kr=0,nfh*nrootu-1                                       2d4s21
               do iv2=ibot,nvirt(isbv2)-1                                2d4s21
                sum=0d0                                                  2d4s21
                do it=0,ntype-1                                          12d7s21
                 jvmat=ivmat(isb)+nvirt(ksbv)*kr+nn*(i1m+irefo(isn)*it)  12d7s21
                 ixv=i3x(ksbv,isbv1,isn)+nvirt(ksbv)*(iv2                12d8s21
     $               +nvirt(isbv2)*(i1m+irefo(isn)*iv1-ilm+nhere*it))   12d8s21
                 do kv=0,nvirt(ksbv)-1                                    2d4s21
                  sum=sum+bc(jvmat+kv)*bc(ixv+kv)                        12d8s21
                 end do                                                  2d4s21
                end do                                                   12d7s21
                itri=((iv2*(iv2-1))/2)+iv1                               2d4s21
                irec=iv1+nvirt(isbv1)*iv2                                2d4s21
                iadv=ioffdnon+itri+isw*(irec-itri)+nvv*kr                2d4s21
                gd(iadv)=gd(iadv)+sum                                    12d7s21
               end do                                                    2d4s21
              end do                                                     2d4s21
             end do                                                      2d4s21
             i10=1                                                       2d4s21
            end do                                                       2d4s21
            call ilimts(irefo(isn),nvirt(isbv2),mynprocg,mynowprog,il,   8d16s21
     $           ih,i1s,i1e,i2s,i2e)                                    8d16s21
            nhere=ih+1-il                                                12d7s21
            ilm=il-1                                                     12d7s21
            i10=i1s                                                      2d4s21
            i1n=irefo(isn)                                               2d4s21
            do i2=i2s,i2e                                                2d4s21
             iv2=i2-1                                                    2d4s21
             itops=iv2-1                                                 2d4s21
             itopn=nvirt(isbv1)-1                                        2d4s21
             itop=itops+isw*(itopn-itops)                                2d4s21
             if(i2.eq.i2e)i1n=i1e                                        2d4s21
             do i1=i10,i1n                                               2d4s21
              i1m=i1-1                                                   2d4s21
              do kr=0,nfh*nrootu-1                                       2d4s21
               do iv1=0,itop                                             2d4s21
                sum=0d0                                                  2d4s21
                do it=0,ntype-1                                          12d7s21
                 itu=i3xk(1,it+1)                                        12d9s21
                 pp=phss(i3xk(2,it+1))                                   12d9s21
                 jvmat=ivmat(isb)+nvirt(ksbv)*kr+nn*(i1m+irefo(isn)*it)  12d7s21
                 ixv=i3x(ksbv,isbv2,isn)+nvirt(ksbv)*(iv1                12d8s21
     $                +nvirt(isbv1)*(i1m+irefo(isn)*iv2-ilm+nhere*itu))  12d9s21
                 do kv=0,nvirt(ksbv)-1                                    2d4s21
                  sum=sum+bc(jvmat+kv)*bc(ixv+kv)*pp                     12d9s21
                 end do                                                  2d4s21
                end do                                                   12d7s21
                itri=((iv2*(iv2-1))/2)+iv1                               2d4s21
                irec=iv1+nvirt(isbv1)*iv2                                2d4s21
                iadv=ioffdnon+itri+isw*(irec-itri)+nvv*kr                2d4s21
                gd(iadv)=gd(iadv)+sum                                    12d9s21
               end do                                                    2d4s21
              end do                                                     2d4s21
             end do                                                      2d4s21
             i10=1                                                       2d4s21
            end do                                                       2d4s21
            ioffdnon=ioffdnon+nvv*nrootu*nfh                             2d4s21
           end if                                                        12d7s21
          end do
         else                                                            12d9s21
c     still need to update ioffdnon!
          do isbv1=1,nsymb                                               12d9s21
           isbv2=multh(isbv1,isbv12)                                     12d9s21
           if(isbv2.ge.isbv1)then                                        12d9s21
            if(isbv1.eq.isbv2)then                                       12d9s21
             ioffdnon=ioffdnon+nvirt(isbv1)*nrootu*nfdat(2,1,isb)        12d9s21
             nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       12d9s21
            else                                                         12d9s21
             nvv=nvirt(isbv1)*nvirt(isbv2)                               12d9s21
            end if                                                       12d9s21
            ioffdnon=ioffdnon+nvv*nrootu*nfh                             12d9s21
           end if                                                        12d9s21
          end do                                                         12d9s21
         end if                                                          12d7s21
        end do                                                           12d7s21
       end if                                                           2d17s22
       ibcoff=ibcvmat                                                   12s15s21
      end do                                                            12d18s20
      ibcoff=ircv                                                       1d30s21
      return
      end
