c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcdsbk(ihsdiagk,nff1,iff1,nff22,nfdat,gd,              8d26s21
     $     nsymb,mdon,mdoop,nec,multh,isymbra,isymket,nvirt,ncsf,ncsf2, 8d16s21
     $     irel,ism,idoubo,irefo,nbasdws,ixw1,ixw2,norb,nrootu,ixmtf,   8d16s21
     $     phase1,lxmt,isymop,n2e,ixmt,i2eop,phase2,ldebug,             8d16s21
     $     maxbx,sr2,iff22,ff22,bc,ibc)                                 11d10s22
      implicit real*8 (a-h,o-z)                                         12d18s20
c
c     to do ...
c     consolidate densities, and perhaps njhere in density calculation
c     did I swap iarg and jarg in 4th gandc4?
c
      integer*8 ihsdiagk(mdoop,nsymb,2),i18,i28,                        8d26s21
     $   i38,i48,i1c,i1o,j2c,j2o,itestc,itesto,ipack8,last8(2),iff22(*),12d21s22
     $     gandcc,gandco,gandcb                                         12d21s22
      integer*1 nab1(2),nab2(2)                                         12d18s20
      external second                                                   2d18s21
      logical l3x,lprt,lnew,ldebug,lchoice                              3d17s21
      equivalence (ipack8,ipack4)                                       12d18s20
      dimension nff1(mdoop,nsymb,2),iff1(*),nff22(mdoop,2,nsymb),       8d13s21
     $     nfdat(5,4,*),gd(*),multh(8,8),nvirt(*),ncsf(*),ncsf2(4,*),   12d18s20
     $     irel(*),ism(*),irefo(*),itest(32,3),ioxx(2),                 12d21s22
     $     nab4(2,3),ipack4(2),nl(4),npre(4),mpre(4),                   8d13s21
     $     id1visv(8,8),nd1visv(8,8),nokdc(8,8,4),nok3v(4),             12d22s20
     $     idhvnotv(4),ndhvnotv(4),id1vnotv(4,8,8),nd1vnotv(4,8,8),     12d21s20
     $     id3vnotv3(4,2),nd3vnotv3(4,2),loff(4),idkeep(2),ndkeep(2),   1d4s21
     $     mdkeep(2),keep(2),idhvnotvf(4),ibmat(8),ivmat(8),            2d4s21
     $     mdhvnotv(4),md3vnotv3(4,2),nok4f(4),nok4(4),nok3vf(4),       12d22s20
     $     ixmtf(*),nbasdws(*),nok33f(4,2),nok33(4,2),ixmt(8,*),
     $     isymop(*),i2eop(2,3),idoubo(*),ff22(*)                       8d16s21
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data loopx/11810000/
      data icall/0/
      save icall
      icall=icall+1
      loop=0
      nrootm=nrootu-1                                                   1d4s21
      idoit=0                                                           3d1s21
      mdoo=mdoop-1                                                      8d13s21
      last8(1)=-1                                                       2d8s21
      ircv=ibcoff                                                       1d29s21
      iacc=ircv+mynprocg                                                1d30s21
      ivs=iacc+mynprocg                                                 1d30s21
      igg=ivs+maxbx                                                     1d30s21
      ibcoff=igg+maxbx                                                  1d30s21
      nacc=0                                                            1d30s21
      call enough('hcdsbk.  1',bc,ibc)
      itransgg=0                                                        1d30s21
      norbx=norb+1                                                      12d18s20
      norbxx=norbx+1                                                    12d18s20
      norbxxx=norbxx+1                                                  12d18s20
      do jsb=1,nsymb                                                    12d18s20
       ksbv=multh(jsb,isymket)                                          8d13s21
c
c     let us form bvrkn=[(vv"|nv')+p(k)(vv'|nv")]Vv'v"rk,               2d3s21
c     v ne v' and v"                                                    2d3s21
c     and space for vvrkn=Vvrj*Djkn                                     2d4s21
c     gandc4(Vsv,Vdv'v"), idx=1 (iv'|vv"), idx=2 (iv"|vv')
c
       ioffvd=1                                                         2d3s21
       ibcbmat=ibcoff                                                   2d3s21
       if(n2e.gt.0)then                                                 8d13s21
        ibcvmat=ibcoff                                                   3d2s21
        nbmat=0                                                         8d27s21
        do isb=1,nsymb                                                   3d2s21
         isbv12=multh(isb,isymbra)                                       8d13s21
         isn=multh(isbv12,ksbv)                                         8d26s21
         ivmat(isb)=ibcoff                                               2d3s21
         if(min(irefo(isn),nvirt(ksbv)).gt.0)then                       8d26s21
          nftrip=nfdat(2,2,isb)+nfdat(2,3,isb)+nfdat(2,4,isb)             2d3s21
          nfh=nfdat(2,1,isb)+nftrip                                       2d3s21
          nn=nvirt(ksbv)*nrootu*nfh                                     8d26s21
          nbmat=nbmat+nn*irefo(isn)                                     8d27s21
          ibcoff=ivmat(isb)+nn*irefo(isn)                                3d2s21
         end if                                                          3d2s21
        end do                                                           3d2s21
        call enough('hcdsbk.  2',bc,ibc)
        do i=ibcbmat,ibcoff-1                                            3d2s21
         bc(i)=0d0                                                       3d2s21
        end do                                                           3d2s21
       end if                                                           8d13s21
       do nclo1p=mdon+1,mdoop                                           8d13s21
        if(min(nff1(nclo1p,jsb,1),nvirt(ksbv)).gt.0)then                8d26s21
         nclo1=nclo1p-1                                                  12d12s20
         jarg=nclo1p-mdon                                                12d12s20
         nopen1=nec-2*nclo1                                              12d12s20
         nggk=nvirt(ksbv)*nrootu                                        8d13s21
         nggg=nff1(nclo1p,jsb,1)*ncsf(jarg)                             12d29s20
         call ilimts(1,nff1(nclo1p,jsb,1),mynprocg,mynowprog,           3d2s21
     $        i1l,i1h,i11s,i11e,i12s,i12e)                              12d19s20
         i1l=1+ncsf(jarg)*(i1l-1)                                       3d2s21
         i1h=ncsf(jarg)*i1h                                             3d2s21
         ibcgg=ibcoff                                                   1d30s21
         i18=1
         i28=nggk                                                       8d13s21
         i38=1                                                          12d18s20
         i48=nggg                                                       12d29s20
         call ddi_iget(bc,ibc,ihsdiagk(nclo1p,jsb,1),i18,i28,i38,i48,   11d15s22
     $        bc(ivs),ibc(ircv),nrcv)                                   11d15s22
         itransvs=0                                                     1d29s21
         ioffdnon=1                                                     12d21s20
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
          do nclo2p=max(mdon+1,nclo1p-2),min(mdoop,nclo1p+3)            8d13s21
           if(nff22(nclo2p,1,isb).gt.0)then                             12d18s20
            nclo2=nclo2p-1                                              12d18s20
            iarg=nclo2p-mdon                                            12d18s20
            iargp=iarg+1                                                12d18s20
            nopen2=nec-2*nclo2p                                         12d18s20
            nopen2p=nopen2+2                                            12d18s20
            idhvisv=ibcoff                                              12d18s20
            lgoal=multh(jsb,isb)                                        12d22s20
            lgoal3=multh(ksbv,isbv12)                                   12d22s20
            do isc=1,nsymb                                              12d18s20
             do l=1,4                                                   12d18s20
              nnl=ncsf(jarg)*nfdat(2,l,isb)*irefo(isc)                  12d19s20
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
              nnn=nn*irefo(jscdv)                                       12d18s20
              do l=1,4                                                  12d19s20
               id1vnotv(l,isd,isc)=ibcoff                               12d19s20
               nd1vnotv(l,isd,isc)=id1vnotv(l,isd,isc)+nnn*ncsf(jarg)   12d19s20
     $              *nfdat(2,l,isb)                                     12d19s20
               ibcoff=nd1vnotv(l,isd,isc)+nnn                           12d21s20
              end do                                                    12d19s20
             end do                                                     12d18s20
            end do                                                      12d18s20
            do ipass=1,1                                                3d1s21
             idkeep(ipass)=ibcoff                                       1d4s21
             ndkeep(ipass)=idkeep(ipass)+irefo(lgoal3)*nnt*ncsf(jarg)    8d13s21
             mdkeep(ipass)=ndkeep(ipass)+irefo(lgoal3)*nnt               8d13s21
             ibcoff=mdkeep(ipass)+irefo(lgoal3)*nnt                      8d13s21
             call enough('hcdsbk.  3',bc,ibc)
            end do                                                      1d4s21
            do l=1,4
             nnl=ncsf(jarg)*nfdat(2,l,isb)*irefo(lgoal3)                 8d13s21
             id3vnotv3(l,1)=ibcoff                                      12d22s20
             nd3vnotv3(l,1)=id3vnotv3(l,1)+nnl                          12d22s20
             id3vnotv3(l,2)=nd3vnotv3(l,1)+irefo(lgoal3)                 8d13s21
             nd3vnotv3(l,2)=id3vnotv3(l,2)+nnl                          12d22s20
             md3vnotv3(l,1)=nd3vnotv3(l,2)+irefo(lgoal3)                 8d13s21
             md3vnotv3(l,2)=md3vnotv3(l,1)+nfdat(2,l,isb)               12d22s20
             ibcoff=md3vnotv3(l,2)+nfdat(2,l,isb)                       12d22s20
            end do                                                      12d22s20
            ibcb4=ibcoff-1                                              12d18s20
            if1o=nff1(nclo1p,jsb,2)                                      12d14s20
            jvs=ivs                                                     12d18s20
            do if1=1,nff1(nclo1p,jsb,1)                                  12d14s20
             ist=1+ncsf(jarg)*(if1-1)                                   1d5s21
             ien=ncsf(jarg)*if1                                         1d5s21
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
               if(nl(l).gt.ncsf2(l,iarg))lchoice=.true.                 1d13s23
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
                   if((btest(i1o,i).and..not.btest(j2c,i)).or.
     $                (btest(i1c,i).and.btest(j2o,i)))then                           10d14s22
                    nab4(1,1)=i
                   else                                                           10d14s22
                    nab4(2,1)=i
                   end if
                  end if                                                          10d14s22
                 end do                                                           10d14s22
c     i1o is rest,a,x
c     j2o is rest,x,xx
c     so nab1(2) is xx,nab1(1) is a.
c     a*xx should have sym isymbra*isymket
c     rest*a*ksbv=isymket, rest*a*jsb*isymket=isymket
c     rest*isbv12=isymbra, rest*isb*isymbra=isymbra,
c     or rest*a*jsb=1=rest*isb or a*jsb=isb or a=jsb*isb
c     a*xx?isymbra*isymket
c     xx=isbv1*isbv12=ksbv*isb*isymbra=jsbv*isymket*isb*isymbra
c     thus a*xx=isymbra*isymket as required.
c     if j is ket, then x is jsbv and xx is jsbv*ksbv12=lgoal
c     if j is bra, then x is ksbv and xx is ksbv*isbv12=lgoal
                call gandc(i1c,i1o,j2c,j2o,nopen1,nopen2p,jarg,iarg,    12d18s20
     $               ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,       12d18s20
     $               ncsfmid1,bc,ibc)                                   11d14s22
                lsa=ism(nab1(1))
                lga=irel(nab1(1))-1                                     12d18s20
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
                iargo=1                                                 2d25s21
                iprod=ibcoff                                            3d19s21
                if(lchoice)then                                         3d19s21
                 ibcoff=iprod+ncsf(jarg)*ncsf(iarg)                     3d19s21
                 call enough('hcdsbk.  4',bc,ibc)
                 call prodn(iwpb1,iwpk1,ncsf(jarg),ncsf(iarg),ncsfmid1, 3d19s21
     $                bc(iprod),bc,ibc,1d0,0d0)                         2d13s23
                end if                                                  3d19s21
                jprod=iprod                                             3d19s21
                do l=1,4                                                12d18s20
                 if(nl(l).gt.0)then                                     12d18s20
                  nnl=ncsf(jarg)*nfdat(2,l,isb)                         12d19s20
                  iad1=jvcv+iff22(jvcv+5+l)                             8d16s21
                  iad2=iad1+nl(l)                                       3d19s21
                  itmp=ibcoff                                           12d18s20
                  ibcoff=itmp+ncsf(jarg)*nl(l)                          12d18s20
                  call enough('hcdsbk.  5',bc,ibc)
                  if(lchoice)then                                       3d19s21
                   call dgemm('n','n',ncsf(jarg),nl(l),ncsf2(l,iarg),   3d19s21
     $                1d0,bc(jprod),ncsf(jarg),ff22(iad2),ncsf2(l,iarg),8d16s21
     $                 0d0,bc(itmp),ncsf(jarg),                         3d19s21
     d' hcdsbk.  1')
                  else                                                  3d19s21
                   call xtimesn2(ncsf(jarg),ncsf(iarg),ncsfmid1,nl(l),   2d25s21
     $                 iwpb1,iwpk1,ff22(iad2),ncsf2(l,iarg),bc(itmp),   8d16s21
     $                 ncsf(jarg),1d0,0d0,iargo,ncsf2(l,iarg),bc,ibc)   11d10s22
                  end if                                                3d19s21
c     jsbv=jsb*isymbra, ksbv=jsb*isymket                                8d13s21
c     isbv12=isb*isymbra, ksbv12=jsb*isymket                            8d13s21
c     isbv(1or2)=ksbv12*jsbv=jsb*isymbra*isb*isymket
c     isbv(1or2)=isbv12*ksbv=jsb*isymket*isb*isymbra                    8d13s21
                  jden=idhvnotv(l)+nnl*lga                              12d21s20
                  ibc(ndhvnotv(l)+lga)=1                                12d21s20
                  jtmp=itmp                                             12d18s20
                  do iii=0,nl(l)-1                                      12d18s20
                   ip=iff22(iad1+iii)-1                                 8d16s21
                   jjden=jden+ncsf(jarg)*ip                             12d18s20
                   ibc(mdhvnotv(l)+ip)=1                                12d21s20
                   do j=0,ncsf(jarg)-1                                  12d18s20
                    bc(jjden+j)=bc(jjden+j)+bc(jtmp+j)                  12d18s20
                   end do                                               12d18s20
                   jtmp=jtmp+ncsf(jarg)                                 12d18s20
                  end do                                                12d18s20
                  if(n2e.gt.0)then                                      8d15s21
                   do i=1,nok-1                                           12d19s20
                    js=ism(itest(i,2))                                           11d13s20
                    jg=irel(itest(i,2))-1                                        11d13s20
                    icolj=jg*(irefo(js)+1)                               8d13s21
                    nn=irefo(js)*irefo(js)                               8d13s21
                    icolj=icolj+nn*lga                                    12d18s20
                    jdenj=id1vnotv(l,js,js)+nnl*icolj                    12d19s20
                    ibc(nd1vnotv(l,js,js)+icolj)=1                       12d19s20
                    if(itest(i,3).eq.2)then                                  12d14s20
                     nn=irefo(lsa)*irefo(js)                             12d18s20
                     icolk=lga+irefo(lsa)*(jg+irefo(js)*jg)             8d16s21
                     jdenk=id1vnotv(l,lsa,js)+nnl*icolk                 8d16s21
                     ibc(nd1vnotv(l,lsa,js)+icolk)=1                    8d16s21
                     jtmp=itmp                                              12d18s20
                     do iii=0,nl(l)-1                                         12d18s20
                      ii=iff22(iad1+iii)-1                              8d16s21
                      jdj=jdenj+ncsf(jarg)*ii                             12d18s20
                      jdk=jdenk+ncsf(jarg)*ii                             12d18s20
                      do j=0,ncsf(jarg)-1                                   12d18s20
                       bc(jdj+j)=bc(jdj+j)+2d0*bc(jtmp+j)                 12d18s20
                       bc(jdk+j)=bc(jdk+j)-bc(jtmp+j)                     12d18s20
                      end do                                                12d18s20
                      jtmp=jtmp+ncsf(jarg)                                  12d18s20
                     end do                                                 12d18s20
                    else                                                     12d14s20
                     jtmp=itmp                                              12d18s20
                     do iii=0,nl(l)-1                                         12d18s20
                      ii=iff22(iad1+iii)-1                              8d16s21
                      jdj=jdenj+ncsf(jarg)*ii                             12d18s20
                      do j=0,ncsf(jarg)-1                                   12d18s20
                       bc(jdj+j)=bc(jdj+j)+bc(jtmp+j)                     12d18s20
                      end do                                                12d18s20
                      jtmp=jtmp+ncsf(jarg)                                  12d18s20
                     end do                                                 12d18s20
                    end if                                                   12d14s20
                   end do                                                    12d14s20
                  end if                                                8d15s21
                  ibcoff=itmp                                           12d18s20
                 end if                                                 12d18s20
                 jprod=jprod+ncsf2(l,iarg)*ncsf(jarg)                   3d19s21
                 iargo=iargo+ncsf2(l,iarg)                              2d25s21
                end do                                                  12d18s20
                ibcoff=iprod                                            3d19s21
                if(n2e.gt.0)then                                        8d15s21
                 do i=1,nok                                              12d19s20
                  if(itest(i,3).eq.1)then                                  12d14s20
                   itestc=i1c                                              12d14s20
                   itesto=i1o                                              12d14s20
                   nopenk=nopen1                                                11d13s20
c
c     anihilate common
c
                   if(btest(itestc,itest(i,2)))then                             11d13s20
                    itestc=ibclr(itestc,itest(i,2))                             11d13s20
                    itesto=ibset(itesto,itest(i,2))                             11d13s20
                    karg=jarg-1                                                11d13s20
                    nopenk=nopenk+1                                             11d13s20
                   else                                                         11d13s20
                    itesto=ibclr(itesto,itest(i,2))                             11d13s20
                    karg=jarg                                                  11d13s20
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
                    call gandc(i1c,i1o,itestc,itesto,nopen1,nopenk,        12d14s20
     $                jarg,karg,ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1, 12d28s20
     $                iwpk1,ncsfmid1,bc,ibc)                            11d14s22
                    call gandc(itestc,itesto,j2c,j2o,nopenk,nopen2p,        12d14s20
     $                karg,iarg,ncsf,norbxx,ixw1,ixw2,nnot2,nab2,iwpb2, 12d28s20
     $                iwpk2,ncsfmid2,bc,ibc)                            11d14s22
c     bra(i) is rest acx, sym(rest)*sym(a)*sym(c)*jsbv=isymbra
c                         sym(rest)*sym(a)*sym(c)*jsb=1
c     test is rest ax,xx,
c     ket(j) is rest c x,xx, sym(rest)*sym(c)*ksbv12=isymket
c                            sym(rest)*sym(c)*isb=1
c                            so sym(a)*jsb*isb=1
c
c     nab1(1)=a,nab1(2)=xx    (a|xx)*(a|c)
c            lsa              lsa    lsc lsb
c     nab2(1)=a,nab2(2)=c
c            lsc       lsb
c
                    if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                     if(max(nab2(1),nab2(2)).le.norb)then                2d26s21
                      lsa=ism(nab2(1))                                      12d14s20
                      lga=irel(nab2(1))-1                                   12d14s20
                      lsb=ism(nab2(2))                                      12d14s20
                      lgb=irel(nab2(2))-1                                   12d14s20
                      lsc=ism(nab1(1))                                    12d18s20
                      lgc=irel(nab1(1))-1                                 12d18s20
                      nn=irefo(lsa)*irefo(lsb)                           12d18s20
                      icol=lga+irefo(lsa)*(lgb+irefo(lsb)*lgc)           12d18s20
                      iprod=ibcoff                                       3d19s21
                      if(lchoice)then                                    3d19s21
                       call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),     3d19s21
     $                  ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod,11d10s22
     $                     bc,ibc)                                      11d10s22
                      end if                                             3d19s21
                      iargo=1                                            2d25s21
                      jprod=iprod                                        3d19s21
                      do l=1,4
                       if(nl(l).gt.0)then                                     12d18s20
                        nnl=ncsf(jarg)*nfdat(2,l,isb)                         12d19s20
                        iad1=jvcv+iff22(jvcv+5+l)                       8d16s21
                        iad2=iad1+nl(l)                                  3d19s21
                        itmp=ibcoff                                           12d18s20
                        ibcoff=itmp+ncsf(jarg)*nl(l)                          12d18s20
                        call enough('hcdsbk.  6',bc,ibc)
                        if(lchoice)then                                  3d19s21
                         call dgemm('n','n',ncsf(jarg),nl(l),            3d19s21
     $                      ncsf2(l,iarg),1d0,bc(jprod),ncsf(jarg),     3d19s21
     $                      ff22(iad2),ncsf2(l,iarg),0d0,bc(itmp),      8d16s21
     $                      ncsf(jarg),                                 3d19s21
     d' hcdsbk.  2')
                        else                                             3d19s21
                         call genmatn2(ncsf(jarg),ncsf(karg),ncsf(iarg),  2d25s21
     $                  ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,      2d25s21
     $                  ff22(iad2),ncsf2(l,iarg),nl(l),bc(itmp),0d0,    8d16s21
     $                  iargo,ncsf2(l,iarg),1d0,bc,ibc)                 11d10s22
                        end if                                           3d19s21
                        jden=id1vnotv(l,lsa,lsb)+nnl*icol                     12d18s20
                        ibc(nd1vnotv(l,lsa,lsb)+icol)=1                       12d18s20
                        jtmp=itmp                                              12d18s20
                        mden=mdhvnotv(l)                                 12d21s20
                        do iii=0,nl(l)-1                                         12d18s20
                         ii=iff22(iad1+iii)-1                           8d16s21
                         jdh=jden+ncsf(jarg)*ii                             12d18s20
                         ibc(mden+ii)=1                                  12d19s20
                         do j=0,ncsf(jarg)-1                                   12d18s20
                          bc(jdh+j)=bc(jdh+j)+bc(jtmp+j)                       12d18s20
                         end do                                                12d18s20
                         jtmp=jtmp+ncsf(jarg)                                  12d18s20
                        end do                                                 12d18s20
                        ibcoff=itmp                                      12d19s20
                       end if                                              12d19s20
                       iargo=iargo+ncsf2(l,iarg)                         2d25s21
                       jprod=jprod+ncsf(jarg)*ncsf2(l,iarg)              3d19s21
                      end do                                             12d19s20
                      ibcoff=iprod                                       3d19s21
                     end if                                                  12d14s20
                    end if                                                   12d14s20
                   end if                                                12d19s20
                  end if                                                 12d19s20
                 end do                                                  12d18s20
                end if                                                  8d15s21
                else
                 nnot=0
                 if(ndifs.eq.4.and.ndifb.eq.4)then                                           10d14s2
                  nnot=4                                                          10d14s22
                  ioxx(1)=1                                                     10d17s22
                  ioxx(2)=1                                                     10d17s22
                  do i=1,norbxx                                                     10d17s22
                   if(btest(gandcb,i))then                                         10d14s22
                    if((btest(j2c,i).and.btest(i1o,i)).or.                         10d17s22
     $                 (btest(j2o,i).and..not.btest(i1c,i)))then                   10d14s22
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
     $                 ((btest(i1c,i).and..not.btest(j2o,i)).or.                 10d14s22
     $                 (btest(j2c,i).and..not.btest(i1o,i))))then                     10d14s22
                     if(btest(j2c,i))iswap=1                                        10d17s22
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
                   if(btest(i1c,nab4(1,2)).and.                         10d26s22
     $                  .not.btest(i1c,nab4(1,1)))nbt=1                 10d26s22
                  else                                                             10d17s22
                   nbt=0                                                           10d17s22
                   if(btest(j2c,nab4(2,2)).and.                         10d26s22
     $                  .not.btest(j2c,nab4(2,1)))nbt=1                 10d26s22
                  end if                                                           10d17s22
                  if(nbt.ne.0)then                                                 10d17s22
                   nab4(1,1)=nab4(1,2)                                             10d17s22
                   nab4(2,1)=nab4(2,2)                                             10d17s22
                  end if                                                           10d17s22
                 else if(ndifs.eq.0.and.ndifd.eq.2)then                            10d14s22
                  nnot=3
                  do i=1,norbxx
                   if(btest(gandcb,i))then                                         10d14s22
                    if(btest(i1c,i))then
                     nab4(1,1)=i
                     nab4(1,2)=i
                    else                                                           10d14s22
                     nab4(2,1)=i
                     nab4(2,2)=i
                    end if
                   end if                                                          10d14s22
                  end do
                 end if                                                            10d14s22
                 ipssx=0                                                 10d21s22
                 if(nnot.eq.3)then                                          12d8s20
                  ipssx=1                                                   12d8s20
                 else if(nnot.eq.4)then                                  10d21s22
                  ipssx=3                                                12d18s20
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
                 itestc=i1c                                              12d8s20
                 itesto=i1o                                              12d8s20
                 if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                  itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                  itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                  nopenk=nopen1+1                                              11d13s20
                  karg=jarg-1                                             12d8s20
                 else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                  itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                  nopenk=nopen1-1                                              11d13s20
                  karg=jarg
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
                 nqq=karg+mdon-1                                        4d20s21
                 if(nqq.ge.mdon.and.nqq.le.mdoo)then                    4d20s21
                  call gandc(i1c,i1o,itestc,itesto,nopen1,nopenk,        12d18s20
     $         jarg,karg,ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,   12d18s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
                  call gandc(itestc,itesto,j2c,j2o,nopenk,nopen2p,         12d8s20
     $         karg,iarg,ncsf,norbxx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,    11d13s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
c     bra(i) is rest abx, sym(rest)*sym(a)*sym(b)*jsbv=isymbra
c                         sym(rest)*sym(a)*sym(b)*jsb=1
c     test is rest bx,xx,    or ax,xx,    or bx,c      or ax,c
c              (a|xx)(b|c)   (b|xx)(a|c)  (a|c)(b|xx)  (b|c)(a|xx)
c               l     l l     l     l l    l l  l       l l  l
c               s     s s     s     s s    s s  s       s s  s
c               c     a b     c     a b    a b  c       a b  c
c     ket(j) is rest c x,xx, sym(rest)*sym(c)*ksbv12=isymket
c                            sym(rest)*sym(c)*isb=1
c                            so sym(a)*sym(b)*sym(c)*jsb*isb=1
c                            or sym(c)=sym(a)*sym(b)*jsb*isb
c     xx*jsbv=ksbv12 or xx*jsb*isymbra=isb*isymket or
c     xx=isb*jsb*isymbra*isymket.
c
                  if(nnot1.eq.2.and.nnot2.eq.2)then                      12d19s20
                   if(max(nab1(1),nab1(2)).gt.norb)then                  12d19s20
                    lsa=ism(nab2(1))                                     12d19s20
                    lga=irel(nab2(1))-1                                  12d19s20
                    lsb=ism(nab2(2))                                     12d19s20
                    lgb=irel(nab2(2))-1                                  12d19s20
                    lsc=ism(nab1(1))                                     12d19s20
                    lgc=irel(nab1(1))-1                                  12d19s20
                   else                                                  12d19s20
                    lsa=ism(nab1(1))                                     12d19s20
                    lga=irel(nab1(1))-1                                  12d19s20
                    lsb=ism(nab1(2))                                     12d19s20
                    lgb=irel(nab1(2))-1                                  12d19s20
                    lsc=ism(nab2(1))                                     12d19s20
                    lgc=irel(nab2(1))-1                                  12d19s20
                   end if                                                12d19s20
                   icol=lga+irefo(lsa)*(lgb+irefo(lsb)*lgc)             12d19s20
                   iprod=ibcoff                                          3d19s21
                   if(lchoice)then                                       3d19s21
                    call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),        3d19s21
     $                 ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod, 11d10s22
     $                  bc,ibc)                                         11d10s22
                   end if                                                3d19s21
                   iargo=1                                               2d25s21
                   jprod=iprod                                           3d19s21
                   do l=1,4
                    if(nl(l).gt.0)then                                     12d18s20
                     nnl=ncsf(jarg)*nfdat(2,l,isb)                         12d19s20
                     iad1=jvcv+iff22(jvcv+l+5)                          8d16s21
                     iad2=iad1+nl(l)                                     3d19s21
                     itmp=ibcoff                                           12d18s20
                     ibcoff=itmp+ncsf(jarg)*nl(l)                          12d18s20
                     call enough('hcdsbk.  7',bc,ibc)
                     if(lchoice)then                                     3d19s21
                      call dgemm('n','n',ncsf(jarg),nl(l),ncsf2(l,iarg), 3d19s21
     $                1d0,bc(jprod),ncsf(jarg),ff22(iad2),ncsf2(l,iarg),8d16s21
     $                  0d0,bc(itmp),ncsf(jarg),                        3d19s21
     d' hcdsbk.  3')
                     else                                                3d19s21
                      call genmatn2(ncsf(jarg),ncsf(karg),ncsf(iarg),     2d25s21
     $                 ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,       2d25s21
     $                 ff22(iad2),ncsf2(l,iarg),nl(l),bc(itmp),0d0,     8d16s21
     $                 iargo,ncsf2(l,iarg),1d0,bc,ibc)                  11d10s22
                     end if                                              3d19s21
                     jden=id1vnotv(l,lsa,lsb)+nnl*icol                     12d18s20
                     ibc(nd1vnotv(l,lsa,lsb)+icol)=1                       12d18s20
                     jtmp=itmp                                              12d18s20
                     mden=mdhvnotv(l)                                    12d21s20
                     do iii=0,nl(l)-1                                         12d18s20
                      ii=iff22(iad1+iii)-1                              8d16s21
                      jdh=jden+ncsf(jarg)*ii                             12d18s20
                      ibc(mden+ii)=1                                     12d19s20
                      do j=0,ncsf(jarg)-1                                   12d18s20
                       bc(jdh+j)=bc(jdh+j)+bc(jtmp+j)                       12d18s20
                      end do                                                12d18s20
                      jtmp=jtmp+ncsf(jarg)                                  12d18s20
                     end do                                                 12d18s20
                     ibcoff=itmp                                         12d19s20
                    end if                                               12d19s20
                    iargo=iargo+ncsf2(l,iarg)                            2d25s21
                    jprod=jprod+ncsf2(l,iarg)*ncsf(jarg)                 3d19s21
                   end do                                                12d19s20
                   ibcoff=iprod                                          3d19s21
                   if(ipss.eq.2)go to 3                                  12d18s20
                  end if                                                 12d18s20
                 end if                                                 4d20s21
                end do                                                  12d18s20
    3           continue                                                12d18s20
               end if                                                   12d18s20
               end if                                                   12d21s22
               if(n2e.gt.0)then                                         8d15s21
                j2o=ibclr(j2o,norbx)                                      12d18s20
                j2o=ibset(j2o,norbxxx)                                     12d18s20
               gandcc=ieor(i1c,j2c)                                     12d21s22
               gandco=ieor(i1o,j2o)                                     12d21s22
               gandcb=ior(gandcc,gandco)                                12d21s22
               ndifb=popcnt(gandcb)                                     12d21s22
               if(ndifb.le.4)then                                       12d21s22
                ndifd=popcnt(gandcc)                                     10d13s22
                ndifs=popcnt(gandco)                                     10d13s22
                 nnot=0
                 if(ndifs.eq.4.and.ndifb.eq.4)then                                           10d14s2
                  nnot=4                                                          10d14s22
                  ioxx(1)=1                                                     10d17s22
                  ioxx(2)=1                                                     10d17s22
                  do i=1,norbxxx                                                     10d17s22
                   if(btest(gandcb,i))then                                         10d14s22
                    if((btest(j2c,i).and.btest(i1o,i)).or.                         10d17s22
     $                 (btest(j2o,i).and..not.btest(i1c,i)))then                   10d14s22
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
     $                 ((btest(i1c,i).and..not.btest(j2o,i)).or.                 10d14s22
     $                 (btest(j2c,i).and..not.btest(i1o,i))))then                     10d14s22
                     if(btest(j2c,i))iswap=1                                        10d17s22
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
                   if(btest(i1c,nab4(1,2)).and.                         10d26s22
     $                  .not.btest(i1c,nab4(1,1)))nbt=1                 10d26s22
                  else                                                             10d17s22
                   nbt=0                                                           10d17s22
                   if(btest(j2c,nab4(2,2)).and.                         10d26s22
     $                  .not.btest(j2c,nab4(2,1)))nbt=1                 10d26s22
                  end if                                                           10d17s22
                  if(nbt.ne.0)then                                                 10d17s22
                   nab4(1,1)=nab4(1,2)                                             10d17s22
                   nab4(2,1)=nab4(2,2)                                             10d17s22
                  end if                                                           10d17s22
                 else if(ndifs.eq.0.and.ndifd.eq.2)then                            10d14s22
                  nnot=3
                  do i=1,norbxxx
                   if(btest(gandcb,i))then                                         10d14s22
                    if(btest(i1c,i))then
                     nab4(1,1)=i
                     nab4(1,2)=i
                    else                                                           10d14s22
                     nab4(2,1)=i
                     nab4(2,2)=i
                    end if
                   end if                                                          10d14s22
                  end do
                 end if                                                            10d14s22
                 ipssx=0                                                 10d21s22
                 if(nnot.eq.3)then                                          12d8s20
                  ipssx=1                                                   12d8s20
                 else if(nnot.eq.4)then                                 12d21s22
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
                  itestc=i1c                                              12d8s20
                  itesto=i1o                                              12d8s20
                  if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                   itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                   itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopen1+1                                              11d13s20
                   karg=jarg-1                                             12d8s20
                  else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                   itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopen1-1                                              11d13s20
                   karg=jarg
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
                  nqq=karg+mdon-1                                        4d20s21
                  if(nqq.ge.mdon.and.nqq.le.mdoo)then                    4d20s21
                   call gandc(i1c,i1o,itestc,itesto,nopen1,nopenk,         12d18s20
     $         jarg,karg,ncsf,norbxxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,   12d18s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
                   call gandc(itestc,itesto,j2c,j2o,nopenk,nopen2p,         12d8s20
     $         karg,iarg,ncsf,norbxxx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,    11d13s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
                   if(nnot1.eq.2.and.nnot2.eq.2)then                       12d19s20
                    if(nab1(2).eq.norbxx)then                              1d25s21
                     iprod=ibcoff                                         3d19s21
                     if(lchoice)then                                      3d19s21
                      call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),       3d19s21
     $               ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod,bc,11d10s22
     $                    ibc)                                          11d10s22
                     end if                                               3d19s21
                     lsa=ism(nab1(1))                                       12d19s20
                     lga=irel(nab1(1))-1                                    12d19s20
                     idx=1                                                 12d19s20
                     jprod=iprod                                          3d19s21
                     iargo=1                                               2d25s21
                     do l=1,4                                               12d19s20
                      if(nl(l).gt.0)then                                     12d18s20
                       nnl=ncsf(jarg)*nfdat(2,l,isb)                         12d19s20
                       iad1=jvcv+iff22(jvcv+5+l)                        8d16s21
                       iad2=iad1+nl(l)                                    3d19s21
                       itmp=ibcoff                                           12d18s20
                       ibcoff=itmp+ncsf(jarg)*nl(l)                          12d18s20
                       call enough('hcdsbk.  8',bc,ibc)
                       if(lchoice)then                                    3d19s21
                        call dgemm('n','n',ncsf(jarg),nl(l),             4d20s21
     $                     ncsf2(l,iarg),                               4d20s21
     $            1d0,bc(jprod),ncsf(jarg),ff22(iad2),ncsf2(l,iarg),0d0,8d16s21
     $                 bc(itmp),ncsf(jarg),                             12d18s20
     d' hcdsbk.  4')
                       else                                               3d19s21
                        call genmatn2(ncsf(jarg),ncsf(karg),ncsf(iarg),          12d18s20
     $                      ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,         2d25s21
     $                   ff22(iad2),ncsf2(l,iarg),nl(l),bc(itmp),0d0,   8d16s21
     $                   iargo,ncsf2(l,iarg),1d0,bc,ibc)                11d10s22
                       end if                                             3d19s21
                       jden=id3vnotv3(l,idx)+nnl*lga                        12d22s20
                       ibc(nd3vnotv3(l,idx)+lga)=1                          12d22s20
                       mden=md3vnotv3(l,idx)                                12d22s20
                       jtmp=itmp                                              12d18s20
                       do iii=0,nl(l)-1                                         12d18s20
                        ii=iff22(iad1+iii)-1                            8d16s21
                        jdh=jden+ncsf(jarg)*ii                              12d18s20
                        ibc(mden+ii)=1                                      12d19s20
                        do j=0,ncsf(jarg)-1                                   12d18s20
                         bc(jdh+j)=bc(jdh+j)+bc(jtmp+j)                       12d18s20
                        end do                                                12d18s20
                        jtmp=jtmp+ncsf(jarg)                                  12d18s20
                       end do                                                 12d18s20
                       ibcoff=itmp                                          12d19s20
                      end if                                                12d19s20
                      jprod=jprod+ncsf(jarg)*ncsf2(l,iarg)                3d19s21
                      iargo=iargo+ncsf2(l,iarg)                            2d25s21
                     end do                                                 12d19s20
                     ibcoff=iprod                                         3d19s21
                    end if                                                 2d25s21
                    if(ipss.eq.2)go to 4                                   12d18s20
                   end if                                                  12d18s20
                  end if                                                 4d20s21
                 end do                                                   12d18s20
    4            continue                                                 12d18s20
                end if                                                    1d6s21
               end if                                                   8d15s21
              end if                                                    3d2s21
              jvcv=jvcv+nspace                                          3d19s21
             end do                                                     12d18s20
c
             if(njhere.gt.0)then                                        3d2s21
              do l=1,4
               nok4f(l)=0                                               11d17s22
               nok4(l)=0                                                11d17s22
               do isc=1,nsymb                                           11d17s22
                do isd=1,nsymb                                          11d17s22
                 nokdc(isd,isc,l)=0                                     11d17s22
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
                 nnl=ncsf(jarg)*nfdat(2,l,isb)                            12d19s20
                 jdhvnotv=idhvnotv(l)                                    12d21s20
                 jdhvnotvf=idhvnotvf(l)                                  1d6s21
                 do i=0,irefo(lgoal)-1                                   12d21s20
                  if(ibc(ndhvnotv(l)+i).ne.0)then                        12d21s20
                   ibc(ndhvnotv(l)+nok4(l))=i                            12d21s20
                   nok4(l)=nok4(l)+1                                     12d21s20
                   iad=idhvnotv(l)+nnl*i
                   do k=0,nok4f(l)-1                                     12d21s20
                    iad=idhvnotv(l)+ncsf(jarg)*ibc(mdhvnotv(l)+k)+nnl*i  1d6s21
                    do j=0,ncsf(jarg)-1                                  1d6s21
                     bc(jdhvnotvf+j)=bc(iad+j)                           1d6s21
                    end do                                               12d21s20
                    jdhvnotvf=jdhvnotvf+ncsf(jarg)                       1d6s21
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
                   nnn=nn*irefo(iscdv)                                       12d18s20
                   nokdc(isd,isc,l)=0                                     12d21s20
                   jd1vnotv=id1vnotv(l,isd,isc)                           12d21s20
                   do i=0,nnn-1
                    if(ibc(nd1vnotv(l,isd,isc)+i).ne.0)then
                     icol=i                                             12d29s20
                     ibc(nd1vnotv(l,isd,isc)+nokdc(isd,isc,l))=icol      12d29s20
                     nokdc(isd,isc,l)=nokdc(isd,isc,l)+1                 12d21s20
                     do k=0,nok4f(l)-1                                    12d21s20
                      iad=id1vnotv(l,isd,isc)+i11s-1                      12d21s20
     $                   +ncsf(jarg)*ibc(mdhvnotv(l)+k)+nnl*i           12d21s20
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
              call enough('hcdsbk.  9',bc,ibc)
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
              if(n2e.gt.0)then                                          8d15s21
               do i=0,irefo(lgoal3)-1                                     8d13s21
                do l=1,4                                                 1d4s21
                 if(nfdat(2,l,isb).gt.0)then                             1d4s21
                  if(ibc(nd3vnotv3(l,ipass)+i).ne.0)then                 1d4s21
                   do k=0,nfdat(2,l,isb)-1                               1d4s21
                    if(ibc(md3vnotv3(l,ipass)+k).ne.0)then               1d4s21
                     ibc(ndkeep(ipass)+nok)=i                            1d5s21
                     ibc(mdkeep(ipass)+nok)=k+loff(l)                    1d5s21
                     iad=id3vnotv3(l,ipass)+ncsf(jarg)*(k                1d4s21
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
              end if                                                    8d15s21
              if(nok.gt.0)then                                          3d2s21
               nnk=nvirt(ksbv)*nrootu                                    8d13s21
               itmpk=ibcoff                                             2d3s21
               itmpd=itmpk+nnk*nok                                       2d3s21
               ibcoff=itmpd+njhere*nok                                  3d2s21
               call enough('hcdsbk. 10',bc,ibc)
               jvss=jvs+nnk*(i11s-1)                                     3d2s21
               call dgemm('n','n',nnk,nok,njhere,1d0,                   8d13s21
     $               bc(jvss),nnk,bc(idkeep(ipass)),njhere,0d0,         8d13s21
     $               bc(itmpk),nnk,                                     8d13s21
     d' hcdsbk.  5')
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
                 call enough('hcdsbk. 11',bc,ibc)
                 fact=0d0                                                  12d21s20
                 if(min(nok4(1),nok4f(1)).gt.0)then                      2d25s21
                  nok=nok4(1)                                            2d25s21
                  nokf=nok4f(1)                                          2d25s21
                  itmp=ibcoff                                              12d21s20
                  ibcoff=itmp+nok*nvirt(ksbv)                           8d13s21
                  call enough('hcdsbk. 12',bc,ibc)
                  jtmp=itmp                                                12d21s20
                  iosym=multh(ksbv,isymop(lxmt))                        8d18s21
                  do iv=0,nvirt(ksbv)-1                                 8d13s21
                   ivp=iv+idoubo(ksbv)+irefo(ksbv)                      8d13s21
                   do i=0,nok-1                                            12d21s20
                    iadh=ixmtf(iosym)+ibc(ndhvnotv(1)+i)+idoubo(iosym)  8d13s21
     $                  +nbasdws(iosym)*ivp                             8d13s21
                    bc(jtmp+i)=bc(iadh)*phase1                          8d16s21
                   end do                                                  12d21s20
                   jtmp=jtmp+nok                                           12d21s20
                  end do                                                   12d21s20
                  call dgemm('n','n',nrow,nvirt(ksbv),nok,sr2,          8d13s21
     $              bc(idhvnotv(1)),nrow,bc(itmp),nok,fact,bc(intden),  2d25s21
     $                nrow,                                             2d25s21
     d' hcdsbk.  6')
                  fact=1d0                                                 12d21s20
                  ibcoff=itmp                                              12d21s20
                 end if                                                    12d21s20
                 nok=0                                                     12d21s20
                 do isc=1,nsymb
                  iscv=multh(isc,ksbv)                                  8d13s21
                  do isd=1,nsymb                                        8d16s21
                   iscdv=multh(iscv,isd)                                    12d19s20
                   nn=irefo(isd)*irefo(isc)                                 12d18s20
                   nnn=nn*irefo(iscdv)                                       12d18s20
                   if(min(nokdc(isd,isc,1),nok4f(1)).gt.0.and.n2e.gt.0) 8d15s21
     $                  then                                            8d15s21
                    nok=nokdc(isd,isc,1)                                 2d25s21
                    nokf=nok4f(1)                                        2d25s21
                    ncol=nok*nokf                                            12d21s20
c     (d|c)(scdv|ksbv)
                    itmp=ibcoff                                            12d21s20
                    ibcoff=itmp+nok*nvirt(ksbv)                            12d21s20
                    call enough('hcdsbk. 13',bc,ibc)
                    do iz=itmp,ibcoff-1                                 8d15s21
                     bc(iz)=0d0                                         8d15s21
                    end do                                              8d15s21
                    irdirc=irefo(isd)*irefo(isc)                        8d15s21
                    do ii2e=1,n2e                                       8d15s21
                     if(multh(isd,isc).eq.isymop(i2eop(1,ii2e)).and.    8d15s21
     $                  multh(iscdv,ksbv).eq.isymop(i2eop(2,ii2e)))then 8d15s21
                      do i=0,nok-1                                      8d15s21
                       jtmp=itmp+i                                      8d15s21
                       idv=ibc(nd1vnotv(1,isd,isc)+i)/irdirc            8d15s21
                       idc=ibc(nd1vnotv(1,isd,isc)+i)-irdirc*idv        8d18s21
                       ic=idc/irefo(isd)                                8d18s21
                       id=idc-irefo(isd)*ic                             8d18s21
                       idc=ixmt(isd,i2eop(1,ii2e))+id+idoubo(isd)       8d15s21
     $                      +nbasdws(isd)*(ic+idoubo(isc))              8d15s21
                       ff=bc(idc)*phase2                                8d15s21
                       idvk=ixmt(iscdv,i2eop(2,ii2e))+idv+idoubo(iscdv) 8d15s21
     $                      +nbasdws(iscdv)*(idoubo(ksbv)+irefo(ksbv))  8d15s21
                       do kv=0,nvirt(ksbv)-1                            8d15s21
                        bc(jtmp+kv*nok)=bc(jtmp+kv*nok)                 8d15s21
     $                       +ff*bc(idvk+nbasdws(iscdv)*kv)             8d15s21
                       end do                                           8d15s21
                      end do                                            8d15s21
                     end if                                             8d15s21
                    end do                                              8d15s21
                    call dgemm('n','n',nrow,nvirt(ksbv),nok,sr2,         2d25s21
     $           bc(id1vnotv(1,isd,isc)),nrow,bc(itmp),nok,fact,        2d25s21
     $                bc(intden),nrow,                                  12d21s20
     d' hcdsbk.  7')
                    fact=1d0                                                 12d21s20
                    ibcoff=itmp                                              12d21s20
                   end if                                                  12d21s20
                  end do                                                   12d21s20
                 end do                                                    12d21s20
                 if(fact.gt.0.5d0)then                                     12d21s20
                  itrans=ibcoff                                            12d21s20
                  ibcoff=itrans+nvirt(ksbv)*nrow                           12d21s20
                  call enough('hcdsbk. 14',bc,ibc)
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
                    call enough('hcdsbk. 15',bc,ibc)
                    fact=0d0                                              12d21s20
                    if(nok4(l).gt.0)then                                 1d1s21
c     Gr jv j = D j k i H0iv2 V jv iv2 r k
                     iosym=multh(isbvu,isymop(lxmt))                    8d16s21
                     itmp=ibcoff                                          12d21s20
                     ibcoff=itmp+nok4(l)*nvirt(isbvu)                     12d21s20
                     call enough('hcdsbk. 16',bc,ibc)
                     jtmp=itmp                                            12d21s20
                     do iv=0,nvirt(isbvu)-1                               12d21s20
                      ivp=iv+idoubo(isbvu)+irefo(isbvu)                 8d16s21
                      ih0=ixmtf(isbvu)+ivp+nbasdws(isbvu)*idoubo(iosym) 8d16s21
                      do i=0,nok4(l)-1                                    12d21s20
                       ii=ibc(ndhvnotv(l)+i)                              12d22s20
                       bc(jtmp+i)=bc(ih0+nbasdws(isbvu)*ii)*phase1      8d16s21
                      end do                                              12d21s20
                      jtmp=jtmp+nok4(l)                                   12d21s20
                     end do                                               12d21s20
                     call dgemm('n','n',nrow,nvirt(isbvu),nok4(l),factt,  12d21s20
     $                  bc(idhvnotv(l)),nrow,bc(itmp),nok4(l),fact,     12d23s20
     $                  bc(intden),nrow,                                12d21s20
     d' hcdsbk.  8')
                     ibcoff=itmp                                          12d21s20
                     fact=1d0                                             12d21s20
                    end if                                                12d21s20
                    if(n2e.gt.0)then                                    8d16s21
                     do isc=1,nsymb                                        12d22s20
                      iscv=multh(isc,lgoal)                                12d22s20
                      do isd=1,nsymb                                    8d16s21
                       irdrc=irefo(isd)*irefo(isc)                      8d16s21
                       iscdv=multh(iscv,isd)                               12d22s20
                       if(nokdc(isd,isc,l).gt.0)then                       12d22s20
                        nn=irefo(isc)*irefo(isd)                          12d22s20
                        nnn=nn*irefo(iscdv)                                12d22s20
                        itmp=ibcoff                                        12d22s20
                        ibcoff=itmp+nokdc(isd,isc,l)*nvirt(isbvu)          12d22s20
                        call enough('hcdsbk. 17',bc,ibc)
                        do iz=itmp,ibcoff-1                              8d16s21
                         bc(iz)=0d0                                      8d16s21
                        end do                                           8d16s21
                        nhit=0
                        do ii2e=1,n2e                                   8d16s21
                         if(multh(isd,isc).eq.isymop(i2eop(1,ii2e)).and.8d16s21
     $                  multh(iscdv,isbvu).eq.isymop(i2eop(2,ii2e)))then8d16s21
                          nhit=1
                          do i=0,nokdc(isd,isc,l)-1                     8d16s21
                           idv=ibc(nd1vnotv(l,isd,isc)+i)/irdrc         8d16s21
                           idc=ibc(nd1vnotv(l,isd,isc)+i)-irdrc*idv     8d16s21
                           ic=idc/irefo(isd)                            8d16s21
                           id=idc-irefo(isd)*ic                         8d16s21
                           idcx=ixmt(isd,i2eop(1,ii2e))+id+idoubo(isd)   8d16s21
     $                          +nbasdws(isd)*(ic+idoubo(isc))          8d18s21
                           idvk=ixmt(iscdv,i2eop(2,ii2e))+idv           8d16s21
     $                          +idoubo(iscdv)                          8d16s21
     $                      +nbasdws(iscdv)*(idoubo(isbvu)+irefo(isbvu))8d16s21
                           ff=bc(idcx)*phase2                            8d16s21
                           jtmp=itmp+i                                     3d23s21
                           do iv=0,nvirt(isbvu)-1                       8d16s21
                            bc(jtmp+iv*nokdc(isd,isc,l))=               8d16s21
     $                      bc(jtmp+iv*nokdc(isd,isc,l))                8d16s21
     $                           +ff*bc(idvk+iv*nbasdws(iscdv))         8d16s21
                           end do                                       8d16s21
                          end do                                        8d16s21
                         end if                                         8d16s21
                        end do                                          8d16s21
                        if(nhit.gt.0)then                               8d16s21
                         call dgemm('n','n',nrow,nvirt(isbvu),
     $                   nokdc(isd,isc,l),factt,bc(id1vnotv(l,isd,isc)),12d22s20
     $                    nrow,bc(itmp),nokdc(isd,isc,l),fact,          12d22s20
     $                    bc(intden),nrow,                              12d22s20
     d' hcdsbk.  9')
                         fact=1d0                                           12d22s20
                        end if                                          8d16s21
                        ibcoff=itmp                                        12d22s20
                       end if                                              12d22s20
                      end do                                               12d22s20
                     end do                                                12d22s20
                    end if                                              8d16s21
                    if(fact.gt.0.5d0)then                                 12d22s20
                     itrans=ibcoff                                        12d22s20
                     ibcoff=itrans+nrow*nvirt(isbvu)                      12d22s20
                     call enough('hcdsbk. 18',bc,ibc)
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
                     call enough('hcdsbk. 19',bc,ibc)
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
                      call enough('hcdsbk. 20',bc,ibc)
                      call dgemm('n','n',mrow,mcol,njhere,1d0,          3d2s21
     $                     bc(itrans),mrow,bc(itmpsv),njhere,0d0,       3d2s21
     $                     bc(itmpdg),mrow,                             3d2s21
     d' hcdsbk. 10')
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
                      call enough('hcdsbk. 21',bc,ibc)
                      call dgemm('n','n',mcol,mrow,njhere,1d0,          3d2s21
     $                     bc(itrans),mcol,bc(itmpsv),njhere,0d0,       3d2s21
     $                     bc(itmpgv),mcol,                             3d2s21
     d' hcdsbk. 11')
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
             jvs=jvs+ncsf(jarg)*nvirt(ksbv)*nrootu                      8d13s21
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
       if(n2e.gt.0)then                                                 8d16s21
        ioffdnon=1                                                       2d4s21
        call dws_gsumf(bc(ibcvmat),nbmat)                                3d2s21
        do isb=1,nsymb                                                   2d4s21
         nfh=nfdat(2,1,isb)+nfdat(2,2,isb)+nfdat(2,3,isb)+nfdat(2,4,isb) 2d4s21
         isbv12=multh(isb,isymbra)                                      8d16s21
         isn=multh(isbv12,ksbv)                                         8d27s21
         if(min(irefo(isn),nvirt(ksbv)).gt.0)then                        2d25s21
          nn=nvirt(ksbv)*nrootu*nfh                                       2d4s21
          ipsr=ibcoff                                                    2d4s21
          ibcoff=ipsr+nrootu*nfh                                         2d4s21
          call enough('hcdsbk. 22',bc,ibc)
          do i=0,nfdat(2,1,isb)*nrootu-1                                 2d4s21
           bc(ipsr+i)=1d0                                                2d4s21
          end do                                                         2d4s21
          do i=nfdat(2,1,isb)*nrootu,nfh*nrootu-1                        2d5s21
           bc(ipsr+i)=-1d0                                               2d4s21
          end do                                                         2d4s21
c
c     form Gv'v"rk=[(vv"|nv')+p(k)(vv'|nv")]vvrkn
c
          do isbv1=1,nsymb                                               2d4s21
           isbv2=multh(isbv1,isbv12)                                     2d4s21
           if(isbv2.ge.isbv1)then                                        2d4s21
            call ilimts(irefo(isn),nvirt(isbv1),mynprocg,mynowprog,il,
     $          ih,i1s,i1e,i2s,i2e)                                         2d4s21
            if(isbv1.eq.isbv2)then                                       2d4s21
             do ii2e=1,n2e                                              8d16s21
              if(multh(ksbv,isbv2).eq.isymop(i2eop(1,ii2e)).and.        8d16s21
     $           multh(isn,isbv1).eq.isymop(i2eop(2,ii2e)))then         8d16s21
               i10=i1s                                                     2d4s21
               i1n=irefo(isn)                                              2d4s21
               do i2=i2s,i2e                                              2d4s21
                iv1=i2-1                                                8d16s21
                iv1p=iv1+idoubo(isbv1)+irefo(isbv1)                     8d16s21
                if(i2.eq.i2e)i1n=i1e                                      2d4s21
                do i1=i10,i1n                                             2d4s21
                 i1m=i1-1
                 ixu=ixmt(isbv1,i2eop(2,ii2e))+iv1p+nbasdws(isbv1)*     8d16s21
     $                (i1m+idoubo(isn))                                 8d16s21
                 ff=bc(ixu)*phase2*sr2                                  8d16s21
                 ixv=ixmt(isbv2,i2eop(1,ii2e))+iv1p                     8d18s21
     $                +nbasdws(isbv1)*(idoubo(ksbv)+irefo(ksbv))        8d18s21
                 do k=0,nfdat(2,1,isb)-1                                  2d4s21
                  do ir=0,nrootm                                          2d4s21
                   sum=0d0                                                2d4s21
                   jvmat=ivmat(isb)+nvirt(ksbv)*(ir+nrootu*(k+nfh*i1m))   2d4s21
                   iadv=ioffdnon+iv1+nvirt(isbv1)*(ir+nrootu*k)           2d4s21
                   do kv=0,nvirt(ksbv)-1                                    2d4s21
                    sum=sum+bc(jvmat+kv)*bc(ixv+nbasdws(isbv1)*kv)      8d16s21
                   end do                                                 2d4s21
                   sum=sum*ff                                           8d16s21
                   gd(iadv)=gd(iadv)+sum                                  2d4s21
                  end do                                                  2d4s21
                 end do                                                   2d4s21
                end do                                                    2d4s21
                i10=1                                                     2d4s21
               end do                                                     2d4s21
              end if                                                    8d16s21
             end do                                                     8d16s21
             ioffdnon=ioffdnon+nvirt(isbv1)*nrootu*nfdat(2,1,isb)        2d4s21
             nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       2d4s21
             isw=0                                                       2d4s21
            else                                                         2d4s21
             nvv=nvirt(isbv1)*nvirt(isbv2)                               2d4s21
             isw=1                                                       2d4s21
            end if                                                       2d4s21
            do ii2e=1,n2e                                               8d16s21
             if(multh(ksbv,isbv2).eq.isymop(i2eop(1,ii2e)).and.         8d16s21
     $           multh(isn,isbv1).eq.isymop(i2eop(2,ii2e)))then         8d16s21
              i10=i1s                                                      2d4s21
              i1n=irefo(isn)                                               2d4s21
              do i2=i2s,i2e                                               2d4s21
               iv1=i2-1                                                 8d16s21
               iv1p=iv1+idoubo(isbv1)+irefo(isbv1)                      8d16s21
               ibots=i2                                                   2d4s21
               ibotn=0                                                    2d4s21
               ibot=ibots+isw*(ibotn-ibots)                               2d4s21
               if(i2.eq.i2e)i1n=i1e                                       2d4s21
               do i1=i10,i1n                                              2d4s21
                i1m=i1-1                                                8d16s21
                ixv=ixmt(isbv1,i2eop(2,ii2e))+iv1p+nbasdws(isbv1)       8d16s21
     $               *(i1m+idoubo(isn))                                 8d16s21
                ff=bc(ixv)*phase2                                       8d16s21
                do kr=0,nfh*nrootu-1                                      2d4s21
                 jvmat=ivmat(isb)+nvirt(ksbv)*kr+nn*i1m                   2d4s21
                 do iv2=ibot,nvirt(isbv2)-1                                2d4s21
                  ixv=ixmt(isbv2,i2eop(1,ii2e))+iv2+idoubo(isbv2)       8d16s21
     $           +irefo(isbv2)+nbasdws(isbv2)*(idoubo(ksbv)+irefo(ksbv))8d16s21
                  sum=0d0                                                 2d4s21
                  itri=((iv2*(iv2-1))/2)+iv1                              2d4s21
                  irec=iv1+nvirt(isbv1)*iv2                               2d4s21
                  iadv=ioffdnon+itri+isw*(irec-itri)+nvv*kr               2d4s21
                  do kv=0,nvirt(ksbv)-1                                    2d4s21
                   sum=sum+bc(jvmat+kv)*bc(ixv+kv*nbasdws(isbv2))       8d16s21
                  end do                                                  2d4s21
                  gd(iadv)=gd(iadv)+sum*ff                              8d16s21
                 end do                                                   2d4s21
                end do                                                    2d4s21
               end do                                                     2d4s21
               i10=1                                                      2d4s21
              end do                                                      2d4s21
             end if                                                     8d16s21
            end do                                                      8d16s21
            call ilimts(irefo(isn),nvirt(isbv2),mynprocg,mynowprog,il,  8d16s21
     $           ih,i1s,i1e,i2s,i2e)                                    8d16s21
            do ii2e=1,n2e                                               8d16s21
             if(multh(isn,isbv2).eq.isymop(i2eop(1,ii2e)).and.          8d16s21
     $          multh(ksbv,isbv1).eq.isymop(i2eop(2,ii2e)))then         8d16s21
              i10=i1s                                                      2d4s21
              i1n=irefo(isn)                                               2d4s21
              do i2=i2s,i2e                                               2d4s21
               iv2=i2-1                                                   2d4s21
               itops=iv2-1                                                2d4s21
               itopn=nvirt(isbv1)-1                                       2d4s21
               itop=itops+isw*(itopn-itops)                               2d4s21
               if(i2.eq.i2e)i1n=i1e                                       2d4s21
               iv2p=iv2+idoubo(isbv2)+irefo(isbv2)                      8d16s21
               do i1=i10,i1n                                              2d4s21
                i1m=i1-1                                                  2d4s21
                ixv=ixmt(isbv2,i2eop(1,ii2e))+iv2p+nbasdws(isbv2)       8d16s21
     $               *(i1m+idoubo(isn))                                 8d16s21
                ff=bc(ixv)*phase2                                       8d16s21
                do kr=0,nfh*nrootu-1                                      2d4s21
                 jvmat=ivmat(isb)+nvirt(ksbv)*kr+nn*i1m                   2d4s21
                 do iv1=0,itop                                            2d4s21
                  sum=0d0                                                 2d4s21
                  itri=((iv2*(iv2-1))/2)+iv1                              2d4s21
                  irec=iv1+nvirt(isbv1)*iv2                               2d4s21
                  iadv=ioffdnon+itri+isw*(irec-itri)+nvv*kr               2d4s21
                  ixv=ixmt(isbv1,i2eop(2,ii2e))+iv1+idoubo(isbv1)       8d16s21
     $           +irefo(isbv1)+nbasdws(isbv1)*(idoubo(ksbv)+irefo(ksbv))8d16s21
                  do kv=0,nvirt(ksbv)-1                                    2d4s21
                   sum=sum+bc(jvmat+kv)*bc(ixv+kv*nbasdws(isbv1))       8d16s21
                  end do                                                  2d4s21
                  gd(iadv)=gd(iadv)+ff*sum*bc(ipsr+kr)                  8d16s21
                 end do                                                   2d4s21
                end do                                                    2d4s21
               end do                                                     2d4s21
               i10=1                                                      2d4s21
              end do                                                      2d4s21
             end if                                                     8d16s21
            end do                                                      8d16s21
            ioffdnon=ioffdnon+nvv*nrootu*nfh                             2d4s21
           end if                                                        2d4s21
          end do                                                         2d4s21
         else                                                            2d5s21
          do isbv1=1,nsymb                                               2d5s21
           isbv2=multh(isbv1,isbv12)                                     2d5s21
           if(isbv2.ge.isbv1)then                                        2d5s21
            if(isbv1.eq.isbv2)then                                       2d5s21
             ioffdnon=ioffdnon+nvirt(isbv1)*nrootu*nfdat(2,1,isb)        2d5s21
             nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       2d5s21
            else                                                         2d5s21
             nvv=nvirt(isbv1)*nvirt(isbv2)                               2d5s21
            end if                                                       2d5s21
            ioffdnon=ioffdnon+nvv*nrootu*nfh                             2d5s21
           end if                                                        2d5s21
          end do                                                         2d4s21
         end if                                                          2d4s21
        end do                                                           2d4s21
       end if                                                           8d16s21
       ibcoff=ibcbmat                                                   2d3s21
      end do                                                            12d18s20
      ibcoff=ircv                                                       1d30s21
      return
      end
