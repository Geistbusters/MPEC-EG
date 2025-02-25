c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcdsd1(ihsdiag,nff1,iff1,nff22,nfdat,vd,nsymb,mdon,    3d10s21
     $     mdoo,nec,multh,isymmrci,nvirt,ncsf,ncsf2,irel,ism,irefo,     12d18s20
     $     ixw1,ixw2,norb,nrootu,iden,nh0av,maxbx,sr2,iff22,ff22,bc,ibc,1d27s23
     $     iroff)                                                       1d27s23
      implicit real*8 (a-h,o-z)                                         12d18s20
c
c     to do ...
c     consolidate densities, and perhaps njhere in density calculation
c     did I swap iarg and jarg in 4th gandc4?
c
      integer*8 ihsdiag(mdoo+1,nsymb,2),i18,i28,i38,i48,i1c,i1o,j2c,    12d18s20
     $     j2o,itestc,itesto,ipack8,last8(2),iff22(*),gandcc,gandco,    2d6s23
     $     gandcb                                                       2d6s23
      integer*1 nab1(2),nab2(2)                                         12d18s20
      external second                                                   2d18s21
      logical l3x,lprt,lnew,ldebug                                      1d25s21
      equivalence (ipack8,ipack4)                                       12d18s20
      dimension nff1(mdoo+1,nsymb,2),iff1(*),nff22(mdoo+1,2,nsymb),     12d18s20
     $     nfdat(5,4,*),multh(8,8),nvirt(*),ncsf(*),ncsf2(4,*),         3d11s21
     $     irel(*),ism(*),irefo(*),vd(*),itest(32,3),                   3d11s21
     $     nab4(2,3),ipack4(2),nl(4),npre(4),mpre(4),                   3d11s21
     $     id1visv(8,8),nd1visv(8,8),nokdc(8,8,4),nok3v(4),             12d22s20
     $     idhvnotv(4),ndhvnotv(4),id1vnotv(4,8,8),nd1vnotv(4,8,8),     12d21s20
     $     id3vnotv3(4,2),nd3vnotv3(4,2),loff(4),idkeep(2),ndkeep(2),   1d4s21
     $     mdkeep(2),keep(2),idhvnotvf(4),ibmat(8),ivmat(8),            2d4s21
     $     mdhvnotv(4),md3vnotv3(4,2),nok4f(4),nok4(4),nok3vf(4),       12d22s20
     $     iden(*),nh0av(*),nok33f(4,2),nok33(4,2),tdends(6),ff22(*)    3d21s22
      include "common.store"
      common/kmfind/invk1(2,8,8,8,2)                                    6d30s18
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      nrootm=nrootu-1                                                   1d4s21
      idoit=0                                                           3d1s21
      lnew=.true.
      last8(1)=-1                                                       2d8s21
      ircv=ibcoff                                                       1d29s21
      ivs=ircv+mynprocg                                                 3d11s21
      ibcoff=ivs+maxbx                                                  6d17s21
      nacc=0                                                            1d30s21
      call enough('hcdsd1.  1',bc,ibc)
      itransgg=0                                                        1d30s21
      nsing=0                                                           12d23s20
      norbx=norb+1                                                      12d18s20
      norbxx=norbx+1                                                    12d18s20
      norbxxx=norbxx+1                                                  12d18s20
      do jsb=1,nsymb                                                    12d18s20
       jsbv=multh(jsb,isymmrci)                                         12d18s20
       do nclo1p=mdon+1,mdoo+1                                          12d18s20
        if(min(nff1(nclo1p,jsb,1),nvirt(jsbv)).gt.0)then                12d18s20
         nclo1=nclo1p-1                                                  12d12s20
         jarg=nclo1p-mdon                                                12d12s20
         nopen1=nec-2*nclo1                                              12d12s20
         ngg=nvirt(jsbv)*nrootu                                         12d15s20
         nggg=nff1(nclo1p,jsb,1)*ncsf(jarg)                             12d29s20
         nsing=nsing+nggg*nvirt(jsbv)                                   12d29s20
         ncolt=nggg*ngg                                                 12d29s20
          call ilimts(1,nff1(nclo1p,jsb,1),mynprocg,mynowprog,          3d2s21
     $        i1l,i1h,i11s,i11e,i12s,i12e)                              12d19s20
          i1l=1+ncsf(jarg)*(i1l-1)                                      3d2s21
          i1h=ncsf(jarg)*i1h                                            3d2s21
         ibcgg=ibcoff                                                   1d30s21
         i18=1
         i28=ngg                                                        12d18s20
         i38=1                                                          12d18s20
         i48=nggg                                                       12d29s20
         call ddi_iget(bc,ibc,ihsdiag(nclo1p,jsb,1),i18,i28,i38,i48,    11d15s22
     $        bc(ivs),ibc(ircv),nrcv)                                   11d15s22
         itransvs=0                                                     1d29s21
         ioffdnon=1                                                     12d21s20
         do isb=1,nsymb                                                 12d18s20
          nfh=nfdat(2,1,isb)+nfdat(2,2,isb)+nfdat(2,3,isb)              2d3s21
     $         +nfdat(2,4,isb)                                          2d3s21
          isbv12=multh(isb,isymmrci)                                    12d18s20
          loff(1)=0                                                     1d4s21
          do l=2,4                                                      1d4s21
           lm=l-1                                                       1d4s21
           loff(l)=loff(lm)+nfdat(2,lm,isb)                             1d4s21
          end do                                                        1d4s21
          nnt=loff(4)+nfdat(2,4,isb)                                    1d4s21
          do nclo2p=max(mdon+1,nclo1p-1),min(mdoo+1,nclo1p+2)           3d11s21
           if(nff22(nclo2p,1,isb).gt.0)then                             12d18s20
            nclo2=nclo2p-1                                              12d18s20
            iarg=nclo2p-mdon                                            12d18s20
            iargp=iarg+1                                                12d18s20
            nopen2=nec-2*nclo2p                                         12d18s20
            nopen2p=nopen2+2                                            12d18s20
            idhvisv=ibcoff                                              12d18s20
            nnj=ncsf(jarg)*nfdat(2,1,isb)                               12d18s20
            lgoal=multh(jsb,isb)                                        12d22s20
            lgoal3=multh(jsbv,isbv12)                                   12d22s20
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
            end do                                                      12d18s20
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
             if(ntest.ne.nclo1)then
              write(6,*)('ntest = '),ntest,if1
              call dcbit(i1c,32,'dorb')
              call dws_synca
              call dws_finalize
              stop
             end if
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
             if(popcnt(i1o).ne.popcnt(iff1(if1o)))then
              write(6,*)('i1o popcnt failure! ')
              call dcbit(iff1(if1o),32,'iff1')
              call dcbit(i1o,64,'i1o')
              stop 'hcdsd1'
             end if
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
              ipack8=iff22(jvcv)                                        3d21s22
              j2c=ipack4(1)                                             12d18s20
              if(popcnt(j2c).ne.popcnt(ipack4(1)))then                  1d30s23
               write(6,*)('j2c popcnt failure! ')
               call dcbit(ipack4(1),32,'ipack4(1)')
               call dcbit(j2c,64,'j2c')
               stop 'hcdsd1'
              end if                                                    1d30s23
              j2o=ipack4(2)                                             12d18s20
              if(popcnt(j2o).ne.popcnt(ipack4(2)))then
               write(6,*)('j2o popcnt failure! ')
               call dcbit(ipack4(2),32,'ipack4(2)')
               call dcbit(j2o,64,'j2o')
               stop 'hcdsd1'
              end if
              nclo=popcnt(ipack4(1))                                    12d18s20
              nspace=iff22(jvcv+1)                                      3d21s22
              do l=1,4                                                  3d19s21
               nl(l)=iff22(jvcv+1+l)                                    3d21s22
              end do                                                    3d19s21
              j2o=ibset(j2o,norbx)                                      12d18s20
              j2o=ibset(j2o,norbxx)                                     12d18s20
              if(njhere.gt.0)then                                       2d26s21
               gandcc=ieor(i1c,j2c)                                     10d14s22
               gandco=ieor(i1o,j2o)                                     10d14s22
               gandcb=ior(gandcc,gandco)                                10d21s22
               ndifb=popcnt(gandcb)                                     10d21s22
               if(ndifb.le.2)then                                       2d6s23
                ndifs=popcnt(gandco)                                     10d14s22
                if(ndifs.eq.2.and.ndifb.eq.2)then                       2d6s23
                 do i=1,norbxx
                  if(btest(gandco,i))then                                          10d14s22
                   if((btest(i1o,i).and..not.btest(j2c,i)).or.
     $                 (btest(i1c,i).and.btest(j2o,i)))then                           10d14s22
                    nab4(1,1)=i
                   else                                                           10d14s22
                    nab4(2,1)=i
                   end if
                  end if                                                          10d14s22
                 end do                                                           10d14s22
                 call gandc(i1c,i1o,j2c,j2o,nopen1,nopen2p,jarg,iarg,    12d18s20
     $               ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,       12d18s20
     $               ncsfmid1,bc,ibc)                                   11d14s22
                 if(nab1(2).ne.norbxx)then                               12d19s20
                  write(6,*)('but expecting virt to be norbxx = '),
     $                norbxx
                  stop
                 end if                                                  12d19s20
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
                 do l=1,4                                                12d18s20
                  if(nl(l).gt.0)then                                     12d18s20
                   nnl=ncsf(jarg)*nfdat(2,l,isb)                         12d19s20
                   iad1=jvcv+iff22(jvcv+5+l)                             3d21s22
                   iad2=iad1+nl(l)                                       3d19s21
                   itmp=ibcoff                                           12d18s20
                   ibcoff=itmp+ncsf(jarg)*nl(l)                          12d18s20
                   call enough('hcdsd1.  2',bc,ibc)
                   call xtimesn2(ncsf(jarg),ncsf(iarg),ncsfmid1,nl(l),   2d25s21
     $                 iwpb1,iwpk1,ff22(iad2),ncsf2(l,iarg),bc(itmp),   3d21s22
     $                 ncsf(jarg),1d0,0d0,iargo,ncsf2(l,iarg),bc,ibc)   11d10s22
c     jsbv=jsb*isymmrci
c     isbv12=isb*isymmrci
c     isbv(1or2)=isbv12*jsbv=jsb*isymmrci*isb*isymmrci=jsb*isb
                   lgoal=multh(jsb,isb)                                  12d21s20
                   if(lsa.ne.lgoal)then
                    write(6,*)('whoa cowboy... lsa vs lgoal'),lsa,lgoal
                    write(6,*)('nab1: '),nab1
                    call dws_synca
                    call dws_finalize
                    stop
                   end if
                   jden=idhvnotv(l)+nnl*lga                              12d21s20
                   ibc(ndhvnotv(l)+lga)=1                                12d21s20
                   jtmp=itmp                                             12d18s20
                   do iii=0,nl(l)-1                                      12d18s20
                    ip=iff22(iad1+iii)-1                                 3d21s22
                    jjden=jden+ncsf(jarg)*ip                             12d18s20
                    ibc(mdhvnotv(l)+ip)=1                                12d21s20
                    do j=0,ncsf(jarg)-1                                  12d18s20
                     bc(jjden+j)=bc(jjden+j)+bc(jtmp+j)                  12d18s20
                    end do                                               12d18s20
                    jtmp=jtmp+ncsf(jarg)                                 12d18s20
                   end do                                                12d18s20
                   ibcoff=itmp                                           12d18s20
                  end if                                                 12d18s20
                  iargo=iargo+ncsf2(l,iarg)                              2d25s21
                 end do                                                  12d18s20
                end if                                                  2d6s23
               end if                                                   12d18s20
              end if                                                    3d2s21
              jvcv=jvcv+nspace                                          3d19s21
             end do                                                     12d18s20
             if(njhere.gt.0)then                                        3d2s21
              do l=1,4
               if(nfdat(2,l,isb).gt.0)then                               12d19s20
                nok4f(l)=0                                               12d21s20
                do i=0,nfdat(2,l,isb)-1                                  12d21s20
                 if(ibc(mdhvnotv(l)+i).ne.0)then                         12d21s20
                  ibc(mdhvnotv(l)+nok4f(l))=i                            12d21s20
                  nok4f(l)=nok4f(l)+1                                    12d21s20
                 end if                                                  12d21s20
                end do                                                   12d21s20
                if(nok4f(l).gt.0)then                                    12d21s20
                 nnl=ncsf(jarg)*nfdat(2,l,isb)                            12d19s20
                 lgoal=multh(jsbv,isbv12)                                12d21s20
                 nok4(l)=0                                               12d21s20
                 jdhvnotv=idhvnotv(l)                                    12d21s20
                 do i=0,irefo(lgoal)-1                                   12d21s20
                  if(ibc(ndhvnotv(l)+i).ne.0)then                        12d21s20
                   ibc(ndhvnotv(l)+nok4(l))=i                            12d21s20
                   nok4(l)=nok4(l)+1                                     12d21s20
                   iad=idhvnotv(l)+nnl*i
                   do k=0,nok4f(l)-1                                     12d21s20
                    iad=idhvnotv(l)+ncsf(jarg)*ibc(mdhvnotv(l)+k)+nnl*i  1d6s21
                    iad=iad+i11s-1                                       1d6s21
                    do j=0,njhere-1                                      12d21s20
                     bc(jdhvnotv+j)=bc(iad+j)                            12d21s20
                    end do                                               12d21s20
                    jdhvnotv=jdhvnotv+njhere                             12d21s20
                   end do                                                12d21s20
   10              format(i5,5x,20i2)
                  end if
                 end do
                end if
               end if                                                    12d21s20
              end do                                                     12d21s20
             end if                                                     3d2s21
c
c     I would like to consolidate the ls so I only have to make 1 pass
c     through the 3x ints.
c
             if(itransvs.eq.0)then                                      1d29s21
              call ddi_done(ibc(ircv),nrcv)                                  1d29s21
              itransvs=1                                                1d29s21
              itmp=ibcoff                                                    1d27s21
              ibcoff=itmp+ngg                                                1d27s21
              call enough('hcdsd1.  3',bc,ibc)
              jjvs=ivs                                                        1d27s21
              nggm=ngg-1                                                     1d27s21
              do i=0,nggg-1                                                  1d27s21
               jvs0=jjvs                                                      1d27s21
               do iv=0,nvirt(jsbv)-1                                         1d27s21
                iad=itmp+iv                                                  1d27s21
                do ir=0,nrootm                                               1d27s21
                 bc(iad+ir*nvirt(jsbv))=bc(jjvs+ir)                           1d27s21
                end do                                                       1d27s21
                jjvs=jjvs+nrootu                                               1d27s21
               end do                                                        1d27s21
               do j=0,nggm                                                   1d27s21
                bc(jvs0+j)=bc(itmp+j)                                        1d27s21
               end do                                                        1d27s21
              end do                                                         1d27s21
              ibcoff=itmp                                                    1d27s21
             end if                                                     1d29s21
             if(njhere.gt.0)then                                        3d2s21
              do isbv1=1,nsymb                                           12d21s20
               isbv2=multh(isbv1,isbv12)                                 12d21s20
               if(isbv1.le.isbv2)then                                    12d21s20
                if(isbv12.eq.1.and.jsbv.eq.isbv1.and.                   11d15s22
     $              nfdat(2,1,isb).gt.0)then                            11d15s22
                 if(min(nok4(1),nok4f(1)).gt.0)then                      2d25s21
                  nok=nok4(1)                                            2d25s21
                  nokf=nok4f(1)                                          2d25s21
c     itmp is compressed h0av
c     intden j,d,v=sum o D j,d,o x H o,v
c     trans v,j,d=sumo D j,d,o x H o,v
c     gd v,r,d = sum j Vs v,r,j x sum o D j,d,o x H o,v
c     H o,v -> sum j Vs v,r,j x sum d D j,d,o x Vd v,r,d
c            = sum d ( sum j Vs v,r,j x D j,d,o ) Vd v,r,d
c            = sum d T v,r,d,o Vd v,r,d
                  mrow=nvirt(jsbv)*nrootu                               3d10s21
                  mcol=nok*nokf                                         3d10s21
                  itmp=ibcoff                                           3d10s21
                  ibcoff=itmp+mrow*mcol                                 3d10s21
                  call enough('hcdsd1.  4',bc,ibc)
                  iads=jvs+mrow*(i11s-1)                                3d10s21
                  call dgemm('n','n',mrow,mcol,njhere,sr2,bc(iads),mrow,1d30s23
     $                 bc(idhvnotv(1)),njhere,0d0,bc(itmp),mrow,        3d10s21
     d'hcdsd1.  1')
                  do io=0,nok-1                                         3d10s21
                   iio=ibc(ndhvnotv(1)+io)                              3d10s21
                   do ir=0,nrootm                                       3d10s21
                    iru=ir*iroff                                        1d27s23
                    jden=iden(jsbv)+irefo(jsbv)+nh0av(jsbv)*(iio        3d11s21
     $                   +nh0av(jsbv)*iru)                              1d27s23
                    do if=0,nokf-1                                      3d10s21
                     iadvd=joffdnon+nvirt(jsbv)*(ir                     3d10s21
     $                    +nrootu*ibc(mdhvnotv(1)+if))                  3d10s21
                     jtmp=itmp+nvirt(jsbv)*(ir+nrootu*(if+nokf*io))     3d10s21
                     do iv=0,nvirt(jsbv)-1                              3d10s21
                      bc(jden+iv)=bc(jden+iv)+vd(iadvd+iv)*bc(jtmp+iv)  3d10s21
                     end do                                             3d10s21
                    end do                                              3d10s21
                   end do                                               3d10s21
                  end do                                                3d10s21
                  ibcoff=itmp                                           3d10s21
                 end if                                                    12d21s20
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
                  if(jsbv.eq.isbv1.or.jsbv.eq.isbv2)then                   12d22s20
                   if(isbv1.eq.jsbv)then                                    12d21s20
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
                   if(min(nok4f(l),nok4(l),nvirt(isbvu)).gt.0)then      3d11s21
                    mrow=nvirt(jsbv)*nrootu                             3d11s21
                    mcol=nok4f(l)*nok4(l)                               3d11s21
                    itmp=ibcoff                                         3d11s21
                    ibcoff=itmp+mrow*mcol                               3d11s21
                    call enough('hcdsd1.  5',bc,ibc)
                    jjvs=jvs+mrow*(i11s-1)                              3d11s21
                    call dgemm('n','n',mrow,mcol,njhere,factt,          3d11s21
     $                   bc(jjvs),mrow,bc(idhvnotv(l)),njhere,0d0,      3d11s21
     $                   bc(itmp),mrow,                                 3d11s21
     d'hcdsd1.  2')
                    if(isbv1.eq.isbv2)then                               12d22s20
                     do ir=0,nrootm                                     3d11s21
                      iru=ir*iroff                                      1d27s23
                      do io=0,nok4(l)-1                                 3d11s21
                       iio=ibc(ndhvnotv(l)+io)                          3d11s21
                       jden=iden(isbvu)+irefo(isbvu)+nh0av(isbvu)*(     3d11s21
     $                      iio+nh0av(isbvu)*iru)                       1d27s23
                       do k=0,nok4f(l)-1                                3d11s21
                        kk=ibc(mdhvnotv(l)+k)                           3d11s21
                        iadv=joffdnon+nvv*(ir+nrootu*kk)                3d11s21
                        iprod=itmp+nvirt(isbv2)*(ir+nrootu*(k           3d11s21
     $                       +nok4f(l)*io))                             3d11s21
                        do iv2=0,nvirt(isbv2)-1                         3d11s21
                         do iv1=0,iv2-1                                 3d11s21
                          bc(jden+iv1)=bc(jden+iv1)                     3d11s21
     $                         +bc(iprod+iv2)*vd(iadv+iv1)*tf2          3d11s21
                          bc(jden+iv2)=bc(jden+iv2)+bc(iprod+iv1)       3d11s21
     $                         *vd(iadv+iv1)                            3d11s21
                         end do                                         3d11s21
                         iadv=iadv+iv2                                  3d11s21
                        end do                                          3d11s21
                       end do                                           3d11s21
                      end do                                            3d11s21
                     end do                                             3d11s21
                    else if(isbvu.eq.isbv1)then                          12d22s20
                     do ir=0,nrootm                                     3d11s21
                      iru=ir*iroff                                      1d27s23
                      do io=0,nok4(l)-1                                 3d11s21
                       iio=ibc(ndhvnotv(l)+io)                          3d11s21
                       jden=iden(isbvu)+irefo(isbvu)+nh0av(isbvu)*(     3d11s21
     $                      iio+nh0av(isbvu)*iru)                       1d27s23
                       do k=0,nok4f(l)-1                                3d11s21
                        kk=ibc(mdhvnotv(l)+k)                           3d11s21
                        iadv=joffdnon+nvv*(ir+nrootu*kk)                3d11s21
                        iprod=itmp+nvirt(jsbv)*(ir+nrootu*(k            3d11s21
     $                       +nok4f(l)*io))                             3d11s21
                        do iv2=0,nvirt(isbv2)-1                         3d11s21
                         do iv1=0,nvirt(isbv1)-1                        3d11s21
                          bc(jden+iv1)=bc(jden+iv1)                     3d11s21
     $                         +bc(iprod+iv2)*vd(iadv+iv1)              3d11s21
                         end do                                         3d11s21
                         iadv=iadv+nvirt(isbv1)                         3d11s21
                        end do                                          3d11s21
                       end do                                           3d11s21
                      end do                                            3d11s21
                     end do                                             3d11s21
                    else                                                 12d22s20
                     do ir=0,nrootm                                     3d11s21
                      iru=ir*iroff                                      1d27s23
                      do io=0,nok4(l)-1                                 3d11s21
                       iio=ibc(ndhvnotv(l)+io)                          3d11s21
                       jden=iden(isbvu)+irefo(isbvu)+nh0av(isbvu)*(     3d11s21
     $                      iio+nh0av(isbvu)*iru)                       1d27s23
                       do k=0,nok4f(l)-1                                3d11s21
                        kk=ibc(mdhvnotv(l)+k)                           3d11s21
                        iadv=joffdnon+nvv*(ir+nrootu*kk)                3d11s21
                        iprod=itmp+nvirt(jsbv)*(ir+nrootu*(k            3d11s21
     $                       +nok4f(l)*io))                             3d11s21
                        do iv2=0,nvirt(isbv2)-1                         3d11s21
                         do iv1=0,nvirt(isbv1)-1                        3d11s21
                          bc(jden+iv2)=bc(jden+iv2)                     3d11s21
     $                         +bc(iprod+iv1)*vd(iadv+iv1)              3d11s21
                         end do                                         3d11s21
                         iadv=iadv+nvirt(isbv1)                         3d11s21
                        end do                                          3d11s21
                       end do                                           3d11s21
                      end do                                            3d11s21
                     end do                                             3d11s21
                    end if
                    ibcoff=itmp                                         3d11s21
                   end if                                                 12d21s20
    5              continue                                               12d21s20
                  end if                                                 12d22s20
                  joffdnon=joffdnon+nvv*nfdat(2,l,isb)*nrootu            12d22s20
                 end if                                                  12d21s20
                end do                                                   12d21s20
               end if                                                     12d21s20
              end do                                                     12d21s20
             end if                                                     3d2s21
             jvs=jvs+ncsf(jarg)*nvirt(jsbv)*nrootu                      12d24s20
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
      end do                                                            12d18s20
      ibcoff=ircv                                                       1d30s21
      return
      end
