c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hccsfdx(vec,ncsft,nff0,mff0,ncsf,nrootz,mdon,mdoop,    1d27s23
     $     ixw1,ixw2,nec,ism,irel,norb,nbasdws,idoubo,                  4d27s23
     $     nsymb,iden,multh,bc,ibc,nh0av,ld2e,i2e)                      4d25s23
      implicit real*8 (a-h,o-z)                                         7d11s19
c
c     one particle density ... sum over roots
c
      external second                                                   8d1s19
      integer*8 ipack,itesta,itestb,i0kc,i0ko,j0bc,j0bo                     5d7s21
      integer*2 ipack2(4)                                               12d1s19
      integer*1 nab1(2),nab2(2)                                         5d7s21
      equivalence (ipack,ipack2)                                        12d1s19
      logical ldebug,lpr,ld2e,lqy                                           4d25s23
      dimension vec(ncsft,nrootz),mff0(*),nff0(mdoop,3),ncsf(*),
     $     nother(2),ism(*),irel(*),nab4(2,3),multh(8,8),               4d27s23
     $     iden(*),nbasdws(*),idoubo(*),itest(64,2),nh0av(*),i2e(*)     4d25s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      data icall/0/
      data loopx/1000/
      include "common.store"                                            7d11s19
      include "common.print"                                            1d3s20
      save
      write(6,*)('hi, my name is hccsfdx ...'),icall
      write(6,*)('irel,ism: ')
      do i=1,norb
       write(6,*)i,irel(i),ism(i)
      end do
      igoal=i2e(1)+56
      write(6,*)('goal at start: '),bc(igoal),igoal
      loop=0
      ibcoffo=ibcoff                                                    11d2s22
      ldebug=.false.                                                    5d12s21
       icall=icall+1
      lpr=icall.eq.1
      if(ldebug.or.lpr)then                                                    1d3s20
       write(6,*)('Hi, my name is hccsfdx '),ncsft,nrootz,
     $     mdon,nec,mdoop
       write(6,*)('call '),icall
       write(6,*)('ism: '),(ism(i),i=1,norb)
       write(6,*)('rel: '),(irel(i),i=1,norb)
       do nclop=mdon+1,mdoop
        iarg=nclop-mdon
        iff=nff0(nclop,2)
        ipk=nff0(nclop,3)
        nclo=nclop-1                                                    7d22s21
        nopen=nec-nclo*2                                                7d27s22
        write(6,*)('for nclo '),nclo,nff0(nclop,1),
     $       iff,ipk
        do if=1,nff0(nclop,1)
         iff0x=iff
         i0kc=mff0(iff)
         if(popcnt(mff0(iff)).ne.nclo)then
          write(6,*)('we want nclo = '),nclo,popcnt(mff0(iff))
          write(6,*)('we''ve got ')
          call dcbit(i0kc,norb,'kc')
          write(6,*)('this is for if = '),if,iff,
     $         i0kc,norb
          call dws_synca
          call dws_finalize
          stop
         end if
         iff=iff+1
         i0ko=mff0(iff)
         if(popcnt(mff0(iff)).ne.nopen)then
          write(6,*)('we want nopen = '),nopen,popcnt(i0ko),
     $         popcnt(mff0(iff))
          write(6,*)('we''ve got ')
          call dcbit(i0ko,norb,'ko')
          write(6,*)('this is for if = '),if
          call dws_synca
          call dws_finalize
          stop
         end if
         iff=iff+1
         sz=0d0                                                         7d16s21
         do ir=1,nrootz                                                 7d16s21
          do i=0,ncsf(iarg)-1                                           7d16s21
           sz=sz+vec(ipk+i,ir)**2                                       7d16s21
          end do                                                        7d16s21
         end do                                                         7d16s21
         sum=sz
         sz=sqrt(sz/dfloat(nrootz*ncsf(iarg)))                          7d16s21
         if(sz.gt.1d-10)then
          write(6,*)('ket '),nclop,if,ipk,sum
          call dcbit(i0kc,norb,'c')
          call dcbit(i0ko,norb,'o')
 1493     continue
          call prntm2(vec(ipk,1),ncsf(iarg),nrootz,ncsft)
         end if                                                         7d16s21
         ipk=ipk+ncsf(iarg)
        end do
       end do
      end if                                                            1d3s20
      igoalr=1
      mdoo=mdoop-1                                                      5d10s21
      loopit=0                                                          1d27s23
      do nclokp=mdon+1,mdoop                                             5d7s21
       nclok=nclokp-1                                                     5d7s21
       iargk=nclokp-mdon                                                  7d11s19
       nopenk=nec-2*nclok                                                 7d11s19
       iff0=nff0(nclokp,2)                                               5d7s21
       ipk=nff0(nclokp,3)                                                 5d7s21
       do if=1,nff0(nclokp,1)                                            5d7s21
        i0kc=mff0(iff0)                                                  5d7s21
        if(popcnt(i0kc).ne.popcnt(mff0(iff0)))then
         write(6,*)('popcnt error on i0kc: ')
         call dcbit(mff0(iff0),32,'mff0')
         call dcbit(i0kc,64,'i0kc')
           stop 'hccsfdx'
        end if
        iff0=iff0+1                                                     5d7s21
        i0ko=mff0(iff0)                                                  5d7s21
        if(popcnt(i0ko).ne.popcnt(mff0(iff0)))then                      1d30s23
         write(6,*)('popcnt error on i0ko: ')
         call dcbit(mff0(iff0),32,'mff0')
         call dcbit(i0ko,64,'i0ko')
           stop 'hccsfdx'
        end if                                                          1d30s23
        iff0=iff0+1                                                     5d7s21
        do nclobp=max(mdon+1,nclokp-2),min(mdoop,nclokp+2)              5d18s23
         nclob=nclobp-1                                                  6d11s19
         nopenb=nec-2*nclob
         jargb=nclobp-mdon                                               6d11s19
         jff0=nff0(nclobp,2)                                               5d7s21
         jpb=nff0(nclobp,3)                                                 5d7s21
         do jf=1,nff0(nclobp,1)                                         5d7s21
          j0bc=mff0(jff0)                                                  5d7s21
          if(popcnt(j0bc).ne.popcnt(mff0(jff0)))then                    1d30s23
           write(6,*)('popcnt error on j0bc')
           call dcbit(mff0(jff0),32,'mff0')
           call dcbit(j0bc,64,'j0bc')
           stop 'hccsfdx'
          end if                                                        1d30s23
          jff0=jff0+1                                                     5d7s21
          j0bo=mff0(jff0)                                                  5d7s21
          if(popcnt(j0bo).ne.popcnt(mff0(jff0)))then                    1d30s23
           write(6,*)('popcnt error on j0bo')
           call dcbit(mff0(jff0),32,'mff0')
           call dcbit(j0bo,64,'j0bo')
           stop 'hccsfdx'
          end if                                                        1d30s23
          jff0=jff0+1                                                     5d7s21
          call gandc4(j0bc,j0bo,i0kc,i0ko,nopenb,nopenk,norb,nnot,nab4, 11d14s22
     $         bc,ibc)                                                  11d14s22
          if(nnot.gt.0)then                                               11d4s19
           if(mod(loopit,mynprocg).eq.mynowprog)then                      12d9s19
            loopit=mynowprog                                            9d10s24
            ibcsav=ibcoff                                                  12d9s19
            nn=ncsf(jargb)*ncsf(iargk)                                     12d9s19
           if(nnot.eq.1)then                                            4d13s21
             sum=0d0                                                    5d20s21
             do ir=1,nrootz                                             11d2s22
              do j=0,ncsf(jargb)-1                                      11d2s22
               sum=sum+vec(ipk+j,ir)**2                                 11d2s22
               if(abs(sum).gt.1d-10)
     $              write(6,*)('sum '),vec(ipk+j,ir),j,ir,sum
              end do                                                    11d2s22
             end do                                                     11d2s22
             do i=1,norb                                                5d21s21
              is=ism(i)                                                 5d20s21
              ig=irel(i)-1                                              4d25s23
              jmat=iden(is)+ig*(nh0av(is)+1)                            4d25s23
              if(btest(i0kc,i))then                                     5d20s21
                orig=bc(igoal)
               bc(jmat)=bc(jmat)+2d0*sum                                11d2s22
                if(abs(orig-bc(igoal)).gt.1d-14)write(6,*)('goal ch'),
     $               orig,sum,bc(igoal),nclokp,if,i
              end if                                                    5d20s21
              if(btest(i0ko,i))then                                     5d20s21
                orig=bc(igoal)
               bc(jmat)=bc(jmat)+sum                                    11d2s22
                if(abs(orig-bc(igoal)).gt.1d-14)write(6,*)('goal oh'),
     $               orig,sum,bc(igoal),nclokp,if,i
              end if                                                    5d20s21
             end do                                                     5d20s21
             if(ld2e)then                                               4d25s23
              if(abs(sum).gt.1d-10)then
               call dcbit(i0kc,norb,'i0kc')
               call dcbit(i0ko,norb,'i0ko')
              end if
              do i1=1,norb                                              5d21s21
               if(btest(i0kc,i1))then                                   5d20s21
                jsa=ism(i1)                                             5d20s21
                jga=irel(i1)-1                                          5d18s23
                do i2=1,norb                                            5d21s21
                 if(btest(i0kc,i2))then                                 5d20s21
                  jsb=ism(i2)                                           5d20s21
                  jsab=multh(jsa,jsb)                                   5d20s21
                  jgb=irel(i2)-1                                        5d18s23
                  jxint=igetint(i2e,jsa,jsa,jsb,jsb,jga,jga,jgb,jgb,nn, 4d25s23
     $                 .false.,bc,ibc)                                  4d25s23
                  orig=bc(igoal)
                  bc(jxint)=bc(jxint)+2d0*sum                           4d25s23
                  if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)('goal 2a'),
     $                 orig,sum,bc(igoal)
                  kxint=igetint(i2e,jsa,jsb,jsb,jsa,jga,jgb,jgb,jga,nn, 4d25s23
     $                 .false.,bc,ibc)                                  4d25s23
                  orig=bc(igoal)
                  bc(kxint)=bc(kxint)-sum                               4d25s23
                  if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)('goal 2b'),
     $                 orig,sum,bc(igoal)
                 end if                                                 5d20s21
                 if(btest(i0ko,i2))then                                 5d20s21
                  jsb=ism(i2)                                           5d20s21
                  jsab=multh(jsa,jsb)                                   5d20s21
                  jgb=irel(i2)-1                                        4d25s23
                  lqy=abs(sum).gt.1d-10
                  jxint=igetint(i2e,jsa,jsa,jsb,jsb,jga,jga,jgb,jgb,nn, 4d25s23
     $                 .false.,bc,ibc)                                  4d25s23
                  orig=bc(igoal)
                  bc(jxint)=bc(jxint)+2d0*sum                           4d25s23
                  if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)('goal 2c'),
     $                 orig,sum,bc(igoal),i2,i1
                  kxint=igetint(i2e,jsa,jsb,jsb,jsa,jga,jgb,jgb,jga,nn, 4d25s23
     $                 lqy,bc,ibc)                                      5d22s23
                  orig=bc(igoal)
                  bc(kxint)=bc(kxint)-sum                               4d25s23
                  if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)('goal 2d'),
     $                 orig,sum,bc(igoal),jga,jsa,jgb,jsb
                 end if                                                 5d20s21
                end do                                                  5d20s21
               end if
              end do                                                    5d20s21
              do i1=1,norb                                              5d21s21
               if(btest(i0ko,i1))then                                   5d20s21
                jsa=ism(i1)                                             5d20s21
                jga=irel(i1)-1                                          4d25s23
                do i2=1,norb                                            5d21s21
                 if(btest(i0ko,i2).and.i1.ne.i2)then                    5d20s21
                  jsb=ism(i2)                                           5d20s21
                  jsab=multh(jsa,jsb)                                   5d20s21
                  jgb=irel(i2)-1                                        4d25s23
                  ixint=igetint(i2e,jsa,jsa,jsb,jsb,jga,jga,jgb,jgb,nn, 4d25s23
     $                 .false.,bc,ibc)                                   4d25s23
                  orig=bc(igoal)
                  bc(ixint)=bc(ixint)+0.5d0*sum                         4d25s23
                  if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)('goal 2e'),
     $                 orig,sum,bc(igoal),i2,i1,jga,jsa,jgb,jsb
                 end if                                                 5d20s21
                end do                                                  5d20s21
               end if
              end do                                                    5d20s21
              imat=ibcoff                                               4d27s23
              ibcoff=imat+ncsf(iargk)*ncsf(iargk)                       4d27s23
              call enough('hccsfdx.mat1',bc,ibc)                        4d27s23
              do i1=1,norb-1                                             5d7s21
               if(btest(i0ko,i1))then                                     5d7s21
                do i2=i1+1,norb                                          5d7s21
                 if(btest(i0ko,i2))then                                   5d7s21
                  itesta=i0kc                                              5d7s21
                  itestb=i0ko                                              5d7s21
                  nopenkk=nopenk-2                                              11d13s20
                  karg=iargk+1                                                 11d13s20
                  nab1(1)=i2                                             5d7s21
                  nab1(2)=i1                                             5d7s21
                  nab2(1)=nab1(2)                                           12d6s20
                  nab2(2)=nab1(1)                                           12d6s20
                  jsa=ism(nab1(1))                                              11d13s20
                  jga=irel(nab1(1))-1                                   4d27s23
                  jsb=ism(nab1(2))                                      5d20s21
                  jgb=irel(nab1(2))-1                                   4d27s23
                  jsc=ism(nab2(1))                                      5d20s21
                  jgc=irel(nab2(1))-1                                   4d27s23
                  jsd=ism(nab2(2))                                      5d20s21
                  jgd=irel(nab2(2))-1                                   4d27s23
                  nqq=karg+mdon-1
                  if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                   itestb=ibclr(itestb,i2)                                5d7s21
                   itestb=ibclr(itestb,i1)                               5d7s21
                   itesta=ibset(itesta,i1)                               5d7s21
                   call enough('hccsfdx.mat1',bc,ibc)                   4d27s23
                   call gandc(j0bc,j0bo,itesta,itestb,nopenb,              5d7s21
     $               nopenkk,jargb,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,   12d6s20
     $               iwpb1,iwpk1,ncsfmid1,bc,ibc)                       11d14s22
                   call gandc(itesta,itestb,i0kc,i0ko,nopenkk,nopenk,        5d7s21
     $         karg,iargk,ncsf,norb,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,    11d13s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
                   call genmat(ncsf(jargb),ncsf(karg),ncsf(iargk),
     $                  ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod,11d10s22
     $                  bc,ibc)                                         11d10s22
                   do i=0,ncsf(iargk)*ncsf(iargk)-1                           4d13s21
                    bc(imat+i)=bc(iprod+i)                                  4d13s21
                   end do                                                   4d13s21
                   ibcoff=iprod                                                11d13s20
                  else                                                        11d13s20
                   do i=0,ncsf(iargk)*ncsf(iargk)-1                           4d13s21
                    bc(imat+i)=0d0                                          4d13s21
                   end do                                                   4d13s21
                  end if                                                      11d13s20
                  do i=0,ncsf(iargk)-1                                       4d13s21
                   ii=imat+i*(ncsf(iargk)+1)                                 4d13s21
                   bc(ii)=bc(ii)-1d0                                        4d13s21
                  end do                                                    4d13s21
                  jsab=multh(jsa,jsb)                                   5d20s21
                  jscd=multh(jsc,jsd)                                   5d20s21
                  ixint=igetint(i2e,jsa,jsb,jsc,jsd,jga,jgb,jgc,jgd,    4d27s23
     $                 nn,.false.,bc,ibc)                               4d27s23
                  itmp=ibcoff                                           4d27s23
                  ibcoff=itmp+ncsf(jargb)*nrootz                        4d27s23
                  call enough('hccsfdx.tmp1',bc,ibc)                    4d27s23
                  call dgemm('n','n',ncsf(jargb),nrootz,ncsf(iargk),    4d27s23
     $                   1d0,bc(imat),ncsf(jargb),vec(ipk,1),ncsft,     5d18s23
     $                   0d0,bc(itmp),ncsf(jargb),                      4d27s23
     d' hccsfdx.  1')
                  sum=0d0                                               4d27s23
                  do i=1,nrootz                                         4d27s23
                   iad1=itmp+ncsf(jargb)*(i-1)                          4d27s23
                   do j=0,ncsf(jargb)-1                                 4d27s23
                    sum=sum+bc(iad1+j)*vec(ipk+j,i)                      4d27s23
                    if(ixint.eq.igoal)write(6,*)bc(iad1+j),vec(ipk+j,i),
     $                   sum,j
                   end do                                               4d27s23
                  orig=bc(igoal)
                   bc(ixint)=bc(ixint)+sum                              4d27s23
                   if(abs(orig-bc(igoal)).gt.1d-10)then
                    write(6,*)('goal 2f'),
     $                 orig,sum,bc(igoal)
                    write(6,*)('imat ')
                    call prntm2(bc(imat),ncsf(iargk),ncsf(iargk),
     $                   ncsf(iargk))
                    call dcbit(j0bc,norb,'j0bc')
                    call dcbit(j0bo,norb,'j0bo')
                    call dcbit(itesta,norb,'itesta')
                    call dcbit(itestb,norb,'itestb')
                    call dcbit(i0kc,norb,'i0kc')
                    call dcbit(i0ko,norb,'i0ko')
                   end if
                  end do                                                4d27s23
                  ibcoff=itmp                                           4d27s23
                 end if                                                  5d7s21
                end do                                                       11d13s20
               end if                                                    5d7s21
              end do                                                        11d13s20
              ibcoff=imat                                                 4d13s21
             end if                                                     5d20s21
            else if(nnot.eq.2)then                                       4d13s21
             call gandc(j0bc,j0bo,i0kc,i0ko,nopenb,                         5d7s21
     $          nopenk,jargb,iargk,ncsf,norb,ixw1,ixw2,nnot1,nab1,iwpb1, 11d16s20
     $            iwpk1,ncsfmid1,bc,ibc)                                11d14s22
             idvtmp=ibcoff                                              5d20s21
             ibcoff=idvtmp+ncsf(jargb)*nrootz                           5d20s21
             call enough('hccsfbk.  3',bc,ibc)                          1d27s23
             call xtimesn(ncsf(jargb),ncsf(iargk),ncsfmid1,nrootz,iwpb1,5d20s21
     $            iwpk1,vec(ipk,1),ncsft,bc(idvtmp),ncsf(jargb),1d0,    11d2s22
     $            0d0,bc,ibc)                                           11d15s22
             sum=0d0                                                    4d27s23
             do ir=1,nrootz                                             4d27s23
              iad1=idvtmp+ncsf(jargb)*(ir-1)                            4d27s23
              do i=0,ncsf(jargb)-1                                      4d27s23
               sum=sum+bc(iad1+i)*vec(jpb+i,ir)                         5d18s23
              end do                                                    4d27s23
             end do                                                     4d27s23
             jsb=ism(nab1(1))                                           5d12s21
             jsk=ism(nab1(2))                                           5d12s21
             jsbk=multh(jsb,jsk)
             if(jsbk.ne.1)then                                          1d27s23
              write(6,*)('wait a sec ... jsbk = '),jsbk,
     $            (' ne isymop = '),1
              write(6,*)('jsb,jsk '),jsb,jsk,nab1(1),nab1(2)            8d16s22
              stop                                                      5d13s21
             end if                                                     5d13s21
             jgb=irel(nab1(1))-1                                        4d25s23
             jgk=irel(nab1(2))-1                                        4d25s23
             if(ld2e)then                                               4d27s23
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
              jsa=ism(nab4(2,1))                                             11d13s20
              jsb=ism(nab4(1,1))                                             11d13s20
              jga=irel(nab4(2,1))-1                                     4d27s23
              jgb=irel(nab4(1,1))-1                                     4d27s23
              do i=1,nok                                                    11d13s20
               js=ism(itest(i,2))                                           11d13s20
               jg=irel(itest(i,2))-1                                    4d27s23
               ixint=igetint(i2e,jsa,jsb,js,js,jga,jgb,jg,jg,nn,.false., 4d27s23
     $              bc,ibc)                                             4d27s23
               if(itest(i,1).eq.2)then                                  4d27s23
                  orig=bc(igoal)
                bc(ixint)=bc(ixint)+2d0*sum                             4d27s23
                if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)('goal 2g'),
     $               orig,sum,bc(igoal)
                ixint=igetint(i2e,jsa,js,js,jsb,jga,jg,jg,jgb,nn,       4d27s23
     $               .false.,bc,ibc)                                    4d27s23
                  orig=bc(igoal)
                bc(ixint)=bc(ixint)-sum                                 4d27s23
                if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)('goal 2h'),
     $               orig,sum,bc(igoal)
               else                                                     4d27s23
                  orig=bc(igoal)
                bc(ixint)=bc(ixint)+sum                                 4d27s23
                if(abs(orig-bc(igoal)).gt.1d-10)write(6,*)('goal 2i'),
     $               orig,sum,bc(igoal)
               end if                                                   4d27s23
              end do                                                    5d20s21
             end if                                                     5d21s21
             iad=iden(jsb)+jgb+nh0av(jsb)*jgk                           4d25s23
             orig=bc(igoal)
             bc(iad)=bc(iad)+sum                                        4d27s23
             if(abs(orig-bc(igoal)).gt.1d-10)then
              write(6,*)('veck: ')
              call prntm2(vec(ipk,1),ncsf(iargk),nrootz,ncsft)
              write(6,*)('dvtmp ')
              call prntm2(bc(idvtmp),ncsf(jargb),nrootz,ncsf(jargb))
              write(6,*)('vecb: ')
              call prntm2(vec(jpb,1),ncsf(jargb),nrootz,ncsf(jargb))
              write(6,*)('goalsum '),
     $            orig,sum,bc(igoal),nclobp,jf,nclokp,if
             end if
             ibcoff=idvtmp                                              5d20s21
             if(ld2e)then                                               4d27s23
              do i=1,nok                                                    11d13s20
               if(itest(i,1).eq.1)then                                      11d13s20
                itesta=j0bc                                               5d7s21
                itestb=j0bo                                               5d7s21
                nopenkk=nopenb                                                11d13s20
c
c     anihilate common
c
                if(btest(itesta,itest(i,2)))then                             11d13s20
                 itesta=ibclr(itesta,itest(i,2))                             11d13s20
                 itestb=ibset(itestb,itest(i,2))                             11d13s20
                 karg=jargb-1                                                11d13s20
                 nopenkk=nopenkk+1                                             11d13s20
                else                                                         11d13s20
                 itestb=ibclr(itestb,itest(i,2))                             11d13s20
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
                nqq=karg+mdon-1
                if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                 call gandc(j0bc,j0bo,itesta,itestb,nopenb,                5d7s21
     $              nopenkk,jargb,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,    11d16s20
     $              iwpb1,iwpk1,ncsfmid1,bc,ibc,bc,ibc)                 11d14s22
                 call gandc(itesta,itestb,i0kc,i0ko,nopenkk,nopenk,          5d7s21
     $         karg,iargk,ncsf,norb,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,    11d13s20
     $         ncsfmid2,bc,ibc,bc,ibc)                                  11d14s22
                 if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                  jsa=ism(nab1(1))                                      5d20s21
                  jga=irel(nab1(1))-1                                   4d27s23
                  jsb=ism(nab1(2))                                      5d20s21
                  jgb=irel(nab1(2))-1                                   4d27s23
                  jsc=ism(nab2(1))                                      5d20s21
                  jgc=irel(nab2(1))-1                                   4d27s23
                  jsd=ism(nab2(2))                                      5d20s21
                  jgd=irel(nab2(2))-1                                   4d27s23
                  jsab=multh(jsa,jsb)                                   5d20s21
                  sum=0d0                                               5d20s21
                  ixint=igetint(i2e,jsa,jsb,jsc,jsd,jga,jgb,jgc,jgd,nn, 4d27s23
     $                 .false.,bc,ibc)                                  4d27s23
                  itmpgb=ibcoff                                         5d20s21
                  ibcoff=itmpgb+ncsf(jargb)*nrootz                      5d20s21
                  call enough('hccsfbk.  4',bc,ibc)
                  do iqq=itmpgb,ibcoff-1                                5d20s21
                   bc(iqq)=0d0                                          5d20s21
                  end do                                                5d20s21
                  call genmatn(ncsf(jargb),ncsf(karg),ncsf(iargk),      5d20s21
     $                ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,        5d20s21
     $                 vec(ipk,1),ncsft,nrootz,bc(itmpgb),0d0,bc,ibc)   4d27s23
                  jtmpgb=itmpgb                                         5d20s21
                  sum=0d0                                               4d27s23
                  do ir=1,nrootz
                   do j=0,ncsf(jargb)-1
                    sum=sum+bc(jtmpgb+j)*vec(jpb+j,ir)                  4d27s23
                   end do                                               5d20s21
                   jtmpgb=jtmpgb+ncsf(jargb)                            5d20s21
                  end do                                                5d20s21
                  orig=bc(igoal)
                  bc(ixint)=bc(ixint)+sum                               4d27s23
                  if(abs(orig-bc(igoal)).gt.1d-10)then
                   write(6,*)('goal 2j'),
     $               orig,sum,bc(igoal)
                   call dcbit(j0bc,norb,'j0bc')
                   call dcbit(j0bo,norb,'j0bo')
                   call dcbit(i0kc,norb,'i0kc')
                   call dcbit(i0ko,norb,'i0ko')
                  end if
                  ibcoff=itmpgb                                         5d20s21
                 else if(max(nnot1,nnot2).lt.2)then                        11d13s20
                  write(6,*)('c:expecting 2s, but got otherwise')
                  call dcbit(itesta,norb,'itesta')
                  call dcbit(itestb,norb,'itestb')
                  stop
                 end if
                end if                                                    11d13s20
               end if                                                       11d13s20
              end do                                                        11d13s20
             end if                                                     5d20s21
            else                                                        4d27s23
             if(ld2e)then                                               4d27s23
              ipssx=0                                                   4d27s23
              if(nnot.eq.3)then                                         2d6s23
               ipssx=1                                                  2d6s23
              else if(nnot.eq.4)then                                    2d6s23
               ipssx=3                                                  2d6s23
              end if                                                    2d6s23
              do ipss=1,ipssx                                             1d22s21
               if(ipss.eq.1)then                                          1d22s21
                iu1=1                                                     1d22s21
                iu2=1                                                     1d22s21
               else if(ipss.eq.2)then                                     1d22s21
                iu1=1                                                     1d22s21
                iu2=2                                                     1d22s21
               else                                                       1d22s21
                iu1=2                                                     1d22s21
                iu2=1                                                     1d22s21
               end if                                                     12d8s20
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
                write(6,*)('iu1 = '),iu1
                write(6,*)('nab4 '),nab4
                call dcbit(j0bc,norb,'j0bc')                            7d27s22
                call dcbit(j0bo,norb,'j0bo')                            7d27s22
                call dcbit(i0kc,norb,'i0kc')                            7d27s22
                call dcbit(i0ko,norb,'i0ko')                            7d27s22
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
               nqq=karg+mdon-1                                            1d26s21
               if(nqq.ge.mdon.and.nqq.le.mdoo)then                        1d26s21
                call gandc(j0bc,j0bo,itesta,itestb,nopenb,                 5d10s21
     $              nopenkk,                                             5d7s21
     $         jargb,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,    11d13s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
                call gandc(itesta,itestb,i0kc,i0ko,nopenkk,nopenk,           5d10s21
     $         karg,iargk,ncsf,norb,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,    11d13s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
                if(nnot1.eq.2.and.nnot2.eq.2)then                         1d22s21
                 jsa=ism(nab1(1))                                       5d20s21
                 jga=irel(nab1(1))-1                                    4d27s23
                 jsb=ism(nab1(2))                                       5d20s21
                 jgb=irel(nab1(2))-1                                    4d27s23
                 jsc=ism(nab2(1))                                       5d20s21
                 jgc=irel(nab2(1))-1                                    4d27s23
                 jsd=ism(nab2(2))                                       5d20s21
                 jgd=irel(nab2(2))-1                                    4d27s23
                 jsab=multh(jsa,jsb)                                    5d20s21
                 ixint=igetint(i2e,jsa,jsb,jsc,jsd,jga,jgb,jgc,jgd,nn,  4d27s23
     $                .false.,bc,ibc)                                   4d27s23
                 itmpgb=ibcoff                                          5d20s21
                 ibcoff=itmpgb+ncsf(jargb)*nrootz                       5d20s21
                 call enough('hccsfbk.  5',bc,ibc)
                 do i=itmpgb,ibcoff-1                                   5d20s21
                  bc(i)=0d0                                             5d20s21
                 end do                                                 5d20s21
                 call genmatn(ncsf(jargb),ncsf(karg),ncsf(iargk),       5d20s21
     $                 ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,       5d20s21
     $                 vec(ipk,1),ncsft,nrootz,bc(itmpgb),0d0,bc,ibc)   4d27s23
                 jtmpgb=itmpgb
                 sum=0d0                                                4d27s23
                 do ir=1,nrootz
                  do j=0,ncsf(jargb)-1
                   orig=sum
                   sum=sum+bc(jtmpgb+j)*vec(jpb+j,ir)                   4d27s23
                   write(6,*)orig,bc(jtmpgb+j),vec(jpb+j,ir),sum,j
                  end do                                                5d20s21
                  jtmpgb=jtmpgb+ncsf(jargb)                             5d20s21
                 end do                                                 5d20s21
                 sum0=sum
                 if(nab1(1).eq.nab2(1).and.nab1(2).eq.nab2(2))then      8d5s21
                  sum=sum*0.5d0                                         8d5s21
                 end if
                  orig=bc(igoal)
                 bc(ixint)=bc(ixint)+sum                                4d27s23
                 if(abs(orig-bc(igoal)).gt.1d-10)then
                  write(6,*)('goal 2k'),
     $               orig,sum,sum0,bc(igoal),jga,jsa,jgb,jsb,jgc,jsc,
     $                 jgd,jsd
                  write(6,*)('vec(ipk')
                  call prntm2(vec(ipk,1),ncsf(iargk),nrootz,ncsft)
                 end if
                 ibcoff=itmpgb                                          5d20s21
                 if(ipss.eq.2)go to 22                                    1d22s21
                end if                                                    1d26s21
               end if                                                     1d22s21
              end do                                                      1d22s21
   22         continue                                                    1d22s21
             end if                                                     4d27s23
            end if                                                       4d13s21
            ibcoff=ibcsav                                                 12d9s19
           end if                                                         12d9s19
           loopit=loopit+1                                                12d9s19
          end if                                                          12d9s19
          jpb=jpb+ncsf(jargb)
         end do                                                         5d7s21
        end do                                                          5d7s21
        ipk=ipk+ncsf(iargk)
       end do                                                           6d11s19
      end do
      ibcoff=ibcoffo                                                    11d2s22
      write(6,*)('goal at end: '),bc(igoal)
      return
      end                                                               7d11s19
