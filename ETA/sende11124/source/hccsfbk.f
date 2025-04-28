c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hccsfbk(gb,ncsftb,veck,ncsftk,nff0b,nff0k,ncsf,mff0b,  5d12s21
     $     mff0k,nrootz,mdon,mdoop,ixw1,ixw2,nec,ism,irel,irefo,        5d12s21
     $     norb,lxmt,ixmt,multh,phase,nbasdws,idoubo,isymop,n2e,i2eop,  5d20s21
     $     phase2,nsymb,shift,ixmtf,bc,ibc)                             11d10s22
      implicit real*8 (a-h,o-z)                                         7d11s19
c
      external second                                                   8d1s19
      integer*8 ipack,itesta,itestb,i0kc,i0ko,j0bc,j0bo,gandcc,gandco,  2d6s23
     $     gandcb                                                       2d6s23
      integer*2 ipack2(4)                                               12d1s19
      integer*1 nab1(2),nab2(2)                                         5d7s21
      equivalence (ipack,ipack2)                                        12d1s19
      logical ldebug,lpr,lquery
      dimension gb(ncsftb,*),veck(ncsftk,nrootz),mff0k(*),mff0b(*),     5d12s21
     $     nff0b(mdoop,3),nff0k(mdoop,3),ncsf(*),nother(2),ism(*),             5d7s21
     $     irel(*),irefo(*),nab4(2,3),multh(8,8),ixmt(8,*),nbasdws(*),  5d20s21
     $     idoubo(*),i2eop(2,3),isymop(*),itest(64,2),ixmtf(8),ioxx(2)  2d6s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      common/fnd2cm/inv(2,8,8,8)                                        4d9s18
      data icall/0/
      include "common.store"                                            7d11s19
      include "common.print"                                            1d3s20
      save
      ldebug=.false.                                                    5d12s21
       icall=icall+1
      lpr=icall.eq.-3
      if(ldebug.or.lpr)then                                                    1d3s20
       write(6,*)('Hi, my name is hccsfbk '),ncsftb,ncsftk,nrootz,
     $     mdon,nec,phase,phase2,n2e,loc(phase)
       ldebug=lpr
       write(6,*)('call '),icall
       write(6,*)('ism: '),(ism(i),i=1,norb)
       write(6,*)('rel: '),(irel(i),i=1,norb)
       n0virt=9
       write(6,*)('n0virt= '),n0virt
       do nclop=mdon+1,mdoop
        iarg=nclop-mdon
        iff=nff0k(nclop,2)
        ipk=nff0k(nclop,3)
        nclo=nclop-1                                                    7d22s21
        nopen=nec-nclo*2                                                7d27s22
        write(6,*)('for nclo '),nclo,nff0k(nclop,1),nff0b(nclop,1),
     $       iff,ipk,loc(nff0k(nclop,3))
        do if=1,nff0k(nclop,1)
         iff0x=iff
         i0kc=mff0k(iff)
         if(popcnt(i0kc).ne.nclo)then
          write(6,*)('we want nclo = '),nclo
          write(6,*)('we''ve got ')
          call dcbit(i0kc,norb,'kc')
          write(6,*)('this is for if = '),if,iff,
     $         i0kc,norb,loc(mff0k(iff))
          call dws_synca
          call dws_finalize
          stop
         end if
         iff=iff+1
         i0ko=mff0k(iff)
         if(popcnt(i0ko).ne.nopen)then
          write(6,*)('we want nopen = '),nopen
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
           sz=sz+veck(ipk+i,ir)**2                                      7d16s21
          end do                                                        7d16s21
         end do                                                         7d16s21
         sz=sqrt(sz/dfloat(nrootz*ncsf(iarg)))                          7d16s21
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
          call prntm2(veck(ipk,1),ncsf(iarg),nrootz,ncsftk)
         end if                                                         7d16s21
         ipk=ipk+ncsf(iarg)
        end do
       end do
       do isb=1,nsymb                                                   5d20s21
        write(6,*)('for symmetry block '),isb,ixmt(isb,lxmt)
        isa=multh(isb,isymop(lxmt))                                     5d20s21
        call prntm2(bc(ixmt(isb,lxmt)),nbasdws(isb),nbasdws(isa),       5d20s21
     $       nbasdws(isb))
       end do
      end if                                                            1d3s20
      shift=0d0                                                         5d20s21
      do isb=1,nsymb                                                    5d20s21
       ixmtf(isb)=ibcoff                                                5d20s21
       isa=multh(isb,isymop(lxmt))                                      5d21s21
       ibcoff=ixmtf(isb)+nbasdws(isa)*nbasdws(isb)                      5d21s21
       call enough('hccsfbk.  1',bc,ibc)
       do i=0,nbasdws(isa)*nbasdws(isb)-1                               5d21s21
        bc(ixmtf(isb)+i)=bc(ixmt(isb,lxmt)+i)                           5d20s21
       end do                                                           5d20s21
      end do                                                            5d20s21
      ibctop=ibcoff                                                     7d9s21
      if(isymop(lxmt).eq.1)then                                         5d20s21
       do isb=1,nsymb                                                   5d20s21
        do i=0,idoubo(isb)-1                                            5d20s21
         jmat=ixmt(isb,lxmt)+i*(nbasdws(isb)+1)                         5d20s21
         shift=shift+2d0*bc(jmat)*phase                                 5d20s21
        end do
       end do                                                           5d20s21
      end if                                                            5d20s21
      do i2e=1,n2e                                                      5d20s21
       if(isymop(i2eop(1,i2e)).eq.1)then                                5d20s21
        do isb=1,nsymb                                                  5d20s21
         if(idoubo(isb).gt.0)then
          do isa=1,nsymb                                                5d20s21
           do ib=0,idoubo(isb)-1                                        5d20s21
            jmatb=ixmt(isb,i2eop(1,i2e))+ib*(nbasdws(isb)+1)            5d20s21
            do ia=0,idoubo(isa)-1                                       5d20s21
             jmata=ixmt(isa,i2eop(1,i2e))+ia*(nbasdws(isa)+1)           5d20s21
             shift=shift+2d0*bc(jmatb)*bc(jmata)*phase2                 5d20s21
            end do                                                      5d20s21
            do ia=idoubo(isa),nbasdws(isa)-1                            5d20s21
             do iap=idoubo(isa),nbasdws(isa)-1                          5d20s21
              jmata=ixmt(isa,i2eop(1,i2e))+iap+nbasdws(isa)*ia          5d20s21
              jxmtf=ixmtf(isb)+iap+nbasdws(isa)*ia                      5d20s21
              bc(jxmtf)=bc(jxmtf)+2d0*bc(jmata)*bc(jmatb)*phase2        5d20s21
             end do                                                     5d20s21
            end do                                                      5d20s21
           end do                                                       5d20s21
          end do                                                        5d20s21
         end if                                                         5d20s21
        end do                                                          5d20s21
       end if                                                           5d20s21
       do isb=1,nsymb                                                   5d20s21
        isa=multh(isb,isymop(i2eop(1,i2e)))                             5d20s21
        do ib=0,idoubo(isb)-1                                           5d20s21
         do ia=0,idoubo(isa)-1                                          5d20s21
          jmat1=ixmt(isa,i2eop(1,i2e))+ia+nbasdws(isa)*ib               5d21s21
          jmat2=ixmt(isb,i2eop(1,i2e))+ib+nbasdws(isb)*ia               5d21s21
          shift=shift-bc(jmat1)*bc(jmat2)*phase2                        5d21s21
         end do                                                         5d20s21
         do ia=idoubo(isa),nbasdws(isa)-1                               5d20s21
          jmata=ixmt(isb,i2eop(1,i2e))+ib+nbasdws(isb)*ia               5d21s21
          do iap=idoubo(isa),nbasdws(isa)-1                             5d20s21
           jmatb=ixmt(isa,i2eop(1,i2e))+iap+nbasdws(isa)*ib             5d20s21
           jxmtf=ixmtf(isa)+iap+nbasdws(isa)*ia                         5d20s21
           bc(jxmtf)=bc(jxmtf)-bc(jmata)*bc(jmatb)*phase2               5d20s21
          end do                                                        5d20s21
         end do                                                         5d20s21
        end do                                                          5d20s21
       end do                                                           5d20s21
      end do                                                            5d20s21
      if(ldebug)then                                                    5d21s21
       write(6,*)('shift = '),shift                                     5d21s21
       write(6,*)('1-e matrix with doubles part folded in')             5d21s21
       do isb=1,nsymb                                                   5d21s21
        isa=multh(isb,isymop(lxmt))                                      5d21s21
        write(6,*)('for symmetry block '),isb,isa,ixmtf(isb)
        call prntm2(bc(ixmtf(isb)),nbasdws(isb),nbasdws(isa),
     $       nbasdws(isb))
       end do                                                           5d21s21
      end if                                                            5d21s21
      igoal=112
      igoalr=1
      mdoo=mdoop-1                                                      5d10s21
      do nclokp=mdon+1,mdoop                                             5d7s21
       nclok=nclokp-1                                                     5d7s21
       iargk=nclokp-mdon                                                  7d11s19
       nopenk=nec-2*nclok                                                 7d11s19
       iff0=nff0k(nclokp,2)                                               5d7s21
       ipk=nff0k(nclokp,3)                                                 5d7s21
       do if=1,nff0k(nclokp,1)                                            5d7s21
        lquery=nclokp.eq.-4.and.if.eq.18449
        i0kc=mff0k(iff0)                                                  5d7s21
        iff0=iff0+1                                                     5d7s21
        i0ko=mff0k(iff0)                                                  5d7s21
        iff0=iff0+1                                                     5d7s21
        do nclobp=max(mdon+1,nclokp-2),min(mdoop,nclokp+2)                5d13s21
         nclob=nclobp-1                                                  6d11s19
         nopenb=nec-2*nclob
         jargb=nclobp-mdon                                               6d11s19
         jff0=nff0b(nclobp,2)                                               5d7s21
         jpb=nff0b(nclobp,3)                                                 5d7s21
         do jf=1,nff0b(nclobp,1)                                         5d7s21
          j0bc=mff0b(jff0)                                                  5d7s21
          jff0=jff0+1                                                     5d7s21
          j0bo=mff0b(jff0)                                                  5d7s21
          jff0=jff0+1                                                     5d7s21
          gandcc=ieor(j0bc,i0kc)                                        2d6s23
          gandco=ieor(j0bo,i0ko)                                        2d6s23
          gandcb=ior(gandcc,gandco)                                     2d6s23
          ndifb=popcnt(gandcb)                                          2d6s23
          nnot=0
          if(ndifb.le.4)then                                            2d6s23
           if(mod(loopit,mynprocg).eq.mynowprog)then                      12d9s19
            loopit=mynowprog                                            9d10s24
            ndifs=popcnt(gandco)                                        2d6s23
            ndifd=popcnt(gandcc)                                        2d6s23
            ibcsav=ibcoff                                                  12d9s19
            nn=ncsf(jargb)*ncsf(iargk)                                     12d9s19
            if(ndifs.eq.0.and.ndifd.eq.0)then                           2d6s23
             if(lquery)write(6,*)('in diag block for query ')
             sum=0d0                                                    5d20s21
             if(isymop(lxmt).eq.1)then                                  5d20s21
              do i=1,norb                                               5d21s21
               is=ism(i)                                                5d20s21
               ig=irel(i)-1+idoubo(is)                                  5d20s21
               jmat=ixmtf(is)+ig*(nbasdws(is)+1)                        5d20s21
               if(btest(i0kc,i))then                                    5d20s21
                orig=sum
                sum=sum+2d0*bc(jmat)*phase                              5d20s21
               end if                                                   5d20s21
               if(btest(i0ko,i))then                                    5d20s21
                orig=sum
                sum=sum+bc(jmat)*phase                                  5d20s21
               end if                                                   5d20s21
              end do                                                    5d20s21
             end if                                                     5d20s21
             sum0=sum
             if(lquery)write(6,*)('sum0 '),sum0
             if(n2e.gt.0)then                                           5d20s21
              do i1=1,norb                                              5d21s21
               if(btest(i0kc,i1))then                                   5d20s21
                jsa=ism(i1)                                             5d20s21
                jga=irel(i1)-1+idoubo(jsa)                              5d20s21
                do i2=1,norb                                            5d21s21
                 if(btest(i0kc,i2))then                                 5d20s21
                  jsb=ism(i2)                                           5d20s21
                  jsab=multh(jsa,jsb)                                   5d20s21
                  jgb=irel(i2)-1+idoubo(jsb)                            5d20s21
                  do i2e=1,n2e                                          5d20s21
                   if(isymop(i2eop(1,i2e)).eq.1)then                    5d20s21
                    jmat1=ixmt(jsa,i2eop(1,i2e))+jga+nbasdws(jsa)*jga   5d20s21
                    jmat2=ixmt(jsb,i2eop(1,i2e))+jgb+nbasdws(jsb)*jgb   5d20s21
                    fact=bc(jmat1)*bc(jmat2)*phase2                     5d20s21
                    orig=sum
                    sum=sum+2d0*fact                                    5d20s21
                   end if                                               5d20s21
                   if(isymop(i2eop(1,i2e)).eq.jsab)then                 5d20s21
                    jmat1=ixmt(jsa,i2eop(1,i2e))+jga+nbasdws(jsa)*jgb   5d20s21
                    jmat2=ixmt(jsb,i2eop(1,i2e))+jgb+nbasdws(jsb)*jga   5d20s21
                    fact=bc(jmat1)*bc(jmat2)*phase2                     5d21s21
                    orig=sum
                    sum=sum-fact                                        5d20s21
                   end if                                               5d20s21
                  end do                                                5d20s21
                 end if                                                 5d20s21
                 if(btest(i0ko,i2))then                                 5d20s21
                  jsb=ism(i2)                                           5d20s21
                  jsab=multh(jsa,jsb)                                   5d20s21
                  jgb=irel(i2)-1+idoubo(jsb)                            5d20s21
                  do i2e=1,n2e                                          5d20s21
                   if(isymop(i2eop(1,i2e)).eq.1)then                    5d20s21
                    jmat1=ixmt(jsa,i2eop(1,i2e))+jga+nbasdws(jsa)*jga   5d20s21
                    jmat2=ixmt(jsb,i2eop(1,i2e))+jgb+nbasdws(jsb)*jgb   5d20s21
                    fact=bc(jmat1)*bc(jmat2)*phase2                     5d20s21
                    orig=sum
                    sum=sum+2d0*fact                                    5d20s21
                   end if                                               5d20s21
                   if(isymop(i2eop(1,i2e)).eq.jsab)then                 5d20s21
                    jmat1=ixmt(jsa,i2eop(1,i2e))+jga+nbasdws(jsa)*jgb   5d21s21
                    jmat2=ixmt(jsb,i2eop(1,i2e))+jgb+nbasdws(jsb)*jga   5d21s21
                    fact=bc(jmat1)*bc(jmat2)*phase2                     5d21s21
                    orig=sum
                    sum=sum-fact                                        5d20s21
                   end if                                               5d20s21
                  end do                                                5d20s21
                 end if                                                 5d20s21
                end do                                                  5d20s21
               end if
              end do                                                    5d20s21
              do i1=1,norb                                              5d21s21
               if(btest(i0ko,i1))then                                   5d20s21
                jsa=ism(i1)                                             5d20s21
                jga=irel(i1)-1+idoubo(jsa)                              5d20s21
                do i2=1,norb                                            5d21s21
                 if(btest(i0ko,i2).and.i1.ne.i2)then                    5d20s21
                  jsb=ism(i2)                                           5d20s21
                  jsab=multh(jsa,jsb)                                   5d20s21
                  jgb=irel(i2)-1+idoubo(jsb)                            5d20s21
                  do i2e=1,n2e                                          5d20s21
                   if(isymop(i2eop(1,i2e)).eq.1)then                    5d20s21
                    jmat1=ixmt(jsa,i2eop(1,i2e))+jga+nbasdws(jsa)*jga   5d20s21
                    jmat2=ixmt(jsb,i2eop(1,i2e))+jgb+nbasdws(jsb)*jgb   5d21s21
                    fact=bc(jmat1)*bc(jmat2)*phase2                     5d20s21
                    orig=sum
                    sum=sum+0.5d0*fact                                  5d20s21
                   end if                                               5d20s21
                  end do                                                5d20s21
                 end if                                                 5d20s21
                end do                                                  5d20s21
               end if
              end do                                                    5d20s21
             end if                                                     5d20s21
             if(lquery)write(6,*)('sum after n2e: '),sum
             do ir=1,nrootz                                             5d20s21
              do j=0,ncsf(jargb)-1                                      5d20s21
               gb(jpb+j,ir)=gb(jpb+j,ir)+sum*veck(ipk+j,ir)             5d20s21
              end do                                                    5d20s21
             end do                                                     5d20s21
             if(lquery)then
              write(6,*)('gb so far: '),jpb
              call prntm2(gb(jpb,1),ncsf(jargb),nrootz,ncsftb)
             end if
             imat=ibcoff                                                5d12s21
             ibcoff=imat+nn                                             5d12s21
             call enough('hccsfbk.  2',bc,ibc)
             if(n2e.gt.0)then                                           5d20s21
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
                  jga=irel(nab1(1))-1+idoubo(jsa)                       5d20s21
                  jsb=ism(nab1(2))                                      5d20s21
                  jgb=irel(nab1(2))-1+idoubo(jsb)                       5d20s21
                  jsc=ism(nab2(1))                                      5d20s21
                  jgc=irel(nab2(1))-1+idoubo(jsc)                       5d20s21
                  jsd=ism(nab2(2))                                      5d20s21
                  jgd=irel(nab2(2))-1+idoubo(jsd)                       5d20s21
                  nqq=karg+mdon-1
                  if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                   itestb=ibclr(itestb,i2)                                5d7s21
                   itestb=ibclr(itestb,i1)                               5d7s21
                   itesta=ibset(itesta,i1)                               5d7s21
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
                  if(lquery)then
                   write(6,*)('mat for '),i1,i2
                   call prntm2(bc(imat),ncsf(iargk),ncsf(iargk),
     $                  ncsf(iargk))
                  end if
                  do i2e=1,n2e                                          5d20s21
                   if(isymop(i2eop(1,i2e)).eq.jsab)then                 5d20s21
                    jmat1=ixmt(jsa,i2eop(1,i2e))+jga+nbasdws(jsa)*jgb   5d20s21
                    jmat2=ixmt(jsc,i2eop(2,i2e))+jgc+nbasdws(jsc)*jgd   5d21s21
                    fact=bc(jmat1)*bc(jmat2)*phase2                     5d20s21
                    if(lquery)write(6,*)('fact '),fact,jga,jsa,jgb,jsb,
     $                   jgc,jsc,jgd,jsd
                    call dgemm('n','n',ncsf(jargb),nrootz,ncsf(iargk),  5d20s21
     $                   fact,bc(imat),ncsf(jargb),veck(ipk,1),ncsftk,  5d21s21
     $                   1d0,gb(jpb,1),ncsftb,                          5d20s21
     d' hccsfbk.  1')
             if(lquery)then
              write(6,*)('gb so far: '),jpb
              call prntm2(gb(jpb,1),ncsf(jargb),nrootz,ncsftb)
             end if
                   end if                                               5d20s21
                  end do                                                5d20s21
                 end if                                                  5d7s21
                end do                                                       11d13s20
               end if                                                    5d7s21
              end do                                                        11d13s20
             end if                                                     5d20s21
             ibcoff=imat                                                 4d13s21
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
             call gandc(j0bc,j0bo,i0kc,i0ko,nopenb,                         5d7s21
     $          nopenk,jargb,iargk,ncsf,norb,ixw1,ixw2,nnot1,nab1,iwpb1, 11d16s20
     $            iwpk1,ncsfmid1,bc,ibc,bc,ibc)                         11d14s22
             idvtmp=ibcoff                                              5d20s21
             ibcoff=idvtmp+ncsf(jargb)*nrootz                           5d20s21
             call enough('hccsfbk.  3',bc,ibc)
             call xtimesn(ncsf(jargb),ncsf(iargk),ncsfmid1,nrootz,iwpb1,5d20s21
     $            iwpk1,veck(ipk,1),ncsftk,bc(idvtmp),ncsf(jargb),1d0,  5d20s21
     $            0d0,bc,ibc)                                           11d10s22
             jsb=ism(nab1(1))                                           5d12s21
             jsk=ism(nab1(2))                                           5d12s21
             jsbk=multh(jsb,jsk)
             if(jsbk.ne.isymop(lxmt))then                               5d20s21
              write(6,*)('wait a sec ... jsbk = '),jsbk,
     $            (' ne isymop = '),isymop(lxmt)                        5d20s21
              write(6,*)('jsb,jsk '),jsb,jsk,nab1(1),nab1(2)            8d16s22
              stop                                                      5d13s21
             end if                                                     5d13s21
             jgb=irel(nab1(1))-1+idoubo(jsb)                            5d13s21
             jgk=irel(nab1(2))-1+idoubo(jsk)                            5d13s21
             iad=ixmtf(jsb)+jgb+nbasdws(jsb)*jgk                        5d20s21
             fact=bc(iad)*phase                                         5d13s21
             sum=fact                                                   5d20s21
             sum0=sum                                                   5d21s21
             if(n2e.gt.0)then                                           5d20s21
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
              jga=irel(nab4(2,1))-1+idoubo(jsa)                         5d20s21
              jgb=irel(nab4(1,1))-1+idoubo(jsb)                         5d20s21
              do i=1,nok                                                    11d13s20
               js=ism(itest(i,2))                                           11d13s20
               jg=irel(itest(i,2))-1+idoubo(js)                         5d20s21
               do i2e=1,n2e                                             5d20s21
                if(isymop(i2eop(1,i2e)).eq.1)then                       5d20s21
                 jmat1=ixmt(jsa,i2eop(1,i2e))+jga+nbasdws(jsa)*jgb      5d20s21
                 jmat2=ixmt(js,i2eop(2,i2e))+jg+nbasdws(js)*jg          5d21s21
                 fact=bc(jmat1)*bc(jmat2)*phase2                        5d20s21
                 if(itest(i,1).eq.2)then                                      11d13s20
                  orig=sum
                  sum=sum+fact*2d0                                      5d20s21
                 else                                                   5d20s21
                  orig=sum
                  sum=sum+fact                                          5d20s21
                 end if                                                 5d20s21
                end if                                                  5d20s21
                if(isymop(i2eop(1,i2e)).eq.multh(jsa,js)                5d20s21
     $               .and.itest(i,1).eq.2)then                          5d20s21
                 jmat1=ixmt(jsa,i2eop(1,i2e))+jga+nbasdws(jsa)*jg       5d20s21
                 jmat2=ixmt(js,i2eop(2,i2e))+jg+nbasdws(js)*jgb         5d21s21
                 fact=bc(jmat1)*bc(jmat2)*phase2                        5d20s21
                 orig=sum
                 sum=sum-fact                                           5d21s21
                end if                                                  5d20s21
               end do                                                   5d20s21
              end do                                                    5d20s21
             end if                                                     5d21s21
             jdvtmp=idvtmp                                              5d20s21
             do ir=1,nrootz                                             5d20s21
              do j=0,ncsf(jargb)-1                                      5d20s21
               gb(jpb+j,ir)=gb(jpb+j,ir)+sum*bc(jdvtmp+j)               5d20s21
              end do                                                    5d20s21
              jdvtmp=jdvtmp+ncsf(jargb)                                 5d20s21
             end do                                                     5d20s21
             ibcoff=idvtmp                                              5d20s21
             if(n2e.gt.0)then                                           5d21s21
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
                  jga=irel(nab1(1))-1+idoubo(jsa)                       5d20s21
                  jsb=ism(nab1(2))                                      5d20s21
                  jgb=irel(nab1(2))-1+idoubo(jsb)                       5d20s21
                  jsc=ism(nab2(1))                                      5d20s21
                  jgc=irel(nab2(1))-1+idoubo(jsc)                       5d20s21
                  jsd=ism(nab2(2))                                      5d20s21
                  jgd=irel(nab2(2))-1+idoubo(jsd)                       5d20s21
                  jsab=multh(jsa,jsb)                                   5d20s21
                  sum=0d0                                               5d20s21
                  do i2e=1,n2e                                          5d20s21
                   if(jsab.eq.isymop(i2eop(1,i2e)))then                 5d20s21
                    jmat1=ixmt(jsa,i2eop(1,i2e))+jga+nbasdws(jsa)*jgb   5d20s21
                    jmat2=ixmt(jsc,i2eop(2,i2e))+jgc+nbasdws(jsc)*jgd   5d21s21
                    fact=bc(jmat1)*bc(jmat2)*phase2                     5d20s21
                    orig=sum
                    sum=sum+fact                                        5d20s21
                   end if                                               5d20s21
                  end do                                                5d20s21
                  itmpgb=ibcoff                                         5d20s21
                  ibcoff=itmpgb+ncsf(jargb)*nrootz                      5d20s21
                  call enough('hccsfbk.  4',bc,ibc)
                  do iqq=itmpgb,ibcoff-1                                5d20s21
                   bc(iqq)=0d0                                          5d20s21
                  end do                                                5d20s21
                  call genmatn(ncsf(jargb),ncsf(karg),ncsf(iargk),      5d20s21
     $                ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,        5d20s21
     $                 veck(ipk,1),ncsftk,nrootz,bc(itmpgb),0d0,bc,ibc) 11d10s22
                  jtmpgb=itmpgb                                         5d20s21
                  do ir=1,nrootz
                   do j=0,ncsf(jargb)-1
                    gb(jpb+j,ir)=gb(jpb+j,ir)+bc(jtmpgb+j)*sum          5d21s21
                   end do                                               5d20s21
                   jtmpgb=jtmpgb+ncsf(jargb)                            5d20s21
                  end do                                                5d20s21
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
            else                                                         4d13s21
c     j is bra and i is ket
             ipssx=0                                                    2d6s23
             if(n2e.gt.0)then                                           2d6s23
              nnot=0                                                     2d6s23
              if(ndifs.eq.4.and.ndifb.eq.4)then                         2d6s23
               nnot=4                                                   2d6s23
               ioxx(1)=1                                                2d6s23
               ioxx(2)=1                                                2d6s23
               do i=1,norb                                              2d6s23
                if(btest(gandcb,i))then                                 2d6s23
                 if((btest(i0kc,i).and.btest(j0bo,i)).or.               2d6s23
     $               (btest(i0ko,i).and..not.btest(j0bc,i)))then        2d6s23
                  nab4(2,ioxx(2))=i                                     2d6s23
                  ioxx(2)=ioxx(2)+1                                     2d6s23
                 else                                                   2d6s23
                  nab4(1,ioxx(1))=i                                     2d6s23
                  ioxx(1)=ioxx(1)+1                                     2d6s23
                 end if                                                 2d6s23
                end if                                                  2d6s23
               end do                                                   2d6s23
              else if(ndifb.eq.3)then                                   2d6s23
               nnot=3                                                   2d6s23
               ioxx(1)=1                                                2d6s23
               ioxx(2)=1                                                2d6s23
               iswap=0                                                  2d6s23
               do i=1,norb                                              2d6s23
                if(btest(gandcb,i))then                                 2d6s23
                 if(btest(gandcc,i).and.                                2d6s23
     $        ((btest(j0bc,i).and..not.btest(i0ko,i)).or.               2d6s23
     $         (btest(i0kc,i).and..not.btest(j0bo,i))))then             2d6s23
                  if(btest(i0kc,i))iswap=1                              2d6s23
                  nab4(1,1)=i                                           2d6s23
                  nab4(1,2)=i                                           2d6s23
                 else                                                   2d6s23
                  nab4(2,ioxx(2))=i                                     2d6s23
                  ioxx(2)=ioxx(2)+1                                     2d6s23
                 end if                                                 2d6s23
                end if                                                  2d6s23
               end do                                                   2d6s23
               if(iswap.ne.0)then                                       2d6s23
                icpy=nab4(1,1)                                          2d6s23
                nab4(1,1)=nab4(2,1)                                     2d6s23
                nab4(2,1)=icpy                                          2d6s23
                icpy=nab4(1,2)                                          2d6s23
                nab4(1,2)=nab4(2,2)                                     2d6s23
                nab4(2,2)=icpy                                          2d6s23
                nbt=0                                                   2d6s23
                if(btest(j0bc,nab4(1,2)).and..not.btest(j0bc,nab4(1,1)))2d6s23
     $       nbt=1                                                      2d6s23
               else                                                     2d6s23
                nbt=0                                                   2d6s23
                if(btest(i0kc,nab4(2,2)).and..not.btest(i0kc,nab4(2,1)))2d6s23
     $              nbt=1                                               2d6s23
               end if                                                   2d6s23
               if(nbt.ne.0)then                                         2d6s23
                nab4(1,1)=nab4(1,2)                                     2d6s23
                nab4(2,1)=nab4(2,2)                                     2d6s23
               end if                                                   2d6s23
              else if(ndifs.eq.0.and.ndifd.eq.2)then                    2d6s23
               nnot=3                                                   2d6s23
               do i=1,norb                                              2d6s23
                if(btest(gandcb,i))then                                 2d6s23
                 if(btest(j0bc,i))then                                  2d6s23
                  nab4(1,1)=i                                           2d6s23
                  nab4(1,2)=i                                           2d6s23
                 else                                                   2d6s23
                  nab4(2,1)=i                                           2d6s23
                  nab4(2,2)=i                                           2d6s23
                 end if                                                 2d6s23
                end if                                                  2d6s23
               end do                                                   2d6s23
              end if                                                    2d6s23
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
                 jga=irel(nab1(1))-1+idoubo(jsa)                        5d21s21
                 jsb=ism(nab1(2))                                       5d20s21
                 jgb=irel(nab1(2))-1+idoubo(jsb)                        5d21s21
                 jsc=ism(nab2(1))                                       5d20s21
                 jgc=irel(nab2(1))-1+idoubo(jsc)                        5d21s21
                 jsd=ism(nab2(2))                                       5d20s21
                 jgd=irel(nab2(2))-1+idoubo(jsd)                        5d21s21
                 jsab=multh(jsa,jsb)                                    5d20s21
                 sum=0d0                                                5d20s21
                 do i2e=1,n2e                                           5d20s21
                  if(jsab.eq.isymop(i2eop(1,i2e)))then                  5d20s21
                   jmat1=ixmt(jsa,i2eop(1,i2e))+jga+nbasdws(jsa)*jgb    5d20s21
                   jmat2=ixmt(jsc,i2eop(2,i2e))+jgc+nbasdws(jsc)*jgd    5d21s21
                   fact=bc(jmat1)*bc(jmat2)*phase2                      5d20s21
                   orig=sum
                   sum=sum+fact                                         5d20s21
                  end if                                                5d20s21
                 end do                                                 5d20s21
                 sum0=sum                                               8d5s21
                 if(nab1(1).eq.nab2(1).and.nab1(2).eq.nab2(2))then      8d5s21
                  sum=sum*0.5d0                                         8d5s21
                 end if
                 itmpgb=ibcoff                                          5d20s21
                 ibcoff=itmpgb+ncsf(jargb)*nrootz                       5d20s21
                 call enough('hccsfbk.  5',bc,ibc)
                 do i=itmpgb,ibcoff-1                                   5d20s21
                  bc(i)=0d0                                             5d20s21
                 end do                                                 5d20s21
                 call genmatn(ncsf(jargb),ncsf(karg),ncsf(iargk),       5d20s21
     $                 ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,       5d20s21
     $                 veck(ipk,1),ncsftk,nrootz,bc(itmpgb),0d0,bc,ibc) 11d10s22
                 jtmpgb=itmpgb
                 do ir=1,nrootz
                  do j=0,ncsf(jargb)-1
                   gb(jpb+j,ir)=gb(jpb+j,ir)+bc(jtmpgb+j)*sum           5d21s21
                  end do                                                5d20s21
                  jtmpgb=jtmpgb+ncsf(jargb)                             5d20s21
                 end do                                                 5d20s21
                 ibcoff=itmpgb                                          5d20s21
                 if(ipss.eq.2)go to 22                                    1d22s21
                end if                                                    1d26s21
               end if                                                     1d22s21
              end do                                                      1d22s21
   22         continue                                                    1d22s21
             end if                                                     5d20s21
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
      ibcoff=ibctop                                                     7d9s21
      if(lpr)then
       ipb=1                                                            7d16s21
       rsum=0d0
       do nclop=mdon+1,mdoop
        iarg=nclop-mdon
        iff=nff0b(nclop,2)
        ipb=nff0b(nclop,3)
        do if=1,nff0b(nclop,1)
         i0bc=mff0b(iff)
         iff=iff+1
         i0bo=mff0b(iff)
         iff=iff+1
         sz=0d0
         do ir=1,nrootz
          do i=0,ncsf(iarg)-1
           sz=sz+gb(ipb+i,ir)**2
          end do
         end do
         sz=sqrt(sz/dfloat(nrootz*ncsf(iarg)))
         if(sz.gt.1d-10)then
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
           write(6,*)('bra '),nclop,if,ipb,ncsf(iarg),ipb-igoal,igoalr
           call prntm2(gb(ipb,1),ncsf(iarg),nrootz,ncsftb)               7d16s21
           do i=0,ncsf(iarg)-1
            rsum=rsum+gb(ipb+i,1)**2
           end do
           write(6,*)('rsum so far: '),rsum
          end if
         end if
         ipb=ipb+ncsf(iarg)
        end do
       end do
      end if
      return
      end                                                               7d11s19
