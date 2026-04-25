c mpec2.1 version theta. For copyright and Disclaimers, see start.f
      subroutine hccsfbkden(gb,ncsftb,veck,ncsftk,nff0b,nff0k,ncsf,
     $     mff0b,mff0k,nrootz,mdon,mdoop,ixw1,ixw2,nec,ism,irel,irefo,  3d10s25
     $     norb,iden,multh,nbasdws,idoubo,nsymb,bc,ibc)                 3d10s25
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
     $     irel(*),irefo(*),nab4(2,3),multh(8,8),nbasdws(*),            3d10s25
     $     idoubo(*),i2eop(2,3),itest(64,2),ixmtf(8),ioxx(2),iden(*)    3d12s25
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      common/fnd2cm/inv(2,8,8,8)                                        4d9s18
      data icall/0/
      include "common.store"                                            7d11s19
      include "common.print"                                            1d3s20
      save
      igoal=iden(1)+8-1+28*(3-1)
      write(6,*)('hi, my name is hccsfbkden '),loc(bc),loc(ibc)
      ldebug=.false.                                                    5d12s21
       icall=icall+1
      lpr=icall.eq.-3
      if(ldebug.or.lpr)then                                                    1d3s20
       write(6,*)('Hi, my name is hccsfbkden '),ncsftb,ncsftk,nrootz,
     $     mdon,nec
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
      end if                                                            1d3s20
      shift=0d0                                                         5d20s21
      ibctop=ibcoff                                                     7d9s21
      if(mynowprog.eq.0)then                                            3d12s25
       do isb=1,nsymb                                                    3d12s25
        do ir=0,nrootz-1                                                 3d12s25
         do i=0,idoubo(isb)-1                                            3d12s25
          jmat=iden(isb)+i+nbasdws(isb)*(i+nbasdws(isb)*ir)              3d12s25
          bc(jmat)=2d0                                                   3d12s25
         end do                                                          3d12s25
        end do                                                           3d12s25
       end do                                                            3d12s25
      end if
      if(ldebug)then                                                    5d21s21
       write(6,*)('shift = '),shift                                     5d21s21
       write(6,*)('1-e matrix with doubles part folded in')             5d21s21
      end if                                                            5d21s21
      igoalr=1
      mdoo=mdoop-1                                                      5d10s21
      do nclokp=mdon+1,mdoop                                             5d7s21
       nclok=nclokp-1                                                     5d7s21
       iargk=nclokp-mdon                                                  7d11s19
       nopenk=nec-2*nclok                                                 7d11s19
       iff0=nff0k(nclokp,2)                                               5d7s21
       ipk=nff0k(nclokp,3)                                                 5d7s21
       do if=1,nff0k(nclokp,1)                                            5d7s21
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
             irsum=ibcoff                                               3d12s25
             ibcoff=irsum+nrootz                                        3d12s25
             call enough('hccsfbkden.irsum',bc,ibc)                     3d12s25
             do iz=irsum,ibcoff-1                                       3d12s25
              bc(iz)=0d0                                                3d12s25
             end do                                                     3d12s25
             do ir=1,nrootz                                             3d12s25
              jrsum=irsum+ir-1                                          3d12s25
              do j=0,ncsf(jargb)-1                                      3d12s25
               bc(jrsum)=bc(jrsum)+veck(ipk+j,ir)*gb(jpb+j,ir)          3d12s25
              end do                                                    3d12s25
             end do                                                     3d12s25
             sum=0d0                                                    5d20s21
             do i=1,norb                                                5d21s21
              is=ism(i)                                                 5d20s21
              ig=irel(i)-1+idoubo(is)                                   5d20s21
              jmat=iden(is)+ig*(nbasdws(is)+1)                          3d12s25
              nnn=nbasdws(is)*nbasdws(is)                               3d12s25
              if(btest(i0kc,i))then                                     5d20s21
               do ir=0,nrootz-1                                         3d12s25
                bc(jmat+ir*nnn)=bc(jmat+ir*nnn)+2d0*bc(irsum+ir)        3d12s25
               end do                                                   3d12s25
              end if                                                    5d20s21
              if(btest(i0ko,i))then                                     5d20s21
               do ir=0,nrootz-1                                         3d12s25
                bc(jmat+ir*nnn)=bc(jmat+ir*nnn)+bc(irsum+ir)            3d12s25
               end do                                                   3d12s25
              end if                                                    5d20s21
             end do                                                     5d20s21
             ibcoff=irsum                                               3d12s25
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
             if(jsbk.ne.1)then                                          3d10s25
              write(6,*)('wait a sec ... jsbk = '),jsbk,
     $            (' ne 1 ')                                            3d10s25
              write(6,*)('jsb,jsk '),jsb,jsk,nab1(1),nab1(2)            8d16s22
              stop                                                      5d13s21
             end if                                                     5d13s21
             jgb=irel(nab1(1))-1+idoubo(jsb)                            5d13s21
             jgk=irel(nab1(2))-1+idoubo(jsk)                            5d13s21
             iad=iden(jsb)+jgb+nbasdws(jsb)*jgk                         53d12s25
             nnn=nbasdws(jsb)*nbasdws(jsb)                              3d12s25
             jdvtmp=idvtmp                                              5d20s21
             do ir=1,nrootz                                             5d20s21
              sumd=0d0                                                   3d12s25
              do j=0,ncsf(jargb)-1                                      3d12s25
               sumd=sumd+gb(jpb+j,ir)*bc(jdvtmp+j)                      3d12s25
              end do                                                    3d12s25
              bc(iad)=bc(iad)+sumd                                      3d12s25
              iad=iad+nnn                                               3d12s25
              jdvtmp=jdvtmp+ncsf(jargb)                                 5d20s21
             end do                                                     5d20s21
             ibcoff=idvtmp                                              5d20s21
            else                                                         4d13s21
c     j is bra and i is ket
             ipssx=0                                                    2d6s23
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
      return
      end                                                               7d11s19
