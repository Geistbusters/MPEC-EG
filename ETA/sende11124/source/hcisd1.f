c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcisd1(nhand,nff1,iff1,nff0,iff0,ncsfv,ncsf,           3d3s21
     $     mdon,mdoo,nsymb,multh,ixw1,ixw2,iden,nh0av,nvirt,vec,bc,ibc, 1d27s23
     $     iroff)                                                       1d27s23
      implicit real*8 (a-h,o-z)
c
c     1-e density
c
      integer*1 nab1(2),nab2(2)
      integer*8 nhand(mdoo+1,nsymb,2),i1,in,i1c,i1o,i0c,i0o,itesta,     11d26s20
     $     itestb,gandcc,gandco,gandcb                                  2d6s23
      logical ldebug                                                    1d25s21
      dimension nff1(mdoo+1,nsymb,2),iff1(*),nff0(mdoo+1,3),            11d25s20
     $     iff0(*),vec(ncsfv,*),ncsf(*),multh(8,8),nab4(2,3),           3d4s21
     $     iden(*),nvirt(*),iden1x(8,8),nden1x(8,8),nh0av(*)            3d3s21
      include "common.store"                                            11d25s20
      include "common.mrci"                                             11d25s20
      common/kmfind/invk1(2,8,8,8,2)
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      ldebug=.false.                                                    12d6s21
      if(ldebug)then
       write(6,*)('hi, I am the new and improved hcisd1!')                3d21s22
       write(6,*)('in hcisd1, nroot = '),nroot,norb                     3d21s22
       write(6,*)('vec ')
       call prntm2(vec,ncsfv,nroot,ncsfv)
       write(6,*)('f',iffi=1,100000)
       do nclo=mdon,mdoo
        nclop=nclo+1
        jarg=nclop-mdon
        jff0=nff0(nclop,2)
        jvs=nff0(nclop,3)
        write(6,*)('for nclo '),nclo,jarg,jff0,jvs,nff0(nclop,1)
        do i=1,nff0(nclop,1)
         i0c=iff0(jff0)
         jff0=jff0+1
         i0o=iff0(jff0)
         jff0=jff0+1
         call dcbit(i0c,norb,'i0c')
         call dcbit(i0o,norb,'i0o')
         call prntm2(vec(jvs,1),ncsf(iarg),nroot,ncsfv)
        end do
       end do
      end if
      igoal=iden(1)+2
      ibcoffo=ibcoff                                                    12d8s20
      norbx=norb+1                                                      11d25s20
      i1=1                                                              11d25s20
      ifirsttime=0                                                      11d26s20
      do isb=1,nsymb                                                    11d25s20
       isbv=multh(isb,isymmrci)
       if(ldebug)write(6,*)('for symmetry '),isb,isbv                   1d25s21
       do nclo1p=mdon+1,mdoo+1                                           11d25s20
        if(min(nvirt(isbv),nff1(nclo1p,isb,1)).gt.0)then                  11d25s20
         nclo1=nclo1p-1                                                  11d25s20
         if(ldebug)write(6,*)('for nclo1 = '),nclo1                     1d25s21
         iarg=nclo1p-mdon                                                11d25s20
         ncolt=nff1(nclo1p,isb,1)*ncsf(iarg)                            11d10s20
         call ilimts(1,ncolt,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)  12d7s20
         ncolh=ih+1-il                                                  12d7s20
         if(ncolh.gt.0)then                                             1d30s23
          ighere=ibcoff                                                   11d25s20
          ibcoff=ighere+nvirt(isbv)*ncolh*nroot                          11d10s20
          i1c=nvirt(isbv)*nroot                                          11d10s20
          i1o=il                                                         12d7s20
          in=ih                                                          12d7s20
          call ddi_get(bc,ibc,nhand(nclo1p,isb,1),i1,i1c,i1o,in,         11d15s22
     $        bc(ighere))                                               11d15s22
          if(ldebug)then
           write(6,*)('my S vectors '),i1,i1c,i1o,in,nhand(nclo1p,isb,1)
           call prntm2(bc(ighere),nvirt(isbv)*nroot,ncolt,
     $         nvirt(isbv)*nroot)
          end if
c
c     ighere is nroot,nvirt;if1,ncsf (=vecs)
c     we compute idenh as ncsf,nroot,irefo
c     we want nvirt,irefo,nroot
          jghere=ighere                                                  12d2s20
          nn=ncsf(iarg)*nroot                                            11d26s20
          nnm=nn-1                                                       11d26s20
          idenh=ibcoff                                                   11d25s20
          ndenh=idenh+nn*irefo(isbv)                                     11d28s20
          ibcoff=ndenh+irefo(isbv)                                       11d28s20
          call enough('hcisd1.  1',bc,ibc)
          ibctop=ibcoff-1                                                11d26s20
          jff1=nff1(nclo1p,isb,2)                                        11d25s20
          ifcol=1                                                        12d7s20
          do if1=1,nff1(nclo1p,isb,1)                                    11d26s20
           itcol=ifcol+ncsf(iarg)-1                                      11d10s20
           if(ifcol.le.ih.and.itcol.ge.il)then                           12d7s20
            nnlow=max(0,il-ifcol)                                        12d7s20
            nnhgh=min(ncsf(iarg)-1,ih-ifcol)                             12d12s20
            nnuse=nnhgh+1-nnlow                                          12d7s20
            nnuser=nnuse*nroot                                           11d10s20
            do i=idenh,ibctop                                             11d26s20
             bc(i)=0d0                                                    11d26s20
            end do                                                        11d26s20
            i1c=iff1(jff1)                                                11d25s20
            ntest=popcnt(i1c)
            if(ntest.ne.nclo1)then
             write(6,*)('ntest = '),ntest,if1
             call dcbit(i1c,32,'dorb')
             call dws_synca
             call dws_finalize
             stop
            end if
            jff1=jff1+1                                                   11d25s20
            i1o=iff1(jff1)                                                11d25s20
            i1o=ibset(i1o,norbx)                                          11d25s20
            nopen1=popcnt(i1o)                                           3d4s21
            jff1=jff1+1                                                   11d25s20
            ngot=0
            do nclop=max(mdon+1,nclo1p-1),min(mdoo+1,nclo1p+1)           3d3s21
             nclo=nclop-1                                                  11d25s20
             jarg=nclop-mdon                                               11d25s20
             jff0=nff0(nclop,2)                                           11d25s20
             jvs=nff0(nclop,3)                                            11d26s20
             do jf=1,nff0(nclop,1)                                        11d26s20
              i0c=iff0(jff0)                                              11d25s20
              jff0=jff0+1                                                 11d25s20
              i0o=iff0(jff0)                                              11d25s20
              nopen=popcnt(i0o)                                          3d3s21
              jff0=jff0+1                                                 11d25s20
              gandcc=ieor(i1c,i0c)                                      2d6s23
              gandco=ieor(i1o,i0o)                                      2d6s23
              gandcb=ior(gandcc,gandco)                                 2d6s23
              ndifb=popcnt(gandcb)                                      10d21s22
              if(ndifb.le.2)then                                        10d21s22
               ndifs=popcnt(gandco)                                     2d6s23
               ndifd=popcnt(gandcc)                                     2d6s23
               if(ndifs.eq.2.and.ndifb.eq.2)then                        2d6s23
                do i=1,norbx                                            2d6s23
                 if(btest(gandco,i))then                                2d6s23
                  if((btest(i1o,i).and..not.btest(i0c,i)).or.           2d6s23
     $                 (btest(i1c,i).and.btest(i0o,i)))then             2d6s23
                   nab4(1,1)=i                                          2d6s23
                  else                                                  2d6s23
                   nab4(2,1)=i                                          2d6s23
                  end if                                                2d6s23
                 end if                                                 2d6s23
                end do                                                  2d6s23
                call gandc(i1c,i1o,i0c,i0o,nopen1,nopen,iarg,jarg,ncsf,    11d25s20
     $            norbx,ixw1,ixw2,nnot1,nab1,iwpb,iwpk,ncsfmid,bc,ibc)  11d14s22
                jgb=irel(nab4(2,1))-1                                      11d27s20
                jdenh=idenh+nn*jgb                                         11d26s20
                call enough('hcisd1.  2',bc,ibc)
                call xtimesn(ncsf(iarg),ncsf(jarg),ncsfmid,                11d26s20
     $            nroot,iwpb,iwpk,vec(jvs,1),ncsfv,bc(jdenh),ncsf(iarg),3d3s21
     $             1d0,1d0,bc,ibc)                                      11d10s22
               end if                                                   2d6s23
              end if                                                      11d25s20
              jvs=jvs+ncsf(jarg)                                          11d26s20
             end do                                                       11d25s20
            end do                                                        11d25s20
            if(irefo(isbv).gt.0)then                                     3d3s21
             itmph=ibcoff                                                3d4s21
             itmpv=itmph+irefo(isbv)*nnuse                               3d4s21
             ibcoff=itmpv+nvirt(isbv)*nnuse                              3d4s21
             call enough('hcisd1.  3',bc,ibc)
             jden=iden(isbv)+irefo(isbv)                                 3d4s21
             do ir=0,nroot-1                                             3d4s21
              do i=0,irefo(isbv)-1                                       3d4s21
               iad1=idenh+nnlow+ncsf(iarg)*(ir+nroot*i)                  3d4s21
               iad2=itmph+nnuse*i                                        3d4s21
               do j=0,nnuse-1                                            3d4s21
                bc(iad2+j)=bc(iad1+j)                                    3d4s21
               end do                                                    3d4s21
              end do                                                     3d4s21
              do j=0,nnuse-1                                             3d4s21
               iad1=jghere+ir+nroot*nvirt(isbv)*j                        3d4s21
               iad2=itmpv+nvirt(isbv)*j                                  3d4s21
               do iv=0,nvirt(isbv)-1                                     3d4s21
                bc(iad2+iv)=bc(iad1+iv*nroot)                            3d4s21
               end do                                                    3d4s21
              end do                                                     3d4s21
              call dgemm('n','n',nvirt(isbv),irefo(isbv),nnuse,1d0,      3d4s21
     $            bc(itmpv),nvirt(isbv),bc(itmph),nnuse,1d0,            3d4s21
     $            bc(jden),nh0av(isbv),                                 3d4s21
     d' hcisd1.  1')
              jden=jden+iroff*nh0av(isbv)*nh0av(isbv)                    1d27s23
             end do                                                      3d4s21
             ibcoff=itmph                                                3d4s21
            end if                                                       3d3s21
            jghere=jghere+nnuser*nvirt(isbv)                             3d4s21
           else                                                          12d8s20
            jff1=jff1+2                                                  12d8s20
           end if                                                        12d7s20
           ifcol=itcol+1                                                 12d7s20
          end do                                                         11d25s20
          ibcoff=ighere                                                  11d25s20
         end if                                                         1d30s23
        end if                                                          11d25s20
       end do                                                            11d25s20
      end do                                                            11d26s20
      return
      end
