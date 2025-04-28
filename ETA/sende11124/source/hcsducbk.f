c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcsducbk(ihsdiagb,nff1,iff1,ihddiagk,                  8d24s21
     $     nff2,iff2,nsymb,mdon,mdoop,nec,multh,isymbra,                8d24s21
     $     isymket,nvirt,                                               8d24s21
     $     ncsf,ncsf2,irel,ism,irefo,ixw1,ixw2,norb,nrootu,             7d29s21
     $     ixmtf,phase1,lxmt,isymop,n2e,ixmt,i2eop,phase2,ldebug,maxbx, 7d30s21
     $     maxbxd,sr2,idoubo,nbasdws,bc,ibc)                            11d14s22
      implicit real*8 (a-h,o-z)                                         12d18s20
c
c     to do ...
c     consolidate densities, and perhaps njhere in density calculation
c     did I swap iarg and jarg in 4th gandc4?
c
      integer*8 ihsdiagb(mdoop,nsymb),i18,i28,i38,i48,j1c,j1o,i2c,      8d21s21
     $     i2o,itestc,itesto,ipack8,last8(2),                           8d24s21
     $     j2o,k2o,ihddiagk(mdoop,nsymb)                                8d24s21
      integer*1 nab1(2),nab2(2)                                         12d18s20
      external second                                                   2d18s21
      logical l3x,lprt,lnew,ldebug,lchoice                              3d17s21
      equivalence (ipack8,ipack4)                                       12d18s20
      dimension nff1(mdoop,nsymb,2),iff1(*),nff2(mdoop,nsymb),iff2(*),  8d24s21
     $     multh(8,8),nvirt(*),ncsf(*),ncsf2(*),ixmt(8,*),              7d29s21
     $     irel(*),ism(*),irefo(*),isymop(*),itest(32,4),               7d29s21
     $     nab4(2,3),ipack4(2),nl(4),npre(4),mpre(4),i2eop(2,*),           12d21s20
     $     id1visv(8,8),nd1visv(8,8),nokdc(8,8,4),nok3v(4),             12d22s20
     $     idhvnotv(4),ndhvnotv(4),id1vnotv(4,8,8),nd1vnotv(4,8,8),     12d21s20
     $     id3vnotv3(4,2),nd3vnotv3(4,2),loff(4),idkeep(2),ndkeep(2),   1d4s21
     $     mdkeep(2),keep(2),idhvnotvf(4),ibmat(8),ivmat(8),            2d4s21
     $     mdhvnotv(4),md3vnotv3(4,2),nok4f(4),nok4(4),nok3vf(4),       12d22s20
     $     ixmtf(*),idoubo(*),nok33f(4,2),nok33(4,2),tdends(6),          6d9s21
     $     nbasdws(*),nfdat(2),ivtmp(2),jvtmp(2),igdst(2),jgdst(2),      6d13s21
     $     itmpx(4,8),jtmpx(4,8),ntmpx(4,8),mtmpx(4),mfdat(2)           8d24s21
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data loopx/1800000/
      nrootm=nrootu-1                                                   7d29s21
      idoit=0                                                           3d1s21
      mdoo=mdoop-1                                                      7d27s21
      last8(1)=-1                                                       2d8s21
      ircv=ibcoff                                                       1d29s21
      iacc=ircv+mynprocg                                                1d30s21
      igg=iacc+mynprocg                                                 1d30s21
      ibcoff=igg+maxbx                                                  1d30s21
      nacc=0                                                            1d30s21
      call enough(' hcsducbk.  1',bc,ibc)
      itransgg=0                                                        1d30s21
      nsing=0                                                           12d23s20
      norbx=norb+1                                                      12d18s20
      norbxx=norbx+1                                                    12d18s20
      norbxxx=norbxx+1                                                  12d18s20
      jff2=1                                                            6d10s21
      idoit=0                                                           6d10s21
      loop=0
      do isb=1,nsymb                                                    6d10s21
       isbv12=multh(isb,isymket)                                        8d24s21
c
c     let us form bvrkn=[(vv"|nv')+p(k)(vv'|nv")]Vv'v"rk,               2d3s21
c     v ne v' and v"                                                    2d3s21
c     and space for vvrkn=Vvrj*Djkn                                     2d4s21
c     gandc4(Vsv,Vdv'v"), idx=1 (iv'|vv"), idx=2 (iv"|vv')
c
       nvisv=0                                                          6d9s21
       nvnotv=0                                                         6d9s21
       do isbv1=1,nsymb                                                 6d9s21
        isbv2=multh(isbv1,isbv12)                                       6d9s21
        if(isbv2.ge.isbv1)then                                          6d9s21
         if(isbv1.eq.isbv2)then                                         6d9s21
          nvisv=nvisv+nvirt(isbv1)                                      6d9s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         6d9s21
         else                                                           6d9s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 6d9s21
         end if                                                         6d9s21
         nvnotv=nvnotv+nvv                                              6d9s21
        end if                                                          6d9s21
       end do                                                           6d9s21
       do nclo2=mdon,mdoo                                               6d11s21
        nclo2p=nclo2+1                                                  6d11s21
        if(nff2(nclo2p,isb).gt.0)then                                   6d11s21
         iarg=nclo2p-mdon                                               6d11s21
         nfh=ncsf(iarg)*nff2(nclo2p,isb)                                6d11s21
         ibcbmat=ibcoff                                                   2d3s21
         if(n2e.gt.0)then                                               7d29s21
          nbmat=0                                                          3d2s21
          do jsb=1,nsymb                                                    12d18s20
           jsbv=multh(jsb,isymbra)                                      8d24s21
           isn=multh(isbv12,jsbv)                                          2d3s21
           ibmat(jsb)=ibcoff                                               2d3s21
           if(min(irefo(isn),nvirt(jsbv)).gt.0)then                        2d3s21
            nn=nvirt(jsbv)*nrootu*nfh                                    7d29s21
            ibcoff=ibmat(jsb)+nn*irefo(isn)                                3d2s21
            nbmat=nbmat+nn*irefo(isn)                                      3d2s21
           end if                                                          3d2s21
          end do                                                           3d2s21
          call enough(' hcsducbk.  2',bc,ibc)
          do i=ibcbmat,ibcoff-1                                            3d2s21
           bc(i)=0d0                                                       3d2s21
          end do                                                           3d2s21
         end if                                                         7d29s21
         mrowk=nrootu*(ncsf2(iarg)*nvisv+ncsf(iarg)*nvnotv)              7d29s21
         mfdat(1)=ncsf2(iarg)                                           7d29s21
         mfdat(2)=ncsf(iarg)-ncsf2(iarg)                                7d29s21
         nfdat(1)=mfdat(1)*nff2(nclo2p,isb)                             6d16s21
         nfdat(2)=mfdat(2)*nff2(nclo2p,isb)                             6d16s21
         ivd=ibcoff                                                     6d9s21
         igd=ivd+mrowk*nff2(nclo2p,isb)                                 8d24s21
         ibcoff=igd+mrowk*nff2(nclo2p,isb)                              8d24s21
         call enough(' hcsducbk.  3',bc,ibc)
         i18=1                                                          6d9s21
         i28=mrowk                                                      7d29s21
         i38=nff2(nclo2p,isb)                                           6d11s21
         if(min(mrowk,nff2(nclo2p,isb)).gt.0)then                       8d24s21
          call ddi_get(bc,ibc,ihddiagk(nclo2p,isb),i18,i28,i18,i38,     11d15s22
     $        bc(igd))                                                  11d15s22
         end if                                                         8d24s21
c
c     ivd is root, ncsf,nvv,iff2.
c     transpose it to nvv,root,ncsf*iff2 to match cotracted case. :)
c
         ntp=ncsf(iarg)-ncsf2(iarg)                                     7d29s21
         szd=0d0
         do iff=0,nff2(nclo2p,isb)-1                                    6d16s21
          jvd=ivd+mrowk*iff                                             7d29s21
          jvd0=jvd
          jtmp=igd+mrowk*iff                                            7d29s21
          do isbv1=1,nsymb                                               6d9s21
           isbv2=multh(isbv1,isbv12)                                     6d10s21
           if(isbv2.ge.isbv1)then                                        6d9s21
            if(isbv1.eq.isbv2)then                                       6d9s21
             do iv=0,nvirt(isbv1)-1                                      6d9s21
              do i=0,ncsf2(iarg)-1                                      7d29s21
               do ir=0,nrootm                                           7d29s21
                kvd=jvd+iv+nvirt(isbv1)*(ir+nrootu*i)                   7d29s21
                bc(kvd)=bc(jtmp+ir)                                     6d9s21
                szd=szd+bc(kvd)**2
               end do                                                   6d9s21
               jtmp=jtmp+nrootu                                         7d29s21
              end do                                                    6d9s21
             end do                                                     6d9s21
             nhere=nvirt(isbv1)*ncsf2(iarg)*nrootu                      7d29s21
             jvd=jvd+nhere                                               6d16s21
             nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       6d9s21
            else                                                         6d9s21
             nvv=nvirt(isbv1)*nvirt(isbv2)                               6d9s21
            end if                                                       6d9s21
            do ivv=0,nvv-1                                              6d17s21
             do i=0,ncsf(iarg)-1                                        6d17s21
              if(i.lt.ncsf2(iarg))then                                  7d29s21
               l=1                                                      6d17s21
              else                                                      6d17s21
               l=2                                                      6d17s21
              end if                                                    6d17s21
              do ir=0,nrootm                                            7d29s21
               kvd=jvd+ivv+nvv*(ir+nrootu*i)                            6d17s21
               bc(kvd)=bc(jtmp+ir)                                      6d9s21
               szd=szd+bc(kvd)**2
              end do                                                    6d9s21
              jtmp=jtmp+nrootu                                          7d29s21
             end do                                                     6d17s21
            end do                                                      6d17s21
            jvd=jvd+nvv*nrootu*ncsf(iarg)                               7d29s21
           end if                                                        6d9s21
          end do                                                         6d9s21
         end do                                                         6d16s21
         ibcoff=igd                                                     8d24s21
         szd=sqrt(szd/dfloat(mrowk*nff2(nclo2p,isb)))
         if(n2e.gt.0)then                                               7d29s21
          do jsb=1,nsymb                                                    12d18s20
           jsbv=multh(jsb,isymbra)                                       7d29s21
           isn=multh(isbv12,jsbv)                                          2d3s21
           if(min(irefo(isn),nvirt(jsbv)).gt.0)then                        2d3s21
            nnn=nfh*nrootu                                              7d29s21
            nn=nnn*nvirt(jsbv)                                             6d9s21
            ioffvd=ivd                                                   6d9s21
            do isbv1=1,nsymb                                               2d3s21
             isbv2=multh(isbv1,isbv12)                                     2d3s21
             if(isbv2.ge.isbv1)then                                        2d3s21
              call ilimts(irefo(isn),nvirt(isbv1),mynprocg,mynowprog,    6d9s21
     $            il,ih,i1s,i1e,i2s,i2e)                                6d9s21
               nrow=nvirt(jsbv)*nvirt(isbv2)                             7d29s21
              nhere=ih+1-il                                              7d29s21
              ii3=ibcoff                                                 7d29s21
              ibcoff=ii3+nrow*nhere                                      7d29s21
              call enough(' hcsducbk.  4',bc,ibc)
              do iz=ii3,ibcoff-1                                         7d29s21
               bc(iz)=0d0                                                7d29s21
              end do                                                     7d29s21
              do ii2e=1,n2e                                              7d29s21
               if(multh(jsbv,isbv2).eq.isymop(i2eop(1,ii2e)).and.         7d29s21
     $           multh(isn,isbv1).eq.isymop(i2eop(2,ii2e)))then         7d29s21
                itmpr=ibcoff                                             7d29s21
                ibcoff=itmpr+nrow                                        7d29s21
                call enough(' hcsducbk.  5',bc,ibc)
                jtmpr=itmpr                                              7d29s21
                do iv2=0,nvirt(isbv2)-1                                  7d29s21
                 iad=ixmt(jsbv,i2eop(1,ii2e))+idoubo(jsbv)+irefo(jsbv)   7d29s21
     $               +nbasdws(jsbv)*(iv2+idoubo(isbv2)+irefo(isbv2))    7d29s21
                 do jv=0,nvirt(jsbv)-1                                  8d2s21
                  bc(jtmpr+jv)=bc(iad+jv)*phase2                         7d29s21
                 end do                                                  7d29s21
                 jtmpr=jtmpr+nvirt(jsbv)                                8d2s21
                end do                                                   7d29s21
                i1n=irefo(isn)                                           7d29s21
                i10=i1s                                                  7d29s21
                jj3=ii3                                                  7d29s21
                do i2=i2s,i2e                                            7d29s21
                 if(i2.eq.i2e)i1n=i1e                                    7d29s21
                 iad=ixmt(isn,i2eop(2,ii2e))-1+idoubo(isn)               7d29s21
     $               +nbasdws(isn)*(i2-1+idoubo(isbv1)+irefo(isbv1))    7d29s21
                 do i1=i10,i1n                                           7d29s21
                  do j=0,nrow-1                                          7d29s21
                   bc(jj3+j)=bc(jj3+j)+bc(iad+i1)*bc(itmpr+j)
                  end do                                                 7d29s21
                  jj3=jj3+nrow                                           7d29s21
                 end do                                                  7d29s21
                 i10=1                                                   7d29s21
                end do                                                   7d29s21
                ibcoff=itmpr                                             7d29s21
               end if                                                    7d29s21
              end do                                                     7d29s21
              if(isbv1.eq.isbv2)then                                       2d3s21
               if(jsbv.eq.isbv1)then                                       2d3s21
                i10=i1s                                                    2d3s21
                i1n=irefo(isn)                                             2d3s21
                jj3=ii3                                                  7d29s21
                do i2=i2s,i2e                                              2d3s21
                 iv=i2-1                                                   2d3s21
                 if(i2.eq.i2e)i1n=i1e                                      2d3s21
                 do i1=i10,i1n                                             2d3s21
                  jb=ibmat(jsb)+nn*(i1-1)                                6d11s21
                  irow=jj3+nvirt(jsbv)*iv                               8d2s21
                  do jv=0,nvirt(jsbv)-1                                    2d4s21
                   xint=bc(irow+jv)*sr2                                 8d2s21
                   do ir=0,nrootm                                       7d29s21
                    do iff=0,nff2(nclo2p,isb)-1                          6d17s21
                     joffvd=ioffvd+mrowk*iff                            7d29s21
                     do k=0,ncsf2(iarg)-1                               7d29s21
                      iadv=joffvd+iv+nvirt(isbv2)*(ir+nrootu*k)         7d29s21
                      bc(jb+k)=bc(jb+k)+xint*bc(iadv)                       2d3s21
                     end do                                                 2d3s21
                     jb=jb+ncsf(iarg)                                    6d17s21
                    end do                                               6d17s21
                   end do                                                6d11s21
                  end do                                                   2d3s21
                  jj3=jj3+nrow                                             2d3s21
                 end do                                                    2d3s21
                 i10=1                                                     2d3s21
                end do                                                     2d3s21
               else                                                      7d29s21
                i10=i1s                                                    2d3s21
                i1n=irefo(isn)                                             2d3s21
                jj3=ii3                                                 7d29s21
                do i2=i2s,i2e                                              2d3s21
                 iv=i2-1                                                   2d3s21
                 if(i2.eq.i2e)i1n=i1e                                      2d3s21
                 do i1=i10,i1n                                             2d3s21
                  jb=ibmat(jsb)+nn*(i1-1)                                6d11s21
                  do jv=0,nvirt(jsbv)-1                                    2d3s21
                   irow=jj3+jv+nvirt(jsbv)*iv                           7d29s21
                   xint=bc(irow)*sr2                                       2d4s21
                   do ir=0,nrootm                                       7d29s21
                    do iff=0,nff2(nclo2p,isb)-1                          6d17s21
                     joffvd=ioffvd+mrowk*iff                            7d29s21
                     do k=0,ncsf2(iarg)-1                               7d29s21
                      iadv=joffvd+iv+nvirt(isbv2)*(ir+nrootu*k)           6d9s21
                      orig=bc(jb+k)
                      bc(jb+k)=bc(jb+k)+xint*bc(iadv)                     6d9s21
                     end do                                               6d9s21
                     jb=jb+ncsf(iarg)                                    6d17s21
                    end do                                               6d17s21
                   end do                                                  2d3s21
                  end do                                                   2d3s21
                  jj3=jj3+nrow                                             2d3s21
                 end do                                                    2d3s21
                 i10=1                                                     2d3s21
                end do                                                     2d3s21
               end if                                                      2d3s21
               ioffvd=ioffvd+nvirt(isbv2)*nrootu*ncsf2(iarg)            7d29s21
               nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       2d3s21
               isw=0                                                       2d3s21
              else                                                         2d3s21
               nvv=nvirt(isbv1)*nvirt(isbv2)                               2d3s21
               isw=1                                                       2d3s21
              end if                                                       2d3s21
              if(jsbv.eq.isbv2)then                                        2d3s21
               jsw=0                                                       2d3s21
              else                                                         2d3s21
               jsw=1                                                       2d3s21
              end if                                                       2d3s21
              i10=i1s                                                      2d3s21
              i1n=irefo(isn)                                               2d3s21
              jj3=ii3                                                   7d29s21
              if(jsbv.eq.isbv2)then                                        2d3s21
               do i2=i2s,i2e                                                2d3s21
                iv1=i2-1                                                   2d3s21
                ibots=i2                                                   2d3s21
                ibotn=0                                                    2d3s21
                ibot=ibots+isw*(ibotn-ibots)                               2d3s21
                if(i2.eq.i2e)i1n=i1e                                        2d3s21
                do i1=i10,i1n                                               2d3s21
                 do iv2=ibot,nvirt(isbv2)-1                                2d3s21
                  itri=((iv2*(iv2-1))/2)+iv1                               2d3s21
                  irec=iv1+nvirt(isbv1)*iv2                                2d3s21
                  ivv=itri+isw*(irec-itri)                                 2d3s21
                  jb=ibmat(jsb)+nn*(i1-1)                                6d11s21
                  irow=jj3+nvirt(jsbv)*iv2                              8d2s21
                  do jv=0,nvirt(jsbv)-1                                    2d4s21
                   do ir=0,nrootm                                       7d29s21
                    do iff=0,nff2(nclo2p,isb)-1                          6d17s21
                     joffvd=ioffvd+mrowk*iff                            7d29s21
                     do k=0,ncsf(iarg)-1                                 6d17s21
                      iadv=joffvd+ivv+nvv*(ir+nrootu*k)                 7d29s21
                      bc(jb+k)=bc(jb+k)+bc(irow+jv)*bc(iadv)            8d2s21
                     end do                                                 2d3s21
                     jb=jb+ncsf(iarg)                                    6d17s21
                    end do                                               6d17s21
                   end do                                                  2d3s21
                  end do                                                   2d3s21
                 end do                                                    2d3s21
                 jj3=jj3+nrow                                           7d29s21
                end do                                                     2d3s21
                i10=1                                                      2d3s21
               end do                                                      2d3s21
              else                                                      7d29s21
               do i2=i2s,i2e                                                2d3s21
                iv1=i2-1                                                   2d3s21
                ibots=i2                                                   2d3s21
                ibotn=0                                                    2d3s21
                ibot=ibots+isw*(ibotn-ibots)                               2d3s21
                if(i2.eq.i2e)i1n=i1e                                        2d3s21
                do i1=i10,i1n                                               2d3s21
                 do iv2=ibot,nvirt(isbv2)-1                                2d3s21
                  itri=((iv2*(iv2-1))/2)+iv1                               2d3s21
                  irec=iv1+nvirt(isbv1)*iv2                                2d3s21
                  ivv=itri+isw*(irec-itri)                                 2d3s21
                  jb=ibmat(jsb)+nn*(i1-1)                                6d11s21
                  irow=jj3+nvirt(jsbv)*iv2                              8d2s21
                  do jv=0,nvirt(jsbv)-1                                    2d3s21
                   do ir=0,nrootm                                       7d29s21
                    do iff=0,nff2(nclo2p,isb)-1                          6d17s21
                     joffvd=ioffvd+iff*mrowk                            7d29s21
                     do k=0,ncsf(iarg)-1                                 6d17s21
                      iadv=joffvd+ivv+nvv*(ir+nrootu*k)                  6d17s21
                      bc(jb+k)=bc(jb+k)+bc(irow+jv)*bc(iadv)            8d2s21
                     end do                                                 2d3s21
                     jb=jb+ncsf(iarg)                                    6d17s21
                    end do                                               6d17s21
                   end do                                                  2d3s21
                  end do                                                   2d3s21
                 end do                                                    2d3s21
                 jj3=jj3+nrow                                           7d29s21
                end do                                                     2d3s21
                i10=1                                                      2d3s21
               end do                                                      2d3s21
              end if                                                       2d3s21
              ibcoff=ii3                                                7d29s21
              call ilimts(irefo(isn),nvirt(isbv2),mynprocg,mynowprog,il, 6d9s21
     $            ih,i1s,i1e,i2s,i2e)                                   6d9s21
              nhere=ih+1-il                                             7d29s21
               nrow=nvirt(jsbv)*nvirt(isbv1)                            7d29s21
              ii3=ibcoff                                                7d29s21
              ibcoff=ii3+nrow*nhere                                     7d29s21
              call enough(' hcsducbk.  6',bc,ibc)
              do iz=ii3,ibcoff-1                                        7d29s21
               bc(iz)=0d0                                               7d29s21
              end do                                                    7d29s21
              do ii2e=1,n2e                                             7d29s21
               if(multh(jsbv,isbv1).eq.isymop(i2eop(1,ii2e)).and.       7d29s21
     $              multh(isn,isbv2).eq.isymop(i2eop(2,ii2e)))then        7d29s21
                itmpr=ibcoff                                            7d29s21
                ibcoff=itmpr+nrow                                       7d29s21
                call enough(' hcsducbk.  7',bc,ibc)
                jtmpr=itmpr                                             7d29s21
                do iv1=0,nvirt(isbv1)-1                                 7d29s21
                 iad=ixmt(jsbv,i2eop(1,ii2e))+idoubo(jsbv)+irefo(jsbv)  7d29s21
     $                +nbasdws(jsbv)*(iv1+idoubo(isbv1)+irefo(isbv1))   7d29s21
                 do jv=0,nvirt(jsbv)-1                                  8d2s21
                  bc(jtmpr+jv)=phase2*bc(iad+jv)                        7d29s21
                 end do                                                 7d29s21
                 jtmpr=jtmpr+nvirt(jsbv)                                8d2s21
                end do                                                  7d29s21
                i10=i1s                                                 7d29s21
                i1n=irefo(isn)                                          7d29s21
                jj3=ii3                                                 7d29s21
                do i2=i2s,i2e                                           7d29s21
                 if(i2.eq.i2e)i1n=i1e                                   7d29s21
                 iad=ixmt(isn,i2eop(2,ii2e))-1+idoubo(isn)              7d29s21
     $               +nbasdws(isn)*(i2-1+idoubo(isbv2)+irefo(isbv2))    7d29s21
                 do i1=i10,i1n                                          7d29s21
                  do j=0,nrow-1                                         7d29s21
                   bc(jj3+j)=bc(jj3+j)+bc(iad+i1)*bc(itmpr+j)           7d29s21
                  end do                                                7d29s21
                  jj3=jj3+nrow                                          7d29s21
                 end do                                                 7d29s21
                 i10=1                                                  7d29s21
                end do                                                  7d29s21
                ibcoff=itmpr                                            7d29s21
               end if                                                   7d29s21
              end do                                                    7d29s21
              jj3=ii3                                                   7d29s21
              i10=i1s                                                      2d3s21
              i1n=irefo(isn)                                               2d3s21
              if(jsbv.eq.isbv1)then                                        2d4s21
               do i2=i2s,i2e                                                2d3s21
                iv2=i2-1                                                   2d3s21
                itops=iv2-1                                                2d3s21
                itopn=nvirt(isbv1)-1                                       2d3s21
                itop=itops+isw*(itopn-itops)                               2d3s21
                if(i2.eq.i2e)i1n=i1e                                        2d3s21
                do i1=i10,i1n                                               2d3s21
                 do iv1=0,itop                                             2d3s21
                  itri=((iv2*(iv2-1))/2)+iv1                               2d3s21
                  irec=iv1+nvirt(isbv1)*iv2                                2d3s21
                  ivv=itri+isw*(irec-itri)                                 2d3s21
                  jb=ibmat(jsb)+nn*(i1-1)                                6d11s21
                  irow=jj3+nvirt(jsbv)*iv1                              8d2s21
                  do jv=0,nvirt(jsbv)-1                                    2d4s21
                   do ir=0,nrootm                                          2d3s21
                    do iff=0,nff2(nclo2p,isb)-1                          6d17s21
                     joffvd=ioffvd+iff*mrowk                            7d29s21
                     do k=0,ncsf2(iarg)-1                                7d29s21
                      iadv=joffvd+ivv+nvv*(ir+nrootu*k)                  6d17s21
                      bc(jb+k)=bc(jb+k)+bc(irow+jv)*bc(iadv)            8d2s21
                     end do                                                 2d3s21
                     do k=ncsf2(iarg),ncsf(iarg)-1                       7d29s21
                      iadv=joffvd+ivv+nvv*(ir+nrootu*k)                 7d29s21
                      bc(jb+k)=bc(jb+k)-bc(irow+jv)*bc(iadv)            8d2s21
                     end do                                                 2d3s21
                     jb=jb+ncsf(iarg)                                    6d17s21
                    end do                                               6d17s21
                   end do                                                  2d3s21
                  end do                                                   2d3s21
                 end do                                                    2d3s21
                 jj3=jj3+nrow                                           7d29s21
                end do                                                     2d3s21
                i10=1                                                      2d3s21
               end do                                                      2d3s21
              else                                                      7d29s21
               do i2=i2s,i2e                                                2d3s21
                iv2=i2-1                                                   2d3s21
                itops=iv2-1                                                2d3s21
                itopn=nvirt(isbv1)-1                                       2d3s21
                itop=itops+isw*(itopn-itops)                               2d3s21
                if(i2.eq.i2e)i1n=i1e                                        2d3s21
                do i1=i10,i1n                                               2d3s21
                 do iv1=0,itop                                             2d3s21
                  itri=((iv2*(iv2-1))/2)+iv1                               2d3s21
                  irec=iv1+nvirt(isbv1)*iv2                                2d3s21
                  ivv=itri+isw*(irec-itri)                                 2d3s21
                  jb=ibmat(jsb)+nn*(i1-1)                                6d11s21
                  irow=jj3+nvirt(jsbv)*iv1                              8d2s21
                  do jv=0,nvirt(jsbv)-1                                    2d3s21
                   do ir=0,nrootm                                          2d3s21
                    do iff=0,nff2(nclo2p,isb)-1                          6d17s21
                     joffvd=ioffvd+mrowk*iff                            7d29s21
                     do k=0,ncsf2(iarg)-1                                7d29s21
                      iadv=joffvd+ivv+nvv*(ir+nrootu*k)                  6d17s21
                      bc(jb+k)=bc(jb+k)+bc(irow+jv)*bc(iadv)            8d2s21
                     end do                                                 2d3s21
                     do k=ncsf2(iarg),ncsf(iarg)-1                       7d29s21
                      iadv=joffvd+ivv+nvv*(ir+nrootu*k)                  6d17s21
                      bc(jb+k)=bc(jb+k)-bc(irow+jv)*bc(iadv)            8d2s21
                     end do                                                 2d3s21
                     jb=jb+ncsf(iarg)                                    6d17s21
                    end do                                               6d17s21
                   end do                                                  2d3s21
                  end do                                                   2d3s21
                 end do                                                    2d3s21
                 jj3=jj3+nrow                                           7d29s21
                end do                                                     2d3s21
                i10=1                                                      2d3s21
               end do                                                      2d3s21
              end if                                                       2d3s21
              ibcoff=ii3                                                7d29s21
              ioffvd=ioffvd+nvv*nrootu*ncsf(iarg)                        6d17s21
             end if                                                      6d9s21
            end do                                                       6d9s21
            itmp=ibcoff                                                    2d3s21
            ibcoff=itmp+nn*irefo(isn)                                      2d3s21
            call enough(' hcsducbk.  8',bc,ibc)
            jbmat=ibmat(jsb)                                               2d3s21
            do n=0,irefo(isn)-1                                            2d3s21
             do jv=0,nvirt(jsbv)-1                                         2d3s21
              do ir=0,nrootm                                               2d3s21
               do k=0,nfh-1                                                2d3s21
                ka=k/ncsf(iarg)
                kb=k-ka*ncsf(iarg)
                iad=itmp+jv+nvirt(jsbv)*(ir+nrootu*(k+nfh*n))              2d3s21
                bc(iad)=bc(jbmat+k)                                        2d3s21
               end do                                                      2d3s21
               jbmat=jbmat+nfh                                             2d3s21
              end do                                                       2d3s21
             end do                                                        2d3s21
            end do                                                         2d3s21
            do i=0,nn*irefo(isn)-1                                         2d3s21
             bc(ibmat(jsb)+i)=bc(itmp+i)                                   2d3s21
            end do                                                         2d3s21
            ibcoff=itmp                                                    6d10s21
           end if                                                          2d3s21
          end do                                                           2d3s21
          call dws_gsumf(bc(ibcbmat),nbmat)                                3d2s21
         end if                                                         7d29s21
         nopen2=nec-2*nclo2p                                            6d10s21
         nopen2p=nopen2+2                                               6d10s21
         iarg=nclo2p-mdon                                               6d10s21
         mrowk=nrootu*(ncsf2(iarg)*nvisv+ncsf(iarg)*nvnotv)              7d29s21
         jff2top=jff2                                                   6d10s21
         if(isbv12.eq.1)then                                            6d10s21
          ntop=nclo2+3                                                  6d10s21
         else                                                           6d10s21
          ntop=nclo2+2                                                  6d10s21
         end if                                                         6d10s21
         do jsb=1,nsymb                                                 6d10s21
          jsbv=multh(jsb,isymbra)                                       7d29s21
          isv=multh(jsbv,isbv12)                                        7d29s21
          do nclo1=max(mdon,nclo2-2),min(mdoo,ntop)                      6d10s21
           nclo1p=nclo1+1                                               6d10s21
           if(nff1(nclo1p,jsb,1).gt.0)then                              6d10s21
            jarg=nclo1p-mdon                                            6d10s21
            ncol=nff1(nclo1p,jsb,1)*ncsf(jarg)                          6d10s21
            i38=ncol                                                    6d10s21
            nrowb=nrootu*nvirt(jsbv)                                    7d29s21
            igs=ibcoff                                                  8d24s21
            ibcoff=igs+nrowb*ncol                                       8d24s21
            call enough(' hcsducbk.  9',bc,ibc)
            do iz=igs,igs+nrowb*ncol-1                                  7d29s21
             bc(iz)=0d0                                                 7d29s21
            end do                                                      7d29s21
            nopen1=nec-2*nclo1                                          6d10s21
            jff2=jff2top
            do if2=1,nff2(nclo2p,isb)                                      6d10s21
             jvd=ivd+mrowk*(if2-1)                                      7d29s21
             i2c=iff2(jff2)                                                6d10s21
             jff2=jff2+1                                                   6d10s21
             i2o=iff2(jff2)                                                6d10s21
             jff2=jff2+1                                                   6d10s21
             ntest=popcnt(i2c)
             if(ntest.ne.nclo2)then
              write(6,*)('oh no!!! '),ntest,nclo2
              call dws_synca
              call dws_finalize
              stop
             end if
             j2o=ibset(i2o,norbx)                                       6d10s21
             j2o=ibset(j2o,norbxx)                                      6d10s21
             do i=1,norb                                                6d10s21
              itest(i,2)=0                                              6d10s21
             end do                                                     6d10s21
             do i=1,norb                                                6d10s21
              if(btest(i2c,i))then                                      6d10s21
               itest(i,2)=2                                             6d10s21
              end if                                                    6d10s21
              if(btest(j2o,i))then                                      6d10s21
               itest(i,2)=1                                             6d10s21
              end if                                                    6d10s21
             end do                                                     6d10s21
             itest(norbx,2)=1                                           6d10s21
             itest(norbxx,2)=1                                          6d10s21
             k2o=ibset(i2o,norbxx)                                      6d10s21
             k2o=ibset(k2o,norbxxx)                                     6d10s21
             jff1=nff1(nclo1p,jsb,2)                                    6d10s21
             do if1=1,nff1(nclo1p,jsb,1)                                6d10s21
              jgs=igs+nrowb*ncsf(jarg)*(if1-1)                          7d29s21
              j1c=iff1(jff1)                                            6d10s21
              jff1=jff1+1                                               6d10s21
              j1o=iff1(jff1)                                            6d10s21
              j1o=ibset(j1o,norbx)                                      6d10s21
              jff1=jff1+1                                               6d10s21
              idoit=idoit+1                                             6d10s21
              if(mod(idoit,mynprocg).eq.mynowprog)then                  6d10s21
               call gandc4(j1c,j1o,i2c,j2o,nopen1,nopen2p,norbxx,nnot,    12d18s20
     $              nab4,bc,ibc)                                        11d14s22
               if(nnot.eq.2)then                                         6d10s21
                do i=1,norbxx                                              12d19s20
                 itest(i,1)=0                                                 11d25s20
                end do                                                        11d25s20
                do i=1,norb                                                   11d25s20
                 if(btest(j1c,i))then                                         11d25s20
                  itest(i,1)=2                                                11d25s20
                 end if                                                       11d25s20
                end do                                                        11d25s20
                do i=1,norb                                                   11d25s20
                 if(btest(j1o,i))then                                         11d25s20
                  itest(i,1)=1                                                11d25s20
                 end if                                                       11d25s20
                end do                                                        11d25s20
                itest(norbx,1)=1                                           12d19s20
                nok=0
                do i=1,norbxx                                            6d10s21
                 ixn=min(itest(i,1),itest(i,2))
                 if(ixn.gt.0)then                                             11d13s20
                  nok=nok+1                                                   11d13s20
                  itest(nok,3)=ixn                                            11d13s20
                  itest(nok,4)=i                                        6d15s21
                 end if                                                       11d13s20
                end do                                                        11d13s20
                call gandc(j1c,j1o,i2c,j2o,nopen1,nopen2p,jarg,iarg,     6d10s21
     $               ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,       12d18s20
     $               ncsfmid1,bc,ibc)                                   11d14s22
c     if ket is doubles,
c     i.e. norbx=multh(jsb,isymbra),
c     norbxx=multh(norbx,multh(isb,isymket))=jsb*isb*isymbra*isymket
c     =jsbv*isbv12
c     if ket is singles,
c     norbx=multh(jsb,isymket),
c     norbxx=multh(norbx,multh(isb,isymbra))=jsb*isb*isymbra*isymket,
c     i.e. the same.
c
                lsa=ism(nab1(1))                                         6d10s21
                lga=irel(nab1(1))-1+idoubo(lsa)                         7d29s21
                iintxsd=ibcoff                                            6d10s21
                iintxds=iintxsd+nvirt(isv)                              7d29s21
                ibcoff=iintxds+nvirt(isv)                               7d29s21
                call enough(' hcsducbk. 10',bc,ibc)
                if(multh(lsa,isv).eq.isymop(lxmt))then                  7d29s21
                 iprod=ibcoff                                             6d10s21
                 ibcoff=iprod+ncsf(jarg)*ncsf(iarg)                       6d10s21
                 call enough(' hcsducbk. 11',bc,ibc)
                 call prodn(iwpb1,iwpk1,ncsf(jarg),ncsf(iarg),ncsfmid1,   6d10s21
     $                bc(iprod),bc,ibc,1d0,0d0)                         2d13s23
                 ih0sd=ixmtf(lsa)+lga                                   7d29s21
     $                +nbasdws(lsa)*(idoubo(isv)+irefo(isv))            7d29s21
                 ih0ds=ixmtf(isv)+idoubo(isv)+irefo(isv)                7d29s21
     $                +nbasdws(isv)*lga                                 7d29s21
                 do jv=0,nvirt(isv)-1                                   7d29s21
                  bc(iintxds+jv)=bc(ih0ds+jv)*phase1                    7d29s21
                 end do                                                  6d10s21
                else                                                    7d29s21
                 do iz=iintxsd,ibcoff-1                                 7d29s21
                  bc(iz)=0d0                                            7d29s21
                 end do                                                 7d29s21
                end if                                                  7d29s21
                if(n2e.gt.0)then                                        7d29s21
                 do i=1,nok-1                                            6d10s21
                  js=ism(itest(i,4))                                     6d15s21
                  jg=irel(itest(i,4))-1+idoubo(js)                       7d29s21
                  ff=phase2                                             7d29s21
                  if(itest(i,3).eq.2)ff=ff*2d0                          7d29s21
                  do ii2e=1,n2e                                         7d29s21
                   if(isymop(i2eop(1,ii2e)).eq.1.and.                   7d29s21
     $                multh(lsa,isv).eq.isymop(i2eop(2,ii2e)))then      7d29s21
                    ifjj=ixmt(js,i2eop(1,ii2e))+jg*(nbasdws(js)+1)      7d29s21
                    fjj=bc(ifjj)*ff                                     7d29s21
                    iadsd=ixmt(lsa,i2eop(2,ii2e))+lga
     $                   +nbasdws(lsa)*(idoubo(isv)+irefo(isv))         7d29s21
                    iadds=ixmt(isv,i2eop(2,ii2e))+idoubo(isv)+irefo(isv)7d29s21
     $                   +nbasdws(isv)*lga                              7d29s21
                    do jv=0,nvirt(isv)-1                                 7d29s21
                     bc(iintxds+jv)=bc(iintxds+jv)+fjj*bc(iadds+jv)     7d29s21
                    end do                                              7d29s21
                   end if                                               7d29s21
                  end do                                                7d29s21
                  if(itest(i,3).eq.2)then                                  12d14s20
                   do ii2e=1,n2e                                        7d29s21
                    if(multh(js,lsa).eq.isymop(i2eop(1,ii2e)).and.      7d29s21
     $                 multh(js,isv).eq.isymop(i2eop(2,ii2e)))then      7d29s21
                     ifkk=ixmt(lsa,i2eop(1,ii2e))+lga+nbasdws(lsa)*jg   7d30s21
                     fkk=bc(ifkk)*phase2                                7d29s21
                     ifkkb=ixmt(js,i2eop(1,ii2e))+jg+nbasdws(js)*lga    7d30s21
                     fkkb=bc(ifkkb)*phase2                                7d29s21
                     iads=ixmt(isv,i2eop(2,ii2e))+idoubo(isv)           7d30s21
     $                    +irefo(isv)+nbasdws(isv)*jg                   7d29s21
                     iasd=ixmt(js,i2eop(2,ii2e))+jg                     7d30s21
     $                    +nbasdws(js)*(idoubo(isv)+irefo(isv))         7d29s21
                     do jv=0,nvirt(isv)-1                               7d30s21
                      bc(iintxds+jv)=bc(iintxds+jv)-fkkb*bc(iads+jv)    7d30s21
                     end do                                             7d29s21
                    end if                                              7d29s21
                   end do                                               7d29s21
                  end if                                                 6d10s21
                 end do                                                  6d10s21
                end if                                                  7d29s21
                idelta=0
                call hcdsuc1sd(bc(iprod),ncsf(jarg),ncsf(iarg),           6d10s21
     $               ncsf2(iarg),bc(iintxds),bc(jgs),bc(jvd),isv,jsbv,  8d2s21
     $               isbv12,nvirt,multh,nrootu,nsymb,sr2,               8d2s21
     $               .false.,bc,ibc,idelta)                             1d27s23
                ibcoff=iintxsd                                          7d30s21
                if(n2e.gt.0)then                                        7d29s21
                 do i=1,nok                                              12d19s20
                  if(itest(i,3).eq.1)then                                  12d14s20
                   itestc=j1c                                              12d14s20
                   itesto=j1o                                              12d14s20
                   nopenk=nopen1                                                11d13s20
c
c     anihilate common
c
                   if(btest(itestc,itest(i,4)))then                      6d15s21
                    itestc=ibclr(itestc,itest(i,4))                      6d15s21
                    itesto=ibset(itesto,itest(i,4))                      6d15s21
                    karg=jarg-1                                                11d13s20
                    nopenk=nopenk+1                                             11d13s20
                   else                                                         11d13s20
                    itesto=ibclr(itesto,itest(i,4))                      6d15s21
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
                    call gandc(j1c,j1o,itestc,itesto,nopen1,nopenk,        12d14s20
     $                jarg,karg,ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1, 12d28s20
     $                iwpk1,ncsfmid1,bc,ibc)                            11d14s22
                    call gandc(itestc,itesto,i2c,j2o,nopenk,nopen2p,        12d14s20
     $                karg,iarg,ncsf,norbxx,ixw1,ixw2,nnot2,nab2,iwpb2, 12d28s20
     $                iwpk2,ncsfmid2,bc,ibc)                            11d14s22
                    if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                     if(max(nab2(1),nab2(2)).le.norb)then                6d11s21
                      lsa=ism(nab2(1))                                      12d14s20
                      lga=irel(nab2(1))-1+idoubo(lsa)                   7d29s21
                      lsb=ism(nab2(2))                                      12d14s20
                      lgb=irel(nab2(2))-1+idoubo(lsb)                   7d29s21
                      lsc=ism(nab1(1))                                    12d18s20
                      lgc=irel(nab1(1))-1+idoubo(lsc)                   7d29s21
                      iintxds=ibcoff                                    7d29s21
                      iintxsd=iintxds+nvirt(isv)                        7d29s21
                      ibcoff=iintxsd+nvirt(isv)                         7d29s21
                      call enough(' hcsducbk. 12',bc,ibc)
                      do iz=iintxds,ibcoff-1                            7d29s21
                       bc(iz)=0d0                                       7d29s21
                      end do                                            7d29s21
                      do ii2e=1,n2e                                     7d29s21
                       if(multh(lsa,lsb).eq.isymop(i2eop(1,ii2e)).and.  7d29s21
     $                    multh(lsc,isv).eq.isymop(i2eop(2,ii2e)))then  7d29s21
                        iab=ixmt(lsa,i2eop(1,ii2e))+lga+nbasdws(lsa)*lgb7d29s21
                        iba=ixmt(lsb,i2eop(1,ii2e))+lgb+nbasdws(lsb)*lga7d30s21
                        icv=ixmt(lsc,i2eop(2,ii2e))+lgc
     $                       +nbasdws(lsc)*(idoubo(isv)+irefo(isv))     7d29s21
                        ivc=ixmt(isv,i2eop(2,ii2e))+idoubo(isv)
     $                       +irefo(isv)+nbasdws(isv)*lgc               7d29s21
                        fba=bc(iba)*phase2                              7d29s21
                        do jv=0,nvirt(isv)-1                            7d29s21
                         bc(iintxds+jv)=bc(iintxds+jv)+fba*bc(ivc+jv)   7d30s21
                        end do                                          7d29s21
                       end if                                           7d29s21
                      end do                                            7d29s21
                      call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),      6d11s21
     $                  ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod,11d10s22
     $                     bc,ibc)                                      11d10s22
                      call hcdsuc1sd(bc(iprod),ncsf(jarg),ncsf(iarg),           6d10s21
     $               ncsf2(iarg),bc(iintxds),bc(jgs),bc(jvd),isv,jsbv,  8d2s21
     $               isbv12,nvirt,multh,nrootu,nsymb,sr2,               8d2s21
     $               loop.eq.-8940,bc,ibc,idelta)                              12d13s22
                      ibcoff=iintxds                                    7d30s21
                     end if                                              6d11s21
                    end if                                               6d11s21
                   end if                                                6d11s21
                  end if                                                 6d11s21
                 end do                                                  6d11s21
                end if                                                  7d29s21
               else if(nnot.gt.2.and.n2e.gt.0)then                      7d29s21
                if(nnot.eq.3)then                                          12d8s20
                 ipssx=1                                                   12d8s20
                else                                                       12d8s20
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
                 itestc=j1c                                              12d8s20
                 itesto=j1o                                              12d8s20
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
                  call gandc(j1c,j1o,itestc,itesto,nopen1,nopenk,        12d18s20
     $         jarg,karg,ncsf,norbxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,   12d18s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
                  call gandc(itestc,itesto,i2c,j2o,nopenk,nopen2p,         12d8s20
     $         karg,iarg,ncsf,norbxx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,    11d13s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
                  if(nnot1.eq.2.and.nnot2.eq.2)then                      12d19s20
                   if(max(nab1(1),nab1(2)).gt.norb)then                  12d19s20
                    lsa=ism(nab2(1))                                     12d19s20
                    lga=irel(nab2(1))-1+idoubo(lsa)                     7d29s21
                    lsb=ism(nab2(2))                                     12d19s20
                    lgb=irel(nab2(2))-1+idoubo(lsb)                     7d29s21
                    lsc=ism(nab1(1))                                     12d19s20
                    lgc=irel(nab1(1))-1+idoubo(lsc)                     7d29s21
                   else                                                  12d19s20
                    lsa=ism(nab1(1))                                     12d19s20
                    lga=irel(nab1(1))-1+idoubo(lsa)                     7d29s21
                    lsb=ism(nab1(2))                                     12d19s20
                    lgb=irel(nab1(2))-1+idoubo(lsb)                     7d29s21
                    lsc=ism(nab2(1))                                     12d19s20
                    lgc=irel(nab2(1))-1+idoubo(lsc)                     7d29s21
                   end if                                                12d19s20
                   iintxds=ibcoff                                       7d29s21
                   iintxsd=iintxds+nvirt(isv)                           7d29s21
                   ibcoff=iintxsd+nvirt(isv)                            7d29s21
                   call enough(' hcsducbk. 13',bc,ibc)
                   do iz=iintxds,ibcoff-1                               7d29s21
                    bc(iz)=0d0                                          7d29s21
                   end do                                               7d29s21
                   do ii2e=1,n2e                                        7d29s21
                    if(multh(lsa,lsb).eq.isymop(i2eop(1,ii2e)).and.     7d29s21
     $                 multh(lsc,isv).eq.isymop(i2eop(2,ii2e)))then     7d29s21
                     iab=ixmt(lsa,i2eop(1,ii2e))+lga+nbasdws(lsa)*lgb   7d29s21
                     fab=bc(iab)*phase2                                 7d29s21
                     iba=ixmt(lsb,i2eop(1,ii2e))+lgb+nbasdws(lsb)*lga   7d29s21
                     fba=bc(iba)*phase2                                 7d29s21
                     icv=ixmt(lsc,i2eop(2,ii2e))+lgc+nbasdws(lsc)       7d29s21
     $                    *(idoubo(isv)+irefo(isv))                     7d29s21
                     ivc=ixmt(isv,i2eop(2,ii2e))+idoubo(isv)+irefo(isv) 7d29s21
     $                    +nbasdws(isv)*lgc                             7d29s21
                     do jv=0,nvirt(isv)-1                               7d29s21
                      orig=bc(iintxds+jv)
                      bc(iintxds+jv)=bc(iintxds+jv)+fab                 7d29s21
     $                     *bc(icv+nbasdws(lsc)*jv)                     7d29s21
                     end do                                             7d29s21
                    end if                                              7d29s21
                   end do                                               7d29s21
                   call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),        3d19s21
     $                 ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod, 11d10s22
     $                  bc,ibc)                                         11d10s22
                   idelta=0
                   call hcdsuc1sd(bc(iprod),ncsf(jarg),ncsf(iarg),           6d10s21
     $               ncsf2(iarg),bc(iintxds),bc(jgs),bc(jvd),isv,jsbv,  8d2s21
     $               isbv12,nvirt,multh,nrootu,nsymb,sr2,               8d2s21
     $               loop.eq.-8940,bc,ibc,idelta)                              12d13s22
                   ibcoff=iintxds                                       7d30s21
                   if(ipss.eq.2)go to 3                                  12d18s20
                  end if                                                 12d18s20
                 end if                                                 4d20s21
                end do                                                  12d18s20
    3           continue                                                12d18s20
               end if                                                    6d10s21
              end if                                                    6d10s21
              idoit=idoit+1                                             6d10s21
              if(mod(idoit,mynprocg).eq.mynowprog.and.n2e.gt.0)then     7d29s21
               call gandc4(j1c,j1o,i2c,k2o,nopen1,nopen2p,norbxxx,nnot,    12d18s20
     $              nab4,bc,ibc)                                        11d14s22
ccc
               if(nnot.gt.1.and.min(n2e,nvirt(jsbv)).gt.0)then          8d24s21
                if(nnot.eq.3)then                                          12d8s20
                 ipssx=1                                                   12d8s20
                else                                                       12d8s20
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
                 itestc=j1c                                             6d11s21
                 itesto=j1o                                             6d11s21
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
                  call gandc(j1c,j1o,itestc,itesto,nopen1,nopenk,       6d11s21
     $         jarg,karg,ncsf,norbxxx,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,   12d18s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
                  call gandc(itestc,itesto,i2c,k2o,nopenk,nopen2p,      6d11s21
     $         karg,iarg,ncsf,norbxxx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,    11d13s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
                  if(nnot1.eq.2.and.nnot2.eq.2)then                       12d19s20
                   if(nab1(2).eq.norbxx)then                              1d25s21
                    call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),       3d19s21
     $               ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod,bc,11d10s22
     $                  ibc)                                            11d10s22
                    lsa=ism(nab1(1))                                       12d19s20
                    lga=irel(nab1(1))-1                                    12d19s20
                    nn=nvirt(jsbv)*nrootu                               6d11s21
                    itmpb=ibcoff                                        6d11s21
                    itmpd=itmpb+nn*ncsf(iarg)                           6d11s21
                    ibcoff=itmpd+ncsf(jarg)*ncsf(iarg)                  6d11s21
                    call enough(' hcsducbk. 14',bc,ibc)
                    do i=0,ncsf(iarg)-1                                 6d11s21
                     do j=0,ncsf(jarg)-1                                6d11s21
                      ji=iprod+j+ncsf(jarg)*i                           6d11s21
                      ij=itmpd+i+ncsf(iarg)*j                           6d11s21
                      bc(ij)=bc(ji)                                     6d11s21
                     end do                                             6d11s21
                    end do                                              6d11s21
                    jbmat=ibmat(jsb)+nn*(ncsf(iarg)*(if2-1)+nfh*lga)    6d11s21
                    call dgemm('n','n',nn,ncsf(jarg),ncsf(iarg),1d0,
     $                   bc(jbmat),nn,bc(itmpd),ncsf(iarg),1d0,bc(jgs), 6d11s21
     $                   nn,                                            12d13s22
     $'hcsducbka')                                                      12d13s22
                    ibcoff=iprod                                         3d19s21
                   end if                                                 2d25s21
                   if(ipss.eq.2)go to 4                                   12d18s20
                  end if                                                  12d18s20
                 end if                                                 4d20s21
                end do                                                   12d18s20
    4           continue                                                 12d18s20
               end if                                                   6d11s21
ccc
              end if                                                    6d10s21
             end do                                                     6d10s21
            end do                                                         6d10s21
            ivs=ibcoff                                                  8d24s21
            ibcoff=ivs+ncol*nrowb                                       8d24s21
            call enough(' hcsducbk. 15',bc,ibc)
            do icol=0,ncol-1                                            6d10s21
             jgs=igs+icol*nrowb                                         7d29s21
             jvs=ivs+icol*nrowb                                         7d29s21
             do jv=0,nvirt(jsbv)-1                                      6d10s21
              do ir=0,nrootm                                            6d10s21
               ij=jvs+ir+nrootu*jv                                      6d10s21
               ji=jgs+jv+nvirt(jsbv)*ir                                 6d10s21
               bc(ij)=bc(ji)                                            6d10s21
              end do                                                    6d10s21
             end do                                                     6d10s21
            end do                                                      6d10s21
            i28=nrowb                                                   7d29s21
            i38=ncol                                                    6d10s21
            call ddi_acc(bc,ibc,ihsdiagb(nclo1p,jsb),i18,i28,i18,i38,   11d15s22
     $           bc(ivs))                                               11d15s22
            ibcoff=igs                                                  8d24s21
           end if                                                       6d10s21
          end do                                                        6d10s21
         end do                                                         6d10s21
         ibcoff=ibcbmat                                                 6d11s21
        end if                                                          6d11s21
       end do                                                           6d11s21
      end do                                                            12d18s20
      ibcoff=ircv                                                       1d30s21
      return
      end
