c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcdsuc(ihsdiag,nff1,iff1,ihddiag,nff2,iff2,nsymb,mdon, 6d9s21
     $     mdoo,nec,multh,isymmrci,nvirt,ncsf,ncsf2,irel,ism,irefo,     12d18s20
     $     ismult,ixw1,ixw2,norb,nrootu,ih0av,nh0av,ionex,i3x,i3x3,     12d21s20
     $     nlocald,ldebug,maxbx,maxbxd,tdends,tovr,ionext,sr2,nwiacc,bc,11d10s22
     $     ibc)                                                         11d10s22
      implicit real*8 (a-h,o-z)                                         12d18s20
c
c     to do ...
c     consolidate densities, and perhaps njhere in density calculation
c     did I swap iarg and jarg in 4th gandc4?
c
      integer*8 ihsdiag(mdoo+1,nsymb,2),i18,i28,i38,i48,j1c,j1o,i2c,    12d18s20
     $     i2o,itestc,itesto,ipack8,last8(2),ihddiag(mdoo+1,nsymb,2),   6d10s21
     $     j2o,k2o                                                      6d10s21
      integer*1 nab1(2),nab2(2)                                         12d18s20
      external second                                                   2d18s21
      logical l3x,lprt,lnew,ldebug,lchoice                              3d17s21
      equivalence (ipack8,ipack4)                                       12d18s20
      dimension nff1(mdoo+1,nsymb,2),iff1(*),nff2(mdoo+1,nsymb),iff2(*),6d9s21
     $     multh(8,8),nvirt(*),ncsf(*),ncsf2(4,*),                      6d9s21
     $     irel(*),ism(*),irefo(*),ionex(*),i3x(*),itest(32,4),         6d15s21
     $     nab4(2,3),ipack4(2),nl(4),npre(4),mpre(4),i3x3(*),           12d21s20
     $     id1visv(8,8),nd1visv(8,8),nokdc(8,8,4),nok3v(4),             12d22s20
     $     idhvnotv(4),ndhvnotv(4),id1vnotv(4,8,8),nd1vnotv(4,8,8),     12d21s20
     $     id3vnotv3(4,2),nd3vnotv3(4,2),loff(4),idkeep(2),ndkeep(2),   1d4s21
     $     mdkeep(2),keep(2),idhvnotvf(4),ibmat(8),ivmat(8),            2d4s21
     $     mdhvnotv(4),md3vnotv3(4,2),nok4f(4),nok4(4),nok3vf(4),       12d22s20
     $     ih0av(*),nh0av(*),nok33f(4,2),nok33(4,2),tdends(6),          6d9s21
     $     ionext(*),nfdat(2),ivtmp(2),jvtmp(2),igdst(2),jgdst(2),      6d13s21
     $     itmpx(4,8),jtmpx(4,8),ntmpx(4,8),mtmpx(4),mfdat(2)           6d16s21
      include "common.store"
      common/kmfind/invk1(2,8,8,8,2)                                    6d30s18
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data loopx/89100000/
      nrootm=nrootu-1                                                   1d4s21
      idoit=0                                                           3d1s21
      nsing=0                                                           6d14s21
      ndoub=0                                                           6d14s21
      do isb=1,nsymb                                                    6d14s21
       isbv12=multh(isb,isymmrci)                                       6d14s21
       nvisv=0                                                          6d14s21
       nvnotv=0                                                         6d14s21
       do isbv1=1,nsymb                                                 6d14s21
        isbv2=multh(isbv1,isbv12)                                       6d14s21
        if(isbv2.ge.isbv1)then                                          6d14s21
         if(isbv2.eq.isbv1)then                                         6d14s21
          nvisv=nvisv+nvirt(isbv1)                                      6d14s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         6d14s21
         else                                                           6d14s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 6d14s21
         end if                                                         6d14s21
         nvnotv=nvnotv+nvv                                              6d14s21
        end if                                                          6d14s21
       end do                                                           6d14s21
       do ii=mdon+1,mdoo+1                                              6d14s21
        if(nff1(ii,isb,1).gt.0)then                                       6d14s21
         iarg=ii-mdon                                                   6d14s21
         nsing=nsing+ncsf(iarg)*nvirt(isbv12)*nff1(ii,isb,1)              6d14s21
        end if                                                          6d14s21
        if(nff2(ii,isb).gt.0)then                                       6d14s21
         iarg=ii-mdon                                                   6d14s21
         ndoub=ndoub+(ncsf2(1,iarg)*nvisv+ncsf(iarg)*nvnotv)            6d14s21
     $         *nff2(ii,isb)                                             6d14s21
        end if                                                          6d14s21
       end do                                                           6d14s21
      end do                                                            6d14s21
      last8(1)=-1                                                       2d8s21
      ircv=ibcoff                                                       1d29s21
      iacc=ircv+mynprocg                                                1d30s21
      igs=iacc+mynprocg                                                 1d30s21
      ibcoff=igs+maxbx                                                  1d30s21
      iaccd=ibcoff                                                      8d11s22
      igd=iaccd+mynprocg                                                8d11s22
      ibcoff=igd+maxbxd                                                 8d11s22
      nacc=0                                                            1d30s21
      naccd=0                                                           8d11s22
      call enough(' hcdsuc.  1',bc,ibc)
      itransgg=0                                                        1d30s21
      loop=0
      nsing=0                                                           12d23s20
      norbx=norb+1                                                      12d18s20
      norbxx=norbx+1                                                    12d18s20
      norbxxx=norbxx+1                                                  12d18s20
      jff2=1                                                            6d10s21
      idoit=0                                                           6d10s21
      loop=0
      do isb=1,nsymb                                                    6d10s21
       isbv12=multh(isb,isymmrci)                                       6d10s21
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
         nbmat=0                                                          3d2s21
         do jsb=1,nsymb                                                    12d18s20
          jsbv=multh(jsb,isymmrci)                                         12d18s20
          isn=multh(isbv12,jsbv)                                          2d3s21
          ibmat(jsb)=ibcoff                                               2d3s21
          if(min(irefo(isn),nvirt(jsbv)).gt.0)then                        2d3s21
           nn=nvirt(jsbv)*nrootu*nfh                                      2d3s21
           ibcoff=ibmat(jsb)+nn*irefo(isn)                                3d2s21
           nbmat=nbmat+nn*irefo(isn)                                      3d2s21
          end if                                                          3d2s21
         end do                                                           3d2s21
         ibcvmat=ibcoff                                                   3d2s21
         do jsb=1,nsymb                                                    12d18s20
          jsbv=multh(jsb,isymmrci)                                         12d18s20
          isn=multh(isbv12,jsbv)                                          2d3s21
          ivmat(jsb)=ibcoff                                               2d3s21
          if(min(irefo(isn),nvirt(jsbv)).gt.0)then                        2d3s21
           nn=nvirt(jsbv)*nrootu*nfh                                      2d3s21
           ibcoff=ivmat(jsb)+nn*irefo(isn)                                3d2s21
          end if                                                          3d2s21
         end do                                                           3d2s21
         call enough(' hcdsuc.  2',bc,ibc)
         do i=ibcbmat,ibcoff-1                                            3d2s21
          bc(i)=0d0                                                       3d2s21
         end do                                                           3d2s21
         mrow=nrootu*(ncsf2(1,iarg)*nvisv+ncsf(iarg)*nvnotv)            6d9s21
         mfdat(1)=ncsf2(1,iarg)                                         6d16s21
         mfdat(2)=ncsf(iarg)-ncsf2(1,iarg)                              6d16s21
         nfdat(1)=mfdat(1)*nff2(nclo2p,isb)                             6d16s21
         nfdat(2)=mfdat(2)*nff2(nclo2p,isb)                             6d16s21
         ivd=ibcoff                                                     6d9s21
         ibcoff=ivd+mrow*nff2(nclo2p,isb)                               8d11s22
         call enough(' hcdsuc.  3',bc,ibc)
         i18=1                                                          6d9s21
         i28=mrow                                                       6d9s21
         i38=nff2(nclo2p,isb)                                           6d11s21
         call ddi_get(bc,ibc,ihddiag(nclo2p,isb,1),i18,i28,i18,i38,     11d15s22
     $        bc(igd))                                                  11d15s22
c
c     ivd is root, ncsf,nvv,iff2.
c     transpose it to nvv,root,ncsf*iff2 to match cotracted case. :)
c
         ntp=ncsf(iarg)-ncsf2(1,iarg)                                   6d10s21
         szd=0d0
         do iff=0,nff2(nclo2p,isb)-1                                    6d16s21
          jvd=ivd+mrow*iff                                              6d16s21
          jvd0=jvd
          jtmp=igd+mrow*iff                                             6d16s21
          do isbv1=1,nsymb                                               6d9s21
           isbv2=multh(isbv1,isbv12)                                     6d10s21
           if(isbv2.ge.isbv1)then                                        6d9s21
            if(isbv1.eq.isbv2)then                                       6d9s21
             do iv=0,nvirt(isbv1)-1                                      6d9s21
              do i=0,ncsf2(1,iarg)-1                                    6d9s21
               do ir=0,nrootm                                           6d9s21
                kvd=jvd+iv+nvirt(isbv1)*(ir+nrootu*i)                   6d16s21
                bc(kvd)=bc(jtmp+ir)                                     6d9s21
                szd=szd+bc(kvd)**2
               end do                                                   6d9s21
               jtmp=jtmp+nrootu                                         6d9s21
              end do                                                    6d9s21
             end do                                                     6d9s21
             nhere=nvirt(isbv1)*ncsf2(1,iarg)*nrootu                     6d9s21
             jvd=jvd+nhere                                               6d16s21
             nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       6d9s21
            else                                                         6d9s21
             nvv=nvirt(isbv1)*nvirt(isbv2)                               6d9s21
            end if                                                       6d9s21
            do ivv=0,nvv-1                                              6d17s21
             do i=0,ncsf(iarg)-1                                        6d17s21
              if(i.lt.ncsf2(1,iarg))then                                6d17s21
               l=1                                                      6d17s21
              else                                                      6d17s21
               l=2                                                      6d17s21
              end if                                                    6d17s21
              do ir=0,nrootm                                            6d17s21
               kvd=jvd+ivv+nvv*(ir+nrootu*i)                            6d17s21
               bc(kvd)=bc(jtmp+ir)                                      6d9s21
               szd=szd+bc(kvd)**2
              end do                                                    6d9s21
              jtmp=jtmp+nrootu                                          6d17s21
             end do                                                     6d17s21
            end do                                                      6d17s21
            jvd=jvd+nvv*nrootu*ncsf(iarg)                               6d17s21
           end if                                                        6d9s21
          end do                                                         6d9s21
         end do                                                         6d16s21
         szd=sqrt(szd/dfloat(mrow*nff2(nclo2p,isb)))
         do jsb=1,nsymb                                                    12d18s20
          jsbv=multh(jsb,isymmrci)                                         12d18s20
          isn=multh(isbv12,jsbv)                                          2d3s21
          if(min(irefo(isn),nvirt(jsbv)).gt.0)then                        2d3s21
           nnn=nfh*nrootu                                                 6d9s21
           nn=nnn*nvirt(jsbv)                                             6d9s21
           ioffvd=ivd                                                   6d9s21
           do isbv1=1,nsymb                                               2d3s21
            isbv2=multh(isbv1,isbv12)                                     2d3s21
            if(isbv2.ge.isbv1)then                                        2d3s21
             call ilimts(irefo(isn),nvirt(isbv1),mynprocg,mynowprog,    6d9s21
     $            il,ih,i1s,i1e,i2s,i2e)                                6d9s21
             i2eu=invk1(1,jsbv,isbv2,isn,2)                               2d3s21
             icase=invk1(2,jsbv,isbv2,isn,2)                              2d3s21
             ii3=i3x(i2eu)                                                2d3s21
             if(isbv1.eq.isbv2)then                                       2d3s21
              if(jsbv.eq.isbv1)then                                       2d3s21
               nrow=(nvirt(jsbv)*(nvirt(jsbv)+1))/2                       2d3s21
               i10=i1s                                                    2d3s21
               i1n=irefo(isn)                                             2d3s21
               do i2=i2s,i2e                                              2d3s21
                iv=i2-1                                                   2d3s21
                if(i2.eq.i2e)i1n=i1e                                      2d3s21
                do i1=i10,i1n                                             2d3s21
                 jb=ibmat(jsb)+nn*(i1-1)                                6d11s21
                 do jv=0,nvirt(jsbv)-1                                    2d4s21
                  ix=max(jv,iv)                                           2d4s21
                  in=min(jv,iv)                                           2d4s21
                  irow=ii3+((ix*(ix+1))/2)+in                             2d4s21
                  xint=bc(irow)*sr2                                       2d4s21
                  do ir=0,nrootm                                        6d11s21
                   do iff=0,nff2(nclo2p,isb)-1                          6d17s21
                    joffvd=ioffvd+mrow*iff                              6d17s21
                    do k=0,ncsf2(1,iarg)-1                              6d17s21
                     iadv=joffvd+iv+nvirt(isbv2)*(ir+nrootu*k)           6d11s21
                     orig=bc(jb+k)
                     bc(jb+k)=bc(jb+k)+xint*bc(iadv)                       2d3s21
                    end do                                                 2d3s21
                    jb=jb+ncsf(iarg)                                    6d17s21
                   end do                                               6d17s21
                  end do                                                6d11s21
                 end do                                                   2d3s21
                 ii3=ii3+nrow                                             2d3s21
                end do                                                    2d3s21
                i10=1                                                     2d3s21
               end do                                                     2d3s21
              else if(icase.eq.1)then                                          2d3s21
               nrow=nvirt(jsbv)*nvirt(isbv1)                              2d3s21
               i10=i1s                                                    2d3s21
               i1n=irefo(isn)                                             2d3s21
               do i2=i2s,i2e                                              2d3s21
                iv=i2-1                                                   2d3s21
                if(i2.eq.i2e)i1n=i1e                                      2d3s21
                do i1=i10,i1n                                             2d3s21
                 jb=ibmat(jsb)+nn*(i1-1)                                6d11s21
                 do jv=0,nvirt(jsbv)-1                                    2d3s21
                  irow=ii3+jv+nvirt(jsbv)*iv                              2d3s21
                  xint=bc(irow)*sr2                                       2d4s21
                  do ir=0,nrootm                                           2d3s21
                   do iff=0,nff2(nclo2p,isb)-1                          6d17s21
                    joffvd=ioffvd+mrow*iff                              6d17s21
                    do k=0,ncsf2(1,iarg)-1                              6d17s21
                     iadv=joffvd+iv+nvirt(isbv2)*(ir+nrootu*k)           6d9s21
                     orig=bc(jb+k)
                     bc(jb+k)=bc(jb+k)+xint*bc(iadv)                     6d9s21
                    end do                                               6d9s21
                    jb=jb+ncsf(iarg)                                    6d17s21
                   end do                                               6d17s21
                  end do                                                  2d3s21
                 end do                                                   2d3s21
                 ii3=ii3+nrow                                             2d3s21
                end do                                                    2d3s21
                i10=1                                                     2d3s21
               end do                                                     2d3s21
              else                                                        2d3s21
               nrow=nvirt(jsbv)*nvirt(isbv1)                              2d3s21
               i10=i1s                                                    2d3s21
               i1n=irefo(isn)                                             2d3s21
               do i2=i2s,i2e                                              2d3s21
                iv=i2-1                                                   2d3s21
                if(i2.eq.i2e)i1n=i1e                                      2d3s21
                do i1=i10,i1n                                             2d3s21
                 jb=ibmat(jsb)+nn*(i1-1)                                6d11s21
                 do jv=0,nvirt(jsbv)-1                                    2d3s21
                  irow=ii3+iv+nvirt(isbv1)*jv                              2d3s21
                  xint=bc(irow)*sr2                                       2d4s21
                  do ir=0,nrootm                                           2d3s21
                   do iff=0,nff2(nclo2p,isb)-1                          6d17s21
                    joffvd=ioffvd+iff*mrow                              6d17s21
                    do k=0,ncsf2(1,iarg)-1                              6d17s21
                     iadv=joffvd+iv+nvirt(isbv2)*(ir+nrootu*k)           6d9s21
                     orig=bc(jb+k)
                     bc(jb+k)=bc(jb+k)+xint*bc(iadv)                     6d9s21
                    end do                                                6d9s21
                    jb=jb+ncsf(iarg)                                    6d17s21
                   end do                                                  2d3s21
                  end do                                                6d17s21
                 end do                                                   2d3s21
                 ii3=ii3+nrow                                             2d3s21
                end do                                                    2d3s21
                i10=1                                                     2d3s21
               end do                                                     2d3s21
              end if                                                      2d3s21
              ioffvd=ioffvd+nvirt(isbv2)*nrootu*ncsf2(1,iarg)           6d17s21
              nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       2d3s21
              isw=0                                                       2d3s21
             else                                                         2d3s21
              nvv=nvirt(isbv1)*nvirt(isbv2)                               2d3s21
              isw=1                                                       2d3s21
             end if                                                       2d3s21
             if(jsbv.eq.isbv2)then                                        2d3s21
              nrow=(nvirt(jsbv)*(nvirt(jsbv)+1))/2                        2d3s21
              jsw=0                                                       2d3s21
             else                                                         2d3s21
              nrow=nvirt(jsbv)*nvirt(isbv2)                               2d3s21
              jsw=1                                                       2d3s21
             end if                                                       2d3s21
             i10=i1s                                                      2d3s21
             i1n=irefo(isn)                                               2d3s21
             ii3=i3x(i2eu)                                                2d3s21
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
                 do jv=0,nvirt(jsbv)-1                                    2d4s21
                  ix=max(iv2,jv)                                          2d4s21
                  in=min(iv2,jv)                                          2d4s21
                  itri=((ix*(ix+1))/2)+in                                 2d4s21
                  irec=jv+nvirt(jsbv)*iv2                                 2d3s21
                  irow=ii3+itri+jsw*(irec-itri)                           2d3s21
                  do ir=0,nrootm                                          2d3s21
                   do iff=0,nff2(nclo2p,isb)-1                          6d17s21
                    joffvd=ioffvd+mrow*iff                              6d17s21
                    do k=0,ncsf(iarg)-1                                 6d17s21
                     iadv=joffvd+ivv+nvv*(ir+nrootu*k)                     2d3s21
                     bc(jb+k)=bc(jb+k)+bc(irow)*bc(iadv)                 6d9s21
                    end do                                                 2d3s21
                    jb=jb+ncsf(iarg)                                    6d17s21
                   end do                                               6d17s21
                  end do                                                  2d3s21
                 end do                                                   2d3s21
                end do                                                    2d3s21
                ii3=ii3+nrow                                              2d3s21
               end do                                                     2d3s21
               i10=1                                                      2d3s21
              end do                                                      2d3s21
             else if(icase.eq.1)then                                      2d3s21
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
                 do jv=0,nvirt(jsbv)-1                                    2d3s21
                  ix=max(jv,iv2)                                          2d3s21
                  in=min(jv,iv2)                                          2d3s21
                  itri=((ix*(ix+1))/2)+in                                 2d3s21
                  irec=jv+nvirt(jsbv)*iv2                                 2d3s21
                  irow=ii3+itri+jsw*(irec-itri)                           2d3s21
                  do ir=0,nrootm                                          2d3s21
                   do iff=0,nff2(nclo2p,isb)-1                          6d17s21
                    joffvd=ioffvd+iff*mrow                              6d17s21
                    do k=0,ncsf(iarg)-1                                 6d17s21
                     iadv=joffvd+ivv+nvv*(ir+nrootu*k)                  6d17s21
                     bc(jb+k)=bc(jb+k)+bc(irow)*bc(iadv)                 6d10s21
                    end do                                                 2d3s21
                    jb=jb+ncsf(iarg)                                    6d17s21
                   end do                                               6d17s21
                  end do                                                  2d3s21
                 end do                                                   2d3s21
                end do                                                    2d3s21
                ii3=ii3+nrow                                              2d3s21
               end do                                                     2d3s21
               i10=1                                                      2d3s21
              end do                                                      2d3s21
             else                                                         2d3s21
              i10=i1s                                                      2d3s21
              i1n=irefo(isn)                                               2d3s21
              ii3=i3x(i2eu)                                               2d3s21
              do i2=i2s,i2e                                                2d3s21
               iv1=i2-1                                                   2d3s21
               ibots=i2                                                   2d3s21
               ibotn=0                                                    2d3s21
               ibot=ibots+isw*(ibotn-ibots)                               2d3s21
               if(i2.eq.i2e)i1n=i1e                                        2d3s21
               do i1=i10,i1n                                               2d3s21
                i1m=i1-1                                                  2d3s21
                do jv=0,nvirt(jsbv)-1                                     2d3s21
                 do iv2=ibot,nvirt(isbv2)-1                               2d3s21
                  itri=((iv2*(iv2-1))/2)+iv1                               2d3s21
                  irec=iv1+nvirt(isbv1)*iv2                                2d3s21
                  ivv=itri+isw*(irec-itri)                                 2d3s21
                  jb=ibmat(jsb)+nnn*(jv+nvirt(jsbv)*i1m)                6d11s21
                  irow=ii3+iv2+nvirt(isbv2)*jv                            2d3s21
                  do ir=0,nrootm                                          2d3s21
                   do iff=0,nff2(nclo2p,isb)-1                          6d17s21
                    joffvd=ioffvd+mrow*iff                              6d17s21
                    do k=0,ncsf(iarg)-1                                 6d17s21
                     iadv=joffvd+ivv+nvv*(ir+nrootu*k)                  6d17s21
                     bc(jb+k)=bc(jb+k)+bc(irow)*bc(iadv)                 6d9s21
                    end do                                                 2d3s21
                    jb=jb+ncsf(iarg)                                    6d17s21
                   end do                                               6d17s21
                  end do                                                  2d3s21
                 end do                                                   2d3s21
                end do                                                    2d3s21
                ii3=ii3+nrow                                              2d3s21
               end do                                                     2d3s21
               i10=1                                                      2d3s21
              end do                                                      2d3s21
             end if                                                       2d3s21
             call ilimts(irefo(isn),nvirt(isbv2),mynprocg,mynowprog,il, 6d9s21
     $            ih,i1s,i1e,i2s,i2e)                                   6d9s21
             i2eu=invk1(1,jsbv,isbv1,isn,2)                               2d3s21
             icase=invk1(2,jsbv,isbv1,isn,2)                              2d3s21
             ii3=i3x(i2eu)                                                2d3s21
             i10=i1s                                                      2d3s21
             i1n=irefo(isn)                                               2d3s21
             if(jsbv.eq.isbv1)then                                        2d3s21
              nrow=(nvirt(jsbv)*(nvirt(jsbv)+1))/2                        2d3s21
              jsw=0                                                       2d3s21
             else                                                         2d3s21
              nrow=nvirt(jsbv)*nvirt(isbv1)                               2d3s21
              jsw=1                                                       2d3s21
             end if                                                       2d3s21
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
                 do jv=0,nvirt(jsbv)-1                                    2d4s21
                  ix=max(jv,iv1)                                          2d4s21
                  in=min(jv,iv1)                                          2d4s21
                  itri=((ix*(ix+1))/2)+in                                 2d4s21
                  irec=jv+nvirt(jsbv)*iv1                                 2d3s21
                  irow=ii3+itri+jsw*(irec-itri)                           2d3s21
                  do ir=0,nrootm                                          2d3s21
                   do iff=0,nff2(nclo2p,isb)-1                          6d17s21
                    joffvd=ioffvd+iff*mrow                              6d17s21
                    do k=0,ncsf2(1,iarg)-1                              6d17s21
                     iadv=joffvd+ivv+nvv*(ir+nrootu*k)                  6d17s21
                     bc(jb+k)=bc(jb+k)+bc(irow)*bc(iadv)                 6d9s21
                    end do                                                 2d3s21
                    do k=ncsf2(1,iarg),ncsf(iarg)-1                     6d17s21
                     iadv=joffvd+ivv+nvv*(ir+nrootu*k)                  6d17s21
                     bc(jb+k)=bc(jb+k)-bc(irow)*bc(iadv)                    2d3s21
                    end do                                                 2d3s21
                    jb=jb+ncsf(iarg)                                    6d17s21
                   end do                                               6d17s21
                  end do                                                  2d3s21
                 end do                                                   2d3s21
                end do                                                    2d3s21
                ii3=ii3+nrow                                              2d3s21
               end do                                                     2d3s21
               i10=1                                                      2d3s21
              end do                                                      2d3s21
             else if(icase.eq.1)then                                      2d4s21
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
                 do jv=0,nvirt(jsbv)-1                                    2d3s21
                  ix=max(jv,iv1)                                          2d3s21
                  in=min(jv,iv1)                                          2d3s21
                  itri=((ix*(ix+1))/2)+in                                 2d3s21
                  irec=jv+nvirt(jsbv)*iv1                                 2d3s21
                  irow=ii3+itri+jsw*(irec-itri)                           2d3s21
                  do ir=0,nrootm                                          2d3s21
                   do iff=0,nff2(nclo2p,isb)-1                          6d17s21
                    joffvd=ioffvd+mrow*iff                              6d17s21
                    do k=0,ncsf2(1,iarg)-1                              6d17s21
                     iadv=joffvd+ivv+nvv*(ir+nrootu*k)                  6d17s21
                     bc(jb+k)=bc(jb+k)+bc(irow)*bc(iadv)                    2d3s21
                    end do                                                 2d3s21
                    do k=ncsf2(1,iarg),ncsf(iarg)-1                     6d17s21
                     iadv=joffvd+ivv+nvv*(ir+nrootu*k)                  6d17s21
                     bc(jb+k)=bc(jb+k)-bc(irow)*bc(iadv)                 6d9s21
                    end do                                                 2d3s21
                    jb=jb+ncsf(iarg)                                    6d17s21
                   end do                                               6d17s21
                  end do                                                  2d3s21
                 end do                                                   2d3s21
                end do                                                    2d3s21
                ii3=ii3+nrow                                              2d3s21
               end do                                                     2d3s21
               i10=1                                                      2d3s21
              end do                                                      2d3s21
             else                                                         2d3s21
              do i2=i2s,i2e                                                2d3s21
               iv2=i2-1                                                   2d3s21
               itops=iv2-1                                                2d3s21
               itopn=nvirt(isbv1)-1                                       2d3s21
               itop=itops+isw*(itopn-itops)                               2d3s21
               if(i2.eq.i2e)i1n=i1e                                        2d3s21
               do i1=i10,i1n                                               2d3s21
                i1m=i1-1                                                  2d3s21
                do jv=0,nvirt(jsbv)-1                                     2d3s21
                 do iv1=0,itop                                            2d3s21
                  itri=((iv2*(iv2-1))/2)+iv1                               2d3s21
                  irec=iv1+nvirt(isbv1)*iv2                                2d3s21
                  ivv=itri+isw*(irec-itri)                                 2d3s21
                  jb=ibmat(jsb)+nnn*(jv+nvirt(jsbv)*i1m)                6d11s21
                  irow=ii3+iv1+nvirt(isbv1)*jv                            2d3s21
                  do ir=0,nrootm                                          2d3s21
                   do iff=0,nff2(nclo2p,isb)-1                          6d17s21
                    joffvd=ioffvd+mrow*iff                              6d17s21
                    do k=0,ncsf2(1,iarg)-1                              6d17s21
                     iadv=joffvd+ivv+nvv*(ir+nrootu*k)                  6d17s21
                     bc(jb+k)=bc(jb+k)+bc(irow)*bc(iadv)                 6d9s21
                    end do                                                 2d3s21
                    do k=ncsf2(1,iarg),ncsf(iarg)-1                     6d17s21
                     iadv=joffvd+ivv+nvv*(ir+nrootu*k)                  6d17s21
                     bc(jb+k)=bc(jb+k)-bc(irow)*bc(iadv)                 6d9s21
                    end do                                                 2d3s21
                    jb=jb+ncsf(iarg)                                    6d17s21
                   end do                                               6d17s21
                  end do                                                  2d3s21
                 end do                                                   2d3s21
                end do                                                    2d3s21
                ii3=ii3+nrow                                              2d3s21
               end do                                                     2d3s21
               i10=1                                                      2d3s21
              end do                                                      2d3s21
             end if                                                       2d3s21
             ioffvd=ioffvd+nvv*nrootu*ncsf(iarg)                        6d17s21
            end if                                                      6d9s21
           end do                                                       6d9s21
           itmp=ibcoff                                                    2d3s21
           ibcoff=itmp+nn*irefo(isn)                                      2d3s21
           call enough(' hcdsuc.  4',bc,ibc)
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
         nopen2=nec-2*nclo2p                                            6d10s21
         nopen2p=nopen2+2                                               6d10s21
         iarg=nclo2p-mdon                                               6d10s21
         mrow=nrootu*(ncsf2(1,iarg)*nvisv+ncsf(iarg)*nvnotv)            6d9s21
         call ddi_done(ibc(iaccd),naccd)                                8d11s22
         do i=0,mrow*nff2(nclo2p,isb)-1
          bc(igd+i)=0d0                                                 6d16s21
         end do                                                         6d16s21
         jff2top=jff2                                                   6d10s21
         if(isbv12.eq.1)then                                            6d10s21
          ntop=nclo2+3                                                  6d10s21
         else                                                           6d10s21
          ntop=nclo2+2                                                  6d10s21
         end if                                                         6d10s21
         do jsb=1,nsymb                                                 6d10s21
          jsbv=multh(jsb,isymmrci)                                      6d10s21
          do nclo1=max(mdon,nclo2-2),min(mdoo,ntop)                      6d10s21
           nclo1p=nclo1+1                                               6d10s21
           if(nff1(nclo1p,jsb,1).gt.0)then                              6d10s21
            jarg=nclo1p-mdon                                            6d10s21
            ncol=nff1(nclo1p,jsb,1)*ncsf(jarg)                          6d10s21
            i38=ncol                                                    6d10s21
            nrow=nrootu*nvirt(jsbv)                                     6d10s21
            i28=nrow                                                    6d10s21
            ivs=ibcoff                                                  6d10s21
            ivst=ivs+nrow*ncol                                          8d11s22
            ibcoff=ivst+nrow*ncol                                       8d11s22
            call enough(' hcdsuc.  5',bc,ibc)
            call ddi_get(bc,ibc,ihsdiag(nclo1p,jsb,1),i18,i28,i18,i38,  11d15s22
     $           bc(ivst))                                              11d15s22
c
c     vectors are stored nroot,nvirt,ncol
c     transpose to nvirt,nroot,ncol
c
            do icol=0,i38-1                                             6d10s21
             ia=icol/ncsf(jarg)
             ib=icol-ia*ncsf(jarg)
             jvst=ivst+icol*nrow                                        8d11s22
             jvs=ivs+icol*nrow                                          6d10s21
             do jv=0,nvirt(jsbv)-1                                      6d10s21
              do ir=0,nrootm                                            6d10s21
               ij=jvst+ir+nrootu*jv                                     8d11s22
               ji=jvs+jv+nvirt(jsbv)*ir                                 6d10s21
               bc(ji)=bc(ij)                                            6d10s21
              end do                                                    6d10s21
             end do                                                     6d10s21
            end do                                                      6d10s21
            ibcoff=ivst                                                 8d11s22
            nopen1=nec-2*nclo1                                          6d10s21
            jff2=jff2top
            do if2=1,nff2(nclo2p,isb)                                      6d10s21
             jvd=ivd+mrow*(if2-1)                                       6d10s21
             jgd=igd+mrow*(if2-1)                                       6d10s21
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
              jvs=ivs+nrow*ncsf(jarg)*(if1-1)                           6d15s21
              jgs=igs+nrow*ncsf(jarg)*(if1-1)                           6d15s21
              j1c=iff1(jff1)                                            6d10s21
              jff1=jff1+1                                               6d10s21
              j1o=iff1(jff1)                                            6d10s21
              j1o=ibset(j1o,norbx)                                      6d10s21
              jff1=jff1+1                                               6d10s21
              idoit=idoit+1                                             6d10s21
              if(mod(idoit,mynprocg).eq.mynowprog)then                  6d10s21
               call gandc4(j1c,j1o,i2c,j2o,nopen1,nopen2p,norbxx,nnot,    12d18s20
     $              nab4,bc,ibc)                                        11d14s22
               loop=loop+1
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
                lsa=ism(nab1(1))                                         6d10s21
                lga=irel(nab1(1))-1                                      6d10s21
                iprod=ibcoff                                             6d10s21
                ibcoff=iprod+ncsf(jarg)*ncsf(iarg)                       6d10s21
                call enough(' hcdsuc.  6',bc,ibc)
                call prodn(iwpb1,iwpk1,ncsf(jarg),ncsf(iarg),ncsfmid1,   6d10s21
     $                bc(iprod),bc,ibc,1d0,0d0)                         2d13s23
                iintx=ibcoff                                            6d10s21
                ibcoff=iintx+nvirt(lsa)                                 6d15s21
                call enough(' hcdsuc.  7',bc,ibc)
                ih0=ih0av(lsa)+irefo(lsa)+nh0av(lsa)*lga                6d15s21
                do jv=0,nvirt(lsa)-1                                    6d15s21
                 bc(iintx+jv)=bc(ih0+jv)                                6d10s21
                end do                                                  6d10s21
                do i=1,nok-1                                            6d10s21
                 js=ism(itest(i,4))                                     6d15s21
                 jg=irel(itest(i,4))-1                                  6d15s21
                 i2eu=invk1(1,js,js,lsa,2)                              6d10s21
                 ncolj=(irefo(js)*(irefo(js)+1))/2                      6d10s21
                 icolj=((jg*(jg+1))/2)+jg+ncolj*lga                     6d10s21
                 jint=ionext(i2eu)+nvirt(lsa)*icolj                     6d15s21
                 if(itest(i,3).eq.2)then                                  12d14s20
                  if(lsa.ne.js)then                                     6d10s21
                   i2eu=invk1(1,js,lsa,js,2)                            6d10s21
                   icase=invk1(2,js,lsa,js,2)                           6d10s21
                   if(icase.eq.1)then                                   6d10s21
                    icolk=jg+irefo(js)*(lga+irefo(lsa)*jg)              6d10s21
                   else                                                 6d10s21
                    icolk=lga+irefo(lsa)*(jg+irefo(js)*jg)              6d15s21
                   end if                                               6d10s21
                  else                                                  6d10s21
                   ix=max(jg,lga)                                       6d10s21
                   in=min(jg,lga)                                       6d10s21
                   icolk=((ix*(ix+1))/2)+in+ncolj*jg                    6d10s21
                  end if                                                6d10s21
                  kint=ionext(i2eu)+nvirt(lsa)*icolk                    6d15s21
                  do jv=0,nvirt(lsa)-1                                  6d15s21
                   bc(iintx+jv)=bc(iintx+jv)+2d0*bc(jint+jv)-bc(kint+jv)6d10s21
                  end do                                                6d10s21
                 else                                                   6d10s21
                  do jv=0,nvirt(lsa)-1                                  6d15s21
                   bc(iintx+jv)=bc(iintx+jv)+bc(jint+jv)                6d10s21
                  end do                                                6d10s21
                 end if                                                 6d10s21
                end do                                                  6d10s21
                if(itransgg.eq.0)then                                   8d11s22
                 call ddi_done(ibc(iacc),nacc)                          8d11s22
                 do iz=0,nrow*ncol-1                                    8d11s22
                  bc(igs+iz)=0d0                                        8d11s22
                 end do                                                 8d11s22
                 itransgg=1                                             8d11s22
                end if                                                  8d11s22
                call hcdsuc1(bc(iprod),ncsf(jarg),ncsf(iarg),           6d10s21
     $               ncsf2(1,iarg),bc(iintx),bc(jvs),bc(jgs),bc(jvd),   6d10s21
     $               bc(jgd),lsa,jsbv,isbv12,nvirt,multh,nrootu,nsymb,  6d10s21
     $               sr2,mrow,loop.eq.-74718,loop,loopx,bc,ibc)         12d13s22
                ibcoff=iprod                                            6d11s21
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
                     if(nab2(1).gt.nab2(2))then                          12d18s20
                      icpy=nab2(1)                                       12d18s20
                      nab2(1)=nab2(2)                                    12d18s20
                      nab2(2)=icpy                                       12d18s20
                     end if                                              12d18s20
                     lsa=ism(nab2(1))                                      12d14s20
                     lga=irel(nab2(1))-1                                   12d14s20
                     lsb=ism(nab2(2))                                      12d14s20
                     lgb=irel(nab2(2))-1                                   12d14s20
                     lsc=ism(nab1(1))                                    12d18s20
                     lgc=irel(nab1(1))-1                                 12d18s20
                     lsab=multh(lsa,lsb)                                6d11s21
                     lsd=multh(lsab,lsc)                                6d11s21
                     i2eu=invk1(1,lsa,lsb,lsc,2)                        6d11s21
                     icase=invk1(2,lsa,lsb,lsc,2)                        6d11s21
                     if(lsa.eq.lsb)then                                  12d18s20
                      nn=(irefo(lsa)*(irefo(lsa)+1))/2                   12d18s20
                      icol=((lgb*(lgb+1))/2)+lga+nn*lgc                  12d18s20
                     else                                                12d18s20
                      if(icase.eq.1)then
                       icol=lga+irefo(lsa)*(lgb+irefo(lsb)*lgc)           12d18s20
                      else                                              6d11s21
                       icol=lgb+irefo(lsb)*(lga+irefo(lsa)*lgc)           12d18s20
                      end if                                            6d11s21
                     end if                                              12d18s20
                     call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),      6d11s21
     $                  ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod,11d10s22
     $                    bc,ibc)                                       11d10s22
                     iintx=ionext(i2eu)+nvirt(lsd)*icol                 6d11s21
                     if(itransgg.eq.0)then                                   8d11s22
                      call ddi_done(ibc(iacc),nacc)                          8d11s22
                      do iz=0,nrow*ncol-1                                    8d11s22
                       bc(igs+iz)=0d0                                        8d11s22
                      end do                                                 8d11s22
                      itransgg=1                                             8d11s22
                     end if                                                  8d11s22
                     call hcdsuc1(bc(iprod),ncsf(jarg),ncsf(iarg),           6d10s21
     $               ncsf2(1,iarg),bc(iintx),bc(jvs),bc(jgs),bc(jvd),   6d10s21
     $               bc(jgd),lsd,jsbv,isbv12,nvirt,multh,nrootu,nsymb,  6d10s21
     $               sr2,mrow,loop.eq.-74718,loop,loopx,bc,ibc)         12d13s22
                     ibcoff=iprod                                            6d11s21
                    end if                                              6d11s21
                   end if                                               6d11s21
                  end if                                                6d11s21
                 end if                                                 6d11s21
                end do                                                  6d11s21
               else if(nnot.gt.2)then                                    6d10s21
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
                    if(nab2(1).gt.nab2(2))then                           12d19s20
                     icpy=nab2(1)                                        12d19s20
                     nab2(1)=nab2(2)                                     12d19s20
                     nab2(2)=icpy                                        12d19s20
                    end if                                               12d19s20
                    lsa=ism(nab2(1))                                     12d19s20
                    lga=irel(nab2(1))-1                                  12d19s20
                    lsb=ism(nab2(2))                                     12d19s20
                    lgb=irel(nab2(2))-1                                  12d19s20
                    lsc=ism(nab1(1))                                     12d19s20
                    lgc=irel(nab1(1))-1                                  12d19s20
                   else                                                  12d19s20
                    if(nab1(1).gt.nab1(2))then                           12d19s20
                     icpy=nab1(1)                                        12d19s20
                     nab1(1)=nab1(2)                                     12d19s20
                     nab1(2)=icpy                                        12d19s20
                    end if                                               12d19s20
                    lsa=ism(nab1(1))                                     12d19s20
                    lga=irel(nab1(1))-1                                  12d19s20
                    lsb=ism(nab1(2))                                     12d19s20
                    lgb=irel(nab1(2))-1                                  12d19s20
                    lsc=ism(nab2(1))                                     12d19s20
                    lgc=irel(nab2(1))-1                                  12d19s20
                   end if                                                12d19s20
                   lsab=multh(lsa,lsb)                                  6d11s21
                   lsd=multh(lsab,lsc)                                  6d11s21
                   i2eu=invk1(1,lsa,lsb,lsc,2)                          6d11s21
                   icase=invk1(2,lsa,lsb,lsc,2)                         6d11s21
                   if(lsa.eq.lsb)then                                    12d19s20
                    nn=(irefo(lsa)*(irefo(lsa)+1))/2                     12d19s20
                    icol=((lgb*(lgb+1))/2)+lga+nn*lgc                    12d19s20
                   else                                                  12d19s20
                    if(icase.eq.1)then                                  6d11s21
                     icol=lga+irefo(lsa)*(lgb+irefo(lsb)*lgc)             12d19s20
                    else                                                6d11s21
                     icol=lgb+irefo(lsb)*(lga+irefo(lsa)*lgc)             12d19s20
                    end if                                              6d11s21
                   end if                                                12d19s20
                   call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),        3d19s21
     $                 ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod, 11d10s22
     $                  bc,ibc)                                         11d10s22
                   iintx=ionext(i2eu)+nvirt(lsd)*icol                   6d11s21
                   if(itransgg.eq.0)then                                   8d11s22
                    call ddi_done(ibc(iacc),nacc)                          8d11s22
                    do iz=0,nrow*ncol-1                                    8d11s22
                     bc(igs+iz)=0d0                                        8d11s22
                    end do                                                 8d11s22
                    itransgg=1                                             8d11s22
                   end if                                                  8d11s22
                   call hcdsuc1(bc(iprod),ncsf(jarg),ncsf(iarg),           6d10s21
     $               ncsf2(1,iarg),bc(iintx),bc(jvs),bc(jgs),bc(jvd),   6d10s21
     $               bc(jgd),lsd,jsbv,isbv12,nvirt,multh,nrootu,nsymb,  6d10s21
     $               sr2,mrow,loop.eq.-79366,loop,loopx,bc,ibc)         12d13s22
                   ibcoff=iprod                                          3d19s21
                   if(ipss.eq.2)go to 3                                  12d18s20
                  end if                                                 12d18s20
                 end if                                                 4d20s21
                end do                                                  12d18s20
    3           continue                                                12d18s20
               end if                                                    6d10s21
              end if                                                    6d10s21
              idoit=idoit+1                                             6d10s21
              if(mod(idoit,mynprocg).eq.mynowprog)then                  6d10s21
               call gandc4(j1c,j1o,i2c,k2o,nopen1,nopen2p,norbxxx,nnot,    12d18s20
     $              nab4,bc,ibc)                                        11d14s22
ccc
               if(nnot.gt.1)then                                         12d18s20
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
                    nn=nvirt(jsbv)*nrootu                               6d11s21
                    if(nn.gt.0)then                                     8d19s21
                     call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),       3d19s21
     $               ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod,bc,11d10s22
     $                   ibc)                                           11d10s22
                     lsa=ism(nab1(1))                                       12d19s20
                     lga=irel(nab1(1))-1                                    12d19s20
                     itmpb=ibcoff                                        6d11s21
                     itmpd=itmpb+nn*ncsf(iarg)                           6d11s21
                     ibcoff=itmpd+ncsf(jarg)*ncsf(iarg)                  6d11s21
                     call enough(' hcdsuc.  8',bc,ibc)
                     call dgemm('n','n',nn,ncsf(iarg),ncsf(jarg),1d0,    6d11s21
     $                   bc(jvs),nn,bc(iprod),ncsf(jarg),0d0,bc(itmpb), 6d11s21
     $                   nn,                                            12d13s22
     d'hcdsuca')
                     jvmat=ivmat(jsb)+nn*(ncsf(iarg)*(if2-1)+nfh*lga)    6d11s21
                     do i=0,nn*ncsf(iarg)-1                              6d11s21
                      bc(jvmat+i)=bc(jvmat+i)+bc(itmpb+i)                6d11s21
                     end do                                              6d11s21
                     do i=0,ncsf(iarg)-1                                 6d11s21
                      do j=0,ncsf(jarg)-1                                6d11s21
                       ji=iprod+j+ncsf(jarg)*i                           6d11s21
                       ij=itmpd+i+ncsf(iarg)*j                           6d11s21
                       bc(ij)=bc(ji)                                     6d11s21
                      end do                                             6d11s21
                     end do                                              6d11s21
                     jbmat=ibmat(jsb)+nn*(ncsf(iarg)*(if2-1)+nfh*lga)    6d11s21
                     if(itransgg.eq.0)then                                   8d11s22
                      call ddi_done(ibc(iacc),nacc)                          8d11s22
                      do iz=0,nrow*ncol-1                                    8d11s22
                       bc(igs+iz)=0d0                                        8d11s22
                      end do                                                 8d11s22
                      itransgg=1                                             8d11s22
                     end if                                                  8d11s22
                     call dgemm('n','n',nn,ncsf(jarg),ncsf(iarg),1d0,
     $                   bc(jbmat),nn,bc(itmpd),ncsf(iarg),1d0,bc(jgs), 6d11s21
     $                   nn
     d,'hcdsucb'                    )                                            6d11s21
                     ibcoff=iprod                                         3d19s21
                    end if                                              8d19s21
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
            if(itransgg.ne.0)then                                       8d11s22
             do icol=0,ncol-1                                            6d10s21
              jgs=igs+icol*nrow                                          6d10s21
              jvs=ivs+icol*nrow                                          6d10s21
              do jv=0,nvirt(jsbv)-1                                      6d10s21
               do ir=0,nrootm                                            6d10s21
                ij=jvs+ir+nrootu*jv                                      6d10s21
                ji=jgs+jv+nvirt(jsbv)*ir                                 6d10s21
                bc(ij)=bc(ji)                                            6d10s21
               end do                                                    6d10s21
              end do                                                     6d10s21
             end do                                                      6d10s21
             i28=nrow                                                    6d10s21
             i38=ncol                                                    6d10s21
             do icol=0,nrow*ncol-1                                       8d11s22
              bc(igs+icol)=bc(ivs+icol)                                  8d11s22
             end do                                                      8d11s22
             if(naccd.gt.0)call ddi_done(ibc(iaccd),naccd)              8d22s22
             call ddi_iacc(bc,ibc,ihsdiag(nclo1p,jsb,2),i18,i28,i18,i38,11d15s22
     $            bc(igs),ibc(iacc),nacc)                               8d11s22
             nwiacc=nwiacc+nrow*ncol                                     8d11s22
             itransgg=0                                                  8d11s22
            end if                                                      8d11s22
            ibcoff=ivs                                                  6d10s21
           end if                                                       6d10s21
          end do                                                        6d10s21
         end do                                                         6d10s21
         call dws_gsumf(bc(ibcvmat),nbmat)                                3d2s21
c
c     we need separate copy of igd here, for it will come out with
c     all singlet pairs before triplet pairs.
c
         igdst(1)=ibcoff                                                6d11s21
         igdst(2)=igdst(1)+nrootu*nfdat(1)*(nvisv+nvnotv)               6d11s21
         ibcoff=igdst(2)+nrootu*nfdat(2)*nvnotv                         6d11s21
         call enough(' hcdsuc.  9',bc,ibc)
         do i=igdst(1),ibcoff-1                                         6d11s21
          bc(i)=0d0                                                     6d11s21
         end do                                                         6d11s21
         do jsb=1,nsymb                                                   2d4s21
          jsbv=multh(isymmrci,jsb)                                      6d17s21
          isn=multh(isbv12,jsbv)                                          2d4s21
          if(min(irefo(isn),nvirt(jsbv)).gt.0)then                        2d25s21
c
c     let us re-order vmat so that singlet pairs occur before triplet   6d11s21
c     pairs                                                             6d11s21
c
           nnn=nvirt(jsbv)*nrootu                                       6d11s21
           do l=1,2                                                     6d11s21
            nn=nnn*nfdat(l)                                             6d11s21
            ivtmp(l)=ibcoff                                             6d11s21
            ibcoff=ivtmp(l)+nn*irefo(isn)                               6d11s21
            jvtmp(l)=ivtmp(l)                                           6d11s21
           end do                                                       6d11s21
           call enough(' hcdsuc. 10',bc,ibc)
           jvmat=ivmat(jsb)                                             6d11s21
           do in=0,irefo(isn)-1                                         6d11s21
            do if2=0,nff2(nclo2p,isb)-1                                 6d11s21
             do i=0,ncsf2(1,iarg)-1                                     6d11s21
              do j=0,nnn-1                                              6d11s21
               bc(jvtmp(1)+j)=bc(jvmat+j)                               6d11s21
              end do                                                    6d11s21
              jvtmp(1)=jvtmp(1)+nnn                                     6d11s21
              jvmat=jvmat+nnn                                           6d11s21
             end do                                                     6d11s21
             do i=ncsf2(1,iarg),ncsf(iarg)-1                            6d11s21
              do j=0,nnn-1                                              6d11s21
               bc(jvtmp(2)+j)=bc(jvmat+j)                               6d11s21
              end do                                                    6d11s21
              jvtmp(2)=jvtmp(2)+nnn                                     6d16s21
              jvmat=jvmat+nnn                                           6d11s21
             end do                                                     6d11s21
            end do                                                      6d11s21
           end do                                                       6d11s21
c
c     form Gv'v"rk=[(vv"|nv')+p(k)(vv'|nv")]vvrkn
c
           jgdst(1)=igdst(1)                                            6d11s21
           jgdst(2)=igdst(2)                                            6d11s21
           do isbv1=1,nsymb                                               2d4s21
            isbv2=multh(isbv1,isbv12)                                     2d4s21
            if(isbv2.ge.isbv1)then                                        2d4s21
             call ilimts(irefo(isn),nvirt(isbv1),mynprocg,mynowprog,il,
     $           ih,i1s,i1e,i2s,i2e)                                         2d4s21
             i2eu=invk1(1,jsbv,isbv2,isn,2)                               2d4s21
             icase=invk1(2,jsbv,isbv2,isn,2)                              2d4s21
             if(isbv1.eq.isbv2)then                                       2d4s21
              i10=i1s                                                     2d4s21
              i1n=irefo(isn)                                              2d4s21
              iint=i3x(i2eu)                                              2d4s21
              if(jsbv.eq.isbv2)then                                       2d4s21
               nrow=(nvirt(jsbv)*(nvirt(jsbv)+1))/2                       2d4s21
               do i2=i2s,i2e                                              2d4s21
                iv1=i2-1                                                  2d4s21
                if(i2.eq.i2e)i1n=i1e                                      2d4s21
                do i1=i10,i1n                                             2d4s21
                 i1m=i1-1                                                 2d4s21
                 do k=0,nfdat(1)-1                                        6d11s21
                  do ir=0,nrootm                                          2d4s21
                   sum=0d0                                                2d4s21
                   jvmat=ivtmp(1)+nvirt(jsbv)*(ir+nrootu*(k             6d11s21
     $                  +nfdat(1)*i1m))                                 6d11s21
                   iadv=jgdst(1)+iv1+nvirt(isbv1)*(ir+nrootu*k)         6d11s21
                   do jv=0,nvirt(jsbv)-1                                    2d4s21
                    ix=max(jv,iv1)                                          2d4s21
                    in=min(jv,iv1)                                          2d4s21
                    irow=iint+((ix*(ix+1))/2)+in                            2d4s21
                    sum=sum+bc(jvmat+jv)*bc(irow)                         2d4s21
                   end do                                                 2d4s21
                   sum=sum*sr2                                            2d4s21
                   bc(iadv)=bc(iadv)+sum                                6d11s21
                  end do                                                  2d4s21
                 end do                                                   2d4s21
                 iint=iint+nrow                                           2d4s21
                end do                                                    2d4s21
                i10=1                                                     2d4s21
               end do                                                     2d4s21
              else if(icase.eq.1)then                                     2d4s21
               nrow=nvirt(jsbv)*nvirt(isbv1)                              2d4s21
               do i2=i2s,i2e                                              2d4s21
                iv1=i2-1                                                  2d4s21
                if(i2.eq.i2e)i1n=i1e                                      2d4s21
                do i1=i10,i1n                                             2d4s21
                 i1m=i1-1                                                 2d4s21
                 do k=0,nfdat(1)-1                                      6d11s21
                  do ir=0,nrootm                                          2d4s21
                   sum=0d0                                                2d4s21
                   jvmat=ivtmp(1)+nvirt(jsbv)*(ir+nrootu*(k             6d11s21
     $                  +nfdat(1)*i1m))                                 6d11s21
                   irow=iint+iv1*nvirt(jsbv)                              2d4s21
                   iadv=jgdst(1)+iv1+nvirt(isbv1)*(ir+nrootu*k)         6d11s21
                   do jv=0,nvirt(jsbv)-1                                    2d4s21
                    sum=sum+bc(jvmat+jv)*bc(irow+jv)                      2d4s21
                   end do                                                 2d4s21
                   sum=sum*sr2                                            2d4s21
                   bc(iadv)=bc(iadv)+sum                                  2d4s21
                  end do                                                  2d4s21
                 end do                                                   2d4s21
                 iint=iint+nrow                                           2d4s21
                end do                                                    2d4s21
                i10=1                                                     2d4s21
               end do                                                     2d4s21
              else                                                        2d4s21
               nrow=nvirt(jsbv)*nvirt(isbv1)                              2d4s21
               do i2=i2s,i2e                                              2d4s21
                iv1=i2-1                                                  2d4s21
                if(i2.eq.i2e)i1n=i1e                                      2d4s21
                do i1=i10,i1n                                             2d4s21
                 i1m=i1-1                                                 2d4s21
                 do k=0,nfdat(1)-1                                      6d11s21
                  do ir=0,nrootm                                          2d4s21
                   sum=0d0                                                2d4s21
                   jvmat=ivtmp(1)+nvirt(jsbv)*(ir+nrootu*(k             6d11s21
     $                  +nfdat(1)*i1m))                                 6d11s21
                   irow=iint+iv1                                          2d4s21
                   iadv=jgdst(1)+iv1+nvirt(isbv1)*(ir+nrootu*k)         6d11s21
                   do jv=0,nvirt(jsbv)-1                                    2d4s21
                    sum=sum+bc(jvmat+jv)*bc(irow+jv*nvirt(isbv1))         2d4s21
                   end do                                                 2d4s21
                   sum=sum*sr2                                            2d4s21
                   bc(iadv)=bc(iadv)+sum                                6d11s21
                  end do                                                  2d4s21
                 end do                                                   2d4s21
                 iint=iint+nrow                                           2d4s21
                end do                                                    2d4s21
                i10=1                                                     2d4s21
               end do                                                     2d4s21
              end if                                                      2d4s21
              jgdst(1)=jgdst(1)+nvirt(isbv1)*nrootu*nfdat(1)             6d11s21
              nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       2d4s21
              isw=0                                                       2d4s21
             else                                                         2d4s21
              nvv=nvirt(isbv1)*nvirt(isbv2)                               2d4s21
              isw=1                                                       2d4s21
             end if                                                       2d4s21
             i10=i1s                                                      2d4s21
             i1n=irefo(isn)                                               2d4s21
             iint=i3x(i2eu)                                               2d4s21
             if(jsbv.eq.isbv2)then                                        2d4s21
              nrow=(nvirt(jsbv)*(nvirt(jsbv)+1))/2                        2d4s21
              do i2=i2s,i2e                                               2d4s21
               iv1=i2-1                                                   2d4s21
               ibots=i2                                                   2d4s21
               ibotn=0                                                    2d4s21
               ibot=ibots+isw*(ibotn-ibots)                               2d4s21
               if(i2.eq.i2e)i1n=i1e                                       2d4s21
               do i1=i10,i1n                                              2d4s21
                i1m=i1-1                                                  2d4s21
                do l=1,2                                                6d11s21
                 nn=nnn*nfdat(l)                                        6d16s21
                 do kr=0,nfdat(l)*nrootu-1                              6d11s21
                  do iv2=ibot,nvirt(isbv2)-1                                2d4s21
                   jvmat=ivtmp(l)+nvirt(jsbv)*kr+nn*i1m                 6d16s21
                   sum=0d0                                                 2d4s21
                   itri=((iv2*(iv2-1))/2)+iv1                              2d4s21
                   irec=iv1+nvirt(isbv1)*iv2                               2d4s21
                   iadv=jgdst(l)+itri+isw*(irec-itri)+nvv*kr            6d11s21
                   do jv=0,nvirt(jsbv)-1                                    2d4s21
                    ix=max(jv,iv2)                                         2d4s21
                    in=min(jv,iv2)                                         2d4s21
                    irow=iint+((ix*(ix+1))/2)+in                           2d4s21
                    orig=sum
                    sum=sum+bc(jvmat+jv)*bc(irow)                          2d4s21
                   end do                                                  2d4s21
                   bc(iadv)=bc(iadv)+sum                                   2d4s21
                  end do                                                   2d4s21
                 end do                                                    2d4s21
                end do                                                  6d11s21
                iint=iint+nrow                                            2d4s21
               end do                                                     2d4s21
               i10=1                                                      2d4s21
              end do                                                      2d4s21
             else if(icase.eq.1)then                                      2d4s21
              nrow=nvirt(jsbv)*nvirt(isbv2)                               2d4s21
              do i2=i2s,i2e                                               2d4s21
               iv1=i2-1                                                   2d4s21
               ibots=i2                                                   2d4s21
               ibotn=0                                                    2d4s21
               ibot=ibots+isw*(ibotn-ibots)                               2d4s21
               if(i2.eq.i2e)i1n=i1e                                       2d4s21
               do i1=i10,i1n                                              2d4s21
                i1m=i1-1                                                  2d4s21
                do l=1,2                                                6d11s21
                 nn=nnn*nfdat(l)                                        6d16s21
                 do kr=0,nfdat(l)*nrootu-1                              6d11s21
                  do iv2=ibot,nvirt(isbv2)-1                                2d4s21
                   jvmat=ivtmp(l)+nvirt(jsbv)*kr+nn*i1m                 6d16s21
                   sum=0d0                                                 2d4s21
                   irow=iint+nvirt(jsbv)*iv2                               2d4s21
                   itri=((iv2*(iv2-1))/2)+iv1                              2d4s21
                   irec=iv1+nvirt(isbv1)*iv2                               2d4s21
                   iadv=jgdst(l)+itri+isw*(irec-itri)+nvv*kr            6d11s21
                   do jv=0,nvirt(jsbv)-1                                    2d4s21
                    sum=sum+bc(jvmat+jv)*bc(irow+jv)                       2d4s21
                   end do                                                  2d4s21
                   bc(iadv)=bc(iadv)+sum                                   2d4s21
                  end do                                                   2d4s21
                 end do                                                    2d4s21
                end do                                                  6d11s21
                iint=iint+nrow                                            2d4s21
               end do                                                     2d4s21
               i10=1                                                      2d4s21
              end do                                                      2d4s21
             else                                                         2d4s21
              nrow=nvirt(jsbv)*nvirt(isbv2)                               2d4s21
              do i2=i2s,i2e                                               2d4s21
               iv1=i2-1                                                   2d4s21
               ibots=i2                                                   2d4s21
               ibotn=0                                                    2d4s21
               ibot=ibots+isw*(ibotn-ibots)                               2d4s21
               if(i2.eq.i2e)i1n=i1e                                       2d4s21
               do i1=i10,i1n                                              2d4s21
                i1m=i1-1                                                  2d4s21
                do l=1,2                                                6d11s21
                 nn=nnn*nfdat(l)                                        6d16s21
                 do kr=0,nfdat(l)*nrootu-1                              6d11s21
                  jvmat=ivtmp(l)+nvirt(jsbv)*kr+nn*i1m                  6d11s21
                  do jv=0,nvirt(jsbv)-1                                    2d4s21
                   irow=iint+nvirt(isbv2)*jv                               2d4s21
                   do iv2=ibot,nvirt(isbv2)-1                                2d4s21
                    itri=((iv2*(iv2-1))/2)+iv1                              2d4s21
                    irec=iv1+nvirt(isbv1)*iv2                               2d4s21
                    iadv=jgdst(l)+itri+isw*(irec-itri)+nvv*kr           6d11s21
                    bc(iadv)=bc(iadv)+bc(jvmat+jv)*bc(irow+iv2)            2d4s21
                   end do                                                  2d4s21
                  end do                                                   2d4s21
                 end do                                                    2d4s21
                end do                                                  6d11s21
                iint=iint+nrow                                            2d4s21
               end do                                                     2d4s21
               i10=1                                                      2d4s21
              end do                                                      2d4s21
             end if                                                       2d4s21
             call ilimts(irefo(isn),nvirt(isbv2),mynprocg,mynowprog,il,
     $            ih,i1s,i1e,i2s,i2e)                                         2d4s21
             i2eu=invk1(1,jsbv,isbv1,isn,2)                               2d4s21
             icase=invk1(2,jsbv,isbv1,isn,2)                              2d4s21
             i10=i1s                                                      2d4s21
             i1n=irefo(isn)                                               2d4s21
             iint=i3x(i2eu)                                               2d4s21
             if(jsbv.eq.isbv1)then                                        2d4s21
              nrow=(nvirt(jsbv)*(nvirt(jsbv)+1))/2                        2d4s21
              do i2=i2s,i2e                                               2d4s21
               iv2=i2-1                                                   2d4s21
               itops=iv2-1                                                2d4s21
               itopn=nvirt(isbv1)-1                                       2d4s21
               itop=itops+isw*(itopn-itops)                               2d4s21
               if(i2.eq.i2e)i1n=i1e                                       2d4s21
               do i1=i10,i1n                                              2d4s21
                i1m=i1-1                                                  2d4s21
                psr=1d0                                                 6d11s21
                do l=1,2                                                6d11s21
                 nn=nnn*nfdat(l)                                        6d16s21
                 do kr=0,nfdat(l)*nrootu-1                              6d11s21
                  do iv1=0,itop                                            2d4s21
                   jvmat=ivtmp(l)+nvirt(jsbv)*kr+nn*i1m                 6d16s21
                   sum=0d0                                                 2d4s21
                   itri=((iv2*(iv2-1))/2)+iv1                              2d4s21
                   irec=iv1+nvirt(isbv1)*iv2                               2d4s21
                   iadv=jgdst(l)+itri+isw*(irec-itri)+nvv*kr            6d11s21
                   do jv=0,nvirt(jsbv)-1                                    2d4s21
                    ix=max(jv,iv1)                                         2d4s21
                    in=min(jv,iv1)                                         2d4s21
                    irow=iint+((ix*(ix+1))/2)+in                           2d4s21
                    orig=sum
                    sum=sum+bc(jvmat+jv)*bc(irow)                          2d4s21
                   end do                                                  2d4s21
                   bc(iadv)=bc(iadv)+sum*psr                            6d11s21
                  end do                                                   2d4s21
                 end do                                                    2d4s21
                 psr=-1d0                                               6d11s21
                end do                                                  6d11s21
                iint=iint+nrow                                            2d4s21
               end do                                                     2d4s21
               i10=1                                                      2d4s21
              end do                                                      2d4s21
             else if(icase.eq.1)then                                      2d4s21
              nrow=nvirt(jsbv)*nvirt(isbv1)                               2d4s21
              do i2=i2s,i2e                                               2d4s21
               iv2=i2-1                                                   2d4s21
               itops=iv2-1                                                2d4s21
               itopn=nvirt(isbv1)-1                                       2d4s21
               itop=itops+isw*(itopn-itops)                               2d4s21
               if(i2.eq.i2e)i1n=i1e                                       2d4s21
               do i1=i10,i1n                                              2d4s21
                i1m=i1-1                                                  2d4s21
                psr=1d0                                                 6d11s21
                do l=1,2                                                6d11s21
                 nn=nnn*nfdat(l)                                        6d16s21
                 do kr=0,nfdat(l)*nrootu-1                              6d11s21
                  do iv1=0,itop                                            2d4s21
                   jvmat=ivtmp(l)+nvirt(jsbv)*kr+nn*i1m                 6d16s21
                   sum=0d0                                                 2d4s21
                   itri=((iv2*(iv2-1))/2)+iv1                              2d4s21
                   irec=iv1+nvirt(isbv1)*iv2                               2d4s21
                   iadv=jgdst(l)+itri+isw*(irec-itri)+nvv*kr            6d11s21
                   irow=iint+nvirt(jsbv)*iv1                               2d4s21
                   do jv=0,nvirt(jsbv)-1                                    2d4s21
                    sum=sum+bc(jvmat+jv)*bc(irow+jv)                       2d4s21
                   end do                                                  2d4s21
                   bc(iadv)=bc(iadv)+sum*psr                            6d11s21
                  end do                                                   2d4s21
                 end do                                                    2d4s21
                 psr=-1d0                                               6d11s21
                end do                                                  6d11s21
                iint=iint+nrow                                            2d4s21
               end do                                                     2d4s21
               i10=1                                                      2d4s21
              end do                                                      2d4s21
             else                                                         2d4s21
              nrow=nvirt(jsbv)*nvirt(isbv1)                               2d4s21
              do i2=i2s,i2e                                               2d4s21
               iv2=i2-1                                                   2d4s21
               itops=iv2-1                                                2d4s21
               itopn=nvirt(isbv1)-1                                       2d4s21
               itop=itops+isw*(itopn-itops)                               2d4s21
               if(i2.eq.i2e)i1n=i1e                                       2d4s21
               do i1=i10,i1n                                              2d4s21
                i1m=i1-1                                                  2d4s21
                psr=1d0                                                 6d16s21
                do l=1,2                                                6d11s21
                 nn=nnn*nfdat(l)                                        6d16s21
                 do kr=0,nfdat(l)*nrootu-1                              6d11s21
                  jvmat=ivtmp(l)+nvirt(jsbv)*kr+nn*i1m                  6d11s21
                  do jv=0,nvirt(jsbv)-1                                    2d4s21
                   irow=iint+nvirt(isbv1)*jv                               2d4s21
                   fact=bc(jvmat+jv)*psr                                6d11s21
                   do iv1=0,itop                                           2d4s21
                    itri=((iv2*(iv2-1))/2)+iv1                              2d4s21
                    irec=iv1+nvirt(isbv1)*iv2                               2d4s21
                    iadv=jgdst(l)+itri+isw*(irec-itri)+nvv*kr           6d11s21
                    bc(iadv)=bc(iadv)+fact*bc(irow+iv1)                 6d11s21
                   end do                                                  2d4s21
                  end do                                                   2d4s21
                 end do                                                    2d4s21
                 psr=-1d0                                               6d11s21
                end do                                                  6d11s21
                iint=iint+nrow                                            2d4s21
               end do                                                     2d4s21
               i10=1                                                      2d4s21
              end do                                                      2d4s21
             end if                                                       2d4s21
             jgdst(1)=jgdst(1)+nvv*nrootu*nfdat(1)                      6d11s21
             jgdst(2)=jgdst(2)+nvv*nrootu*nfdat(2)                      6d11s21
            end if                                                        2d4s21
           end do                                                         2d4s21
           ibcoff=ivtmp(1)                                              6d11s21
          end if                                                        6d11s21
         end do                                                         6d11s21
         jgdst(1)=igdst(1)                                              6d11s21
         jgdst(2)=igdst(2)                                              6d11s21
c     (nvv,r,ncsf,icol)>(nvv,r,ncsf),icol
         jgd=igd                                                        6d16s21
         do isbv1=1,nsymb                                               6d16s21
          isbv2=multh(isbv1,isbv12)                                     6d16s21
          if(isbv2.ge.isbv1)then                                        6d16s21
           if(isbv1.eq.isbv2)then                                       6d16s21
            nn=ncsf2(1,iarg)*nrootu*nvirt(isbv1)                        6d16s21
            do iff=0,nff2(nclo2p,isb)-1                                 6d16s21
             lgd=jgd+mrow*iff                                           6d16s21
             do irv=0,nn-1                                              6d16s21
              orig=bc(lgd+irv)
              bc(lgd+irv)=bc(lgd+irv)+bc(jgdst(1)+irv)                  6d16s21
             end do                                                     6d16s21
             jgdst(1)=jgdst(1)+nn                                       6d16s21
            end do                                                      6d16s21
            jgd=jgd+nn                                                  6d16s21
            nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       6d16s21
           else                                                         6d16s21
            nvv=nvirt(isbv1)*nvirt(isbv2)                               6d16s21
           end if                                                       6d16s21
           nvvu=nvv*nrootu                                              6d16s21
           do iff=0,nff2(nclo2p,isb)-1                                  6d16s21
            lgd=jgd+mrow*iff                                            6d16s21
            do irv=0,mfdat(1)*nvvu-1                                    6d16s21
             orig=bc(lgd+irv)
             bc(lgd+irv)=bc(lgd+irv)+bc(jgdst(1)+irv)                   6d16s21
            end do                                                      6d16s21
            lgd=lgd+mfdat(1)*nvvu                                       6d16s21
            jgdst(1)=jgdst(1)+mfdat(1)*nvvu                             6d16s21
            do irv=0,mfdat(2)*nvvu-1                                    6d16s21
             orig=bc(lgd+irv)
             bc(lgd+irv)=bc(lgd+irv)+bc(jgdst(2)+irv)                   6d16s21
            end do                                                      6d16s21
            lgd=lgd+mfdat(2)*nvvu                                       6d16s21
            jgdst(2)=jgdst(2)+mfdat(2)*nvvu                             6d16s21
           end do                                                       6d16s21
           jgd=jgd+nvvu*ncsf(iarg)                                      6d16s21
          end if                                                        6d16s21
         end do                                                         6d16s21
         i28=mrow                                                       6d10s21
         i38=nff2(nclo2p,isb)                                           6d10s21
c     (nvv,r,ncsf),icol
         do icol=0,nff2(nclo2p,isb)-1                                   6d10s21
          jgd=igd+mrow*icol                                             6d10s21
          jvd=ivd+mrow*icol                                             6d10s21
          do isbv1=1,nsymb                                              6d10s21
           isbv2=multh(isbv1,isbv12)                                    6d10s21
           if(isbv2.ge.isbv1)then                                       6d10s21
            if(isbv1.eq.isbv2)then                                      6d10s21
             do i=0,ncsf2(1,iarg)-1                                     6d10s21
              do iv=0,nvirt(isbv1)-1                                    6d10s21
               do j=0,nrootm                                            6d10s21
                jvi=jvd+j+nrootu*(i+ncsf2(1,iarg)*iv)                   6d15s21
                ivj=jgd+iv+nvirt(isbv1)*(j+nrootu*i)                    6d10s21
                bc(jvi)=bc(ivj)                                         6d10s21
               end do                                                   6d10s21
              end do                                                    6d10s21
             end do                                                     6d10s21
             nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                      6d10s21
             nhere=nrootu*nvirt(isbv1)*ncsf2(1,iarg)                    6d10s21
             jgd=jgd+nhere                                              6d10s21
             jvd=jvd+nhere                                              6d10s21
            else                                                        6d10s21
             nvv=nvirt(isbv1)*nvirt(isbv2)                              6d10s21
            end if                                                      6d10s21
            do i=0,ncsf(iarg)-1                                         6d10s21
             do ivv=0,nvv-1                                             6d10s21
              do j=0,nrootm                                             6d10s21
               jvi=jvd+j+nrootu*(i+ncsf(iarg)*ivv)                      6d15s21
               ivj=jgd+ivv+nvv*(j+nrootu*i)                             6d10s21
               bc(jvi)=bc(ivj)                                          6d10s21
              end do                                                    6d10s21
             end do                                                     6d10s21
            end do                                                      6d10s21
            nhere=nrootu*nvv*ncsf(iarg)                                  6d10s21
            jgd=jgd+nhere                                               6d10s21
            jvd=jvd+nhere                                               6d10s21
           end if                                                       6d10s21
          end do                                                        6d10s21
         end do                                                         6d10s21
         do i=0,mrow*nff2(nclo2p,isb)-1                                 8d11s22
          bc(igd+i)=bc(ivd+i)                                           8d11s22
         end do                                                         8d11s22
c
c     it seems if the message is too large, it can get lost ...         8d22s22
c
         if(nacc.gt.0)call ddi_done(ibc(iacc),nacc)                     8d22s22
         call ddi_iacc(bc,ibc,ihddiag(nclo2p,isb,2),i18,i28,i18,i38,    11d15s22
     $        bc(igd),ibc(iaccd),naccd)                                 11d15s22
         nwiacc=nwiacc+mrow*nff2(nclo2p,isb)                            8d11s22
        end if                                                          6d11s21
       end do                                                           6d11s21
      end do                                                            12d18s20
      if(nacc.gt.0)call ddi_done(ibc(iacc),nacc)                        8d31s22
      if(naccd.gt.0)call ddi_done(ibc(iaccd),naccd)                     8d31s22
      ibcoff=ircv                                                       1d30s21
      if(ldebug)then                                                    6d13s21
       write(6,*)('all done in hcdsuc')
       call dws_synca                                                    1d24s21
       nsing=0                                                          6d14s21
       ndoub=0                                                          6d14s21
       do isb=1,nsymb                                                   6d14s21
        isbv12=multh(isb,isymmrci)                                      6d14s21
        nvisv=0                                                         6d14s21
        nvnotv=0                                                        6d14s21
        do isbv1=1,nsymb                                                6d14s21
         isbv2=multh(isbv1,isbv12)                                      6d14s21
         if(isbv2.ge.isbv1)then                                         6d14s21
          if(isbv2.eq.isbv1)then                                        6d14s21
           nvisv=nvisv+nvirt(isbv1)                                     6d14s21
           nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        6d14s21
          else                                                          6d14s21
           nvv=nvirt(isbv1)*nvirt(isbv2)                                6d14s21
          end if                                                        6d14s21
          nvnotv=nvnotv+nvv                                             6d14s21
         end if                                                         6d14s21
        end do                                                          6d14s21
        do ii=mdon+1,mdoo+1                                             6d14s21
         if(nff1(ii,isb,1).gt.0)then                                      6d14s21
          iarg=ii-mdon                                                  6d14s21
          nsing=nsing+ncsf(iarg)*nvirt(isbv12)*nff1(ii,isb,1)             6d14s21
         end if                                                         6d14s21
         if(nff2(ii,isb).gt.0)then                                      6d14s21
          iarg=ii-mdon                                                  6d14s21
          ndoub=ndoub+(ncsf2(1,iarg)*nvisv+ncsf(iarg)*nvnotv)           6d14s21
     $         *nff2(ii,isb)                                            6d14s21
         end if                                                         6d14s21
        end do                                                          6d14s21
       end do                                                           6d14s21
       write(6,*)('what we have for hcds: '),ndoub
       igdmaster=ibcoff                                                 6d15s21
       ibcoff=igdmaster+ndoub*nrootu
       jgdmaster=igdmaster
       call enough('hcdsuc.11',bc,ibc)
       do isb=1,nsymb                                                    12d7s20
        isbv12=multh(isb,isymmrci)                                         12d7s20
        do l=1,4                                                        6d11s21
         mtmpx(l)=0                                                     6d11s21
        end do                                                          6d11s21
        do ii=mdon+1,mdoo+1                                             6d11s21
         if(nff2(ii,isb).gt.0)then                                      6d11s21
          iarg=ii-mdon                                                  6d11s21
          do l=1,4                                                      6d11s21
           mtmpx(l)=mtmpx(l)+ncsf2(l,iarg)*nff2(ii,isb)                 6d11s21
          end do                                                        6d11s21
         end if                                                         6d11s21
        end do                                                          6d11s21
        nvisv=0                                                         6d8s21
        nvnotv=0                                                        6d8s21
        ibcsrt=ibcoff                                                   6d11s21
        nvisv=0                                                         6d12s21
        do isbv1=1,nsymb                                                6d8s21
         isbv2=multh(isbv1,isbv12)                                      6d8s21
         if(isbv2.ge.isbv1)then                                         6d8s21
          nvisv0=0                                                      6d11s21
          if(isbv1.eq.isbv2)then                                        6d8s21
           nvisv=nvisv+nvirt(isbv1)                                     6d8s21
           nvisv0=nvirt(isbv1)                                          6d11s21
           nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        6d9s21
          else                                                          6d8s21
           nvv=nvirt(isbv1)*nvirt(isbv2)                                6d8s21
          end if                                                        6d8s21
          do l=1,4                                                      6d11s21
           if(l.eq.1)then                                               6d11s21
            nuse=nvv+nvisv0                                             6d11s21
           else                                                         6d11s21
            nuse=nvv                                                    6d11s21
           end if                                                       6d11s21
           ntmpx(l,isbv1)=nuse*mtmpx(l)                                 6d11s21
           itmpx(l,isbv1)=jgdmaster                                     6d15s21
           jtmpx(l,isbv1)=jgdmaster                                     6d15s21
           jgdmaster=itmpx(l,isbv1)+nrootu*ntmpx(l,isbv1)                   6d11s21
          end do                                                        6d11s21
          nvnotv=nvnotv+nvv                                             6d8s21
         end if                                                         6d8s21
        end do                                                          6d8s21
        call enough(' hcdsuc. 12',bc,ibc)
        xnan=-2d0                                                       6d11s21
        do i=ibcsrt,ibcoff-1                                            6d11s21
         bc(i)=xnan                                                     6d11s21
        end do                                                          6d11s21
        do ii=1,mdoo+1                                                   12d7s20
         if(nff2(ii,isb).gt.0)then                                      6d8s21
          iarg=ii-mdon                                                    12d7s20
          nrow=nrootu*(nvisv*ncsf2(1,iarg)+nvnotv*ncsf(iarg))
          itmp=ibcoff                                                   6d8s21
          ibcoff=itmp+nrow*nff2(ii,isb)                                 6d8s21
          call enough(' hcdsuc. 13',bc,ibc)
          do i=itmp,ibcoff-1
           bc(i)=xnan
          end do
          i2o=1                                                         6d8s21
          i2c=nrow                                                      6d8s21
          j2o=1                                                         6d8s21
          k2o=nff2(ii,isb)                                              6d8s21
          call ddi_get(bc,ibc,ihddiag(ii,isb,2),i2o,i2c,j2o,k2o,        11d15s22
     $         bc(itmp))                                                11d15s22
          jtmp=itmp                                                     6d11s21
          do if2=0,nff2(ii,isb)-1                                       6d11s21
           do isbv1=1,nsymb                                              6d11s21
            isbv2=multh(isbv1,isbv12)                                    6d11s21
            if(isbv2.ge.isbv1)then                                       6d11s21
             if(isbv1.eq.isbv2)then                                      6d11s21
              nvvs=(nvirt(isbv1)*(nvirt(isbv1)+1))/2                    6d11s21
              nvvt=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                    6d11s21
              do iv=0,nvirt(isbv1)-1                                    6d11s21
               ivv=((iv*(iv+1))/2)+iv                                   6d11s21
               do i=0,ncsf2(1,iarg)-1                                    6d11s21
                do ir=0,nrootu-1                                         6d11s21
                 iad=jtmpx(1,isbv1)+i+mtmpx(1)*(ivv+nvvs*ir)            6d11s21
                 bc(iad)=bc(jtmp+ir)                                    6d11s21
                end do                                                  6d11s21
                jtmp=jtmp+nrootu                                         6d11s21
               end do                                                   6d11s21
              end do                                                    6d11s21
              isw=0                                                     6d11s21
             else                                                       6d11s21
              nvvs=nvirt(isbv1)*nvirt(isbv2)                            6d11s21
              nvvt=nvvs                                                 6d11s21
              isw=1                                                     6d11s21
             end if                                                     6d11s21
             do iv2=0,nvirt(isbv2)-1                                    6d11s21
              itop=iv2-1                                                6d11s21
              itop=itop+isw*(nvirt(isbv1)-1-itop)                       6d11s21
              do iv1=0,itop                                             6d11s21
               nvv=nvvs                                                   6d11s21
               ip=+1                                                      6d11s21
               do l=1,4                                                   6d11s21
                if(ntmpx(l,isbv1).gt.0)then                               6d11s21
                itri=((iv2*(iv2+ip))/2)+iv1                              6d11s21
                irec=iv1+nvirt(isbv1)*iv2                                6d11s21
                ivv=itri+isw*(irec-itri)                                 6d11s21
                 do i=0,ncsf2(l,iarg)-1                                   6d11s21
                  do ir=0,nrootu-1                                       6d11s21
                   iad=jtmpx(l,isbv1)+i+mtmpx(l)*(ivv+nvv*ir)           6d11s21
                   bc(iad)=bc(jtmp+ir)                                  6d11s21
                  end do                                                6d11s21
                  jtmp=jtmp+nrootu                                       6d11s21
                 end do                                                 6d11s21
                end if                                                   6d11s21
                nvv=nvvt                                                  6d11s21
                ip=-1                                                     6d11s21
               end do                                                   6d11s21
              end do                                                    6d11s21
             end do                                                     6d11s21
             do l=1,4
              jtmpx(l,isbv1)=jtmpx(l,isbv1)+ncsf2(l,iarg)               6d11s21
             end do
            end if                                                       6d11s21
           end do                                                       6d11s21
          end do                                                        6d11s21
          ibcoff=itmp                                                   6d8s21
         end if                                                         6d8s21
        end do
       end do                                                            12d7s20
       write(6,*)('the whole kit and kaboodle')
       call prntm2(bc(igdmaster),ndoub,nrootu,ndoub)
       write(6,*)('what we have for gs ')
       igg=ibcoff                                                        12d23s20
       ibcoff=igg+nsing*nrootu
       call enough('hcdsuc.14',bc,ibc)
       jgg=igg
       do isb=1,nsymb                                                    12d22s20
        isbv=multh(isb,isymmrci)                                         12d22s20
        do nclo1=mdon,mdoo                                               12d22s20
         nclo1p=nclo1+1                                                  12d22s20
         if(min(nff1(nclo1p,isb,1),nvirt(isbv)).gt.0)then                12d22s20
          jarg=nclo1p-mdon
          ngg=nvirt(isbv)*nrootu                                         12d15s20
          ncol=nff1(nclo1p,isb,1)*ncsf(jarg)                              12d18s20
          igs=ibcoff                                                     12d22s20
          ibcoff=igs+ngg*ncol                                            12d22s20
          call enough(' hcdsuc. 15',bc,ibc)
          i18=1
          i28=ngg                                                        12d18s20
          i38=1                                                          12d18s20
          i48=ncol
          call ddi_get(bc,ibc,ihsdiag(nclo1p,isb,2),i18,i28,i38,i48,    11d15s22
     $         bc(igs))                                                 11d15s22
          jgs=igs
          do i=0,nff1(nclo1p,isb,1)*ncsf(jarg)-1
           do iv=0,nvirt(isbv)-1                                         12d23s20
            do ir=0,nrootu-1
             bc(jgg+ir*nsing)=bc(jgs)
             jgs=jgs+1
            end do
            jgg=jgg+1
           end do
          end do
          ibcoff=igs                                                     12d22s20
         end if                                                          12d22s20
        end do                                                           12d22s20
       end do                                                            12d22s20
       write(6,*)('jgg '),jgg-igg,nsing
       call prntm2(bc(igg),nsing,nrootu,nsing)
       ibcoff=igg                                                       1d25s21
      end if                                                            6d13s21
      return
      end
