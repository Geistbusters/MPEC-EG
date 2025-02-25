c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcds12(ihsdiag,nff1,iff1,nff22,nfdat,gd,vd,nsymb,mdon,    9d12s23
     $     mdoo,nec,multh,isymmrci,nvirt,ncsf,ncsf2,irel,ism,irefo,     12d18s20
     $     ixw1,ixw2,norb,nrootu,ih0av,nh0av,i3x,mdoub,ndoub,ldebug,    9d12s23
     $     maxbx,ionext,sr2,idenh0,id1x,id3x,idcont,idorth,veco,bc,ibc, 9d15s23
     $     nsdlk1,isblk1,nfdatd,ibufvcvds,idvmt,nder,igqqq)                   4d29s24
      implicit real*8 (a-h,o-z)                                         12d18s20
c
c     to do ...
c     consolidate densities, and perhaps njhere in density calculation
c     did I swap iarg and jarg in 4th gandc4?
c
      integer*8 ihsdiag(mdoo+1,nsymb,2),i18,i28,i38,i48,i1c,i1o,j2c,    12d18s20
     $     j2o,itestc,itesto,ipack8,last8(2),gandcc,gandco,gandcb       10d21s22
      integer*1 nab1(2),nab2(2)                                         12d18s20
      external second                                                   2d18s21
      logical l3x,lprt,lnew,ldebug,lchoice                              3d17s21
      equivalence (ipack8,ipack4)                                       12d18s20
      dimension nff1(mdoo+1,nsymb,2),iff1(*),nff22(mdoo+1,2,nsymb),     12d18s20
     $     nfdat(5,4,*),multh(8,8),nvirt(*),ncsf(*),ncsf2(4,*),         9d12s23
     $     irel(*),ism(*),irefo(*),i3x(*),vd(*),itest(32,3),            9d12s23
     $     nab4(2,3),ipack4(2),nl(4),npre(4),mpre(4),idenh0(*),                   9d12s23
     $     id1visv(8,8),nd1visv(8,8),nokdc(8,8,4),nok3v(4),             12d22s20
     $     idhvnotv(4),ndhvnotv(4),id1vnotv(4,8,8),nd1vnotv(4,8,8),     12d21s20
     $     id3vnotv3(4,2),nd3vnotv3(4,2),loff(4),idkeep(2),ndkeep(2),   1d4s21
     $     mdkeep(2),keep(2),idhvnotvf(4),ibmat(8),ivmat(8),            2d4s21
     $     mdhvnotv(4),md3vnotv3(4,2),nok4f(4),nok4(4),nok3vf(4),       12d22s20
     $     ih0av(*),nh0av(*),nok33f(4,2),nok33(4,2),veco(*),            9d13s23
     $     ionext(*),ioxx(2),ichvnotv(4),ic1vnotv(4,8,8),id1x(*),       9d12s23
     $     id3x(*),idcont(*),gd(*),idorth(4,*),isblk1(4,*),             9d18s23
     $     idchvnotv(8,4),iddch(8,4),nfdatd(8,*),sumdc(4),sumdh(4),     9d19s23
     $     id1xvisv(8,8),id1xvnotv(4,8,8,8),iddc1x(8,8,8,4),            9d26s23
     $     iddc3x(4),ibmatdc(8,4)                                       9d26s23
      include "common.store"                                            9d12s23
      common/kmfind/invk1(2,8,8,8,2)                                    6d30s18
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      common/cpucom/tovrx,top(10),tso(11)                                10d25s22
      data loopx/10000000/
      igoul=22713902
      xnan=0d0
      fctr=1d0                                                          4d29s24
      ncomp=0
      nucomp=0
      nrootm=nrootu-1                                                   1d4s21
      idoit=0                                                           3d1s21
      last8(1)=-1                                                       2d8s21
      ircv=ibcoff                                                       1d29s21
      iacc=ircv+mynprocg                                                1d30s21
      ivs=iacc+mynprocg                                                 1d30s21
      igg=ivs+maxbx                                                     1d30s21
      ibcoff=igg+maxbx                                                  1d30s21
      nacc=0                                                            1d30s21
      call enough('hcds.  1',bc,ibc)
      itransgg=0                                                        1d30s21
      loop=0
      nsing=0                                                           12d23s20
      norbx=norb+1                                                      12d18s20
      norbxx=norbx+1                                                    12d18s20
      norbxxx=norbxx+1                                                  12d18s20
      sumhh=0d0
      sumcc=0d0
      igoal=22740976
      do jsb=1,nsymb                                                    12d18s20
       jsbv=multh(jsb,isymmrci)                                         12d18s20
c
c     let us form bvrkn=[(vv"|nv')+p(k)(vv'|nv")]Vv'v"rk,               2d3s21
c     v ne v' and v"                                                    2d3s21
c     and space for vvrkn=Vvrj*Djkn                                     2d4s21
c     gandc4(Vsv,Vdv'v"), idx=1 (iv'|vv"), idx=2 (iv"|vv')
c
       ioffvd=1                                                         2d3s21
       ibcbmat=ibcoff                                                   2d3s21
       nbmat=0                                                          3d2s21
       do isb=1,nsymb                                                   3d2s21
        isbv12=multh(isb,isymmrci)                                      2d3s21
        isn=multh(isbv12,jsbv)                                          2d3s21
        ibmat(isb)=ibcoff                                               2d3s21
        if(min(irefo(isn),nvirt(jsbv)).gt.0)then                        2d3s21
         nftrip=nfdat(2,2,isb)+nfdat(2,3,isb)+nfdat(2,4,isb)             2d3s21
         nfh=nfdat(2,1,isb)+nftrip                                       2d3s21
         nn=nvirt(jsbv)*nrootu*nfh                                      2d3s21
         ibcoff=ibmat(isb)+nn*irefo(isn)                                3d2s21
         nbmat=nbmat+nn*irefo(isn)                                      3d2s21
        end if                                                          3d2s21
       end do                                                           3d2s21
       ibcvmat=ibcoff                                                   3d2s21
       do isb=1,nsymb                                                   3d2s21
        isbv12=multh(isb,isymmrci)                                      2d3s21
        isn=multh(isbv12,jsbv)                                          2d3s21
        ivmat(isb)=ibcoff                                               2d3s21
        if(min(irefo(isn),nvirt(jsbv)).gt.0)then                        2d3s21
         nftrip=nfdat(2,2,isb)+nfdat(2,3,isb)+nfdat(2,4,isb)             2d3s21
         nfh=nfdat(2,1,isb)+nftrip                                       2d3s21
         nn=nvirt(jsbv)*nrootu*nfh                                      2d3s21
         ibcoff=ivmat(isb)+nn*irefo(isn)                                3d2s21
        end if                                                          3d2s21
       end do                                                           3d2s21
       call enough('hcds.  2',bc,ibc)
       do i=ibcbmat,ibcoff-1                                            3d2s21
        bc(i)=0d0                                                       3d2s21
       end do                                                           3d2s21
       do isb=1,nsymb                                                   2d3s21
        isbv12=multh(isb,isymmrci)                                      2d3s21
        isn=multh(isbv12,jsbv)                                          2d3s21
        nftrip=nfdat(2,2,isb)+nfdat(2,3,isb)+nfdat(2,4,isb)             2d3s21
        nfh=nfdat(2,1,isb)+nftrip                                       2d3s21
        if(min(irefo(isn),nvirt(jsbv)).gt.0)then                        2d3s21
         nn=nvirt(jsbv)*nrootu*nfh                                      2d3s21
         nnn=nfh*nrootu                                                 2d3s21
         do isbv1=1,nsymb                                               2d3s21
          isbv2=multh(isbv1,isbv12)                                     2d3s21
          if(isbv2.ge.isbv1)then                                        2d3s21
           call ilimts(irefo(isn),nvirt(isbv1),mynprocg,mynowprog,il,ih,2d3s21
     $         i1s,i1e,i2s,i2e)                                         2d3s21
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
               jb=ibmat(isb)+nn*(i1-1)                                  2d3s21
               do jv=0,nvirt(jsbv)-1                                    2d4s21
                ix=max(jv,iv)                                           2d4s21
                in=min(jv,iv)                                           2d4s21
                irow=ii3+((ix*(ix+1))/2)+in                             2d4s21
                xint=bc(irow)*sr2                                       2d4s21
                do ir=0,nrootm                                           2d3s21
                 do k=0,nfdat(2,1,isb)-1                                2d3s21
                  iadv=ioffvd+iv+nvirt(isbv2)*(ir+nrootu*k)             2d3s21
                  bc(jb+k)=bc(jb+k)+xint*vd(iadv)                       2d3s21
                 end do                                                 2d3s21
                 jb=jb+nfh                                              2d3s21
                end do                                                  2d3s21
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
               jb=ibmat(isb)+nn*(i1-1)                                  2d3s21
               do jv=0,nvirt(jsbv)-1                                    2d3s21
                irow=ii3+jv+nvirt(jsbv)*iv                              2d3s21
                xint=bc(irow)*sr2                                       2d4s21
                do ir=0,nrootm                                           2d3s21
                 do k=0,nfdat(2,1,isb)-1                                2d3s21
                  iadv=ioffvd+iv+nvirt(isbv2)*(ir+nrootu*k)             2d3s21
                  bc(jb+k)=bc(jb+k)+xint*vd(iadv)                       2d3s21
                 end do                                                 2d3s21
                 jb=jb+nfh                                              2d3s21
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
               jb=ibmat(isb)+nn*(i1-1)                                  2d3s21
               do jv=0,nvirt(jsbv)-1                                    2d3s21
                irow=ii3+iv+nvirt(isbv1)*jv                              2d3s21
                xint=bc(irow)*sr2                                       2d4s21
                do ir=0,nrootm                                           2d3s21
                 do k=0,nfdat(2,1,isb)-1                                2d3s21
                  iadv=ioffvd+iv+nvirt(isbv2)*(ir+nrootu*k)             2d3s21
                  bc(jb+k)=bc(jb+k)+xint*vd(iadv)                       2d3s21
                 end do                                                 2d3s21
                 jb=jb+nfh                                              2d3s21
                end do                                                  2d3s21
               end do                                                   2d3s21
               ii3=ii3+nrow                                             2d3s21
              end do                                                    2d3s21
              i10=1                                                     2d3s21
             end do                                                     2d3s21
            end if                                                      2d3s21
            ioffvd=ioffvd+nvirt(isbv2)*nrootu*nfdat(2,1,isb)            2d3s21
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
               jb=ibmat(isb)+nn*(i1-1)                                   2d3s21
               do jv=0,nvirt(jsbv)-1                                    2d4s21
                ix=max(iv2,jv)                                          2d4s21
                in=min(iv2,jv)                                          2d4s21
                itri=((ix*(ix+1))/2)+in                                 2d4s21
                irec=jv+nvirt(jsbv)*iv2                                 2d3s21
                irow=ii3+itri+jsw*(irec-itri)                           2d3s21
                do ir=0,nrootm                                          2d3s21
                 do k=0,nfh-1                                           2d3s21
                  iadv=ioffvd+ivv+nvv*(ir+nrootu*k)                     2d3s21
                  bc(jb+k)=bc(jb+k)+bc(irow)*vd(iadv)                   2d3s21
                 end do                                                 2d3s21
                 jb=jb+nfh                                              2d3s21
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
               jb=ibmat(isb)+nn*(i1-1)                                   2d3s21
               do jv=0,nvirt(jsbv)-1                                    2d3s21
                ix=max(jv,iv2)                                          2d3s21
                in=min(jv,iv2)                                          2d3s21
                itri=((ix*(ix+1))/2)+in                                 2d3s21
                irec=jv+nvirt(jsbv)*iv2                                 2d3s21
                irow=ii3+itri+jsw*(irec-itri)                           2d3s21
                do ir=0,nrootm                                          2d3s21
                 do k=0,nfh-1                                           2d3s21
                  iadv=ioffvd+ivv+nvv*(ir+nrootu*k)                     2d3s21
                  bc(jb+k)=bc(jb+k)+bc(irow)*vd(iadv)                   2d3s21
                 end do                                                 2d3s21
                 jb=jb+nfh                                              2d3s21
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
                jb=ibmat(isb)+nnn*(jv+nvirt(jsbv)*i1m)                  2d4s21
                irow=ii3+iv2+nvirt(isbv2)*jv                            2d3s21
                do ir=0,nrootm                                          2d3s21
                 do k=0,nfh-1                                           2d3s21
                  iadv=ioffvd+ivv+nvv*(ir+nrootu*k)                     2d3s21
                  bc(jb+k)=bc(jb+k)+bc(irow)*vd(iadv)                   2d3s21
                 end do                                                 2d3s21
                 jb=jb+nfh                                              2d3s21
                end do                                                  2d3s21
               end do                                                   2d3s21
              end do                                                    2d3s21
              ii3=ii3+nrow                                              2d3s21
             end do                                                     2d3s21
             i10=1                                                      2d3s21
            end do                                                      2d3s21
           end if                                                       2d3s21
           call ilimts(irefo(isn),nvirt(isbv2),mynprocg,mynowprog,il,ih,2d3s21
     $         i1s,i1e,i2s,i2e)                                         2d3s21
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
               jb=ibmat(isb)+nn*(i1-1)                                   2d3s21
               do jv=0,nvirt(jsbv)-1                                    2d4s21
                ix=max(jv,iv1)                                          2d4s21
                in=min(jv,iv1)                                          2d4s21
                itri=((ix*(ix+1))/2)+in                                 2d4s21
                irec=jv+nvirt(jsbv)*iv1                                 2d3s21
                irow=ii3+itri+jsw*(irec-itri)                           2d3s21
                do ir=0,nrootm                                          2d3s21
                 do k=0,nfdat(2,1,isb)-1                                2d3s21
                  iadv=ioffvd+ivv+nvv*(ir+nrootu*k)                     2d3s21
                  bc(jb+k)=bc(jb+k)+bc(irow)*vd(iadv)                    2d3s21
                 end do                                                 2d3s21
                 do k=nfdat(2,1,isb),nfh-1                              2d3s21
                  iadv=ioffvd+ivv+nvv*(ir+nrootu*k)                     2d3s21
                  bc(jb+k)=bc(jb+k)-bc(irow)*vd(iadv)                    2d3s21
                 end do                                                 2d3s21
                 jb=jb+nfh                                              2d3s21
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
               jb=ibmat(isb)+nn*(i1-1)                                   2d3s21
               do jv=0,nvirt(jsbv)-1                                    2d3s21
                ix=max(jv,iv1)                                          2d3s21
                in=min(jv,iv1)                                          2d3s21
                itri=((ix*(ix+1))/2)+in                                 2d3s21
                irec=jv+nvirt(jsbv)*iv1                                 2d3s21
                irow=ii3+itri+jsw*(irec-itri)                           2d3s21
                do ir=0,nrootm                                          2d3s21
                 do k=0,nfdat(2,1,isb)-1                                2d3s21
                  iadv=ioffvd+ivv+nvv*(ir+nrootu*k)                     2d3s21
                  bc(jb+k)=bc(jb+k)+bc(irow)*vd(iadv)                    2d3s21
                 end do                                                 2d3s21
                 do k=nfdat(2,1,isb),nfh-1                              2d3s21
                  iadv=ioffvd+ivv+nvv*(ir+nrootu*k)                     2d3s21
                  bc(jb+k)=bc(jb+k)-bc(irow)*vd(iadv)                    2d3s21
                 end do                                                 2d3s21
                 jb=jb+nfh                                              2d3s21
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
                jb=ibmat(isb)+nnn*(jv+nvirt(jsbv)*i1m)                  2d4s21
                irow=ii3+iv1+nvirt(isbv1)*jv                            2d3s21
                do ir=0,nrootm                                          2d3s21
                 do k=0,nfdat(2,1,isb)-1                                2d3s21
                  iadv=ioffvd+ivv+nvv*(ir+nrootu*k)                     2d3s21
                  bc(jb+k)=bc(jb+k)+bc(irow)*vd(iadv)                   2d3s21
                 end do                                                 2d3s21
                 do k=nfdat(2,1,isb),nfh-1                              2d3s21
                  iadv=ioffvd+ivv+nvv*(ir+nrootu*k)                     2d3s21
                  bc(jb+k)=bc(jb+k)-bc(irow)*vd(iadv)                   2d3s21
                 end do                                                 2d3s21
                 jb=jb+nfh                                              2d3s21
                end do                                                  2d3s21
               end do                                                   2d3s21
              end do                                                    2d3s21
              ii3=ii3+nrow                                              2d3s21
             end do                                                     2d3s21
             i10=1                                                      2d3s21
            end do                                                      2d3s21
           end if                                                       2d3s21
           ioffvd=ioffvd+nvv*nrootu*nfh                                 2d3s21
          end if                                                        2d3s21
         end do                                                         2d3s21
         itmp=ibcoff                                                    2d3s21
         ibcoff=itmp+nn*irefo(isn)                                      2d3s21
         call enough('hcds.  3',bc,ibc)
         jbmat=ibmat(isb)                                               2d3s21
         do n=0,irefo(isn)-1                                            2d3s21
          do jv=0,nvirt(jsbv)-1                                         2d3s21
           do ir=0,nrootm                                               2d3s21
            do k=0,nfh-1                                                2d3s21
             iad=itmp+jv+nvirt(jsbv)*(ir+nrootu*(k+nfh*n))              2d3s21
             bc(iad)=bc(jbmat+k)                                        2d3s21
            end do                                                      2d3s21
            jbmat=jbmat+nfh                                             2d3s21
           end do                                                       2d3s21
          end do                                                        2d3s21
         end do                                                         2d3s21
         do i=0,nn*irefo(isn)-1                                         2d3s21
          bc(ibmat(isb)+i)=bc(itmp+i)                                   2d3s21
         end do                                                         2d3s21
         ibcoff=itmp                                                    6d10s21
        else                                                            2d3s21
         do isbv1=1,nsymb                                               2d3s21
          isbv2=multh(isbv1,isbv12)                                     2d3s21
          if(isbv2.ge.isbv1)then                                        2d3s21
           if(isbv1.eq.isbv2)then                                       2d3s21
            ioffvd=ioffvd+nvirt(isbv1)*nrootu*nfdat(2,1,isb)            2d3s21
            nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       2d3s21
           else                                                         2d3s21
            nvv=nvirt(isbv1)*nvirt(isbv2)                               2d3s21
           end if                                                       2d3s21
           ioffvd=ioffvd+nvv*nrootu*nfh                                 2d3s21
          end if                                                        2d3s21
         end do                                                         2d3s21
        end if                                                          2d3s21
       end do                                                           2d3s21
       call dws_gsumf(bc(ibcbmat),nbmat)                                3d2s21
c
c     for dcont, take v,r,nf,irefo order to v,nf,irefo,r,l order.
c
       do isb=1,nsymb                                                   9d26s23
        isbv12=multh(isb,isymmrci)                                      2d3s21
        isn=multh(isbv12,jsbv)                                          2d3s21
        if(min(irefo(isn),nvirt(jsbv)).gt.0)then                        2d3s21
         nftrip=nfdat(2,2,isb)+nfdat(2,3,isb)+nfdat(2,4,isb)             2d3s21
         nfh=nfdat(2,1,isb)+nftrip                                       2d3s21
         nn=nvirt(jsbv)*nrootu*nfh                                      2d3s21
         ibdcoff=0                                                      9d26s23
         do l=1,4                                                       9d26s23
          ibmatdc(isb,l)=ibcoff                                         9d26s23
          ibcoff=ibmatdc(isb,l)+nvirt(jsbv)*nrootu*nfdat(2,l,isb)       9d26s23
     $         *irefo(isn)                                              9d26s23
          do n=0,irefo(isn)-1                                            9d26s23
           do k=0,nfdat(2,l,isb)-1                                      9d26s23
            kp=k+ibdcoff                                                9d26s23
            do ir=0,nrootu-1                                             9d26s23
             iad1=ibmat(isb)+nvirt(jsbv)*(ir+nrootu*(kp+nfh*n))           9d26s23
             iad2=ibmatdc(isb,l)+nvirt(jsbv)*(k+nfdat(2,l,isb)*(n       9d26s23
     $            +irefo(isn)*ir))                                      9d26s23
             do iv=0,nvirt(jsbv)-1                                       9d26s23
              bc(iad2+iv)=bc(iad1+iv)                                    9d26s23
             end do                                                      9d26s23
            end do                                                       9d26s23
           end do                                                        9d26s23
          end do                                                         9d26s23
          ibdcoff=ibdcoff+nfdat(2,l,isb)                                9d26s23
         end do                                                         9d26s23
        end if                                                          3d2s21
       end do                                                           9d26s23
       do nclo1p=mdon+1,mdoo+1                                          12d18s20
        if(min(nff1(nclo1p,jsb,1),nvirt(jsbv)).gt.0)then                12d18s20
         nclo1=nclo1p-1                                                  12d12s20
         jarg=nclo1p-mdon                                                12d12s20
         nopen1=nec-2*nclo1                                              12d12s20
         ngg=nvirt(jsbv)*nrootu                                         12d15s20
         nggg=nff1(nclo1p,jsb,1)*ncsf(jarg)                             12d29s20
         nsing=nsing+nggg*nvirt(jsbv)                                   12d29s20
         ncolt=nggg*ngg                                                 12d29s20
         x2=dfloat(nff1(nclo1p,jsb,1))/dfloat(mynprocg)                 3d2s21
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
          idhvisv=ibcoff                                                10d27s22
          nnj=ncsf(jarg)*nfdat(2,1,isb)                                 10d27s22
          lgoal=multh(jsb,isb)                                          10d27s22
          lgoal3=multh(jsbv,isbv12)                                     10d27s22
          do isc=1,nsymb                                                10d27s22
           do l=1,4                                                     10d27s22
            nnl=ncsf(jarg)*nfdat(2,l,isb)*irefo(isc)                    10d27s22
            if(isc.eq.lgoal)then                                        10d27s22
             idhvnotvf(l)=ibcoff                                        10d27s22
             idhvnotv(l)=idhvnotvf(l)+nnl                               1d6s21
             ndhvnotv(l)=idhvnotv(l)+nnl                                12d21s20
             mdhvnotv(l)=ndhvnotv(l)+irefo(isc)                         12d31s20
             ibcoff=mdhvnotv(l)+nfdat(2,l,isb)                          12d21s20
            end if                                                      12d21s20
           end do                                                       12d18s20
           iscv=multh(isc,jsbv)                                         12d18s20
           jscv=multh(isc,lgoal)                                        12d21s20
           do isd=1,isc                                                 12d18s20
            iscdv=multh(iscv,isd)                                       12d18s20
            jscdv=multh(jscv,isd)                                       12d21s20
            if(isc.eq.isd)then                                          12d18s20
             nn=(irefo(isd)*(irefo(isd)+1))/2                           12d18s20
            else                                                        12d18s20
             nn=irefo(isd)*irefo(isc)                                   12d18s20
            end if                                                      12d18s20
            nnn=nn*irefo(jscdv)                                         12d18s20
            do l=1,4                                                    12d19s20
             id1vnotv(l,isd,isc)=ibcoff                                 12d19s20
             nd1vnotv(l,isd,isc)=id1vnotv(l,isd,isc)+nnn*ncsf(jarg)     12d19s20
     $            *nfdat(2,l,isb)                                       12d19s20
             ibcoff=nd1vnotv(l,isd,isc)+nnn                             12d21s20
            end do                                                      12d19s20
           end do                                                       12d18s20
          end do                                                        12d18s20
          do ipass=1,1                                                  3d1s21
           idkeep(ipass)=ibcoff                                         1d4s21
           ndkeep(ipass)=idkeep(ipass)+irefo(lgoal3)*nnt*ncsf(jarg)     1d4s21
           mdkeep(ipass)=ndkeep(ipass)+irefo(lgoal3)*nnt                1d5s21
           ibcoff=mdkeep(ipass)+irefo(lgoal3)*nnt                       1d5s21
           call enough('hcds.  4',bc,ibc)
          end do                                                        1d4s21
          do l=1,4
           nnl=ncsf(jarg)*nfdat(2,l,isb)*irefo(lgoal3)                  12d22s20
           id3vnotv3(l,1)=ibcoff                                        12d22s20
           nd3vnotv3(l,1)=id3vnotv3(l,1)+nnl                            12d22s20
           id3vnotv3(l,2)=nd3vnotv3(l,1)+irefo(lgoal3)                  12d22s20
           nd3vnotv3(l,2)=id3vnotv3(l,2)+nnl                            12d22s20
           md3vnotv3(l,1)=nd3vnotv3(l,2)+irefo(lgoal3)                  12d22s20
           md3vnotv3(l,2)=md3vnotv3(l,1)+nfdat(2,l,isb)                 12d22s20
           ibcoff=md3vnotv3(l,2)+nfdat(2,l,isb)                         12d22s20
          end do                                                        12d22s20
          ibcb4=ibcoff-1                                                12d18s20
c
c     for idcont ...
c     visv: isbv1=isbv2=jsbv, i.e. isbv12=1=multh(isb,isymmrci)
c
          idchvisv=ibcoff                                               9d18s23
          if(isbv12.eq.1)then                                           9d18s23
           ibcoff=idchvisv+nvirt(jsbv)*nrootu*nfdat(2,1,isb)*irefo(jsbv)9d18s23
c recall           lgoal=multh(jsb,isb)= isbv1(or 2) that isn't jsbv    9d20s23
c so integral has args isd,isc,q,lgoal, so q=isd*isc*lgoal
           do isc=1,nsymb                                               9d20s23
            jscv=multh(isc,lgoal)                                       9d20s23
            do isd=1,isc                                                9d20s23
             jscdv=multh(jscv,isd)                                      9d20s23
             if(isc.eq.isd)then                                         9d20s23
              nn=(irefo(isd)*(irefo(isd)+1))/2                          9d20s23
             else                                                       9d20s23
              nn=irefo(isd)*irefo(isc)                                  9d20s23
             end if                                                     9d20s23
             nnn=nn*irefo(multh(isd,multh(isc,jsbv)))                   9d20s23
             id1xvisv(isd,isc)=ibcoff                                   9d20s23
             ibcoff=id1xvisv(isd,isc)+nnn*nvirt(jsbv)*nrootu*           9d22s23
     $            nfdat(2,1,isb)                                        9d22s23
            end do                                                      9d20s23
           end do                                                       9d20s23
           call enough('hcds12.dchvisv',bc,ibc)                         9d18s23
           loffdnon=ioffdnon                                            9d18s23
           do isbv1=1,jsbv                                              9d18s23
            if(isbv1.eq.jsbv)then                                       9d18s23
             do k=0,nfdat(2,1,isb)-1                                    9d18s23
              do ia=0,irefo(jsbv)-1                                      9d18s23
               iadh=ih0av(jsbv)+irefo(jsbv)+nh0av(jsbv)*ia              9d18s23
               do ir=0,nrootu-1                                         9d18s23
                jdchvisv=idchvisv+nvirt(jsbv)*(k+nfdat(2,1,isb)*(ia     9d18s23
     $               +irefo(jsbv)*ir))                                  9d18s23
                iadl=loffdnon+nvirt(jsbv)*(ir+nrootu*k)                 9d18s23
                do iv=0,nvirt(jsbv)-1                                   9d18s23
                 bc(jdchvisv+iv)=bc(iadh+iv)*vd(iadl+iv)                9d18s23
                end do                                                  9d18s23
               end do                                                   9d18s23
              end do                                                    9d18s23
             end do                                                     9d18s23
c     dchvisv is v,k,a,r
             nvn=nvirt(jsbv)*nfdat(2,1,isb)                             9d20s23
c     dchvisv is v,k,b,r
c     js,js,lsa ... = jsb*isb = one not jsbv
             do isc=1,nsymb                                             9d20s23
              jscv=multh(isc,lgoal)                                     9d20s23
              do isd=1,isc                                              9d19s23
               jscdv=multh(isd,multh(isc,jsbv))                         9d20s23
               if(isc.eq.isd)then                                         9d20s23
                nn=(irefo(isd)*(irefo(isd)+1))/2                          9d20s23
               else                                                       9d20s23
                nn=irefo(isd)*irefo(isc)                                  9d20s23
               end if                                                     9d20s23
               i2eu=invk1(1,isd,isc,jscdv,2)                            9d20s23
               icase=invk1(2,isd,isc,jscdv,2)                           9d25s23
               nnn=nn*irefo(jscdv)                                      9d20s23
               do k=0,nfdat(2,1,isb)-1                                  9d20s23
                do ir=0,nrootu-1                                        9d20s23
                 iadl=loffdnon+nvirt(jsbv)*(ir+nrootu*k)                9d20s23
                 int1x=ionext(i2eu)                                     9d20s23
                 do iabc=0,nnn-1                                        9d20s23
                  jd1xvisv=id1xvisv(isd,isc)+nvirt(jsbv)*(k              9d20s23
     $                +nfdat(2,1,isb)*(iabc+nnn*ir))                    9d21s23
                  do iv=0,nvirt(jsbv)-1                                 9d20s23
                   bc(jd1xvisv+iv)=bc(int1x+iv)*vd(iadl+iv)             9d20s23
                  end do                                                9d20s23
                  int1x=int1x+nvirt(jsbv)                               9d20s23
                 end do                                                 9d20s23
                end do                                                  9d20s23
               end do                                                   9d20s23
              end do                                                    9d20s23
             end do                                                     9d20s23
            end if                                                      9d18s23
            loffdnon=loffdnon+nvirt(isbv1)*nrootu*nfdat(2,1,isb)        9d18s23
            nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       9d18s23
            do l=1,4
             loffdnon=loffdnon+nvv*nrootu*nfdat(2,l,isb)                9d18s23
            end do                                                      9d18s23
           end do                                                       9d18s23
          end if                                                        9d18s23
c
c     vnotv: one of isbv1 or isbv2 is jsbv
c
          loffdnon=ioffdnon                                             9d18s23
          do isbv1=1,nsymb                                              9d18s23
           isbv2=multh(isbv1,isbv12)                                    9d18s23
           if(isbv2.ge.isbv1)then                                       9d18s23
            if(isbv1.eq.isbv2)then                                      9d18s23
             loffdnon=loffdnon+nvirt(isbv1)*nrootu*nfdat(2,1,isb)           9d18s23
             nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                      9d18s23
            else                                                        9d18s23
             nvv=nvirt(isbv1)*nvirt(isbv2)                              9d18s23
            end if                                                      9d18s23
            factt=1d0                                                   9d18s23
            isbvu=0                                                     9d18s23
            if(isbv1.eq.jsbv)then                                       9d18s23
             isbvu=isbv2                                                9d18s23
            else if(isbv2.eq.jsbv)then                                  9d18s23
             isbvu=isbv1                                                9d18s23
            end if                                                      9d18s23
            do l=1,4                                                    9d18s23
             if(isbvu.ne.0.and.nfdat(2,l,isb).ne.0)then                 6d24s24
              if(isbv1.eq.isbv2)then                                    9d18s23
               isbvu=isbv1                                              9d18s23
               nrw=nvirt(isbv1)*nfdat(2,l,isb)*nrootu                   9d18s23
               nnp=nvirt(isbv1)*nfdat(2,l,isb)                          9d18s23
               idchvnotv(isbv1,l)=ibcoff                                9d18s23
               ibcoff=idchvnotv(isbv1,l)+nrw*irefo(isbv1)               9d20s23
               do isc=1,nsymb                                           9d20s23
                do isd=1,isc                                            9d20s23
                 if(isc.eq.isd)then                                     9d20s23
                  ncol=(irefo(isc)*(irefo(isc)+1))/2                    9d20s23
                 else                                                   9d20s23
                  ncol=irefo(isc)*irefo(isd)                            9d20s23
                 end if                                                 9d20s23
                 jscdv=multh(isd,multh(isc,isbvu))                      9d20s23
                 nrcol=ncol*irefo(jscdv)                                9d20s23
                 id1xvnotv(l,isd,isc,jscdv)=ibcoff                      9d25s23
                 ibcoff=id1xvnotv(l,isd,isc,jscdv)+nrw*nrcol            9d25s23
                end do                                                  9d20s23
               end do                                                   9d20s23
               itmp=ibcoff                                              9d20s23
               itmpb=itmp+nrw*nvirt(isbv1)                              9d20s23
               ibcoff=itmpb+nrw*nvirt(isbv1)                            9d20s23
               itmp2=ibcoff                                             9d18s23
               ibcoff=itmp2+nrw*irefo(isbv1)                            9d18s23
               call enough('hcds12.idchvnotv',bc,ibc)                   9d18s23
               do iz=idchvnotv(isbv1,l),itmp-1                          7d17s24
                bc(iz)=xnan                                             7d17s24
               end do                                                   7d17s24
               do iz=itmp,ibcoff-1                                      9d18s23
                bc(iz)=0d0                                              9d18s23
               end do                                                   9d18s23
               iad=loffdnon
               do k=0,nfdat(2,l,isb)-1                                  9d18s23
                do ir=0,nrootu-1                                        9d18s23
                 do iv2=1,nvirt(isbv2)-1                                9d18s23
                  iad2=itmp+nvirt(isbv1)*(k+nfdat(2,l,isb)*(ir          9d18s23
     $                 +nrootu*iv2))                                    9d18s23
                  do iv1=0,iv2-1                                        9d18s23
                   bc(iad2+iv1)=vd(iad+iv1)                             9d18s23
                  end do                                                9d18s23
                  iad=iad+iv2
                 end do                                                 9d18s23
                end do                                                  9d18s23
               end do                                                   9d18s23
               iad=ih0av(isbv1)+irefo(isbv1)                            9d18s23
               if(irefo(isbv1).gt.0)then                                9d20s23
                call dgemm('n','n',nrw,irefo(isbv1),nvirt(isbv1),1d0,    9d18s23
     $              bc(itmp),nrw,bc(iad),nh0av(isbv1),0d0,              9d18s23
     $              bc(itmp2),nrw,'hcds12.chvnotva')                    9d19s23
               end if                                                   9d20s23
               do iz=itmpb,itmp2-1                                      9d20s23
                bc(iz)=0d0                                              9d18s23
               end do                                                   9d18s23
               iad=loffdnon
               nnt=nvirt(isbv1)*nfdat(2,l,isb)*nrootu                   9d18s23
               do k=0,nfdat(2,l,isb)-1                                  9d18s23
                do ir=0,nrootu-1                                        9d18s23
                 do iv2=1,nvirt(isbv2)-1                                9d18s23
                  iad2=itmpb+iv2+nvirt(isbv1)*(k+nfdat(2,l,isb)*ir)     9d20s23
                  do iv1=0,iv2-1                                        9d18s23
                   bc(iad2)=vd(iad+iv1)                                 9d18s23
                   iad2=iad2+nnt                                        9d18s23
                  end do                                                9d18s23
                  iad=iad+iv2                                           9d18s23
                 end do                                                 9d18s23
                end do                                                  9d18s23
               end do                                                   9d18s23
c     tmp is v'krv
               if(irefo(isbv1).gt.0)then                                9d20s23
                iad=ih0av(isbv1)+irefo(isbv1)                            9d18s23
                call dgemm('n','n',nrw,irefo(isbv1),nvirt(isbv1),factt,  9d18s23
     $              bc(itmpb),nrw,bc(iad),nh0av(isbv1),1d0,             9d20s23
     $              bc(itmp2),nrw,'hcds12.chvnotvb')                    9d18s23
               end if                                                   9d20s23
c     tmp2 is v'krb want dchvnotv v'kbr
               sz=0d0
               do ib=0,irefo(isbv1)-1                                   9d18s23
                do ir=0,nrootu-1                                        9d18s23
                 iadrb=itmp2+nnp*(ir+nrootu*ib)                         9d18s23
                 iadbr=idchvnotv(isbv1,l)+nnp*(ib+irefo(isbv1)*ir)      9d18s23
                 do k=0,nnp-1                                           9d18s23
                  bc(iadbr+k)=bc(iadrb+k)                               9d18s23
                  sz=sz+bc(iadbr+k)**2
                 end do                                                 9d18s23
                end do                                                  9d18s23
               end do                                                   9d18s23
               ibcoff=itmp2                                             9d20s23
               do isc=1,nsymb                                           9d20s23
                do isd=1,isc                                            9d20s23
                 if(isc.eq.isd)then                                     9d20s23
                  ncol=(irefo(isc)*(irefo(isc)+1))/2
                 else                                                   9d20s23
                  ncol=irefo(isc)*irefo(isd)                            9d20s23
                 end if                                                 9d20s23
                 jscdv=multh(isd,multh(isc,isbv1))                      9d20s23
                 nrcol=ncol*irefo(jscdv)                                9d20s23
                 if(min(nrw,nrcol,nvirt(isbv1)).gt.0)then               9d20s23
                  itmp2=ibcoff                                           9d20s23
                  ibcoff=itmp2+nrw*nrcol                                 9d20s23
                  call enough('hcds12.nrw',bc,ibc)                       9d20s23
                  i2eu=invk1(1,isd,isc,jscdv,2)                          9d20s23
                  call dgemm('n','n',nrw,nrcol,nvirt(isbv1),1d0,         9d20s23
     $              bc(itmp),nrw,bc(ionext(i2eu)),nvirt(isbv1),0d0,     9d20s23
     $              bc(itmp2),nrw,'hcds12.c1xvnotva')                    9d19s23
                  call dgemm('n','n',nrw,nrcol,nvirt(isbv1),factt,       9d20s23
     $              bc(itmpb),nrw,bc(ionext(i2eu)),nvirt(isbv1),1d0,    9d26s23
     $              bc(itmp2),nrw,'hcds12.c1xvnotvb')                    9d19s23
c     tmp2 is v'krabc want dchvnotv v'kabcr
                  do ib=0,nrcol-1                                       9d20s23
                   do ir=0,nrootu-1                                        9d18s23
                    iadrb=itmp2+nnp*(ir+nrootu*ib)                         9d18s23
                    iadbr=id1xvnotv(l,isd,isc,jscdv)+nnp*(ib+nrcol*ir)  9d25s23
                    do k=0,nnp-1                                           9d18s23
                     bc(iadbr+k)=bc(iadrb+k)                               9d18s23
                    end do                                                 9d18s23
                   end do                                                  9d18s23
                  end do                                                   9d18s23
                  ibcoff=itmp2                                          9d20s23
                 end if                                                 9d20s23
                end do                                                  9d20s23
               end do                                                   9d20s23
               ibcoff=itmp                                              9d18s23
              else if(isbv1.eq.jsbv)then                                9d18s23
               isbvu=isbv2                                              9d18s23
               nnp=nvirt(isbv1)*nfdat(2,l,isb)                          9d18s23
               nrw=nvirt(jsbv)*nfdat(2,l,isb)*nrootu                    9d18s23
               idchvnotv(isbvu,l)=ibcoff                                9d18s23
               ibcoff=idchvnotv(isbvu,l)+nrw*irefo(isbvu)               9d18s23
               do isc=1,nsymb                                           9d20s23
                do isd=1,isc                                            9d20s23
                 if(isc.eq.isd)then                                     9d20s23
                  ncol=(irefo(isc)*(irefo(isc)+1))/2                    9d20s23
                 else                                                   9d20s23
                  ncol=irefo(isc)*irefo(isd)                            9d20s23
                 end if                                                 9d20s23
                 jsdc=multh(isd,multh(isc,isbvu))                       9d20s23
                 nrcol=ncol*irefo(jsdc)                                 9d20s23
                 id1xvnotv(l,isd,isc,jsdc)=ibcoff                       9d25s23
                 ibcoff=id1xvnotv(l,isd,isc,jsdc)+nrw*nrcol             9d25s23
                end do                                                  9d20s23
               end do                                                   9d20s23
               itmp=ibcoff                                              9d18s23
               ibcoff=itmp+nrw*nvirt(isbvu)                             9d18s23
               call enough('hcds12.idchvnotv',bc,ibc)                   9d18s23
               do iz=idchvnotv(isbvu,l),itmp-1                          7d17s24
                bc(iz)=xnan                                             7d17s24
               end do                                                   7d17s24
               iad=loffdnon
               do k=0,nfdat(2,l,isb)-1                                  9d18s23
                do ir=0,nrootu-1                                        9d18s23
                 do iv2=0,nvirt(isbv2)-1                                9d18s23
                  iad2=itmp+nvirt(isbv1)*(k+nfdat(2,l,isb)*(ir          9d18s23
     $                 +nrootu*iv2))                                    9d18s23
                  do iv1=0,nvirt(isbv1)-1                               9d18s23
                   bc(iad2+iv1)=vd(iad+iv1)                             9d18s23
                  end do                                                9d18s23
                  iad=iad+nvirt(isbv1)                                  9d18s23
                 end do                                                 9d18s23
                end do                                                  9d18s23
               end do                                                   9d18s23
c     tmp is v,k,r,v'
               iad=ih0av(isbvu)+irefo(isbvu)                            9d18s23
               if(irefo(isbvu).gt.0)then                                9d20s23
                itmp2=ibcoff                                            9d26s23
                ibcoff=itmp2+nrw*irefo(isbvu)                            9d18s23
                call enough('hcds12.tmp2',bc,ibc)                        9d18s23
                call dgemm('n','n',nrw,irefo(isbvu),nvirt(isbvu),1d0,    9d18s23
     $              bc(itmp),nrw,bc(iad),nh0av(isbvu),0d0,              9d18s23
     $              bc(itmp2),nrw,'hcds12.chvnotvc')                    9d19s23
c     tmp2 is v,k,r,b want idchvnotv to be v,k,b,r
                do ib=0,irefo(isbvu)-1                                   9d18s23
                 do ir=0,nrootu-1                                        9d18s23
                  iadrb=itmp2+nnp*(ir+nrootu*ib)                         9d18s23
                  iadbr=idchvnotv(isbvu,l)+nnp*(ib+irefo(isbvu)*ir)      9d18s23
                  do kv=0,nnp-1                                          9d18s23
                   bc(iadbr+kv)=bc(iadrb+kv)                             9d18s23
                  end do                                                 9d18s23
                 end do                                                  9d18s23
                end do                                                   9d18s23
                ibcoff=itmp2                                             9d20s23
               end if                                                   9d20s23
               do isc=1,nsymb                                           9d20s23
                do isd=1,isc                                            9d20s23
                 if(isc.eq.isd)then                                     9d20s23
                  ncol=(irefo(isc)*(irefo(isc)+1))/2                    9d20s23
                 else                                                   9d20s23
                  ncol=irefo(isc)*irefo(isd)                            9d20s23
                 end if                                                 9d20s23
                 jscdv=multh(isd,multh(isc,isbvu))                      9d20s23
                 nrcol=ncol*irefo(jscdv)                                9d20s23
                 if(min(nrcol,nrw,nvirt(isbvu)).gt.0)then               7d11s24
                  itmp2=ibcoff                                          9d20s23
                  ibcoff=itmp2+nrw*nrcol                                9d20s23
                  call enough('hcds12.1xa',bc,ibc)                      9d20s23
                  i2eu=invk1(1,isd,isc,jscdv,2)                         9d20s23
                  call dgemm('n','n',nrw,nrcol,nvirt(isbvu),1d0,        9d20s23
     $                 bc(itmp),nrw,bc(ionext(i2eu)),nvirt(isbvu),0d0,  9d25s23
     $                 bc(itmp2),nrw,'hcds12.1xa')                      9d20s23
                  do ib=0,nrcol-1                                       9d20s23
                   do ir=0,nrootu-1                                     9d20s23
                    iadrb=itmp2+nnp*(ir+nrootu*ib)                      9d20s23
                    iadbr=id1xvnotv(l,isd,isc,jscdv)+nnp*(ib+nrcol*ir)  9d25s23
                    do kv=0,nnp-1                                       9d20s23
                     bc(iadbr+kv)=bc(iadrb+kv)                          9d20s23
                    end do                                              9d20s23
                   end do                                               9d20s23
                  end do                                                9d20s23
                  ibcoff=itmp2                                          9d20s23
                 end if                                                 9d20s23
                end do                                                  9d20s23
               end do                                                   9d20s23
               ibcoff=itmp                                              9d18s23
              else                                                      9d18s23
               isbvu=isbv1                                              9d18s23
               nrw=nvirt(jsbv)*nfdat(2,l,isb)*nrootu                    9d18s23
               idchvnotv(isbvu,l)=ibcoff                                9d18s23
               ibcoff=idchvnotv(isbvu,l)+nrw*irefo(isbvu)               9d18s23
               do isc=1,nsymb                                           9d20s23
                do isd=1,isc                                            9d20s23
                 if(isc.eq.isd)then                                     9d20s23
                  ncol=(irefo(isc)*(irefo(isc)+1))/2                    9d20s23
                 else                                                   9d20s23
                  ncol=irefo(isc)*irefo(isd)                            9d20s23
                 end if                                                 9d20s23
                 jsdc=multh(isd,multh(isc,isbvu))                       9d20s23
                 nrcol=ncol*irefo(jsdc)                                 9d20s23
                 id1xvnotv(l,isd,isc,jsdc)=ibcoff                       9d25s23
                 ibcoff=id1xvnotv(l,isd,isc,jsdc)+nrw*nrcol             9d25s23
                end do                                                  9d20s23
               end do                                                   9d20s23
               itmp=ibcoff                                              9d18s23
               ibcoff=itmp+nrw*nvirt(isbvu)                             9d18s23
               call enough('hcds12.idchvnotv',bc,ibc)                   9d18s23
               do iz=idchvnotv(isbvu,l),itmp-1                          7d17s24
                bc(iz)=xnan                                             7d17s24
               end do                                                   7d17s24
               iad=loffdnon                                             9d18s23
               nnp=nvirt(isbv2)*nfdat(2,l,isb)                          9d18s23
               nnt=nnp*nrootu                                           9d18s23
               do k=0,nfdat(2,l,isb)-1                                  9d18s23
                do ir=0,nrootu-1                                        9d18s23
                 do iv2=0,nvirt(isbv2)-1                                9d18s23
                  iad2=itmp+iv2+nvirt(isbv2)*(k+nfdat(2,l,isb)*ir)      9d18s23
                  do iv1=0,nvirt(isbv1)-1                               9d18s23
                   bc(iad2)=vd(iad+iv1)                                 9d18s23
                   iad2=iad2+nnt                                        9d18s23
                  end do                                                9d18s23
                  iad=iad+nvirt(isbv1)                                  9d18s23
                 end do                                                 9d18s23
                end do                                                  9d18s23
               end do                                                   9d18s23
c     tmp is v',k,r,v
               iad=ih0av(isbvu)+irefo(isbvu)                            9d18s23
               if(irefo(isbvu).gt.0)then                                9d20s23
                itmp2=ibcoff                                             9d18s23
                ibcoff=itmp2+nrw*nvirt(isbvu)                            9d18s23
                call enough('hcds12.tmp2',bc,ibc)                        9d18s23
                szv=0d0
                szh=0d0
                do i=0,nrw*nvirt(isbvu)-1
                 szv=szv+bc(itmp+i)**2
                end do
                do j=0,irefo(isbvu)-1
                 do i=0,nvirt(isbvu)-1
                  ij=iad+i+nh0av(isbvu)*j
                  szh=szh+bc(ij)**2
                 end do
                end do
                call dgemm('n','n',nrw,irefo(isbvu),nvirt(isbvu),factt,  9d18s23
     $              bc(itmp),nrw,bc(iad),nh0av(isbvu),0d0,              9d18s23
     $              bc(itmp2),nrw,'hcds12.chvnotvd')                    9d18s23
c     tmp2 is  v',k,r,b want dchvnotv to be v',k,b,r
                do ib=0,irefo(isbvu)-1                                   9d18s23
                 do ir=0,nrootu-1                                        9d18s23
                  iadrb=itmp2+nnp*(ir+nrootu*ib)                         9d18s23
                  iadbr=idchvnotv(isbvu,l)+nnp*(ib+irefo(isbvu)*ir)      9d18s23
                  do kv=0,nnp-1                                          9d18s23
                   bc(iadbr+kv)=bc(iadrb+kv)                             9d18s23
                  end do                                                 9d18s23
                 end do                                                  9d18s23
                end do                                                   9d18s23
                ibcoff=itmp2                                             9d20s23
               end if                                                   9d20s23
               do isc=1,nsymb                                           9d20s23
                do isd=1,isc                                            9d20s23
                 if(isc.eq.isd)then                                     9d20s23
                  ncol=(irefo(isc)*(irefo(isc)+1))/2                    9d20s23
                 else                                                   9d20s23
                  ncol=irefo(isc)*irefo(isd)                            9d20s23
                 end if                                                 9d20s23
                 jscdv=multh(isd,multh(isc,isbvu))                      9d20s23
                 nrcol=ncol*irefo(jscdv)                                9d20s23
                 if(min(nrcol,nrw,nvirt(isbvu)).gt.0)then               8d5s24
                  itmp2=ibcoff                                          9d20s23
                  ibcoff=itmp2+nrw*nrcol                                9d20s23
                  call enough('hcds12.1xb',bc,ibc)                      9d20s23
                  i2eu=invk1(1,isd,isc,jscdv,2)                         9d20s23
                  call dgemm('n','n',nrw,nrcol,nvirt(isbvu),factt,        9d20s23
     $                 bc(itmp),nrw,bc(ionext(i2eu)),nvirt(isbvu),0d0,  9d25s23
     $                 bc(itmp2),nrw,'hcds12.1xb')                      9d20s23
                  do ib=0,nrcol-1                                       9d20s23
                   do ir=0,nrootu-1                                     9d20s23
                    iadrb=itmp2+nnp*(ir+nrootu*ib)                      9d20s23
                    iadbr=id1xvnotv(l,isd,isc,jscdv)+nnp*(ib+nrcol*ir)  9d25s23
                    do kv=0,nnp-1                                       9d20s23
                     bc(iadbr+kv)=bc(iadrb+kv)                          9d20s23
                    end do                                              9d20s23
                   end do                                               9d20s23
                  end do                                                9d20s23
                  ibcoff=itmp2                                          9d20s23
                 end if                                                 9d20s23
                end do                                                  9d20s23
               end do                                                   9d20s23
               ibcoff=itmp                                              9d18s23
              end if                                                    9d18s23
             end if                                                     9d18s23
             factt=-1d0                                                 9d18s23
             loffdnon=loffdnon+nvv*nrootu*nfdat(2,l,isb)                    9d18s23
            end do                                                      9d18s23
           end if                                                       9d18s23
          end do                                                        9d18s23
          if1o=nff1(nclo1p,jsb,2)                                        12d14s20
          jvs=ivs                                                       12d18s20
          do if1=1,nff1(nclo1p,jsb,1)                                    12d14s20
           ist=1+ncsf(jarg)*(if1-1)                                     1d5s21
           ien=ncsf(jarg)*if1                                           1d5s21
           istu=max(ist,i1l)                                            1d5s21
           ienu=min(ien,i1h)                                            1d5s21
           if(ienu.ge.istu)then                                         1d5s21
            i11s=istu-ist+1                                             1d5s21
            njhere=ienu+1-istu                                          1d5s21
           else                                                         1d5s21
            njhere=0                                                    1d5s21
           end if                                                       1d5s21
           do l=1,4
            sumdc(l)=0d0
            sumdh(l)=0d0                                                9d19s23
           end do
           if(itransvs.eq.0)then                                        1d29s21
            call ddi_done(ibc(ircv),nrcv)                                    1d29s21
            itransvs=1                                                  1d29s21
            itmp=ibcoff                                                    1d27s21
            ibcoff=itmp+ngg                                                1d27s21
            call enough('hcds. 10',bc,ibc)
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
c
c     factor of 2 because we only do ds densities rather than ds and sd.
c
              bc(jvs0+j)=bc(itmp+j)*2d0                                 9d27s23
             end do                                                        1d27s21
            end do                                                         1d27s21
            ibcoff=itmp                                                    1d27s21
           end if                                                       1d29s21
           if(njhere.gt.0)then                                          9d18s23
            itmp1=ibcoff                                                9d18s23
            ibcoff=itmp1+nvirt(jsbv)*nrootu*njhere                      9d18s23
            call enough('hcds12.tmp1',bc,ibc)                           9d18s23
            nnt=njhere*nvirt(jsbv)                                      9d18s23
            do ir=0,nrootu-1                                            9d18s23
             do j=0,njhere-1                                            9d19s23
              iad=jvs+nvirt(jsbv)*(ir+nrootu*(j+i11s-1))                9d18s23
              iad1=itmp1+j+nnt*ir                                       9d18s23
              do iv=0,nvirt(jsbv)-1                                     9d18s23
               bc(iad1)=bc(iad+iv)                                      9d18s23
               iad1=iad1+njhere                                         9d19s23
              end do                                                    9d18s23
             end do                                                     9d18s23
            end do                                                      9d18s23
            sz1=sz
c     tmp1 is j,v,r
            do isbv1=1,nsymb                                            9d18s23
             isbv2=multh(isbv1,isbv12)                                  9d18s23
             if(isbv2.ge.isbv1)then                                     9d18s23
              isbvu=0                                                   9d18s23
              if(isbv1.eq.jsbv)then                                     9d18s23
               isbvu=isbv2                                              9d18s23
              else if(isbv2.eq.jsbv)then                                9d18s23
               isbvu=isbv1                                              9d18s23
              end if                                                    9d18s23
              if(isbvu.ne.0)then                                        9d20s23
               nn=njhere*irefo(isbvu)*nrootu                              9d18s23
c recall           lgoal=multh(jsb,isb)= isbv1(or 2) that isn't jsbv    9d20s23
c so integral has args isd,isc,q,lgoal, so q=isd*isc*lgoal
               do l=1,4                                                 9d18s23
                iddch(isbvu,l)=ibcoff                                   9d18s23
                ibcoff=iddch(isbvu,l)+nn*nfdat(2,l,isb)                 9d18s23
                isn=multh(isbv12,jsbv)                                  9d26s23
                if(irefo(isn).gt.0)then                                 9d26s23
                 iddc3x(l)=ibcoff                                       9d26s23
                 ibcoff=iddc3x(l)+njhere*nfdat(2,l,isb)*nrootu          9d26s23
     $                 *irefo(isn)                                      9d26s23
                end if                                                  9d26s23
                do isc=1,nsymb                                          9d20s23
                 do isd=1,isc                                           9d20s23
                  if(isc.eq.isd)then
                   ncol=(irefo(isc)*(irefo(isc)+1))/2                   9d20s23
                  else                                                  9d20s23
                   ncol=irefo(isc)*irefo(isd)                           9d20s23
                  end if                                                9d20s23
                  jsdc=multh(isd,multh(isc,isbvu))                      9d20s23
                  nrcol=ncol*irefo(jsdc)                                9d20s23
                  iddc1x(isd,isc,jsdc,l)=ibcoff                         9d20s23
                  ibcoff=iddc1x(isd,isc,jsdc,l)+njhere*nrcol*nrootu     9d21s23
     $                 *nfdat(2,l,isb)                                  9d21s23
                 end do                                                 9d20s23
                end do                                                  9d20s23
               end do                                                   9d18s23
               call enough('hcds12.be4',bc,ibc)                         9d18s23
               do iz=iddch(isbvu,1),ibcoff-1
                bc(iz)=xnan
               end do
               ff=0d0                                                   9d18s23
c                    jdenj=id1vnotv(l,js,js)+nnl*icolj                    12d19s20
               if(isbv12.eq.1.and.isbv1.eq.jsbv)then                    9d18s23
                if(irefo(isbvu).gt.0)then                               9d20s23
                 ncl=nfdat(2,1,isb)*irefo(isbvu)                         9d18s23
                 jtmp1=itmp1                                                9d18s23
                 jdchvisv=idchvisv                                          9d18s23
                 jans=iddch(jsbv,1)                                         9d18s23
c     dchvisv is v,k,b,r
c     js,js,lsa ... = jsb*isb = one not jsbv
                 do ir=0,nrootu-1                                           9d18s23
                  call dgemm('n','n',njhere,ncl,nvirt(jsbv),sr2,         9d18s23
     $                bc(jtmp1),njhere,bc(jdchvisv),nvirt(jsbv),0d0,    9d18s23
     $                bc(jans),njhere,'hcds12.jansa')                    9d18s23
                  jtmp1=jtmp1+njhere*nvirt(jsbv)                            9d18s23
                  jdchvisv=jdchvisv+nvirt(jsbv)*ncl                         9d18s23
                  jans=jans+njhere*ncl                                   9d18s23
                 end do                                                     9d18s23
                end if                                                  9d20s23
                do isc=1,nsymb                                          9d20s23
                 do isd=1,isc                                           9d20s23
                  if(isc.eq.isd)then                                    9d20s23
                   ncol=(irefo(isc)*(irefo(isc)+1))/2                   9d20s23
                  else                                                  9d20s23
                   ncol=irefo(isc)*irefo(isd)                           9d20s23
                  end if                                                9d20s23
                  jsdc=multh(isd,multh(isc,isbvu))                      9d20s23
                  nrcol=ncol*irefo(jsdc)                                9d20s23
                  if(min(nrcol,nvirt(jsbv)).gt.0)then                   9d20s23
                   ncl=nfdat(2,1,isb)*nrcol                             9d20s23
                   jtmp1=itmp1                                           9d20s23
                   jd1xvisv=id1xvisv(isd,isc)                            9d20s23
                   jans=iddc1x(isd,isc,jsdc,1)                           9d20s23
                   do ir=0,nrootu-1                                      9d20s23
                    call dgemm('n','n',njhere,ncl,nvirt(jsbv),sr2,       9d20s23
     $                  bc(jtmp1),njhere,bc(jd1xvisv),nvirt(jsbv),0d0,  9d20s23
     $                  bc(jans),njhere,'hcds12.jand1x')                9d20s23
                    jtmp1=jtmp1+njhere*nvirt(jsbv)                       9d20s23
                    jans=jans+njhere*ncl                                 9d20s23
                    jd1xvisv=jd1xvisv+nvirt(jsbv)*ncl                    9d20s23
                   end do                                                9d20s23
                  end if                                                9d20s23
                 end do                                                 9d20s23
                end do                                                  9d20s23
                ff=1d0                                                  9d18s23
               end if                                                      9d18s23
               do l=1,4                                                 9d18s23
                if(min(nfdat(2,l,isb),irefo(isbvu)).gt.0)then           6d24s24
                 ncl=nfdat(2,l,isb)*irefo(isbvu)                         9d18s23
                 jtmp1=itmp1                                             9d18s23
                 jdchvnotv=idchvnotv(isbvu,l)                            9d18s23
c     v'kbr
                 jans=iddch(isbvu,l)                                     9d18s23
                 do ir=0,nrootu-1                                           9d18s23
                  call dgemm('n','n',njhere,ncl,nvirt(jsbv),1d0,         9d18s23
     $                bc(jtmp1),njhere,bc(jdchvnotv),nvirt(jsbv),ff,    9d18s23
     $                bc(jans),njhere,'hcds12.jansb')                    9d18s23
c     so idcch is jkbr
                  jtmp1=jtmp1+njhere*nvirt(jsbv)                            9d18s23
                  jdchvnotv=jdchvnotv+nvirt(jsbv)*ncl                    9d18s23
                  jans=jans+njhere*ncl                                   9d18s23
                 end do                                                  9d20s23
                end if                                                  9d20s23
                isn=multh(isbv12,jsbv)                                  9d26s23
                if(min(nvirt(jsbv),irefo(isn),nfdat(2,l,isb)).gt.0)then 9d26s23
c     ibmatdc(isb,l) is nvirt(jsbv),nfdat,irefo(isn),nrootu
                 jbmatdc=ibmatdc(isb,l)                                 9d26s23
                 jddc3x=iddc3x(l)                                       9d26s23
                 ncl=nfdat(2,l,isb)*irefo(isn)                          9d26s23
                 jtmp1=itmp1                                            9d26s23
                 do ir=0,nrootu-1                                       9d26s23
                  call dgemm('n','n',njhere,ncl,nvirt(jsbv),1d0,        9d26s23
     $                 bc(jtmp1),njhere,bc(jbmatdc),nvirt(jsbv),0d0,    9d26s23
     $                 bc(jddc3x),njhere,'hcds12.jddc3x')               9d26s23
                  jtmp1=jtmp1+njhere*nvirt(jsbv)                        9d26s23
                  jbmatdc=jbmatdc+nvirt(jsbv)*ncl                       9d26s23
                  jddc3x=jddc3x+njhere*ncl                              9d26s23
                 end do                                                 9d26s23
                end if                                                  9d26s23
                do isc=1,nsymb                                          9d20s23
                 do isd=1,isc                                           9d20s23
                  if(isc.eq.isd)then                                    9d20s23
                   ncol=(irefo(isc)*(irefo(isc)+1))/2                   9d20s23
                  else                                                  9d20s23
                   ncol=irefo(isc)*irefo(isd)                           9d20s23
                  end if                                                9d20s23
                  jsdc=multh(isd,multh(isc,isbvu))                      9d20s23
                  nrcol=ncol*irefo(jsdc)                                9d20s23
                  if(min(nrcol,nfdat(2,l,isb)).gt.0)then                9d20s23
                   jtmp1=itmp1                                           9d20s23
                   jd1xvnotv=id1xvnotv(l,isd,isc,jsdc)                  9d25s23
                   jans=iddc1x(isd,isc,jsdc,l)                          9d20s23
                   ncl=nrcol*nfdat(2,l,isb)                             9d20s23
                   do ir=0,nrootu-1                                      9d20s23
                    call dgemm('n','n',njhere,ncl,nvirt(jsbv),1d0,      9d20s23
     $                   bc(jtmp1),njhere,bc(jd1xvnotv),nvirt(jsbv),    9d20s23
     $                   ff,bc(jans),njhere,'hcds12.jans1x')            9d20s23
                    jtmp1=jtmp1+njhere*nvirt(jsbv)                      9d20s23
                    jd1xvnotv=jd1xvnotv+nvirt(jsbv)*ncl                 9d20s23
                    jans=jans+njhere*ncol                               9d20s23
                   end do                                               9d20s23
                  end if                                                9d20s23
                 end do                                                 9d20s23
                end do                                                     9d18s23
                ff=0d0                                                  9d20s23
               end do                                                   9d18s23
              end if                                                    9d18s23
             end if                                                     9d18s23
            end do                                                      9d18s23
           end if                                                       9d18s23
           idoit=idoit+1                                                3d1s21
           joffdnon=ioffdnon                                            12d21s20
           do i=idhvisv,ibcb4                                           12d18s20
            bc(i)=0d0                                                   12d18s20
           end do                                                       12d18s20
           if1oc=if1o                                                   10d27s22
           do nclo2p=max(mdon+1,nclo1p-2),min(mdoo+1,nclo1p+3)           12d18s20
            if(nff22(nclo2p,1,isb).gt.0)then                             12d18s20
             if1o=if1oc                                                 10d27s22
             nclo2=nclo2p-1                                              12d18s20
             iarg=nclo2p-mdon                                            12d18s20
             nopen2=nec-2*nclo2p                                         12d18s20
             nopen2p=nopen2+2                                            12d18s20
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
             jdcont=idcont(isb)+nfdatd(isb,nclo2p)                      9d18s23
             do if2=1,nff22(nclo2p,1,isb)                               12d18s20
              ipack8=ibc(jvcv)                                          12d18s20
              j2c=ipack4(1)                                             12d18s20
              j2o=ipack4(2)                                             12d18s20
              nclo=popcnt(ipack4(1))                                    12d18s20
              nspace=ibc(jvcv+1)                                        3d19s21
              lchoice=.false.                                           3d18s21
              ndcont=0                                                  9d18s23
              do l=1,4                                                   3d17s21
               nl(l)=ibc(jvcv+1+l)                                      3d19s21
               ndcont=ndcont+nl(l)*ncsf2(l,iarg)                        9d18s23
               if(nl(l).gt.ncsf2(l,iarg))lchoice=.true.                 1d12s23
              end do                                                     3d17s21
              ndcont=ndcont*nrootu                                      9d18s23
              lchoice=.true.
              j2o=ibset(j2o,norbx)                                      12d18s20
              j2o=ibset(j2o,norbxx)                                     12d18s20
              if(njhere.gt.0)then                                       2d26s21
               gandcc=ieor(i1c,j2c)                                     10d13s22
               gandco=ieor(i1o,j2o)                                     10d13s22
               gandcb=ior(gandcc,gandco)                                10d21s22
               ndifb=popcnt(gandcb)                                     10d21s22
               if(ndifb.le.4)then                                       10d21s22
                ndifd=popcnt(gandcc)                                     10d13s22
                ndifs=popcnt(gandco)                                     10d13s22
                iprod=-1
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
                  call enough('hcds.  5',bc,ibc)
                  call prodn(iwpb1,iwpk1,ncsf(jarg),ncsf(iarg),ncsfmid1, 3d19s21
     $                bc(iprod),bc,ibc,1d0,0d0)                         2d13s23
                 end if                                                  3d19s21
                 jprod=iprod                                             3d19s21
                 ldcont=jdcont                                          9d18s23
                 do l=1,4                                                12d18s20
                  if(nl(l).gt.0)then                                     12d18s20
                   nnl=ncsf(jarg)*nfdat(2,l,isb)                         12d19s20
                   iad1=jvcv+ibc(jvcv+5+l)                               3d19s21
                   iad2=iad1+nl(l)                                       3d19s21
                   itmp=ibcoff                                           12d18s20
                   ibcoff=itmp+ncsf(jarg)*nl(l)                          12d18s20
                   itmpdc1=ibcoff                                       9d18s23
                   itmpdc2=itmpdc1+ncsf(jarg)*nl(l)*nrootu              9d18s23
                   ibcoff=itmpdc2+ncsf2(l,iarg)*nl(l)*nrootu            9d18s23
                   call enough('hcds12.tmpdc1',bc,ibcoff)               9d18s23
                   do iz=itmpdc1,itmpdc2-1                              9d18s23
                    bc(iz)=0d0                                          9d18s23
                   end do                                               9d18s23
                   nnt=nl(l)*nrootu                                     9d18s23
                   sz=0d0
c     k,r,j
                   if(sz.gt.1d-10)write(6,*)('tmpddc1 '),sz,lcoice
                   if(lchoice)then                                       3d19s21
                    call dgemm('n','n',ncsf(jarg),nl(l),ncsf2(l,iarg),   3d19s21
     $                 1d0,bc(jprod),ncsf(jarg),bc(iad2),ncsf2(l,iarg), 3d19s21
     $                 0d0,bc(itmp),ncsf(jarg),                         3d19s21
     d' hcds12.  1a')
                   else                                                  3d19s21
                    call xtimesn2(ncsf(jarg),ncsf(iarg),ncsfmid1,nl(l),   2d25s21
     $                 iwpb1,iwpk1,bc(iad2),ncsf2(l,iarg),bc(itmp),     2d25s21
     $                 ncsf(jarg),1d0,0d0,iargo,ncsf2(l,iarg),bc,ibc)   11d10s22
                   end if                                                3d19s21
c     k,r,i
                   call hcds12c(nrootu,nl(l),ibc(iad1),iddch(lsa,l),    9d20s23
     $                  njhere,nfdat(2,l,isb),lga,irefo(lsa),itmpdc1,   9d20s23
     $                  itmpdc2,i11s,nnt,ncsf2(l,iarg),ncsf(jarg),      9d20s23
     $             bc(jprod),ldcont,iad2,sumdc(l),fctr,bc,ibc,igqqq,'a')     4d29s24
c     jsbv=jsb*isymmrci
c     isbv12=isb*isymmrci
c     isbv(1or2)=isbv12*jsbv=jsb*isymmrci*isb*isymmrci=jsb*isb
                   lgoal=multh(jsb,isb)                                  12d21s20
                   jden=idhvnotv(l)+nnl*lga                              12d21s20
                   ibc(ndhvnotv(l)+lga)=1                                12d21s20
                   jtmp=itmp                                             12d18s20
                   do iii=0,nl(l)-1                                      12d18s20
                    ip=ibc(iad1+iii)-1                                   12d18s20
                    jjden=jden+ncsf(jarg)*ip                             12d18s20
                    ibc(mdhvnotv(l)+ip)=1                                12d21s20
                    do j=0,ncsf(jarg)-1                                  12d18s20
                     bc(jjden+j)=bc(jjden+j)+bc(jtmp+j)                  12d18s20
                    end do                                               12d18s20
                    jtmp=jtmp+ncsf(jarg)                                 12d18s20
                   end do                                                12d18s20
                   do i=1,nok-1                                           12d19s20
                    js=ism(itest(i,2))                                           11d13s20
                    jg=irel(itest(i,2))-1                                        11d13s20
                    icolj=((jg*(jg+1))/2)+jg                                12d14s20
                    nn=(irefo(js)*(irefo(js)+1))/2                        12d18s20
                    icolj=icolj+nn*lga                                    12d18s20
                    jdenj=id1vnotv(l,js,js)+nnl*icolj                    12d19s20
                    nxxj=nn*irefo(lsa)                                  9d20s23
                    ibc(nd1vnotv(l,js,js)+icolj)=1                       12d19s20
                    if(itest(i,3).eq.2)then                                  12d14s20
                     if(lsa.ne.js)then                                    12d18s20
                      nn=irefo(lsa)*irefo(js)                             12d18s20
                      nxxk=nn*irefo(js)                                 9d21s23
                      if(itest(i,2).lt.nab1(1))then                       12d18s20
                       icolk=jg+irefo(js)*lga                             12d18s20
                       icolk=icolk+nn*jg                                    12d18s20
                       jdenk=id1vnotv(l,js,lsa)+nnl*icolk                12d19s20
                       ibc(nd1vnotv(l,js,lsa)+icolk)=1                   12d19s20
                       iddk=iddc1x(js,lsa,js,l)                          9d20s23
                       icase=invk1(2,js,lsa,js,2)
                       if(icase.eq.1)then                               9d25s23
                        icoldc=icolk                                    9d25s23
                       else                                             9d25s23
                        icoldc=lga+irefo(lsa)*(jg+irefo(js)*jg)         9d25s23
                       end if                                           9d25s23
                      else                                                12d18s20
                       icolk=lga+irefo(lsa)*jg                            12d18s20
                       icolk=icolk+nn*jg                                    12d18s20
                       jdenk=id1vnotv(l,lsa,js)+nnl*icolk                    12d18s20
                       ibc(nd1vnotv(l,lsa,js)+icolk)=1                   12d19s20
                       iddk=iddc1x(lsa,js,js,l)                         9d22s23
                       icase=invk1(2,lsa,js,js,2)
                       if(icase.eq.1)then                               9d25s23
                        icoldc=icolk                                    9d25s23
                       else                                             9d25s23
                        icoldc=jg+irefo(js)*(lga+irefo(lsa)*jg)         9d25s23
                       end if                                           9d25s23
                      end if                                              12d18s20
                     else                                                 12d18s20
                      ix=max(jg,lga)                                      12d18s20
                      in=min(jg,lga)                                      12d18s20
                      icolk=((ix*(ix+1))/2)+in+nn*jg                      12d18s20
                      jdenk=id1vnotv(l,js,js)+nnl*icolk                  12d19s20
                      ibc(nd1vnotv(l,js,js)+icolk)=1                     12d19s20
                      iddk=iddc1x(js,js,lsa,l)                           9d20s23
                      nxxk=nxxj                                         9d22s23
                      icoldc=icolk                                      9d25s23
                     end if                                               12d18s20
                     call hcds12c(nrootu,nl(l),ibc(iad1),                9d20s23
     $                  iddc1x(js,js,lsa,l),njhere,nfdat(2,l,isb),icolj,9d20s23
     $                   nxxj,itmpdc1,itmpdc2,i11s,nnt,ncsf2(l,iarg),    9d20s23
     $               ncsf(jarg),bc(jprod),ldcont,iad2,sumdc(l),2d0*fctr,4d29s24
     $                    bc,ibc,igqqq,'b')                                       9d20s23
                     call hcds12c(nrootu,nl(l),ibc(iad1),                9d20s23
     $                   iddk,njhere,nfdat(2,l,isb),icoldc,             9d25s23
     $                   nxxk,itmpdc1,itmpdc2,i11s,nnt,ncsf2(l,iarg),    9d20s23
     $                  ncsf(jarg),bc(jprod),ldcont,iad2,sumdc(l),-fctr,4d29s24
     $                    bc,ibc,igqqq,'c')                                       9d20s23
                     jtmp=itmp                                              12d18s20
                     do iii=0,nl(l)-1                                         12d18s20
                      ii=ibc(iad1+iii)-1                                      12d18s20
                      jdj=jdenj+ncsf(jarg)*ii                             12d18s20
                      jdk=jdenk+ncsf(jarg)*ii                             12d18s20
                      do j=0,ncsf(jarg)-1                                   12d18s20
                       bc(jdj+j)=bc(jdj+j)+2d0*bc(jtmp+j)                 12d18s20
                       bc(jdk+j)=bc(jdk+j)-bc(jtmp+j)                     12d18s20
                      end do                                                12d18s20
                      jtmp=jtmp+ncsf(jarg)                                  12d18s20
                     end do                                                 12d18s20
                    else                                                     12d14s20
                     call hcds12c(nrootu,nl(l),ibc(iad1),                9d20s23
     $                  iddc1x(js,js,lsa,l),njhere,nfdat(2,l,isb),icolj,9d20s23
     $                   nxxj,itmpdc1,itmpdc2,i11s,nnt,ncsf2(l,iarg),    9d20s23
     $                   ncsf(jarg),bc(jprod),ldcont,iad2,sumdc(l),fctr,4d29s24
     $                    bc,ibc,igqqq,'d')                                       9d20s23
                     jtmp=itmp                                              12d18s20
                     do iii=0,nl(l)-1                                         12d18s20
                      ii=ibc(iad1+iii)-1                                      12d18s20
                      jdj=jdenj+ncsf(jarg)*ii                             12d18s20
                      do j=0,ncsf(jarg)-1                                   12d18s20
                       bc(jdj+j)=bc(jdj+j)+bc(jtmp+j)                     12d18s20
                      end do                                                12d18s20
                      jtmp=jtmp+ncsf(jarg)                                  12d18s20
                     end do                                                 12d18s20
                    end if                                                   12d14s20
                   end do                                                    12d14s20
                   ibcoff=itmp                                           12d18s20
                   ldcont=ldcont+nrootu*nl(l)*ncsf2(l,iarg)             9d18s23
                  end if                                                 12d18s20
                  jprod=jprod+ncsf2(l,iarg)*ncsf(jarg)                   3d19s21
                  iargo=iargo+ncsf2(l,iarg)                              2d25s21
                 end do                                                  12d18s20
                 ibcoff=iprod                                            3d19s21
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
                    if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                     if(max(nab2(1),nab2(2)).le.norb)then                2d26s21
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
                      if(lsa.eq.lsb)then                                  12d18s20
                       nn=(irefo(lsa)*(irefo(lsa)+1))/2                   12d18s20
                       icol=((lgb*(lgb+1))/2)+lga+nn*lgc                  12d18s20
                       icoldc=icol                                      9d25s23
                      else                                                12d18s20
                       nn=irefo(lsa)*irefo(lsb)                           12d18s20
                       icol=lga+irefo(lsa)*(lgb+irefo(lsb)*lgc)           12d18s20
                       icase=invk1(2,lsa,lsb,lsc,2)                     9d25s23
                       icoldc=icol                                      9d25s23
                       if(icase.eq.2)icoldc=lgb+irefo(lsb)*(lga         9d25s23
     $                      +irefo(lsa)*lgc)                            9d25s23
                      end if                                              12d18s20
                      nxx=nn*irefo(lsc)                                 9d20s23
                      iprod=ibcoff                                       3d19s21
                      if(lchoice)then                                    3d19s21
                       call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),     3d19s21
     $                  ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod,11d10s22
     $                     bc,ibc)                                      11d10s22
                      end if                                             3d19s21
                      iargo=1                                            2d25s21
                      jprod=iprod                                        3d19s21
                      ldcont=jdcont                                     9d20s23
                      do l=1,4
                       if(nl(l).gt.0)then                                     12d18s20
                        nnl=ncsf(jarg)*nfdat(2,l,isb)                         12d19s20
                        iad1=jvcv+ibc(jvcv+5+l)                          3d19s21
                        iad2=iad1+nl(l)                                  3d19s21
                        itmp=ibcoff                                           12d18s20
                        ibcoff=itmp+ncsf(jarg)*nl(l)                          12d18s20
                        itmpdc1=ibcoff                                       9d18s23
                        itmpdc2=itmpdc1+ncsf(jarg)*nl(l)*nrootu              9d18s23
                        ibcoff=itmpdc2+ncsf2(l,iarg)*nl(l)*nrootu            9d18s23
                        call enough('hcds.  7',bc,ibc)
                        nnt=nl(l)*nrootu                                     9d18s23
                        call hcds12c(nrootu,nl(l),ibc(iad1),            9d20s23
     $                       iddc1x(lsa,lsb,lsc,l),                     9d20s23
     $                  njhere,nfdat(2,l,isb),icoldc,nxx,itmpdc1,       9d25s23
     $                  itmpdc2,i11s,nnt,ncsf2(l,iarg),ncsf(jarg),      9d20s23
     $            bc(jprod),ldcont,iad2,sumdc(l),fctr,bc,ibc,igqqq,'e')     4d29s24
                        if(lchoice)then                                  3d19s21
                         call dgemm('n','n',ncsf(jarg),nl(l),            3d19s21
     $                      ncsf2(l,iarg),1d0,bc(jprod),ncsf(jarg),     3d19s21
     $                      bc(iad2),ncsf2(l,iarg),0d0,bc(itmp),        3d19s21
     $                      ncsf(jarg),                                 3d19s21
     d' hcds.  2')
                        else                                             3d19s21
                         call genmatn2(ncsf(jarg),ncsf(karg),ncsf(iarg),  2d25s21
     $                  ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,      2d25s21
     $                  bc(iad2),ncsf2(l,iarg),nl(l),bc(itmp),0d0,      2d25s21
     $                  iargo,ncsf2(l,iarg),1d0,bc,ibc)                 11d10s22
                        end if                                           3d19s21
                        jden=id1vnotv(l,lsa,lsb)+nnl*icol                     12d18s20
                        ibc(nd1vnotv(l,lsa,lsb)+icol)=1                       12d18s20
                        jtmp=itmp                                              12d18s20
                        mden=mdhvnotv(l)                                 12d21s20
                        do iii=0,nl(l)-1                                         12d18s20
                         ii=ibc(iad1+iii)-1                                      12d18s20
                         jdh=jden+ncsf(jarg)*ii                             12d18s20
                         ibc(mden+ii)=1                                  12d19s20
                         do j=0,ncsf(jarg)-1                                   12d18s20
                          bc(jdh+j)=bc(jdh+j)+bc(jtmp+j)                       12d18s20
                         end do                                                12d18s20
                         jtmp=jtmp+ncsf(jarg)                                  12d18s20
                        end do                                                 12d18s20
                        ibcoff=itmp                                      12d19s20
                        ldcont=ldcont+nrootu*nl(l)*ncsf2(l,iarg)        9d20s23
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
     $         ncsfmid1,bc,ibc)                                         12d10s22
                   call gandc(itestc,itesto,j2c,j2o,nopenk,nopen2p,         12d8s20
     $         karg,iarg,ncsf,norbxx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,    11d13s20
     $         ncsfmid2,bc,ibc)                                         12d10s22
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
                    if(lsa.eq.lsb)then                                    12d19s20
                     nn=(irefo(lsa)*(irefo(lsa)+1))/2                     12d19s20
                     icol=((lgb*(lgb+1))/2)+lga+nn*lgc                    12d19s20
                     icoldc=icol                                        9d25s23
                    else                                                  12d19s20
                     nn=irefo(lsa)*irefo(lsb)                           9d20s23
                     icol=lga+irefo(lsa)*(lgb+irefo(lsb)*lgc)             12d19s20
                     icase=invk1(2,lsa,lsb,lsc,2)                       9d25s23
                     icoldc=icol                                        9d25s23
                     if(icase.eq.2)then                                 9d25s23
                      icoldc=lgb+irefo(lsb)*(lga+irefo(lsa)*lgc)        9d25s23
                     end if                                             9d25s23
                    end if                                                12d19s20
                    nxx=nn*irefo(lsc)                                   9d20s23
                    iprod=ibcoff                                          3d19s21
                    if(lchoice)then                                       3d19s21
                     call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),        3d19s21
     $                 ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod, 11d10s22
     $                   bc,ibc)                                        11d10s22
                    end if                                                3d19s21
                    iargo=1                                               2d25s21
                    jprod=iprod                                           3d19s21
                    ldcont=jdcont                                       9d20s23
                    do l=1,4
                     if(nl(l).gt.0)then                                     12d18s20
                      nnl=ncsf(jarg)*nfdat(2,l,isb)                         12d19s20
                      iad1=jvcv+ibc(jvcv+l+5)                             3d19s21
                      iad2=iad1+nl(l)                                     3d19s21
                      itmp=ibcoff                                           12d18s20
                      ibcoff=itmp+ncsf(jarg)*nl(l)                          12d18s20
                      itmpdc1=ibcoff                                       9d18s23
                      itmpdc2=itmpdc1+ncsf(jarg)*nl(l)*nrootu              9d18s23
                      ibcoff=itmpdc2+ncsf2(l,iarg)*nl(l)*nrootu            9d18s23
                      call enough('hcds.  8',bc,ibc)
                      nnt=nl(l)*nrootu                                     9d18s23
c     nclo1p,isb,if1,isbv1,nclo2p,if2,ipss,l
                      call hcds12c(nrootu,nl(l),ibc(iad1),              9d20s23
     $                       iddc1x(lsa,lsb,lsc,l),                     9d20s23
     $                  njhere,nfdat(2,l,isb),icoldc,nxx,itmpdc1,       9d25s23
     $                  itmpdc2,i11s,nnt,ncsf2(l,iarg),ncsf(jarg),      9d20s23
     $             bc(jprod),ldcont,iad2,sumdc(l),fctr,bc,ibc,igqqq,'f')     4d29s24
                      if(lchoice)then                                     3d19s21
                       call dgemm('n','n',ncsf(jarg),nl(l),             10d26s22
     $                     ncsf2(l,iarg),1d0,bc(jprod),ncsf(jarg),      10d26s22
     $                     bc(iad2),ncsf2(l,iarg),0d0,bc(itmp),         10d26s22
     $                     ncsf(jarg),                                  10d26s22
     d' hcds.  3')
                      else                                                3d19s21
                       call genmatn2(ncsf(jarg),ncsf(karg),ncsf(iarg),     2d25s21
     $                 ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,       2d25s21
     $                 bc(iad2),ncsf2(l,iarg),nl(l),bc(itmp),0d0,       2d25s21
     $                 iargo,ncsf2(l,iarg),1d0,bc,ibc)                  11d10s22
                      end if                                              3d19s21
                      jden=id1vnotv(l,lsa,lsb)+nnl*icol                     12d18s20
                      ibc(nd1vnotv(l,lsa,lsb)+icol)=1                       12d18s20
                      jtmp=itmp                                              12d18s20
                      mden=mdhvnotv(l)                                    12d21s20
                      do iii=0,nl(l)-1                                         12d18s20
                       ii=ibc(iad1+iii)-1                                      12d18s20
                       jdh=jden+ncsf(jarg)*ii                             12d18s20
                       ibc(mden+ii)=1                                     12d19s20
                       do j=0,ncsf(jarg)-1                                   12d18s20
                        bc(jdh+j)=bc(jdh+j)+bc(jtmp+j)                       12d18s20
                       end do                                                12d18s20
                       jtmp=jtmp+ncsf(jarg)                                  12d18s20
                      end do                                                 12d18s20
                      ibcoff=itmp                                         12d19s20
                      ldcont=ldcont+nrootu*nl(l)*ncsf2(l,iarg)          9d20s23
                     end if                                               12d19s20
                     iargo=iargo+ncsf2(l,iarg)                            2d25s21
                     jprod=jprod+ncsf2(l,iarg)*ncsf(jarg)                 3d19s21
                    end do                                                12d19s20
                    ibcoff=iprod                                          3d19s21
                    if(ipss.eq.2)go to 3                                  12d18s20
                   end if                                                 12d18s20
                  end if                                                 4d20s21
                 end do                                                  12d18s20
    3            continue                                                12d18s20
                end if                                                   12d18s20
               end if                                                   10d13s22
               j2o=ibclr(j2o,norbx)                                      12d18s20
               j2o=ibset(j2o,norbxxx)                                     12d18s20
               gandco=ieor(i1o,j2o)                                     10d13s22
               gandcb=ior(gandcc,gandco)                                10d21s22
               ndifb=popcnt(gandcb)                                     10d21s22
               if(ndifb.le.4)then                                       10d21s22
                ndifs=popcnt(gandco)                                     10d13s22
                ndifd=popcnt(gandcc)                                     10d13s22
                nnot=0                                                  10d21s22
                if(ndifs.eq.4.and.ndifb.eq.4)then                                           10d14s22
                 nnot=4                                                          10d14s22
                 ioxx(1)=1                                                     10d17s22
                 ioxx(2)=1                                                     10d17s22
                 do i=1,norbxxx                                                     10d17s22
                  if(btest(gandcb,i))then                                         10d14s22
                   if((btest(j2c,i).and.btest(i1o,i)).or.                         10d17s22
     $      (btest(j2o,i).and..not.btest(i1c,i)))then                   10d14s22
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
     $                ((btest(i1c,i).and..not.btest(j2o,i)).or.                 10d14s22
     $         (btest(j2c,i).and..not.btest(i1o,i))))then                     10d14s22
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
                  if(btest(i1c,nab4(1,2)).and.                          10d26s22
     $                  .not.btest(i1c,nab4(1,1)))nbt=1                 10d26s22
                 else                                                             10d17s22
                  nbt=0                                                           10d17s22
                  if(btest(j2c,nab4(2,2)).and.
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
                ipssx=0                                                 10d26s22
                if(nnot.eq.3)then                                          12d8s20
                 ipssx=1                                                   12d8s20
                else if(nnot.eq.4)then                                  10d26s22
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
     $         ncsfmid1,bc,ibc)                                         12d10s22
                  call gandc(itestc,itesto,j2c,j2o,nopenk,nopen2p,         12d8s20
     $         karg,iarg,ncsf,norbxxx,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,    11d13s20
     $         ncsfmid2,bc,ibc)                                         12d10s22
                  if(nnot1.eq.2.and.nnot2.eq.2)then                       12d19s20
                   if(nab1(2).eq.norbxx)then                              1d25s21
                    iprod=ibcoff                                         3d19s21
                    if(lchoice)then                                      3d19s21
                     call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),       3d19s21
     $               ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod,bc,11d10s22
     $                   ibc)                                           11d10s22
                    end if                                               3d19s21
                    lsa=ism(nab1(1))                                       12d19s20
                    lga=irel(nab1(1))-1                                    12d19s20
                    idx=1                                                 12d19s20
                    jprod=iprod                                          3d19s21
                    iargo=1                                               2d25s21
                    ldcont=jdcont                                       9d26s23
                    do l=1,4                                               12d19s20
                     if(nl(l).gt.0)then                                     12d18s20
                      nnl=ncsf(jarg)*nfdat(2,l,isb)                         12d19s20
                      iad1=jvcv+ibc(jvcv+5+l)                            3d19s21
                      iad2=iad1+nl(l)                                    3d19s21
                      itmp=ibcoff                                           12d18s20
                      ibcoff=itmp+ncsf(jarg)*nl(l)                          12d18s20
                      itmpdc1=ibcoff                                       9d18s23
                      itmpdc2=itmpdc1+ncsf(jarg)*nl(l)*nrootu              9d18s23
                      ibcoff=itmpdc2+ncsf2(l,iarg)*nl(l)*nrootu            9d18s23
                      call enough('hcds12.tmpdc13x',bc,ibcoff)               9d18s23
                      do iz=itmpdc1,itmpdc2-1                              9d18s23
                       bc(iz)=0d0                                          9d18s23
                      end do                                               9d18s23
                      nnt=nl(l)*nrootu                                     9d18s23
                      call hcds12c(nrootu,nl(l),ibc(iad1),iddc3x(l),    9d20s23
     $                  njhere,nfdat(2,l,isb),lga,irefo(lsa),itmpdc1,   9d20s23
     $                  itmpdc2,i11s,nnt,ncsf2(l,iarg),ncsf(jarg),      9d20s23
     $             bc(jprod),ldcont,iad2,sumdc(l),fctr,bc,ibc,igqqq,'g')     4d29s24
                      ldcont=ldcont+nrootu*nl(l)*ncsf2(l,iarg)          9d20s23
                      if(lchoice)then                                    3d19s21
                       call dgemm('n','n',ncsf(jarg),nl(l),             4d20s21
     $                     ncsf2(l,iarg),                               4d20s21
     $              1d0,bc(jprod),ncsf(jarg),bc(iad2),ncsf2(l,iarg),0d0,12d19s20
     $                 bc(itmp),ncsf(jarg),                             12d18s20
     d' hcds.  4')
                      else                                               3d19s21
                       call genmatn2(ncsf(jarg),ncsf(karg),ncsf(iarg),          12d18s20
     $                      ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,         2d25s21
     $                   bc(iad2),ncsf2(l,iarg),nl(l),bc(itmp),0d0,     2d25s21
     $                   iargo,ncsf2(l,iarg),1d0,bc,ibc)                11d10s22
                      end if                                             3d19s21
                      jden=id3vnotv3(l,idx)+nnl*lga                        12d22s20
                      ibc(nd3vnotv3(l,idx)+lga)=1                          12d22s20
                      mden=md3vnotv3(l,idx)                                12d22s20
                      jtmp=itmp                                              12d18s20
                      do iii=0,nl(l)-1                                         12d18s20
                       ii=ibc(iad1+iii)-1                                      12d18s20
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
    4           continue                                                 12d18s20
               end if                                                    10d13s22
              end if                                                    3d2s21
              jvcv=jvcv+nspace                                          3d19s21
              jdcont=jdcont+ndcont                                      9d18s23
c     end if2 loop
             end do                                                     12d18s20
            end if                                                      10d27s22
c     end nclo2p loop
           end do                                                       10d27s22
c
           if(njhere.gt.0)then                                          3d2s21
            do l=1,4
             if(nfdat(2,l,isb).gt.0)then                                12d19s20
              n1=ncsf(jarg)*nfdat(2,l,isb)
              nok4f(l)=0                                                12d21s20
              do i=0,nfdat(2,l,isb)-1                                   12d21s20
               if(ibc(mdhvnotv(l)+i).ne.0)then                          12d21s20
                ibc(mdhvnotv(l)+nok4f(l))=i                             12d21s20
                nok4f(l)=nok4f(l)+1                                     12d21s20
               end if                                                   12d21s20
              end do                                                    12d21s20
              if(nok4f(l).gt.0)then                                     12d21s20
               nnl=ncsf(jarg)*nfdat(2,l,isb)                             12d19s20
               lgoal=multh(jsbv,isbv12)                                 12d21s20
               nok4(l)=0                                                12d21s20
               jdhvnotv=idhvnotv(l)                                     12d21s20
               jdhvnotvf=idhvnotvf(l)                                   1d6s21
               do i=0,irefo(lgoal)-1                                    12d21s20
                if(ibc(ndhvnotv(l)+i).ne.0)then                         12d21s20
                 ibc(ndhvnotv(l)+nok4(l))=i                             12d21s20
                 nok4(l)=nok4(l)+1                                      12d21s20
                 iad=idhvnotv(l)+nnl*i
                 sz=0d0                                                 10d27s22
                 jdhvnotvf0=jdhvnotvf                                   10d27s22
                 do k=0,nok4f(l)-1                                      12d21s20
                  iad=idhvnotv(l)+ncsf(jarg)*ibc(mdhvnotv(l)+k)+nnl*i   1d6s21
                  do j=0,ncsf(jarg)-1                                   1d6s21
                   bc(jdhvnotvf+j)=bc(iad+j)                            1d6s21
                   sz=sz+bc(iad+j)**2                                   10d27s22
                  end do                                                12d21s20
                  jdhvnotvf=jdhvnotvf+ncsf(jarg)                        1d6s21
                  iad=iad+i11s-1                                        1d6s21
                  do j=0,njhere-1                                       12d21s20
                   bc(jdhvnotv+j)=bc(iad+j)                             12d21s20
                  end do                                                12d21s20
                  jdhvnotv=jdhvnotv+njhere                              12d21s20
                 end do                                                 12d21s20
                 sz=sqrt(sz/dfloat(ncsf(jarg)*nok4f(l)))                10d27s22
                 if(sz.lt.1d-14)then                                    10d27s22
                  nok4(l)=nok4(l)-1                                     10d27s22
                  jdhvnotvf=jdhvnotvf0                                  10d27s22
                 end if                                                 10d27s22
   10            format(i5,5x,20i2)
                end if
               end do
               if(nok4(l).gt.0)then
                 ichvnotv(l)=idhvnotv(l)                                11d9s22
               end if
               do isc=1,nsymb                                           12d21s20
                iscv=multh(isc,lgoal)                                   12d21s20
                do isd=1,isc                                               12d18s20
                 iscdv=multh(iscv,isd)                                     12d18s20
                 if(isc.eq.isd)then                                        12d18s20
                  nn=(irefo(isd)*(irefo(isd)+1))/2                         12d18s20
                 else                                                      12d18s20
                  nn=irefo(isd)*irefo(isc)                                 12d18s20
                 end if                                                    12d18s20
                 nnn=nn*irefo(iscdv)                                       12d18s20
                 nokdc(isd,isc,l)=0                                     12d21s20
                 jd1vnotv=id1vnotv(l,isd,isc)                           12d21s20
                 icase=invk1(2,isd,isc,iscdv,2)                         12d29s20
                 do i=0,nnn-1
                  if(ibc(nd1vnotv(l,isd,isc)+i).ne.0)then
                   if(icase.eq.1)then                                   12d29s20
                    icol=i                                              12d29s20
                   else                                                 12d29s20
                    i3=i/nn                                             12d29s20
                    left=i-nn*i3                                        12d29s20
                    ic=left/irefo(isd)                                  12d29s20
                    id=left-ic*irefo(isd)                               12d29s20
                    icol=ic+irefo(isc)*(id+irefo(isd)*i3)               12d29s20
                   end if                                               12d29s20
                   ibc(nd1vnotv(l,isd,isc)+nokdc(isd,isc,l))=icol       12d29s20
                   jd1vnotv0=jd1vnotv                                   10d27s22
                   sz=0d0                                               10d27s22
                   do k=0,nok4f(l)-1                                    12d21s20
                    iad=id1vnotv(l,isd,isc)+i11s-1                      12d21s20
     $                   +ncsf(jarg)*ibc(mdhvnotv(l)+k)+nnl*i           12d21s20
                    do j=0,njhere-1                                     12d21s20
                     bc(jd1vnotv+j)=bc(iad+j)                           12d21s20
                     sz=sz+bc(iad+j)**2                                 10d27s22
                    end do                                              12d21s20
                    jd1vnotv=jd1vnotv+njhere                            12d21s20
                   end do                                               12d21s20
                   sz=sqrt(sz/dfloat(nok4f(l)*njhere))                  10d27s22
                   if(sz.le.1d-14)then                                  10d27s22
                    jd1vnotv=jd1vnotv0                                  10d27s22
                   else                                                 10d27s22
                    nokdc(isd,isc,l)=nokdc(isd,isc,l)+1                  12d21s20
                   end if                                               10d27s22
                  end if
                 end do
                 if(nokdc(isd,isc,l).gt.0)then
                   ic1vnotv(l,isd,isc)=id1vnotv(l,isd,isc)              11d9s22
                 end if
                end do
               end do                                                   12d21s20
              end if
             end if                                                     12d21s20
            end do                                                      12d21s20
           end if                                                       3d2s21
           if(itransgg.eq.0)then                                        1d30s21
            do i=0,ncolt-1                                              1d30s21
             bc(igg+i)=0d0                                              1d30s21
            end do                                                      1d30s21
            itransgg=1                                                  1d30s21
           end if                                                       1d30s21
           if(njhere.gt.0)then                                          3d2s21
            ipass=1                                                     3d2s21
            jdkeep=idkeep(ipass)                                        1d4s21
            nok=0                                                       1d4s21
            do i=0,irefo(lgoal3)-1                                       1d4s21
             do l=1,4                                                   1d4s21
              if(nfdat(2,l,isb).gt.0)then                               1d4s21
               if(ibc(nd3vnotv3(l,ipass)+i).ne.0)then                   1d4s21
                do k=0,nfdat(2,l,isb)-1                                 1d4s21
                 if(ibc(md3vnotv3(l,ipass)+k).ne.0)then                 1d4s21
                  ibc(ndkeep(ipass)+nok)=i                              1d5s21
                  ibc(mdkeep(ipass)+nok)=k+loff(l)                      1d5s21
                  iad=id3vnotv3(l,ipass)+ncsf(jarg)*(k                  1d4s21
     $                   +nfdat(2,l,isb)*i)                             1d4s21
     $                   +i11s-1                                        3d2s21
                  ibc(ibcoff+nok)=l
                  do j=0,njhere-1                                       3d2s21
                   bc(jdkeep+j)=bc(iad+j)                               1d4s21
                  end do                                                1d4s21
                  jdkeep=jdkeep+njhere                                  3d2s21
                  nok=nok+1                                             1d4s21
                 end if                                                 1d4s21
                end do                                                  1d4s21
               end if                                                   1d4s21
              end if                                                    1d4s21
             end do                                                     1d4s21
            end do                                                       1d4s21
            keep(ipass)=nok                                             1d4s21
            if(nok.gt.0)then                                            3d2s21
             nn=nvirt(jsbv)*nrootu                                      2d3s21
             itmpb=ibcoff                                               2d3s21
             itmpd=itmpb+nn*nok                                         2d3s21
             ibcoff=itmpd+njhere*nok                                    3d2s21
             call enough('hcds. 11',bc,ibc)
             jvss=jvs+nn*(i11s-1)                                       3d2s21
             call dgemm('n','n',nn,nok,njhere,1d0,                      3d2s21
     $               bc(jvss),nn,bc(idkeep(ipass)),njhere,0d0,          3d2s21
     $               bc(itmpb),nn,                                      2d4s21
     d' hcds.  5')
             jtmpb=itmpb                                                2d3s21
             do i=0,nok-1                                               2d3s21
              do j=0,njhere-1                                           3d2s21
               ji=idkeep(ipass)+j+njhere*i                              3d2s21
               ij=itmpd+i+nok*j                                         2d3s21
               bc(ij)=bc(ji)                                            2d3s21
              end do                                                    2d3s21
              in=ibc(ndkeep(ipass)+i)                                   2d3s21
              k=ibc(mdkeep(ipass)+i)                                    2d3s21
              jvmat=ivmat(isb)+nn*(k+nfh*in)                            2d4s21
              do j=0,nn-1                                               2d3s21
               bc(jvmat+j)=bc(jvmat+j)+bc(jtmpb+j)                      2d4s21
              end do                                                    2d3s21
              jtmpb=jtmpb+nn                                            2d3s21
             end do                                                     2d3s21
             jggs=jvss+igg-ivs                                          2d3s21
             ibcoff=itmpb                                               2d3s21
            end if                                                      2d3s21
            itmpsvd=ibcoff                                              9d13s23
            ibcoff=itmpsvd+nvirt(jsbv)*njhere*nrootu                    9d13s23
            call enough('hcds12.tmpsvd',bc,ibc)                         9d13s23
            do ir=0,nrootu-1                                            9d13s23
             do j=0,njhere-1                                            9d13s23
              iad1=itmpsvd+nvirt(jsbv)*(j+njhere*ir)                    9d13s23
              iad2=jvs+nvirt(jsbv)*(ir+nrootu*(j+i11s-1))               9d13s23
              do iv=0,nvirt(jsbv)-1                                     9d13s23
               bc(iad1+iv)=bc(iad2+iv)                                  9d13s23
              end do                                                    9d13s23
             end do                                                     9d13s23
            end do                                                      9d13s23
            sumcc=sumcc+sumdc(1)+sumdc(2)+sumdc(3)+sumdc(4)
            do isbv1=1,nsymb                                            12d21s20
             isbv2=multh(isbv1,isbv12)                                  12d21s20
             if(isbv1.le.isbv2)then                                     12d21s20
              if(isbv12.eq.1.and.jsbv.eq.isbv1.and.                     11d15s22
     $            nfdat(2,1,isb).gt.0)then                              11d15s22
               nokf=nok4f(1)                                            2d25s21
               nrow=njhere*nokf                                          12d21s20
               intden=ibcoff                                             12d21s20
               ibcoff=intden+nrow*nvirt(jsbv)                            12d21s20
               call enough('hcds12. 12',bc,ibc)
               fact=0d0                                                  12d21s20
               if(nok4f(1).gt.0)then                                    9d18s23
                nok=nok4(1)                                             2d25s21
                nokf=nok4f(1)                                           2d25s21
                itmp=ibcoff                                             9d15s23
                ivvtmp=itmp+nok*nvirt(jsbv)                             9d15s23
                ibcoff=ivvtmp+nvirt(jsbv)*nrootu*njhere*nokf              9d15s23
                call enough('hcds12. 13',bc,ibc)
                do ir=0,nrootu-1                                        9d15s23
                 do k=0,nokf-1                                           9d15s23
                  iadvd=joffdnon+nvirt(jsbv)*(ir                        12d21s20
     $                   +nrootu*ibc(mdhvnotv(1)+k))                    9d15s23
                  do j=0,njhere-1                                        9d15s23
                   jvvtmp=ivvtmp+nvirt(jsbv)*(j+njhere*(k+nokf*ir))     9d15s23
                   iads=jvs+nvirt(jsbv)*(ir+nrootu*(i11s+j-1))          1d28s21
                   do iv=0,nvirt(jsbv)-1                                 9d15s23
                    bc(jvvtmp+iv)=bc(iads+iv)*vd(iadvd+iv)              9d15s23
                   end do                                               9d15s23
                  end do                                                9d15s23
                 end do                                                 9d15s23
                end do                                                  9d15s23
                nmul=njhere*nokf                                        9d15s23
                if(nok.gt.0)then                                        9d15s23
                 do ir=0,nrootu-1                                       9d15s23
                  jvvtmp=ivvtmp+nvirt(jsbv)*njhere*nokf*ir              9d15s23
                  call dgemm('n','n',nvirt(jsbv),nok,nmul,sr2,          9d15s23
     $                 bc(jvvtmp),nvirt(jsbv),bc(ichvnotv(1)),nmul,0d0, 9d15s23
     $                 bc(itmp),nvirt(jsbv),'hcds12.tmp')               9d15s23
                  do i=0,nok-1                                          9d15s23
                   ii=ibc(ndhvnotv(1)+i)                                9d15s23
                   jtmp=itmp+nvirt(jsbv)*i                              9d15s23
                   jd=idenh0(jsbv)+irefo(jsbv)+nh0av(jsbv)*(ii+         9d15s23
     $                  nh0av(jsbv)*ir)                                 9d15s23
                   do iv=0,nvirt(jsbv)-1                                9d15s23
                    bc(jd+iv)=bc(jd+iv)+bc(jtmp+iv)                     9d15s23
                   end do                                               9d15s23
                  end do                                                9d15s23
                 end do                                                 9d15s23
                end if                                                  9d15s23
                do isc=1,nsymb                                          9d15s23
                 iscv=multh(isc,lgoal)                                  9d15s23
                 do isd=1,isc                                           9d15s23
                  iscdv=multh(iscv,isd)                                 9d15s23
                  if(nokdc(isd,isc,1).gt.0)then                         9d15s23
                   if(isd.eq.isc)then                                   12d22s20
                    nn=(irefo(isc)*(irefo(isc)+1))/2                    12d22s20
                   else                                                 12d22s20
                    nn=irefo(isc)*irefo(isd)                            12d22s20
                   end if                                               12d22s20
                   nnn=nn*irefo(iscdv)                                  12d22s20
                   i2eu=invk1(1,isd,isc,iscdv,2)                        12d22s20
                   itmpp=ibcoff                                          12d22s20
                   ibcoff=itmpp+nokdc(isd,isc,1)*nvirt(jsbv)            12d22s20
                   call enough('hcds12. 18',bc,ibc)
c     vvtmp is v,j,k,root
                   do ir=0,nrootu-1                                     9d15s23
                    jvvtmp=ivvtmp+nvirt(jsbv)*nmul*ir                   9d15s23
                    call dgemm('n','n',nvirt(jsbv),nokdc(isd,isc,1),    9d15s23
     $                 nmul,sr2,bc(jvvtmp),nvirt(jsbv),                 9d15s23
     $                 bc(ic1vnotv(1,isd,isc)),nmul,                    9d15s23
     $                 0d0,bc(itmpp),nvirt(jsbv),'hcds12.tmp')          9d15s23
                    jj=id1x(i2eu)+nnn*nvirt(jsbv)*ir                    9d15s23
                    jtmpp=itmpp                                         9d15s23
                    do i=0,nokdc(isd,isc,1)-1                           9d15s23
                     jjj=jj+nvirt(jsbv)*ibc(nd1vnotv(1,isd,isc)+i)      9d15s23
                     ixx=ionext(i2eu)
     $                  +nvirt(jsbv)*ibc(nd1vnotv(1,isd,isc)+i)
                     do iv=0,nvirt(jsbv)-1                              9d15s23
                      bc(jjj+iv)=bc(jjj+iv)+bc(jtmpp+iv)                9d15s23
                     end do                                             9d15s23
                     jtmpp=jtmpp+nvirt(jsbv)                            9d15s23
                    end do                                              9d15s23
                   end do                                               9d15s23
                   ibcoff=itmpp                                         9d15s23
                  end if                                                9d15s23
                 end do                                                 9d15s23
                end do                                                  9d15s23
                ibcoff=ivvtmp                                           9d15s23
                if(nok.gt.0)then                                        9d18s23
                 jtmp=itmp                                                12d21s20
c
c     visv case
c     ichvnotv=idhvnotv:ncsf,ncont,irefo (or compressed)
c     ichvnotv_*h0av=intden:ncsf,ncont,v
c     trans=intdenT=v,ncont,ncsf
c     gd=itransVs: gd_v,root,ncont=sum_ncsf v,cont,ncsf*v,root,ncsf
c     so den_=sum_v,
                 do iv=0,nvirt(jsbv)-1                                    12d21s20
                  ivp=iv+irefo(jsbv)                                      12d21s20
                  do i=0,nok-1                                            12d21s20
                   iadh=ih0av(jsbv)+ibc(ndhvnotv(1)+i)+nh0av(jsbv)*ivp   2d25s21
                   bc(jtmp+i)=bc(iadh)                                    12d21s20
                  end do                                                  12d21s20
                  jtmp=jtmp+nok                                           12d21s20
                 end do                                                   12d21s20
                 if(ichvnotv(1).lt.0)then                                11d8s22
                  nusedi=ibc(-ichvnotv(1))/2                             11d8s22
                  if(nusedi*2.ne.ibc(-ichvnotv(1)))nusedi=nusedi+1       11d8s22
                  if(nusedi.gt.0)then                                    11d8s22
                   icmp1=1-ichvnotv(1)                                   11d8s22
                   icmp2=icmp1+nok                                       11d9s22
                   icmp3=icmp2+nusedi                                    11d8s22
                   call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),      11d8s22
     $                 bc(itmp),nok,nok,bc(intden),nrow,                11d9s22
     $                 nrow,nvirt(jsbv),fact,sr2)
                   ncomp=ncomp+1
                  end if                                                 11d8s22
                 else                                                    11d8s22
                  call dgemm('n','n',nrow,nvirt(jsbv),nok,sr2,           11d9s22
     $                 bc(ichvnotv(1)),nrow,bc(itmp),nok,fact,           11d9s22
     $                bc(intden),nrow,'hcds12.intden')                  9d18s23
                  nucomp=nucomp+1
                 end if                                                  11d8s22
                 fact=1d0                                                 12d21s20
                end if                                                  9d18s23
                ibcoff=itmp                                              12d21s20
               end if                                                    12d21s20
               nok=0                                                     12d21s20
               do isc=1,nsymb
                iscv=multh(isc,jsbv)                                       12d18s20
                do isd=1,isc                                              12d19s20
                 iscdv=multh(iscv,isd)                                    12d19s20
                 if(isc.eq.isd)then                                        12d18s20
                  nn=(irefo(isd)*(irefo(isd)+1))/2                         12d18s20
                 else                                                      12d18s20
                  nn=irefo(isd)*irefo(isc)                                 12d18s20
                 end if                                                    12d18s20
                 nnn=nn*irefo(iscdv)                                       12d18s20
                 if(min(nokdc(isd,isc,1),nok4f(1)).gt.0)then            2d25s21
                  nok=nokdc(isd,isc,1)                                  2d25s21
                  nokf=nok4f(1)                                         2d25s21
                  ncol=nok*nokf                                            12d21s20
                  i2eu=invk1(1,isd,isc,iscdv,2)                          12d21s20
                  icase=invk1(2,isd,isc,iscdv,2)                          12d21s20
                  itmp=ibcoff                                            12d21s20
                  ibcoff=itmp+nok*nvirt(jsbv)                            12d21s20
                  call enough('hcds. 14',bc,ibc)
                  do i=0,nok-1                                          12d21s20
                   jtmp=itmp+i                                          3d23s21
                   iad=ionext(i2eu)                                     3d23s21
     $                    +nvirt(jsbv)*ibc(nd1vnotv(1,isd,isc)+i)       3d23s21
                   do iv=0,nvirt(jsbv)-1                                  12d21s20
                    bc(jtmp+iv*nok)=bc(iad+iv)                          3d23s21
                   end do                                                12d21s20
                  end do                                                 12d21s20
                  if(ic1vnotv(1,isd,isc).lt.0)then                      11d9s22
                   nusedi=ibc(-ic1vnotv(1,isd,isc))/2                   11d9s22
                   if(nusedi*2.ne.ibc(-ic1vnotv(1,isd,isc)))            11d9s22
     $                  nusedi=nusedi+1                                 11d9s22
                   if(nusedi.gt.0)then                                    11d8s22
                    icmp1=1-ic1vnotv(1,isd,isc)                         11d9s22
                    icmp2=icmp1+nok                                       11d9s22
                    icmp3=icmp2+nusedi                                    11d8s22
                    call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),      11d8s22
     $                 bc(itmp),nok,nok,bc(intden),nrow,                11d9s22
     $                 nrow,nvirt(jsbv),fact,sr2)
                    ncomp=ncomp+1
                   end if                                                 11d8s22
                  else                                                    11d8s22
                   call dgemm('n','n',nrow,nvirt(jsbv),nok,sr2,           11d9s22
     $                bc(ic1vnotv(1,isd,isc)),nrow,bc(itmp),nok,fact,           11d9s22
     $                bc(intden),nrow,                                  11d9s22
     d                'hcds.8')
                   nucomp=nucomp+1
                  end if                                                  11d8s22
                  fact=1d0                                                 12d21s20
                  ibcoff=itmp                                              12d21s20
                 end if                                                  12d21s20
                end do                                                   12d21s20
               end do                                                    12d21s20
               if(fact.gt.0.5d0)then                                     12d21s20
                itrans=ibcoff                                            12d21s20
                ibcoff=itrans+nvirt(jsbv)*nrow                           12d21s20
                call enough('hcds. 15',bc,ibc)
                do i=0,nvirt(jsbv)-1                                     12d21s20
                 do j=0,nrow-1                                           12d21s20
                  ji=intden+j+nrow*i                                    12d29s20
                  ij=itrans+i+nvirt(jsbv)*j                              12d21s20
                  bc(ij)=bc(ji)                                          12d21s20
                 end do                                                  12d21s20
                end do                                                   12d21s20
                sumt=0d0
                do if=0,nokf-1                                           12d21s20
                 do ir=0,nrootu-1                                        12d21s20
                  iadvd=joffdnon+nvirt(jsbv)*(ir                        12d21s20
     $                   +nrootu*ibc(mdhvnotv(1)+if))                    2d25s21
                  do j=0,njhere-1                                        12d21s20
                   iad=itrans+nvirt(jsbv)*(j+njhere*if)                  12d21s20
                   iads=jvs+nvirt(jsbv)*(ir+nrootu*(i11s+j-1))          1d28s21
                   iadg=iads+igg-ivs                                     12d21s20
                   do iv=0,nvirt(jsbv)-1                                  12d21s20
                    bc(iadg+iv)=bc(iadg+iv)+bc(iad+iv)*vd(iadvd+iv)     1d27s21
                    prod=fctr*bc(iad+iv)*bc(iads+iv)                    4d29s24
                    gd(iadvd+iv)=gd(iadvd+iv)+prod                      9d19s23
                   end do                                                12d21s20
                  end do                                                 12d21s20
                 end do                                                  12d21s20
                end do                                                   12d21s20
               end if                                                    12d21s20
               ibcoff=intden                                             12d21s20
              end if                                                     12d21s20
              if(isbv12.eq.1)then                                       12d22s20
               joffdnon=joffdnon+nfdat(2,1,isb)*nvirt(isbv1)*nrootu     12d21s20
               nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                    12d21s20
               isw=0                                                    12d21s20
              else                                                      12d21s20
               nvv=nvirt(isbv1)*nvirt(isbv2)                            12d21s20
               isw=1                                                    12d21s20
              end if                                                    12d21s20
              do l=1,4                                                  12d21s20
               if(nfdat(2,l,isb).gt.0)then                              12d21s20
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
                  tf2=-1d0                                              12d29s20
                  factt=tf                                              12d21s20
                 end if                                                 12d21s20
                 if(nok4f(l).gt.0.and.nvirt(isbvu).gt.0)then            3d2s21
c
                  ivvtmp=ibcoff                                         9d13s23
                  itmpdv=ivvtmp+nvirt(isbvu)*nok4f(l)*njhere*nrootu     9d13s23
                  ibcoff=itmpdv+nvv*nok4f(l)                            9d13s23
                  ivvtmp2=ibcoff                                        9d13s23
                  nvn=nvirt(isbvu)*nok4f(l)                             8d7s24
                  ibcoff=ivvtmp2+nvn*njhere                             9d13s23
                  call enough('hcds12.vvtmp',bc,ibc)                    9d13s23
                  if(isbv1.eq.isbv2)then                                9d13s23
                   nmul=njhere*nok4f(l)                                 9d14s23
                   nvnz=nvn*njhere*nrootu                               9d13s23
                   do iz=0,nvnz-1                                       9d13s23
                    bc(ivvtmp+iz)=0d0                                   9d13s23
                   end do                                               9d13s23
                   mrow=nvirt(isbv1)*nok4f(l)                           9d15s23
                   do ir=0,nrootu-1
                    do iz=0,nvirt(isbv2)*mrow-1                         9d15s23
                     bc(itmpdv+iz)=0d0                                  9d15s23
                    end do                                              9d15s23
                    jtmpsvd=itmpsvd+nvirt(jsbv)*njhere*ir               9d13s23
                    do ivd=0,nvirt(isbv2)-2                             9d14s23
                     do k=0,nok4f(l)-1                                    9d13s23
                      kk=ibc(mdhvnotv(l)+k)                               9d13s23
                      iad1=joffdnon+ivd+nvv*(ir+nrootu*kk)              9d14s23
                      iad2=itmpdv+ivd+nvirt(isbv1)*(k+nok4f(l)*(ivd+1)) 9d14s23
                      do iv2=ivd+1,nvirt(isbv2)-1                       9d14s23
                       itri=((iv2*(iv2-1))/2)                           9d14s23
                       bc(iad2)=vd(iad1+itri)*tf2                       9d14s23
                       iad2=iad2+nvn                                    9d14s23
                      end do                                            9d13s23
                     end do                                             9d13s23
                    end do                                              9d13s23
                    jvvtmp=ivvtmp+mrow*njhere*ir                        9d13s23
                    call dgemm('n','n',mrow,njhere,nvirt(jsbv),1d0,     9d13s23
     $                   bc(itmpdv),mrow,bc(jtmpsvd),nvirt(jsbv),0d0,   9d13s23
     $                   bc(jvvtmp),mrow,'hcds12.jvvtmpa')               9d13s23
                    do iz=0,nvirt(isbv2)*mrow-1                         9d15s23
                     bc(itmpdv+iz)=0d0                                  9d15s23
                    end do                                              9d15s23
                    do ivd=0,nvirt(isbv2)-1                             9d14s23
                     do k=0,nok4f(l)-1                                    9d13s23
                      kk=ibc(mdhvnotv(l)+k)                               9d13s23
                      iad1=joffdnon+((ivd*(ivd-1))/2)+nvv*(ir+nrootu*kk)9d13s23
                      iad2=itmpdv+ivd+nvirt(isbv2)*k                    9d14s23
                      do iv1=0,ivd-1                                    9d13s23
                       bc(iad2)=vd(iad1+iv1)                            9d13s23
                       iad2=iad2+nvn                                    9d14s23
                      end do                                            9d13s23
                     end do                                             9d13s23
                    end do                                              9d13s23
                    call dgemm('n','n',mrow,njhere,nvirt(jsbv),1d0,     9d13s23
     $                   bc(itmpdv),mrow,bc(jtmpsvd),nvirt(jsbv),1d0,   9d13s23
     $                   bc(jvvtmp),mrow,'hcds12.jvvtmpb')               9d13s23
                   end do                                               9d13s23
                  else                                                  9d14s23
                   if(isbvu.eq.isbv1)then                               9d14s23
                    do ir=0,nrootu-1
                     jtmpsvd=itmpsvd+nvirt(jsbv)*njhere*ir               9d13s23
                     do k=0,nok4f(l)-1                                    9d13s23
                      kk=ibc(mdhvnotv(l)+k)                               9d13s23
                      do iv2=0,nvirt(isbv2)-1                            9d13s23
                       iad1=joffdnon+nvirt(isbv1)*(iv2+nvirt(isbv2)*(ir  9d13s23
     $                     +nrootu*kk))                                 9d13s23
                       iad2=itmpdv+nvirt(isbv1)*(k+nok4f(l)*iv2)         9d14s23
                       do iv1=0,nvirt(isbv1)-1                           9d13s23
                        bc(iad2+iv1)=vd(iad1+iv1)                          9d13s23
                       end do                                            9d13s23
                      end do                                             9d13s23
                     end do                                              9d13s23
                     mrow=nvirt(isbv1)*nok4f(l)                          9d13s23
                     jvvtmp=ivvtmp+mrow*njhere*ir                        9d13s23
                     call dgemm('n','n',mrow,njhere,nvirt(jsbv),1d0,     9d13s23
     $                   bc(itmpdv),mrow,bc(jtmpsvd),nvirt(jsbv),0d0,   9d13s23
     $                   bc(jvvtmp),mrow,'hcds12.jvvtmpc')               9d13s23
                    end do                                               9d13s23
                   else                                                  9d13s23
                    mrow=nvirt(isbv2)*nok4f(l)                           9d14s23
                    do ir=0,nrootu-1
                     jtmpsvd=itmpsvd+nvirt(jsbv)*njhere*ir               9d13s23
                     do k=0,nok4f(l)-1                                    9d13s23
                      kk=ibc(mdhvnotv(l)+k)                               9d13s23
                      do iv2=0,nvirt(isbv2)-1                            9d13s23
                       iad1=joffdnon+nvirt(isbv1)*(iv2+nvirt(isbv2)*(ir  9d13s23
     $                     +nrootu*kk))                                 9d13s23
                       iad2=itmpdv+iv2+nvirt(isbv2)*k                    9d14s23
                       do iv1=0,nvirt(isbv1)-1                           9d13s23
                        bc(iad2)=vd(iad1+iv1)                            9d13s23
                        iad2=iad2+mrow                                   9d14s23
                       end do                                            9d13s23
                      end do                                             9d13s23
                     end do                                              9d13s23
                     jvvtmp=ivvtmp+mrow*njhere*ir                        9d13s23
                     call dgemm('n','n',mrow,njhere,nvirt(jsbv),1d0,     9d13s23
     $                   bc(itmpdv),mrow,bc(jtmpsvd),nvirt(jsbv),0d0,   9d13s23
     $                   bc(jvvtmp),mrow,'hcds12.jvvtmpd')               9d13s23
                    end do                                               9d13s23
                   end if                                               9d14s23
                  end if                                                9d13s23
                  do ir=0,nrootu-1                                      9d13s23
                   jvvtmp=ivvtmp+nvn*njhere*ir                          9d13s23
                   do k=0,nok4f(l)-1                                     9d13s23
                    do j=0,njhere-1                                      9d13s23
                     jk=ivvtmp2+nvirt(isbvu)*(j+njhere*k)               9d13s23
                     kj=jvvtmp+nvirt(isbvu)*(k+nok4f(l)*j)              9d13s23
                     do iv=0,nvirt(isbvu)-1                              9d13s23
                      bc(jk+iv)=bc(kj+iv)                               9d13s23
                     end do                                             9d13s23
                    end do                                              9d13s23
                   end do                                               9d13s23
                   do k=0,nvn*njhere-1                                  9d13s23
                    bc(jvvtmp+k)=bc(ivvtmp2+k)                          9d13s23
                   end do
                  end do                                                9d13s23
                  nmul=njhere*nok4f(l)                                  9d13s23
                  if(nok4(l).gt.0)then                                  9d13s23
                   ivdprod=ibcoff                                        9d13s23
                   ibcoff=ivdprod+nvirt(isbvu)*nok4(l)                   9d13s23
                   call enough('hcds12.vdprod',bc,ibc)                   9d13s23
                   do ir=0,nrootu-1                                      9d13s23
                    jvvtmp=ivvtmp+nvirt(isbvu)*nmul*ir                   9d13s23
                    call dgemm('n','n',                                  9d13s23
     $                   nvirt(isbvu),nok4(l),nmul,                     9d13s23
     $                   factt,                                         9d13s23
     $                   bc(jvvtmp),nvirt(isbvu),                       9d13s23
     $                   bc(ichvnotv(l)),nmul,                          9d13s23
     $                   0d0,                                           9d13s23
     $                   bc(ivdprod),nvirt(isbvu),                      9d13s23
     $                   'hcds12.vdprod')                               9d13s23
                    do i=0,nok4(l)-1                                     9d13s23
                     ii=ibc(ndhvnotv(l)+i)                               9d13s23
                     iad1=idenh0(isbvu)+irefo(isbvu)+nh0av(isbvu)*(ii    9d13s23
     $                    +nh0av(isbvu)*ir)                              9d13s23
                     iad2=ivdprod+nvirt(isbvu)*i                         9d13s23
                     iad3=ih0av(isbvu)+irefo(isbvu)+nh0av(isbvu)*ii
                     do iv=0,nvirt(isbvu)-1                              9d13s23
                      bc(iad1+iv)=bc(iad1+iv)+bc(iad2+iv)                9d13s23
                     end do                                              9d13s23
                    end do                                               9d13s23
                   end do                                                9d13s23
                  end if                                                9d13s23
                  do isc=1,nsymb                                        12d22s20
                   iscv=multh(isc,lgoal)                                12d22s20
                   do isd=1,isc                                         12d22s20
                    iscdv=multh(iscv,isd)                               12d22s20
                    if(nokdc(isd,isc,l).gt.0)then                       12d22s20
                     if(isd.eq.isc)then                                 12d22s20
                      nn=(irefo(isc)*(irefo(isc)+1))/2                  12d22s20
                     else                                               12d22s20
                      nn=irefo(isc)*irefo(isd)                          12d22s20
                     end if                                             12d22s20
                     nnn=nn*irefo(iscdv)                                12d22s20
                     i2eu=invk1(1,isd,isc,iscdv,2)                      12d22s20
                     itmp=ibcoff                                        12d22s20
                     ibcoff=itmp+nokdc(isd,isc,l)*nvirt(isbvu)          12d22s20
                     call enough('hcds12. 18',bc,ibc)
c     vvtmp is v,j,k,root
                     do ir=0,nrootu-1                                   9d15s23
                      jvvtmp=ivvtmp+nvirt(isbvu)*nmul*ir                9d15s23
                      call dgemm('n','n',nvirt(isbvu),nokdc(isd,isc,l), 9d15s23
     $                     nmul,factt,bc(jvvtmp),nvirt(isbvu),          9d15s23
     $                     bc(ic1vnotv(l,isd,isc)),nmul,                9d15s23
     $                     0d0,bc(itmp),nvirt(isbvu),'hcds12.tmp')      9d15s23
                      jj=id1x(i2eu)+nnn*nvirt(isbvu)*ir                 9d15s23
                      jtmp=itmp                                         9d15s23
                      do i=0,nokdc(isd,isc,l)-1                         9d15s23
                       jjj=jj+nvirt(isbvu)*ibc(nd1vnotv(l,isd,isc)+i)   9d15s23
                       ixx=ionext(i2eu)
     $                      +nvirt(isbvu)*ibc(nd1vnotv(l,isd,isc)+i)
                       do iv=0,nvirt(isbvu)-1                            9d15s23
                        bc(jjj+iv)=bc(jjj+iv)+bc(jtmp+iv)               9d15s23
                       end do                                           9d15s23
                       jtmp=jtmp+nvirt(isbvu)                           9d15s23
                      end do                                            9d15s23
                     end do                                             9d15s23
                     ibcoff=itmp                                        9d15s23
                    end if                                              9d15s23
                   end do                                               9d15s23
                  end do                                                9d15s23
                  ibcoff=ivvtmp                                         9d13s23
                  nrow=njhere*nok4f(l)                                  12d21s20
                  intden=ibcoff                                         12d21s20
                  ibcoff=intden+nrow*nvirt(isbvu)                       12d21s20
                  call enough('hcds. 16',bc,ibc)
                  fact=0d0                                              12d21s20
                  if(nok4(l).gt.0)then                                  1d1s21
c     Gr jv j = D j k i H0iv2 V jv iv2 r k
                   itmp=ibcoff                                          12d21s20
                   ibcoff=itmp+nok4(l)*nvirt(isbvu)                     12d21s20
                   call enough('hcds. 17',bc,ibc)
                   jtmp=itmp                                            12d21s20
                   do iv=0,nvirt(isbvu)-1                               12d21s20
                    ivp=iv+irefo(isbvu)                                 12d21s20
                    ih0=ih0av(isbvu)+nh0av(isbvu)*ivp                   12d21s20
                    do i=0,nok4(l)-1                                    12d21s20
                     ii=ibc(ndhvnotv(l)+i)                              12d22s20
                     bc(jtmp+i)=bc(ih0+ii)                              12d21s20
                    end do                                              12d21s20
                    jtmp=jtmp+nok4(l)                                   12d21s20
                   end do                                               12d21s20
                   if(ichvnotv(l).lt.0)then
                    nusedi=ibc(-ichvnotv(l))/2                             11d8s22
                    if(nusedi*2.ne.ibc(-ichvnotv(l)))nusedi=nusedi+1       11d8s22
                    if(nusedi.gt.0)then                                    11d8s22
                     icmp1=1-ichvnotv(l)                                   11d8s22
                     icmp2=icmp1+nok4(l)                                11d9s22
                     icmp3=icmp2+nusedi                                    11d8s22
                     call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),      11d8s22
     $                 bc(itmp),nok4(l),nok4(l),bc(intden),nrow,        11d9s22
     $                 nrow,nvirt(isbvu),fact,factt)                    11d9s22
                     ncomp=ncomp+1
                    end if                                                 11d8s22
                   else
                    call dgemm('n','n',nrow,nvirt(isbvu),nok4(l),factt,  11d9s22
     $                bc(ichvnotv(l)),nrow,bc(itmp),nok4(l),fact,       11d9s22
     $                bc(intden),nrow,'hcds12.intdenb')                 9d18s23
                    nucomp=nucomp+1
                   end if
                   ibcoff=itmp                                          12d21s20
                   fact=1d0                                             12d21s20
                  end if                                                12d21s20
                  do isc=1,nsymb                                        12d22s20
                   iscv=multh(isc,lgoal)                                12d22s20
                   do isd=1,isc                                         12d22s20
                    iscdv=multh(iscv,isd)                               12d22s20
                    if(nokdc(isd,isc,l).gt.0)then                       12d22s20
                     if(isd.eq.isc)then                                 12d22s20
                      nn=(irefo(isc)*(irefo(isc)+1))/2                  12d22s20
                     else                                               12d22s20
                      nn=irefo(isc)*irefo(isd)                          12d22s20
                     end if                                             12d22s20
                     nnn=nn*irefo(iscdv)                                12d22s20
                     i2eu=invk1(1,isd,isc,iscdv,2)                      12d22s20
                     icase=invk1(2,isd,isc,iscdv,2)                      12d22s20
                     itmp=ibcoff                                        12d22s20
                     ibcoff=itmp+nokdc(isd,isc,l)*nvirt(isbvu)          12d22s20
                     call enough('hcds. 18',bc,ibc)
                     do i=0,nokdc(isd,isc,l)-1                          12d22s20
                      jtmp=itmp+i                                       3d23s21
                      iad=ionext(i2eu)                                  3d23s21
     $                       +nvirt(isbvu)*ibc(nd1vnotv(l,isd,isc)+i)   3d23s21
                      do iv=0,nvirt(isbvu)-1                             12d22s20
                       bc(jtmp+iv*nokdc(isd,isc,l))=bc(iad+iv)          3d23s21
                      end do                                            12d22s20
                     end do                                             12d22s20
                     if(ic1vnotv(l,isd,isc).lt.0)then                      11d9s22
                      nusedi=ibc(-ic1vnotv(l,isd,isc))/2                   11d9s22
                      if(nusedi*2.ne.ibc(-ic1vnotv(l,isd,isc)))            11d9s22
     $                  nusedi=nusedi+1                                 11d9s22
                      if(nusedi.gt.0)then                                    11d8s22
                       icmp1=1-ic1vnotv(l,isd,isc)                         11d9s22
                       icmp2=icmp1+nokdc(isd,isc,l)                     11d9s22
                       icmp3=icmp2+nusedi                                    11d8s22
                       call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),      11d8s22
     $                 bc(itmp),nokdc(isd,isc,l),nokdc(isd,isc,l),
     $                     bc(intden),nrow,nrow,nvirt(isbvu),fact,factt)11d9s22
                       ncomp=ncomp+1
                      end if                                                 11d8s22
                     else                                                    11d8s22
                      call dgemm('n','n',nrow,nvirt(isbvu),             11d9s22
     $                     nokdc(isd,isc,l),factt,                      11d9s22
     $                     bc(ic1vnotv(l,isd,isc)),nrow,bc(itmp),       11d9s22
     $                     nokdc(isd,isc,l),fact,bc(intden),nrow,       11d9s22
     d                'hcds.10')
                      nucomp=nucomp+1
                     end if                                                  11d8s22
                     ibcoff=itmp                                        12d22s20
                     fact=1d0                                           12d22s20
                    end if                                              12d22s20
                   end do                                               12d22s20
                  end do                                                12d22s20
                  if(fact.gt.0.5d0)then                                 12d22s20
                   itrans=ibcoff                                        12d22s20
                   ibcoff=itrans+nrow*nvirt(isbvu)                      12d22s20
                   call enough('hcds. 19',bc,ibc)
                   if(isbv1.eq.isbv2)then                               12d29s20
                    do i=0,nvirt(isbvu)-1                                12d22s20
                     do j=0,nrow-1                                       12d22s20
                      ji=intden+j+nrow*i                                 12d22s20
c     vjk
                      ij=itrans+i+nvirt(isbvu)*j                         12d22s20
                      bc(ij)=bc(ji)                                      12d22s20
                     end do                                              12d22s20
                    end do                                               12d22s20
                   else                                                 12d29s20
                    do iv=0,nvirt(isbvu)-1                              12d29s20
                     do k=0,nok4f(l)-1                                  12d29s20
                      do j=0,njhere-1                                   12d29s20
                       jki=intden+j+njhere*(k+nok4f(l)*iv)              12d29s20
                       kij=itrans+k+nok4f(l)*(iv+nvirt(isbvu)*j)        12d29s20
c     kvj
                       bc(kij)=bc(jki)                                  12d29s20
                      end do                                            12d29s20
                     end do                                             12d29s20
                    end do                                              12d29s20
                   end if                                               12d29s20
                   ncol=nrootu*nvirt(jsbv)                              12d22s20
                   itmpsv=ibcoff                                         12d22s20
                   itmpsg=itmpsv+njhere*ncol                              12d22s20
                   ibcoff=itmpsg+njhere*ncol                             12d22s20
                   call enough('hcds. 20',bc,ibc)
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
                    do i=itmpsg,ibcoff-1                                 12d22s20
                     bc(i)=0d0                                           12d22s20
                    end do                                               12d22s20
                    do ir=0,nrootm                                      1d27s21
                     do j=0,njhere-1                                     12d22s20
                      iad=jvs+nvirt(jsbv)*(ir+nrootu*(j+i11s-1))        1d27s21
                      jad=itmpsv+nvirt(jsbv)*(ir+nrootu*j)              1d27s21
                      do jv=0,nvirt(jsbv)-1                              12d22s20
                       bc(jad+jv)=bc(iad+jv)                            1d28s21
                      end do                                            12d22s20
                     end do                                             12d22s20
                    end do                                              12d22s20
                    sumh=0d0
                    do k=0,nok4f(l)-1                                   12d22s20
                     kk=ibc(mdhvnotv(l)+k)                              12d22s20
                     do ir=0,nrootm                                     10d26s22
                      iadv=joffdnon+nvv*(ir+nrootu*kk)                  9d14s23
                      do j=0,njhere-1                                    12d22s20
                       iadd=itrans+nvirt(isbvu)*(j+njhere*k)            12d22s20
                       iadsv=itmpsv+nvirt(jsbv)*(ir+nrootu*j)            12d22s20
                       iadvu=iadv                                       9d14s23
                       do iv2=0,nvirt(isbv2)-1                          12d22s20
                        do iv1=0,iv2-1                                  12d22s20
                         gd(iadvu+iv1)=gd(iadvu+iv1)                    9d14s23
     $                        +fctr*bc(iadd+iv1)*bc(iadsv+iv2)*tf2      4d29s24
                         gd(iadvu+iv1)=gd(iadvu+iv1)                    9d14s23
     $                        +fctr*bc(iadd+iv2)*bc(iadsv+iv1)          4d29s24
                        end do                                          12d22s20
                        iadvu=iadvu+iv2                                 9d14s23
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
                    call enough('hcds. 21',bc,ibc)
c     vu,j,k
                    call dgemm('n','n',mrow,mcol,njhere,1d0,            3d2s21
     $                     bc(itrans),mrow,bc(itmpsv),njhere,0d0,       3d2s21
     $                     bc(itmpdg),mrow,                             3d2s21
     d' hcds. 12')
                    do ir=0,nrootm                                      10d26s22
                     do k=0,nok4f(l)-1                                  12d22s20
                      kk=ibc(mdhvnotv(l)+k)                             12d22s20
                      do iv2=0,nvirt(isbv2)-1                             12d22s20
                       do iv1=0,nvirt(isbv1)-1                           12d22s20
                        iad1=itmpdg+k+nok4f(l)*(iv1+nvirt(isbv1)*(iv2   1d27s21
     $                        +nvirt(isbv2)*ir))                        1d27s21
                        iad2=joffdnon+iv1+nvirt(isbv1)*(iv2             12d22s20
     $                        +nvirt(isbv2)*(ir+nrootu*kk))             12d22s20
                        gd(iad2)=gd(iad2)+fctr*bc(iad1)                 4d29s24
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
                    call enough('hcds. 22',bc,ibc)
                    call dgemm('n','n',mcol,mrow,njhere,1d0,            3d2s21
     $                     bc(itrans),mcol,bc(itmpsv),njhere,0d0,       3d2s21
     $                     bc(itmpgv),mcol,                             3d2s21
     d' hcds. 14')
                    do ir=0,nrootm                                      1d28s21
                     do k=0,nok4f(l)-1                                  12d22s20
                      kk=ibc(mdhvnotv(l)+k)                             12d22s20
                      do iv1=0,nvirt(isbv1)-1                             12d22s20
                       do iv2=0,nvirt(isbv2)-1                           12d22s20
                        iad1=itmpgv+k+nok4f(l)*(iv2+nvirt(isbv2)        12d22s20
     $                       *(iv1+nvirt(isbv1)*ir))                    1d28s21
                        iad2=joffdnon+iv1+nvirt(isbv1)*(iv2             12d22s20
     $                       +nvirt(isbv2)*(ir+nrootu*kk))              12d22s20
                        gd(iad2)=gd(iad2)+fctr*bc(iad1)                 4d29s24
                       end do                                           12d22s20
                      end do                                            12d22s20
                     end do                                             12d22s20
                    end do                                              12d22s20
                   end if                                               12d22s20
                  end if                                                12d22s20
                  ibcoff=intden                                         12d22s20
                 end if                                                 12d21s20
    5            continue                                               12d21s20
                end if                                                  12d22s20
                joffdnon=joffdnon+nvv*nfdat(2,l,isb)*nrootu             12d22s20
               end if                                                   12d21s20
              end do                                                    12d21s20
             end if                                                     12d21s20
            end do                                                      12d18s20
            ibcoff=itmpsvd                                              9d13s23
           end if                                                       12d18s20
           jvs=jvs+ncsf(jarg)*nvirt(jsbv)*nrootu                        10d27s22
           if(njhere.gt.0)ibcoff=itmp1                                  9d18s23
c     end if1 loop
          end do                                                        12d18s20
          ibcoff=idhvisv                                                10d27s22
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
c     end isb loop
         end do                                                         12d18s20
c
c     igs is nvirt,nroot, ...
c     transpose it to nroot,nvirt,...
c
         itmp=ibcoff                                                    1d27s21
         ibcoff=itmp+ngg                                                1d27s21
         call enough('hcds. 23',bc,ibc)
         jgg=igg                                                        1d27s21
         do i=0,nggg-1                                                  1d27s21
          do ir=0,nrootm                                                1d28s21
           do iv=0,nvirt(jsbv)-1                                         1d27s21
            ji=jgg+iv+nvirt(jsbv)*ir                                    1d27s21
            ij=itmp+ir+nrootu*iv                                        1d27s21
            bc(ij)=bc(ji)                                               1d27s21
           end do                                                       1d27s21
          end do                                                        1d27s21
          do j=0,nggm                                                   1d27s21
           bc(jgg+j)=bc(itmp+j)                                         1d27s21
          end do                                                        1d27s21
          jgg=jgg+ngg                                                   1d27s21
         end do                                                         1d27s21
         ibcoff=itmp                                                    1d27s21
         last8(1)=ihsdiag(nclo1p,jsb,2)                                 2d8s21
         last8(2)=i48                                                   2d8s21
         itransgg=0                                                     1d30s21
         ibcoff=ibcgg                                                   1d30s21
        end if                                                          12d18s20
c     nclo1p
       end do                                                           12d18s20
       ioffdnon=1                                                       2d4s21
       call dws_gsumf(bc(ibcvmat),nbmat)                                3d2s21
       do isb=1,nsymb                                                   2d4s21
        nfh=nfdat(2,1,isb)+nfdat(2,2,isb)+nfdat(2,3,isb)+nfdat(2,4,isb) 2d4s21
        isbv12=multh(isb,isymmrci)                                      2d4s21
        isn=multh(isbv12,jsbv)                                          2d4s21
        if(min(irefo(isn),nvirt(jsbv)).gt.0)then                        2d25s21
         nn=nvirt(jsbv)*nrootu*nfh                                       2d4s21
         ipsr=ibcoff                                                    2d4s21
         ibcoff=ipsr+nrootu*nfh                                         2d4s21
         call enough('hcds. 24',bc,ibc)
         do i=0,nfdat(2,1,isb)*nrootu-1                                 2d4s21
          bc(ipsr+i)=1d0                                                2d4s21
         end do                                                         2d4s21
         do i=nfdat(2,1,isb)*nrootu,nfh*nrootu-1                        2d5s21
          bc(ipsr+i)=-1d0                                               2d4s21
         end do                                                         2d4s21
c
c     form Gv'v"rk=[(vv"|nv')+p(k)(vv'|nv")]vvrkn
c     and densities
c     Dvv"nv'r=sum_k vvrkn*Vv'v"rk and
c     Dvv'nv"r=sum_k vvrkn*Vv'v"rk*p(k)
c
         do isbv1=1,nsymb                                               2d4s21
          isbv2=multh(isbv1,isbv12)                                     2d4s21
          if(isbv2.ge.isbv1)then                                        2d4s21
           call ilimts(irefo(isn),nvirt(isbv1),mynprocg,mynowprog,il,ih,2d4s21
     $         i1s,i1e,i2s,i2e)                                         2d4s21
           nhere=ih+1-il                                                9d26s23
           i2eu=invk1(1,jsbv,isbv2,isn,2)                               2d4s21
           icase=invk1(2,jsbv,isbv2,isn,2)                              2d4s21
           if(isbv1.eq.isbv2)then                                       2d4s21
            i10=i1s                                                     2d4s21
            i1n=irefo(isn)                                              2d4s21
            iint=i3x(i2eu)                                              2d4s21
            iden=id3x(i2eu)                                             9d26s23
            if(jsbv.eq.isbv2)then                                       2d4s21
             nrow=(nvirt(jsbv)*(nvirt(jsbv)+1))/2                       2d4s21
             nrowh=nrow*nhere                                           9d26s23
             do i2=i2s,i2e                                              2d4s21
              iv1=i2-1                                                  2d4s21
              if(i2.eq.i2e)i1n=i1e                                      2d4s21
              do i1=i10,i1n                                             2d4s21
               i1m=i1-1                                                 2d4s21
               do k=0,nfdat(2,1,isb)-1                                  2d4s21
                do ir=0,nrootm                                          2d4s21
                 sum=0d0                                                2d4s21
                 jvmat=ivmat(isb)+nvirt(jsbv)*(ir+nrootu*(k+nfh*i1m))   2d4s21
                 iadv=ioffdnon+iv1+nvirt(isbv1)*(ir+nrootu*k)           2d4s21
                 do jv=0,nvirt(jsbv)-1                                    2d4s21
                  ix=max(jv,iv1)                                          2d4s21
                  in=min(jv,iv1)                                          2d4s21
                  irow=iint+((ix*(ix+1))/2)+in                            2d4s21
                  irowd=iden+((ix*(ix+1))/2)+in+ir*nrowh                9d26s23
                  bc(irowd)=bc(irowd)+bc(jvmat+jv)*vd(iadv)*sr2         9d26s23
                  sum=sum+bc(jvmat+jv)*bc(irow)                         2d4s21
                 end do                                                 2d4s21
                 sum=sum*sr2*fctr                                       4d29s24
                 gd(iadv)=gd(iadv)+sum                                  2d4s21
                end do                                                  2d4s21
               end do                                                   2d4s21
               iint=iint+nrow                                           2d4s21
               iden=iden+nrow                                           9d26s23
              end do                                                    2d4s21
              i10=1                                                     2d4s21
             end do                                                     2d4s21
            else if(icase.eq.1)then                                     2d4s21
             nrow=nvirt(jsbv)*nvirt(isbv1)                              2d4s21
             nrowh=nrow*nhere                                           9d26s23
             do i2=i2s,i2e                                              2d4s21
              iv1=i2-1                                                  2d4s21
              if(i2.eq.i2e)i1n=i1e                                      2d4s21
              do i1=i10,i1n                                             2d4s21
               i1m=i1-1                                                 2d4s21
               do k=0,nfdat(2,1,isb)-1                                  2d4s21
                do ir=0,nrootm                                          2d4s21
                 sum=0d0                                                2d4s21
                 jvmat=ivmat(isb)+nvirt(jsbv)*(ir+nrootu*(k+nfh*i1m))   2d4s21
                 irow=iint+iv1*nvirt(jsbv)                              2d4s21
                 iadv=ioffdnon+iv1+nvirt(isbv1)*(ir+nrootu*k)           2d4s21
                 irowd=iden+iv1*nvirt(jsbv)+nrowh*ir                    9d26s23
                 do jv=0,nvirt(jsbv)-1                                    2d4s21
                  sum=sum+bc(jvmat+jv)*bc(irow+jv)                      2d4s21
                  bc(irowd+jv)=bc(irowd+jv)+bc(jvmat+jv)*vd(iadv)*sr2   9d26s23
                 end do                                                 2d4s21
                 sum=sum*sr2*fctr                                       4d29s24
                 gd(iadv)=gd(iadv)+sum                                  2d4s21
                end do                                                  2d4s21
               end do                                                   2d4s21
               iint=iint+nrow                                           2d4s21
               iden=iden+nrow                                           9d26s23
              end do                                                    2d4s21
              i10=1                                                     2d4s21
             end do                                                     2d4s21
            else                                                        2d4s21
             nrow=nvirt(jsbv)*nvirt(isbv1)                              2d4s21
             nrowh=nrow*nhere                                           9d26s23
             do i2=i2s,i2e                                              2d4s21
              iv1=i2-1                                                  2d4s21
              if(i2.eq.i2e)i1n=i1e                                      2d4s21
              do i1=i10,i1n                                             2d4s21
               i1m=i1-1                                                 2d4s21
               do k=0,nfdat(2,1,isb)-1                                  2d4s21
                do ir=0,nrootm                                          2d4s21
                 sum=0d0                                                2d4s21
                 jvmat=ivmat(isb)+nvirt(jsbv)*(ir+nrootu*(k+nfh*i1m))   2d4s21
                 irow=iint+iv1                                          2d4s21
                 iadv=ioffdnon+iv1+nvirt(isbv1)*(ir+nrootu*k)           2d4s21
                 irowd=iden+iv1+nrowh*ir                                9d26s23
                 do jv=0,nvirt(jsbv)-1                                    2d4s21
                  sum=sum+bc(jvmat+jv)*bc(irow+jv*nvirt(isbv1))         2d4s21
                  bc(irowd+jv*nvirt(isbv1))=bc(irowd+jv*nvirt(isbv1))   9d26s23
     $                 +bc(jvmat+jv)*vd(iadv)*sr2                       9d26s23
                 end do                                                 2d4s21
                 sum=sum*sr2*fctr                                       4d29s24
                 gd(iadv)=gd(iadv)+sum                                  2d4s21
                end do                                                  2d4s21
               end do                                                   2d4s21
               iint=iint+nrow                                           2d4s21
               iden=iden+nrow                                           9d26s23
              end do                                                    2d4s21
              i10=1                                                     2d4s21
             end do                                                     2d4s21
            end if                                                      2d4s21
            ioffdnon=ioffdnon+nvirt(isbv1)*nrootu*nfdat(2,1,isb)        2d4s21
            nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       2d4s21
            isw=0                                                       2d4s21
           else                                                         2d4s21
            nvv=nvirt(isbv1)*nvirt(isbv2)                               2d4s21
            isw=1                                                       2d4s21
           end if                                                       2d4s21
           i10=i1s                                                      2d4s21
           i1n=irefo(isn)                                               2d4s21
           iint=i3x(i2eu)                                               2d4s21
           iden=id3x(i2eu)                                              9d26s23
           if(jsbv.eq.isbv2)then                                        2d4s21
            nrow=(nvirt(jsbv)*(nvirt(jsbv)+1))/2                        2d4s21
            nrowh=nrow*nhere                                            9d26s23
            do i2=i2s,i2e                                               2d4s21
             iv1=i2-1                                                   2d4s21
             ibots=i2                                                   2d4s21
             ibotn=0                                                    2d4s21
             ibot=ibots+isw*(ibotn-ibots)                               2d4s21
             if(i2.eq.i2e)i1n=i1e                                       2d4s21
             do i1=i10,i1n                                              2d4s21
              i1m=i1-1                                                  2d4s21
              do k=0,nfh-1                                              9d26s23
               do ir=0,nrootu-1                                         9d26s23
                kr=ir+nrootu*k                                          9d26s23
                jvmat=ivmat(isb)+nvirt(jsbv)*kr+nn*i1m                   2d4s21
                do iv2=ibot,nvirt(isbv2)-1                                2d4s21
                 sum=0d0                                                 2d4s21
                 itri=((iv2*(iv2-1))/2)+iv1                              2d4s21
                 irec=iv1+nvirt(isbv1)*iv2                               2d4s21
                 iadv=ioffdnon+itri+isw*(irec-itri)+nvv*kr               2d4s21
                 do jv=0,nvirt(jsbv)-1                                    2d4s21
                  ix=max(jv,iv2)                                         2d4s21
                  in=min(jv,iv2)                                         2d4s21
                  irow=iint+((ix*(ix+1))/2)+in                           2d4s21
                  sum=sum+bc(jvmat+jv)*bc(irow)                          2d4s21
                  irowd=iden+((ix*(ix+1))/2)+in+ir*nrowh                 9d26s23
                  bc(irowd)=bc(irowd)+bc(jvmat+jv)*vd(iadv)              9d26s23
                 end do                                                  2d4s21
                 sum=sum*fctr                                           4d29s24
                 gd(iadv)=gd(iadv)+sum                                   2d4s21
                end do                                                  9d26s23
               end do                                                   2d4s21
              end do                                                    2d4s21
              iint=iint+nrow                                            2d4s21
              iden=iden+nrow                                            9d26s23
             end do                                                     2d4s21
             i10=1                                                      2d4s21
            end do                                                      2d4s21
           else if(icase.eq.1)then                                      2d4s21
            nrow=nvirt(jsbv)*nvirt(isbv2)                               2d4s21
            nrowh=nrow*nhere                                            9d26s23
            do i2=i2s,i2e                                               2d4s21
             iv1=i2-1                                                   2d4s21
             ibots=i2                                                   2d4s21
             ibotn=0                                                    2d4s21
             ibot=ibots+isw*(ibotn-ibots)                               2d4s21
             if(i2.eq.i2e)i1n=i1e                                       2d4s21
             do i1=i10,i1n                                              2d4s21
              i1m=i1-1                                                  2d4s21
              do k=0,nfh-1                                              9d26s23
               do ir=0,nrootu-1                                         9d26s23
                kr=ir+nrootu*k                                          9d26s23
                jvmat=ivmat(isb)+nvirt(jsbv)*kr+nn*i1m                   2d4s21
                do iv2=ibot,nvirt(isbv2)-1                                2d4s21
                 sum=0d0                                                 2d4s21
                 irow=iint+nvirt(jsbv)*iv2                               2d4s21
                 itri=((iv2*(iv2-1))/2)+iv1                              2d4s21
                 irec=iv1+nvirt(isbv1)*iv2                               2d4s21
                 iadv=ioffdnon+itri+isw*(irec-itri)+nvv*kr               2d4s21
                 irowd=iden+nvirt(jsbv)*iv2+nrowh*ir                    9d26s23
                 do jv=0,nvirt(jsbv)-1                                    2d4s21
                  sum=sum+bc(jvmat+jv)*bc(irow+jv)                       2d4s21
                  bc(irowd+jv)=bc(irowd+jv)+bc(jvmat+jv)*vd(iadv)       9d26s23
                 end do                                                  2d4s21
                 gd(iadv)=gd(iadv)+sum*fctr                             4d29s24
                end do                                                  9d26s23
               end do                                                   2d4s21
              end do                                                    2d4s21
              iint=iint+nrow                                            2d4s21
              iden=iden+nrow                                            9d26s23
             end do                                                     2d4s21
             i10=1                                                      2d4s21
            end do                                                      2d4s21
           else                                                         2d4s21
            nrow=nvirt(jsbv)*nvirt(isbv2)                               2d4s21
            nrowh=nrow*nhere                                            9d26s23
            do i2=i2s,i2e                                               2d4s21
             iv1=i2-1                                                   2d4s21
             ibots=i2                                                   2d4s21
             ibotn=0                                                    2d4s21
             ibot=ibots+isw*(ibotn-ibots)                               2d4s21
             if(i2.eq.i2e)i1n=i1e                                       2d4s21
             do i1=i10,i1n                                              2d4s21
              i1m=i1-1                                                  2d4s21
              do k=0,nfh-1                                              9d26s23
               do ir=0,nrootu-1                                         9d26s23
                kr=ir+nrootu*k                                          9d26s23
                jvmat=ivmat(isb)+nvirt(jsbv)*kr+nn*i1m                   2d4s21
                do jv=0,nvirt(jsbv)-1                                    2d4s21
                 irow=iint+nvirt(isbv2)*jv                               2d4s21
                 irowd=iden+nvirt(isbv2)*jv+nrowh*ir                    9d26s23
                 do iv2=ibot,nvirt(isbv2)-1                                2d4s21
                  itri=((iv2*(iv2-1))/2)+iv1                              2d4s21
                  irec=iv1+nvirt(isbv1)*iv2                               2d4s21
                  iadv=ioffdnon+itri+isw*(irec-itri)+nvv*kr               2d4s21
                  gd(iadv)=gd(iadv)+fctr*bc(jvmat+jv)*bc(irow+iv2)      4d29s24
                  bc(irowd+iv2)=bc(irowd+iv2)+bc(jvmat+jv)*vd(iadv)     9d26s23
                 end do                                                  2d4s21
                end do                                                   2d4s21
               end do                                                    2d4s21
              end do                                                    9d26s23
              iint=iint+nrow                                            2d4s21
              iden=iden+nrow                                            9d26s23
             end do                                                     2d4s21
             i10=1                                                      2d4s21
            end do                                                      2d4s21
           end if                                                       2d4s21
           call ilimts(irefo(isn),nvirt(isbv2),mynprocg,mynowprog,il,ih,2d4s21
     $         i1s,i1e,i2s,i2e)                                         2d4s21
           nhere=ih+1-il                                                9d26s23
           i2eu=invk1(1,jsbv,isbv1,isn,2)                               2d4s21
           icase=invk1(2,jsbv,isbv1,isn,2)                              2d4s21
           i10=i1s                                                      2d4s21
           i1n=irefo(isn)                                               2d4s21
           iint=i3x(i2eu)                                               2d4s21
           iden=id3x(i2eu)                                              9d26s23
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
              do k=0,nfh-1                                              9d26s23
               do ir=0,nrootu-1                                         9d26s23
                kr=ir+nrootu*k                                          9d26s23
                jvmat=ivmat(isb)+nvirt(jsbv)*kr+nn*i1m                   2d4s21
                do iv1=0,itop                                            2d4s21
                 sum=0d0                                                 2d4s21
                 itri=((iv2*(iv2-1))/2)+iv1                              2d4s21
                 irec=iv1+nvirt(isbv1)*iv2                               2d4s21
                 iadv=ioffdnon+itri+isw*(irec-itri)+nvv*kr               2d4s21
                 do jv=0,nvirt(jsbv)-1                                    2d4s21
                  ix=max(jv,iv1)                                         2d4s21
                  in=min(jv,iv1)                                         2d4s21
                  irow=iint+((ix*(ix+1))/2)+in                           2d4s21
                  irowd=iden+((ix*(ix+1))/2)+in+nrowh*ir                9d26s23
                  sum=sum+bc(jvmat+jv)*bc(irow)                          2d4s21
                  bc(irowd)=bc(irowd)+bc(jvmat+jv)*vd(iadv)*bc(ipsr+kr) 9d26s23
                 end do                                                  2d4s21
                 gd(iadv)=gd(iadv)+fctr*sum*bc(ipsr+kr)                 4d29s24
                end do                                                   2d4s21
               end do                                                   9d26s23
              end do                                                    2d4s21
              iint=iint+nrow                                            2d4s21
              iden=iden+nrow                                            9d26s23
             end do                                                     2d4s21
             i10=1                                                      2d4s21
            end do                                                      2d4s21
           else if(icase.eq.1)then                                      2d4s21
            nrow=nvirt(jsbv)*nvirt(isbv1)                               2d4s21
            nrowh=nrow*nhere                                            9d26s23
            do i2=i2s,i2e                                               2d4s21
             iv2=i2-1                                                   2d4s21
             itops=iv2-1                                                2d4s21
             itopn=nvirt(isbv1)-1                                       2d4s21
             itop=itops+isw*(itopn-itops)                               2d4s21
             if(i2.eq.i2e)i1n=i1e                                       2d4s21
             do i1=i10,i1n                                              2d4s21
              i1m=i1-1                                                  2d4s21
              do k=0,nfh-1                                              9d26s23
               do ir=0,nrootu-1                                         9d26s23
                kr=ir+nrootu*k                                          9d26s23
                jvmat=ivmat(isb)+nvirt(jsbv)*kr+nn*i1m                   2d4s21
                do iv1=0,itop                                            2d4s21
                 sum=0d0                                                 2d4s21
                 itri=((iv2*(iv2-1))/2)+iv1                              2d4s21
                 irec=iv1+nvirt(isbv1)*iv2                               2d4s21
                 iadv=ioffdnon+itri+isw*(irec-itri)+nvv*kr               2d4s21
                 irow=iint+nvirt(jsbv)*iv1                               2d4s21
                 irowd=iden+nvirt(jsbv)*iv1+nrowh*ir                    9d26s23
                 do jv=0,nvirt(jsbv)-1                                    2d4s21
                  sum=sum+bc(jvmat+jv)*bc(irow+jv)                       2d4s21
                  bc(irowd+jv)=bc(irowd+jv)+bc(jvmat+jv)*vd(iadv)       9d27s23
     $                 *bc(ipsr+kr)                                     9d27s23
                 end do                                                  2d4s21
                 gd(iadv)=gd(iadv)+fctr*sum*bc(ipsr+kr)                 4d29s24
                end do                                                   2d4s21
               end do                                                   9d26s23
              end do                                                    2d4s21
              iint=iint+nrow                                            2d4s21
              iden=iden+nrow                                            9d27s23
             end do                                                     2d4s21
             i10=1                                                      2d4s21
            end do                                                      2d4s21
           else                                                         2d4s21
            nrow=nvirt(jsbv)*nvirt(isbv1)                               2d4s21
            nrowh=nrow*nhere                                            9d27s23
            do i2=i2s,i2e                                               2d4s21
             iv2=i2-1                                                   2d4s21
             itops=iv2-1                                                2d4s21
             itopn=nvirt(isbv1)-1                                       2d4s21
             itop=itops+isw*(itopn-itops)                               2d4s21
             if(i2.eq.i2e)i1n=i1e                                       2d4s21
             do i1=i10,i1n                                              2d4s21
              i1m=i1-1                                                  2d4s21
              do k=0,nfh-1                                              9d26s23
               do ir=0,nrootu-1                                         9d26s23
                kr=ir+nrootu*k                                          9d26s23
                jvmat=ivmat(isb)+nvirt(jsbv)*kr+nn*i1m                   2d4s21
                do jv=0,nvirt(jsbv)-1                                    2d4s21
                 irow=iint+nvirt(isbv1)*jv                               2d4s21
                 irowd=iden+nvirt(isbv1)*jv+nrowh*ir                    9d27s23
                 fact=bc(jvmat+jv)*bc(ipsr+kr)*fctr                     4d29s24
                 do iv1=0,itop                                           2d4s21
                  itri=((iv2*(iv2-1))/2)+iv1                              2d4s21
                  irec=iv1+nvirt(isbv1)*iv2                               2d4s21
                  iadv=ioffdnon+itri+isw*(irec-itri)+nvv*kr               2d4s21
                  gd(iadv)=gd(iadv)+fact*bc(irow+iv1)                    2d4s21
                  bc(irowd+iv1)=bc(irowd+iv1)+fact*vd(iadv)             9d27s23
                 end do                                                  2d4s21
                end do                                                   2d4s21
               end do                                                    2d4s21
              end do                                                    9d27s23
              iint=iint+nrow                                            2d4s21
              iden=iden+nrow                                            9d27s23
             end do                                                     2d4s21
             i10=1                                                      2d4s21
            end do                                                      2d4s21
           end if                                                       2d4s21
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
       ibcoff=ibcbmat                                                   2d3s21
      end do                                                            12d18s20
      ibcoff=ircv                                                       1d30s21
      return
      end
