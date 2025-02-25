c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dhf(nbtot,ih0,nbb,ivecs,ilsz,idbly,idoub,iacto,isinfo,
     $     icas,idorel4c,eposcut,nnbb,idorel,ascale,nbaslarge,nbassmall,
     $     isnorm,iblstor,ibsstor,nsqbas,natom,ngaus,ibdat,nbasis,scopy,3d17s20
     $     potdws,myguess,ivecr,lprint,bc,ibc)                          11d9s22
      implicit real*8 (a-h,o-z)
      real*16 fl                                                        2d26s20
      character*19 ostng                                                4d20s18
      character*1 extrapq                                               3d17s20
      include "common.store"
      logical myguess,lprint                                            10d5s22
      dimension idbly(*),idoub(*),iacto(*),isinfo(3,*),ivecs(*),
     $     irtyp(5),ioooo(2),scopy(*)                                   3d16s20
      COMMON/FACT16/FL(922),NCALL                                       2d26s20
      common/singcm/iuse,nff
      if(lprint)then                                                    10d5s22
       write(6,*)('Hi, my name is dhf')
       write(6,*)('nbaslarge = '),nbaslarge                              3d17s20
       write(6,*)('nbassmall = '),nbassmall                              3d17s20
      end if                                                            10d5s22
      isb=1
      ibcoffo=ibcoff                                                    3d1s20
         ih0i=ih0+nbtot*nbtot
         do i=0,nbtot-1                                                 2d22s20
          do j=0,i-1                                                    2d21s20
           ji=j+nbtot*i                                                 2d21s20
           ij=i+nbtot*j                                                 2d21s20
           bc(ih0+ji)=bc(ih0+ij)                                        2d21s20
           bc(ih0i+ji)=-bc(ih0i+ij)                                     2d21s20
          end do                                                        2d21s20
         end do                                                         2d21s20
         ih0ao=ibcoff                                                   2d24s20
         idwst1=ih0ao+nbb                                               2d24s20
         ibcoff=idwst1+nbb                                              2d21s20
         call enough('dhf.  1',bc,ibc)
         do i=0,nbb-1                                                   2d24s20
          bc(ih0ao+i)=bc(ih0+i)                                         2d24s20
         end do                                                         2d24s20
         idwst1i=idwst1+nbtot*nbtot                                     2d21s20
         do ipass=1,2                                                   2d21s20
          call dgemm('n','n',nbtot,nbtot,nbtot,1d0,bc(ih0),nbtot,       2d21s20
     $         bc(ivecs(isb)),nbtot,0d0,bc(idwst1),nbtot,               2d21s20
     d' dhf.  1')
          do i=0,nbtot-1                                                2d21s20
           do j=0,nbtot-1                                               2d21s20
            ji=idwst1+j+nbtot*i                                         2d21s20
            ij=ih0+i+nbtot*j                                            2d21s20
            bc(ij)=bc(ji)                                               2d21s20
           end do                                                       2d21s20
          end do                                                        2d21s20
          call dgemm('n','n',nbtot,nbtot,nbtot,1d0,bc(ih0i),nbtot,      2d21s20
     $         bc(ivecs(isb)),nbtot,0d0,bc(idwst1),nbtot,               2d21s20
     d' dhf.  2')
          do i=0,nbtot-1                                                2d21s20
           do j=0,nbtot-1                                               2d21s20
            ji=idwst1+j+nbtot*i                                         2d21s20
            ij=ih0i+i+nbtot*j                                           2d21s20
            bc(ij)=bc(ji)                                               2d21s20
           end do                                                       2d21s20
          end do                                                        2d21s20
         end do                                                         2d21s20
         if(icas.eq.0)then                                              2d25s20
          nelec=idbly(1)                                                10d5s22
          nstate=1                                                      2d25s20
          iocc=ibcoff                                                   2d25s20
          ibcoff=iocc+nelec                                             2d25s20
          call enough('dhf.  2',bc,ibc)
          do i=0,nelec-1                                                2d25s20
           ibc(iocc+i)=i+1                                              2d25s20
          end do                                                        2d25s20
          nox=nelec
         else                                                           2d25s20
          write(6,*)('how many electrons are there? '),idoub(1),        2d26s20
     $         isinfo(3,1)                                              2d26s20
          nelec=idoub(1)*2+isinfo(3,1)                                  2d26s20
          nox=(idoub(1)+iacto(1))*2                                     2d26s20
          write(6,*)('I think '),nelec
          itop=iacto(1)*2+1                                             2d26s20
          ibot1=isinfo(3,1)+1                                           2d26s20
          ibot2=itop-ibot1+1                                            2d26s20
          if(min(itop,ibot1,ibot2).lt.1)then                            2d26s20
           write(6,*)('factorial failure !!! '),itop,ibot1,ibot2
           call dws_sync
           call dws_finalize
          end if                                                        2d26s20
          xgot=fl(itop)-fl(ibot1)-fl(ibot2)                             2d26s20
          xgote=exp(xgot)                                                2d26s20
          write(6,*)('no. states? '),xgot,xgote
          nstate=nint(xgote)                                            2d26s20
          write(6,*)('setting iocc to ibcoff = '),ibcoff
          iocc=ibcoff                                                   2d26s20
          ibcoff=iocc+nstate*nelec                                      2d26s20
          call enough('dhf.  3',bc,ibc)
          istate=0                                                      2d26s20
          do i=1,iacto(1)*2-isinfo(3,1)+1                               2d26s20
           jocc=iocc+istate*nelec                                       2d26s20
           do j=1,idoub(1)*2                                            2d26s20
            ibc(jocc)=j                                                 2d26s20
            jocc=jocc+1
           end do                                                       2d26s20
           ibc(jocc)=i+idoub(1)*2                                       2d26s20
           jocc=jocc+1                                                  2d26s20
           do j=idoub(1)*2+2,nelec                                      2d26s20
            ibc(jocc)=0                                                 2d26s20
            jocc=jocc+1                                                 2d26s20
           end do                                                       2d26s20
           jocc=iocc+istate*nelec                                       2d26s20
           istate=istate+1                                              2d26s20
          end do                                                        2d26s20
          iiso=0                                                         2d26s20
          nso=istate                                                    2d26s20
          do ie=2,isinfo(3,1)                                           2d26s20
           nox=(idoub(1)+iacto(1))*2-(isinfo(3,1)-ie)                   2d26s20
           iel=idoub(1)*2+ie-1                                          2d26s20
           isn=istate                                                   2d26s20
           ise=iiso+nso-1                                                2d26s20
           do i=iiso,ise                                                 2d26s20
            jocc=iocc+nelec*i                                           2d26s20
            iol=ibc(jocc+iel-1)                                           2d26s20
            do j=iol+1,nox                                              2d26s20
             kocc=iocc+nelec*istate                                     2d26s20
             do l=0,nelec-1                                             2d26s20
              ibc(kocc+l)=ibc(jocc+l)                                   2d26s20
             end do                                                     2d26s20
             ibc(kocc+iel)=j                                            2d26s20
             istate=istate+1                                            2d26s20
            end do                                                      2d26s20
            nso=istate-isn                                              2d26s20
            iiso=isn                                                     2d26s20
           end do
          end do                                                        2d26s20
          write(6,*)('what we have for iiso and nso '),iiso,nso
          if(nso.ne.nstate)then
           write(6,*)('ooops, nso ne nstate '),nso,nstate
           call dws_sync
           call dws_finalize
           stop
          end if
          if(iiso.ne.0)then                                             2d26s20
           nmove=nso*nelec                                              2d26s20
           jocc=iocc+iiso*nelec                                         2d26s20
           do i=0,nmove-1                                               2d26s20
            ibc(iocc+i)=ibc(jocc+i)                                     2d26s20
           end do                                                       2d26s20
          end if                                                        2d26s20
          if(idorel4c.gt.0)then                                         2d26s20
           write(6,*)('set nstate to min nstate '),idorel4c             2d26s20
           nstate=min(nstate,idorel4c)                                  2d26s20
          end if                                                        2d26s20
          write(6,*)('states to average over: ')
          do i=1,nstate
           jocc=iocc+nelec*(i-1)
           write(6,*)i,(ibc(jocc+j),j=0,nelec-1)                        2d26s20
          end do                                                        2d26s20
         end if                                                         2d25s20
         ieigz=ibcoff                                                   2d22s20
         ivecz=ieigz+nbtot                                              2d22s20
         iveczi=ivecz+nnbb                                              3d16s20
         iden=ivecz+nbb                                                 2d23s20
         ibcoff=iden+nbb
         ifock=ibcoff                                                   2d24s20
         ibcoff=ifock+nbb                                               2d24s20
         call enough('dhf.  4',bc,ibc)
         do i=0,nbb-1                                                   2d24s20
          bc(ifock+i)=bc(ih0+i)                                         2d24s20
         end do                                                         2d24s20
         call chj(bc(ifock),nbtot,bc(ieigz),bc(ivecz),eposcut,nes,iden, 2d24s20
     $        nelec,bc(ivecs(isb)),nstate,ibc(iocc),1,scopy,nbaslarge,  3d16s20
     $       nbassmall,bc,ibc)                                          11d9s22
         ifocki=ifock+nnbb                                              2d24s20
         ih0aoi=ih0ao+nnbb                                              2d24s20
         iscf=0                                                         2d24s20
         extrapq=' '                                                    3d17s20
         ndiis=100                                                      3d16s20
         iscfdiis=0                                                     3d17s20
         iresidiis=ibcoff                                               3d16s20
         idendiis=iresidiis+ndiis*nnbb                                  3d16s20
         ibcoff=idendiis+ndiis*nnbb                                     3d16s20
         ideno=ibcoff                                                   2d24s20
         ibcoff=ideno+nbb                                               2d24s20
         iveco=ibcoff                                                   3d16s20
         ivecoi=iveco+nbtot*nox                                         3d16s20
         ibcoff=ivecoi+nbtot*nox                                        3d16s20
         call enough('dhf.  5',bc,ibc)
         idenoi=ideno+nnbb                                              2d24s20
         ideni=iden+nnbb                                                2d24s20
         if(lprint)write(6,*)('idorel = '),idorel                       10d5s22
         irtyp(1)=1                                                     2d25s20
         irtyp(4)=0                                                     2d25s20
         irtyp(5)=0                                                     2d25s20
         if(iabs(idorel).eq.1)then                                      2d26s20
          irtyp(2)=0                                                     2d25s20
          irtyp(3)=0                                                     2d25s20
         else if(mod(iabs(idorel),2).eq.0)then                          2d26s20
          irtyp(2)=1                                                     2d25s20
          irtyp(3)=0                                                     2d25s20
         else                                                           2d26s20
          irtyp(2)=1                                                     2d25s20
          irtyp(3)=1                                                     2d25s20
         end if                                                         2d26s20
         if(lprint)write(6,2822)irtyp                                   10d5s22
 2822    format('2e integral type factors ',/1x,'llll',5x,i1,/1x,'llss',2d25s20
     $        5x,i1,/1x,'ssss',5x,i1,/1x,'1/2rij',3x,i1,/1x,'1/2rij^3', 11d17s22
     $        x,i1)                                                     11d17s22
         if(lprint)then                                                 10d5s22
          write(6,*)('starting SCF iterations ...')
          write(6,2821)                                                  2d25s20
         end if                                                         10d5s22
 2821    format('iter.    energy',8x,'  e change den change diis')      10d5s22
         idiis=0                                                        3d16s20
 2820    continue                                                       2d24s20
         iscf=iscf+1                                                    2d24s20
         if(iscf.gt.100)then
          call dws_sync
          call dws_finalize
          stop
         end if
         do i=0,nbb-1                                                   2d24s20
          bc(ideno+i)=bc(iden+i)                                        2d24s20
         end do                                                         2d24s20
         nbb1=nint(sqrt(dfloat(nbb/2)))                                 11d17s22
         call paraeri4c(natom,ngaus,ibdat,nbasis,bc(ifock),bc(iden),    2d24s20
     $       ibc(iblstor),ibc(ibsstor),nsqbas,0,ascale,                 2d25s20
     $            nbaslarge,nbassmall,ilsz,isnorm,nbtot,irtyp,bc,ibc)   11d9s22
         do i=0,nbb-1                                                   2d24s20
          bc(ifock+i)=bc(ifock+i)+bc(ih0ao+i)                           2d24s20
         end do                                                         2d24s20
c
         egy=potdws                                                     3d17s20
         do i=0,nnbb-1
          avgr=0.5d0*(bc(ih0ao+i)+bc(ifock+i))                          2d24s20
          avgi=0.5d0*(bc(ih0aoi+i)+bc(ifocki+i))                        2d24s20
          egy=egy+bc(ideno+i)*avgr-bc(idenoi+i)*avgi                    2d24s20
         end do
         do i=0,nbtot*nox-1
          bc(iveco+i)=bc(ivecz+i)                                       3d16s20
          bc(ivecoi+i)=bc(iveczi+i)                                     3d16s20
         end do                                                         3d16s20
         do ipass=1,2                                                   2d21s20
          call dgemm('n','n',nbtot,nbtot,nbtot,1d0,bc(ifock),nbtot,       2d21s20
     $         bc(ivecs(isb)),nbtot,0d0,bc(idwst1),nbtot,               2d21s20
     d' dhf.  3')
          do i=0,nbtot-1                                                2d21s20
           do j=0,nbtot-1                                               2d21s20
            ji=idwst1+j+nbtot*i                                         2d21s20
            ij=ifock+i+nbtot*j                                            2d21s20
            bc(ij)=bc(ji)                                               2d21s20
           end do                                                       2d21s20
          end do                                                        2d21s20
          call dgemm('n','n',nbtot,nbtot,nbtot,1d0,bc(ifocki),nbtot,      2d21s20
     $         bc(ivecs(isb)),nbtot,0d0,bc(idwst1),nbtot,               2d21s20
     d' dhf.  4')
          do i=0,nbtot-1                                                2d21s20
           do j=0,nbtot-1                                               2d21s20
            ji=idwst1+j+nbtot*i                                         2d21s20
            ij=ifocki+i+nbtot*j                                           2d21s20
            bc(ij)=bc(ji)                                               2d21s20
           end do                                                       2d21s20
          end do                                                        2d21s20
         end do                                                         2d21s20
         call chj(bc(ifock),nbtot,bc(ieigz),bc(ivecz),eposcut,nes,iden, 2d24s20
     $        nelec,bc(ivecs(1)),nstate,ibc(iocc),1,scopy,nbaslarge,    3d16s20
     $        nbassmall,bc,ibc)                                         11d9s22
         tdiff=0d0                                                       10d5s22
         iterm=idiis                                                    10d5s22
         jdendiis=idendiis+nnbb*iterm                                   10d5s22
         extrapq=' '                                                    3d17s20
         do i=0,nbtot-1                                                 3d16s20
          do j=0,i-1                                                    3d16s20
           iad=ideno+j+nbtot*i                                          3d16s20
           jad=iden+j+nbtot*i                                           3d16s20
           diff=bc(jad)-bc(iad)                                         3d16s20
           tdiff=tdiff+diff**2                                          10d5s22
           bc(iad)=bc(jad)                                              10d5s22
           bc(jdendiis)=bc(jad)                                         10d5s22
           jdendiis=jdendiis+1                                          3d16s20
           iad=iad+nnbb                                                 3d16s20
           jad=jad+nnbb                                                 3d16s20
           diff=bc(jad)-bc(iad)                                         3d16s20
           tdiff=tdiff+diff**2                                          10d5s22
           bc(iad)=bc(jad)                                              10d5s22
           bc(jdendiis)=bc(jad)                                         3d16s20
           jdendiis=jdendiis+1                                          3d16s20
          end do                                                        3d16s20
          iad=ideno+i+nbtot*i                                           3d16s20
          jad=iden+i+nbtot*i                                            3d16s20
          diff=bc(jad)-bc(iad)                                          3d16s20
          tdiff=tdiff+diff**2                                           10d5s22
          bc(iad)=bc(jad)                                               10d5s22
          bc(jdendiis)=bc(jad)                                          3d16s20
          jdendiis=jdendiis+1                                           3d16s20
         end do                                                         3d16s20
         dden=sqrt(tdiff/dfloat(nbb))                                   10d5s22
         if(idiis.eq.0)then                                              10d5s22
          jresidiis=iresidiis+nnbb                                      10d5s22
          jresidiis0=iresidiis                                          10d5s22
          do i=0,nbtot-1                                                 3d16s20
           do j=0,i-1                                                    3d16s20
            jad=iden+j+nbtot*i                                           3d16s20
            bc(jresidiis)=bc(jad)                                         10d5s22
            jresidiis=jresidiis+1                                       10d5s22
            bc(jresidiis0)=0d0                                          10d5s22
            jresidiis0=jresidiis0+1                                       10d5s22
            jad=jad+nnbb                                                 3d16s20
            bc(jresidiis)=bc(jad)                                         3d16s20
            jresidiis=jresidiis+1                                       10d5s22
            bc(jresidiis0)=0d0                                          10d5s22
            jresidiis0=jresidiis0+1                                     10d5s22
           end do                                                        3d16s20
           jad=iden+i+nbtot*i                                            3d16s20
           bc(jresidiis)=bc(jad)                                          3d16s20
           jresidiis=jresidiis+1                                        10d5s22
           bc(jresidiis0)=0d0                                           10d5s22
           jresidiis0=jresidiis0+1                                        10d5s22
          end do                                                         3d16s20
         end if                                                         10d5s22
         if(idiis.ge.1)then                                             10d5s22
          itern=iterm-1                                                 10d5s22
          ifmat=ibcoff                                                  10d5s22
          ivec2=ifmat+nnbb*iterm                                        10d5s22
          ibcoff=ivec2+nnbb                                             10d5s22
          call enough('dhf.  6',bc,ibc)
          do i=0,iterm-1                                                10d5s22
           jfmat=ifmat+nnbb*i                                           10d5s22
           iad1=idendiis+i*nnbb                                         10d5s22
           iad3=iresidiis+i*nnbb                                        10d5s22
           do j=0,nnbb-1                                                10d5s22
            bc(jfmat+j)=bc(iad1+j)-bc(iad3+j)                           10d5s22
           end do                                                       10d5s22
          end do                                                        10d5s22
          iad1=idendiis+nnbb*iterm                                      5d19s22
          iad3=iresidiis+nnbb*iterm                                     10d5s22
          do j=0,nnbb-1                                                 10d5s22
           bc(ivec2+j)=bc(iad3+j)-bc(iad1+j)                            10d5s22
          end do                                                        10d5s22
          do i=0,iterm-1                                                10d5s22
           jfmat=ifmat+nnbb*i                                           10d5s22
           do j=0,nnbb-1                                                10d5s22
            bc(jfmat+j)=bc(jfmat+j)+bc(ivec2+j)                         10d5s22
           end do                                                       10d5s22
          end do                                                        10d5s22
          iwgt=ibcoff                                                   10d5s22
          icoef=iwgt+nnbb                                               10d5s22
          iscr=icoef+iterm                                              10d5s22
          ibcoff=iscr+2*iterm+iterm*iterm+2*iterm*nnbb+nnbb             10d5s22
          call enough('dhf.  7',bc,ibc)
          do i=iwgt,icoef-1                                             10d5s22
           bc(i)=1d0                                                    10d5s22
          end do                                                        10d5s22
          iuse=0                                                        10d5s22
          call lsqfit2(bc(ifmat),nnbb,iterm,bc(ivec2),nnbb,1,           10d5s22
     $         nnbb,bc(icoef),iterm,bc(iscr),bc(iwgt),0,rmsdiis,bc,ibc) 3d27s23
          sum=0d0                                                       10d5s22
          do i=0,iterm-1                                                10d5s22
           sum=sum+bc(icoef+i)                                          10d5s22
          end do                                                        10d5s22
          cn=1d0-sum                                                    10d5s22
          if(cn.lt.0.25d0)then                                          10d5s22
           idiis=-1                                                     10d5s22
          else                                                          10d5s22
           idiisp=idiis+1                                               10d5s22
           bc(icoef+iterm)=cn                                            10d5s22
           jresidiis=iresidiis+nnbb*idiisp                              10d5s22
           call dgemm('n','n',nnbb,1,idiisp,1d0,bc(idendiis),nnbb,      10d5s22
     $         bc(icoef),idiisp,0d0,bc(jresidiis),nnbb,                 10d5s22
     d'dhf.  1')
           call dws_bcast(bc(jresidiis),nnbb)                            10d5s22
           do i=0,nbtot-1                                                 3d16s20
            do j=0,i-1                                                    3d16s20
             jad=iden+j+nbtot*i                                           3d16s20
             iad=iden+i+nbtot*j                                          10d5s22
             bc(jad)=bc(jresidiis)                                       10d5s22
             bc(iad)=bc(jresidiis)                                       10d5s22
             jresidiis=jresidiis+1                                       10d5s22
             jad=jad+nnbb                                                 3d16s20
             iad=iad+nnbb                                                10d5s22
             bc(jad)=bc(jresidiis)                                       10d5s22
             bc(iad)=-bc(jresidiis)                                       10d5s22
             jresidiis=jresidiis+1                                       10d5s22
            end do                                                        3d16s20
            jad=iden+i+nbtot*i                                            3d16s20
            bc(jad)=bc(jresidiis)                                        10d5s22
            jresidiis=jresidiis+1                                        10d5s22
           end do                                                         3d16s20
          end if                                                        10d5s22
          ibcoff=ifmat                                                  10d5s22
         end if                                                         10d5s22
         idiis=idiis+1                                                  10d5s22
         if(iscf.gt.1)then                                              2d25s20
 1552     format(i4,1x,a19,2es9.1,2x,i4,a1)                                  3d17s20
 1551     format(f19.12)                                                  4d20s18
          write(ostng,1551)egy                                          2d25s20
          change=egy-eold                                               2d25s20
          if(change.ne.0d0)then                                         2d25s20
           ndig=nint(-log10(abs(change)))+1                             2d25s20
           ist=7+ndig                                                   2d25s20
           do k=ist,19                                                  2d25s20
            ostng(k:k)=' '                                              2d25s20
           end do                                                       2d25s20
          end if                                                        2d25s20
          if(lprint)write(6,1552)iscf,ostng,change,dden,idiis           10d5s22
  310     continue                                                      3d17s20
         end if                                                         3d16s20
         eold=egy
         if(dden.gt.1d-10)go to 2820                                     2d24s20
         if(lprint)then                                                 10d5s22
          write(6,*)('SCF calculations converged!')                      2d25s20
          write(6,*)('final orbitals ')
          write(6,*)('Fock eigenvalues ')
          call prntm2(bc(ieigz),1,nox,1)                                 3d17s20
          write(6,*)('real part ')
          call prntm2(bc(ivecz),nbtot,nox,nbtot)                         2d26s20
          write(6,*)('imag part ')
          call prntm2(bc(ivecz+nnbb),nbtot,nox,nbtot)                    2d26s20
         end if                                                         10d5s22
         if(idorel.lt.0)then                                            11d17s22
          if(lprint)                                                    11d17s22
     $        write(6,*)('Let us look at Gaunt or Breit corrections')   11d17s22
          nxx=nox*nox                                                   11d17s22
          if(lprint)write(6,*)('to do this we will build H matrix ...') 11d17s22
          ih0ms=ibcoff                                                  2d27s20
          itmps=ih0ms+nxx*2                                             2d27s20
          itmpt=itmps+nbb                                               2d27s20
          ibcoff=itmpt+nbb                                              2d27s20
          call enough('dhf.8',bc,ibc)
          do i=0,nbb-1                                                  2d27s20
           bc(itmps+i)=bc(ih0ao+i)                                      2d27s20
          end do                                                        2d27s20
          nrow=nbtot                                                    2d27s20
          ff=1d0
          do ipass=1,2                                                  2d27s20
           fm=-ff
           call dgemm('n','n',nrow,nox,nbtot,1d0,bc(itmps),nrow,        2d27s20
     $          bc(ivecz),nbtot,0d0,bc(itmpt),nrow)                     2d27s20
           call dgemm('n','n',nrow,nox,nbtot,fm,bc(itmps+nnbb),nrow,    2d28s20
     $          bc(ivecz+nnbb),nbtot,1d0,bc(itmpt),nrow)                2d27s20
           call dgemm('n','n',nrow,nox,nbtot,1d0,bc(itmps+nnbb),nrow,   2d27s20
     $          bc(ivecz),nbtot,0d0,bc(itmpt+nnbb),nrow)                2d27s20
           call dgemm('n','n',nrow,nox,nbtot,ff,bc(itmps),nrow,         2d28s20
     $          bc(ivecz+nnbb),nbtot,1d0,bc(itmpt+nnbb),nrow)           2d27s20
           do i=0,nox-1                                                 2d27s20
            do j=0,nrow-1                                               2d27s20
             jir=itmpt+j+nrow*i                                         2d27s20
             jii=jir+nnbb                                               2d27s20
             ijr=itmps+i+nox*j                                          2d27s20
             iji=ijr+nnbb                                               2d27s20
             bc(ijr)=bc(jir)                                            2d27s20
             bc(iji)=+bc(jii)                                           2d28s20
            end do                                                      2d27s20
           end do                                                       2d27s20
           nrow=nox                                                     2d27s20
           ff=-1d0                                                      2d28s20
          end do                                                        2d27s20
          do i=0,nox-1                                                  2d27s20
           do j=0,nox-1                                                 2d27s20
            ji1=itmps+j+nox*i                                           2d28s20
            ji1i=ji1+nnbb                                               2d27s20
            ji2=ih0ms+j+nox*i                                           2d27s20
            ji2i=ji2+nxx                                                2d27s20
            bc(ji2)=bc(ji1)                                             2d27s20
            bc(ji2i)=bc(ji1i)                                           2d27s20
           end do                                                       2d27s20
          end do                                                        2d27s20
          ibcoff=itmps                                                  2d27s20
          ioooo(1)=ibcoff                                               2d27s20
          ioooo(2)=ioooo(1)+nxx*nxx                                     2d28s20
          ibcoff=ioooo(2)+nxx*nxx                                       2d28s20
          call enough('dhf.9',bc,ibc)
          npass=3                                                       11d17s22
          do ipass=1,npass                                              11d17s22
           if(ipass.eq.1)then                                           11d17s22
            if(lprint)                                                  11d17s22
     $          write(6,*)('pass 1, check what we got from SCF calns')  11d17s22
           else if(ipass.eq.2)then                                      11d17s22
            if(lprint)write(6,*)('pass 2: Gaunt approximation to Breit')11d17s22
            irtyp(4)=2
           else                                                         11d17s22
            if(lprint)write(6,*)('pass 3: Full Breit')                            11d17s22
            irtyp(4)=1
            irtyp(5)=1
           end if                                                       11d17s22
           if(lprint)write(6,2822)irtyp                                             2d25s20
           do i=ioooo(1),ibcoff-1                                        2d27s20
            bc(i)=0d0                                                    2d27s20
           end do                                                        2d27s20
           call paraeri4ct(natom,ngaus,ibdat,nbasis,bc(ioooo(1)),nox,    2d27s20
     $         bc(ivecz),                                               2d27s20
     $       ibc(iblstor),ibc(ibsstor),0,ascale,                        2d28s20
     $            nbaslarge,nbassmall,ilsz,isnorm,nbtot,irtyp,bc,ibc)   11d17s22
           nxx2=nxx*nxx
           hr=potdws                                                    11d17s22
           hi=0d0
           sumaa=0d0                                                    3d14s20
           sumab=0d0                                                    3d14s20
           do jv=0,nelec-1                                              2d28s20
            iv=jv                                                       11d17s22
            ihr=ih0ms+iv*(nox+1)                                        2d28s20
            ihi=ihr+nxx                                                 2d28s20
            hr=hr+bc(ihr)                                               11d17s22
            hi=hi+bc(ihi)                                               11d17s22
            do kv=0,nelec-1                                             2d28s20
             lv=kv                                                      11d17s22
             jadd=lv+nox*(lv+nox*(iv+nox*iv))                           2d28s20
             kadd=lv+nox*(iv+nox*(iv+nox*lv))                           2d28s20
             term=0.5d0*(bc(ioooo(1)+jadd)-bc(ioooo(1)+kadd))           11d17s22
             hr=hr+term                                                 11d17s22
             hi=hi+0.5d0*(bc(ioooo(2)+jadd)-bc(ioooo(2)+kadd))          11d17s22
            end do                                                      2d28s20
           end do                                                       2d28s20
           if(lprint)then
            if(ipass.eq.1)then
             erefr=hr
             write(6,447)hr                                             12d2s22
  447        format('>variational Coulomb energy: ',f16.10)              12d2s22
            else if(ipass.eq.2)then
             gaunt=hr-erefr
             write(6,446)hr,gaunt                                       12d2s22
  446        format('>Gaunt energy: ',f16.10,' Gaunt correction: ',     12d2s22
     $            f16.10)                                               12d2s22
            else                                                         11d17s22
             breit=hr-erefr                                              11d17s22
             write(6,445)hr,breit                                       12d2s22
  445        format('>Breit energy: ',f16.10,' Breit correction: ',     12d2s22
     $            f16.10)                                               12d2s22
            end if                                                       11d17s22
           end if
          end do                                                        11d17s22
         end if                                                         11d17s22
         ibcoff=ibcoffo
         return
         end
