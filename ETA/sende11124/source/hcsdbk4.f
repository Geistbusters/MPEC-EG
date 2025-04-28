c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcsdbk4(ihsdiagb,nff1,iff1,i2sb,i2smb,mdoobp,ncsfb,    12d9s21
     $     isymbra,nff22,iff22,ff22,nfdat,ncsfk,ncsfk2,vd,i2sk,i2smk,   12d9s21
     $     mdookp,isymket,nec,mdon,nsymb,multh,irw1,irw2,nvirt,         12d9s21
     $     nrootu,irorip,isopt,ism,irel,irefo,norb,ih0n,                12d9s21
     $     nh0,ionex,i3x,iifmx,ntype,maxbx,sr2,l2e,nwiacc,bc,ibc)       1d26s23
      implicit real*8 (a-h,o-z)                                         12d18s20
c
c     to do ...
c     consolidate densities, and perhaps njhere in density calculation
c     did I swap iarg and jarg in 4th gandc4?
c
      integer*8 ihsdiagb(mdoobp,nsymb),i18,i28,ipack,ipackc,            12d9s21
     $   i38,i48,i1c,i1o,j2c,j2o,itestc,itesto,ipack8,last8(2),iff22(*),2d7s23
     $     gandcc,gandco,gandcb                                         2d7s23
      integer*1 nab1(2),nab2(2),imap(64),icode(64),ipackc1(8)           12d9s21
      integer*2 ipack2(4)                                               12d9s21
      external second                                                   2d18s21
      logical l3x,lprt,lnew,ldebug,lchoice                              3d17s21
      equivalence (ipack8,ipack4),(ipackc,ipackc1),(ipack,ipack2)       12d9s21
      dimension nff1(mdoobp,nsymb,2),iff1(*),nff22(mdookp,2,nsymb),       8d13s21
     $     nfdat(5,4,*),multh(8,8),nvirt(*),ncsfb(*),ncsfk(*),
     $     irel(*),ism(*),irefo(*),vd(*),itest(32,3),ncsfk2(4,*),         8d26s21
     $     nab4(2,3),ipack4(2),nl(4),npre(4),mpre(4),                   8d13s21
     $     id1visv(8,8),nd1visv(8,8),nokdc(8,8,4),nok3v(4),             12d22s20
     $     idhvnotv(4),ndhvnotv(4),id1vnotv(4,8,8),nd1vnotv(4,8,8),     12d21s20
     $     id3vnotv3(4),nd3vnotv3(4),loff(4),idkeep(2),ndkeep(2),       12d9s21
     $     mdkeep(2),keep(2),idhvnotvf(4),ibmat(8),ivmat(8),            2d4s21
     $     mdhvnotv(4),md3vnotv3(4),nok4f(4),nok4(4),nok3vf(4),         12d9s21
     $     nok33f(4,2),nok33(4,2),ff22(*),isopt(*),ih0n(*),nh0(*),      12d9s21
     $     ionex(8,8,8),i3x(8,8,8),iifmx(*),iifo(2,32),phss(2),         12d9s21
     $     iwpb1(4),iwpk1(4),iwpb2(4),iwpk2(4),mcsf(2),ipack2a(4),      12d9s21
     $     imy(4),igya(4),isy(4),ioxx(2)                                2d7s23
      include "common.store"
      data phss/1d0,-1d0/                                               12d9s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data loopx/11810000/
      common/paddcm/npadddi                                             6d7s22
      loop=0
      nrootm=nrootu-1                                                   1d4s21
      idoit=0                                                           3d1s21
      mdoob=mdoobp-1                                                      8d13s21
      irori=irorip-1                                                    9d6s21
      if(isopt(3).ne.0)then
       ieoro=1-irori                                                      9d6s21
      else                                                              10d14s21
       ieoro=irori                                                      9d6s21
      end if                                                            9d6s21
      do i=1,16                                                         12d13s21
       iifo(2,i)=1                                                      12d13s21
       iifo(1,i)=0                                                      12d13s21
      end do                                                            12d13s21
      do i=1,16                                                         12d9s21
       if(iifmx(i).ge.0)then                                            12d9s21
        im=i-1                                                          12d9s21
        ids=im/8
        im=im-ids*8
        ics=im/4
        im=im-ics*4
        ibs=im/2
        ias=im-ibs*2
        isgoal=ics+2*(ibs+2*(ias+2*ids))                                12d13s21
        isgoal=isgoal+1
        isgoal2=1-ics+2*(1-ibs+2*(1-ias+2*(1-ids)))                     12d13s21
        isgoal2=isgoal2+1
        iifo(2,iifmx(i)+1)=1
        if(iifmx(isgoal).ge.0)then
         iifo(1,iifmx(i)+1)=iifmx(isgoal)
        else if(iifmx(isgoal2).ge.0)then
         iifo(1,iifmx(i)+1)=iifmx(isgoal2)
         if(isopt(4).ne.0)iifo(2,iifmx(i)+1)=2
        end if                                                          12d9s21
       end if                                                           12d9s21
      end do                                                            12d9s21
      igoal=0
      igoal1=2
      igoal2=2
      igoal3=2
      igoal4=2
      igoal5=2
      igoal6=2
      igoal7=2
      last8(1)=-1                                                       2d8s21
c
c     maxbx is for input wavefcns ...
c     but here bra has ket roots, so maxbx may not be large enough.     6d7s22
c
      maxbxb=0                                                          6d7s22
      do isb=1,nsymb                                                    6d7s22
       isbv=multh(isb,isymbra)                                          6d7s22
       do ii=mdon+1,mdoobp                                              6d30s23
        iarg=ii-mdon                                                    6d7s22
        maxbxb=max(maxbxb,ncsfb(iarg)*nff1(ii,isb,1)*nvirt(isbv)*nrootu)6d8s22
       end do                                                           6d7s22
      end do                                                            6d7s22
      maxbxb=maxbxb+npadddi                                             6d7s22
      ircv=ibcoff                                                       1d29s21
      iacc=ircv+mynprocg                                                1d30s21
      ivs=iacc+mynprocg                                                 1d30s21
      igg=ivs+maxbx                                                     1d30s21
      ibcoff=igg+maxbxb                                                 6d7s22
      nacc=0                                                            1d30s21
      call enough('hcsdbk4.  1',bc,ibc)
      itransgg=0                                                        1d30s21
      loop=0
      nsing=0                                                           12d23s20
      norbx=norb+1                                                      12d18s20
      norbxx=norbx+1                                                    12d18s20
      norbxxx=norbxx+1                                                  12d18s20
      isymbk=multh(isymbra,isymket)                                     8d13s21
      do jsb=1,nsymb                                                    12d18s20
       jsbv=multh(jsb,isymbra)                                          8d13s21
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
        isbv12=multh(isb,isymket)                                       8d13s21
        isn=multh(isopt(1),multh(isbv12,jsbv))                          12d9s21
        ibmat(isb)=ibcoff                                               2d3s21
        if(min(irefo(isn),nvirt(jsbv)).gt.0)then                        2d3s21
         nftrip=nfdat(2,2,isb)+nfdat(2,3,isb)+nfdat(2,4,isb)             2d3s21
         nfh=nfdat(2,1,isb)+nftrip                                       2d3s21
         nn=nvirt(jsbv)*nrootu*nfh*ntype                                      2d3s21
         ibcoff=ibmat(isb)+nn*irefo(isn)                                3d2s21
         nbmat=nbmat+nn*irefo(isn)                                      3d2s21
        end if                                                          3d2s21
       end do                                                           3d2s21
       call enough('hcsdbk4.  2',bc,ibc)
       do i=ibcbmat,ibcoff-1                                            3d2s21
        bc(i)=0d0                                                       3d2s21
       end do                                                           3d2s21
       if(l2e.eq.0)then                                                 2d17s22
        ixint=ibcoff                                                     8d13s21
        ibcoff=ixint+nvirt(jsbv)                                         8d13s21
        call enough('hcsdbk4.  3',bc,ibc)
        do isb=1,nsymb                                                   2d3s21
         isbv12=multh(isb,isymket)                                       8d13s21
         isn=multh(isopt(1),multh(isbv12,jsbv))                          12d9s21
         nftrip=nfdat(2,2,isb)+nfdat(2,3,isb)+nfdat(2,4,isb)             2d3s21
         nfh=nfdat(2,1,isb)+nftrip                                       2d3s21
         if(min(irefo(isn),nvirt(jsbv)).gt.0)then                        2d3s21
          nn=nvirt(jsbv)*nrootu*nfh                                      12d9s21
          nnn=nfh*nrootu                                                 2d3s21
          do isbv1=1,nsymb                                               2d3s21
           isbv2=multh(isbv1,isbv12)                                     2d3s21
           if(isbv2.ge.isbv1)then                                        2d3s21
            call ilimts(irefo(isn),nvirt(isbv1),mynprocg,mynowprog,il,   12d9s21
     $          ih,i1s,i1e,i2s,i2e)                                     8d13s21
            nhere=ih+1-il                                                12d9s21
            if(isbv1.eq.isbv2)then                                       2d3s21
             i10=i1s                                                     2d3s21
             i1n=irefo(isn)                                              2d3s21
             do i2=i2s,i2e                                               2d3s21
              iv=i2-1                                                    2d3s21
              if(i2.eq.i2e)i1n=i1e                                       2d3s21
              do i1=i10,i1n                                              2d3s21
               i1m=i1-1
c                    ito=i3x(isb,isc,isd,it)+ibv+nvirt(isb)*(iav          10d6s21
c     $                  +nvirt(isa)*icoldc)                             10d6s21
               do it=0,ntype-1                                           12d9s21
                ixint=i3x(jsbv,isbv1,isn)+nvirt(jsbv)*(iv+nvirt(isbv1)    12d9s21
     $             *(i1+irefo(isn)*iv-il+nhere*it))                     12d9s21
                jb=ibmat(isb)+nn*(i1m+irefo(isn)*it)                      12d9s21
                do jv=0,nvirt(jsbv)-1                                     2d4s21
                 xxint=bc(ixint+jv)*sr2                                  12d13s21
                 do ir=0,nrootm                                           2d3s21
                  do k=0,nfdat(2,1,isb)-1                                 2d3s21
                   iadv=ioffvd+iv+nvirt(isbv2)*(ir+nrootu*k)              2d3s21
                   bc(jb+k)=bc(jb+k)+xxint*vd(iadv)                      12d13s21
                  end do                                                 2d3s21
                  jb=jb+nfh                                              2d3s21
                 end do                                                  2d3s21
                end do                                                   2d3s21
               end do                                                    12d9s21
              end do                                                     2d3s21
              i10=1                                                      2d3s21
             end do                                                      2d3s21
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
            do i2=i2s,i2e                                                2d3s21
             iv1=i2-1                                                    2d3s21
             ibots=i2                                                    2d3s21
             ibotn=0                                                     2d3s21
             ibot=ibots+isw*(ibotn-ibots)                                2d3s21
             if(i2.eq.i2e)i1n=i1e                                        2d3s21
             do i1=i10,i1n                                               2d3s21
              i1m=i1-1
              do iv2=ibot,nvirt(isbv2)-1                                 2d3s21
               itri=((iv2*(iv2-1))/2)+iv1                                2d3s21
               irec=iv1+nvirt(isbv1)*iv2                                 2d3s21
               ivv=itri+isw*(irec-itri)                                  2d3s21
               do it=0,ntype-1
                jb=ibmat(isb)+nn*(i1m+irefo(isn)*it)                      12d9s21
c     ito=i3x(isb,isc,isd,it)+ibv+nvirt(isb)*(iav          10d6s21
c                  +nvirt(isa)*icoldc)                             10d6s21
c i.e.
c     (v"v|v'n)Vv'v"
                ixint=i3x(jsbv,isbv1,isn)+nvirt(jsbv)*(iv2+nvirt(isbv2)   12d9s21
     $             *(i1+irefo(isn)*iv1-il+nhere*it))                    12d9s21
                do jv=0,nvirt(jsbv)-1                                     2d4s21
                 irow=jv+ixint                                            8d13s21
                 do ir=0,nrootm                                           2d3s21
                  do k=0,nfh-1                                            2d3s21
                   iadv=ioffvd+ivv+nvv*(ir+nrootu*k)                      2d3s21
                   bc(jb+k)=bc(jb+k)+bc(irow)*vd(iadv)                    2d3s21
                  end do                                                  2d3s21
                  jb=jb+nfh                                               2d3s21
                 end do                                                   2d3s21
                end do                                                    2d3s21
               end do                                                     12d9s21
              end do                                                      2d3s21
             end do                                                       2d3s21
             i10=1                                                        2d3s21
            end do                                                        2d3s21
            call ilimts(irefo(isn),nvirt(isbv2),mynprocg,mynowprog,il,   8d13s21
     $           ih,i1s,i1e,i2s,i2e)                                    8d13s21
            nhere=ih+1-il                                                 12d9s21
            i10=i1s                                                      2d3s21
            i1n=irefo(isn)                                               2d3s21
            do i2=i2s,i2e                                                2d3s21
             iv2=i2-1                                                    2d3s21
             itops=iv2-1                                                 2d3s21
             itopn=nvirt(isbv1)-1                                        2d3s21
             itop=itops+isw*(itopn-itops)                                2d3s21
             if(i2.eq.i2e)i1n=i1e                                        2d3s21
             do i1=i10,i1n                                               2d3s21
              i1m=i1-1
              do iv1=0,itop                                              2d3s21
               itri=((iv2*(iv2-1))/2)+iv1                                2d3s21
               irec=iv1+nvirt(isbv1)*iv2                                 2d3s21
               ivv=itri+isw*(irec-itri)                                  2d3s21
               do it=0,ntype-1                                            12d9s21
                itp=it+1                                                 12d13s21
                jb=ibmat(isb)+nn*(i1m+irefo(isn)*it)                      12d9s21
                ixint=i3x(jsbv,isbv2,isn)+nvirt(jsbv)*(iv1+nvirt(isbv1)   12d9s21
     $             *(i1+irefo(isn)*iv2-il+nhere*iifo(1,itp)))           12d13s21
                do jv=0,nvirt(jsbv)-1                                    2d4s21
                 irow=ixint+jv                                           8d13s21
                 xint=-bc(ixint+jv)*phss(iifo(2,itp))                    12d13s21
                 do ir=0,nrootm                                          2d3s21
                  do k=0,nfh-1                                            12d9s21
                   iadv=ioffvd+ivv+nvv*(ir+nrootu*k)                     2d3s21
                   bc(jb+k)=bc(jb+k)+xint*vd(iadv)                        12d9s21
                  end do                                                 2d3s21
                  jb=jb+nfh                                              2d3s21
                 end do                                                  2d3s21
                end do                                                   2d3s21
               end do                                                    2d3s21
              end do                                                      12d9s21
             end do                                                      2d3s21
             i10=1                                                       2d3s21
            end do                                                       2d3s21
            ioffvd=ioffvd+nvv*nrootu*nfh                                 2d3s21
           end if                                                        2d3s21
          end do                                                         2d3s21
          itmp=ibcoff                                                    2d3s21
          ibcoff=itmp+nn*irefo(isn)*ntype                                12d9s21
          call enough('hcsdbk4.  4',bc,ibc)
          jbmat=ibmat(isb)                                               2d3s21
          do n=0,irefo(isn)*ntype-1                                      12d9s21
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
          do i=0,nn*irefo(isn)*ntype-1                                   12d9s21
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
       end if                                                           2d17s22
       do nclo1p=mdon+1,mdoobp                                           8d13s21
        if(min(nff1(nclo1p,jsb,1),nvirt(jsbv)).gt.0)then                8d26s21
         nclo1=nclo1p-1                                                  12d12s20
         jarg=nclo1p-mdon                                                12d12s20
         nopen1=nec-2*nclo1                                              12d12s20
         nggb=nvirt(jsbv)*nrootu                                        8d13s21
         nggg=nff1(nclo1p,jsb,1)*ncsfb(jarg)                             12d29s20
         ncolt=nggg*nggb                                                8d13s21
         x2=dfloat(nff1(nclo1p,jsb,1))/dfloat(mynprocg)                 3d2s21
         call ilimts(1,nff1(nclo1p,jsb,1),mynprocg,mynowprog,           3d2s21
     $        i1l,i1h,i11s,i11e,i12s,i12e)                              12d19s20
         i1l=1+ncsfb(jarg)*(i1l-1)                                      3d2s21
         i1h=ncsfb(jarg)*i1h                                            3d2s21
         ibcgg=ibcoff                                                   1d30s21
         i18=1
         i38=1                                                          12d18s20
         i48=nggg                                                       12d29s20
         itransvs=0                                                     1d29s21
         koffdnon=1                                                     8d16s21
         do isb=1,nsymb                                                 12d18s20
          nfh=nfdat(2,1,isb)+nfdat(2,2,isb)+nfdat(2,3,isb)              2d3s21
     $         +nfdat(2,4,isb)                                          2d3s21
          ksbv12=multh(isb,isymket)                                     8d13s21
          loff(1)=0                                                     1d4s21
          do l=2,4                                                      1d4s21
           lm=l-1                                                       1d4s21
           loff(l)=loff(lm)+nfdat(2,lm,isb)                             1d4s21
          end do                                                        1d4s21
          nnt=loff(4)+nfdat(2,4,isb)                                    1d4s21
          do nclo2p=max(mdon+1,nclo1p-2),min(mdookp,nclo1p+3)            8d13s21
           if(nff22(nclo2p,1,isb).gt.0)then                             12d18s20
            nclo2=nclo2p-1                                              12d18s20
            iarg=nclo2p-mdon                                            12d18s20
            iargp=iarg+1                                                12d18s20
            nopen2=nec-2*nclo2p                                         12d18s20
            nopen2p=nopen2+2                                            12d18s20
            idhvisv=ibcoff                                              12d18s20
            nnj=ncsfb(jarg)*nfdat(2,1,isb)                               12d18s20
            lgoal=multh(jsb,isb)                                        12d22s20
            lgoal3=multh(isopt(1),multh(jsbv,ksbv12))                   12d9s21
            do isc=1,nsymb                                              12d18s20
             do l=1,4                                                   12d18s20
              nnl=ncsfb(jarg)*nfdat(2,l,isb)*irefo(isc)                 12d9s21
              if(isc.eq.lgoal)then                                      12d21s20
               idhvnotvf(l)=ibcoff                                      1d6s21
               idhvnotv(l)=idhvnotvf(l)+nnl                             1d6s21
               ndhvnotv(l)=idhvnotv(l)+nnl                              12d21s20
               mdhvnotv(l)=ndhvnotv(l)+irefo(isc)                       12d31s20
               ibcoff=mdhvnotv(l)+nfdat(2,l,isb)                        12d21s20
              end if                                                    12d21s20
             end do                                                     12d18s20
             iscv=multh(isc,jsbv)                                       12d18s20
             jscv=multh(isc,lgoal)                                      12d21s20
             do isd=1,nsymb                                             8d16s21
              iscdv=multh(iscv,isd)                                     12d18s20
              jscdv=multh(jscv,isd)                                     12d21s20
              nn=irefo(isd)*irefo(isc)                                  8d13s21
              nnn=nn*irefo(jscdv)*ntype                                 12d9s21
              do l=1,4                                                  12d19s20
               id1vnotv(l,isd,isc)=ibcoff                               12d19s20
               nd1vnotv(l,isd,isc)=id1vnotv(l,isd,isc)+nnn*ncsfb(jarg)   12d19s20
     $              *nfdat(2,l,isb)                                     12d19s20
               ibcoff=nd1vnotv(l,isd,isc)+nnn                           12d21s20
              end do                                                    12d19s20
             end do                                                     12d18s20
            end do                                                      12d18s20
            do ipass=1,1                                                3d1s21
             idkeep(ipass)=ibcoff                                       1d4s21
             ndkeep(ipass)=idkeep(ipass)+irefo(lgoal3)*nnt*ncsfb(jarg)  12d11s21
     $            *ntype                                                12d11s21
             mdkeep(ipass)=ndkeep(ipass)+irefo(lgoal3)*nnt*ntype        12d11s21
             ibcoff=mdkeep(ipass)+irefo(lgoal3)*nnt*ntype               12d13s21
            end do                                                      1d4s21
            do l=1,4
             nnl=ntype*ncsfb(jarg)*nfdat(2,l,isb)*irefo(lgoal3)         12d9s21
             id3vnotv3(l)=ibcoff                                        12d22s20
             nd3vnotv3(l)=id3vnotv3(l)+nnl                              12d22s20
             md3vnotv3(l)=nd3vnotv3(l)+irefo(lgoal3)*ntype              12d9s21
             ibcoff=md3vnotv3(l)+nfdat(2,l,isb)*ntype                   12d9s21
            end do                                                      12d22s20
            call enough('hcsdbk4.  5',bc,ibc)
            ibcb4=ibcoff-1                                              12d18s20
            if1o=nff1(nclo1p,jsb,2)                                      12d14s20
            jgs=igg                                                     8d13s21
            do if1=1,nff1(nclo1p,jsb,1)                                  12d14s20
             ist=1+ncsfb(jarg)*(if1-1)                                   1d5s21
             ien=ncsfb(jarg)*if1                                         1d5s21
             istu=max(ist,i1l)                                          1d5s21
             ienu=min(ien,i1h)                                          1d5s21
             if(ienu.ge.istu)then                                       1d5s21
              i11s=istu-ist+1                                           1d5s21
              njhere=ienu+1-istu                                        1d5s21
             else                                                       1d5s21
              njhere=0                                                  1d5s21
             end if                                                     1d5s21
             idoit=idoit+1                                              3d1s21
             loffdnon=koffdnon                                          8d16s21
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
               if(nl(l).gt.ncsfk2(l,iarg))lchoice=.true.                1d13s23
              end do                                                     3d17s21
              j2o=ibset(j2o,norbx)                                      12d18s20
              j2o=ibset(j2o,norbxx)                                     12d18s20
              if(njhere.gt.0)then                                       2d26s21
               gandcc=ieor(i1c,j2c)                                      10d13s22
               gandco=ieor(i1o,j2o)                                      10d13s22
               gandcb=ior(gandcc,gandco)                                 10d20s22
               ndifb=popcnt(gandcb)                                      10d20s22
               if(ndifb.le.4)then                                        10d20s22
                ndifs=popcnt(gandco)                                      10d13s22
                ndifd=popcnt(gandcc)                                      10d13s22
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
                 iunit=ibcoff                                            12d6s21
                 ibcoff=iunit+ncsfk(iarg)*ncsfk(iarg)                    12d9s21
                 call enough('hcsdbk4.  6',bc,ibc)
                 do iz=iunit,ibcoff-1                                    12d6s21
                  bc(iz)=0d0                                             12d6s21
                 end do                                                  12d6s21
                 do iz=0,ncsfk(iarg)-1                                   12d6s21
                  iad=iunit+iz*(ncsfk(iarg)+1)                           12d6s21
                  bc(iad)=1d0                                            12d6s21
                 end do                                                  12d6s21
                 call gandcr(i1c,i1o,j2c,j2o,nopen1,nopen2p,norbxx,     2d7s23
     $                nnot1,nab1,icode,imap,nx,irw1,irw2,iwpb1,iwpk1,   2d7s23
     $                bc,ibc)                                           2d7s23
                 call gencup(i2sb,i2smb,i2sk,i2smk,nopen1,nopen2p,nab1,  12d6s21
     $               iwpb1,iwpk1,ioutg,imatg,ntypeg,mcsf,bc(iunit),     12d6s21
     $               ncsfk(iarg),ncsfk(iarg),bc,ibc)                    11d14s22
                 do ixi=0,ntypeg-1                                        12d6s21
                  ipackc=ibc(ioutg+ixi)                                   12d6s21
                  imatu=imatg+ncsfb(jarg)*ncsfk(iarg)*ixi                12d7s21
                  if(ipackc1(1).gt.0)then                                   9d10s21
                   imy(3)=0                                              12d6s21
                   ipack2a(3)=ipackc1(1)
                  else                                                   12d6s21
                   imy(3)=1                                              12d6s21
                   ipack2a(3)=-ipackc1(1)                                12d6s21
                  end if                                                 12d6s21
                  if(ipackc1(2).gt.0)then                                12d6s21
                   imy(4)=0                                              12d6s21
                  else                                                   12d6s21
                   imy(4)=1                                              12d6s21
                  end if                                                 12d6s21
                  isy(3)=ism(ipack2a(3))                                 12d6s21
                  igya(3)=irel(ipack2a(3))-1                             12d6s21
                  phsh=1d0                                               12d6s21
                  isy(4)=multh(isy(3),isopt(1))                          12d7s21
                  if(ih0n(isy(3)).gt.0)then                              12d13s21
                   if(imy(3).ne.0.and.isopt(4).ne.0)phsh=-phsh           12d13s21
                  else                                                   12d6s21
                   if(isopt(2).ne.0)phsh=-phsh                           12d6s21
                   if(imy(4).ne.0.and.isopt(4).ne.0)phsh=-phsh           12d13s21
                  end if                                                 12d6s21
                  jmat=imatu                                             12d9s21
                  do l=1,4                                                12d18s20
                   if(nl(l).gt.0)then                                     12d18s20
                    nnl=ncsfb(jarg)*nfdat(2,l,isb)                       12d6s21
                    iad1=jvcv+iff22(jvcv+5+l)                             8d16s21
                    iad2=iad1+nl(l)                                       3d19s21
                    itmp=ibcoff                                           12d18s20
                    ibcoff=itmp+ncsfb(jarg)*nl(l)                        12d6s21
                    call enough('hcsdbk4.  7',bc,ibc)
                    call dgemm('n','n',ncsfb(jarg),nl(l),ncsfk2(l,iarg), 12d6s21
     $               1d0,bc(jmat),ncsfb(jarg),ff22(iad2),ncsfk2(l,iarg),8d16s21
     $                 0d0,bc(itmp),ncsfb(jarg),                        12d9s21
     d' hcsdbk4.  1')
c     jsbv=jsb*isymbra, ksbv=jsb*isymket                                8d13s21
c     isbv12=isb*isymbra, ksbv12=jsb*isymket                            8d13s21
c     isbv(1or2)=ksbv12*jsbv=jsb*isymbra*isb*isymket
c     isbv(1or2)=isbv12*ksbv=jsb*isymket*isb*isymbra                    8d13s21
                    jden=idhvnotv(l)+nnl*igya(3)                         12d6s21
                    ibc(ndhvnotv(l)+igya(3))=1                           12d6s21
                    jtmp=itmp                                             12d18s20
                    do iii=0,nl(l)-1                                      12d18s20
                     ip=iff22(iad1+iii)-1                                 8d16s21
                     jjden=jden+ncsfb(jarg)*ip                             12d18s20
                     ibc(mdhvnotv(l)+ip)=1                                12d21s20
                     do j=0,ncsfb(jarg)-1                                  12d18s20
                      bc(jjden+j)=bc(jjden+j)+bc(jtmp+j)*phsh            12d6s21
                     end do                                               12d18s20
                     jtmp=jtmp+ncsfb(jarg)                               12d6s21
                    end do                                                12d18s20
                    do i=1,norb                                          12d6s21
                     if(btest(j2c,i).and.btest(i1c,i).and.l2e.eq.0)then  2d17s22
                      js=ism(i)                                               9d10s21
                      jg=irel(i)-1                                            9d10s21
                      nn=irefo(js)*irefo(js)                             12d6s21
                      do jpass=0,1                                       11d22s21
c
c     (3v|ii)=(ii|3v)=+/-(ii|v3)                                                           12d6s21
c
                       ltest=jpass+2*(jpass+2*(imy(4)+2*imy(3)))         12d9s21
                       jtest=1-jpass+2*(1-jpass+2*(1-imy(4)              12d9s21
     $                      +2*(1-imy(3))))                              12d9s21
                       itestp=ltest+1                                         10d12s21
                       jtestp=jtest+1                                    11d22s21
                       iuse=-1                                           11d22s21
                       phsj=1d0                                          11d22s21
                       if(isopt(2).ne.0)phsj=-phsj                       12d9s21
                       if(iifmx(itestp).ge.0)then                        11d22s21
                        iuse=iifmx(itestp)                               11d22s21
                       else if(iifmx(jtestp).ge.0)then                   11d22s21
                        iuse=iifmx(jtestp)                               11d22s21
                        if(isopt(4).ne.0)phsj=-phsj                      11d22s21
                       end if                                            11d22s21
                       if(iuse.ge.0)then                                 11d22s21
                        icolj=jg+irefo(js)*(jg+irefo(js)*(igya(3)+       12d6s21
     $                     irefo(isy(3))*iuse))                         12d6s21
                        jdenj=id1vnotv(l,js,js)+nnl*icolj                    12d19s20
                        ibc(nd1vnotv(l,js,js)+icolj)=1                       12d19s20
                        jtmp=itmp                                              12d18s20
                        do iii=0,nl(l)-1                                         12d18s20
                         ii=iff22(iad1+iii)-1                              8d16s21
                         jdj=jdenj+ncsfb(jarg)*ii                             12d18s20
                         do j=0,ncsfb(jarg)-1                                   12d18s20
                          bc(jdj+j)=bc(jdj+j)+bc(jtmp+j)*phsj            12d6s21
                         end do                                                12d18s20
                         jtmp=jtmp+ncsfb(jarg)                           12d7s21
                        end do                                                 12d18s20
                       end if                                            12d6s21
c
c     (3i|iv)=+/-(i3|vi)
c
                       ltest=jpass+2*(imy(3)+2*(imy(4)+2*jpass))         12d6s21
                       jtest=1-jpass+2*(1-imy(3)+2*(1-imy(4)             12d6s21
     $                     +2*(1-jpass)))                               12d6s21
                       itestp=ltest+1                                    12d6s21
                       jtestp=jtest+1                                    12d6s21
                       iuse=-1                                           12d6s21
                       phsk=-1d0                                         12d6s21
                       if(isopt(2).ne.0)phsk=-phsk                       12d6s21
                       if(iifmx(itestp).ge.0)then                        12d6s21
                        iuse=iifmx(itestp)                               12d6s21
                       else if(iifmx(jtestp).ge.0)then                   12d6s21
                        iuse=iifmx(jtestp)                               12d6s21
                        if(isopt(4).ne.0)phsk=-phsk                      12d6s21
                       end if                                            12d6s21
                       if(iuse.ge.0)then                                 12d6s21
                        icolk=jg+irefo(js)*(igya(3)+irefo(isy(3))*(jg+    12d6s21
     $                     irefo(js)*iuse))                             12d6s21
                        jdenj=id1vnotv(l,js,isy(3))+nnl*icolk                    12d19s20
                        ibc(nd1vnotv(l,js,isy(3))+icolk)=1               12d6s21
                        jtmp=itmp                                              12d18s20
                        do iii=0,nl(l)-1                                         12d18s20
                         ii=iff22(iad1+iii)-1                              8d16s21
                         jdj=jdenj+ncsfb(jarg)*ii                             12d18s20
                         do j=0,ncsfb(jarg)-1                                   12d18s20
                          bc(jdj+j)=bc(jdj+j)+bc(jtmp+j)*phsk            12d6s21
                         end do                                                12d18s20
                         jtmp=jtmp+ncsfb(jarg)                           12d7s21
                        end do                                                 12d18s20
                       end if                                            12d6s21
                      end do                                             12d6s21
                     end if                                              12d6s21
                    end do                                               12d6s21
                    if(btest(i1c,ipack2a(3)).and.l2e.eq.0)then           2d17s22
c
c     x term
c         --   --         --
c     (3v|33)=(33|3v)=+/-(33|v3)
c
                     ltest=1-imy(3)+2*(1-imy(3)+2*(imy(4)+2*imy(3)))     12d6s21
                     jtest=imy(3)+2*(imy(3)+2*(1-imy(4)+2*(1-imy(3))))   12d6s21
                     itestp=ltest+1                                      12d6s21
                     jtestp=jtest+1                                      12d6s21
                     phsj=1d0                                            12d6s21
                     if(isopt(2).ne.0)phsj=-phsj                         12d9s21
                     iuse=-1                                             12d6s21
                     if(iifmx(itestp).ge.0)then                          12d6s21
                      iuse=iifmx(itestp)                                 12d6s21
                     else if(iifmx(jtestp).ge.0)then                     12d6s21
                      iuse=iifmx(jtestp)                                 12d6s21
                      if(isopt(4).ne.0)phsj=-phsj                        12d7s21
                     end if                                              12d6s21
                     if(iuse.ge.0)then                                   12d6s21
                      icolj=igya(3)+irefo(isy(3))*(igya(3)+irefo(isy(3)) 12d6s21
     $                   *(igya(3)+irefo(isy(3))*iuse))                 12d6s21
                      jdenj=id1vnotv(l,isy(3),isy(3))+nnl*icolj          12d6s21
                      ibc(nd1vnotv(l,isy(3),isy(3))+icolj)=1             12d6s21
                      jtmp=itmp                                              12d18s20
                      do iii=0,nl(l)-1                                         12d18s20
                       ii=iff22(iad1+iii)-1                              8d16s21
                       jdj=jdenj+ncsfb(jarg)*ii                             12d18s20
                       do j=0,ncsfb(jarg)-1                                   12d18s20
                        bc(jdj+j)=bc(jdj+j)+bc(jtmp+j)*phsj              12d6s21
                       end do                                                12d18s20
                       jtmp=jtmp+ncsfb(jarg)                             12d7s21
                      end do                                                 12d18s20
                     end if                                              12d6s21
c
c       - -       -   -
c     (33|3v)=+/-(33|v3)
c
                     ltest=1-imy(3)+2*(imy(3)+2*(imy(4)+2*(1-imy(3))))   12d6s21
                     jtest=imy(3)+2*(1-imy(3)+2*(1-imy(4)+2*imy(3)))     12d6s21
                     itestp=ltest+1                                      12d6s21
                     jtestp=jtest+1                                      12d6s21
                     phsk=-1d0                                            12d6s21
                     if(isopt(2).ne.0)phsk=-phsk                         12d9s21
                     iuse=-1                                             12d6s21
                     if(iifmx(itestp).ge.0)then                          12d6s21
                      iuse=iifmx(itestp)                                 12d6s21
                     else if(iifmx(jtestp).ge.0)then                     12d6s21
                      iuse=iifmx(jtestp)                                 12d6s21
                      if(isopt(4).ne.0)phsk=-phsk                        12d6s21
                     end if                                              12d6s21
                     if(iuse.ge.0)then                                   12d6s21
                      icolk=igya(3)+irefo(isy(3))*(igya(3)+irefo(isy(3)) 12d6s21
     $                   *(igya(3)+irefo(isy(3))*iuse))                 12d6s21
                      jdenj=id1vnotv(l,isy(3),isy(3))+nnl*icolk          12d6s21
                      ibc(nd1vnotv(l,isy(3),isy(3))+icolk)=1             12d6s21
                      jtmp=itmp                                              12d18s20
                      do iii=0,nl(l)-1                                         12d18s20
                       ii=iff22(iad1+iii)-1                              8d16s21
                       jdj=jdenj+ncsfb(jarg)*ii                             12d18s20
                       do j=0,ncsfb(jarg)-1                                   12d18s20
                        bc(jdj+j)=bc(jdj+j)+bc(jtmp+j)*phsk              12d6s21
                       end do                                                12d18s20
                       jtmp=jtmp+ncsfb(jarg)                             12d7s21
                      end do                                                 12d18s20
                     end if                                              12d6s21
                    end if                                               12d6s21
                    ibcoff=itmp                                          12d6s21
                   end if                                                12d6s21
                   jmat=jmat+ncsfb(jarg)*ncsfk2(l,iarg)                  12d6s21
                  end do                                                 12d7s21
                 end do
                 ibcoff=ioutg                                            12d6s21
                 do i=1,norb                                             12d9s21
                  if(btest(i1o,i).and.btest(j2o,i).and.l2e.eq.0)then     2d17s22
                   itestc=j2c                                            12d6s21
                   itesto=j2o                                            12d6s21
                   nopenk=nopen2p                                        12d6s21
c
c     anihilate common
c
                   if(btest(itestc,i))then                               12d9s21
                    itestc=ibclr(itestc,i)                               12d9s21
                    itesto=ibset(itesto,i)                               12d9s21
                    nopenk=nopenk+1                                             11d13s20
                   else                                                         11d13s20
                    itesto=ibclr(itesto,i)                               12d9s21
                    nopenk=nopenk-1                                             11d13s20
                   end if                                                       11d13s20
c
c     create bra
c
                   if(btest(itesto,nab4(1,1)))then                            12d9s21
                    itestc=ibset(itestc,nab4(1,1))                            12d9s21
                    itesto=ibclr(itesto,nab4(1,1))                            12d9s21
                    nopenk=nopenk-1                                             11d13s20
                   else                                                         11d13s20
                    itesto=ibset(itesto,nab4(1,1))                            12d9s21
                    nopenk=nopenk+1                                             11d13s20
                   end if                                                       11d13s20
                   call gandcr(i1c,i1o,itestc,itesto,nopen1,nopenk,      12d6s21
     $            norbxx,nnot1,nab1,icode,imap,nx1,irw1,irw2,iwpb1,     12d6s21
     $               iwpk1,bc,ibc)                                      11d14s22
                   call gandcr(itestc,itesto,j2c,j2o,nopenk,nopen2p,     12d6s21
     $            norbxx,nnot2,nab2,icode,imap,nx1,irw1,irw2,iwpb2,     12d6s21
     $                  iwpk2,bc,ibc)                                   11d14s22
                   if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                    call spinloop(i2sb,i2smb,i2sk,i2smk,nopen1,nopen2p,  12d6s21
     $               nopenk,ncsfb(jarg),ncsfk(iarg),itype,imatx,ntypeq, 11d22s21
     $                  nab1,iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(iunit),    11d22s21
     $                  ncsfk(iarg),ieoro,bc,ibc)                       11d14s22
                    do ixi=0,ntypeq-1                                     12d6s21
                     ipack=ibc(itype+ixi)                                12d7s21
                     imatu=imatx+ncsfk(iarg)*ncsfb(jarg)*ixi             12d7s21
                     do j=1,4                                               9d9s21
                      if(ipack2(j).gt.0)then                                9d9s21
                       ipack2a(j)=ipack2(j)                              12d6s21
                       imy(j)=0                                             10d13s21
                      else                                                  9d9s21
                       ipack2a(j)=-ipack2(j)                             12d6s21
                       imy(j)=1                                             10d13s21
                      end if                                                9d9s21
                      if(ipack2a(j).le.norb)then                         12d6s21
                       isy(j)=ism(ipack2a(j))                                9d9s21
                       igya(j)=irel(ipack2a(j))-1                             9d9s21
                      end if                                             12d6s21
                     end do                                                 9d9s21
                     if(ipack2a(2).gt.norb)then                          12d11s21
c
c     (1v|34)=(34|1v)=+/-(43|v1)
c
                      ltest=imy(4)+2*(imy(3)+2*(imy(2)+2*imy(1)))        12d11s21
                      jtest=1-imy(4)+2*(1-imy(3)+2*(1-imy(2)             12d11s21
     $                     +2*(1-imy(1))))                               12d11s21
                      itestp=ltest+1                                     12d7s21
                      jtestp=jtest+1                                     12d7s21
                      juse=-1                                            12d7s21
                      phsj=1d0                                           12d7s21
                      if(isopt(2).ne.0)phsj=-phsj                        12d9s21
                      if(iifmx(itestp).ge.0)then                         12d7s21
                       juse=iifmx(itestp)                                12d7s21
                      else if(iifmx(jtestp).ge.0)then                    12d7s21
                       juse=iifmx(jtestp)                                12d7s21
                       if(isopt(4).ne.0)phsj=-phsj                       12d7s21
                      end if                                             12d7s21
c
c     (14|3v)=+/-(41|v3)
c
                      ltest=imy(4)+2*(imy(1)+2*(imy(2)+2*imy(3)))        12d11s21
                      jtest=1-imy(4)+2*(1-imy(1)+2*(1-imy(2)             12d11s21
     $                    +2*(1-imy(3))))                               12d11s21
                      itestp=ltest+1                                     12d7s21
                      jtestp=jtest+1                                     12d7s21
                      kuse=-1                                            12d7s21
                      phsk=-1d0                                           12d7s21
                      if(isopt(2).ne.0)phsk=-phsk                        12d11s21
                      if(iifmx(itestp).ge.0)then                         12d7s21
                       kuse=iifmx(itestp)                                12d7s21
                      else if(iifmx(jtestp).ge.0)then                    12d7s21
                       kuse=iifmx(jtestp)                                12d7s21
                       if(isopt(4).ne.0)phsk=-phsk                       12d7s21
                      end if                                             12d7s21
                      if(max(juse,kuse).ge.0)then                        12d7s21
                       jcol=igya(4)+irefo(isy(4))*(igya(3)+irefo(isy(3)) 12d11s21
     $                     *(igya(1)+irefo(isy(1))*juse))               12d11s21
                       njkcol=irefo(isy(4))*irefo(isy(3))*irefo(isy(1))
     $                     *ntype
                       kcol=igya(4)+irefo(isy(4))*(igya(1)+irefo(isy(1)) 12d11s21
     $                     *(igya(3)+irefo(isy(3))*kuse))               12d11s21
                       jmat=imatu                                        12d11s21
                       do l=1,4
                        if(nl(l).gt.0)then                                     12d18s20
                         nnl=ncsfb(jarg)*nfdat(2,l,isb)                  12d13s21
                         iad1=jvcv+iff22(jvcv+5+l)                        12d7s21
                         iad2=iad1+nl(l)                                  3d19s21
                         itmp=ibcoff                                           12d18s20
                         ibcoff=itmp+ncsfb(jarg)*nl(l)                   12d11s21
                         call enough('hcsdbk4.  8',bc,ibc)
                         call dgemm('n','n',ncsfb(jarg),nl(l),           12d11s21
     $                      ncsfk2(l,iarg),1d0,bc(jmat),ncsfb(jarg),    12d11s21
     $                      ff22(iad2),ncsfk2(l,iarg),0d0,bc(itmp),     12d11s21
     $                      ncsfb(jarg),                                12d11s21
     d' hcsdbk4.  2')
                         if(juse.ge.0)then                                12d7s21
                          jden=id1vnotv(l,isy(4),isy(3))+nnl*jcol        12d11s21
                          ibc(nd1vnotv(l,isy(4),isy(3))+jcol)=1          12d11s21
                          jtmp=itmp                                              12d18s20
                          mden=mdhvnotv(l)                                 12d21s20
                          do iii=0,nl(l)-1                                         12d18s20
                           ii=iff22(iad1+iii)-1                           8d16s21
                           jdh=jden+ncsfb(jarg)*ii                       12d11s21
                           ibc(mden+ii)=1                                  12d19s20
                           do j=0,ncsfb(jarg)-1                          12d11s21
                            bc(jdh+j)=bc(jdh+j)+bc(jtmp+j)*phsj           12d7s21
                           end do                                                12d18s20
                           jtmp=jtmp+ncsfb(jarg)                         12d11s21
                          end do                                                 12d18s20
                         end if                                           12d7s21
                         if(kuse.ge.0)then                                12d7s21
                          jden=id1vnotv(l,isy(4),isy(1))+nnl*kcol        12d11s21
                          ibc(nd1vnotv(l,isy(4),isy(1))+kcol)=1          12d11s21
                          jtmp=itmp                                              12d18s20
                          mden=mdhvnotv(l)                                 12d21s20
                          do iii=0,nl(l)-1                                         12d18s20
                           ii=iff22(iad1+iii)-1                           8d16s21
                           jdh=jden+ncsfb(jarg)*ii                       12d11s21
                           ibc(mden+ii)=1                                  12d19s20
                           do j=0,ncsfb(jarg)-1                          12d11s21
                            orig=bc(jdh+j)
                            bc(jdh+j)=bc(jdh+j)+bc(jtmp+j)*phsk           12d7s21
                           end do                                                12d18s20
                           jtmp=jtmp+ncsfb(jarg)                         12d11s21
                          end do                                                 12d18s20
                         end if                                           12d7s21
                         ibcoff=itmp                                      12d19s20
                        end if                                              12d19s20
                        jmat=jmat+ncsfb(jarg)*ncsfk2(l,iarg)             12d11s21
                       end do                                             12d19s20
                      end if                                             12d7s21
                     else                                                12d6s21
                      write(6,*)('don''t know how to handle ipack2! '),  12d6s21
     $                    ipack2                                        12d6s21
                      stop 'hcsdbk4'                                     12d6s21
                     end if                                              12d6s21
                    end do                                               12d6s21
                    ibcoff=itype                                         12d6s21
                   end if                                                12d6s21
                  end if                                                 12d7s21
                 end do                                                  12d18s20
                 ibcoff=iunit                                            12d7s21
                else if(l2e.eq.0)then                                   2d8s23
                 nnot=0                                                 2d7s23
                 if(ndifs.eq.4.and.ndifb.eq.4)then                      2d7s23
                  nnot=4                                                2d7s23
                  ioxx(1)=1                                             2d7s23
                  ioxx(2)=1                                             2d7s23
                  do i=1,norbxx                                         2d7s23
                   if(btest(gandcb,i))then                              2d7s23
                    if((btest(j2c,i).and.btest(i1o,i)).or.              2d8s23
     $                (btest(j2o,i).and..not.btest(i1c,i)))then         2d8s23
                     nab4(2,ioxx(2))=i                                  2d7s23
                     ioxx(2)=ioxx(2)+1                                  2d7s23
                    else                                                2d7s23
                     nab4(1,ioxx(1))=i                                  2d7s23
                     ioxx(1)=ioxx(1)+1                                  2d7s23
                    end if                                              2d7s23
                   end if                                               2d7s23
                  end do                                                2d7s23
                 else if(ndifb.eq.3)then                                2d7s23
                  nnot=3                                                2d7s23
                  ioxx(1)=1                                             2d7s23
                  ioxx(2)=1                                             2d7s23
                  iswap=0                                               2d7s23
                  do i=1,norbxx                                         2d7s23
                   if(btest(gandcb,i))then                              2d7s23
                    if(btest(gandcc,i).and.                             2d7s23
     $        ((btest(i1c,i).and..not.btest(j2o,i)).or.                 2d7s23
     $         (btest(j2c,i).and..not.btest(i1o,i))))then               2d7s23
                     if(btest(j2c,i))iswap=1                            2d7s23
                     nab4(1,1)=i                                        2d7s23
                     nab4(1,2)=i                                        2d7s23
                    else                                                2d7s23
                     nab4(2,ioxx(2))=i                                  2d7s23
                     ioxx(2)=ioxx(2)+1                                  2d7s23
                    end if                                              2d7s23
                   end if                                               2d7s23
                  end do                                                2d7s23
                  if(iswap.ne.0)then                                    2d7s23
                   icpy=nab4(1,1)                                       2d7s23
                   nab4(1,1)=nab4(2,1)                                  2d7s23
                   nab4(2,1)=icpy                                       2d7s23
                   icpy=nab4(1,2)                                       2d7s23
                   nab4(1,2)=nab4(2,2)                                  2d7s23
                   nab4(2,2)=icpy                                       2d7s23
                   nbt=0                                                2d7s23
                   if(btest(j2c,nab4(2,2)).and.                         2d7s23
     $                  .not.btest(j2c,nab4(2,1)))nbt=1                 2d7s23
                  else                                                  2d7s23
                   nbt=0                                                2d7s23
                   if(btest(i1c,nab4(1,2)).and.                         2d8s23
     $                  .not.btest(i1c,nab4(1,1)))nbt=1                 2d8s23
                  end if                                                2d7s23
                  if(nbt.ne.0)then                                      2d7s23
                   nab4(1,1)=nab4(1,2)                                  2d7s23
                   nab4(2,1)=nab4(2,2)                                  2d7s23
                  end if                                                2d7s23
                 else if(ndifs.eq.0.and.ndifd.eq.2)then                 2d7s23
                  nnot=3                                                2d7s23
                  do i=1,norbxx                                         2d7s23
                   if(btest(gandcb,i))then                              2d7s23
                    if(btest(i1c,i))then                                2d8s23
                     nab4(1,1)=i                                        2d7s23
                     nab4(1,2)=i                                        2d7s23
                    else                                                2d7s23
                     nab4(2,1)=i                                        2d7s23
                     nab4(2,2)=i                                        2d7s23
                    end if                                              2d7s23
                   end if                                               2d7s23
                  end do                                                2d7s23
                 end if                                                 2d7s23
                 if(nnot.ne.0)then                                      2d7s23
                  iunit=ibcoff                                            12d6s21
                  ibcoff=iunit+ncsfk(iarg)*ncsfk(iarg)                    12d6s21
                  call enough('hcsdbk4.  9',bc,ibc)
                  do iz=iunit,ibcoff-1                                    12d6s21
                   bc(iz)=0d0                                             12d6s21
                  end do                                                  12d6s21
                  do iz=0,ncsfk(iarg)-1                                   12d6s21
                   iad=iunit+iz*(ncsfk(iarg)+1)                           12d6s21
                   bc(iad)=1d0                                            12d6s21
                  end do                                                  12d6s21
                  iu1=1
                  iu2=1
                  itestc=i1c                                              12d7s21
                  itesto=i1o                                              12d7s21
                  if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                   itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                   itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopen1+1                                              11d13s20
                  else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                   itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopen1-1                                              11d13s20
                  end if                                                        11d13s20
                  if(btest(itesto,nab4(2,iu2)))then                         12d8s20
                   itestc=ibset(itestc,nab4(2,iu2))                         12d8s20
                   itesto=ibclr(itesto,nab4(2,iu2))                             11d27s20
                   nopenk=nopenk-1                                              11d13s20
                  else                                                          11d13s20
                   itesto=ibset(itesto,nab4(2,iu2))                         12d8s20
                   nopenk=nopenk+1                                              11d13s20
                  end if                                                        11d13s20
                  call gandcr(i1c,i1o,itestc,itesto,nopen1,nopenk,        12d6s21
     $            norbxx,nnot1,nab1,icode,imap,nx1,irw1,irw2,iwpb1,     12d6s21
     $               iwpk1,bc,ibc)                                      11d14s22
                  call gandcr(itestc,itesto,j2c,j2o,nopenk,nopen2p,       12d6s21
     $            norbxx,nnot2,nab2,icode,imap,nx1,irw1,irw2,iwpb2,     12d6s21
     $                  iwpk2,bc,ibc)                                   11d14s22
                  if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                   call spinloop(i2sb,i2smb,i2sk,i2smk,nopen1,nopen2p,    12d6s21
     $               nopenk,ncsfb(jarg),ncsfk(iarg),itype,imatx,ntypeq, 11d22s21
     $                  nab1,iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(iunit),    11d22s21
     $                  ncsfk(iarg),ieoro,bc,ibc)                       11d14s22
                   do ixi=0,ntypeq-1                                      12d6s21
                    ipack=ibc(itype+ixi)                                  12d6s21
                    imatu=imatx+ncsfb(jarg)*ncsfk(iarg)*ixi               12d7s21
                    do j=1,4                                               9d9s21
                     if(ipack2(j).gt.0)then                                9d9s21
                      ipack2a(j)=ipack2(j)                                12d6s21
                      imy(j)=0                                             10d13s21
                     else                                                  9d9s21
                      ipack2a(j)=-ipack2(j)                               12d6s21
                      imy(j)=1                                             10d13s21
                     end if                                                9d9s21
                     if(ipack2a(j).le.norb)then                           12d6s21
                      isy(j)=ism(ipack2a(j))                                9d9s21
                      igya(j)=irel(ipack2a(j))-1                             9d9s21
                     end if                                               12d6s21
                    end do                                                 9d9s21
                    if(ipack2a(4).gt.norb)then                            12d11s21
c
c     (12|3v)=+/-(21|v3)
c
                     ltest=imy(2)+2*(imy(1)+2*(imy(4)+2*imy(3)))          12d11s21
                     jtest=1-imy(2)+2*(1-imy(1)+2*(1-imy(4)             2d7s23
     $                    +2*(1-imy(3))))                               2d7s23
                     itestp=ltest+1                                       12d7s21
                     jtestp=jtest+1                                       12d7s21
                     juse=-1                                              12d7s21
                     phsj=1d0                                             12d7s21
                     if(isopt(2).ne.0)phsj=-phsj                          12d11s21
                     if(iifmx(itestp).ge.0)then                           12d7s21
                      juse=iifmx(itestp)                                  12d7s21
                     else if(iifmx(jtestp).ge.0)then                      12d7s21
                      juse=iifmx(jtestp)                                  12d7s21
                      if(isopt(4).ne.0)phsj=-phsj                         12d7s21
                     end if                                               12d7s21
                     kuse=-1                                              12d7s21
                     if(nnot.eq.4)then                                    12d7s21
c
c     (1v|32)=(32|1v)=+/-(23|v1)                                        12d11s21
c
                      ltest=imy(2)+2*(imy(3)+2*(imy(4)+2*imy(1)))         12d11s21
                      jtest=1-imy(2)+2*(1-imy(3)+2*(1-imy(4)              12d11s21
     $                    +2*(1-imy(1))))                               12d11s21
                      itestp=ltest+1                                      12d7s21
                      jtestp=jtest+1                                      12d7s21
                      phsk=-1d0                                           12d7s21
                      if(isopt(2).ne.0)phsk=-phsk                         12d11s21
                      if(iifmx(itestp).ge.0)then                          12d7s21
                       kuse=iifmx(itestp)                                 12d7s21
                      else if(iifmx(jtestp).ge.0)then                     12d7s21
                       kuse=iifmx(jtestp)                                 12d7s21
                       if(isopt(4).ne.0)phsk=-phsk                        12d7s21
                      end if                                              12d7s21
                     end if                                               12d7s21
                     if(max(juse,kuse).ge.0)then                          12d7s21
                      jcol=igya(2)+irefo(isy(2))*(igya(1)+irefo(isy(1))   12d11s21
     $                     *(igya(3)+irefo(isy(3))*juse))               12d11s21
                      kcol=igya(2)+irefo(isy(2))*(igya(3)+irefo(isy(3))   12d11s21
     $                     *(igya(1)+irefo(isy(1))*kuse))               12d11s21
                      jmat=imatu                                          12d11s21
                      do l=1,4
                       if(nl(l).gt.0)then                                     12d18s20
                        nnl=ncsfb(jarg)*nfdat(2,l,isb)                         12d19s20
                        iad1=jvcv+iff22(jvcv+5+l)                         12d7s21
                        iad2=iad1+nl(l)                                   12d7s21
                        itmp=ibcoff                                           12d18s20
                        ibcoff=itmp+ncsfb(jarg)*nl(l)                          12d18s20
                        call enough('hcsdbk4. 10',bc,ibc)
                        call dgemm('n','n',ncsfb(jarg),nl(l),             12d7s21
     $                      ncsfk2(l,iarg),1d0,bc(jmat),ncsfb(jarg),     3d19s21
     $                      ff22(iad2),ncsfk2(l,iarg),0d0,bc(itmp),      8d16s21
     $                      ncsfb(jarg),                                 3d19s21
     d' hcsdbk4.  3')
                        if(juse.ge.0)then                                 12d7s21
                         jden=id1vnotv(l,isy(2),isy(1))+nnl*jcol          12d11s21
                         ibc(nd1vnotv(l,isy(2),isy(1))+jcol)=1            12d11s21
                         jtmp=itmp                                              12d18s20
                         mden=mdhvnotv(l)                                 12d21s20
                         do iii=0,nl(l)-1                                         12d18s20
                          ii=iff22(iad1+iii)-1                            12d7s21
                          jdh=jden+ncsfb(jarg)*ii                         12d11s21
                          ibc(mden+ii)=1                                  12d19s20
                          do j=0,ncsfb(jarg)-1                            12d11s21
                           bc(jdh+j)=bc(jdh+j)+bc(jtmp+j)*phsj            12d7s21
                          end do                                                12d18s20
                          jtmp=jtmp+ncsfb(jarg)                           12d11s21
                         end do                                                 12d18s20
                        end if                                            12d7s21
                        if(kuse.ge.0)then                                 12d7s21
                         jden=id1vnotv(l,isy(2),isy(3))+nnl*kcol          12d11s21
                         ibc(nd1vnotv(l,isy(2),isy(3))+kcol)=1            12d11s21
                         jtmp=itmp                                              12d18s20
                         mden=mdhvnotv(l)                                 12d21s20
                         do iii=0,nl(l)-1                                         12d18s20
                          ii=iff22(iad1+iii)-1                            12d7s21
                          jdh=jden+ncsfb(jarg)*ii                         12d11s21
                          ibc(mden+ii)=1                                  12d19s20
                          do j=0,ncsfb(jarg)-1                            12d11s21
                           bc(jdh+j)=bc(jdh+j)+bc(jtmp+j)*phsk            12d7s21
                          end do                                                12d18s20
                          jtmp=jtmp+ncsfb(jarg)                           12d11s21
                         end do                                                 12d18s20
                        end if                                            12d7s21
                        ibcoff=itmp                                       12d7s21
                       end if                                              12d19s20
                       jmat=jmat+ncsfb(jarg)*ncsfk2(l,iarg)               12d7s21
                      end do                                              12d7s21
                     end if                                               12d7s21
                    else                                                  12d7s21
                     write(6,*)('don''t know how to handle ipack2! '),    12d7s21
     $                    ipack2                                        12d6s21
                     stop 'hcsdbk4'                                       12d7s21
                    end if                                                12d7s21
                   end do                                                 12d7s21
                   ibcoff=itype                                           12d7s21
                  end if                                                  12d7s21
                 end if                                                 2d7s23
                end if                                                   12d18s20
               end if                                                   2d7s23
               if(l2e.eq.0)then                                         2d7s23
                j2o=ibclr(j2o,norbx)                                      12d18s20
                j2o=ibset(j2o,norbxxx)                                     12d18s20
                gandcc=ieor(i1c,j2c)                                      10d13s22
                gandco=ieor(i1o,j2o)                                      10d13s22
                gandcb=ior(gandcc,gandco)                                 10d20s22
                ndifb=popcnt(gandcb)                                      10d20s22
                if(ndifb.le.4)then                                        10d20s22
                 ndifs=popcnt(gandco)                                      10d13s22
                 ndifd=popcnt(gandcc)                                      10d13s22
                 nnot=0                                                 2d7s23
                 if(ndifs.eq.4.and.ndifb.eq.4)then                      2d7s23
                  nnot=4                                                2d7s23
                  ioxx(1)=1                                             2d7s23
                  ioxx(2)=1                                             2d7s23
                  do i=1,norbxx                                         2d7s23
                   if(btest(gandcb,i))then                              2d7s23
                    if((btest(j2c,i).and.btest(i1o,i)).or.              2d8s23
     $                (btest(j2o,i).and..not.btest(i1c,i)))then         2d8s23
                     nab4(2,ioxx(2))=i                                  2d7s23
                     ioxx(2)=ioxx(2)+1                                  2d7s23
                    else                                                2d7s23
                     nab4(1,ioxx(1))=i                                  2d7s23
                     ioxx(1)=ioxx(1)+1                                  2d7s23
                    end if                                              2d7s23
                   end if                                               2d7s23
                  end do                                                2d7s23
                 else if(ndifb.eq.3)then                                2d7s23
                  nnot=3                                                2d7s23
                  ioxx(1)=1                                             2d7s23
                  ioxx(2)=1                                             2d7s23
                  iswap=0                                               2d7s23
                  do i=1,norbxx                                         2d7s23
                   if(btest(gandcb,i))then                              2d7s23
                    if(btest(gandcc,i).and.                             2d7s23
     $        ((btest(i1c,i).and..not.btest(j2o,i)).or.                 2d8s23
     $         (btest(j2c,i).and..not.btest(i1o,i))))then               2d8s23
                     if(btest(j2c,i))iswap=1                            2d7s23
                     nab4(1,1)=i                                        2d7s23
                     nab4(1,2)=i                                        2d7s23
                    else                                                2d7s23
                     nab4(2,ioxx(2))=i                                  2d7s23
                     ioxx(2)=ioxx(2)+1                                  2d7s23
                    end if                                              2d7s23
                   end if                                               2d7s23
                  end do                                                2d7s23
                  if(iswap.ne.0)then                                    2d7s23
                   icpy=nab4(1,1)                                       2d7s23
                   nab4(1,1)=nab4(2,1)                                  2d7s23
                   nab4(2,1)=icpy                                       2d7s23
                   icpy=nab4(1,2)                                       2d7s23
                   nab4(1,2)=nab4(2,2)                                  2d7s23
                   nab4(2,2)=icpy                                       2d7s23
                   nbt=0                                                2d7s23
                   if(btest(j2c,nab4(2,2)).and.                         2d7s23
     $                  .not.btest(j2c,nab4(2,1)))nbt=1                 2d7s23
                  else                                                  2d7s23
                   nbt=0                                                2d7s23
                   if(btest(i1c,nab4(1,2)).and.                         2d7s23
     $                  .not.btest(i1c,nab4(1,1)))nbt=1                 2d7s23
                  end if                                                2d7s23
                  if(nbt.ne.0)then                                      2d7s23
                   nab4(1,1)=nab4(1,2)                                  2d7s23
                   nab4(2,1)=nab4(2,2)                                  2d7s23
                  end if                                                2d7s23
                 else if(ndifs.eq.0.and.ndifd.eq.2)then                 2d7s23
                  nnot=3                                                2d7s23
                  do i=1,norbxxx                                        2d7s23
                   if(btest(gandcb,i))then                              2d7s23
                    if(btest(i1c,i))then                                2d7s23
                     nab4(1,1)=i                                        2d7s23
                     nab4(1,2)=i                                        2d7s23
                    else                                                2d7s23
                     nab4(2,1)=i                                        2d7s23
                     nab4(2,2)=i                                        2d7s23
                    end if                                              2d7s23
                   end if                                               2d7s23
                  end do                                                2d7s23
                 end if                                                 2d7s23
                 if(nnot.eq.4)then                                      2d7s23
                  iunit=ibcoff                                            12d6s21
                  ibcoff=iunit+ncsfk(iarg)*ncsfk(iarg)                    12d6s21
                  call enough('hcsdbk4. 11',bc,ibc)
                  do iz=iunit,ibcoff-1                                    12d6s21
                   bc(iz)=0d0                                             12d6s21
                  end do                                                  12d6s21
                  do iz=0,ncsfk(iarg)-1                                   12d6s21
                   iad=iunit+iz*(ncsfk(iarg)+1)                           12d6s21
                   bc(iad)=1d0                                            12d6s21
                  end do                                                  12d6s21
                  iu1=1
                  iu2=1
                  itestc=i1c                                              12d11s21
                  itesto=i1o                                              12d11s21
                  if(btest(itestc,nab4(1,iu1)))then                             11d27s20
                   itestc=ibclr(itestc,nab4(1,iu1))                             11d27s20
                   itesto=ibset(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopen1+1                                        12d11s21
                  else if(btest(itesto,nab4(1,iu1)))then                        11d27s20
                   itesto=ibclr(itesto,nab4(1,iu1))                             11d27s20
                   nopenk=nopen1-1                                              11d13s20
                  end if                                                        11d13s20
                  if(btest(itesto,nab4(2,iu2)))then                         12d8s20
                   itestc=ibset(itestc,nab4(2,iu2))                         12d8s20
                   itesto=ibclr(itesto,nab4(2,iu2))                             11d27s20
                   nopenk=nopenk-1                                              11d13s20
                  else                                                          11d13s20
                   itesto=ibset(itesto,nab4(2,iu2))                         12d8s20
                   nopenk=nopenk+1                                              11d13s20
                  end if                                                        11d13s20
                  call gandcr(i1c,i1o,itestc,itesto,nopen1,nopenk,        12d6s21
     $            norbxxx,nnot1,nab1,icode,imap,nx1,irw1,irw2,iwpb1,     12d6s21
     $               iwpk1,bc,ibc)                                      11d14s22
                  call gandcr(itestc,itesto,j2c,j2o,nopenk,nopen2p,       12d7s21
     $              norbxxx,nnot2,nab2,icode,imap,nx1,irw1,irw2,iwpb2,     12d6s21
     $                  iwpk2,bc,ibc)                                   11d14s22
                  if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                   call spinloop(i2sb,i2smb,i2sk,i2smk,nopen1,nopen2p,    12d6s21
     $               nopenk,ncsfb(jarg),ncsfk(iarg),itype,imatx,ntypeq, 11d22s21
     $                  nab1,iwpb1,iwpk1,nab2,iwpb2,iwpk2,bc(iunit),    11d22s21
     $                  ncsfk(iarg),ieoro,bc,ibc)                       11d14s22
                   do ixi=0,ntypeq-1                                       12d6s21
                    ipack=ibc(itype+ixi)                                   12d6s21
                    imatu=imatx+ncsfk(iarg)*ncsfb(jarg)*ixi                12d7s21
                    do j=1,4                                               9d9s21
                     if(ipack2(j).gt.0)then                                9d9s21
                      ipack2a(j)=ipack2(j)                                12d6s21
                      imy(j)=0                                              10d13s21
                     else                                                  9d9s21
                      ipack2a(j)=-ipack2(j)                               12d6s21
                      imy(j)=1                                             10d13s21
                     end if                                                9d9s21
                     if(ipack2a(j).le.norb)then                           12d6s21
                      isy(j)=ism(ipack2a(j))                                9d9s21
                      igya(j)=irel(ipack2a(j))-1                             9d9s21
                     end if                                               12d6s21
                    end do                                                 9d9s21
                    if(ipack2a(1).le.norb)then                            12d7s21
c
c     (v"v|v'n)Vv'v" is bmat
c but 3x is stored badc, indexed b,c,d
c                ba d c
c       d  sd'   sd'  d      sd   d'
c     (1v'|vv")=(vv"|1v') & (vv'|1v")
c      12  34    34  12      32  14
c     (1v'|vv")=(vv"|1v')=+/-(v"v|v'1)
c      12  34    34  12       4 3 2 1
c
                     ltest=imy(4)+2*(imy(3)+2*(imy(2)+2*imy(1)))          12d7s21
                     jtest=1-imy(4)+2*(1-imy(3)+2*(1-imy(2)             2d7s23
     $                    +2*(1-imy(1))))                               2d7s23
                     itestp=ltest+1                                       12d7s21
                     jtestp=jtest+1                                       12d7s21
                     iuse=-1                                              12d7s21
                     phs=1d0                                              12d7s21
                     if(isopt(2).ne.0)phs=-phs                            12d11s21
                     if(iifmx(itestp).ge.0)then                           12d7s21
                      iuse=iifmx(itestp)                                  12d7s21
                     else if(iifmx(jtestp).ge.0)then                      12d7s21
                      iuse=iifmx(jtestp)                                  12d7s21
                      if(isopt(4).ne.0)phs=-phs                           12d7s21
                     end if                                               12d7s21
c     (1v"|vv')=(vv'|1v")=+/-(v'v|v"1)
c      14  32    32  14       2 3 4 1
                     kuse=-1                                              12d9s21
                     phsk=-1d0                                            12d9s21
                     if(isopt(2).ne.0)phsk=-phsk                          12d11s21
                     ltest=imy(2)+2*(imy(3)+2*(imy(4)+2*imy(1)))          12d11s21
                     jtest=1-imy(2)+2*(1-imy(3)+2*(1-imy(4)             2d7s23
     $                    +2*(1-imy(1))))                               2d7s23
                     itestp=ltest+1                                       12d7s21
                     jtestp=jtest+1                                       12d7s21
                     if(iifmx(itestp).ge.0)then                           12d7s21
                      kuse=iifmx(itestp)                                  12d7s21
                     else if(iifmx(jtestp).ge.0)then                      12d7s21
                      kuse=iifmx(jtestp)                                  12d7s21
                      if(isopt(4).ne.0)phsk=-phsk                         12d9s21
                     end if                                               12d7s21
                     if(iuse.ge.0)then                                    12d7s21
                      icol=igya(1)+irefo(isy(1))*iuse                     12d11s21
                      jmat=imatu                                          12d7s21
                      do l=1,4                                               12d19s20
                       if(nl(l).gt.0)then                                     12d18s20
                        nnl=ncsfb(jarg)*nfdat(2,l,isb)                    12d11s21
                        iad1=jvcv+iff22(jvcv+5+l)                         8d16s21
                        iad2=iad1+nl(l)                                    3d19s21
                        itmp=ibcoff                                           12d18s20
                        ibcoff=itmp+ncsfb(jarg)*nl(l)                     12d11s21
                        call enough('hcsdbk4. 12',bc,ibc)
                        call dgemm('n','n',ncsfb(jarg),nl(l),             12d11s21
     $                     ncsfk2(l,iarg),                              12d11s21
     $            1d0,bc(jmat),ncsfb(jarg),ff22(iad2),ncsfk2(l,iarg),   12d11s21
     $                     0d0,bc(itmp),ncsfb(jarg),                    12d11s21
     d' hcsdbk4.  4')
                        jden=id3vnotv3(l)+nnl*icol                        12d7s21
                        ibc(nd3vnotv3(l)+icol)=1                          12d22s20
                        mden=md3vnotv3(l)                                 12d7s21
                        jtmp=itmp                                              12d18s20
                        do iii=0,nl(l)-1                                         12d18s20
                         ii=iff22(iad1+iii)-1                             12d7s21
                         jdh=jden+ncsfb(jarg)*ii                          12d11s21
                         ibc(mden+ii)=1                                      12d19s20
                         do j=0,ncsfb(jarg)-1                             12d11s21
                          bc(jdh+j)=bc(jdh+j)+bc(jtmp+j)*phs                       12d18s20
                         end do                                                12d18s20
                         jtmp=jtmp+ncsfb(jarg)                            12d11s21
                        end do                                                 12d18s20
                        ibcoff=itmp                                          12d19s20
                       end if                                                12d19s20
                       jmat=jmat+ncsfb(jarg)*ncsfk2(l,iarg)               12d11s21
                      end do                                               12d7s21
                     end if                                               12d7s21
                    else                                                  12d7s21
                     write(6,*)('I don''t know how to deal with '),
     $                  ipack2
                     stop 'hcsdbk4'
                    end if                                                12d7s21
                   end do                                                 12d7s21
                   ibcoff=iunit                                           12d7s21
                  end if
                 end if                                                 2d7s23
                end if                                                  2d7s23
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
                 nokdc(isd,isc,l)=0                                     12d21s20
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
                 nnl=ncsfb(jarg)*nfdat(2,l,isb)                            12d19s20
                 nok4(l)=0                                               12d21s20
                 jdhvnotv=idhvnotv(l)                                    12d21s20
                 jdhvnotvf=idhvnotvf(l)                                  1d6s21
                 do i=0,irefo(lgoal)-1                                   12d21s20
                  if(ibc(ndhvnotv(l)+i).ne.0)then                        12d21s20
                   ibc(ndhvnotv(l)+nok4(l))=i                            12d21s20
                   nok4(l)=nok4(l)+1                                     12d21s20
                   iad=idhvnotv(l)+nnl*i
                   do k=0,nok4f(l)-1                                     12d21s20
                    iad=idhvnotv(l)+ncsfb(jarg)*ibc(mdhvnotv(l)+k)+nnl*i  1d6s21
                    do j=0,ncsfb(jarg)-1                                  1d6s21
                     bc(jdhvnotvf+j)=bc(iad+j)                           1d6s21
                    end do                                               12d21s20
                    jdhvnotvf=jdhvnotvf+ncsfb(jarg)                       1d6s21
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
                   nnn=nn*irefo(iscdv)*ntype                            12d9s21
                   nokdc(isd,isc,l)=0                                     12d21s20
                   jd1vnotv=id1vnotv(l,isd,isc)                           12d21s20
                   do i=0,nnn-1
                    if(ibc(nd1vnotv(l,isd,isc)+i).ne.0.and.l2e.eq.0)then2d17s22
                     icol=i                                             12d29s20
                     ibc(nd1vnotv(l,isd,isc)+nokdc(isd,isc,l))=icol      12d29s20
                     nokdc(isd,isc,l)=nokdc(isd,isc,l)+1                 12d21s20
                     do k=0,nok4f(l)-1                                    12d21s20
                      iad=id1vnotv(l,isd,isc)+i11s-1                      12d21s20
     $                   +ncsfb(jarg)*ibc(mdhvnotv(l)+k)+nnl*i           12d21s20
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
             if(itransgg.eq.0)then                                      1d30s21
              call ddi_done(ibc(iacc),nacc)                             1d30s21
              do i=0,ncolt-1                                            1d30s21
               bc(igg+i)=0d0                                            1d30s21
              end do                                                    1d30s21
              itransgg=1                                                1d30s21
             end if                                                     1d30s21
             nok=0                                                      2d25s22
             if(njhere.gt.0.and.l2e.eq.0)then                           2d17s22
              ipass=1                                                   3d2s21
              jdkeep=idkeep(ipass)                                      1d4s21
              do i=0,irefo(lgoal3)*ntype-1                              12d9s21
               do l=1,4                                                 1d4s21
                if(nfdat(2,l,isb).gt.0)then                             1d4s21
                 if(ibc(nd3vnotv3(l)+i).ne.0)then                       12d9s21
                  do k=0,nfdat(2,l,isb)-1                               1d4s21
                   if(ibc(md3vnotv3(l)+k).ne.0)then                     12d9s21
                    ibc(ndkeep(ipass)+nok)=i                            1d5s21
                    ibc(mdkeep(ipass)+nok)=k+loff(l)                    1d5s21
                    iad=id3vnotv3(l)+ncsfb(jarg)*(k                     12d9s21
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
             end if                                                     8d15s21
             if(nok.gt.0.and.njhere.gt.0)then                           12s15s21
              nnb=nvirt(jsbv)*nrootu                                    8d13s21
              itmpb=ibcoff                                              8d26s21
              itmpd=itmpb+nnb*nok                                       2d3s21
              ibcoff=itmpd+njhere*nok                                   3d2s21
              call enough('hcsdbk4. 13',bc,ibc)
              jtmpb=itmpb                                               2d3s21
              do i=0,nok-1                                              2d3s21
               do j=0,njhere-1                                          3d2s21
                ji=idkeep(ipass)+j+njhere*i                             3d2s21
                ij=itmpd+i+nok*j                                        2d3s21
                bc(ij)=bc(ji)                                           2d3s21
               end do                                                   2d3s21
               in=ibc(ndkeep(ipass)+i)                                  2d3s21
               k=ibc(mdkeep(ipass)+i)                                   2d3s21
               jbmat=ibmat(isb)+nnb*(k+nfh*in)                          8d13s21
               do j=0,nnb-1                                             8d13s21
                bc(jtmpb+j)=bc(jbmat+j)                                 2d3s21
               end do                                                   2d3s21
               jtmpb=jtmpb+nnb                                          8d13s21
              end do                                                    2d3s21
              jggs=jgs+nnb*(i11s-1)                                     8d13s21
              call dgemm('n','n',nnb,njhere,nok,1d0,bc(itmpb),nnb,      8d13s21
     $               bc(itmpd),nok,1d0,bc(jggs),nnb,                    8d13s21
     d' hcsdbk4.  5')
              ibcoff=itmpb                                              8d26s21
             end if                                                     2d3s21
             do isbv1=1,nsymb                                           12d21s20
              ksbv2=multh(isbv1,ksbv12)                                 12d21s20
c     D is ket, S is bra
              ksbv1=isbv1                                               8d16s21
              if(ksbv1.le.ksbv2.and.njhere.gt.0)then                    12s15s21
               if(ksbv12.eq.1.and.jsbv.eq.ksbv1)then                    1d6s21
                nokf=nok4f(1)                                           2d25s21
                nrow=njhere*nokf                                          12d21s20
                intden=ibcoff                                             12d21s20
                ibcoff=intden+nrow*nvirt(jsbv)                          8d13s21
                call enough('hcsdbk4. 14',bc,ibc)
                fact=0d0                                                  12d21s20
                if(min(nok4(1),nok4f(1)).gt.0)then                      2d25s21
                 nok=nok4(1)                                            2d25s21
                 nokf=nok4f(1)                                          2d25s21
                 itmp=ibcoff                                              12d21s20
                 ibcoff=itmp+nok*nvirt(jsbv)                              12d21s20
                 call enough('hcsdbk4. 15',bc,ibc)
                 jtmp=itmp                                                12d21s20
                 iosym=multh(jsbv,isopt(1))                             12d9s21
c
c     the v here is for the double
c     isymbra=jsb*jsbv=jsb*jv=isb*io*jv
c     isymket=isb*isbv12=isb*jv*iv
c     isymbra*isymket=isopt(1)
c     isopt(1)=isb*io*jv*isb*jv*iv=io*iv
c
                 if(ih0n(iosym).gt.0)then                               2d24s22
                  do iv=0,nvirt(jsbv)-1                                 2d24s22
                   ivp=iv+irefo(jsbv)                                   2d24s22
                   do i=0,nok-1                                         2d24s22
                    iadh=ih0n(iosym)+ibc(ndhvnotv(1)+i)                 2d24s22
     $                   +nh0(iosym)*ivp                                2d24s22
                    bc(jtmp+i)=bc(iadh)                                 2d24s22
                   end do                                               2d24s22
                   jtmp=jtmp+nok                                        2d24s22
                  end do                                                2d24s22
                 else if(ih0n(jsbv).gt.0)then                           2d24s22
                  do iv=0,nvirt(jsbv)-1                                 2d24s22
                   ivp=iv+irefo(jsbv)                                   2d24s22
                   do i=0,nok-1                                         2d24s22
                    iadh=ih0n(jsbv)+ivp+nh0(jsbv)*ibc(ndhvnotv(1)+i)    2d24s22
                    bc(jtmp+i)=bc(iadh)                                 2d24s22
                   end do                                               2d24s22
                   jtmp=jtmp+nok                                        2d24s22
                  end do                                                2d24s22
                 end if                                                 12d7s21
                 call dgemm('n','n',nrow,nvirt(jsbv),nok,sr2,           2d25s21
     $              bc(idhvnotv(1)),nrow,bc(itmp),nok,fact,bc(intden),  2d25s21
     $                nrow,                                             2d25s21
     d'hcsdbk4.  1')
                 fact=1d0                                                 12d21s20
                 ibcoff=itmp                                              12d21s20
                end if                                                    12d21s20
                nok=0                                                     12d21s20
                do isc=1,nsymb
                 iscv=multh(isc,lgoal)                                  12d11s21
                 do isd=1,nsymb                                         8d16s21
                  iscdv=multh(iscv,isd)                                    12d19s20
                  nn=irefo(isd)*irefo(isc)                                 12d18s20
                  nnn=nn*irefo(iscdv)                                       12d18s20
                  if(min(nokdc(isd,isc,1),nok4f(1)).gt.0)then           2d25s21
                   nok=nokdc(isd,isc,1)                                 2d25s21
                   nokf=nok4f(1)                                        2d25s21
                   ncol=nok*nokf                                            12d21s20
                   irdirc=irefo(isd)*irefo(isc)                          8d16s21
                   itmp=ibcoff                                            12d21s20
                   ibcoff=itmp+nok*nvirt(jsbv)                            12d21s20
                   call enough('hcsdbk4. 16',bc,ibc)
                   do iz=itmp,ibcoff-1                                  8d16s21
                    bc(iz)=0d0                                          8d16s21
                   end do                                               8d16s21
                   do i=0,nok-1                                         8d16s21
                    idt=ibc(nd1vnotv(1,isd,isc)+i)/nnn                  12d7s21
                    left=ibc(nd1vnotv(1,isd,isc)+i)-idt*nnn             12d7s21
                    idv=left/irdirc                                     12d7s21
                    idc=left-irdirc*idv                                 12d7s21
                    ic=idc/irefo(isd)                                   12d7s21
                    id=idc-irefo(isd)*ic                                12d7s21
                    idd=ionex(isc,isd,iscdv)+ic+irefo(isc)*(id+         12d8s21
     $                    irefo(isd)*(idv+irefo(iscdv)*nvirt(jsbv)*idt))12d13s21
                    jtmp=itmp+i                                         8d16s21
                    do jv=0,nvirt(jsbv)-1                               8d16s21
                     bc(jtmp+jv*nok)=bc(idd+jv*nnn)                     12d9s21
                    end do                                              8d16s21
                   end do                                               8d16s21
                   call dgemm('n','n',nrow,nvirt(jsbv),nok,sr2,         2d25s21
     $           bc(id1vnotv(1,isd,isc)),nrow,bc(itmp),nok,fact,        2d25s21
     $                bc(intden),nrow,                                  12d21s20
     d'hcsdbk4.  2')
                   fact=1d0                                                 12d21s20
                   ibcoff=itmp                                              12d21s20
                  end if
                 end do                                                    12d21s20
                end do                                                  12d9s21
                if(fact.gt.0.5d0)then                                     12d21s20
                 itrans=ibcoff                                            12d21s20
                 ibcoff=itrans+nvirt(jsbv)*nrow                           12d21s20
                 call enough('hcsdbk4. 17',bc,ibc)
                 do i=0,nvirt(jsbv)-1                                     12d21s20
                  do j=0,nrow-1                                           12d21s20
                   ji=intden+j+nrow*i                                   12d29s20
                   ij=itrans+i+nvirt(jsbv)*j                              12d21s20
                   bc(ij)=bc(ji)                                          12d21s20
                  end do                                                  12d21s20
                 end do                                                   12d21s20
                 do if=0,nokf-1                                           12d21s20
                  do ir=0,nrootu-1                                        12d21s20
                   iadvd=loffdnon+nvirt(jsbv)*(ir                       8d16s21
     $                  +nrootu*ibc(mdhvnotv(1)+if))                    2d25s21
                   do j=0,njhere-1                                        12d21s20
                    iad=itrans+nvirt(jsbv)*(j+njhere*if)                  12d21s20
                    iadg=jgs+nvirt(jsbv)*(ir+nrootu*(i11s+j-1))         1d28s21
                    do iv=0,nvirt(jsbv)-1                                  12d21s20
                     bc(iadg+iv)=bc(iadg+iv)+bc(iad+iv)*vd(iadvd+iv)    1d27s21
                    end do                                                12d21s20
                   end do                                                 12d21s20
                  end do                                                  12d21s20
                 end do                                                   12d21s20
                end if                                                    12d21s20
                ibcoff=intden                                             12d21s20
               end if                                                     12d21s20
               if(ksbv12.eq.1)then                                      12d22s20
                loffdnon=loffdnon+nfdat(2,1,isb)*nvirt(ksbv1)*nrootu    8d16s21
                nvv=(nvirt(ksbv1)*(nvirt(ksbv1)-1))/2                    12d21s20
                isw=0                                                   12d21s20
               else                                                      12d21s20
                nvv=nvirt(ksbv1)*nvirt(ksbv2)                            12d21s20
                isw=1                                                   12d21s20
               end if                                                    12d21s20
               do l=1,4                                                 12d21s20
                if(nfdat(2,l,isb).gt.0)then                             12d21s20
                 if(jsbv.eq.ksbv1.or.jsbv.eq.ksbv2)then                   12d22s20
                  if(ksbv1.eq.jsbv)then                                    12d21s20
                   isbvu=ksbv2                                             12d21s20
                   tf=1d0                                                  12d21s20
                  else                                                     12d21s20
                   isbvu=ksbv1                                             12d21s20
                   tf=-1d0                                                 12d21s20
                  end if                                                   12d21s20
                  if(l.eq.1)then                                         12d21s20
                   factt=1d0                                             12d21s20
                   tf2=1d0
                  else                                                   12d21s20
                   tf2=-1d0                                             12d29s20
                   factt=tf                                              12d21s20
                  end if                                                 12d21s20
                  if(nok4f(l).gt.0.and.nvirt(isbvu).gt.0)then           3d2s21
                   nrow=njhere*nok4f(l)                                  12d21s20
                   intden=ibcoff                                         12d21s20
                   ibcoff=intden+nrow*nvirt(isbvu)                       12d21s20
                   call enough('hcsdbk4. 18',bc,ibc)
                   fact=0d0                                              12d21s20
                   if(nok4(l).gt.0)then                                 1d1s21
c     Gr jv j = D j k i H0iv2 V jv iv2 r k
                    itmp=ibcoff                                          12d21s20
                    ibcoff=itmp+nok4(l)*nvirt(isbvu)                     12d21s20
                    call enough('hcsdbk4. 19',bc,ibc)
                    jtmp=itmp                                            12d21s20
                    iosym=multh(isbvu,isopt(1))                         12d9s21
                    if(ih0n(iosym).gt.0)then                            2d24s22
                     do iv=0,nvirt(isbvu)-1                             2d24s22
                      ivp=iv+irefo(isbvu)                               2d24s22
                      ih0=ih0n(iosym)+nh0(iosym)*ivp                    2d24s22
                      do i=0,nok4(l)-1                                  2d24s22
                       ii=ibc(ndhvnotv(l)+i)                            2d24s22
                       bc(jtmp+i)=bc(ih0+ii)                            2d24s22
                      end do                                            2d24s22
                      jtmp=jtmp+nok4(l)                                 2d24s22
                     end do                                             2d24s22
                    else if(ih0n(isbvu).gt.0)then                       2d24s22
                     do iv=0,nvirt(isbvu)-1                             2d24s22
                      ivp=iv+irefo(isbvu)                               2d24s22
                      ih0=ih0n(isbvu)+ivp                               2d24s22
                      do i=0,nok4(l)-1                                  2d24s22
                       ii=ibc(ndhvnotv(l)+i)                            2d24s22
                       bc(jtmp+i)=bc(ih0+nh0(isbvu)*ii)                 2d24s22
                      end do                                            2d24s22
                      jtmp=jtmp+nok4(l)                                 2d24s22
                     end do                                             12d7s21
                    end if                                              12d7s21
                    call dgemm('n','n',nrow,nvirt(isbvu),nok4(l),factt,  12d21s20
     $                  bc(idhvnotv(l)),nrow,bc(itmp),nok4(l),fact,     12d23s20
     $                  bc(intden),nrow,                                12d21s20
     d'hcsdbk4.  3')
                    ibcoff=itmp                                          12d21s20
                    fact=1d0                                             12d21s20
                   end if                                                12d21s20
                   do isc=1,nsymb                                        12d22s20
                    iscv=multh(isc,lgoal)                                12d22s20
                    do isd=1,nsymb                                      8d16s21
                     iscdv=multh(iscv,isd)                               12d22s20
                     if(nokdc(isd,isc,l).gt.0)then                       12d22s20
                      irdrc=irefo(isd)*irefo(isc)                       8d16s21
                      nn=irefo(isc)*irefo(isd)                          12d22s20
                      nnn=nn*irefo(iscdv)                                12d22s20
                      itmp=ibcoff                                        12d22s20
                      ibcoff=itmp+nokdc(isd,isc,l)*nvirt(isbvu)          12d22s20
                      call enough('hcsdbk4. 20',bc,ibc)
                      do iz=itmp,ibcoff-1                               8d16s21
                       bc(iz)=0d0                                       8d16s21
                      end do                                            8d16s21
                      do i=0,nokdc(isd,isc,l)-1                         12d22s20
                       idt=ibc(nd1vnotv(l,isd,isc)+i)/nnn               12d7s21
                       left=ibc(nd1vnotv(l,isd,isc)+i)-nnn*idt          12d7s21
                       idv=left/irdrc                                   12d7s21
                       idc=left-irdrc*idv                               12d7s21
                       ic=idc/irefo(isd)                                12d7s21
                       id=idc-irefo(isd)*ic                             12d7s21
                       jtmp=itmp+i                                      3d23s21
                       idvk=ionex(isc,isd,iscdv)+ic+irefo(isc)*(id      12d8s21
     $            +irefo(isd)*(idv+irefo(iscdv)*nvirt(isbvu)*idt))      12d8s21
                       do kv=0,nvirt(isbvu)-1                           8d16s21
                        bc(jtmp+kv*nokdc(isd,isc,l))=bc(idvk+nnn*kv)          8d16s21
                       end do                                           8d16s21
                      end do                                            8d16s21
                      call dgemm('n','n',nrow,nvirt(isbvu),
     $                   nokdc(isd,isc,l),factt,bc(id1vnotv(l,isd,isc)),12d22s20
     $                    nrow,bc(itmp),nokdc(isd,isc,l),fact,          12d22s20
     $                    bc(intden),nrow,                              12d22s20
     d'hcsdbk4.  4')
                      ibcoff=itmp                                        12d22s20
                      fact=1d0                                           12d22s20
                     end if                                              8d16s21
                    end do                                               12d22s20
                   end do                                                12d22s20
                   if(fact.gt.0.5d0)then                                 12d22s20
                    itrans=ibcoff                                        12d22s20
                    ibcoff=itrans+nrow*nvirt(isbvu)                      12d22s20
                    call enough('hcsdbk4. 21',bc,ibc)
                    if(ksbv1.eq.ksbv2)then                              12d29s20
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
                    ncol=nrootu*nvirt(jsbv)                              12d22s20
                    itmpsv=ibcoff                                         12d22s20
                    itmpsg=itmpsv+njhere*ncol                              12d22s20
                    ibcoff=itmpsg+njhere*ncol                             12d22s20
                    call enough('hcsdbk4. 22',bc,ibc)
                    if(ksbv1.eq.ksbv2)then                               12d22s20
                     do i=itmpsg,ibcoff-1                                 12d22s20
                      bc(i)=0d0                                           12d22s20
                     end do                                               12d22s20
                     do k=0,nok4f(l)-1                                   12d22s20
                      kk=ibc(mdhvnotv(l)+k)                              12d22s20
                      do j=0,njhere-1                                    12d22s20
                       do ir=0,nrootu-1                                  12d22s20
                        iadv=loffdnon+nvv*(ir+nrootu*kk)                 12d22s20
                        iadd=itrans+nvirt(isbvu)*(j+njhere*k)            12d22s20
                        iadsg=itmpsg+nvirt(jsbv)*(ir+nrootu*j)            12d22s20
                        do iv2=0,nvirt(ksbv2)-1                          12d22s20
                         do iv1=0,iv2-1                                  12d22s20
                          bc(iadsg+iv2)=bc(iadsg+iv2)                    12d22s20
     $                        +bc(iadd+iv1)*vd(iadv+iv1)*tf2            12d29s20
                          bc(iadsg+iv1)=bc(iadsg+iv1)                   12d29s20
     $                         +bc(iadd+iv2)*vd(iadv+iv1)                12d31s20
                         end do                                          12d22s20
                         iadv=iadv+iv2                                   12d22s20
                        end do                                           12d22s20
                       end do                                            12d22s20
                      end do                                             12d22s20
                     end do                                              12d22s20
                     do ir=0,nrootm                                     1d27s21
                      do j=0,njhere-1                                     12d22s20
                       iadg=itmpsg+nvirt(jsbv)*(ir+nrootu*j)            1d27s21
                       iad=jgs+nvirt(jsbv)*(ir+nrootu*(j+i11s-1))       1d27s21
                       do iv=0,nvirt(jsbv)-1                              12d22s20
                        bc(iad+iv)=bc(iad+iv)+bc(iadg+iv)               1d27s21
                       end do                                            12d22s20
                      end do                                             12d22s20
                     end do                                              12d22s20
                    else if(isbvu.eq.ksbv1)then                          12d22s20
c     r v2 j, j f v1, v1v2 r f
c     f v1 r v2, f v1 j, j r v2
                     mrow=nvirt(ksbv1)*nok4f(l)                          12d22s20
                     mcol=nvirt(ksbv2)*nrootu                            12d22s20
                     itmpdv=ibcoff                                        12d22s20
                     itmpdg=itmpdv+mrow*mcol                             12d22s20
                     ibcoff=itmpdg+mrow*mcol                             12d22s20
                     call enough('hcsdbk4. 23',bc,ibc)
                     do k=0,nok4f(l)-1                                   12d22s20
                      kk=ibc(mdhvnotv(l)+k)                              12d22s20
                      do ir=0,nrootu-1                                   12d22s20
                       do iv2=0,nvirt(ksbv2)-1                           12d22s20
                        do iv1=0,nvirt(ksbv1)-1                          12d22s20
                         iad1=loffdnon+iv1                              8d16s21
     $                  +nvirt(ksbv1)*(iv2+nvirt(ksbv2)*(ir+nrootu*kk)) 12d22s20
                         iad2=itmpdv+k+nok4f(l)*(iv1+nvirt(ksbv1)        12d22s20
     $                       *(iv2+nvirt(ksbv2)*ir))                    1d29s21
                         bc(iad2)=vd(iad1)                               12d22s20
                        end do                                           12d22s20
                       end do                                            12d22s20
                      end do                                             12d22s20
                     end do                                              12d22s20
                     call dgemm('n','n',njhere,mcol,mrow,1d0,           3d2s21
     $                     bc(intden),njhere,bc(itmpdv),mrow,0d0,       3d2s21
     $                     bc(itmpsg),njhere,                           3d2s21
     d'hcsdbk4.  5')
                     do ir=0,nrootm                                     1d27s21
                      do jv=0,nvirt(jsbv)-1                              12d22s20
                       iad1=itmpsg+njhere*(jv+nvirt(jsbv)*ir)           1d27s21
                       do j=0,njhere-1                                   12d22s20
                        iad2=jgs+jv+nvirt(jsbv)*(ir+nrootu*(i11s+j-1))  1d27s21
                        bc(iad2)=bc(iad2)+bc(iad1+j)                       12d22s20
                       end do                                            12d22s20
                      end do                                             12d22s20
                     end do                                              12d22s20
                    else                                                 12d22s20
c     r v1 j, j f v2, v1v2 r f
c     f v1 r v2, f v2 j, j r v1
                     mrow=nrootu*nvirt(ksbv1)                            12d22s20
                     mcol=nok4f(l)*nvirt(ksbv2)                          12d22s20
                     itmpdv=ibcoff                                       12d22s20
                     itmpgv=itmpdv+mrow*mcol                             12d22s20
                     ibcoff=itmpgv+mrow*mcol                             12d22s20
                     call enough('hcsdbk4. 24',bc,ibc)
                     do k=0,nok4f(l)-1                                   12d22s20
                      kk=ibc(mdhvnotv(l)+k)                              12d22s20
                      do ir=0,nrootu-1                                   12d22s20
                       do iv2=0,nvirt(ksbv2)-1                            12d22s20
                        do iv1=0,nvirt(ksbv1)-1                          12d22s20
                         iad1=loffdnon+iv1+nvirt(ksbv1)*(iv2            12d22s20
     $                        +nvirt(ksbv2)*(ir+nrootu*kk))             12d22s20
                         iad2=itmpdv+k+nok4f(l)*(iv2+nvirt(ksbv2)        12d22s20
     $                       *(iv1+nvirt(ksbv1)*ir))                    1d27s21
                         bc(iad2)=vd(iad1)                               12d22s20
                        end do                                           12d22s20
                       end do                                            12d22s20
                      end do                                             12d22s20
                     end do                                              12d22s20
                     call dgemm('n','n',njhere,mrow,mcol,1d0,           3d2s21
     $                     bc(intden),njhere,bc(itmpdv),mcol,0d0,       3d2s21
     $                     bc(itmpsg),njhere,                           3d2s21
     d'hcsdbk4.  6')
                     do i=0,ncol-1                                       12d22s20
                      do j=0,njhere-1                                    12d22s20
                       ji=itmpsg+j+njhere*i                              12d22s20
                       ij=jgs+i+ncol*(j+i11s-1)                          12d22s20
                       bc(ij)=bc(ij)+bc(ji)                              12d22s20
                      end do                                             12d22s20
                     end do                                              12d22s20
                    end if                                               12d22s20
                   end if                                                12d22s20
                   ibcoff=intden                                         12d22s20
                  end if                                                 12d21s20
                 end if                                                 12d22s20
                 loffdnon=loffdnon+nvv*nfdat(2,l,isb)*nrootu            12d22s20
                end if                                                  12d21s20
               end do                                                   12d21s20
              end if                                                     12d21s20
             end do                                                     12d21s20
             jgs=jgs+ncsfb(jarg)*nvirt(jsbv)*nrootu                      12d24s20
            end do                                                      12d18s20
            ibcoff=idhvisv                                              12d24s20
           end if                                                       12d18s20
          end do                                                        12d18s20
          do isbv1=1,nsymb                                              12d21s20
           ksbv1=isbv1                                                  8d16s21
           ksbv2=multh(ksbv1,ksbv12)                                    12d21s20
           if(ksbv1.le.ksbv2)then                                       12d21s20
            if(ksbv1.eq.ksbv2)then                                      12d21s20
             nvvs=(nvirt(ksbv1)*(nvirt(ksbv1)+1))/2                     12d21s20
             nvvt=(nvirt(ksbv1)*(nvirt(ksbv1)-1))/2                     12d21s20
            else                                                        12d21s20
             nvvs=nvirt(ksbv1)*nvirt(ksbv2)                             12d21s20
             nvvt=nvvs                                                  12d21s20
            end if                                                      12d21s20
            koffdnon=koffdnon+(nvvs*nfdat(2,1,isb)                      8d16s21
     $           +nvvt*(nfdat(2,2,isb)+nfdat(2,3,isb)+nfdat(2,4,isb)))  12d21s20
     $             *nrootu                                              12d21s20
           end if                                                       12d21s20
          end do                                                        12d21s20
         end do                                                         12d18s20
c
c     igs is nvirt,nroot, ...
c     transpose it to nroot,nvirt,...
c
         itmp=ibcoff                                                    1d27s21
         ibcoff=itmp+nggb                                               8d16s21
         nggm=nggb-1                                                    8d16s21
         call enough('hcsdbk4. 25',bc,ibc)
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
          jgg=jgg+nggb                                                  8d16s21
         end do                                                         1d27s21
         ibcoff=itmp                                                    1d27s21
         i28=nggb                                                       8d19s21
         call ddi_iacc(bc,ibc,ihsdiagb(nclo1p,jsb),i18,i28,i38,i48,     11d15s22
     $        bc(igg),ibc(iacc),nacc)                                   11d15s22
         nwiacc=nwiacc+i28*(i48+1-i38)                                  1d26s23
         last8(1)=ihsdiagb(nclo1p,jsb)                                  8d21s21
         last8(2)=i48                                                   2d8s21
         itransgg=0                                                     1d30s21
         ibcoff=ibcgg                                                   1d30s21
        end if                                                          12d18s20
       end do                                                           12d18s20
       ibcoff=ibcbmat                                                   2d3s21
      end do                                                            12d18s20
      ibcoff=ircv                                                       1d30s21
      return
      end
