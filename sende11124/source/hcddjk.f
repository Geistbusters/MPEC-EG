c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcddjk(nff22,nfdat,gd,vd,nsymb,mdon,mdoo,nec,multh,    1d8s21
     $     isymmrci,nvirt,ncsf,ncsf2,irel,ism,irefo,ismult,ixw1,ixw2,   1d8s21
     $     norb,nrootu,ih0av,nh0av,ioooo,jmats,kmats,ndoub,mdoub,shift, 2d18s21
     $     tdendd,tovr,sr2,srh,iden1e,idenhvv,idenj,idenk,bc,ibc,       5d17s23
     $     idenput,idenorg)                                             5d17s23
      implicit real*8 (a-h,o-z)                                         1d8s21
c                                                                       1d8s21
c     0 and 2 virt contribution to gd=hdd*vd                            1d8s21
c                                                                       1d8s21
      logical lprt,lchoice                                              3d19s21
      integer*1 nab1(2),nab2(2),nab1b(2),nab2b(2)                       2d22s21
      integer*8 ipack8,itc,ito,jtc,jto,itestc,itesto,itesta,itestb,     1d17s21
     $     jtesta,jtestb,gandcc,gandco,gandcb,igoal,                    5d17s23
     $     idenput(nsymb,nsymb,2)                                       5d17s23
      external second                                                   2d18s21
      equivalence (ipack8,ipack4)                                       1d8s21
      dimension nff22(mdoo+1,2,nsymb),nfdat(5,4,*),vd(*),gd(*),         1d8s21
     $     multh(8,8),nvirt(*),ncsf(*),ncsf2(4,*),irel(*),ism(*),       1d8s21
     $     irefo(*),ih0av(*),nh0av(*),ioooo(*),jmats(*),kmats(*),       1d8s21
     $     ipack4(2),nprei(4),mprei(4),nli(4),nprej(4),                 1d12s21
     $     mprej(4),nlj(4),nab4(2,3),idenj(4,8,8,8),ndenj(4,8),         11d7s22
     $     mdenj(4,8),idenk(4,4,8,8,8),ndenk(4,8),mdenk(4,8),loope(6),  11d7s22
     $     idenhvv(4,8,*),mdenhvv(4),itest(32,3),iden1e(4,*),mden1e(4), 11d7s22
     $     ldenj(4,8),ldenk(4,8),iden1ef(4),idenhvvf(4),ivects(4),      11d3s22
     $     idenjf(4,8),idenkf(4,4,8),ndenjf(4,8),ndenkf(4,4,8),
     $     nokj(4,8),nokk(4,4,8),isorb(32),idorb(32),idenjn(4,8),       2d23s21
     $     idenkn(4,4,8),iden1en(4),idenhvvn(4),tdendd(*),ioxx(2)       11d1s22
      include "common.store"                                            1d8s21
      common/kmfind/invk1(2,8,8,8,2)                                    6d30s18
      common/fnd2cm/inv(2,8,8,8)                                        9d2s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      nnes=0
      ndelta2=0
      nsamex=0
      ntotcalc=0
      nnzznn=0
      mplanb=0
      ibcoffo=ibcoff
      norbx=norb+1                                                      1d8s21
      norbxx=norbx+1                                                    1d8s21
      norbxxx=norbxx+1                                                  1d8s21
      ioffdnon=1                                                        1d8s21
      nrootm=nrootu-1                                                   1d11s21
      ioffd=1                                                           1d11s21
      do isb=1,nsymb                                                    1d8s21
       isbv12=multh(isb,isymmrci)                                       1d8s21
       joffdnon=1                                                       1d8s21
       do jsb=1,isb                                                     1d8s21
        ijsb=multh(jsb,isb)                                             1d9s21
        jsbv12=multh(jsb,isymmrci)                                      1d8s21
        ibc0=ibcoff                                                     1d8s21
        do l=1,4                                                        1d8s21
         if(nfdat(2,l,isb).gt.0)then                                    2d24s21
          do ll=1,4                                                     1d8s21
           if(nfdat(2,ll,jsb).gt.0)then                                 2d24s21
            nll=nfdat(2,l,isb)*nfdat(2,ll,jsb)                          2d24s21
            if(l.eq.ll)then                                             1d11s21
             iden1en(l)=ibcoff                                          2d24s21
             idenhvvn(l)=iden1en(l)+nll                                 2d24s21
             ibcoff=idenhvvn(l)+nll                                     2d24s21
            end if                                                      1d11s21
            do lsa=1,nsymb                                              1d8s21
             lsb=multh(lsa,ijsb)                                        1d9s21
             if(lsb.ge.lsa.and.l.eq.ll)then                             1d11s21
              if(lsb.eq.lsa)then                                        1d9s21
               nn=(irefo(lsa)*(irefo(lsa)+1))/2                         1d9s21
              else                                                      1d9s21
               nn=irefo(lsa)*irefo(lsb)                                 1d9s21
              end if                                                    1d9s21
              ndenjf(l,lsa)=ibcoff                                      2d25s21
              ibcoff=ndenjf(l,lsa)+nn                                    1d11s21
              idenjn(l,lsa)=ibcoff                                      1d11s21
              ibcoff=idenjn(l,lsa)+nn*nll                               2d24s21
             end if                                                     1d9s21
             nn=irefo(lsa)*irefo(lsb)                                   1d9s21
             ndenkf(l,ll,lsa)=ibcoff                                    2d25s21
             ibcoff=ndenkf(l,ll,lsa)+nn                                     1d11s21
             idenkn(l,ll,lsa)=ibcoff                                    1d9s21
             ibcoff=idenkn(l,ll,lsa)+nn*nll                             2d24s21
            end do                                                      1d8s21
           end if                                                       1d8s21
          end do                                                        1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
        call enough('hcddjk.  1',bc,ibc)
        nwds=ibcoff-ibc0                                                1d8s21
        do i=ibc0,ibcoff-1                                              1d8s21
         bc(i)=0d0                                                      1d8s21
        end do                                                          1d8s21
        ibuff=ibcoff                                                    5d17s23
        ibcoff=ibuff+idenput(isb,jsb,2)                                 5d17s23
        call enough('hcddjk.buff',bc,ibc)                               5d17s23
        ibsft=ibuff-idenorg                                             5d17s23
        itc=1                                                           5d17s23
        ito=idenput(isb,jsb,2)                                          5d17s23
        call ddi_get(bc,ibc,idenput(isb,jsb,1),itc,itc,itc,ito,         5d17s23
     $       bc(ibuff))                                                 5d17s23
        do l=1,4                                                        1d8s21
         if(nfdat(2,l,isb).gt.0)then                                    1d8s21
          do ll=1,4                                                     1d8s21
           if(nfdat(2,ll,jsb).gt.0)then                                 1d8s21
            nll=nfdat(2,l,isb)*nfdat(2,ll,jsb)                          1d8s21
            do lsa=1,nsymb                                              1d8s21
             lsb=multh(lsa,ijsb)                                        1d9s21
             if(lsb.ge.lsa.and.l.eq.ll.and.                             2d23s21
     $            min(irefo(lsa),irefo(lsb),nll).gt.0)then              2d23s21
              if(lsa.eq.lsb)then                                        11d7s22
               nn=(irefo(lsa)*(irefo(lsa)+1))/2                         11d7s22
              else                                                      11d7s22
               nn=irefo(lsa)*irefo(lsb)                                 11d7s22
              end if                                                    11d7s22
              itrial=ibcoff                                             11d7s22
              ibcoff=itrial+nll*nn                                      11d7s22
              call enough('hcddjk.x',bc,ibc)
              if(idenj(l,lsa,isb,jsb).lt.0)then                         11d7s22
               idenu=-idenj(l,lsa,isb,jsb)+ibsft                        5d17s23
               nusedi=ibc(idenu)/2                                      5d18s23
               if(2*nusedi.ne.ibc(idenu))nusedi=nusedi+1                5d17s23
               if(nusedi.eq.0)then
                do iz=itrial,ibcoff-1
                 bc(iz)=0d0
                end do
               else
                icmp1=1+idenu                                           5d17s23
                icmp2=icmp1+nn                                           11d7s22
                icmp3=icmp2+nusedi                                       11d7s22
                call uncompxu(ibc(icmp1),ibc(icmp2),bc(icmp3),          11d7s22
     $               bc(itrial),nll,nn)                                 11d7s22
               end if
              else                                                      11d7s22
               idenu=idenj(l,lsa,isb,jsb)+ibsft                         5d17s23
               do i=0,nll*nn-1                                          11d8s22
                bc(itrial+i)=bc(idenu+i)                                5d17s23
               end do                                                   11d8s22
              end if                                                    11d7s22
              idenjf(l,lsa)=itrial                                      11d7s22
             end if                                                     1d9s21
             nn=irefo(lsa)*irefo(lsb)                                   1d9s21
             if(min(nn,nll).gt.0)then                                   2d24s21
              itrial=ibcoff                                             11d7s22
              ibcoff=itrial+nll*nn                                      11d7s22
              call enough('hcddjk.y',bc,ibc)
              if(idenk(l,ll,lsa,isb,jsb).lt.0)then                         11d7s22
               idenu=-idenk(l,ll,lsa,isb,jsb)+ibsft                     5d17s23
               nusedi=ibc(idenu)/2                                      5d17s23
               if(2*nusedi.ne.ibc(idenu))nusedi=nusedi+1                5d17s23
               if(nusedi.eq.0)then
                do iz=itrial,ibcoff-1
                 bc(iz)=0d0
                end do
               else
                icmp1=1+idenu                                           5d17s23
                icmp2=icmp1+nn                                           11d7s22
                icmp3=icmp2+nusedi                                       11d7s22
                call uncompxu(ibc(icmp1),ibc(icmp2),bc(icmp3),          11d7s22
     $               bc(itrial),nll,nn)                                 11d7s22
               end if
              else                                                      11d7s22
               idenu=idenk(l,ll,lsa,isb,jsb)+ibsft                      5d17s23
               do i=0,nll*nn-1                                          11d8s22
                bc(itrial+i)=bc(idenu+i)                                5d17s23
               end do                                                   11d8s22
              end if                                                    11d7s22
              idenkf(l,ll,lsa)=itrial                                   11d7s22
             end if                                                     2d24s21
            end do                                                      1d8s21
           end if                                                       1d8s21
          end do                                                        1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
        nall=0                                                          5d10s21
        do ll=1,4                                                       1d8s21
         if(nfdat(2,ll,jsb).gt.0)then                                   1d8s21
          do l=1,4                                                      1d8s21
           if(nfdat(2,l,isb).gt.0)then                                  1d8s21
            do lsa=1,nsymb                                              1d8s21
             lsb=multh(lsa,ijsb)                                        1d9s21
             nokk(l,ll,lsa)=0                                           1d15s21
             if(min(irefo(lsa),irefo(lsb)).gt.0)then                    1d13s21
              nnn=nfdat(2,l,isb)*nfdat(2,ll,jsb)                        1d9s21
              if(lsb.ge.lsa.and.l.eq.ll)then                            1d13s21
               if(lsa.eq.lsb)then                                        1d8s21
                nn=(irefo(lsa)*(irefo(lsa)+1))/2                         1d8s21
               else                                                      1d8s21
                nn=irefo(lsa)*irefo(lsb)                                 1d8s21
               end if                                                    1d8s21
               jden=idenjf(l,lsa)                                        1d11s21
               jdent=idenjf(l,lsa)                                       1d11s21
               nokj(l,lsa)=0                                             1d12s21
               kden=idenjf(l,lsa)                                        1d12s21
               do icol=0,nn-1                                            1d9s21
                sz=0d0                                                    1d8s21
                do i=0,nnn-1                                             1d9s21
                 sz=sz+bc(jden+i)**2                                    1d13s21
                end do                                                    1d8s21
                sz=sqrt(sz/dfloat(nnn))                                  1d9s21
                if(sz.gt.1d-10)then                                       1d8s21
                 ibc(ndenjf(l,lsa)+nokj(l,lsa))=icol                     1d12s21
                 nokj(l,lsa)=nokj(l,lsa)+1                               1d12s21
                 do i=0,nnn-1                                            1d12s21
                  bc(kden+i)=bc(jden+i)                                  1d12s21
                 end do                                                  1d12s21
                 kden=kden+nnn                                           1d12s21
                 nall=nall+1                                            5d10s21
                end if                                                   1d13s21
                jden=jden+nnn                                            1d9s21
               end do                                                    1d9s21
              end if                                                    1d8s21
              nn=irefo(lsa)*irefo(lsb)                                  1d8s21
              jden=idenkf(l,ll,lsa)                                     1d9s21
              jdent=idenkf(ll,l,lsb)                                    1d10s21
              kden=idenkf(l,ll,lsa)                                     1d12s21
              do icol=0,nn-1                                            1d9s21
               ib=icol/irefo(lsa)
               ia=icol-irefo(lsa)*ib
               icolt=ib+irefo(lsb)*ia                                   1d10s21
               jdent=idenkf(ll,l,lsb)+nnn*icolt                         1d10s21
               sz=0d0                                                    1d8s21
               do i=0,nnn-1                                             1d9s21
                sz=sz+bc(jden+i)**2                                     1d9s21
               end do                                                    1d8s21
               sz=sqrt(sz/dfloat(nnn))                                  1d9s21
               if(sz.gt.1d-10)then                                       1d8s21
                ibc(ndenkf(l,ll,lsa)+nokk(l,ll,lsa))=icol               1d12s21
                nokk(l,ll,lsa)=nokk(l,ll,lsa)+1                         1d12s21
                do i=0,nnn-1                                            1d12s21
                 bc(kden+i)=bc(jden+i)                                  1d12s21
                end do                                                  1d12s21
                kden=kden+nnn                                           1d12s21
                nall=nall+1                                             5d10s21
               end if                                                    1d8s21
               jden=jden+nnn                                            1d9s21
              end do                                                    1d9s21
             end if                                                     1d8s21
            end do                                                      1d8s21
           end if                                                       1d8s21
          end do                                                        1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
        if(nall.gt.0)then                                               5d10s21
         ioff=ioffdnon                                                   1d12s21
         do isbv1=1,nsymb                                                1d12s21
          isbv2=multh(isbv1,isbv12)                                      1d12s21
          if(isbv2.ge.isbv1)then                                         1d12s21
           if(isbv12.eq.1)then                                           1d12s21
            isw=0                                                        1d12s21
            mvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        1d12s21
            ioffp=ioff+nvirt(isbv1)*nrootu*nfdat(2,1,isb)
           else                                                          1d12s21
            isw=1                                                        1d12s21
            mvv=nvirt(isbv1)*nvirt(isbv2)                                1d12s21
            ioffp=ioff
           end if                                                        1d12s21
           mmvv=mvv*nrootu                                               1d12s21
           joff=joffdnon                                                 1d12s21
           do jsbv1=1,nsymb                                              1d12s21
            jsbv2=multh(jsbv1,jsbv12)                                    1d12s21
            if(jsbv2.ge.jsbv1)then                                       1d12s21
             if(jsbv12.eq.1)then                                         1d12s21
              jsw=0                                                      1d12s21
              nvv=(nvirt(jsbv1)*(nvirt(jsbv2)-1))/2                      1d12s21
              joffp=joff+nvirt(jsbv1)*nrootu*nfdat(2,1,jsb)              1d13s21
             else                                                        1d12s21
              jsw=1                                                      1d12s21
              nvv=nvirt(jsbv1)*nvirt(jsbv2)                              1d12s21
              joffp=joff                                                 1d13s21
             end if                                                      1d12s21
             nnvv=nvv*nrootu                                             1d12s21
c     Gvvrj =djiab (ab|vv") Vvv"ri
c     Gvv'rj=djiab (ab|vv") Vv'v"ri
c     Gvv'rj=djiab (ab|vv") Vv"v'ri
c     Gv'vrj=djiab (ab|vv") Vv'v"ri
c     Gv'vrj=djiab (ab|vv") Vv"v'ri
             if(jsbv2.eq.isbv1.or.jsbv2.eq.isbv2.or.jsbv1.eq.isbv1.or.   1d13s21
     $         jsbv1.eq.isbv2)then                                      1d13s21
              do imatch=1,4                                              1d13s21
               kase=0                                                    1d14s21
               if(jsbv2.eq.isbv1.and.imatch.eq.1)then                    1d13s21
                isbvc=jsbv2                                               1d13s21
                isbl=jsbv1                                                1d13s21
                isbr=isbv2                                                1d13s21
                itf=0                                                     1d13s21
                jtf=0                                                     1d13s21
                ibl=1                                                     1d13s21
                ibr=1                                                     1d13s21
                kase=1
               end if                                                    1d13s21
               if(jsbv2.eq.isbv2.and.imatch.eq.2)then                    1d13s21
                isbvc=jsbv2                                               1d13s21
                isbl=jsbv1                                                1d13s21
                isbr=isbv1                                                1d13s21
                itf=1                                                     1d13s21
                jtf=0                                                    1d15s21
                ibl=1                                                     1d13s21
                ibr=2                                                     1d13s21
                kase=2
               end if                                                    1d13s21
               if(jsbv1.eq.isbv1.and.imatch.eq.3)then                    1d14s21
                isbvc=jsbv1                                               1d13s21
                isbl=jsbv2                                                1d13s21
                isbr=isbv2                                                1d13s21
                jtf=1                                                     1d13s21
                itf=0                                                    1d15s21
                ibl=2                                                     1d13s21
                ibr=1                                                     1d13s21
                kase=3
               end if                                                    1d13s21
               if(jsbv1.eq.isbv2.and.imatch.eq.4)then                    1d13s21
                isbvc=jsbv1                                               1d13s21
                isbl=jsbv2                                                1d13s21
                isbr=isbv1                                                1d13s21
                itf=1                                                     1d13s21
                jtf=1                                                     1d13s21
                ibl=2                                                     1d13s21
                ibr=2                                                     1d13s21
                kase=4                                                    1d13s21
               end if                                                     1d13s21
               iioffp=ioffp                                               1d12s21
               tfi=1d0                                                   1d14s21
               do li=1,4                                                 1d14s21
                if(min(kase,nfdat(2,li,isb)).gt.0)then                   1d14s21
                 jjoffp=joffp                                            1d14s21
                 tfj=1d0                                                 1d14s21
                 do lj=1,4                                               1d14s21
                  if(nfdat(2,lj,jsb).gt.0)then                           1d14s21
                   tf=tfi*tfj
                   call ilimts(nvirt(isbl),nvirt(isbr),mynprocg,            1d13s21
     $                mynowprog,il,ih,i1s,i1e,i2s,i2e)                  1d12s21
                   nhere=ih+1-il                                            1d12s21
                   nhere2=i2e+1-i2s                                         1d12s21
                   if(min(nhere,nvirt(isbvc)).gt.0.and.isb.ne.jsb)then   3d19s21
                    if(min(lj,li).eq.1.and.max(lj,li).gt.1)then                   1d11s21
                     phs=-1d0                                                   1d11s21
                    else                                                        1d11s21
                     phs=+1d0                                                   1d11s21
                    end if                                                      1d11s21
                    mn=nfdat(2,li,isb)*nfdat(2,lj,jsb)                   1d14s21
                    intden=ibcoff                                            1d12s21
                    ibcoff=intden+mn*nhere                                   1d12s21
                    call enough('hcddjk. 22',bc,ibc)
                    fact=0d0                                                1d13s21
                    do lsa=1,nsymb                                              1d8s21
                     lsb=multh(lsa,ijsb)                                        1d9s21
                     if(min(irefo(lsa),irefo(lsb)).gt.0)then             1d14s21
                      i2eu=invk1(1,lsb,lsa,isbl,1)                              1d13s21
                      if(nokk(li,lj,lsa).gt.0)then                              1d12s21
                       itmpi=ibcoff                                         1d12s21
                       ibcoff=itmpi+nokk(li,lj,lsa)*nhere                       1d12s21
                       call enough('hcddjk. 23',bc,ibc)
                       kint=kmats(i2eu)                                     1d12s21
                       do j=0,nokk(li,lj,lsa)-1                                1d12s21
                        jj=ibc(ndenkf(li,lj,lsa)+j)                            1d12s21
                        jb=jj/irefo(lsa)                                   1d12s21
                        ja=jj-irefo(lsa)*jb                                1d12s21
                        kk=jb+irefo(lsb)*ja                                1d12s21
                        do i=0,nhere-1                                       1d12s21
                         ij=kint+i+nhere*kk                                1d12s21
                         ji=itmpi+j+nokk(li,lj,lsa)*i                           1d12s21
                         bc(ji)=bc(ij)                                      1d12s21
                        end do                                              1d12s21
                       end do                                              1d12s21
                       call dgemm('n','n',mn,nhere,nokk(li,lj,lsa),1d0,   1d14s21
     $               bc(idenkf(li,lj,lsa)),mn,bc(itmpi),nokk(li,lj,lsa),1d14s21
     $                    fact,bc(intden),mn,                           1d14s21
     d' hcddjk. 14')
                       fact=1d0                                             1d12s21
                       ibcoff=itmpi                                         1d13s21
                      end if                                                1d12s21
                     end if                                                 1d13s21
                    end do                                                  1d13s21
                    if(fact.gt.0.5d0)then                                    1d12s21
                     nrow=nfdat(2,lj,jsb)*nvirt(isbl)                        1d13s21
                     nmul=nfdat(2,li,isb)*nhere2
                     itmp1=ibcoff                                            1d12s21
                     ibcoff=itmp1+nrow*nmul                                 1d12s21
                     call enough('hcddjk. 24',bc,ibc)
                     do i=itmp1,ibcoff-1                                    1d12s21
                      bc(i)=0d0                                             1d12s21
                     end do                                                 1d12s21
                     i10=i1s                                                1d12s21
                     i1n=nvirt(isbl)                                        1d13s21
                     jntden=intden                                          1d12s21
                     do i2=i2s,i2e                                          1d12s21
                      i2n=i2-i2s                                            1d12s21
                      if(i2.eq.i2e)i1n=i1e                                  1d12s21
                      do i1=i10,i1n                                         1d12s21
                       i1m=i1-1                                             1d12s21
                       do j=0,nfdat(2,lj,jsb)-1                              1d12s21
                        do i=0,nfdat(2,li,isb)-1                             1d12s21
                         iad=itmp1+j+nfdat(2,lj,jsb)*(i1m+nvirt(isbl)*(  1d14s21
     $                       i2n+nhere2*i))                             1d14s21
                         bc(iad)=bc(jntden+i)                               1d13s21
                        end do                                              1d12s21
                        jntden=jntden+nfdat(2,li,isb)                        1d12s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                      i10=1                                                 1d12s21
                     end do                                                 1d12s21
                     ncol=nvirt(isbvc)*nrootu                               1d12s21
                     itmp2=ibcoff                                           1d12s21
                     ibcoff=itmp2+nmul*ncol                                 1d12s21
                     call enough('hcddjk. 25',bc,ibc)
                     nx=nhere2*nfdat(2,li,isb)*nrootu                        1d12s21
                     do i=itmp2,ibcoff-1                                    1d12s21
                      bc(i)=0d0                                             1d12s21
                     end do                                                 1d12s21
                     if(ibr.eq.1)then                                       1d13s21
                      do i2=i2s,i2e                                          1d12s21
                       i2n=i2-i2s
                       iv2=i2-1                                              1d12s21
                       ntop=iv2-1                                            1d12s21
                       ntop=ntop+isw*(nvirt(isbvc)-1-ntop)                   1d12s21
                       do i=0,nfdat(2,li,isb)-1                               1d12s21
                        do ir=0,nrootm                                       1d12s21
                         iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                         jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,li,isb)*ir)        1d13s21
                         do iv1=0,ntop                                       1d12s21
                          irec=iv1+nvirt(isbv1)*iv2                          1d13s21
                          itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                          itri=itri+isw*(irec-itri)                          1d13s21
                          bc(jtmp2+iv1*nx)=vd(iuse+itri)                     1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     else                                                   1d13s21
                      do i2=i2s,i2e                                          1d12s21
                       i2n=i2-i2s
                       iv1=i2-1                                              1d12s21
                       nbot=iv1+1                                            1d12s21
                       nbot=nbot+isw*(0-nbot)                               1d13s21
                       do i=0,nfdat(2,li,isb)-1                               1d12s21
                        do ir=0,nrootm                                       1d12s21
                         iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                         jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,li,isb)*ir)        1d13s21
                         do iv2=nbot,nvirt(isbv2)-1                         1d13s21
                          irec=iv1+nvirt(isbv1)*iv2                          1d13s21
                          itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                          itri=itri+isw*(irec-itri)                          1d13s21
                          bc(jtmp2+iv2*nx)=vd(iuse+itri)                     1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     end if                                                 1d13s21
                     if(li.eq.1.and.isbv12.eq.1)then                     1d15s21
                      iioff=iioffp-nvirt(isbv1)*nrootu*nfdat(2,1,isb)    1d15s21
                      do i2=i2s,i2e                                      1d15s21
                       i2n=i2-i2s                                        1d15s21
                       iv1=i2-1                                          1d15s21
                       do i=0,nfdat(2,1,isb)-1                           1d15s21
                        do ir=0,nrootm                                   1d15s21
                         iuse=iioff+iv1+nvirt(isbv1)*(ir+nrootu*i)       1d15s21
                         jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,1,isb)*(ir    1d15s21
     $                       +nrootu*iv1))                              1d15s21
                         bc(jtmp2)=vd(iuse)*srh                          1d15s21
                        end do                                           1d15s21
                       end do                                            1d15s21
                      end do                                             1d15s21
                     end if                                              1d15s21
                     iprod=ibcoff                                           1d12s21
                     ibcoff=iprod+nrow*ncol                                 1d12s21
                     call enough('hcddjk. 26',bc,ibc)
                     tff=tf*phs                                          1d14s21
                     call dgemm('n','n',nrow,ncol,nmul,tff,              1d14s21
     $                 bc(itmp1),nrow,bc(itmp2),nmul,0d0,                1d12s21
     $                bc(iprod),nrow,                                   1d12s21
     d' hcddjk. 15')
                     if(ibl.eq.1)then                                       1d13s21
                      do iv2=0,nvirt(isbvc)-1                                1d12s21
                       ntop=iv2-1                                            1d12s21
                       ntop=ntop+jsw*(nvirt(jsbv1)-1-ntop)                   1d13s21
                       do ir=0,nrootm                                        1d12s21
                        do iv1=0,ntop                                        1d12s21
                         irec=iv1+nvirt(jsbv1)*iv2                           1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                         itri=itri+jsw*(irec-itri)                          1d13s21
                         iad=jjoffp+itri+nvv*ir                            1d13s21
                         jprod=iprod+nfdat(2,lj,jsb)*(iv1+nvirt(jsbv1)*( 1d14s21
     $                       ir+nrootu*iv2))                            1d14s21
                         do j=0,nfdat(2,lj,jsb)-1                             1d13s21
                          gd(iad+j*nnvv)=gd(iad+j*nnvv)+bc(jprod+j)          1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     else                                                   1d13s21
                      do iv1=0,nvirt(isbvc)-1                                1d12s21
                       nbot=iv1+1                                           1d13s21
                       nbot=nbot+jsw*(0-nbot)                               1d13s21
                       do ir=0,nrootm                                        1d12s21
                        do iv2=nbot,nvirt(jsbv2)-1                          1d13s21
                         irec=iv1+nvirt(jsbv1)*iv2                           1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                         itri=itri+jsw*(irec-itri)                          1d13s21
                         iad=jjoffp+itri+nvv*ir                            1d13s21
                         jprod=iprod+nfdat(2,lj,jsb)*(iv2+nvirt(jsbv2)*( 1d14s21
     $                       ir+nrootu*iv1))                            1d14s21
                         do j=0,nfdat(2,lj,jsb)-1                             1d13s21
                          gd(iad+j*nnvv)=gd(iad+j*nnvv)+bc(jprod+j)          1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     end if                                                 1d13s21
                     if(lj.eq.1.and.jsbv12.eq.1)then                     1d14s21
                      nv=nvirt(jsbv1)*nrootu                             1d14s21
                      do iv2=0,nvirt(isbvc)-1                            1d14s21
                       do ir=0,nrootm                                    1d14s21
                        iad=joff+iv2+nvirt(jsbv1)*ir                     1d14s21
                        jprod=iprod+nfdat(2,lj,jsb)*(iv2+nvirt(jsbv1)    1d14s21
     $                      *(ir+nrootu*iv2))                           1d14s21
                        do j=0,nfdat(2,lj,jsb)-1                         1d14s21
                         gd(iad+j*nv)=gd(iad+j*nv)+bc(jprod+j)*srh       1d14s21
                        end do                                           1d14s21
                       end do                                            1d14s21
                      end do                                             1d14s21
                     end if                                              1d14s21
                    end if                                                   1d12s21
                    ibcoff=intden                                            1d12s21
                   end if                                                1d14s21
                   call ilimts(nvirt(isbr),nvirt(isbl),mynprocg,            1d13s21
     $                   mynowprog,il,ih,i1s,i1e,i2s,i2e)                  1d12s21
                   nhere=ih+1-il                                            1d12s21
                   nhere2=i2e+1-i2s                                         1d12s21
                   if(min(nhere,nvirt(isbvc)).gt.0)then                 3d19s21
                    if(min(lj,li).eq.1.and.max(lj,li).gt.1)then                   1d11s21
                     phs=-1d0                                                   1d11s21
                    else                                                        1d11s21
                     phs=+1d0                                                   1d11s21
                    end if                                                      1d11s21
                    mn=nfdat(2,li,isb)*nfdat(2,lj,jsb)                  1d15s21
                    intden=ibcoff                                            1d12s21
                    ibcoff=intden+mn*nhere                                   1d12s21
                    if(ibcoff.lt.0)write(6,*)('enough6 '),intden,ibcoff,
     $                   mn,nhere
                    call enough('hcddjk. 27',bc,ibc)
                    fact=0d0                                                1d13s21
                    do lsa=1,nsymb                                              1d8s21
                     lsb=multh(lsa,ijsb)                                        1d9s21
                     if(min(irefo(lsa),irefo(lsb)).gt.0)then            1d14s21
                      i2eu=invk1(1,lsa,lsb,isbr,1)                      1d14s21
                      if(nokk(li,lj,lsa).gt.0)then                      1d14s21
                       itmpi=ibcoff                                         1d12s21
                       ibcoff=itmpi+nokk(li,lj,lsa)*nhere                       1d12s21
                       call enough('hcddjk. 28',bc,ibc)
                       kint=kmats(i2eu)                                 1d15s21
                       do j=0,nokk(li,lj,lsa)-1                                 1d12s21
                        jj=ibc(ndenkf(li,lj,lsa)+j)                              1d12s21
                        do i=0,nhere-1                                       1d12s21
                         ij=kint+i+nhere*jj                                 1d12s21
                         ji=itmpi+j+nokk(li,lj,lsa)*i                   1d14s21
                         bc(ji)=bc(ij)                                      1d12s21
                        end do                                              1d12s21
                       end do                                               1d12s21
                       call dgemm('n','n',mn,nhere,nokk(li,lj,lsa),1d0, 1d14s21
     $                  bc(idenkf(li,lj,lsa)),mn,bc(itmpi),             1d14s21
     $                       nokk(li,lj,lsa),fact,bc(intden),mn,        1d14s21
     d' hcddjk. 16')
                       fact=1d0                                             1d12s21
                       ibcoff=itmpi                                         1d13s21
                      end if                                                1d12s21
                     end if                                                 1d13s21
                    end do                                                  1d13s21
                    if(fact.gt.0.5d0)then                                    1d12s21
c        intden
c       isbr  isbl   ibr ibl
c     1 isbv2 jsbv1   1   1
c     2 isbv1 jsbv1   2   1
c     3 isbv2 jsbv2   1   2
c     4 isbv1 jsbv2   2   2
                     nrow=nfdat(2,li,isb)*nvirt(isbr)                       1d13s21
                     nmul=nfdat(2,lj,jsb)*nhere2                            1d13s21
                     itmp1=ibcoff                                            1d12s21
                     ibcoff=itmp1+nrow*nmul                                 1d12s21
                     call enough('hcddjk. 29',bc,ibc)
                     do i=itmp1,ibcoff-1                                    1d12s21
                      bc(i)=0d0                                             1d12s21
                     end do                                                 1d12s21
                     i10=i1s                                                1d12s21
                     i1n=nvirt(isbr)                                        1d13s21
                     jntden=intden                                          1d12s21
                     do i2=i2s,i2e                                          1d12s21
                      i2n=i2-i2s                                            1d12s21
                      if(i2.eq.i2e)i1n=i1e                                  1d12s21
                      do i1=i10,i1n                                         1d12s21
                       i1m=i1-1                                             1d12s21
                       do j=0,nfdat(2,lj,jsb)-1                             1d12s21
                        iad=itmp1+nfdat(2,li,isb)*(i1m+nvirt(isbr)*(i2n     1d13s21
     $                  +nhere2*j))                                     1d12s21
                        do i=0,nfdat(2,li,isb)-1                              1d12s21
                         bc(iad+i)=bc(jntden+i)                            1d13s21
                        end do                                              1d12s21
                        jntden=jntden+nfdat(2,li,isb)                        1d12s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                      i10=1                                                 1d12s21
                     end do                                                 1d12s21
                     ncol=nvirt(isbvc)*nrootu                               1d12s21
                     itmp2=ibcoff                                           1d12s21
                     ibcoff=itmp2+nmul*ncol                                 1d12s21
                     call enough('hcddjk. 30',bc,ibc)
                     nx=nhere2*nfdat(2,lj,jsb)*nrootu                        1d12s21
                     do i=itmp2,ibcoff-1                                    1d12s21
                      bc(i)=0d0                                             1d12s21
                     end do                                                 1d12s21
                     if(ibl.eq.2)then                                      1d13s21
                      do i2=i2s,i2e                                          1d12s21
                       i2n=i2-i2s
                       iv2=i2-1                                              1d12s21
                       ntop=iv2-1                                            1d12s21
                       ntop=ntop+jsw*(nvirt(isbvc)-1-ntop)                   1d12s21
                       do j=0,nfdat(2,lj,jsb)-1                               1d12s21
                        do ir=0,nrootm                                       1d12s21
                         iuse=jjoffp+nvv*(ir+nrootu*j)                    1d13s21
                         jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,lj,jsb)*ir)        1d13s21
                         do iv1=0,ntop                                       1d12s21
                          irec=iv1+nvirt(jsbv1)*iv2                          1d13s21
                          itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                          itri=itri+jsw*(irec-itri)                          1d13s21
                          bc(jtmp2+iv1*nx)=vd(iuse+itri)                     1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     else                                                   1d13s21
                      do i2=i2s,i2e                                          1d12s21
                       i2n=i2-i2s
                       iv1=i2-1                                              1d12s21
                       nbot=iv1+1                                            1d12s21
                       nbot=nbot+jsw*(0-nbot)                               1d13s21
                       do j=0,nfdat(2,lj,jsb)-1                               1d12s21
                        do ir=0,nrootm                                       1d12s21
                         iuse=jjoffp+nvv*(ir+nrootu*j)                    1d13s21
                         jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,lj,jsb)*ir)        1d13s21
                         do iv2=nbot,nvirt(jsbv2)-1                         1d13s21
                          irec=iv1+nvirt(jsbv1)*iv2                          1d13s21
                          itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                          itri=itri+jsw*(irec-itri)                          1d13s21
                          bc(jtmp2+iv2*nx)=vd(iuse+itri)                   1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     end if                                                 1d13s21
                     if(lj.eq.1.and.jsbv12.eq.1)then                     1d15s21
                      jjoff=jjoffp-nvirt(jsbv1)*nrootu*nfdat(2,1,jsb)    1d15s21
                      do i2=i2s,i2e                                      1d15s21
                       i2n=i2-i2s                                        1d15s21
                       iv1=i2-1                                          1d15s21
                       do j=0,nfdat(2,1,jsb)-1                           1d15s21
                        do ir=0,nrootm                                   1d15s21
                         iuse=jjoff+iv1+nvirt(jsbv1)*(ir+nrootu*j)       1d15s21
                         jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,1,jsb)*(ir    1d15s21
     $                        +nrootu*iv1))                              1d15s21
                         bc(jtmp2)=vd(iuse)*srh                          1d15s21
                        end do                                           1d15s21
                       end do                                            1d15s21
                      end do                                             1d15s21
                     end if                                              1d15s21
                     iprod=ibcoff                                           1d12s21
                     ibcoff=iprod+nrow*ncol                                 1d12s21
                     call enough('hcddjk. 31',bc,ibc)
                     tfu=tf*phs                                         1d15s21
                     call dgemm('n','n',nrow,ncol,nmul,tfu,             1d15s21
     $                 bc(itmp1),nrow,bc(itmp2),nmul,0d0,                1d12s21
     $                bc(iprod),nrow,                                   1d12s21
     d' hcddjk. 17')
c     prod
c       isbr  isbvc  ibr ibl
c     1 isbv2 isbv1   1   1
c     2 isbv1 isbv2   2   1
c     3 isbv2 isbv1   1   2
c     4 isbv1 isbv2   2   2
                     if(ibr.eq.2)then                                       1d13s21
                      do iv2=0,nvirt(isbvc)-1                                1d12s21
                       ntop=iv2-1                                            1d12s21
                       ntop=ntop+isw*(nvirt(isbv1)-1-ntop)                   1d13s21
                       do ir=0,nrootm                                        1d12s21
                        do iv1=0,ntop                                        1d12s21
                         irec=iv1+nvirt(isbv1)*iv2                           1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                         itri=itri+isw*(irec-itri)                         1d13s21
                         iad=iioffp+itri+mvv*ir                               1d13s21
                         jprod=iprod+nfdat(2,li,isb)*(iv1+nvirt(isbv1)* 1d14s21
     $                        (ir+nrootu*iv2))                          1d14s21
                         do i=0,nfdat(2,li,isb)-1                             1d13s21
                          gd(iad+i*mmvv)=gd(iad+i*mmvv)+bc(jprod+i)          1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     else                                                   1d13s21
                      do iv1=0,nvirt(isbvc)-1                                1d12s21
                       nbot=iv1+1                                           1d13s21
                       nbot=nbot+isw*(0-nbot)                               1d13s21
                       do ir=0,nrootm                                        1d12s21
                        do iv2=nbot,nvirt(isbv2)-1                          1d13s21
                         irec=iv1+nvirt(isbv1)*iv2                           1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                         itri=itri+isw*(irec-itri)                         1d13s21
                         iad=iioffp+itri+mvv*ir                            1d13s21
                         jprod=iprod+nfdat(2,li,isb)*(iv2+nvirt(isbv2)* 1d14s21
     $                        (ir+nrootu*iv1))                          1d14s21
                         do i=0,nfdat(2,li,isb)-1                             1d13s21
                          gd(iad+i*mmvv)=gd(iad+i*mmvv)+bc(jprod+i)          1d13s21
                         end do                                              1d12s21
                        end do                                               1d12s21
                       end do                                                1d12s21
                      end do                                                 1d12s21
                     end if                                                 1d13s21
                     if(li.eq.1.and.isbv12.eq.1)then                    1d14s21
                      nv=nvirt(isbv1)*nrootu                            1d14s21
                      do iv2=0,nvirt(isbvc)-1                           1d14s21
                       do ir=0,nrootm                                   1d14s21
                        iad=ioff+iv2+nvirt(isbv1)*ir                    1d14s21
                        jprod=iprod+nfdat(2,li,isb)*(iv2+nvirt(isbv1)   1d14s21
     $                       *(ir+nrootu*iv2))                          1d14s21
                        do i=0,nfdat(2,li,isb)-1                        1d14s21
                         gd(iad+i*nv)=gd(iad+i*nv)+bc(jprod+i)*srh      1d18s21
                        end do                                          1d14s21
                       end do                                           1d14s21
                      end do                                            1d14s21
                     end if                                             1d14s21
                    end if                                                   1d12s21
                    ibcoff=intden                                            1d12s21
                   end if                                                   1d13s21
                   jjoffp=jjoffp+nnvv*nfdat(2,lj,jsb)                     1d14s21
                  end if                                                  1d14s21
                  if(jtf.ne.0)tfj=-1d0                                   1d15s21
                 end do                                                  1d14s21
                end if                                                   1d14s21
                iioffp=iioffp+mmvv*nfdat(2,li,isb)                       1d14s21
                if(kase.ne.0)then                                       1d18s23
                 if(itf.ne.0)tfi=-1d0                                     1d14s21
                end if                                                  1d18s23
               end do                                                    1d14s21
               iioffp=ioffp                                               1d12s21
               jjoffp=joffp                                              1d13s21
               tf=1d0                                                     1d13s21
               do l=1,4                                                   1d12s21
                if(min(kase,nfdat(2,l,isb),nfdat(2,l,jsb)).gt.0)then     1d14s21
                 mn=nfdat(2,l,isb)*nfdat(2,l,jsb)                         1d12s21
                 call ilimts(nvirt(isbl),nvirt(isbr),mynprocg,            1d13s21
     $                mynowprog,il,ih,i1s,i1e,i2s,i2e)                  1d12s21
                 nhere=ih+1-il                                            1d12s21
                 nhere2=i2e+1-i2s                                         1d12s21
                 if(min(nhere,nvirt(isbvc)).gt.0)then                    3d19s21
                  intden=ibcoff                                            1d12s21
                  ibcoff=intden+mn*nhere                                   1d12s21
                  call enough('hcddjk. 32',bc,ibc)
                  fact=0d0                                                1d13s21
                  do lsa=1,nsymb                                              1d8s21
                   lsb=multh(lsa,ijsb)                                        1d9s21
                   if(min(irefo(lsa),irefo(lsb)).gt.0.and.lsb.ge.lsa)
     $                  then                                            1d12s21
                    i2eu=inv(1,lsa,lsb,isbl)                              1d13s21
                    icase=inv(2,lsa,lsb,isbl)                             1d13s21
                    if(nokj(l,lsa).gt.0)then                              1d12s21
                     itmpi=ibcoff                                         1d12s21
                     ibcoff=itmpi+nokj(l,lsa)*nhere                       1d12s21
                     call enough('hcddjk. 33',bc,ibc)
                     jint=jmats(i2eu)                                     1d12s21
                     if(icase.eq.1)then                                   1d12s21
                      do j=0,nokj(l,lsa)-1                                 1d12s21
                       jj=ibc(ndenjf(l,lsa)+j)                              1d12s21
                       do i=0,nhere-1                                       1d12s21
                        ij=jint+i+nhere*jj                                 1d12s21
                        ji=itmpi+j+nokj(l,lsa)*i                           1d12s21
                        bc(ji)=bc(ij)                                      1d12s21
                       end do                                              1d12s21
                      end do                                               1d12s21
                     else                                                 1d12s21
                      do j=0,nokj(l,lsa)-1                                1d12s21
                       jj=ibc(ndenjf(l,lsa)+j)                            1d12s21
                       jb=jj/irefo(lsa)                                   1d12s21
                       ja=jj-irefo(lsa)*jb                                1d12s21
                       kk=jb+irefo(lsb)*ja                                1d12s21
                       do i=0,nhere-1                                       1d12s21
                        ij=jint+i+nhere*kk                                1d12s21
                        ji=itmpi+j+nokj(l,lsa)*i                           1d12s21
                        bc(ji)=bc(ij)                                      1d12s21
                       end do                                              1d12s21
                      end do                                              1d12s21
                     end if                                               1d12s21
                     call dgemm('n','n',mn,nhere,nokj(l,lsa),1d0,         1d12s21
     $                  bc(idenjf(l,lsa)),mn,bc(itmpi),nokj(l,lsa),fact,1d12s21
     $                  bc(intden),mn,                                  1d12s21
     d' hcddjk. 18')
                     fact=1d0                                             1d12s21
                     ibcoff=itmpi                                         1d13s21
                    end if                                                1d12s21
                   end if                                                 1d13s21
                  end do                                                  1d13s21
                  if(fact.gt.0.5d0)then                                    1d12s21
                   nrow=nfdat(2,l,jsb)*nvirt(isbl)                        1d13s21
                   nmul=nfdat(2,l,isb)*nhere2
                   itmp1=ibcoff                                            1d12s21
                   ibcoff=itmp1+nrow*nmul                                 1d12s21
                   call enough('hcddjk. 34',bc,ibc)
                   do i=itmp1,ibcoff-1                                    1d12s21
                    bc(i)=0d0                                             1d12s21
                   end do                                                 1d12s21
c     intden
c       isbl  isbr  ibl ibr
c     1 jsbv1 isbv2  1   1
c     2 jsbv1 isbv1  1   2
c     3 jsbv2 isbv2  2   1
c     4 jsbv2 isbv1  2   2
                   i10=i1s                                                1d12s21
                   i1n=nvirt(isbl)                                        1d13s21
                   jntden=intden                                          1d12s21
                   do i2=i2s,i2e                                          1d12s21
                    i2n=i2-i2s                                            1d12s21
                    if(i2.eq.i2e)i1n=i1e                                  1d12s21
                    do i1=i10,i1n                                         1d12s21
                     i1m=i1-1                                             1d12s21
                     do j=0,nfdat(2,l,jsb)-1                              1d12s21
                      do i=0,nfdat(2,l,isb)-1                             1d12s21
                       iad=itmp1+j+nfdat(2,l,jsb)*(i1m+nvirt(isbl)*(i2n   1d13s21
     $                  +nhere2*i))                                     1d12s21
                       bc(iad)=bc(jntden+i)                               1d13s21
                      end do                                              1d12s21
                      jntden=jntden+nfdat(2,l,isb)                        1d12s21
                     end do                                               1d12s21
                    end do                                                1d12s21
                    i10=1                                                 1d12s21
                   end do                                                 1d12s21
                   ncol=nvirt(isbvc)*nrootu                               1d12s21
                   itmp2=ibcoff                                           1d12s21
                   ibcoff=itmp2+nmul*ncol                                 1d12s21
                   call enough('hcddjk. 35',bc,ibc)
                   nx=nhere2*nfdat(2,l,isb)*nrootu                        1d12s21
                   do i=itmp2,ibcoff-1                                    1d12s21
                    bc(i)=0d0                                             1d12s21
                   end do                                                 1d12s21
                   if(ibr.eq.1)then                                       1d13s21
                    do i2=i2s,i2e                                          1d12s21
                     i2n=i2-i2s
                     iv2=i2-1                                              1d12s21
                     ntop=iv2-1                                            1d12s21
                     ntop=ntop+isw*(nvirt(isbvc)-1-ntop)                   1d12s21
                     do i=0,nfdat(2,l,isb)-1                               1d12s21
                      do ir=0,nrootm                                       1d12s21
                       iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                       jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,l,isb)*ir)        1d13s21
                       do iv1=0,ntop                                       1d12s21
                        irec=iv1+nvirt(isbv1)*iv2                          1d13s21
                        itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                        itri=itri+isw*(irec-itri)                          1d13s21
                        bc(jtmp2+iv1*nx)=vd(iuse+itri)                     1d13s21
                       end do                                              1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                    end do                                                 1d12s21
                   else                                                   1d13s21
                    do i2=i2s,i2e                                          1d12s21
                     i2n=i2-i2s
                     iv1=i2-1                                              1d12s21
                     nbot=iv1+1                                            1d12s21
                     nbot=nbot+isw*(0-nbot)                               1d13s21
                     do i=0,nfdat(2,l,isb)-1                               1d12s21
                      do ir=0,nrootm                                       1d12s21
                       iuse=iioffp+mvv*(ir+nrootu*i)                       1d13s21
                       jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,l,isb)*ir)        1d13s21
                       do iv2=nbot,nvirt(isbv2)-1                         1d13s21
                        irec=iv1+nvirt(isbv1)*iv2                          1d13s21
                        itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                        itri=itri+isw*(irec-itri)                          1d13s21
                        bc(jtmp2+iv2*nx)=vd(iuse+itri)                     1d13s21
                       end do                                              1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                    end do                                                 1d12s21
                   end if                                                 1d13s21
                   if(l.eq.1.and.isbv12.eq.1)then                        1d15s21
                    iioff=iioffp-nvirt(isbv1)*nrootu*nfdat(2,1,isb)      1d15s21
                    do i2=i2s,i2e                                        1d15s21
                     iv1=i2-1                                            1d15s21
                     i2n=i2-i2s                                          1d15s21
                     do i=0,nfdat(2,1,isb)-1                              1d15s21
                      do ir=0,nrootm                                     1d15s21
                       iuse=iioff+iv1+nvirt(isbv1)*(ir+nrootu*i)         1d15s21
                       jtmp2=itmp2+i2n+nhere2*(i+nfdat(2,1,isb)*(ir      1d15s21
     $                     +nrootu*iv1))                                1d15s21
                       bc(jtmp2)=vd(iuse)*srh                            1d15s21
                      end do                                             1d15s21
                     end do                                              1d15s21
                    end do                                               1d15s21
                   end if                                                1d15s21
                   iprod=ibcoff                                           1d12s21
                   ibcoff=iprod+nrow*ncol                                 1d12s21
                   call enough('hcddjk. 36',bc,ibc)
                   call dgemm('n','n',nrow,ncol,nmul,tf,                  1d13s21
     $                 bc(itmp1),nrow,bc(itmp2),nmul,0d0,                1d12s21
     $                bc(iprod),nrow,                                   1d12s21
     d' hcddjk. 19')
c        prod
c       isbl  isbvc ibl ibr
c     1 jsbv1 jsbv2  1   1
c     2 jsbv1 jsbv2  1   2
c     3 jsbv2 jsbv1  2   1
c     4 jsbv2 jsbv1  2   2
                   if(ibl.eq.1)then                                       1d13s21
                    do iv2=0,nvirt(isbvc)-1                                1d12s21
                     ntop=iv2-1                                            1d12s21
                     ntop=ntop+jsw*(nvirt(jsbv1)-1-ntop)                   1d13s21
                     do ir=0,nrootm                                        1d12s21
                      do iv1=0,ntop                                        1d12s21
                       irec=iv1+nvirt(jsbv1)*iv2                           1d13s21
                       itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                       itri=itri+jsw*(irec-itri)                          1d13s21
                       iad=jjoffp+itri+nvv*ir                            1d13s21
                       jprod=iprod+nfdat(2,l,jsb)*(iv1+nvirt(jsbv1)*(ir    1d13s21
     $                     +nrootu*iv2))                                  1d12s21
                       do j=0,nfdat(2,l,jsb)-1                             1d13s21
                        gd(iad+j*nnvv)=gd(iad+j*nnvv)+bc(jprod+j)          1d13s21
                       end do                                              1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                    end do                                                 1d12s21
                   else                                                   1d13s21
                    do iv1=0,nvirt(isbvc)-1                                1d12s21
                     nbot=iv1+1                                           1d13s21
                     nbot=nbot+jsw*(0-nbot)                               1d13s21
                     do ir=0,nrootm                                        1d12s21
                      do iv2=nbot,nvirt(jsbv2)-1                          1d13s21
                       irec=iv1+nvirt(jsbv1)*iv2                           1d13s21
                       itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                       itri=itri+jsw*(irec-itri)                          1d13s21
                       iad=jjoffp+itri+nvv*ir                            1d13s21
                       jprod=iprod+nfdat(2,l,jsb)*(iv2+nvirt(jsbv2)*(ir   1d13s21
     $                     +nrootu*iv1))                                1d13s21
                       do j=0,nfdat(2,l,jsb)-1                             1d13s21
                        gd(iad+j*nnvv)=gd(iad+j*nnvv)+bc(jprod+j)          1d13s21
                       end do                                              1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                    end do                                                 1d12s21
                   end if                                                 1d13s21
                   if(l.eq.1.and.jsbv12.eq.1)then                        1d15s21
                    jjoff=jjoffp-nvirt(isbvc)*nrootu*nfdat(2,1,jsb)      1d15s21
                    nv=nrootu*nvirt(isbvc)                               1d15s21
                    do iv1=0,nvirt(isbvc)-1                              1d15s21
                     do ir=0,nrootm                                      1d15s21
                      iad=jjoff+iv1+nvirt(isbvc)*ir                      1d15s21
                      jprod=iprod+nfdat(2,l,jsb)*(iv1+nvirt(jsbv2)*(ir   1d15s21
     $                    +nrootu*iv1))                                 1d15s21
                      do j=0,nfdat(2,l,jsb)-1                            1d15s21
                       gd(iad+j*nv)=gd(iad+j*nv)+bc(jprod+j)*srh         1d15s21
                      end do                                             1d15s21
                     end do                                              1d15s21
                    end do                                               1d15s21
                   end if                                                1d15s21
                  end if                                                   1d12s21
                  ibcoff=intden                                            1d12s21
                 end if                                                   1d13s21
                 if(isb.ne.jsb)then
                  call ilimts(nvirt(isbr),nvirt(isbl),mynprocg,            1d13s21
     $                mynowprog,il,ih,i1s,i1e,i2s,i2e)                  1d12s21
                  nhere=ih+1-il                                            1d12s21
                  nhere2=i2e+1-i2s                                         1d12s21
                  if(min(nhere,nvirt(isbvc)).gt.0)then                   3d19s21
                   intden=ibcoff                                            1d12s21
                   ibcoff=intden+mn*nhere                                   1d12s21
                   call enough('hcddjk. 37',bc,ibc)
                   fact=0d0                                                1d13s21
                   do lsa=1,nsymb                                              1d8s21
                    lsb=multh(lsa,ijsb)                                        1d9s21
                    if(min(irefo(lsa),irefo(lsb)).gt.0.and.              1d13s21
     $                  lsb.ge.lsa)then                                 1d13s21
                     i2eu=inv(1,lsa,lsb,isbr)                             1d13s21
                     icase=inv(2,lsa,lsb,isbr)                            1d13s21
                     if(nokj(l,lsa).gt.0)then                              1d12s21
                      itmpi=ibcoff                                         1d12s21
                      ibcoff=itmpi+nokj(l,lsa)*nhere                       1d12s21
                      call enough('hcddjk. 38',bc,ibc)
                      jint=jmats(i2eu)                                     1d12s21
                      if(icase.eq.1)then                                   1d12s21
                       do j=0,nokj(l,lsa)-1                                 1d12s21
                        jj=ibc(ndenjf(l,lsa)+j)                              1d12s21
                        do i=0,nhere-1                                       1d12s21
                         ij=jint+i+nhere*jj                                 1d12s21
                         ji=itmpi+j+nokj(l,lsa)*i                           1d12s21
                         bc(ji)=bc(ij)                                      1d12s21
                        end do                                              1d12s21
                       end do                                               1d12s21
                      else                                                 1d12s21
                       do j=0,nokj(l,lsa)-1                                1d12s21
                        jj=ibc(ndenjf(l,lsa)+j)                            1d12s21
                        jb=jj/irefo(lsa)                                   1d12s21
                        ja=jj-irefo(lsa)*jb                                1d12s21
                        kk=jb+irefo(lsb)*ja                                1d12s21
                        do i=0,nhere-1                                       1d12s21
                         ij=jint+i+nhere*kk                                1d12s21
                         ji=itmpi+j+nokj(l,lsa)*i                           1d12s21
                         bc(ji)=bc(ij)                                      1d12s21
                        end do                                              1d12s21
                       end do                                              1d12s21
                      end if                                               1d12s21
                      call dgemm('n','n',mn,nhere,nokj(l,lsa),1d0,         1d12s21
     $                  bc(idenjf(l,lsa)),mn,bc(itmpi),nokj(l,lsa),fact,1d12s21
     $                  bc(intden),mn,                                  1d12s21
     d' hcddjk. 20')
                      fact=1d0                                             1d12s21
                      ibcoff=itmpi                                         1d13s21
                     end if                                                1d12s21
                    end if                                                 1d13s21
                   end do                                                  1d13s21
                   if(fact.gt.0.5d0)then                                    1d12s21
c     intden
c       isbr  isbl   ibr ibl
c     1 isbv2 jsbv1   1   1
c     2 isbv1 jsbv1   2   1
c     3 isbv2 jsbv2   1   2
c     4 isbv1 jsbv2   2   2
                    nrow=nfdat(2,l,isb)*nvirt(isbr)                       1d13s21
                    nmul=nfdat(2,l,jsb)*nhere2                            1d13s21
                    itmp1=ibcoff                                            1d12s21
                    ibcoff=itmp1+nrow*nmul                                 1d12s21
                    call enough('hcddjk. 39',bc,ibc)
                    do i=itmp1,ibcoff-1                                    1d12s21
                     bc(i)=0d0                                             1d12s21
                    end do                                                 1d12s21
                    i10=i1s                                                1d12s21
                    i1n=nvirt(isbr)                                        1d13s21
                    jntden=intden                                          1d12s21
                    do i2=i2s,i2e                                          1d12s21
                     i2n=i2-i2s                                            1d12s21
                     if(i2.eq.i2e)i1n=i1e                                  1d12s21
                     do i1=i10,i1n                                         1d12s21
                      i1m=i1-1                                             1d12s21
                      do j=0,nfdat(2,l,jsb)-1                             1d12s21
                       iad=itmp1+nfdat(2,l,isb)*(i1m+nvirt(isbr)*(i2n     1d13s21
     $                  +nhere2*j))                                     1d12s21
                       do i=0,nfdat(2,l,isb)-1                              1d12s21
                        bc(iad+i)=bc(jntden+i)                            1d13s21
                       end do                                              1d12s21
                       jntden=jntden+nfdat(2,l,isb)                        1d12s21
                      end do                                               1d12s21
                     end do                                                1d12s21
                     i10=1                                                 1d12s21
                    end do                                                 1d12s21
                    ncol=nvirt(isbvc)*nrootu                               1d12s21
                    itmp2=ibcoff                                           1d12s21
                    ibcoff=itmp2+nmul*ncol                                 1d12s21
                    call enough('hcddjk. 40',bc,ibc)
                    nx=nhere2*nfdat(2,l,jsb)*nrootu                        1d12s21
                    do i=itmp2,ibcoff-1                                    1d12s21
                     bc(i)=0d0                                             1d12s21
                    end do                                                 1d12s21
                    if(ibl.eq.2)then                                      1d13s21
                     do i2=i2s,i2e                                          1d12s21
                      i2n=i2-i2s
                      iv2=i2-1                                              1d12s21
                      ntop=iv2-1                                            1d12s21
                      ntop=ntop+jsw*(nvirt(isbvc)-1-ntop)                   1d12s21
                      do j=0,nfdat(2,l,jsb)-1                               1d12s21
                       do ir=0,nrootm                                       1d12s21
                        iuse=jjoffp+nvv*(ir+nrootu*j)                    1d13s21
                        jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,l,jsb)*ir)        1d13s21
                        do iv1=0,ntop                                       1d12s21
                         irec=iv1+nvirt(jsbv1)*iv2                          1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                         itri=itri+jsw*(irec-itri)                          1d13s21
                         bc(jtmp2+iv1*nx)=vd(iuse+itri)                     1d13s21
                        end do                                              1d12s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                     end do                                                 1d12s21
                    else                                                   1d13s21
                     do i2=i2s,i2e                                          1d12s21
                      i2n=i2-i2s
                      iv1=i2-1                                              1d12s21
                      nbot=iv1+1                                            1d12s21
                      nbot=nbot+jsw*(0-nbot)                               1d13s21
                      do j=0,nfdat(2,l,jsb)-1                               1d12s21
                       do ir=0,nrootm                                       1d12s21
                        iuse=jjoffp+nvv*(ir+nrootu*j)                    1d13s21
                        jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,l,jsb)*ir)        1d13s21
                        do iv2=nbot,nvirt(jsbv2)-1                         1d13s21
                         irec=iv1+nvirt(jsbv1)*iv2                          1d13s21
                         itri=((iv2*(iv2-1))/2)+iv1                         1d12s21
                         itri=itri+jsw*(irec-itri)                          1d13s21
                         bc(jtmp2+iv2*nx)=vd(iuse+itri)                   1d13s21
                        end do                                              1d12s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                     end do                                                 1d12s21
                    end if                                                 1d13s21
                    if(l.eq.1.and.jsbv12.eq.1)then                       1d15s21
                     jjoff=jjoffp-nvirt(jsbv1)*nrootu*nfdat(2,1,jsb)     1d15s21
                     do i2=i2s,i2e                                       1d15s21
                      iv1=i2-1                                           1d15s21
                      i2n=i2-i2s                                         1d15s21
                      do j=0,nfdat(2,1,jsb)-1                             1d15s21
                       do ir=0,nrootm                                    1d15s21
                        iuse=jjoff+iv1+nvirt(jsbv1)*(ir+nrootu*j)        1d15s21
                        jtmp2=itmp2+i2n+nhere2*(j+nfdat(2,1,jsb)*(ir     1d15s21
     $                      +nrootu*iv1))                               1d15s21
                        bc(jtmp2)=vd(iuse)*srh                           1d15s21
                       end do                                            1d15s21
                      end do                                             1d15s21
                     end do                                              1d15s21
                    end if                                               1d15s21
                    iprod=ibcoff                                           1d12s21
                    ibcoff=iprod+nrow*ncol                                 1d12s21
                    call enough('hcddjk. 41',bc,ibc)
                    call dgemm('n','n',nrow,ncol,nmul,tf,                 1d13s21
     $                 bc(itmp1),nrow,bc(itmp2),nmul,0d0,                1d12s21
     $                bc(iprod),nrow,                                   1d12s21
     d' hcddjk. 21')
c     prod
c       isbr  isbvc  ibr ibl
c     1 isbv2 isbv1   1   1
c     2 isbv1 isbv2   2   1
c     3 isbv2 isbv1   1   2
c     4 isbv1 isbv2   2   2
                    if(ibr.eq.2)then                                       1d13s21
                     do iv2=0,nvirt(isbvc)-1                                1d12s21
                      ntop=iv2-1                                            1d12s21
                      ntop=ntop+isw*(nvirt(isbv1)-1-ntop)                   1d13s21
                      do ir=0,nrootm                                        1d12s21
                       do iv1=0,ntop                                        1d12s21
                        irec=iv1+nvirt(isbv1)*iv2                           1d13s21
                        itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                        itri=itri+isw*(irec-itri)                         1d13s21
                        iad=iioffp+itri+mvv*ir                               1d13s21
                        jprod=iprod+nfdat(2,l,isb)*(iv1+nvirt(isbv1)*(ir    1d13s21
     $                     +nrootu*iv2))                                  1d12s21
                        do i=0,nfdat(2,l,isb)-1                             1d13s21
                         gd(iad+i*mmvv)=gd(iad+i*mmvv)+bc(jprod+i)          1d13s21
                        end do                                              1d12s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                     end do                                                 1d12s21
                    else                                                   1d13s21
                     do iv1=0,nvirt(isbvc)-1                                1d12s21
                      nbot=iv1+1                                           1d13s21
                      nbot=nbot+isw*(0-nbot)                               1d13s21
                      do ir=0,nrootm                                        1d12s21
                       do iv2=nbot,nvirt(isbv2)-1                          1d13s21
                        irec=iv1+nvirt(isbv1)*iv2                           1d13s21
                        itri=((iv2*(iv2-1))/2)+iv1                          1d12s21
                        itri=itri+isw*(irec-itri)                         1d13s21
                        iad=iioffp+itri+mvv*ir                            1d13s21
                        jprod=iprod+nfdat(2,l,isb)*(iv2+nvirt(isbv2)*(ir    1d13s21
     $                      +nrootu*iv1))                                  1d12s21
                        do i=0,nfdat(2,l,isb)-1                             1d13s21
                         gd(iad+i*mmvv)=gd(iad+i*mmvv)+bc(jprod+i)          1d13s21
                        end do                                              1d12s21
                       end do                                               1d12s21
                      end do                                                1d12s21
                     end do                                                 1d12s21
                    end if                                                 1d13s21
                    if(l.eq.1.and.isbv12.eq.1)then                        1d15s21
                     iioff=iioffp-nvirt(isbvc)*nrootu*nfdat(2,1,isb)      1d15s21
                     nv=nrootu*nvirt(isbvc)                               1d15s21
                     do iv1=0,nvirt(isbvc)-1                              1d15s21
                      do ir=0,nrootm                                      1d15s21
                       iad=iioff+iv1+nvirt(isbvc)*ir                      1d15s21
                       jprod=iprod+nfdat(2,l,isb)*(iv1+nvirt(isbv2)*(ir   1d15s21
     $                    +nrootu*iv1))                                 1d15s21
                       do i=0,nfdat(2,l,isb)-1                            1d15s21
                        gd(iad+i*nv)=gd(iad+i*nv)+bc(jprod+i)*srh        1d15s21
                       end do                                             1d15s21
                      end do                                              1d15s21
                     end do                                               1d15s21
                    end if                                                1d15s21
                   end if                                                   1d12s21
                   ibcoff=intden                                            1d12s21
                  end if                                                   1d13s21
                 end if                                                    1d13s21
                end if                                                      1d13s21
                if(kase.ne.0)then                                       1d18s23
                 if(mod(jtf+itf,2).ne.0)then                              1d15s21
                  tf=-1d0                                                 1d14s21
                 end if                                                   1d14s21
                end if                                                  1d18s23
                iioffp=iioffp+mmvv*nfdat(2,l,isb)                          1d13s21
                jjoffp=jjoffp+nnvv*nfdat(2,l,jsb)                            1d13s21
               end do                                                      1d13s21
              end do                                                     1d13s21
             end if                                                      1d13s21
             if(jsbv12.eq.1)joff=joff                                    1d12s21
     $           +nvirt(jsbv1)*nrootu*nfdat(2,1,jsb)                    1d12s21
             do l=1,4                                                    1d12s21
              joff=joff+nnvv*nfdat(2,l,jsb)                              1d12s21
             end do                                                      1d12s21
            end if                                                       1d12s21
           end do                                                        1d12s21
           if(isbv12.eq.1)ioff=ioff+nvirt(isbv1)*nrootu*nfdat(2,1,isb)   1d12s21
           do l=1,4                                                      1d12s21
            ioff=ioff+mmvv*nfdat(2,l,isb)                                1d12s21
           end do                                                        1d12s21
          end if                                                         1d12s21
         end do                                                          1d12s21
         if(ijsb.eq.1)then                                               1d11s21
          ioff=ioffdnon                                                  1d11s21
          do isbv1=1,nsymb                                               1d11s21
           isbv2=multh(isbv1,isbv12)                                     1d11s21
           if(isbv2.ge.isbv1)then                                        1d11s21
            nvv=0                                                        4d28s21
            if(isbv1.eq.isbv2)then                                       1d11s21
             ioffs=ioff                                                  1d12s21
             nnn=nvirt(isbv1)*nrootu                                     1d11s21
             call ilimts(nvirt(isbv1),nrootu,mynprocg,mynowprog,il,ih,   1d11s21
     $           i1s,i1e,i2s,i2e)                                       1d11s21
             nhere=ih+1-il                                               1d11s21
             if(nhere.gt.0)then                                          1d11s21
              itmp=ibcoff                                                1d11s21
              ibcoff=itmp+nhere*nfdat(2,1,isb)                           1d11s21
              call enough('hcddjk. 42',bc,ibc)
              iden1u=iden1e(1,isb)+ibsft                                5d18s23
              call dgemm('n','n',nhere,nfdat(2,1,isb),nfdat(2,1,isb),    1d11s21
     $           1d0,vd(ioff+il-1),nnn,bc(iden1u),nfdat(2,1,isb),       5d18s23
     $            0d0,bc(itmp),nhere,                                   1d11s21
     d' hcddjk. 22')
              do i=0,nfdat(2,1,isb)-1                                    1d11s21
               i10=i1s                                                   1d11s21
               i1n=nvirt(isbv1)                                          1d11s21
               jtmp=itmp+nhere*i                                         1d11s21
               do i2=i2s,i2e                                             1d11s21
                if(i2.eq.i2e)i1n=i1e                                     1d11s21
                iad=ioff-1+nvirt(isbv1)*(i2-1+nrootu*i)                  1d11s21
                do i1=i10,i1n                                            1d11s21
                 gd(iad+i1)=gd(iad+i1)+bc(jtmp)                          1d11s21
                 jtmp=jtmp+1                                             1d11s21
                end do                                                   1d11s21
                i10=1                                                    1d11s21
               end do                                                    1d11s21
              end do                                                     1d11s21
              ibcoff=itmp                                                1d11s21
             end if                                                      1d11s21
             if(nnn.gt.0)then                                            1d25s21
              ioff=ioff+nnn*nfdat(2,1,isb)                                1d11s21
              nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       1d11s21
              isw=0                                                       1d11s21
              itmp=ibcoff                                                 1d18s21
              ibcoff=itmp+nnn*nfdat(2,1,isb)                              1d18s21
              call enough('hcddjk. 43',bc,ibc)
              idenhuu=idenhvv(1,isb,jsb)+ibsft                          5d17s23
              call dgemm('n','n',nnn,nfdat(2,1,isb),nfdat(2,1,isb),       1d18s21
     $          2d0,vd(ioffs),nnn,bc(idenhuu),nfdat(2,1,isb),           5d18s23
     $              0d0,bc(itmp),nnn,                                   1d18s21
     d' hcddjk. 23')
              do i=0,nfdat(2,1,isb)-1                                     1d18s21
               do ir=0,nrootm                                             1d18s21
                do iv=mynowprog,nvirt(isbv1)-1,mynprocg                   1d18s21
                 iad=ioffs+iv+nvirt(isbv1)*(ir+nrootu*i)                  1d18s21
                 ih0=ih0av(isbv1)+(iv+irefo(isbv1))*(nh0av(isbv1)+1)      1d18s21
                 jtmp=itmp+iv+nvirt(isbv1)*(ir+nrootu*i)                  1d18s21
                 gd(iad)=gd(iad)+bc(ih0)*bc(jtmp)                         1d18s21
                end do                                                    1d18s21
               end do                                                     1d18s21
              end do                                                      1d18s21
              ibcoff=itmp                                                 1d18s21
             end if                                                      1d25s21
            else                                                         1d11s21
             nvv=nvirt(isbv1)*nvirt(isbv2)                               1d11s21
             isw=1                                                       1d11s21
            end if                                                       1d11s21
            nnn=nvv*nrootu                                               1d11s21
            tf=1d0                                                       1d11s21
            do l=1,4                                                     1d11s21
             if(nfdat(2,l,isb).gt.0)then                                 1d11s21
              call ilimts(nvv,nrootu,mynprocg,mynowprog,il,ih,            1d11s21
     $            i1s,i1e,i2s,i2e)                                       1d11s21
              nhere=ih+1-il                                               1d11s21
              if(nhere.gt.0)then                                          1d11s21
               itmp=ibcoff                                                1d11s21
               ibcoff=itmp+nhere*nfdat(2,l,isb)                           1d11s21
               if(isbv12.eq.1.and.l.eq.1)then                            1d12s21
                nrow=nvirt(isbv1)*nrootu                                 1d12s21
                ibcoff=max(ibcoff,itmp+nrow*nfdat(2,l,isb))              1d12s21
               end if                                                    1d12s21
               call enough('hcddjk. 44',bc,ibc)
               iden1u=iden1e(l,isb)+ibsft                               5d18s23
               call dgemm('n','n',nhere,nfdat(2,l,isb),nfdat(2,l,isb),    1d11s21
     $           1d0,vd(ioff+il-1),nnn,bc(iden1u),nfdat(2,l,isb),       5d18s23
     $            0d0,bc(itmp),nhere,                                   1d11s21
     d' hcddjk. 24')
               do i=0,nfdat(2,l,isb)-1                                    1d11s21
                i10=i1s                                                   1d11s21
                i1n=nvv                                                  1d11s21
                jtmp=itmp+nhere*i                                         1d11s21
                do i2=i2s,i2e                                             1d11s21
                 if(i2.eq.i2e)i1n=i1e                                     1d11s21
                 iad=ioff-1+nvv*(i2-1+nrootu*i)                          1d11s21
                 do i1=i10,i1n                                            1d11s21
                  gd(iad+i1)=gd(iad+i1)+bc(jtmp)                          1d11s21
                  jtmp=jtmp+1                                             1d11s21
                 end do                                                   1d11s21
                 i10=1                                                    1d11s21
                end do                                                    1d11s21
               end do                                                     1d11s21
c
c     Gvv'ri=dij Hvv" Vv'v"rj ... if sym v = sym v' = sym v" we get this
c     if v"=v' > Gvv'ri=dij Hvv' Vv'v'rj.
c     if v=v' > Gvvri dij Hvv" Vvv"rj
               if(isbv12.eq.1.and.l.eq.1.and.nrow.gt.0)then              4d28s21
                idenhuu=idenhvv(1,isb,jsb)+ibsft                        5d18s23
                call dgemm('n','n',nrow,nfdat(2,1,isb),nfdat(2,1,isb),
     $         sr2,vd(ioffs),nrow,bc(idenhuu),nfdat(2,1,isb),           5d18s23
     $              0d0,bc(itmp),nrow,                                  1d12s21
     d' hcddjk. 25')
                do i=0,nfdat(2,l,isb)-1                                   1d11s21
                 i10=i1s                                                   1d11s21
                 i1n=nvv                                                   1d11s21
                 do i2=i2s,i2e                                             1d11s21
                  ir=i2-1                                                 1d11s21
                  if(i2.eq.i2e)i1n=i1e                                    1d11s21
                  do i1=i10,i1n                                           1d11s21
                   i1m=i1-1                                               1d11s21
c     i1m=((iv2*(iv2-1))/2)+iv1,
c     2*i1m ge iv2*(iv2-1)=(iv2-1/2)^2-1/4
                   iv2=i1m/nvirt(isbv1)                                   1d11s21
                   iv1=i1m-nvirt(isbv1)*iv2                               1d11s21
                   try=2d0*dfloat(i1m)+0.25d0                             1d11s21
                   iv2t=sqrt(try)+0.5d0                                   1d11s21
                   iv1t=i1m-((iv2t*(iv2t-1))/2)                           1d11s21
                   iv1=iv1t+isw*(iv1-iv1t)
                   iv2=iv2t+isw*(iv2-iv2t)
c     Gvv'ri=dij Hvv' Vv'v'rj
                   ih0=ih0av(isbv1)+iv1+irefo(isbv1)+nh0av(isbv1)*(      1d12s21
     $                irefo(isbv1)+iv2)                                 1d12s21
                   iad=ioff+i1m+nvv*(ir+nrootu*i)                        1d12s21
                   jtmp=itmp+iv2+nvirt(isbv1)*(ir+nrootu*i)              1d12s21
                   ktmp=itmp+iv1+nvirt(isbv1)*(ir+nrootu*i)              1d12s21
                   gd(iad)=gd(iad)+(bc(ktmp)+bc(jtmp))*bc(ih0)           1d12s21
                  end do                                                 1d12s21
                  i10=1                                                  4d21s21
                 end do                                                  1d12s21
                end do                                                   1d12s21
               end if                                                    1d12s21
c     Gvv'ri=dij Hvv" Vv"v'rj, Gvv'ri=dij Hv'v" Vvv"rj
c
               idenhuu=idenhvv(l,isb,jsb)+ibsft                         5d18s23
               call dgemm('n','n',nhere,nfdat(2,l,isb),nfdat(2,l,isb),   1d11s21
     $       tf,vd(ioff+il-1),nnn,bc(idenhuu),nfdat(2,l,isb),           5d18s23
     $             0d0,bc(itmp),nhere,                                  1d11s21
     d' hcddjk. 26')
               if(l.eq.1.and.isbv12.eq.1)then                            1d18s21
                do i=0,nfdat(2,l,isb)-1                                   1d11s21
                 jtmp=itmp+nhere*i                                        1d11s21
                 i10=i1s                                                   1d11s21
                 i1n=nvv                                                   1d11s21
                 do i2=i2s,i2e                                             1d11s21
                  ir=i2-1                                                 1d11s21
                  if(i2.eq.i2e)i1n=i1e                                    1d11s21
                  ltmp=jtmp                                              1d18s21
                  do i1=i10,i1n                                           1d11s21
                   i1m=i1-1                                               1d11s21
c     i1m=((iv2*(iv2-1))/2)+iv1,
c     2*i1m ge iv2*(iv2-1)=(iv2-1/2)^2-1/4
                   iv2=i1m/nvirt(isbv1)                                   1d11s21
                   iv1=i1m-nvirt(isbv1)*iv2                               1d11s21
                   try=2d0*dfloat(i1m)+0.25d0                             1d11s21
                   iv2t=sqrt(try)+0.5d0                                   1d11s21
                   iv1t=i1m-((iv2t*(iv2t-1))/2)                           1d11s21
                   iv1=iv1t+isw*(iv1-iv1t)
                   iv2=iv2t+isw*(iv2-iv2t)
c     if v=v' > Gvvri dij Hvv" Vvv"rj
                   ih0=ih0av(isbv1)+irefo(isbv1)+iv1+nh0av(isbv1)*(      1d18s21
     $                irefo(isbv1)+iv2)                                 1d18s21
                   iad=ioffs+iv1+nvirt(isbv1)*(ir+nrootu*i)              1d18s21
                   gd(iad)=gd(iad)+bc(ih0)*bc(jtmp)*sr2                  1d18s21
                   jtmp=jtmp+1                                           1d18s21
                  end do                                                 1d18s21
                  if(isbv1.eq.isbv2)then                                 1d18s21
                   jtmp=ltmp                                             1d18s21
                   do i1=i10,i1n                                           1d11s21
                    i1m=i1-1                                               1d11s21
c     i1m=((iv2*(iv2-1))/2)+iv1,
c     2*i1m ge iv2*(iv2-1)=(iv2-1/2)^2-1/4
                    iv2=i1m/nvirt(isbv1)                                   1d11s21
                    iv1=i1m-nvirt(isbv1)*iv2                               1d11s21
                    try=2d0*dfloat(i1m)+0.25d0                             1d11s21
                    iv2t=sqrt(try)+0.5d0                                   1d11s21
                    iv1t=i1m-((iv2t*(iv2t-1))/2)                           1d11s21
                    iv1=iv1t+isw*(iv1-iv1t)
                    iv2=iv2t+isw*(iv2-iv2t)
c     if v=v' > Gvvri dij Hvv" Vv"vrj
                    ih0=ih0av(isbv1)+irefo(isbv1)+iv1+nh0av(isbv1)*(      1d18s21
     $                irefo(isbv1)+iv2)                                 1d18s21
                    iad=ioffs+iv2+nvirt(isbv1)*(ir+nrootu*i)              1d18s21
                    gd(iad)=gd(iad)+bc(ih0)*bc(jtmp)*sr2                  1d18s21
                    jtmp=jtmp+1                                           1d18s21
                   end do                                                 1d18s21
                  end if                                                 1d18s21
                  i10=1                                                  4d21s21
                 end do                                                  1d18s21
                end do                                                   1d18s21
               end if                                                    1d18s21
               do i=0,nfdat(2,l,isb)-1                                   1d11s21
                jtmp=itmp+nhere*i                                        1d11s21
                i10=i1s                                                   1d11s21
                i1n=nvv                                                   1d11s21
                do i2=i2s,i2e                                             1d11s21
                 ir=i2-1                                                 1d11s21
                 if(i2.eq.i2e)i1n=i1e                                    1d11s21
                 do i1=i10,i1n                                           1d11s21
                  i1m=i1-1                                               1d11s21
c     i1m=((iv2*(iv2-1))/2)+iv1,
c     2*i1m ge iv2*(iv2-1)=(iv2-1/2)^2-1/4
                  iv2=i1m/nvirt(isbv1)                                   1d11s21
                  iv1=i1m-nvirt(isbv1)*iv2                               1d11s21
                  try=2d0*dfloat(i1m)+0.25d0                             1d11s21
                  iv2t=sqrt(try)+0.5d0                                   1d11s21
                  iv1t=i1m-((iv2t*(iv2t-1))/2)                           1d11s21
                  iv1=iv1t+isw*(iv1-iv1t)
                  iv2=iv2t+isw*(iv2-iv2t)
c     Gvv'ri=dij Hvv" Vv"v'rj
                  ntop=iv2-1                                             1d11s21
                  ntop=ntop+isw*(nvirt(isbv1)-1-ntop)                    1d11s21
                  ih0=ih0av(isbv1)+irefo(isbv1)+nh0av(isbv1)*(           1d12s21
     $                irefo(isbv1)+iv1)                                 1d12s21
                  iad=ioff+nvv*(ir+nrootu*i)                             1d11s21
                  do iv=0,ntop                                           1d11s21
                   irec=iv+nvirt(isbv1)*iv2                              1d11s21
                   itri=((iv2*(iv2-1))/2)+iv                             1d18s21
                   irow=itri+isw*(irec-itri)                             1d11s21
                   gd(iad+irow)=gd(iad+irow)+bc(ih0+iv)*bc(jtmp)         1d11s21
                  end do                                                 1d11s21
c     Gvv'ri=dij Hv'v" Vvv"rj
                  nbot=iv1+1                                             1d11s21
                  nbot=nbot+isw*(0-nbot)                                 1d11s21
                  ih0=ih0av(isbv2)+irefo(isbv2)+nh0av(isbv2)*            1d12s21
     $                (irefo(isbv2)+iv2)                                1d12s21
                  do iv=nbot,nvirt(isbv2)-1                              1d11s21
                   irec=iv1+nvirt(isbv1)*iv                              1d11s21
                   itri=((iv*(iv-1))/2)+iv1                              1d18s21
                   irow=itri+isw*(irec-itri)                             1d11s21
                   gd(iad+irow)=gd(iad+irow)+bc(ih0+iv)*bc(jtmp)
                  end do                                                 1d11s21
                  if(isbv1.eq.isbv2)then
c     Gvv'ri=dij Hvv" Vv'v"rj
                   ntop=iv1-1                                             1d18s21
                   ntop=ntop+isw*(nvirt(isbv1)-1-ntop)                    1d18s21
                   ih0=ih0av(isbv1)+irefo(isbv1)+nh0av(isbv1)*(           1d12s21
     $                irefo(isbv1)+iv2)                                 1d18s21
                   iad=ioff+nvv*(ir+nrootu*i)                             1d11s21
                   do iv=0,ntop                                           1d11s21
                    irec=iv+nvirt(isbv1)*iv1                              1d11s21
                    itri=((iv1*(iv1-1))/2)+iv                             1d18s21
                    irow=itri+isw*(irec-itri)                             1d11s21
                    gd(iad+irow)=gd(iad+irow)+bc(ih0+iv)*bc(jtmp)*tf      1d18s21
                   end do                                                 1d11s21
c     Gvv'ri=dij Hv'v" Vv"vrj                                           1d18s21
                   nbot=iv2+1                                             1d18s21
                   nbot=nbot+isw*(0-nbot)                                 1d18s21
                   ih0=ih0av(isbv2)+irefo(isbv2)+nh0av(isbv2)*            1d18s21
     $                  (irefo(isbv2)+iv1)                                1d18s21
                   do iv=nbot,nvirt(isbv2)-1                              1d18s21
                    irec=iv2+nvirt(isbv1)*iv                              1d18s21
                    itri=((iv*(iv-1))/2)+iv2                              1d18s21
                    irow=itri+isw*(irec-itri)                             1d18s21
                    gd(iad+irow)=gd(iad+irow)+bc(ih0+iv)*bc(jtmp)*tf      1d18s21
                   end do                                                 1d11s21
                  end if                                                 1d18s21
                  jtmp=jtmp+1                                            1d11s21
                 end do                                                  1d11s21
                 i10=1                                                   1d11s21
                end do                                                   1d11s21
               end do                                                    1d11s21
               ibcoff=itmp                                               1d12s21
              end if                                                     1d11s21
              ioff=ioff+nnn*nfdat(2,l,isb)                               1d11s21
              ioffd=ioffd+nvv*nfdat(2,l,isb)                             1d11s21
             end if                                                      1d11s21
             tf=-1d0                                                     1d11s21
            end do                                                       1d11s21
           end if                                                        1d11s21
          end do                                                         1d11s21
         end if                                                          1d11s21
        end if                                                          5d10s21
        ibcoff=ibc0                                                     1d8s21
        do jsbv1=1,nsymb                                                1d8s21
         jsbv2=multh(jsbv1,jsbv12)                                      1d8s21
         if(jsbv2.ge.jsbv1)then                                         1d8s21
          if(jsbv12.eq.1)then                                           1d8s21
           joffdnon=joffdnon+nrootu*nvirt(jsbv1)*nfdat(2,1,jsb)         1d8s21
           nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                        1d8s21
          else                                                          1d8s21
           nvv=nvirt(jsbv1)*nvirt(jsbv2)                                1d8s21
          end if                                                        1d8s21
          nvv=nvv*nrootu                                                1d8s21
          do l=1,4                                                      1d8s21
           joffdnon=joffdnon+nvv*nfdat(2,l,jsb)                         1d8s21
          end do                                                        1d8s21
         end if                                                         1d8s21
        end do                                                          1d8s21
       end do                                                           1d8s21
       do isbv1=1,nsymb                                                 1d8s21
        isbv2=multh(isbv1,isbv12)                                       1d8s21
        if(isbv2.ge.isbv1)then                                          1d8s21
         if(isbv12.eq.1)then                                            1d8s21
          ioffdnon=ioffdnon+nvirt(isbv1)*nrootu*nfdat(2,1,isb)          1d8s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         1d8s21
         else                                                           1d8s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 1d8s21
         end if                                                         1d8s21
         nvv=nvv*nrootu                                                 1d8s21
         do l=1,4                                                       1d8s21
          ioffdnon=ioffdnon+nvv*nfdat(2,l,isb)                          1d8s21
         end do                                                         1d8s21
        end if                                                          1d8s21
       end do                                                           1d8s21
      end do                                                            1d8s21
      return                                                            1d8s21
      end                                                               1d8s21
