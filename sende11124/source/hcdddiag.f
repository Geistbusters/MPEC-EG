c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcdddiag(nff2,nsymb,mdon,mdoo,                         11d25s20
     $     multh,isymmrci,nvirt,ncsf,nec,ncsf2,nfdat,                   11d19s20
     $     irel,ism,ih0av,nh0av,ioooo,shift,kmats,jmats,irefo,          11d19s20
     $     ismult,gb,ixw1,ixw2,norb,lprt,bc,ibc)                        11d10s22
      implicit real*8 (a-h,o-z)                                         10d19s20
      integer*8 itesta,itestb,jtesta,jtestb,ipack8                      11d20s20
      integer*1 idorb(64),isorb(64),nab(2),iorb(64),imap(64),icode(64)    10d19s20
      logical ldebug,lbail,lprt                                         1d25s21
      dimension nff2(mdoo+1,nsymb,*),multh(8,8),nvirt(*),               11d25s20
     $     ncsf(*),ncsf2(4,*),nfdat(5,4,*),irel(*),ism(*),ih0av(*),     11d19s20
     $     nh0av(*),ioooo(*),jvcv(8),iovr(4,8),ipack4(2),npre(4),nl(4), 11d20s20
     $     kmats(*),jmats(*),ioffb(4),irefo(*),gb(*),iaddg(8,8),mpre(4),11d22s20
     $     iham1to4(4,8),idenvisv(8),idenvnotvj(4,8,8),                 11d23s20
     $     idenvnotvk(4,8,8),imatl(4),itrans(4),ntril(4),kvcvl(4),      11d23s20
     $     switch(5),ntrill(4,8)                                        1d22s21
      equivalence (ipack8,ipack4)                                       11d20s20
      include "common.store"
      common/fnd2cm/inv(2,8,8,8)                                        4d9s18
      common/kmfind/invk1(2,8,8,8,2)                                    6d30s18
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      ldebug=.false.                                                    10d19s20
      norbx=norb+1                                                      11d19s20
      norbxx=norbx+1                                                    11d19s20
      loop=0
      ismultm=ismult-1                                                  10d19s20
      ibcoffo=ibcoff                                                    11d22s20
      iham1=ibcoff                                                      11d23s20
      ntri11=(nfdat(2,1,isymmrci)*(nfdat(2,1,isymmrci)+1))/2             11d24s20
      ibcoff=iham1+ntri11                                                11d23s20
      if(lprt)write(6,*)('in hcddd '),ibcoff,ntri11                     1d25s21
      do isb=1,nsymb                                                    10d19s20
       jvcv(isb)=nfdat(5,1,isb)                                         11d20s20
       idenvisv(isb)=ibcoff                                             11d23s20
       ibcoff=idenvisv(isb)+irefo(isb)*ntri11                            11d23s20
       do l=1,4
        iham1to4(l,isb)=ibcoff                                          11d23s20
        ntrill(l,isb)=(nfdat(2,l,isb)*(nfdat(2,l,isb)+1))/2             11d25s20
        ibcoff=iham1to4(l,isb)+ntrill(l,isb)                            11d25s20
        do jsb=1,nsymb                                                  11d23s20
         idenvnotvj(l,jsb,isb)=ibcoff                                   11d23s20
         idenvnotvk(l,jsb,isb)=idenvnotvj(l,jsb,isb)+irefo(jsb)         11d24s20
     $        *ntrill(l,isb)                                            11d25s20
         ibcoff=idenvnotvk(l,jsb,isb)+irefo(jsb)*ntrill(l,isb)          11d25s20
        end do                                                          11d23s20
       end do                                                           11d23s20
      end do                                                            10d19s20
      igoal=1
      do i=ibcoffo,ibcoff-1                                             11d22s20
       bc(i)=0d0                                                        11d22s20
      end do                                                            11d22s20
      ntogether=ibcoff-ibcoffo                                          2d19s21
      idoit=0                                                           2d19s21
      jgoal=ibcoffo+1109
      ibclast=ibcoff                                                    10d17s20
      do nclok=mdon,mdoo                                                   10d19s20
       if(ldebug)write(6,*)('for '),nclok,(' closed shells ')
       nclokp=nclok+1                                                         10d19s20
       iargk=nclokp-mdon                                                    10d19s20
       iargkp=iargk+1                                                   2d1s21
       nopenk=nec-2*nclokp
       nopenkp=nopenk+2                                                   10d19s20
       norbt=nopenkp+nclok                                              2d1s21
       do isb=1,nsymb                                                   10d19s20
        do if=1,nff2(nclokp,isb,7)                                      11d20s20
         if(ldebug)write(6,*)('for if = '),if
         if(mod(idoit,mynprocg).eq.mynowprog)then                       2d19s21
          ipack8=ibc(jvcv(isb))                                          11d20s20
          nclot=popcnt(ipack4(1))                                        11d20s20
          if(ldebug)then
           call dcbit(ipack4(1),norb,'closed')
           call dcbit(ipack4(2),norb,'open')
          end if
          if(nclot.ne.nclok)then                                         11d20s20
           write(6,*)('for isb = '),isb,(' if = '),if,(' nclot = '),
     $         nclot,(' is ne nclok = '),nclok
           write(6,*)ibc(jvcv(isb)),jvcv(isb)
           call dws_synca
           call dws_finalize
           stop
          end if
          ii=1                                                           11d21s20
          do i=1,norb                                                    11d21s20
           if(btest(ipack4(1),i))then                                    11d21s20
            idorb(ii)=i                                                  11d21s20
            ii=ii+1                                                      11d21s20
           end if                                                        11d21s20
          end do                                                         11d21s20
          ii=1                                                           11d21s20
          do i=1,norb                                                    11d21s20
           if(btest(ipack4(2),i))then                                    11d21s20
            isorb(ii)=i                                                  11d21s20
            ii=ii+1                                                      11d21s20
           end if                                                        11d21s20
          end do                                                         11d21s20
          sumc=0d0                                                         9d23s19
          do i=1,nclok                                                   11d23s20
           ig=irel(idorb(i))-1                                           11d23s20
           is=ism(idorb(i))                                              11d23s20
           iad=ih0av(is)+ig*(nh0av(is)+1)                                  7d15s19
           h0=bc(iad)
           do j=1,nclok                                                  11d23s20
            jg=irel(idorb(j))-1                                          11d23s20
            js=ism(idorb(j))                                             11d23s20
            xj=getint(ioooo,is,is,js,js,ig,ig,jg,jg,bc,ibc)              11d15s22
            xk=getint(ioooo,is,js,js,is,ig,jg,jg,ig,bc,ibc)              11d15s22
            sumc=sumc+2d0*xj-xk                                            7d15s19
           end do                                                          7d15s19
           sumc=sumc+h0*2d0                                                7d15s19
          end do                                                           7d15s19
          sumo=0d0                                                         7d15s19
          do i=1,nopenk                                                  11d23s20
           ig=irel(isorb(i))-1                                           11d23s20
           is=ism(isorb(i))                                              11d23s20
           iad=ih0av(is)+ig*(nh0av(is)+1)                                  7d15s19
           h0=bc(iad)
           sumo=sumo+h0                                                    7d15s19
           do j=1,nopenk                                                 11d23s20
            if(i.ne.j)then                                                 7d15s19
             jg=irel(isorb(j))-1                                         11d23s20
             js=ism(isorb(j))                                            11d23s20
             xj=getint(ioooo,is,is,js,js,ig,ig,jg,jg,bc,ibc)             11d15s22
             sumo=sumo+xj*0.5d0                                            7d15s19
            end if                                                         7d15s19
           end do                                                          7d15s19
          end do                                                           7d15s19
          sum=sumc+sumo+shift
          do i=1,nclok                                                   11d23s20
           is=ism(idorb(i))                                              11d23s20
           ig=irel(idorb(i))-1                                           11d23s20
           do j=1,nopenk                                                 11d23s20
            js=ism(isorb(j))                                             11d23s20
            jg=irel(isorb(j))-1                                          11d23s20
            xj=getint(ioooo,is,is,js,js,ig,ig,jg,jg,bc,ibc)              11d15s22
            xk=getint(ioooo,is,js,js,is,ig,jg,jg,ig,bc,ibc)              11d15s22
            sum=sum+xj*2d0-xk                                              7d15s19
           end do                                                          7d15s19
          end do                                                           7d15s19
          if(ldebug)write(6,*)('what we have for sum: '),sum
          if(isb.eq.isymmrci)then                                           11d23s20
c
c     v is v
c
           nspace=ibc(jvcv(isb)+1)                                       3d19s21
           n1=ibc(jvcv(isb)+2)                                           3d19s21
           if(n1.gt.0)then                                               11d23s20
            ntri1=(n1*(n1+1))/2                                          11d23s20
            imat1=ibcoff                                                  11d23s20
            itmpt=imat1+ntri1                                            11d23s20
            ibcoff=itmpt+n1*ncsf2(1,iargk)                               11d23s20
            call enough('hcdddiag.  1',bc,ibc)
            do i=imat1,itmpt-1                                           11d23s20
             bc(i)=0d0                                                   11d23s20
            end do                                                       11d23s20
            iad1=jvcv(isb)+ibc(jvcv(isb)+6)                              3d19s21
            kvcv=iad1+n1                                                 3d19s21
            do i=0,n1-1                                                  11d23s20
             do j=0,ncsf2(1,iargk)-1                                     11d23s20
              ji=kvcv+j+ncsf2(1,iargk)*i                                 11d23s20
              ij=itmpt+i+n1*j                                            11d23s20
              bc(ij)=bc(ji)                                              11d23s20
             end do                                                      11d23s20
            end do                                                       11d23s20
            call mxmt(bc(itmpt),n1,bc(kvcv),ncsf2(1,iargk),bc(imat1))    11d23s20
            ji=imat1                                                     11d23s20
            do i=0,n1-1                                                  11d23s20
             ii=ibc(iad1+i)                                              11d23s20
             iii=((ii*(ii-1))/2)+iham1-1                                 11d23s20
             do j=0,i                                                    11d23s20
              jj=ibc(iad1+j)                                             11d23s20
              jjii=iii+jj                                                11d23s20
              bc(jjii)=bc(jjii)+bc(ji+j)*sum                             11d23s20
             end do                                                      11d23s20
             ji=ji+i+1                                                   11d23s20
            end do                                                       11d23s20
            do j=1,nclok                                                 11d23s20
             js=ism(idorb(j))                                             8d5s19
             jg=irel(idorb(j))-1                                          8d5s19
             jden=idenvisv(js)+ntri11*jg                                 11d24s20
             ji=imat1                                                     11d23s20
             do i=0,n1-1                                                  11d23s20
              ii=ibc(iad1+i)                                              11d23s20
              iiij=((ii*(ii-1))/2)+jden-1                                11d23s20
              do k=0,i                                                   11d23s20
               jj=ibc(iad1+k)                                            11d23s20
               jjii=iiij+jj                                               11d23s20
               bc(jjii)=bc(jjii)+bc(ji+k)*4d0                            11d23s20
              end do                                                     11d23s20
              ji=ji+i+1                                                  11d23s20
             end do                                                      11d23s20
            end do                                                       11d23s20
            do j=1,nopenk                                                11d23s20
             js=ism(isorb(j))                                            8d2s19
             jg=irel(isorb(j))-1                                         8d2s19
             jden=idenvisv(js)+ntri11*jg                                 11d24s20
             ji=imat1                                                    11d23s20
             do i=0,n1-1                                                 11d23s20
              ii=ibc(iad1+i)                                             11d23s20
              iiij=((ii*(ii-1))/2)+jden-1                                11d23s20
              do k=0,i                                                   11d23s20
               jj=ibc(iad1+k)                                            11d23s20
               jjii=iiij+jj                                               11d23s20
               bc(jjii)=bc(jjii)+bc(ji+k)*2d0                            11d23s20
              end do                                                     11d23s20
              ji=ji+i+1                                                  11d23s20
             end do                                                      11d23s20
            end do                                                       11d23s20
            if(nopenk.gt.0)then                                          11d28s20
             nn=ncsf2(1,iargk)*ncsf2(1,iargk)                             11d23s20
             itmp=ibcoff                                                   11d22s20
             ibcoff=itmp+nn                                               11d23s20
             call enough('hcdddiag.  2',bc,ibc)
             do i=itmp,ibcoff-1                                          2d1s21
              bc(i)=0d0                                                  1d31s21
             end do                                                      1d31s21
             itestb=0                                                        11d19s20
             do i=1,nopenk                                                  11d19s20
              itestb=ibset(itestb,i)                                        11d19s20
             end do                                                          11d19s20
             itesta=0                                                       11d19s20
             do i=1,nclok+1                                                    11d19s20
              itesta=ibset(itesta,i+nopenk)                                 11d19s20
             end do                                                          11d19s20
             do i1=1,nopenk-1                                              11d22s20
              jsa=ism(isorb(i1))                                           11d22s20
              jga=irel(isorb(i1))-1                                        11d22s20
              do i2=i1+1,nopenk                                           11d23s20
               jsb=ism(isorb(i2))                                          11d22s20
               jgb=irel(isorb(i2))-1                                       11d22s20
               xint=getint(ioooo,jsa,jsb,jsa,jsb,jga,jgb,jga,jgb,bc,ibc) 11d15s22
               kout=ibcoff                                               2d1s21
               ibcoff=kout+nn                                            2d1s21
               call enough('hcdddiag.  3',bc,ibc)
               jtesta=itesta                                                11d19s20
               jtestb=itestb                                                11d19s20
               nopenj=nopenk-2                                              11d19s20
               jarg=iargkp+1                                                11d19s20
               nqq=jarg+mdon-1
               if(nqq.ge.mdon.and.nqq.le.mdoo)then                          11d19s20
                jtestb=ibclr(jtestb,i2)                                      11d19s20
                jtestb=ibclr(jtestb,i1)                                      11d19s20
                jtesta=ibset(jtesta,i1)                                      11d19s20
                call gandc(itesta,itestb,jtesta,jtestb,nopenk,nopenj,        11d19s20
     $         iargkp,jarg,ncsf,norbt,ixw1,ixw2,nnot,nab,iwpb,iwpk,      11d19s20
     $         ncsfmid,bc,ibc)                                          11d14s22
                itmp1=ibcoff                                                11d19s20
                itmp2=itmp1+ncsf2(1,iargk)*ncsf(jarg)                       11d19s20
                ibcoff=itmp2+ncsf2(1,iargk)*ncsf(jarg)                         11d13s20
                call enough('hcdddiag.  4',bc,ibc)
                call prodn(iwpb,iwpk,ncsf2(1,iargk),ncsf(jarg),ncsfmid,     12d4s20
     $           bc(itmp1),bc,ibc,1d0,0d0)                              2d13s23
                do k=0,ncsf(jarg)-1                                         11d13s20
                 do j=0,ncsf2(1,iargk)-1                                       11d13s20
                  jk=itmp1+j+ncsf2(1,iargk)*k                                  11d13s20
                  kj=itmp2+k+ncsf(jarg)*j                                   11d13s20
                  bc(kj)=bc(jk)                                             11d13s20
                 end do                                                     11d13s20
                end do                                                      11d13s20
                call dgemm('n','n',ncsf2(1,iargk),ncsf2(1,iargk),
     $               ncsf(jarg),1d0,bc(itmp1),ncsf2(1,iargk),bc(itmp2),
     $           ncsf(jarg),0d0,bc(kout),ncsf2(1,iargk),                11d23s20
     d' hcdddiag.  1')
                ibcoff=itmp1                                                11d19s20
               else                                                         11d19s20
                do i=0,ncsf2(1,iargk)*ncsf2(1,iargk)-1                           11d19s20
                 bc(kout+i)=0d0                                             11d23s20
                end do                                                      11d19s20
               end if                                                       11d19s20
               do i=0,ncsf2(1,iargk)-1                                      11d19s20
                ii=kout+i*(ncsf2(1,iargk)+1)                                11d23s20
                bc(ii)=bc(ii)-1d0                                           11d19s20
               end do                                                       11d22s20
               do i=0,nn-1                                               2d1s21
                bc(itmp+i)=bc(itmp+i)+xint*bc(kout+i)                    2d1s21
               end do                                                    2d1s21
               ibcoff=kout                                               2d1s21
              end do                                                       11d22s20
             end do                                                        11d22s20
             itmp2=ibcoff                                                 11d23s20
             ibcoff=itmp2+ncsf2(1,iargk)*n1                               11d23s20
             call enough('hcdddiag.  5',bc,ibc)
             call dgemm('n','n',ncsf2(1,iargk),n1,ncsf2(1,iargk),1d0,     11d23s20
     $           bc(itmp),ncsf2(1,iargk),bc(kvcv),ncsf2(1,iargk),0d0,   11d23s20
     $           bc(itmp2),ncsf2(1,iargk),                              11d23s20
     d' hcdddiag.  2')
             do i=0,ntri1-1                                               11d23s20
              bc(imat1+i)=0d0                                             11d23s20
             end do                                                       11d23s20
             call mxmt(bc(itmpt),n1,bc(itmp2),ncsf2(1,iargk),bc(imat1))   11d23s20
             ji=imat1                                                     11d23s20
             do i=0,n1-1                                                  11d23s20
              ii=ibc(iad1+i)                                              11d23s20
              iii=((ii*(ii-1))/2)+iham1-1                                 11d23s20
              do j=0,i                                                    11d23s20
               jj=ibc(iad1+j)                                             11d23s20
               jjii=iii+jj                                                11d23s20
               bc(jjii)=bc(jjii)+bc(ji+j)                                 11d23s20
              end do                                                      11d23s20
              ji=ji+i+1                                                   11d23s20
             end do                                                       11d23s20
            end if                                                       11d28s20
            ibcoff=imat1                                                 11d23s20
           end if                                                        11d23s20
          end if                                                         11d23s20
         end if                                                         2d19s21
c
c     v not v
c
         nspace=ibc(jvcv(isb)+1)                                        3d19s21
         if(mod(idoit,mynprocg).eq.mynowprog)then                       2d19s21
          ibcb4=ibcoff                                                   11d23s20
          if(ldebug)then
           do l=1,4
            nl(l)=ibc(jvcv(isb)+1+l)                                     3d19s21
            if(nl(l).gt.0)then
             write(6,*)('l = '),l,nl(l)
            end if
           end do
          end if
          do l=1,4
           nl(l)=ibc(jvcv(isb)+1+l)                                      3d19s21
           if(nl(l).gt.0)then
            iadu=jvcv(isb)+ibc(jvcv(isb)+5+l)                            3d19s21
            iad2=iadu+nl(l)                                              3d19s21
            kvcvl(l)=iad2                                                11d23s20
            itrans(l)=ibcoff                                                  11d20s20
            ibcoff=itrans(l)+nl(l)*ncsf2(l,iargk)                          11d20s20
            call enough('hcdddiag.  6',bc,ibc)
            do i=0,nl(l)-1                                                    11d20s20
             do j=0,ncsf2(l,iargk)-1                                        11d20s20
              ji=iad2+j+ncsf2(l,iargk)*i                                    11d20s20
              ij=itrans(l)+i+nl(l)*j                                     11d23s20
              bc(ij)=bc(ji)                                               11d20s20
             end do                                                        11d20s20
            end do                                                         11d20s20
            ntril(l)=(nl(l)*(nl(l)+1))/2                                 11d23s20
            imatl(l)=ibcoff                                              11d23s20
            ibcoff=imatl(l)+ntril(l)                                     11d23s20
            call enough('hcdddiag.  7',bc,ibc)
            do i=imatl(l),ibcoff-1                                       11d23s20
             bc(i)=0d0                                                   11d23s20
            end do                                                       11d23s20
            call mxmt(bc(itrans(l)),nl(l),bc(kvcvl(l)),ncsf2(l,iargk),   11d23s20
     $           bc(imatl(l)))                                          11d23s20
            ji=imatl(l)                                                  11d23s20
            do i=0,nl(l)-1                                               11d23s20
             ii=ibc(iadu+i)                                              11d23s20
             iii=((ii*(ii-1))/2)+iham1to4(l,isb)-1                       11d23s20
             do j=0,i                                                    11d23s20
              jj=ibc(iadu+j)                                             11d23s20
              jjii=iii+jj                                                11d23s20
              bc(jjii)=bc(jjii)+bc(ji+j)*sum                             11d23s20
             end do                                                      11d23s20
             ji=ji+i+1                                                   11d23s20
            end do                                                       11d23s20
            do j=1,nclok                                                 11d23s20
             js=ism(idorb(j))                                            8d5s19
             jg=irel(idorb(j))-1                                         8d5s19
             jden=idenvnotvj(l,js,isb)+ntrill(l,isb)*jg                  11d25s20
             kden=idenvnotvk(l,js,isb)+ntrill(l,isb)*jg                  11d25s20
             ji=imatl(l)                                                 11d23s20
             do i=0,nl(l)-1                                                 11d23s20
              ii=ibc(iadu+i)                                             11d23s20
              iiij=((ii*(ii-1))/2)+jden-1                                11d23s20
              iiik=((ii*(ii-1))/2)+kden-1                                11d23s20
              do k=0,i                                                   11d23s20
               jj=ibc(iadu+k)                                            11d23s20
               jjii=iiij+jj                                               11d23s20
               bc(jjii)=bc(jjii)+bc(ji+k)*2d0                            11d23s20
               jjii=iiik+jj                                               11d23s20
               bc(jjii)=bc(jjii)-bc(ji+k)                                11d23s20
              end do                                                     11d23s20
              ji=ji+i+1                                                  11d23s20
             end do                                                      11d23s20
            end do                                                       11d23s20
            do j=1,nopenk                                                11d23s20
             js=ism(isorb(j))                                            8d2s19
             jg=irel(isorb(j))-1                                         8d2s19
             jden=idenvnotvj(l,js,isb)+ntrill(l,isb)*jg                  11d25s20
             ji=imatl(l)                                                    11d23s20
             do i=0,nl(l)-1                                                 11d23s20
              ii=ibc(iadu+i)                                             11d23s20
              iiij=((ii*(ii-1))/2)+jden-1                                11d23s20
              do k=0,i                                                   11d23s20
               jj=ibc(iadu+k)                                            11d23s20
               jjii=iiij+jj                                               11d23s20
               bc(jjii)=bc(jjii)+bc(ji+k)                                11d23s20
              end do                                                     11d23s20
              ji=ji+i+1                                                  11d23s20
             end do                                                      11d23s20
            end do                                                       11d23s20
           end if                                                        11d23s20
          end do                                                         11d23s20
          nnl=ncsf(iargk)*ncsf(iargk)                                    11d23s20
          nll=ncsf2(1,iargk)**2+ncsf2(2,iargk)**2+ncsf2(3,iargk)**2      11d23s20
     $         +ncsf2(4,iargk)**2                                        11d23s20
          itmp=ibcoff                                                    11d23s20
          itmpo=itmp+nll                                                 11d23s20
          ibcoff=itmpo+nll                                               11d23s20
          call enough('hcdddiag.  8',bc,ibc)
          do i=itmp,itmpo-1                                              11d23s20
           bc(i)=0d0                                                     11d23s20
          end do                                                         11d23s20
          itestb=0                                                        11d19s20
          do i=1,nopenkp                                                  11d19s20
           itestb=ibset(itestb,i)                                        11d19s20
          end do                                                          11d19s20
          itesta=0                                                       11d19s20
          do i=1,nclok                                                    11d19s20
           itesta=ibset(itesta,i+nopenkp)                                11d19s20
          end do                                                          11d19s20
          do i1=1,nopenk                                                 11d23s20
           jsa=ism(isorb(i1))                                            11d22s20
           jga=irel(isorb(i1))-1                                         11d22s20
           do i2=i1+1,nopenkp                                            11d19s20
            jtmpo=itmpo                                                  11d23s20
            kout=ibcoff                                                  2d1s21
            ibcoff=kout+nnl                                              2d1s21
            call enough('hcdddiag.  9',bc,ibc)
            jtesta=itesta                                                11d19s20
            jtestb=itestb                                                11d19s20
            nopenj=nopenkp-2                                             11d19s20
            jarg=iargk+1                                                 11d19s20
            nqq=jarg+mdon-1
            if(nqq.ge.mdon.and.nqq.le.mdoo)then                          11d19s20
             jtestb=ibclr(jtestb,i2)                                     11d19s20
             jtestb=ibclr(jtestb,i1)                                     11d19s20
             jtesta=ibset(jtesta,i1)                                     11d19s20
             call gandc(itesta,itestb,jtesta,jtestb,nopenkp,nopenj,      11d19s20
     $         iargk,jarg,ncsf,norbt,ixw1,ixw2,nnot,nab,iwpb,iwpk,      11d19s20
     $         ncsfmid,bc,ibc)                                          11d14s22
             itmp1=ibcoff                                                11d19s20
             itmp2=itmp1+ncsf(iargk)*ncsf(jarg)                          11d13s20
             ibcoff=itmp2+ncsf(iargk)*ncsf(jarg)                         11d13s20
             call enough('hcdddiag. 10',bc,ibc)
             call prodn(iwpb,iwpk,ncsf(iargk),ncsf(jarg),ncsfmid,        12d4s20
     $            bc(itmp1),bc,ibc,1d0,0d0)                             2d13s23
             do k=0,ncsf(jarg)-1                                         11d13s20
              do j=0,ncsf(iargk)-1                                       11d13s20
               jk=itmp1+j+ncsf(iargk)*k                                  11d13s20
               kj=itmp2+k+ncsf(jarg)*j                                   11d13s20
               bc(kj)=bc(jk)                                             11d13s20
              end do                                                     11d13s20
             end do                                                      11d13s20
             call dgemm('n','n',ncsf(iargk),ncsf(iargk),ncsf(jarg),1d0,  11d19s20
     $            bc(itmp1),ncsf(iargk),bc(itmp2),ncsf(jarg),0d0,        11d13s20
     $            bc(kout),ncsf(iargk),                                 11d23s20
     d' hcdddiag.  3')
             ibcoff=itmp1                                                11d19s20
            else                                                         11d19s20
             do i=0,ncsf(iargk)*ncsf(iargk)-1                            11d19s20
              bc(kout+i)=0d0                                             11d23s20
             end do                                                      11d19s20
            end if                                                       11d19s20
            do i=0,ncsf(iargk)-1                                         11d19s20
             ii=kout+i*(ncsf(iargk)+1)                                   11d23s20
             bc(ii)=bc(ii)-1d0                                           11d19s20
            end do                                                       11d19s20
            loff=0                                                       11d23s20
            do l=1,4                                                     11d23s20
             do i=0,ncsf2(l,iargk)-1                                     11d23s20
              ip=i+loff                                                  11d23s20
              do j=0,ncsf2(l,iargk)-1                                    11d23s20
               jp=j+loff                                                 11d23s20
               ifrom=kout+jp+ncsf(iargk)*ip                              2d1s21
               ito=jtmpo+j+ncsf2(l,iargk)*i                              11d23s20
               bc(ito)=bc(ifrom)                                         11d23s20
              end do                                                     11d23s20
             end do                                                      11d23s20
             loff=loff+ncsf2(l,iargk)                                    11d23s20
             jtmpo=jtmpo+ncsf2(l,iargk)*ncsf2(l,iargk)                   11d23s20
            end do                                                       11d23s20
            if(i2.le.nopenk)then                                         11d23s20
             jsb=ism(isorb(i2))                                          11d22s20
             jgb=irel(isorb(i2))-1                                       11d22s20
             xint=getint(ioooo,jsa,jsb,jsb,jsa,jga,jgb,jgb,jga,bc,ibc)   11d15s22
             do i=0,nll-1                                                11d23s20
              bc(itmp+i)=bc(itmp+i)+xint*bc(itmpo+i)                     11d23s20
             end do                                                      11d23s20
            else if(i2.eq.nopenk+1)then                                  11d23s20
             jtmpo=itmpo                                                 11d23s20
             do l=1,4                                                    11d23s20
              if(min(nl(l),ncsf2(l,iargk)).gt.0)then                     1d22s21
               itmpl=ibcoff                                              11d23s20
               itmpt=itmpl+ncsf2(l,iargk)*nl(l)                          11d23s20
               ibcoff=itmpt+ntril(l)                                     11d23s20
               call enough('hcdddiag. 11',bc,ibc)
               call dgemm('n','n',ncsf2(l,iargk),nl(l),ncsf2(l,iargk),   11d23s20
     $              1d0,bc(jtmpo),ncsf2(l,iargk),bc(kvcvl(l)),
     $              ncsf2(l,iargk),0d0,bc(itmpl),ncsf2(l,iargk),        11d23s20
     d' hcdddiag.  4')
               do i=0,ntril(l)-1                                         11d23s20
                bc(itmpt+i)=0d0                                          11d23s20
               end do                                                    11d23s20
               call mxmt(bc(itrans(l)),nl(l),bc(itmpl),ncsf2(l,iargk),   11d23s20
     $              bc(itmpt))                                          11d23s20
               ji=itmpt                                                  11d23s20
               iadu=jvcv(isb)+ibc(jvcv(isb)+5+l)                         3d19s21
               idenvk=idenvnotvk(l,jsa,isb)+ntrill(l,isb)*jga            11d25s20
               do i=0,nl(l)-1                                              11d23s20
                ii=ibc(iadu+i)                                             11d23s20
                iii=((ii*(ii-1))/2)-1+idenvk                             11d24s20
                do j=0,i                                                   11d23s20
                 jj=ibc(iadu+j)                                            11d23s20
                 jjii=iii+jj                                               11d23s20
                 bc(jjii)=bc(jjii)+bc(ji+j)                              11d23s20
                end do                                                     11d23s20
                ji=ji+i+1                                                  11d23s20
               end do                                                      11d23s20
               ibcoff=itmpl                                              11d23s20
              end if                                                     11d23s20
              jtmpo=jtmpo+ncsf2(l,iargk)**2                              12d28s20
             end do                                                      11d23s20
            end if                                                       11d23s20
            ibcoff=kout                                                  2d1s21
           end do
          end do
          jtmp=itmp                                                      11d23s20
          do l=1,4                                                       11d23s20
           if(min(nl(l),ncsf2(l,iargk)).gt.0)then                        11d23s20
            itmpl=ibcoff                                                 11d23s20
            itmpt=itmpl+ncsf2(l,iargk)*nl(l)                             11d23s20
            ibcoff=itmpt+ntril(l)                                        11d23s20
            call enough('hcdddiag. 12',bc,ibc)
            call dgemm('n','n',ncsf2(l,iargk),nl(l),ncsf2(l,iargk),      11d23s20
     $          1d0,bc(jtmp),ncsf2(l,iargk),bc(kvcvl(l)),               11d23s20
     $           ncsf2(l,iargk),0d0,bc(itmpl),ncsf2(l,iargk),           11d23s20
     d' hcdddiag.  5')
            do i=0,ntril(l)-1                                            11d23s20
             bc(itmpt+i)=0d0                                             11d23s20
            end do                                                       11d23s20
            call mxmt(bc(itrans(l)),nl(l),bc(itmpl),ncsf2(l,iargk),      11d23s20
     $          bc(itmpt))                                              11d23s20
            ji=itmpt                                                     11d23s20
            iadu=jvcv(isb)+ibc(jvcv(isb)+5+l)                            3d19s21
            do i=0,nl(l)-1                                               11d23s20
             ii=ibc(iadu+i)                                              11d23s20
             iii=((ii*(ii-1))/2)-1+iham1to4(l,isb)                       11d24s20
             do j=0,i                                                    11d23s20
              jj=ibc(iadu+j)                                             11d23s20
              jjii=iii+jj                                                11d23s20
              bc(jjii)=bc(jjii)+bc(ji+j)                                 11d23s20
             end do                                                      11d23s20
             ji=ji+i+1                                                   11d23s20
            end do                                                       11d23s20
            ibcoff=itmpl                                                 11d23s20
           end if                                                        11d23s20
           jtmp=jtmp+ncsf2(l,iargk)**2                                   11d23s20
          end do                                                         11d23s20
          ibcoff=ibcb4                                                   11d23s20
         end if                                                         2d19s21
         idoit=idoit+1                                                  2d19s21
         jvcv(isb)=jvcv(isb)+nspace                                     3d19s21
        end do                                                          11d23s20
       end do                                                           11d23s20
      end do                                                            11d23s20
      call dws_gsumf(bc(ibcoffo),ntogether)                             2d19s21
      if(nfdat(3,1,isymmrci).gt.0)then                                  11d17s21
       call sandtd(bc(iham1),nfdat(2,1,isymmrci),nfdat(3,1,isymmrci),    11d23s20
     $     1,bc(nfdat(4,1,isymmrci)),bc,ibc)                            11d10s22
       do jsb=1,nsymb                                                    11d23s20
        if(irefo(jsb).gt.0)then                                          11d23s20
         call sandtd(bc(idenvisv(jsb)),nfdat(2,1,isymmrci),              11d23s20
     $      nfdat(3,1,isymmrci),irefo(jsb),bc(nfdat(4,1,isymmrci)),bc,  11d10s22
     $       ibc)                                                       11d10s22
        end if                                                           11d23s20
       end do                                                            11d23s20
      end if                                                            11d17s21
      do isb=1,nsymb                                                    11d23s20
       do l=1,4                                                          11d23s20
        if(nfdat(3,l,isb).gt.0)then                                      11d23s20
         call sandtd(bc(iham1to4(l,isb)),nfdat(2,l,isb),nfdat(3,l,isb), 11d23s20
     $       1,bc(nfdat(4,l,isb)),bc,ibc)                               11d10s22
         do jsb=1,nsymb                                                 11d23s20
          if(irefo(jsb).gt.0)then                                       11d23s20
           call sandtd(bc(idenvnotvj(l,jsb,isb)),nfdat(2,l,isb),        11d23s20
     $        nfdat(3,l,isb),irefo(jsb),bc(nfdat(4,l,isb)),bc,ibc)      11d10s22
           call sandtd(bc(idenvnotvk(l,jsb,isb)),nfdat(2,l,isb),        11d23s20
     $        nfdat(3,l,isb),irefo(jsb),bc(nfdat(4,l,isb)),bc,ibc)      11d10s22
          end if                                                        11d23s20
         end do                                                         11d23s20
        end if                                                           11d23s20
       end do                                                            11d23s20
      end do                                                            11d23s20
      ioff=1                                                            11d23s20
      do isb=1,nsymb                                                    11d23s20
       isbv12=multh(isb,isymmrci)                                       11d23s20
       ncdtt=nfdat(3,1,isb)+nfdat(3,2,isb)+nfdat(3,3,isb)+nfdat(3,4,isb)11d23s20
       do isbv1=1,nsymb                                                 11d23s20
        isbv2=multh(isbv12,isbv1)                                       11d23s20
        if(isbv1.le.isbv2)then                                          11d23s20
         call ilimts(nvirt(isbv1),nvirt(isbv1),mynprocg,mynowprog,il,     11d23s20
     $          ih,i1s,i1e,i2s,i2e)                                     10d19s20
         nhere=ih+1-il                                                  11d23s20
         nv=0                                                           11d23s20
         iv0=nvirt(isbv1)                                                11d23s20
         do iv=i2s,i2e                                                  11d23s20
          ivm=iv-1                                                      11d23s20
          irow=iv+nvirt(isbv1)*ivm                                       11d23s20
          if(irow.ge.il.and.irow.le.ih)then                             11d23s20
           iv0=min(iv0,ivm)                                             11d23s20
           nv=nv+1                                                      11d23s20
          end if                                                        11d23s20
         end do                                                         11d23s20
         nvm=nv-1                                                        11d23s20
         iadh1=ih0av(isbv1)+irefo(isbv1)+nh0av(isbv1)*irefo(isbv1)      11d23s20
         nadh1=nh0av(isbv1)+1                                           11d23s20
         if(isbv1.eq.isbv2)then                                         11d23s20
          ioffw=ioff                                                     11d23s20
          do i=0,nfdat(3,1,isb)-1                                       11d23s20
           do iv=mynowprog,nvirt(isbv1)-1,mynprocg                      11d23s20
            iad=iadh1+iv*nadh1                                          11d23s20
            gb(ioffw+iv)=gb(ioffw+iv)+bc(iham1+i)+2d0*bc(iad)
           end do                                                        11d23s20
           ioffw=ioffw+nvirt(isbv1)                                     11d23s20
          end do                                                         11d23s20
          ioffw=ioff+iv0                                                 11d23s20
          if(min(nfdat(3,1,isb),nv).gt.0)then                           11d17s21
           iquerry=ibcoff
           ibcoff=iquerry+nfdat(3,1,isb)*nv
           call enough('hcdddiag. 13',bc,ibc)
           do i=iquerry,ibcoff-1
            bc(i)=0d0
           end do
           do jsb=1,nsymb                                                 11d23s20
            if(irefo(jsb).gt.0)then                                       11d23s20
             i2euj=inv(1,jsb,jsb,isbv1)                                    11d23s20
             i2euk=invk1(1,jsb,jsb,isbv1,1)                                11d23s20
             itmpjk=ibcoff                                                11d23s20
             ibcoff=itmpjk+nv*irefo(jsb)                                    11d23s20
             call enough('hcdddiag. 14',bc,ibc)
             do ig=0,irefo(jsb)-1                                         11d23s20
              icolj=((ig*(ig+1))/2)+ig                                    11d23s20
              iadj=jmats(i2euj)+nhere*icolj-il+1                         11d23s20
              icolk=ig*(irefo(jsb)+1)                                     11d23s20
              iadk=kmats(i2euk)+nhere*icolk-il+1                         11d23s20
              ito=itmpjk+nv*ig                                            11d23s20
              do jv=0,nvm                                                 11d23s20
               iv=jv+iv0                                                11d24s20
               irow=iv+nvirt(isbv1)*iv                                   11d23s20
               bc(ito+jv)=bc(iadj+irow)
     $              -0.5d0*bc(iadk+irow)
              end do                                                      11d23s20
             end do                                                       11d23s20
             call dgemm('n','n',nv,nfdat(3,1,isb),irefo(jsb),1d0,        11d23s20
     $          bc(itmpjk),nv,bc(idenvisv(jsb)),irefo(jsb),1d0,         11d23s20
     $           gb(ioffw),nvirt(isbv1),                                  11d23s20
     d' hcdddiag.  6')
             ibcoff=itmpjk                                                11d23s20
            end if                                                        11d23s20
           end do                                                         11d23s20
          end if                                                        11d23s20
          ioff=ioff+nvirt(isbv1)*nfdat(3,1,isb)                         11d23s20
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         11d23s20
          isw=0                                                         11d23s20
          iadh2=iadh1                                                   11d23s20
          nadh2=nadh1                                                   11d23s20
          nherej=nhere                                                  11d23s20
          jl=il                                                         11d25s20
          nvj=nv                                                        11d23s20
          iv0j=iv0                                                      11d23s20
          nvmj=nvm                                                      11d23s20
         else                                                           11d23s20
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 11d23s20
          isw=1                                                         11d23s20
          iadh2=ih0av(isbv2)+irefo(isbv2)+nh0av(isbv2)*irefo(isbv2)      11d23s20
          nadh2=nh0av(isbv2)+1                                           11d23s20
          call ilimts(nvirt(isbv2),nvirt(isbv2),mynprocg,mynowprog,jl,   11d23s20
     $          jh,j1s,j1e,j2s,j2e)                                     10d19s20
          nherej=jh+1-jl                                                  11d23s20
          nvj=0                                                           11d23s20
          iv0j=nvirt(isbv2)                                                11d23s20
          do iv=j2s,j2e                                                  11d23s20
           ivm=iv-1                                                      11d23s20
           irow=iv+nvirt(isbv2)*ivm                                       11d23s20
           if(irow.ge.jl.and.irow.le.jh)then                             11d23s20
            iv0j=min(iv0j,ivm)                                             11d23s20
            nvj=nvj+1                                                      11d23s20
           end if                                                        11d23s20
          end do                                                         11d23s20
          nvmj=nvj-1                                                        11d23s20
         end if                                                         11d23s20
         if(nvv.gt.0)then                                               11d23s20
          ioffw=ioff                                                     11d23s20
          do l=1,4                                                       11d23s20
           do i=0,nfdat(3,l,isb)-1                                       11d23s20
            ih0p=iham1to4(l,isb)+i                                       11d23s20
            do iv2=mynowprog,nvirt(isbv2)-1,mynprocg                     11d23s20
             iv2m=iv2-1                                                  11d23s20
             itop=iv2m+isw*(nvirt(isbv1)-1-iv2m)                         11d23s20
             ih02=iadh2+iv2*nadh2                                        11d23s20
             do iv1=0,itop                                               11d23s20
              itri=((iv2*iv2m)/2)+iv1                                   11d24s20
              irec=iv1+nvirt(isbv1)*iv2                                  11d23s20
              itri=itri+isw*(irec-itri)                                 11d23s20
              ih01=iadh1+iv1*nadh1                                       11d23s20
              gb(ioffw+itri)=gb(ioffw+itri)+bc(ih0p)+bc(ih02)
     $             +bc(ih01)
             end do                                                      11d23s20
            end do
            ioffw=ioffw+nvv                                              11d23s20
           end do                                                        11d23s20
          end do                                                         11d23s20
          itmpp=ibcoff                                                   11d23s20
          itmppj=itmpp+nv*ncdtt                                          11d23s20
          ibcoff=itmppj+nvj*ncdtt                                        11d23s20
          call enough('hcdddiag. 15',bc,ibc)
          do i=itmpp,ibcoff-1                                            11d23s20
           bc(i)=0d0                                                     11d23s20
          end do                                                         11d23s20
          do jsb=1,nsymb                                                 11d23s20
           if(irefo(jsb).gt.0)then                                       11d23s20
            if(nv.gt.0)then                                             11d23s20
             i2euj=inv(1,jsb,jsb,isbv1)                                   11d23s20
             i2euk=invk1(1,jsb,jsb,isbv1,1)                               11d23s20
             itmpj=ibcoff                                                 11d23s20
             itmpk=itmpj+nv*irefo(jsb)                                    11d23s20
             ibcoff=itmpk+nv*irefo(jsb)                                   11d23s20
             call enough('hcdddiag. 16',bc,ibc)
             jtmpj=itmpj                                                  11d23s20
             jtmpk=itmpk                                                  11d23s20
             do ig=0,irefo(jsb)-1                                         11d23s20
              icolj=((ig*(ig+1))/2)+ig                                    11d23s20
              iadj=jmats(i2euj)+icolj*nhere+1-il                          11d23s20
              icolk=ig*(irefo(jsb)+1)                                     11d23s20
              iadk=kmats(i2euk)+icolk*nhere+1-il                          11d23s20
              do jv=0,nvm                                                 11d23s20
               iv=jv+iv0                                                11d24s20
               irow=iv*(nvirt(isbv1)+1)                                   11d23s20
               bc(jtmpj+jv)=bc(iadj+irow)
               bc(jtmpk+jv)=bc(iadk+irow)
              end do                                                      11d23s20
              jtmpj=jtmpj+nv                                              11d23s20
              jtmpk=jtmpk+nv                                              11d23s20
             end do                                                       11d23s20
             jtmpp=itmpp                                                  11d23s20
             do l=1,4                                                     11d23s20
              if(nfdat(3,l,isb).gt.0)then                                11d23s20
               call dgemm('n','n',nv,nfdat(3,l,isb),irefo(jsb),1d0,        11d23s20
     $           bc(itmpj),nv,bc(idenvnotvj(l,jsb,isb)),irefo(jsb),1d0, 11d23s20
     $           bc(jtmpp),nv,                                          11d23s20
     d' hcdddiag.  7')
               call dgemm('n','n',nv,nfdat(3,l,isb),irefo(jsb),1d0,        11d23s20
     $           bc(itmpk),nv,bc(idenvnotvk(l,jsb,isb)),irefo(jsb),1d0, 11d23s20
     $           bc(jtmpp),nv,                                          11d23s20
     d' hcdddiag.  8')
              end if                                                     11d23s20
              jtmpp=jtmpp+nv*nfdat(3,l,isb)                               11d23s20
             end do                                                       11d23s20
             ibcoff=itmpj                                                 11d23s20
            end if                                                      11d23s20
c
            if(nvj.gt.0)then                                            11d23s20
             i2euj=inv(1,jsb,jsb,isbv2)                                   11d23s20
             i2euk=invk1(1,jsb,jsb,isbv2,1)                               11d23s20
             itmpj=ibcoff                                                 11d23s20
             itmpk=itmpj+nvj*irefo(jsb)                                    11d23s20
             ibcoff=itmpk+nvj*irefo(jsb)                                   11d23s20
             call enough('hcdddiag. 17',bc,ibc)
             jtmpj=itmpj                                                  11d23s20
             jtmpk=itmpk                                                  11d23s20
             do ig=0,irefo(jsb)-1                                         11d23s20
              icolj=((ig*(ig+1))/2)+ig                                    11d23s20
              iadj=jmats(i2euj)+icolj*nherej+1-jl                         11d23s20
              icolk=ig*(irefo(jsb)+1)                                     11d23s20
              iadk=kmats(i2euk)+icolk*nherej+1-jl                         11d23s20
              do jv=0,nvmj                                                 11d23s20
               iv=jv+iv0j                                               11d24s20
               irow=iv*(nvirt(isbv2)+1)                                   11d23s20
               bc(jtmpj+jv)=bc(iadj+irow)
               bc(jtmpk+jv)=bc(iadk+irow)
              end do                                                      11d23s20
              jtmpj=jtmpj+nvj                                              11d23s20
              jtmpk=jtmpk+nvj                                              11d23s20
             end do                                                       11d23s20
             jtmppj=itmppj                                                  11d23s20
             do l=1,4                                                     11d23s20
              if(nfdat(3,l,isb).gt.0)then                                11d23s20
               call dgemm('n','n',nvj,nfdat(3,l,isb),irefo(jsb),1d0,        11d23s20
     $           bc(itmpj),nvj,bc(idenvnotvj(l,jsb,isb)),irefo(jsb),1d0, 11d23s20
     $           bc(jtmppj),nvj,                                          11d23s20
     d' hcdddiag.  9')
               call dgemm('n','n',nvj,nfdat(3,l,isb),irefo(jsb),1d0,        11d23s20
     $           bc(itmpk),nvj,bc(idenvnotvk(l,jsb,isb)),irefo(jsb),1d0, 11d23s20
     $           bc(jtmppj),nvj,                                          11d23s20
     d' hcdddiag. 10')
              end if                                                     11d23s20
              jtmppj=jtmppj+nvj*nfdat(3,l,isb)                               11d23s20
             end do                                                       11d23s20
             ibcoff=itmpj                                                 11d23s20
            end if                                                      11d23s20
           end if
          end do                                                         11d23s20
          do i=0,ncdtt-1                                                 11d23s20
           do jv2=0,nvmj                                                 11d23s20
            iv2=iv0j+jv2                                                11d24s20
            iv2m=iv2-1                                                   11d23s20
            itop=iv2m+isw*(nvirt(isbv1)-1-iv2m)                         11d23s20
            jtmppj=itmppj+jv2+nvj*i                                      11d23s20
            do iv1=0,itop                                               11d23s20
             itri=((iv2*iv2m)/2)+iv1                                    11d23s20
             irec=iv1+nvirt(isbv1)*iv2                                   11d23s20
             itri=itri+isw*(irec-itri)                                   11d23s20
             gb(ioff+itri)=gb(ioff+itri)+bc(jtmppj)                     11d23s20
            end do                                                       11d23s20
           end do                                                        11d23s20
           do iv2=0,nvirt(isbv2)-1                                      11d23s20
            iv2m=iv2-1                                                  11d23s20
            itop=min(iv2m,iv0+nvm)+isw*(iv0+nvm-min(iv2m,iv0+nvm))      11d25s20
            do iv1=iv0,itop                                             11d24s20
             itri=((iv2*iv2m)/2)+iv1                                    11d24s20
             irec=iv1+nvirt(isbv1)*iv2                                   11d23s20
             itri=itri+isw*(irec-itri)                                   11d23s20
             jtmpp=itmpp+iv1-iv0+nv*i                                   11d24s20
             gb(ioff+itri)=gb(ioff+itri)+bc(jtmpp)                      11d23s20
            end do                                                       11d23s20
           end do                                                        11d23s20
           ioff=ioff+nvv                                                 11d23s20
          end do                                                         11d23s20
          ibcoff=itmpp                                                   11d23s20
         end if                                                         11d23s20
        end if                                                           11d23s20
       end do                                                           11d23s20
      end do                                                            11d23s20
      ibcoff=ibcoffo                                                    11d23s20
      ioff=ioff-1                                                       11d26s20
      call dws_gsumf(gb,ioff)                                           11d26s20
      if(lprt)then                                                      1d25s21
       write(6,*)('all done!'),ioff
       call prntm2(gb,1,ioff,1)
       ioff=1
       do isb=1,nsymb
        isbv12=multh(isb,isymmrci)
        do isbv1=1,nsymb
         isbv2=multh(isbv1,isbv12)
         if(isbv2.ge.isbv1)then
          write(6,*)('for symmetries '),isb,isbv1,isbv2
          if(isbv1.eq.isbv2)then                                        1d22s21
           nvvs=(nvirt(isbv1)*(nvirt(isbv1)+1))/2                       1d22s21
           nvvt=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       1d22s21
           isw=0                                                        1d22s21
          else                                                          1d22s21
           nvvs=nvirt(isbv1)*nvirt(isbv2)                               1d22s21
           nvvt=nvvs                                                    1d22s21
           isw=1                                                        1d22s21
          end if                                                        1d22s21
          ip=1
          nvv=nvvs
          do l=1,4
           if(min(nvv,nfdat(3,l,isb)).gt.0)then
            itmp=ibcoff                                                 1d22s21
            ibcoff=itmp+nvv*nfdat(3,l,isb)
            call enough('hcdddiag. 18',bc,ibc)
            if(l.eq.1.and.isbv1.eq.isbv2)then
             do i=0,nfdat(3,l,isb)-1
              do iv=0,nvirt(isbv1)-1
               ivv=((iv*(iv+1))/2)+iv
               ito=itmp+i+nfdat(3,l,isb)*ivv
               bc(ito)=gb(ioff)
               ioff=ioff+1
              end do
             end do
             do i=0,nfdat(3,l,isb)-1                                    1d22s21
              do iv2=0,nvirt(isbv1)-1
               do iv1=0,iv2-1
                ivv=((iv2*(iv2+1))/2)+iv1
                ito=itmp+i+nfdat(3,l,isb)*ivv
                bc(ito)=gb(ioff)
                ioff=ioff+1
               end do
              end do
             end do
            else
             do i=0,nfdat(3,l,isb)-1
              do ivv=0,nvv-1
               ito=itmp+i+nfdat(3,l,isb)*ivv
               bc(ito)=gb(ioff)
               ioff=ioff+1
              end do
             end do
            end if
            if(l.eq.1)then
             write(6,*)('singlet coupled pairs ')
             call prntm2(bc(itmp),nfdat(3,l,isb),nvvs,nfdat(3,l,isb))
            else
             write(6,*)('triplet coupled pairs for spin '),l
             call prntm2(bc(itmp),nfdat(3,l,isb),nvvt,nfdat(3,l,isb))
            end if
            ibcoff=itmp
           end if
           ip=-1
           nvv=nvvt
          end do
         end if
        end do
       end do
       write(6,*)('ioff = '),ioff
      end if
      return                                                            10d19s20
      end                                                               10d19s20
