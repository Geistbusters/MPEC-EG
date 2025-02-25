c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hccsfd(vecx,ncsft,ibasis,ncsf,iptr,nfcn,iden1,iden2,   12d13s19
     $     nrootz,mdon,isorb,idorb,icsfpd,nec,wgt,istype,ism,irel,irefo,4d13s21
     $     ixw1,ixw2,iptrbit,mdoo,mysym,norb,bc,ibc)                    11d10s22
      implicit real*8 (a-h,o-z)                                         7d11s19
c
c     compute 1 and 2 electron densities
c     note on istype: if one is dealing with 2-e integrals with the
c     symmetry (ab|cd)=(ba|cd)=(ab|dc)=(ba|dc)=(cd|ab)= etc, then
c     istype should be zero.
c     if one is dealing with 2-e integrals with the symmetry
c     (ab|cd)=-(ba|cd)=-(ab|dc)=(ba|dc)=(cd|ab)=- etc, then istype
c     should be non-zero. An example of the latter is the Lz matrices
c     involved in the calculation of the Lz^2 operator.
c     pfr=+1 for symmetric operators, like two-electron integrals, and
c         -1 for skew symmetric operators, like lz
c
      external second                                                   8d1s19
      integer*1 isorb(*),idorb(*),icode(64),imap(64),nab(2),isame(64),  4d13s21
     $     nab1(2),nab2(2)                                              4d12s21
      integer*8 ipack,itesta,itestb                                     4d13s21
      integer*2 ipack2(4)                                               12d1s19
      equivalence (ipack,ipack2)                                        12d1s19
      logical ldebug                                                    1d3s20
      dimension vecx(ncsft,nrootz),iptr(4,*),wgt(*),iden1(*),           12d13s19
     $     iden2(*),ibasis(3,*),icsfpd(*),ncsf(*),nother(2),ism(*),     12d28s19
     $     irel(*),irefo(*),iptrbit(2,mdoo+1,*),nab4(2,3),itest(64,2)   4d13s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      common/fnd2cm/inv(2,8,8,8)                                        4d9s18
      include "common.store"                                            7d11s19
      include "common.print"                                            1d3s20
      if(iprtr(6).eq.0)then                                             1d3s20
       ldebug=.false.                                                   1d3s20
      else                                                              1d3s20
       ldebug=.true.                                                    1d3s20
      end if                                                            1d3s20
      if(ldebug)then                                                    1d3s20
      write(6,*)('Hi, my name is hccsfd '),ncsft,nfcn,nrootz,mdon,nec
      write(6,*)('mysym = '),mysym
      write(6,*)('input vector ')
      call prntm2(vecx,ncsft,nrootz,ncsft)
      write(6,*)('root weights: '),(wgt(i),i=1,nrootz)
      write(6,*)('what''s under ibasis: ')
      do i=1,nfcn
       write(6,*)i,(ibasis(j,i),j=1,3)
      end do
      end if                                                            1d3s20
      ips=0                                                             7d11s19
      loopit=0                                                          7d11s19
      do if=1,nfcn                                                      7d11s19
       nclo=ibasis(1,if)                                                7d11s19
       nclop=nclo+1                                                     7d11s19
       iarg=nclop-mdon                                                  7d11s19
       nopen=nec-2*nclo                                                 7d11s19
       ip=ips+1                                                         7d11s19
       iic=iptrbit(1,nclop,mysym)+ibasis(2,if)-1                        4d13s21
       iio=iptrbit(2,nclop,mysym)+ibasis(3,if)-1                        4d13s21
       if(mynowprog.eq.0)then                                           12d14s19
        ww=0d0                                                          12d17s19
        do ir=1,nrootz                                                  12d17s19
         do ix=0,ncsf(iarg)-1                                           12d17s19
          ww=ww+vecx(ip+ix,ir)*vecx(ip+ix,ir)*wgt(ir)                   12d17s19
         end do                                                         12d17s19
        end do                                                          12d17s19
        do i=1,norb                                                     8d2s22
         if(btest(ibc(iic),i))then                                      8d2s22
          is=ism(i)                                                     8d2s22
          ig=irel(i)-1                                                  8d2s22
          jden1=iden1(is)+ig*(irefo(is)+1)                                12d14s19
          bc(jden1)=bc(jden1)+2d0*ww                                     12d17s19
          irow=((ig*(ig+1))/2)+ig                                        12d17s19
          nrow=(irefo(is)*(irefo(is)+1))/2                               12d17s19
          do j=1,norb                                                   8d2s22
           if(btest(ibc(iic),j))then                                    8d2s22
            js=ism(j)                                                   8d2s22
            jg=irel(j)-1                                                8d2s22
            i2eu=inv(1,is,is,js)                                          12d17s19
            icol=((jg*(jg+1))/2)+jg                                       12d17s19
            jden2j=iden2(i2eu)+irow+nrow*icol                             12d17s19
            bc(jden2j)=bc(jden2j)+2d0*ww                                  12d17s19
            i2euk=inv(1,is,js,js)                                         12d17s19
            if(is.eq.js)then                                              12d17s19
             ix=max(ig,jg)                                                12d17s19
             in=min(ig,jg)                                                12d17s19
             ircowk=((ix*(ix+1))/2)+in                                     12d17s19
             jden2k=iden2(i2euk)+ircowk*(nrow+1)                          12d17s19
            else                                                          12d17s19
             icase=inv(2,is,js,js)                                        12d17s19
             nrowk=irefo(is)*irefo(js)                                    12d17s19
             if(icase.eq.1.or.icase.eq.3)then                                           12d17s19
              irowk=ig+irefo(is)*jg                                       12d17s19
             else
              irowk=jg+irefo(js)*ig                                       12d17s19
             end if                                                       12d17s19
             if(icase.eq.1.or.icase.eq.4)then                             12d17s19
              icolk=jg+irefo(js)*ig                                       12d17s19
             else                                                         12d17s19
              icolk=ig+irefo(is)*jg                                       12d17s19
             end if                                                       12d17s19
             jden2k=iden2(i2euk)+irowk+nrowk*icolk                        12d17s19
            end if                                                        12d17s19
            bc(jden2k)=bc(jden2k)-ww                                      12d27s19
           end if                                                       8d2s22
          end do                                                           12d14s19
         end if                                                         8d2s22
        end do
        do i=1,norb                                                     8d2s22
         if(btest(ibc(iio),i))then                                      8d2s22
          is=ism(i)                                                     8d2s22
          ig=irel(i)-1                                                  8d2s22
          jden1=iden1(is)+ig*(irefo(is)+1)                                12d14s19
          bc(jden1)=bc(jden1)+ww                                         12d17s19
          irow=((ig*(ig+1))/2)+ig                                        12d17s19
          nrow=(irefo(is)*(irefo(is)+1))/2                               12d17s19
          do j=1,norb                                                   8d2s22
           if(btest(ibc(iio),j).and.i.ne.j)then                         8d2s22
            js=ism(j)                                                   8d2s22
            jg=irel(j)-1                                                8d2s22
            i2eu=inv(1,is,is,js)                                         12d17s19
            icol=((jg*(jg+1))/2)+jg                                      12d17s19
            jden2j=iden2(i2eu)+irow+nrow*icol                            12d17s19
            bc(jden2j)=bc(jden2j)+ww*0.5d0                               12d17s19
           end if                                                        12d17s19
          end do                                                         12d17s19
          do j=1,norb                                                   8d2s22
           if(btest(ibc(iic),j))then                                    8d2s22
            js=ism(j)                                                   8d2s22
            jg=irel(j)-1                                                8d2s22
            i2eu=inv(1,is,is,js)                                          12d17s19
            icol=((jg*(jg+1))/2)+jg                                       12d17s19
            jden2j=iden2(i2eu)+irow+nrow*icol                             12d17s19
            bc(jden2j)=bc(jden2j)+ww*2d0                                  12d17s19
            i2euk=inv(1,is,js,js)                                         12d17s19
            if(is.eq.js)then                                              12d17s19
             ix=max(ig,jg)                                                12d17s19
             in=min(ig,jg)                                                12d17s19
             irc=((ix*(ix+1))/2)+in                                       12d17s19
             jden2k=iden2(i2euk)+irc*(nrow+1)                             12d17s19
            else                                                          12d17s19
             icase=inv(2,is,js,js)                                        12d17s19
             nrowk=irefo(is)*irefo(js)                                    12d17s19
             if(icase.eq.1.or.icase.eq.3)then                             12d17s19
              irowk=ig+irefo(is)*jg                                       12d17s19
             else                                                         12d17s19
              irowk=jg+irefo(js)*ig                                       12d17s19
             end if
             if(icase.eq.1.or.icase.eq.4)then                             12d17s19
              icolk=jg+irefo(js)*ig                                       12d17s19
             else                                                         12d17s19
              icolk=ig+irefo(is)*jg                                       12d17s19
             end if                                                       12d17s19
             jden2k=iden2(i2euk)+irowk+nrowk*icolk                        12d17s19
            end if                                                        12d17s19
            bc(jden2k)=bc(jden2k)-ww                                      12d27s19
           end if                                                       8d2s22
          end do                                                         12d17s19
         end if                                                         8d2s22
        end do                                                           12d14s19
       end if                                                           12d14s19
       jps=ips                                                          12d11s19
       do jf=if,nfcn                                                    12d9s19
        ncloj=ibasis(1,jf)                                              7d11s19
        if(iabs(ncloj-nclo).le.2)then                                   12d9s19
         if(jf.eq.if)then                                               12d13s19
          fact=1d0                                                      12d13s19
         else                                                           12d13s19
          fact=2d0                                                      12d13s19
         end if                                                         12d13s19
         if(ldebug)write(6,*)('for jf,if: '),jf,if,ncloj,mdon,intmul
         nclopj=ncloj+1                                                  6d11s19
         nopenj=nec-2*ncloj
         jarg=ncloj+1-mdon                                               6d11s19
         jjo=iptrbit(2,nclopj,mysym)+ibasis(3,jf)-1                     4d13s21
         jjc=iptrbit(1,nclopj,mysym)+ibasis(2,jf)-1                     4d13s21
         jp=jps+1
         call gandc4(ibc(jjc),ibc(jjo),ibc(iic),ibc(iio),nopenj,nopen,  11d16s20
     $        norb,nnot,nab4,bc,ibc)                                    11d14s22
         if(ldebug)write(6,*)('nnot from gandc4 '),nnot
         if(nnot.gt.0)then                                               11d4s19
          if(mod(loopit,mynprocg).eq.mynowprog)then                      12d9s19
           loopit=mynowprog                                             9d10s24
           ibcsav=ibcoff                                                  12d9s19
           nn=ncsf(jarg)*ncsf(iarg)                                     12d9s19
           if(nnot.eq.1)then                                            4d13s21
            imat=ibcoff                                                 4d13s21
            ibcoff=imat+ncsf(iarg)*ncsf(iarg)                           4d13s21
            call enough('hccsfd.  1',bc,ibc)
            do i1=1,norb-1                                              8d2s22
             if(btest(ibc(iio),i1))then                                 8d2s22
              do i2=i1+1,norb                                           8d2s22
               if(btest(ibc(iio),i2))then                               8d2s22
                nab1(2)=i1                                              8d2s22
                nab1(1)=i2                                              8d2s22
                itesta=ibc(iic)                                             11d13s20
                itestb=ibc(iio)                                             11d13s20
                nopenk=nopen-2                                              11d13s20
                karg=iarg+1                                                 11d13s20
                nab2(1)=nab1(2)                                           12d6s20
                nab2(2)=nab1(1)                                           12d6s20
                jsa=ism(nab1(1))                                              11d13s20
                jga=irel(nab1(1))-1                                           11d13s20
                jsb=ism(nab1(2))                                              11d13s20
                jgb=irel(nab1(2))-1                                           11d13s20
                jsc=ism(nab2(1))                                              11d13s20
                jgc=irel(nab2(1))-1                                           11d13s20
                jsd=ism(nab2(2))                                              11d13s20
                jgd=irel(nab2(2))-1                                           11d13s20
                nqq=karg+mdon-1
                if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                 itestb=ibclr(itestb,nab1(1))                           8d2s22
                 itestb=ibclr(itestb,nab1(2))                           8d2s22
                 itesta=ibset(itesta,nab1(2))                           8d2s22
                 call gandc(ibc(jjc),ibc(jjo),itesta,itestb,nopenj,       12d6s20
     $               nopenk,jarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,   12d6s20
     $               iwpb1,iwpk1,ncsfmid1,bc,ibc)                       11d14s22
                 call gandc(itesta,itestb,ibc(iic),ibc(iio),nopenk,     8d2s22
     $                nopen,karg,iarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,   8d2s22
     $                iwpb2,iwpk2,ncsfmid2,bc,ibc)                      11d14s22
                 call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),ncsfmid1,   4d13s21
     $               iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod,bc,ibc)     11d10s22
                 do i=0,ncsf(iarg)*ncsf(iarg)-1                           4d13s21
                  bc(imat+i)=bc(iprod+i)                                  4d13s21
                 end do                                                   4d13s21
                 ibcoff=iprod                                                11d13s20
                else                                                        11d13s20
                 do i=0,ncsf(iarg)*ncsf(iarg)-1                           4d13s21
                  bc(imat+i)=0d0                                          4d13s21
                 end do                                                   4d13s21
                end if                                                      11d13s20
                do i=0,ncsf(iarg)-1                                       4d13s21
                 ii=imat+i*(ncsf(iarg)+1)                                 4d13s21
                 bc(ii)=bc(ii)-1d0                                        4d13s21
                end do                                                    4d13s21
                i2eu=inv(1,jsa,jsb,jsc)                                      12d13s19
                icase=inv(2,jsa,jsb,jsc)                                      12d13s19
                pz=1d0                                                      12d22s19
                if(jsa.eq.jsb)then                                          12d13s19
                 ix=max(jga,jgb)                                            12d13s19
                 in=min(jga,jgb)                                            12d13s19
                 irow=((ix*(ix+1))/2)+in                                    12d13s19
                 nrow=(irefo(jsa)*(irefo(jsa)+1))/2                         12d13s19
                 ix=max(jgc,jgd)                                            12d13s19
                 in=min(jgc,jgd)                                            12d13s19
                 icol=((ix*(ix+1))/2)+in                                    12d13s19
                else                                                        12d13s19
                 nrow=irefo(jsa)*irefo(jsb)                                 12d13s19
                 if(icase.ge.3.and.istype.ne.0)pz=-1d0                      12d27s19
                 if(icase.eq.1)then                                         12d13s19
                  irow=jga+irefo(jsa)*jgb                                   12d13s19
                  icol=jgc+irefo(jsc)*jgd                                   12d13s19
                 else if(icase.eq.2)then                                    12d13s19
                  irow=jgb+irefo(jsb)*jga                                   12d13s19
                  icol=jgd+irefo(jsd)*jgc                                   12d13s19
                 else if(icase.eq.3)then                                    12d13s19
                  irow=jga+irefo(jsa)*jgb                                   12d13s19
                  icol=jgd+irefo(jsd)*jgc                                   12d13s19
                 else                                                       12d13s19
                  irow=jgb+irefo(jsb)*jga                                   12d13s19
                  icol=jgc+irefo(jsc)*jgd                                   12d13s19
                 end if                                                     12d13s19
                end if                                                      12d13s19
                jden2=iden2(i2eu)+irow+nrow*icol                            12d17s19
                ww=0d0                                                      12d17s19
                do ir=1,nrootz                                              12d13s19
                 wwx=wgt(ir)*fact                                           12d27s19
                 do i=0,ncsf(iarg)-1                                        12d13s19
                  wfact=vecx(ip+i,ir)*wwx                                    12d13s19
                  jout=imat+ncsf(jarg)*i                                    12d13s19
                  do j=0,ncsf(jarg)-1                                       12d13s19
                   ww=ww+wfact*vecx(jp+j,ir)*bc(jout+j)                     12d17s19
                  end do                                                    12d13s19
                 end do                                                     12d13s19
                end do                                                      12d13s19
                orig=bc(jden2)
                bc(jden2)=bc(jden2)+ww*pz                                   12d22s19
               end if                                                   8d2s22
              end do                                                       11d13s20
             end if                                                     8d2s22
            end do                                                        11d13s20
            ibcoff=imat                                                 4d13s21
           else if(nnot.eq.2)then                                       4d13s21
            call gandc(ibc(jjc),ibc(jjo),ibc(iic),ibc(iio),nopenj,      11d16s20
     $            nopen,jarg,iarg,ncsf,norb,ixw1,ixw2,nnot1,nab1,iwpb1, 11d16s20
     $            iwpk1,ncsfmid1,bc,ibc)                                11d14s22
            call gandc(ibc(iic),ibc(iio),ibc(jjc),ibc(jjo),nopen,       11d16s20
     $            nopenj,iarg,jarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,iwpb2, 11d16s20
     $            iwpk2,ncsfmid2,bc,ibc)                                11d14s22
            iprod=ibcoff                                                  11d13s20
            ibcoff=iprod+ncsf(iarg)*ncsf(jarg)                          4d13s21
            call enough('hccsfd.  2',bc,ibc)
            call prodn(iwpb1,iwpk1,ncsf(jarg),ncsf(iarg),ncsfmid1,      4d13s21
     $           bc(iprod),bc,ibc,1d0,0d0)                              2d13s23
            do i=1,norb                                                   11d13s20
             itest(i,1)=0                                                   11d13s20
             itest(i,2)=0                                                    11d13s20
             if(btest(ibc(iic),i))itest(i,1)=2                          8d2s22
             if(btest(ibc(iio),i))itest(i,1)=1                          8d2s22
             if(btest(ibc(jjc),i))itest(i,2)=2                          8d2s22
             if(btest(ibc(jjo),i))itest(i,2)=1                          8d2s22
            end do                                                        11d13s20
            nok=0                                                         11d13s20
            do i=1,norb                                                   11d13s20
             ixn=min(itest(i,1),itest(i,2))
             if(ixn.gt.0)then                                             11d13s20
              nok=nok+1                                                   11d13s20
              itest(nok,1)=ixn                                            11d13s20
              itest(nok,2)=i                                              11d13s20
             end if                                                       11d13s20
            end do                                                        11d13s20
            is=ism(nab4(1,1))                                             11d13s20
            jga=irel(nab4(2,1))-1                                       11d27s20
            jgb=irel(nab4(1,1))-1                                       11d27s20
            jden1=iden1(is)+jga+irefo(is)*jgb                           12d13s19
            ww=0d0                                                      12d17s19
            do ir=1,nrootz                                              12d17s19
             do i=0,ncsf(iarg)-1                                        12d17s19
              jout=iprod+ncsf(jarg)*i                                    12d17s19
              wfact=fact*wgt(ir)*vecx(ip+i,ir)                          12d17s19
              do j=0,ncsf(jarg)-1                                       12d17s19
               ww=ww+wfact*vecx(jp+j,ir)*bc(jout+j)
              end do                                                    12d17s19
             end do                                                     12d17s19
            end do                                                      12d17s19
            bc(jden1)=bc(jden1)+ww                                      12d17s19
            ix=max(jga,jgb)                                             12d17s19
            in=min(jga,jgb)                                             12d17s19
            icolj=((ix*(ix+1))/2)+in                                    12d17s19
            do i=1,nok                                                    11d13s20
             js=ism(itest(i,2))                                           11d13s20
             jg=irel(itest(i,2))-1                                        11d13s20
             i2eu=inv(1,js,js,is)                                       4d13s21
             irow=((jg*(jg+1))/2)+jg                                    4d13s21
             nrow=(irefo(js)*(irefo(js)+1))/2                           4d13s21
             jden2j=iden2(i2eu)+irow+nrow*icolj                         12d17s19
             if(itest(i,1).eq.2)then                                      11d13s20
              bc(jden2j)=bc(jden2j)+ww*2d0                               12d17s19
              i2eu=inv(1,js,is,is)                                       4d13s21
              icase=inv(2,js,is,is)                                      4d13s21
              if(js.eq.is)then                                           4d13s21
               ix=max(jg,jgb)                                            4d13s21
               in=min(jg,jgb)                                            4d13s21
               irow=((ix*(ix+1))/2)+in                                   4d13s21
               nrow=(irefo(js)*(irefo(js)+1))/2                          4d13s21
               ix=max(jga,jg)                                            4d13s21
               in=min(jga,jg)                                            4d13s21
               icol=((ix*(ix+1))/2)+in                                   12d13s19
              else                                                       12d13s19
               nrow=irefo(js)*irefo(is)                                  4d13s21
               if(icase.eq.1)then                                        12d13s19
                irow=jg+irefo(js)*jgb                                    4d13s21
                icol=jga+irefo(is)*jg                                    4d13s21
               else if(icase.eq.2)then                                    12d13s19
                irow=jgb+irefo(is)*jg                                    4d13s21
                icol=jg+irefo(js)*jga                                    4d13s21
               else if(icase.eq.3)then                                    12d13s19
                irow=jg+irefo(js)*jgb                                    4d13s21
                icol=jg+irefo(js)*jga                                    4d13s21
               else                                                      12d13s19
                irow=jgb+irefo(is)*jg                                    4d13s21
                icol=jga+irefo(is)*jg                                    4d13s21
               end if                                                    12d13s19
              end if                                                     12d13s19
              jden2k=iden2(i2eu)+irow+nrow*icol                           12d13s19
              bc(jden2k)=bc(jden2k)-ww                                   12d27s19
             else                                                         11d13s20
              bc(jden2j)=bc(jden2j)+ww                                  4d13s21
             end if                                                       11d13s20
            end do                                                        11d13s20
            ibcoff=iprod                                                4d13s21
            do i=1,nok                                                    11d13s20
             if(itest(i,1).eq.1)then                                      11d13s20
              itesta=ibc(jjc)                                              11d13s20
              itestb=ibc(jjo)                                              11d13s20
              nopenk=nopenj                                                11d13s20
c
c     anihilate common
c
              if(btest(itesta,itest(i,2)))then                             11d13s20
               itesta=ibclr(itesta,itest(i,2))                             11d13s20
               itestb=ibset(itestb,itest(i,2))                             11d13s20
               karg=jarg-1                                                11d13s20
               nopenk=nopenk+1                                             11d13s20
              else                                                         11d13s20
               itestb=ibclr(itestb,itest(i,2))                             11d13s20
               karg=jarg                                                  11d13s20
               nopenk=nopenk-1                                             11d13s20
              end if                                                       11d13s20
c
c     create ket
c
              if(btest(itestb,nab4(2,1)))then                           11d27s20
               itesta=ibset(itesta,nab4(2,1))                           11d27s20
               itestb=ibclr(itestb,nab4(2,1))                           11d27s20
               karg=karg+1                                                11d13s20
               nopenk=nopenk-1                                             11d13s20
              else                                                         11d13s20
               itestb=ibset(itestb,nab4(2,1))                           11d27s20
               nopenk=nopenk+1                                             11d13s20
              end if                                                       11d13s20
              nqq=karg+mdon-1
              if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
               call gandc(ibc(jjc),ibc(jjo),itesta,itestb,nopenj,       11d16s20
     $              nopenk,jarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,    11d16s20
     $              iwpb1,iwpk1,ncsfmid1,bc,ibc)                        11d14s22
               call gandc(itesta,itestb,ibc(iic),ibc(iio),nopenk,nopen,      11d13s20
     $         karg,iarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,    11d13s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
               if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),ncsfmid1,   4d13s21
     $               iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,imat,bc,ibc)      11d10s22
                jsa=ism(nab1(1))                                              11d13s20
                jga=irel(nab1(1))-1                                           11d13s20
                jsb=ism(nab1(2))                                              11d13s20
                jgb=irel(nab1(2))-1                                           11d13s20
                jsc=ism(nab2(1))                                              11d13s20
                jgc=irel(nab2(1))-1                                           11d13s20
                jsd=ism(nab2(2))                                              11d13s20
                jgd=irel(nab2(2))-1                                           11d13s20
                i2eu=inv(1,jsa,jsb,jsc)                                      12d13s19
                icase=inv(2,jsa,jsb,jsc)                                      12d13s19
                pz=1d0                                                      12d22s19
                if(jsa.eq.jsb)then                                          12d13s19
                 ix=max(jga,jgb)                                            12d13s19
                 in=min(jga,jgb)                                            12d13s19
                 irow=((ix*(ix+1))/2)+in                                    12d13s19
                 nrow=(irefo(jsa)*(irefo(jsa)+1))/2                         12d13s19
                 ix=max(jgc,jgd)                                            12d13s19
                 in=min(jgc,jgd)                                            12d13s19
                 icol=((ix*(ix+1))/2)+in                                    12d13s19
                else                                                        12d13s19
                 nrow=irefo(jsa)*irefo(jsb)                                 12d13s19
                 if(icase.ge.3.and.istype.ne.0)pz=-1d0                      12d27s19
                 if(icase.eq.1)then                                         12d13s19
                  irow=jga+irefo(jsa)*jgb                                   12d13s19
                  icol=jgc+irefo(jsc)*jgd                                   12d13s19
                 else if(icase.eq.2)then                                    12d13s19
                  irow=jgb+irefo(jsb)*jga                                   12d13s19
                  icol=jgd+irefo(jsd)*jgc                                   12d13s19
                 else if(icase.eq.3)then                                    12d13s19
                  irow=jga+irefo(jsa)*jgb                                   12d13s19
                  icol=jgd+irefo(jsd)*jgc                                   12d13s19
                 else                                                       12d13s19
                  irow=jgb+irefo(jsb)*jga                                   12d13s19
                  icol=jgc+irefo(jsc)*jgd                                   12d13s19
                 end if                                                     12d13s19
                end if                                                      12d13s19
                jden2=iden2(i2eu)+irow+nrow*icol                            12d17s19
                ww=0d0                                                      12d17s19
                do ir=1,nrootz                                              12d13s19
                 wwx=wgt(ir)*fact                                           12d27s19
                 do ii=0,ncsf(iarg)-1                                        12d13s19
                  wfact=vecx(ip+ii,ir)*wwx                                    12d13s19
                  jout=imat+ncsf(jarg)*ii                                    12d13s19
                  do j=0,ncsf(jarg)-1                                       12d13s19
                   ww=ww+wfact*vecx(jp+j,ir)*bc(jout+j)                     12d17s19
                  end do                                                    12d13s19
                 end do                                                     12d13s19
                end do                                                      12d13s19
                bc(jden2)=bc(jden2)+ww*pz                                   12d22s19
                ibcoff=imat                                             4d13s21
               else if(max(nnot1,nnot2).lt.2)then                        11d13s20
                write(6,*)('c:expecting 2s, but got otherwise')
                call dcbit(itesta,norb,'itesta')
                call dcbit(itestb,norb,'itestb')
                stop
               end if
              end if                                                    11d13s20
             end if                                                       11d13s20
            end do                                                        11d13s20
           else                                                         4d13s21
c     j is bra and i is ket
            if(nnot.eq.3)then                                           1d22s21
             ipssx=1                                                    1d22s21
            else                                                        1d22s21
             ipssx=3                                                    1d22s21
            end if                                                      1d22s21
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
             itesta=ibc(jjc)                                               11d13s20
             itestb=ibc(jjo)                                               11d13s20
             if(btest(itesta,nab4(1,iu1)))then                          1d22s21
              itesta=ibclr(itesta,nab4(1,iu1))                          1d22s21
              itestb=ibset(itestb,nab4(1,iu1))                          1d22s21
              nopenk=nopenj+1                                           1d22s21
              karg=jarg-1                                               1d22s21
             else if(btest(itestb,nab4(1,iu1)))then                     1d22s21
              itestb=ibclr(itestb,nab4(1,iu1))                          1d22s21
              nopenk=nopenj-1                                              11d13s20
              karg=jarg                                                  11d13s20
             else                                                          11d13s20
              write(6,*)('bit not set for nab4(1,1) = '),nab4(1,iu1)      11d27s20
              stop 'nab4(1,1)'                                          11d27s20
             end if                                                        11d13s20
             if(btest(itestb,nab4(2,iu2)))then                          1d22s21
              itesta=ibset(itesta,nab4(2,iu2))                          1d22s21
              itestb=ibclr(itestb,nab4(2,iu2))                          1d22s21
              nopenk=nopenk-1                                              11d13s20
              karg=karg+1                                                  11d13s20
             else if(btest(itesta,nab4(2,iu2)))then                     1d22s21
              write(6,*)('already double in nab4(2,1) = '),nab4(2,iu2)  1d22s21
              stop 'nab4(2,1)'                                          11d27s20
             else                                                          11d13s20
              itestb=ibset(itestb,nab4(2,iu2))                          1d22s21
              nopenk=nopenk+1                                              11d13s20
             end if                                                        11d13s20
             nqq=karg+mdon-1                                            1d26s21
             if(nqq.ge.mdon.and.nqq.le.mdoo)then                        1d26s21
              call gandc(ibc(jjc),ibc(jjo),itesta,itestb,nopenj,nopenk,     11d13s20
     $         jarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,    11d13s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
              call gandc(itesta,itestb,ibc(iic),ibc(iio),nopenk,nopen,      11d13s20
     $         karg,iarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,    11d13s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
              if(nnot1.eq.2.and.nnot2.eq.2)then                         1d22s21
               call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),ncsfmid1,    4d13s21
     $               iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,imat,bc,ibc)      11d10s22
               jsa=ism(nab1(1))                                              11d13s20
               jga=irel(nab1(1))-1                                           11d13s20
               jsb=ism(nab1(2))                                              11d13s20
               jgb=irel(nab1(2))-1                                           11d13s20
               jsc=ism(nab2(1))                                              11d13s20
               jgc=irel(nab2(1))-1                                           11d13s20
               jsd=ism(nab2(2))                                              11d13s20
               jgd=irel(nab2(2))-1                                           11d13s20
               i2eu=inv(1,jsa,jsb,jsc)                                      12d13s19
               icase=inv(2,jsa,jsb,jsc)                                      12d13s19
               pz=1d0                                                      12d22s19
               if(jsa.eq.jsb)then                                          12d13s19
                ix=max(jga,jgb)                                            12d13s19
                in=min(jga,jgb)                                            12d13s19
                irow=((ix*(ix+1))/2)+in                                    12d13s19
                nrow=(irefo(jsa)*(irefo(jsa)+1))/2                         12d13s19
                ix=max(jgc,jgd)                                            12d13s19
                in=min(jgc,jgd)                                            12d13s19
                icol=((ix*(ix+1))/2)+in                                    12d13s19
               else                                                        12d13s19
                nrow=irefo(jsa)*irefo(jsb)                                 12d13s19
                if(icase.ge.3.and.istype.ne.0)pz=-1d0                      12d27s19
                if(icase.eq.1)then                                         12d13s19
                 irow=jga+irefo(jsa)*jgb                                   12d13s19
                 icol=jgc+irefo(jsc)*jgd                                   12d13s19
                else if(icase.eq.2)then                                    12d13s19
                 irow=jgb+irefo(jsb)*jga                                   12d13s19
                 icol=jgd+irefo(jsd)*jgc                                   12d13s19
                else if(icase.eq.3)then                                    12d13s19
                 irow=jga+irefo(jsa)*jgb                                   12d13s19
                 icol=jgd+irefo(jsd)*jgc                                   12d13s19
                else                                                       12d13s19
                 irow=jgb+irefo(jsb)*jga                                   12d13s19
                 icol=jgc+irefo(jsc)*jgd                                   12d13s19
                end if                                                     12d13s19
               end if                                                      12d13s19
               jden2=iden2(i2eu)+irow+nrow*icol                            12d17s19
               ww=0d0                                                      12d17s19
               xfact=fact                                               4d13s21
               if(nab1(1).eq.nab2(1).and.nab1(2).eq.nab2(2))              11d16s20
     $              xfact=fact*0.5d0                                    4d13s21
               if(istype.ne.0)pz=-pz                                    4d13s21
               do ir=1,nrootz                                              12d13s19
                wwx=wgt(ir)*xfact                                       4d13s21
                do i=0,ncsf(iarg)-1                                        12d13s19
                 wfact=vecx(ip+i,ir)*wwx                                    12d13s19
                 jout=imat+ncsf(jarg)*i                                    12d13s19
                 do j=0,ncsf(jarg)-1                                       12d13s19
                  ww=ww+wfact*vecx(jp+j,ir)*bc(jout+j)                     12d17s19
                 end do                                                    12d13s19
                end do                                                     12d13s19
               end do                                                      12d13s19
               bc(jden2)=bc(jden2)+ww*pz                                   12d22s19
               ibcoff=imat                                                  11d13s20
               if(ipss.eq.2)go to 22                                    1d22s21
              end if                                                    1d26s21
             end if                                                     1d22s21
            end do                                                      1d22s21
   22       continue                                                    1d22s21
           end if                                                       4d13s21
           ibcoff=ibcsav                                                 12d9s19
          end if                                                         12d9s19
          loopit=loopit+1                                                12d9s19
         end if                                                          12d9s19
        end if                                                          12d9s19
        jps=jps+ncsf(jarg)
       end do                                                           6d11s19
       ips=ips+ncsf(iarg)
      end do
      nwds=ncsft*nrootz                                                 7d11s19
      if(ldebug)write(6,*)('all done in hccsfd')
      return
      end                                                               7d11s19
