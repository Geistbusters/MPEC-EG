c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine cont2ecsfb(vc,nct2,nrtot,ncol,mdon,ncsf,               11d19s20
     $     icsfpd,nec,thrs,nkkk,iscp,ivcc,nct2c,ncsf2,                  11d19s20
     $     lprint,ilcol,ihcol,iff2,nff2,mdoo,isb,nsymb,nfdat,           11d19s20
     $     needcv,nff22,bc,ibc,nder,idvmt,ntackon)                      8d7s24
      implicit real*8 (a-h,o-z)
      real*16 fl                                                        5d27s19
      integer*8 iff2(*),ipack8,idvmt(2),i18,il8,ih8                     9d1s23
      logical ldebug,lprint                                             12d12s19
      equivalence(ipack8,ipack4)                                        11d20s20
      character*6 line6                                                 9d1s23
      dimension vc(nct2,ncol*nrtot),ncsf(*),icsfpd(*),nff22(mdoo+1,2),  12d18s20
     $     itcp(3),nkkk(3,2),ivcc(4),nct2c(4),kkk(4),ncsf2(4,*),        11d7s19
     $     kkkt(4),n2ch(4),nff2(mdoo+1,nsymb,*),nfdat(5,4),noffs(4,2),  11d20s20
     $     jvcc(4),ipack4(2),nll(4),mpre(4),mpred(4)                    8d30s23
      COMMON/FACT16/FL(922),NCALL                                       5d27s19
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      include "common.store"
      include "common.basis"                                            8d1s19
      include "common.input"
      include "common.mrci"                                             8d1s19
      include "common.print"                                            1d5s20
      if(iprtr(9).eq.0)then                                             1d5s20
       ldebug=.false.                                                    12d12s19
      else                                                              1d5s20
       ldebug=.true.                                                    1d5s20
      end if                                                            1d5s20
      if(ldebug)write(6,*)('in cont2ecsfb for symmetry '),isb           11d20s20
      do i=1,3                                                          8d2s19
       nkkk(i,1)=0                                                      10d20s20
       nkkk(i,2)=0                                                      8d2s19
      end do                                                            8d2s19
      do l=1,4                                                          11d19s20
       do i=1,5                                                         11d19s20
        nfdat(i,l)=0                                                    11d19s20
       end do                                                           11d19s20
      end do                                                            11d19s20
      ndim=ncol*nrtot
      keep=0                                                            8d1s19
      thrsk=smallest
      if(ldebug)then                                                    12d12s19
       write(6,*)('kill threshold: '),thrsk
       write(6,*)('keep threshold: '),thrs
       write(6,*)('nder: '),nder                                        8d30s23
      end if                                                            12d12s19
      if(ldebug)write(6,*)('normalizing ivcc keepers ')
      ndimh=(ihcol+1-ilcol)*nrtot                                       10d17s20
      ndimh0=(ihcol+1-ilcol)                                            8d30s23
      nderp=nder+1                                                      8d30s23
      do i=1,4                                                          8d15s19
       kkk(i)=0                                                         8d15s19
       if(nct2c(i).gt.0)then                                            10d17s20
        iad=ivcc(i)                                                      8d15s19
        jad=ivcc(i)                                                      8d15s19
        if(nder.gt.0)then                                               8d30s23
c
c     in cont2ecsfa, we pretended ders were just extra roots, so        8d30s23
c     ivcc is csf,der+1,nroot,col ...                                      8d30s23
c     for cont2ecsfb, we want ivcc to be csf,root,col,der+1
c
         nnd=ndimh*nderp                                                8d30s23
         itmp=ibcoff                                                    8d30s23
         ibcoff=itmp+nnd                                                8d30s23
         call enough('cont2ecsfb.tmp',bc,ibc)                           8d30s23
         do k=0,nct2c(i)-1                                              8d30s23
          do icl=0,ndimh0-1                                             8d30s23
           do idr=0,nder                                                8d30s23
            do irt=0,nrtot-1                                              8d30s23
             iad1=ivcc(i)+k+nct2c(i)*(idr+nderp*(irt+nrtot*icl))        8d30s23
             iad2=itmp+irt+nrtot*(icl+ndimh0*idr)                       8d30s23
             bc(iad2)=bc(iad1)                                          8d30s23
            end do                                                      8d30s23
           end do                                                       8d30s23
          end do                                                        8d30s23
          do j=0,nnd-1                                                  8d30s23
           iad1=ivcc(i)+k+nct2c(i)*j                                    8d30s23
           bc(iad1)=bc(itmp+j)                                          8d30s23
          end do                                                        8d30s23
         end do                                                         8d30s23
         ibcoff=itmp                                                    8d30s23
        end if                                                          8d30s23
        if(ldebug)call prntm2(bc(ivcc(i)),nct2c(i),ndimh*nderp,         8d30s23
     $       nct2c(i))                                                  8d30s23
        do j=1,ndimh                                                      8d15s19
         sz=0d0                                                          8d15s19
         vx=0d0                                                         12d23s20
         do k=0,nct2c(i)-1                                               8d15s19
          xvx=abs(bc(iad+k))                                            12d23s20
          vx=max(vx,xvx)                                                12d23s20
          sz=sz+bc(iad+k)**2                                             8d15s19
         end do                                                          8d15s19
         if(sz.gt.thrsk)then                                             8d15s19
          if(lcrenorm.ne.0)then                                          1d22s20
           szx=1d0/sqrt(sz)                                               8d15s19
           if(nder.gt.0)then                                            8d30s23
            do idr=1,nder                                               8d30s23
             iaddd=ivcc(i)+nct2c(i)*(j-1+ndimh*idr)                     8d30s23
             jaddd=ivcc(i)+nct2c(i)*(kkk(i)+ndimh*idr)                  8d30s23
             dsz=0d0                                                    8d30s23
             do k=0,nct2c(i)-1                                          8d30s23
              dsz=dsz+bc(iad+k)*bc(iaddd+k)                             8d30s23
             end do                                                     8d30s23
             dsz=dsz*2d0                                                8d30s23
             dszx=-dsz*0.5d0*(szx**3)                                   8d30s23
c
c     v=u*szx, v'=u'*szx+u*szx'
c
             do k=0,nct2c(i)-1                                          8d30s23
              bc(jaddd+k)=bc(iaddd+k)*szx+bc(iad+k)*dszx                8d30s23
             end do                                                     8d30s23
            end do                                                      8d30s23
           end if                                                       8d30s23
          else                                                           1d22s20
           do idr=1,nder                                                8d30s23
            iaddd=ivcc(i)+nct2c(i)*(j-1+ndimh*idr)                      8d30s23
            jaddd=ivcc(i)+nct2c(i)*(kkk(i)+ndimh*idr)                   8d30s23
            do k=0,nct2c(i)-1                                           8d30s23
             bc(jaddd+k)=bc(iaddd+k)                                    8d30s23
            end do                                                      8d30s23
           end do                                                       8d30s23
           szx=1d0
          end if                                                         1d22s20
          do k=0,nct2c(i)-1                                              8d15s19
           bc(jad+k)=bc(iad+k)*szx                                       8d15s19
          end do                                                         8d15s19
          jad=jad+nct2c(i)                                               8d15s19
          kkk(i)=kkk(i)+1                                                8d15s19
          if(ldebug)write(6,*)i,kkk(i),j,sz
         end if                                                          8d15s19
         iad=iad+nct2c(i)                                                8d15s19
        end do                                                           8d15s19
        if(nder.gt.0)then                                               8d30s23
         do idr=1,nder                                                  8d30s23
          do j=0,kkk(i)-1                                               8d30s23
           iaddo=ivcc(i)+nct2c(i)*(j+ndimh*idr)                         8d30s23
           iaddn=ivcc(i)+nct2c(i)*(j+kkk(i)*idr)                        8d30s23
           do k=0,nct2c(i)-1                                            8d30s23
            bc(iaddn+k)=bc(iaddo+k)                                     8d30s23
           end do                                                       8d30s23
          end do                                                        8d30s23
         end do                                                         8d30s23
        end if                                                          8d30s23
       end if                                                           10d17s20
      end do                                                            8d15s19
      if(ldebug)then                                                    1d25s21
       write(6,*)('number of keepers = '),kkk
       do l=1,4
        if(nct2c(l).gt.0)then
         write(6,*)('for l = '),l
         call prntm2(bc(ivcc(l)),nct2c(l),kkk(l)*(1+nder),nct2c(l))     8d30s23
        end if
       end do
      end if                                                            1d25s21
c
      nfcn=0                                                            11d17s20
      do nclo=mdon,mdoo                                                 11d17s20
       nclop=nclo+1                                                     11d17s20
       nfcn=nfcn+nff2(nclop,isb,1)                                      11d17s20
      end do                                                            11d17s20
      if(ldebug)write(6,*)('no. fcns of this sym '),nfcn                1d25s21
      nct2cn=0                                                          11d17s20
      nls=0                                                             11d19s20
      do l=1,4                                                          11d17s20
       if(nct2c(l).gt.0)then                                            11d20s20
        nct2cn=nct2cn+nfcn                                              11d20s20
        nls=nls+1                                                       11d19s20
       end if                                                           11d20s20
      end do                                                            11d17s20
      ifkeep=ibcoff                                                     11d17s20
      ibcoff=ifkeep+nct2cn*mynprocg                                     11d17s20
      call enough('cont2ecsfb.  1',bc,ibc)
      do i=ifkeep,ibcoff-1                                              11d17s20
       bc(i)=0d0                                                        11d17s20
      end do                                                            11d17s20
      jfkeep=ifkeep+nct2cn*mynowprog                                    11d17s20
      do l=1,4
       if(kkk(l).gt.0)then
        ioffv=0                                                         11d17s20
        ikeep=ibcoff
        ibcoff=ikeep+kkk(l)
        call enough('cont2ecsfb.  2',bc,ibc)
        do k=0,kkk(l)-1
         ibc(ikeep+k)=0
        end do
        do nclo=mdon,mdoo
         nclop=nclo+1                                                   11d17s20
         iarg=nclop-mdon                                                11d17s20
         do if=1,nff2(nclop,isb,1)
          iff=nff2(nclop,isb,2)+if-1
          do j=0,kkk(l)-1
           kvcc=ivcc(l)+ioffv+nct2c(l)*j                                11d17s20
           if(ncsf2(l,iarg).gt.0)then                                   1d26s21
            sz=0d0                                                        11d17s20
            do k=0,ncsf2(l,iarg)-1                                       11d17s20
             sz=sz+bc(kvcc+k)**2                                         11d17s20
            end do                                                       11d17s20
            sz=sqrt(sz/dfloat(ncsf2(l,iarg)))                            11d17s20
            if(sz.gt.1d-10)then
             if(ldebug)call dcbit(iff2(iff),32+norb,'orb')               1d25s21
             bc(jfkeep)=bc(jfkeep)+1d0                                         11d17s20
             ibc(ikeep+j)=ibc(ikeep+j)+ncsf2(l,iarg)
            end if
           end if                                                       1d26s21
          end do
          ioffv=ioffv+ncsf2(l,iarg)
          jfkeep=jfkeep+1                                               11d17s20
         end do
        end do
        if(ldebug)write(6,*)(ibc(ikeep+j),j=0,kkk(l)-1)                 1d25s21
        ibcoff=ikeep                                                    11d17s20
       else                                                             4d28s21
        if(nct2c(l).gt.0)jfkeep=jfkeep+nfcn                             4d28s21
       end if
      end do
      call dws_gsumf(bc(ifkeep),nct2cn*mynprocg)                        11d17s20
      ikkk=ibcoff                                                       10d17s20
      ibcoff=ikkk+4*mynprocg                                            10d17s20
      call enough('cont2ecsfb.  4',bc,ibc)
      do i=ikkk,ibcoff-1                                                10d17s20
       bc(i)=0d0                                                        10d17s20
      end do                                                            10d17s20
      jkkk=ikkk-1+4*mynowprog                                           10d17s20
      do l=1,4                                                          10d17s20
       bc(jkkk+l)=dfloat(kkk(l))                                        10d17s20
      end do                                                            10d17s20
      call dws_gsumf(bc(ikkk),4*mynprocg)                               10d17s20
      do l=1,4                                                          10d17s20
       lm=l-1                                                           10d17s20
       sum=0d0                                                          10d17s20
       sumpre=0d0                                                       11d17s20
       do ip=0,mynprocg-1                                               10d17s20
        jkkk=ikkk+lm+4*ip                                               10d17s20
        sum=sum+bc(jkkk)                                                10d17s20
        if(ip.lt.mynowprog)sumpre=sumpre+bc(jkkk)                       11d17s20
       end do                                                           10d17s20
       ibc(ikkk+lm)=nint(sumpre)                                        11d17s20
       kkkt(l)=nint(sum)                                                10d17s20
      end do                                                            10d17s20
      if(ldebug)then                                                    1d25s21
       write(6,*)('number of keepers across all procs '),(kkkt(l),l=1,4) 10d17s20
       write(6,*)('number before me: '),(ibc(ikkk+lm),lm=0,3)            11d17s20
      end if                                                            1d25s21
      jfkeep=ifkeep                                                     11d17s20
      needcv=0                                                          11d17s20
      needcvd=0                                                         8d30s23
      ioff4=1
      do nclo=mdon,mdoo                                                 11d20s20
       nclop=nclo+1                                                     11d20s20
       iarg=nclop-mdon                                                  11d20s20
       do if=0,nff2(nclop,isb,1)-1                                      11d20s20
        kfkeep=jfkeep                                                   11d20s20
        sum=0d0                                                         11d20s20
        nvec=0                                                          11d20s20
        ioff4=ioff4+ncsf2(4,iarg)
        do ip=0,mynprocg-1                                              11d20s20
         do l=1,4                                                       11d20s20
          if(nct2c(l).gt.0)then                                         11d20s20
           nvec=nvec+nint(bc(kfkeep))*ncsf2(l,iarg)                     11d28s20
           sum=sum+bc(kfkeep)                                             11d20s20
           kfkeep=kfkeep+nfcn                                             11d20s20
          end if                                                        11d20s20
         end do                                                         11d20s20
        end do                                                          11d20s20
        if(sum.gt.1d-8)then
         nhere=nint(sum)                                                11d20s20
c     got
c     1: bit codes for closed&open
c     2: nhere
c     nhere: ks
c     nhere: ls
c     want
c     1: bit codes
c     2: nvec
c     3-7: nl
c     8-12: npre
         nspace=2+8+nvec+nhere                                          3d18s21
         needcv=needcv+nspace                                           11d20s20
         needcvd=needcvd+nvec                                           8d30s23
        end if                                                          11d20s20
        jfkeep=jfkeep+1                                                 11d20s20
       end do                                                           11d20s20
      end do                                                            11d20s20
      ivcv=ibcoff                                                       11d17s20
      ibcoff=ivcv+needcv                                                11d17s20
      call enough('cont2ecsfb.  5',bc,ibc)
      ndorth=0                                                          9d5s23
      do l=1,4                                                          9d5s23
       ndorth=ndorth+kkkt(l)**2                                         9d5s23
      end do                                                            9d5s23
      iorth=ibcoff                                                      9d5s23
      ibcoff=iorth+ndorth                                               9d5s23
      ndorth0=ndorth                                                    10d17s23
      ndorth=ndorth*nder                                                9d5s23
      idorth=ibcoff                                                     9d5s23
      ivcvd=idorth+ndorth                                               9d5s23
      ibcoff=ivcvd+needcvd*nder                                         8d30s23
      call enough('cont2ecsfb.vcvd',bc,ibc)                             8d30s23
      do i=ivcv,ibcoff-1                                                11d17s20
       bc(i)=0d0                                                        11d17s20
      end do                                                            11d17s20
      jvcv=ivcv                                                         11d17s20
      jvcvd=ivcvd                                                       8d30s23
      jfkeep=ifkeep                                                     11d17s20
      igoal=ivcvd
      do l=1,4                                                          11d20s20
       jvcc(l)=ivcc(l)                                                  11d20s20
      end do                                                            11d20s20
      nl=0                                                              11d20s20
      do nclo=mdon,mdoo                                                 11d20s20
       nclop=nclo+1                                                     11d20s20
       iarg=nclop-mdon                                                  11d20s20
       nff2(nclop,isb,7)=0                                              11d20s20
       do if=0,nff2(nclop,isb,1)-1                                      11d20s20
        kfkeep=jfkeep                                                   11d20s20
        sum=0d0                                                         11d20s20
        do l=1,4                                                        11d20s20
         noffs(l,1)=0                                                   11d20s20
         noffs(l,2)=0                                                   11d20s20
        end do                                                          11d20s20
        nvec=0                                                          11d20s20
        nme=0                                                           11d20s20
        do l=1,4                                                        3d18s21
         nll(l)=0                                                       3d18s21
        end do                                                          3d18s21
        do ip=0,mynprocg-1                                              11d20s20
         do l=1,4                                                       11d20s20
          if(nct2c(l).gt.0)then                                         11d20s20
           if(ip.lt.mynowprog)then                                      11d20s20
            noffs(l,1)=noffs(l,1)+nint(bc(kfkeep))*ncsf2(l,iarg)        11d28s20
            noffs(l,2)=noffs(l,2)+nint(bc(kfkeep))                      11d20s20
           end if                                                       11d20s20
           nvec=nvec+nint(bc(kfkeep))*ncsf2(l,iarg)                     11d28s20
           nll(l)=nll(l)+nint(bc(kfkeep))                               3d18s21
           sum=sum+bc(kfkeep)                                             11d20s20
           if(ip.eq.mynowprog)nme=nme+nint(bc(kfkeep))                  11d20s20
           kfkeep=kfkeep+nfcn                                             11d20s20
          end if                                                        11d20s20
         end do                                                         11d20s20
        end do
        if(sum.gt.1d-10)then                                            11d20s20
         nl=nl+1                                                        11d20s20
         nhere=nint(sum)                                                11d20s20
         iff=nff2(nclop,isb,2)+if                                       11d20s20
         npre=10                                                        3d18s21
         npred=0                                                        8d30s23
         do l=1,4
          mpre(l)=npre                                                  3d18s21
          npre=npre+nll(l)*(ncsf2(l,iarg)+1)                            3d18s21
          mpred(l)=npred                                                8d30s23
          npred=npred+nll(l)*ncsf2(l,iarg)                              8d30s23
         end do                                                         3d18s21
         if(mynowprog.eq.0)then                                         11d20s20
          ibc(jvcv)=iff2(iff)                                           11d20s20
          ibc(jvcv+1)=nvec+nhere+10                                     3d18s21
          npre=10                                                       3d18s21
          do l=1,4                                                      3d18s21
           ijvcv=jvcv+1+l                                               3d18s21
           ibc(ijvcv)=nll(l)                                            3d18s21
           ijvcv=ijvcv+4                                                3d18s21
           ibc(ijvcv)=npre                                              3d18s21
           npre=npre+nll(l)*(ncsf2(l,iarg)+1)                           3d18s21
          end do                                                        3d18s21
         end if                                                         11d20s20
         nff2(nclop,isb,7)=nff2(nclop,isb,7)+1                          11d20s20
         mme=0                                                          11d20s20
         do l=1,4                                                       11d20s20
          jvcvsto1=jvcv+mpre(l)+noffs(l,2)                              3d18s21
          jvcvsto2=jvcv+mpre(l)+nll(l)+noffs(l,2)*ncsf2(l,iarg)         3d18s21
          lvcvd=jvcvd+mpred(l)+noffs(l,2)*ncsf2(l,iarg)                 8d30s23
          joff=1+ibc(ikkk+l-1)
          do j=0,kkk(l)-1
           kvcc=jvcc(l)+nct2c(l)*j                                      11d20s20
           sz=0d0                                                        11d17s20
           if(ncsf2(l,iarg).gt.0)then                                   1d26s21
            do k=0,ncsf2(l,iarg)-1                                       11d17s20
             sz=sz+bc(kvcc+k)**2                                         11d17s20
            end do                                                       11d17s20
            sz=sqrt(sz/dfloat(ncsf2(l,iarg)))                            11d17s20
            if(sz.gt.1d-10)then
             mme=mme+1                                                   11d20s20
             ibc(jvcvsto1)=j+joff                                       3d18s21
             do k=0,ncsf2(l,iarg)-1                                      11d17s20
              orig=bc(igoal)
              bc(jvcvsto2+k)=bc(kvcc+k)                                 3d18s21
             end do                                                      11d17s20
             llvcvd=lvcvd                                                8d30s23
             do idr=1,nder                                              8d30s23
              kvcc=jvcc(l)+nct2c(l)*(j+kkk(l)*idr)                      8d30s23
              do k=0,ncsf2(l,iarg)-1                                    8d30s23
               orig=bc(igoal)
               bc(llvcvd+k)=bc(kvcc+k)                                   8d30s23
              end do                                                    8d30s23
              llvcvd=llvcvd+needcvd                                       8d30s23
             end do                                                     8d30s23
             lvcvd=lvcvd+ncsf2(l,iarg)                                  8d30s23
             jvcvsto1=jvcvsto1+1                                        3d18s21
             jvcvsto2=jvcvsto2+ncsf2(l,iarg)                            3d18s21
            end if                                                       11d17s20
           end if                                                       1d26s21
          end do                                                        11d17s20
         end do                                                         11d20s20
         nspace=10+nvec+nhere                                           3d18s21
         jvcv=jvcv+nspace                                               11d17s20
         jvcvd=jvcvd+nvec                                               8d30s23
         if(nme.ne.mme)then                                             11d20s20
          write(6,*)('wtf? ')
          do l=1,4
           write(6,*)l,jvcc(l)-ivcc(l)
          end do
          call dws_synca
          call dws_finalize
          stop
         end if
        end if                                                          11d20s20
        do l=1,4                                                        11d20s20
         jvcc(l)=jvcc(l)+ncsf2(l,iarg)                                  11d20s20
        end do                                                          11d20s20
        jfkeep=jfkeep+1                                                 11d20s20
       end do                                                           11d20s20
      end do                                                            11d20s20
      nfdat(1,1)=nl                                                     11d20s20
      call dws_gbor(bc(ivcv),needcv+nder*needcvd+ndorth0+ndorth)        10d17s23
      do ii=mdon+1,mdoo+1                                               3d1s21
       nff22(ii,1)=0                                                    3d1s21
      end do                                                            3d1s21
      jdorth=idorth                                                     9d5s23
      jorth=iorth                                                       9d5s23
      do l=1,4                                                          11d17s20
       if(kkkt(l).gt.0)then                                             11d17s20
        jvcv=ivcv                                                         11d17s20
        jvcvd=ivcvd                                                     8d30s23
        if(ldebug)write(6,*)('computing overlap for l = '),l
        iovr=ibcoff                                                     11d17s20
        ibcoff=iovr+kkkt(l)*kkkt(l)*nderp                               8d30s23
        call enough('cont2ecsfb.  6',bc,ibc)
        do i=iovr,ibcoff-1
         bc(i)=0d0
        end do                                                          11d17s20
        if(ldebug)write(6,*)('no. fcns for this = '),nfdat(1,1)         1d25s21
        do ii=mdon+1,mdoo+1                                             12d18s20
         nff22(ii,1)=0                                                  12d18s20
         nff22(ii,2)=ibcoff                                             12d18s20
        end do                                                          12d18s20
        do if=1,nfdat(1,1)                                              11d20s20
         ipack8=ibc(jvcv)                                               11d20s20
         nclo=popcnt(ipack4(1))                                         11d20s20
         nclop=nclo+1                                                   11d17s20
         nff22(nclop,1)=nff22(nclop,1)+1                                12d18s20
         nff22(nclop,2)=min(nff22(nclop,2),jvcv-ivcv)                   12d18s20
         iarg=nclop-mdon                                                11d17s20
         nvecp=ibc(jvcv+1)                                               3d18s21
         nvecpd=0                                                       8d30s23
         nvecpl=0                                                       8d30s23
         do ll=1,4                                                      8d30s23
          nvecpd=nvecpd+ibc(jvcv+1+ll)*ncsf2(ll,iarg)                   8d30s23
         end do                                                         8d30s23
         do ll=1,l-1                                                    8d30s23
          nvecpl=nvecpl+ibc(jvcv+1+ll)*ncsf2(ll,iarg)                   8d30s23
         end do                                                         8d30s23
         nl=ibc(jvcv+1+l)                                               3d18s21
         iad1=jvcv+ibc(jvcv+5+l)                                        3d18s21
         iad2=iad1+nl                                                   3d18s21
         if(nl.gt.0)then
          itrans=ibcoff                                                  11d20s20
          ibcoff=itrans+nl*ncsf2(l,iarg)                                 11d20s20
          call enough('cont2ecsfb.  7',bc,ibc)
          do i=0,nl-1                                                    11d20s20
           do j=0,ncsf2(l,iarg)-1                                        11d20s20
            ji=iad2+j+ncsf2(l,iarg)*i                                    11d20s20
            ij=itrans+i+nl*j                                             11d20s20
            bc(ij)=bc(ji)                                               11d20s20
           end do                                                        11d20s20
          end do                                                         11d20s20
          iprod=ibcoff                                                   11d20s20
          ibcoff=iprod+nl*nl                                             11d20s20
          call dgemm('n','n',nl,nl,ncsf2(l,iarg),1d0,bc(itrans),nl,      11d20s20
     $        bc(iad2),ncsf2(l,iarg),0d0,bc(iprod),nl,                  11d20s20
     d' cont2ecsfb.  1')
          do i=0,nl-1                                                    11d20s20
           ii=ibc(iad1+i)-1                                              11d20s20
           do j=0,nl-1                                                   11d20s20
            jj=ibc(iad1+j)-1                                             11d20s20
            iiad=iovr+jj+kkkt(l)*ii                                      11d20s20
            iad=iprod+j+nl*i                                             11d20s20
            bc(iiad)=bc(iiad)+bc(iad)                                    11d20s20
           end do                                                        11d20s20
          end do                                                         11d20s20
          lvcvd=jvcvd+nvecpl                                            8d30s23
          do idr=1,nder                                                 8d30s23
           call dgemm('n','n',nl,nl,ncsf2(l,iarg),1d0,bc(itrans),nl,    8d30s23
     $        bc(lvcvd),ncsf2(l,iarg),0d0,bc(iprod),nl,                 8d30s23
     d' cont2ecsfb.  1')
           do i=0,nl-1                                                    11d20s20
            ii=ibc(iad1+i)-1                                              11d20s20
            do j=0,nl-1                                                   11d20s20
             jj=ibc(iad1+j)-1                                             11d20s20
             iiad=iovr+jj+kkkt(l)*(ii+kkkt(l)*idr)                      8d30s23
             iad=iprod+j+nl*i                                             11d20s20
             bc(iiad)=bc(iiad)+bc(iad)                                    11d20s20
             iiad=iovr+ii+kkkt(l)*(jj+kkkt(l)*idr)                      8d31s23
             bc(iiad)=bc(iiad)+bc(iad)                                  8d31s23
            end do                                                        11d20s20
           end do                                                         11d20s20
           lvcvd=lvcvd+needcvd                                          8d30s23
          end do                                                        8d30s23
          ibcoff=itrans                                                  11d20s20
         end if                                                         11d20s20
         jvcv=jvcv+nvecp                                                3d18s21
         jvcvd=jvcvd+nvecpd                                             8d30s23
        end do                                                          11d20s20
        if(ldebug)then                                                  8d31s23
         jovr=iovr                                                      8d31s23
         do idr=0,nder                                                  8d31s23
          jovr=jovr+kkkt(l)*kkkt(l)                                     8d31s23
         end do                                                         8d31s23
        end if
        call orthovrb(bc(iovr),kkkt(l),nok,thrs,ivec,bc,ibc,isb,l,nder) 8d30s23
        if(ldebug)then                                                  1d25s21
         write(6,*)('orthogonalization matrix ')
         call prntm2(bc(ivec),kkkt(l),nok,kkkt(l))
         do idr=0,nder
          write(6,*)('for derivative '),idr
          do i=0,nok-1
           jvec=ivec+kkkt(l)*(i+nok*idr)
           do j=0,kkkt(l)-1
            if(abs(bc(jvec+j)).gt.1d-8)then
             write(line6,1777)j,i
 1777        format(2i3)
             do k=4,6
              if(line6(k:k).eq.' ')line6(k:k)='0'
             end do
             write(6,1778)line6,bc(jvec+j)
 1778        format(a6,x,es22.14)
            end if
           end do
          end do
         end do
        end if                                                          1d25s21
        nfdat(2,l)=kkkt(l)                                              11d18s20
        nfdat(3,l)=nok                                                  11d18s20
        do i=0,kkkt(l)*nok-1                                            9d5s23
         bc(jorth+i)=bc(ivec+i)                                         9d5s23
        end do                                                          11d18s20
        nfdat(4,l)=jorth                                                9d5s23
        jorth=jorth+kkkt(l)*nok                                         9d5s23
        ibcoff=iovr                                                     9d5s23
        jvec=ivec+kkkt(l)*nok                                           9d5s23
        do i=0,kkkt(l)*nok*nder-1                                       9d5s23
         bc(jdorth+i)=bc(jvec+i)                                        9d5s23
        end do                                                          9d5s23
        jdorth=jdorth+kkkt(l)*kkkt(l)*nder                              9d5s23
       end if                                                           11d20s20
      end do                                                            11d20s20
      ibcoff=jorth                                                      9d5s23
      nfdat(5,1)=ivcv                                                   11d20s20
      ntackon=0                                                         8d7s24
      if(nder.gt.0)then                                                 9d1s23
c
c     memory mapping so far ...
c     vcv,for l=1,4 ortho, (l=1,4,for der dortho),vcvd                  9d5s23
c     but rest of code is expecting
c     vcv, for l=1,4 ortho,
c     so let us store vcvd and dortho in distributed memory ...
c
       need2save=needcvd*nder+ndorth                                    9d5s23
       idvmt(2)=need2save                                               9d1s23
       call ilimts(1,need2save,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)9d1s23
       i18=1                                                            9d1s23
       il8=il                                                           9d1s23
       ih8=ih                                                           9d1s23
       nhere=ih+1-il                                                    10d17s23
       iorgn=idorth+il-1                                                9d5s23
       ntackon=nhere                                                    8d7s24
       do i=0,nhere-1                                                   8d7s24
        bc(ibcoff+i)=bc(iorgn+i)                                        8d7s24
       end do                                                           8d7s24
       ibcoff=ibcoff+nhere                                              8d7s24
      end if                                                            9d1s23
      if(ldebug)then                                                    1d25s21
       write(6,*)('what we have for nff22: ')
       nsum=0                                                            12d18s20
       do i=mdon+1,mdoo+1
        write(6,*)i,nff22(i,1),nff22(i,2)
        nsum=nsum+nff22(i,1)                                             12d18s20
       end do
       write(6,*)('nsum vs. nfdat: '),nsum,nfdat(1,1)
      end if                                                            1d25s21
      return
      end
