c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine mover2qs(vecin,nin,veco,no,ipoint,nps,nroot)
      implicit real*8 (a-h,o-z)
      dimension vecin(*),veco(no,*),ipoint(*)                           10d30s24
      j=1
      do l=1,nroot
       do i=1,no
        veco(i,l)=0d0
        do k=1,nps
         if(ipoint(k).eq.i)go to 1
        end do
        veco(i,l)=vecin(j)
        j=j+1
    1   continue
       end do
      end do
      return
      end
      subroutine movback2ps(vecin,veco,ncsft,nps,ipoint)
      implicit real*8 (a-h,o-z)
      dimension vecin(ncsft),veco(nps),ipoint(nps)
      do i=1,nps
       write(6,*)i,ipoint(i),veco(i)
       vecin(ipoint(i))=veco(i)
      end do
      return
      end
      subroutine hccsflr(gl,ndl,vecr,ndr,ibasisl,ibasisr,ncsf,nec,nfcnl,9d4s24
     $     nfcnr,ih0a,i2e,nrootz,multh,mdon,mdoo,iptrbl,iptrbr,ixw1,    9d4s24
     $     ixw2,bc,ibc)                                                 9d4s24
      implicit real*8 (a-h,o-z)                                         9d4s24
      integer*1 nab1(2),nab2(2)                                         9d4s24
      integer*8 gandcc,gandco,gandcb,itesta,itestb                      9d4s24
      include "common.hf"                                               3d31s20
      include "common.store"                                            9d4s24
      include "common.mrci"                                             6d19s19
      include "common.print"                                            1d13s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      dimension vecr(ndr,*),gl(ndl,*),ih0a(*),i2e(*),ibasisl(3,*),      9d4s24
     $     ibasisr(3,*),ncsf(*),multh(8,8),iptrbl(2,*),iptrbr(2,*),     9d4s24
     $     ioxx(2),nab4(2,3),itest(64,2)                                9d4s24
      ibcoffo=ibcoff                                                    4d8s20
      do i=1,nrootz                                                     7d11s19
       do j=1,ndl                                                       3d25s21
        gl(j,i)=0d0                                                     7d11s19
       end do                                                           7d11s19
      end do                                                            7d11s19
      ioff=1                                                            9d4s24
      loopit=0                                                          9d4s24
      do if=1,nfcnr                                                     9d4s24
       nclor=ibasisr(1,if)                                              9d4s24
       nopenr=nec-2*nclor                                               9d4s24
       nclorp=nclor+1                                                   9d4s24
       iarg=nclorp-mdon                                                 9d4s24
       iic=iptrbr(1,nclorp)+ibasisr(2,if)-1                             9d4s24
       iio=iptrbr(2,nclorp)+ibasisr(3,if)-1                             9d4s24
       joff=1                                                           9d4s24
       do jf=1,nfcnl                                                    9d4s24
        nclol=ibasisl(1,jf)                                              9d4s24
        nclolp=nclol+1                                                   9d4s24
        jarg=nclolp-mdon                                                 9d4s24
        if(mod(loopit,mynprocg).eq.mynowprog)then                       9d4s24
         loopit=mynowprog                                               9d6s24
         nopenl=nec-2*nclol                                             9d4s24
         nn=ncsf(jarg)*ncsf(iarg)                                        9d4s24
         jjc=iptrbl(1,nclolp)+ibasisl(2,jf)-1                            9d4s24
         jjo=iptrbl(2,nclolp)+ibasisl(3,jf)-1                            9d4s24
         gandcc=ieor(ibc(iic),ibc(jjc))                                 9d4s24
         gandco=ieor(ibc(iio),ibc(jjo))                                 9d4s24
         gandcb=ior(gandcc,gandco)                                       9d4s24
         ndifb=popcnt(gandcb)                                            9d4s24
         if(ndifb.le.4)then                                              9d4s24
          ndifs=popcnt(gandco)                                           9d4s24
          ndifd=popcnt(gandcc)                                           9d4s24
          if(ndifs.eq.0.and.ndifd.eq.0)then                              9d4s24
           sum=0d0                                                       9d4s24
           do i=1,norb                                                   9d4s24
            is=ism(i)                                                    9d4s24
            ig=irel(i)-1                                                 9d4s24
            jmat=ih0a(is)+ig*(irefo(is)+1)                               9d4s24
            if(btest(ibc(iic),i))then                                    9d4s24
             sum=sum+2d0*bc(jmat)                                        9d4s24
            end if                                                       9d4s24
            if(btest(ibc(iio),i))then                                    9d4s24
             sum=sum+bc(jmat)                                            9d4s24
            end if                                                       9d4s24
           end do                                                        9d4s24
           do i1=1,norb                                                  9d4s24
            if(btest(ibc(iic),i1))then                                   9d4s24
             jsa=ism(i1)                                                  9d4s24
             jga=irel(i1)-1                                               9d4s24
             do i2=1,norb                                                 9d4s24
              if(btest(ibc(iic),i2).or.btest(ibc(iio),i2))then           9d4s24
               jsb=ism(i2)                                                9d4s24
               jgb=irel(i2)-1                                             9d4s24
               xint=getint(i2e,jsa,jsa,jsb,jsb,jga,jga,jgb,jgb,bc,ibc)    9d4s24
               sum=sum+2d0*xint                                           9d4s24
               xint=getint(i2e,jsa,jsb,jsb,jsa,jga,jgb,jgb,jga,bc,ibc)    9d4s24
               sum=sum-xint                                               9d4s24
              end if                                                      9d4s24
             end do                                                      9d4s24
            end if                                                       9d4s24
            if(btest(ibc(iio),i1))then                                   5d20s21
             jsa=ism(i1)                                                 9d4s24
             jga=irel(i1)-1                                              9d4s24
             do i2=1,norb                                                9d4s24
              if(btest(ibc(iio),i2).and.i1.ne.i2)then                    9d4s24
               jsb=ism(i2)                                               9d4s24
               jgb=irel(i2)-1                                            9d4s24
               xint=getint(i2e,jsa,jsa,jsb,jsb,jga,jga,jgb,jgb,bc,ibc)   9d4s24
               sum=sum+0.5d0*xint                                       9d4s24
              end if                                                     9d4s24
             end do                                                      9d4s24
            end if                                                       9d4s24
           end do                                                        9d4s24
           do ir=1,nrootz                                                9d4s24
            do j=0,ncsf(jarg)-1                                          9d4s24
             gl(joff+j,ir)=gl(joff+j,ir)+sum*vecr(ioff+j,ir)             5d20s21
            end do                                                       9d4s24
           end do                                                        9d4s24
           imat=ibcoff                                                  9d4s24
           ibcoff=imat+nn                                               9d4s24
           call enough('hccsflr.1',bc,ibc)                              9d4s24
           do i1=1,norb-1                                               9d4s24
            if(btest(ibc(iio),i1))then                                  9d4s24
             do i2=i1+1,norb                                            9d4s24
              if(btest(ibc(iio),i2))then                                9d4s24
               itesta=ibc(iic)                                          9d4s24
               itestb=ibc(iio)                                          9d4s24
               nopenkk=nopenr-2                                         9d4s24
               karg=iarg+1                                              9d4s24
               nab1(1)=i2                                               9d4s24
               nab1(2)=i1                                               9d4s24
               nab2(1)=nab1(2)                                           12d6s20
               nab2(2)=nab1(1)                                           12d6s20
               jsa=ism(nab1(1))                                              11d13s20
               jga=irel(nab1(1))-1                                      9d4s24
               jsb=ism(nab1(2))                                         9d4s24
               jgb=irel(nab1(2))-1                                      9d4s24
               jsc=ism(nab2(1))                                         9d4s24
               jgc=irel(nab2(1))-1                                      9d4s24
               jsd=ism(nab2(2))                                         9d4s24
               jgd=irel(nab2(2))-1                                      9d4s24
               nqq=karg+mdon-1                                          9d5s24
               if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                itestb=ibclr(itestb,i2)                                 9d4s24
                itestb=ibclr(itestb,i1)                                 9d4s24
                itesta=ibset(itesta,i1)                                 9d4s24
                call gandc(ibc(jjc),ibc(jjo),itesta,itestb,nopenl,      9d4s24
     $               nopenkk,jarg,iarg,ncsf,norb,ixw1,ixw2,nnot1,nab1,  9d4s24
     $               iwpb1,iwpk1,ncsfmid1,bc,ibc)                       9d4s24
                call gandc(itesta,itestb,ibc(iic),ibc(iio),nopenkk,     9d4s24
     $               nopenr,karg,iarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,   9d4s24
     $               iwpb2,iwpk2,ncsfmid2,bc,ibc)                       9d4s24
                call genmat(ncsf(jarg),ncsf(karg),ncsf(iarg),           9d4s24
     $               ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,iprod,   9d4s24
     $               bc,ibc)                                            9d4s24
                do i=0,ncsf(iarg)*ncsf(iarg)-1                          9d4s24
                 bc(imat+i)=bc(iprod+i)                                  4d13s21
                end do                                                   4d13s21
                ibcoff=iprod                                                11d13s20
               else                                                        11d13s20
                do i=0,ncsf(iarg)*ncsf(iarg)-1                          9d4s24
                 bc(imat+i)=0d0                                          4d13s21
                end do                                                   4d13s21
               end if                                                   9d4s24
               do i=0,ncsf(iarg)-1                                       4d13s21
                ii=imat+i*(ncsf(iarg)+1)                                9d4s24
                bc(ii)=bc(ii)-1d0                                        4d13s21
               end do                                                    4d13s21
               xint=getint(i2e,jsa,jsb,jsc,jsd,jga,jgb,jgc,jgd,bc,ibc)  9d4s24
               call dgemm('n','n',ncsf(jarg),nrootz,ncsf(iarg),         9d4s24
     $                   xint,bc(imat),ncsf(jarg),vecr(ioff,1),ndr,     9d4s24
     $                   1d0,gl(joff,1),ndl,'hccsflr.1')                9d4s24
              end if                                                    9d4s24
             end do                                                       11d13s20
            end if                                                      9d4s24
           end do                                                        11d13s20
           ibcoff=imat                                                  9d4s24
          else if(ndifs.eq.2.and.ndifb.eq.2)then                         9d4s24
           do i=1,norb                                                  9d4s24
            if(btest(gandco,i))then                                     9d4s24
             if((btest(ibc(jjo),i).and..not.btest(ibc(iic),i)).or.      9d4s24
     $          (btest(ibc(jjc),i).and.btest(ibc(iio),i)))then          9d4s24
              nab4(1,1)=i                                               9d4s24
             else                                                       9d4s24
              nab4(2,1)=i                                               9d4s24
             end if                                                     9d4s24
            end if                                                      9d4s24
           end do                                                       9d4s24
           call gandc(ibc(jjc),ibc(jjo),ibc(iic),ibc(iio),nopenl,       9d4s24
     $        nopenr,jarg,iarg,ncsf,norb,ixw1,ixw2,nnot1,nab1,iwpb1,    9d4s24
     $        iwpk1,ncsfmid1,bc,ibc)                                    9d4s24
           idvtmp=ibcoff                                                9d4s24
           ibcoff=idvtmp+ncsf(jarg)*nrootz                              9d4s24
           call enough('hccsflr.  2',bc,ibc)
           call xtimesn(ncsf(jarg),ncsf(iarg),ncsfmid1,nrootz,iwpb1,    9d4s24
     $        iwpk1,vecr(ioff,1),ndr,bc(idvtmp),ncsf(jarg),1d0,         9d4s24
     $        0d0,bc,ibc)                                               9d4s24
           jsb=ism(nab1(1))                                             9d4s24
           jsk=ism(nab1(2))                                             9d4s24
           jgb=irel(nab1(1))-1                                          9d4s24
           jgk=irel(nab1(2))-1                                          9d4s24
           iad=ih0a(jsb)+jgb+irefo(jsb)*jgk                             9d4s24
           sum=bc(iad)                                                  9d4s24
           nok=0                                                         11d13s20
           do i=1,norb                                                   11d13s20
            itest(i,1)=0                                                   11d13s20
            itest(i,2)=0                                                    11d13s20
            if(btest(ibc(iic),i))itest(i,1)=2                           9d4s24
            if(btest(ibc(iio),i))itest(i,1)=1                           9d4s24
            if(btest(ibc(jjc),i))itest(i,2)=2                           9d4s24
            if(btest(ibc(jjo),i))itest(i,2)=1                           9d4s24
            ixn=min(itest(i,1),itest(i,2))
            if(ixn.gt.0)then                                             11d13s20
             nok=nok+1                                                   11d13s20
             itest(nok,1)=ixn                                            11d13s20
             itest(nok,2)=i                                              11d13s20
            end if                                                       11d13s20
           end do                                                        11d13s20
           jsa=ism(nab4(2,1))                                             11d13s20
           jsb=ism(nab4(1,1))                                             11d13s20
           jga=irel(nab4(2,1))-1                                        9d4s24
           jgb=irel(nab4(1,1))-1                                        9d4s24
           do i=1,nok                                                    11d13s20
            js=ism(itest(i,2))                                           11d13s20
            jg=irel(itest(i,2))-1                                       9d4s24
            xint=getint(i2e,jsa,jsb,js,js,jga,jgb,jg,jg,bc,ibc)         9d4s24
            if(itest(i,1).eq.2)then                                     9d4s24
             sum=sum+xint*2d0                                           9d4s24
            else                                                        9d4s24
             sum=sum+xint                                               9d4s24
            end if                                                      9d4s24
            if(itest(i,1).eq.2)then                                     9d4s24
             xint=getint(i2e,jsa,js,js,jsb,jga,jg,jg,jgb,bc,ibc)         9d4s24
             sum=sum-xint                                               9d4s24
            end if                                                      9d4s24
           end do                                                       9d4s24
           jdvtmp=idvtmp                                                 9d4s24
           do ir=1,nrootz                                                9d4s24
            do j=0,ncsf(jarg)-1                                          9d4s24
             gl(joff+j,ir)=gl(joff+j,ir)+sum*bc(jdvtmp+j)                9d4s24
            end do                                                       9d4s24
            jdvtmp=jdvtmp+ncsf(jarg)                                     9d4s24
           end do                                                        9d4s24
           ibcoff=idvtmp                                                 9d4s24
           do i=1,nok                                                    9d4s24
            if(itest(i,1).eq.1)then                                      9d4s24
             itesta=ibc(jjc)                                             9d4s24
             itestb=ibc(jjo)                                             9d4s24
             nopenkk=nopenl                                              9d4s24
c
c     anihilate common
c
             if(btest(itesta,itest(i,2)))then                            9d4s24
              itesta=ibclr(itesta,itest(i,2))                            9d4s24
              itestb=ibset(itestb,itest(i,2))                            9d4s24
              karg=jarg-1                                                9d4s24
              nopenkk=nopenkk+1                                          9d4s24
             else                                                        9d4s24
              itestb=ibclr(itestb,itest(i,2))                            9d4s24
              karg=jarg                                                  9d4s24
              nopenkk=nopenkk-1                                          9d4s24
             end if                                                      9d4s24
c
c     create ket
c
             if(btest(itestb,nab4(2,1)))then                             9d4s24
              itesta=ibset(itesta,nab4(2,1))                             9d4s24
              itestb=ibclr(itestb,nab4(2,1))                             9d4s24
              karg=karg+1                                                9d4s24
              nopenkk=nopenkk-1                                          9d4s24
             else                                                        9d4s24
              itestb=ibset(itestb,nab4(2,1))                             9d4s24
              nopenkk=nopenkk+1                                          9d4s24
             end if                                                      9d4s24
             nqq=karg+mdon-1                                             9d4s24
             if(nqq.ge.mdon.and.nqq.le.mdoo)then                         9d4s24
              call gandc(ibc(jjc),ibc(jjo),itesta,itestb,nopenl,         9d4s24
     $           nopenkk,jarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,      9d4s24
     $           iwpb1,iwpk1,ncsfmid1,bc,ibc,bc,ibc)                    9d4s24
              call gandc(itesta,itestb,ibc(iic),ibc(iio),nopenkk,nopenr, 9d4s24
     $            karg,iarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2, 9d4s24
     $            ncsfmid2,bc,ibc,bc,ibc)                               9d4s24
              if(nnot1.eq.2.and.nnot2.eq.2)then                          9d4s24
               jsa=ism(nab1(1))                                          9d4s24
               jga=irel(nab1(1))-1                                       9d4s24
               jsb=ism(nab1(2))                                          9d4s24
               jgb=irel(nab1(2))-1                                       9d4s24
               jsc=ism(nab2(1))                                          9d4s24
               jgc=irel(nab2(1))-1                                       9d4s24
               jsd=ism(nab2(2))                                          9d4s24
               jgd=irel(nab2(2))-1                                       9d4s24
               xint=getint(i2e,jsa,jsb,jsc,jsd,jga,jgb,jgc,jgd,bc,ibc)   9d4s24
               itmpgb=ibcoff                                             9d4s24
               ibcoff=itmpgb+ncsf(jarg)*nrootz                           9d4s24
               call enough('hccsflr.  3',bc,ibc)
               do iqq=itmpgb,ibcoff-1                                    9d4s24
                bc(iqq)=0d0                                              9d4s24
               end do                                                    9d4s24
               call genmatn(ncsf(jarg),ncsf(karg),ncsf(iarg),            9d4s24
     $             ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,           9d4s24
     $     vecr(ioff,1),ndr,nrootz,bc(itmpgb),0d0,bc,ibc)               9d4s24
               jtmpgb=itmpgb                                             9d4s24
               do ir=1,nrootz
                do j=0,ncsf(jarg)-1                                      9d4s24
                 gl(joff+j,ir)=gl(joff+j,ir)+bc(jtmpgb+j)*xint           9d4s24
                end do                                                   9d4s24
                jtmpgb=jtmpgb+ncsf(jarg)                                 9d4s24
               end do                                                    9d4s24
               ibcoff=itmpgb                                             9d4s24
              else if(max(nnot1,nnot2).lt.2)then                         9d4s24
               write(6,*)('c:expecting 2s, but got otherwise')           9d4s24
               call dcbit(itesta,norb,'itesta')                          9d4s24
               call dcbit(itestb,norb,'itestb')                          9d4s24
               stop 'hccsflr'                                            9d4s24
              end if                                                     9d4s24
             end if                                                      9d4s24
            end if                                                       11d13s20
           end do                                                        11d13s20
          else                                                           9d4s24
           nnot=0                                                        9d4s24
           ipssx=0                                                       9d4s24
           if(ndifs.eq.4.and.ndifb.eq.4)then                             9d4s24
            nnot=4                                                       9d4s24
            ioxx(1)=1                                                    9d4s24
            ioxx(2)=1                                                    9d4s24
            do i=1,norb                                                  9d4s24
             if(btest(gandcb,i))then                                     9d4s24
              if((btest(ibc(iic),i).and.btest(ibc(jjo),i)).or.          9d4s24
     $               (btest(ibc(iio),i).and..not.btest(ibc(jjc),i)))then9d4s24
               nab4(2,ioxx(2))=i                                         9d4s24
               ioxx(2)=ioxx(2)+1                                         9d4s24
              else                                                       9d4s24
               nab4(1,ioxx(1))=i                                         9d4s24
               ioxx(1)=ioxx(1)+1                                         9d4s24
              end if                                                     9d4s24
             end if                                                      9d4s24
            end do                                                       9d4s24
           else if(ndifb.eq.3)then                                       9d4s24
            nnot=3                                                       9d4s24
            ioxx(1)=1                                                    9d4s24
            ioxx(2)=1                                                    9d4s24
            iswap=0                                                      9d4s24
            do i=1,norb                                                  9d4s24
             if(btest(gandcb,i))then                                     9d4s24
              if(btest(gandcc,i).and.                                    9d4s24
     $        ((btest(ibc(jjc),i).and..not.btest(ibc(iio),i)).or.       9d4s24
     $         (btest(ibc(iic),i).and..not.btest(ibc(jjo),i))))then     9d4s24
               if(btest(ibc(iic),i))iswap=1                              2d6s23
               nab4(1,1)=i                                               9d4s24
               nab4(1,2)=i                                               9d4s24
              else                                                       9d4s24
               nab4(2,ioxx(2))=i                                         9d4s24
               ioxx(2)=ioxx(2)+1                                         9d4s24
              end if                                                     9d4s24
             end if                                                      9d4s24
            end do                                                       9d4s24
            if(iswap.ne.0)then                                           9d4s24
             icpy=nab4(1,1)                                              9d4s24
             nab4(1,1)=nab4(2,1)                                         9d4s24
             nab4(2,1)=icpy                                              9d4s24
             icpy=nab4(1,2)                                              9d4s24
             nab4(1,2)=nab4(2,2)                                         9d4s24
             nab4(2,2)=icpy                                              9d4s24
             nbt=0                                                       9d4s24
             if(btest(ibc(jjc),nab4(1,2)).and.                           9d4s24
     $           .not.btest(ibc(jjc),nab4(1,1)))nbt=1                   9d4s24
            else                                                         9d4s24
             nbt=0                                                       9d4s24
             if(btest(ibc(iic),nab4(2,2)).and.                           9d4s24
     $            .not.btest(ibc(iic),nab4(2,1)))nbt=1                  9d4s24
            end if                                                       9d4s24
            if(nbt.ne.0)then                                             9d4s24
             nab4(1,1)=nab4(1,2)                                         9d4s24
             nab4(2,1)=nab4(2,2)                                         9d4s24
            end if                                                       9d4s24
           else if(ndifs.eq.0.and.ndifd.eq.2)then                        9d4s24
            nnot=3                                                       9d4s24
            do i=1,norb                                                  9d4s24
             if(btest(gandcb,i))then                                     9d4s24
              if(btest(ibc(jjc),i))then                                  2d6s23
               nab4(1,1)=i                                               9d4s24
               nab4(1,2)=i                                               9d4s24
              else                                                       9d4s24
               nab4(2,1)=i                                               9d4s24
               nab4(2,2)=i                                               9d4s24
              end if                                                     9d4s24
             end if                                                      9d4s24
            end do                                                       9d4s24
           end if                                                        9d4s24
           if(nnot.eq.3)then                                             9d4s24
            ipssx=1                                                      9d4s24
           else if(nnot.eq.4)then                                        9d4s24
            ipssx=3                                                      9d4s24
           end if                                                        9d4s24
           do ipss=1,ipssx                                              9d4s24
            if(ipss.eq.1)then                                           9d4s24
             iu1=1                                                      9d4s24
             iu2=1                                                      9d4s24
            else if(ipss.eq.2)then                                      9d4s24
             iu1=1                                                      9d4s24
             iu2=2                                                      9d4s24
            else                                                        9d4s24
             iu1=2                                                      9d4s24
             iu2=1                                                      9d4s24
            end if                                                      9d4s24
            itesta=ibc(jjc)                                             9d4s24
            itestb=ibc(jjo)                                             9d4s24
            if(btest(itesta,nab4(1,iu1)))then                           9d4s24
             itesta=ibclr(itesta,nab4(1,iu1))                           9d4s24
             itestb=ibset(itestb,nab4(1,iu1))                           9d4s24
             nopenkk=nopenl+1                                           9d4s24
             karg=jarg-1                                                9d4s24
            else if(btest(itestb,nab4(1,iu1)))then                      9d4s24
             itestb=ibclr(itestb,nab4(1,iu1))                           9d4s24
             nopenkk=nopenl-1                                              11d13s20
             karg=jarg                                                  9d4s24
            end if                                                        11d13s20
            if(btest(itestb,nab4(2,iu2)))then                           9d4s24
             itesta=ibset(itesta,nab4(2,iu2))                           9d4s24
             itestb=ibclr(itestb,nab4(2,iu2))                           9d4s24
             nopenkk=nopenkk-1                                          9d4s24
             karg=karg+1                                                9d4s24
            else if(btest(itesta,nab4(2,iu2)))then                      9d4s24
             write(6,*)('already double in nab4(2,1) = '),nab4(2,iu2)   9d4s24
             stop 'nab4(2,1)'                                           9d4s24
            else                                                        9d4s24
             itestb=ibset(itestb,nab4(2,iu2))                           9d4s24
             nopenkk=nopenkk+1                                          9d4s24
            end if                                                      9d4s24
            nqq=karg+mdon-1                                             9d4s24
            if(nqq.ge.mdon.and.nqq.le.mdoo)then                         9d4s24
             call gandc(ibc(jjc),ibc(jjo),itesta,itestb,nopenl,         9d4s24
     $        nopenkk,jarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,iwpb1,   9d4s24
     $           iwpk1,ncsfmid1,bc,ibc)                                 9d4s24
             call gandc(itesta,itestb,ibc(iic),ibc(iio),nopenkk,nopenr, 9d4s24
     $        karg,iarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,     9d4s24
     $        ncsfmid2,bc,ibc)                                          9d4s24
             if(nnot1.eq.2.and.nnot2.eq.2)then                          9d4s24
              jsa=ism(nab1(1))                                          9d4s24
              jga=irel(nab1(1))-1                                       9d4s24
              jsb=ism(nab1(2))                                          9d4s24
              jgb=irel(nab1(2))-1                                       9d4s24
              jsc=ism(nab2(1))                                          9d4s24
              jgc=irel(nab2(1))-1                                       9d4s24
              jsd=ism(nab2(2))                                          9d4s24
              jgd=irel(nab2(2))-1                                       9d4s24
              xint=getint(i2e,jsa,jsb,jsc,jsd,jga,jgb,jgc,jgd,bc,ibc)   9d4s24
              if(nab1(1).eq.nab2(1).and.nab1(2).eq.nab2(2))then         9d4s24
               xint=xint*0.5d0                                          9d4s24
              end if
              itmpgb=ibcoff                                             9d4s24
              ibcoff=itmpgb+ncsf(jarg)*nrootz                           9d4s24
              call enough('hccsflr.  4',bc,ibc)                         9d4s24
              do i=itmpgb,ibcoff-1                                      9d4s24
               bc(i)=0d0                                                9d4s24
              end do                                                    9d4s24
              call genmatn(ncsf(jarg),ncsf(karg),ncsf(iarg),            9d4s24
     $             ncsfmid1,iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,           9d4s24
     $             vecr(ioff,1),ndr,nrootz,bc(itmpgb),0d0,bc,ibc)       9d4s24
              jtmpgb=itmpgb
              do ir=1,nrootz
               do j=0,ncsf(jarg)-1
                gl(joff+j,ir)=gl(joff+j,ir)+bc(jtmpgb+j)*xint           9d4s24
               end do                                                   9d4s24
               jtmpgb=jtmpgb+ncsf(jarg)                                 9d4s24
              end do                                                    9d4s24
              ibcoff=itmpgb                                             9d4s24
              if(ipss.eq.2)go to 22                                     9d4s24
             end if                                                     9d4s24
            end if                                                      9d4s24
           end do                                                       9d4s24
   22      continue                                                     9d4s24
          end if                                                         9d4s24
         end if                                                          9d4s24
        end if                                                          9d4s24
        joff=joff+ncsf(jarg)
        loopit=loopit+1                                                 9d4s24
       end do                                                           9d4s24
       ioff=ioff+ncsf(iarg)
      end do                                                            9d4s24
      call dws_gsumf(gl,ndl*nrootz)                                     9d4s24
      ibcoff=ibcoffo                                                    9d4s24
      return                                                            9d4s24
      end                                                               9d4s24
      subroutine hccsf(vecx,gx,ncsft,ibasis,ncsf,iptr,nfcn,ih0a,i2e,    7d11s19
     $     hdig,nrootz,mdon,isorb,idorb,icsfpd,nec,chc,mysym,iptrbit,   11d16s20
     $     ixw1,ixw2,mdoo,bc,ibc)                                       11d10s22
      implicit real*8 (a-h,o-z)                                         7d11s19
      external second                                                   8d1s19
      logical ltest,lnew                                                11d16s20
      integer*1 isorb(*),idorb(*),icode(64),imap(64),nab(2),isame(64),  11d16s20
     $     nab1(2),nab2(2)                                              11d16s20
      integer*8 ipack,itesta,itestb                                     11d16s20
      integer*2 ipack2(4)                                               12d1s19
      equivalence (ipack,ipack2)                                        12d1s19
      dimension vecx(ncsft,nrootz),gx(ncsft,nrootz),iptr(4,*),ih0a(*),  7d11s19
     $    i2e(*),hdig(*),ibasis(3,*),icsfpd(*),ncsf(*),chc(*),nother(2),11d16s20
     $     iptrbit(2,mdoo+1,*),nab4(2,3),itest(64,2)                    11d16s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d15s19
     $     mynnode                                                      5d15s19
      include "common.store"                                            7d11s19
      include "common.mrci"                                             6d19s19
      igoul=67301746
      ltest=.false.                                                     11d16s20
      lnew=.true.                                                       11d16s20
      if(ltest)write(6,*)('hi, my name is hccsf')
      ngcode=0
      do i=1,nrootz                                                     7d11s19
       do j=1,ncsft                                                     7d11s19
        gx(j,i)=0d0                                                     7d11s19
       end do                                                           7d11s19
      end do                                                            7d11s19
      ips=0                                                             7d11s19
      loopit=0                                                          7d11s19
      do if=1,nfcn                                                      7d11s19
       nclo=ibasis(1,if)                                                7d11s19
       nclop=nclo+1                                                     7d11s19
       iarg=nclop-mdon                                                  7d11s19
       if(mynowprog.eq.0)then                                           7d11s19
        do i=1,nrootz                                                   7d11s19
         do j=1,ncsf(iarg)                                              7d11s19
          gx(j+ips,i)=gx(j+ips,i)+hdig(if)*vecx(j+ips,i)                7d11s19
         end do                                                         7d11s19
        end do                                                          7d11s19
       end if                                                           7d11s19
       nopen=nec-2*nclo                                                 7d11s19
       ip=ips+1                                                         7d11s19
       ic=iptr(2,nclop)+nclo*(ibasis(2,if)-1)                           7d11s19
       io=iptr(4,nclop)+nopen*(ibasis(3,if)-1)                          12d9s19
       iic=iptrbit(1,nclop,mysym)+ibasis(2,if)-1                        11d16s20
       iio=iptrbit(2,nclop,mysym)+ibasis(3,if)-1                        11d16s20
       jps=ips                                                          12d11s19
       do jf=if,nfcn                                                    12d9s19
        ncloj=ibasis(1,jf)                                              7d11s19
        if(iabs(ncloj-nclo).le.2)then                                   12d9s19
         nclopj=ncloj+1                                                  6d11s19
         nopenj=nec-2*ncloj
         jarg=ncloj+1-mdon                                               6d11s19
         jo=iptr(4,nclopj)+nopenj*(ibasis(3,jf)-1)                       7d11s19
         jc=iptr(2,nclopj)+ncloj*(ibasis(2,jf)-1)                        7d11s19
         jp=jps+1
         jjo=iptrbit(2,nclopj,mysym)+ibasis(3,jf)-1                     11d16s20
         jjc=iptrbit(1,nclopj,mysym)+ibasis(2,jf)-1                     11d16s20
         if(lnew.or.ltest)then                                          11d16s20
          call gandc4(ibc(jjc),ibc(jjo),ibc(iic),ibc(iio),nopenj,nopen,  11d16s20
     $        norb,nnot4,nab4,bc,ibc)                                   11d14s22
          if(lnew)nnotu=nnot4                                           11d16s20
         end if                                                         11d16s20
         if(nnotu.gt.0)then                                             11d16s20
          if(mod(loopit,mynprocg).eq.mynowprog)then                      12d9s19
           loopit=mynowprog                                             9d10s24
           if(nnotu.gt.1.and.nab4(1,1).eq.0)then                        11d16s20
            write(6,*)('gandc4 failure! '),nnot4                        11d16s20
            write(6,*)nab4                                              11d16s20
            stop                                                        11d16s20
           end if                                                       11d16s20
           if(ltest)then
            write(6,*)('nnotu = '),nnotu
            write(6,*)('bra c:o '),(idorb(jc+i),i=0,ncloj-1),(':'),
     $         (isorb(jo+i),i=0,nopenj-1)
            write(6,*)('ket c:o '),(idorb(ic+i),i=0,nclo-1),(':'),
     $         (isorb(io+i),i=0,nopen-1)
           end if
           ibcsav=ibcoff                                                  12d9s19
           if(ltest.or.lnew)then                                        11d16s20
            jhtmpgg=ibcoff                                              11d16s20
            ibcoff=jhtmpgg+ncsf(jarg)*nrootz                            11d16s20
            if(if.ne.jf)then                                            11d16s20
             ihtmpgg=ibcoff                                             11d16s20
             ibcoff=ihtmpgg+ncsf(iarg)*nrootz                           11d16s20
            end if                                                      11d16s20
            call enough('hccsf.  1',bc,ibc)
            if(nnot4.eq.1)then                                             11d13s20
             do i=0,nrootz-1                                            11d16s20
              do j=0,ncsf(iarg)-1                                       11d16s20
               ji=jhtmpgg+j+ncsf(jarg)*i                                   11d13s20
               bc(ji)=0d0                                                  11d13s20
              end do                                                       11d13s20
             end do                                                        11d13s20
             if(ltest.and.lnew)write(6,*)('ncsf: '),ncsf(jarg)
             do i1=0,nopen-2                                               11d13s20
              do i2=i1+1,nopen-1                                           11d13s20
               if(ltest.and.lnew)then
                write(6,*)('anihilate '),isorb(io+i2)
                write(6,*)('create '),isorb(io+i1)
               end if
               itesta=ibc(iic)                                             11d13s20
               itestb=ibc(iio)                                             11d13s20
               nopenk=nopen-2                                              11d13s20
               karg=iarg+1                                                 11d13s20
               nab1(1)=isorb(io+i2)                                     12d6s20
               nab1(2)=isorb(io+i1)                                     12d6s20
               nab2(1)=nab1(2)                                          12d6s20
               nab2(2)=nab1(1)                                          12d6s20
               if(ltest.and.lnew)write(6,*)('nnot1, nab1 '),
     $             nnot1,nab1,ncsfmid1,ncsf(karg)
               jsa=ism(nab1(1))                                              11d13s20
               jga=irel(nab1(1))-1                                           11d13s20
               jsb=ism(nab1(2))                                              11d13s20
               jgb=irel(nab1(2))-1                                           11d13s20
               jsc=ism(nab2(1))                                              11d13s20
               jgc=irel(nab2(1))-1                                           11d13s20
               jsd=ism(nab2(2))                                              11d13s20
               jgd=irel(nab2(2))-1                                           11d13s20
               xint=getint(i2e,jsa,jsb,jsc,jsd,jga,jgb,jgc,jgd,bc,ibc)  11d15s22
               nqq=karg+mdon-1
               if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                itestb=ibclr(itestb,isorb(io+i2))                                     11d13s20
                itestb=ibclr(itestb,isorb(io+i1))                                     11d13s20
                itesta=ibset(itesta,isorb(io+i1))                           11d13s20
                call gandc(ibc(jjc),ibc(jjo),itesta,itestb,nopenj,      12d6s20
     $               nopenk,jarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,   12d6s20
     $               iwpb1,iwpk1,ncsfmid1,bc,ibc)                       11d14s22
                call gandc(itesta,itestb,ibc(iic),ibc(iio),nopenk,nopen, 11d16s20
     $         karg,iarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,    11d13s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
                iprod=ibcoff                                               11d13s20
                ibcoff=iprod+ncsf(jarg)*nrootz                           11d16s20
                call enough('hccsf.  2',bc,ibc)
                call genmatn(ncsf(jarg),ncsf(karg),ncsf(iarg),ncsfmid1, 12d4s20
     $               iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,vecx(ip,1),ncsft, 12d4s20
     $               nrootz,bc(iprod),0d0,bc,ibc)                       11d10s22
                if(ltest.and.lnew)then
                 write(6,*)('coupling matrix '),xint
                 call prntm2(bc(iprod),ncsf(jarg),nrootz,ncsf(jarg))
                 write(6,*)('subtract unit matrix ')
                end if
                do ir=0,nrootz-1                                        11d16s20
                 irp=ir+1                                               11d16s20
                 do i=0,ncsf(jarg)-1                                         11d13s20
                  ii=iprod+i+ncsf(jarg)*ir                              11d16s20
                  bc(ii)=bc(ii)-vecx(ip+i,irp)                          11d16s20
                 end do                                                 11d16s20
                end do                                                      11d13s20
                if(ltest.and.lnew)then
                 call prntm2(bc(iprod),ncsf(jarg),nrootz,ncsf(jarg))    11d16s20
                end if
                do i=0,ncsf(jarg)*nrootz-1                              11d16s20
                 bc(jhtmpgg+i)=bc(jhtmpgg+i)+xint*bc(iprod+i)               11d13s20
                end do                                                      11d13s20
                ibcoff=iprod                                                11d13s20
               else                                                        11d13s20
                if(ltest.and.lnew)then
                 write(6,*)
     $            ('karg is out of range, so just subtract unit matrix')
                 write(6,*)('nqq '),nqq,mdon,mdoo
                end if
                do ir=0,nrootz-1                                        11d16s20
                 irp=ir+1                                               11d16s20
                 do i=0,ncsf(jarg)-1                                        11d13s20
                  ii=jhtmpgg+i+ncsf(jarg)*ir                            11d16s20
                  bc(ii)=bc(ii)-xint*vecx(ip+i,irp)                     11d16s20
                 end do                                                 11d16s20
                end do                                                     11d13s20
               end if                                                      11d13s20
              end do                                                       11d13s20
             end do                                                        11d13s20
            else if(nnot4.eq.2)then                                        11d13s20
             call gandc(ibc(jjc),ibc(jjo),ibc(iic),ibc(iio),nopenj,     11d16s20
     $            nopen,jarg,iarg,ncsf,norb,ixw1,ixw2,nnot1,nab1,iwpb1, 11d16s20
     $            iwpk1,ncsfmid1,bc,ibc)                                11d14s22
             call gandc(ibc(iic),ibc(iio),ibc(jjc),ibc(jjo),nopen,      11d16s20
     $            nopenj,iarg,jarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,iwpb2, 11d16s20
     $            iwpk2,ncsfmid2,bc,ibc)                                11d14s22
             iprod=ibcoff                                                  11d13s20
             jprod=iprod+ncsf(iarg)*nrootz                              11d16s20
             ibcoff=jprod+ncsf(jarg)*nrootz                             12d4s20
             call enough('hccsf.  3',bc,ibc)
             call xtimesn(ncsf(jarg),ncsf(iarg),ncsfmid1,nrootz,iwpb1,  12d4s20
     $            iwpk1,vecx(ip,1),ncsft,bc(jprod),ncsf(jarg),1d0,0d0,  11d10s22
     $            bc,ibc)                                               11d10s22
             call xtimesn(ncsf(iarg),ncsf(jarg),ncsfmid2,nrootz,iwpb2,  12d4s20
     $            iwpk2,vecx(jp,1),ncsft,bc(iprod),ncsf(iarg),1d0,0d0,  11d10s22
     $            bc,ibc)                                               11d10s22
             do i=1,norb                                                   11d13s20
              itest(i,1)=0                                                   11d13s20
              itest(i,2)=0                                                    11d13s20
             end do                                                        11d13s20
             do i=0,nclo-1                                                 11d13s20
              itest(idorb(ic+i),1)=2                                        11d13s20
             end do                                                        11d13s20
             do i=0,nopen-1                                                 11d13s20
              itest(isorb(io+i),1)=1                                        11d13s20
             end do                                                        11d13s20
             do i=0,ncloj-1                                                 11d13s20
              itest(idorb(jc+i),2)=2                                         11d13s20
             end do                                                        11d13s20
             do i=0,nopenj-1                                                 11d13s20
              itest(isorb(jo+i),2)=1                                         11d13s20
             end do                                                        11d13s20
             if(ltest.and.lnew)write(6,*)('sames? ')
             nok=0                                                         11d13s20
             do i=1,norb                                                   11d13s20
              ixn=min(itest(i,1),itest(i,2))
              if(ltest.and.lnew)write(6,*)i,itest(i,1),itest(i,2),ixn                        11d13s2
              if(ixn.gt.0)then                                             11d13s20
               nok=nok+1                                                   11d13s20
               itest(nok,1)=ixn                                            11d13s20
               itest(nok,2)=i                                              11d13s20
              end if                                                       11d13s20
             end do                                                        11d13s20
             is=ism(nab4(1,1))                                             11d13s20
             jga=irel(nab4(2,1))-1                                      11d27s20
             jgb=irel(nab4(1,1))-1                                      11d27s20
             ih0u=ih0a(is)+jga+irefo(is)*jgb                               11d13s20
             sum=bc(ih0u)                                                  11d13s20
             if(ltest.and.lnew)write(6,*)('h0 '),sum
             do i=1,nok                                                    11d13s20
              js=ism(itest(i,2))                                           11d13s20
              jg=irel(itest(i,2))-1                                        11d13s20
              xint=getint(i2e,is,is,js,js,jga,jgb,jg,jg,bc,ibc)         11d15s22
              if(itest(i,1).eq.2)then                                      11d13s20
               xintk=getint(i2e,is,js,is,js,jga,jg,jgb,jg,bc,ibc)       11d15s22
               if(ltest.and.lnew)write(6,*)('jk '),xint,xintk
               sum=sum+2d0*xint-xintk                                      11d13s20
              else                                                         11d13s20
               if(ltest.and.lnew)write(6,*)('j '),xint
               sum=sum+xint                                                11d13s20
              end if                                                       11d13s20
             end do                                                        11d13s20
             if(ltest.and.lnew)then
              write(6,*)('sum '),sum
              write(6,*)('iprod ')
              call prntm2(bc(iprod),ncsf(iarg),nrootz,ncsf(iarg))
              write(6,*)('jprod ')
              call prntm2(bc(jprod),ncsf(jarg),nrootz,ncsf(jarg))
             end if
             do i=0,ncsf(iarg)*nrootz-1                                 11d16s20
              bc(ihtmpgg+i)=sum*bc(iprod+i)                                11d13s20
             end do                                                        11d13s20
             do j=0,ncsf(jarg)*nrootz-1                                 11d16s20
              bc(jhtmpgg+j)=sum*bc(jprod+j)                                11d13s20
             end do                                                        11d13s20
             ibcoff=iprod                                                11d13s20
             do i=1,nok                                                    11d13s20
              if(itest(i,1).eq.1)then                                      11d13s20
               itesta=ibc(jjc)                                              11d13s20
               itestb=ibc(jjo)                                              11d13s20
               nopenk=nopenj                                                11d13s20
c
c     anihilate common
c
               if(ltest.and.lnew)write(6,*)('anhihilate '),itest(i,2),
     $             loopit
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
               if(ltest.and.lnew)write(6,*)('create '),nab4(2,1)        11d27s20
               if(btest(itestb,nab4(2,1)))then                          11d27s20
                itesta=ibset(itesta,nab4(2,1))                          11d27s20
                itestb=ibclr(itestb,nab4(2,1))                          11d27s20
                karg=karg+1                                                11d13s20
                nopenk=nopenk-1                                             11d13s20
               else                                                         11d13s20
                itestb=ibset(itestb,nab4(2,1))                          11d27s20
                nopenk=nopenk+1                                             11d13s20
               end if                                                       11d13s20
               nqq=karg+mdon-1
               if(nqq.ge.mdon.and.nqq.le.mdoo)then                       11d13s20
                call gandc(ibc(jjc),ibc(jjo),itesta,itestb,nopenj,      11d16s20
     $              nopenk,jarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,    11d16s20
     $              iwpb1,iwpk1,ncsfmid1,bc,ibc)                        11d14s22
                call gandc(itesta,itestb,ibc(iic),ibc(iio),nopenk,nopen,      11d13s20
     $         karg,iarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,    11d13s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
                if(ltest.and.lnew)then
                 write(6,*)('nnot1, nab1 '),nnot1,nab1
                 write(6,*)('nnot2, nab2 '),nnot2,nab2
                 write(6,*)('jarg,karg,iarg '),jarg,karg,iarg
                 write(6,*)ncsf(jarg),ncsf(karg),ncsf(iarg),ncsfmid1,
     $              ncsfmid2
                 write(6,*)nopenj,nopenk,nopen
                 write(6,*)('iwpb1b '),iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2
                end if
                if(nnot1.eq.2.and.nnot2.eq.2)then                         11d13s20
                 jprod=ibcoff                                           11d16s20
                 ibcoff=jprod+ncsf(jarg)*nrootz                         12d4s20
                 call enough('hccsf.  4',bc,ibc)
                 call genmatn(ncsf(jarg),ncsf(karg),ncsf(iarg),ncsfmid1,12d4s20
     $                iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,vecx(ip,1),ncsft,12d4s20
     $                nrootz,bc(jprod),0d0,bc,ibc)                      11d10s22
                 jsa=ism(nab1(1))                                              11d13s20
                 jga=irel(nab1(1))-1                                           11d13s20
                 jsb=ism(nab1(2))                                              11d13s20
                 jgb=irel(nab1(2))-1                                           11d13s20
                 jsc=ism(nab2(1))                                              11d13s20
                 jgc=irel(nab2(1))-1                                           11d13s20
                 jsd=ism(nab2(2))                                              11d13s20
                 jgd=irel(nab2(2))-1                                           11d13s20
                 xint=getint(i2e,jsa,jsb,jsc,jsd,jga,jgb,jgc,jgd,bc,ibc)11d15s22
                 if(ltest.and.lnew)then
                  write(6,*)('coupling matrix '),xint
                  call prntm2(bc(jprod),ncsf(jarg),nrootz,ncsf(jarg))
                 end if
                 do j=0,ncsf(jarg)*nrootz-1                             11d16s20
                  bc(jhtmpgg+j)=bc(jhtmpgg+j)+xint*bc(jprod+j)                11d13s20
                 end do                                                       11d13s20
                 ibcoff=jprod                                                 11d13s20
                 call gandc(ibc(iic),ibc(iio),itesta,itestb,nopen,      11d16s20
     $              nopenk,iarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,    11d16s20
     $              iwpb1,iwpk1,ncsfmid1,bc,ibc)                        11d14s22
                 call gandc(itesta,itestb,ibc(jjc),ibc(jjo),nopenk,     11d16s20
     $                nopenj,karg,jarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,  11d16s20
     $                iwpb2,iwpk2,ncsfmid2,bc,ibc)                      11d14s22
                 iprod=ibcoff                                           11d16s20
                 ibcoff=iprod+ncsf(iarg)*nrootz                         12d4s20
                 call enough('hccsf.  5',bc,ibc)
                 call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),ncsfmid1,12d4s20
     $                iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,vecx(jp,1),ncsft,12d4s20
     $                nrootz,bc(iprod),0d0,bc,ibc)                      11d10s22
                 if(ltest.and.lnew)then
                  write(6,*)('coupling matrix '),xint
                  call prntm2(bc(iprod),ncsf(iarg),nrootz,ncsf(iarg))
                 end if
                 do j=0,ncsf(iarg)*nrootz-1                             11d16s20
                  bc(ihtmpgg+j)=bc(ihtmpgg+j)+xint*bc(iprod+j)          11d16s20
                 end do                                                       11d13s20
                 ibcoff=iprod                                           11d16s20
                else if(max(nnot1,nnot2).lt.2)then                        11d13s20
                 write(6,*)('c:expecting 2s, but got otherwise')
                 call dcbit(itesta,norb,'itesta')
                 call dcbit(itestb,norb,'itestb')
                 stop
                end if
               end if                                                    11d13s20
              end if                                                       11d13s20
             end do                                                        11d13s20
            else if(nnot4.ge.3)then                                        11d13s20
c     j is bra and i is ket
             if(nnot4.eq.3)then                                         1d22s21
              ipssx=1                                                   1d22s21
             else                                                       1d22s21
              ipssx=3                                                   1d22s21
             end if                                                     1d22s21
             fact=0d0                                                   1d22s21
             do ipss=1,ipssx                                            1d22s21
              if(ipss.eq.1)then                                         1d22s21
               iu1=1                                                    1d22s21
               iu2=1                                                    1d22s21
              else if(ipss.eq.2)then                                    1d22s21
               iu1=1                                                    1d22s21
               iu2=2                                                    1d22s21
              else                                                      1d22s21
               iu1=2                                                    1d22s21
               iu2=1                                                    1d22s21
              end if                                                    12d8s20
              itesta=ibc(jjc)                                               11d13s20
              itestb=ibc(jjo)                                               11d13s20
              if(btest(itesta,nab4(1,iu1)))then                         1d22s21
               if(ltest.and.lnew)
     $             write(6,*)('anihilate '),nab4(1,iu1)
               itesta=ibclr(itesta,nab4(1,iu1))                         1d22s21
               itestb=ibset(itestb,nab4(1,iu1))                         1d22s21
               nopenk=nopenj+1                                          1d22s21
               karg=jarg-1                                              1d22s21
              else if(btest(itestb,nab4(1,iu1)))then                    1d22s21
               if(ltest.and.lnew)
     $              write(6,*)('anihilate '),nab4(1,iu1)
               itestb=ibclr(itestb,nab4(1,iu1))                         1d22s21
               nopenk=nopenj-1                                              11d13s20
               karg=jarg                                                  11d13s20
              else                                                          11d13s20
               write(6,*)('bit not set for nab4(1,1) = '),nab4(1,iu1)      11d27s20
               stop 'nab4(1,1)'                                          11d27s20
              end if                                                        11d13s20
              if(btest(itestb,nab4(2,iu2)))then                         1d22s21
               if(ltest.and.lnew)
     $             write(6,*)('create '),nab4(2,iu2)
               itesta=ibset(itesta,nab4(2,iu2))                         1d22s21
               itestb=ibclr(itestb,nab4(2,iu2))                         1d22s21
               nopenk=nopenk-1                                              11d13s20
               karg=karg+1                                                  11d13s20
              else if(btest(itesta,nab4(2,iu2)))then                    1d22s21
               write(6,*)('already double in nab4(2,1) = '),nab4(2,iu2) 1d22s21
               stop 'nab4(2,1)'                                          11d27s20
              else                                                          11d13s20
               if(ltest.and.lnew)
     $              write(6,*)('create '),nab4(2,iu2)
               itestb=ibset(itestb,nab4(2,iu2))                         1d22s21
               nopenk=nopenk+1                                              11d13s20
              end if                                                        11d13s20
              nqq=karg+mdon-1                                           1d26s21
              if(nqq.ge.mdon.and.nqq.le.mdoo)then                       1d26s21
               call gandc(ibc(jjc),ibc(jjo),itesta,itestb,nopenj,nopenk,     11d13s20
     $         jarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,    11d13s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
               call gandc(itesta,itestb,ibc(iic),ibc(iio),nopenk,nopen,      11d13s20
     $         karg,iarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,iwpb2,iwpk2,    11d13s20
     $         ncsfmid2,bc,ibc)                                         11d14s22
               if(ltest.and.lnew)then
                write(6,*)('nnot1, nab1 '),nnot1,nab1
                write(6,*)('nnot2, nab2 '),nnot2,nab2
               end if
               if(nnot1.eq.2.and.nnot2.eq.2)then                         1d22s21
                jprod=ibcoff                                               11d16s20
                ibcoff=jprod+ncsf(jarg)*nrootz                             12d4s20
                call enough('hccsf.  6',bc,ibc)
                call genmatn(ncsf(jarg),ncsf(karg),ncsf(iarg),ncsfmid1,    12d4s20
     $            iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,vecx(ip,1),ncsft,    12d4s20
     $            nrootz,bc(jprod),0d0,bc,ibc)                          11d10s22
                jsa=ism(nab1(1))                                              11d13s20
                jga=irel(nab1(1))-1                                           11d13s20
                jsb=ism(nab1(2))                                              11d13s20
                jgb=irel(nab1(2))-1                                           11d13s20
                jsc=ism(nab2(1))                                              11d13s20
                jgc=irel(nab2(1))-1                                           11d13s20
                jsd=ism(nab2(2))                                              11d13s20
                jgd=irel(nab2(2))-1                                           11d13s20
                xint=getint(i2e,jsa,jsb,jsc,jsd,jga,jgb,jgc,jgd,bc,ibc) 11d15s22
                if(ltest.and.lnew)then
                 write(6,*)('coupling matrix '),xint
                 call prntm2(bc(jprod),ncsf(jarg),nrootz,ncsf(jarg))
                end if
                if(nab1(1).eq.nab2(1).and.nab1(2).eq.nab2(2))              11d16s20
     $              xint=xint*0.5d0                                       11d16s20
                do j=0,ncsf(jarg)*nrootz-1                                 11d16s20
                 bc(jhtmpgg+j)=bc(jhtmpgg+j)*fact+bc(jprod+j)*xint       1d22s21
                end do                                                        11d13s20
                ibcoff=jprod                                                  11d13s20
                call gandc(ibc(iic),ibc(iio),itesta,itestb,nopen,nopenk,   11d16s20
     $         iarg,karg,ncsf,norb,ixw1,ixw2,nnot1,nab1,iwpb1,iwpk1,    11d13s20
     $         ncsfmid1,bc,ibc)                                         11d14s22
                call gandc(itesta,itestb,ibc(jjc),ibc(jjo),nopenk,      1d26s21
     $               nopenj,karg,jarg,ncsf,norb,ixw1,ixw2,nnot2,nab2,   1d26s21
     $               iwpb2,iwpk2,ncsfmid2,bc,ibc)                       11d14s22
                iprod=ibcoff                                               11d16s20
                ibcoff=iprod+ncsf(iarg)*nrootz                             12d4s20
                call enough('hccsf.  7',bc,ibc)
                call genmatn(ncsf(iarg),ncsf(karg),ncsf(jarg),ncsfmid1,    12d4s20
     $            iwpb1,iwpk1,ncsfmid2,iwpb2,iwpk2,vecx(jp,1),ncsft,    12d4s20
     $            nrootz,bc(iprod),0d0,bc,ibc)                          11d10s22
                if(ltest.and.lnew)then
                 write(6,*)('coupling matrix '),xint
                 call prntm2(bc(iprod),ncsf(iarg),nrootz,ncsf(iarg))
                end if
                do i=0,ncsf(iarg)*nrootz-1                                 11d16s20
                 bc(ihtmpgg+i)=bc(ihtmpgg+i)*fact+bc(iprod+i)*xint       1d22s21
                end do                                                        11d13s20
                ibcoff=iprod                                                  11d13s20
                if(ipss.eq.2)go to 22                                    1d22s21
                fact=1d0                                                 1d22s21
               end if                                                   1d26s21
              end if                                                    1d22s21
             end do                                                     1d22s21
   22        continue                                                   1d22s21
            end if                                                         11d13s20
             itmp=ihtmpgg
             if(if.ne.jf)then                                           11d16s20
              iout=ibcoff                                               11d16s20
              ibcoff=iout+ncsf(iarg)*ncsf(jarg)                         11d16s20
             end if                                                     11d16s20
           end if                                                       11d16s20
            jht=jhtmpgg                                                 11d16s20
            do ir=0,nrootz-1                                            11d16s20
             irp=ir+1                                                   11d16s20
             do j=0,ncsf(jarg)-1                                        11d16s20
              gx(jp+j,irp)=gx(jp+j,irp)+bc(jht+j)                       11d16s20
             end do                                                     11d16s20
             jht=jht+ncsf(jarg)                                         11d16s20
            end do                                                      11d16s20
           if(jf.ne.if)then                                             12d9s19
             iht=ihtmpgg                                                 11d16s20
             do ir=0,nrootz-1                                            11d16s20
              irp=ir+1                                                   11d16s20
              do i=0,ncsf(iarg)-1                                        11d16s20
               gx(ip+i,irp)=gx(ip+i,irp)+bc(iht+i)                       11d16s20
              end do                                                     11d16s20
              iht=iht+ncsf(iarg)                                        11d16s20
             end do                                                      11d16s20
           end if                                                       12d9s19
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
      call dws_gsumf(gx,nwds)
      do i=1,nrootz                                                     7d12s19
       sum=0d0                                                          7d12s19
       do j=1,ncsft                                                     7d12s19
        sum=sum+gx(j,i)*vecx(j,i)                                       7d12s19
       end do                                                           7d12s19
       chc(i)=sum                                                       7d12s19
      end do                                                            7d12s19
      return
      end                                                               7d11s19
