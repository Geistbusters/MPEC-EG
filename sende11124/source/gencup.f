c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine gencup(i2sb,i2msb,i2sk,i2msk,nopenb,nopenk,nab,iwpb,   9d17s21
     $     iwpk,iout,imat,ntype,mcsf,vec,nvec,nroot,bc,ibc)             11d14s22
      implicit real*8 (a-h,o-z)                                         9d17s21
      logical ldebug,lcompb,lcompk                                      3d25s22
      integer*1 nab(2),ipackc1(8)                                       9d17s21
      integer*4 ipack48(2)                                              9d17s21
      integer*8 ipack8,ipackc                                           9d17s21
      equivalence (ipack8,ipack48),(ipackc,ipackc1)                     9d17s21
      dimension iwpb(4),iwpk(4),mcsf(2),vec(*)                          10d19s21
      include "common.store"                                            9d17s21
      data icall/0/                                                     9d17s21
      save icall                                                        9d17s21
      icall=icall+1                                                     9d17s21
      ldebug=icall.eq.-5453
      if(ldebug)then
       write(6,*)('addresses: '),ibc(iwpb(1)),ibc(iwpk(1))
       write(6,*)('iwpb '),iwpb
       write(6,*)('iwpk '),iwpk
       write(6,*)('nopenb,nopenk '),nopenb,nopenk
       write(6,*)('nab '),nab
       write(6,*)('bra spin: '),i2sb,i2msb                               9d17s21
       write(6,*)('ket spin: '),i2sk,i2msk                               9d17s21
      end if
      iso=mod(i2sb,2)                                                   9d17s21
      islo=max(iso,i2sb-2,i2sk-2)                                       9d17s21
      ishi=min(i2sb+2,i2sk+2)                                           9d17s21
      if(iwpb(4).eq.1)then                                              9d17s21
       ishi=min(ishi,nopenb)                                            9d17s21
      else                                                              9d17s21
       ishi=min(ishi,nopenb+2)                                          9d17s21
      end if                                                            9d17s21
      if(iwpk(4).eq.1)then                                              9d17s21
       ishi=min(ishi,nopenk)                                            9d17s21
      else                                                              9d17s21
       ishi=min(ishi,nopenk+2)                                          9d17s21
      end if                                                            9d17s21
      mxmt=4                                                            9d17s21
      iout=ibcoff                                                       9d17s21
      imat=iout+mxmt                                                    9d17s21
      ntype=0
      do is=islo,ishi,2                                                 9d17s21
       mslo=max(-is,i2msb-2,i2msk-2)                                    9d17s21
       mshi=min(is,i2msb+2,i2msk+2)                                     9d17s21
       do ms=mslo,mshi,2                                                9d17s21
        if(ldebug)then
         write(6,*)('for intermediate spin '),is,ms
         write(6,*)('for bra ')
        end if
        ff=0.5d0                                                        9d21s21
        isigb=1                                                         9d21s21
        if(iwpb(3).eq.0)then                                            9d17s21
         ioffb=(i2sb-iso)*(i2sb-iso+2)                                  9d21s21
         ioffb=ioffb/8                                                  9d21s21
         if(i2msb.lt.0)then                                             9d21s21
          ioffb=ioffb+((mod(i2sb,2)-i2msb)/2)                           9d21s21
          idelm=-(i2msb-ms)/2                                             9d17s21
          if(is.ne.i2sb)ff=-ff                                          9d21s21
          isigb=-1                                                      9d21s21
         else                                                           9d21s21
          ioffb=ioffb+((mod(i2sb,2)+i2msb)/2)                           9d21s21
          idelm=(i2msb-ms)/2                                             9d17s21
         end if                                                         9d21s21
         idels=(i2sb-is)/2                                              9d17s21
         idels=idels+1                                                  9d17s21
         idelm=idelm+1                                                  9d17s21
         ioff=ibc(iwpb(1))+idels+3*(idelm+3*ioffb)                      9d17s21
         if(ldebug)write(6,*)('normal '),ioffb,idels,idelm
        else                                                            9d17s21
         ioffb=(is-iso)*(is-iso+2)                                      9d21s21
         ioffb=ioffb/8                                                  9d21s21
         if(ms.lt.0)then                                                9d21s21
          if(is.ne.i2sb)ff=-ff                                          9d21s21
          ioffb=ioffb+((mod(is,2)-ms)/2)                                9d21s21
          idelm=-(ms-i2msb)/2                                           9d21s21
          isigb=-1                                                      9d21s21
         else                                                           9d21s21
          ioffb=ioffb+((mod(is,2)+ms)/2)                                9d21s21
          idelm=(ms-i2msb)/2                                             9d17s21
         end if                                                         9d21s21
         idels=(is-i2sb)/2                                              9d17s21
         idels=idels+1                                                  9d17s21
         idelm=idelm+1                                                  9d17s21
         ioff=ibc(iwpb(1))+idels+3*(idelm+3*ioffb)                      9d17s21
         if(ldebug)write(6,*)('transposed '),idels,idelm,ioffb
        end if                                                          9d17s21
        if(ldebug)write(6,*)('address of my data '),ioff,ibc(ioff)                9d17s21
        idata=ibc(ioff)                                                 9d17s21
        ipack8=ibc(idata)                                                9d17s21
        if(ldebug)write(6,*)('ipack8 at '),idata,('is '),ipack8,ipack48
        ncsfb=ipack48(1)                                                9d17s21
        ncsfk=ipack48(2)                                                9d17s21
        if(ldebug)write(6,*)('ncsfb,ncsfk '),ipack48                              9d17s21
        idata=idata+1                                                   9d17s21
        nxm=ibc(idata)                                                  3d25s22
        itoc=idata+1                                                    3d25s22
        idata=itoc+nxm                                                  3d25s22
        nt=ibc(idata)                                                   9d17s21
        if(ldebug)write(6,*)('number of matrices '),nt
        if(nt.lt.0.or.nt.gt.10)then
         write(6,*)('bad nt'),nt
         call dws_synca
         call dws_finalize
         stop
        end if
        nn=ipack48(1)*ipack48(2)*2+1                                    9d17s21
        imats=ibc(itoc+iwpb(2))                                         3d25s22
        if(ldebug)write(6,*)('start of data '),imats
        ipack8=ibc(imats)                                               9d17s21
        isigk=1                                                         9d21s21
        if(ldebug)write(6,*)('for ket ')
        if(iwpk(3).ne.0)then                                            9d17s21
         ioffk=(i2sk-iso)*(i2sk-iso+2)                                  9d21s21
         ioffk=ioffk/8                                                  9d21s21
         if(i2msk.lt.0)then                                             9d21s21
          ioffk=ioffk+((mod(i2sk,2)-i2msk)/2)                           9d21s21
          idelm=-(i2msk-ms)/2                                             9d17s21
          if(i2sk.ne.is)ff=-ff                                          9d21s21
          isigk=-1                                                      9d21s21
         else                                                           9d21s21
          ioffk=ioffk+((mod(i2sk,2)+i2msk)/2)                           9d21s21
          idelm=(i2msk-ms)/2                                             9d17s21
         end if                                                         9d21s21
         idels=(i2sk-is)/2                                              9d17s21
         idels=idels+1                                                  9d17s21
         idelm=idelm+1                                                  9d17s21
         joff=ibc(iwpk(1))+idels+3*(idelm+3*ioffk)                      9d17s21
         if(ldebug)write(6,*)('transposed '),ioffk,idels,idelm
        else                                                            9d17s21
         ioffk=(is-iso)*(is-iso+2)                                      9d21s21
         ioffk=ioffk/8                                                  9d21s21
         if(ms.lt.0)then                                                9d21s21
          ioffk=ioffk+((mod(is,2)-ms)/2)                                9d21s21
          idelm=-(ms-i2msk)/2                                             9d17s21
          if(is.ne.i2sk)ff=-ff                                          9d21s21
          isigk=-isigk                                                  9d21s21
         else                                                           9d21s21
          ioffk=ioffk+((mod(is,2)+ms)/2)                                9d21s21
          idelm=(ms-i2msk)/2                                             9d17s21
         end if                                                         9d21s21
         idels=(is-i2sk)/2                                              9d17s21
         idels=idels+1                                                  9d17s21
         idelm=idelm+1                                                  9d17s21
         joff=ibc(iwpk(1))+idels+3*(idelm+3*ioffk)                      9d17s21
         if(ldebug)write(6,*)('normal '),ioffk,idels,idelm
        end if                                                          9d17s21
        if(ldebug)write(6,*)('address of my data '),joff,ibc(joff)                9d17s21
        jdata=ibc(joff)                                                 9d17s21
        ipack8=ibc(jdata)                                                9d17s21
        if(ldebug)write(6,*)('ipack8 at '),jdata,('is '),ipack8,ipack48
        mcsfb=ipack48(1)                                                9d17s21
        mcsfk=ipack48(2)                                                9d17s21
        if(ldebug)write(6,*)('mcsfb,mcsfk '),ipack48                              9d17s21
        if(ntype.eq.0)then                                              10d19s21
         if(iwpb(3).eq.0)then                                           10d19s21
          n1=ncsfb                                                      10d19s21
          n2=ncsfk                                                      10d19s21
         else                                                           10d19s21
          n1=ncsfk                                                      10d19s21
          n2=ncsfb                                                      10d19s21
         end if                                                         10d19s21
         if(iwpk(3).eq.0)then                                           10d19s21
          n3=mcsfb                                                      10d19s21
          n4=mcsfk                                                      10d19s21
         else                                                           10d19s21
          n3=mcsfk                                                      10d19s21
          n4=mcsfb                                                      10d19s21
         end if                                                         10d19s21
         if(n2.ne.n3)then
          write(6,*)('yikes!! middle dims don''t match! '),n1,n2,n3,n4
          stop
         end if
         if(nroot.gt.0.and.ldebug)then                                             10d19s21
          write(6,*)('multiply ket by vectors ')
          call prntm2(vec,n4,nroot,nvec)                                10d19s21
         end if                                                         10d19s21
         if(nroot.gt.0)then                                             3d28s22
          ibcoff=imat+n1*nroot*2                                         10d19s21
         else                                                           3d28s22
          ibcoff=imat+n1*n4*2                                           3d28s22
         end if                                                         3d28s22
        end if                                                          10d19s21
        jdata=jdata+1                                                   9d17s21
        nxmj=ibc(jdata)                                                 3d25s22
        jtoc=jdata+1                                                    3d25s22
        jdata=jtoc+nxmj                                                 3d25s22
        mt=ibc(jdata)                                                   9d17s21
        if(ldebug)write(6,*)('number of matrices '),mt
        if(mt.lt.0.or.mt.gt.10)then
         write(6,*)('bad mt'),mt
         call dws_synca
         call dws_finalize
         stop
        end if
        mm=ipack48(1)*ipack48(2)*2+1                                    9d17s21
        jmats=ibc(jtoc+iwpk(2))                                         3d25s22
        if(ldebug)write(6,*)('start of data '),jmats,mm,iwpk(2)
        do jt=0,mt-1                                                    9d17s21
         ipack8=ibc(jmats)                                              9d17s21
         if(iabs(ipack48(1)).gt.iabs(ipack48(2)))then                   9d17s21
          ick1=ipack48(1)                                                9d17s21
          ick2=ipack48(2)                                                9d17s21
         else                                                           9d17s21
          ick2=ipack48(1)                                                9d17s21
          ick1=ipack48(2)                                                9d17s21
         end if                                                         9d17s21
         ick1=ick1*isigk                                                9d21s21
         ick2=ick2*isigk                                                9d21s21
         jmats=jmats+1                                                  3d25s22
         lcompb=ibc(jmats).ne.0                                         3d25s22
         if(lcompb)then                                                 3d25s22
          ipack8=ibc(jmats+1)                                           3d25s22
          nused=ipack48(1)                                              3d25s22
          noff=ipack48(2)                                               3d25s22
          nusedi=nused/2                                                3d25s22
          if(2*nusedi.ne.nused)nusedi=nusedi+1                          3d25s22
         end if                                                         3d25s22
         if(iwpk(3).eq.0)then                                           9d17s21
          nc=mcsfb                                                      9d17s21
          nd=mcsfk                                                      9d17s21
          if(lcompb)then                                                3d25s22
           icdm=jmats+2+noff                                            3d25s22
          else                                                          3d25s22
           icdm=jmats+1                                                  9d17s21
          end if                                                        3d25s22
         else                                                           9d17s21
          nc=mcsfk                                                      9d17s21
          nd=mcsfb                                                      9d17s21
          if(lcompb)then                                                3d25s22
           icdm=jmats+2                                                  9d17s21
          else                                                          3d25s22
           icdm=jmats+1+nc*nd                                            9d17s21
          end if                                                        3d25s22
         end if                                                         9d17s21
         ibctoper=ibcoff                                                3d28s22
         if(nroot.gt.0)then                                             10d19s21
          itmpv=ibcoff                                                  10d19s21
          ibcoff=itmpv+nc*nroot                                         10d19s21
          call enough('gencup.  1',bc,ibc)
          if(ldebug)then
           write(6,*)('multiply ')
           call prntm2(vec,nd,nroot,nvec)
          end if
          if(lcompb)then                                                3d25s22
           icmp1=icdm                                                    3d25s22
           icmp2=icmp1+nd                                               3d25s22
           icmp3=icmp2+nusedi                                            3d25s22
           call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),vec,nvec,    3d25s22
     $          nd,bc(itmpv),nc,nc,nroot,0d0,1d0)                       3d25s22
          else                                                          3d25s22
           call dgemm('n','n',nc,nroot,nd,1d0,bc(icdm),nc,vec,nvec,0d0,  10d19s21
     $         bc(itmpv),nc,                                            10d19s21
     d' gencup.  1')
          end if                                                        3d25s22
          if(ldebug)then
           write(6,*)('to obtain ')
           call prntm2(bc(itmpv),nc,nroot,nc)
          end if
         else if(lcompb)then                                            3d25s22
          itmpcmp=ibcoff                                                3d28s22
          ibcoff=itmpcmp+nc*nd                                          3d28s22
          call enough('gencup.  2',bc,ibc)
          icmp1=icdm                                                    3d25s22
          icmp2=icmp1+nd                                                3d25s22
          icmp3=icmp2+nusedi                                            3d25s22
          call uncompxu(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itmpcmp),nc, 3d28s22
     $         nd)                                                      3d28s22
         else                                                           3d28s22
          itmpcmp=icdm                                                  3d28s22
         end if                                                         10d19s21
         kmats=imats                                                    9d17s21
         do it=0,nt-1                                                   9d17s21
          ipack8=ibc(kmats)                                             9d17s21
          if(iabs(ipack48(1)).gt.iabs(ipack48(2)))then                   9d17s21
           icb1=ipack48(1)                                                9d17s21
           icb2=ipack48(2)                                                9d17s21
          else                                                           9d17s21
           icb2=ipack48(1)                                                9d17s21
           icb1=ipack48(2)                                                9d17s21
          end if                                                         9d17s21
          kmats=kmats+1                                                 3d25s22
          lcompk=ibc(kmats).ne.0                                        3d25s22
          if(lcompk)then                                                3d25s22
           ipack8=ibc(kmats+1)                                          3d25s22
           nusedk=ipack48(1)                                             3d25s22
           noffk=ipack48(2)                                              3d25s22
           nusedki=nusedk/2                                             3d25s22
           if(2*nusedki.ne.nusedk)nusedki=nusedki+1                     3d25s22
          end if                                                        3d25s22
          icb1=icb1*isigb                                               9d21s21
          icb2=icb2*isigb                                               9d21s21
          if(ldebug)write(6,*)('code '),icb1,icb2,ick1,ick2,
     $         ipack48,ipack8,kmats
          if(ick1*icb1.gt.0)then                                        9d17s21
           ipack48(1)=icb2                                              9d17s21
           ipack48(2)=ick2                                              9d17s21
           if(ldebug)write(6,*)('keep this one '),ipack8,
     $          (ibc(iout+kt),kt=0,ntype-1)
           mymat=-1                                                     9d17s21
           do kt=0,ntype-1                                              9d17s21
            if(ibc(iout+kt).eq.ipack8)then                              9d17s21
             mymat=kt                                                   9d17s21
            end if                                                      9d17s21
           end do                                                       9d17s21
           fact=1d0
           if(mymat.lt.0)then                                           9d17s21
            fact=0d0                                                    9d17s21
            ibc(iout+ntype)=ipack8                                      9d17s21
            mymat=ntype                                                 9d17s21
            ntype=ntype+1                                               9d17s21
            if(ntype.gt.mxmt)then                                       9d17s21
             write(6,*)('toooo many types in gencup! ')
             do kt=0,ntype-1
              ipack8=ibc(iout+kt)
              write(6,*)kt+1,ipack48
             end do
             stop 'gencup'
            end if                                                      9d17s21
           end if                                                       9d17s21
           if(iwpb(3).eq.0)then                                         9d17s21
            na=ncsfb                                                    9d17s21
            nb=ncsfk                                                    9d17s21
            if(lcompk)then                                              3d25s22
             iabm=kmats+2+noffk                                         3d25s22
            else                                                        3d25s22
             iabm=kmats+1                                                9d17s21
            end if                                                      3d25s22
           else                                                         9d17s21
            na=ncsfk                                                    9d17s21
            nb=ncsfb                                                    9d17s21
            if(lcompk)then                                              3d25s22
             iabm=kmats+2                                                9d17s21
            else                                                        3d25s22
             iabm=kmats+1+na*nb                                          9d17s21
            end if                                                      3d25s22
           end if                                                       9d17s21
           if(lcompk)then                                               3d25s22
            kcmp1=iabm                                                  3d25s22
            kcmp2=kcmp1+nb                                              3d25s22
            kcmp3=kcmp2+nusedki                                         3d25s22
           end if                                                       3d25s22
           if(nb.ne.nc)then                                             9d17s21
            write(6,*)('inter dimensions do not match!!! '),nb,nc
            stop 'gencup'
           end if                                                       9d17s21
           if(nroot.eq.0)then                                           10d19s21
            nsz=na*nd                                                    9d17s21
            imatu=imat+mymat*nsz                                         9d17s21
            if(ldebug)then
             write(6,*)('multiplyb '),fact,ff,imatu,kcmp1
             call prntm2(bc(itmpcmp),nb,nd,nb)                          3d28s22
            end if
            if(lcompk)then                                              3d25s22
             call compxtimes(ibc(kcmp1),ibc(kcmp2),bc(kcmp3),           3d28s22
     $           bc(itmpcmp),nb,nb,bc(imatu),na,na,nd,fact,ff)          3d28s22
             if(ldebug)then
              write(6,*)('itmpcmp after '),itmpcmp,imatu,imatu+na*nd
              call prntm2(bc(itmpcmp),nb,nd,nb)                          3d28s22
             end if
            else                                                        3d25s22
             if(ldebug)then                                             3d28s22
              write(6,*)('abm ')                                        3d28s22
              call prntm2(bc(iabm),na,nb,na)                            3d28s22
             end if                                                     3d28s22
             call dgemm('n','n',na,nd,nb,ff,bc(iabm),na,bc(itmpcmp),nb, 3d28s22
     $          fact,bc(imatu),na,                                      9d17s21
     d' gencup.  2')
             if(ldebug)then
              write(6,*)('itmpcmp after '),itmpcmp,imatu,imatu+na*nd
              call prntm2(bc(itmpcmp),nb,nd,nb)                          3d28s22
             end if
            end if                                                      3d25s22
            if(ldebug)then
             write(6,*)('to obtain ')
             call prntm2(bc(imatu),na,nd,na)
            end if
           else                                                         10d19s21
            nsz=na*nroot                                                10d19s21
            imatu=imat+mymat*nsz                                         9d17s21
            if(ldebug)then
             write(6,*)('multiplyc '),iwpb(3),fact,ff,imatu
             call prntm2(bc(itmpv),nb,nroot,nb)
            end if
            if(lcompk)then                                              3d25s22
             call compxtimes(ibc(kcmp1),ibc(kcmp2),bc(kcmp3),bc(itmpv),  3d25s22
     $           nb,nb,bc(imatu),na,na,nroot,fact,ff)                   3d25s22
            else                                                        3d25s22
             call dgemm('n','n',na,nroot,nb,ff,bc(iabm),na,bc(itmpv),nb, 10d19s21
     $          fact,bc(imatu),na,                                      9d17s21
     d' gencup.  3')
            end if                                                      3d25s22
            if(ldebug)then
             write(6,*)('to obtain ')
             call prntm2(bc(imatu),na,nroot,na)
            end if
           end if                                                       10d19s21
          end if                                                        9d17s21
          if(lcompk)then                                                3d25s22
           kmats=kmats+noffk*2-ncsfb+ncsfk+2                            3d25s22
          else                                                          3d25s22
           kmats=kmats+nn                                                9d17s21
          end if                                                        3d25s22
         end do                                                         9d17s21
         ibcoff=ibctoper                                                3d28s22
         if(lcompb)then                                                 3d25s22
          jmats=jmats+noff*2-mcsfb+mcsfk+2                              3d25s22
         else                                                           3d25s22
          jmats=jmats+mm                                                 9d17s21
         end if                                                         3d25s22
        end do                                                          9d17s21
       end do                                                           9d17s21
      end do                                                            9d17s21
      if(ldebug)write(6,*)('my matrices: '),nab,icall
      do i=3,8                                                          9d17s21
       ipackc1(i)=0                                                     9d17s21
      end do                                                            9d17s21
      if(nroot.ne.0)nd=nroot                                            10d19s21
      do i=0,ntype-1
       ipack8=ibc(iout+i)                                               9d17s21
       imatu=imat+i*nsz                                                 9d17s21
       if(ldebug)write(6,*)('type '),ipack48                            10d12s21
       if(ipack48(1).gt.0)then                                          9d17s21
        ipackc1(1)=nab(1)                                               9d17s21
       else                                                             9d17s21
        ipackc1(1)=-nab(1)                                              9d17s21
       end if                                                           9d17s21
       if(ipack48(2).gt.0)then                                          9d17s21
        ipackc1(2)=nab(2)                                               9d17s21
       else                                                             9d17s21
        ipackc1(2)=-nab(2)                                              9d17s21
       end if                                                           9d17s21
       if(ldebug)then
        write(6,*)('mapped to '),ipackc1(1),ipackc1(2),iwpb(3),iwpb(4),
     $      iwpk(3),iwpk(4)                                             9d17s21
        write(6,*)('address: '),imatu,imatu+na*nd
        call prntm2(bc(imatu),na,nd,na)
       end if                                                           10d12s21
       ibc(iout+i)=ipackc                                               9d17s21
      end do
      mcsf(1)=na                                                        10d12s21
      mcsf(2)=nd                                                        10d12s21
      ibcoff=imat+ntype*nsz                                             9d17s21
      return                                                            9d17s21
      end                                                               9d17s21
