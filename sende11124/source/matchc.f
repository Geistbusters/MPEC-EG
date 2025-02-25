c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine matchc(inst,ias,ndet,nopen,isen,ndetj,jas,nopenj,vdet, 9d12s19
     $     vsp,vtp,nrootu,ismult,ncsf2,joffg,nct2c,icol,ndimv,iextra,   11d4s19
     $     vcsf,ndim,lcrenorm,vnorm,irooto,fact3jst,fact3jht,fact3jlt,  9d28s20
     $     fact3jst2,fact3jht2,srh,nrtot,bc,ibc)                        11d14s22
      implicit real*8 (a-h,o-z)
      integer*1 ias(nopen,*),jas(nopenj,ndetj),iorb(64)                 11d4s19
      logical ldebug                                                    12d12s19
      dimension inst(2),vdet(max(1,ndet),nrootu),vsp(ndetj,nrootu),     11d4s19
     $     vtp(ndetj,nrootu),ncsf2(4),joffg(4),nct2c(4),vcsf(ndim,1),   1d22s20
     $     vnorm(*),phstp(4),ivecl(4)                                   12d1s20
      include "common.store"                                            9d12s19
      include "common.print"                                            1d5s20
      data phstp/1d0,-1d0,1d0,-1d0/                                     9d29s20
      if(iprtr(9).eq.0)then                                             1d5s20
       ldebug=.false.                                                    12d12s19
      else                                                              1d5s20
       ldebug=.true.                                                    1d5s20
      end if                                                            1d5s20
      if(ldebug)then                                                    12d12s19
       write(6,*)('in matchc for scenario '),isen
       write(6,*)('iextra = '),iextra
       write(6,*)('icol = '),icol
       write(6,*)('inst: '),inst
       write(6,*)('ndet,ndetj '),ndet,ndetj
      end if                                                            12d12s19
      if(nopenj.eq.0)then                                               9d12s19
       if(ldebug)write(6,*)('there are no open shells ...')             12d12s19
       if(ismult.eq.1)then                                              9d12s19
        if(ldebug)                                                      12d12s19
     $    write(6,*)('straight copy to sp'),joffg(1),nct2c(1),vdet(1,1),12d12s19
     $      vcsf(1,1)                                                   12d12s19
        if(isen.eq.4)then                                               11d4s19
         factr=sqrt(2d0)                                                11d4s19
        else                                                            11d4s19
         factr=1d0                                                      11d4s19
        end if                                                          11d4s19
        do i=0,nrootu-1                                                   9d12s19
         ip=i+1                                                         9d12s19
         iadd=joffg(1)+nct2c(1)*(i+irooto+nrtot*icol)                   4d28s21
         bc(iadd)=bc(iadd)+vcsf(1,ip)*factr                             11d4s19
        end do                                                            9d12s19
       else if(ismult.eq.3)then                                         9d12s19
        if(ldebug)                                                      12d12s19
     $       write(6,*)('straight copy to ls tp '),joffg(2),nct2c(2)    12d12s19
        do i=0,nrootu-1                                                   9d12s19
         ip=i+1                                                         9d12s19
         iadd=joffg(2)+nct2c(2)*(i+irooto+nrtot*icol)                   4d28s21
         bc(iadd)=bc(iadd)+vcsf(1,ip)                                   11d4s19
        end do                                                            9d12s19
       end if                                                           9d12s19
       return                                                           9d12s19
      end if                                                            9d12s19
      ismultm=ismult-1                                                  9d12s19
      if(ndetj.gt.100 000)then                                          9d30s20
       write(6,*)('ndetj is out of line! ')
       call dws_sync
       call dws_finalize
       stop 'matchc'
      end if
      do i=1,nrootu                                                     9d12s19
       do j=1,ndetj                                                     9d12s19
        vsp(j,i)=0d0                                                    9d12s19
        vtp(j,i)=0d0                                                    9d12s19
       end do                                                           9d12s19
      end do                                                            9d12s19
      do idet=1,max(1,ndet)                                             11d4s19
       if(ldebug)write(6,*)('for det '),idet,(ias(j,idet),j=1,nopen)    12d12s19
       do iab=1,2                                                       9d10s19
        if(iab.eq.1)then                                                9d10s19
         if(ldebug)write(6,*)('for alpha-beta excitation ')             12d12s19
        else                                                            9d10s19
         if(ldebug)write(6,*)('for beta-alpha excitation ')             12d12s19
        end if                                                          9d10s19
        if(isen.eq.2)then
         if(ldebug)write(6,*)('nopen,inst: '),nopen,inst                12d12s19
         nb=inst(2)-inst(1)+1                                           11d18s19
         jj=0                                                           9d10s19
         do i=1,min(inst(1)-1,nopen)                                    9d10s19
          jj=jj+1                                                       9d10s19
          iorb(jj)=ias(i,idet)                                          9d10s19
         end do                                                         9d10s19
         if(ldebug)write(6,*)('inserting first ')                       12d12s19
         jj=jj+1                                                        9d10s19
         if(iab.eq.1)then                                               9d10s19
          iorb(jj)=-1                                                   9d10s19
         else                                                           9d10s19
          iorb(jj)=+1                                                   9d10s19
         end if                                                         9d10s19
         do i=inst(1),min(inst(2)-1,nopen)                              9d10s19
          jj=jj+1                                                       9d10s19
          iorb(jj)=ias(i,idet)                                          9d10s19
         end do                                                         9d10s19
         jj=jj+1                                                        9d10s19
         if(ldebug)write(6,*)('inserting second')                       12d12s19
         if(iab.eq.1)then                                               9d10s19
          iorb(jj)=+1                                                   9d10s19
         else                                                           9d10s19
          iorb(jj)=-1                                                   9d10s19
         end if                                                         9d10s19
         do i=inst(2),nopen                                             9d10s19
          jj=jj+1                                                       9d10s19
          iorb(jj)=ias(i,idet)                                          9d10s19
         end do                                                         9d10s19
        else if(isen.eq.3)then                                           9d10s19
         nb=2*nopen-inst(1)-inst(2)+1                                   11d18s19
         if(inst(1).lt.inst(2))then                                     9d10s19
          nb=nb+1                                                       9d16s19
          if(mod(iab+iextra,2).eq.1)then                                  9d16s19
           igoal=1                                                      9d10s19
          else                                                          9d10s19
           igoal=-1                                                      9d10s19
           nb=nb+1                                                      9d13s19
          end if                                                        9d10s19
         else                                                           9d10s19
          if(mod(iab+iextra,2).eq.1)then                                  9d16s19
           igoal=+1                                                     9d16s19
          else                                                          9d10s19
           igoal=-1                                                     9d16s19
           nb=nb+1                                                      9d13s19
          end if                                                        9d10s19
         end if                                                         9d10s19
         if(ldebug)write(6,*)('igoal '),igoal                           12d12s19
         if(ias(inst(1),idet).ne.igoal)go to 2                          9d10s19
         jj=0                                                           9d10s19
         do i=1,min(inst(2)-1,nopen)                                    9d10s19
          if(i.ne.inst(1))then                                          9d10s19
           jj=jj+1                                                      9d10s19
           iorb(jj)=ias(i,idet)                                         9d10s19
          end if                                                        9d10s19
         end do                                                         9d10s19
         jj=jj+1                                                        9d10s19
         iorb(jj)=igoal                                                 9d10s19
         ihit=0                                                         9d10s19
         do i=inst(2),nopen                                             9d10s19
          if(i.ne.inst(1))then                                          9d10s19
           ihit=1                                                       9d10s19
           jj=jj+1                                                      9d10s19
           iorb(jj)=ias(i,idet)                                         9d10s19
          end if                                                        9d10s19
         end do                                                         9d10s19
         if(inst(2).gt.nopen.and.ihit.eq.1)then                         9d10s19
          jj=jj+1                                                        9d10s19
          iorb(jj)=igoal                                                9d10s19
         end if                                                         9d10s19
        else                                                             9d10s19
         if(iab.eq.1.and.ias(inst(1),idet).lt.0.or.                     9d10s19
     $      iab.eq.1.and.ias(inst(2),idet).gt.0.or.                     9d10s19
     $      iab.eq.2.and.ias(inst(1),idet).gt.0.or.                     9d10s19
     $      iab.eq.2.and.ias(inst(2),idet).lt.0)go to 2                 9d10s19
         jj=0                                                           9d10s19
         do i=1,nopen                                                   9d10s19
          if(i.ne.inst(1).and.i.ne.inst(2))then                         9d10s19
           jj=jj+1                                                      9d10s19
           iorb(jj)=ias(i,idet)                                         9d10s19
          end if                                                        9d10s19
         end do                                                         9d10s19
         nb=nopen-inst(2)+nopen-inst(1)                                 11d18s19
        end if                                                           9d10s19
        if(ldebug)then                                                  12d12s19
         write(6,*)('nb for phase: '),nb
         write(6,*)('promoted det: '),(iorb(ij),ij=1,jj)                 9d10s19
        end if                                                          12d12s19
        if(jj.ne.nopenj)then                                            9d10s19
         write(6,*)('jj = '),jj,(' ne nopenj = '),nopenj
         call dws_sync
         call dws_finalize
         stop 'matchc'                                                  9d10s19
        end if                                                          9d10s19
        do i=1,ndetj                                                    9d10s19
         do j=1,nopenj                                                  9d10s19
          if(iorb(j).ne.jas(j,i))go to 1                                9d10s19
         end do                                                         9d10s19
         if(ldebug)write(6,*)('matches det no. '),i                     12d12s19
         if(mod(nb,2).eq.0)then                                         9d12s19
          ff=1d0                                                        9d12s19
         else                                                           9d12s19
          ff=-1d0                                                       9d12s19
         end if                                                         9d12s19
         if(mod(iab+iextra,2).eq.1)then                                 9d16s19
          do j=1,nrootu                                                 9d12s19
           prod=ff*vdet(idet,j)                                         9d12s19
           vsp(i,j)=vsp(i,j)+prod                                       9d12s19
          end do                                                        9d12s19
         else                                                           9d12s19
          do j=1,nrootu                                                 9d12s19
           prod=ff*vdet(idet,j)                                         9d12s19
           vtp(i,j)=vtp(i,j)+prod                                       9d12s19
          end do                                                        9d12s19
         end if                                                         9d12s19
         go to 2                                                        9d10s19
    1    continue                                                       9d10s19
        end do                                                          9d10s19
        write(6,*)('no match found!!! ')                                9d10s19
        call dws_sync                                                   9d10s19
        call dws_finalize                                               9d10s19
        stop 'matchc'                                                   9d10s19
    2   continue                                                        9d10s19
       end do                                                           9d10s19
      end do
      do i=1,nrootu                                                     9d12s19
       do j=1,ndetj                                                     9d12s19
        sum=vsp(j,i)+vtp(j,i)                                           9d12s19
        dif=vsp(j,i)-vtp(j,i)                                           9d12s19
        vsp(j,i)=dif                                                    9d12s19
c                                                                       9d29s20
c     I'm getting wrong phase for scenario's 2 and 4                     9d29s20
c     but this doesn't really matter, because we actually use fcns      9d29s20
c     generated using matchcl, which has correct phase.                 9d29s20
c
        vtp(j,i)=sum*phstp(isen)                                        9d29s20
       end do                                                           9d12s19
      end do                                                            9d12s19
      if(ldebug)then                                                    12d12s19
       write(6,*)('sp dets ')                                            9d12s19
       call prntm2(vsp,ndetj,nrootu,ndetj)                               9d12s19
       write(6,*)('tp dets ')                                            9d12s19
       call prntm2(vtp,ndetj,nrootu,ndetj)                               9d12s19
       write(6,*)('ncsf2: '),ncsf2
      end if                                                            12d12s19
      iaorb=ibcoff                                                      9d12s19
      ivtmp=iaorb+ndetj*nopenj                                          3d13s20
      ibcoff=ivtmp+1                                                    12d1s20
      if(ldebug)write(6,*)('for same spin vectors ')                    12d12s19
      call fetchcsf2(nopenj,ismultm,nfcn,bc(ivtmp),idum,ibc(iaorb),
     $     ismultm,idumxx,bc,ibc)                                       11d15s22
      if(ldebug)write(6,*)('nfcn,ndetj,idumxx '),nfcn,ndetj,idumxx
      itmpt=ibcoff                                                      9d12s19
      itmpv=itmpt+ndetj*nfcn                                            9d12s19
      ibcoff=itmpv+nrootu*nfcn                                          9d12s19
      ivtmpa=nint(bc(ivtmp))                                            12d1s20
      call enough('matchc.  1',bc,ibc)
      do i=0,nfcn-1                                                     9d12s19
       do j=0,ndetj-1                                                   9d12s19
        ji=ivtmpa+j+ndetj*i                                             12d1s20
        ij=itmpt+i+nfcn*j                                               9d12s19
        bc(ij)=bc(ji)                                                   9d12s19
       end do                                                           9d12s19
      end do                                                            9d12s19
      if(ldebug.and.abs(fact3jst).gt.1d-10.and.min(nfcn,ndetj).gt.0)then10d2s20
       ff=srh*fact3jst
       call dgemm('n','n',nfcn,nrootu,ndetj,ff,bc(itmpt),nfcn,vtp,      9d28s20
     $      ndetj,0d0,bc(itmpv),nfcn,
     d' matchc.  1')
       write(6,*)('same spin tp ')
       call prntm2(bc(itmpv),nfcn,nrootu,nfcn)
      end if                                                            12d12s19
      if(min(nfcn,ndetj).gt.0)then                                      10d2s20
       call dgemm('n','n',nfcn,nrootu,ndetj,1d0,bc(itmpt),nfcn,vsp,     10d8s20
     $     ndetj,0d0,bc(itmpv),nfcn,                                    10d2s20
     d' matchc.  2')
       if(ldebug)then                                                    12d12s19
        write(6,*)('contribution to sp vector ')
        call prntm2(bc(itmpv),nfcn,nrootu,nfcn)
       end if                                                            12d12s19
       do i=0,nrootu-1                                                   9d12s19
        iadd=joffg(1)+nct2c(1)*(i+irooto+nrtot*icol)                    4d28s21
        jadd=itmpv+nfcn*i                                                9d12s19
        do j=0,nfcn-1                                                    9d12s19
         bc(iadd+j)=bc(iadd+j)+bc(jadd+j)                                9d12s19
        end do                                                           9d12s19
       end do                                                            9d12s19
      end if                                                            10d2s20
      if(ldebug)then
       write(6,*)('for high spin vectors ')
       ismultxx=ismult+1                                                  9d28s20
       call fetchcsf2(nopenj,ismultxx,nfcnxx,bc(ivtmp),idum,ibc(iaorb),
     $     ismultm,idum,bc,ibc)                                         11d15s22
       ivtmpa=nint(bc(ivtmp))                                           12d1s20
       call prntm2(bc(ivtmpa),ndetj,nfcnxx,ndetj)                       12d1s20
       if(min(nfcnxx,ndetj).gt.0)then                                   9d29s20
        ff=srh*fact3jht                                                  9d28s20
        call dgemm('t','n',nfcnxx,nrootu,ndetj,ff,bc(ivtmpa),ndetj,vtp, 12d1s20
     $      ndetj,0d0,bc(itmpv),nfcnxx,
     d' matchc.  3')
        write(6,*)('high spin spin tp ')
        call prntm2(bc(itmpv),nfcnxx,nrootu,nfcnxx)
       end if                                                            12d12s19
      end if
      ibcoff=ivtmp                                                      9d12s19
 1676 continue                                                          11d18s19
      nfcnt=0                                                           9d12s19
      ivec=ibcoff                                                       9d12s19
      ismm=ismultm-2                                                    9d12s19
      do l=2,4                                                          9d12s19
       if(ncsf2(l).gt.0)then                                            9d12s19
        if(ldebug)write(6,*)('load vectors for tcp spin '),l-1          12d12s19
        ivecu=ibcoff                                                    9d12s19
        isx=ismultm+2*(l-3)                                             9d12s19
        call fetchcsf2(nopenj,isx,nfcn,bc(ivecu),idum,ibc(iaorb),        9d12s19
     $      ismm,ndetk,bc,ibc)                                          11d15s22
        ivecl(l)=nint(bc(ivecu))                                        12d1s20
        if(nfcn.ne.ncsf2(l))then                                        9d12s19
         write(6,*)('nfcn = '),nfcn,(' doesn''t match ncsf2 = '),
     $       ncsf2(l)
         call dws_sync
         call dws_finalize
         stop
        end if                                                          9d12s19
        nfcnt=nfcnt+nfcn                                                9d12s19
       end if                                                           9d12s19
      end do                                                            9d12s19
      if(nfcnt.gt.0)then                                                9d12s19
       itmpt=ibcoff                                                     9d12s19
       itmp=itmpt+ndetk*nfcnt                                           9d13s19
       iout=itmp+ndetk*nrootu                                           9d13s19
       ibcoff=iout+nfcnt*nrootu                                         9d13s19
       call enough('matchc.  2',bc,ibc)
       loffv=0                                                          12d1s20
       do l=2,4                                                         12d1s20
        if(ncsf2(l).gt.0)then                                           12d1s20
         do i=0,ncsf2(l)-1                                              12d1s20
          ip=i+loffv                                                    12d1s20
          do j=0,ndetk-1                                                12d1s20
           ji=ivecl(l)+j+ndetk*i                                        12d1s20
           ij=itmpt+ip+nfcnt*j                                          12d1s20
           bc(ij)=bc(ji)                                                12d1s20
          end do                                                        12d1s20
         end do                                                         12d1s20
         loffv=loffv+ncsf2(l)                                           12d1s20
        end if                                                          12d1s20
       end do                                                           12d1s20
       call matchcl(inst,ias,ndet,nopen,isen,ndetk,ibc(iaorb),nopenj,   9d12s19
     $      vdet,bc(itmp),bc(itmpt),nrootu,nfcnt,bc(iout),iextra,       1d22s20
     $      lcrenorm,vnorm,fact3jlt,fact3jst2,fact3jht2)                9d28s20
       jout=iout                                                        9d13s19
       do l=2,4                                                         9d13s19
        if(l.eq.2)then                                                  9d28s20
         ff=fact3jlt                                                    9d28s20
        else if(l.eq.3)then                                             9d28s20
         ff=fact3jst2                                                   9d28s20
        else                                                            9d28s20
         ff=fact3jht2                                                   9d28s20
        end if                                                          9d28s20
        do i=0,nrootu-1                                                   9d12s19
         iadd=joffg(l)+nct2c(l)*(i+irooto+nrtot*icol)                   4d28s21
         jadd=jout+nfcnt*i                                              9d13s19
         do j=0,ncsf2(l)-1                                              9d13s19
          xx=bc(jadd+j)*ff                                              9d28s20
          orig=bc(iadd+j)
          bc(iadd+j)=bc(iadd+j)+xx                                      9d28s20
          if(ldebug)write(6,*)l,j,i,xx,orig,bc(iadd+j)
         end do                                                           9d12s19
        end do                                                            9d12s19
        jout=jout+ncsf2(l)                                              9d13s19
       end do                                                           9d13s19
      end if                                                            9d13s19
ccccc
ccc
      ibcoff=iaorb                                                      9d12s19
      return
      end
