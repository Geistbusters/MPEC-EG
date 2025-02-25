c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dotvall(g,ihsg,ihdg,nrootk,v,ihsv,ihdv,nrootb,ncsft,   7d22s21
     $     mff1,mff2,nff1,nff2,mdon,mdoop,isymmrci,multh,nsymb,ncsf,    7d22s21
     $     ncsf2,nvirt,result,iff1,iff2,norb,irefo,nff0,iff0,nfdat,     8d3s21
     $     vdb,gdb,mdoub,ndoub,vdk,n2e,isymop,i2eop,ixmt,phase2,sr2,    8d12s21
     $     idoubo,nbasdws,ioverwrite,icoder,idv4,ndv4,iwket,itu,ntype,  12d3s21
     $     iifmx,bc,ibc)                                                11d10s22
      implicit real*8 (a-h,o-z)                                         7d12s21
      integer*8 ihsg(mdoop,nsymb),ihsv(mdoop,nsymb),i18,i28,i38,i48,    8d21s21
     $     ihdg(mdoop,nsymb),ihdv(mdoop,nsymb)                          8d21s21
      dimension g(ncsft,*),v(ncsft,*),nff1(mdoop,nsymb,2),multh(8,8),     7d12s21
     $     ncsf(*),nvirt(*),result(nrootb,nrootk),iff1(*),irefo(*),     7d19s21
     $     nff0(mdoop,3),iff0(*),nff2(mdoop,nsymb,2),ncsf2(*),iff2(*),  8d3s21
     $     nfdat(5,4,*),vdb(*),gdb(*),vdk(*),isymop(*),i2eop(2,3),      8d12s21
     $     ixmt(8,*),idoubo(*),nbasdws(*)                               8d12s21
      include "common.store"                                            7d12s21
      data loopx/20/
      data icall/0/
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      save icall
      icall=icall+1
      loop=0
      do ik=1,nrootk                                                    7d12s21
       do ib=1,nrootb                                                   7d12s21
        result(ib,ik)=0d0                                               7d12s21
       end do                                                           7d12s21
      end do                                                            7d12s21
      if(mff1.gt.0)then                                                 7d12s21
       i18=1                                                            7d12s21
       do nclop=mdon+1,mdoop                                            7d19s21
        iarg=nclop-mdon                                                 7d19s21
        if(icall.eq.-1)then
        jf0=nff0(nclop,2)                                               7d19s21
        jv=nff0(nclop,3)                                                7d19s21
        do if=1,nff0(nclop,1)                                           7d19s21
         sz=0d0                                                         7d19s21
         do ir=1,nrootk                                                 7d19s21
          do i=0,ncsf(iarg)-1                                            7d19s21
           sz=sz+g(jv+i,ir)**2                                          7d19s21
          end do                                                        7d19s21
         end do                                                         7d19s21
         sz=sqrt(sz/dfloat(nrootk*ncsf(iarg)))
         jf0=jf0+2
         jv=jv+ncsf(iarg)                                               7d19s21
        end do                                                          7d19s21
        end if
        do isb=1,nsymb                                                   7d12s21
         isbv=multh(isb,isymmrci)                                        7d12s21
         if(min(nvirt(isbv),nff1(nclop,isb,1)).gt.0)then                  7d12s21
          ncol=nff1(nclop,isb,1)*ncsf(iarg)                               7d12s21
          call ilimts(1,ncol,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)  7d12s21
          nhere=ih+1-il                                                 7d12s21
          if(nhere.gt.0)then                                            7d12s21
           nrowb=nvirt(isbv)*nrootb                                      7d12s21
           nrowk=nvirt(isbv)*nrootk                                      7d12s21
           itmpv=ibcoff                                                 7d12s21
           itmpg=itmpv+nhere*nrowb                                      7d12s21
           ibcoff=itmpg+nhere*nrowk                                     7d12s21
           call enough('dotvall.  1',bc,ibc)
           i28=nrowk                                                    7d12s21
           i38=il                                                       7d12s21
           i48=ih                                                       7d12s21
           call ddi_get(bc,ibc,ihsg(nclop,isb),i18,i28,i38,i48,         11d15s22
     $          bc(itmpg))                                              11d15s22
           if(icoder.eq.10)then
            itmpxxx=ibcoff
            ibcoff=itmpxxx+ncol*nrowk
            call enough('dotvall.  2',bc,ibc)
            i38=1
            i48=ncol
            call ddi_get(bc,ibc,ihsg(nclop,isb),i18,i28,i38,i48,        11d15s22
     $           bc(itmpxxx))                                           11d15s22
            jff1=nff1(nclop,isb,2)                                       7d16s21
            ipert=ibcoff
            ibcoff=ipert+ncsf(iarg)*nrootk                               11d12s21
            call enough('dotvall.  3',bc,ibc)
            do if=1,nff1(nclop,isb,1)
             jtmpg=itmpxxx+nrowk*ncsf(iarg)*(if-1)
             write(6,*)('for fcn '),nclop,isb,if,jtmpg-itmpg
             call dcbit(iff1(jff1),norb,'closed')                       7d16s21
             call dcbit(iff1(jff1+1),norb,'open')                       7d16s21
             ihead=0
             do iv=1,nvirt(isbv)
              rms=0d0
              do i=0,ncsf(iarg)-1
               iad=jtmpg+nrootk*(iv-1+nvirt(isbv)*i)                     11d12s21
               do ir=0,nrootk-1                                          11d12s21
                rms=rms+bc(iad+ir)**2
                jpert=ipert+i+ncsf(iarg)*ir
                bc(jpert)=bc(iad+ir)
                if(.not.(abs(bc(iad+ir)).lt.1d-14))
     $              write(6,*)ir,i,bc(iad+ir),
     $              iad+ir-itmpxxx,nclop,isb
               end do
              end do
              rms=sqrt(rms/dfloat(nrootk*ncsf(iarg)))                    11d12s21
              if(.not.(rms.lt.1d-14))then
               write(6,*)('for iv = '),iv+irefo(isbv),isbv,iv
               iad=jtmpg+nrootk*(iv-1)                                   11d12s21
               call prntm2(bc(ipert),ncsf(iarg),nrootk,ncsf(iarg))       11d12s21
              end if
             end do
             jff1=jff1+2                                                 7d16s21
            end do
            ibcoff=itmpxxx
           end if
           i28=nrowb                                                    7d12s21
           call ddi_get(bc,ibc,ihsv(nclop,isb),i18,i28,i38,i48,         11d15s22
     $          bc(itmpv))                                              11d15s22
           if(icall.eq.-7)then
            sz=0d0
            do i=0,nrowk*nhere-1
             sz=sz+bc(itmpg+i)**2
            end do
            do i=0,nrowb*nhere-1
             sz=sz+bc(itmpv+i)**2
            end do
            sz=sqrt(sz/dfloat(nhere*(nrowk+nrowb)))
            if(sz.gt.1d-10)then
             write(6,*)('getting vectors from '),ihsv(nclop,isb)
             write(6,*)('getting hg from '),ihsg(nclop,isb)             8d21s21
             call prntm2(bc(itmpg),nrowk,nhere,nrowk)
             call prntm2(bc(itmpv),nrowb,nhere,nrowb)
            end if
           end if
           nn=nhere*nvirt(isbv)                                         7d12s21
           jtmpg=itmpg-1                                                7d12s21
           jtmpv=itmpv-1                                                7d12s21
           do i=0,nn-1                                                  7d12s21
            do ik=1,nrootk                                              7d12s21
             do ib=1,nrootb                                             7d12s21
              result(ib,ik)=result(ib,ik)+bc(jtmpv+ib)*bc(jtmpg+ik)     7d12s21
             end do                                                     7d12s21
            end do                                                      7d12s21
            jtmpg=jtmpg+nrootk                                          7d12s21
            jtmpv=jtmpv+nrootb                                          7d12s21
           end do                                                       7d12s21
           ibcoff=itmpv                                                 7d12s21
          end if                                                        7d12s21
         end if                                                         7d12s21
        end do                                                          7d12s21
       end do                                                           7d12s21
       result1=result(1,1)
       if(icoder.eq.10)then
        write(6,*)('result after singles: ')
        call prntm2(result,nrootb,nrootk,nrootb)
       end if
      end if                                                            7d12s21
      if(mff2.gt.0)then                                                 7d22s21
       i18=1                                                            7d12s21
       do nclop=mdon+1,mdoop                                            7d19s21
        iarg=nclop-mdon                                                 7d19s21
        if(icall.eq.-1)then
        jf0=nff0(nclop,2)                                               7d19s21
        jv=nff0(nclop,3)                                                7d19s21
        do if=1,nff0(nclop,1)                                           7d19s21
         sz=0d0                                                         7d19s21
         do ir=1,nrootk                                                 7d19s21
          do i=0,ncsf(iarg)-1                                            7d19s21
           sz=sz+g(jv+i,ir)**2                                          7d19s21
          end do                                                        7d19s21
         end do                                                         7d19s21
         sz=sqrt(sz/dfloat(nrootk*ncsf(iarg)))
         jf0=jf0+2
         jv=jv+ncsf(iarg)                                               7d19s21
        end do                                                          7d19s21
        end if
        do ipass=1,2                                                    8d20s21
         do isb=1,nsymb                                                   7d12s21
          isbv12=multh(isb,isymmrci)                                     7d22s21
          nvisv=0                                                        7d22s21
          nvnotv=0                                                       7d22s21
          do isbv1=1,nsymb                                               7d22s21
           isbv2=multh(isbv1,isbv12)                                     7d22s21
           if(isbv2.ge.isbv1)then                                        7d22s21
            if(isbv1.eq.isbv2)then                                       7d22s21
             nvisv=nvisv+nvirt(isbv1)                                    7d22s21
             nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       7d22s21
            else                                                         7d22s21
             nvv=nvirt(isbv1)*nvirt(isbv2)                               7d22s21
            end if                                                       7d22s21
            nvnotv=nvnotv+nvv                                            7d22s21
           end if                                                        7d22s21
          end do                                                         7d22s21
          if(min(nvisv+nvnotv,nff2(nclop,isb,1)).gt.0)then               7d22s21
           nvv=nvisv*ncsf2(iarg)+nvnotv*ncsf(iarg)                       7d22s21
           ncol=nff2(nclop,isb,1)                                        7d22s21
           call ilimts(1,ncol,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)  7d12s21
           nhere=ih+1-il                                                 7d12s21
           if(nhere.gt.0)then                                            7d12s21
            nrowb=nvv*nrootb                                             7d22s21
            nrowk=nvv*nrootk                                             7d22s21
            itmpv=ibcoff                                                 7d12s21
            itmpg=itmpv+nhere*nrowb                                      7d12s21
            ibcoff=itmpg+nhere*nrowk                                     7d12s21
            call enough('dotvall.  4',bc,ibc)
            i28=nrowk                                                    7d12s21
            i38=il                                                       7d12s21
            i48=ih                                                       7d12s21
            call ddi_get(bc,ibc,ihdg(nclop,isb),i18,i28,i38,i48,        11d15s22
     $           bc(itmpg))                                             11d15s22
            if(icall.eq.-3)then
             write(6,*)('for externals with isb = '),isb,
     $           ihdg(nclop,isb),i28,i38,i48                            8d21s21
             write(6,*)('ipass = '),ipass
             write(6,*)('hdg from '),ihdg(nclop,isb),nclop,isb,il
             call prntm2(bc(itmpg),nrowk,nhere,nrowk)
             iprt=ibcoff
             ibcoff=iprt+ncsf(iarg)*nrootk
             call enough('dotvall.  5',bc,ibc)
             jff=nff2(nclop,isb,2)
             do if2=il,ih
              jtmpg=itmpg+nrowk*(if2-il)
              jtmpg0=jtmpg                                              8d20s21
              write(6,*)nclop,isb,if2
              call dcbit(iff2(jff),norb,'closed')
              call dcbit(iff2(jff+1),norb,'opened')
              jff=jff+2
              do isbv1=1,nsymb                                           7d23s21
               isbv2=multh(isbv1,isbv12)                                 7d23s21
               if(isbv2.ge.isbv1)then                                    7d23s21
                if(isbv2.eq.isbv1)then                                   7d23s21
                 if(ipass.eq.2)then                                     8d20s21
                  do iv=0,nvirt(isbv1)-1                                  7d23s21
                   sz=0d0
                   do i=0,ncsf2(iarg)-1                                   7d23s21
                    do ir=0,nrootk-1
                     sz=sz+bc(jtmpg+ir)**2
                     jprt=iprt+i+ncsf2(iarg)*ir
                     bc(jprt)=bc(jtmpg+ir)
                     if(abs(bc(jtmpg+ir)).gt.1d-10)
     $                    write(6,*)bc(jtmpg+ir),
     $                  nclop,isb,jtmpg+ir-itmpg
                    end do
                    jtmpg=jtmpg+nrootk
                   end do                                                 7d23s21
                   sz=sqrt(sz/dfloat(nrootk*ncsf2(iarg)))
                   if(sz.gt.1d-10)then
                    write(6,*)('for visv '),iv,('s'),isbv1
                    call prntm2(bc(iprt),ncsf2(iarg),nrootk,ncsf2(iarg))
                   end if
                  end do                                                  7d23s21
                 else
                  jtmpg=jtmpg+nrootk*ncsf2(iarg)*nvirt(isbv1)           8d20s21
                 end if                                                 8d20s21
                isw=0                                                   7d23s21
               else                                                     7d23s21
                isw=1                                                   7d23s21
               end if                                                   7d23s21
               do iv2=0,nvirt(isbv2)-1                                  7d23s21
                itop=(iv2+isw*(nvirt(isbv1)-iv2))-1                     7d23s21
                do iv1=0,itop                                           7d23s21
                 itri=((iv2*(iv2-1))/2)+iv1                             7d23s21
                 irec=iv1+nvirt(isbv1)*iv2
                 icol=itri+isw*(irec-itri)
                 if(ipass.eq.1)then
                  sz=0d0
                  do i=0,ncsf(iarg)-1
                   do ir=0,nrootk-1
                    sz=sz+bc(jtmpg+ir)**2
                    if(abs(bc(jtmpg+ir)).gt.1d-10)
     $                   write(6,*)bc(jtmpg+ir),
     $                  nclop,isb,jtmpg+ir-itmpg,jtmpg+ir-jtmpg0
                     jprt=iprt+i+ncsf(iarg)*ir
                    bc(jprt)=bc(jtmpg+ir)
                   end do
                   jtmpg=jtmpg+nrootk
                  end do
                  sz=sqrt(sz/dfloat(nrootk*ncsf(iarg)))
                  if(sz.gt.1d-10)then
                   write(6,*)('for vnotv '),iv1,('s'),isbv1,(' and '),
     $                 iv2,('s'),isbv2
                   call prntm2(bc(iprt),ncsf(iarg),nrootk,ncsf(iarg))
                  end if
                 else                                                   8d20s21
                  jtmpg=jtmpg+nrootk*ncsf(iarg)                         8d20s21
                 end if                                                 8d20s21
                end do
               end do
              end if                                                    7d23s21
             end do                                                     7d23s21
            end do
           end if
           i28=nrowb                                                    7d12s21
           call ddi_get(bc,ibc,ihdv(nclop,isb),i18,i28,i38,i48,         11d15s22
     $          bc(itmpv))                                              11d15s22
           nn=nhere*nvv                                                 7d22s21
           jtmpg=itmpg-1                                                7d12s21
           jtmpv=itmpv-1                                                7d12s21
           if(ipass.eq.1)then                                           8d20s21
            if(icall.eq.-3)then
             write(6,*)('result b4 doubles: '),nvv
             call prntm2(result,nrootb,nrootk,nrootb)
             write(6,*)('getting hdv from '),ihdv(nclop,isb)
             call prntm2(bc(itmpv),nrowb,nhere,nrowb)
            end if
            do i=0,nn-1                                                  7d12s21
             do ik=1,nrootk                                              7d12s21
              do ib=1,nrootb                                             7d12s21
               prod=bc(jtmpv+ib)*bc(jtmpg+ik)
               result(ib,ik)=result(ib,ik)+bc(jtmpv+ib)*bc(jtmpg+ik)     7d12s21
              end do                                                     7d12s21
             end do                                                      7d12s21
             jtmpg=jtmpg+nrootk                                          7d12s21
             jtmpv=jtmpv+nrootb                                          7d12s21
            end do                                                       7d12s21
           end if                                                       8d20s21
           ibcoff=itmpv                                                 7d12s21
          end if                                                        7d12s21
         end if                                                         7d12s21
        end do                                                          7d12s21
       end do                                                           7d12s21
       end do                                                           8d20s21
      else if(mff2.lt.0)then                                            8d3s21
       nwds=mdoub*nrootk                                                8d6s21
       ibctoper=ibcoff                                                  11d21s22
       call dws_gsumf(gdb,nwds)                                         8d6s21
       if(icoder.eq.10)then
        write(6,*)('in non-orthogonal basis ')
        ioffg=1
        do isb=1,nsymb
         isbv12=multh(isb,isymmrci)
         do isbv1=1,nsymb
          isbv2=multh(isbv1,isbv12)
          if(isbv2.ge.isbv1)then
           if(isbv1.eq.isbv2)then
            if(nfdat(2,1,isb).gt.0)then
             nrow=nvirt(isbv1)*nrootk                                      8d9s21
             sz=0d0
             do is=0,nrow*nfdat(2,1,isb)-1
              sz=sz+gdb(ioffg+is)**2
             end do
             sz=sqrt(sz/dfloat(nrow*nfdat(2,1,isb)))
             if(.not.(sz.lt.1d-14).and.nrow*nfdat(2,1,isb).gt.0)then
              write(6,*)('visv for syms '),isb,isbv1,ioffg
              call prntm2(gdb(ioffg),nrow,nfdat(2,1,isb),nrow)
             end if
             ioffg=ioffg+nrow*nfdat(2,1,isb)
            end if
            nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2
           else
            nvv=nvirt(isbv1)*nvirt(isbv2)
           end if
           nrow=nvv*nrootk
           do l=1,4
            if(nfdat(2,l,isb).gt.0)then
             sz=0d0
             do is=0,nrow*nfdat(2,l,isb)-1
              sz=sz+gdb(ioffg+is)**2
             end do
             sz=sqrt(sz/dfloat(nrow*nfdat(2,l,isb)))
             if(.not.(sz.lt.1d-14).and.nrow*nfdat(2,l,isb).gt.0)then
              write(6,*)('vnotv for spins '),l,('syms '),isb,isbv1,
     $             isbv2,
     $            ioffg
              call prntm2(gdb(ioffg),nrow,nfdat(2,l,isb),nrow)
             end if
             ioffg=ioffg+nrow*nfdat(2,l,isb)
            end if
           end do
          end if
         end do
        end do
        write(6,*)('ioffg at end '),ioffg,ioffg/mdoub
       end if                                                           8d13s21
       itmpg=ibcoff                                                     8d12s21
       ibcoff=itmpg+ndoub*nrootk                                        8d3s21
       call enough('dotvall.  6',bc,ibc)
       call tofrob(bc(itmpg),gdb,nrootk,nfdat,nvirt,nsymb,multh,        8d3s21
     $      isymmrci,2,ndoub,mdoub,iff2,bc,ibc)                         11d10s22
       if(ioverwrite.ne.0)then                                          8d21s21
        jtmpg=itmpg-1                                                   8d21s21
        do ic=1,nrootk*ndoub                                            8d21s21
         gdb(ic)=bc(jtmpg+ic)                                           8d21s21
        end do                                                          8d21s21
       end if                                                           8d21s21
       if(n2e.gt.0)then                                                 8d12s21
        phasez=0d0
        nwds=ndoub*nrootk                                               8d12s21
        igtmp=ibcoff                                                    8d12s21
        ibcoff=igtmp+nwds                                               8d12s21
        call enough('dotvall.  7',bc,ibc)
        do iz=igtmp,ibcoff-1                                            8d12s21
         bc(iz)=0d0                                                     8d12s21
        end do                                                          8d12s21
        call do4v(bc(itmpg),bc(igtmp),nwds,vdk,nvirt,nsymb,multh,       8d13s21
     $       isymmrci,nfdat,nrootk,n2e,isymop,i2eop,ixmt,phase2,sr2,    8d13s21
     $       idoubo,irefo,nbasdws,icall.eq.-7,bc,ibc)                   11d14s22
        ibcoff=igtmp                                                    8d12s21
       else if(icoder.ne.0)then
       end if                                                           8d12s21
       ivoff=1                                                          8d12s21
       igoff=itmpg                                                      8d3s21
       idoit=0                                                          8d3s21
       do isb=1,nsymb                                                   8d3s21
        isbv12=multh(isb,isymmrci)                                       8d3s21
        do isbv1=1,nsymb                                                8d3s21
         isbv2=multh(isbv1,isbv12)                                      8d3s21
         if(isbv2.ge.isbv1)then                                         8d3s21
          if(isbv1.eq.isbv2)then                                        8d3s21
           do if=1,nfdat(3,1,isb)                                       8d3s21
            if(mod(idoit,mynprocg).eq.mynowprog)then                    8d3s21
             nhit=0
             do irk=1,nrootk                                            8d3s21
              iadk=igoff+nvirt(isbv1)*(irk-1)                           8d3s21
              do irb=1,nrootb                                           8d3s21
               iadb=ivoff+nvirt(isbv1)*(irb-1)                          8d3s21
               do iv=0,nvirt(isbv1)-1                                   8d3s21
                result(irb,irk)=result(irb,irk)+bc(iadk+iv)*vdb(iadb+iv)8d12s21
               end do                                                   8d3s21
              end do                                                    8d3s21
             end do                                                     8d3s21
            end if                                                      8d3s21
            idoit=idoit+1                                               8d3s21
            ivoff=ivoff+nvirt(isbv1)*nrootb                             8d3s21
            igoff=igoff+nvirt(isbv1)*nrootk                             8d3s21
           end do                                                       8d3s21
           nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                        8d3s21
           isw=0
          else                                                          8d3s21
           nvv=nvirt(isbv1)*nvirt(isbv2)                                8d3s21
           isw=1
          end if                                                        8d3s21
          do l=1,4                                                      8d3s21
           do if=1,nfdat(3,l,isb)                                       8d3s21
            if(mod(idoit,mynprocg).eq.mynowprog)then                    8d3s21
             do irk=1,nrootk                                            8d3s21
              iadk=igoff+nvv*(irk-1)                                    8d3s21
              do irb=1,nrootb                                           8d3s21
               iadb=ivoff+nvv*(irb-1)                                   8d3s21
               do iv=0,nvv-1                                            8d3s21
                result(irb,irk)=result(irb,irk)+bc(iadk+iv)*vdb(iadb+iv)8d12s21
               end do                                                   8d3s21
              end do                                                    8d3s21
             end do                                                     8d3s21
            end if                                                      8d3s21
            idoit=idoit+1                                               8d3s21
            ivoff=ivoff+nvv*nrootb                                      8d3s21
            igoff=igoff+nvv*nrootk                                      8d3s21
           end do                                                       8d3s21
          end do                                                        8d3s21
         end if                                                         8d3s21
        end do                                                          8d3s21
       end do                                                           8d3s21
       ibcoff=ibctoper                                                  11d21s22
       result2=result(1,1)
       if(icoder.eq.10)then
        write(6,*)('result after contracted doubles: ')
        call prntm2(result,nrootb,nrootk,nrootb)
       end if
      end if                                                            8d25s21
      if(icoder.eq.10)then
        rsum=0d0
        do ii=mdon+1,mdoop
         if(nff0(ii,1).gt.0)then
          iarg=ii-mdon
          ivoff=nff0(ii,3)
          jff=nff0(ii,2)
          do if=1,nff0(ii,1)
           sz=0d0
           do ir=1,nrootk
            do i=0,ncsf(iarg)-1
             sz=sz+g(ivoff+i,ir)**2
            end do
           end do
           sz=sqrt(sz/dfloat(nrootk*ncsf(iarg)))
           if(sz.gt.-1d-14)write(6,*)('ivoff: '),ivoff,ii
           if(sz.gt.-1d-14)call dcbit(iff0(jff),8,'closed')
           jff=jff+1
           if(sz.gt.-1d-14)call dcbit(iff0(jff),8,'open')
           jff=jff+1
           if(sz.gt.1d-14)then
            call prntm2(g(ivoff,1),ncsf(iarg),nrootk,ncsft)
            do i=0,ncsf(iarg)-1
             rsum=rsum+g(ivoff+i,1)**2
            end do
            write(6,*)('rsum so far: '),rsum
           end if
           ivoff=ivoff+ncsf(iarg)
          end do
         end if
        end do
       end if                                                           8d5s21
      call ilimts(ncsft,1,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)     7d12s21
      nhere=ih+1-il                                                     7d12s21
      if(nhere.gt.0)then                                                7d12s21
       do i=il,ih                                                       7d12s21
        do ik=1,nrootk                                                  7d12s21
         do ib=1,nrootb                                                 7d12s21
          result(ib,ik)=result(ib,ik)+v(i,ib)*g(i,ik)                   7d12s21
         end do                                                         7d12s21
        end do                                                          7d12s21
       end do                                                           7d12s21
      end if                                                            7d12s21
      nwds=nrootb*nrootk                                                7d12s21
      result3=result(1,1)
      call dws_gsumf(result,nwds)                                       7d12s21
      if(icall.eq.-1)then
       write(6,*)('result after global sum: ')
       call prntm2(result,nrootb,nrootk,nrootb)
      end if
      return                                                            7d12s21
      end                                                               7d12s21
