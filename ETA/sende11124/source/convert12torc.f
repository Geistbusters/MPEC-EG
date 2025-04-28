c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine convert12torc(iff22,nff22,nfdat,vd,ncsf,ncsf2,isymmrci,8d13s21
     $     nvirt,nroot,mdon,mdoop,nsymb,multh,vx2,nsing,vnew,nff20,ism, 8d3s21
     $     ismv,irel,irelv,iff20,norb,irefo,kff20,veci,ncsfti,nff0,     8d3s21
     $     iff0,ff22,rsum,ihsdiag,nff1,iff1,bc,ibc)                     11d14s22
      implicit real*8 (a-h,o-z)                                         8d3s21
      logical lgot                                                      8d5s21
      integer*8 iff22(*),ipack8,ihsdiag(mdoop,*),i18,i28,i38                        8d13s21
      dimension ipack4(2),nff22(mdoop,2,*),nfdat(5,4,*),vd(*),ncsf(*),  8d3s21
     $     ncsf2(4,*),nvirt(*),multh(8,8),vx2(nroot),vnew(nsing,*),     8d3s21
     $     nff20(mdoop,3),ism(*),ismv(*),irel(*),irelv(*),irefo(*),     8d3s21
     $     iptb(200,8),kff20(*),veci(ncsfti,*),nff0(mdoop,3),iff0(*),   8d3s21
     $     nl(4),ff22(*),ioffl(4),rsum(nroot,nroot),nff1(mdoop,nsymb,*),8d13s21
     $     iff1(*)                                                      8d13s21
      equivalence (ipack8,ipack4)                                       8d3s21
      include "common.store"                                            8d3s21
      data loopx/2692410/
      data icall/0/
      save icall
      icall=icall+1
      i18=1                                                             8d16s21
      do ir=1,nroot
       vx2(ir)=0d0
      end do
      do ii=mdon+1,mdoop
       iarg=ii-mdon
       do isb=1,nsymb                                                   8d16s21
        isbv=multh(isb,isymmrci)
        if(min(nff1(ii,isb,1),nvirt(isbv)).gt.0)then
         ncol=ncsf(iarg)*nff1(ii,isb,1)
         nrow=nroot*nvirt(isbv)
         itmp=ibcoff
         ibcoff=itmp+ncol*nrow
         call enough('convert12torc.  1',bc,ibc)
         i28=nrow
         i38=ncol
         call ddi_get(bc,ibc,ihsdiag(ii,isb),i18,i28,i18,i38,bc(itmp))  11d15s22
         do icol=0,ncol-1
          do iv=0,nvirt(isbv)-1
           jtmp=itmp+nroot*iv+nrow*icol-1
           do ir=1,nroot
            if(abs(bc(jtmp+ir)).gt.vx2(ir))vx2(ir)=abs(bc(jtmp+ir))
           end do
          end do
         end do
         ibcoff=itmp
        end if
       end do                                                           8d16s21
      end do
      igoal=66897+1
      irgoal=2
      loop=0
      write(6,*)('setting ismv and irelv in convert12torc ')
      do i=1,norb
       ismv(i)=ism(i)
       irelv(i)=irel(i)
       write(6,*)('orb '),i,ism(i),irel(i)
      end do
      next=norb+1
      do isb=1,nsymb                                                    7d22s21
       if(nvirt(isb).gt.200)stop 'convert12torc iptb'
       do iv=1,nvirt(isb)                                               7d15s21
        ismv(next)=isb                                                  7d15s21
        irelv(next)=iv+irefo(isb)                                       7d15s21
        iptb(iv,isb)=next                                               7d15s21
        write(6,*)('orb '),next,isb,iv+irefo(isb)
        next=next+1
       end do
      end do
      ibcoff=iff20                                                      7d15s21
      ivoff=1                                                           7d15s21
      ifoff=0                                                           7d15s21
      jff20=1                                                           7d15s21
      do ii=mdon+1,mdoop                                                7d23s21
       nff20(ii,1)=0                                                    7d23s21
       nff20(ii,2)=-1
      end do                                                            7d23s21
      icntr=0
      do ii=mdon+1,mdoop                                                7d16s21
       write(6,*)('for nclo = '),ii-1,nff20(ii,1)
       iarg=ii-mdon
       write(6,*)('ncsf = '),ncsf(iarg)
       nff20(ii,1)=nff20(ii,1)+nff0(ii,1)                               7d23s21
       if(nff20(ii,2).lt.0)then
        nff20(ii,2)=jff20                                                7d16s21
        nff20(ii,3)=ivoff                                                7d16s21
       end if
       jf=nff0(ii,2)                                                    7d19s21
       jv=nff0(ii,3)                                                    7d19s21
       do if=1,nff0(ii,1)                                               7d19s21
        write(6,*)('for internal fcn '),if
        kff20(jff20)=iff0(jf)                                           7d19s21
        jff20=jff20+1                                                   7d19s21
        write(6,*)('new fcn '),jff20/2
        jf=jf+1                                                         7d19s21
        kff20(jff20)=iff0(jf)                                           7d19s21
        call dcbit(iff0(jf-1),norb,'c')
        call dcbit(iff0(jf),norb,'o')
        jff20=jff20+1                                                   7d19s21
        jf=jf+1                                                         7d19s21
        jv0=jv
        do i=1,ncsf(iarg)                                               7d19s21
         do ir=1,nroot                                                   7d19s21
          vnew(ivoff,ir)=veci(jv,ir)                                      7d19s21
         end do                                                         7d19s21
         ivoff=ivoff+1                                                  7d19s21
         jv=jv+1                                                        7d19s21
        end do                                                          7d19s21
       end do
       do isb=1,nsymb                                                   8d13s21
        isbv=multh(isb,isymmrci)                                        8d13s21
        if(min(nff1(ii,isb,1),nvirt(isbv)).gt.0)then
         ncol=ncsf(iarg)*nff1(ii,isb,1)
         nrow=nroot*nvirt(isbv)
         ibcoff=jff20+nff1(ii,isb,1)*nvirt(isbv)+iff20-1                12d23s21
         itmp=ibcoff
         ibcoff=itmp+nrow*ncol
         call enough('convert12torc.  2',bc,ibc)
         i28=nrow
         i38=ncol
         call ddi_get(bc,ibc,ihsdiag(ii,isb),i18,i28,i18,i38,bc(itmp))  11d15s22
         nff20(ii,1)=nff20(ii,1)+nff1(ii,isb,1)*nvirt(isbv)             7d16s21
         jff1=nff1(ii,isb,2)
         nkeep=0
         write(6,*)('for singles symmetry '),isb,isbv
         do if1=1,nff1(ii,isb,1)
          jff1p=jff1+1
          call dcbit(iff1(jff1),norb,'c')
          call dcbit(iff1(jff1p),norb,'o')
          do iv=1,nvirt(isbv)
           kff20(jff20)=iff1(jff1)
           jff20=jff20+1
           kff20(jff20)=ibset(iff1(jff1p),iptb(iv,isbv))                7d15s21
           write(6,*)('for fcn,v '),if1,iv,iptb(iv,isbv)
           write(6,*)('new fcn no. '),jff20/2
           call dcbit(kff20(jff20-1),31,'c')
           call dcbit(kff20(jff20),31,'o')
           jff20=jff20+1                                                7d15s21
           do i=0,ncsf(iarg)-1                                          7d15s21
            iad=itmp+nrow*(i+ncsf(iarg)*(if1-1))-1+nroot*(iv-1)         7d16s21
            do ir=1,nroot                                               7d15s21
             vnew(ivoff,ir)=bc(iad+ir)
            end do                                                      7d15s21
            ivoff=ivoff+1                                               7d15s21
           end do                                                       7d15s21
          end do                                                        7d15s21
          jff1=jff1p+1                                                  7d15s21
         end do
         call ddi_put(bc,ibc,ihsdiag(ii,isb),i18,i28,i18,i38,bc(itmp))  11d15s22
         ibcoff=itmp
        end if
       end do                                                           8d13s21
       do ipass=1,2
        ioffvd=1                                                          8d3s21
        write(6,*)('vd for pass '),ipass
        do isb=1,nsymb                                                    8d3s21
         isbv12=multh(isb,isymmrci)                                       8d3s21
         do isbv1=1,nsymb                                               8d3s21
          isbv2=multh(isbv1,isbv12)                                     8d3s21
          if(isbv2.ge.isbv1)then                                        8d3s21
           if(isbv1.eq.isbv2)then                                       8d3s21
            if(nff22(ii,1,isb).gt.0.and.ipass.eq.2)then                 8d3s21
             jvcv=nfdat(5,1,isb)+nff22(ii,2,isb)                        8d3s21
             do if=1,nff22(ii,1,isb)                                    8d3s21
              ipack8=iff22(jvcv)                                        8d3s21
              nspace=iff22(jvcv+1)                                      8d3s21
              nff20(ii+1,1)=nff20(ii+1,1)+nvirt(isbv1)                  8d3s21
              if(nff20(ii+1,2).lt.0)then                                8d3s21
               nff20(ii+1,2)=jff20                                      8d3s21
               nff20(ii+1,3)=ivoff                                      8d3s21
              end if                                                    8d3s21
              nl(1)=iff22(jvcv+2)
              iad1=jvcv+iff22(jvcv+6)                                   8d3s21
c     iad2 is ncsf2xnl
              iad2=iad1+nl(1)                                           8d3s21
              sz=0d0
              do iz=0,nvirt(isbv1)*nroot*nfdat(2,1,isb)-1
               sz=sz+vd(ioffvd+iz)**2
              end do
              sz=sqrt(sz/dfloat(nvirt(isbv1)*nroot*nfdat(2,1,isb)))
              ivoff0=ivoff
              do iv=1,nvirt(isbv1)                                      8d3s21
               write(6,*)('visv '),iv
               kff20(jff20)=ibset(ipack4(1),iptb(iv,isbv1))             8d3s21
               jff20=jff20+1                                            8d3s21
               kff20(jff20)=ipack4(2)                                   8d3s21
               write(6,*)('new fcn no. '),jff20/2
               call dcbit(kff20(jff20-1),31,'c')
               call dcbit(kff20(jff20),31,'o')
               jff20=jff20+1                                            8d3s21
               do ir=1,nroot                                            8d3s21
                do i=0,ncsf2(1,iarg)-1                                  8d3s21
                 vnew(ivoff+i,ir)=0d0                                   8d3s21
                end do                                                  8d3s21
                jad2=iad2                                               8d3s21
                do iii=0,nl(1)-1
                 iiii=iff22(iad1+iii)-1
                 iad=ioffvd+iv-1+nvirt(isbv1)*(ir-1+nroot*iiii)
                 do i=0,ncsf2(1,iarg)-1                                 8d3s21
                  vnew(ivoff+i,ir)=vnew(ivoff+i,ir)+ff22(jad2+i)        8d3s21
     $                 *vd(iad)                                         8d3s21
                 end do                                                 8d3s21
                 jad2=jad2+ncsf2(1,iarg)                                8d3s21
                end do                                                  8d3s21
               end do                                                   8d3s21
               ivoff=ivoff+ncsf2(1,iarg)                                8d3s21
              end do                                                    8d3s21
              nxnv=nvirt(isbv1)*ncsf2(1,iarg)
              jvcv=jvcv+nspace                                          8d3s21
             end do                                                     8d3s21
            end if                                                      8d3s21
            nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       8d3s21
            ioffvd=ioffvd+nroot*nvirt(isbv1)*nfdat(2,1,isb)             8d3s21
            isw=0                                                       8d3s21
           else                                                         8d3s21
            nvv=nvirt(isbv1)*nvirt(isbv2)                               8d3s21
            isw=1                                                       8d3s21
           end if                                                       8d3s21
           nft=0
           do l=1,4
            ioffl(l)=ioffvd                                             8d3s21
            if(nfdat(2,l,isb).gt.0)then                                 8d3s21
             ioffvd=ioffvd+nroot*nvv*nfdat(2,l,isb)
             nft=nft+nfdat(2,l,isb)
            end if                                                      8d3s21
           end do                                                       8d3s21
           if(nff22(ii,1,isb).gt.0.and.ipass.eq.1)then                  8d3s21
            nff20(ii,1)=nff20(ii,1)+nff22(ii,1,isb)*nvv                 8d3s21
            jvcv=nfdat(5,1,isb)+nff22(ii,2,isb)                         8d3s21
            do if=1,nff22(ii,1,isb)                                     8d3s21
             ipack8=iff22(jvcv)                                         8d3s21
             nspace=iff22(jvcv+1)                                       8d3s21
             ivoff0=ivoff
             do iv2=0,nvirt(isbv2)-1                                    8d3s21
              iv2p=iv2+1                                                8d3s21
              itop=(iv2+isw*(nvirt(isbv1)-iv2))-1                       8d3s21
              do iv1=0,itop                                             8d3s21
               write(6,*)('vnotv '),iv1,iv2
               iv1p=iv1+1
               irec=iv1+nvirt(isbv1)*iv2                                8d3s21
               itri=((iv2*(iv2-1))/2)+iv1                               8d3s21
               irow=itri+isw*(irec-itri)                                8d3s21
               kff20(jff20)=ipack4(1)                                   8d3s21
               jff20=jff20+1                                            8d3s21
               kff20(jff20)=ibset(ipack4(2),iptb(iv1p,isbv1))            8d3s21
               kff20(jff20)=ibset(kff20(jff20),iptb(iv2p,isbv2))         8d3s21
               write(6,*)('new fcn no. '),jff20/2
               call dcbit(kff20(jff20-1),31,'c')
               call dcbit(kff20(jff20),31,'o')
               jff20=jff20+1                                            8d3s21
               do ir=1,nroot                                            8d3s21
                do i=0,ncsf(iarg)-1                                     8d3s21
                 vnew(ivoff+i,ir)=0d0                                   8d3s21
                end do                                                  8d3s21
                jvoff=ivoff                                             8d3s21
                do l=1,4                                                 8d3s21
                 nl(l)=iff22(jvcv+1+l)                                   8d3s21
                 iad1=jvcv+iff22(jvcv+5+l)                              8d3s21
                 iad2=iad1+nl(l)                                         8d3s21
                 jad2=iad2                                              8d3s21
                 do iii=0,nl(l)-1                                          8d3s21
                  iiii=iff22(iad1+iii)-1                                8d3s21
                  iad=ioffl(l)+irow+nvv*(ir-1+nroot*iiii)               8d3s21
                  do i=0,ncsf2(l,iarg)-1                                8d3s21
                   vnew(jvoff+i,ir)=vnew(jvoff+i,ir)+ff22(jad2+i)       8d3s21
     $                    *vd(iad)                                        8d3s21
                  end do                                                8d3s21
                  jad2=jad2+ncsf2(l,iarg)                               8d3s21
                 end do                                                 8d3s21
                 jvoff=jvoff+ncsf2(l,iarg)                              8d3s21
                end do                                                  8d3s21
               end do                                                   8d3s21
               ivoff=ivoff+ncsf(iarg)                                   8d3s21
              end do                                                    8d3s21
             end do                                                     8d3s21
             nxnv=nvv*ncsf(iarg)
             jvcv=jvcv+nspace                                           8d3s21
            end do                                                      8d3s21
           end if                                                        8d3s21
          end if                                                        8d3s21
         end do                                                         8d3s21
        end do
       end do
      end do
      ibcoff=iff20+(jff20/2)
      return                                                            8d3s21
      end                                                               8d3s21
c mpec2.1 version zeta copyright u.s. government
      subroutine convertto(ncsf,nroot,mdon,mdoop,vx2,vnew,nff20,        5d18s23
     $     ism,ismv,irel,irelv,iff20,norb,kff20,veci,ncsfti,nff0,       5d18s23
     $     iff0,bc,ibc)                                                 5d18s23
      implicit real*8 (a-h,o-z)                                         8d3s21
      logical lgot                                                      8d5s21
      integer*8 ipack8,i18,i28,i38                                      5d18s23
      dimension ipack4(2),ncsf(*),vx2(nroot),vnew(ncsfti,*),            5d18s23
     $     nff20(mdoop,3),ism(*),ismv(*),irel(*),irelv(*),              5d18s23
     $     iptb(200,8),kff20(*),veci(ncsfti,*),nff0(mdoop,3),iff0(*),   8d3s21
     $     nl(4),ioffl(4)                                               5d18s23
      equivalence (ipack8,ipack4)                                       8d3s21
      include "common.store"                                            8d3s21
      data loopx/2692410/
      data icall/0/
      save icall
      icall=icall+1
      write(6,*)('starting vectors ')
      call prntm2(veci,ncsfti,nroot,ncsfti)
      i18=1                                                             8d16s21
      do ir=1,nroot
       vx2(ir)=0d0
      end do
      igoal=66897+1
      irgoal=2
      loop=0
      write(6,*)('setting ismv and irelv in convertto ')
      do i=1,norb
       ismv(i)=ism(i)
       irelv(i)=irel(i)
       write(6,*)('orb '),i,ism(i),irel(i)
      end do
      ibcoff=iff20                                                      7d15s21
      ivoff=1                                                           7d15s21
      ifoff=0                                                           7d15s21
      jff20=1                                                           7d15s21
      do ii=mdon+1,mdoop                                                7d23s21
       nff20(ii,1)=0                                                    7d23s21
       nff20(ii,2)=-1
      end do                                                            7d23s21
      do ii=mdon+1,mdoop                                                7d16s21
       write(6,*)('for nclo = '),ii-1,nff20(ii,1)
       iarg=ii-mdon
       write(6,*)('ncsf = '),ncsf(iarg)
       nff20(ii,1)=nff20(ii,1)+nff0(ii,1)                               7d23s21
       if(nff20(ii,2).lt.0)then
        nff20(ii,2)=jff20                                                7d16s21
        nff20(ii,3)=ivoff                                                7d16s21
       end if
       jf=nff0(ii,2)                                                    7d19s21
       jv=nff0(ii,3)                                                    7d19s21
       do if=1,nff0(ii,1)                                               7d19s21
        write(6,*)('for internal fcn '),if
        kff20(jff20)=iff0(jf)                                           7d19s21
        jff20=jff20+1                                                   7d19s21
        write(6,*)('new fcn '),jff20/2
        jf=jf+1                                                         7d19s21
        kff20(jff20)=iff0(jf)                                           7d19s21
        call dcbit(iff0(jf-1),norb,'c')
        call dcbit(iff0(jf),norb,'o')
        jff20=jff20+1                                                   7d19s21
        jf=jf+1                                                         7d19s21
        jv0=jv
        do i=1,ncsf(iarg)                                               7d19s21
         do ir=1,nroot                                                   7d19s21
          vnew(ivoff,ir)=veci(jv,ir)                                      7d19s21
         end do                                                         7d19s21
         ivoff=ivoff+1                                                  7d19s21
         jv=jv+1                                                        7d19s21
        end do                                                          7d19s21
       end do
      end do
      ibcoff=iff20+(jff20/2)
      write(6,*)('returning vectors ')
      call prntm2(vnew,ncsfti,nroot,ncsfti)
      return                                                            8d3s21
      end                                                               8d3s21
