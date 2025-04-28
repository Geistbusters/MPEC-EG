c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine convert2torc(iff22,nff22,nfdat,vd,ncsf,ncsf2,isymmrci, 8d3s21
     $     nvirt,nroot,mdon,mdoop,nsymb,multh,vx2,nsing,vnew,nff20,ism, 8d3s21
     $     ismv,irel,irelv,iff20,norb,irefo,kff20,veci,ncsfti,nff0,     8d3s21
     $     iff0,ff22,rsum,bc,ibc)                                       11d14s22
      implicit real*8 (a-h,o-z)                                         8d3s21
      logical lgot                                                      8d5s21
      integer*8 iff22(*),ipack8                                         8d3s21
      dimension ipack4(2),nff22(mdoop,2,*),nfdat(5,4,*),vd(*),ncsf(*),  8d3s21
     $     ncsf2(4,*),nvirt(*),multh(8,8),vx2(nroot),vnew(nsing,*),     8d3s21
     $     nff20(mdoop,3),ism(*),ismv(*),irel(*),irelv(*),irefo(*),     8d3s21
     $     iptb(200,8),kff20(*),veci(ncsfti,*),nff0(mdoop,3),iff0(*),   8d3s21
     $     nl(4),ff22(*),ioffl(4),rsum(nroot,nroot)                                       8d3s21
      equivalence (ipack8,ipack4)                                       8d3s21
      include "common.store"                                            8d3s21
      data loopx/2692410/
      write(6,*)('isymmrci in convert2torc: '),isymmrci,loc(kff20)
      do ir=1,nroot
       vx2(ir)=0d0
       do k=1,ncsfti
        if(abs(veci(k,ir)).gt.vx2(ir))vx2(ir)=abs(veci(k,ir))
       end do
       vx2(ir)=0.99d0*vx2(ir)
      end do
      igoal=66897+1
      irgoal=2
      loop=0
      do i=1,norb
       ismv(i)=ism(i)
       irelv(i)=irel(i)
      end do
      next=norb+1
      do isb=1,nsymb                                                    7d22s21
       if(nvirt(isb).gt.200)stop 'convert2tor iptb'
       do iv=1,nvirt(isb)                                               7d15s21
        ismv(next)=isb                                                  7d15s21
        irelv(next)=iv+irefo(isb)                                       7d15s21
        iptb(iv,isb)=next                                               7d15s21
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
       write(6,*)('for ii = '),ii
       write(6,*)(nff20(jj,1),jj=mdon+1,mdoop)
       iarg=ii-mdon
       nff20(ii,1)=nff20(ii,1)+nff0(ii,1)                               7d23s21
       write(6,*)('nff0: '),nff0(ii,1)
       if(nff20(ii,2).lt.0)then
        nff20(ii,2)=jff20                                                7d16s21
        nff20(ii,3)=ivoff                                                7d16s21
       end if
       jf=nff0(ii,2)                                                    7d19s21
       jv=nff0(ii,3)                                                    7d19s21
       do if=1,nff0(ii,1)                                               7d19s21
        write(6,*)('for internal fcn '),jff20/2,if,loc(kff20(jff20)),
     $       kff20(2)
        kff20(jff20)=iff0(jf)                                           7d19s21
        call dcbit(kff20(jff20),31,'closed')
        jff20=jff20+1                                                   7d19s21
        jf=jf+1                                                         7d19s21
        kff20(jff20)=iff0(jf)                                           7d19s21
        call dcbit(kff20(jff20),31,'open')
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
       do ipass=1,2
        ioffvd=1                                                          8d3s21
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
              lgot=sz.gt.1d-10
              if(lgot)then                                              8d5s21
               write(6,*)('contracting visv vd '),isb,isbv1,if,ioffvd
               call prntm2(vd(ioffvd),nvirt(isbv1)*nroot,nfdat(2,1,isb),
     $             nvirt(isbv1)*nroot)
               write(6,*)('cols '),(iff22(iad1+iii),iii=0,nl(1)-1)
               last=iff22(iad1)
               do iii=1,nl(1)-1
                if(iff22(iad1+iii).le.last.or.
     $              iff22(iad1+iii).gt.nfdat(2,1,isb))then
                 write(6,*)('bad col!!! '),last,iad1+iii
                 call dws_synca
                 call dws_finalize
                 stop
                end if
                last=iff22(iad1+iii)
               end do
               write(6,*)('with block vectors '),if,ii,isb,iad2,1,iarg,
     $              jvcv,nfdat(5,1,isb),nff22(ii,2,isb)
               call prntm2(ff22(iad2),ncsf2(1,iarg),nl(1),ncsf2(1,iarg))
              end if
              ivoff0=ivoff
              do iv=1,nvirt(isbv1)                                      8d3s21
               write(6,*)('for visv fcn '),jff20/2,iv,if,isbv1,
     $              isb,ii,loc(kff20(jff20)),kff20(2)
               kff20(jff20)=ibset(ipack4(1),iptb(iv,isbv1))             8d3s21
               call dcbit(kff20(jff20),31,'closed')
               jff20=jff20+1                                            8d3s21
               kff20(jff20)=ipack4(2)                                   8d3s21
               call dcbit(kff20(jff20),31,'open')
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
                  orig=vnew(igoal,irgoal)
                  vnew(ivoff+i,ir)=vnew(ivoff+i,ir)+ff22(jad2+i)        8d3s21
     $                 *vd(iad)                                         8d3s21
                  if(abs(orig-vnew(igoal,irgoal)).gt.1d-10)write(6,*)
     $                 ('vnew visv: '),orig,ff22(jad2+i),vd(iad),
     $                 vnew(igoal,irgoal)
                 end do                                                 8d3s21
                 jad2=jad2+ncsf2(1,iarg)                                8d3s21
                end do                                                  8d3s21
               end do                                                   8d3s21
               ivoff=ivoff+ncsf2(1,iarg)                                8d3s21
              end do                                                    8d3s21
              nxnv=nvirt(isbv1)*ncsf2(1,iarg)
              sz=0d0
              do ir=1,nroot
               do i=0,nxnv-1
                sz=sz+vnew(ivoff0+i,ir)**2
               end do
              end do
              sz=sqrt(sz/dfloat(nroot*nxnv))
              if(sz.gt.1d-10)then
               write(6,*)('to yield '),vnew(igoal,irgoal),ivoff0-igoal
               call prntm2(vnew(ivoff0,1),nxnv,nroot,nsing)
               if(.not.lgot)write(6,*)('but wait - vd was zero!!!')
              end if
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
             lgot=.false.                                               8d5s21
             do l=1,4
              if(nfdat(2,l,isb).gt.0)then
               sz=0d0                                                   8d5s21
               do iz=0,nvv*nroot*nfdat(2,l,isb)-1
                sz=sz+vd(ioffl(l)+iz)**2
               end do
               sz=sqrt(sz/dfloat(nvv*nroot*nfdat(2,l,isb)))
               if(sz.gt.1d-10)then
                lgot=.true.
                write(6,*)('transforming vnotv vd '),ioffl(l),l,if,
     $              isb,isbv1,isbv2
                call prntm2(vd(ioffl(l)),nvv*nroot,nfdat(2,l,isb),
     $             nvv*nroot)
                nl(l)=iff22(jvcv+1+l)                                     8d3s21
                iad1=jvcv+iff22(jvcv+5+l)                                8d3s21
                write(6,*)('cols '),(iff22(iad1+iii),iii=0,nl(l)-1)
                last=iff22(iad1)
                do iii=1,nl(l)-1
                 if(iff22(iad1+iii).le.last.or.
     $              iff22(iad1+iii).gt.nfdat(2,l,isb))then
                  write(6,*)('bad col!!! '),iad1+iii,iad1
                  call dws_synca
                  call dws_finalize
                  stop
                 end if
                 last=iff22(iad1+iii)
                end do
                iad2=iad1+nl(l)                                           8d3s21
                write(6,*)('by block vectors '),l,if,ii,isb,iad2
                call prntm2(ff22(iad2),ncsf2(l,iarg),nl(l),
     $               ncsf2(l,iarg))
               end if                                                   8d5s21
              end if
             end do
             ivoff0=ivoff
             do iv2=0,nvirt(isbv2)-1                                    8d3s21
              iv2p=iv2+1                                                8d3s21
              itop=(iv2+isw*(nvirt(isbv1)-iv2))-1                       8d3s21
              do iv1=0,itop                                             8d3s21
               iv1p=iv1+1
               irec=iv1+nvirt(isbv1)*iv2                                8d3s21
               itri=((iv2*(iv2-1))/2)+iv1                               8d3s21
               irow=itri+isw*(irec-itri)                                8d3s21
               kff20(jff20)=ipack4(1)                                   8d3s21
               write(6,*)('for vnotv fcn '),jff20/2,iv1,iv2,if,isbv1,
     $              isbv2,isb,ii,loc(kff20(jff20)),kff20(2)
               call dcbit(kff20(jff20),31,'closed')
               jff20=jff20+1                                            8d3s21
               kff20(jff20)=ibset(ipack4(2),iptb(iv1p,isbv1))            8d3s21
               kff20(jff20)=ibset(kff20(jff20),iptb(iv2p,isbv2))         8d3s21
               call dcbit(kff20(jff20),31,'open')
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
                   orig=vnew(igoal,irgoal)
                   vnew(jvoff+i,ir)=vnew(jvoff+i,ir)+ff22(jad2+i)       8d3s21
     $                    *vd(iad)                                        8d3s21
                   if(abs(orig-vnew(igoal,irgoal)).gt.1d-10)write(6,*)
     $                  ('vnew vnotv: '),orig,ff22(jad2+i),vd(iad),
     $                  vnew(igoal,irgoal)
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
             sz=0d0
             do ir=1,nroot
              do i=0,nxnv-1
               sz=sz+vnew(ivoff0+i,ir)**2
              end do
             end do
             sz=sqrt(sz/dfloat(nxnv*nroot))
             if(sz.gt.1d-10)then
              write(6,*)('to yield'),vnew(igoal,irgoal),ivoff0-igoal
              call prntm2(vnew(ivoff0,1),nxnv,nroot,nsing)
              if(.not.lgot)write(6,*)('wait - vd was zero!')
             end if
             jvcv=jvcv+nspace                                           8d3s21
            end do                                                      8d3s21
           end if                                                        8d3s21
          end if                                                        8d3s21
         end do                                                         8d3s21
        end do
       end do
      end do
      write(6,*)('jff20 at end'),jff20,kff20(2)
      write(6,*)('ioffvd at end '),ioffvd
      write(6,*)('ivoff at end is '),ivoff,nsing
      ibcoff=iff20+(jff20/2)
      return                                                            8d3s21
      end                                                               8d3s21
