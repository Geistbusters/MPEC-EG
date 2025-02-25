c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine deconvert12torc(iff22,nff22,nfdat,ncsf,ncsf2,           8d9s21
     $     isymmrci,nvirt,nroot,mdon,mdoop,nsymb,multh,ntot,gunc,       8d9s21
     $     ff22,mdoub,nff0,gc,ndoub,nff1,bc,ibc)                        11d10s22
      implicit real*8 (a-h,o-z)                                         8d3s21
      logical lgot                                                      8d5s21
      integer*8 iff22(*),ipack8                                         8d3s21
      dimension ipack4(2),nff22(mdoop,2,*),nfdat(5,4,*),ncsf(*),        8d9s21
     $     ncsf2(4,*),nvirt(*),multh(8,8),gunc(ntot,*),                 8d9s21
     $     nff20(mdoop,3),nff0(mdoop,3),nff1(mdoop,nsymb),              8d16s21
     $     nl(4),ff22(*),ioffl(4),gc(*)                                 8d16s21
      equivalence (ipack8,ipack4)                                       8d3s21
      include "common.store"                                            8d3s21
      write(6,*)('isymmrci in deconvert12torc: '),isymmrci,ndoub,
     $     loc(ndoub)
      igoal=177+1
      do iz=0,mdoub*nroot                                               8d9s21
       gc(iz)=0d0                                                       8d9s21
      end do                                                            8d9s21
      igoff=1                                                           7d15s21
      ifoff=0                                                           7d15s21
      jff20=1                                                           7d15s21
      icntr=0
      do ii=mdon+1,mdoop                                                7d16s21
       write(6,*)('for ii = '),ii
       iarg=ii-mdon
       igoff=igoff+nff0(ii,1)*ncsf(iarg)
       do isb=1,nsymb                                                   8d13s21
        isbv=multh(isb,isymmrci)                                        8d13s21
        igoff=igoff+nff1(ii,isb)*ncsf(iarg)*nvirt(isbv)                 8d13s21
       end do                                                           8d13s21
       do ipass=1,2
        ioffgd=1
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
              nl(1)=iff22(jvcv+2)
              iad1=jvcv+iff22(jvcv+6)                                   8d3s21
c     iad2 is ncsf2xnl
              iad2=iad1+nl(1)                                           8d3s21
              nxnv=nvirt(isbv1)*ncsf2(1,iarg)
              sz=0d0
              do ir=1,nroot
               do i=0,nxnv-1
                sz=sz+gunc(igoff+i,ir)**2
               end do
              end do
              sz=sqrt(sz/dfloat(nroot*nxnv))
              lgot=sz.gt.1d-14
              if(lgot)then                                              8d5s21
               last=iff22(iad1)
               do iii=1,nl(1)-1
                if(iff22(iad1+iii).le.last.or.
     $              iff22(iad1+iii).gt.nfdat(2,1,isb))then
                 write(6,*)('bad col!!! ')
                 call dws_synca
                 call dws_finalize
                 stop
                end if
                last=iff22(iad1+iii)
               end do
              end if
              do iv=1,nvirt(isbv1)                                      8d3s21
               do iii=0,nl(1)-1                                         8d9s21
                kad2=iad2+ncsf2(1,iarg)*iii                             8d9s21
                do ir=1,nroot                                           8d9s21
                 jgc=ioffgd+iv-1+nvirt(isbv1)*(ir-1+nroot*                 8d9s21
     $                (iff22(iad1+iii)-1))                              8d9s21
                 do i=0,ncsf2(1,iarg)-1                                  8d9s21
                  orig=gc(igoal)
                  gc(jgc)=gc(jgc)+ff22(kad2+i)*gunc(igoff+i,ir)         8d9s21
                  if(abs(orig-gc(igoal)).gt.1d-14)write(6,*)
     $                 ('goalvisv '),
     $                 orig,ff22(kad2+i),gunc(igoff+i,ir),gc(igoal),
     $                 isb,isbv1,iv,ii,igoff+i,ir,ioffgd,igoal-ioffgd
                 end do                                                 8d9s21
                end do                                                  8d9s21
               end do                                                   8d9s21
               igoff=igoff+ncsf2(1,iarg)                                8d3s21
              end do
              jvcv=jvcv+nspace                                          8d3s21
             end do                                                     8d3s21
            end if                                                      8d3s21
            nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                       8d3s21
            ioffgd=ioffgd+nroot*nvirt(isbv1)*nfdat(2,1,isb)             8d3s21
           else                                                         8d3s21
            nvv=nvirt(isbv1)*nvirt(isbv2)                               8d3s21
           end if                                                       8d3s21
           nft=0
           do l=1,4
            ioffl(l)=ioffgd                                             8d9s21
            if(nfdat(2,l,isb).gt.0)then                                 8d3s21
             ioffgd=ioffgd+nroot*nvv*nfdat(2,l,isb)                     8d9s21
             nft=nft+nfdat(2,l,isb)
            end if                                                      8d3s21
           end do                                                       8d3s21
           if(nff22(ii,1,isb).gt.0.and.ipass.eq.1)then                  8d3s21
            jvcv=nfdat(5,1,isb)+nff22(ii,2,isb)                         8d3s21
            do if=1,nff22(ii,1,isb)                                     8d3s21
             ipack8=iff22(jvcv)                                         8d3s21
             nspace=iff22(jvcv+1)                                       8d3s21
             do ivv=0,nvv-1                                             8d9s21
              do l=1,4
               if(nfdat(2,l,isb).gt.0)then
                nl(l)=iff22(jvcv+1+l)                                     8d3s21
                iad1=jvcv+iff22(jvcv+5+l)                                8d3s21
                iad2=iad1+nl(l)                                           8d3s21
                do iii=0,nl(l)-1
                 do ir=1,nroot                                          8d9s21
                  jgc=ioffl(l)+ivv+nvv*(ir-1+nroot*(iff22(iad1+iii)-1)) 8d9s21
                  do i=0,ncsf2(l,iarg)-1
                   orig=gc(igoal)
                   gc(jgc)=gc(jgc)+ff22(iad2+i)*gunc(i+igoff,ir)        8d9s21
                  if(abs(orig-gc(igoal)).gt.1d-14)write(6,*)
     $                 ('goalvnotv '),
     $                 orig,ff22(iad2+i),gunc(igoff+i,ir),gc(igoal),
     $                 isb,isbv1,isbv2,l,ivv,ii,igoff+i,ir
                  end do                                                8d9s21
                 end do                                                 8d9s21
                 iad2=iad2+ncsf2(l,iarg)                                8d9s21
                end do                                                  8d9s21
               end if                                                   8d5s21
               igoff=igoff+ncsf2(l,iarg)                                8d9s21
              end do                                                    8d9s21
             end do
             jvcv=jvcv+nspace                                           8d3s21
            end do                                                      8d3s21
           end if                                                        8d3s21
          end if                                                        8d3s21
         end do                                                         8d3s21
        end do
       end do
      end do
      write(6,*)('doubles vector in non-orthogonal basis: '),gc(igoal)
      ioffgd=1                                                          8d9s21
      do isb=1,nsymb                                                    8d9s21
       isbv12=multh(isb,isymmrci)                                       8d9s21
       do isbv1=1,nsymb                                                 8d9s21
        isbv2=multh(isbv1,isbv12)                                       8d9s21
        if(isbv1.le.isbv2)then                                          8d9s21
         if(isbv1.eq.isbv2)then
          sz=0d0
          nrow=nvirt(isbv1)*nroot                                       8d9s21
          do is=0,nrow*nfdat(2,1,isb)-1
           sz=sz+gc(ioffgd+is)**2
          end do
          sz=sqrt(sz/dfloat(nrow*nfdat(2,1,isb)))
          if(sz.gt.1d-14)then
           write(6,*)('visv for syms '),isb,isbv1,ioffgd
           call prntm2(gc(ioffgd),nrow,nfdat(2,1,isb),nrow)              8d9s21
          end if
          ioffgd=ioffgd+nrow*nfdat(2,1,isb)                             8d9s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         8d9s21
         else                                                           8d9s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 8d9s21
         end if
         nrow=nvv*nroot                                                 8d9s21
         do l=1,4                                                       8d9s21
          if(nfdat(2,l,isb).gt.0)then                                   8d9s21
           sz=0d0
           do is=0,nrow*nfdat(2,l,isb)-1
            sz=sz+gc(ioffgd+is)**2
           end do
           sz=sqrt(sz/dfloat(nrow*nfdat(2,l,isb)))
           if(sz.gt.1d-14)then
            write(6,*)('vnotv for spins '),l,('syms '),isb,isbv1,isbv2, 8d9s21
     $           ioffgd
            call prntm2(gc(ioffgd),nrow,nfdat(2,l,isb),nrow)             8d9s21
           end if
           ioffgd=ioffgd+nrow*nfdat(2,l,isb)                            8d9s21
          end if                                                        8d9s21
         end do                                                         8d9s21
        end if                                                          8d9s21
       end do                                                           8d9s21
      end do                                                            8d9s21
      itmpg=ibcoff                                                      8d12s21
      ibcoff=itmpg+ndoub*nroot                                          8d12s21
      write(6,*)('enough for tmpg '),itmpg,ndoub,nroot,ibcoff
      call enough('deconvert12torc.  1',bc,ibc)
      call tofrob(bc(itmpg),gc,nroot,nfdat,nvirt,nsymb,multh,isymmrci,2,
     $     ndoub,mdoub,iff22,bc,ibc)                                    11d10s22
      write(6,*)('gd in orthogonal basis ')                             8d12s21
      jtmpg=itmpg                                                       8d12s21
      do isb=1,nsymb                                                    8d12s21
       isbv12=multh(isb,isymmrci)
       do isbv1=1,nsymb                                                 8d12s21
        isbv2=multh(isbv1,isbv12)                                       8d12s21
        if(isbv2.ge.isbv1)then                                          8d12s21
         if(isbv1.eq.isbv2)then                                         8d12s21
          nrow=nroot*nvirt(isbv1)                                       8d12s21
          sz=0d0                                                        8d12s21
          do iz=0,nrow*nfdat(3,1,isb)-1                                 8d12s21
           sz=sz+bc(jtmpg+iz)**2                                        8d12s21
          end do                                                        8d12s21
          sz=sqrt(sz/dfloat(nrow*nfdat(3,1,isb)))                       8d12s21
          if(sz.gt.1d-14)then                                           8d12s21
           write(6,*)('visv for syms '),isb,isbv1,jtmpg-itmpg
           call prntm2(bc(jtmpg),nrow,nfdat(3,1,isb),nrow)              8d12s21
          end if                                                        8d12s21
          jtmpg=jtmpg+nrow*nfdat(3,1,isb)                               8d12s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         8d12s21
         else                                                           8d12s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 8d12s21
         end if                                                         8d12s21
         nrow=nvv*nroot                                                 8d12s21
         do l=1,4                                                       8d12s21
          if(nfdat(3,l,isb).gt.0)then                                   8d12s21
           sz=0d0                                                       8d12s21
           do iz=0,nrow*nfdat(3,l,isb)-1                                8d12s21
            sz=sz+bc(jtmpg+iz)**2                                       8d12s21
           end do                                                       8d12s21
           sz=sqrt(sz/dfloat(nrow*nfdat(3,l,isb)))                      8d12s21
           if(sz.gt.1d-14)then                                          8d12s21
            write(6,*)('vnotv for spins '),l,('symmetries '),isb,isbv1,
     $          isbv2,jtmpg-itmpg
            call prntm2(bc(jtmpg),nrow,nfdat(3,l,isb),nrow)             8d12s21
           end if                                                       8d12s21
           jtmpg=jtmpg+nrow*nfdat(3,l,isb)                              8d12s211
          end if                                                        8d12s21
         end do                                                         8d12s21
        end if                                                          8d12s21
       end do                                                           8d12s21
      end do                                                            8d12s21
      ibcoff=itmpg                                                      8d12s21
      return                                                            8d3s21
      end                                                               8d3s21
