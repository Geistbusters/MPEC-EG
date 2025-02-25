c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine sotest2c(iff22,nff22,nfdat,vd,ncsf,ncsf2,isymmrci,     8d3s21
     $     nvirt,test,nroot,mdon,mdoop,nsymb,multh,ndoubx,vx2,ff22,bc,  11d10s22
     $     ibc)                                                         11d10s22
      implicit real*8 (a-h,o-z)                                         8d3s21
      integer*8 iff22(*),ipack8                                         8d5s21
      dimension nff22(mdoop,2,*),nfdat(5,4,*),vd(*),ncsf(*),            8d5s21
     $     ncsf2(4,*),nvirt(*),test(nroot,nroot),multh(8,8),vx2(nroot), 8d5s21
     $     ff22(*),ivcx(4),ipack4(2)                                    8d5s21
      equivalence (ipack8,ipack4)                                       8d5s21
      include "common.store"                                            8d5s21
      data icall/0/
      save icall
      icall=icall+1
       thrx=0.99d0
      do isb=1,nsymb                                                    8d5s21
       ibc0=ibcoff                                                      8d5s21
       do l=1,4                                                         8d5s21
        ivcx(l)=ibcoff                                                  8d5s21
        ibcoff=ivcx(l)+nfdat(3,l,isb)                                   8d5s21
       end do                                                           8d5s21
       call enough('sotest2c.  1',bc,ibc)
       do iz=ibc0,ibcoff-1                                              8d5s21
        bc(iz)=0d0                                                      8d5s21
       end do                                                           8d5s21
       do ii=mdon+1,mdoop                                               8d5s21
        if(nff22(ii,1,isb).gt.0)then                                    8d5s21
         iarg=ii-mdon                                                   8d5s21
         jvcv=nfdat(5,1,isb)+nff22(ii,2,isb)                            8d5s21
         do if=1,nff22(ii,1,isb)                                        8d5s21
          nspace=iff22(jvcv+1)                                          8d5s21
          do l=1,4                                                      8d5s21
           nl=iff22(jvcv+1+l)                                           8d5s21
           iad1=jvcv+iff22(jvcv+5+l)                                    8d5s21
           iad2=iad1+nl                                                 8d5s21
           do iii=0,nl-1                                                8d5s21
            iiii=ivcx(l)+iff22(iad1+iii)-1                              8d5s21
            do j=0,ncsf2(l,iarg)-1                                      8d5s21
             tt=abs(ff22(iad2+j))
             if(tt.gt.bc(iiii))bc(iiii)=tt                              8d5s21
            end do                                                      8d5s21
            iad2=iad2+ncsf2(l,iarg)                                     8d5s21
           end do                                                       8d5s21
          end do                                                        8d5s21
          jvcv=jvcv+nspace                                              8d5s21
         end do
        end if                                                          8d5s21
       end do                                                           8d5s21
       if(icall.eq.-1)then
       do l=1,4
        if(nfdat(3,l,isb).gt.0)then
         write(6,*)('take transformation matrix for l = '),l,('from')
         call prntm2(ff22(nfdat(4,l,isb)),nfdat(2,l,isb),
     $        nfdat(3,l,isb),nfdat(2,l,isb))
         do i=0,nfdat(3,l,isb)-1
          do j=0,nfdat(2,l,isb)-1
           iad=nfdat(4,l,isb)+j+nfdat(2,l,isb)*i
           if(j.ne.i)then
            ff22(iad)=0d0
           else                                                         8d6s21
            ff22(iad)=1d0                                               8d6s21
           end if
          end do
         end do
         write(6,*)('to ')
         call prntm2(ff22(nfdat(4,l,isb)),nfdat(2,l,isb),
     $        nfdat(3,l,isb),nfdat(2,l,isb))
        end if
       end do
       end if
       do l=1,4
        if(nfdat(3,l,isb).gt.0)then
         do i=0,nfdat(3,l,isb)-1
          bc(ivcx(l)+i)=bc(ivcx(l)+i)*thrx                              8d5s21
         end do
        end if
       end do
       do ii=mdon+1,mdoop                                               8d5s21
        if(nff22(ii,1,isb).gt.0)then                                    8d5s21
         iarg=ii-mdon                                                   8d5s21
         jvcv=nfdat(5,1,isb)+nff22(ii,2,isb)                            8d5s21
         do if=1,nff22(ii,1,isb)                                        8d5s21
          ipack8=iff22(jvcv)                                            8d5s21
          nspace=iff22(jvcv+1)                                          8d5s21
          if(icall.eq.-1)then
           write(6,*)('for isb,ii,if: '),isb,ii,if
           write(6,*)('we are at jvcv = '),jvcv,loc(iff22(jvcv))
           write(6,*)('nspace = '),nspace
          do l=1,4                                                      8d5s21
           nl=iff22(jvcv+1+l)                                           8d5s21
           write(6,*)('nl for l = '),l,(' is '),nl
           if(min(ncsf2(l,iarg),nl).gt.0)then
            iad1=jvcv+iff22(jvcv+5+l)                                    8d5s21
            write(6,*)('cols at '),iad1,(' are '),
     $           (iff22(iad1+iii),iii=0,nl-1)
            iad2=iad1+nl                                                 8d5s21
            write(6,*)('while vectors at '),iad2,(' are')
            iad20=iad2                                                  11d23s21
            call prntm2(ff22(iad2),ncsf2(l,iarg),nl,ncsf2(l,iarg))
            do iii=0,nl-1                                                8d5s21
             iiii=ivcx(l)+iff22(iad1+iii)-1                              8d5s21
             do j=0,ncsf2(l,iarg)-1                                      8d5s21
              tt=abs(ff22(iad2+j))
              if(tt.gt.bc(iiii))then                                     8d5s21
               ff22(iad2+j)=1d0
               write(6,*)('keep '),tt,j,iff22(iad1+iii),l,if,ii,
     $             isb
               call dcbit(ipack4,10,'closed')
               call dcbit(ipack4(2),10,'open')
              else
               ff22(iad2+j)=0d0
               write(6,*)('zeroing '),iad2+j
              end if
             end do                                                      8d5s21
             iad2=iad2+ncsf2(l,iarg)                                     8d5s21
            end do                                                       8d5s21
            call prntm2(ff22(iad20),ncsf2(l,iarg),nl,ncsf2(l,iarg))     11d23s21
           end if
          end do                                                        8d5s21
          end if
          jvcv=jvcv+nspace                                              8d5s21
         end do
        end if                                                          8d5s21
       end do                                                           8d5s21
       ibcoff=ibc0                                                      8d5s21
      end do                                                            8d5s21
      do ir=1,nroot
       vx2(ir)=0d0
      end do
      ndoubx=0                                                          8d3s21
      ioffvd=1                                                          8d3s21
      do isb=1,nsymb                                                    8d3s21
       isbv12=multh(isb,isymmrci)                                       8d3s21
       nvisv=0                                                          8d3s21
       nvnotv=0                                                         8d3s21
       do isbv1=1,nsymb                                                 8d3s21
        isbv2=multh(isbv1,isbv12)                                       8d3s21
        if(isbv2.ge.isbv1)then                                          8d3s21
         if(isbv1.eq.isbv2)then                                         8d3s21
          nvisv=nvisv+nvirt(isbv1)                                      8d3s21
          do k=0,nfdat(3,1,isb)-1                                       8d3s21
           do ir=1,nroot                                                8d3s21
            iad=ioffvd+nvirt(isbv1)*(ir-1+nroot*k)                      8d3s21
            do iv=0,nvirt(isbv1)-1                                      8d3s21
             if(vd(iad+iv).ne.vd(iad+iv))then
              write(6,*)('vd is NaN!!! '),iv,ir,k,isbv1,isb,iad+iv
              call dws_synca
              call dws_finalize
              stop
             end if
             if(abs(vd(iad+iv)).gt.vx2(ir))vx2(ir)=abs(vd(iad+iv))      8d3s21
            end do                                                      8d3s21
            do jr=1,nroot                                               8d3s21
             jad=ioffvd+nvirt(isbv1)*(jr-1+nroot*k)                     8d3s21
             do iv=0,nvirt(isbv1)-1                                     8d3s21
              test(jr,ir)=test(jr,ir)+vd(iad+iv)*vd(jad+iv)             8d3s21
             end do                                                     8d3s21
            end do                                                      8d3s21
           end do                                                       8d3s21
          end do                                                        8d3s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         8d3s21
          ioffvd=ioffvd+nvirt(isbv1)*nroot*nfdat(3,1,isb)               8d3s21
         else                                                           8d3s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 8d3s21
         end if                                                         8d3s21
         nvnotv=nvnotv+nvv                                              8d3s21
         do l=1,4                                                       8d3s21
          if(nfdat(3,l,isb).gt.0)then                                   8d3s21
           do k=0,nfdat(3,l,isb)-1                                      8d3s21
            do ir=1,nroot                                                8d3s21
             iad=ioffvd+nvv*(ir-1+nroot*k)                              8d3s21
             do iv=0,nvv-1                                              8d3s21
              if(vd(iad+iv).ne.vd(iad+iv))then
               write(6,*)('vd is NaNb!!! '),iv,ir,k,l,isbv1,isbv2,isb,
     $              iad+iv,loc(vd(iad+iv))
               call dws_synca
               call dws_finalize
               stop
              end if
              if(abs(vd(iad+iv)).gt.vx2(ir))vx2(ir)=abs(vd(iad+iv))      8d3s21
             end do                                                      8d3s21
             do jr=1,nroot                                               8d3s21
              jad=ioffvd+nvv*(jr-1+nroot*k)                             8d3s21
              do ivv=0,nvv-1                                             8d3s21
               test(jr,ir)=test(jr,ir)+vd(iad+ivv)*vd(jad+ivv)          8d3s21
              end do                                                     8d3s21
             end do                                                      8d3s21
            end do                                                       8d3s21
           end do                                                       8d3s21
           ioffvd=ioffvd+nvv*nroot*nfdat(3,l,isb)                       8d3s21
          end if                                                        8d3s21
         end do                                                         8d3s21
        end if                                                          8d3s21
       end do                                                           8d3s21
       do ii=mdon+1,mdoop                                               8d3s21
        if(nff22(ii,1,isb).gt.0)then                                    8d3s21
         iarg=ii-mdon                                                   8d3s21
         jvcv=nfdat(5,1,isb)+nfdat(ii,2,isb)                            8d3s21
         ndoubx=ndoubx+ncsf2(1,iarg)*nvisv*nff22(ii,1,isb)              8d3s21
     $        +ncsf(iarg)*nff22(ii,1,isb)*nvnotv                        8d3s21
        end if                                                          8d3s21
       end do                                                           8d3s21
      end do                                                            8d3s21
      do i=1,nroot
      end do
      ioffvd=1                                                          8d3s21
      igotone=0                                                         12d30s21
      do isb=1,nsymb                                                    8d3s21
       isbv12=multh(isb,isymmrci)                                       8d3s21
       do isbv1=1,nsymb                                                 8d3s21
        isbv2=multh(isbv1,isbv12)                                       8d3s21
        if(isbv2.ge.isbv1)then                                          8d3s21
         if(isbv1.eq.isbv2)then                                         8d3s21
          do k=0,nfdat(3,1,isb)-1                                       8d3s21
           do ir=1,nroot                                                8d3s21
            iad=ioffvd+nvirt(isbv1)*(ir-1+nroot*k)                      8d3s21
            do iv=0,nvirt(isbv1)-1                                      8d3s21
            end do                                                      8d3s21
           end do                                                       8d3s21
          end do                                                        8d3s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         8d3s21
          ioffvd=ioffvd+nvirt(isbv1)*nroot*nfdat(3,1,isb)               8d3s21
         else                                                           8d3s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 8d3s21
         end if                                                         8d3s21
         nvnotv=nvnotv+nvv                                              8d3s21
         do l=1,4                                                       8d3s21
          if(nfdat(3,l,isb).gt.0)then                                   8d3s21
           do k=0,nfdat(3,l,isb)-1                                      8d3s21
            do ir=1,nroot                                                8d3s21
             iad=ioffvd+nvv*(ir-1+nroot*k)                              8d3s21
             do iv=0,nvv-1                                              8d3s21
             end do                                                      8d3s21
            end do                                                       8d3s21
           end do                                                       8d3s21
           ioffvd=ioffvd+nvv*nroot*nfdat(3,l,isb)                       8d3s21
          end if                                                        8d3s21
         end do                                                         8d3s21
        end if                                                          8d3s21
       end do                                                           8d3s21
      end do                                                            8d3s21
      return                                                            8d3s21
      end                                                               8d3s21
