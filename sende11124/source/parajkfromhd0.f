c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine parajkfromhd0(nocc,multh,jmats,kmats,ioooo,ionex,iter, 7d26s16
     $     idwsdeb,nvirtc,iprop,itrans,i4od,ionexd,nbasdwsc,isblkder,   8d3s16
     $     isblkxder,nsblkder,nsblkxder,isblkkder,nsblkkder,jmatd,kmatd,8d4s16
     $     i3x,i4od2b,ionexd2,isblkxder1,nsblkxder1,isblkder1,nsblkder1,5d4s22
     $     idxdom,igoal,bc,ibc)                                         11d10s22
      implicit real*8 (a-h,o-z)
c
c     like jkfromh, except for operators of symmetry iprop.
c     this version gives the contribution to the integrals from the
c     undifferentiated ao integrals.
c     by suitable choice of inputs, this will do the first ders of      7d27s16
c     4o and onex and the pure second ders of 4o and onex.              7d27s16
c     if nsblkkder gt 0, we will do first ders of j and k.              8d4s16
c
      external second
      logical ltest                                                     3d29s13
      integer*8 ifiddle8
      include "common.store"
      include "common.hf"
      dimension nocc(1),jmats(1),kmats(1),jmatd(1),kmatd(1)             8d22s16
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension multh(8,8),itrans(1),nbasdwsc(8),ioooo(1),i3x(idbk),    7d26s16
     $     ionex(idbk),nvirtc(8),i4od(idbk),ionexd(idbk),               8d19s16
     $     isblkder(4,idbk),isblkxder(4,idbk),isblkkder(4,idbk),        8d31s16
     $     i4od2b(idbk),ionexd2(idbk),isblkxder1(4,idbk),                9d12s16
     $     isblkder1(4,idbk)                                            9d12s16
      data icall/0/                                                     5d9s22
      save icall                                                        5d9s22
      icall=icall+1                                                     5d9s22
      call second(time1)                                                11d27s12
      kgoal=2176320
      if(idwsdeb.gt.10)then
       write(6,*)('in parajkfromhd0 for iprop = '),iprop,nsblkxder,
     $      ibcoff
       write(6,*)('for icall = '),icall
       write(6,*)('idwsdeb = '),idwsdeb
       write(6,*)('what we have for trans: ')
       do isb=1,nsymb
        isk=multh(isb,iprop)
        if(min(nbasdwsc(isb),nbasdwsc(isk)).gt.0)then                   11d28s22
         write(6,*)('for symmetry block '),isb,isk,itrans(isb),
     $       nbasdwsc(isb),nbasdwsc(isk)
         call prntm2(bc(itrans(isb)),nbasdwsc(isb),nbasdwsc(isk),
     $      nbasdwsc(isb))
        end if                                                          11d28s22
       end do
      end if                                                            5d19s22
      do is=1,nsblkxder
       i4od(is)=ibcoff
       nrow=nocc(isblkxder(1,is))*nocc(isblkxder(2,is))
       ncol=nocc(isblkxder(3,is))*nocc(isblkxder(4,is))
       nn=nrow*ncol
       if(nn.gt.0)then
        ibcoff=i4od(is)+nrow*ncol
        call enough('parajkfromhd0.  1',bc,ibc)
        do i=0,nrow*ncol-1
         bc(i4od(is)+i)=0d0
        end do
c
c     dxor contribution
c
        is4=multh(isblkxder(4,is),iprop)
        if(nocc(is4).eq.0)go to 1010                                    4d18s16
        igot=0
        do is2=1,nsdlk
         if(isblk(1,is2).eq.isblkxder(1,is).and.
     $      isblk(2,is2).eq.isblkxder(2,is).and.
     $      isblk(3,is2).eq.isblkxder(3,is).and.
     $      isblk(4,is2).eq.is4)then
          igot=1
          nocc1=nocc(isblk(1,is2))
          nocc2=nocc(isblk(2,is2))
          if(isblk(1,is2).eq.isblk(2,is2))then
           nr4o=(nocc1*(nocc1+1))/2
          else
           nr4o=nocc1*nocc2
          end if
          call ilimts(nocc(isblk(3,is2)),nocc(isblk(4,is2)),mynprocg,   3d15s16
     $         mynowprog,ilo,iho,io1s,io1e,io2s,io2e)                   3d15s16
          noccm=nocc(isblkxder(4,is))
          nocc3=nocc(isblk(3,is2))
          nocc4=nocc(isblk(4,is2))
          n3=nbasdwsc(isblkxder(4,is))                                  3d24s16
          n4=nbasdwsc(isblk(4,is2))                                     3d24s16
          do im=0,noccm-1
           i10=io1s
           i1n=nocc3
           ii0=ioooo(is2)
           do i2=io2s,io2e
            if(i2.eq.io2e)i1n=io1e
            i2m=i2-1
            iad2=itrans(isblk(4,is2))+i2m+n4*im
            do i1=i10,i1n
             i1m=i1-1
             iad1=i4od(is)+nrow*(i1m+nocc3*im)
             if(isblk(1,is2).eq.isblk(2,is2))then
              do i3=0,nocc2-1
               do i4=0,nocc1-1
                in=min(i3,i4)
                ix=max(i3,i4)
                ii=ii0+((ix*(ix+1))/2)+in
                term=bc(ii)*bc(iad2)
                orig=bc(iad1)
                bc(iad1)=bc(iad1)+bc(ii)*bc(iad2)
                iad1=iad1+1
               end do
              end do
              ii0=ii0+nr4o
             else
              do i3=0,nocc2-1
               do i4=0,nocc1-1
                bc(iad1)=bc(iad1)+bc(ii0)*bc(iad2)
                iad1=iad1+1
                ii0=ii0+1
               end do
              end do
             end if
            end do
            i10=1
           end do
          end do
          go to 1010
         end if
         if(isblk(1,is2).eq.isblkxder(1,is).and.
     $      isblk(2,is2).eq.isblkxder(2,is).and.
     $      isblk(4,is2).eq.isblkxder(3,is).and.
     $      isblk(3,is2).eq.is4)then
          igot=1
          nocc1=nocc(isblk(1,is2))
          nocc2=nocc(isblk(2,is2))
          if(isblk(1,is2).eq.isblk(2,is2))then
           nr4o=(nocc1*(nocc1+1))/2
          else
           nr4o=nocc1*nocc2
          end if
          call ilimts(nocc(isblk(3,is2)),nocc(isblk(4,is2)),mynprocg,   3d15s16
     $         mynowprog,ilo,iho,io1s,io1e,io2s,io2e)                   3d15s16
          noccm=nocc(isblkxder(4,is))
          nocc3=nocc(isblk(3,is2))
          nocc4=nocc(isblk(4,is2))
          n3=nbasdwsc(isblkxder(4,is))                                  3d24s16
          n4=nbasdwsc(isblk(3,is2))                                     3d24s16
          do im=0,noccm-1
           i10=io1s
           i1n=nocc3
           ii0=ioooo(is2)
           do i2=io2s,io2e
            if(i2.eq.io2e)i1n=io1e
            i2m=i2-1
            iad10=i4od(is)+nrow*(i2m+nocc4*im)
            do i1=i10,i1n
             i1m=i1-1
             iad2=itrans(isblk(3,is2))+i1m+n4*im
             iad1=iad10
             if(isblk(1,is2).eq.isblk(2,is2))then
              do i3=0,nocc2-1
               do i4=0,nocc1-1
                in=min(i3,i4)
                ix=max(i3,i4)
                ii=ii0+((ix*(ix+1))/2)+in
                bc(iad1)=bc(iad1)+bc(ii)*bc(iad2)
                iad1=iad1+1
               end do
              end do
              ii0=ii0+nr4o
             else
              do i3=0,nocc2-1
               do i4=0,nocc1-1
                bc(iad1)=bc(iad1)+bc(ii0)*bc(iad2)
                iad1=iad1+1
                ii0=ii0+1
               end do
              end do
             end if
            end do
            i10=1
           end do
          end do
          go to 1010
         end if
         if(isblk(2,is2).eq.isblkxder(1,is).and.
     $      isblk(1,is2).eq.isblkxder(2,is).and.
     $      isblk(3,is2).eq.isblkxder(3,is).and.
     $      isblk(4,is2).eq.is4)then
          igot=1
          nocc1=nocc(isblk(1,is2))
          nocc2=nocc(isblk(2,is2))
          if(isblk(1,is2).eq.isblk(2,is2))then
           nr4o=(nocc1*(nocc1+1))/2
          else
           nr4o=nocc1*nocc2
          end if
          call ilimts(nocc(isblk(3,is2)),nocc(isblk(4,is2)),mynprocg,   3d15s16
     $         mynowprog,ilo,iho,io1s,io1e,io2s,io2e)                   3d15s16
          noccm=nocc(isblkxder(4,is))
          nocc3=nocc(isblk(3,is2))
          nocc4=nocc(isblk(4,is2))
          n3=nbasdwsc(isblkxder(4,is))                                  3d24s16
          n4=nbasdwsc(isblk(4,is2))                                     3d24s16
          do im=0,noccm-1
           i10=io1s
           i1n=nocc3
           ii0=ioooo(is2)
           do i2=io2s,io2e
            if(i2.eq.io2e)i1n=io1e
            i2m=i2-1
            iad2=itrans(isblk(4,is2))+i2m+n4*im
            do i1=i10,i1n
             i1m=i1-1
             iad1=i4od(is)+nrow*(i1m+nocc3*im)
             do i3=0,nocc2-1
              do i4=0,nocc1-1
               jad1=iad1+i3+nocc2*i4                                    4d18s16
               bc(jad1)=bc(jad1)+bc(ii0)*bc(iad2)                       4d18s16
               ii0=ii0+1
              end do
             end do
            end do
            i10=1
           end do
          end do
          go to 1010
         end if
        end do
        if(igot.eq.0.and.nocc(is4).gt.0)then                            7d11s22
         write(6,*)('could not find oooo match for (oooo)'' ')
         write(6,*)('want ')
         write(6,12)(isblkxder(j,is),j=1,3),is4
         write(6,*)('got ')
         do is2=1,nsdlk
          write(6,12)(isblk(j,is2),j=1,4),is2
         end do
         call dws_sync
         call dws_finalize
         stop
        end if
 1010   continue
        igot=0
        do is2=1,nsdlk1
         if(isblk1(1,is2).eq.isblkxder(1,is).and.
     $      isblk1(2,is2).eq.isblkxder(2,is).and.
     $      isblk1(3,is2).eq.isblkxder(3,is).and.
     $      isblk1(4,is2).eq.is4)then
          igot=1
          nocc1=nocc(isblk1(1,is2))
          nocc2=nocc(isblk1(2,is2))
          if(isblk1(1,is2).eq.isblk1(2,is2))then
           nr1x=(nocc1*(nocc1+1))/2
          else
           nr1x=nocc1*nocc2
          end if
          nocc3=nocc(isblk1(3,is2))
          n4v=nvirtc(isblk1(4,is2))
          call ilimts(nocc3,n4v,mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,    3d15s16
     $         ix2s,ix2e)                                               3d15s16
          noccm=nocc(isblkxder(4,is))
          n3=nbasdwsc(isblkxder(4,is))                                  3d24s16
          n4=nbasdwsc(isblk1(4,is2))                                    3d24s16
          do im=0,noccm-1
           i10=ix1s
           i1n=nocc3
           ii0=ionex(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2m=i2-1
            i2p=i2m+nocc(isblk1(4,is2))
             iad2=itrans(isblk1(4,is2))+i2p+n4*im
            do i1=i10,i1n
             i1m=i1-1
             iad1=i4od(is)+nrow*(i1m+nocc3*im)
             if(isblk1(1,is2).eq.isblk1(2,is2))then
              do i3=0,nocc2-1
               do i4=0,nocc1-1
                in=min(i3,i4)
                ix=max(i3,i4)
                ii=ii0+((ix*(ix+1))/2)+in
                bc(iad1)=bc(iad1)+bc(ii)*bc(iad2)
                iad1=iad1+1
               end do
              end do
              ii0=ii0+nr1x
             else
              do i3=0,nocc2-1
               do i4=0,nocc1-1
                bc(iad1)=bc(iad1)+bc(ii0)*bc(iad2)
                ii0=ii0+1
                iad1=iad1+1
               end do
              end do
             end if
            end do
            i10=1
           end do
          end do
          go to 1011
         end if
         if(isblk1(2,is2).eq.isblkxder(1,is).and.                       4d18s16
     $      isblk1(1,is2).eq.isblkxder(2,is).and.                       4d18s16
     $      isblk1(3,is2).eq.isblkxder(3,is).and.
     $      isblk1(4,is2).eq.is4)then
          igot=1
          nocc1=nocc(isblk1(1,is2))
          nocc2=nocc(isblk1(2,is2))
          if(isblk1(1,is2).eq.isblk1(2,is2))then
           nr1x=(nocc1*(nocc1+1))/2
          else
           nr1x=nocc1*nocc2
          end if
          nocc3=nocc(isblk1(3,is2))
          n4v=nvirtc(isblk1(4,is2))
          call ilimts(nocc3,n4v,mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,    3d15s16
     $         ix2s,ix2e)                                               3d15s16
          noccm=nocc(isblkxder(4,is))
          n3=nbasdwsc(isblkxder(4,is))                                  3d24s16
          n4=nbasdwsc(isblk1(4,is2))                                    3d24s16
          do im=0,noccm-1
           i10=ix1s
           i1n=nocc3
           ii0=ionex(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2m=i2-1
            i2p=i2m+nocc(isblk1(4,is2))
             iad2=itrans(isblk1(4,is2))+i2p+n4*im
            do i1=i10,i1n
             i1m=i1-1
             iad1=i4od(is)+nrow*(i1m+nocc3*im)
             do i3=0,nocc2-1
              do i4=0,nocc1-1
               jad1=iad1+i3+nocc2*i4                                    4d18s16
               bc(jad1)=bc(jad1)+bc(ii0)*bc(iad2)
               ii0=ii0+1
              end do
             end do
            end do
            i10=1
           end do
          end do
          go to 1011
         end if
        end do
        if(igot.eq.0.and.nvirt(is4).gt.0)then                           7d11s22
         write(6,*)('could not find onex match for (oooo)'' ')
         write(6,*)('want ')
         write(6,12)(isblkxder(j,is),j=1,3),is4
         write(6,*)('got ')
         do is2=1,nsdlk1
          write(6,12)(isblk1(j,is2),j=1,4),is2
         end do
         call dws_sync
         call dws_finalize
         stop
        end if
 1011   continue
        if(idwsdeb.gt.10)then
         write(6,12)(isblkxder(j,is),j=1,4)
          write(6,*)('final oooo'' '),i4od(is)
         if(mynprocg.eq.1)then
          call prntm2(bc(i4od(is)),nrow,ncol,nrow)
         else
          do i=0,nrow*ncol-1
           bc(ibcoff+i)=bc(i4od(is)+i)
          end do
          call dws_gsumf(bc(ibcoff),nrow*ncol)
          call prntm2(bc(ibcoff),nrow,ncol,nrow)
         end if
        end if
       end if
      end do
      if(nsblkkder.eq.-10)return                                        5d16s22
      do is=1,nsblkxder
       ionexd(is)=ibcoff
       nrow=nocc(isblkxder(1,is))*nocc(isblkxder(2,is))
       ncol=nocc(isblkxder(3,is))*nvirtc(isblkxder(4,is))
       nn=nrow*ncol
       if(nn.gt.0)then
        ibcoff=ionexd(is)+nn
        call enough('parajkfromhd0.  3',bc,ibc)
        do i=0,nn-1
         bc(ionexd(is)+i)=0d0
        end do
c
c     xor contribution to o''oox: ooox and K
c
        is1=multh(isblkxder(1,is),iprop)
        if(nocc(is1).eq.0)go to 1101                                    4d18s16
        igot=0
        do is2=1,nsdlk1
         if(isblk1(1,is2).eq.is1.and.isblk1(2,is2).eq.isblkxder(2,is)
     $        .and.isblk1(3,is2).eq.isblkxder(3,is).and.
     $        isblk1(4,is2).eq.isblkxder(4,is))then
          call ilimts(nocc(isblk1(3,is2)),nvirtc(isblk1(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          if(isblk1(1,is2).eq.isblk1(2,is2))then
           nrowx=(nocc(isblk1(1,is2))*(nocc(isblk1(1,is2))+1))/2
          else
           nrowx=nocc(isblk1(1,is2))*nocc(isblk1(2,is2))
          end if
          nm=nbasdwsc(isblk1(1,is2))                                    3d24s16
          do im=0,nocc(isblkxder(1,is))-1
           iad2=itrans(isblk1(1,is2))+nm*im
           i10=ix1s
           i1n=nocc(isblk1(3,is2))
           ii0=ionex(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2m=i2-1
            do i1=i10,i1n
             i1m=i1-1
             if(isblk1(2,is2).ne.isblk1(1,is2))then
              do i3=0,nocc(isblk1(2,is2))-1
               iad1=ionexd(is)+im+nocc(isblkxder(1,is))*(i3
     $               +nocc(isblk1(2,is2))*(i1m+nocc(isblk1(3,is2))*i2m))
               do i4=0,nocc(isblk1(1,is2))-1
                bc(iad1)=bc(iad1)+bc(iad2+i4)*bc(ii0)
                ii0=ii0+1
               end do
              end do
             else
              do i3=0,nocc(isblk1(2,is2))-1
               iad1=ionexd(is)+im+nocc(isblkxder(1,is))*(i3
     $               +nocc(isblk1(2,is2))*(i1m+nocc(isblk1(3,is2))*i2m))
               do i4=0,nocc(isblk1(1,is2))-1
                ix=max(i3,i4)
                in=min(i3,i4)
                ii=((ix*(ix+1))/2)+in+ii0
                bc(iad1)=bc(iad1)+bc(iad2+i4)*bc(ii)
               end do
              end do
              ii0=ii0+nrowx
             end if
            end do
            i10=1
           end do
          end do
          igot=1
         else if(isblk1(2,is2).eq.is1.and.
     $         isblk1(1,is2).eq.isblkxder(2,is)
     $        .and.isblk1(3,is2).eq.isblkxder(3,is).and.
     $        isblk1(4,is2).eq.isblkxder(4,is))then
          call ilimts(nocc(isblk1(3,is2)),nvirtc(isblk1(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          nm=nbasdwsc(isblk1(2,is2))                                    3d24s16
          do im=0,nocc(isblkxder(1,is))-1
           iad2=itrans(isblk1(2,is2))+nm*im
           i10=ix1s
           i1n=nocc(isblk1(3,is2))
           ii0=ionex(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2m=i2-1
            do i1=i10,i1n
             i1m=i1-1
             do i3=0,nocc(isblk1(2,is2))-1
              iad2=itrans(isblk1(2,is2))+nm*im+i3
              do i4=0,nocc(isblk1(1,is2))-1
              iad1=ionexd(is)+im+nocc(isblkxder(1,is))*(i4
     $              +nocc(isblk1(1,is2))*(i1m+nocc(isblk1(3,is2))*i2m))
               bc(iad1)=bc(iad1)+bc(iad2)*bc(ii0)
               ii0=ii0+1
              end do
             end do
            end do
            i10=1
           end do
          end do
          igot=1
         end if
         if(igot.ne.0)go to 1101
        end do
        write(6,*)('could not find ooox for o''oox: '),is1,
     $       (isblkxder(j,is),j=2,4)
        call dws_sync
        call dws_finalize
        stop
 1101   continue
        igot=0
        do is2=1,nsdlkk
         if(isblkk(1,is2).eq.isblkxder(2,is)
     $        .and.isblkk(2,is2).eq.isblkxder(3,is)
     $        .and.isblkk(3,is2).eq.isblkxder(4,is).and.
     $        isblkk(4,is2).eq.is1)then
          call ilimts(nvirtc(isblkk(3,is2)),nvirtc(isblkk(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          nm=nbasdwsc(isblkk(4,is2))                                    3d24s16
          do im=0,nocc(isblkxder(1,is))-1
           i10=ix1s
           i1n=nvirtc(isblkk(3,is2))
           ii0=kmats(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2p=i2-1+nocc(isblkk(4,is2))
            iad2=itrans(isblkk(4,is2))+i2p+nm*im
            do i1=i10,i1n
             i1m=i1-1
             do i3=0,nocc(isblkk(2,is2))-1
              do i4=0,nocc(isblkk(1,is2))-1
               iad1=ionexd(is)+im+nocc(isblkxder(1,is))*(i4
     $               +nocc(isblkk(1,is2))*(i3+nocc(isblkk(2,is2))*i1m))
               bc(iad1)=bc(iad1)+bc(iad2)*bc(ii0)
               ii0=ii0+1
              end do
             end do
            end do
            i10=1
           end do
          end do
          igot=1
         else if(isblkk(1,is2).eq.isblkxder(3,is)
     $        .and.isblkk(2,is2).eq.isblkxder(2,is)
     $        .and.isblkk(3,is2).eq.is1.and.
     $        isblkk(4,is2).eq.isblkxder(4,is))then
          call ilimts(nvirtc(isblkk(3,is2)),nvirtc(isblkk(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          nm=nbasdwsc(isblkk(3,is2))                                    3d24s16
          do im=0,nocc(isblkxder(1,is))-1
           i10=ix1s
           i1n=nvirtc(isblkk(3,is2))
           ii0=kmats(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2m=i2-1
            do i1=i10,i1n
             i1p=i1-1+nocc(isblkk(3,is2))
             iad2=itrans(isblkk(3,is2))+i1p+nm*im
             do i3=0,nocc(isblkk(2,is2))-1
              do i4=0,nocc(isblkk(1,is2))-1
               iad1=ionexd(is)+im+nocc(isblkxder(1,is))*(i3
     $               +nocc(isblkk(2,is2))*(i4+nocc(isblkk(1,is2))*i2m))
               bc(iad1)=bc(iad1)+bc(iad2)*bc(ii0)
               ii0=ii0+1
              end do
             end do
            end do
            i10=1
           end do
          end do
          igot=1
         end if
         if(igot.ne.0)go to 1102
        end do
        write(6,*)('could not find kmat for o''oox: '),
     $       isblkxder(2,is),isblkxder(3,is),isblkxder(4,is),is1
        do is2=1,nsdlkk
         write(6,12)(isblkk(j,is2),j=1,4),is2
        end do
        call dws_sync
        call dws_finalize
        stop
 1102   continue
c
c     xor contribution to oo''ox: ooox and K
c
        is1=multh(isblkxder(2,is),iprop)
        if(nocc(is1).eq.0)go to 1103                                    4d18s16
        igot=0
        do is2=1,nsdlk1
         if(isblk1(1,is2).eq.isblkxder(1,is).and.isblk1(2,is2).eq.is1   3d21s16
     $        .and.isblk1(3,is2).eq.isblkxder(3,is).and.
     $        isblk1(4,is2).eq.isblkxder(4,is))then
          call ilimts(nocc(isblk1(3,is2)),nvirtc(isblk1(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          if(isblk1(1,is2).eq.isblk1(2,is2))then
           nrowx=(nocc(isblk1(1,is2))*(nocc(isblk1(1,is2))+1))/2
          else
           nrowx=nocc(isblk1(1,is2))*nocc(isblk1(2,is2))
          end if
          nm=nbasdwsc(isblk1(2,is2))                                    3d24s16
          do im=0,nocc(isblkxder(2,is))-1
           i10=ix1s
           i1n=nocc(isblk1(3,is2))
           ii0=ionex(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2m=i2-1
            do i1=i10,i1n
             i1m=i1-1
             if(isblk1(2,is2).ne.isblk1(1,is2))then
              do i3=0,nocc(isblk1(2,is2))-1
               iad2=itrans(isblk1(2,is2))+nm*im+i3
               iad1=ionexd(is)+nocc(isblkxder(1,is))*(im
     $             +nocc(isblkxder(2,is))*(i1m+nocc(isblk1(3,is2))*i2m))
               do i4=0,nocc(isblk1(1,is2))-1
                bc(iad1+i4)=bc(iad1+i4)+bc(iad2)*bc(ii0)
                ii0=ii0+1
               end do
              end do
             else
              do i3=0,nocc(isblk1(2,is2))-1
               iad1=ionexd(is)+nocc(isblkxder(1,is))*(im
     $             +nocc(isblkxder(2,is))*(i1m+nocc(isblk1(3,is2))*i2m))
               iad2=itrans(isblk1(2,is2))+nm*im+i3
               do i4=0,nocc(isblk1(1,is2))-1
                ix=max(i3,i4)
                in=min(i3,i4)
                ii=((ix*(ix+1))/2)+in+ii0
                bc(iad1+i4)=bc(iad1+i4)+bc(iad2)*bc(ii)
               end do
              end do
              ii0=ii0+nrowx
             end if
            end do
            i10=1
           end do
          end do
          igot=1
         else if(isblk1(2,is2).eq.isblkxder(1,is).and.                  3d21s16
     $         isblk1(1,is2).eq.is1                                     3d21s16
     $        .and.isblk1(3,is2).eq.isblkxder(3,is).and.
     $        isblk1(4,is2).eq.isblkxder(4,is))then
          call ilimts(nocc(isblk1(3,is2)),nvirtc(isblk1(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          nm=nbasdwsc(isblk1(1,is2))                                    3d24s16
          do im=0,nocc(isblkxder(2,is))-1
           iad2=itrans(isblk1(1,is2))+nm*im
           i10=ix1s
           i1n=nocc(isblk1(3,is2))
           ii0=ionex(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2m=i2-1
            do i1=i10,i1n
             i1m=i1-1
             do i3=0,nocc(isblk1(2,is2))-1
              iad1=ionexd(is)+i3+nocc(isblkxder(1,is))*(im
     $             +nocc(isblkxder(2,is))*(i1m+nocc(isblk1(3,is2))*i2m))
              do i4=0,nocc(isblk1(1,is2))-1
               bc(iad1)=bc(iad1)+bc(iad2+i4)*bc(ii0)
               ii0=ii0+1
              end do
             end do
            end do
            i10=1
           end do
          end do
          igot=1
         end if
         if(igot.ne.0)go to 1103
        end do
        write(6,*)('could not find ooox for oo''ox: '),isblkxder(1,is),
     $       is1,(isblkxder(j,is),j=3,4)
        do is2=1,nsdlk1
         write(6,12)(isblk1(j,is2),j=1,4),is2
        end do
        call dws_sync
        call dws_finalize
        stop
 1103   continue
c k_{ad}^{cb}=(ab|dc)
        igot=0
        do is2=1,nsdlkk
         if(isblkk(1,is2).eq.isblkxder(1,is)
     $        .and.isblkk(2,is2).eq.isblkxder(3,is)
     $        .and.isblkk(3,is2).eq.isblkxder(4,is).and.
     $        isblkk(4,is2).eq.is1)then
          call ilimts(nvirtc(isblkk(3,is2)),nvirtc(isblkk(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          nm=nbasdwsc(isblkk(4,is2))                                    3d24s16
          do im=0,nocc(isblkxder(2,is))-1
           i10=ix1s
           i1n=nvirtc(isblkk(3,is2))
           ii0=kmats(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2p=i2-1+nocc(isblkk(4,is2))
            iad2=itrans(isblkk(4,is2))+i2p+nm*im
            do i1=i10,i1n
             i1m=i1-1
             do i3=0,nocc(isblkk(2,is2))-1
              iad1=ionexd(is)+nocc(isblkxder(1,is))*(im
     $              +nocc(isblkxder(2,is))*(i3+nocc(isblkk(2,is2))*i1m))
              do i4=0,nocc(isblkk(1,is2))-1
               bc(iad1+i4)=bc(iad1+i4)+bc(iad2)*bc(ii0)
               ii0=ii0+1
              end do
             end do
            end do
            i10=1
           end do
          end do
          igot=1
         else if(isblkk(1,is2).eq.isblkxder(3,is)
     $        .and.isblkk(2,is2).eq.isblkxder(1,is)
     $        .and.isblkk(3,is2).eq.is1.and.
     $        isblkk(4,is2).eq.isblkxder(4,is))then
          call ilimts(nvirtc(isblkk(3,is2)),nvirtc(isblkk(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          nm=nbasdwsc(isblkk(3,is2))                                    3d24s16
          do im=0,nocc(isblkxder(2,is))-1
           i10=ix1s
           i1n=nvirtc(isblkk(3,is2))
           ii0=kmats(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2m=i2-1
            do i1=i10,i1n
             i1p=i1-1+nocc(isblkk(3,is2))
             iad2=itrans(isblkk(3,is2))+i1p+nm*im
             do i3=0,nocc(isblkk(2,is2))-1
              do i4=0,nocc(isblkk(1,is2))-1
               iad1=ionexd(is)+i3+nocc(isblkk(2,is2))*(im
     $              +nocc(isblkxder(2,is))*(i4+nocc(isblkk(1,is2))*i2m))
               bc(iad1)=bc(iad1)+bc(iad2)*bc(ii0)
               ii0=ii0+1
              end do
             end do
            end do
            i10=1
           end do
          end do
          igot=1
         end if
         if(igot.ne.0)go to 1104
        end do
        write(6,*)('could not find kmat for oo''ox: '),
     $       isblkxder(1,is),is1,isblkxder(3,is),isblkxder(4,is)
        call dws_sync
        call dws_finalize
        stop
 1104   continue
c
c     xor contribution to ooo''ox: ooox and J
c
        is1=multh(isblkxder(3,is),iprop)
        igot=0
        if(nocc(is1).eq.0)go to 1105                                    4d19s22
        do is2=1,nsdlk1
         if(isblk1(1,is2).eq.isblkxder(1,is).and.                       3d21s16
     $        isblk1(2,is2).eq.isblkxder(2,is)                          3d21s16
     $        .and.isblk1(3,is2).eq.is1.and.                            3d21s16
     $        isblk1(4,is2).eq.isblkxder(4,is))then
          call ilimts(nocc(isblk1(3,is2)),nvirtc(isblk1(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          if(isblk1(1,is2).eq.isblk1(2,is2))then
           nrowx=(nocc(isblk1(1,is2))*(nocc(isblk1(1,is2))+1))/2
          else
           nrowx=nocc(isblk1(1,is2))*nocc(isblk1(2,is2))
          end if
          nm=nbasdwsc(isblk1(3,is2))                                    3d24s16
          do im=0,nocc(isblkxder(3,is))-1
           i10=ix1s
           i1n=nocc(isblk1(3,is2))
           ii0=ionex(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2m=i2-1
            do i1=i10,i1n
             i1m=i1-1
             iad2=itrans(isblk1(3,is2))+nm*im+i1m
             if(isblk1(2,is2).ne.isblk1(1,is2))then
              do i3=0,nocc(isblk1(2,is2))-1
               iad1=ionexd(is)+nocc(isblk1(1,is2))*(i3
     $              +nocc(isblk1(2,is2))*(im+nocc(isblkxder(3,is))*i2m))
               do i4=0,nocc(isblk1(1,is2))-1
                bc(iad1+i4)=bc(iad1+i4)+bc(iad2)*bc(ii0)
                ii0=ii0+1
               end do
              end do
             else
              do i3=0,nocc(isblk1(2,is2))-1
               iad1=ionexd(is)+nocc(isblk1(1,is2))*(i3
     $              +nocc(isblk1(2,is2))*(im+nocc(isblkxder(3,is))*i2m))
               do i4=0,nocc(isblk1(1,is2))-1
                ix=max(i3,i4)
                in=min(i3,i4)
                ii=((ix*(ix+1))/2)+in+ii0
                bc(iad1+i4)=bc(iad1+i4)+bc(iad2)*bc(ii)
               end do
              end do
              ii0=ii0+nrowx
             end if
            end do
            i10=1
           end do
          end do
          igot=1
          go to 1105
         end if
         if(isblk1(2,is2).eq.isblkxder(1,is).and.                       4d18s16
     $        isblk1(1,is2).eq.isblkxder(2,is)                          4d18s16
     $        .and.isblk1(3,is2).eq.is1.and.                            3d21s16
     $        isblk1(4,is2).eq.isblkxder(4,is))then
          call ilimts(nocc(isblk1(3,is2)),nvirtc(isblk1(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          if(isblk1(1,is2).eq.isblk1(2,is2))then
           nrowx=(nocc(isblk1(1,is2))*(nocc(isblk1(1,is2))+1))/2
          else
           nrowx=nocc(isblk1(1,is2))*nocc(isblk1(2,is2))
          end if
          nm=nbasdwsc(isblk1(3,is2))                                    3d24s16
          do im=0,nocc(isblkxder(3,is))-1
           i10=ix1s
           i1n=nocc(isblk1(3,is2))
           ii0=ionex(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2m=i2-1
            do i1=i10,i1n
             i1m=i1-1
             iad2=itrans(isblk1(3,is2))+nm*im+i1m
             do i3=0,nocc(isblk1(2,is2))-1
              iad1=ionexd(is)+i3+nocc(isblk1(2,is2))*(                  4d18s16
     $              +nocc(isblk1(1,is2))*(im+nocc(isblkxder(3,is))*i2m))4d18s16
              do i4=0,nocc(isblk1(1,is2))-1
               jad1=iad1+i4*nocc(isblk1(2,is2))                         4d18s16
               bc(jad1)=bc(jad1)+bc(iad2)*bc(ii0)                       4d18s16
               ii0=ii0+1
              end do
             end do
            end do
            i10=1
           end do
          end do
          igot=1
          go to 1105
         end if
        end do
        write(6,*)('could not find ooox for ooo''x: '),isblkxder(1,is),
     $       isblkxder(2,is),is1,isblkxder(4,is)
        call dws_sync
        call dws_finalize
        stop
 1105   continue
        igot=0
        do is2=1,nsdlk
         if(isblk(1,is2).eq.isblkxder(1,is)
     $        .and.isblk(2,is2).eq.isblkxder(2,is)
     $        .and.isblk(3,is2).eq.is1.and.
     $        isblk(4,is2).eq.isblkxder(4,is))then
          call ilimts(nvirtc(isblk(3,is2)),nvirtc(isblk(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          if(isblk(1,is2).eq.isblk(2,is2))then                          3d22s16
           nrowj=(nocc(isblk(1,is2))*(nocc(isblk(1,is2))+1))/2
          else
           nrowj=nocc(isblk(1,is2))*nocc(isblk(2,is2))
          end if
          nm=nbasdwsc(isblk(3,is2))                                     3d24s16
          do im=0,nocc(isblkxder(3,is))-1
           i10=ix1s
           i1n=nvirtc(isblk(3,is2))
           ii0=jmats(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2m=i2-1
            do i1=i10,i1n
             i1p=i1-1+nocc(isblk(3,is2))
             iad2=itrans(isblk(3,is2))+i1p+nm*im
             if(isblk(1,is2).ne.isblk(2,is2))then                       3d21s16
              do i3=0,nocc(isblk(2,is2))-1
               iad1=ionexd(is)+nocc(isblk(1,is2))*(i3
     $               +nocc(isblk(2,is2))*(im+nocc(isblkxder(3,is))*i2m))
               do i4=0,nocc(isblk(1,is2))-1
                bc(iad1+i4)=bc(iad1+i4)+bc(iad2)*bc(ii0)
                ii0=ii0+1
               end do
              end do
             else
              do i3=0,nocc(isblk(2,is2))-1
               iad1=ionexd(is)+nocc(isblk(1,is2))*(i3
     $               +nocc(isblk(2,is2))*(im+nocc(isblkxder(3,is))*i2m))
               do i4=0,nocc(isblk(1,is2))-1
                ix=max(i3,i4)
                in=min(i3,i4)
                ii=ii0+((ix*(ix+1))/2)+in
                bc(iad1+i4)=bc(iad1+i4)+bc(iad2)*bc(ii)
               end do
              end do
              ii0=ii0+nrowj
             end if
            end do
            i10=1
           end do
          end do
          igot=1
         else if(isblk(1,is2).eq.isblkxder(1,is)
     $        .and.isblk(2,is2).eq.isblkxder(2,is)
     $        .and.isblk(4,is2).eq.is1.and.
     $        isblk(3,is2).eq.isblkxder(4,is))then
          call ilimts(nvirtc(isblk(3,is2)),nvirtc(isblk(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          nm=nbasdwsc(isblk(4,is2))                                     3d24s16
          do im=0,nocc(isblkxder(3,is))-1
           i10=ix1s
           i1n=nvirtc(isblk(3,is2))
           ii0=jmats(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2p=i2-1+nocc(isblk(4,is2))
            iad2=itrans(isblk(4,is2))+i2p+nm*im
            do i1=i10,i1n
             i1m=i1-1
             do i3=0,nocc(isblk(2,is2))-1
              iad1=ionexd(is)+nocc(isblk(1,is2))*(i3
     $              +nocc(isblk(2,is2))*(im+nocc(isblkxder(3,is))*i1m))
              do i4=0,nocc(isblk(1,is2))-1
               bc(iad1+i4)=bc(iad1+i4)+bc(iad2)*bc(ii0)
               ii0=ii0+1
              end do
             end do
            end do
            i10=1
           end do
          end do
          igot=1
         else if(isblk(2,is2).eq.isblkxder(1,is)                        4d18s16
     $        .and.isblk(1,is2).eq.isblkxder(2,is)                      4d18s16
     $        .and.isblk(3,is2).eq.is1.and.
     $        isblk(4,is2).eq.isblkxder(4,is))then
          call ilimts(nvirtc(isblk(3,is2)),nvirtc(isblk(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          if(isblk(1,is2).eq.isblk(2,is2))then                          3d22s16
           nrowj=(nocc(isblk(1,is2))*(nocc(isblk(1,is2))+1))/2
          else
           nrowj=nocc(isblk(1,is2))*nocc(isblk(2,is2))
          end if
          nm=nbasdwsc(isblk(3,is2))                                     3d24s16
          do im=0,nocc(isblkxder(3,is))-1
           i10=ix1s
           i1n=nvirtc(isblk(3,is2))
           ii0=jmats(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2m=i2-1
            do i1=i10,i1n
             i1p=i1-1+nocc(isblk(3,is2))
             iad2=itrans(isblk(3,is2))+i1p+nm*im
             do i3=0,nocc(isblk(2,is2))-1
              iad1=ionexd(is)+i3+nocc(isblk(2,is2))*(                   4d18s16
     $              +nocc(isblk(1,is2))*(im+nocc(isblkxder(3,is))*i2m)) 4d18s16
              do i4=0,nocc(isblk(1,is2))-1
               jad1=iad1+i4*nocc(isblk(2,is2))                          4d18s16
               bc(jad1)=bc(jad1)+bc(iad2)*bc(ii0)                       4d18s16
               ii0=ii0+1
              end do
             end do
            end do
            i10=1
           end do
          end do
          igot=1
         else if(isblk(1,is2).eq.isblkxder(2,is)                        4d18s16
     $        .and.isblk(2,is2).eq.isblkxder(1,is)                      4d18s16
     $        .and.isblk(4,is2).eq.is1.and.
     $        isblk(3,is2).eq.isblkxder(4,is))then
          call ilimts(nvirtc(isblk(3,is2)),nvirtc(isblk(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          nm=nbasdwsc(isblk(4,is2))                                     3d24s16
          do im=0,nocc(isblkxder(3,is))-1
           i10=ix1s
           i1n=nvirtc(isblk(3,is2))
           ii0=jmats(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2p=i2-1+nocc(isblk(4,is2))
            iad2=itrans(isblk(4,is2))+i2p+nm*im
            do i1=i10,i1n
             i1m=i1-1
             do i3=0,nocc(isblk(2,is2))-1
              iad1=ionexd(is)+i3+nocc(isblk(2,is2))*(                   4d18s16
     $              +nocc(isblk(1,is2))*(im+nocc(isblkxder(3,is))*i1m)) 4d18s16
              do i4=0,nocc(isblk(1,is2))-1
               jad1=iad1+i4*nocc(isblk(2,is2))                          4d18s16
               bc(jad1)=bc(jad1)+bc(iad2)*bc(ii0)                       4d18s16
               ii0=ii0+1
              end do
             end do
            end do
            i10=1
           end do
          end do
          igot=1
         end if
         if(igot.ne.0)go to 1106
        end do
        write(6,*)('could not find jmat for ooo''x: '),
     $       isblkxder(1,is),isblkxder(2,is),is1,isblkxder(4,is)
        do is2=1,nsdlk
         write(6,12)(isblk(j,is2),j=1,4),is2
        end do
        call dws_sync
        call dws_finalize
        stop
 1106   continue
c
c     xor contribution to ooox'': oooo and ooox
c
        is1=multh(isblkxder(4,is),iprop)
        if(nocc(is1).eq.0)go to 1107                                    4d18s16
        igot=0
        do is2=1,nsdlk
         if(isblk(1,is2).eq.isblkxder(1,is)
     $        .and.isblk(2,is2).eq.isblkxder(2,is)
     $        .and.isblk(3,is2).eq.isblkxder(3,is).and.
     $        isblk(4,is2).eq.is1)then
          call ilimts(nocc(isblk(3,is2)),nocc(isblk(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          if(isblk(1,is2).eq.isblk(2,is2))then
           nrowj=(nocc(isblk(1,is2))*(nocc(isblk(1,is2))+1))/2
          else
           nrowj=nocc(isblk(1,is2))*nocc(isblk(2,is2))
          end if
          nm=nbasdwsc(isblk(4,is2))                                     3d24s16
          do im=0,nvirtc(isblkxder(4,is))-1
           imp=im+nocc(isblkxder(4,is))                                 3d21s16
           i10=ix1s
           i1n=nocc(isblk(3,is2))
           ii0=ioooo(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2m=i2-1
            iad2=itrans(isblk(4,is2))+i2m+nm*imp
            do i1=i10,i1n
             i1m=i1-1
             if(isblk(1,is2).ne.isblk(2,is2))then
              do i3=0,nocc(isblk(2,is2))-1
               iad1=ionexd(is)+nocc(isblk(1,is2))*(i3
     $               +nocc(isblk(2,is2))*(i1m+nocc(isblk(3,is2))*im))
               do i4=0,nocc(isblk(1,is2))-1
                bc(iad1+i4)=bc(iad1+i4)+bc(iad2)*bc(ii0)
                ii0=ii0+1
               end do
              end do
             else
              do i3=0,nocc(isblk(2,is2))-1
               iad1=ionexd(is)+nocc(isblk(1,is2))*(i3
     $               +nocc(isblk(2,is2))*(i1m+nocc(isblk(3,is2))*im))
               do i4=0,nocc(isblk(1,is2))-1
                ix=max(i3,i4)
                in=min(i3,i4)
                ii=ii0+((ix*(ix+1))/2)+in
                bc(iad1+i4)=bc(iad1+i4)+bc(iad2)*bc(ii)
               end do
              end do
              ii0=ii0+nrowj
             end if
            end do
            i10=1
           end do
          end do
          igot=1
         else if(isblk(1,is2).eq.isblkxder(1,is)
     $        .and.isblk(2,is2).eq.isblkxder(2,is)
     $        .and.isblk(4,is2).eq.isblkxder(3,is).and.
     $        isblk(3,is2).eq.is1)then
          call ilimts(nocc(isblk(3,is2)),nocc(isblk(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          if(isblk(1,is2).eq.isblk(2,is2))then
           nrowj=(nocc(isblk(1,is2))*(nocc(isblk(1,is2))+1))/2
          else
           nrowj=nocc(isblk(1,is2))*nocc(isblk(2,is2))
          end if
          nm=nbasdwsc(isblk(3,is2))                                     3d24s16
          do im=0,nvirtc(isblkxder(4,is))-1                             3d21s16
           imp=im+nocc(isblkxder(4,is))                                 3d21s16
           i10=ix1s
           i1n=nocc(isblk(3,is2))
           ii0=ioooo(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2m=i2-1
            do i1=i10,i1n
             i1m=i1-1
             iad2=itrans(isblk(3,is2))+i1m+nm*imp                       3d21s16
             if(isblk(1,is2).ne.isblk(2,is2))then
              do i3=0,nocc(isblk(2,is2))-1
               iad1=ionexd(is)+nocc(isblk(1,is2))*(i3
     $               +nocc(isblk(2,is2))*(i2m+nocc(isblk(4,is2))*im))
               do i4=0,nocc(isblk(1,is2))-1
                bc(iad1+i4)=bc(iad1+i4)+bc(iad2)*bc(ii0)
                ii0=ii0+1
               end do
              end do
             else
              do i3=0,nocc(isblk(2,is2))-1
               iad1=ionexd(is)+nocc(isblk(1,is2))*(i3
     $               +nocc(isblk(2,is2))*(i2m+nocc(isblk(4,is2))*im))
               do i4=0,nocc(isblk(1,is2))-1
                ix=max(i3,i4)
                in=min(i3,i4)
                ii=ii0+((ix*(ix+1))/2)+in
                bc(iad1+i4)=bc(iad1+i4)+bc(iad2)*bc(ii)
               end do
              end do
              ii0=ii0+nrowj
             end if
            end do
            i10=1
           end do
          end do
          igot=1
         else if(isblk(2,is2).eq.isblkxder(1,is)                        4d18s16
     $        .and.isblk(1,is2).eq.isblkxder(2,is)                      4d18s16
     $        .and.isblk(3,is2).eq.isblkxder(3,is).and.
     $        isblk(4,is2).eq.is1)then
          call ilimts(nocc(isblk(3,is2)),nocc(isblk(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          if(isblk(1,is2).eq.isblk(2,is2))then
           nrowj=(nocc(isblk(1,is2))*(nocc(isblk(1,is2))+1))/2
          else
           nrowj=nocc(isblk(1,is2))*nocc(isblk(2,is2))
          end if
          nm=nbasdwsc(isblk(4,is2))                                     3d24s16
          do im=0,nvirtc(isblkxder(4,is))-1
           imp=im+nocc(isblkxder(4,is))                                 3d21s16
           i10=ix1s
           i1n=nocc(isblk(3,is2))
           ii0=ioooo(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2m=i2-1
            iad2=itrans(isblk(4,is2))+i2m+nm*imp
            do i1=i10,i1n
             i1m=i1-1
             do i3=0,nocc(isblk(2,is2))-1
              iad1=ionexd(is)+i3+nocc(isblk(2,is2))*(                   4d18s16
     $               +nocc(isblk(1,is2))*(i1m+nocc(isblk(3,is2))*im))   4d18s16
              do i4=0,nocc(isblk(1,is2))-1
               ix=max(i3,i4)
               in=min(i3,i4)
               ii=ii0+((ix*(ix+1))/2)+in
               jad1=iad1+i4*nocc(isblk(2,is2))                          4d18s16
               bc(jad1)=bc(jad1)+bc(iad2)*bc(ii)                        4d18s16
              end do
             end do
             ii0=ii0+nrowj
            end do
            i10=1
           end do
          end do
          igot=1
         else if(isblk(2,is2).eq.isblkxder(1,is)                        4d18s16
     $        .and.isblk(1,is2).eq.isblkxder(2,is)                      4d18s16
     $        .and.isblk(4,is2).eq.isblkxder(3,is).and.
     $        isblk(3,is2).eq.is1)then
          call ilimts(nocc(isblk(3,is2)),nocc(isblk(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          if(isblk(1,is2).eq.isblk(2,is2))then
           nrowj=(nocc(isblk(1,is2))*(nocc(isblk(1,is2))+1))/2
          else
           nrowj=nocc(isblk(1,is2))*nocc(isblk(2,is2))
          end if
          nm=nbasdwsc(isblk(3,is2))                                     3d24s16
          do im=0,nvirtc(isblkxder(4,is))-1                             3d21s16
           imp=im+nocc(isblkxder(4,is))                                 3d21s16
           i10=ix1s
           i1n=nocc(isblk(3,is2))
           ii0=ioooo(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2m=i2-1
            do i1=i10,i1n
             i1m=i1-1
             iad2=itrans(isblk(3,is2))+i1m+nm*imp                       3d21s16
             do i3=0,nocc(isblk(2,is2))-1
              iad1=ionexd(is)+i3+nocc(isblk(2,is2))*(                   4d18s16
     $              +nocc(isblk(1,is2))*(i2m+nocc(isblk(4,is2))*im))    4d18s16
              do i4=0,nocc(isblk(1,is2))-1                              4d18s16
               jad1=iad1+i4*nocc(isblk(2,is2))                          4d18s16
               bc(jad1)=bc(jad1)+bc(iad2)*bc(ii0)                       4d18s16
               ii0=ii0+1
              end do
             end do
            end do
            i10=1
           end do
          end do
          igot=1
         end if
         if(igot.ne.0)go to 1107
        end do
        write(6,*)('could not find oooo for ooox'': '),
     $       (isblkxder(j,is),j=1,3),is1
        do is2=1,nsdlk
         write(6,12)(isblk(j,is2),j=1,4),is2
        end do
        call dws_sync
        call dws_finalize
        stop
 1107   continue
        igot=0
        do is2=1,nsdlk1
         if(isblk1(1,is2).eq.isblkxder(1,is).and.                       3d21s16
     $        isblk1(2,is2).eq.isblkxder(2,is)                          3d21s16
     $        .and.isblk1(3,is2).eq.isblkxder(3,is).and.                3d21s16
     $        isblk1(4,is2).eq.is1)then                                 3d21s16
          call ilimts(nocc(isblk1(3,is2)),nvirtc(isblk1(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          if(isblk1(1,is2).eq.isblk1(2,is2))then
           nrowx=(nocc(isblk1(1,is2))*(nocc(isblk1(1,is2))+1))/2
          else
           nrowx=nocc(isblk1(1,is2))*nocc(isblk1(2,is2))
          end if
          nm=nbasdwsc(isblk1(4,is2))                                    3d24s16
          do im=0,nvirtc(isblkxder(4,is))-1
           imp=im+nocc(isblkxder(4,is))
           i10=ix1s
           i1n=nocc(isblk1(3,is2))
           ii0=ionex(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2p=i2-1+nocc(isblk1(4,is2))
            iad2=itrans(isblk1(4,is2))+nm*imp+i2p
            do i1=i10,i1n
             i1m=i1-1
              if(isblk1(2,is2).ne.isblk1(1,is2))then
              do i3=0,nocc(isblk1(2,is2))-1
               iad1=ionexd(is)+nocc(isblkxder(1,is))*(i3
     $               +nocc(isblk1(2,is2))*(i1m+nocc(isblk1(3,is2))*im))
               do i4=0,nocc(isblk1(1,is2))-1
                bc(iad1+i4)=bc(iad1+i4)+bc(iad2)*bc(ii0)
                ii0=ii0+1
               end do
              end do
             else
              do i3=0,nocc(isblk1(2,is2))-1
               iad1=ionexd(is)+nocc(isblkxder(1,is))*(i3
     $               +nocc(isblk1(2,is2))*(i1m+nocc(isblk1(3,is2))*im))
               do i4=0,nocc(isblk1(1,is2))-1
                ix=max(i3,i4)
                in=min(i3,i4)
                ii=((ix*(ix+1))/2)+in+ii0
                bc(iad1+i4)=bc(iad1+i4)+bc(iad2)*bc(ii)
               end do
              end do
              ii0=ii0+nrowx
             end if
            end do
            i10=1
           end do
          end do
          igot=1
          go to 1108
         end if
         if(isblk1(2,is2).eq.isblkxder(1,is).and.                       4d18s16
     $        isblk1(1,is2).eq.isblkxder(2,is)                          4d18s16
     $        .and.isblk1(3,is2).eq.isblkxder(3,is).and.                3d21s16
     $        isblk1(4,is2).eq.is1)then                                 3d21s16
          call ilimts(nocc(isblk1(3,is2)),nvirtc(isblk1(4,is2)),
     $         mynprocg,mynowprog,ilx,ihx,ix1s,ix1e,ix2s,ix2e)
          if(isblk1(1,is2).eq.isblk1(2,is2))then
           nrowx=(nocc(isblk1(1,is2))*(nocc(isblk1(1,is2))+1))/2
          else
           nrowx=nocc(isblk1(1,is2))*nocc(isblk1(2,is2))
          end if
          nm=nbasdwsc(isblk1(4,is2))                                    3d24s16
          do im=0,nvirtc(isblkxder(4,is))-1
           imp=im+nocc(isblkxder(4,is))
           i10=ix1s
           i1n=nocc(isblk1(3,is2))
           ii0=ionex(is2)
           do i2=ix2s,ix2e
            if(i2.eq.ix2e)i1n=ix1e
            i2p=i2-1+nocc(isblk1(4,is2))
            iad2=itrans(isblk1(4,is2))+nm*imp+i2p
            do i1=i10,i1n
             i1m=i1-1
             do i3=0,nocc(isblk1(2,is2))-1
              iad1=ionexd(is)+i3+nocc(isblkxder(1,is))*(                4d18s16
     $               +nocc(isblk1(1,is2))*(i1m+nocc(isblk1(3,is2))*im)) 4d18s16
              do i4=0,nocc(isblk1(1,is2))-1
               jad1=iad1+i4*nocc(isblkxder(1,is))                       4d18s16
               bc(jad1)=bc(jad1)+bc(iad2)*bc(ii0)                       4d18s16
               ii0=ii0+1
              end do
             end do
            end do
            i10=1
           end do
          end do
          igot=1
          go to 1108
         end if
        end do
        if(nvirt(is1).gt.0)then                                         7d11s22
         write(6,*)('could not find ooox for ooox'': '),
     $       (isblkxder(j,is),j=1,3),is1
         call dws_sync
         call dws_finalize
         stop
        end if                                                          7d11s22
 1108   continue
       end if
      end do
      if(idwsdeb.gt.10)then
       write(6,*)('for (ooox)'' tau contribution not gsummed ')
       do is=1,nsblkxder
        nrow=nocc(isblkxder(1,is))*nocc(isblkxder(2,is))
        ncol=nocc(isblkxder(3,is))*nvirtc(isblkxder(4,is))
        if(min(nrow,ncol).gt.0)then
         write(6,12)(isblkxder(j,is),j=1,4)
         write(6,*)('ionexd(is) '),ionexd(is),ionexd(is)+nrow*ncol
         if(mynprocg.eq.1)then
          call prntm2(bc(ionexd(is)),nrow,ncol,nrow)
          if(nsymb.eq.1)then
           call printa(bc(ionexd(is)),nocc,0,nocc,0,nocc,0,nvirtc,nocc,
     $          bc(ibcoff))
          end if
         else
          do i=0,nrow*ncol-1
           bc(ibcoff+i)=bc(ionexd(is)+i)
          end do
          call dws_gsumf(bc(ibcoff),nrow*ncol)
          call prntm2(bc(ibcoff),nrow,ncol,nrow)
         end if
        end if
        itmp=ibcoff
        ibcoff=itmp+nrow*ncol
        call enough('parajkfromhd0.  4',bc,ibc)
        do i=0,nrow*ncol-1
         bc(itmp+i)=bc(ionexd(is)+i)
        end do
        write(6,*)('gsumf '),nrow,ncol,nrow*ncol
        call dws_gsumf(bc(itmp),nrow*ncol)
        ibcoff=itmp
       end do
       write(6,*)('for (oooo)'' tau contribution not gsummed ')
       do is=1,nsblkxder
        nrow=nocc(isblkxder(1,is))*nocc(isblkxder(2,is))
        ncol=nocc(isblkxder(3,is))*nocc(isblkxder(4,is))
        if(min(nrow,ncol).gt.0)then                                     5d4s22
         write(6,*)('i4od(is) '),i4od(is),i4od(is)+nrow*ncol
         write(6,12)(isblkxder(j,is),j=1,4)                              5d16s22
         if(mynprocg.eq.1)then
          call prntm2(bc(i4od(is)),nrow,ncol,nrow)
          if(nsymb.eq.1)then
           call printa(bc(i4od(is)),nocc,0,nocc,0,nocc,0,nocc,0,
     $          bc(ibcoff))
          end if
         else
          do i=0,nrow*ncol-1                                            5d31s22
           bc(ibcoff+i)=bc(i4od(is)+i)                                  5d31s22
          end do                                                        5d31s22
          call dws_gsumf(bc(ibcoff),nrow*ncol)                          5d31s22
          call prntm2(bc(ibcoff),nrow,ncol,nrow)                        5d31s22
         end if
         itmp=ibcoff
         ibcoff=itmp+nrow*ncol
         call enough('parajkfromhd0.  5',bc,ibc)
         do i=0,nrow*ncol-1
          bc(itmp+i)=bc(i4od(is)+i)
         end do
         write(6,*)('gsumf2 '),nrow,ncol,nrow*ncol
         call dws_gsumf(bc(itmp),nrow*ncol)
        end if                                                           5d4s22
        ibcoff=itmp
       end do
      end if                                                            3d21s16
      if(nsblkkder.lt.0)return                                          8d4s16
c
c     for j and k.                                                      8d4s16
c     since contributions for a single matrix element are comming from  8d4s16
c     all over the place, need to store all and do a global sum at end  8d4s16
c
   12 format('integral type ',4i2,5x,i5)                                3d4s13
      do isb=1,nsblkkder                                                8d4s16
       nrow=nocc(isblkkder(1,isb))*nocc(isblkkder(2,isb))               8d22s16
       ncol=nbasdwsc(isblkkder(3,isb))*nbasdwsc(isblkkder(4,isb))       8d4s16
       nall=nrow*ncol                                                   8d4s16
       if(nall.gt.0)then                                                11d30s16
       call ilimts(nvirtc(isblkkder(3,isb)),nvirtc(isblkkder(4,isb)),   8d22s16
     $      mynprocg,mynowprog,il,ih,i1s,i1e,is,i2e)                    8d22s16
       jmatd(isb)=ibcoff
       kmatd(isb)=jmatd(isb)+nrow*(ih+1-il)                             8d22s16
       ibcoff=kmatd(isb)+nrow*(ih+1-il)                                 8d22s16
       jmatda=ibcoff
       ibcoff=jmatda+nall
       call enough('parajkfromhd0.  6',bc,ibc)
       do i=0,nall-1                                                    8d4s16
        bc(jmatda+i)=0d0                                                8d4s16
       end do                                                           8d4s16
c
c     xor contribution to j'=(o'o|vv): (oo|vv) and (vo|vv)
c     to go from a to a', use trans(symof(a)) a,a'
c
       isw=multh(isblkkder(1,isb),iprop)                                8d4s16
c
c     in this naming scheme, isblkkder(1,isb) is a' and isw is a.
c
       igotit=0                                                         8d22s16
       do is=1,nsdlk
        if(isblk(1,is).eq.isw.and.isblk(2,is).eq.isblkkder(2,isb).and.  8d4s16
     $     isblk(3,is).eq.isblkkder(3,isb).and.
     $       isblk(4,is).eq.isblkkder(4,isb))then
         igotit=1                                                       8d22s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  8d4s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d19s16
         i10=i1s                                                        8d19s16
         i1n=nvirtc(isblk(3,is))                                        8d19s16
         ii=jmats(is)                                                   8d19s16
         do i2=i2s,i2e                                                  8d19s16
          if(i2.eq.i2e)i1n=i1e                                          8d19s16
          i2m=i2-1                                                      8d19s16
          do i1=i10,i1n                                                 8d19s16
           i1m=i1-1                                                     8d19s16
           if(isblk(1,is).eq.isblk(2,is))then                           8d19s16
            do i3=0,nocc(isblk(1,is))-1                                 8d19s16
             do i4=0,i3-1                                               8d19s16
              do im=0,nocc(isblkkder(1,isb))-1                          8d19s16
               iad1=itrans(isw)+i4+nbasdwsc(isw)*im                         8d19s16
               iad2=jmatda+im+nocc(isblkkder(1,isb))*(i3                8d19s16
     $              +nocc(isblk(2,is))*(i1m+nvirtc(isblk(3,is))*i2m))   8d19s16
               bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                        8d19s16
               iad1=itrans(isw)+i3+nbasdwsc(isw)*im                         8d19s16
               iad2=jmatda+im+nocc(isblkkder(1,isb))*(i4+               8d19s16
     $              nocc(isblk(2,is))*(i1m+nvirtc(isblk(3,is))*i2m))    8d19s16
               bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                        8d19s16
              end do                                                    8d19s16
              ii=ii+1                                                   8d19s16
             end do                                                     8d19s16
             do im=0,nocc(isblkkder(1,isb))-1                           8d19s16
              iad1=itrans(isw)+i3+nbasdwsc(isw)*im                          8d19s16
              iad2=jmatda+im+nocc(isblkkder(1,isb))*(i3+                8d19s16
     $             nocc(isblk(2,is))*(i1m+nvirtc(isblk(3,is))*i2m))     8d19s16
              bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                         8d19s16
             end do                                                     8d19s16
             ii=ii+1                                                    8d19s16
            end do                                                      8d19s16
           else                                                         8d19s16
            do i3=0,nocc(isblk(2,is))-1                                 8d19s16
             do i4=0,nocc(isblk(1,is))-1                                8d19s16
              do im=0,nocc(isblkkder(1,isb))-1                          8d19s16
               iad1=itrans(isw)+i4+nbasdwsc(isw)*im                         8d19s16
               iad2=jmatda+im+nocc(isblkkder(1,isb))*(i3+               8d19s16
     $              nocc(isblk(2,is))*(i1m+nvirtc(isblk(3,is))*i2m))    8d19s16
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                        8d19s16
              end do                                                    8d19s16
              ii=ii+1                                                   8d19s16
             end do                                                     8d19s16
            end do                                                      8d19s16
           end if                                                       8d19s16
          end do                                                        8d19s16
          i10=1                                                         8d19s16
         end do                                                         8d19s16
         go to 300                                                      8d22s16
        else if(isblk(2,is).eq.isw.and.isblk(1,is).eq.isblkkder(2,isb)  8d22s16
     $        .and.isblk(3,is).eq.isblkkder(3,isb).and.                 8d22s16
     $       isblk(4,is).eq.isblkkder(4,isb))then                       8d22s16
         igotit=1                                                       8d22s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  8d4s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d19s16
         i10=i1s                                                        8d19s16
         i1n=nvirtc(isblk(3,is))                                        8d19s16
         ii=jmats(is)                                                   8d19s16
         do i2=i2s,i2e                                                  8d19s16
          if(i2.eq.i2e)i1n=i1e                                          8d19s16
          i2m=i2-1                                                      8d19s16
          do i1=i10,i1n                                                 8d19s16
           i1m=i1-1                                                     8d19s16
           do i3=0,nocc(isblk(2,is))-1                                  8d22s16
            do i4=0,nocc(isblk(1,is))-1                                 8d22s16
             do im=0,nocc(isblkkder(1,isb))-1                           8d22s16
              iad1=itrans(isw)+i3+nbasdwsc(isw)*im                      8d22s16
              iad2=jmatda+im+nocc(isblkkder(1,isb))*(i4+                8d22s16
     $             nocc(isblk(1,is))*(i1m+nvirtc(isblk(3,is))*i2m))     8d22s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d22s16
             end do                                                     8d22s16
             ii=ii+1                                                    8d22s16
            end do                                                      8d22s16
           end do                                                       8d22s16
          end do                                                        8d19s16
          i10=1                                                         8d19s16
         end do                                                         8d19s16
         go to 300                                                      8d22s16
        else if(isblk(1,is).eq.isw.and.isblk(2,is).eq.isblkkder(2,isb)  8d22s16
     $        .and.isblk(4,is).eq.isblkkder(3,isb).and.                 8d22s16
     $       isblk(3,is).eq.isblkkder(4,isb))then
         igotit=1                                                       8d22s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  8d4s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d19s16
         i10=i1s                                                        8d19s16
         i1n=nvirtc(isblk(3,is))                                        8d19s16
         ii=jmats(is)                                                   8d19s16
         do i2=i2s,i2e                                                  8d19s16
          if(i2.eq.i2e)i1n=i1e                                          8d19s16
          i2m=i2-1                                                      8d19s16
          do i1=i10,i1n                                                 8d19s16
           i1m=i1-1                                                     8d19s16
           do i3=0,nocc(isblk(2,is))-1                                  8d22s16
            do i4=0,nocc(isblk(1,is))-1                                 8d22s16
             do im=0,nocc(isblkkder(1,isb))-1                           8d22s16
              iad1=itrans(isw)+i4+nbasdwsc(isw)*im                      8d22s16
              iad2=jmatda+im+nocc(isblkkder(1,isb))*(i3+                8d22s16
     $             nocc(isblk(2,is))*(i2m+nvirtc(isblk(4,is))*i1m))     8d22s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d22s16
             end do                                                     8d22s16
             ii=ii+1                                                    8d22s16
            end do                                                      8d22s16
           end do                                                       8d22s16
          end do                                                        8d19s16
          i10=1                                                         8d19s16
         end do                                                         8d19s16
         go to 300                                                      8d22s16
        else if(isblk(2,is).eq.isw.and.isblk(1,is).eq.isblkkder(2,isb)  8d22s16
     $        .and.isblk(4,is).eq.isblkkder(3,isb).and.                 8d22s16
     $       isblk(3,is).eq.isblkkder(4,isb))then
         igotit=1                                                       8d22s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  8d4s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d19s16
         i10=i1s                                                        8d19s16
         i1n=nvirtc(isblk(3,is))                                        8d19s16
         ii=jmats(is)                                                   8d19s16
         do i2=i2s,i2e                                                  8d19s16
          if(i2.eq.i2e)i1n=i1e                                          8d19s16
          i2m=i2-1                                                      8d19s16
          do i1=i10,i1n                                                 8d19s16
           i1m=i1-1                                                     8d19s16
           do i3=0,nocc(isblk(2,is))-1                                  8d22s16
            do i4=0,nocc(isblk(1,is))-1                                 8d22s16
             do im=0,nocc(isblkkder(1,isb))-1                           8d22s16
              iad1=itrans(isw)+i3+nbasdwsc(isw)*im                      8d22s16
              iad2=jmatda+im+nocc(isblkkder(1,isb))*(i4+                8d22s16
     $             nocc(isblk(1,is))*(i2m+nvirtc(isblk(4,is))*i1m))     8d22s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d22s16
             end do                                                     8d22s16
             ii=ii+1                                                    8d22s16
            end do                                                      8d22s16
           end do                                                       8d22s16
          end do                                                        8d19s16
          i10=1                                                         8d19s16
         end do                                                         8d19s16
         go to 300                                                      8d22s16
        end if
       end do
       write(6,*)('could not find jmats for (o''o|vv)'),isw,
     $      (isblkkder(j,isb),j=2,4)
       write(6,*)('got '),isblkkder(1,isb),iprop
       do i=1,nsdlk
        write(6,12)(isblk(j,i),j=1,4)
       end do
       call dws_sync
       call dws_finalize
       stop
  300  continue                                                         8d22s16
       if(nocc(isblkkder(2,isb)).eq.0)go to 301                         8d24s16
       do is=1,nsdlk1                                                   8d19s16
        if(isblk1(4,is).eq.isw.and.isblk1(3,is).eq.isblkkder(2,isb).and.8d19s16
     $     isblk1(1,is).eq.isblkkder(3,isb).and.                        8d19s16
     $       isblk1(2,is).eq.isblkkder(4,isb))then                      8d19s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  8d19s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d19s16
         i10=i1s                                                        8d19s16
         i1n=nocc(isblk1(3,is))                                         8d19s16
         ii=i3x(is)                                                     8d19s16
         do i2=i2s,i2e                                                  8d19s16
          i2p=i2-1+nocc(isblk1(4,is))                                   8d29s16
          if(i2.eq.i2e)i1n=i1e                                          8d19s16
          do i1=i10,i1n                                                 8d19s16
           i1m=i1-1                                                     8d19s16
           if(isblk1(1,is).eq.isblk1(2,is))then                           8d19s16
            do i3=0,nvirtc(isblk1(1,is))-1                                 8d19s16
             do i4=0,i3-1                                               8d19s16
              do im=0,nocc(isblkkder(1,isb))-1                          8d19s16
               iad1=itrans(isw)+i2p+nbasdwsc(isw)*im                         8d19s16
               iad2=jmatda+im+nocc(isblkkder(1,isb))*(i1m                8d19s16
     $              +nocc(isblk1(3,is))*(i4+nvirtc(isblk1(1,is))*i3))   8d19s16
               bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                        8d19s16
               iad2=jmatda+im+nocc(isblkkder(1,isb))*(i1m+               8d19s16
     $              nocc(isblk1(3,is))*(i3+nvirtc(isblk1(1,is))*i4))    8d19s16
               bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                        8d19s16
              end do                                                    8d19s16
              ii=ii+1                                                   8d19s16
             end do                                                     8d19s16
             do im=0,nocc(isblkkder(1,isb))-1                           8d19s16
              iad1=itrans(isw)+i2p+nbasdwsc(isw)*im                          8d19s16
              iad2=jmatda+im+nocc(isblkkder(1,isb))*(i1m+                8d19s16
     $             nocc(isblk1(3,is))*(i3+nvirtc(isblk1(1,is))*i3))     8d19s16
              bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                         8d19s16
             end do                                                     8d19s16
             ii=ii+1                                                    8d19s16
            end do                                                      8d19s16
           else                                                         8d19s16
            do i3=0,nvirtc(isblk1(2,is))-1                                 8d19s16
             do i4=0,nvirtc(isblk1(1,is))-1                                8d19s16
              do im=0,nocc(isblkkder(1,isb))-1                          8d19s16
               iad1=itrans(isw)+i2p+nbasdwsc(isw)*im                         8d19s16
               iad2=jmatda+im+nocc(isblkkder(1,isb))*(i1m+               8d19s16
     $              nocc(isblk1(3,is))*(i4+nvirtc(isblk1(1,is))*i3))    8d19s16
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                        8d19s16
              end do                                                    8d19s16
              ii=ii+1                                                   8d19s16
             end do                                                     8d19s16
            end do                                                      8d19s16
           end if                                                       8d19s16
          end do                                                        8d19s16
          i10=1                                                         8d19s16
         end do                                                         8d19s16
         go to 301                                                      8d22s16
        else if(isblk1(4,is).eq.isw.and.isblk1(3,is).eq.isblkkder(2,isb)8d19s16
     $     .and.isblk1(2,is).eq.isblkkder(3,isb).and.                   8d22s16
     $       isblk1(1,is).eq.isblkkder(4,isb))then                      8d22s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  8d19s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d19s16
         i10=i1s                                                        8d19s16
         i1n=nocc(isblk1(3,is))                                         8d19s16
         ii=i3x(is)                                                     8d19s16
         do i2=i2s,i2e                                                  8d19s16
          i2p=i2-1+nocc(isblk1(4,is))                                   8d29s16
          if(i2.eq.i2e)i1n=i1e                                          8d19s16
          do i1=i10,i1n                                                 8d19s16
           i1m=i1-1                                                     8d19s16
           do i3=0,nvirtc(isblk1(2,is))-1                               8d22s16
            do i4=0,nvirtc(isblk1(1,is))-1                              8d22s16
             do im=0,nocc(isblkkder(1,isb))-1                           8d22s16
              iad1=itrans(isw)+i2p+nbasdwsc(isw)*im                     8d22s16
              iad2=jmatda+im+nocc(isblkkder(1,isb))*(i1m+               8d22s16
     $             nocc(isblk1(3,is))*(i3+nvirtc(isblk1(2,is))*i4))     8d22s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d22s16
             end do                                                     8d22s16
             ii=ii+1                                                    8d22s16
            end do                                                      8d22s16
           end do                                                       8d22s16
          end do                                                        8d19s16
          i10=1                                                         8d19s16
         end do                                                         8d19s16
         go to 301                                                      8d22s16
        end if
       end do                                                           8d19s16
       if(min(nbasdwsc(isblkkder(3,isb)),nbasdwsc(isblkkder(4,isb)),    11d28s22
     $      nbasdwsc(isw)).gt.0)then                                    11d28s22
        write(6,*)('could not find x3 for (o''o|vv)'),isw,
     $      (isblkkder(j,isb),j=2,4)
        write(6,*)('want: '),isblkkder(3,isb),isblkkder(4,isb),
     $      isblkkder(2,isb),isw
        write(6,*)('got ')
        do i=1,nsdlk1
         write(6,12)(isblk1(j,i),j=1,4)
        end do
        call dws_sync
        call dws_finalize
        stop
       end if                                                           11d28s22
  301  continue                                                         8d22s16
c
c     xor contribution to j'=(oo'|vv): (oo|vv) and (ov|vv)
c
       isw=multh(isblkkder(2,isb),iprop)                                8d19s16
       if(nocc(isw).eq.0)go to 302                                      8d24s16
       do is=1,nsdlk
        if(isblk(1,is).eq.isblkkder(1,isb).and.isblk(2,is).eq.isw.and.  8d19s16
     $     isblk(3,is).eq.isblkkder(3,isb).and.
     $       isblk(4,is).eq.isblkkder(4,isb))then
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  8d4s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d19s16
         i10=i1s                                                        8d19s16
         i1n=nvirtc(isblk(3,is))                                        8d19s16
         ii=jmats(is)                                                   8d19s16
         do i2=i2s,i2e                                                  8d19s16
          if(i2.eq.i2e)i1n=i1e                                          8d19s16
          i2m=i2-1                                                      8d19s16
          do i1=i10,i1n                                                 8d19s16
           i1m=i1-1                                                     8d19s16
           if(isblk(1,is).eq.isblk(2,is))then                           8d19s16
            do i3=0,nocc(isblk(1,is))-1                                 8d19s16
             do i4=0,i3-1                                               8d19s16
              do im=0,nocc(isblkkder(2,isb))-1                          8d19s16
               iad1=itrans(isw)+i4+nbasdwsc(isw)*im                         8d19s16
               iad2=jmatda+i3+nocc(isblk(1,is))*(im                     8d25s16
     $            +nocc(isblkkder(2,isb))*(i1m+nvirtc(isblk(3,is))*i2m))8d25s16
               bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                        8d19s16
               iad1=itrans(isw)+i3+nbasdwsc(isw)*im                         8d19s16
               iad2=jmatda+i4+nocc(isblk(1,is))*(im+                    8d25s16
     $             nocc(isblkkder(2,isb))*(i1m+nvirtc(isblk(3,is))*i2m))8d25s16
               bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                        8d19s16
              end do                                                    8d19s16
              ii=ii+1                                                   8d19s16
             end do                                                     8d19s16
             do im=0,nocc(isblkkder(2,isb))-1                           8d19s16
              iad1=itrans(isw)+i3+nbasdwsc(isw)*im                          8d19s16
              iad2=jmatda+i3+nocc(isblk(1,is))*(im+                     8d25s16
     $             nocc(isblkkder(2,isb))*(i1m+nvirtc(isblk(3,is))*i2m))8d25s16
              bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                         8d19s16
             end do                                                     8d19s16
             ii=ii+1                                                    8d19s16
            end do                                                      8d19s16
           else                                                         8d19s16
            do i3=0,nocc(isblk(2,is))-1                                 8d19s16
             do i4=0,nocc(isblk(1,is))-1                                8d19s16
              do im=0,nocc(isblkkder(2,isb))-1                          8d19s16
               iad1=itrans(isw)+i3+nbasdwsc(isw)*im                         8d19s16
               iad2=jmatda+i4+nocc(isblk(1,is))*(im+                    8d25s16
     $             nocc(isblkkder(2,isb))*(i1m+nvirtc(isblk(3,is))*i2m))8d25s16
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                        8d19s16
              end do                                                    8d19s16
              ii=ii+1                                                   8d19s16
             end do                                                     8d19s16
            end do                                                      8d19s16
           end if                                                       8d19s16
          end do                                                        8d19s16
          i10=1                                                         8d19s16
         end do                                                         8d19s16
         go to 302                                                      8d22s16
        else if(isblk(2,is).eq.isblkkder(1,isb).and.                    8d22s16
     $        isblk(1,is).eq.isw.and.                                   8d22s16
     $     isblk(3,is).eq.isblkkder(3,isb).and.                         8d22s16
     $       isblk(4,is).eq.isblkkder(4,isb))then                       8d22s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  8d4s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d19s16
         i10=i1s                                                        8d19s16
         i1n=nvirtc(isblk(3,is))                                        8d19s16
         ii=jmats(is)                                                   8d19s16
         do i2=i2s,i2e                                                  8d19s16
          if(i2.eq.i2e)i1n=i1e                                          8d19s16
          i2m=i2-1                                                      8d19s16
          do i1=i10,i1n                                                 8d19s16
           i1m=i1-1                                                     8d19s16
           do i3=0,nocc(isblk(2,is))-1                                  8d22s16
            do i4=0,nocc(isblk(1,is))-1                                 8d22s16
             do im=0,nocc(isblkkder(2,isb))-1                           8d22s16
              iad1=itrans(isw)+i4+nbasdwsc(isw)*im                      8d22s16
              iad2=jmatda+i3+nocc(isblk(2,is))*(im+                     8d25s16
     $             nocc(isblkkder(2,isb))*(i1m+nvirtc(isblk(3,is))*i2m))8d25s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d22s16
             end do                                                     8d22s16
             ii=ii+1                                                    8d22s16
            end do                                                      8d22s16
           end do                                                       8d22s16
          end do                                                        8d19s16
          i10=1                                                         8d19s16
         end do                                                         8d19s16
         go to 302                                                      8d22s16
        else if(isblk(2,is).eq.isblkkder(1,isb).and.                    8d22s16
     $        isblk(1,is).eq.isw.and.                                   8d22s16
     $     isblk(4,is).eq.isblkkder(3,isb).and.                         8d22s16
     $       isblk(3,is).eq.isblkkder(4,isb))then                       8d22s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  8d4s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d19s16
         i10=i1s                                                        8d19s16
         i1n=nvirtc(isblk(3,is))                                        8d19s16
         ii=jmats(is)                                                   8d19s16
         do i2=i2s,i2e                                                  8d19s16
          if(i2.eq.i2e)i1n=i1e                                          8d19s16
          i2m=i2-1                                                      8d19s16
          do i1=i10,i1n                                                 8d19s16
           i1m=i1-1                                                     8d19s16
           do i3=0,nocc(isblk(2,is))-1                                  8d22s16
            do i4=0,nocc(isblk(1,is))-1                                 8d22s16
             do im=0,nocc(isblkkder(2,isb))-1                           8d22s16
              iad1=itrans(isw)+i4+nbasdwsc(isw)*im                      8d22s16
              iad2=jmatda+i3+nocc(isblk(2,is))*(im+                     8d25s16
     $             nocc(isblkkder(2,isb))*(i2m+nvirtc(isblk(4,is))*i1m))8d25s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d22s16
             end do                                                     8d22s16
             ii=ii+1                                                    8d22s16
            end do                                                      8d22s16
           end do                                                       8d22s16
          end do                                                        8d19s16
          i10=1                                                         8d19s16
         end do                                                         8d19s16
         go to 302                                                      8d22s16
        else if(isblk(1,is).eq.isblkkder(1,isb).and.                    8d22s16
     $        isblk(2,is).eq.isw.and.                                   8d22s16
     $     isblk(4,is).eq.isblkkder(3,isb).and.                         8d22s16
     $       isblk(3,is).eq.isblkkder(4,isb))then                       8d22s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  8d4s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d19s16
         i10=i1s                                                        8d19s16
         i1n=nvirtc(isblk(3,is))                                        8d19s16
         ii=jmats(is)                                                   8d19s16
         do i2=i2s,i2e                                                  8d19s16
          if(i2.eq.i2e)i1n=i1e                                          8d19s16
          i2m=i2-1                                                      8d19s16
          do i1=i10,i1n                                                 8d19s16
           i1m=i1-1                                                     8d19s16
           do i3=0,nocc(isblk(2,is))-1                                  8d22s16
            do i4=0,nocc(isblk(1,is))-1                                 8d22s16
             do im=0,nocc(isblkkder(2,isb))-1                           8d22s16
              iad1=itrans(isw)+i3+nbasdwsc(isw)*im                      8d22s16
              iad2=jmatda+i4+nocc(isblk(1,is))*(im+                     8d25s16
     $             nocc(isblkkder(2,isb))*(i2m+nvirtc(isblk(4,is))*i1m))8d25s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d22s16
             end do                                                     8d22s16
             ii=ii+1                                                    8d22s16
            end do                                                      8d22s16
           end do                                                       8d22s16
          end do                                                        8d19s16
          i10=1                                                         8d19s16
         end do                                                         8d19s16
         go to 302                                                      8d22s16
        end if
       end do
       write(6,*)('need jhelp for (oo''|vv) '),isblkkder(1,isb),isw,
     $      (isblkkder(j,isb),j=3,4)
       do i=1,nsdlk
        write(6,12)(isblk(j,i),j=1,4)
       end do
       call dws_sync
       call dws_finalize
       stop
  302  continue                                                         8d22s16
       if(nocc(isblkkder(1,isb)).eq.0)go to 303                         8d24s16
       do is=1,nsdlk1                                                   8d19s16
        if(isblk1(4,is).eq.isw.and.isblk1(3,is).eq.isblkkder(1,isb).and.8d19s16
     $     isblk1(1,is).eq.isblkkder(3,isb).and.                        8d19s16
     $       isblk1(2,is).eq.isblkkder(4,isb))then                      8d19s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  8d19s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d19s16
         i10=i1s                                                        8d19s16
         i1n=nocc(isblk1(3,is))                                         8d19s16
         ii=i3x(is)                                                     8d19s16
         do i2=i2s,i2e                                                  8d19s16
          i2p=i2-1+nocc(isblk1(4,is))                                     8d19s16
          if(i2.eq.i2e)i1n=i1e                                          8d19s16
          do i1=i10,i1n                                                 8d19s16
           i1m=i1-1
           if(isblk1(1,is).eq.isblk1(2,is))then                           8d19s16
            do i3=0,nvirtc(isblk1(1,is))-1                                 8d19s16
             do i4=0,i3-1                                               8d19s16
              do im=0,nocc(isblkkder(2,isb))-1                          8d19s16
               iad1=itrans(isw)+i2p+nbasdwsc(isw)*im                         8d19s16
               iad2=jmatda+i1m+nocc(isblk1(3,is))*(im                   8d19s16
     $             +nocc(isblkkder(2,isb))*(i4+nvirtc(isblk1(1,is))*i3))8d25s16
               bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                        8d19s16
               iad2=jmatda+i1m+nocc(isblk1(3,is))*(im+                  8d19s16
     $              nocc(isblkkder(2,isb))*(i3+nvirtc(isblk1(1,is))*i4))8d25s16
               bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                        8d19s16
              end do                                                    8d19s16
              ii=ii+1                                                   8d19s16
             end do                                                     8d19s16
             do im=0,nocc(isblkkder(2,isb))-1                           8d19s16
              iad1=itrans(isw)+i2p+nbasdwsc(isw)*im                          8d19s16
              iad2=jmatda+i1m+nocc(isblk1(3,is))*(im+                   8d19s16
     $             nocc(isblkkder(2,isb))*(i3+nvirtc(isblk1(1,is))*i3)) 8d25s16
              bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                         8d19s16
             end do                                                     8d19s16
             ii=ii+1                                                    8d19s16
            end do                                                      8d19s16
           else                                                         8d19s16
            do i3=0,nvirtc(isblk1(2,is))-1                                 8d19s16
             do i4=0,nvirtc(isblk1(1,is))-1                                8d19s16
              do im=0,nocc(isblkkder(2,isb))-1                          8d19s16
               iad1=itrans(isw)+i2p+nbasdwsc(isw)*im                         8d19s16
               iad2=jmatda+i1m+nocc(isblk1(3,is))*(im+                  8d19s16
     $              nocc(isblkkder(2,isb))*(i4+nvirtc(isblk1(1,is))*i3))8d25s16
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                        8d19s16
              end do                                                    8d19s16
              ii=ii+1                                                   8d19s16
             end do                                                     8d19s16
            end do                                                      8d19s16
           end if                                                       8d19s16
          end do                                                        8d19s16
          i10=1                                                         8d19s16
         end do                                                         8d19s16
         go to 303                                                      8d22s16
        else if(isblk1(4,is).eq.isw.and.isblk1(3,is).eq.isblkkder(1,isb)8d22s16
     $     .and.isblk1(2,is).eq.isblkkder(3,isb).and.                   8d22s16
     $       isblk1(1,is).eq.isblkkder(4,isb))then                      8d22s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  8d19s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d19s16
         i10=i1s                                                        8d19s16
         i1n=nocc(isblk1(3,is))                                         8d19s16
         ii=i3x(is)                                                     8d19s16
         do i2=i2s,i2e                                                  8d19s16
          i2p=i2-1+nocc(isblk1(4,is))                                     8d19s16
          if(i2.eq.i2e)i1n=i1e                                          8d19s16
          do i1=i10,i1n                                                 8d19s16
           i1m=i1-1
           do i3=0,nvirtc(isblk1(2,is))-1                                 8d19s16
            do i4=0,nvirtc(isblk1(1,is))-1                                8d19s16
             do im=0,nocc(isblkkder(2,isb))-1                           8d22s16
              iad1=itrans(isw)+i2p+nbasdwsc(isw)*im                     8d22s16
              iad2=jmatda+i1m+nocc(isblk1(3,is))*(im+                   8d22s16
     $             nocc(isblkkder(2,isb))*(i3+nvirtc(isblk1(2,is))*i4)) 8d25s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d22s16
             end do                                                     8d22s16
             ii=ii+1                                                    8d22s16
            end do                                                      8d22s16
           end do                                                       8d22s16
          end do                                                        8d19s16
          i10=1                                                         8d19s16
         end do                                                         8d19s16
         go to 303                                                      8d22s16
        end if
       end do                                                           8d19s16
       if(min(nocc(isblkkder(1,isb)),nvirtc(isw)).gt.0)then             11d28s22
        write(6,*)('could not find i3x for (oo''|vv) '),
     $      isblkkder(1,isb),isw,(isblkkder(j,isb),j=3,4)
        write(6,*)('want: '),isblkkder(3,isb),isblkkder(4,isb),
     $      isblkkder(1,isb),isw
        write(6,*)('got: ')
        do i=1,nsdlk1
         write(6,12)(isblk1(j,i),j=1,4)
        end do
        call dws_sync
        call dws_finalize
        stop
       end if                                                           11d28s22
  303  continue
c
c     xor contribution to j'=(oo|v'v): (oo|ov) and (oo|vv)
c
       isw=multh(isblkkder(3,isb),iprop)                                8d19s16
       if(nocc(isw).eq.0)go to 304                                      8d24s16
       do is=1,nsdlk1                                                   8d19s16
        if(isblk1(3,is).eq.isw.and.isblk1(4,is).eq.isblkkder(4,isb).and.8d19s16
     $     isblk1(1,is).eq.isblkkder(1,isb).and.                        8d19s16
     $       isblk1(2,is).eq.isblkkder(2,isb))then                      8d19s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  8d19s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d19s16
         i10=i1s                                                        8d19s16
         i1n=nocc(isblk1(3,is))                                         8d19s16
         ii=ionex(is)                                                   8d19s16
         do i2=i2s,i2e                                                  8d19s16
          if(i2.eq.i2e)i1n=i1e                                          8d19s16
          i2m=i2-1                                                      8d19s16
          do i1=i10,i1n                                                 8d19s16
           i1m=i1-1
           if(isblk1(1,is).eq.isblk1(2,is))then                           8d19s16
            do i3=0,nocc(isblk1(1,is))-1                                 8d19s16
             do i4=0,i3-1                                               8d19s16
              do im=0,nvirtc(isblkkder(3,isb))-1                        8d19s16
               imp=im+nocc(isblkkder(3,isb))                            8d19s16
               iad1=itrans(isw)+i1m+nbasdwsc(isw)*imp                   8d19s16
               iad2=jmatda+i4+nocc(isblk1(1,is))*(i3                    8d19s16
     $            +nocc(isblk1(2,is))*(im+nvirtc(isblkkder(3,isb))*i2m))            8d19s16
               bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                        8d19s16
               iad2=jmatda+i3+nocc(isblk1(1,is))*(i4+                   8d19s16
     $             nocc(isblk1(2,is))*(im+nvirtc(isblkkder(3,isb))*i2m))    8d19s16
               bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                        8d19s16
              end do                                                    8d19s16
              ii=ii+1                                                   8d19s16
             end do                                                     8d19s16
             do im=0,nvirtc(isblkkder(3,isb))-1                           8d19s16
              imp=im+nocc(isblkkder(3,isb))                             8d19s16
              iad1=itrans(isw)+i1m+nbasdwsc(isw)*imp                         8d19s16
              iad2=jmatda+i3+nocc(isblk1(1,is))*(i3+                    8d19s16
     $             nocc(isblk1(1,is))*(im+nvirtc(isblkkder(3,isb))*i2m))     8d19s16
              bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                         8d19s16
             end do                                                     8d19s16
             ii=ii+1                                                    8d19s16
            end do                                                      8d19s16
           else                                                         8d19s16
            do i3=0,nocc(isblk1(2,is))-1                                 8d19s16
             do i4=0,nocc(isblk1(1,is))-1                                8d19s16
              do im=0,nvirtc(isblkkder(3,isb))-1                        8d19s16
               imp=im+nocc(isblkkder(3,isb))                            8d19s16
               iad1=itrans(isw)+i1m+nbasdwsc(isw)*imp                         8d19s16
               iad2=jmatda+i4+nocc(isblk1(1,is))*(i3+                   8d19s16
     $             nocc(isblk1(2,is))*(im+nvirtc(isblkkder(3,isb))*i2m))    8d19s16
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                        8d19s16
              end do                                                    8d19s16
              ii=ii+1                                                   8d19s16
             end do                                                     8d19s16
            end do                                                      8d19s16
           end if                                                       8d19s16
          end do                                                        8d19s16
          i10=1                                                         8d19s16
         end do                                                         8d19s16
         go to 304                                                      8d22s16
        else if(isblk1(3,is).eq.isw.and.isblk1(4,is).eq.isblkkder(4,isb)8d24s16
     $     .and.isblk1(2,is).eq.isblkkder(1,isb).and.                   8d24s16
     $       isblk1(1,is).eq.isblkkder(2,isb))then                      8d19s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  8d24s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         8d24s16
         ii=ionex(is)                                                   8d24s16
         do i2=i2s,i2e                                                  8d24s16
          if(i2.eq.i2e)i1n=i1e                                          8d24s16
          i2m=i2-1                                                      8d24s16
          do i1=i10,i1n                                                 8d24s16
           i1m=i1-1                                                     8d24s16
           do i3=0,nocc(isblk1(2,is))-1                                 8d24s16
            do i4=0,nocc(isblk1(1,is))-1                                8d24s16
             do im=0,nvirtc(isblkkder(3,isb))-1                         8d24s16
              imp=im+nocc(isblkkder(3,isb))                             8d24s16
              iad1=itrans(isw)+i1m+nbasdwsc(isw)*imp                    8d24s16
              iad2=jmatda+i3+nocc(isblk1(2,is))*(i4+                    8d24s16
     $             nocc(isblk1(1,is))*(im+nvirtc(isblkkder(3,isb))*i2m))11d30s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d24s16
             end do                                                     8d24s16
             ii=ii+1                                                    8d24s16
            end do                                                      8d24s16
           end do                                                       8d24s16
          end do                                                        8d24s16
          i10=1                                                         8d24s16
         end do                                                         8d24s16
         go to 304                                                      8d24s16
        end if
       end do                                                           8d19s16
       write(6,*)('could not find onex for (oo|v''v)'),(isblkkder(j,isb)
     $      ,j=1,2),isw,isblkkder(4,isb)
       call dws_sync
       call dws_finalize
       stop
  304  continue                                                         8d22s16
       if(nocc(isblkkder(1,isb))*nocc(isblkkder(2,isb)).eq.0)go to 305  8d24s16
       do is=1,nsdlk
        if(isblk(1,is).eq.isblkkder(1,isb).and.                         8d19s16
     $       isblk(2,is).eq.isblkkder(2,isb).and.isblk(3,is).eq.isw.and.8d19s16
     $       isblk(4,is).eq.isblkkder(4,isb))then                       8d19s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  8d4s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d19s16
         i10=i1s                                                        8d19s16
         i1n=nvirtc(isblk(3,is))                                        8d19s16
         ii=jmats(is)                                                   8d19s16
         do i2=i2s,i2e                                                  8d19s16
          if(i2.eq.i2e)i1n=i1e                                          8d19s16
          i2m=i2-1                                                      8d19s16
          do i1=i10,i1n                                                 8d19s16
           i1p=i1-1+nocc(isblk(3,is))                                   8d19s16
           if(isblk(1,is).eq.isblk(2,is))then                           8d19s16
            do i3=0,nocc(isblk(1,is))-1                                 8d19s16
             do i4=0,i3-1                                               8d19s16
              do im=0,nvirtc(isblkkder(3,isb))-1                        8d19s16
               imp=im+nocc(isblkkder(3,isb))                            8d19s16
               iad1=itrans(isw)+i1p+nbasdwsc(isw)*imp                   8d19s16
               iad2=jmatda+i3+nocc(isblk(1,is))*(i4                     8d19s16
     $             +nocc(isblk(1,is))*(im+nvirtc(isblkkder(3,isb))*i2m))   8d19s16
               bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                        8d19s16
               iad2=jmatda+i4+nocc(isblk(1,is))*(i3+                    8d30s16
     $              nocc(isblk(1,is))*(im+nvirtc(isblkkder(3,isb))*i2m))8d19s16
               bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                        8d19s16
              end do                                                    8d19s16
              ii=ii+1                                                   8d19s16
             end do                                                     8d19s16
             do im=0,nvirtc(isblkkder(3,isb))-1                         8d19s16
              imp=im+nocc(isblkkder(3,isb))                             8d19s16
              iad1=itrans(isw)+i1p+nbasdwsc(isw)*imp                          8d19s16
              iad2=jmatda+i3+nocc(isblk(1,is))*(i3+                     8d19s16
     $             nocc(isblk(1,is))*(im+nvirtc(isblkkder(3,isb))*i2m))  8d19s16
              bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                         8d19s16
             end do                                                     8d19s16
             ii=ii+1                                                    8d19s16
            end do                                                      8d19s16
           else                                                         8d19s16
            do i3=0,nocc(isblk(2,is))-1                                 8d19s16
             do i4=0,nocc(isblk(1,is))-1                                8d19s16
              do im=0,nvirtc(isblkkder(3,isb))-1                          8d19s16
               imp=im+nocc(isblkkder(3,isb))
               iad1=itrans(isw)+i1p+nbasdwsc(isw)*imp                         8d19s16
               iad2=jmatda+i4+nocc(isblk(1,is))*(i3+                    8d19s16
     $              nocc(isblk(2,is))*(im+nvirtc(isblkkder(3,isb))*i2m))8d19s16
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                        8d19s16
              end do                                                    8d19s16
              ii=ii+1                                                   8d19s16
             end do                                                     8d19s16
            end do                                                      8d19s16
           end if                                                       8d19s16
          end do                                                        8d19s16
          i10=1                                                         8d19s16
         end do                                                         8d19s16
         go to 305                                                      8d22s16
        else if(isblk(2,is).eq.isblkkder(1,isb).and.                    8d24s16
     $       isblk(1,is).eq.isblkkder(2,isb).and.isblk(3,is).eq.isw.and.8d24s16
     $       isblk(4,is).eq.isblkkder(4,isb))then                       8d24s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  8d24s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        8d24s16
         ii=jmats(is)                                                   8d24s16
         do i2=i2s,i2e                                                  8d24s16
          if(i2.eq.i2e)i1n=i1e                                          8d24s16
          i2m=i2-1                                                      8d24s16
          do i1=i10,i1n                                                 8d24s16
           i1p=i1-1+nocc(isblk(3,is))                                   8d24s16
           do i3=0,nocc(isblk(2,is))-1                                  8d24s16
            do i4=0,nocc(isblk(1,is))-1                                 8d24s16
             do im=0,nvirtc(isblkkder(3,isb))-1                         8d24s16
              imp=im+nocc(isblkkder(3,isb))                             8d24s16
              iad1=itrans(isw)+i1p+nbasdwsc(isw)*imp                    8d24s16
              iad2=jmatda+i3+nocc(isblk(2,is))*(i4+                     8d24s16
     $             nocc(isblk(1,is))*(im+nvirtc(isblkkder(3,isb))*i2m)) 8d24s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d24s16
             end do                                                     8d24s16
             ii=ii+1                                                    8d24s16
            end do                                                      8d24s16
           end do                                                       8d24s16
          end do                                                        8d24s16
          i10=1                                                         8d24s16
         end do                                                         8d24s16
         go to 305                                                      8d24s16
        else if(isblk(2,is).eq.isblkkder(1,isb).and.                    8d24s16
     $       isblk(1,is).eq.isblkkder(2,isb).and.isblk(4,is).eq.isw.and.8d24s16
     $       isblk(3,is).eq.isblkkder(4,isb))then                       8d24s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  8d24s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        8d24s16
         ii=jmats(is)                                                   8d24s16
         do i2=i2s,i2e                                                  8d24s16
          if(i2.eq.i2e)i1n=i1e                                          8d24s16
          i2p=i2-1+nocc(isblk(4,is))                                    8d24s16
          do i1=i10,i1n                                                 8d24s16
           i1m=i1-1                                                     8d24s16
           do i3=0,nocc(isblk(2,is))-1                                  8d24s16
            do i4=0,nocc(isblk(1,is))-1                                 8d24s16
             do im=0,nvirtc(isblkkder(3,isb))-1                         8d24s16
              imp=im+nocc(isblkkder(3,isb))                             8d24s16
              iad1=itrans(isw)+i2p+nbasdwsc(isw)*imp                    8d24s16
              iad2=jmatda+i3+nocc(isblk(2,is))*(i4+                     8d24s16
     $             nocc(isblk(1,is))*(im+nvirtc(isblkkder(3,isb))*i1m)) 8d24s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d24s16
             end do                                                     8d24s16
             ii=ii+1                                                    8d24s16
            end do                                                      8d24s16
           end do                                                       8d24s16
          end do                                                        8d24s16
          i10=1                                                         8d24s16
         end do                                                         8d24s16
         go to 305                                                      8d24s16
        else if(isblk(1,is).eq.isblkkder(1,isb).and.                    8d24s16
     $       isblk(2,is).eq.isblkkder(2,isb).and.isblk(4,is).eq.isw.and.8d24s16
     $       isblk(3,is).eq.isblkkder(4,isb))then                       8d24s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  8d24s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        8d24s16
         ii=jmats(is)                                                   8d24s16
         do i2=i2s,i2e                                                  8d24s16
          if(i2.eq.i2e)i1n=i1e                                          8d24s16
          i2p=i2-1+nocc(isblk(4,is))                                    8d24s16
          do i1=i10,i1n                                                 8d24s16
           i1m=i1-1                                                     8d24s16
           do i3=0,nocc(isblk(2,is))-1                                  8d24s16
            do i4=0,nocc(isblk(1,is))-1                                 8d24s16
             do im=0,nvirtc(isblkkder(3,isb))-1                         8d24s16
              imp=im+nocc(isblkkder(3,isb))                             8d24s16
              iad1=itrans(isw)+i2p+nbasdwsc(isw)*imp                    8d24s16
              iad2=jmatda+i4+nocc(isblk(1,is))*(i3+                     8d24s16
     $             nocc(isblk(2,is))*(im+nvirtc(isblkkder(3,isb))*i1m)) 8d24s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d24s16
             end do                                                     8d24s16
             ii=ii+1                                                    8d24s16
            end do                                                      8d24s16
           end do                                                       8d24s16
          end do                                                        8d24s16
          i10=1                                                         8d24s16
         end do                                                         8d24s16
         go to 305                                                      8d24s16
        end if
       end do
       write(6,*)('need jhelp for (oo|v''v)'),(isblkkder(j,isb),j=1,2),
     $      isw,isblkkder(4,isb)
       call dws_sync
       call dws_finalize
       stop
  305  continue
c
c     xor contribution to j'=(oo|vv'): (oo|vo) and (oo|vv)
c
       isw=multh(isblkkder(4,isb),iprop)                                8d19s16
       if(nocc(isw).eq.0)go to 306                                      8d24s16
       do is=1,nsdlk1                                                   8d19s16
        if(isblk1(3,is).eq.isw.and.isblk1(4,is).eq.isblkkder(3,isb).and.8d19s16
     $     isblk1(1,is).eq.isblkkder(1,isb).and.                        8d19s16
     $       isblk1(2,is).eq.isblkkder(2,isb))then                      8d19s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  8d19s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d19s16
         i10=i1s                                                        8d19s16
         i1n=nocc(isblk1(3,is))                                         8d19s16
         ii=ionex(is)                                                   8d19s16
         do i2=i2s,i2e                                                  8d19s16
          if(i2.eq.i2e)i1n=i1e                                          8d19s16
          i2m=i2-1                                                      8d19s16
          do i1=i10,i1n                                                 8d19s16
           i1m=i1-1
           if(isblk1(1,is).eq.isblk1(2,is))then                           8d19s16
            do i3=0,nocc(isblk1(1,is))-1                                 8d19s16
             do i4=0,i3-1                                               8d19s16
              do im=0,nvirtc(isblkkder(4,isb))-1                        8d19s16
               imp=im+nocc(isblkkder(4,isb))                            8d19s16
               iad1=itrans(isw)+i1m+nbasdwsc(isw)*imp                   8d19s16
               iad2=jmatda+i4+nocc(isblk1(1,is))*(i3                    8d19s16
     $             +nocc(isblk1(2,is))*(i2m+nvirtc(isblk1(4,is))*im))   8d22s16
               bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                        8d19s16
               iad2=jmatda+i3+nocc(isblk1(1,is))*(i4+                   8d19s16
     $             nocc(isblk1(2,is))*(i2m+nvirtc(isblk1(4,is))*im))    8d19s16
               bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                        8d19s16
              end do                                                    8d19s16
              ii=ii+1                                                   8d19s16
             end do                                                     8d19s16
             do im=0,nvirtc(isblkkder(4,isb))-1                           8d19s16
              imp=im+nocc(isblkkder(4,isb))                             8d19s16
              iad1=itrans(isw)+i1m+nbasdwsc(isw)*imp                         8d19s16
              iad2=jmatda+i3+nocc(isblk1(1,is))*(i3+                    8d19s16
     $             nocc(isblk1(1,is))*(i2m+nvirtc(isblk1(4,is))*im))    8d25s16
              bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                         8d19s16
             end do                                                     8d19s16
             ii=ii+1                                                    8d19s16
            end do                                                      8d19s16
           else                                                         8d19s16
            do i3=0,nocc(isblk1(2,is))-1                                 8d19s16
             do i4=0,nocc(isblk1(1,is))-1                                8d19s16
              do im=0,nvirtc(isblkkder(4,isb))-1                        8d19s16
               imp=im+nocc(isblkkder(4,isb))                            8d19s16
               iad1=itrans(isw)+i1m+nbasdwsc(isw)*imp                         8d19s16
               iad2=jmatda+i4+nocc(isblk1(1,is))*(i3+                   8d19s16
     $              nocc(isblk1(2,is))*(i2m+nvirtc(isblk1(4,is))*im))   8d25s16
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                        8d19s16
              end do                                                    8d19s16
              ii=ii+1                                                   8d19s16
             end do                                                     8d19s16
            end do                                                      8d19s16
           end if                                                       8d19s16
          end do                                                        8d19s16
          i10=1                                                         8d19s16
         end do                                                         8d19s16
         go to 306                                                      8d22s16
        else if(isblk1(3,is).eq.isw.and.isblk1(4,is).eq.isblkkder(3,isb)8d24s16
     $     .and.isblk1(2,is).eq.isblkkder(1,isb).and.                   8d24s16
     $       isblk1(1,is).eq.isblkkder(2,isb))then                      8d24s16
         if(iprop.eq.1)then
          write(6,*)('in block 2 for is = '),is
          call dws_sync
          call dws_finalize
          stop
         end if
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  8d24s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         8d24s16
         ii=ionex(is)                                                   8d24s16
         do i2=i2s,i2e                                                  8d24s16
          if(i2.eq.i2e)i1n=i1e                                          8d24s16
          i2m=i2-1                                                      8d24s16
          do i1=i10,i1n                                                 8d24s16
           i1m=i1-1                                                     8d24s16
           do i3=0,nocc(isblk1(2,is))-1                                 8d24s16
            do i4=0,nocc(isblk1(1,is))-1                                8d24s16
             do im=0,nvirtc(isblkkder(4,isb))-1                         8d24s16
              imp=im+nocc(isblkkder(4,isb))                             8d24s16
              iad1=itrans(isw)+i1m+nbasdwsc(isw)*imp                    8d24s16
              iad2=jmatda+i3+nocc(isblk1(2,is))*(i4+                    8d24s16
     $             nocc(isblk1(1,is))*(i2m+nvirtc(isblk1(4,is))*im))    8d25s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d24s16
             end do                                                     8d24s16
             ii=ii+1                                                    8d24s16
            end do                                                      8d24s16
           end do                                                       8d24s16
          end do                                                        8d24s16
          i10=1                                                         8d24s16
         end do                                                         8d24s16
         go to 306                                                      8d24s16
        end if
       end do                                                           8d19s16
       write(6,*)('could not fit onex for (oo|vv'')'),
     $      (isblkkder(j,isb),j=1,3),isw
       call dws_sync
       call dws_finalize
       stop
  306  continue                                                         8d22s16
       if(nocc(isblkkder(1,isb))*nocc(isblkkder(2,isb)).eq.0)go to 307  8d24s16
       do is=1,nsdlk
        if(isblk(1,is).eq.isblkkder(1,isb).and.                         8d19s16
     $       isblk(2,is).eq.isblkkder(2,isb).and.isblk(4,is).eq.isw.and.8d19s16
     $       isblk(3,is).eq.isblkkder(3,isb))then                       8d19s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  8d4s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d19s16
         i10=i1s                                                        8d19s16
         i1n=nvirtc(isblk(3,is))                                        8d19s16
         ii=jmats(is)                                                   8d19s16
         do i2=i2s,i2e                                                  8d19s16
          if(i2.eq.i2e)i1n=i1e                                          8d19s16
          i2p=i2-1+nocc(isblk(4,is))                                    8d19s16
          do i1=i10,i1n                                                 8d19s16
           i1m=i1-1
           if(isblk(1,is).eq.isblk(2,is))then                           8d19s16
            do i3=0,nocc(isblk(1,is))-1                                 8d19s16
             do i4=0,i3-1                                               8d19s16
              do im=0,nvirtc(isblkkder(4,isb))-1                        8d19s16
               imp=im+nocc(isblkkder(4,isb))                            8d19s16
               iad1=itrans(isw)+i2p+nbasdwsc(isw)*imp                   8d19s16
               iad2=jmatda+i3+nocc(isblk(1,is))*(i4                     8d19s16
     $             +nocc(isblk(1,is))*(i1m+nvirtc(isblk(3,is))*im))     8d19s16
               bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                        8d19s16
               iad2=jmatda+i4+nocc(isblk(1,is))*(i3+                    8d30s16
     $              nocc(isblk(1,is))*(i1m+nvirtc(isblk(3,is))*im))     8d19s16
               bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                        8d19s16
              end do                                                    8d19s16
              ii=ii+1                                                   8d19s16
             end do                                                     8d19s16
             do im=0,nvirtc(isblkkder(4,isb))-1                         8d19s16
              imp=im+nocc(isblkkder(4,isb))                             8d19s16
              iad1=itrans(isw)+i2p+nbasdwsc(isw)*imp                          8d19s16
              iad2=jmatda+i3+nocc(isblk(1,is))*(i3+                     8d19s16
     $             nocc(isblk(1,is))*(i1m+nvirtc(isblk(3,is))*im))      8d25s16
              bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                         8d19s16
             end do                                                     8d19s16
             ii=ii+1                                                    8d19s16
            end do                                                      8d19s16
           else                                                         8d19s16
            do i3=0,nocc(isblk(2,is))-1                                 8d19s16
             do i4=0,nocc(isblk(1,is))-1                                8d19s16
              do im=0,nvirtc(isblkkder(4,isb))-1                          8d19s16
               imp=im+nocc(isblkkder(4,isb))
               iad1=itrans(isw)+i2p+nbasdwsc(isw)*imp                         8d19s16
               iad2=jmatda+i4+nocc(isblk(1,is))*(i3+                    8d19s16
     $              nocc(isblk(2,is))*(i1m+nvirtc(isblk(3,is))*im))     8d19s16
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                        8d19s16
              end do                                                    8d19s16
              ii=ii+1                                                   8d19s16
             end do                                                     8d19s16
            end do                                                      8d19s16
           end if                                                       8d19s16
          end do                                                        8d19s16
          i10=1                                                         8d19s16
         end do                                                         8d19s16
         go to 307                                                      8d22s16
         else if(isblk(1,is).eq.isblkkder(1,isb).and.                         8d19s16
     $       isblk(2,is).eq.isblkkder(2,isb).and.isblk(3,is).eq.isw.and.8d22s16
     $       isblk(4,is).eq.isblkkder(3,isb))then                       8d22s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  8d4s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d19s16
         i10=i1s                                                        8d19s16
         i1n=nvirtc(isblk(3,is))                                        8d19s16
         ii=jmats(is)                                                   8d19s16
         do i2=i2s,i2e                                                  8d19s16
          if(i2.eq.i2e)i1n=i1e                                          8d19s16
          i2m=i2-1                                                      8d22s16
          do i1=i10,i1n                                                 8d19s16
           i1p=i1-1+nocc(isblk(3,is))                                   8d22s16
           do i3=0,nocc(isblk(2,is))-1                                  8d19s16
            do i4=0,nocc(isblk(1,is))-1                                 8d19s16
             do im=0,nvirtc(isblkkder(4,isb))-1                         8d19s16
              imp=im+nocc(isblkkder(4,isb))                             8d22s16
              iad1=itrans(isw)+i1p+nbasdwsc(isw)*imp                    8d22s16
              iad2=jmatda+i4+nocc(isblk(1,is))*(i3+                     8d22s16
     $             nocc(isblk(2,is))*(i2m+nvirtc(isblk(4,is))*im))      8d19s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d19s16
             end do                                                     8d19s16
             ii=ii+1                                                    8d19s16
            end do                                                      8d19s16
           end do                                                       8d19s16
          end do                                                        8d19s16
          i10=1                                                         8d19s16
         end do                                                         8d19s16
         go to 307                                                      8d22s16
         else if(isblk(2,is).eq.isblkkder(1,isb).and.                   8d24s16
     $       isblk(1,is).eq.isblkkder(2,isb).and.isblk(3,is).eq.isw.and.8d24s16
     $       isblk(4,is).eq.isblkkder(3,isb))then                       8d24s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  8d24s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        8d24s16
         ii=jmats(is)                                                   8d24s16
         do i2=i2s,i2e                                                  8d24s16
          if(i2.eq.i2e)i1n=i1e                                          8d24s16
          i2m=i2-1                                                      8d24s16
          do i1=i10,i1n                                                 8d24s16
           i1p=i1-1+nocc(isblk(3,is))                                   8d24s16
           do i3=0,nocc(isblk(2,is))-1                                  8d24s16
            do i4=0,nocc(isblk(1,is))-1                                 8d24s16
             do im=0,nvirtc(isblkkder(4,isb))-1                         8d24s16
              imp=im+nocc(isblkkder(4,isb))                             8d24s16
              iad1=itrans(isw)+i1p+nbasdwsc(isw)*imp                    8d24s16
              iad2=jmatda+i3+nocc(isblk(2,is))*(i4+                     8d24s16
     $             nocc(isblk(1,is))*(i2m+nvirtc(isblk(4,is))*im))      8d24s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d24s16
             end do                                                     8d24s16
             ii=ii+1                                                    8d24s16
            end do                                                      8d24s16
           end do                                                       8d24s16
          end do                                                        8d24s16
          i10=1                                                         8d24s16
         end do                                                         8d24s16
         go to 307                                                      8d24s16
         else if(isblk(2,is).eq.isblkkder(1,isb).and.                   8d24s16
     $       isblk(1,is).eq.isblkkder(2,isb).and.isblk(4,is).eq.isw.and.8d24s16
     $       isblk(3,is).eq.isblkkder(3,isb))then                       8d24s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  8d24s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        8d24s16
         ii=jmats(is)                                                   8d24s16
         do i2=i2s,i2e                                                  8d24s16
          if(i2.eq.i2e)i1n=i1e                                          8d24s16
          i2p=i2-1+nocc(isblk(4,is))                                    8d24s16
          do i1=i10,i1n                                                 8d24s16
           i1m=i1-1                                                     8d24s16
           do i3=0,nocc(isblk(2,is))-1                                  8d24s16
            do i4=0,nocc(isblk(1,is))-1                                 8d24s16
             do im=0,nvirtc(isblkkder(4,isb))-1                         8d24s16
              imp=im+nocc(isblkkder(4,isb))                             8d24s16
              iad1=itrans(isw)+i2p+nbasdwsc(isw)*imp                    8d24s16
              iad2=jmatda+i3+nocc(isblk(2,is))*(i4+                     8d24s16
     $             nocc(isblk(1,is))*(i1m+nvirtc(isblk(3,is))*im))      8d24s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d24s16
             end do                                                     8d24s16
             ii=ii+1                                                    8d24s16
            end do                                                      8d24s16
           end do                                                       8d24s16
          end do                                                        8d24s16
          i10=1                                                         8d24s16
         end do                                                         8d24s16
         go to 307                                                      8d24s16
        end if
       end do
       write(6,*)('need jhelp for (oo|vv'')'),(isblkkder(j,isb),j=1,3),
     $      isw
       write(6,*)('got: ')
       do i=1,nsdlk
        write(6,12)(isblk(j,i),j=1,4)
       end do
       call dws_sync
       call dws_finalize
       stop
  307  continue                                                         8d22s16
       kmatda=ibcoff
       ibcoff=kmatda+nall
       call enough('parajkfromhd0.  7',bc,ibc)
       do i=0,nall-1                                                    8d4s16
        bc(kmatda+i)=0d0                                                8d4s16
       end do                                                           8d4s16
c
c     xor contribution to k'_nm^ab=(n'b|ma): (ov|ov) and (vv|ov)
c
       isw=multh(isblkkder(1,isb),iprop)
       if(nocc(isw).eq.0)go to 400                                      8d24s16
       do is=1,nsdlkk
        if(isblkk(1,is).eq.isw.and.isblkk(2,is).eq.isblkkder(2,isb).and.
     $     isblkk(3,is).eq.isblkkder(3,isb).and.
     $       isblkk(4,is).eq.isblkkder(4,isb))then
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,8d22s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d22s16
         i10=i1s                                                        8d22s16
         i1n=nvirtc(isblkk(3,is))                                       8d22s16
         ii=kmats(is)                                                   8d22s16
         do i2=i2s,i2e                                                  8d22s16
          if(i2.eq.i2e)i1n=i1e                                          8d22s16
          i2m=i2-1                                                      8d22s16
          do i1=i10,i1n                                                 8d22s16
           i1m=i1-1                                                     8d22s16
           do i3=0,nocc(isblkk(2,is))-1                                 8d22s16
            do i4=0,nocc(isblkk(1,is))-1                                8d22s16
             do im=0,nocc(isblkkder(1,isb))-1                           8d22s16
              iad1=itrans(isw)+i4+nbasdwsc(isw)*im                      8d22s16
              iad2=kmatda+im+nocc(isblkkder(1,isb))*(i3                 8d22s16
     $             +nocc(isblkk(2,is))*(i1m+nvirtc(isblkk(3,is))*i2m))  8d22s16
              bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                         8d22s16
             end do                                                     8d22s16
             ii=ii+1                                                    8d22s16
            end do                                                      8d22s16
           end do                                                       8d22s16
          end do                                                        8d22s16
          i10=1                                                         8d22s16
         end do                                                         8d22s16
         go to 400                                                      8d22s16
        else if(isblkk(2,is).eq.isw.and.isblkk(1,is).eq.isblkkder(2,isb)8d22s16
     $     .and.isblkk(4,is).eq.isblkkder(3,isb).and.                   8d22s16
     $       isblkk(3,is).eq.isblkkder(4,isb))then                      8d22s16
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,8d22s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d22s16
         i10=i1s                                                        8d22s16
         i1n=nvirtc(isblkk(3,is))                                       8d22s16
         ii=kmats(is)                                                   8d22s16
         do i2=i2s,i2e                                                  8d22s16
          if(i2.eq.i2e)i1n=i1e                                          8d22s16
          i2m=i2-1                                                      8d22s16
          do i1=i10,i1n                                                 8d22s16
           i1m=i1-1                                                     8d22s16
           do i3=0,nocc(isblkk(2,is))-1                                 8d22s16
            do i4=0,nocc(isblkk(1,is))-1                                8d22s16
             do im=0,nocc(isblkkder(1,isb))-1                           8d22s16
              iad1=itrans(isw)+i3+nbasdwsc(isw)*im                      8d22s16
              iad2=kmatda+im+nocc(isblkkder(1,isb))*(i4                 8d22s16
     $             +nocc(isblkk(1,is))*(i2m+nvirtc(isblkk(4,is))*i1m))  8d22s16
              bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                         8d22s16
             end do                                                     8d22s16
             ii=ii+1                                                    8d22s16
            end do                                                      8d22s16
           end do                                                       8d22s16
          end do                                                        8d22s16
          i10=1                                                         8d22s16
         end do                                                         8d22s16
         go to 400                                                      8d22s16
        end if
       end do
       write(6,*)('could not find k for (n''b|ma)')
       write(6,*)('want: '),isw,(isblkkder(j,isb),j=2,4)
       write(6,*)('got: ')
       do i=1,nsdlkk
        write(6,12)(isblkk(j,i),j=1,4)
       end do
       call dws_sync
       call dws_finalize
       stop
  400  continue                                                         8d22s16
c     xor contribution to k'_nm^ab=(n'b|ma): (ov|ov) and (vv|ov)
       if(nocc(isblkkder(2,isb)).eq.0)go to 401                         8d24s16
       do is=1,nsdlk1
        if(isblk1(3,is).eq.isblkkder(2,isb).and.
     $     isblk1(4,is).eq.isblkkder(3,isb).and.
     $     isblk1(1,is).eq.isw.and.
     $     isblk1(2,is).eq.isblkkder(4,isb))then
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  8d22s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d22s16
         nhere=ih+1-il
         nrw=nvirtc(isblk1(1,is))*nvirtc(isblk1(2,is))
         i10=i1s                                                        8d22s16
         i1n=nocc(isblk1(3,is))                                         8d22s16
         ii=i3x(is)                                                     8d22s16
         do i2=i2s,i2e                                                  8d22s16
          if(i2.eq.i2e)i1n=i1e                                          8d22s16
          i2m=i2-1                                                      8d22s16
          do i1=i10,i1n                                                 8d22s16
           i1m=i1-1                                                     8d22s16
           if(isblk1(1,is).eq.isblk1(2,is))then                         8d22s16
            do i3=0,nvirtc(isblk1(2,is))-1                              8d22s16
             i3p=i3+nocc(isblk1(2,is))
             do i4=0,i3-1                                               8d22s16
              i4p=i4+nocc(isblk1(2,is))
              do im=0,nocc(isblkkder(1,isb))-1                          8d22s16
               iad1=itrans(isw)+i3p+nbasdwsc(isw)*im                     8d22s16
               iad2=kmatda+im+nocc(isblkkder(1,isb))*(i1m               8d22s16
     $              +nocc(isblk1(3,is))*(i2m+nvirtc(isblk1(4,is))*i4))  8d29s16
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                        8d22s16
               iad1=itrans(isw)+i4p+nbasdwsc(isw)*im                     8d22s16
               iad2=kmatda+im+nocc(isblkkder(1,isb))*(i1m               8d22s16
     $              +nocc(isblk1(3,is))*(i2m+nvirtc(isblk1(4,is))*i3))  8d29s16
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                        8d22s16
              end do                                                    8d22s16
              ii=ii+1                                                   8d22s16
             end do                                                     8d22s16
             do im=0,nocc(isblkkder(1,isb))-1                           8d22s16
              iad1=itrans(isw)+i3p+nbasdwsc(isw)*im                     8d22s16
              iad2=kmatda+im+nocc(isblkkder(1,isb))*(i1m                8d22s16
     $             +nocc(isblk1(3,is))*(i2m+nvirtc(isblk1(4,is))*i3))   8d29s16
   87         format(5i5,2es15.7,i8,es15.7)
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d22s16
             end do                                                     8d22s16
             ii=ii+1                                                    8d22s16
            end do                                                      8d22s16
           else                                                         8d22s16
            do i3=0,nvirtc(isblk1(2,is))-1                              8d22s16
             do i4=0,nvirtc(isblk1(1,is))-1                             8d22s16
              i4p=i4+nocc(isblk1(1,is))
              do im=0,nocc(isblkkder(1,isb))-1                           8d22s16
               iad1=itrans(isw)+i4p+nbasdwsc(isw)*im                    8d22s16
               iad2=kmatda+im+nocc(isblkkder(1,isb))*(i1m               8d22s16
     $             +nocc(isblk1(3,is))*(i2m+nvirtc(isblk1(4,is))*i3))    8d22s16
               bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                         8d22s16
              end do                                                     8d22s16
              ii=ii+1                                                    8d22s16
             end do                                                      8d22s16
            end do                                                       8d22s16
           end if                                                       8d22s16
          end do                                                        8d22s16
          i10=1                                                         8d22s16
         end do                                                         8d22s16
         go to 401                                                      8d22s16
        else if(isblk1(3,is).eq.isblkkder(2,isb).and.
     $     isblk1(4,is).eq.isblkkder(3,isb).and.
     $     isblk1(2,is).eq.isw.and.                                     8d22s16
     $     isblk1(1,is).eq.isblkkder(4,isb))then                        8d22s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  8d22s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d22s16
         i10=i1s                                                        8d22s16
         i1n=nocc(isblk1(3,is))                                         8d22s16
         ii=i3x(is)                                                     8d22s16
         do i2=i2s,i2e                                                  8d22s16
          if(i2.eq.i2e)i1n=i1e                                          8d22s16
          i2m=i2-1                                                      8d22s16
          do i1=i10,i1n                                                 8d22s16
           i1m=i1-1                                                     8d22s16
           do i3=0,nvirtc(isblk1(2,is))-1                               8d22s16
            i3p=i3+nocc(isblk1(2,is))                                   8d22s16
            do i4=0,nvirtc(isblk1(1,is))-1                              8d22s16
             do im=0,nocc(isblkkder(1,isb))-1                           8d22s16
              iad1=itrans(isw)+i3p+nbasdwsc(isw)*im                     8d22s16
              iad2=kmatda+im+nocc(isblkkder(1,isb))*(i1m                8d22s16
     $             +nocc(isblk1(3,is))*(i2m+nvirtc(isblk1(4,is))*i4))   8d22s16
              bc(iad2)=bc(iad2)+bc(ii)*bc(iad1)                         8d22s16
             end do                                                     8d22s16
             ii=ii+1                                                    8d22s16
            end do                                                      8d22s16
           end do                                                       8d22s16
          end do                                                        8d22s16
          i10=1                                                         8d22s16
         end do                                                         8d22s16
         go to 401                                                      8d22s16
        end if
       end do
       if(min(nbasdwsc(isw),nbasdwsc(isblkkder(4,isb)),                 11d28s22
     $      nbasdwsc(isblkkder(3,isb))).gt.0)then                       11d28s22
        write(6,*)('could not find 3x for (n''b|ma)')
        write(6,*)('want: '),isw,isblkkder(4,isb),isblkkder(2,isb),
     $      isblkkder(3,isb)
        write(6,*)('got: ')
        do i=1,nsdlk1
         write(6,12)(isblk1(j,i),j=1,4)
        end do
        call dws_sync
        call dws_finalize
        stop
       end if                                                           11d28s22
  401  continue                                                         8d22s16
c
c     xor contribution to k'=(nb'|ma): (oo|ov) and (ov|ov)
c
       isw=multh(iprop,isblkkder(4,isb))                                8d22s16
       if(nocc(isw)*nocc(isblkkder(1,isb))*nocc(isblkkder(2,isb)).eq.0) 8d24s16
     $      go to 402                                                   8d24s16
       do is=1,nsdlk1
        if(isblk1(1,is).eq.isblkkder(1,isb).and.
     $     isblk1(2,is).eq.isw.and.
     $     isblk1(3,is).eq.isblkkder(2,isb).and.
     $       isblk1(4,is).eq.isblkkder(3,isb))then
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  8d22s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d22s16
         i10=i1s                                                        8d22s16
         i1n=nocc(isblk1(3,is))
         ii=ionex(is)
         do i2=i2s,i2e
          i2m=i2-1
          if(i2.eq.i2e)i1n=i1e
          do i1=i10,i1n
           i1m=i1-1
           if(isblk1(1,is).eq.isblk1(2,is))then
            do i3=0,nocc(isblk1(1,is))-1
             do i4=0,i3-1
              do im=0,nvirtc(isblkkder(4,isb))-1
               imp=im+nocc(isblkkder(4,isb))
               iad1=itrans(isw)+i3+nbasdwsc(isw)*imp
               iad2=kmatda+i4+nocc(isblk1(1,is))*(i1m+nocc(isblk1(3,is))
     $              *(i2m+nvirtc(isblk1(4,is))*im))
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)
               iad1=itrans(isw)+i4+nbasdwsc(isw)*imp
               iad2=kmatda+i3+nocc(isblk1(1,is))*(i1m+nocc(isblk1(3,is))
     $              *(i2m+nvirtc(isblk1(4,is))*im))
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)
              end do
              ii=ii+1
             end do
             do im=0,nvirtc(isblkkder(4,isb))-1
              imp=im+nocc(isblkkder(4,isb))
              iad1=itrans(isw)+i3+nbasdwsc(isw)*imp
              iad2=kmatda+i3+nocc(isblk1(1,is))*(i1m+nocc(isblk1(3,is))
     $             *(i2m+nvirtc(isblk1(4,is))*im))
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)
             end do
             ii=ii+1
            end do
           else
            do i3=0,nocc(isblk1(2,is))-1
             do i4=0,nocc(isblk1(1,is))-1
              do im=0,nvirtc(isblkkder(4,isb))-1
               imp=im+nocc(isblkkder(4,isb))
               iad1=itrans(isw)+i3+nbasdwsc(isw)*imp
               iad2=kmatda+i4+nocc(isblk1(1,is))*(i1m+nocc(isblk1(3,is))8d22s16
     $             *(i2m+nvirtc(isblk1(4,is))*im))                      8d22s16
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 402
        else if(isblk1(2,is).eq.isblkkder(1,isb).and.                   8d22s16
     $     isblk1(1,is).eq.isw.and.                                     8d22s16
     $     isblk1(3,is).eq.isblkkder(2,isb).and.
     $       isblk1(4,is).eq.isblkkder(3,isb))then
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  8d22s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d22s16
         i10=i1s                                                        8d22s16
         i1n=nocc(isblk1(3,is))
         ii=ionex(is)
         do i2=i2s,i2e
          i2m=i2-1
          if(i2.eq.i2e)i1n=i1e
          do i1=i10,i1n
           i1m=i1-1
           do i3=0,nocc(isblk1(2,is))-1
            do i4=0,nocc(isblk1(1,is))-1
             do im=0,nvirtc(isblkkder(4,isb))-1
              imp=im+nocc(isblkkder(4,isb))
              iad1=itrans(isw)+i4+nbasdwsc(isw)*imp                     8d22s16
              iad2=kmatda+i3+nocc(isblk1(2,is))*(i1m+nocc(isblk1(3,is)) 8d22s16
     $             *(i2m+nvirtc(isblk1(4,is))*im))                      8d22s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 402
        end if
       end do
       write(6,*)('could not find onex for (nb''|ma)')
       write(6,*)('want '),isblkkder(1,isb),isw,isblkkder(2,isb),
     $      isblkkder(3,isb)
       do i=1,nsdlk1
        write(6,12)(isblk1(j,i),j=1,4)
       end do
       call dws_sync
       call dws_finalize
       stop
  402  continue
c
c     xor contribution to k'=(nb'|ma): (oo|ov) and (ov|ov)
c
       if(nocc(isblkkder(1,isb))*nocc(isblkkder(2,isb)).eq.0)go to 403  8d24s16
       do is=1,nsdlkk                                                   8d22s16
        if(isblkk(1,is).eq.isblkkder(1,isb).and.                        8d22s16
     $     isblkk(2,is).eq.isblkkder(2,isb).and.                        8d22s16
     $     isblkk(3,is).eq.isblkkder(3,isb).and.                        8d22s16
     $       isblkk(4,is).eq.isw)then                                   8d22s16
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,8d22s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d22s16
         i10=i1s                                                        8d22s16
         i1n=nvirtc(isblkk(3,is))                                       8d22s16
         ii=kmats(is)                                                   8d22s16
         do i2=i2s,i2e                                                  8d22s16
          i2p=i2-1+nocc(isblkk(4,is))                                   8d22s16
          if(i2.eq.i2e)i1n=i1e                                          8d22s16
          do i1=i10,i1n                                                 8d22s16
           i1m=i1-1                                                     8d22s16
           do i3=0,nocc(isblkk(2,is))-1                                 8d22s16
            do i4=0,nocc(isblkk(1,is))-1                                8d22s16
             do im=0,nvirtc(isblkkder(4,isb))-1                         8d22s16
              imp=im+nocc(isblkkder(4,isb))
              iad1=itrans(isw)+i2p+nbasdwsc(isw)*imp                     8d22s16
              iad2=kmatda+i4+nocc(isblkk(1,is))*(i3+nocc(isblkk(2,is))  8d22s16
     $             *(i1m+nvirtc(isblkk(3,is))*im))                      8d22s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d22s16
             end do                                                     8d22s16
             ii=ii+1                                                    8d22s16
            end do                                                      8d22s16
           end do                                                       8d22s16
          end do                                                        8d22s16
          i10=1                                                         8d22s16
         end do                                                         8d22s16
         go to 403                                                      8d22s16
        else if(isblkk(2,is).eq.isblkkder(1,isb).and.                        8d22s16
     $     isblkk(1,is).eq.isblkkder(2,isb).and.                        8d22s16
     $     isblkk(4,is).eq.isblkkder(3,isb).and.                        8d22s16
     $       isblkk(3,is).eq.isw)then                                   8d22s16
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,8d22s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d22s16
         i10=i1s                                                        8d22s16
         i1n=nvirtc(isblkk(3,is))                                       8d22s16
         ii=kmats(is)                                                   8d22s16
         do i2=i2s,i2e                                                  8d22s16
          i2m=i2-1
          if(i2.eq.i2e)i1n=i1e                                          8d22s16
          do i1=i10,i1n                                                 8d22s16
           i1p=i1-1+nocc(isblkk(3,is))
           do i3=0,nocc(isblkk(2,is))-1                                 8d22s16
            do i4=0,nocc(isblkk(1,is))-1                                8d22s16
             do im=0,nvirtc(isblkkder(4,isb))-1                         8d22s16
              imp=im+nocc(isblkkder(4,isb))
              iad1=itrans(isw)+i1p+nbasdwsc(isw)*imp                    8d22s16
              iad2=kmatda+i3+nocc(isblkk(2,is))*(i4+nocc(isblkk(1,is))  8d22s16
     $             *(i2m+nvirtc(isblkk(4,is))*im))                      8d22s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d22s16
             end do                                                     8d22s16
             ii=ii+1                                                    8d22s16
            end do                                                      8d22s16
           end do                                                       8d22s16
          end do                                                        8d22s16
          i10=1                                                         8d22s16
         end do                                                         8d22s16
         go to 403                                                      8d22s16
        end if                                                          8d22s16
       end do                                                           8d22s16
       write(6,*)('could not find k for (nb''|ma)')                     8d22s16
       write(6,*)('want: '),(isblkkder(j,isb),j=1,3),isw                8d22s16
       write(6,*)('got: ')                                              8d22s16
       do i=1,nsdlkk                                                    8d22s16
        write(6,12)(isblkk(j,i),j=1,4)                                  8d22s16
       end do                                                           8d22s16
       call dws_sync                                                    8d22s16
       call dws_finalize                                                8d22s16
       stop                                                             8d22s16
  403  continue                                                         8d22s16
c
c     xor contribution to k'=(nb|m'a): (ov|ov) and (ov|vv)
c
       isw=multh(isblkkder(2,isb),iprop)                                8d22s16
       if(nocc(isw).eq.0)go to 404                                      8d24s16
       do is=1,nsdlkk
        if(isblkk(1,is).eq.isblkkder(1,isb).and.
     $     isblkk(2,is).eq.isw.and.
     $     isblkk(3,is).eq.isblkkder(3,isb).and.                        8d22s16
     $       isblkk(4,is).eq.isblkkder(4,isb))then                      8d22s16
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,8d22s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d22s16
         i10=i1s                                                        8d22s16
         i1n=nvirtc(isblkk(3,is))                                       8d22s16
         ii=kmats(is)                                                   8d22s16
         do i2=i2s,i2e                                                  8d22s16
          if(i2.eq.i2e)i1n=i1e                                          8d22s16
          i2m=i2-1
          do i1=i10,i1n                                                 8d22s16
           i1m=i1-1
           do i3=0,nocc(isblkk(2,is))-1                                 8d29s16
            do i4=0,nocc(isblkk(1,is))-1                                8d29s16
             do im=0,nocc(isblkkder(2,isb))-1                           8d22s16
              iad1=itrans(isw)+i3+nbasdwsc(isw)*im
              iad2=kmatda+i4+nocc(isblkk(1,is))*(im
     $           +nocc(isblkkder(2,isb))*(i1m+nvirtc(isblkk(3,is))*i2m))
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)
             end do                                                     8d22s16
             ii=ii+1
            end do                                                      8d22s16
           end do                                                       8d22s16
          end do                                                        8d22s16
          i10=1                                                         8d22s16
         end do                                                         8d22s16
         go to 404
        else if(isblkk(2,is).eq.isblkkder(1,isb).and.                   8d22s16
     $     isblkk(1,is).eq.isw.and.                                     8d22s16
     $     isblkk(4,is).eq.isblkkder(3,isb).and.                        8d22s16
     $       isblkk(3,is).eq.isblkkder(4,isb))then                      8d22s16
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,8d22s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d22s16
         i10=i1s                                                        8d22s16
         i1n=nvirtc(isblkk(3,is))                                       8d22s16
         ii=kmats(is)                                                   8d22s16
         do i2=i2s,i2e                                                  8d22s16
          if(i2.eq.i2e)i1n=i1e                                          8d22s16
          i2m=i2-1
          do i1=i10,i1n                                                 8d22s16
           i1m=i1-1
           do i3=0,nocc(isblkk(2,is))-1                                 8d29s16
            do i4=0,nocc(isblkk(1,is))-1                                8d29s16
             do im=0,nocc(isblkkder(2,isb))-1                           8d22s16
              iad1=itrans(isw)+i4+nbasdwsc(isw)*im                      8d22s16
              iad2=kmatda+i3+nocc(isblkk(2,is))*(im                     8d22s16
     $           +nocc(isblkkder(2,isb))*(i2m+nvirtc(isblkk(4,is))*i1m))8d22s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)
             end do                                                     8d22s16
             ii=ii+1
            end do                                                      8d22s16
           end do                                                       8d22s16
          end do                                                        8d22s16
          i10=1                                                         8d22s16
         end do                                                         8d22s16
         go to 404
        end if
       end do
       write(6,*)('could not find k for (nb|m''a)')
       write(6,*)('want: '),isblkkder(1,isb),isw,(isblkkder(j,isb),j=3, 8d22s16
     $      4)                                                          8d22s16
       write(6,*)('got: ')                                              8d22s16
       do i=1,nsdlkk                                                    8d22s16
        write(6,12)(isblkk(j,i),j=1,4)                                  8d22s16
       end do                                                           8d22s16
       call dws_sync
       call dws_finalize
       stop
  404  continue
c     xor contribution to k'=(nb|m'a): (ov|ov) and (ov|vv)
       if(nocc(isblkkder(1,isb)).eq.0)go to 405                         8d24s16
       do is=1,nsdlk1
        if(isblk1(1,is).eq.isw.and.isblk1(2,is).eq.isblkkder(3,isb)
     $       .and.isblk1(3,is).eq.isblkkder(1,isb).and.
     $       isblk1(4,is).eq.isblkkder(4,isb))then
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  8d22s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d22s16
         i10=i1s                                                        8d22s16
         i1n=nocc(isblk1(3,is))                                         8d22s16
         ii=i3x(is)                                                     8d22s16
         do i2=i2s,i2e                                                  8d22s16
          if(i2.eq.i2e)i1n=i1e                                          8d22s16
          i2m=i2-1
          do i1=i10,i1n                                                 8d22s16
           i1m=i1-1
           if(isblk1(1,is).eq.isblk1(2,is))then
            do i4=0,nvirtc(isblk1(1,is))-1                              8d22s16
             i4p=i4+nocc(isblk1(1,is))
             do i3=0,i4-1
              i3p=i3+nocc(isblk1(1,is))
              do im=0,nocc(isblkkder(2,isb))-1
               iad1=itrans(isw)+i3p+nbasdwsc(isw)*im
               iad2=kmatda+i1m+nocc(isblk1(3,is))*(im
     $            +nocc(isblkkder(2,isb))*(i4+nvirtc(isblk1(1,is))*i2m))5d4s22
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)
               iad1=itrans(isw)+i4p+nbasdwsc(isw)*im
               iad2=kmatda+i1m+nocc(isblk1(3,is))*(im
     $            +nocc(isblkkder(2,isb))*(i3+nvirtc(isblk1(1,is))*i2m))5d4s22
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)
              end do
              ii=ii+1
             end do
             do im=0,nocc(isblkkder(2,isb))-1
              iad1=itrans(isw)+i4p+nbasdwsc(isw)*im
              iad2=kmatda+i1m+nocc(isblk1(3,is))*(im
     $            +nocc(isblkkder(2,isb))*(i4+nvirtc(isblk1(1,is))*i2m))5d4s22
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)
             end do
             ii=ii+1
            end do
           else
            do i4=0,nvirtc(isblk1(2,is))-1
             do i3=0,nvirtc(isblk1(1,is))-1
              i3p=i3+nocc(isblk1(1,is))
              do im=0,nocc(isblkkder(2,isb))-1
               iad1=itrans(isw)+i3p+nbasdwsc(isw)*im
               iad2=kmatda+i1m+nocc(isblk1(3,is))*(im
     $            +nocc(isblkkder(2,isb))*(i4+nvirtc(isblk1(2,is))*i2m))
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 405
        else if(isblk1(2,is).eq.isw.and.isblk1(1,is).eq.isblkkder(3,isb)8d22s16
     $       .and.isblk1(3,is).eq.isblkkder(1,isb).and.
     $       isblk1(4,is).eq.isblkkder(4,isb))then
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  8d22s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d22s16
         i10=i1s                                                        8d22s16
         i1n=nocc(isblk1(3,is))                                         8d22s16
         ii=i3x(is)                                                     8d22s16
         do i2=i2s,i2e                                                  8d22s16
          if(i2.eq.i2e)i1n=i1e                                          8d22s16
          i2m=i2-1
          do i1=i10,i1n                                                 8d22s16
           i1m=i1-1
           do i4=0,nvirtc(isblk1(2,is))-1
            i4p=i4+nocc(isblk1(2,is))                                   8d22s16
            do i3=0,nvirtc(isblk1(1,is))-1                              8d22s16
             do im=0,nocc(isblkkder(2,isb))-1                           8d22s16
              iad1=itrans(isw)+i4p+nbasdwsc(isw)*im                     8d22s16
              iad2=kmatda+i1m+nocc(isblk1(3,is))*(im                    8d22s16
     $            +nocc(isblkkder(2,isb))*(i3+nvirtc(isblk1(1,is))*i2m))8d22s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d22s16
             end do                                                     8d22s16
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 405
        end if
       end do
       if(min(nbasdwsc(isw),nbasdwsc(isblkkder(3,isb)),                 11d28s22
     $      nbasdwsc(isblkkder(4,isb))).gt.0)then                       11d28s22
        write(6,*)('could not find 3x for (nb|m''a)')
        write(6,*)('want: '),isw,isblkkder(3,isb),isblkkder(1,isb),
     $      isblkkder(4,isb)
        write(6,*)('got ')
        do i=1,nsdlk1
         write(6,12)(isblk1(j,i),j=1,4)
        end do
        call dws_sync
        call dws_finalize
        stop
       end if                                                           11d28s22
  405  continue
c
c     xor contribution to k'=(nb|ma'): (ov|oo) and (ov|ov)
c
       isw=multh(isblkkder(3,isb),iprop)
       if(nocc(isw).eq.0)go to 406                                      8d24s16
       if(nocc(isblkkder(1,isb)).eq.0)go to 406                         8d24s16
       do is=1,nsdlk1
        if(isblk1(1,is).eq.isblkkder(2,isb).and.
     $     isblk1(2,is).eq.isw.and.
     $     isblk1(3,is).eq.isblkkder(1,isb).and.
     $       isblk1(4,is).eq.isblkkder(4,isb))then                      8d22s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  8d22s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d22s16
         i10=i1s                                                        8d22s16
         i1n=nocc(isblk1(3,is))                                         8d22s16
         ii=ionex(is)                                                   8d22s16
         do i2=i2s,i2e                                                  8d22s16
          if(i2.eq.i2e)i1n=i1e                                          8d22s16
          i2m=i2-1
          do i1=i10,i1n                                                 8d22s16
           i1m=i1-1
           if(isblk1(1,is).eq.isblk1(2,is))then
            do i3=0,nocc(isblk1(1,is))-1
             do i4=0,i3-1
              do im=0,nvirtc(isblkkder(3,isb))-1
               imp=im+nocc(isblkkder(3,isb))
               iad1=itrans(isw)+i4+nbasdwsc(isw)*imp
               iad2=kmatda+i1m+nocc(isblk1(3,is))*(i3+nocc(isblk1(1,is))
     $              *(im+nvirtc(isblkkder(3,isb))*i2m))
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)
               iad1=itrans(isw)+i3+nbasdwsc(isw)*imp
               iad2=kmatda+i1m+nocc(isblk1(3,is))*(i4+nocc(isblk1(1,is))
     $              *(im+nvirtc(isblkkder(3,isb))*i2m))
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)
              end do
              ii=ii+1
             end do
             do im=0,nvirtc(isblkkder(3,isb))-1
              imp=im+nocc(isblkkder(3,isb))
              iad1=itrans(isw)+i3+nbasdwsc(isw)*imp
              iad2=kmatda+i1m+nocc(isblk1(3,is))*(i3+nocc(isblk1(1,is))
     $             *(im+nvirtc(isblkkder(3,isb))*i2m))
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)
             end do
             ii=ii+1
            end do
           else
            do i3=0,nocc(isblk1(2,is))-1
             do i4=0,nocc(isblk1(1,is))-1
              do im=0,nvirtc(isblkkder(3,isb))-1
               imp=im+nocc(isblkkder(3,isb))
               iad1=itrans(isw)+i3+nbasdwsc(isw)*imp
               iad2=kmatda+i1m+nocc(isblk1(3,is))*(i4+nocc(isblk1(1,is))
     $              *(im+nvirtc(isblkkder(3,isb))*i2m))
               bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 406
        else if(isblk1(2,is).eq.isblkkder(2,isb).and.                   8d22s16
     $     isblk1(1,is).eq.isw.and.                                     8d22s16
     $     isblk1(3,is).eq.isblkkder(1,isb).and.
     $       isblk1(4,is).eq.isblkkder(4,isb))then                      8d22s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  8d22s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d22s16
         i10=i1s                                                        8d22s16
         i1n=nocc(isblk1(3,is))                                         8d22s16
         ii=ionex(is)                                                   8d22s16
         do i2=i2s,i2e                                                  8d22s16
          if(i2.eq.i2e)i1n=i1e                                          8d22s16
          i2m=i2-1
          do i1=i10,i1n                                                 8d22s16
           i1m=i1-1
           do i3=0,nocc(isblk1(2,is))-1
            do i4=0,nocc(isblk1(1,is))-1
             do im=0,nvirtc(isblkkder(3,isb))-1
              imp=im+nocc(isblkkder(3,isb))
              iad1=itrans(isw)+i4+nbasdwsc(isw)*imp                     8d22s16
              iad2=kmatda+i1m+nocc(isblk1(3,is))*(i3+nocc(isblk1(2,is)) 8d22s16
     $             *(im+nvirtc(isblkkder(3,isb))*i2m))
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 406
        end if
       end do
       write(6,*)('could not find onex for (nb|ma'')')
       write(6,*)('want: '),isblkkder(2,isb),isw,isblkkder(1,isb),
     $      isblkkder(4,isb)
       write(6,*)('got: ')
       do i=1,nsdlk1
        write(6,12)(isblk1(j,i),j=1,4)
       end do
       call dws_sync
       call dws_finalize
       stop
  406  continue
c     xor contribution to k'=(nb|ma'): (ov|oo) and (ov|ov)
       if(nocc(isblkkder(1,isb))*nocc(isblkkder(2,isb)).eq.0)go to 407  8d24s16
       do is=1,nsdlkk
        if(isblkk(1,is).eq.isblkkder(1,isb).and.
     $     isblkk(2,is).eq.isblkkder(2,isb).and.
     $     isblkk(3,is).eq.isw.and.
     $       isblkk(4,is).eq.isblkkder(4,isb))then
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,8d22s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d22s16
         i10=i1s                                                        8d22s16
         i1n=nvirtc(isblkk(3,is))                                         8d22s16
         ii=kmats(is)                                                   8d22s16
         do i2=i2s,i2e                                                  8d22s16
          if(i2.eq.i2e)i1n=i1e                                          8d22s16
          i2m=i2-1
          do i1=i10,i1n                                                 8d22s16
           i1p=i1-1+nocc(isblkk(3,is))                                  8d25s16
           do i4=0,nocc(isblkk(2,is))-1
            do i3=0,nocc(isblkk(1,is))-1
             do im=0,nvirtc(isblkkder(3,isb))-1
              imp=im+nocc(isblkkder(3,isb))
              iad1=itrans(isw)+i1p+nbasdwsc(isw)*imp
              iad2=kmatda+i3+nocc(isblkk(1,is))*(i4                     8d22s16
     $            +nocc(isblkk(2,is))*(im+nvirtc(isblkkder(3,isb))*i2m))
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 407
        else if(isblkk(2,is).eq.isblkkder(1,isb).and.                   8d24s16
     $     isblkk(1,is).eq.isblkkder(2,isb).and.                        8d24s16
     $     isblkk(4,is).eq.isw.and.                                     8d24s16
     $       isblkk(3,is).eq.isblkkder(4,isb))then                      8d24s16
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,8d24s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       8d24s16
         ii=kmats(is)                                                   8d24s16
         do i2=i2s,i2e                                                  8d24s16
          if(i2.eq.i2e)i1n=i1e                                          8d24s16
          i2p=i2-1+nocc(isblkk(4,is))                                   8d25s16
          do i1=i10,i1n                                                 8d24s16
           i1m=i1-1                                                     8d25s16
           do i4=0,nocc(isblkk(2,is))-1                                 8d24s16
            do i3=0,nocc(isblkk(1,is))-1                                8d24s16
             do im=0,nvirtc(isblkkder(3,isb))-1                         8d24s16
              imp=im+nocc(isblkkder(3,isb))                             8d24s16
              iad1=itrans(isw)+i2p+nbasdwsc(isw)*imp                    8d25s16
              iad2=kmatda+i4+nocc(isblkk(2,is))*(i3                     8d24s16
     $            +nocc(isblkk(1,is))*(im+nvirtc(isblkkder(3,isb))*i1m))8d25s16
              bc(iad2)=bc(iad2)+bc(iad1)*bc(ii)                         8d24s16
             end do                                                     8d24s16
             ii=ii+1                                                    8d24s16
            end do                                                      8d24s16
           end do                                                       8d24s16
          end do                                                        8d24s16
          i10=1                                                         8d24s16
         end do                                                         8d24s16
         go to 407                                                      8d24s16
        end if
       end do
       write(6,*)('could not find kmat for (nb|ma'')'),
     $      isblkkder(1,isb),isblkkder(2,isb),isw,isblkkder(4,isb)
       call dws_sync
       call dws_finalize
       stop
  407  continue
       nsum=nall*2                                                      8d22s16
       call dws_gsumf(bc(jmatda),nsum)                                  8d22s16
       if(idwsdeb.gt.10)then
      end if
       call ilimts(nvirtc(isblkkder(3,isb)),nvirtc(isblkkder(4,isb)),   8d22s16
     $      mynprocg,mynowprog,il,ih,i1s,i1e,is,i2e)                    8d22s16
       do i12=il,ih                                                     8d22s16
        i12m=i12-1                                                      8d22s16
        i12k=i12-il                                                     8d22s16
        iad1=jmatd(isb)+i12k*nrow                                       8d22s16
        iad2=jmatda+i12m*nrow                                           8d22s16
        do i34=0,nrow-1                                                 8d22s16
         bc(iad1+i34)=bc(iad2+i34)                                      8d22s16
        end do                                                          8d22s16
        iad1=kmatd(isb)+i12k*nrow                                       8d22s16
        iad2=kmatda+i12m*nrow                                           8d22s16
        do i34=0,nrow-1                                                 8d22s16
         bc(iad1+i34)=bc(iad2+i34)                                      8d22s16
        end do                                                          8d22s16
       end do
       if(idwsdeb.gt.10)then                                            9d16s16
        write(6,12)(isblkkder(j,isb),j=1,4)
        ncol=ih+1-il
       end if                                                           9d16s16
       ibcoff=jmatda                                                    8d22s16
       end if                                                           11d30s16
      end do                                                            8d4s16
c
c     now do mixed 2nd der for 4o and 3ox. note, that in contrast to
c     j and k, 4o and 3ox ders are not global summed yet. also note
c     that 4o ders are only for one of the 8 symmetries, however for
c     the mixed 2nd, we lose the 8 fold symmetry, so we will explicitly
c     symmetrize this part, then divide by the appropriate symmetry
c     number.
c     nope. I decided I'd store the part needing to be symmetrized
c     separately from this that does not.
c
      do isb=1,nsblkder1                                                9d2s16
       nrow=nocc(isblkder1 (1,isb))*nocc(isblkder1 (2,isb))
       ncol=nocc(isblkder1 (3,isb))*nocc(isblkder1 (4,isb))
       i4od2b(isb)=ibcoff                                               9d16s16
       ibcoff=ibcoff+nrow*ncol                                          9d16s16
       call enough('parajkfromhd0.  8',bc,ibc)
       do i=0,nrow*ncol-1
        bc(i4od2b(isb)+i)=0d0
       end do
c
c     (o'o'|oo): contributions from (oo|oo), (vo|oo), (ov|oo), and (vv|oo)
c
       isw1=multh(isblkder1 (1,isb),iprop)                               9d2s16
       isw2=multh(isblkder1 (2,isb),iprop)                               9d2s16
       if(nocc(isw1)*nocc(isw2).gt.0)then                               11d30s16
       do is=1,nsdlk                                                    9d2s16
        if(isblk(1,is).eq.isw1.and.isblk(2,is).eq.isw2.and.             9d2s16
     $     isblk(3,is).eq.isblkder1 (3,isb).and.                         9d2s16
     $     isblk(4,is).eq.isblkder1 (4,isb))then                         9d2s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           if(isblk(1,is).eq.isblk(2,is))then                           9d2s16
            do i4=0,nocc(isblk(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d2s16
               iad2=itrans(isw2)+i3+nbasdwsc(isw2)*m2                   9d2s16
               jad2=itrans(isw2)+i4+nbasdwsc(isw2)*m2                   9d2s16
               fa=2d0*bc(ii)*bc(iad2)                                   9d16s16
               fb=2d0*bc(ii)*bc(jad2)                                   9d16s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i4+nbasdwsc(isw1)*m1                  9d2s16
                jad1=itrans(isw1)+i3+nbasdwsc(isw1)*m1                  9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2           9d2s16
     $             +nocc(isblkder1 (2,isb))*(i1m+nocc(isblk(3,is))*i2m))9d2s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa+bc(jad1)*fb               9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m2=0,nocc(isblkder1 (2,isb))-1                           9d2s16
              iad2=itrans(isw2)+i4+nbasdwsc(isw2)*m2                    9d2s16
              fa=2d0*bc(iad2)*bc(ii)                                    9d16s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d2s16
               iad1=itrans(isw1)+i4+nbasdwsc(isw1)*m1                   9d2s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2            9d2s16
     $             +nocc(isblkder1 (2,isb))*(i1m+nocc(isblk(3,is))*i2m))9d2s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk(2,is))-1                                 9d2s16
             do i3=0,nocc(isblk(1,is))-1                                9d2s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d2s16
               iad2=itrans(isw2)+i4+nbasdwsc(isw2)*m2                   9d2s16
               fa=2d0*bc(ii)*bc(iad2)                                   9d16s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i3+nbasdwsc(isw1)*m1                  9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2           9d2s16
     $             +nocc(isblkder1 (2,isb))*(i1m+nocc(isblk(3,is))*i2m))9d2s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
             ii=ii+1                                                    9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 500                                                      9d2s16
        end if                                                          9d2s16
        if(isblk(2,is).eq.isw1.and.isblk(1,is).eq.isw2.and.             9d8s16
     $     isblk(3,is).eq.isblkder1 (3,isb).and.                         9d2s16
     $     isblk(4,is).eq.isblkder1 (4,isb))then                         9d2s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m2=0,nocc(isblkder1 (2,isb))-1                           9d8s16
              iad2=itrans(isw2)+i3+nbasdwsc(isw2)*m2                    9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d8s16
               iad1=itrans(isw1)+i4+nbasdwsc(isw1)*m1                   9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $             +nocc(isblkder1 (2,isb))*(i1m+nocc(isblk(3,is))*i2m))9d2s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 500                                                      9d2s16
        end if                                                          9d2s16
        if(isblk(2,is).eq.isw1.and.isblk(1,is).eq.isw2.and.             9d8s16
     $     isblk(4,is).eq.isblkder1 (3,isb).and.                         9d8s16
     $     isblk(3,is).eq.isblkder1 (4,isb))then                         9d8s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m2=0,nocc(isblkder1 (2,isb))-1                           9d8s16
              iad2=itrans(isw2)+i3+nbasdwsc(isw2)*m2                    9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d8s16
               iad1=itrans(isw1)+i4+nbasdwsc(isw1)*m1                   9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $             +nocc(isblkder1 (2,isb))*(i2m+nocc(isblk(4,is))*i1m))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 500                                                      9d2s16
        end if                                                          9d2s16
        if(isblk(1,is).eq.isw1.and.isblk(2,is).eq.isw2.and.             9d8s16
     $     isblk(4,is).eq.isblkder1 (3,isb).and.                         9d8s16
     $     isblk(3,is).eq.isblkder1 (4,isb))then                         9d8s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m2=0,nocc(isblkder1 (2,isb))-1                           9d8s16
              iad2=itrans(isw2)+i4+nbasdwsc(isw2)*m2                    9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d8s16
               iad1=itrans(isw1)+i3+nbasdwsc(isw1)*m1                   9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $             +nocc(isblkder1 (2,isb))*(i2m+nocc(isblk(4,is))*i1m))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 500                                                      9d2s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       write(6,*)('no 4o for '),(isblkder1 (j,isb),j=1,4)                9d2s16
       write(6,*)('need '),isw1,isw2,(isblkder(j,isb),j=3,4)
       write(6,*)('got '),nocc(isw1),nocc(isw2)
       do is=1,nsdlk
        write(6,12)(isblk(j,is),j=1,4),is
       end do
       call dws_sync                                                    9d2s16
       call dws_finalize                                                9d2s16
       stop                                                             9d2s16
  500  continue                                                         9d2s16
       end if                                                           11d30s16
c
c     (o'o'|oo): (vo|oo) part
c
       if(nocc(isw2).gt.0)then                                          11d30s16
       do is=1,nsdlk1                                                   9d2s16
        if(isblk1(4,is).eq.isw1.and.isblk1(3,is).eq.isw2.and.           9d2s16
     $     isblk1(1,is).eq.isblkder1 (3,isb).and.                        9d2s16
     $     isblk1(2,is).eq.isblkder1 (4,isb))then                        9d2s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           if(isblk1(1,is).eq.isblk1(2,is))then                           9d2s16
            do i4=0,nocc(isblk1(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d2s16
               iad2=itrans(isw2)+i1m+nbasdwsc(isw2)*m2                  9d2s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                 9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2           9d2s16
     $              +nocc(isblkder1 (2,isb))*(i3+nocc(isblk1(1,is))*i4)) 9d2s16
                jad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2           9d2s16
     $              +nocc(isblkder1 (2,isb))*(i4+nocc(isblk1(1,is))*i3)) 9d2s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
                bc(jad3)=bc(jad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m2=0,nocc(isblkder1 (2,isb))-1                           9d2s16
              iad2=itrans(isw2)+i1m+nbasdwsc(isw2)*m2                   9d2s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d2s16
               iad1=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                  9d2s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2            9d2s16
     $              +nocc(isblkder1 (2,isb))*(i4+nocc(isblk1(1,is))*i4)) 9d2s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk1(2,is))-1                                9d2s16
             do i3=0,nocc(isblk1(1,is))-1                               9d2s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d2s16
               iad2=itrans(isw2)+i1m+nbasdwsc(isw2)*m2                  9d2s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                 9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2           9d2s16
     $              +nocc(isblkder1 (2,isb))*(i3+nocc(isblk1(1,is))*i4)) 9d2s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 501                                                      9d2s16
        end if                                                          9d2s16
        if(isblk1(4,is).eq.isw1.and.isblk1(3,is).eq.isw2.and.           9d2s16
     $     isblk1(2,is).eq.isblkder1 (3,isb).and.                        9d8s16
     $     isblk1(1,is).eq.isblkder1 (4,isb))then                        9d8s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk1(2,is))-1                                 9d8s16
            do i3=0,nocc(isblk1(1,is))-1                                9d8s16
             do m2=0,nocc(isblkder1 (2,isb))-1                           9d8s16
              iad2=itrans(isw2)+i1m+nbasdwsc(isw2)*m2                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d8s16
               iad1=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                  9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $              +nocc(isblkder1 (2,isb))*(i4+nocc(isblk1(2,is))*i3)) 9d2s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 501                                                      9d2s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       if(min(nbasdwsc(isw1),nocc(isw2)).gt.0)then                      11d28s22
        write(6,*)('no onex (vo|oo) for (o''o''|oo) '),                  9d2s16
     $      (isblkder1 (j,isb),j=1,4)                                    9d2s16
        write(6,*)('what we have for sblk1: ')
        write(6,*)('nsdlk1: '),nsdlk1
        do is=1,nsdlk1                                                   6d13s22
         write(6,*)is,(isblk1(j,is),j=1,4)
        end do                                                           6d13s22
        call dws_sync                                                    9d2s16
        call dws_finalize                                                9d2s16
        stop                                                             9d2s16
       end if                                                           11d28s22o
  501  continue                                                         9d2s16
       end if
c
c     (o'o'|oo): (ov|oo) part
c
       if(min(nocc(isw1),nbasdwsc(isw2)).gt.0)then                      11d28s22
       do is=1,nsdlk1                                                   9d2s16
        if(isblk1(4,is).eq.isw2.and.isblk1(3,is).eq.isw1.and.           9d2s16
     $     isblk1(1,is).eq.isblkder1 (3,isb).and.                        9d2s16
     $     isblk1(2,is).eq.isblkder1 (4,isb))then                        9d2s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d8s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d8s16
           if(isblk1(1,is).eq.isblk1(2,is))then                           9d2s16
            do i4=0,nocc(isblk1(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d2s16
               iad2=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                  9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i1m+nbasdwsc(isw1)*m1                 9d8s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2           9d2s16
     $              +nocc(isblkder1 (2,isb))*(i3+nocc(isblk1(1,is))*i4)) 9d2s16
                jad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2           9d2s16
     $              +nocc(isblkder1 (2,isb))*(i4+nocc(isblk1(1,is))*i3)) 9d2s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
                bc(jad3)=bc(jad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m2=0,nocc(isblkder1 (2,isb))-1                           9d2s16
              iad2=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                   9d8s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d2s16
               iad1=itrans(isw1)+i1m+nbasdwsc(isw1)*m1                  9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2            9d2s16
     $              +nocc(isblkder1 (2,isb))*(i4+nocc(isblk1(1,is))*i4)) 9d2s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk1(2,is))-1                                9d2s16
             do i3=0,nocc(isblk1(1,is))-1                               9d2s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d2s16
               iad2=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                  9d9s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i1m+nbasdwsc(isw1)*m1                 9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2           9d2s16
     $              +nocc(isblkder1 (2,isb))*(i3+nocc(isblk1(1,is))*i4)) 9d2s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 502                                                      9d2s16
        end if                                                          9d2s16
        if(isblk1(4,is).eq.isw2.and.isblk1(3,is).eq.isw1.and.           9d2s16
     $     isblk1(2,is).eq.isblkder1 (3,isb).and.                        9d8s16
     $     isblk1(1,is).eq.isblkder1 (4,isb))then                        9d8s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d8s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d8s16
           do i4=0,nocc(isblk1(2,is))-1                                 9d8s16
            do i3=0,nocc(isblk1(1,is))-1                                9d8s16
             do m2=0,nocc(isblkder1 (2,isb))-1                           9d8s16
              iad2=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d8s16
               iad1=itrans(isw1)+i1m+nbasdwsc(isw1)*m1                  9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $              +nocc(isblkder1 (2,isb))*(i4+nocc(isblk1(2,is))*i3)) 9d2s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 502                                                      9d2s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       write(6,*)('no onex (ov|oo)  for (o''o''|oo) '),                 9d2s16
     $      (isblkder1 (j,isb),j=1,4)                                    9d2s16
       call dws_sync                                                    9d2s16
       call dws_finalize                                                9d2s16
       stop                                                             9d2s16
  502  continue                                                         9d2s16
       end if                                                           11d30s16
c
c     (o'o'|oo): (vv|oo) part
c
       do is=1,nsdlk                                                    9d2s16
        if(isblk(4,is).eq.isw2.and.isblk(3,is).eq.isw1.and.             9d2s16
     $     isblk(1,is).eq.isblkder1 (3,isb).and.                         9d2s16
     $     isblk(2,is).eq.isblkder1 (4,isb))then                         9d2s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d2s16
         ii=jmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk(4,is))                                    9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblk(3,is))                                   9d2s16
           if(isblk(1,is).eq.isblk(2,is))then                           9d2s16
            do i4=0,nocc(isblk(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d2s16
               iad2=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                  9d2s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i1p+nbasdwsc(isw1)*m1                 9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2           9d2s16
     $              +nocc(isblkder1 (2,isb))*(i3+nocc(isblk(1,is))*i4))  9d2s16
                jad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2           9d2s16
     $              +nocc(isblkder1 (2,isb))*(i4+nocc(isblk(1,is))*i3))  9d2s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
                bc(jad3)=bc(jad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m2=0,nocc(isblkder1 (2,isb))-1                           9d2s16
              iad2=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                   9d2s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d2s16
               iad1=itrans(isw1)+i1p+nbasdwsc(isw1)*m1                  9d2s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2            9d2s16
     $              +nocc(isblkder1 (2,isb))*(i4+nocc(isblk(1,is))*i4))  9d2s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk(2,is))-1                                 9d2s16
             do i3=0,nocc(isblk(1,is))-1                                9d2s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d2s16
               iad2=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                  9d2s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i1p+nbasdwsc(isw1)*m1                 9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2           9d2s16
     $              +nocc(isblkder1 (2,isb))*(i3+nocc(isblk(1,is))*i4))  9d2s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 503                                                      9d2s16
        end if                                                          9d2s16
        if(isblk(4,is).eq.isw2.and.isblk(3,is).eq.isw1.and.             9d2s16
     $     isblk(2,is).eq.isblkder1 (3,isb).and.                         9d8s16
     $       isblk(1,is).eq.isblkder1 (4,isb))then                         9d8s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d2s16
         ii=jmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk(4,is))                                    9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblk(3,is))                                   9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m2=0,nocc(isblkder1 (2,isb))-1                           9d8s16
              iad2=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d8s16
               iad1=itrans(isw1)+i1p+nbasdwsc(isw1)*m1                  9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $              +nocc(isblkder1 (2,isb))*(i4+nocc(isblk(2,is))*i3))  9d2s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 503                                                      9d2s16
        end if                                                          9d2s16
        if(isblk(3,is).eq.isw2.and.isblk(4,is).eq.isw1.and.             9d8s16
     $     isblk(2,is).eq.isblkder1 (3,isb).and.                         9d8s16
     $     isblk(1,is).eq.isblkder1 (4,isb))then                         9d8s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d2s16
         ii=jmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk(4,is))                                    9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblk(3,is))                                   9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m2=0,nocc(isblkder1 (2,isb))-1                           9d8s16
              iad2=itrans(isw2)+i1p+nbasdwsc(isw2)*m2                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d8s16
               iad1=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                  9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $              +nocc(isblkder1 (2,isb))*(i4+nocc(isblk(2,is))*i3))  9d2s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 503                                                      9d2s16
        end if                                                          9d2s16
        if(isblk(3,is).eq.isw2.and.isblk(4,is).eq.isw1.and.             9d8s16
     $     isblk(1,is).eq.isblkder1 (3,isb).and.                         9d8s16
     $     isblk(2,is).eq.isblkder1 (4,isb))then                         9d8s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d2s16
         ii=jmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk(4,is))                                    9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblk(3,is))                                   9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m2=0,nocc(isblkder1 (2,isb))-1                           9d8s16
              iad2=itrans(isw2)+i1p+nbasdwsc(isw2)*m2                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d8s16
               iad1=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                  9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $              +nocc(isblkder1 (2,isb))*(i3+nocc(isblk(1,is))*i4))  9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 503                                                      9d2s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       if(min(nbasdwsc(isw1),nbasdwsc(isw2)).gt.0)then                  11d28s22
        write(6,*)('no jmats (vv|oo)  for (o''o''|oo) '),                 9d2s16
     $      (isblkder1 (j,isb),j=1,4),isw1,isw2                                    9d2s16
        call dws_sync                                                    9d2s16
        call dws_finalize                                                9d2s16
        stop                                                             9d2s16
       end if                                                           11d28s22
  503  continue                                                         9d2s16
c
c     (o'o|o'o), (oo|oo) part
c
       isw1=multh(isblkder1 (1,isb),iprop)                               9d2s16
       isw3=multh(isblkder1 (3,isb),iprop)                               9d2s16
       if(nocc(isw1)*nocc(isw3).ne.0)then                               11d30s16
       do is=1,nsdlk                                                    9d2s16
        if(isblk(1,is).eq.isw1.and.isblk(2,is).eq.isblkder1 (2,isb).and. 9d2s16
     $     isblk(3,is).eq.isw3.and.                                     9d2s16
     $     isblk(4,is).eq.isblkder1 (4,isb))then                         9d2s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           if(isblk(1,is).eq.isblk(2,is))then                           9d2s16
            do i4=0,nocc(isblk(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d2s16
               iad2=itrans(isw3)+i1m+nbasdwsc(isw3)*m3                  9d2s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i4+nbasdwsc(isw1)*m1                  9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i3           9d2s16
     $              +nocc(isblk(2,is))*(m3+nocc(isblkder1 (3,isb))*i2m)) 9d2s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
                iad1=itrans(isw1)+i3+nbasdwsc(isw1)*m1                  9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i4           9d2s16
     $              +nocc(isblk(2,is))*(m3+nocc(isblkder1 (3,isb))*i2m)) 9d2s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d2s16
              iad2=itrans(isw3)+i1m+nbasdwsc(isw3)*m3                   9d2s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d2s16
               iad1=itrans(isw1)+i4+nbasdwsc(isw1)*m1                   9d2s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i4            9d2s16
     $              +nocc(isblk(2,is))*(m3+nocc(isblkder1 (3,isb))*i2m)) 9d7s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk(2,is))-1                                 9d2s16
             do i3=0,nocc(isblk(1,is))-1                                9d2s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d2s16
               iad2=itrans(isw3)+i1m+nbasdwsc(isw3)*m3                  9d2s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i3+nbasdwsc(isw1)*m1                  9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i4           9d2s16
     $              +nocc(isblk(2,is))*(m3+nocc(isblkder1 (3,isb))*i2m)) 9d7s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
             ii=ii+1                                                    9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 504                                                      9d2s16
        end if                                                          9d2s16
        if(isblk(2,is).eq.isw1.and.isblk(1,is).eq.isblkder1 (2,isb).and. 9d8s16
     $     isblk(3,is).eq.isw3.and.                                     9d2s16
     $     isblk(4,is).eq.isblkder1 (4,isb))then                         9d2s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d8s16
              iad2=itrans(isw3)+i1m+nbasdwsc(isw3)*m3                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d8s16
               iad1=itrans(isw1)+i4+nbasdwsc(isw1)*m1                   9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i3            9d8s16
     $              +nocc(isblk(1,is))*(m3+nocc(isblkder1 (3,isb))*i2m)) 9d7s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 504                                                      9d2s16
        end if                                                          9d2s16
        if(isblk(2,is).eq.isw1.and.isblk(1,is).eq.isblkder1 (2,isb).and. 9d8s16
     $     isblk(4,is).eq.isw3.and.                                     9d8s16
     $     isblk(3,is).eq.isblkder1 (4,isb))then                         9d8s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d8s16
              iad2=itrans(isw3)+i2m+nbasdwsc(isw3)*m3                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d8s16
               iad1=itrans(isw1)+i4+nbasdwsc(isw1)*m1                   9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i3            9d8s16
     $              +nocc(isblk(1,is))*(m3+nocc(isblkder1 (3,isb))*i1m)) 9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 504                                                      9d2s16
        end if                                                          9d2s16
        if(isblk(1,is).eq.isw1.and.isblk(2,is).eq.isblkder1 (2,isb).and. 9d8s16
     $     isblk(4,is).eq.isw3.and.                                     9d8s16
     $     isblk(3,is).eq.isblkder1 (4,isb))then                         9d8s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d8s16
              iad2=itrans(isw3)+i2m+nbasdwsc(isw3)*m3                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d8s16
               iad1=itrans(isw1)+i3+nbasdwsc(isw1)*m1                   9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i4            9d8s16
     $              +nocc(isblk(2,is))*(m3+nocc(isblkder1 (3,isb))*i1m)) 9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 504                                                      9d2s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       write(6,*)('no 4o for (o''o|o''o)'),(isblkder1 (j,isb),j=1,4)     9d2s16
       call dws_sync                                                    9d2s16
       call dws_finalize                                                9d2s16
       stop                                                             9d2s16
  504  continue                                                         9d2s16
       end if
c
c     (o'o|o'o): (vo|oo) part
c
       if(nocc(isw3).ne.0)then                                          11d30s16
       do is=1,nsdlk1                                                   9d2s16
        if(isblk1(4,is).eq.isw1.and.isblk1(3,is).eq.isblkder1 (2,isb)    9d2s16
     $       .and.isblk1(1,is).eq.isw3.and.                             9d2s16
     $     isblk1(2,is).eq.isblkder1 (4,isb))then                        9d2s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           if(isblk1(1,is).eq.isblk1(2,is))then                           9d2s16
            do i4=0,nocc(isblk1(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d2s16
               iad2=itrans(isw3)+i3+nbasdwsc(isw3)*m3                   9d2s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               jad2=itrans(isw3)+i4+nbasdwsc(isw3)*m3                   9d2s16
               fb=2d0*bc(ii)*bc(jad2)                                       9d2s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                 9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i1m          9d2s16
     $              +nocc(isblk1(3,is))*(m3+nocc(isblkder1 (3,isb))*i4)) 9d2s16
                jad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i1m          9d2s16
     $              +nocc(isblk1(3,is))*(m3+nocc(isblkder1 (3,isb))*i3)) 9d2s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
                bc(jad3)=bc(jad3)+bc(iad1)*fb                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d2s16
              iad2=itrans(isw3)+i4+nbasdwsc(isw3)*m3                    9d2s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d2s16
               iad1=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                  9d2s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i1m           9d2s16
     $              +nocc(isblk1(3,is))*(m3+nocc(isblkder1 (3,isb))*i4)) 9d2s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk1(2,is))-1                                9d2s16
             do i3=0,nocc(isblk1(1,is))-1                               9d2s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d2s16
               iad2=itrans(isw3)+i3+nbasdwsc(isw3)*m3                   9d2s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                 9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i1m          9d2s16
     $              +nocc(isblk1(3,is))*(m3+nocc(isblkder1 (3,isb))*i4)) 9d2s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 505                                                      9d2s16
        end if                                                          9d2s16
        if(isblk1(4,is).eq.isw1.and.isblk1(3,is).eq.isblkder1 (2,isb)    9d2s16
     $       .and.isblk1(2,is).eq.isw3.and.                             9d8s16
     $     isblk1(1,is).eq.isblkder1 (4,isb))then                        9d8s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk1(2,is))-1                                 9d8s16
            do i3=0,nocc(isblk1(1,is))-1                                9d8s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d8s16
              iad2=itrans(isw3)+i4+nbasdwsc(isw3)*m3                    9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d8s16
               iad1=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                  9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i1m           9d8s16
     $              +nocc(isblk1(3,is))*(m3+nocc(isblkder1 (3,isb))*i3)) 9d2s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 505                                                      9d2s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       if(nbasdwsc(isw1).gt.0)then                                      11d28s22
        write(6,*)('no onex (vo|oo) for (o''o|o''o) '),                  9d2s16
     $      (isblkder1 (j,isb),j=1,4),isw1
        call dws_sync                                                    9d2s16
        call dws_finalize                                                9d2s16
        stop                                                             9d2s16
       end if                                                           11d28s22
  505  continue                                                         9d2s16
       end if                                                           11d30s16
c
c     (o'o|o'o): (oo|vo) part
c
       if(min(nocc(isw1),nbasdwsc(isw3)).gt.0)then                      11d28s22
       do is=1,nsdlk1                                                   9d2s16
        if(isblk1(4,is).eq.isw3.and.isblk1(3,is).eq.isblkder1 (4,isb)    9d2s16
     $       .and.isblk1(1,is).eq.isw1.and.                             9d2s16
     $     isblk1(2,is).eq.isblkder1 (2,isb))then                        9d2s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1
           if(isblk1(1,is).eq.isblk1(2,is))then                           9d2s16
            do i4=0,nocc(isblk1(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d2s16
               iad2=itrans(isw3)+i2p+nbasdwsc(isw3)*m3                  9d2s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i3+nbasdwsc(isw1)*m1                  9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i4           9d2s16
     $             +nocc(isblk1(2,is))*(m3+nocc(isblkder1 (3,isb))*i1m)) 9d2s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
                iad1=itrans(isw1)+i4+nbasdwsc(isw1)*m1                  9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i3           9d2s16
     $             +nocc(isblk1(2,is))*(m3+nocc(isblkder1 (3,isb))*i1m))9d2s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d2s16
              iad2=itrans(isw3)+i2p+nbasdwsc(isw3)*m3                   9d2s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d2s16
               iad1=itrans(isw1)+i4+nbasdwsc(isw1)*m1                   9d2s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i4            9d2s16
     $             +nocc(isblk1(2,is))*(m3+nocc(isblkder1 (3,isb))*i1m))9d2s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk1(2,is))-1                                9d2s16
             do i3=0,nocc(isblk1(1,is))-1                               9d2s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d2s16
               iad2=itrans(isw3)+i2p+nbasdwsc(isw3)*m3                  9d2s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i3+nbasdwsc(isw1)*m1                  9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i4           9d2s16
     $             +nocc(isblk1(2,is))*(m3+nocc(isblkder1 (3,isb))*i1m))9d2s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 506                                                      9d2s16
        end if                                                          9d2s16
        if(isblk1(4,is).eq.isw3.and.isblk1(3,is).eq.isblkder1 (4,isb)    9d2s16
     $       .and.isblk1(2,is).eq.isw1.and.                             9d8s16
     $     isblk1(1,is).eq.isblkder1 (2,isb))then                        9d8s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1
           do i4=0,nocc(isblk1(2,is))-1                                 9d8s16
            do i3=0,nocc(isblk1(1,is))-1                                9d8s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d8s16
              iad2=itrans(isw3)+i2p+nbasdwsc(isw3)*m3                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d8s16
               iad1=itrans(isw1)+i4+nbasdwsc(isw1)*m1                   9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i3            9d8s16
     $             +nocc(isblk1(1,is))*(m3+nocc(isblkder1 (3,isb))*i1m))9d2s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 506                                                      9d2s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       write(6,*)('no onex (oo|vo)  for (o''o|o''o) '),                 9d2s16
     $      (isblkder1 (j,isb),j=1,4)                                    9d2s16
       call dws_sync                                                    9d2s16
       call dws_finalize                                                9d2s16
       stop                                                             9d2s16
  506  continue                                                         9d2s16
       end if                                                           11d30s16
c
c     (o'o|o'o): (vo|vo) part
c     recall k_nm^ab=(nb|ma)
c
       if(min(nbasdwsc(isw1),nbasdwsc(isw3)).eq.0)go to 507             11d28s22
       do is=1,nsdlkk                                                   9d2s16
        if(isblkk(4,is).eq.isw1.and.isblkk(3,is).eq.isw3                9d2s16
     $       .and.isblkk(2,is).eq.isblkder1 (4,isb).and.                         9d2s16
     $     isblkk(1,is).eq.isblkder1 (2,isb))then                         9d2s16
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d2s16
         ii=kmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblkk(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblkk(3,is))                                   9d2s16
           do i4=0,nocc(isblkk(2,is))-1                                  9d2s16
            do i3=0,nocc(isblkk(1,is))-1                                 9d2s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d2s16
              iad2=itrans(isw3)+i1p+nbasdwsc(isw3)*m3                   9d2s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d2s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d2s16
               iad1=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                  9d2s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i3            9d2s16
     $             +nocc(isblkk(1,is))*(m3+nocc(isblkder1 (3,isb))*i4))  9d2s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           end do                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 507                                                      9d2s16
        end if                                                          9d2s16
        if(isblkk(3,is).eq.isw1.and.isblkk(4,is).eq.isw3                9d8s16
     $       .and.isblkk(1,is).eq.isblkder1 (4,isb).and.                 9d8s16
     $     isblkk(2,is).eq.isblkder1 (2,isb))then                        9d8s16
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d2s16
         ii=kmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblkk(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblkk(3,is))                                   9d2s16
           do i4=0,nocc(isblkk(2,is))-1                                  9d2s16
            do i3=0,nocc(isblkk(1,is))-1                                 9d2s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d2s16
              iad2=itrans(isw3)+i2p+nbasdwsc(isw3)*m3                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d2s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d2s16
               iad1=itrans(isw1)+i1p+nbasdwsc(isw1)*m1                  9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i4            9d8s16
     $             +nocc(isblkk(2,is))*(m3+nocc(isblkder1 (3,isb))*i3))  9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           end do                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 507                                                      9d2s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       write(6,*)('no kmats  for (o''o|o''o) '),                        9d2s16
     $      (isblkder1 (j,isb),j=1,4)                                    9d2s16
       call dws_sync                                                    9d2s16
       call dws_finalize                                                9d2s16
       stop                                                             9d2s16
  507  continue                                                         9d2s16
c
c     (o'o|oo'), (oo|oo) part
c
       isw1=multh(isblkder1 (1,isb),iprop)                               9d2s16
       isw4=multh(isblkder1 (4,isb),iprop)                               9d2s16
       do is=1,nsdlk                                                    9d2s16
        if(isblk(1,is).eq.isw1.and.isblk(2,is).eq.isblkder1 (2,isb).and. 9d2s16
     $     isblk(3,is).eq.isblkder1 (3,isb).and.                         9d2s16
     $     isblk(4,is).eq.isw4)then                                     9d2s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           if(isblk(1,is).eq.isblk(2,is))then                           9d2s16
            do i4=0,nocc(isblk(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d2s16
               iad2=itrans(isw4)+i2m+nbasdwsc(isw4)*m4                  9d2s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i4+nbasdwsc(isw1)*m1                  9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i3           9d2s16
     $              +nocc(isblk(2,is))*(i1m+nocc(isblk(3,is))*m4))      9d2s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
                iad1=itrans(isw1)+i3+nbasdwsc(isw1)*m1                  9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i4           9d2s16
     $              +nocc(isblk(2,is))*(i1m+nocc(isblk(3,is))*m4))      9d2s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d2s16
              iad2=itrans(isw4)+i2m+nbasdwsc(isw4)*m4                   9d2s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d2s16
               iad1=itrans(isw1)+i4+nbasdwsc(isw1)*m1                   9d2s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i4            9d2s16
     $              +nocc(isblk(2,is))*(i1m+nocc(isblk(3,is))*m4))      9d7s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk(2,is))-1                                 9d2s16
             do i3=0,nocc(isblk(1,is))-1                                9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d2s16
               iad2=itrans(isw4)+i2m+nbasdwsc(isw4)*m4                  9d2s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i3+nbasdwsc(isw1)*m1                  9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i4           9d2s16
     $              +nocc(isblk(2,is))*(i1m+nocc(isblk(3,is))*m4))      9d7s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
             ii=ii+1                                                    9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 508                                                      9d2s16
        end if                                                          9d2s16
        if(isblk(2,is).eq.isw1.and.isblk(1,is).eq.isblkder1 (2,isb).and. 9d8s16
     $     isblk(3,is).eq.isblkder1 (3,isb).and.                         9d2s16
     $     isblk(4,is).eq.isw4)then                                     9d2s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i2m+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d8s16
               iad1=itrans(isw1)+i4+nbasdwsc(isw1)*m1                   9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i3            9d8s16
     $              +nocc(isblk(1,is))*(i1m+nocc(isblk(3,is))*m4))      9d7s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 508                                                      9d2s16
        end if                                                          9d2s16
        if(isblk(2,is).eq.isw1.and.isblk(1,is).eq.isblkder1 (2,isb).and. 9d8s16
     $     isblk(4,is).eq.isblkder1 (3,isb).and.                         9d2s16
     $     isblk(3,is).eq.isw4)then                                     9d2s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i1m+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d8s16
               iad1=itrans(isw1)+i4+nbasdwsc(isw1)*m1                   9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i3            9d8s16
     $              +nocc(isblk(1,is))*(i2m+nocc(isblk(4,is))*m4))      9d9s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 508                                                      9d2s16
        end if                                                          9d2s16
        if(isblk(1,is).eq.isw1.and.isblk(2,is).eq.isblkder1 (2,isb).and. 9d8s16
     $     isblk(4,is).eq.isblkder1 (3,isb).and.                         9d2s16
     $     isblk(3,is).eq.isw4)then                                     9d2s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i1m+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d8s16
               iad1=itrans(isw1)+i3+nbasdwsc(isw1)*m1                   9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i4            9d8s16
     $              +nocc(isblk(2,is))*(i2m+nocc(isblk(4,is))*m4))      9d9s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 508                                                      9d2s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       write(6,*)('no 4o for (o''o|oo'')'),(isblkder1 (j,isb),j=1,4)     9d2s16
       call dws_sync                                                    9d2s16
       call dws_finalize                                                9d2s16
       stop                                                             9d2s16
  508  continue                                                         9d2s16
c
c     (o'o|oo'): (vo|oo) part
c
       if(min(nbasdwsc(isw1),nocc(isw4)).eq.0)go to 509                 11d28s22
       do is=1,nsdlk1                                                   9d2s16
        if(isblk1(4,is).eq.isw1.and.isblk1(3,is).eq.isblkder1 (2,isb)    9d2s16
     $       .and.isblk1(1,is).eq.isblkder1 (3,isb).and.                 9d2s16
     $     isblk1(2,is).eq.isw4)then                                    9d2s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           if(isblk1(1,is).eq.isblk1(2,is))then                           9d2s16
            do i4=0,nocc(isblk1(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d2s16
               iad2=itrans(isw4)+i3+nbasdwsc(isw4)*m4                   9d2s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               jad2=itrans(isw4)+i4+nbasdwsc(isw4)*m4                   9d2s16
               fb=2d0*bc(ii)*bc(jad2)                                       9d2s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                 9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i1m         9d9s16
     $              +nocc(isblk1(3,is))*(i4+nocc(isblk1(1,is))*m4))     9d9s16
                jad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i1m         9d9s16
     $              +nocc(isblk1(3,is))*(i3+nocc(isblk1(1,is))*m4))     9d9s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
                bc(jad3)=bc(jad3)+bc(iad1)*fb                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d2s16
              iad2=itrans(isw4)+i4+nbasdwsc(isw4)*m4                    9d2s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d2s16
               iad1=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                  9d2s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i1m          9d9s16
     $              +nocc(isblk1(3,is))*(i4+nocc(isblk1(1,is))*m4))     9d9s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk1(2,is))-1                                9d2s16
             do i3=0,nocc(isblk1(1,is))-1                               9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d2s16
               iad2=itrans(isw4)+i4+nbasdwsc(isw4)*m4                   9d2s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                 9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i1m          9d2s16
     $              +nocc(isblk1(3,is))*(i3+nocc(isblk1(1,is))*m4))     9d2s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 509                                                      9d2s16
        end if                                                          9d2s16
        if(isblk1(4,is).eq.isw1.and.isblk1(3,is).eq.isblkder1 (2,isb)    9d2s16
     $       .and.isblk1(2,is).eq.isblkder1 (3,isb).and.                 9d8s16
     $     isblk1(1,is).eq.isw4)then                                    9d8s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk1(2,is))-1                                 9d8s16
            do i3=0,nocc(isblk1(1,is))-1                                9d8s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i3+nbasdwsc(isw4)*m4                    9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d8s16
               iad1=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                  9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i1m           9d8s16
     $              +nocc(isblk1(3,is))*(i4+nocc(isblk1(2,is))*m4))     9d2s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 509                                                      9d2s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       if(nbasdwsc(isw1).gt.0)then                                      11d28s22
        write(6,*)('no onex (vo|oo) for (o''o|oo'') '),                  9d2s16
     $      (isblkder1 (j,isb),j=1,4),isw1                                    9d2s16
        call dws_sync                                                    9d2s16
        call dws_finalize                                                9d2s16
        stop                                                             9d2s16
       end if                                                           11d28s22
  509  continue                                                         9d2s16
c
c     (o'o|oo'): (oo|ov) part
c
       do is=1,nsdlk1                                                   9d2s16
        if(isblk1(4,is).eq.isw4.and.isblk1(3,is).eq.isblkder1 (3,isb)    9d7s16
     $       .and.isblk1(1,is).eq.isw1.and.                             9d2s16
     $     isblk1(2,is).eq.isblkder1 (2,isb))then                        9d2s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1
           if(isblk1(1,is).eq.isblk1(2,is))then                           9d2s16
            do i4=0,nocc(isblk1(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d7s16
               iad2=itrans(isw4)+i2p+nbasdwsc(isw4)*m4                  9d7s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i3+nbasdwsc(isw1)*m1                  9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i4           9d2s16
     $              +nocc(isblk1(2,is))*(i1m+nocc(isblk1(3,is))*m4))    9d7s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
                iad1=itrans(isw1)+i4+nbasdwsc(isw1)*m1                  9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i3           9d2s16
     $              +nocc(isblk1(2,is))*(i1m+nocc(isblk1(3,is))*m4))    9d7s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d7s16
              iad2=itrans(isw4)+i2p+nbasdwsc(isw4)*m4                   9d7s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d2s16
               iad1=itrans(isw1)+i4+nbasdwsc(isw1)*m1                   9d2s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i4            9d2s16
     $              +nocc(isblk1(2,is))*(i1m+nocc(isblk1(3,is))*m4))    9d9s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk1(2,is))-1                                9d2s16
             do i3=0,nocc(isblk1(1,is))-1                               9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d7s16
               iad2=itrans(isw4)+i2p+nbasdwsc(isw4)*m4                  9d7s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m1=0,nocc(isblkder1 (1,isb))-1                         9d2s16
                iad1=itrans(isw1)+i3+nbasdwsc(isw1)*m1                  9d2s16
                iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i4           9d2s16
     $              +nocc(isblk1(2,is))*(i1m+nocc(isblk1(3,is))*m4))    9d9s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 510                                                      9d2s16
        end if                                                          9d2s16
        if(isblk1(4,is).eq.isw4.and.isblk1(3,is).eq.isblkder1 (3,isb)    9d7s16
     $       .and.isblk1(2,is).eq.isw1.and.                             9d8s16
     $     isblk1(1,is).eq.isblkder1 (2,isb))then                        9d8s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1
           do i4=0,nocc(isblk1(2,is))-1                                 9d8s16
            do i3=0,nocc(isblk1(1,is))-1                                9d8s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i2p+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d8s16
               iad1=itrans(isw1)+i4+nbasdwsc(isw1)*m1                   9d8s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i3            9d8s16
     $              +nocc(isblk1(1,is))*(i1m+nocc(isblk1(3,is))*m4))    9d9s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 510                                                      9d2s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       if(min(nocc(isw1),nbasdwsc(isw4)).gt.0)then                      11d28s22
        write(6,*)('no onex (oo|ov)  for (o''o|oo'') '),                 9d7s16
     $      (isblkder1 (j,isb),j=1,4),isw1,isw4                              11d28s22
        call dws_sync                                                    9d2s16
        call dws_finalize                                                9d2s16
        stop                                                             9d2s16
       end if                                                           11d28s22
  510  continue                                                         9d2s16
c
c     (o'o|oo'): (vo|ov) part
c     recall k_nm^ab=(nb|ma)
c
       if(min(nbasdwsc(isw1),nbasdwsc(isw4)).eq.0)go to 511             11d28s22
       do is=1,nsdlkk                                                   9d2s16
        if(isblkk(4,is).eq.isw1.and.isblkk(3,is).eq.isw4                9d7s16
     $       .and.isblkk(2,is).eq.isblkder1 (3,isb).and.                 9d7s16
     $     isblkk(1,is).eq.isblkder1 (2,isb))then                         9d2s16
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d2s16
         ii=kmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblkk(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblkk(3,is))                                   9d2s16
           do i4=0,nocc(isblkk(2,is))-1                                  9d2s16
            do i3=0,nocc(isblkk(1,is))-1                                 9d2s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d7s16
              iad2=itrans(isw4)+i1p+nbasdwsc(isw4)*m4                   9d7s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d2s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d2s16
               iad1=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                  9d2s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i3            9d2s16
     $             +nocc(isblkk(1,is))*(i4+nocc(isblkk(2,is))*m4))      9d7s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           end do                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 511                                                      9d7s16
        end if                                                          9d2s16
        if(isblkk(3,is).eq.isw1.and.isblkk(4,is).eq.isw4                9d8s16
     $       .and.isblkk(1,is).eq.isblkder1 (3,isb).and.                 9d8s16
     $     isblkk(2,is).eq.isblkder1 (2,isb))then                         9d2s16
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d2s16
         ii=kmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblkk(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblkk(3,is))                                   9d2s16
           do i4=0,nocc(isblkk(2,is))-1                                  9d2s16
            do i3=0,nocc(isblkk(1,is))-1                                 9d2s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d7s16
              iad2=itrans(isw4)+i2p+nbasdwsc(isw4)*m4                   9d7s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d2s16
              do m1=0,nocc(isblkder1 (1,isb))-1                          9d2s16
               iad1=itrans(isw1)+i1p+nbasdwsc(isw1)*m1                  9d2s16
               iad3=i4od2b(isb)+m1+nocc(isblkder1 (1,isb))*(i4            9d2s16
     $             +nocc(isblkk(2,is))*(i3+nocc(isblkk(1,is))*m4))      9d7s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           end do                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 511                                                      9d7s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       write(6,*)('no kmats  for (o''o|oo'') '),                        9d7s16
     $      (isblkder1 (j,isb),j=1,4)                                    9d2s16
       call dws_sync                                                    9d2s16
       call dws_finalize                                                9d2s16
       stop                                                             9d2s16
  511  continue                                                         9d7s16
c
c     (oo'|o'o), (oo|oo) part                                           9d7s16
c
       isw2=multh(isblkder1 (2,isb),iprop)                               9d7s16
       isw3=multh(isblkder1 (3,isb),iprop)                               9d7s16
       do is=1,nsdlk                                                    9d2s16
        if(isblk(1,is).eq.isblkder1 (1,isb).and.isblk(2,is).eq.isw2.and. 9d7s16
     $     isblk(3,is).eq.isw3.and.isblk(4,is).eq.isblkder1 (4,isb))then
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           if(isblk(1,is).eq.isblk(2,is))then                           9d2s16
            do i4=0,nocc(isblk(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d2s16
               iad2=itrans(isw3)+i1m+nbasdwsc(isw3)*m3                  9d7s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m2=0,nocc(isblkder1 (2,isb))-1                         9d7s16
                iad1=itrans(isw2)+i4+nbasdwsc(isw2)*m2                  9d7s16
                iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(m2           9d7s16
     $        +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i2m))9d7s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
                iad1=itrans(isw2)+i3+nbasdwsc(isw2)*m2                  9d7s16
                iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(m2           9d7s16
     $        +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i2m))9d7s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d7s16
              iad2=itrans(isw3)+i1m+nbasdwsc(isw3)*m3                   9d7s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d7s16
               iad1=itrans(isw2)+i4+nbasdwsc(isw2)*m2                   9d7s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(m2            9d7s16
     $        +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i2m))9d7s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk(2,is))-1                                 9d2s16
             do i3=0,nocc(isblk(1,is))-1                                9d2s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d7s16
               iad2=itrans(isw3)+i1m+nbasdwsc(isw3)*m3                  9d7s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m2=0,nocc(isblkder1 (2,isb))-1                         9d7s16
                iad1=itrans(isw2)+i4+nbasdwsc(isw2)*m2                  9d7s16
                iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(m2           9d7s16
     $        +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i2m))9d7s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
             ii=ii+1                                                    9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 512                                                      9d7s16
        end if                                                          9d2s16
        if(isblk(2,is).eq.isblkder1 (1,isb).and.isblk(1,is).eq.isw2.and. 9d7s16
     $     isblk(3,is).eq.isw3.and.isblk(4,is).eq.isblkder1 (4,isb))then
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d8s16
              iad2=itrans(isw3)+i1m+nbasdwsc(isw3)*m3                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d8s16
               iad1=itrans(isw2)+i3+nbasdwsc(isw2)*m2                   9d8s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (2,isb))*(m2            9d8s16
     $        +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i2m))9d7s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 512                                                      9d7s16
        end if                                                          9d2s16
        if(isblk(2,is).eq.isblkder1 (1,isb).and.isblk(1,is).eq.isw2.and. 9d7s16
     $     isblk(4,is).eq.isw3.and.isblk(3,is).eq.isblkder1 (4,isb))then 9d8s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d8s16
              iad2=itrans(isw3)+i2m+nbasdwsc(isw3)*m3                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d8s16
               iad1=itrans(isw2)+i3+nbasdwsc(isw2)*m2                   9d8s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (2,isb))*(m2            9d8s16
     $        +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i1m))9d7s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 512                                                      9d7s16
        end if                                                          9d2s16
        if(isblk(1,is).eq.isblkder1 (1,isb).and.isblk(2,is).eq.isw2.and. 9d7s16
     $     isblk(4,is).eq.isw3.and.isblk(3,is).eq.isblkder1 (4,isb))then 9d8s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d8s16
              iad2=itrans(isw3)+i2m+nbasdwsc(isw3)*m3                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d8s16
               iad1=itrans(isw2)+i4+nbasdwsc(isw2)*m2                   9d8s16
               iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $        +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i1m))9d7s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 512                                                      9d7s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       write(6,*)('no 4o for (oo''|o''o)'),(isblkder1 (j,isb),j=1,4)     9d7s16
       call dws_sync                                                    9d2s16
       call dws_finalize                                                9d2s16
       stop                                                             9d2s16
  512  continue                                                         9d7s16
c
c     (oo'|o'o): (ov|oo) part
c
       do is=1,nsdlk1                                                   9d2s16
        if(isblk1(4,is).eq.isw2.and.isblk1(3,is).eq.isblkder1 (1,isb)    9d8s16
     $       .and.isblk1(1,is).eq.isw3.and.                             9d8s16
     $     isblk1(2,is).eq.isblkder1 (4,isb))then                        9d8s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           if(isblk1(1,is).eq.isblk1(2,is))then                           9d2s16
            do i4=0,nocc(isblk1(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d8s16
               iad2=itrans(isw3)+i3+nbasdwsc(isw3)*m3                   9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               jad2=itrans(isw3)+i4+nbasdwsc(isw3)*m3                   9d8s16
               fb=2d0*bc(ii)*bc(jad2)                                       9d2s16
               do m2=0,nocc(isblkder1 (2,isb))-1                         9d8s16
                iad1=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                 9d8s16
                iad3=i4od2b(isb)+i1m+nocc(isblkder1 (1,isb))*(m2          9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i4))9d8s16
                jad3=i4od2b(isb)+i1m+nocc(isblkder1 (1,isb))*(m2          9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i3))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
                bc(jad3)=bc(jad3)+bc(iad1)*fb                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d8s16
              iad2=itrans(isw3)+i4+nbasdwsc(isw3)*m3                    9d8s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d8s16
               iad1=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                  9d8s16
               iad3=i4od2b(isb)+i1m+nocc(isblkder1 (1,isb))*(m2           9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk1(2,is))-1                                9d2s16
             do i3=0,nocc(isblk1(1,is))-1                               9d2s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d8s16
               iad2=itrans(isw3)+i3+nbasdwsc(isw3)*m3                   9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m2=0,nocc(isblkder1 (2,isb))-1                         9d8s16
                iad1=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                 9d8s16
                iad3=i4od2b(isb)+i1m+nocc(isblkder1 (1,isb))*(m2          9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i4))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 513                                                      9d8s16
        end if                                                          9d2s16
        if(isblk1(4,is).eq.isw2.and.isblk1(3,is).eq.isblkder1 (1,isb)    9d8s16
     $       .and.isblk1(2,is).eq.isw3.and.                             9d8s16
     $     isblk1(1,is).eq.isblkder1 (4,isb))then                        9d8s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk1(2,is))-1                                 9d8s16
            do i3=0,nocc(isblk1(1,is))-1                                9d8s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d8s16
              iad2=itrans(isw3)+i4+nbasdwsc(isw3)*m3                    9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d8s16
               iad1=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                  9d8s16
               iad3=i4od2b(isb)+i1m+nocc(isblkder1 (1,isb))*(m2           9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i3))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 513                                                      9d8s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       if(min(nbasdwsc(isw2),nocc(isw3)).gt.0)then                      11d28s22
        write(6,*)('no onex (ov|oo) for (oo''|o''o) '),                  9d8s16
     $      (isblkder1 (j,isb),j=1,4),isw3,isw2                                    9d2s16
        call dws_sync                                                    9d2s16
        call dws_finalize                                                9d2s16
        stop                                                             9d2s16
       end if                                                           11d28s22
  513  continue                                                         9d8s16
c
c     (oo'|o'o): (oo|vo) part
c
       do is=1,nsdlk1                                                   9d2s16
        if(isblk1(4,is).eq.isw3.and.isblk1(3,is).eq.isblkder1 (4,isb)    9d8s16
     $       .and.isblk1(1,is).eq.isblkder1 (1,isb).and.                 9d8s16
     $     isblk1(2,is).eq.isw2)then                                    9d8s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1
           if(isblk1(1,is).eq.isblk1(2,is))then                           9d2s16
            do i4=0,nocc(isblk1(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d8s16
               iad2=itrans(isw3)+i2p+nbasdwsc(isw3)*m3                  9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m2=0,nocc(isblkder1 (2,isb))-1                         9d8s16
                iad1=itrans(isw2)+i3+nbasdwsc(isw2)*m2                  9d8s16
                iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(m2           9d8s16
     $        +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i1m))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
                iad1=itrans(isw2)+i4+nbasdwsc(isw2)*m2                  9d8s16
                iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(m2           9d8s16
     $        +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i1m))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d8s16
              iad2=itrans(isw3)+i2p+nbasdwsc(isw3)*m3                   9d8s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d8s16
               iad1=itrans(isw2)+i4+nbasdwsc(isw2)*m2                   9d8s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $        +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i1m))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk1(2,is))-1                                9d2s16
             do i3=0,nocc(isblk1(1,is))-1                               9d2s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d8s16
               iad2=itrans(isw3)+i2p+nbasdwsc(isw3)*m3                  9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m2=0,nocc(isblkder1 (2,isb))-1                         9d8s16
                iad1=itrans(isw2)+i4+nbasdwsc(isw2)*m2                  9d8s16
                iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(m2           9d8s16
     $        +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i1m))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 514                                                      9d8s16
        end if                                                          9d2s16
        if(isblk1(4,is).eq.isw3.and.isblk1(3,is).eq.isblkder1 (4,isb)    9d8s16
     $       .and.isblk1(2,is).eq.isblkder1 (1,isb).and.                 9d8s16
     $     isblk1(1,is).eq.isw2)then                                    9d8s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1
           do i4=0,nocc(isblk1(2,is))-1                                 9d8s16
            do i3=0,nocc(isblk1(1,is))-1                                9d8s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d8s16
              iad2=itrans(isw3)+i2p+nbasdwsc(isw3)*m3                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d8s16
               iad1=itrans(isw2)+i3+nbasdwsc(isw2)*m2                   9d8s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $        +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i1m))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 514                                                      9d8s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       if(min(nbasdwsc(isw3),nocc(isw2)).gt.0)then                      11d28s22
        write(6,*)('no onex (oo|vo)  for (oo''|o''o) '),                 9d8s16
     $      (isblkder1 (j,isb),j=1,4),isw2,isw3                                    9d2s16
        call dws_sync                                                    9d2s16
        call dws_finalize                                                9d2s16
        stop                                                             9d2s16
       end if                                                           11d28s22
  514  continue                                                         9d8s16
c
c     (oo'|o'o): (ov|vo) part
c     recall k_nm^ab=(nb|ma)
c
       if(min(nbasdwsc(isw2),nbasdwsc(isw3)).eq.0)go to 515             11d28s22
       do is=1,nsdlkk                                                   9d2s16
        if(isblkk(4,is).eq.isw2.and.isblkk(3,is).eq.isw3                9d8s16
     $       .and.isblkk(1,is).eq.isblkder1 (1,isb).and.                 9d8s16
     $     isblkk(2,is).eq.isblkder1 (4,isb))then                        9d8s16
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d2s16
         ii=kmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblkk(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblkk(3,is))                                   9d2s16
           do i4=0,nocc(isblkk(2,is))-1                                  9d2s16
            do i3=0,nocc(isblkk(1,is))-1                                 9d2s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d8s16
              iad2=itrans(isw3)+i1p+nbasdwsc(isw3)*m3                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d2s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d8s16
               iad1=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                  9d8s16
               iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i4))9d9s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           end do                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 515                                                      9d7s16
        end if                                                          9d2s16
        if(isblkk(3,is).eq.isw2.and.isblkk(4,is).eq.isw3                9d8s16
     $       .and.isblkk(2,is).eq.isblkder1 (1,isb).and.                 9d8s16
     $     isblkk(1,is).eq.isblkder1 (4,isb))then                        9d8s16
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d2s16
         ii=kmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblkk(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblkk(3,is))                                   9d2s16
           do i4=0,nocc(isblkk(2,is))-1                                  9d2s16
            do i3=0,nocc(isblkk(1,is))-1                                 9d2s16
             do m3=0,nocc(isblkder1 (3,isb))-1                           9d8s16
              iad2=itrans(isw3)+i2p+nbasdwsc(isw3)*m3                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d2s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d8s16
               iad1=itrans(isw2)+i1p+nbasdwsc(isw2)*m2                  9d8s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*i3))9d9s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           end do                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 515                                                      9d7s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       write(6,*)('no kmats  for (oo''|o''o) '),                        9d8s16
     $      (isblkder1 (j,isb),j=1,4)                                    9d2s16
       call dws_sync                                                    9d2s16
       call dws_finalize                                                9d2s16
       stop                                                             9d2s16
  515  continue                                                         9d8s16
c
c     (oo'|oo'), (oo|oo) part                                           9d8s16
c
       isw2=multh(isblkder1 (2,isb),iprop)                               9d7s16
       isw4=multh(isblkder1 (4,isb),iprop)                               9d8s16
       do is=1,nsdlk                                                    9d2s16
        if(isblk(1,is).eq.isblkder1 (1,isb).and.isblk(2,is).eq.isw2.and. 9d7s16
     $     isblk(3,is).eq.isblkder1 (3,isb).and.isblk(4,is).eq.isw4)then 9d8s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           if(isblk(1,is).eq.isblk(2,is))then                           9d2s16
            do i4=0,nocc(isblk(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d8s16
               iad2=itrans(isw4)+i2m+nbasdwsc(isw4)*m4                  9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m2=0,nocc(isblkder1 (2,isb))-1                         9d7s16
                iad1=itrans(isw2)+i4+nbasdwsc(isw2)*m2                  9d7s16
                iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(m2           9d7s16
     $        +nocc(isblkder1 (2,isb))*(i1m+nocc(isblkder1 (3,isb))*m4))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
                iad1=itrans(isw2)+i3+nbasdwsc(isw2)*m2                  9d7s16
                iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(m2           9d7s16
     $        +nocc(isblkder1 (2,isb))*(i1m+nocc(isblkder1 (3,isb))*m4))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i2m+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d7s16
               iad1=itrans(isw2)+i4+nbasdwsc(isw2)*m2                   9d7s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(m2            9d7s16
     $        +nocc(isblkder1 (2,isb))*(i1m+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk(2,is))-1                                 9d2s16
             do i3=0,nocc(isblk(1,is))-1                                9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d8s16
               iad2=itrans(isw4)+i2m+nbasdwsc(isw4)*m4                  9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m2=0,nocc(isblkder1 (2,isb))-1                         9d7s16
                iad1=itrans(isw2)+i4+nbasdwsc(isw2)*m2                  9d7s16
                iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(m2           9d7s16
     $        +nocc(isblkder1 (2,isb))*(i1m+nocc(isblkder1 (3,isb))*m4))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
             ii=ii+1                                                    9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 516                                                      9d8s16
        end if                                                          9d2s16
        if(isblk(2,is).eq.isblkder1 (1,isb).and.isblk(1,is).eq.isw2.and. 9d8s16
     $     isblk(3,is).eq.isblkder1 (3,isb).and.isblk(4,is).eq.isw4)then 9d8s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i2m+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d8s16
               iad1=itrans(isw2)+i3+nbasdwsc(isw2)*m2                   9d8s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $        +nocc(isblkder1 (2,isb))*(i1m+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 516                                                      9d8s16
        end if                                                          9d2s16
        if(isblk(2,is).eq.isblkder1 (1,isb).and.isblk(1,is).eq.isw2.and. 9d8s16
     $     isblk(4,is).eq.isblkder1 (3,isb).and.isblk(3,is).eq.isw4)then 9d8s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i1m+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d8s16
               iad1=itrans(isw2)+i3+nbasdwsc(isw2)*m2                   9d8s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $        +nocc(isblkder1 (2,isb))*(i2m+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 516                                                      9d8s16
        end if                                                          9d2s16
        if(isblk(1,is).eq.isblkder1 (1,isb).and.isblk(2,is).eq.isw2.and. 9d8s16
     $     isblk(4,is).eq.isblkder1 (3,isb).and.isblk(3,is).eq.isw4)then 9d8s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i1m+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d8s16
               iad1=itrans(isw2)+i4+nbasdwsc(isw2)*m2                   9d8s16
               iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $        +nocc(isblkder1 (2,isb))*(i2m+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 516                                                      9d8s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       write(6,*)('no 4o for (oo''|oo'')'),(isblkder1 (j,isb),j=1,4)     9d8s16
       call dws_sync                                                    9d2s16
       call dws_finalize                                                9d2s16
       stop                                                             9d2s16
  516  continue                                                         9d8s16
c
c     (oo'|oo'): (ov|oo) part
c
       do is=1,nsdlk1                                                   9d2s16
        if(isblk1(4,is).eq.isw2.and.isblk1(3,is).eq.isblkder1 (1,isb)    9d8s16
     $       .and.isblk1(1,is).eq.isblkder1 (3,isb).and.                 9d8s16
     $     isblk1(2,is).eq.isw4)then                                    9d8s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           if(isblk1(1,is).eq.isblk1(2,is))then                           9d2s16
            do i4=0,nocc(isblk1(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d8s16
               iad2=itrans(isw4)+i3+nbasdwsc(isw4)*m4                   9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               jad2=itrans(isw4)+i4+nbasdwsc(isw4)*m4                   9d8s16
               fb=2d0*bc(ii)*bc(jad2)                                       9d2s16
               do m2=0,nocc(isblkder1 (2,isb))-1                         9d8s16
                iad1=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                 9d8s16
                iad3=i4od2b(isb)+i1m+nocc(isblkder1 (1,isb))*(m2          9d8s16
     $         +nocc(isblkder1 (2,isb))*(i4+nocc(isblkder1 (3,isb))*m4))9d8s16
                jad3=i4od2b(isb)+i1m+nocc(isblkder1 (1,isb))*(m2          9d8s16
     $         +nocc(isblkder1 (2,isb))*(i3+nocc(isblkder1 (3,isb))*m4))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
                bc(jad3)=bc(jad3)+bc(iad1)*fb                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i4+nbasdwsc(isw4)*m4                    9d8s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d8s16
               iad1=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                  9d8s16
               iad3=i4od2b(isb)+i1m+nocc(isblkder1 (1,isb))*(m2           9d8s16
     $         +nocc(isblkder1 (2,isb))*(i4+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk1(2,is))-1                                9d2s16
             do i3=0,nocc(isblk1(1,is))-1                               9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d8s16
               iad2=itrans(isw4)+i4+nbasdwsc(isw4)*m4                   9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m2=0,nocc(isblkder1 (2,isb))-1                         9d8s16
                iad1=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                 9d8s16
                iad3=i4od2b(isb)+i1m+nocc(isblkder1 (1,isb))*(m2          9d8s16
     $         +nocc(isblkder1 (2,isb))*(i3+nocc(isblkder1 (3,isb))*m4))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 517                                                      9d8s16
        end if                                                          9d2s16
        if(isblk1(4,is).eq.isw2.and.isblk1(3,is).eq.isblkder1 (1,isb)    9d8s16
     $       .and.isblk1(2,is).eq.isblkder1 (3,isb).and.                 9d8s16
     $     isblk1(1,is).eq.isw4)then                                    9d8s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk1(2,is))-1                                 9d8s16
            do i3=0,nocc(isblk1(1,is))-1                                9d8s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i3+nbasdwsc(isw4)*m4                    9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d8s16
               iad1=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                  9d8s16
               iad3=i4od2b(isb)+i1m+nocc(isblkder1 (1,isb))*(m2           9d8s16
     $         +nocc(isblkder1 (2,isb))*(i4+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 517                                                      9d8s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       if(min(nocc(isw4),nbasdwsc(isw2)).gt.0)then                      11d28s22
        write(6,*)('no onex (ov|oo) for (oo''|oo'') '),                  9d8s16
     $      (isblkder1 (j,isb),j=1,4),isw4,isw2                                    9d2s16
        call dws_sync                                                    9d2s16
        call dws_finalize                                                9d2s16
        stop                                                             9d2s16
       end if                                                           11d28s22
  517  continue                                                         9d8s16
c
c     (oo'|oo'): (oo|ov) part
c
       do is=1,nsdlk1                                                   9d2s16
        if(isblk1(4,is).eq.isw4.and.isblk1(3,is).eq.isblkder1 (3,isb)    9d8s16
     $       .and.isblk1(1,is).eq.isblkder1 (1,isb).and.                 9d8s16
     $     isblk1(2,is).eq.isw2)then                                    9d8s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1
           if(isblk1(1,is).eq.isblk1(2,is))then                           9d2s16
            do i4=0,nocc(isblk1(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d8s16
               iad2=itrans(isw4)+i2p+nbasdwsc(isw4)*m4                  9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m2=0,nocc(isblkder1 (2,isb))-1                         9d8s16
                iad1=itrans(isw2)+i3+nbasdwsc(isw2)*m2                  9d8s16
                iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(m2           9d8s16
     $        +nocc(isblkder1 (2,isb))*(i1m+nocc(isblkder1 (3,isb))*m4))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
                iad1=itrans(isw2)+i4+nbasdwsc(isw2)*m2                  9d8s16
                iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(m2           9d8s16
     $        +nocc(isblkder1 (2,isb))*(i1m+nocc(isblkder1 (3,isb))*m4))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i2p+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d8s16
               iad1=itrans(isw2)+i4+nbasdwsc(isw2)*m2                   9d8s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $        +nocc(isblkder1 (2,isb))*(i1m+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk1(2,is))-1                                9d2s16
             do i3=0,nocc(isblk1(1,is))-1                               9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d8s16
               iad2=itrans(isw4)+i2p+nbasdwsc(isw4)*m4                  9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m2=0,nocc(isblkder1 (2,isb))-1                         9d8s16
                iad1=itrans(isw2)+i4+nbasdwsc(isw2)*m2                  9d8s16
                iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(m2           9d8s16
     $        +nocc(isblkder1 (2,isb))*(i1m+nocc(isblkder1 (3,isb))*m4))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 518                                                      9d8s16
        end if                                                          9d2s16
        if(isblk1(4,is).eq.isw4.and.isblk1(3,is).eq.isblkder1 (3,isb)    9d8s16
     $       .and.isblk1(2,is).eq.isblkder1 (1,isb).and.                 9d8s16
     $     isblk1(1,is).eq.isw2)then                                    9d8s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1
           do i4=0,nocc(isblk1(2,is))-1                                 9d8s16
            do i3=0,nocc(isblk1(1,is))-1                                9d8s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i2p+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d8s16
               iad1=itrans(isw2)+i3+nbasdwsc(isw2)*m2                   9d8s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $        +nocc(isblkder1 (2,isb))*(i1m+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 518                                                      9d8s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       if(min(nocc(isw2),nbasdwsc(isw4)).gt.0)then                      11d28s22
        write(6,*)('no onex (oo|ov)  for (oo''|oo'') '),                 9d8s16
     $      (isblkder1 (j,isb),j=1,4),isw2,isw4                                    9d2s16
        call dws_sync                                                    9d2s16
        call dws_finalize                                                9d2s16
        stop                                                             9d2s16
       end if                                                           11d28s22
  518  continue                                                         9d8s16
c
c     (oo'|oo'): (ov|ov) part
c     recall k_nm^ab=(nb|ma)
c
       do is=1,nsdlkk                                                   9d2s16
        if(isblkk(4,is).eq.isw2.and.isblkk(3,is).eq.isw4                9d8s16
     $       .and.isblkk(1,is).eq.isblkder1 (1,isb).and.                 9d8s16
     $     isblkk(2,is).eq.isblkder1 (3,isb))then                        9d8s16
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d2s16
         ii=kmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblkk(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblkk(3,is))                                   9d2s16
           do i4=0,nocc(isblkk(2,is))-1                                  9d2s16
            do i3=0,nocc(isblkk(1,is))-1                                 9d2s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i1p+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d2s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d8s16
               iad1=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                  9d8s16
               iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $             +nocc(isblkder1 (2,isb))*(i4+nocc(isblkk(2,is))*m4)) 9d9s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           end do                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 519                                                      9d8s16
        end if                                                          9d2s16
        if(isblkk(3,is).eq.isw2.and.isblkk(4,is).eq.isw4                9d8s16
     $       .and.isblkk(2,is).eq.isblkder1 (1,isb).and.                 9d8s16
     $     isblkk(1,is).eq.isblkder1 (3,isb))then                        9d8s16
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d2s16
         ii=kmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblkk(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblkk(3,is))                                   9d2s16
           do i4=0,nocc(isblkk(2,is))-1                                  9d2s16
            do i3=0,nocc(isblkk(1,is))-1                                 9d2s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i2p+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d2s16
              do m2=0,nocc(isblkder1 (2,isb))-1                          9d8s16
               iad1=itrans(isw2)+i1p+nbasdwsc(isw2)*m2                  9d8s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(m2            9d8s16
     $             +nocc(isblkder1 (2,isb))*(i3+nocc(isblkk(1,is))*m4)) 9d9s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           end do                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 519                                                      9d8s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       if(min(nbasdwsc(isw2),nbasdwsc(isw4)).gt.0)then                  11d28s22
        write(6,*)('no kmats  for (oo''|oo'') '),                        9d8s16
     $      (isblkder1 (j,isb),j=1,4)                                    9d2s16
        call dws_sync                                                    9d2s16
        call dws_finalize                                                9d2s16
        stop                                                             9d2s16
       end if                                                           11d28s22
  519  continue                                                         9d8s16
c
c     (oo|o'o'): (oo|oo) part
c
       isw3=multh(isblkder1 (3,isb),iprop)                               9d8s16
       isw4=multh(isblkder1 (4,isb),iprop)                               9d8s16
       if(nocc(isw3)*nocc(isw4).ne.0)then                               11d30s16
       do is=1,nsdlk                                                    9d2s16
        if(isblk(1,is).eq.isblkder1 (1,isb).and.                         9d8s16
     $       isblk(2,is).eq.isblkder1 (2,isb).and.                       9d8s16
     $     isblk(3,is).eq.isw3.and.                                     9d8s16
     $     isblk(4,is).eq.isw4)then                                     9d8s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           if(isblk(1,is).eq.isblk(2,is))then                           9d2s16
            do i4=0,nocc(isblk(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d2s16
               iad2=itrans(isw4)+i2m+nbasdwsc(isw4)*m4                  9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m3=0,nocc(isblkder1 (3,isb))-1                         9d8s16
                iad1=itrans(isw3)+i1m+nbasdwsc(isw3)*m3                 9d8s16
                iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(i4           9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
                iad4=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(i3           9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
                term=bc(iad1)*fa                                        9d8s16
                bc(iad3)=bc(iad3)+term                                  9d8s16
                bc(iad4)=bc(iad4)+term                                  9d8s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i2m+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d8s16
               iad1=itrans(isw3)+i1m+nbasdwsc(isw3)*m3                  9d8s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(i4            9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk(2,is))-1                                 9d2s16
             do i3=0,nocc(isblk(1,is))-1                                9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d8s16
               iad2=itrans(isw4)+i2m+nbasdwsc(isw4)*m4                  9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m3=0,nocc(isblkder1 (3,isb))-1                         9d8s16
                iad1=itrans(isw3)+i1m+nbasdwsc(isw3)*m3                 9d8s16
                iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(i4           9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
             ii=ii+1                                                    9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 520                                                      9d8s16
        end if                                                          9d2s16
        if(isblk(2,is).eq.isblkder1 (1,isb).and.                         9d8s16
     $       isblk(1,is).eq.isblkder1 (2,isb).and.                       9d8s16
     $     isblk(3,is).eq.isw3.and.                                     9d8s16
     $     isblk(4,is).eq.isw4)then                                     9d8s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i2m+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d8s16
               iad1=itrans(isw3)+i1m+nbasdwsc(isw3)*m3                  9d8s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(i3            9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 520                                                      9d8s16
        end if                                                          9d2s16
        if(isblk(2,is).eq.isblkder1 (1,isb).and.                         9d8s16
     $       isblk(1,is).eq.isblkder1 (2,isb).and.                       9d8s16
     $     isblk(4,is).eq.isw3.and.                                     9d8s16
     $     isblk(3,is).eq.isw4)then                                     9d8s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i1m+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d8s16
               iad1=itrans(isw3)+i2m+nbasdwsc(isw3)*m3                  9d8s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(i3            9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 520                                                      9d8s16
        end if                                                          9d2s16
        if(isblk(1,is).eq.isblkder1 (1,isb).and.                         9d8s16
     $       isblk(2,is).eq.isblkder1 (2,isb).and.                       9d8s16
     $     isblk(4,is).eq.isw3.and.                                     9d8s16
     $     isblk(3,is).eq.isw4)then                                     9d8s16
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d2s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i1m+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d8s16
               iad1=itrans(isw3)+i2m+nbasdwsc(isw3)*m3                  9d8s16
               iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(i4            9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 520                                                      9d8s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       write(6,*)('no 4o for (oo|o''o'')'),(isblkder1 (j,isb),j=1,4)     9d8s16
       call dws_sync                                                    9d2s16
       call dws_finalize                                                9d2s16
       stop                                                             9d2s16
  520  continue                                                         9d8s16
       end if
c
c     (oo|o'o'): (oo|vo) part
c
       if(min(nocc(isw4),nbasdwsc(isw3)).gt.0)then                      11d28s22
       do is=1,nsdlk1                                                   9d2s16
        if(isblk1(4,is).eq.isw3.and.isblk1(3,is).eq.isw4.and.           9d8s16
     $     isblk1(1,is).eq.isblkder1 (1,isb).and.                        9d8s16
     $     isblk1(2,is).eq.isblkder1 (2,isb))then                        9d8s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           if(isblk1(1,is).eq.isblk1(2,is))then                           9d2s16
            do i4=0,nocc(isblk1(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d8s16
               iad2=itrans(isw4)+i1m+nbasdwsc(isw4)*m4                  9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m3=0,nocc(isblkder1 (3,isb))-1                         9d8s16
                iad1=itrans(isw3)+i2p+nbasdwsc(isw3)*m3                 9d8s16
                iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(i4           9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
                iad4=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(i3           9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
                term=bc(iad1)*fa                                        9d8s16
                bc(iad3)=bc(iad3)+term                                  9d8s16
                bc(iad4)=bc(iad4)+term                                  9d8s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i1m+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d8s16
               iad1=itrans(isw3)+i2p+nbasdwsc(isw3)*m3                  9d2s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(i4            9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk1(2,is))-1                                9d2s16
             do i3=0,nocc(isblk1(1,is))-1                               9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d8s16
               iad2=itrans(isw4)+i1m+nbasdwsc(isw4)*m4                  9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m3=0,nocc(isblkder1 (3,isb))-1                         9d8s16
                iad1=itrans(isw3)+i2p+nbasdwsc(isw3)*m3                 9d8s16
                iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(i4           9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 521                                                      9d8s16
        end if                                                          9d2s16
        if(isblk1(4,is).eq.isw3.and.isblk1(3,is).eq.isw4.and.           9d8s16
     $     isblk1(2,is).eq.isblkder1 (1,isb).and.                        9d8s16
     $     isblk1(1,is).eq.isblkder1 (2,isb))then                        9d8s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d2s16
           do i4=0,nocc(isblk1(2,is))-1                                 9d8s16
            do i3=0,nocc(isblk1(1,is))-1                                9d8s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i1m+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d8s16
               iad1=itrans(isw3)+i2p+nbasdwsc(isw3)*m3                  9d8s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(i3            9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 521                                                      9d8s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       write(6,*)('no onex (oo|vo) for (oo|o''o'') '),                  9d8s16
     $      (isblkder1 (j,isb),j=1,4)                                    9d2s16
       call dws_sync                                                    9d2s16
       call dws_finalize                                                9d2s16
       stop                                                             9d2s16
  521  continue                                                         9d8s16
       end if                                                           11d30s16
c
c     (oo|o'o'): (oo|ov) part
c
       if(min(nocc(isw3),nbasdwsc(isw4)).gt.0)then                      11d28s22
       do is=1,nsdlk1                                                   9d2s16
        if(isblk1(4,is).eq.isw4.and.isblk1(3,is).eq.isw3.and.           9d8s16
     $     isblk1(1,is).eq.isblkder1 (1,isb).and.                        9d8s16
     $     isblk1(2,is).eq.isblkder1 (2,isb))then                        9d8s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d8s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d8s16
           if(isblk1(1,is).eq.isblk1(2,is))then                           9d2s16
            do i4=0,nocc(isblk1(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d8s16
               iad2=itrans(isw4)+i2p+nbasdwsc(isw4)*m4                  9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m3=0,nocc(isblkder1 (3,isb))-1                         9d8s16
                iad1=itrans(isw3)+i1m+nbasdwsc(isw3)*m3                 9d8s16
                iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(i4           9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
                jad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(i3           9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
                bc(jad3)=bc(jad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i2p+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d8s16
               iad1=itrans(isw3)+i1m+nbasdwsc(isw3)*m3                  9d8s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(i4            9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk1(2,is))-1                                9d2s16
             do i3=0,nocc(isblk1(1,is))-1                               9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d8s16
               iad2=itrans(isw4)+i2p+nbasdwsc(isw4)*m4                  9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m3=0,nocc(isblkder1 (3,isb))-1                         9d8s16
                iad1=itrans(isw3)+i1m+nbasdwsc(isw3)*m3                 9d8s16
                iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(i4           9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 522                                                      9d2s16
        end if                                                          9d2s16
        if(isblk1(4,is).eq.isw4.and.isblk1(3,is).eq.isw3.and.           9d8s16
     $     isblk1(2,is).eq.isblkder1 (1,isb).and.                        9d8s16
     $     isblk1(1,is).eq.isblkder1 (2,isb))then                        9d8s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d2s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d8s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d8s16
           do i4=0,nocc(isblk1(2,is))-1                                 9d8s16
            do i3=0,nocc(isblk1(1,is))-1                                9d8s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i2p+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d8s16
               iad1=itrans(isw3)+i1m+nbasdwsc(isw3)*m3                  9d8s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(i3            9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 522                                                      9d2s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       write(6,*)('no onex (oo|ov)  for (oo|o''o'') '),                 9d8s16
     $      (isblkder1 (j,isb),j=1,4)                                    9d2s16
       call dws_sync                                                    9d2s16
       call dws_finalize                                                9d2s16
       stop                                                             9d2s16
  522  continue                                                         9d2s16
       end if                                                           11d30s16
c
c     (oo|o'o'): (oo|vv) part
c
       if(min(nbasdwsc(isw4),nbasdwsc(isw3)).eq.0)go to 523             11d28s22
       do is=1,nsdlk                                                    9d2s16
        if(isblk(4,is).eq.isw4.and.isblk(3,is).eq.isw3.and.             9d8s16
     $     isblk(1,is).eq.isblkder1 (1,isb).and.                         9d8s16
     $     isblk(2,is).eq.isblkder1 (2,isb))then                         9d8s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d2s16
         ii=jmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk(4,is))                                    9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblk(3,is))                                   9d2s16
           if(isblk(1,is).eq.isblk(2,is))then                           9d2s16
            do i4=0,nocc(isblk(1,is))-1                                 9d2s16
             do i3=0,i4-1                                               9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d8s16
               iad2=itrans(isw4)+i2p+nbasdwsc(isw4)*m4                  9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m3=0,nocc(isblkder1 (3,isb))-1                         9d8s16
                iad1=itrans(isw3)+i1p+nbasdwsc(isw3)*m3                 9d8s16
                iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(i3           9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
                jad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(i4           9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
                bc(jad3)=bc(jad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i2p+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(iad2)*bc(ii)                                        9d2s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d2s16
               iad1=itrans(isw3)+i1p+nbasdwsc(isw3)*m3                  9d2s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(i4            9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d2s16
              end do                                                    9d2s16
             end do                                                     9d2s16
             ii=ii+1                                                    9d2s16
            end do                                                      9d2s16
           else                                                         9d2s16
            do i4=0,nocc(isblk(2,is))-1                                 9d2s16
             do i3=0,nocc(isblk(1,is))-1                                9d2s16
              do m4=0,nocc(isblkder1 (4,isb))-1                          9d8s16
               iad2=itrans(isw4)+i2p+nbasdwsc(isw4)*m4                  9d8s16
               fa=2d0*bc(ii)*bc(iad2)                                       9d2s16
               do m3=0,nocc(isblkder1 (3,isb))-1                         9d8s16
                iad1=itrans(isw3)+i1p+nbasdwsc(isw3)*m3                 9d8s16
                iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(i4           9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
                bc(iad3)=bc(iad3)+bc(iad1)*fa                           9d2s16
               end do                                                   9d2s16
              end do                                                    9d2s16
              ii=ii+1                                                   9d2s16
             end do                                                     9d2s16
            end do                                                      9d2s16
           end if                                                       9d2s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 523                                                      9d8s16
        end if                                                          9d2s16
        if(isblk(4,is).eq.isw4.and.isblk(3,is).eq.isw3.and.             9d8s16
     $     isblk(2,is).eq.isblkder1 (1,isb).and.                         9d8s16
     $     isblk(1,is).eq.isblkder1 (2,isb))then                         9d8s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d2s16
         ii=jmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk(4,is))                                    9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblk(3,is))                                   9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i2p+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d8s16
               iad1=itrans(isw3)+i1p+nbasdwsc(isw3)*m3                  9d8s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(i3            9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 523                                                      9d8s16
        end if                                                          9d2s16
        if(isblk(3,is).eq.isw4.and.isblk(4,is).eq.isw3.and.             9d8s16
     $     isblk(2,is).eq.isblkder1 (1,isb).and.                         9d8s16
     $     isblk(1,is).eq.isblkder1 (2,isb))then                         9d8s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d2s16
         ii=jmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk(4,is))                                    9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblk(3,is))                                   9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i1p+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d8s16
               iad1=itrans(isw3)+i2p+nbasdwsc(isw3)*m3                  9d8s16
               iad3=i4od2b(isb)+i4+nocc(isblkder1 (1,isb))*(i3            9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 523                                                      9d8s16
        end if                                                          9d2s16
        if(isblk(3,is).eq.isw4.and.isblk(4,is).eq.isw3.and.             9d8s16
     $     isblk(1,is).eq.isblkder1 (1,isb).and.                         9d8s16
     $     isblk(2,is).eq.isblkder1 (2,isb))then                         9d8s16
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d2s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d2s16
         ii=jmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk(4,is))                                    9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblk(3,is))                                   9d2s16
           do i4=0,nocc(isblk(2,is))-1                                  9d8s16
            do i3=0,nocc(isblk(1,is))-1                                 9d8s16
             do m4=0,nocc(isblkder1 (4,isb))-1                           9d8s16
              iad2=itrans(isw4)+i1p+nbasdwsc(isw4)*m4                   9d8s16
              fa=2d0*bc(ii)*bc(iad2)                                        9d8s16
              do m3=0,nocc(isblkder1 (3,isb))-1                          9d8s16
               iad1=itrans(isw3)+i2p+nbasdwsc(isw3)*m3                  9d8s16
               iad3=i4od2b(isb)+i3+nocc(isblkder1 (1,isb))*(i4            9d8s16
     $         +nocc(isblkder1 (2,isb))*(m3+nocc(isblkder1 (3,isb))*m4))9d8s16
               bc(iad3)=bc(iad3)+bc(iad1)*fa                            9d8s16
              end do                                                    9d8s16
             end do                                                     9d8s16
             ii=ii+1                                                    9d8s16
            end do                                                      9d8s16
           end do                                                       9d8s16
          end do                                                        9d2s16
          i10=1                                                         9d2s16
         end do                                                         9d2s16
         go to 523                                                      9d8s16
        end if                                                          9d2s16
       end do                                                           9d2s16
       write(6,*)('no jmats (oo|vv)  for (oo|o''o'') '),                9d8s16
     $      (isblkder1 (j,isb),j=1,4)                                    9d2s16
       call dws_sync                                                    9d2s16
       call dws_finalize                                                9d2s16
       stop                                                             9d2s16
  523  continue                                                         9d8s16
       if(idwsdeb.gt.10)then
        write(6,*)('mixed 4od2 after xor parts')
        write(6,12)(isblkder1(i,isb),i=1,4)
        write(6,*)('not yet global summed ')
        call prntm2(bc(i4od2b(isb)),nrow,ncol,nrow)
        itmp=ibcoff
        ibcoff=itmp+nrow*ncol
        do i=0,nrow*ncol-1
         bc(itmp+i)=bc(i4od2b(isb)+i)                                   10d31s16
        end do
        call dws_gsumf(bc(itmp),nrow*ncol)
        write(6,*)('gsummed: ')
        if(nsymb.eq.1)then
         call printa(bc(itmp),nocc(isblkder1 (1,isb)),0,
     $        nocc(isblkder1 (2,isb)),0,nocc(isblkder1 (3,isb)),0,
     $        nocc(isblkder1 (4,isb)),0,bc(ibcoff))
        else
         call prntm2(bc(itmp),nrow,ncol,nrow)
        end if
        ibcoff=itmp
       end if
      end do                                                            9d2s16
      do isb=1,nsblkxder1
       nrow=nocc(isblkxder1(1,isb))*nocc(isblkxder1(2,isb))
       ncol=nocc(isblkxder1(3,isb))*nvirtc(isblkxder1(4,isb))
       if(idwsdeb.gt.10.and.min(nrow,ncol).gt.0)then
       write(6,*)('starting onexd2 ooox'),ionexd2(isb)
       write(6,12)(isblkxder1(j,isb),j=1,4)
       call prntm2(bc(ionexd2(isb)),nrow,ncol,nrow)
       itmp=ibcoff
       ibcoff=itmp+nrow*ncol
       call enough('parajkfromhd0.  9',bc,ibc)
       do i=0,nrow*ncol-1
        bc(itmp+i)=bc(ionexd2(isb)+i)
       end do
       call dws_gsumf(bc(itmp),nrow*ncol)
       write(6,*)('after global sum: ')
       call prntm2(bc(itmp),nrow,ncol,nrow)
        if(nsymb.eq.1)then
         call printa(bc(itmp),nocc(isblkxder1 (1,isb)),0,
     $        nocc(isblkxder1 (2,isb)),0,nocc(isblkxder1 (3,isb)),0,
     $        nvirt(isblkxder1(4,isb)),nocc(isblkxder1 (4,isb)),
     $        bc(ibcoff))
        end if
       ibcoff=itmp
       end if
       isw1=multh(isblkxder1(1,isb),iprop)
       isw2=multh(isblkxder1(2,isb),iprop)
       isw3=multh(isblkxder1(3,isb),iprop)
       isw4=multh(isblkxder1(4,isb),iprop)
c
c     (o'o'|ov): (oo|ov) part ok
c
       if(min(nocc(isw1),nocc(isw2)).eq.0)go to 600                     11d28s22
       do is=1,nsdlk1
        if(isblk1(1,is).eq.isw1.and.isblk1(2,is).eq.isw2.and.
     $     isblk1(3,is).eq.isblkxder1(3,isb).and.
     $       isblk1(4,is).eq.isblkxder1(4,isb))then
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d12s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           if(isblk1(1,is).eq.isblk1(2,is))then
            do i4=0,nocc(isblk1(1,is))-1
             do i3=0,i4-1
              do m2=0,nocc(isblkxder1(2,isb))-1
               iad1=itrans(isw2)+i3+nbasdwsc(isw2)*m2
               fa=2d0*bc(ii)*bc(iad1)
               jad1=itrans(isw2)+i4+nbasdwsc(isw2)*m2
               fb=2d0*bc(ii)*bc(jad1)
               do m1=0,nocc(isblkxder1(1,isb))-1
                iad2=itrans(isw1)+i4+nbasdwsc(isw1)*m1
                jad2=itrans(isw1)+i3+nbasdwsc(isw1)*m1
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $        *(m2+nocc(isblkxder1(2,isb))*(i1m+nocc(isblk1(3,is))*i2m))
                bc(iad3)=bc(iad3)+fa*bc(iad2)+fb*bc(jad2)
               end do
              end do
              ii=ii+1
             end do
             do m2=0,nocc(isblkxder1(2,isb))-1
              iad1=itrans(isw2)+i4+nbasdwsc(isw2)*m2
              fa=2d0*bc(ii)*bc(iad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i4+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $        *(m2+nocc(isblkxder1(2,isb))*(i1m+nocc(isblk1(3,is))*i2m))
               bc(iad3)=bc(iad3)+fa*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           else
            do i4=0,nocc(isblk1(2,is))-1
             do i3=0,nocc(isblk1(1,is))-1
              do m2=0,nocc(isblkxder1(2,isb))-1
               jad1=itrans(isw2)+i4+nbasdwsc(isw2)*m2
               fb=2d0*bc(ii)*bc(jad1)
               do m1=0,nocc(isblkxder1(1,isb))-1
                jad2=itrans(isw1)+i3+nbasdwsc(isw1)*m1
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $        *(m2+nocc(isblkxder1(2,isb))*(i1m+nocc(isblk1(3,is))*i2m))
                bc(iad3)=bc(iad3)+fb*bc(jad2)
               end do
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 600
        end if                                                          9d12s16
        if(isblk1(2,is).eq.isw1.and.isblk1(1,is).eq.isw2.and.
     $     isblk1(3,is).eq.isblkxder1(3,isb).and.
     $       isblk1(4,is).eq.isblkxder1(4,isb))then
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d12s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           do i4=0,nocc(isblk1(2,is))-1
            do i3=0,nocc(isblk1(1,is))-1
             do m2=0,nocc(isblkxder1(2,isb))-1
              jad1=itrans(isw2)+i3+nbasdwsc(isw2)*m2
              fb=2d0*bc(ii)*bc(jad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               jad2=itrans(isw1)+i4+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $        *(m2+nocc(isblkxder1(2,isb))*(i1m+nocc(isblk1(3,is))*i2m))
                bc(iad3)=bc(iad3)+fb*bc(jad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 600
        end if                                                          9d12s16
       end do
       write(6,*)('for isb = '),isb
       write(6,*)('could not find ionex for (o''o''|ov) '),
     $      (isblkxder1(i,isb),i=1,4)
       write(6,*)('what we have for sblk1: ')
       write(6,*)('nsdlk1: '),nsdlk1
       do is=1,nsdlk1                                                   6d13s22
        write(6,*)is,(isblk1(j,is),j=1,4)
       end do                                                           6d13s22
       call dws_sync
       call dws_finalize
       stop
  600  continue
c
c     (o'o'|ov): (vo|ov) part ok
c     recall k_nm^ab=(nb|ma)
c
       if(min(nocc(isw2),nvirtc(isw1)).eq.0)go to 606                   2d24s23
       do is=1,nsdlkk
        if(isblkk(1,is).eq.isw2.and.isblkk(2,is).eq.isblkxder1(3,isb)
     $       .and.isblkk(3,is).eq.isblkxder1(4,isb).and.
     $       isblkk(4,is).eq.isw1)then
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d13s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d13s16
         ii=kmats(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblkk(4,is))
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d13s16
           do i4=0,nocc(isblkk(2,is))-1
            do i3=0,nocc(isblkk(1,is))-1
             do m2=0,nocc(isblkxder1(2,isb))-1
              iad1=itrans(isw2)+i3+nbasdwsc(isw2)*m2                    9d13s16
              fa=2d0*bc(ii)*bc(iad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i2p+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(m2         9d13s16
     $         +nocc(isblkxder1(2,isb))*(i4+nocc(isblkk(2,is))*i1m))    9d13s16
               bc(iad3)=bc(iad3)+bc(iad2)*fa                            9d13s16
              end do
             end do
             ii=ii+1
            end do
           end do
          end do                                                        9d13s16
          i10=1                                                         9d13s16
         end do                                                         9d13s16
         go to 606
        end if
        if(isblkk(2,is).eq.isw2.and.isblkk(1,is).eq.isblkxder1(3,isb)
     $       .and.isblkk(4,is).eq.isblkxder1(4,isb).and.
     $       isblkk(3,is).eq.isw1)then
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d13s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d13s16
         ii=kmats(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblkk(3,is))                                  9d13s16
           do i4=0,nocc(isblkk(2,is))-1
            do i3=0,nocc(isblkk(1,is))-1
             do m2=0,nocc(isblkxder1(2,isb))-1
              iad1=itrans(isw2)+i4+nbasdwsc(isw2)*m2                    9d13s16
              fa=2d0*bc(ii)*bc(iad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i1p+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(m2         9d13s16
     $         +nocc(isblkxder1(2,isb))*(i3+nocc(isblkk(1,is))*i2m))    9d13s16
               bc(iad3)=bc(iad3)+bc(iad2)*fa                            9d13s16
              end do
             end do
             ii=ii+1
            end do
           end do
          end do                                                        9d13s16
          i10=1                                                         9d13s16
         end do                                                         9d13s16
         go to 606
        end if
       end do
       write(6,*)('could not find K for (vo|ov) part of (o''o''|ov)')
       call dws_sync
       call dws_finalize
       stop
  606  continue
c
c     (o'o'|ov): (ov|ov) part ok
c     recall k_nm^ab=(nb|ma)
c
       if(min(nocc(isw1),nvirtc(isw2)).eq.0)go to 607                   2d24s23
       do is=1,nsdlkk
        if(isblkk(1,is).eq.isw1.and.isblkk(2,is).eq.isblkxder1(3,isb)
     $       .and.isblkk(3,is).eq.isblkxder1(4,isb).and.
     $       isblkk(4,is).eq.isw2)then
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d13s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d13s16
         ii=kmats(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblkk(4,is))
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d13s16
           do i4=0,nocc(isblkk(2,is))-1
            do i3=0,nocc(isblkk(1,is))-1
             do m2=0,nocc(isblkxder1(2,isb))-1
              iad1=itrans(isw2)+i2p+nbasdwsc(isw2)*m2                    9d13s16
              fa=2d0*bc(ii)*bc(iad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i3+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(m2         9d13s16
     $         +nocc(isblkxder1(2,isb))*(i4+nocc(isblkk(2,is))*i1m))    9d13s16
               bc(iad3)=bc(iad3)+bc(iad2)*fa                            9d13s16
              end do
             end do
             ii=ii+1
            end do
           end do
          end do                                                        9d13s16
          i10=1                                                         9d13s16
         end do                                                         9d13s16
         go to 607
        end if
        if(isblkk(2,is).eq.isw1.and.isblkk(1,is).eq.isblkxder1(3,isb)
     $       .and.isblkk(4,is).eq.isblkxder1(4,isb).and.
     $       isblkk(3,is).eq.isw2)then
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d13s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d13s16
         ii=kmats(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblkk(3,is))                                  9d13s16
           do i4=0,nocc(isblkk(2,is))-1
            do i3=0,nocc(isblkk(1,is))-1
             do m2=0,nocc(isblkxder1(2,isb))-1
              iad1=itrans(isw2)+i1p+nbasdwsc(isw2)*m2                    9d13s16
              fa=2d0*bc(ii)*bc(iad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i4+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(m2         9d13s16
     $         +nocc(isblkxder1(2,isb))*(i3+nocc(isblkk(1,is))*i2m))    9d13s16
               bc(iad3)=bc(iad3)+bc(iad2)*fa                            9d13s16
              end do
             end do
             ii=ii+1
            end do
           end do
          end do                                                        9d13s16
          i10=1                                                         9d13s16
         end do                                                         9d13s16
         go to 607
        end if
       end do
       write(6,*)('could not find K for (ov|ov) part of (o''o''|ov)')
       call dws_sync
       call dws_finalize
       stop
  607  continue
c
c     (o'o'|ov): (vv|ov) part ok
c
       if(min(nvirtc(isw1),nvirtc(isw2)).eq.0)go to 608                 2d24s23
       do is=1,nsdlk1
        if(isblk1(1,is).eq.isw1.and.isblk1(2,is).eq.isw2.and.
     $     isblk1(3,is).eq.isblkxder1(3,isb).and.
     $       isblk1(4,is).eq.isblkxder1(4,isb))then
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d12s16
         ii=i3x(is)                                                     9d13s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           if(isblk1(1,is).eq.isblk1(2,is))then
            do i4=0,nvirtc(isblk1(1,is))-1
             i4p=i4+nocc(isblk1(1,is))                                  9d13s16
             do i3=0,i4-1
              i3p=i3+nocc(isblk1(1,is))                                 9d13s16
              do m2=0,nocc(isblkxder1(2,isb))-1
               iad1=itrans(isw2)+i3p+nbasdwsc(isw2)*m2                  9d13s16
               fa=2d0*bc(ii)*bc(iad1)
               jad1=itrans(isw2)+i4p+nbasdwsc(isw2)*m2                  9d13s16
               fb=2d0*bc(ii)*bc(jad1)
               do m1=0,nocc(isblkxder1(1,isb))-1
                iad2=itrans(isw1)+i4p+nbasdwsc(isw1)*m1                 9d13s16
                jad2=itrans(isw1)+i3p+nbasdwsc(isw1)*m1                 9d13s16
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $        *(m2+nocc(isblkxder1(2,isb))*(i1m+nocc(isblk1(3,is))*i2m))
                bc(iad3)=bc(iad3)+fa*bc(iad2)+fb*bc(jad2)
               end do
              end do
              ii=ii+1
             end do
             do m2=0,nocc(isblkxder1(2,isb))-1
              iad1=itrans(isw2)+i4p+nbasdwsc(isw2)*m2
              fa=2d0*bc(ii)*bc(iad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i4p+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $        *(m2+nocc(isblkxder1(2,isb))*(i1m+nocc(isblk1(3,is))*i2m))
               bc(iad3)=bc(iad3)+fa*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           else
            do i4=0,nvirtc(isblk1(2,is))-1
             i4p=i4+nocc(isblk1(2,is))                                  9d13s16
             do i3=0,nvirtc(isblk1(1,is))-1
              i3p=i3+nocc(isblk1(1,is))                                 9d13s16
              do m2=0,nocc(isblkxder1(2,isb))-1
               jad1=itrans(isw2)+i4p+nbasdwsc(isw2)*m2                  9d13s16
               fb=2d0*bc(ii)*bc(jad1)
               do m1=0,nocc(isblkxder1(1,isb))-1
                jad2=itrans(isw1)+i3p+nbasdwsc(isw1)*m1                 9d13s16
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $        *(m2+nocc(isblkxder1(2,isb))*(i1m+nocc(isblk1(3,is))*i2m))
                bc(iad3)=bc(iad3)+fb*bc(jad2)
               end do
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 608
        end if                                                          9d12s16
        if(isblk1(2,is).eq.isw1.and.isblk1(1,is).eq.isw2.and.           9d13s16
     $     isblk1(3,is).eq.isblkxder1(3,isb).and.
     $       isblk1(4,is).eq.isblkxder1(4,isb))then
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d12s16
         ii=i3x(is)                                                     9d13s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           do i4=0,nvirtc(isblk1(2,is))-1
            i4p=i4+nocc(isblk1(2,is))                                   9d13s16
            do i3=0,nvirtc(isblk1(1,is))-1
             i3p=i3+nocc(isblk1(1,is))                                  9d13s16
             do m2=0,nocc(isblkxder1(2,isb))-1
              jad1=itrans(isw2)+i3p+nbasdwsc(isw2)*m2                   9d13s16
              fb=2d0*bc(ii)*bc(jad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               jad2=itrans(isw1)+i4p+nbasdwsc(isw1)*m1                  9d13s16
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $        *(m2+nocc(isblkxder1(2,isb))*(i1m+nocc(isblk1(3,is))*i2m))
               bc(iad3)=bc(iad3)+fb*bc(jad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 608
        end if                                                          9d12s16
       end do
       write(6,*)('could not find (vv|ov) for (o''o''|ov) ')
       call dws_sync
       call dws_finalize
       stop
  608  continue
c
c     (o'o|o'v): (oo|ov) part ok
c
      if(min(nocc(isw1),nocc(isw3)).gt.0)then                           2d24s23
       do is=1,nsdlk1
        if(isblk1(1,is).eq.isw1.and.isblk1(2,is).eq.isblkxder1(2,isb)
     $     .and.isblk1(3,is).eq.isw3.and.
     $       isblk1(4,is).eq.isblkxder1(4,isb))then
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d12s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           if(isblk1(1,is).eq.isblk1(2,is))then
            do i4=0,nocc(isblk1(1,is))-1
             do i3=0,i4-1
              do m3=0,nocc(isblkxder1(3,isb))-1
               iad1=itrans(isw3)+i1m+nbasdwsc(isw3)*m3
               fa=2d0*bc(ii)*bc(iad1)
               do m1=0,nocc(isblkxder1(1,isb))-1
                iad2=itrans(isw1)+i4+nbasdwsc(isw1)*m1
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $         *(i3+nocc(isblk1(1,is))*(m3+nocc(isblkxder1(3,isb))*i2m))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
                iad2=itrans(isw1)+i3+nbasdwsc(isw1)*m1
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $         *(i4+nocc(isblk1(1,is))*(m3+nocc(isblkxder1(3,isb))*i2m))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+i1m+nbasdwsc(isw3)*m3
              fa=2d0*bc(ii)*bc(iad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i4+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $         *(i4+nocc(isblk1(1,is))*(m3+nocc(isblkxder1(3,isb))*i2m))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           else
            do i4=0,nocc(isblk1(2,is))-1
             do i3=0,nocc(isblk1(1,is))-1
              do m3=0,nocc(isblkxder1(3,isb))-1
               jad1=itrans(isw3)+i1m+nbasdwsc(isw3)*m3
               fb=2d0*bc(ii)*bc(jad1)
               do m1=0,nocc(isblkxder1(1,isb))-1
                jad2=itrans(isw1)+i3+nbasdwsc(isw1)*m1
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $         *(i4+nocc(isblk1(2,is))*(m3+nocc(isblkxder1(3,isb))*i2m))
                bc(iad3)=bc(iad3)+fb*bc(jad2)
               end do
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 601
        end if                                                          9d12s16
        if(isblk1(2,is).eq.isw1.and.isblk1(1,is).eq.isblkxder1(2,isb)
     $     .and.isblk1(3,is).eq.isw3.and.
     $       isblk1(4,is).eq.isblkxder1(4,isb))then
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d12s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           do i4=0,nocc(isblk1(2,is))-1
            do i3=0,nocc(isblk1(1,is))-1
             do m3=0,nocc(isblkxder1(3,isb))-1
              jad1=itrans(isw3)+i1m+nbasdwsc(isw3)*m3
              fb=2d0*bc(ii)*bc(jad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               jad2=itrans(isw1)+i4+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $         *(i3+nocc(isblk1(1,is))*(m3+nocc(isblkxder1(3,isb))*i2m))
               bc(iad3)=bc(iad3)+fb*bc(jad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 601
        end if                                                          9d12s16
       end do
       write(6,*)('could not find ionex for (o''o|o''v) '),
     $      (isblkxder1(i,isb),i=1,4)
       call dws_sync
       call dws_finalize
       stop
  601  continue
       end if                                                           11d30s16
c
c     (o'o|o'v): (vo|ov) part ok
c     recall k_nm^ab=(nb|ma)
c
       if(min(nocc(isw3),nvirtc(isw1)).ne.0)then                        2d24s23
       do is=1,nsdlkk
        if(isblkk(1,is).eq.isblkxder1(2,isb).and.isblkk(2,is).eq.isw3   9d13s16
     $       .and.isblkk(3,is).eq.isblkxder1(4,isb).and.
     $       isblkk(4,is).eq.isw1)then
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d13s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d13s16
         ii=kmats(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblkk(4,is))
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d13s16
           do i4=0,nocc(isblkk(2,is))-1
            do i3=0,nocc(isblkk(1,is))-1
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+i4+nbasdwsc(isw3)*m3                    9d13s16
              fa=2d0*bc(ii)*bc(iad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i2p+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i3         9d13s16
     $         +nocc(isblkk(1,is))*(m3+nocc(isblkxder1(3,isb))*i1m))    9d13s16
               bc(iad3)=bc(iad3)+bc(iad2)*fa                            9d13s16
              end do
             end do
             ii=ii+1
            end do
           end do
          end do                                                        9d13s16
          i10=1                                                         9d13s16
         end do                                                         9d13s16
         go to 609
        end if                                                          9d13s16
        if(isblkk(2,is).eq.isblkxder1(2,isb).and.isblkk(1,is).eq.isw3   9d13s16
     $       .and.isblkk(4,is).eq.isblkxder1(4,isb).and.
     $       isblkk(3,is).eq.isw1)then
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d13s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d13s16
         ii=kmats(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d13s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblkk(3,is))                                  9d13s16
           do i4=0,nocc(isblkk(2,is))-1
            do i3=0,nocc(isblkk(1,is))-1
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+i3+nbasdwsc(isw3)*m3                    9d13s16
              fa=2d0*bc(ii)*bc(iad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i1p+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i4         9d13s16
     $         +nocc(isblkk(2,is))*(m3+nocc(isblkxder1(3,isb))*i2m))    9d13s16
               bc(iad3)=bc(iad3)+bc(iad2)*fa                            9d13s16
              end do
             end do
             ii=ii+1
            end do
           end do
          end do                                                        9d13s16
          i10=1                                                         9d13s16
         end do                                                         9d13s16
         go to 609
        end if                                                          9d13s16
       end do
       write(6,*)('could not find (vo|ov) for (o''o|o''v)')
       call dws_sync
       call dws_finalize
       stop
  609  continue
       end if                                                           11d30s16
c
c     (o'o|o'v): (oo|vv) part ok
c
       if(min(nocc(isw1),nvirtc(isw3)).ne.0)then                        2d24s23
       do is=1,nsdlk
        if(isblk(1,is).eq.isw1.and.isblk(2,is).eq.isblkxder1(2,isb)     9d13s16
     $       .and.isblk(3,is).eq.isw3.and.
     $       isblk(4,is).eq.isblkxder1(4,isb))then
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d13s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d13s16
         ii=jmats(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblk(3,is))                                   4d20s22
           if(isblk(1,is).eq.isblk(2,is))then
            do i4=0,nocc(isblk(1,is))-1
             do i3=0,i4-1
              do m3=0,nocc(isblkxder1(3,isb))-1
               iad1=itrans(isw3)+i1p+nbasdwsc(isw3)*m3
               fa=2d0*bc(iad1)*bc(ii)
               do m1=0,nocc(isblkxder1(1,isb))-1
                iad2=itrans(isw1)+i4+nbasdwsc(isw1)*m1                  2d24s23
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i3
     $              +nocc(isblk(1,is))*(m3+nocc(isblkxder1(3,isb))*i2m))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
                iad2=itrans(isw1)+i3+nbasdwsc(isw1)*m1
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i4
     $              +nocc(isblk(1,is))*(m3+nocc(isblkxder1(3,isb))*i2m))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+i1p+nbasdwsc(isw3)*m3
              fa=2d0*bc(iad1)*bc(ii)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i4+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i4
     $              +nocc(isblk(1,is))*(m3+nocc(isblkxder1(3,isb))*i2m))
               bc(iad3)=bc(iad3)+fa*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           else
            do i4=0,nocc(isblk(2,is))-1
             do i3=0,nocc(isblk(1,is))-1
              do m3=0,nocc(isblkxder1(3,isb))-1
               iad1=itrans(isw3)+i1p+nbasdwsc(isw3)*m3
               fa=2d0*bc(iad1)*bc(ii)
               do m1=0,nocc(isblkxder1(1,isb))-1
                iad2=itrans(isw1)+i3+nbasdwsc(isw1)*m1
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i4
     $              +nocc(isblk(2,is))*(m3+nocc(isblkxder1(3,isb))*i2m))
                bc(iad3)=bc(iad3)+bc(iad2)*fa
               end do
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 610
        end if
        if(isblk(2,is).eq.isw1.and.isblk(1,is).eq.isblkxder1(2,isb)     9d13s16
     $       .and.isblk(3,is).eq.isw3.and.
     $       isblk(4,is).eq.isblkxder1(4,isb))then
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d13s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d13s16
         ii=jmats(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblk(3,is))                                   4d20s22
           do i4=0,nocc(isblk(2,is))-1
            do i3=0,nocc(isblk(1,is))-1
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+i1p+nbasdwsc(isw3)*m3
              fa=2d0*bc(iad1)*bc(ii)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i4+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i3
     $              +nocc(isblk(1,is))*(m3+nocc(isblkxder1(3,isb))*i2m))
               bc(iad3)=bc(iad3)+bc(iad2)*fa
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 610
        end if
        if(isblk(2,is).eq.isw1.and.isblk(1,is).eq.isblkxder1(2,isb)     9d13s16
     $       .and.isblk(4,is).eq.isw3.and.
     $       isblk(3,is).eq.isblkxder1(4,isb))then
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d13s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d13s16
         ii=jmats(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk(4,is))                                    4d20s22
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1
           do i4=0,nocc(isblk(2,is))-1
            do i3=0,nocc(isblk(1,is))-1
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+i2p+nbasdwsc(isw3)*m3
              fa=2d0*bc(iad1)*bc(ii)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i4+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i3
     $              +nocc(isblk(1,is))*(m3+nocc(isblkxder1(3,isb))*i1m))
               bc(iad3)=bc(iad3)+bc(iad2)*fa
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 610
        end if
        if(isblk(1,is).eq.isw1.and.isblk(2,is).eq.isblkxder1(2,isb)     9d13s16
     $       .and.isblk(4,is).eq.isw3.and.
     $       isblk(3,is).eq.isblkxder1(4,isb))then
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d13s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d13s16
         ii=jmats(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk(4,is))                                    4d20s22
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1
           do i4=0,nocc(isblk(2,is))-1
            do i3=0,nocc(isblk(1,is))-1
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+i2p+nbasdwsc(isw3)*m3
              fa=2d0*bc(iad1)*bc(ii)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i3+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i4
     $              +nocc(isblk(2,is))*(m3+nocc(isblkxder1(3,isb))*i1m))
               bc(iad3)=bc(iad3)+bc(iad2)*fa
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 610
        end if
       end do
       write(6,*)('could not find (oo|vv) for (o''o|o''v) ')
       call dws_sync
       call dws_finalize
       stop
  610  continue
       end if                                                           11d30s16
c
c     (o'o|o'v): (vo|vv) part ok
c
       if(min(nvirtc(isw3),nvirtc(isw1)).eq.0)go to 611                 2d24s23
       do is=1,nsdlk1
        if(isblk1(1,is).eq.isw3.and.isblk1(2,is).eq.isblkxder1(4,isb)
     $     .and.isblk1(3,is).eq.isblkxder1(2,isb).and.
     $       isblk1(4,is).eq.isw1)then
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d12s16
         ii=i3x(is)                                                     9d13s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d13s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           if(isblk1(1,is).eq.isblk1(2,is))then
            do i4=0,nvirtc(isblk1(1,is))-1
             i4p=i4+nocc(isblk1(1,is))                                  9d13s16
             do i3=0,i4-1
              i3p=i3+nocc(isblk1(1,is))                                 9d13s16
              do m3=0,nocc(isblkxder1(3,isb))-1
               iad1=itrans(isw3)+i3p+nbasdwsc(isw3)*m3                  9d13s16
               fa=2d0*bc(ii)*bc(iad1)
               jad1=itrans(isw3)+i4p+nbasdwsc(isw3)*m3                  9d13s16
               fb=2d0*bc(ii)*bc(jad1)
               do m1=0,nocc(isblkxder1(1,isb))-1
                iad2=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                 9d13s16
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $        *(i1m+nocc(isblk1(3,is))*(m3+nocc(isblkxder1(3,isb))*i4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $        *(i1m+nocc(isblk1(3,is))*(m3+nocc(isblkxder1(3,isb))*i3))
                bc(iad3)=bc(iad3)+fb*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+i4p+nbasdwsc(isw3)*m3
              fa=2d0*bc(ii)*bc(iad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i2p+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $        *(i1m+nocc(isblk1(3,is))*(m3+nocc(isblkxder1(3,isb))*i4))
               bc(iad3)=bc(iad3)+fa*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           else
            do i4=0,nvirtc(isblk1(2,is))-1
             do i3=0,nvirtc(isblk1(1,is))-1
              i3p=i3+nocc(isblk1(1,is))                                 9d13s16
              do m3=0,nocc(isblkxder1(3,isb))-1
               jad1=itrans(isw3)+i3p+nbasdwsc(isw3)*m3                  9d13s16
               fb=2d0*bc(ii)*bc(jad1)
               do m1=0,nocc(isblkxder1(1,isb))-1
                jad2=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                 9d13s16
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $        *(i1m+nocc(isblk1(3,is))*(m3+nocc(isblkxder1(3,isb))*i4))
                bc(iad3)=bc(iad3)+fb*bc(jad2)
               end do
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 611
        end if
        if(isblk1(2,is).eq.isw3.and.isblk1(1,is).eq.isblkxder1(4,isb)
     $     .and.isblk1(3,is).eq.isblkxder1(2,isb).and.
     $       isblk1(4,is).eq.isw1)then
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d12s16
         ii=i3x(is)                                                     9d13s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d13s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           do i4=0,nvirtc(isblk1(2,is))-1
            i4p=i4+nocc(isblk1(2,is))                                   9d13s16
            do i3=0,nvirtc(isblk1(1,is))-1
             do m3=0,nocc(isblkxder1(3,isb))-1
              jad1=itrans(isw3)+i4p+nbasdwsc(isw3)*m3                   9d13s16
              fb=2d0*bc(ii)*bc(jad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               jad2=itrans(isw1)+i2p+nbasdwsc(isw1)*m1                  9d13s16
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $        *(i1m+nocc(isblk1(3,is))*(m3+nocc(isblkxder1(3,isb))*i3))
               bc(iad3)=bc(iad3)+fb*bc(jad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 611
        end if
       end do
       write(6,*)('could not find (vo|vv) for (o''o|o''v)')
       call dws_sync
       call dws_finalize
       stop
  611  continue
c
c     (o'o|ov'): (oo|oo) part ok
c
       if(min(nocc(isw1),nvirtc(isw4)).eq.0)go to 602                   2d24s23
       do is=1,nsdlk
        if(isblk(1,is).eq.isw1.and.isblk(2,is).eq.isblkxder1(2,isb)
     $     .and.isblk(3,is).eq.isblkxder1(3,isb).and.
     $       isblk(4,is).eq.isw4)then
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d12s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           if(isblk(1,is).eq.isblk(2,is))then
            do i4=0,nocc(isblk(1,is))-1
             do i3=0,i4-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))                           4d20s22
               iad1=itrans(isw4)+i2m+nbasdwsc(isw4)*m4p
               fa=2d0*bc(ii)*bc(iad1)
               do m1=0,nocc(isblkxder1(1,isb))-1
                iad2=itrans(isw1)+i4+nbasdwsc(isw1)*m1
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $         *(i3+nocc(isblk(1,is))*(i1m+nocc(isblk(3,is))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
                iad2=itrans(isw1)+i3+nbasdwsc(isw1)*m1
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $         *(i4+nocc(isblk(1,is))*(i1m+nocc(isblk(3,is))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+i2m+nbasdwsc(isw4)*m4p
              fa=2d0*bc(ii)*bc(iad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i4+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $         *(i4+nocc(isblk(1,is))*(i1m+nocc(isblk(3,is))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           else
            do i4=0,nocc(isblk(2,is))-1
             do i3=0,nocc(isblk(1,is))-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               jad1=itrans(isw4)+i2m+nbasdwsc(isw4)*m4p
               fb=2d0*bc(ii)*bc(jad1)
               do m1=0,nocc(isblkxder1(1,isb))-1
                jad2=itrans(isw1)+i3+nbasdwsc(isw1)*m1
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $         *(i4+nocc(isblk(2,is))*(i1m+nocc(isblk(3,is))*m4))
                bc(iad3)=bc(iad3)+fb*bc(jad2)
               end do
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 602
        end if                                                          9d12s16
        if(isblk(2,is).eq.isw1.and.isblk(1,is).eq.isblkxder1(2,isb)
     $     .and.isblk(3,is).eq.isblkxder1(3,isb).and.
     $       isblk(4,is).eq.isw4)then
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d12s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           do i4=0,nocc(isblk(2,is))-1
            do i3=0,nocc(isblk(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              jad1=itrans(isw4)+i2m+nbasdwsc(isw4)*m4p
              fb=2d0*bc(ii)*bc(jad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               jad2=itrans(isw1)+i4+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $         *(i3+nocc(isblk(1,is))*(i1m+nocc(isblk(3,is))*m4))
               bc(iad3)=bc(iad3)+fb*bc(jad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 602
        end if                                                          9d12s16
        if(isblk(2,is).eq.isw1.and.isblk(1,is).eq.isblkxder1(2,isb)
     $     .and.isblk(4,is).eq.isblkxder1(3,isb).and.
     $       isblk(3,is).eq.isw4)then
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d12s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           do i4=0,nocc(isblk(2,is))-1
            do i3=0,nocc(isblk(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              jad1=itrans(isw4)+i1m+nbasdwsc(isw4)*m4p
              fb=2d0*bc(ii)*bc(jad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               jad2=itrans(isw1)+i4+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $         *(i3+nocc(isblk(1,is))*(i2m+nocc(isblk(4,is))*m4))
               bc(iad3)=bc(iad3)+fb*bc(jad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 602
        end if                                                          9d12s16
        if(isblk(1,is).eq.isw1.and.isblk(2,is).eq.isblkxder1(2,isb)
     $     .and.isblk(4,is).eq.isblkxder1(3,isb).and.
     $       isblk(3,is).eq.isw4)then
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d12s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           do i4=0,nocc(isblk(2,is))-1
            do i3=0,nocc(isblk(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              jad1=itrans(isw4)+i1m+nbasdwsc(isw4)*m4p
              fb=2d0*bc(ii)*bc(jad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               jad2=itrans(isw1)+i3+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))
     $         *(i4+nocc(isblk(2,is))*(i2m+nocc(isblk(4,is))*m4))
               bc(iad3)=bc(iad3)+fb*bc(jad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 602
        end if                                                          9d12s16
       end do
       write(6,*)('could not find oooo for (o''o|ov'') '),
     $      (isblkxder1(i,isb),i=1,4)
       call dws_sync
       call dws_finalize
       stop
  602  continue
c
c     (o'o|ov'): (vo|oo) part ok
c
       if(min(nocc(isw4),nvirtc(isw1)).eq.0)go to 612                   2d24s23
       do is=1,nsdlk1
        if(isblk1(1,is).eq.isblkxder1(3,isb).and.                       9d14s16
     $       isblk1(2,is).eq.isw4                                       9d14s16
     $       .and.isblk1(3,is).eq.isblkxder1(2,isb).and.                9d14s16
     $       isblk1(4,is).eq.isw1)then                                  9d14s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d14s16
         ii=ionex(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d14s16
           if(isblk1(1,is).eq.isblk1(2,is))then                         9d14s16
            do i4=0,nocc(isblk1(1,is))-1
             do i3=0,i4-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               iad1=itrans(isw4)+i3+nbasdwsc(isw4)*m4p
               fa=2d0*bc(ii)*bc(iad1)
               jad1=itrans(isw4)+i4+nbasdwsc(isw4)*m4p
               fb=2d0*bc(ii)*bc(jad1)
               do m1=0,nocc(isblkxder1(1,isb))-1
                iad2=itrans(isw1)+i2p+nbasdwsc(isw1)*m1
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i1m
     $               +nocc(isblk1(3,is))*(i4+nocc(isblk1(1,is))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i1m
     $               +nocc(isblk1(3,is))*(i3+nocc(isblk1(1,is))*m4))
                bc(iad3)=bc(iad3)+fb*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              jad1=itrans(isw4)+i4+nbasdwsc(isw4)*m4p
              fb=2d0*bc(ii)*bc(jad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i2p+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i1m
     $               +nocc(isblk1(3,is))*(i4+nocc(isblk1(1,is))*m4))
               bc(iad3)=bc(iad3)+fb*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           else
            do i4=0,nocc(isblk1(2,is))-1
             do i3=0,nocc(isblk1(1,is))-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               iad1=itrans(isw4)+i4+nbasdwsc(isw4)*m4p
               fa=2d0*bc(ii)*bc(iad1)
               do m1=0,nocc(isblkxder1(1,isb))-1
                iad2=itrans(isw1)+i2p+nbasdwsc(isw1)*m1
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i1m
     $               +nocc(isblk1(3,is))*(i3+nocc(isblk1(1,is))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 612
        end if
        if(isblk1(2,is).eq.isblkxder1(3,isb).and.                       9d14s16
     $       isblk1(1,is).eq.isw4                                       9d14s16
     $       .and.isblk1(3,is).eq.isblkxder1(2,isb).and.                9d14s16
     $       isblk1(4,is).eq.isw1)then                                  9d14s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d14s16
         ii=ionex(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d14s16
           do i4=0,nocc(isblk1(2,is))-1
            do i3=0,nocc(isblk1(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+i3+nbasdwsc(isw4)*m4p
              fa=2d0*bc(ii)*bc(iad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i2p+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i1m
     $               +nocc(isblk1(3,is))*(i4+nocc(isblk1(2,is))*m4))
               bc(iad3)=bc(iad3)+fa*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 612
        end if
       end do
       write(6,*)('could not find (vo|oo) for (o''o|ov'')')
       call dws_sync
       call dws_finalize
       stop
  612  continue
c
c     (o'o|ov'): (oo|ov) part ok
c
       if(min(nocc(isw1),nvirtc(isw4)).eq.0)go to 613                   2d24s23
       do is=1,nsdlk1
        if(isblk1(1,is).eq.isw1.and.                                    9d14s16
     $       isblk1(2,is).eq.isblkxder1(2,isb)                          9d14s16
     $       .and.isblk1(3,is).eq.isblkxder1(3,isb).and.                9d14s16
     $       isblk1(4,is).eq.isw4)then                                  9d14s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d14s16
         ii=ionex(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d14s16
           if(isblk1(1,is).eq.isblk1(2,is))then                         9d14s16
            do i4=0,nocc(isblk1(1,is))-1
             do i3=0,i4-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               iad1=itrans(isw4)+i2p+nbasdwsc(isw4)*m4p
               fa=2d0*bc(ii)*bc(iad1)
               do m1=0,nocc(isblkxder1(1,isb))-1
                iad2=itrans(isw1)+i3+nbasdwsc(isw1)*m1
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i4
     $               +nocc(isblk1(1,is))*(i1m+nocc(isblk1(3,is))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
                iad2=itrans(isw1)+i4+nbasdwsc(isw1)*m1
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i3
     $               +nocc(isblk1(1,is))*(i1m+nocc(isblk1(3,is))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              jad1=itrans(isw4)+i2p+nbasdwsc(isw4)*m4p
              fb=2d0*bc(ii)*bc(jad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i4+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i4
     $               +nocc(isblk1(1,is))*(i1m+nocc(isblk1(3,is))*m4))
               bc(iad3)=bc(iad3)+fb*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           else
            do i4=0,nocc(isblk1(2,is))-1
             do i3=0,nocc(isblk1(1,is))-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               iad1=itrans(isw4)+i2p+nbasdwsc(isw4)*m4p
               fa=2d0*bc(ii)*bc(iad1)
               do m1=0,nocc(isblkxder1(1,isb))-1
                iad2=itrans(isw1)+i3+nbasdwsc(isw1)*m1
                iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i4
     $               +nocc(isblk1(2,is))*(i1m+nocc(isblk1(3,is))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 613
        end if
        if(isblk1(2,is).eq.isw1.and.                                    9d14s16
     $       isblk1(1,is).eq.isblkxder1(2,isb)                          9d14s16
     $       .and.isblk1(3,is).eq.isblkxder1(3,isb).and.                9d14s16
     $       isblk1(4,is).eq.isw4)then                                  9d14s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d14s16
         ii=ionex(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d14s16
           do i4=0,nocc(isblk1(2,is))-1
            do i3=0,nocc(isblk1(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+i2p+nbasdwsc(isw4)*m4p
              fa=2d0*bc(ii)*bc(iad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i4+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i3
     $               +nocc(isblk1(1,is))*(i1m+nocc(isblk1(3,is))*m4))
               bc(iad3)=bc(iad3)+fa*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 613
        end if
       end do
       write(6,*)('could not find (oo|ov) for (o''o|ov'')')
       call dws_sync
       call dws_finalize
       stop
  613  continue
c
c     (o'o|ov'): (vo|ov) part ok
c     recall k_nm^ab=(nb|ma)
c
       if(min(nvirtc(isw1),nvirtc(isw4)).eq.0)go to 614                 2d24s23
       do is=1,nsdlkk
        if(isblkk(1,is).eq.isblkxder1(2,isb).and.                       9d14s16
     $       isblkk(2,is).eq.isblkxder1(3,isb)                          9d14s16
     $       .and.isblkk(3,is).eq.isw4.and.                             9d14s16
     $       isblkk(4,is).eq.isw1)then
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d13s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d13s16
         ii=kmats(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblkk(4,is))
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblkk(3,is))                                  9d14s16
           do i4=0,nocc(isblkk(2,is))-1
            do i3=0,nocc(isblkk(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+i1p+nbasdwsc(isw4)*m4p                  9d14s16
              fa=2d0*bc(ii)*bc(iad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i2p+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i3         9d13s16
     $         +nocc(isblkk(1,is))*(i4+nocc(isblkk(2,is))*m4))          9d14s16
               bc(iad3)=bc(iad3)+bc(iad2)*fa                            9d13s16
              end do
             end do
             ii=ii+1
            end do
           end do
          end do                                                        9d13s16
          i10=1                                                         9d13s16
         end do                                                         9d13s16
         go to 614
        end if
        if(isblkk(2,is).eq.isblkxder1(2,isb).and.                       9d14s16
     $       isblkk(1,is).eq.isblkxder1(3,isb)                          9d14s16
     $       .and.isblkk(4,is).eq.isw4.and.                             9d14s16
     $       isblkk(3,is).eq.isw1)then
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d13s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d13s16
         ii=kmats(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblkk(4,is))
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblkk(3,is))                                  9d14s16
           do i4=0,nocc(isblkk(2,is))-1
            do i3=0,nocc(isblkk(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+i2p+nbasdwsc(isw4)*m4p                  9d14s16
              fa=2d0*bc(ii)*bc(iad1)
              do m1=0,nocc(isblkxder1(1,isb))-1
               iad2=itrans(isw1)+i1p+nbasdwsc(isw1)*m1
               iad3=ionexd2(isb)+m1+nocc(isblkxder1(1,isb))*(i4         9d13s16
     $         +nocc(isblkk(2,is))*(i3+nocc(isblkk(1,is))*m4))          9d14s16
               bc(iad3)=bc(iad3)+bc(iad2)*fa                            9d13s16
              end do
             end do
             ii=ii+1
            end do
           end do
          end do                                                        9d13s16
          i10=1                                                         9d13s16
         end do                                                         9d13s16
         go to 614
        end if
       end do
       write(6,*)('could not find (vo|ov) for (o''o|ov'')')
       call dws_sync
       call dws_finalize
       stop
  614  continue
c
c     (oo'|o'v): (oo|ov) part ok
c
       if(min(nocc(isw2),nocc(isw3)).ne.0)then                          2d24s23
       do is=1,nsdlk1
        if(isblk1(1,is).eq.isblkxder1(1,isb).and.isblk1(2,is).eq.isw2   9d12s16
     $     .and.isblk1(3,is).eq.isw3.and.
     $       isblk1(4,is).eq.isblkxder1(4,isb))then
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d12s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           if(isblk1(1,is).eq.isblk1(2,is))then
            do i4=0,nocc(isblk1(1,is))-1
             do i3=0,i4-1
              do m3=0,nocc(isblkxder1(3,isb))-1
               iad1=itrans(isw3)+i1m+nbasdwsc(isw3)*m3
               fa=2d0*bc(ii)*bc(iad1)
               do m2=0,nocc(isblkxder1(2,isb))-1
                iad2=itrans(isw2)+i4+nbasdwsc(isw2)*m2
                iad3=ionexd2(isb)+i3+nocc(isblk1(1,is))
     $    *(m2+nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i2m))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
                iad2=itrans(isw2)+i3+nbasdwsc(isw2)*m2
                iad3=ionexd2(isb)+i4+nocc(isblk1(1,is))
     $    *(m2+nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i2m))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+i1m+nbasdwsc(isw3)*m3
              fa=2d0*bc(ii)*bc(iad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               iad2=itrans(isw2)+i4+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i4+nocc(isblk1(1,is))
     $    *(m2+nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i2m))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           else
            do i4=0,nocc(isblk1(2,is))-1
             do i3=0,nocc(isblk1(1,is))-1
              do m3=0,nocc(isblkxder1(3,isb))-1
               jad1=itrans(isw3)+i1m+nbasdwsc(isw3)*m3
               fb=2d0*bc(ii)*bc(jad1)
               do m2=0,nocc(isblkxder1(2,isb))-1
                jad2=itrans(isw2)+i4+nbasdwsc(isw2)*m2
                iad3=ionexd2(isb)+i3+nocc(isblk1(1,is))
     $    *(m2+nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i2m))
                bc(iad3)=bc(iad3)+fb*bc(jad2)
               end do
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 603
        end if                                                          9d12s16
        if(isblk1(2,is).eq.isblkxder1(1,isb).and.isblk1(1,is).eq.isw2   9d12s16
     $     .and.isblk1(3,is).eq.isw3.and.
     $       isblk1(4,is).eq.isblkxder1(4,isb))then
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d12s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           do i4=0,nocc(isblk1(2,is))-1
            do i3=0,nocc(isblk1(1,is))-1
             do m3=0,nocc(isblkxder1(3,isb))-1
              jad1=itrans(isw3)+i1m+nbasdwsc(isw3)*m3
              fb=2d0*bc(ii)*bc(jad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               jad2=itrans(isw2)+i3+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i4+nocc(isblk1(2,is))
     $    *(m2+nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i2m))
               bc(iad3)=bc(iad3)+fb*bc(jad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 603
        end if                                                          9d12s16
       end do
       write(6,*)('could not find ionex for (oo''|o''v) '),
     $      (isblkxder1(i,isb),i=1,4)
       call dws_sync
       call dws_finalize
       stop
  603  continue
       end if                                                           11d30s16
c
c     (oo'|o'v): (ov|ov) part ok
c     recall k_nm^ab=(nb|ma)
c
       if(min(nocc(isw3),nvirtc(isw2)).ne.0)then                        2d24s23
       do is=1,nsdlkk
        if(isblkk(1,is).eq.isblkxder1(1,isb).and.                       9d14s16
     $       isblkk(2,is).eq.isw3                                       9d14s16
     $       .and.isblkk(3,is).eq.isblkxder1(4,isb).and.                9d14s16
     $       isblkk(4,is).eq.isw2)then                                  9d14s16
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d13s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d13s16
         ii=kmats(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblkk(4,is))
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d14s16
           do i4=0,nocc(isblkk(2,is))-1
            do i3=0,nocc(isblkk(1,is))-1
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+i4+nbasdwsc(isw3)*m3                    9d14s16
              fa=2d0*bc(ii)*bc(iad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               iad2=itrans(isw2)+i2p+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i3+nocc(isblkk(1,is))*(m2              9d14s16
     $        +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i1m))9d14s16
               bc(iad3)=bc(iad3)+bc(iad2)*fa                            9d13s16
              end do
             end do
             ii=ii+1
            end do
           end do
          end do                                                        9d13s16
          i10=1                                                         9d13s16
         end do                                                         9d13s16
         go to 615
        end if
        if(isblkk(2,is).eq.isblkxder1(1,isb).and.                       9d14s16
     $       isblkk(1,is).eq.isw3                                       9d14s16
     $       .and.isblkk(4,is).eq.isblkxder1(4,isb).and.                9d14s16
     $       isblkk(3,is).eq.isw2)then                                  9d14s16
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d13s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d13s16
         ii=kmats(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2m=i2-1                                                      9d14s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblkk(3,is))                                  9d14s16
           do i4=0,nocc(isblkk(2,is))-1
            do i3=0,nocc(isblkk(1,is))-1
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+i3+nbasdwsc(isw3)*m3                    9d14s16
              fa=2d0*bc(ii)*bc(iad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               iad2=itrans(isw2)+i1p+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i4+nocc(isblkk(2,is))*(m2              9d14s16
     $        +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i2m))9d14s16
               bc(iad3)=bc(iad3)+bc(iad2)*fa                            9d13s16
              end do
             end do
             ii=ii+1
            end do
           end do
          end do                                                        9d13s16
          i10=1                                                         9d13s16
         end do                                                         9d13s16
         go to 615
        end if
       end do
       write(6,*)('could not find (ov|ov) for (oo''|o''v)')
       call dws_sync
       call dws_finalize
       stop
  615  continue
       end if                                                           11d30s16
c
c     (oo'|o'v): (oo|vv) part ok
c
       if(min(nocc(isw2),nvirtc(isw3)).ne.0)then                        2d24s23
       do is=1,nsdlk
        if(isblk(1,is).eq.isblkxder1(1,isb).and.
     $     isblk(2,is).eq.isw2.and.
     $     isblk(3,is).eq.isw3.and.
     $     isblk(4,is).eq.isblkxder1(4,isb))then
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d14s16
         ii=jmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblk(3,is))
           if(isblk(1,is).eq.isblk(2,is))then
            do i4=0,nocc(isblk(1,is))-1
             do i3=0,i4-1
              do m3=0,nocc(isblkxder1(3,isb))-1
               iad1=itrans(isw3)+i1p+nbasdwsc(isw3)*m3
               fa=2d0*bc(ii)*bc(iad1)
               do m2=0,nocc(isblkxder1(2,isb))-1
                iad2=itrans(isw2)+i3+nbasdwsc(isw2)*m2
                iad3=ionexd2(isb)+i4+nocc(isblk(1,is))*(m2
     $        +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i2m))
                bc(iad3)=bc(iad3)+bc(iad2)*fa
                iad2=itrans(isw2)+i4+nbasdwsc(isw2)*m2
                iad3=ionexd2(isb)+i3+nocc(isblk(1,is))*(m2
     $        +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i2m))
                bc(iad3)=bc(iad3)+bc(iad2)*fa
               end do
              end do
              ii=ii+1
             end do
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+i1p+nbasdwsc(isw3)*m3
              fa=2d0*bc(ii)*bc(iad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               iad2=itrans(isw2)+i4+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i4+nocc(isblk(1,is))*(m2
     $        +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i2m))
               bc(iad3)=bc(iad3)+bc(iad2)*fa
              end do
             end do
             ii=ii+1
            end do
           else
            do i4=0,nocc(isblk(2,is))-1
             do i3=0,nocc(isblk(1,is))-1
              do m3=0,nocc(isblkxder1(3,isb))-1
               iad1=itrans(isw3)+i1p+nbasdwsc(isw3)*m3
               fa=2d0*bc(ii)*bc(iad1)
               do m2=0,nocc(isblkxder1(2,isb))-1
                iad2=itrans(isw2)+i4+nbasdwsc(isw2)*m2
                iad3=ionexd2(isb)+i3+nocc(isblk(1,is))*(m2
     $        +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i2m))
                bc(iad3)=bc(iad3)+bc(iad2)*fa
               end do
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 616
        end if
        if(isblk(2,is).eq.isblkxder1(1,isb).and.
     $     isblk(1,is).eq.isw2.and.
     $     isblk(3,is).eq.isw3.and.
     $     isblk(4,is).eq.isblkxder1(4,isb))then
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d14s16
         ii=jmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblk(3,is))
           do i4=0,nocc(isblk(2,is))-1
            do i3=0,nocc(isblk(1,is))-1
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+i1p+nbasdwsc(isw3)*m3
              fa=2d0*bc(ii)*bc(iad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               iad2=itrans(isw2)+i3+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i4+nocc(isblk(2,is))*(m2
     $        +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i2m))
               bc(iad3)=bc(iad3)+bc(iad2)*fa
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 616
        end if
        if(isblk(2,is).eq.isblkxder1(1,isb).and.
     $     isblk(1,is).eq.isw2.and.
     $     isblk(4,is).eq.isw3.and.
     $     isblk(3,is).eq.isblkxder1(4,isb))then
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d14s16
         ii=jmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2p=i2-1+nocc(isblk(4,is))                                    9d14s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1
           do i4=0,nocc(isblk(2,is))-1
            do i3=0,nocc(isblk(1,is))-1
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+i2p+nbasdwsc(isw3)*m3
              fa=2d0*bc(ii)*bc(iad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               iad2=itrans(isw2)+i3+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i4+nocc(isblk(2,is))*(m2
     $        +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i1m))
               bc(iad3)=bc(iad3)+bc(iad2)*fa
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 616
        end if
        if(isblk(1,is).eq.isblkxder1(1,isb).and.
     $     isblk(2,is).eq.isw2.and.
     $     isblk(4,is).eq.isw3.and.
     $     isblk(3,is).eq.isblkxder1(4,isb))then
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d14s16
         ii=jmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2p=i2-1+nocc(isblk(4,is))                                    9d14s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1
           do i4=0,nocc(isblk(2,is))-1
            do i3=0,nocc(isblk(1,is))-1
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+i2p+nbasdwsc(isw3)*m3
              fa=2d0*bc(ii)*bc(iad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               iad2=itrans(isw2)+i4+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i3+nocc(isblk(1,is))*(m2
     $        +nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i1m))
               bc(iad3)=bc(iad3)+bc(iad2)*fa
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 616
        end if
       end do
       write(6,*)('could not find (oo|vv) for (oo''|o''v)')
       call dws_sync
       call dws_finalize
       stop
  616  continue
       end if                                                           11d30s16
c
c     (oo'|o'v): (ov|vv) part ok
c
       if(min(nvirt(isw3),nvirtc(isw2)).eq.0)go to 617                  2d24s23
       do is=1,nsdlk1
        if(isblk1(1,is).eq.isw3.and.
     $     isblk1(2,is).eq.isblkxder1(4,isb).and.
     $     isblk1(3,is).eq.isblkxder1(1,isb).and.
     $       isblk1(4,is).eq.isw2)then
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d14s16
         ii=i3x(is)                                                     9d14s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d14s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1
           if(isblk1(1,is).eq.isblk1(2,is))then
            do i4=0,nvirtc(isblk1(1,is))-1
             i4p=i4+nocc(isblk1(1,is))
             do i3=0,i4-1
              i3p=i3+nocc(isblk1(1,is))
              do m3=0,nocc(isblkxder1(3,isb))-1
               iad1=itrans(isw3)+i3p+nbasdwsc(isw3)*m3
               fa=2d0*bc(ii)*bc(iad1)                                       9d14s16
               jad1=itrans(isw3)+i4p+nbasdwsc(isw3)*m3
               fb=2d0*bc(ii)*bc(jad1)                                       9d14s16
               do m2=0,nocc(isblkxder1(2,isb))-1
                iad2=itrans(isw2)+i2p+nbasdwsc(isw2)*m2
                iad3=ionexd2(isb)+i1m+nocc(isblk1(3,is))*(m2+
     $          nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
                iad3=ionexd2(isb)+i1m+nocc(isblk1(3,is))*(m2+
     $          nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i3))
                bc(iad3)=bc(iad3)+fb*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
             do m3=0,nocc(isblkxder1(3,isb))-1
              jad1=itrans(isw3)+i4p+nbasdwsc(isw3)*m3
              fb=2d0*bc(ii)*bc(jad1)                                        9d14s16
              do m2=0,nocc(isblkxder1(2,isb))-1
               iad2=itrans(isw2)+i2p+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i1m+nocc(isblk1(3,is))*(m2+
     $          nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i4))
               bc(iad3)=bc(iad3)+fb*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           else
            do i4=0,nvirtc(isblk1(2,is))-1
             do i3=0,nvirtc(isblk1(1,is))-1
              i3p=i3+nocc(isblk1(1,is))
              do m3=0,nocc(isblkxder1(3,isb))-1
               iad1=itrans(isw3)+i3p+nbasdwsc(isw3)*m3
               fa=2d0*bc(ii)*bc(iad1)
               do m2=0,nocc(isblkxder1(2,isb))-1
                iad2=itrans(isw2)+i2p+nbasdwsc(isw2)*m2
                iad3=ionexd2(isb)+i1m+nocc(isblk1(3,is))*(m2+
     $          nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
            end do
           end if
          end do                                                        9d14s16
          i10=1
         end do                                                         9d14s16
         go to 617
        end if
        if(isblk1(2,is).eq.isw3.and.
     $     isblk1(1,is).eq.isblkxder1(4,isb).and.
     $     isblk1(3,is).eq.isblkxder1(1,isb).and.
     $       isblk1(4,is).eq.isw2)then
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d14s16
         ii=i3x(is)                                                     9d14s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                   9d14s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1
           do i4=0,nvirtc(isblk1(2,is))-1
            i4p=i4+nocc(isblk1(2,is))
            do i3=0,nvirtc(isblk1(1,is))-1
             do m3=0,nocc(isblkxder1(3,isb))-1
              iad1=itrans(isw3)+i4p+nbasdwsc(isw3)*m3
              fa=2d0*bc(ii)*bc(iad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               iad2=itrans(isw2)+i2p+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i1m+nocc(isblk1(3,is))*(m2+
     $          nocc(isblkxder1(2,isb))*(m3+nocc(isblkxder1(3,isb))*i3))
               bc(iad3)=bc(iad3)+fa*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do                                                        9d14s16
          i10=1
         end do                                                         9d14s16
         go to 617
        end if
       end do
       write(6,*)('could not find (ov|vv) for (oo''|o''v)')
       call dws_sync
       call dws_finalize
       stop
  617  continue
c
c     (oo'|ov'): (oo|oo) part ok
c
       if(min(nocc(isw2),nocc(isw4)).ne.0)then                          2d24s23
       do is=1,nsdlk
        if(isblk(1,is).eq.isblkxder1(1,isb).and.isblk(2,is).eq.isw2
     $     .and.isblk(3,is).eq.isblkxder1(3,isb).and.
     $       isblk(4,is).eq.isw4)then
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d12s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           if(isblk(1,is).eq.isblk(2,is))then
            do i4=0,nocc(isblk(1,is))-1
             do i3=0,i4-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               iad1=itrans(isw4)+i2m+nbasdwsc(isw4)*m4p
               fa=2d0*bc(ii)*bc(iad1)
               do m2=0,nocc(isblkxder1(2,isb))-1
                iad2=itrans(isw2)+i4+nbasdwsc(isw2)*m2
                iad3=ionexd2(isb)+i3+nocc(isblk(1,is))
     $         *(m2+nocc(isblkxder1(2,isb))*(i1m+nocc(isblk(3,is))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
                iad2=itrans(isw2)+i3+nbasdwsc(isw2)*m2
                iad3=ionexd2(isb)+i4+nocc(isblk(1,is))
     $         *(m2+nocc(isblkxder1(2,isb))*(i1m+nocc(isblk(3,is))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+i2m+nbasdwsc(isw4)*m4p
              fa=2d0*bc(ii)*bc(iad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               iad2=itrans(isw2)+i4+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i4+nocc(isblk(1,is))
     $         *(m2+nocc(isblkxder1(2,isb))*(i1m+nocc(isblk(3,is))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           else
            do i4=0,nocc(isblk(2,is))-1
             do i3=0,nocc(isblk(1,is))-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               jad1=itrans(isw4)+i2m+nbasdwsc(isw4)*m4p
               fb=2d0*bc(ii)*bc(jad1)
               do m2=0,nocc(isblkxder1(2,isb))-1
                jad2=itrans(isw2)+i4+nbasdwsc(isw2)*m2
                iad3=ionexd2(isb)+i3+nocc(isblk(1,is))
     $         *(m2+nocc(isblkxder1(2,isb))*(i1m+nocc(isblk(3,is))*m4))
                bc(iad3)=bc(iad3)+fb*bc(jad2)
               end do
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 604
        end if                                                          9d12s16
        if(isblk(2,is).eq.isblkxder1(1,isb).and.isblk(1,is).eq.isw2
     $     .and.isblk(3,is).eq.isblkxder1(3,isb).and.
     $       isblk(4,is).eq.isw4)then
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d12s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           do i4=0,nocc(isblk(2,is))-1
            do i3=0,nocc(isblk(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              jad1=itrans(isw4)+i2m+nbasdwsc(isw4)*m4p
              fb=2d0*bc(ii)*bc(jad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               jad2=itrans(isw2)+i3+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i4+nocc(isblk(2,is))
     $         *(m2+nocc(isblkxder1(2,isb))*(i1m+nocc(isblk(3,is))*m4))
               bc(iad3)=bc(iad3)+fb*bc(jad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 604
        end if                                                          9d12s16
        if(isblk(2,is).eq.isblkxder1(1,isb).and.isblk(1,is).eq.isw2
     $     .and.isblk(4,is).eq.isblkxder1(3,isb).and.
     $       isblk(3,is).eq.isw4)then
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d12s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           do i4=0,nocc(isblk(2,is))-1
            do i3=0,nocc(isblk(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              jad1=itrans(isw4)+i1m+nbasdwsc(isw4)*m4p
              fb=2d0*bc(ii)*bc(jad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               jad2=itrans(isw2)+i3+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i4+nocc(isblk(2,is))
     $         *(m2+nocc(isblkxder1(2,isb))*(i2m+nocc(isblk(4,is))*m4))
               bc(iad3)=bc(iad3)+fb*bc(jad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 604
        end if                                                          9d12s16
        if(isblk(1,is).eq.isblkxder1(1,isb).and.isblk(2,is).eq.isw2
     $     .and.isblk(4,is).eq.isblkxder1(3,isb).and.
     $       isblk(3,is).eq.isw4)then
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d12s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           do i4=0,nocc(isblk(2,is))-1
            do i3=0,nocc(isblk(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              jad1=itrans(isw4)+i1m+nbasdwsc(isw4)*m4p
              fb=2d0*bc(ii)*bc(jad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               jad2=itrans(isw2)+i4+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i3+nocc(isblk(1,is))
     $         *(m2+nocc(isblkxder1(2,isb))*(i2m+nocc(isblk(4,is))*m4))
               bc(iad3)=bc(iad3)+fb*bc(jad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 604
        end if                                                          9d12s16
       end do
       write(6,*)('could not find oooo for (oo''|ov'') '),
     $      (isblkxder1(i,isb),i=1,4)
       call dws_sync
       call dws_finalize
       stop
  604  continue
       end if
c
c     (oo'|ov'): (ov|oo) part ok
c
       if(min(nocc(isw4),nbasdws(isw2)).ne.0)then                       11d28s22
       do is=1,nsdlk1
        if(isblk1(1,is).eq.isblkxder1(3,isb).and.isblk1(2,is).eq.isw4   9d12s16
     $     .and.isblk1(3,is).eq.isblkxder1(1,isb).and.                  9d14s16
     $       isblk1(4,is).eq.isw2)then                                  9d14s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d12s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1
           if(isblk1(1,is).eq.isblk1(2,is))then
            do i4=0,nocc(isblk1(1,is))-1
             do i3=0,i4-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               iad1=itrans(isw4)+i3+nbasdwsc(isw4)*m4p
               fa=2d0*bc(ii)*bc(iad1)
               jad1=itrans(isw4)+i4+nbasdwsc(isw4)*m4p
               fb=2d0*bc(ii)*bc(jad1)
               do m2=0,nocc(isblkxder1(2,isb))-1
                iad2=itrans(isw2)+i2p+nbasdwsc(isw2)*m2
                iad3=ionexd2(isb)+i1m+nocc(isblk1(3,is))*(m2+
     $               nocc(isblkxder1(2,isb))*(i4+nocc(isblk1(1,is))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
                iad3=ionexd2(isb)+i1m+nocc(isblk1(3,is))*(m2+
     $               nocc(isblkxder1(2,isb))*(i3+nocc(isblk1(1,is))*m4))
                bc(iad3)=bc(iad3)+fb*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+i4+nbasdwsc(isw4)*m4p
              fa=2d0*bc(ii)*bc(iad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               iad2=itrans(isw2)+i2p+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i1m+nocc(isblk1(3,is))*(m2+
     $              nocc(isblkxder1(2,isb))*(i4+nocc(isblk1(1,is))*m4))
               bc(iad3)=bc(iad3)+fa*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           else
            do i4=0,nocc(isblk1(2,is))-1
             do i3=0,nocc(isblk1(1,is))-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               iad1=itrans(isw4)+i4+nbasdwsc(isw4)*m4p
               fa=2d0*bc(ii)*bc(iad1)
               do m2=0,nocc(isblkxder1(2,isb))-1
                iad2=itrans(isw2)+i2p+nbasdwsc(isw2)*m2
                iad3=ionexd2(isb)+i1m+nocc(isblk1(3,is))*(m2+
     $               nocc(isblkxder1(2,isb))*(i3+nocc(isblk1(1,is))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 618
        end if
        if(isblk1(2,is).eq.isblkxder1(3,isb).and.isblk1(1,is).eq.isw4   9d12s16
     $     .and.isblk1(3,is).eq.isblkxder1(1,isb).and.                  9d14s16
     $       isblk1(4,is).eq.isw2)then                                  9d14s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d12s16
         ii=ionex(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2p=i2-1+nocc(isblk1(4,is))                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1
           do i4=0,nocc(isblk1(2,is))-1
            do i3=0,nocc(isblk1(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+i3+nbasdwsc(isw4)*m4p
              fa=2d0*bc(ii)*bc(iad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               iad2=itrans(isw2)+i2p+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i1m+nocc(isblk1(3,is))*(m2+
     $              nocc(isblkxder1(2,isb))*(i4+nocc(isblk1(2,is))*m4))
               bc(iad3)=bc(iad3)+fa*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 618
        end if
       end do
       write(6,*)('could not find (ov|oo) for (oo''|ov'')')
       call dws_sync
       call dws_finalize
       stop
  618  continue
       end if
c
c     (oo'|ov'): (oo|ov) part ok
c
       if(min(nocc(isw2),nbasdwsc(isw4)).ne.0)then                      11d28s22
       do is=1,nsdlk1
        if(isblk1(1,is).eq.isblkxder1(1,isb).and.                                    9d14s16
     $       isblk1(2,is).eq.isw2
     $       .and.isblk1(3,is).eq.isblkxder1(3,isb).and.                9d14s16
     $       isblk1(4,is).eq.isw4)then                                  9d14s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d14s16
         ii=ionex(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d14s16
           if(isblk1(1,is).eq.isblk1(2,is))then                         9d14s16
            do i4=0,nocc(isblk1(1,is))-1
             do i3=0,i4-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               iad1=itrans(isw4)+i2p+nbasdwsc(isw4)*m4p
               fa=2d0*bc(ii)*bc(iad1)
               do m2=0,nocc(isblkxder1(2,isb))-1
                iad2=itrans(isw2)+i3+nbasdwsc(isw2)*m2
                iad3=ionexd2(isb)+i4+nocc(isblk1(1,is))*(m2
     $             +nocc(isblkxder1(2,isb))*(i1m+nocc(isblk1(3,is))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
                iad2=itrans(isw2)+i4+nbasdwsc(isw2)*m2
                iad3=ionexd2(isb)+i3+nocc(isblk1(1,is))*(m2
     $             +nocc(isblkxder1(2,isb))*(i1m+nocc(isblk1(3,is))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              jad1=itrans(isw4)+i2p+nbasdwsc(isw4)*m4p
              fb=2d0*bc(ii)*bc(jad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               iad2=itrans(isw2)+i4+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i4+nocc(isblk1(1,is))*(m2
     $             +nocc(isblkxder1(2,isb))*(i1m+nocc(isblk1(3,is))*m4))
               bc(iad3)=bc(iad3)+fb*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           else
            do i4=0,nocc(isblk1(2,is))-1
             do i3=0,nocc(isblk1(1,is))-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               iad1=itrans(isw4)+i2p+nbasdwsc(isw4)*m4p
               fa=2d0*bc(ii)*bc(iad1)
               do m2=0,nocc(isblkxder1(2,isb))-1
                iad2=itrans(isw2)+i4+nbasdwsc(isw2)*m2
                iad3=ionexd2(isb)+i3+nocc(isblk1(1,is))*(m2
     $             +nocc(isblkxder1(2,isb))*(i1m+nocc(isblk1(3,is))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 619
        end if
        if(isblk1(2,is).eq.isblkxder1(1,isb).and.                                    9d14s16
     $       isblk1(1,is).eq.isw2
     $       .and.isblk1(3,is).eq.isblkxder1(3,isb).and.                9d14s16
     $       isblk1(4,is).eq.isw4)then                                  9d14s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d14s16
         ii=ionex(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d14s16
           do i4=0,nocc(isblk1(2,is))-1
            do i3=0,nocc(isblk1(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+i2p+nbasdwsc(isw4)*m4p
              fa=2d0*bc(ii)*bc(iad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               iad2=itrans(isw2)+i3+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i4+nocc(isblk1(2,is))*(m2
     $             +nocc(isblkxder1(2,isb))*(i1m+nocc(isblk1(3,is))*m4))
               bc(iad3)=bc(iad3)+fa*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 619
        end if
       end do
       write(6,*)('could not find (oo|ov) for (oo''|ov'')')
       call dws_sync
       call dws_finalize
       stop
  619  continue
       end if
c
c     (oo'|ov'): (ov|ov) part ok
c     recall k_nm^ab=(nb|ma)
c
       if(min(nvirtc(isw2),nvirtc(isw4)).eq.0)go to 620                 2d24s23
       do is=1,nsdlkk
        if(isblkk(1,is).eq.isblkxder1(1,isb).and.                       9d14s16
     $       isblkk(2,is).eq.isblkxder1(3,isb)                          9d14s16
     $       .and.isblkk(3,is).eq.isw4.and.                             9d14s16
     $       isblkk(4,is).eq.isw2)then
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d13s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d13s16
         ii=kmats(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblkk(4,is))
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblkk(3,is))                                  9d14s16
           do i4=0,nocc(isblkk(2,is))-1
            do i3=0,nocc(isblkk(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+i1p+nbasdwsc(isw4)*m4p                  9d14s16
              fa=2d0*bc(ii)*bc(iad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               iad2=itrans(isw2)+i2p+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i3+nocc(isblkk(1,is))*(m2              9d16s16
     $         +nocc(isblkxder1(2,isb))*(i4+nocc(isblkk(2,is))*m4))     9d16s16
               bc(iad3)=bc(iad3)+bc(iad2)*fa                            9d13s16
              end do
             end do
             ii=ii+1
            end do
           end do
          end do                                                        9d13s16
          i10=1                                                         9d13s16
         end do                                                         9d13s16
         go to 620
        end if
        if(isblkk(2,is).eq.isblkxder1(1,isb).and.                       9d14s16
     $       isblkk(1,is).eq.isblkxder1(3,isb)                          9d14s16
     $       .and.isblkk(4,is).eq.isw4.and.                             9d14s16
     $       isblkk(3,is).eq.isw2)then
         call ilimts(nvirtc(isblkk(3,is)),nvirtc(isblkk(4,is)),mynprocg,9d13s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblkk(3,is))                                       9d13s16
         ii=kmats(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblkk(4,is))
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblkk(3,is))                                  9d14s16
           do i4=0,nocc(isblkk(2,is))-1
            do i3=0,nocc(isblkk(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+i2p+nbasdwsc(isw4)*m4p                  9d14s16
              fa=2d0*bc(ii)*bc(iad1)
              do m2=0,nocc(isblkxder1(2,isb))-1
               iad2=itrans(isw2)+i1p+nbasdwsc(isw2)*m2
               iad3=ionexd2(isb)+i4+nocc(isblkk(2,is))*(m2              9d16s16
     $         +nocc(isblkxder1(2,isb))*(i3+nocc(isblkk(1,is))*m4))     9d16s16
               bc(iad3)=bc(iad3)+bc(iad2)*fa                            9d13s16
              end do
             end do
             ii=ii+1
            end do
           end do
          end do                                                        9d13s16
          i10=1                                                         9d13s16
         end do                                                         9d13s16
         go to 620
        end if
       end do
       write(6,*)('could not find (ov|ov) for (oo''|ov'')')
       call dws_sync
       call dws_finalize
       stop
  620  continue
c
c     (oo|o'v'): (oo|oo) part ok
c
       if(min(nocc(isw3),nocc(isw4)).ne.0)then                          2d24s23
       do is=1,nsdlk
        if(isblk(1,is).eq.isblkxder1(1,isb).and.isblk(2,is).eq.
     $       isblkxder1(2,isb).and.isblk(3,is).eq.isw3.and.
     $       isblk(4,is).eq.isw4)then
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d12s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           if(isblk(1,is).eq.isblk(2,is))then
            do i4=0,nocc(isblk(1,is))-1
             do i3=0,i4-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               iad1=itrans(isw4)+i2m+nbasdwsc(isw4)*m4p
               fa=2d0*bc(ii)*bc(iad1)
               do m3=0,nocc(isblkxder1(3,isb))-1
                iad2=itrans(isw3)+i1m+nbasdwsc(isw3)*m3
                iad3=ionexd2(isb)+i3+nocc(isblk(1,is))
     $         *(i4+nocc(isblk(1,is))*(m3+nocc(isblkxder1(3,isb))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
                iad3=ionexd2(isb)+i4+nocc(isblk(1,is))
     $         *(i3+nocc(isblk(2,is))*(m3+nocc(isblkxder1(3,isb))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+i2m+nbasdwsc(isw4)*m4p
              fa=2d0*bc(ii)*bc(iad1)
              do m3=0,nocc(isblkxder1(3,isb))-1
               iad2=itrans(isw3)+i1m+nbasdwsc(isw3)*m3
               iad3=ionexd2(isb)+i4+nocc(isblk(1,is))
     $         *(i4+nocc(isblk(2,is))*(m3+nocc(isblkxder1(3,isb))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           else
            do i4=0,nocc(isblk(2,is))-1
             do i3=0,nocc(isblk(1,is))-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               jad1=itrans(isw4)+i2m+nbasdwsc(isw4)*m4p
               fb=2d0*bc(ii)*bc(jad1)
               do m3=0,nocc(isblkxder1(3,isb))-1
                jad2=itrans(isw3)+i1m+nbasdwsc(isw3)*m3
                iad3=ionexd2(isb)+i3+nocc(isblk(1,is))
     $         *(i4+nocc(isblk(2,is))*(m3+nocc(isblkxder1(3,isb))*m4))
                bc(iad3)=bc(iad3)+fb*bc(jad2)
               end do
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 605
        end if                                                          9d12s16
        if(isblk(2,is).eq.isblkxder1(1,isb).and.isblk(1,is).eq.
     $       isblkxder1(2,isb).and.isblk(3,is).eq.isw3.and.
     $       isblk(4,is).eq.isw4)then
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d12s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           do i4=0,nocc(isblk(2,is))-1
            do i3=0,nocc(isblk(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              jad1=itrans(isw4)+i2m+nbasdwsc(isw4)*m4p
              fb=2d0*bc(ii)*bc(jad1)
              do m3=0,nocc(isblkxder1(3,isb))-1
               jad2=itrans(isw3)+i1m+nbasdwsc(isw3)*m3
               iad3=ionexd2(isb)+i4+nocc(isblk(2,is))
     $         *(i3+nocc(isblk(1,is))*(m3+nocc(isblkxder1(3,isb))*m4))
               bc(iad3)=bc(iad3)+fb*bc(jad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 605
        end if                                                          9d12s16
        if(isblk(2,is).eq.isblkxder1(1,isb).and.isblk(1,is).eq.
     $       isblkxder1(2,isb).and.isblk(4,is).eq.isw3.and.
     $       isblk(3,is).eq.isw4)then
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d12s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           do i4=0,nocc(isblk(2,is))-1
            do i3=0,nocc(isblk(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              jad1=itrans(isw4)+i1m+nbasdwsc(isw4)*m4p
              fb=2d0*bc(ii)*bc(jad1)
              do m3=0,nocc(isblkxder1(3,isb))-1
               jad2=itrans(isw3)+i2m+nbasdwsc(isw3)*m3
               iad3=ionexd2(isb)+i4+nocc(isblk(2,is))
     $         *(i3+nocc(isblk(1,is))*(m3+nocc(isblkxder1(3,isb))*m4))
               bc(iad3)=bc(iad3)+fb*bc(jad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 605
        end if                                                          9d12s16
        if(isblk(1,is).eq.isblkxder1(1,isb).and.isblk(2,is).eq.
     $       isblkxder1(2,isb).and.isblk(4,is).eq.isw3.and.
     $       isblk(3,is).eq.isw4)then
         call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,      9d12s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk(3,is))                                          9d12s16
         ii=ioooo(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2m=i2-1                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d12s16
           do i4=0,nocc(isblk(2,is))-1
            do i3=0,nocc(isblk(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              jad1=itrans(isw4)+i1m+nbasdwsc(isw4)*m4p
              fb=2d0*bc(ii)*bc(jad1)
              do m3=0,nocc(isblkxder1(3,isb))-1
               jad2=itrans(isw3)+i2m+nbasdwsc(isw3)*m3
               iad3=ionexd2(isb)+i3+nocc(isblk(1,is))
     $         *(i4+nocc(isblk(2,is))*(m3+nocc(isblkxder1(3,isb))*m4))
               bc(iad3)=bc(iad3)+fb*bc(jad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 605
        end if                                                          9d12s16
       end do
       write(6,*)('could not find oooo for (oo|o''v'') '),
     $      (isblkxder1(i,isb),i=1,4)
       call dws_sync
       call dws_finalize
       stop
  605  continue
       end if                                                           11d30s16
c
c     (oo|o'v'): (oo|vo) part ok
c
       if(min(nocc(isw4),nvirtc(isw3)).ne.0)then                        2d24s23
       do is=1,nsdlk1
        if(isblk1(1,is).eq.isblkxder1(1,isb).and.                                    9d14s16
     $       isblk1(2,is).eq.isblkxder1(2,isb)                          9d16s16
     $       .and.isblk1(3,is).eq.isw4.and.                             9d16s16
     $       isblk1(4,is).eq.isw3)then                                  9d14s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d14s16
         ii=ionex(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d14s16
           if(isblk1(1,is).eq.isblk1(2,is))then                         9d14s16
            do i4=0,nocc(isblk1(1,is))-1
             do i3=0,i4-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               iad1=itrans(isw4)+i1m+nbasdwsc(isw4)*m4p                 9d16s16
               fa=2d0*bc(ii)*bc(iad1)
               do m3=0,nocc(isblkxder1(3,isb))-1
                iad2=itrans(isw3)+i2p+nbasdwsc(isw3)*m3
                iad3=ionexd2(isb)+i4+nocc(isblk1(1,is))*(i3             9d16s16
     $             +nocc(isblk1(1,is))*(m3+nocc(isblkxder1(3,isb))*m4)) 9d16s16
                bc(iad3)=bc(iad3)+fa*bc(iad2)
                iad3=ionexd2(isb)+i3+nocc(isblk1(1,is))*(i4             9d16s16
     $             +nocc(isblk1(1,is))*(m3+nocc(isblkxder1(3,isb))*m4)) 9d16s16
                bc(iad3)=bc(iad3)+fa*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              jad1=itrans(isw4)+i1m+nbasdwsc(isw4)*m4p
              fb=2d0*bc(ii)*bc(jad1)
              do m3=0,nocc(isblkxder1(3,isb))-1
               iad2=itrans(isw3)+i2p+nbasdwsc(isw3)*m3
               iad3=ionexd2(isb)+i4+nocc(isblk1(1,is))*(i4
     $             +nocc(isblk1(1,is))*(m3+nocc(isblkxder1(3,isb))*m4))
               bc(iad3)=bc(iad3)+fb*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           else
            do i4=0,nocc(isblk1(2,is))-1
             do i3=0,nocc(isblk1(1,is))-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               iad1=itrans(isw4)+i1m+nbasdwsc(isw4)*m4p
               fa=2d0*bc(ii)*bc(iad1)
               do m3=0,nocc(isblkxder1(3,isb))-1
                iad2=itrans(isw3)+i2p+nbasdwsc(isw3)*m3
                iad3=ionexd2(isb)+i3+nocc(isblk1(1,is))*(i4
     $             +nocc(isblk1(2,is))*(m3+nocc(isblkxder1(3,isb))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 621
        end if
        if(isblk1(2,is).eq.isblkxder1(1,isb).and.                                    9d14s16
     $       isblk1(1,is).eq.isblkxder1(2,isb)                          9d16s16
     $       .and.isblk1(3,is).eq.isw4.and.                             9d16s16
     $       isblk1(4,is).eq.isw3)then                                  9d14s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d14s16
         ii=ionex(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d14s16
           do i4=0,nocc(isblk1(2,is))-1
            do i3=0,nocc(isblk1(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+i1m+nbasdwsc(isw4)*m4p
              fa=2d0*bc(ii)*bc(iad1)
              do m3=0,nocc(isblkxder1(3,isb))-1
               iad2=itrans(isw3)+i2p+nbasdwsc(isw3)*m3
               iad3=ionexd2(isb)+i4+nocc(isblk1(2,is))*(i3
     $              +nocc(isblk1(1,is))*(m3+nocc(isblkxder1(3,isb))*m4))
               bc(iad3)=bc(iad3)+fa*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 621
        end if
       end do
       write(6,*)('could not find (oo|vo) for (oo|o''v'')')
       call dws_sync
       call dws_finalize
       stop
  621 continue
      end if                                                            11d30s16
c
c     (oo|o'v'): (oo|ov) part ok
c
      if(min(nocc(isw3),nvirtc(isw4)).ne.0)then                         2d24s23
      do is=1,nsdlk1
        if(isblk1(1,is).eq.isblkxder1(1,isb).and.                                    9d14s16
     $       isblk1(2,is).eq.isblkxder1(2,isb)                          9d16s16
     $       .and.isblk1(3,is).eq.isw3.and.                             9d16s16
     $       isblk1(4,is).eq.isw4)then                                  9d14s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d14s16
         ii=ionex(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d14s16
           if(isblk1(1,is).eq.isblk1(2,is))then                         9d14s16
            do i4=0,nocc(isblk1(1,is))-1
             do i3=0,i4-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               iad1=itrans(isw4)+i2p+nbasdwsc(isw4)*m4p                 9d16s16
               fa=2d0*bc(ii)*bc(iad1)
               do m3=0,nocc(isblkxder1(3,isb))-1
                iad2=itrans(isw3)+i1m+nbasdwsc(isw3)*m3
                iad3=ionexd2(isb)+i4+nocc(isblk1(1,is))*(i3             9d16s16
     $             +nocc(isblk1(1,is))*(m3+nocc(isblkxder1(3,isb))*m4)) 9d16s16
                bc(iad3)=bc(iad3)+fa*bc(iad2)
                iad3=ionexd2(isb)+i3+nocc(isblk1(1,is))*(i4             9d16s16
     $             +nocc(isblk1(1,is))*(m3+nocc(isblkxder1(3,isb))*m4)) 9d16s16
                bc(iad3)=bc(iad3)+fa*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              jad1=itrans(isw4)+i2p+nbasdwsc(isw4)*m4p
              fb=2d0*bc(ii)*bc(jad1)
              do m3=0,nocc(isblkxder1(3,isb))-1
               iad2=itrans(isw3)+i1m+nbasdwsc(isw3)*m3
               iad3=ionexd2(isb)+i4+nocc(isblk1(1,is))*(i4
     $             +nocc(isblk1(1,is))*(m3+nocc(isblkxder1(3,isb))*m4))
               bc(iad3)=bc(iad3)+fb*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           else
            do i4=0,nocc(isblk1(2,is))-1
             do i3=0,nocc(isblk1(1,is))-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               iad1=itrans(isw4)+i2p+nbasdwsc(isw4)*m4p
               fa=2d0*bc(ii)*bc(iad1)
               do m3=0,nocc(isblkxder1(3,isb))-1
                iad2=itrans(isw3)+i1m+nbasdwsc(isw3)*m3
                iad3=ionexd2(isb)+i3+nocc(isblk1(1,is))*(i4
     $             +nocc(isblk1(2,is))*(m3+nocc(isblkxder1(3,isb))*m4))
                bc(iad3)=bc(iad3)+fa*bc(iad2)
               end do
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 622
        end if
        if(isblk1(2,is).eq.isblkxder1(1,isb).and.                                    9d14s16
     $       isblk1(1,is).eq.isblkxder1(2,isb)                          9d16s16
     $       .and.isblk1(3,is).eq.isw3.and.                             9d16s16
     $       isblk1(4,is).eq.isw4)then                                  9d14s16
         call ilimts(nocc(isblk1(3,is)),nvirtc(isblk1(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nocc(isblk1(3,is))                                         9d14s16
         ii=ionex(is)                                                   9d13s16
         do i2=i2s,i2e                                                  9d2s16
          i2p=i2-1+nocc(isblk1(4,is))
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          do i1=i10,i1n                                                 9d2s16
           i1m=i1-1                                                     9d14s16
           do i4=0,nocc(isblk1(2,is))-1
            do i3=0,nocc(isblk1(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+i2p+nbasdwsc(isw4)*m4p
              fa=2d0*bc(ii)*bc(iad1)
              do m3=0,nocc(isblkxder1(3,isb))-1
               iad2=itrans(isw3)+i1m+nbasdwsc(isw3)*m3
               iad3=ionexd2(isb)+i4+nocc(isblk1(2,is))*(i3
     $              +nocc(isblk1(1,is))*(m3+nocc(isblkxder1(3,isb))*m4))
               bc(iad3)=bc(iad3)+fa*bc(iad2)
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 622
        end if
       end do
       write(6,*)('could not find (oo|ov) for (oo|o''v'')')
       call dws_sync
       call dws_finalize
       stop
  622 continue
      end if
c
c     (oo|o'v'): (oo|vv) part ok
c
      if(min(nvirtc(isw3),nvirtc(isw4)).eq.0)go to 623                  2d24s23
       do is=1,nsdlk
        if(isblk(1,is).eq.isblkxder1(1,isb).and.
     $     isblk(2,is).eq.isblkxder1(2,isb).and.
     $     isblk(3,is).eq.isw3.and.
     $     isblk(4,is).eq.isw4)then
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d14s16
         ii=jmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2p=i2-1+nocc(isblk(4,is))                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblk(3,is))
           if(isblk(1,is).eq.isblk(2,is))then
            do i4=0,nocc(isblk(1,is))-1
             do i3=0,i4-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               iad1=itrans(isw4)+i2p+nbasdwsc(isw4)*m4p
               fa=2d0*bc(ii)*bc(iad1)
               do m3=0,nocc(isblkxder1(3,isb))-1
                iad2=itrans(isw3)+i1p+nbasdwsc(isw3)*m3
                iad3=ionexd2(isb)+i4+nocc(isblk(1,is))*(i3
     $        +nocc(isblk(1,is))*(m3+nocc(isblkxder1(3,isb))*m4))
                bc(iad3)=bc(iad3)+bc(iad2)*fa
                iad3=ionexd2(isb)+i3+nocc(isblk(1,is))*(i4
     $        +nocc(isblk(1,is))*(m3+nocc(isblkxder1(3,isb))*m4))
                bc(iad3)=bc(iad3)+bc(iad2)*fa
               end do
              end do
              ii=ii+1
             end do
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+i2p+nbasdwsc(isw4)*m4p
              fa=2d0*bc(ii)*bc(iad1)
              do m3=0,nocc(isblkxder1(3,isb))-1
               iad2=itrans(isw3)+i1p+nbasdwsc(isw3)*m3
               iad3=ionexd2(isb)+i4+nocc(isblk(1,is))*(i4
     $        +nocc(isblk(1,is))*(m3+nocc(isblkxder1(3,isb))*m4))
               bc(iad3)=bc(iad3)+bc(iad2)*fa
              end do
             end do
             ii=ii+1
            end do
           else
            do i4=0,nocc(isblk(2,is))-1
             do i3=0,nocc(isblk(1,is))-1
              do m4=0,nvirtc(isblkxder1(4,isb))-1
               m4p=m4+nocc(isblkxder1(4,isb))
               iad1=itrans(isw4)+i2p+nbasdwsc(isw4)*m4p
               fa=2d0*bc(ii)*bc(iad1)
               do m3=0,nocc(isblkxder1(3,isb))-1
                iad2=itrans(isw3)+i1p+nbasdwsc(isw3)*m3
                iad3=ionexd2(isb)+i3+nocc(isblk(1,is))*(i4
     $        +nocc(isblk(2,is))*(m3+nocc(isblkxder1(3,isb))*m4))
                bc(iad3)=bc(iad3)+bc(iad2)*fa
               end do
              end do
              ii=ii+1
             end do
            end do
           end if
          end do
          i10=1
         end do
         go to 623
        end if
        if(isblk(2,is).eq.isblkxder1(1,isb).and.
     $     isblk(1,is).eq.isblkxder1(2,isb).and.
     $     isblk(3,is).eq.isw3.and.
     $     isblk(4,is).eq.isw4)then
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d14s16
         ii=jmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2p=i2-1+nocc(isblk(4,is))                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblk(3,is))
           do i4=0,nocc(isblk(2,is))-1
            do i3=0,nocc(isblk(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+i2p+nbasdwsc(isw4)*m4p
              fa=2d0*bc(ii)*bc(iad1)
              do m3=0,nocc(isblkxder1(3,isb))-1
               iad2=itrans(isw3)+i1p+nbasdwsc(isw3)*m3
               iad3=ionexd2(isb)+i4+nocc(isblk(2,is))*(i3
     $        +nocc(isblk(1,is))*(m3+nocc(isblkxder1(3,isb))*m4))
               bc(iad3)=bc(iad3)+bc(iad2)*fa
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 623
        end if
        if(isblk(2,is).eq.isblkxder1(1,isb).and.
     $     isblk(1,is).eq.isblkxder1(2,isb).and.
     $     isblk(4,is).eq.isw3.and.
     $     isblk(3,is).eq.isw4)then
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d14s16
         ii=jmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2p=i2-1+nocc(isblk(4,is))                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblk(3,is))
           do i4=0,nocc(isblk(2,is))-1
            do i3=0,nocc(isblk(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+i1p+nbasdwsc(isw4)*m4p
              fa=2d0*bc(ii)*bc(iad1)
              do m3=0,nocc(isblkxder1(3,isb))-1
               iad2=itrans(isw3)+i2p+nbasdwsc(isw3)*m3
               iad3=ionexd2(isb)+i4+nocc(isblk(2,is))*(i3
     $        +nocc(isblk(1,is))*(m3+nocc(isblkxder1(3,isb))*m4))
               bc(iad3)=bc(iad3)+bc(iad2)*fa
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 623
        end if
        if(isblk(1,is).eq.isblkxder1(1,isb).and.
     $     isblk(2,is).eq.isblkxder1(2,isb).and.
     $     isblk(4,is).eq.isw3.and.
     $     isblk(3,is).eq.isw4)then
         call ilimts(nvirtc(isblk(3,is)),nvirtc(isblk(4,is)),mynprocg,  9d14s16
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d24s16
         i10=i1s                                                        8d24s16
         i1n=nvirtc(isblk(3,is))                                        9d14s16
         ii=jmats(is)                                                   9d2s16
         do i2=i2s,i2e                                                  9d2s16
          if(i2.eq.i2e)i1n=i1e                                          9d2s16
          i2p=i2-1+nocc(isblk(4,is))                                                      9d12s16
          do i1=i10,i1n                                                 9d2s16
           i1p=i1-1+nocc(isblk(3,is))
           do i4=0,nocc(isblk(2,is))-1
            do i3=0,nocc(isblk(1,is))-1
             do m4=0,nvirtc(isblkxder1(4,isb))-1
              m4p=m4+nocc(isblkxder1(4,isb))
              iad1=itrans(isw4)+i1p+nbasdwsc(isw4)*m4p
              fa=2d0*bc(ii)*bc(iad1)
              do m3=0,nocc(isblkxder1(3,isb))-1
               iad2=itrans(isw3)+i2p+nbasdwsc(isw3)*m3
               iad3=ionexd2(isb)+i3+nocc(isblk(1,is))*(i4
     $        +nocc(isblk(2,is))*(m3+nocc(isblkxder1(3,isb))*m4))
               bc(iad3)=bc(iad3)+bc(iad2)*fa
              end do
             end do
             ii=ii+1
            end do
           end do
          end do
          i10=1
         end do
         go to 623
        end if
       end do
       write(6,*)('could not find (oo|vv) for (oo|o''v'')')
       call dws_sync
       call dws_finalize
       stop
  623 continue
      if(idwsdeb.gt.10.and.min(nrow,ncol).gt.0)then
       write(6,*)('onexd2 ooox before global sum '),ionexd2(isb)
       write(6,12)(isblkxder1(j,isb),j=1,4)
       call prntm2(bc(ionexd2(isb)),nrow,ncol,nrow)
       itmp=ibcoff
       ibcoff=itmp+nrow*ncol
       call enough('parajkfromhd0. 10',bc,ibc)
       do i=0,nrow*ncol-1
        bc(itmp+i)=bc(ionexd2(isb)+i)
       end do
       call dws_gsumf(bc(itmp),nrow*ncol)
       write(6,*)('after global sum ')
       write(6,12)(isblkxder1(j,isb),j=1,4)
       call prntm2(bc(itmp),nrow,ncol,nrow)
        if(nsymb.eq.1)then
         call printa(bc(itmp),nocc(isblkxder1 (1,isb)),0,
     $        nocc(isblkxder1 (2,isb)),0,nocc(isblkxder1 (3,isb)),0,
     $        nvirt(isblkxder1(4,isb)),nocc(isblkxder1 (4,isb)),
     $        bc(ibcoff))
        end if
       ibcoff=itmp
      end if
      end do
      return
      end
