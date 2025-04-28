c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine core(ih0mo,shift,ioooo,noc4,ih0e,shift1,iooooa,ncomp,  8d27s19
     $     idwsdeb,ldeb,bc,ibc)                                         11d14s22
      implicit real*8 (a-h,o-z)
      integer*8 iarg1,iarg2
c
c     take care of doubly occupied orbitals in cas
c     space.
c
      dimension ih0mo(8)                                                8d8s14
      logical ldeb
      include "common.store"
      include "common.hf"
      include "common.cas"
      include "common.print"                                            2d14s20
      dimension noc4(8),ih0e(8),ioooo(1),iooooa(1)                      8d8s14
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
c
c     diagonal shift due to h0
c
      if(iprtr(8).ne.0.or.ldeb)write(6,*)('in core '),mynowprog                  8d27s19
      shift=0d0
      jh0mo=ih0mo(1)                                                    8d8s06
      ih0e(1)=ibcoff                                                    8d8s06
      jh0e=ih0e(1)                                                      8d8s06
      do isb=1,nsymb
       nbasdwsb=nbasdws(isb)                                            2d1s19
       if(iprtr(8).ne.0.or.ldeb)then                                             8d27s19
        write(6,*)('for symmetry block '),isb,jh0mo,loc(bc(jh0mo))
        call prntm2(bc(jh0mo),nbasdwsb,nbasdwsb,nbasdwsb)                8d31s15
       end if                                                           8d27s19
       ih0e(isb)=jh0e                                                   8d8s06
       do i=1,idoub(isb)
        iadd=jh0mo+(i-1)*(nbasdwsb+1)                                   8d31s15
        shift=shift+2d0*bc(iadd)
       end do
       if(mynowprog.eq.0)then                                            8d8s14
        do i2=idoub(isb)+1,noc4(isb)
         do i1=idoub(isb)+1,noc4(isb)
         iadd=jh0mo+i1-1+(i2-1)*nbasdwsb                                8d31s15
         bc(jh0e)=bc(iadd)
          jh0e=jh0e+1
         end do
        end do
       else                                                             8d8s14
        do i2=idoub(isb)+1,noc4(isb)                                    8d8s14
         do i1=idoub(isb)+1,noc4(isb)                                   8d8s14
          bc(jh0e)=0d0                                                  8d8s14
          jh0e=jh0e+1                                                   8d8s14
         end do                                                         8d8s14
        end do                                                          8d8s14
       end if                                                           8d8s14
       na=noc4(isb)-idoub(isb)
       if(iprtr(8).ne.0.or.ldeb)then                                             8d27s19
        write(6,*)('printing ih0e '),noc4(isb),idoub(isb),na
        call prntm2(bc(ih0e(isb)),na,na,na)                              8d8s14
       end if                                                           8d27s19
       jh0mo=jh0mo+nbasdwsb*nbasdwsb                                    8d31s15
      end do
      ibcoff=jh0e
      call enough('core.  1',bc,ibc)
      if(iprtr(8).ne.0.or.ldeb)
     $                write(6,*)('1e part of diagonal shift '),shift    8d27s19
      shift1=shift                                                      11d8s05
c
c     2e integrals
c     recall oooo ints are distributed, so need to do global sum at end 8d8s14
c
      shift2=0d0
      do idws=1,nsdlk
       n1=isblk(1,idws)
       n2=isblk(2,idws)
       n3=isblk(3,idws)
       n4=isblk(4,idws)
   51  format('integral type ',4i2)
       m1=idoub(n1)+iacto(n1)
       m2=idoub(n2)+iacto(n2)
       m3=idoub(n3)+iacto(n3)
       m4=idoub(n4)+iacto(n4)
       if(n1.eq.n2)then
        nn=(m1*(m1+1))/2
       else
        nn=m1*m2
       end if
       if(n3.eq.n4)then
        mm=(m3*(m3+1))/2
       else
        mm=m3*m4
       end if
       call ilimts(m3,m4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)      8d8s14
       nhere=ih+1-il                                                    8d8s14
       if(nn*nhere.gt.0)then                                            8d8s14
 1011  format('try integral type ',4i1)
c
c     2e part of shift
c
        if(n1.eq.n2.and.n3.eq.n4)then
         do i34=1,idoub(n4)
          icol=i34+m3*(i34-1)
          if(icol.ge.il.and.icol.le.ih)then                              8d8s14
           do i12=1,idoub(n1)
            iad1=ioooo(idws)+((i12*(i12+1))/2)-1+nn*(icol-il)             8d8s14
            shift2=shift2+2d0*bc(iad1)
           end do
          end if                                                         8d8s14
         end do
        end if
        if(n1.eq.n3.and.n2.eq.n4)then
         if(n1.eq.n2)then
          do i24=1,idoub(n4)
           do i13=1,idoub(n3)
            icol=i13+m3*(i24-1)
            if(icol.ge.il.and.icol.le.ih)then
             in=min(i13,i24)
             ix=max(i13,i24)
             ii=((ix*(ix-1))/2)+in-1
             iad1=ioooo(idws)+ii+nn*(icol-il)                            8d8s14
             shift2=shift2-bc(iad1)
            end if                                                       8d8s14
           end do
          end do
         else
          do i24=1,idoub(n4)
           do i13=1,idoub(n3)
            icol=i13+m3*(i24-1)                                          8d8s14
            if(icol.ge.il.and.icol.le.ih)then                            8d8s14
             iad1=ioooo(idws)+i13-1+m1*(i24-1+m2*(icol-il))              8d8s14
             shift2=shift2-bc(iad1)*2d0
            end if                                                       8d8s14
           end do
          end do
         end if
        end if
c
c     fold doubly occupied orbitals into effective h0
c
        if(n1.eq.n2.and.n3.eq.n4)then
         do i4=1,iacto(n4)
          i4p=i4+idoub(n4)
          do i3=1,iacto(n3)
           i3p=i3+idoub(n3)
           iad2=ih0e(n4)+i3-1+iacto(n3)*(i4-1)                           8d8s06
           icol=i3p+m3*(i4p-1)                                           8d8s14
           if(icol.ge.il.and.icol.le.ih)then                             8d8s14
            do i12=1,idoub(n1)
             iad1=ioooo(idws)+((i12*(i12+1))/2)-1+nn*(icol-il)           8d8s14
             bc(iad2)=bc(iad2)+bc(iad1)*2d0                               8d8s06
            end do
           end if                                                        8d8s14
          end do
         end do
        end if
        if(n1.eq.n3.and.n2.eq.n4)then
         if(n1.eq.n2)then
          do i4=1,iacto(n4)
           i4p=i4+idoub(n4)
           do i2=1,iacto(n2)
            i2p=i2+idoub(n2)
            iad2=ih0e(n4)+i2-1+iacto(n2)*(i4-1)                           8d8s06
            do i13=1,idoub(n1)
             icol=i13+m3*(i4p-1)                                         8d8s14
             if(icol.ge.il.and.icol.le.ih)then                           8d8s14
              in=min(i13,i2p)
              ix=max(i13,i2p)
              ii=((ix*(ix-1))/2)+in-1
              iad1=ioooo(idws)+ii+nn*(icol-il)                           8d8s14
              bc(iad2)=bc(iad2)-bc(iad1)                                  8d8s06
             end if                                                      8d8s14
            end do
           end do
          end do
         else
          do i4=1,iacto(n4)
           i4p=i4+idoub(n4)
           do i2=1,iacto(n2)
            i2p=i2+idoub(n2)
            iad2=ih0e(n4)+i2-1+iacto(n2)*(i4-1)                           8d8s06
            do i13=1,idoub(n1)
             icol=i13+m3*(i4p-1)                                         8d8s14
             if(icol.ge.il.and.icol.le.ih)then                           8d8s14
              iad1=ioooo(idws)+i13-1+m1*(i2p-1+m2*(icol-il))             8d8s14
              bc(iad2)=bc(iad2)-bc(iad1)                                  8d8s06
             end if                                                      8d8s14
            end do
           end do
          end do
          do i24=1,idoub(n4)                                             8d10s06
           do i3=1,iacto(n3)
            i3p=i3+idoub(n3)                                             8d10s06
            icol=i3p+m3*(i24-1)                                          8d8s14
            if(icol.ge.il.and.icol.le.ih)then                            8d8s14
             do i1=1,iacto(n1)                                            8d10s06
              i1p=i1+idoub(n1)                                            8d10s06
              iad2=ih0e(n3)+i1-1+iacto(n1)*(i3-1)                         8d10s06
              iad1=ioooo(idws)+i1p-1+m1*(i24-1+m2*(icol-il))             8d8s14
              bc(iad2)=bc(iad2)-bc(iad1)                                  8d8s06
             end do
            end if                                                       8d8s14
           end do
          end do
         end if
        end if
       end if                                                           8d11s14
c
c     throw away all 2e ints except for all active
c     and replicate ints on all procs
c
       if(n1.eq.n2)then
        nnn=(iacto(n1)*(iacto(n1)+1))/2
       else
        nnn=iacto(n1)*iacto(n2)
       end if
       if(n3.eq.n4)then
        mmm=(iacto(n3)*(iacto(n3)+1))/2
       else
        mmm=iacto(n3)*iacto(n4)
       end if
       icpy=ibcoff
       ibcoff=icpy+nnn*mmm
       call enough('core.  2',bc,ibc)
       do i=0,nnn*mmm-1                                                 8d8s14
        bc(icpy+i)=0d0                                                  8d8s14
       end do                                                           8d8s14
       jcpy=icpy
       if(n1.eq.n2.and.n3.eq.n4)then
        do i4=1,iacto(n4)
         i4p=i4+idoub(n4)
         do i3=1,i4
          i3p=i3+idoub(n3)
          icol=i3p+m3*(i4p-1)                                           8d8s14
          if(icol.ge.il.and.icol.le.ih)then                             8d8s14
           ii=nn*(icol-il)+ioooo(idws)-1+idoub(n1)                       8d8s14
           do i2=1,iacto(n2)
            i2p=i2+idoub(n2)
            iii=ii+((i2p*(i2p-1))/2)
            do i1=1,i2
             bc(jcpy)=bc(iii+i1)
             jcpy=jcpy+1
            end do
           end do
          else                                                          8d8s14
           jcpy=jcpy+(iacto(n2)*(iacto(n2)+1))/2                        8d8s14
          end if                                                        8d8s14
         end do
        end do
       else
        do i4=1,iacto(n4)
         i4p=i4+idoub(n4)
         do i3=1,iacto(n3)
          i3p=i3+idoub(n3)
          icol=i3p+m3*(i4p-1)                                           8d8s14
          if(icol.ge.il.and.icol.le.ih)then                             8d8s14
           ii=nn*(icol-il)+ioooo(idws)-1+idoub(n1)                      8d8s14
           do i2=1,iacto(n2)
            i2p=i2+idoub(n2)
            iii=ii+m1*(i2p-1)
            do i1=1,iacto(n1)
             bc(jcpy)=bc(iii+i1)
             jcpy=jcpy+1
            end do
           end do
          else                                                          8d8s14
           jcpy=jcpy+iacto(n1)*iacto(n2)                                8d8s14
          end if                                                        8d8s14
         end do
        end do
       end if
       iooooa(idws)=icpy                                                8d8s14
       call dws_gsumf(bc(icpy),nnn*mmm)                                 8d8s14
       if(iprtr(8).ne.0.or.ldeb)then                                             8d27s19
        write(6,*)('final 2e ints '),idws,n1,n2,n3,n4
        call prntm2(bc(icpy),nnn,mmm,nnn)
       end if                                                           8d27s19
      end do
      call dws_gsumf(shift2,1)                                          8d8s14
      do isb=1,nsymb
       nn=iacto(isb)**2                                                 8d8s14
       call dws_gsumf(bc(ih0e(isb)),nn)                                 8d8s14
      end do
      if(iprtr(8).ne.0.or.ldeb)then                                              8d27s19
       write(6,*)('final h0')
       do i=1,nsymb
        if(iacto(i).gt.0)then
         write(6,*)('for symmetry block '),i
        iarg1=iacto(i)
        call prntm2(bc(ih0e(i)),iarg1,iarg1,iarg1)
        if(nsymb.eq.1)call printa(bc(ih0e(i)),iacto(i),idoub(i),1,0,
     $       iacto(i),idoub(i),1,0,bc(ibcoff))
        end if
       end do
      end if                                                            8d27s19
      if(iprtr(8).ne.0.or.ldeb)write(6,*)
     $                          ('2e part of diagonal shift '),shift2   8d27s19
      shift=shift+shift2
      return
      end
