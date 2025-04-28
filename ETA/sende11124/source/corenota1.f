c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine corenota1(ipuse,ih0mo,shift,ioooo,noc4,ih0e,shift1,    6d10s22
     $     iooooa,isblkder,nsblkder,multh,idwsdeb,ldeb,bc,ibc)          11d14s22
      implicit real*8 (a-h,o-z)
      integer*8 iarg1,iarg2
c
c     take care of doubly occupied orbitals in cas
c     space. This version is for perturbations that are not totally     6d10s22
c     symmetric.                                                        6d10s22
c
      dimension ih0mo(*)                                                8d8s14
      logical ldeb                                                      6d22s22
      include "common.store"
      include "common.hf"
      include "common.cas"
      include "common.print"                                            2d14s20
      dimension noc4(*),ih0e(*),ioooo(*),iooooa(*),isblkder(4,*),       6d10s22
     $     multh(8,8)
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
c
c     there is no diagonal shift                                        6d10s22
c
      if(iprtr(8).ne.0)write(6,*)('in corenota1 '),mynowprog                  8d27s19
      shift=0d0
      jh0mo=ih0mo(1)                                                    8d8s06
      ih0e(1)=ibcoff                                                    8d8s06
      jh0e=ih0e(1)                                                      8d8s06
      do isb=1,nsymb
       jsb=multh(isb,ipuse)                                             6d10s22
       if(iprtr(8).ne.0.or.ldeb)then                                    6d22s22
        write(6,*)('for symmetry block '),isb,jsb,jh0mo                 6d10s22
        call prntm2(bc(jh0mo),nbasdws(isb),nbasdws(jsb),nbasdws(isb))   6d10s22
       end if                                                           8d27s19
       ih0e(isb)=jh0e                                                   8d8s06
       if(mynowprog.eq.0)then                                            8d8s14
        do i2=idoub(jsb)+1,noc4(jsb)                                    6d10s22
         do i1=idoub(isb)+1,noc4(isb)
          iadd=jh0mo+i1-1+(i2-1)*nbasdws(isb)                            6d10s22
          bc(jh0e)=bc(iadd)
          jh0e=jh0e+1
         end do
        end do
       else                                                             8d8s14
        do i2=idoub(jsb)+1,noc4(jsb)                                    6d10s22
         do i1=idoub(isb)+1,noc4(isb)                                   8d8s14
          bc(jh0e)=0d0                                                  8d8s14
          jh0e=jh0e+1                                                   8d8s14
         end do                                                         8d8s14
        end do                                                          8d8s14
       end if                                                           8d8s14
       na=noc4(isb)-idoub(isb)
       naj=noc4(jsb)-idoub(jsb)                                         6d10s22
       if(iprtr(8).ne.0.or.ldeb)then                                             8d27s19
        write(6,*)('printing ih0e ')                                    6d10s22
        call prntm2(bc(ih0e(isb)),na,naj,na)                            6d10s22
       end if                                                           8d27s19
       jh0mo=jh0mo+nbasdws(isb)*nbasdws(jsb)                            6d10s22
      end do
      ibcoff=jh0e
      call enough('corenota1.  1',bc,ibc)
      shift1=shift                                                      11d8s05
c
c     2e integrals
c     recall oooo ints are distributed, so need to do global sum at end 8d8s14
c
      shift2=0d0
      do idws=1,nsblkder                                                6d10s22
       n1=isblkder(1,idws)                                              6d10s22
       n2=isblkder(2,idws)                                              6d10s22
       n3=isblkder(3,idws)                                              6d10s22
       n4=isblkder(4,idws)                                              6d10s22
   51  format('integral type ',4i2)
       m1=idoub(n1)+iacto(n1)
       m2=idoub(n2)+iacto(n2)
       m3=idoub(n3)+iacto(n3)
       m4=idoub(n4)+iacto(n4)
       nn=m1*m2                                                         6d10s22
       mm=m3*m4                                                         6d10s22
       call ilimts(m3,m4,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)      8d8s14
       nhere=ih+1-il                                                    8d8s14
       if(nn*nhere.gt.0)then                                            8d8s14
 1011  format('try integral type ',4i1)
c
c     fold doubly occupied orbitals into effective h0
c
        if(n1.eq.n2)then                                                6d10s22
         do i4=1,iacto(n4)
          i4p=i4+idoub(n4)
          do i3=1,iacto(n3)
           i3p=i3+idoub(n3)
           iad2=ih0e(n3)+i3-1+iacto(n3)*(i4-1)                           8d8s06
           iad3=ih0e(n4)+i4-1+iacto(n4)*(i3-1)                          6d10s22
           icol=i3p+m3*(i4p-1)                                           8d8s14
           if(icol.ge.il.and.icol.le.ih)then                             8d8s14
            do i12=0,idoub(n1)-1                                        6d10s22
             iad1=ioooo(idws)+i12*(m1+1)+nn*(icol-il)                   6d10s22
             bc(iad2)=bc(iad2)+bc(iad1)*2d0                               8d8s06
             bc(iad3)=bc(iad3)+bc(iad1)*2d0                             6d10s22
            end do
           end if                                                        8d8s14
          end do
         end do
        end if
        if(n1.eq.n3)then                                                6d10s22
         do i4=1,iacto(n4)
          i4p=i4+idoub(n4)
          do i2=1,iacto(n2)
           i2p=i2+idoub(n2)
           iad2=ih0e(n2)+i2-1+iacto(n2)*(i4-1)                           8d8s06
           do i13=1,idoub(n1)
            icol=i13+m3*(i4p-1)                                         8d8s14
            if(icol.ge.il.and.icol.le.ih)then                           8d8s14
             irow=i13-1+noc4(n1)*(i2p-1)                                6d10s22
             iad1=ioooo(idws)+irow+nn*(icol-il)                         6d10s22
             bc(iad2)=bc(iad2)-bc(iad1)                                  8d8s06
            end if                                                      8d8s14
           end do
          end do
         end do
        end if                                                          6d10s22
        if(n2.eq.n4)then                                                6d10s22
         do i24=1,idoub(n4)                                             8d10s06
          do i3=1,iacto(n3)
           i3p=i3+idoub(n3)                                             8d10s06
           icol=i3p+m3*(i24-1)                                          8d8s14
           if(icol.ge.il.and.icol.le.ih)then                            8d8s14
            do i1=1,iacto(n1)                                            8d10s06
             i1p=i1+idoub(n1)                                            8d10s06
             iad2=ih0e(n1)+i1-1+iacto(n1)*(i3-1)                        6d10s22
             iad1=ioooo(idws)+i1p-1+m1*(i24-1+m2*(icol-il))             8d8s14
             bc(iad2)=bc(iad2)-bc(iad1)                                  8d8s06
            end do
           end if                                                       8d8s14
          end do
         end do
        end if
       end if                                                           8d11s14
c
c     throw away all 2e ints except for all active
c     and replicate ints on all procs
c
       nnn=iacto(n1)*iacto(n2)
       mmm=iacto(n3)*iacto(n4)
       icpy=ibcoff
       ibcoff=icpy+nnn*mmm
       call enough('corenota1.  2',bc,ibc)
       do i=0,nnn*mmm-1                                                 8d8s14
        bc(icpy+i)=0d0                                                  8d8s14
       end do                                                           8d8s14
       jcpy=icpy
       do i4=1,iacto(n4)
        i4p=i4+idoub(n4)
        do i3=1,iacto(n3)
         i3p=i3+idoub(n3)
         icol=i3p+m3*(i4p-1)                                            8d8s14
         if(icol.ge.il.and.icol.le.ih)then                              8d8s14
          ii=nn*(icol-il)+ioooo(idws)-1+idoub(n1)                       8d8s14
          do i2=1,iacto(n2)
           i2p=i2+idoub(n2)
           iii=ii+m1*(i2p-1)
           do i1=1,iacto(n1)
            bc(jcpy)=bc(iii+i1)
            jcpy=jcpy+1
           end do
          end do
         else                                                           8d8s14
          jcpy=jcpy+iacto(n1)*iacto(n2)                                 8d8s14
         end if                                                         8d8s14
        end do
       end do
       iooooa(idws)=icpy                                                8d8s14
       call dws_gsumf(bc(icpy),nnn*mmm)                                 8d8s14
       if(iprtr(8).ne.0.or.ldeb)then                                             8d27s19
        write(6,*)('final 2e ints '),idws,n1,n2,n3,n4,icpy
        call prntm2(bc(icpy),nnn,mmm,nnn)
       end if                                                           8d27s19
      end do
      call dws_gsumf(shift2,1)                                          8d8s14
      do isb=1,nsymb
       jsb=multh(isb,ipuse)                                             6d10s22
       nn=iacto(isb)*iacto(jsb)                                         6d10s22
       call dws_gsumf(bc(ih0e(isb)),nn)                                 8d8s14
      end do
      if(iprtr(8).ne.0.or.ldeb)then                                              8d27s19
       write(6,*)('final h0'),nsymb,ipuse,iacto(1),iacto(2)
       do i=1,nsymb
        j=multh(i,ipuse)                                                6d10s22
         write(6,*)('for symmetry block '),i,j
        if(min(iacto(j),iacto(i)).gt.0)then                             6d10s22
         call prntm2(bc(ih0e(i)),iacto(i),iacto(j),iacto(i))            6d10s22
        end if
       end do
      end if                                                            8d27s19
      if(iprtr(8).ne.0.or.ldeb)
     $     write(6,*)('2e part of diagonal shift '),shift2              6d22s22
      shift=shift+shift2
      return
      end
