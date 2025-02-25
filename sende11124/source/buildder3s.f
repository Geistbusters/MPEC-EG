c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine buildder3s(ionex,ider3,h0mo,noc,nvirtc,nbasdwsc,       1d23s17
     $     idwsdeb,bc,ibc)                                              11d14s22
      implicit real*8 (a-h,o-z)
c
c     compute sum over occupied orbitals of term going into 3rd
c     derivative of the energy wrt orbital rotations.
c
      include "common.store"
      include "common.hf"
      dimension noc(8),nvirtc(8),multh(8,8),ionex(1),h0mo(1),
     $     nbasdwsc(8)
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      ider3=ibcoff                                                      1d23s17
      jder3=ider3                                                       1d23s17
      do isb=1,nsymb                                                    1d23s17
       ibcoff=jder3+nvirtc(isb)*noc(isb)                                1d23s17
       call enough('buildder3s.  1',bc,ibc)
       do i=0,nvirtc(isb)*noc(isb)-1                                    1d23s17
        bc(jder3+i)=0d0                                                 1d23s17
       end do                                                           1d23s17
       do is=1,nsdlk1                                                   1d23s17
        if(isblk1(4,is).eq.isb)then                                     1d23s17
c
c     coulomb term
c
         if(isblk1(3,is).eq.isb)then
          call ilimts(noc(isblk1(3,is)),nvirtc(isblk1(4,is)),           1d23s17
     $         mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)
          nrow=(noc(isblk1(1,is))*(noc(isblk1(1,is))+1))/2
          i1n=noc(isblk1(3,is))
          i10=i1s
          ii=ionex(is)-1
          do i2=i2s,i2e
           i2m=i2-1
           if(i2.eq.i2e)i1n=i1e
           do i1=i10,i1n
            i1m=i1-1
            isum=jder3+i2m+nvirtc(isblk1(4,is))*i1m
            do i34=1,noc(isblk1(1,is))
             iad=ii+(i34*(i34+1))/2
             bc(isum)=bc(isum)-2d0*bc(iad)
            end do
            ii=ii+nrow
           end do
           i10=1
          end do
         end if
c
c     exchange term
c
         if(isblk1(1,is).eq.isb)then
          call ilimts(noc(isblk1(3,is)),nvirtc(isblk1(4,is)),           1d23s17
     $         mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)
          if(isblk1(1,is).eq.isblk1(2,is))then                          1d23s17
           nrow=(noc(isblk1(1,is))*(noc(isblk1(1,is))+1))/2
          else
           nrow=noc(isblk1(1,is))*noc(isblk1(2,is))
          end if
          i1n=noc(isblk1(3,is))
          i10=i1s
          ii=ionex(is)
          do i2=i2s,i2e
           i2m=i2-1
           if(i2.eq.i2e)i1n=i1e
           do i1=i10,i1n
            i1m=i1-1
            if(isblk1(1,is).eq.isblk1(2,is))then
             do iq=0,noc(isb)-1
              ix=max(iq,i1m)
              in=min(iq,i1m)
              isum=jder3+i2m+nvirtc(isblk1(4,is))*iq
              iad=ii+((ix*(ix+1))/2)+in
              bc(isum)=bc(isum)+bc(iad)
             end do
            else
             do iq=0,noc(isb)-1
              isum=jder3+i2m+nvirtc(isblk1(4,is))*iq
              iad=ii+iq+noc(isb)*i1m
              bc(isum)=bc(isum)+bc(iad)
             end do
            end if
            ii=ii+nrow
           end do
           i10=1
          end do
         else if(isblk1(2,is).eq.isb)then
          call ilimts(noc(isblk1(3,is)),nvirtc(isblk1(4,is)),           1d23s17
     $         mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)
          nrow=noc(isblk1(1,is))*noc(isblk1(2,is))
          i1n=noc(isblk1(3,is))
          i10=i1s
          ii=ionex(is)
          do i2=i2s,i2e
           i2m=i2-1
           if(i2.eq.i2e)i1n=i1e
           do i1=i10,i1n
            i1m=i1-1
            do iq=0,noc(isb)-1
             isum=jder3+i2m+nvirtc(isblk1(4,is))*iq
             iad=ii+i1m+noc(isblk1(3,is))*iq
             bc(isum)=bc(isum)+bc(iad)
            end do
            ii=ii+nrow
           end do
           i10=1
          end do
         end if
        end if                                                          1d23s17
       end do                                                           1d23s17
       jder3=jder3+nvirtc(isb)*noc(isb)                                 1d23s17
      end do                                                            1d23s17
      nder3=jder3-ider3                                                 1d23s17
      call dws_gsumf(bc(ider3),nder3)                                   1d23s17
      jder3=ider3
      ioff=1
      do isb=1,nsymb
       do iq=0,noc(isb)-1
        do il=0,nvirtc(isb)-1
         ilp=il+noc(isb)
         isum=jder3+il+nvirtc(isb)*iq
         ih0=ioff+ilp+nbasdwsc(isb)*iq
         bc(isum)=bc(isum)-h0mo(ih0)
        end do
       end do
       if(idwsdeb.gt.10)then
        write(6,*)('der3s matrix for symmetry block '),isb
        call prntm2(bc(jder3),nvirtc(isb),noc(isb),nvirtc(isb))
       end if
       ioff=ioff+nbasdwsc(isb)*nbasdwsc(isb)
       jder3=jder3+noc(isb)*nvirtc(isb)
      end do
      return
      end
