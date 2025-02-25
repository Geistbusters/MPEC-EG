c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine genergy(nsymb,ih0mo,idoub,iact,nbasdws,iden1,isblk,    5d17s16
     $     nsdlk,jdenpt,potnuc,                                         5d19s16
     $     numpro,nowpro,ioooo,ehf,ncomp,idwsdeb,enobreit,lprint,iconv, 11d11s22
     $     bc,ibc)                                                      11d11s22
      implicit real*8 (a-h,o-z)
c
c     compute energy using implied densities for doubly occupied orbs   5d17s16
c     and explicit densities for active orbs.                           5d17s16
c                                                                       5d17s16
      logical lprint,lprintx                                            3d3s20
      dimension iact(8),nbasdws(8),iden1(8),isblk(4,nsdlk),
     $     jdenpt(nsdlk),idoub(8),ioooo(nsdlk),noc(8)
      include "common.store"
      include "common.print"                                            3d3s20
      common/gencm/return(8)
      data icall/0/
      data loopx/10000000/
      save icall
      loop=0
      icall=icall+1
      if(iprtr(16).eq.0)then                                            3d3s20
       lprintx=.false.                                                  3d3s20
      else                                                              3d3s20
       lprintx=.true.                                                   3d3s20
      end if                                                            3d3s20
      if(lprintx)write(6,*)('Hi, my name is genergy'),ih0mo,loop             3d3s20
c
c     using shift, modified h0
c
      eshift1=0d0                                                       4d18s17
      eshift2=0d0                                                       4d18s17
      eshiftc=0d0
      jh0mo=ih0mo                                                       4d18s17
      do isb=1,nsymb                                                    4d18s17
       do i=0,idoub(isb)-1                                              4d18s17
        iad=jh0mo+i*(1+nbasdws(isb))                                    4d18s17
        eshift1=eshift1+2d0*bc(iad)                                     4d18s17
       end do                                                           4d18s17
       if(lprintx)then                                                  3d3s20
        write(6,*)('h0mo for symmetry block '),isb,jh0mo,idoub(isb),loop
        call prntm2(bc(jh0mo),nbasdws(isb),nbasdws(isb),nbasdws(isb))
        write(6,*)('eshift1 so far '),eshift1
        write(6,*)('density '),loop
        call prntm2(bc(iden1(isb)),iact(isb),iact(isb),iact(isb))
       end if                                                           3d3s20
       jh0mo=jh0mo+nbasdws(isb)*nbasdws(isb)                            4d18s17
      end do                                                            4d18s17
      if(lprintx)write(6,*)('eshift1 '),eshift1,nsdlk,loop                   3d3s20
      do is=1,nsdlk                                                     4d18s17
       n3=idoub(isblk(3,is))+iact(isblk(3,is))
       n4=idoub(isblk(4,is))+iact(isblk(4,is))
       call ilimts(n3,n4,numpro,nowpro,il,ih,i1s,i1e,i2s,i2e)
       nhere=ih+1-il
       n1=iact(isblk(1,is))+idoub(isblk(1,is))
       if(isblk(1,is).eq.isblk(2,is))then
        nrow=(n1*(n1+1))/2
       else
        n2=iact(isblk(2,is))+idoub(isblk(2,is))
        nrow=n1*n2
       end if
       if(isblk(1,is).eq.isblk(2,is))then                               4d18s17
        is12=isblk(1,is)                                                4d18s17
        is34=isblk(3,is)                                                4d18s17
        noc34=idoub(is34)+iact(is34)                                    4d19s17
        noc12=idoub(is12)+iact(is12)                                    4d19s17
        call ilimts(noc34,noc34,numpro,nowpro,il,ih,i1s,i1e,i2s,i2e)    4d19s17
        i10=i1s                                                         4d18s17
        i1n=noc34                                                       4d19s17
        i4o=ioooo(is)-1                                                 4d18s17
        nrow=(noc12*(noc12+1))/2                                        4d19s17
        if(lprintx)then
         write(6,*)('4o for symmetry '),(isblk(j,is),j=1,4)
         ncolx=(noc34*(noc34+1))/2
         call prntm2(bc(ioooo(is)),nrow,ncolx,nrow)
        end if
        do i2=i2s,i2e                                                   4d18s17
         if(i2.eq.i2e)i1n=i1e                                           4d18s17
         do i1=i10,i1n                                                  4d18s17
          if(i1.eq.i2.and.i1.le.idoub(is34))then                        4d18s17
           do i12=1,idoub(is12)                                         4d18s17
            iad=i4o+((i12*(i12+1))/2)                                   4d18s17
            eshift2=eshift2+2d0*bc(iad)                                 4d18s17
            eshiftc=eshiftc+2d0*bc(iad)
           end do                                                       4d18s17
          end if                                                        4d18s17
          i4o=i4o+nrow                                                  4d18s17
         end do                                                         4d18s17
         i10=1                                                          4d18s17
        end do                                                          4d18s17
       end if                                                           4d18s17
       if(isblk(1,is).eq.isblk(3,is))then                               4d18s17
        is13=isblk(1,is)                                                4d18s17
        is24=isblk(2,is)                                                4d18s17
        noc13=idoub(is13)+iact(is13)                                    4d19s17
        noc24=idoub(is24)+iact(is24)                                    4d19s17
        call ilimts(noc13,noc24,numpro,nowpro,il,ih,i1s,i1e,i2s,i2e)    4d19s17
        i10=i1s                                                         4d18s17
        i1n=noc13                                                       4d19s17
        i4o=ioooo(is)-1                                                 4d18s17
        if(is13.eq.is24)then                                            4d18s17
         nrow=((noc13*(noc13+1))/2)                                     4d19s17
        else                                                            4d18s17
         nrow=noc13*noc24                                               4d19s17
        end if                                                          4d18s17
        if(lprintx.and.nrow.gt.0)then                                   6d12s24
         write(6,*)('4o for symmetry '),(isblk(j,is),j=1,4)
         call prntm2(bc(ioooo(is)),nrow,nrow,nrow)
        end if
        do i2=i2s,i2e                                                   4d18s17
         if(i2.eq.i2e)i1n=i1e                                           4d18s17
         do i1=i10,i1n                                                  4d18s17
          if(i2.le.idoub(is24).and.i1.le.idoub(is13))then               4d18s17
           if(is13.eq.is24)then                                          4d18s17
            ix=max(i1,i2)                                                4d18s17
            in=min(i1,i2)                                                4d18s17
            ix=((ix*(ix-1))/2)+in+i4o                                   4d18s17
            eshift2=eshift2-bc(ix)                                      4d18s17
           else                                                          4d18s17
            ix=i1+noc13*(i2-1)+i4o                                      4d19s17
            eshift2=eshift2-bc(ix)*2d0                                  4d18s17
           end if                                                        4d18s17
          end if                                                        4d18s17
          i4o=i4o+nrow                                                  4d18s17
         end do                                                         4d18s17
         i10=1                                                          4d18s17
        end do                                                          4d18s17
       end if                                                           4d18s17
      end do                                                            4d18s17
      call dws_gsumf(eshift2,1)                                         4d18s17
      call dws_gsumf(eshiftc,1)                                         4d18s17
      shift=eshift1+eshift2                                             4d18s17
      if(lprintx)then                                                   3d3s20
       write(6,*)('one e contribution to shift: '),eshift1               4d18s17
       write(6,*)('two e contribution to shift: '),eshift2
       write(6,*)('Coulomb part: '),eshiftc
       write(6,*)('Exchange part: '),eshift2-eshiftc
       write(6,*)('total shift: '),shift                                 4d18s17
      end if                                                            3d3s20
      jh0mo=ih0mo                                                       4d18s17
      trac1a=0d0                                                        4d18s17
      trac1b=0d0                                                        4d18s17
      tracx=0d0
      do isb=1,nsymb                                                    4d18s17
       nocsb=idoub(isb)+iact(isb)                                       4d19s17
       do im=0,iact(isb)-1                                              4d18s17
        imp=im+idoub(isb)+1                                             4d19s17
        do in=0,iact(isb)-1                                             4d18s17
         inp=in+idoub(isb)+1                                            4d19s17
         iadd=iden1(isb)+in+iact(isb)*im                                4d18s17
         iadh0=jh0mo+inp-1+nbasdws(isb)*(imp-1)                         4d19s17
         trac1a=trac1a+bc(iadd)*bc(iadh0)                               4d18s17
         sum=0d0                                                        4d18s17
         do is=1,nsdlk                                                  4d18s17
          if(isblk(1,is).eq.isblk(2,is).and.isblk(3,is).eq.isb.and.     4d18s17
     $         idoub(isblk(1,is)).gt.0)then                             4d18s17
           is12=isblk(1,is)                                             4d18s17
           noc12=idoub(is12)+iact(is12)                                 4d19s17
           call ilimts(nocsb,nocsb,numpro,nowpro,il,ih,i1s,i1e,i2s,i2e) 4d19s17
           i10=i1s                                                      4d18s17
           i1n=nocsb                                                    4d19s17
           i4o=ioooo(is)-1                                              4d19s17
           nrow=((noc12*(noc12+1))/2)                                   4d19s17
           do i2=i2s,i2e                                                4d18s17
            if(i2.eq.i2e)i1n=i1e                                        4d18s17
            do i1=i10,i1n                                               4d18s17
             if(i1.eq.inp.and.i2.eq.imp)then                            4d18s17
              do i12=1,idoub(is12)                                      4d18s17
               iad=i4o+((i12*(i12+1))/2)                                4d18s17
               sum=sum+2d0*bc(iad)                                      4d18s17
              end do                                                    4d18s17
             end if                                                     4d18s17
             i4o=i4o+nrow                                               4d18s17
            end do                                                      4d18s17
            i10=1                                                       4d18s17
           end do                                                       4d18s17
          end if                                                        4d18s17
          if(isblk(1,is).eq.isblk(3,is).and.isblk(2,is).eq.isb.and.     4d18s17
     $         idoub(isblk(1,is)).gt.0)then                             4d18s17
           is13=isblk(1,is)                                             4d18s17
           noc13=idoub(is13)+iact(is13)                                 4d19s17
           call ilimts(noc13,nocsb,numpro,nowpro,il,ih,i1s,i1e,i2s,i2e) 4d19s17
           i10=i1s                                                      4d18s17
           i1n=noc13                                                    4d19s17
           i4o=ioooo(is)-1                                              4d18s17
           if(is13.eq.isb)then                                          4d18s17
            nrow=((noc13*(noc13+1))/2)                                  4d19s17
           else                                                         4d18s17
            nrow=noc13*nocsb                                            4d19s17
           end if                                                       4d18s17
           do i2=i2s,i2e                                                4d18s17
            if(i2.eq.i2e)i1n=i1e                                        4d18s17
            do i1=i10,i1n                                               4d18s17
             if(i1.le.idoub(is13).and.i2.eq.imp)then                    4d18s17
              if(is13.eq.isb)then                                       4d18s17
               ix=max(i1,inp)                                           4d18s17
               jn=min(i1,inp)                                           4d19s17
               ii=((ix*(ix-1))/2)+jn-1
               iad=i4o+((ix*(ix-1))/2)+jn
               sum=sum-bc(iad)                                          4d18s17
              else                                                      4d18s17
               iad=i4o+i1+noc13*(inp-1)                                 4d19s17
               sum=sum-bc(iad)                                          4d19s17
              end if                                                    4d18s17
             end if                                                     4d18s17
             i4o=i4o+nrow                                               4d18s17
            end do                                                      4d18s17
            i10=1                                                       4d18s17
           end do                                                       4d18s17
          end if                                                        4d18s17
          if(isblk(1,is).eq.isblk(3,is).and.isblk(1,is).eq.isb.and.     4d20s17
     $         idoub(isblk(2,is)).gt.0.and.isblk(2,is).ne.isb)then      4d20s17
           is24=isblk(2,is)                                             4d20s17
           noc24=idoub(is24)+iact(is24)                                 4d20s17
           call ilimts(nocsb,noc24,numpro,nowpro,il,ih,i1s,i1e,i2s,i2e) 4d20s17
           i10=i1s                                                      4d20s17
           i1n=nocsb                                                    4d20s17
           i4o=ioooo(is)-1                                              4d18s17
           nrow=noc24*nocsb                                             4d20s17
           do i2=i2s,i2e                                                4d18s17
            if(i2.eq.i2e)i1n=i1e                                        4d18s17
            do i1=i10,i1n                                               4d18s17
             if(i2.le.idoub(is24).and.i1.eq.imp)then                    4d20s17
              iad=i4o+inp+nocsb*(i2-1)                                  4d20s17
              sum=sum-bc(iad)                                           4d20s17
             end if                                                     4d20s17
             i4o=i4o+nrow                                               4d20s17
            end do                                                      4d20s17
            i10=1                                                       4d20s17
           end do                                                       4d20s17
          end if                                                        4d20s17
         end do                                                         4d18s17
         trac1b=trac1b+bc(iadd)*sum                                     4d18s17
        end do                                                          4d18s17
       end do                                                           4d18s17
       jh0mo=jh0mo+nbasdws(isb)*nbasdws(isb)                            4d18s17
      end do                                                            4d18s17
      call dws_gsumf(trac1b,1)                                          4d18s17
      if(lprintx)                                                       3d3s20
     $     write(6,*)('trace 1prt density '),trac1a,trac1b,trac1a+trac1b3d3s20
      trac1=0d0
      trac1sv=0d0
      jh0mo=ih0mo                                                       8d13s04
      do 264 isb=1,nsymb                                                8d13s04
       noc(isb)=idoub(isb)+iact(isb)                                    5d17s16
       iarg1=noc(isb)
       iarg2=nbasdws(isb)                                               12d14s15
       do i=1,idoub(isb)                                                5d17s16
        i1=jh0mo+(i-1)*(nbasdws(isb)+1)                                 5d17s16
        trac1sv=trac1sv+2d0*bc(i1)                                          5d17s16
       end do                                                           5d17s16
       do i=1,iact(isb)                                                 5d17s16
        ip=i+idoub(isb)                                                 5d17s16
        i1=jh0mo+idoub(isb)-1+(ip-1)*nbasdws(isb)                       5d17s16
        i2=iden1(isb)-1+(i-1)*iact(isb)                                 5d17s16
        do j=1,iact(isb)                                                5d17s16
         trac1=trac1+bc(i1+j)*bc(i2+j)                                  5d17s16
        end do                                                          5d17s16
       end do                                                           5d17s16
       jh0mo=jh0mo+nbasdws(isb)*nbasdws(isb)                            12d14s15
  264 continue                                                          8d13s04
      if(lprintx)then                                                   3d3s20
       write(6,*)('doub contribution to trace: '),trac1sv
       write(6,*)('active contribution to trace: '),trac1
      end if                                                            3d3s20
      tracx=0d0
      traceshift=0d0                                                    7d20s17
      traceh0=0d0                                                       7d20s17
      thp2=0d0                                                          8d8s17
      thp3=0d0                                                          8d8s17
      sum1=0d0
      sum5=0d0
      sum2=0d0
      sum3=0d0
      sumk1=0d0
      sumk2=0d0
      sumk3=0d0
      sumk4=0d0
      saaaa=0d0                                                         9d5s17
      call dws_sync                                                     10d12s04
      do 265 idws=1,nsdlk                                               8d13s04
       trac2=0d0
       k1=isblk(1,idws)
       k2=isblk(2,idws)
       k3=isblk(3,idws)
       k4=isblk(4,idws)
       if(k1.eq.k3.and.k1.eq.k2)then                                    8d24s04
        xmult=1d0                                                       8d24s04
       else                                                             8d24s04
        xmult=2d0                                                       8d24s04
       end if                                                           8d24s04
       n3=noc(k3)                                                       1d5s12
       n4=noc(k4)                                                       1d5s12
       call ilimts(n3,n4,numpro,nowpro,ilj,ihj,i1s,i1e,i2s,i2e)         8d16s04
       if(k1.eq.k2)then                                                 8d16s04
        nocx=(noc(k1)*(noc(k1)+1))/2                                    8d16s04
       else                                                             8d16s04
        nocx=noc(k1)*noc(k2)                                            8d16s04
       end if                                                           8d16s04
       fact1=2d0                                                        3d30s06
       if(k1.ne.k3)fact1=8d0                                            3d30s06
       fact2=4d0
       if(k1.ne.k3)fact2=8d0
       dws=0d0                                                          3d28s06
       dws2=0d0                                                         3d28s06
       dwso=0d0                                                         3d28s06
       idwso=0
       dwso2=0d0                                                        3d28s06
       dwsd=0d0                                                         3d28s06
       dwsd2=0d0                                                        3d28s06
       dwse=0d0                                                         3d28s06
       dwsf2=0d0                                                        3d28s06
       dwse1=0d0
       dwse2=0d0
       dwse3=0d0
       dwse4=0d0
       dum1=0d0
       dum2a=0d0
       dum2b=0d0
       dum3=0d0                                                         6d13s16
       dum2=0d0                                                         6d13s16
       dum3a=0d0
       dum3b=0d0
       dum4=0d0
       dum5=0d0                                                         5d19s16
       if(idoub(k1)*idoub(k2)*idoub(k3)*idoub(k4).gt.0)then             5d17s16
c
c     2e part of shift
c
        if(k1.eq.k2.and.k3.eq.k4)then
         do i34=1,idoub(k4)
          icol=i34+noc(k3)*(i34-1)
          if(icol.ge.ilj.and.icol.le.ihj)then
           do i12=1,idoub(k1)
            iad1=ioooo(idws)+((i12*(i12+1))/2)-1+nocx*(icol-ilj)
            dum1=dum1+2d0*bc(iad1)
           end do
          end if
         end do
        end if
        if(k1.eq.k3.and.k2.eq.k4)then
         if(k1.eq.k2)then
          do i24=1,idoub(k4)
           do i13=1,idoub(k3)
            icol=i13+noc(k3)*(i24-1)
            if(icol.ge.ilj.and.icol.le.ihj)then
             in=min(i13,i24)
             ix=max(i13,i24)
             ii=((ix*(ix-1))/2)+in-1
             iad1=ioooo(idws)+ii+nocx*(icol-ilj)
             dum1=dum1-bc(iad1)
            end if
           end do
          end do
         else
          do i24=1,idoub(k4)
           do i13=1,idoub(k3)
            icol=i13+noc(k3)*(i24-1)
            if(icol.ge.ilj.and.icol.le.ihj)then
             iad1=ioooo(idws)+i13-1+noc(k1)*(i24-1+noc(k2)*(icol-ilj))
             dum1=dum1-bc(iad1)*2d0
            end if
           end do
          end do
         end if
        end if
       end if                                                           5d17s16
       if(jdenpt(idws).gt.0)then
c
c     jdenpt is addressed by triangle if possible on first 2 and        5d17s16
c     last 2 indices                                                    5d17s16
c                                                                       5d17s16
        if(k1.eq.k2)then                                                5d17s16
         nrow=(iact(k1)*(iact(k1)+1))/2                                 5d17s16
        else                                                            5d17s16
         nrow=iact(k1)*iact(k2)                                         5d17s16
        end if                                                          5d17s16
        i10=i1s
        i1n=noc(k3)                                                     1d5s12
        jj=ioooo(idws)                                                  4d10s17
        do i2=i2s,i2e                                                   1d5s12
         i2m=i2-idoub(k4)                                               5d17s16
         if(i2.eq.i2e)i1n=i1e                                           1d5s12
         do i1=i10,i1n
          i1m=i1-idoub(k3)                                              5d17s16
          if(i1.gt.idoub(k3).and.i2.gt.idoub(k4))then                   5d19s16
           if(k3.eq.k4)then                                             5d17s16
            ix=max(i2m,i1m)                                             5d17s16
            in=min(i2m,i1m)                                             5d17s16
            icol=((ix*(ix-1))/2)+in-1                                   5d17s16
           else                                                         5d17s16
            icol=i1m+iact(k3)*(i2m-1)-1                                 5d17s16
           end if                                                       5d17s16
           jd=jdenpt(idws)+nrow*icol                                    5d17s16
           if(k1.eq.k2)then                                             5d17s16
            do i3=0,iact(k1)-1                                          4d10s17
             do i4=0,iact(k1)-1                                         4d10s17
              ix=max(i3,i4)                                             4d10s17
              in=min(i3,i4)                                             4d10s17
              kd=jd+((ix*(ix+1))/2)+in                                  4d10s17
              ix=ix+idoub(k1)                                           4d10s17
              kk=jj+((ix*(ix+1))/2)+in+idoub(k1)                        4d10s17
              dum5=dum5+bc(kd)*bc(kk)                                   5d17s16
              saaaa=saaaa+bc(kd)*bc(kk)
             end do                                                     5d17s16
            end do                                                      5d17s16
           else                                                         5d17s16
            do i3=1,iact(k2)                                            5d17s16
             i3p=i3+idoub(k2)-1                                         5d17s16
             do i4=1,iact(k1)                                           5d17s16
              i4p=i4+idoub(k1)                                          5d17s16
              kk=jj+i4p+noc(k1)*i3p-1                                   4d10s17
              dum5=dum5+bc(jd)*bc(kk)                                   5d17s16
              saaaa=saaaa+bc(jd)*bc(kk)
              jd=jd+1                                                   5d17s16
             end do                                                     5d17s16
            end do                                                      5d17s16
           end if                                                       5d19s16
          end if                                                        5d17s16
          jj=jj+nocx                                                    5d17s16
         end do
         i10=1
        end do
       end if
       if(idoub(k3)*iact(k1).gt.0.and.k1.eq.k2.and.k3.eq.k4)then        5d19s16
        i10=i1s
        i1n=noc(k3)                                                     5d17s16
        jj=ioooo(idws)-1                                                5d17s16
        do i2=i2s,i2e                                                   5d17s16
         if(i2.eq.i2e)i1n=i1e                                           5d17s16
         do i1=i10,i1n                                                  5d17s16
          if(i1.le.idoub(k3).and.i2.eq.i1)then                          5d17s16
           ii=iden1(k1)                                                 5d17s16
           do i3=1,iact(k2)                                             5d17s16
            i3p=i3+idoub(k2)                                            5d17s16
            do i4=1,iact(k1)                                            5d17s16
             i4p=i4+idoub(k1)                                           5d17s16
             ix=max(i3p,i4p)                                            5d17s16
             in=min(i3p,i4p)                                            5d17s16
             kk=jj+((ix*(ix-1))/2)+in                                   5d17s16
             dum3=dum3+bc(ii)*bc(kk)                                    5d19s16
             ii=ii+1                                                    5d17s16
            end do                                                      5d17s16
           end do                                                       5d17s16
          end if                                                        5d17s16
          jj=jj+nocx                                                    5d17s16
         end do                                                         5d17s16
         i10=1                                                          5d17s16
        end do                                                          5d17s16
       end if
       if(idoub(k1)*iact(k3).gt.0.and.k1.eq.k2.and.k3.eq.k4)then        5d19s16
        i10=i1s
        i1n=noc(k3)                                                     5d17s16
        jj=ioooo(idws)-1                                                5d17s16
        do i2=i2s,i2e                                                   5d17s16
         if(i2.eq.i2e)i1n=i1e                                           5d17s16
         do i1=i10,i1n                                                  5d17s16
          if(i1.gt.idoub(k3).and.i2.gt.idoub(k4))then                   5d17s16
           i1m=i1-idoub(k3)                                             5d17s16
           i2m=i2-idoub(k4)                                             5d17s16
           ii=iden1(k3)+i1m-1+iact(k3)*(i2m-1)                          5d19s16
           value=bc(ii)                                                 5d19s16
           do i3=1,idoub(k2)                                             5d17s16
            kk=jj+((i3*(i3+1))/2)                                        5d19s16
            dum3=dum3+value*bc(kk)                                      5d17s16
           end do                                                       5d17s16
          end if                                                        5d17s16
          jj=jj+nocx                                                    5d17s16
         end do                                                         5d17s16
         i10=1                                                          5d17s16
        end do                                                          5d17s16
       end if
c
c     we store both aabb and bbaa ints, but only the abab part and
c     not baab, abba, and baba
c
       t1=dum2
       fmul=1d0                                                         5d19s16
       if(k1.ne.k2)fmul=4d0                                             6d17s16
       if(idoub(k1)*iact(k2).gt.0.and.k1.eq.k3.and.k2.eq.k4)then        5d19s16
        dum2start=dum2
        i10=i1s
        i1n=noc(k3)                                                     5d17s16
        jj=ioooo(idws)-1                                                5d17s16
        do i2=i2s,i2e                                                   5d17s16
         if(i2.eq.i2e)i1n=i1e                                           5d17s16
         do i1=i10,i1n                                                  5d17s16
          if(i1.le.idoub(k3).and.i2.gt.idoub(k4))then                   5d17s16
           i2m=i2-idoub(k4)                                             5d17s16
           ii=iden1(k2)-1+iact(k4)*(i2m-1)                              5d17s16
           if(k1.eq.k2)then                                             5d19s16
            do i3=1,iact(k2)                                             5d17s16
             i3p=i3+idoub(k2)                                            5d17s16
             ix=max(i1,i3p)                                              5d17s16
             in=min(i1,i3p)                                              5d17s16
             kk=jj+((ix*(ix-1))/2)+in                                    5d17s16
             dum2=dum2-fmul*0.25d0*bc(ii+i3)*bc(kk)                      5d19s16
            end do                                                       5d17s16
           else                                                         5d19s16
            do i3=1,iact(k2)                                             5d17s16
             i3p=i3+idoub(k2)                                            5d17s16
             kk=jj+i1+noc(k1)*(i3p-1)                                   5d19s16
             dum2=dum2-fmul*0.25d0*bc(ii+i3)*bc(kk)                      5d19s16
            end do                                                       5d17s16
           end if                                                       5d19s16
          end if                                                        5d17s16
          jj=jj+nocx                                                    5d17s16
         end do                                                         5d17s16
         i10=1                                                          5d17s16
        end do                                                          5d17s16
       end if
       t1=dum2-t1
       t2=dum2
       if(idoub(k1)*iact(k2).gt.0.and.k1.eq.k4.and.k2.eq.k3)then        5d19s16
        dum2start=dum2
        i10=i1s
        i1n=noc(k3)                                                     5d17s16
        jj=ioooo(idws)-1                                                5d17s16
        do i2=i2s,i2e                                                   5d17s16
         if(i2.eq.i2e)i1n=i1e                                           5d17s16
         do i1=i10,i1n                                                  5d17s16
          if(i1.gt.idoub(k3).and.i2.le.idoub(k4))then                   5d17s16
           i1m=i1-idoub(k3)                                             5d17s16
           ii=iden1(k2)-1+iact(k2)*(i1m-1)                              5d17s16
           if(k1.eq.k2)then                                             5d19s16
            do i3=1,iact(k2)                                             5d17s16
             i3p=i3+idoub(k2)                                            5d17s16
             ix=max(i2,i3p)                                              5d17s16
             in=min(i2,i3p)                                              5d17s16
             kk=jj+((ix*(ix-1))/2)+in                                    5d17s16
             dum2=dum2-fmul*0.25d0*bc(ii+i3)*bc(kk)                      5d19s16
            end do                                                       5d17s16
           else                                                         5d19s16
            do i3=1,iact(k2)                                             5d17s16
             i3p=i3+idoub(k2)                                            5d17s16
             kk=jj+i2+noc(k1)*(i3p-1)                                   5d19s16
             dum2=dum2-fmul*0.25d0*bc(ii+i3)*bc(kk)                      5d19s16
            end do                                                       5d17s16
           end if                                                       5d19s16
          end if                                                        5d17s16
          jj=jj+nocx                                                    5d17s16
         end do                                                         5d17s16
         i10=1                                                          5d17s16
        end do                                                          5d17s16
       end if
       t2=dum2-t2
       t3=dum2
       if(idoub(k2)*iact(k1).gt.0.and.k2.eq.k4.and.k1.eq.k3)then        5d19s16
        dum2start=dum2
        i10=i1s
        i1n=noc(k3)                                                     5d17s16
        jj=ioooo(idws)-1                                                5d17s16
        do i2=i2s,i2e                                                   5d17s16
         if(i2.eq.i2e)i1n=i1e                                           5d17s16
         do i1=i10,i1n                                                  5d17s16
          if(i1.gt.idoub(k3).and.i2.le.idoub(k4))then                   5d17s16
           i1m=i1-idoub(k3)                                             5d17s16
           ii=iden1(k1)-1+iact(k3)*(i1m-1)                              5d17s16
           if(k1.eq.k2)then                                             5d19s16
            do i3=1,iact(k1)                                             5d17s16
             i3p=i3+idoub(k1)                                            5d17s16
             ix=max(i2,i3p)                                              5d17s16
             in=min(i2,i3p)                                              5d17s16
             kk=jj+((ix*(ix-1))/2)+in                                    5d17s16
             dum2=dum2-fmul*0.25d0*bc(ii+i3)*bc(kk)                      5d19s16
            end do                                                       5d17s16
           else                                                         5d19s16
            do i3=1,iact(k1)                                             5d17s16
             i3p=i3+idoub(k1)                                            5d17s16
             kk=jj+i3p+noc(k1)*(i2-1)                                   5d19s16
             dum2=dum2-fmul*0.25d0*bc(ii+i3)*bc(kk)                      5d19s16
20203        format(3i5,2es15.7)
            end do                                                       5d17s16
           end if                                                       5d19s16
          end if                                                        5d17s16
          jj=jj+nocx                                                    5d17s16
         end do                                                         5d17s16
         i10=1                                                          5d17s16
        end do                                                          5d17s16
       end if
       t3=dum2-t3
       t4=dum2
       if(idoub(k2)*iact(k1).gt.0.and.k2.eq.k3.and.k1.eq.k4)then        5d20s16
        dum2start=dum2
        i10=i1s
        i1n=noc(k3)                                                     5d17s16
        jj=ioooo(idws)-1                                                5d17s16
        do i2=i2s,i2e                                                   5d17s16
         if(i2.eq.i2e)i1n=i1e                                           5d17s16
         do i1=i10,i1n                                                  5d17s16
          if(i1.le.idoub(k3).and.i2.gt.idoub(k4))then                   5d17s16
           i2m=i2-idoub(k4)                                             5d17s16
           ii=iden1(k1)-1+iact(k4)*(i2m-1)                              5d17s16
           if(k1.eq.k2)then                                             5d19s16
            do i3=1,iact(k1)                                             5d17s16
             i3p=i3+idoub(k1)                                            5d17s16
             ix=max(i1,i3p)                                              5d17s16
             in=min(i1,i3p)                                              5d17s16
             kk=jj+((ix*(ix-1))/2)+in                                    5d17s16
             dum2=dum2-fmul*0.25d0*bc(ii+i3)*bc(kk)                      5d19s16
            end do                                                       5d17s16
           else                                                         5d19s16
            do i3=1,iact(k1)                                             5d17s16
             i3p=i3+idoub(k1)                                            5d17s16
             kk=jj+i3p+noc(k1)*(i1-1)                                   5d19s16
             dum2=dum2-fmul*0.25d0*bc(ii+i3)*bc(kk)                      5d19s16
            end do                                                      5d19s16
           end if                                                       5d19s16
          end if                                                        5d17s16
          jj=jj+nocx                                                    5d17s16
         end do                                                         5d17s16
         i10=1                                                          5d17s16
        end do                                                          5d17s16
       end if
       trac2=dum1*2d0+dum2*2d0+dum3*2d0+dum4+2d0*dum5+dum2b-dum2a+dum3b     5d19s16
     $      -dum3a                                                      5d19s16
       traceshift=traceshift+dum1                                       7d20s17
       traceh0=traceh0+dum2+dum3                                        7d20s17
       thp2=thp2+dum2                                                   8d8s17
       thp3=thp3+dum3                                                   8d8s17
       sum1=sum1+dum1
       sum5=sum5+dum5
       sum2=sum2+dum2
       sum3=sum3+dum3
       tracx=tracx+trac2
  265 continue                                                          8d13s04
      trac2=tracx
      call dws_sync                                                     10d12s04
      iarg1=1                                                           1d25s05
      call dws_gsumf(sum1,iarg1)
      shift=trac1sv+sum1
      call dws_gsumf(trac2,iarg1)                                        4d12s07
      call dws_gsumf(traceshift,1)                                      7d20s17
      call dws_gsumf(traceh0,1)                                         7d20s17
      call dws_gsumf(thp2,1)                                            8d8s17
      call dws_gsumf(thp3,1)                                            8d8s17
      call dws_gsumf(saaaa,1)
      return(1)=trac1sv                                                 1d5s18
      return(2)=trac1                                                   1d5s18
      return(3)=traceshift                                              1d5s18
      return(4)=traceh0                                                 1d5s18
      return(5)=saaaa                                                   1d5s18
      stry=saaaa+traceh0+traceshift
      return(6)=thp2                                                    1d5s18
      return(7)=thp3                                                    1d5s18
      trac1=trac1+trac1sv                                               5d19s16
      ehf=trac1+0.5d0*trac2+potnuc                                      2d3s17
      if(lprintx)write(6,*)('ehf: '),ehf,trac1,trac2,potnuc             7d5s24
      return(8)=ehf                                                     1d5s18
      if((lprint.and.iconv.eq.1))                                       4d20s18
     $     write(6,*)('>energy from densities '),ehf,trac1,
     $     trac2*0.5d0
      if(idwsdeb.gt.10)then
       write(6,*)('energy parts '),trac1,trac2,trac2*0.5d0,potnuc
      end if
      if(ehf.ne.ehf)then                                                5d7s18
       write(6,*)('in genergy, ehf is nan')                             5d7s18
       write(6,*)('energy parts '),trac1,trac2,trac2*0.5d0
       call dws_sync                                                    5d7s18
       call dws_finalize                                                5d7s18
       stop                                                             5d7s18
      end if                                                            5d7s18
      if(enobreit.ne.0d0.and.lprint)then                                2d2s16
       write(6,*)('>breit correction '),ehf-enobreit                     10d19s15
      end if                                                            10d19s15
      return
      end
