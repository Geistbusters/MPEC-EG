c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine paramp2(kmats,noc,ieigr,nowpro,ehf,numpro,icore,       8d9s24
     $     nbasisp,bc,ibc)                                              8d9s24
      implicit real*8 (a-h,o-z)
      external second
c
c      closed shell mp2
c
      logical lprint                                                    10d13s04
      dimension kmats(*),noc(*),ifeig(8),icore(8),sums(2),telap(2),     8d9s24
     $     nbasisp(*)                                                   8d9s24
      include "common.store"
      include "common.hf"
      include "common.print"                                            8d9s24
      idwsdeb=0
      if(iprtr(33).ne.0)idwsdeb=10                                      8d9s24
      lprint=nowpro.eq.0                                                3d19s12
      if(lprint)then
       writE(6,*)('in para mp2 ')
       write(6,*)('reference energy: '),ehf
       write(6,*)('icore = '),(icore(i),i=1,nsymb)
       write(6,*)('noc = '),(noc(i),i=1,nsymb)
       write(6,*)('nbas = '),(nbasdws(i),i=1,nsymb)
       write(6,*)('nbasisp = '),(nbasisp(i),i=1,nsymb)                  8d9s24
      end if
      call second(time1)
      ifeig(1)=ieigr-1
      if(idwsdeb.ne.0)then                                              8d9s24
       write(6,*)('fock eigenvalues for symmetry 1'),ieigr
       call prntm2(bc(ifeig(1)+1),1,nbasdws(1),1)
      end if                                                            8d9s24
      do i=2,nsymb
       ifeig(i)=ifeig(i-1)+nbasisp(i-1)                                 8d9s24
       if(idwsdeb.ne.0)then                                             8d9s24
        write(6,*)('fock eigenvalues for symmetry'),i,ifeig(i)+1
        call prntm2(bc(ifeig(i)+1),1,nbasdws(i),1)
       end if                                                           8d9s24
      end do
      emp2=0d0
      emp2v=0d0
      sp=0d0                                                            10d12s04
      tp=0d0                                                            10d12s04
      lxxx=0                                                            9d10s07
      do idws=1,nsdlkk
       n1=noc(isblkk(1,idws))*noc(isblkk(2,idws))
       n2=nvirt(isblkk(3,idws))*nvirt(isblkk(4,idws))
       if(n1*n2.gt.0)then
        if(lprint)write(6,1351)(isblkk(j,idws),j=1,4)
 1351   format('mp2 contribution from integral type :',4i3)
        m1=isblkk(1,idws)
        m2=isblkk(2,idws)
        m3=isblkk(3,idws)
        m4=isblkk(4,idws)
        n3=nvirt(m3)
        n4=nvirt(m4)
        call ilimts(n3,n4,
     $        numpro,nowpro,ilj,ihj,i1s,i1e,i2s,i2e)
        nhere=ihj+1-ilj
        if(idebug.ne.0)then                                             8d9s24
         call prntm2(bc(kmats(idws)),n1,nhere,n1)                       8d9s24
        end if                                                          8d9s24
        if(m1.eq.m2)then
         if(idwsdeb.gt.10)then
          write(6,*)('ij in same symmetry block')
         end if
         call ilimts(n3,n4,
     $        numpro,nowpro,ilj,ihj,i1s,i1e,i2s,i2e)
         ii=0
         i10=i1s                                                        3d19s12
         i1n=n3                                                         3d19s12
         do i2=i2s,i2e
          if(i2.eq.i2e)i1n=i1e                                          3d19s12
          i2p=i2+noc(m4)                                                3d19s12
          do i1=i10,i1n
           i1p=i1+noc(m3)                                               3d19s12
           do i=icore(m1)+1,noc(m1)
            do j=icore(m1)+1,noc(m1)
             denom=1d0/(bc(ifeig(m1)+i)+bc(ifeig(m1)+j)
     $             -bc(ifeig(m3)+i1p)-bc(ifeig(m4)+i2p))
             k1=kmats(idws)+j-1+noc(m1)*(i-1+noc(m1)*ii)
             k2=kmats(idws)+i-1+noc(m1)*(j-1+noc(m1)*ii)
             emp2v=emp2v+bc(k1)*(2d0*bc(k1)-bc(k2))*denom
             if(i.eq.3.and.j.eq.3)then
              vec=(2d0*bc(k1)-bc(k2))*denom
             end if
 3315        format(5i5,8es15.7)
             tmp=bc(k1)*denom
             sp=sp+tmp*bc(k1)                                           10d12s04
             tp=tp+tmp*bc(k2)                                           10d12s04
            end do
           end do
           ii=ii+1
          end do
          i10=1                                                         3d22s12
         end do
        else
         if(idwsdeb.gt.10)then
          write(6,*)('ij not same symmetry: '),m1,m2
          write(6,*)('look for matching pair ')
         end if
         do idws2=idws+1,nsdlkk
          if(isblkk(1,idws2).eq.m2.and.isblkk(2,idws2).eq.m1.and.
     $       isblkk(3,idws2).eq.m3.and.isblkk(4,idws2).eq.m4)then
           if(idwsdeb.gt.10)then
            write(6,*)('got a match ')
           end if
           call ilimts(n3,n4,numpro,nowpro,
     $          ilj,ihj,i1s,i1e,i2s,i2e)
           ii=0
           do i2=i2s,i2e
            i2p=i2+noc(m4)                                              3d19s12
            if(i2.eq.i2s)then
             i10=i1s
            else
             i10=1
            end if
            if(i2.eq.i2e)then
             i1n=i1e
            else
             i1n=n3
            end if
            do i1=i10,i1n
             i1p=i1+noc(m3)                                             3d19s12
             do i=icore(m2)+1,noc(m2)
              do j=icore(m1)+1,noc(m1)
               denom=2d0/(bc(ifeig(m2)+i)+bc(ifeig(m1)+j)
     $               -bc(ifeig(m3)+i1p)-bc(ifeig(m4)+i2p))              3d19s12
               k1=kmats(idws)+j-1+noc(m1)*(i-1+noc(m2)*ii)              3d19s12
               k2=kmats(idws2)+i-1+noc(m2)*(j-1+noc(m1)*ii)             3d19s12
               emp2v=emp2v+(bc(k1)*(2d0*bc(k1)-bc(k2))
     $                      +bc(k2)*(2d0*bc(k2)-bc(k1)))*denom
               tmp1=bc(k1)*denom
               tmp2=bc(k2)*denom
               sp=sp+tmp1*bc(k1)+tmp2*bc(k2)                            10d12s04
               tp=tp+tmp1*bc(k2)+tmp2*bc(k1)                            10d12s04
              end do
             end do
             ii=ii+1
            end do
           end do
          end if
         end do
        end if
       end if
      end do
      call dws_sync
      iarg1=1
      call dws_gsumf(emp2v,iarg1)                                       1d25s05
      sums(1)=sp                                                        10d12s04
      sums(2)=tp                                                        10d12s04
      iarg1=2
      call dws_gsumf(sums,iarg1)                                        1d25s05
      emp2v=sums(1)*2d0-sums(2)                                         10d12s04
      sp=0.5d0*(sums(1)+sums(2))                                        10d12s04
      tp=1.5d0*(sums(1)-sums(2))                                        10d12s04
      call second(time2)
      telap(1)=time2-time1
      telap(2)=telap(1)**2
      call dws_gsumf(telap,iarg1)                                       1d25s05
      if(lprint)write(6,*)('total elapsed time '),telap(1)
      tavg=telap(1)/dfloat(numpro)
      tdev=telap(2)/dfloat(numpro)
      tdev=sqrt(abs(tdev-tavg*tavg))
      if(lprint)then
      writE(6,*)('ave time '),tavg,(' rms dev '),tdev
      write(6,*)('valence correlation energy '),emp2v
      write(6,*)('singlet pair contribution '),sp
      write(6,*)('triplet pair contribution '),tp
      writE(6,*)('total energy '),emp2v+ehf
      end if
      return
      end
