c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine orbder(idamat,igmat,noc,nvirtc,ipuse,multh,morb,       3d29s16
     $     propmat,nbasdwsc,dere1,h0mo,ihessd,iflag,bodc,idwsdeb,ider3, 1d24s17
     $     ionex,i3x,bodcfact,bc,ibc,npass)                             3d2s23
      implicit real*8 (a-h,o-z)
c
c     solve for orbital ders given damat and gmat.
c
      dimension idamat(8),igmat(8,8),noc(8),nvirtc(8),multh(8,8),
     $     ipt(8),morb(1),propmat(1),nbasdwsc(1),h0mo(1),ipt2(8),       12d12s16
     $     ihessd(8,8),ipt3(8)                                          12d12s16
      include "common.store"
      include "common.hf"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      ibcoffo=ibcoff
      igoal=2
      kgoal=2
      n2ndoff=0                                                         12d12s16
      n3ndoff=0                                                         12d12s16
      do isb=1,nsymb                                                    12d12s16
       isk=multh(isb,ipuse)                                             12d12s16
       n2ndoff=n2ndoff+nbasdwsc(isb)*nbasdwsc(isk)                      12d12s16
       n3ndoff=n3ndoff+nbasdwsc(isb)*nbasdwsc(isb)                      12d12s16
      end do                                                            12d12s16
      n3ndoff=n3ndoff+n2ndoff                                           12d12s16
      ioff=1                                                            3d29s16
      do isb=1,nsymb                                                    3d29s16
       isk=multh(isb,ipuse)                                             3d29s16
       if(isk.ge.isb)then                                               4d28s22
        ipt2(isb)=ioff                                                   4d7s16
        ioff=ioff+nbasdwsc(isb)*nbasdwsc(isk)                            3d29s16
        if(isk.ne.isb)then                                              4d28s22
         ipt2(isk)=ioff                                                 4d28s22
         ioff=ioff+nbasdwsc(isb)*nbasdwsc(isk)                          4d28s22
        end if                                                          4d28s22
       end if                                                           4d28s22
      end do                                                            3d29s16
      do isb=1,nsymb                                                    12d12s16
       ipt3(isb)=ioff                                                   12d12s16
       ioff=ioff+nbasdwsc(isb)*nbasdwsc(isb)                            12d12s16
      end do                                                            12d12s16
c
c     solve linear equations via conjugate gradient
c
      nun=0                                                             4d4s23
      do isb=1,nsymb
       isk=multh(isb,ipuse)
       ipt(isb)=nun                                                     4d12s23
       nun=nun+nvirtc(isb)*noc(isk)                                     4d4s23
       if(idwsdeb.gt.10)then
        write(6,*)('damat for symmetry block '),isb,isk,idamat(isb),
     $      loc(bc(idamat(isb)))
        write(6,*)('ipt = '),ipt(isb)
        call prntm2(bc(idamat(isb)),nvirtc(isb),noc(isk),nvirtc(isb))
       end if
      end do
      nunm=nun-1
      isoln=ibcoff                                                      4d4s23
      ix=isoln+nun                                                      4d4s23
      ir=ix+nun
      iz=ir+nun
      ip=iz+nun
      igdiag=ip+nun
      iq=igdiag+nun
      ibcoff=iq+nun                                                     4d4s23
      call enough('orbder.  2',bc,ibc)
      do i=0,nunm
       bc(igdiag+i)=0d0
       bc(ir+i)=bc(idamat(1)+i)
      end do
      if(idwsdeb.gt.10)then                                             12d17s16
       write(6,*)('ibcoff = '),ibcoff
       write(6,*)('we we have for rhs ')
       do isb=1,nsymb
        isk=multh(isb,ipuse)
        write(6,*)('for symmetry block  '),isb,isk,idamat(isb)
        call prntm2(bc(idamat(isb)),nvirtc(isb),noc(isk),
     $       nvirtc(isb))
        ioff=idamat(isb)-1
        do i=1,noc(isk)
         write(6,*)(bc(ioff+j),j=1,nvirtc(isb))
         ioff=ioff+nvirtc(isb)
        end do
       end do
      end if                                                            12d17s16
      if(idwsdeb.gt.10)then                                             12d17s16
       write(6,*)('what we have for gmat: ')
       do isb=1,nsymb
        do isk=1,nsymb
         if(igmat(isb,isk).gt.0)then
          isbo=multh(isb,ipuse)
          isko=multh(isk,ipuse)
          write(6,*)('symmetry blocks '),isbo,isko,isk,isb,
     $         igmat(isb,isk)
          call ilimts(nvirtc(isk),nvirtc(isb),mynprocg,mynowprog,       5d2s22
     $        il,ih,i1s,i1e,i2s,i2e)
          nrow=noc(isbo)*noc(isko)
          nhere=ih+1-il
          call prntm2(bc(igmat(isb,isk)),nrow,nhere,nrow)
         end if
        end do
       end do
      end if                                                            12d17s16
      do isb=1,nsymb
       isbv=isb
       if(igmat(isb,isb).gt.0)then                                      3d29s16
        isbo=multh(isb,ipuse)                                           3d29s16
        call ilimts(nvirtc(isbv),nvirtc(isbv),mynprocg,mynowprog,       3d29s16
     $       il,ih,i1s,i1e,i2s,i2e)
        i10=i1s
        i1n=nvirtc(isbv)
        ii=igmat(isb,isb)
        nrow=noc(isbo)*noc(isbo)
        nhere=ih+1-il
        jgdiag=igdiag+ipt(isbv)
        do i2=i2s,i2e
         if(i2.eq.i2e)i1n=i1e
         do i1=i10,i1n
          if(i1.eq.i2)then
           do i34=0,noc(isbo)-1
            iad1=jgdiag+i2-1+nvirtc(isbv)*i34
            iad2=ii+i34*(noc(isbo)+1)
            bc(iad1)=bc(iad2)
           end do
          end if
          ii=ii+nrow
         end do
         i10=1
        end do
       end if
      end do
      call dws_gsumf(bc(igdiag),nun)
      if(idwsdeb.gt.10)then                                             12d17s16
       write(6,*)('diagonal elements of hessian: ')
       call prntm2(bc(igdiag),1,nun,1)
      end if                                                            12d17s16
      do i=0,nunm                                                       4d4s23
       if(abs(bc(igdiag+i)).gt.1d-10)then                               4d4s23
        bc(igdiag+i)=1d0/bc(igdiag+i)                                   4d4s23
       else                                                             4d4s23
        bc(igdiag+i)=0d0                                                4d4s23
       end if                                                           4d4s23
      end do                                                            4d4s23
      xdot=0d0                                                          4d4s23
      dotrz=0d0                                                         4d4s23
      do i=0,nunm                                                       4d4s23
       bc(iz+i)=bc(ir+i)*bc(igdiag+i)                                   4d4s23
       xdot=xdot+bc(ir+i)**2                                            4d4s23
       dotrz=dotrz+bc(iz+i)*bc(ir+i)                                    4d4s23
      end do                                                            4d4s23
      do i=0,nunm                                                       4d4s23
       bc(isoln+i)=0d0                                                  4d4s23
      end do                                                            4d4s23
      scale=1d0/xdot                                                    4d4s23
      macit=0                                                           5d9s22
      iter=0
      itmax=10
      ierror=ibcoff                                                     4d4s23
      ibcoff=ierror+itmax                                               4d4s23
      tol=1d-14
      nrepreive=0                                                       4d19s16
  100 continue                                                          4d4s23
       macit=macit+1                                                     5d9s22
       if(macit.gt.10)then                                               5d9s22
        write(6,*)('too many restarts!!!')
        call dws_synca
        call dws_finalize
        stop 'orbder'
       end if
       do i=0,nunm
        bc(ix+i)=0d0
       end do
       do i=0,nunm
        bc(ip+i)=bc(iz+i)                                                4d3s23
       end do
       iter=0
    1  continue
        iter=iter+1
        if(iter.gt.itmax)then
         if(mynowprog.eq.0.and.idwsdeb.ne.10)then                         7d14s22
          write(6,52)(bc(ierror+izz),izz=0,iter-2)                          5d19s22
         end if
         do i=0,nunm
          bc(isoln+i)=bc(isoln+i)+bc(ix+i)                                4d4s23
         end do                                                           4d4s23
         do i=0,nunm
          bc(iq+i)=0d0
         end do
         do isb=1,nsymb                                                   3d29s16
          isbv=isb                                                        3d29s16
          do isa=1,nsymb                                                  3d29s16
           isav=isa                                                       3d29s16
           if(igmat(isa,isb).gt.0)then                                    3d29s16
            isao=multh(isa,ipuse)                                         3d29s16
            isbo=multh(isb,ipuse)                                         3d29s16
            call ilimts(nvirtc(isbv),nvirtc(isav),mynprocg,mynowprog,      3d29s16
     $         il,ih,i1s,i1e,i2s,i2e)
            i10=i1s
            i1n=nvirtc(isbv)
            ii=igmat(isa,isb)
            jp=isoln+ipt(isav)                                          4d4s23
            jq=iq+ipt(isbv)                                               5d2s22
            do i2=i2s,i2e
             i2m=i2-1
             if(i2.eq.i2e)i1n=i1e
             do i1=i10,i1n
              i1m=i1-1
              do i3=0,noc(isbo)-1
               do i4=0,noc(isao)-1
                iadp=jp+i2m+nvirtc(isav)*i4                               5d2s22
                iadq=jq+i1m+nvirtc(isbv)*i3                               5d2s22
                bc(iadq)=bc(iadq)+bc(ii)*bc(iadp)                       4d4s23
                ii=ii+1
               end do
              end do
             end do
             i10=1
            end do
            if(isa.ne.isb)then                                            3d29s16
             i10=i1s
             i1n=nvirtc(isbv)
             ii=igmat(isa,isb)
             jp=isoln+ipt(isbv)                                         4d4s23
             jq=iq+ipt(isav)                                              5d2s22
             do i2=i2s,i2e
              i2m=i2-1
              if(i2.eq.i2e)i1n=i1e
              do i1=i10,i1n
               i1m=i1-1
               do i3=0,noc(isbo)-1
                do i4=0,noc(isao)-1
                 iadp=jp+i1m+nvirtc(isbv)*i3                              5d2s22
                 iadq=jq+i2m+nvirtc(isav)*i4                              5d2s22
                 bc(iadq)=bc(iadq)+bc(ii)*bc(iadp)
                 ii=ii+1
                end do
               end do
              end do
              i10=1
             end do
            end if                                                        3d29s16
           end if
          end do
         end do
         call dws_gsumf(bc(iq),nun)
         xdot=0d0                                                         4d4s23
         dotrz=0d0                                                        4d4s23
         do i=0,nunm                                                     4d3s23
          bc(ir+i)=bc(idamat(1)+i)-bc(iq+i)                              4d3s23
          bc(iz+i)=bc(ir+i)*bc(igdiag+i)                                 4d3s23
          xdot=xdot+bc(ir+i)**2
          dotrz=dotrz+bc(ir+i)*bc(iz+i)                                  4d3s23
         end do
         xtest=sqrt(xdot*scale)                                          4d4s23
         go to 100                                                        4d4s23
        end if
         if(iter.ne.1)then                                              4d4s23
          beta=dotrznew/dotrz                                           4d4s23
          dotrz=dotrznew                                                4d4s23
          do i=0,nunm                                                     4d3s23
           bc(ip+i)=bc(iz+i)+beta*bc(ip+i)                                4d3s23
          end do                                                          4d3s23
         end if                                                           4d3s23
         do i=0,nunm
          bc(iq+i)=0d0
         end do
         do isb=1,nsymb                                                   3d29s16
          isbv=isb                                                        3d29s16
          do isa=1,nsymb                                                  3d29s16
           isav=isa                                                       3d29s16
           if(igmat(isa,isb).gt.0)then                                    3d29s16
            isao=multh(isa,ipuse)                                         3d29s16
            isbo=multh(isb,ipuse)                                         3d29s16
            call ilimts(nvirtc(isbv),nvirtc(isav),mynprocg,mynowprog,      3d29s16
     $         il,ih,i1s,i1e,i2s,i2e)
            i10=i1s
            i1n=nvirtc(isbv)
            ii=igmat(isa,isb)
            jp=ip+ipt(isav)                                               5d2s22
            jq=iq+ipt(isbv)                                               5d2s22
            do i2=i2s,i2e
             i2m=i2-1
             if(i2.eq.i2e)i1n=i1e
             do i1=i10,i1n
              i1m=i1-1
              do i3=0,noc(isbo)-1
               do i4=0,noc(isao)-1
                iadp=jp+i2m+nvirtc(isav)*i4                               5d2s22
                iadq=jq+i1m+nvirtc(isbv)*i3                               5d2s22
                bc(iadq)=bc(iadq)+bc(ii)*bc(iadp)
                ii=ii+1
               end do
              end do
             end do
             i10=1
            end do
            if(isa.ne.isb)then                                            3d29s16
             i10=i1s
             i1n=nvirtc(isbv)
             ii=igmat(isa,isb)
             jp=ip+ipt(isbv)                                              5d2s22
             jq=iq+ipt(isav)                                              5d2s22
             do i2=i2s,i2e
              i2m=i2-1
              if(i2.eq.i2e)i1n=i1e
              do i1=i10,i1n
               i1m=i1-1
               do i3=0,noc(isbo)-1
                do i4=0,noc(isao)-1
                 iadp=jp+i1m+nvirtc(isbv)*i3                              5d2s22
                 iadq=jq+i2m+nvirtc(isav)*i4                            4d12s23
                 bc(iadq)=bc(iadq)+bc(ii)*bc(iadp)                      4d12s23
                 ii=ii+1
                end do
               end do
              end do
              i10=1
             end do
            end if                                                        3d29s16
           end if
          end do
         end do
         call dws_gsumf(bc(iq),nun)
         alpha=0d0                                                        4d3s23
         do i=0,nunm
          alpha=alpha+bc(ip+i)*bc(iq+i)                                   4d3s23
         end do                                                           4d3s23
         alpha=dotrz/alpha                                                4d3s23
         xdotnew=0d0                                                      4d3s23
         dotrznew=0d0                                                     4d3s23
         do i=0,nunm                                                      4d3s23
          bc(ix+i)=bc(ix+i)+alpha*bc(ip+i)                                4d3s23
          bc(ir+i)=bc(ir+i)-alpha*bc(iq+i)                                4d3s23
          bc(iz+i)=bc(ir+i)*bc(igdiag+i)                                  4d3s23
          xdotnew=xdotnew+bc(ir+i)**2                                     4d3s23
          dotrznew=dotrznew+bc(ir+i)*bc(iz+i)                             4d3s23
         end do                                                           4d3s23
         xtest=sqrt(xdotnew*scale)                                        4d4s23
         bc(ierror+iter-1)=xtest                                          4d3s23
    2    format('after iteration ',i2,' error is ',es10.3)
         if(xtest.gt.tol)go to 1
         if(mynowprog.eq.0.and.idwsdeb.ne.10)then                                            3d3s17
          write(6,*)('cg iterations converged after '),iter,
     $      (' iterations and '),macit-1,(' restarts')
   52   format(10es8.1)                                                 5d19s22
        end if                                                             3d3s17
        do i=0,nunm                                                       4d3s23
         bc(isoln+i)=bc(isoln+i)+bc(ix+i)                               4d4s23
        end do                                                            5d16s22
        ioff=0
        do isb=1,nsymb
         isk=multh(isb,ipuse)
         if(idwsdeb.gt.10)then
          write(6,*)('solution for symmetry '),isb,isk,ioff,
     $      isoln+ioff
          call prntm2(bc(isoln+ioff),nvirtc(isb),noc(isk),
     $      nvirtc(isb))
         end if
       koff=isoln+ioff-1
       do i=1,noc(isk)
        koff=koff+nvirtc(isb)
       end do
       ioff=ioff+nvirtc(isb)*noc(isk)
      end do
      if(iflag.eq.1)then
       jsoln=isoln                                                      4d4s23
       itmp=ibcoff                                                      12d12s16
       if(mynprocg.eq.1)then
        igoal=1893746
       else if(mynowprog.eq.0)then
        igoal=1786756
       else if(mynowprog.eq.1)then
        igoal=1783464
       else if(mynowprog.eq.2)then
        igoal=1779963
       else if(mynowprog.eq.3)then
        igoal=1777835
       else if(mynowprog.eq.4)then
        igoal=1776670
       else
        igoal=1777210
       end if
       igoal=igoal+2
       nwds=0                                                           12d12s16
       npoff=1                                                          12d12s16
       do isa=1,nsymb                                                   12d12s16
        nwds=nwds+noc(isa)*nvirtc(isa)                                  12d12s16
        isk=multh(isa,ipuse)                                            12d12s16
        npoff=npoff+nbasdwsc(isa)*(nbasdwsc(isa)+nbasdwsc(isk))         12d12s16
       end do                                                           12d12s16
       ibcoff=itmp+nwds                                                 12d12s16
       call enough('orbder.  4',bc,ibc)
       do i=0,nwds-1                                                    12d12s16
        bc(itmp+i)=0                                                    12d12s16
       end do                                                           12d12s16
       do iskv=1,nsymb                                                   12d12s16
        isko=multh(iskv,ipuse)                                            12d12s16
        if(noc(isko).gt.0)then                                           12d12s16
c                                                                       1d6s17
c     da^2 contribution to (n|d2|m)                                     1d6s17
c                                                                       1d6s17
         itrans=ibcoff                                                  1d6s17
         ibcoff=itrans+nvirtc(iskv)*noc(isko)                            1d6s17
         call enough('orbder.  5',bc,ibc)
         do i=0,noc(isko)-1                                             1d6s17
          do j=0,nvirtc(iskv)-1                                         1d6s17
           ji=jsoln+j+nvirtc(iskv)*i                                    14d4s23
           ij=itrans+i+noc(isko)*j                                      1d6s17
           bc(ij)=bc(ji)                                                1d6s17
          end do                                                        1d6s17
         end do                                                         1d6s17
         if(idwsdeb.gt.10)then
         write(6,*)('starting (n|d2/dq2|m) matrix '),isko,ipt3(isko)
         call prntm2(propmat(ipt3(isko)),nbasdwsc(isko),nbasdwsc(isko), 1d6s17
     $        nbasdwsc(isko))                                           1d6s17
         if(nsymb.eq.1)call printa(propmat(ipt3(isko)),nbasdwsc,0,1,0,
     $        nbasdwsc,0,1,0,bc(ibcoff))
         end if
         if(idwsdeb.gt.10)then
         write(6,*)('transpose of da ')
         call prntm2(bc(itrans),noc(isko),nvirtc(iskv),noc(isko))
         if(nsymb.eq.1)call printa(bc(itrans),noc,0,1,0,nvirtc,noc,1,0
     $        ,bc(ibcoff))
         write(6,*)('4dt T dt')
         end if
         if(idwsdeb.gt.10)then
         call dgemm('n','n',noc(isko),noc(isko),nvirtc(iskv),-4d0,
     $        bc(itrans),noc(isko),bc(jsoln),nvirtc(iskv),0d0,          4d4s23
     $        bc(ibcoff),noc(isko),
     d' orbder.  1')
         call prntm2(bc(ibcoff),noc(isko),noc(isko),noc(isko))
         end if
c
c     the factor of 4 comes from using half the gradient to generate
c     the first derivative.
c
         if(idwsdeb.gt.10)then
         write(6,*)('starting propmata '),ipt3(isko),isko
         call prntm2(propmat(ipt3(isko)),nbasdwsc(isko),noc(isko),
     $        nbasdwsc(isko))
         if(nsymb.eq.1)call printa(propmat(ipt3(isko)),nbasdwsc,0,
     $        1,0,noc,0,1,0,bc(ibcoff))
         write(6,*)('jsoln: ')
         call prntm2(bc(jsoln),nvirtc(iskv),noc(isko),nvirt(iskv))      4d4s23
         if(nsymb.eq.1)call printa(bc(jsoln),nvirtc,noc,1,0,            4d4s23
     $        noc,0,1,0,bc(ibcoff))
         end if
         call dgemm('n','n',noc(isko),noc(isko),nvirtc(iskv),-4d0,      1d6s17
     $        bc(itrans),noc(isko),bc(jsoln),nvirtc(iskv),1d0,          4d4s23
     $        propmat(ipt3(isko)),nbasdwsc(isko),                       1d6s17
     d' orbder.  2')
         if(idwsdeb.gt.10)then
         write(6,*)('(n|d2/dq2|m) matrix after  daT da')                1d6s17
         call prntm2(propmat(ipt3(isko)),nbasdwsc(isko),nbasdwsc(isko), 1d6s17
     $        nbasdwsc(isko))                                           1d6s17
         if(nsymb.eq.1)call printa(propmat(ipt3(isko)),nbasdwsc,0,1,0,
     $        nbasdwsc,0,1,0,bc(ibcoff))
         end if
         if(idwsdeb.gt.10)then
         write(6,*)('starting (n|d2/dq2|m) matrix '),ipt3(iskv)
         call prntm2(propmat(ipt3(iskv)),nbasdwsc(iskv),nbasdwsc(iskv), 1d6s17
     $        nbasdwsc(iskv))                                           1d6s17
         if(nsymb.eq.1)call printa(propmat(ipt3(iskv)),nbasdwsc,0,1,0,
     $        nbasdwsc,0,1,0,bc(ibcoff))
         end if
         jpt3=ipt3(iskv)+noc(iskv)*(1+nbasdwsc(iskv))                   1d6s17
         if(idwsdeb.gt.10)then
         write(6,*)('jsoln ')
         call prntm2(bc(jsoln),nvirtc(iskv),noc(isko),nvirtc(iskv))
         if(nsymb.eq.1)call printa(bc(jsoln),nvirtc,noc,1,0,noc,0,
     $        1,0,bc(ibcoff))
         write(6,*)('trans')
         call prntm2(bc(itrans),noc(isko),nvirtc(iskv),noc(isko))
         if(nsymb.eq.1)call printa(bc(itrans),noc,0,1,0,nvirtc,noc,1,0,
     $        bc(ibcoff))
         write(6,*)('starting propmatb '),jpt3,iskv,isko
         call prntm2(propmat(jpt3),nvirtc(iskv),noc(iskv),              2d24s23
     $        nbasdwsc(iskv))
         if(nsymb.eq.1)then
          icpyx=ibcoff
          ibcoff=icpyx+nvirtc(1)*noc(1)
          call enough('orbder.cpyx',bc,ibc)
          jcpyx=icpyx
          do i1=0,noc(1)-1
           iad=jpt3+nbasdwsc(1)*i1
           do i2=0,nvirtc(1)-1
            bc(jcpyx+i2)=propmat(iad+i2)
           end do
           jcpyx=jcpyx+nvirtc(1)
          end do
          call printa(bc(icpyx),nvirtc,noc,1,0,noc,0,
     $        1,0,bc(ibcoff))
          ibcoff=icpyx
         end if
         end if
         call dgemm('n','n',nvirtc(iskv),nvirtc(iskv),noc(isko),-4d0,   1d6s17
     $        bc(jsoln),nvirtc(iskv),bc(itrans),noc(isko),1d0,          1d6s17
     $        propmat(jpt3),nbasdwsc(iskv),                             1d6s17
     d' orbder.  3')
         if(idwsdeb.gt.10)then
         write(6,*)('(n|d2/dq2|m) matrix after  da daT')                1d6s17
         call prntm2(propmat(ipt3(iskv)),nbasdwsc(iskv),nbasdwsc(iskv), 1d6s17
     $        nbasdwsc(iskv))                                           1d6s17
         if(nsymb.eq.1)call printa(propmat(ipt3(iskv)),nbasdwsc,0,1,0,
     $        nbasdwsc,0,1,0,bc(ibcoff))
         end if
         ibcoff=itrans
         jtmp=itmp                                                      12d12s16
         igoal=itmp
         do isb=1,nsymb                                                  12d12s16
          if(noc(isb).gt.0)then                                         12d12s16
           nrow=noc(isb)*noc(isko)                                         12d12s16
           call ilimts(nvirtc(iskv),nvirtc(isb),mynprocg,mynowprog,ilh,  12d12s16
     $          ihh,i1sh,i1eh,i2sh,i2eh)                                12d12s16
           nhere=ihh+1-ilh                                              12d12s16
           ii=ihessd(isb,iskv)                                           12d12s16
           i1n=nvirtc(iskv)                                              12d12s16
           i10=i1sh                                                     12d12s16
           do i2=i2sh,i2eh                                              12d12s16
            if(i2.eq.i2eh)i1n=i1eh                                      12d12s16
            i2m=i2-1                                                    12d12s16
            do i1=i10,i1n                                               12d12s16
             i1m=i1-1                                                   12d12s16
             do i4=0,noc(isko)-1                                         12d12s16
              do i3=0,noc(isb)-1                                        12d12s16
               iaddg=jsoln+i1m+nvirtc(iskv)*i4                          12d12s16
               iaddt=jtmp+i2m+nvirtc(isb)*i3                            12d12s16
               bc(iaddt)=bc(iaddt)+bc(ii)*bc(iaddg)                     12d12s16
               ii=ii+1                                                  12d12s16
              end do                                                    12d12s16
             end do                                                     12d12s16
            end do                                                      12d12s16
            i10=1                                                       12d12s16
           end do                                                       12d12s16
           jtmp=jtmp+noc(isb)*nvirtc(isb)                               12d12s16
          end if                                                        12d12s16
         end do                                                          12d12s16
         jsoln=jsoln+nvirtc(iskv)*noc(isko)                              12d12s16
        end if                                                          12d12s16
       end do                                                           12d12s16
       call dws_gsumf(bc(itmp),nwds)                                    12d12s16
         if(idwsdeb.gt.10)then
          jtmp=itmp
          do isa=1,nsymb
           write(6,*)('result of dH*dorb '),jtmp,(' sym '),isa                                 12d12
           call prntm2(bc(jtmp),nvirtc(isa),noc(isa),nvirtc(isa))
           jtmp=jtmp+nvirtc(isa)*noc(isa)                               2d24s23
          end do                                                        2d24s23
       if(nsymb.eq.1)call printa(bc(itmp),nvirtc,noc,1,0,noc,0,1,0,
     $      bc(ibcoff))
       end if
       jtmp=itmp                                                        12d12s16
       npoff=1+n3ndoff                                                  12d12s16
       do isa=1,nsymb                                                   12d12s16
        if(noc(isa).gt.0)then                                           12d12s16
         call der3part(propmat,ider3,ionex,i3x,bc(isoln),noc,nvirtc,
     $      ipuse,multh,isa,ider3p,bc,ibc)                              11d14s22
         if(idwsdeb.gt.10)then
       write(6,*)('now for Cijk da*da contributions ')                  1d23s17
         call prntm2(bc(ider3p),nvirtc(isa),noc(isa),nvirtc(isa))
         if(nsymb.eq.1)call printa(bc(ider3p),nvirtc,noc,1,0,noc,0,1,0,
     $        bc(ibcoff))
         write(6,*)('starting '),npoff
         call prntm2(propmat(npoff),nvirtc(isa),noc(isa),nvirtc(isa))   12d12s16
         if(nsymb.eq.1)call printa(propmat(npoff),nvirtc,noc,1,0,noc,0,
     $        1,0,bc(ibcoff))
         write(6,*)('jtmp ')
         call prntm2(bc(jtmp),nvirtc(isa),noc(isa),nvirtc(isa))
         if(nsymb.eq.1)call printa(bc(jtmp),nvirtc,noc,1,0,noc,0,1,0,
     $        bc(ibcoff))
         end if
         do i=0,nvirtc(isa)*noc(isa)-1                                  12d12s16
          propmat(npoff+i)=propmat(npoff+i)
     $         -2d0*bc(jtmp+i)+2d0*bc(ider3p+i)                         1d26s17
         end do                                                         12d12s16
         ibcoff=ider3p                                                  1d26s17
         if(idwsdeb.gt.10)then
         write(6,*)('result '),npoff,loc(propmat(npoff))
         call prntm2(propmat(npoff),nvirtc(isa),noc(isa),nvirtc(isa))   12d12s16
         if(nsymb.eq.1)call printa(propmat(npoff),nvirtc,noc,1,0,
     $        noc,0,1,0,bc(ibcoff))
         end if
         npoff=npoff+nvirtc(isa)*noc(isa)                               12d12s16
        end if                                                          12d12s16
        jtmp=jtmp+nvirtc(isa)*noc(isa)                                  12d17s16
       end do                                                           12d12s16
       ibcoff=itmp                                                      12d12s16
      end if
      ibcoff=ix                                                         4d4s23
c
c     final contribution of (n|d/dq|m) is the mo matrix times the       3d29s16
c     orbital der matrix                                                3d29s16
c     the other parts of the matrix are skew symmetric, so for the isk lt    3d29s16
c     isb block, we have no contribution yet                            3d29s16
c
      if(ipuse.ne.1)then                                                3d29s16
       do isb=1,nsymb                                                   3d29s16
        isk=multh(isb,ipuse)                                            3d29s16
        if(isk.gt.isb)then                                              3d29s16
         do i=0,nbasdwsc(isk)-1                                         3d29s16
          do j=0,nbasdwsc(isb)-1                                        3d29s16
           ji=j+nbasdwsc(isb)*i                                         3d29s16
           ij=i+nbasdwsc(isk)*j                                         3d29s16
           propmat(ipt2(isk)+ij)=-propmat(ipt2(isb)+ji)                 4d7s16
          end do                                                        3d29s16
         end do                                                         3d29s16
        end if                                                          3d29s16
       end do                                                           3d29s16
      end if                                                            3d29s16
      if(iflag.eq.1)then                                                12d12s16
       ipcopy=ibcoff                                                     12d12s16
       ibcoff=ipcopy+n2ndoff                                             12d12s16
       call enough('orbder.  6',bc,ibc)
       jpcopy=ipcopy-1                                                   12d12s16
       do i=1,n2ndoff                                                    12d12s16
        bc(jpcopy+i)=propmat(i)                                          12d12s16
       end do                                                            12d12s16
      end if                                                            12d12s16
      jsoln=isoln                                                       4d4s23
      derbodc=0d0                                                       12d28s16
      do isb=1,nsymb                                                    3d29s16
       isk=multh(isb,ipuse)                                             3d29s16
       itmp1=ibcoff                                                       3d29s16
       ibcoff=itmp1+nbasdwsc(isb)*nbasdwsc(isk)                         12d12s16
       call enough('orbder.  7',bc,ibc)
       do i=0,nbasdwsc(isb)*nbasdwsc(isk)-1                             3d29s16
        bc(itmp1+i)=0d0                                                 3d29s16
       end do                                                           3d29s16
       if(isb.eq.isk)then
        do icol=0,noc(isb)-1                                             3d29s16
         do irow=0,nvirtc(isb)-1                                          3d29s16
          irowp=irow+noc(isb)                                            4d7s16
          iad1=itmp1+irowp+nbasdwsc(isb)*icol                           4d7s16
          iad2=itmp1+icol+nbasdwsc(isb)*irowp                           4d7s16
          val=bc(jsoln)*2d0                                             4d4s23
          bc(iad1)=val                                                   3d29s16
          bc(iad2)=-val                                                  3d29s16
          jsoln=jsoln+1                                                 4d4s23
         end do                                                          3d29s16
        end do                                                           3d29s16
       else                                                             4d7s16
        igb=isoln+ipt(isb)                                              4d4s23
        igk=isoln+ipt(isk)                                              4d4s23
        do icol=0,noc(isk)-1
         do irow=0,nvirtc(isb)-1
          irowp=irow+noc(isb)
          iad1=itmp1+irowp+nbasdwsc(isb)*icol
          val=bc(igb)*2d0
          bc(iad1)=val
          igb=igb+1
         end do
        end do
        do icol=0,noc(isb)-1
         do irow=0,nvirtc(isk)-1
          irowp=irow+noc(isk)
          iad2=itmp1+icol+nbasdwsc(isb)*irowp
          val=bc(igk)*2d0
          bc(iad2)=-val
          igk=igk+1
         end do
        end do
       end if                                                           4d7s16
       if(iflag.eq.1)then                                               12d12s16
        call dgemm('n','n',nbasdwsc(isk),nbasdwsc(isk),nbasdwsc(isb),    12d12s16
     $     -2d0,bc(jpcopy+ipt2(isk)),nbasdwsc(isk),bc(itmp1),           12d28s16
     $      nbasdwsc(isb),1d0,propmat(ipt3(isk)),nbasdwsc(isk),         12d12s16
     d' orbder.  4')
       end if                                                           12d12s16
       if(iflag.eq.1)then                                               12d12s16
       else                                                             12d12s16
       end if                                                           12d12s16
       do iz=0,nbasdwsc(isb)*nbasdwsc(isk)-1
        propmat(ipt2(isb)+iz)=propmat(ipt2(isb)+iz)-bc(itmp1+iz)        5d9s16
       end do
       if(iflag.ne.1)then                                               12d28s16
c
c     2 is from 2 electrons/orbital
c
        do i=0,noc(isb)-1                                                12d28s16
         iad=ipt2(isb)+i*(nbasdwsc(isb)+1)                               12d28s16
         derbodc=derbodc+2d0*propmat(iad)                               3d1s17
        end do                                                           12d28s16
       end if                                                           12d28s16
       ibcoff=itmp1
      end do                                                            3d29s16
      if(iflag.eq.1)then                                                12d29s16
       do isb=1,nsymb                                                   12d28s16
        do ib=0,noc(isb)-1                                              12d28s16
c
c     j part - this doesn't happen because matrix element is skew
c     symmetric
c
c
c     k part
c
         isk=multh(isb,ipuse)                                           12d29s16
         do ik=0,noc(isk)-1                                             12d29s16
          iad=ipt2(isb)+ib+nbasdwsc(isb)*ik                             12d29s16
          derbodc=derbodc+(propmat(iad)**2)                             3d1s17
         end do                                                         12d29s16
        end do                                                          12d28s16
       end do                                                           12d28s16
c
c     slater's rules are for sum i=1,n sum j=i+1,n 1/rij, but           3d1s17
c     for here we get 2 times that. Slater's rules tell us that         3d1s17
c     the k part is -0.5 sum n,m (m|dn)(n|dm). since (n|dm)=-(m|dn),    3d1s17
c     we get +0.5 sum n,m (m|dn)^2. splitting sums to be just over      3d1s17
c     space orbs and separate sums over spins, we get                    3d1s17
c     sum n,m (m|dn)^2.                                                 3d1s17
c     this is then multiplied by 2                                      3d1s17
c
       derbodc=derbodc*2d0
      end if                                                            12d28s16
c
c     factor of -1/2 is from KE operator
c
      derbodc=-derbodc*0.5d0*bodcfact                                   5d4s22
      bodc=bodc+derbodc                                                 12d29s16
      ibcoff=ibcoffo                                                    3d29s16
      return
      end
