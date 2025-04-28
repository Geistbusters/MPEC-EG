c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine trans2e(jmats,kmats,ionex,i3x,i3xb,irefo,iuncd,ionexb, 4d24s20
     $     jmatt,kmatt,ionexc,kmatd,i3x3,ionexbt,bc,ibc)                11d10s22
      implicit real*8 (a-h,o-z)
c
c     for mrci, transpose order of indices so virts come first
c
      include "common.store"
      include "common.hf"
      dimension jmats(*),kmats(*),ionex(*),irefo(*),i3x(*),i3xb(*),     1d25s20
     $     ionexb(*),jmatt(*),kmatt(*),ionexc(*),kmatd(8,8),i3x3(8),    3d23s21
     $     ionexbt(*)                                                   3d23s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      common/kmfind/invk1(2,8,8,8,2)                                    6d30s18
      do is=1,nsdlk
       call ilimts(nvirt(isblk(3,is)),nvirt(isblk(4,is)),mynprocg,
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)
       nhere=ih+1-il
       if(isblk(1,is).eq.isblk(2,is))then
        nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2
       else
        nrow=irefo(isblk(1,is))*irefo(isblk(2,is))
       end if
       need=nrow*nhere
       if(need.gt.0)then
        jmatt(is)=ibcoff                                                4d24s20
        ibcoff=jmatt(is)+need                                           4d24s20
        call enough('trans2e.  1',bc,ibc)
        do i=0,need-1                                                   4d24s20
         bc(jmatt(is)+i)=bc(jmats(is)+i)                                4d24s20
        end do                                                          4d24s20
        itmp=ibcoff
        ibcoff=itmp+need
        call enough('trans2e.  2',bc,ibc)
        do irow=0,nrow-1
         do icol=0,nhere-1
          ji=jmats(is)+irow+nrow*icol
          ij=itmp+icol+nhere*irow
          bc(ij)=bc(ji)
         end do
        end do
        do i=0,need-1
         bc(jmats(is)+i)=bc(itmp+i)
        end do
        ibcoff=itmp
       end if
      end do
      do is=1,nsdlkk
       call ilimts(nvirt(isblkk(3,is)),nvirt(isblkk(4,is)),mynprocg,
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)
       nhere=ih+1-il
       nrow=irefo(isblkk(1,is))*irefo(isblkk(2,is))
       need=nrow*nhere
       if(need.gt.0)then
        kmatt(is)=ibcoff                                                4d27s20
        ibcoff=kmatt(is)+need                                           4d27s20
        call enough('trans2e.  3',bc,ibc)
        do i=0,need-1                                                   4d27s20
         bc(kmatt(is)+i)=bc(kmats(is)+i)                                4d27s20
        end do                                                          4d27s20
        itmp=ibcoff
        ibcoff=itmp+need
        call enough('trans2e.  4',bc,ibc)
        do irow=0,nrow-1
         do icol=0,nhere-1
          ji=kmats(is)+irow+nrow*icol
          ij=itmp+icol+nhere*irow
          bc(ij)=bc(ji)
         end do
        end do
        do i=0,need-1
         bc(kmats(is)+i)=bc(itmp+i)
        end do
        ibcoff=itmp
       end if
      end do
      ib0=ibcoff                                                        8d28s20
      do isb=1,nsymb                                                    8d19s20
       call ilimts(nvirt(isb),nvirt(isb),mynprocg,mynowprog,            8d19s20
     $              il,ih,i1s,i1e,i2s,i2e)                              5d27s20
       nhere=ih+1-il                                                    8d19s20
       do isa=1,nsymb                                                   8d19s20
        nv=(irefo(isa)*(irefo(isa)+1))/2                                8d19s20
        kmatd(isa,isb)=ibcoff                                           8d19s20
        ibcoff=kmatd(isa,isb)+nv*nvirt(isb)                             8d19s20
        call enough('trans2e.  5',bc,ibc)
        do i=kmatd(isa,isb),ibcoff-1                                    8d29s20
         bc(i)=0d0                                                      8d29s20
        end do                                                          8d29s20
        if(min(nhere,nv).gt.0)then                                      8d19s20
         i2eu=invk1(1,isa,isa,isb,1)                                    8d19s20
         do i2=i2s,i2e                                                   8d19s20
          i2m=i2-1                                                       8d19s20
          irow=i2+nvirt(isb)*i2m                                         8d19s20
          if(irow.ge.il.and.irow.le.ih)then                              8d19s20
           iint=kmats(i2eu)+irow-il                                     8d19s20
           do ia=0,irefo(isa)-1                                         8d19s20
            do ic=0,ia                                                  8d19s20
             irec=ic+irefo(isa)*ia                                      8d19s20
             itri=((ia*(ia+1))/2)+ic                                    8d19s20
             iad1=iint+nhere*irec                                       8d19s20
             iad2=kmatd(isa,isb)+i2m+nvirt(isb)*itri                    8d19s20
             bc(iad2)=bc(iad1)                                          8d19s20
            end do                                                      8d19s20
           end do                                                       8d19s20
          end if                                                        8d19s20
         end do                                                          8d19s20
        end if                                                          8d19s20
       end do                                                           8d19s20
      end do                                                            8d19s20
      do isb=1,nsymb                                                    12d21s20
       i3x3(isb)=ibcoff                                                 12d21s20
       ibcoff=i3x3(isb)+irefo(isb)*nvirt(isb)                           12d21s20
       call enough('trans2e.  6',bc,ibc)
       do i=i3x3(isb),ibcoff-1                                          12d21s20
        bc(i)=0d0                                                       12d21s20
       end do                                                           12d21s20
       if(min(irefo(isb),nvirt(isb)).gt.0)then                          12d21s20
        i2eu=invk1(1,isb,isb,isb,2)                                     12d21s20
        call ilimts(irefo(isb),nvirt(isb),mynprocg,mynowprog,il,ih,i1s, 12d21s20
     $       i1e,i2s,i2e)                                               12d21s20
        iint=i3x(i2eu)                                                  12d21s20
        i10=i1s                                                         12d21s20
        i1n=irefo(isb)                                                  12d21s20
        nrow=(nvirt(isb)*(nvirt(isb)+1))/2                              12d21s20
        do i2=i2s,i2e                                                   12d21s20
         i2m=i2-1                                                       12d21s20
         if(i2.eq.i2e)i1n=i1e                                           12d21s20
         do i1=i10,i1n                                                  12d21s20
          iad=i3x3(isb)+i1-1+irefo(isb)*i2m                             12d21s20
          jad=iint+((i2m*(i2m+1))/2)+i2m                                12d21s20
          bc(iad)=bc(jad)                                               12d21s20
          iint=iint+nrow                                                12d21s20
         end do                                                         12d21s20
         i10=1                                                          12d21s20
        end do                                                          12d21s20
       end if                                                           12d21s20
      end do                                                            12d21s20
      nwds=ibcoff-ib0                                                   8d28s20
      call dws_gsumf(bc(ib0),nwds)                                      8d28s20
      do is=1,nsdlk1                                                    7d13s20
       if(isblk1(1,is).eq.isblk1(2,is))then                             7d13s20
        nrow=(irefo(isblk1(1,is))*(irefo(isblk1(1,is))+1))/2            7d13s20
       else                                                             7d13s20
        nrow=irefo(isblk1(1,is))*irefo(isblk1(2,is))                    7d13s20
       end if                                                           7d13s20
       call ilimts(irefo(isblk1(3,is)),nvirt(isblk1(4,is)),mynprocg,    12d11s19
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            12d12s19
       nhere=ih+1-il                                                    12d11s19
       ionexc(is)=ibcoff                                                7d13s20
       ibcoff=ionexc(is)+nhere*nrow                                     7d13s20
       call enough('trans2e.  7',bc,ibc)
       do i=0,nrow-1                                                    7d13s20
        do j=0,nhere-1                                                  7d13s20
         ji=ionexc(is)+j+nhere*i                                        7d13s20
         ij=ionex(is)+i+nrow*j                                          7d13s20
         bc(ji)=bc(ij)                                                  7d13s20
        end do                                                          7d13s20
       end do                                                           7d13s20
      end do                                                            7d13s20
      if(iuncd.ne.-1492)then                                            6d7s21
       do is=1,nsdlk1                                                    12d11s19
        if(isblk1(1,is).eq.isblk1(2,is))then                            1d25s20
         nrow=(irefo(isblk1(1,is))*(irefo(isblk1(1,is))+1))/2           1d25s20
        else                                                            1d25s20
         nrow=irefo(isblk1(1,is))*irefo(isblk1(2,is))                   1d25s20
        end if                                                          1d25s20
        nwds=nrow*irefo(isblk1(3,is))*nvirt(isblk1(4,is))               1d25s20
        ionexb(is)=ibcoff                                               1d25s20
        ionexbt(is)=ionexb(is)+nwds                                     3d23s21
        ibcoff=ionexbt(is)+nwds                                         3d23s21
        call enough('trans2e.  8',bc,ibc)
        do i=ionexb(is),ibcoff-1                                        1d25s20
         bc(i)=0d0                                                      1d25s20
        end do                                                          1d25s20
        call ilimts(irefo(isblk1(3,is)),nvirt(isblk1(4,is)),mynprocg,    12d11s19
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            12d12s19
        nhere=ih+1-il                                                    12d11s19
        if(nhere.gt.0)then                                              1d25s20
         i10=i1s                                                        1d25s20
         i1n=irefo(isblk1(3,is))                                        1d25s20
         ii=ionex(is)                                                   1d25s20
         do i2=i2s,i2e                                                  1d25s20
          if(i2.eq.i2e)i1n=i1e                                          1d25s20
          do i1=i10,i1n                                                 1d25s20
           iad1=ionexb(is)+nrow*(i1-1+irefo(isblk1(3,is))*(i2-1))       1d25s20
           do i34=0,nrow-1                                              1d25s20
            bc(iad1+i34)=bc(ii+i34)                                     1d25s20
           end do                                                       1d25s20
           ii=ii+nrow                                                   1d25s20
          end do                                                        1d25s20
          i10=1                                                         1d25s20
         end do                                                         1d25s20
        end if                                                          1d25s20
        call dws_gsumf(bc(ionexb(is)),nwds)                             1d25s20
        ncol=nrow*irefo(isblk1(3,is))                                   3d23s21
        do i=0,ncol-1                                                   3d23s21
         do j=0,nvirt(isblk1(4,is))-1                                   3d23s21
          ji=ionexbt(is)+j+nvirt(isblk1(4,is))*i                        3d23s21
          ij=ionexb(is)+i+ncol*j                                        3d23s21
          bc(ji)=bc(ij)                                                 3d23s21
         end do                                                         3d23s21
        end do                                                          3d23s21
        if(nhere.gt.0)then
         if(isblk1(1,is).eq.isblk1(2,is))then                             12d11s19
          nrow=(nvirt(isblk1(1,is))*(nvirt(isblk1(1,is))+1))/2            12d11s19
         else                                                             12d11s19
          nrow=nvirt(isblk1(1,is))*nvirt(isblk1(2,is))                    12d11s19
         end if                                                           12d11s19
         i3xb(is)=ibcoff                                                  12d11s19
         ibcoff=i3xb(is)+nrow*nhere                                       12d11s19
         call enough('trans2e.  9',bc,ibc)
         ii3x=i3x(is)                                                     12d11s19
         ii3xb=i3xb(is)                                                 7d2s20
         i10=i1s                                                          12d11s19
         i1n=irefo(isblk1(3,is))                                          12d11s19
         do i2=i2s,i2e                                                    12d11s19
          if(i2.eq.i2e)i1n=i1e                                            12d11s19
          nhere1=i1n+1-i10                                                12d11s19
          do i1=i10,i1n                                                   12d11s19
           do i34=0,nrow-1                                                12d11s19
            iad1=ii3xb+i1-i10+nhere1*i34                                  12d11s19
            bc(iad1)=bc(ii3x+i34)                                         12d11s19
           end do                                                         12d11s19
           ii3x=ii3x+nrow                                                 12d11s19
          end do                                                          12d11s19
          ii3xb=ii3xb+nhere1*nrow                                         12d11s19
          i10=1                                                           12d11s19
         end do                                                           12d11s19
        end if                                                           12d11s19
       end do                                                            12d11s19
      end if                                                            12d12s19
      return
      end
