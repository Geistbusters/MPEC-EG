c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine secoi(iptoh,nocc,multh,iumat,ioooo,itmat,ih0mo,ih0new, 4d19s13
     $     ioooo2,nbasdwsc,bc,ibc)                                      11d9s22
      implicit real*8 (a-h,o-z)
c
c     compute transformed intetegrals for second order energy
c     we will compute this processors part to the global h0 and 4o ints,
c     then do a global sum.
c
      include "common.store"
      include "common.hf"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension iptoh(8,8,8),nocc(8),multh(8,8),iumat(8),ipart(6),
     $     ioooo(1),itmat(8),ioooo2(1),nbasdwsc(8)
      ipty=ibcoff
      ibcoff=ipty+nsdlk
      i2ebuf=ibcoff
      n2tot=0
      do is=1,nsdlk
       ncol=nocc(isblk(3,is))*nocc(isblk(4,is))
       if(isblk(1,is).eq.isblk(2,is))then
        nrow=(nocc(isblk(1,is))*(nocc(isblk(1,is))+1))/2
       else
        nrow=nocc(isblk(1,is))*nocc(isblk(2,is))
       end if
       if(nrow*ncol.gt.0)then
        ibc(ipty+is-1)=ibcoff
        ibcoff=ibcoff+nrow*ncol
        n2tot=n2tot+nrow*ncol
        call enough('secoi.  1',bc,ibc)
        do i=0,nrow*ncol-1
         bc(ibc(ipty+is-1)+i)=0d0
        end do
        call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,
     $       mynowprog,il,ih,i1,i2,i3,i4)
        nhere=ih+1-il
        if(nhere.gt.0)then
         do i34=0,nhere-1
          iad1=ibc(ipty+is-1)+nrow*(i34+il-1)
          iad2=ioooo(is)+nrow*i34
          do i12=0,nrow-1
           bc(iad1+i12)=-bc(iad2+i12)
          end do
         end do
        end if
       end if
      end do
      jh0mo=ih0mo                                                       4d12s13
      ih0new=ibcoff                                                     4d12s13
      do is=1,nsymb                                                     4d12s13
       if(nocc(is).gt.0)then                                            4d12s13
        ih0u=ibcoff                                                     4d12s13
        ibcoff=ibcoff+nocc(is)*nocc(is)                                 4d12s13
        n2tot=n2tot+nocc(is)*nocc(is)                                   4d12s13
        itmp=ibcoff                                                     4d12s13
        call ilimts(nbasdwsc(is),1,mynprocg,mynowprog,il,ih,i1,i2,i3,i4)8d31s15
        nhere=ih+1-il                                                   4d12s13
        itmp2=itmp+nhere*nocc(is)                                       4d12s13
        ibcoff=itmp2+nhere*nocc(is)                                     4d12s13
        call enough('secoi.  2',bc,ibc)
        do i=0,nocc(is)*nocc(is)-1                                      4d12s13
         bc(ih0u+i)=0d0                                                 4d12s13
        end do                                                          4d12s13
        if(nhere.gt.0)then                                              4d12s13
         kh0mo=jh0mo+il-1                                               4d12s13
         call dgemm('n','n',nhere,nocc(is),nbasdwsc(is),1d0,bc(kh0mo),  8d31s15
     $       nbasdwsc(is),bc(iumat(is)),nbasdwsc(is),0d0,bc(itmp),nhere,8d31s15
     d' secoi.  1')
         do i=1,nocc(is)
          do j=1,nhere
           iad1=itmp+j-1+nhere*(i-1)
           iad2=itmp2+i-1+nocc(is)*(j-1)
           bc(iad2)=bc(iad1)
          end do
         end do
         kumat=iumat(is)+il-1
         call dgemm('n','n',nocc(is),nocc(is),nhere,1d0,bc(itmp2),
     $        nocc(is),bc(kumat),nbasdwsc(is),0d0,bc(ih0u),nocc(is),    8d31s15
     d' secoi.  2')
         ibcoff=itmp                                                    4d12s13
        end if                                                          4d12s13
       end if                                                           4d12s13
       jh0mo=jh0mo+nbasdwsc(is)*nbasdwsc(is)                            8d31s15
      end do                                                            4d12s13
    1 format('integral type ',4i1)
      do isb=1,nsymb
       do isc=1,nsymb
        isbc=multh(isc,isb)
        do isd=1,nsymb
         isa=multh(isd,isbc)
         ii=iptoh(isd,isc,isb)
         if(ii.gt.0)then
          call ilimts(nocc(isd),nbasdwsc(isc),mynprocg,mynowprog,il,ih, 8d31s15
     $         i1s,i1e,i2s,i2e)
          ncol=ih+1-il
          if(isa.eq.isb)then
           nrow=(nbasdwsc(isa)*(nbasdwsc(isa)+1))/2                     8d31s15
          else
           nrow=nbasdwsc(isa)*nbasdwsc(isb)                             8d31s15
          end if
          if(nrow*ncol.gt.0.and.
     $         nocc(isa)*nocc(isb)*nocc(isc)*nocc(isd).ne.0)then        4d11s13
c
c     J type transformation
c
           do i2=1,nocc(isc)                                            4d11s13
            do i1=1,nocc(isd)                                           4d11s13
             icol=i1+nocc(isd)*(i2-1)                                   4d11s13
             if(icol.ge.il.and.icol.le.ih)then                          4d11s13
              itmp2=ibcoff                                              4d11s13
              itmp3=itmp2+nbasdwsc(isa)*nocc(isb)                       8d31s15
              ibcoff=itmp3+nbasdwsc(isa)*nocc(isb)                      8d31s15
              call enough('secoi.  3',bc,ibc)
              if(isa.eq.isb)then                                        4d11s13
               itmp=ibcoff                                              4d11s13
               ibcoff=itmp+nbasdwsc(isa)*nbasdwsc(isa)                  8d31s15
               call enough('secoi.  4',bc,ibc)
               ij=ihcol(ii)+nrow*(icol-il)                              4d11s13
               do i=1,nbasdwsc(isa)                                     8d31s15
                do j=1,i                                                4d11s13
                 iad1=itmp+j-1+nbasdwsc(isa)*(i-1)                      8d31s15
                 bc(iad1)=bc(ij)                                        4d11s13
                 iad1=itmp+i-1+nbasdwsc(isa)*(j-1)                      8d31s15
                 bc(iad1)=bc(ij)                                        4d11s13
                 ij=ij+1                                                4d11s13
                end do                                                  4d11s13
               end do                                                   4d11s13
              else                                                      4d11s13
               itmp=ihcol(ii)+nrow*(icol-il)                            7d18s14
              end if
              call dgemm('n','n',nbasdwsc(isa),nocc(isb),nbasdwsc(isb), 8d31s15
     $               1d0,bc(itmp),nbasdwsc(isa),bc(iumat(isb)),         8d31s15
     $               nbasdwsc(isb),0d0,bc(itmp2),nbasdwsc(isa),         8d31s15
     d' secoi.  3')
              do i=1,nbasdwsc(isa)                                      8d31s15
               do j=1,nocc(isb)                                         4d11s13
                iad1=itmp2+i-1+nbasdwsc(isa)*(j-1)                      8d31s15
                iad2=itmp3+j-1+nocc(isb)*(i-1)                          4d11s13
                bc(iad2)=bc(iad1)                                       4d11s13
               end do                                                   4d11s13
              end do                                                    4d11s13
              call dgemm('n','n',nocc(isb),nocc(isa),nbasdwsc(isa),     8d31s15
     $             1d0,bc(itmp3),nocc(isb),bc(iumat(isa)),nbasdwsc(isa),8d31s15
     $             0d0,bc(itmp2),nocc(isb),                             4d11s13
     d' secoi.  4')
              do is=1,nsdlk                                             4d11s13
               if(isblk(1,is).eq.isblk(2,is))then                       4d11s13
                nsrow=(nocc(isblk(1,is))*(nocc(isblk(1,is))+1))/2       4d11s13
               else                                                     4d11s13
                nsrow=nocc(isblk(1,is))*nocc(isblk(2,is))               4d11s13
               end if                                                   4d11s13
               if(isa.eq.isblk(1,is).and.isb.eq.isblk(2,is).and.        4d11s13
     $            isd.eq.isblk(4,is).and.isc.eq.isblk(3,is))then        4d11s13
                do i=1,nocc(isb)                                        4d11s13
                 ixx=ibc(ipty+is-1)-1+nsrow*(i2-1+nocc(isc)*(i1-1))     4d11s13
                 if(isa.eq.isb)then                                     4d11s13
                  ixx=ixx+((i*(i-1))/2)                                 4d11s13
                  jtop=i                                                4d11s13
                 else                                                   4d11s13
                  jtop=nocc(isa)                                        4d11s13
                  ixx=ixx+nocc(isa)*(i-1)                               4d11s13
                 end if                                                 4d11s13
                 do j=1,jtop                                            4d11s13
                  iad=itmp2+i-1+nocc(isb)*(j-1)                         7d18s14
                  bc(ixx+j)=bc(ixx+j)+bc(iad)                           4d11s13
                 end do                                                 4d11s13
                end do                                                  4d11s13
               else if(isb.eq.isblk(1,is).and.isa.eq.isblk(2,is).and.   4d11s13
     $            isd.eq.isblk(4,is).and.isc.eq.isblk(3,is))then        4d11s13
                do i=1,nocc(isa)                                        4d11s13
                 ixx=ibc(ipty+is-1)-1+nsrow*(i2-1+nocc(isc)*(i1-1))     4d11s13
                 jtop=nocc(isb)                                         4d11s13
                 ixx=ixx+nocc(isb)*(i-1)                                4d11s13
                 do j=1,jtop                                            4d11s13
                  iad=itmp2+j-1+nocc(isb)*(i-1)                         7d18s14
                  bc(ixx+j)=bc(ixx+j)+bc(iad)                           4d11s13
                 end do                                                 4d11s13
                end do                                                  4d11s13
               end if                                                   4d11s13
               if(isa.eq.isblk(3,is).and.isb.eq.isblk(4,is).and.        4d11s13
     $            isd.eq.isblk(2,is).and.isc.eq.isblk(1,is))then        4d11s13
                do i=1,nocc(isb)                                        4d11s13
                 ixx=ibc(ipty+is-1)-1+nsrow*(-1+nocc(isa)*(i-1))        4d11s13
                 if(isa.eq.isb)then                                     4d11s13
                  if(i2.ge.i1)then
                   ixx=ixx+((i2*(i2-1))/2)+i1                           4d12s13
                   do j=1,nocc(isa)                                     4d12s13
                    iad=itmp2+i-1+nocc(isb)*(j-1)                       7d18s14
                    bc(ixx+j*nsrow)=bc(ixx+j*nsrow)+bc(iad)             4d12s13
                   end do                                               4d12s13
                  end if                                                4d12s13
                 else                                                   4d11s13
                  ixx=ixx+i2+nocc(isc)*(i1-1)                           4d11s13
                  do j=1,nocc(isa)                                       4d11s13
                   iad=itmp2+i-1+nocc(isb)*(j-1)                         4d11s13
                   bc(ixx+j*nsrow)=bc(ixx+j*nsrow)+bc(iad)               4d11s13
                  end do                                                 4d11s13
                 end if                                                 4d11s13
                end do                                                  4d11s13
               else if(isb.eq.isblk(3,is).and.isa.eq.isblk(4,is).and.   4d11s13
     $            isd.eq.isblk(2,is).and.isc.eq.isblk(1,is))then        4d11s13
                do i=1,nocc(isa)                                        4d11s13
                 ixx=ibc(ipty+is-1)-1+nsrow*(-1+nocc(isb)*(i-1))        4d11s13
                 ixx=ixx+i2+nocc(isc)*(i1-1)                            4d11s13
                 do j=1,nocc(isb)                                       4d11s13
                  iad=itmp2+j-1+nocc(isb)*(i-1)                         4d11s13
                  bc(ixx+j*nsrow)=bc(ixx+j*nsrow)+bc(iad)               4d11s13
                 end do                                                 4d11s13
                end do                                                  4d11s13
               end if                                                   4d11s13
              end do                                                    4d11s13
              ibcoff=itmp2                                              4d11s13
             end if                                                     4d11s13
            end do
           end do
c
c     K type transformation
c
c     transform second index
c
           itmp=ibcoff                                                  4d11s13
           nhrow=nocc(isa)*nocc(isb)                                    4d11s13
           itmp3=itmp+nhrow*ncol                                        4d11s13
           ibcoff=itmp3+nhrow*nocc(isd)*nocc(isc)                       4d11s13
           call enough('secoi.  5',bc,ibc)
           do icol=1,ncol                                               4d11s13
            itmp2=ibcoff                                                4d11s13
            ibcoff=itmp2+nocc(isa)*nbasdwsc(isb)                        8d31s15
            call enough('secoi.  6',bc,ibc)
            if(isa.eq.isb)then                                          4d11s13
             do i=1,nocc(isa)                                            4d15s13
              do j=1,i                                                   4d15s13
               iarg1=ihcol(ii)+((i*(i-1))/2)+j-1+nrow*(icol-1)          4d11s13
               iarg2=itmp2+j-1+nocc(isa)*(i-1)                          4d15s13
               bc(iarg2)=bc(iarg1)                                      4d15s13
               iarg2=itmp2+i-1+nocc(isa)*(j-1)                          4d15s13
               bc(iarg2)=bc(iarg1)                                      4d15s13
              end do                                                    4d15s13
             end do                                                     4d15s13
             do i=nocc(isa)+1,nbasdwsc(isa)                             8d31s15
              do j=1,nocc(isa)                                          4d15s13
               iarg1=ihcol(ii)+((i*(i-1))/2)+j-1+nrow*(icol-1)           4d15s13
               iarg2=itmp2+j-1+nocc(isa)*(i-1)                          4d15s13
               bc(iarg2)=bc(iarg1)                                      4d15s13
              end do                                                    4d15s13
             end do                                                     4d15s13
            else                                                        4d11s13
             do i=1,nbasdwsc(isb)                                       8d31s15
              do j=1,nocc(isa)                                          4d11s13
               iarg1=itmp2+j-1+nocc(isa)*(i-1)                          4d11s13
               iarg2=ihcol(ii)+j-1+nbasdwsc(isa)*(i-1)+nrow*(icol-1)    8d31s15
               bc(iarg1)=bc(iarg2)                                      4d11s13
              end do                                                    4d11s13
             end do                                                     4d11s13
            end if                                                      4d11s13
            jtmp=itmp+nhrow*(icol-1)                                    4d11s13
            call dgemm('n','n',nocc(isa),nocc(isb),nbasdwsc(isb),1d0,   8d31s15
     $           bc(itmp2),nocc(isa),bc(itmat(isb)),nbasdwsc(isb),      8d31s15
     $           0d0,bc(jtmp),nocc(isa),                                4d11s13
     d' secoi.  5')
            ibcoff=itmp2                                                4d11s13
           end do                                                       4d11s13
c
c     partial transformation of 3th index
c
           do i=0,nhrow*nocc(isd)*nocc(isc)-1                           4d11s13
            bc(itmp3+i)=0d0                                             4d11s13
           end do                                                       4d11s13
           do i3=1,nocc(isc)                                            4d11s13
            i10=i1s                                                     4d11s13
            i1x=nocc(isd)                                               4d11s13
            jtmp=itmp-1                                                 4d11s13
            do i2=i2s,i2e                                               4d11s13
             targ=bc(itmat(isc)+i2-1+nbasdwsc(isc)*(i3-1))              8d31s15
             if(i2.eq.i2e)i1x=i1e                                       4d11s13
             do i1=i10,i1x                                              4d11s13
              iarg1=itmp3-1+nhrow*(i1-1+nocc(isd)*(i3-1))               4d11s13
              do i4=1,nhrow                                             4d11s13
               bc(iarg1+i4)=bc(iarg1+i4)+targ*bc(jtmp+i4)               4d11s13
              end do                                                    4d11s13
              jtmp=jtmp+nhrow                                           4d11s13
             end do                                                     4d11s13
             i10=1                                                      4d11s13
            end do                                                      4d11s13
           end do                                                       4d11s13
           do is=1,nsdlk                                                4d11s13
            if(isblk(1,is).eq.isblk(2,is))then
             nsrow=(nocc(isblk(1,is))*(nocc(isblk(1,is))+1))/2
            else
             nsrow=nocc(isblk(1,is))*nocc(isblk(2,is))
            end if
            nscol=nocc(isblk(3,is))*nocc(isblk(4,is))
            if(isblk(1,is).eq.isa.and.isblk(2,is).eq.isb.and.
     $         isblk(3,is).eq.isc.and.isblk(4,is).eq.isd.and.is.ne.499)
     $           then
             jtmp3=itmp3                                                4d11s13
             if(isblk(1,is).eq.isblk(2,is))then
              do i4=1,nocc(isc)                                          4d11s13
               do i3=1,nocc(isd)                                         4d11s13
                do i2=1,nocc(isb)                                        4d11s13
                 do i1=1,nocc(isa)                                       4d11s13
                  if(i1.ge.i2)then                                      4d11s13
                   iarg=ibc(ipty+is-1)+i2-1+((i1*(i1-1))/2)+nsrow       4d11s13
     $                 *(i4-1+nocc(isc)*(i3-1))                         4d11s13
                   bc(iarg)=bc(iarg)+bc(jtmp3)                          4d15s13
                  end if                                                4d11s13
                  jtmp3=jtmp3+1                                          4d11s13
                 end do                                                  4d11s13
                end do                                                   4d11s13
               end do                                                    4d11s13
              end do                                                     4d11s13
             else                                                       4d11s13
              do i4=1,nocc(isc)                                          4d11s13
               do i3=1,nocc(isd)                                         4d11s13
                do i2=1,nocc(isb)                                        4d11s13
                 do i1=1,nocc(isa)                                       4d11s13
                  iarg=ibc(ipty+is-1)+i1-1+nocc(isa)*(i2-1+nocc(isb)    4d11s13
     $                 *(i4-1+nocc(isc)*(i3-1)))                         4d11s13
                  bc(iarg)=bc(iarg)+bc(jtmp3)                            4d11s13
                  jtmp3=jtmp3+1                                          4d11s13
                 end do                                                  4d11s13
                end do                                                   4d11s13
               end do                                                    4d11s13
              end do                                                     4d11s13
             end if                                                     4d11s13
            end if
            if(isblk(2,is).eq.isa.and.isblk(1,is).eq.isb.and.
     $         isblk(3,is).eq.isc.and.isblk(4,is).eq.isd.and.is.ne.499)
     $           then
             jtmp3=itmp3                                                4d11s13
             if(isblk(1,is).eq.isblk(2,is))then
              do i4=1,nocc(isc)                                          4d11s13
               do i3=1,nocc(isd)                                         4d11s13
                do i2=1,nocc(isb)                                        4d11s13
                 do i1=1,nocc(isa)                                       4d11s13
                  if(i2.ge.i1)then                                      4d11s13
                   iarg=ibc(ipty+is-1)+i1-1+((i2*(i2-1))/2)+nsrow       4d11s13
     $                 *(i4-1+nocc(isc)*(i3-1))                         4d11s13
                   bc(iarg)=bc(iarg)+bc(jtmp3)                          4d15s13
                  end if                                                4d11s13
                  jtmp3=jtmp3+1                                          4d11s13
                 end do                                                  4d11s13
                end do                                                   4d11s13
               end do                                                    4d11s13
              end do                                                     4d11s13
             else                                                       4d11s13
              do i4=1,nocc(isc)                                          4d11s13
               do i3=1,nocc(isd)                                         4d11s13
                do i2=1,nocc(isb)                                        4d11s13
                 do i1=1,nocc(isa)                                       4d11s13
                  iarg=ibc(ipty+is-1)+i2-1+nocc(isb)*(i1-1+nocc(isa)    4d11s13
     $                 *(i4-1+nocc(isc)*(i3-1)))                         4d11s13
                  bc(iarg)=bc(iarg)+bc(jtmp3)                            4d11s13
                  jtmp3=jtmp3+1                                          4d11s13
                 end do                                                  4d11s13
                end do                                                   4d11s13
               end do                                                    4d11s13
              end do                                                     4d11s13
             end if                                                     4d11s13
            end if
            if(isblk(1,is).eq.isa.and.isblk(2,is).eq.isb.and.
     $         isblk(4,is).eq.isc.and.isblk(3,is).eq.isd.and.is.ne.499)
     $           then
             jtmp3=itmp3                                                4d11s13
             if(isblk(1,is).eq.isblk(2,is))then
              do i4=1,nocc(isc)                                          4d11s13
               do i3=1,nocc(isd)                                         4d11s13
                do i2=1,nocc(isb)                                        4d11s13
                 do i1=1,nocc(isa)                                       4d11s13
                  if(i1.ge.i2)then                                      4d11s13
                   iarg=ibc(ipty+is-1)+i2-1+((i1*(i1-1))/2)+nsrow       4d11s13
     $                 *(i3-1+nocc(isd)*(i4-1))                         4d11s13
                   bc(iarg)=bc(iarg)+bc(jtmp3)                          4d15s13
                  end if                                                4d11s13
                  jtmp3=jtmp3+1                                          4d11s13
                 end do                                                  4d11s13
                end do                                                   4d11s13
               end do                                                    4d11s13
              end do                                                     4d11s13
             else                                                       4d11s13
              do i4=1,nocc(isc)                                          4d11s13
               do i3=1,nocc(isd)                                         4d11s13
                do i2=1,nocc(isb)                                        4d11s13
                 do i1=1,nocc(isa)                                       4d11s13
                  iarg=ibc(ipty+is-1)+i1-1+nocc(isa)*(i2-1+nocc(isb)    4d11s13
     $                 *(i3-1+nocc(isd)*(i4-1)))                         4d11s13
                  bc(iarg)=bc(iarg)+bc(jtmp3)                            4d11s13
                  jtmp3=jtmp3+1                                          4d11s13
                 end do                                                  4d11s13
                end do                                                   4d11s13
               end do                                                    4d11s13
              end do                                                     4d11s13
             end if                                                     4d11s13
            end if
            if(isblk(2,is).eq.isa.and.isblk(1,is).eq.isb.and.
     $         isblk(4,is).eq.isc.and.isblk(3,is).eq.isd.and.is.ne.499)
     $           then
             jtmp3=itmp3                                                4d11s13
             if(isblk(1,is).eq.isblk(2,is))then
              do i4=1,nocc(isc)                                          4d11s13
               do i3=1,nocc(isd)                                         4d11s13
                do i2=1,nocc(isb)                                        4d11s13
                 do i1=1,nocc(isa)                                       4d11s13
                  if(i2.ge.i1)then                                      4d11s13
                   iarg=ibc(ipty+is-1)+i1-1+((i2*(i2-1))/2)+nsrow       4d11s13
     $                 *(i3-1+nocc(isd)*(i4-1))                         4d11s13
                   bc(iarg)=bc(iarg)+bc(jtmp3)                          4d15s13
                  end if                                                4d11s13
                  jtmp3=jtmp3+1                                          4d11s13
                 end do                                                  4d11s13
                end do                                                   4d11s13
               end do                                                    4d11s13
              end do                                                     4d11s13
             else                                                       4d11s13
              do i4=1,nocc(isc)                                          4d11s13
               do i3=1,nocc(isd)                                         4d11s13
                do i2=1,nocc(isb)                                        4d11s13
                 do i1=1,nocc(isa)                                       4d11s13
                  iarg=ibc(ipty+is-1)+i2-1+nocc(isb)*(i1-1+nocc(isa)    4d11s13
     $                 *(i3-1+nocc(isd)*(i4-1)))                         4d11s13
                  bc(iarg)=bc(iarg)+bc(jtmp3)                            4d11s13
                  jtmp3=jtmp3+1                                          4d11s13
                 end do                                                  4d11s13
                end do                                                   4d11s13
               end do                                                    4d11s13
              end do                                                     4d11s13
             end if                                                     4d11s13
            end if
           end do                                                       4d11s13
c
c     transform first index
c
           if(isa.ne.isb)then                                           4d16s13
            do icol=1,ncol                                               4d11s13
             itmp2=ibcoff                                                4d11s13
             ibcoff=itmp2+nocc(isb)*nbasdwsc(isa)                       8d31s15
             call enough('secoi.  7',bc,ibc)
             if(isa.eq.isb)then
              do i=1,nocc(isa)                                           4d11s13
               do j=1,i                                                  4d11s13
                iarg1=ihcol(ii)+j-1+((i*(i-1))/2)+nrow*(icol-1)          4d11s13
                iarg2=itmp2+j-1+nocc(isa)*(i-1)
                bc(iarg2)=bc(iarg1)                                      4d11s13
                iarg2=itmp2+i-1+nocc(isa)*(j-1)
                bc(iarg2)=bc(iarg1)                                      4d11s13
               end do                                                    4d11s13
              end do                                                     4d11s13
              do i=nocc(isa)+1,nbasdwsc(isa)                            8d31s15
               do j=1,nocc(isa)                                          4d11s13
                iarg1=ihcol(ii)+j-1+((i*(i-1))/2)+nrow*(icol-1)          4d11s13
                iarg2=itmp2+j-1+nocc(isa)*(i-1)
                bc(iarg2)=bc(iarg1)                                      4d11s13
               end do                                                    4d11s13
              end do                                                     4d11s13
             else                                                        4d11s13
              do i=1,nocc(isb)                                           4d11s13
               do j=1,nbasdwsc(isa)                                     8d31s15
                iarg1=itmp2+i-1+nocc(isb)*(j-1)                         4d18s13
                iarg2=ihcol(ii)+j-1+nbasdwsc(isa)*(i-1)+nrow*(icol-1)   8d31s15
                bc(iarg1)=bc(iarg2)                                      4d11s13
               end do                                                    4d11s13
              end do                                                     4d11s13
             end if                                                      4d11s13
             jtmp=itmp+nhrow*(icol-1)                                    4d11s13
             call dgemm('n','n',nocc(isb),nocc(isa),nbasdwsc(isa),1d0,  8d31s15
     $            bc(itmp2),nocc(isb),bc(itmat(isa)),nbasdwsc(isa),     8d31s15
     $            0d0,bc(jtmp),nocc(isb),                                4d11s13
     d' secoi.  6')
             ibcoff=itmp2                                                4d11s13
            end do                                                       4d11s13
c
c     partial transformation of 3th index
c
            do i=0,nhrow*nocc(isd)*nocc(isc)-1                           4d11s13
             bc(itmp3+i)=0d0                                             4d11s13
            end do                                                       4d11s13
            do i3=1,nocc(isc)                                            4d11s13
             i10=i1s                                                     4d11s13
             i1x=nocc(isd)                                               4d11s13
             jtmp=itmp-1                                                 4d11s13
             do i2=i2s,i2e                                               4d11s13
              targ=bc(itmat(isc)+i2-1+nbasdwsc(isc)*(i3-1))             8d31s15
              if(i2.eq.i2e)i1x=i1e                                       4d11s13
              do i1=i10,i1x                                              4d11s13
               iarg1=itmp3-1+nhrow*(i1-1+nocc(isd)*(i3-1))               4d11s13
               do i4=1,nhrow                                             4d11s13
                bc(iarg1+i4)=bc(iarg1+i4)+targ*bc(jtmp+i4)               4d11s13
               end do                                                    4d11s13
               jtmp=jtmp+nhrow                                           4d11s13
              end do                                                     4d11s13
              i10=1                                                      4d11s13
             end do                                                      4d11s13
            end do                                                       4d11s13
            do is=1,nsdlk                                                4d11s13
             if(isblk(1,is).eq.isblk(2,is))then
              nsrow=(nocc(isblk(1,is))*(nocc(isblk(1,is))+1))/2
             else
              nsrow=nocc(isblk(1,is))*nocc(isblk(2,is))
             end if
             nscol=nocc(isblk(3,is))*nocc(isblk(4,is))
             if(isblk(1,is).eq.isa.and.isblk(2,is).eq.isb.and.
     $         isblk(3,is).eq.isc.and.isblk(4,is).eq.isd)then
              jtmp3=itmp3                                                4d11s13
              if(isblk(1,is).eq.isblk(2,is))then
               do i4=1,nocc(isc)                                          4d11s13
                do i3=1,nocc(isd)                                         4d11s13
                 do i1=1,nocc(isa)                                       4d11s13
                  do i2=1,nocc(isb)                                        4d11s13
                   if(i1.ge.i2)then                                      4d11s13
                    iarg=ibc(ipty+is-1)+i2-1+((i1*(i1-1))/2)+nsrow       4d11s13
     $                  *(i4-1+nocc(isc)*(i3-1))                         4d11s13
                    bc(iarg)=bc(iarg)+bc(jtmp3)                          4d15s13
                   end if                                                4d11s13
                   jtmp3=jtmp3+1                                          4d11s13
                  end do                                                  4d11s13
                 end do                                                   4d11s13
                end do                                                    4d11s13
               end do                                                     4d11s13
              else                                                       4d11s13
               do i4=1,nocc(isc)                                          4d11s13
                do i3=1,nocc(isd)                                         4d11s13
                 do i1=1,nocc(isa)                                       4d11s13
                  do i2=1,nocc(isb)                                        4d11s13
                   iarg=ibc(ipty+is-1)+i1-1+nocc(isa)*(i2-1+nocc(isb)    4d11s13
     $                 *(i4-1+nocc(isc)*(i3-1)))                         4d11s13
                   bc(iarg)=bc(iarg)+bc(jtmp3)                            4d11s13
                   jtmp3=jtmp3+1                                          4d11s13
                  end do                                                  4d11s13
                 end do                                                   4d11s13
                end do                                                    4d11s13
               end do                                                     4d11s13
              end if                                                     4d11s13
             end if
             if(isblk(2,is).eq.isa.and.isblk(1,is).eq.isb.and.
     $            isblk(3,is).eq.isc.and.isblk(4,is).eq.isd)then
              jtmp3=itmp3                                                4d11s13
              if(isblk(1,is).eq.isblk(2,is))then
               do i4=1,nocc(isc)                                          4d11s13
                do i3=1,nocc(isd)                                         4d11s13
                 do i1=1,nocc(isa)                                       4d11s13
                  do i2=1,nocc(isb)                                        4d11s13
                   if(i1.ge.i2)then                                      4d11s13
                    iarg=ibc(ipty+is-1)+i2-1+((i1*(i1-1))/2)+nsrow       4d11s13
     $                 *(i4-1+nocc(isc)*(i3-1))                         4d11s13
                    bc(iarg)=bc(iarg)+bc(jtmp3)                          4d15s13
                   end if                                                4d11s13
                   jtmp3=jtmp3+1                                          4d11s13
                  end do                                                  4d11s13
                 end do                                                   4d11s13
                end do                                                    4d11s13
               end do                                                     4d11s13
              else                                                       4d11s13
               do i4=1,nocc(isc)                                          4d11s13
                do i3=1,nocc(isd)                                         4d11s13
                 do i1=1,nocc(isa)                                       4d11s13
                  do i2=1,nocc(isb)                                        4d11s13
                   iarg=ibc(ipty+is-1)+i2-1+nocc(isb)*(i1-1+nocc(isa)    4d11s13
     $                 *(i4-1+nocc(isc)*(i3-1)))                         4d11s13
                   bc(iarg)=bc(iarg)+bc(jtmp3)                            4d11s13
                   jtmp3=jtmp3+1                                          4d11s13
                  end do                                                  4d11s13
                 end do                                                   4d11s13
                end do                                                    4d11s13
               end do                                                     4d11s13
              end if                                                     4d11s13
             end if
             if(isblk(1,is).eq.isa.and.isblk(2,is).eq.isb.and.
     $         isblk(4,is).eq.isc.and.isblk(3,is).eq.isd)then
              jtmp3=itmp3                                                4d11s13
              if(isblk(1,is).eq.isblk(2,is))then
               do i4=1,nocc(isc)                                          4d11s13
                do i3=1,nocc(isd)                                         4d11s13
                 do i1=1,nocc(isa)                                       4d11s13
                  do i2=1,nocc(isb)                                        4d11s13
                   if(i1.ge.i2)then                                      4d11s13
                    iarg=ibc(ipty+is-1)+i2-1+((i1*(i1-1))/2)+nsrow       4d11s13
     $                 *(i3-1+nocc(isd)*(i4-1))                         4d11s13
                    bc(iarg)=bc(iarg)+bc(jtmp3)                          4d15s13
                   end if                                                4d11s13
                   jtmp3=jtmp3+1                                          4d11s13
                  end do                                                  4d11s13
                 end do                                                   4d11s13
                end do                                                    4d11s13
               end do                                                     4d11s13
              else                                                       4d11s13
               do i4=1,nocc(isc)                                          4d11s13
                do i3=1,nocc(isd)                                         4d11s13
                 do i1=1,nocc(isa)                                       4d11s13
                  do i2=1,nocc(isb)                                        4d11s13
                   iarg=ibc(ipty+is-1)+i1-1+nocc(isa)*(i2-1+nocc(isb)    4d11s13
     $                 *(i3-1+nocc(isd)*(i4-1)))                         4d11s13
                   bc(iarg)=bc(iarg)+bc(jtmp3)                            4d11s13
                   jtmp3=jtmp3+1                                          4d11s13
                  end do                                                  4d11s13
                 end do                                                   4d11s13
                end do                                                    4d11s13
               end do                                                     4d11s13
              end if                                                     4d11s13
             end if
             if(isblk(2,is).eq.isa.and.isblk(1,is).eq.isb.and.
     $         isblk(4,is).eq.isc.and.isblk(3,is).eq.isd)then
              jtmp3=itmp3                                                4d11s13
              if(isblk(1,is).eq.isblk(2,is))then
               do i4=1,nocc(isc)                                          4d11s13
                do i3=1,nocc(isd)                                         4d11s13
                 do i1=1,nocc(isa)                                       4d11s13
                  do i2=1,nocc(isb)                                        4d11s13
                   if(i1.ge.i2)then                                      4d11s13
                    iarg=ibc(ipty+is-1)+i2-1+((i1*(i1-1))/2)+nsrow       4d11s13
     $                 *(i3-1+nocc(isd)*(i4-1))                         4d11s13
                    bc(iarg)=bc(iarg)+bc(jtmp3)                          4d15s13
                   end if                                                4d11s13
                   jtmp3=jtmp3+1                                          4d11s13
                  end do                                                  4d11s13
                 end do                                                   4d11s13
                end do                                                    4d11s13
               end do                                                     4d11s13
              else                                                       4d11s13
               do i4=1,nocc(isc)                                          4d11s13
                do i3=1,nocc(isd)                                         4d11s13
                 do i1=1,nocc(isa)                                       4d11s13
                  do i2=1,nocc(isb)                                        4d11s13
                   iarg=ibc(ipty+is-1)+i2-1+nocc(isb)*(i1-1+nocc(isa)    4d11s13
     $                  *(i3-1+nocc(isd)*(i4-1)))                         4d11s13
                   bc(iarg)=bc(iarg)+bc(jtmp3)                            4d11s13
                   jtmp3=jtmp3+1                                          4d11s13
                  end do                                                  4d11s13
                 end do                                                   4d11s13
                end do                                                    4d11s13
               end do                                                     4d11s13
              end if                                                     4d11s13
             end if
            end do                                                       4d11s13
           end if                                                       4d16s13
           ibcoff=itmp
          end if
         end if
        end do
       end do
      end do
      call dws_gsumf(bc(i2ebuf),n2tot)                                  4d11s13
      do is=1,nsdlk
       ncol=nocc(isblk(3,is))*nocc(isblk(4,is))
       if(isblk(1,is).eq.isblk(2,is))then
        nrow=(nocc(isblk(1,is))*(nocc(isblk(1,is))+1))/2
       else
        nrow=nocc(isblk(1,is))*nocc(isblk(2,is))
       end if
       ioooo2(is)=ibcoff                                                4d19s13
       if(nrow*ncol.gt.0)then
        call ilimts(nocc(isblk(3,is)),nocc(isblk(4,is)),mynprocg,       4d19s13
     $       mynowprog,il,ih,i1,i2,i3,i4)                               4d19s13
        nhere=ih+1-il                                                   4d19s13
        if(nhere.gt.0)then
         ibcoff=ioooo2(is)+nhere*nrow                                    4d19s13
         call enough('secoi.  8',bc,ibc)
         ioff=ibc(ipty+is-1)+nrow*(il-1)                                4d19s13
         do i=0,nhere*nrow-1                                            4d19s13
          bc(ioooo2(is)+i)=bc(ioff+i)                                   4d19s13
         end do                                                         4d19s13
        end if                                                          4d19s13
       end if
      end do
      ih0u=ih0new
      do is=1,nsymb
       if(nocc(is).gt.0)then
        sym=0d0
        do i=1,nocc(is)
         do j=1,i-1
          iad1=ih0u+j-1+nocc(is)*(i-1)
          iad2=ih0u+i-1+nocc(is)*(j-1)
          sym=sym+(bc(iad1)-bc(iad2))**2
         end do
        end do
        isym=(nocc(is)*(nocc(is)+1))/2
        sym=sqrt(sym/dfloat(isym))
        ih0u=ih0u+nocc(is)*nocc(is)
       end if
      end do
      return
      end
