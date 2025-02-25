c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine postden(nsdlk,isblk,nsdlk1,isblk1,nsdlkk,isblkk,irefo, 8d6s24
     $     nvirt,id4o,ioooo,id1x,ionex,jmden,jmats,kmden,kmats,id3x,i3x,8d6s24
     $     ldebug,nsing,ndoub,dot4v,shift,ee,igoal,bc,ibc)              8d6s24
      implicit real*8 (a-h,o-z)                                         8d6s24
      include "common.store"                                            8d6s24
      logical ldebug                                                    8d6s24
      dimension isblk(4,*),irefo(*),id4o(*),ioooo(*),isblk1(4,*),       8d6s24
     $     id1x(*),ionex(*),jmden(*),jmats(*),isblkk(4,*),kmden(*),     8d6s24
     $     kmats(*),id3x(*),i3x(*),nvirt(*)                             8d6s24
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      do ib=1,nsdlk
       if(isblk(1,ib).ne.isblk(2,ib))then                               5d5s23
        do ibb=ib+1,nsdlk                                               5d5s23
         if(isblk(1,ib).eq.isblk(1,ibb).and.isblk(2,ib).eq.isblk(2,ibb) 5d5s23
     $        .and.isblk(3,ib).eq.isblk(4,ibb))then                     5d5s23
          if(min(irefo(isblk(1,ib)),irefo(isblk(2,ib)),
     $        irefo(isblk(3,ib)),irefo(isblk(4,ib))).gt.0)then
           ncol=irefo(isblk(3,ib))*irefo(isblk(4,ib))                   5d5s23
           nn=irefo(isblk(1,ib))*irefo(isblk(2,ib))                     5d5s23
           do i4=0,irefo(isblk(4,ib))-1                                  5d5s23
            do i3=0,irefo(isblk(3,ib))-1                                 5d5s23
             i34=id4o(ib)+nn*(i3+irefo(isblk(3,ib))*i4)                  5d5s23
             i43=id4o(ibb)+nn*(i4+irefo(isblk(4,ib))*i3)                 5d5s23
             do i12=0,nn-1                                               5d5s23
              sum=0.5d0*(bc(i12+i34)+bc(i12+i43))                       7d17s23
              bc(i12+i34)=sum                                            5d5s23
              bc(i12+i43)=sum                                            5d5s23
             end do                                                      5d5s23
            end do                                                       5d5s23
           end do                                                        5d5s23
          end if
         end if                                                         5d5s23
        end do                                                          5d5s23
       end if                                                           5d5s23
      end do
      do ib=1,nsdlk                                                     8d19s14
       n1=isblk(1,ib)                                                   8d19s14
       n2=isblk(2,ib)                                                   8d19s14
       n3=isblk(3,ib)
       n4=isblk(4,ib)                                                   8d19s14
       if(n1.eq.n2)then
        nnn=(irefo(n1)*(irefo(n1)+1))/2
       else
        nnn=irefo(n1)*irefo(n2)
       end if
       if(n3.eq.n4)then
        mmm=(irefo(n3)*(irefo(n3)+1))/2
       else
        mmm=irefo(n3)*irefo(n4)
       end if
       if(nnn*mmm.gt.0)then
        if(n1.eq.n3.and.n2.eq.n4)then                                   7d6s23
         do j=0,nnn-1                                                    8d20s14
          do k=0,j-1                                                     8d20s14
           iad1=id4o(ib)+k+nnn*j                                        8d20s14
           iad2=id4o(ib)+j+nnn*k                                        8d20s14
           avg=0.5d0*(bc(iad1)+bc(iad2))                                 8d20s14
           bc(iad1)=avg                                                  8d20s14
           bc(iad2)=avg                                                  8d20s14
          end do                                                         8d20s14
         end do                                                          8d20s14
        end if                                                           8d20s14
        if(.not.(n1.eq.n2.and.n3.eq.n4.and.n1.eq.n3))then                9d18s14
         do io=ib+1,nsdlk                                                 9d18s14
          if(isblk(1,io).eq.n3.and.isblk(2,io).eq.n4.and.                 9d18s14
     $       isblk(3,io).eq.n1.and.isblk(4,io).eq.n2)then               9d18s14
           do j=0,nnn-1                                                    9d18s14
            do k=0,mmm-1                                                   9d18s14
             iad1=id4o(ib)+j+nnn*k                                       9d18s14
             iad2=id4o(io)+k+mmm*j                                       9d18s14
             avg=0.5d0*(bc(iad1)+bc(iad2))                                 9d18s14
             bc(iad1)=avg                                                  9d18s14
             bc(iad2)=avg                                                  9d18s14
            end do                                                         9d18s14
           end do                                                          9d18s14
          else if(isblk(1,io).eq.n4.and.isblk(2,io).eq.n3.and.                 9d18s14
     $       isblk(3,io).eq.n2.and.isblk(4,io).eq.n1)then               9d18s14
           do i4=0,irefo(n4)-1                                          7d25s23
            do i3=0,irefo(n3)-1                                         7d25s23
             i34=i3+irefo(n3)*i4                                        7d25s23
             i43=i4+irefo(n4)*i3                                        7d25s23
             do i2=0,irefo(n2)-1                                        7d25s23
              do i1=0,irefo(n1)-1                                       7d25s23
               i12=i1+irefo(n1)*i2                                      7d25s23
               i21=i2+irefo(n2)*i1                                      7d25s23
               iad1=id4o(ib)+i12+nnn*i34                                7d25s23
               iad2=id4o(io)+i43+mmm*i21                                7d25s23
               avg=0.5d0*(bc(iad1)+bc(iad2))                            7d25s23
               bc(iad1)=avg                                                  9d18s14
               bc(iad2)=avg                                                  9d18s14
              end do                                                    7d25s23
             end do                                                     7d25s23
            end do                                                         9d18s14
           end do                                                          9d18s14
          end if                                                           9d18s14
         end do                                                            9d18s14
        end if                                                            9d18s14
       end if                                                            9d18s14
      end do                                                            9d18s14
      if(ldebug)write(6,*)('4o '),bc(igoal)
      t4o=0d0                                                           7d17s23
      do is=1,nsdlk                                                     5d2s23
       if(isblk(1,is).eq.isblk(2,is))then                               7d28s22
        nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2              7d28s22
        ncol=(irefo(isblk(3,is))*(irefo(isblk(3,is))+1))/2              7d28s22
       else                                                             7d28s22
        nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                      7d28s22
        ncol=irefo(isblk(3,is))*irefo(isblk(4,is))                       7d28s22
       end if                                                           7d28s22
       if(min(nrow,ncol).gt.0)then                                      5d2s23
        if(ldebug)then                                                  6d20s24
         write(6,*)('4o for '),(isblk(j,is),j=1,4),id4o(is)              5d23s23
         call prntm2(bc(id4o(is)),nrow,ncol,nrow)                        5d2s23
        end if                                                          6d20s24
        call ilimts(ncol,1,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)    7d17s23
        nhere=ih+1-il                                                   7d17s23
        trace=0d0                                                       7d17s23
        if(nhere.gt.0.and.ldebug)then                                   6d20s24
         do icol=il,ih                                                  7d17s23
          icolm=icol-1                                                  7d17s23
          iad1=id4o(is)+nrow*icolm                                      7d17s23
          iad2=ioooo(is)+nrow*icolm                                     7d17s23
          do j=0,nrow-1                                                 7d17s23
           term=bc(iad1+j)*bc(iad2+j)                                   7d17s23
           trace=trace+term                                             7d17s23
          end do                                                        7d17s23
         end do                                                         7d17s23
        end if                                                          7d17s23
        if(ldebug)then
         call dws_gsumf(trace,1)                                         7d17s23
         t4o=t4o+trace                                                   7d17s23
        end if
       end if                                                           5d2s23
      end do                                                            5d2s23
      if(ldebug)then                                                    6d20s24
       write(6,*)('final 4o trace '),t4o                                 7d17s23
       ee=ee+t4o                                                         7d17s23
      end if                                                            6d20s24
      if(nsing.gt.0)then                                                7d20s23
       if(ldebug)then                                                   6d20s24
        write(6,*)('1x '),bc(igoal),nsdlk1
        t1x=0d0                                                           7d17s23
       end if                                                           6d20s24
       do is=1,nsdlk1                                                    7d28s22
        if(isblk1(1,is).eq.isblk1(2,is))then                             7d28s22
         nrow1=(irefo(isblk1(1,is))*(irefo(isblk1(1,is))+1))/2           7d28s22
        else                                                             7d28s22
         nrow1=irefo(isblk1(1,is))*irefo(isblk1(2,is))                   7d28s22
        end if                                                           7d28s22
        ncol=irefo(isblk1(3,is))*nvirt(isblk1(4,is))                     5d18s23
        if(min(nrow1,ncol).gt.0)then
         call dws_gsumf(bc(id1x(is)),nrow1*ncol)                         7d6s23
         call ilimts(nrow1,1,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)   7d17s23
         nhere=ih+1-il                                                   7d17s23
         trace=0d0                                                       7d17s23
         if(nhere.gt.0.and.ldebug)then                                  6d20s24
          do icol=il,ih                                                  7d17s23
           icolm=icol-1                                                  7d17s23
           iad1=id1x(is)+ncol*icolm                                      7d17s23
           iad2=ionex(is)+ncol*icolm                                     7d17s23
           do j=0,ncol-1                                                 7d17s23
            trace=trace+bc(iad1+j)*bc(iad2+j)                            7d17s23
           end do                                                        7d17s23
          end do                                                         7d17s23
         end if                                                          7d17s23
         if(ldebug)then                                                 6d20s24
          call dws_gsumf(trace,1)                                         7d17s23
          t1x=t1x+trace                                                   7d17s23
         end if                                                         6d20s24
        end if
       end do                                                            7d28s22
       if(ldebug)then                                                   6d20s24
        write(6,*)('final trace for onex '),t1x
        ee=ee+t1x                                                        10d18s23
       end if                                                           6d20s24
      end if                                                            7d20s23
      if(max(nsing,ndoub).gt.0.and.ldebug)then                          6d20s24
       write(6,*)('J '),bc(igoal)
       tj=0d0                                                            7d17s23
       do is=1,nsdlk
        if(isblk(1,is).eq.isblk(2,is))then
         nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2
        else
         nrow=irefo(isblk(1,is))*irefo(isblk(2,is))
        end if
        call ilimts(nvirt(isblk(3,is)),nvirt(isblk(4,is)),mynprocg,
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)
        ncol=ih+1-il
        trace=0d0                                                        7d17s23
        if(min(nrow,ncol).gt.0)then
         do i=0,ncol*nrow-1                                              7d17s23
          trace=trace+bc(jmden(is)+i)*bc(jmats(is)+i)                    7d17s23
         end do                                                          7d17s23
        end if
        call dws_gsumf(trace,1)                                          7d17s23
        tj=tj+trace                                                      7d17s23
       end do
       write(6,*)('final tj trace '),tj
       ee=ee+tj                                                          7d17s23
       write(6,*)('K '),bc(igoal),igoal
       tk=0d0
       do is=1,nsdlkk
        nrow=irefo(isblkk(1,is))*irefo(isblkk(2,is))
        call ilimts(nvirt(isblkk(3,is)),nvirt(isblkk(4,is)),mynprocg,
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)
        ncol=ih+1-il
        trace=0d0                                                        7d17s23
        if(min(nrow,ncol).gt.0)then
         write(6,*)('K for '),(isblkk(j,is),j=1,4),kmden(is),('is = '),
     $        is
         call prntm2(bc(kmden(is)),ncol,nrow,ncol)
         do i=0,ncol*nrow-1                                              7d17s23
          trace=trace+bc(kmden(is)+i)*bc(kmats(is)+i)                    7d17s23
         end do                                                          7d17s23
        end if
        call dws_gsumf(trace,1)                                          7d17s23
        tk=tk+trace                                                      7d17s23
       end do
       write(6,*)('final trace for tk '),tk
       ee=ee+tk                                                          7d17s23
      end if                                                            7d20s23
      if(ldebug)then                                                    6d20s24
       write(6,*)('3x '),bc(igoal)
       t3x=0d0                                                           10d18s23
       do is=1,nsdlk1                                                    7d28s22
        if(isblk1(1,is).eq.isblk1(2,is))then                             7d28s22
         nrow3=(nvirt(isblk1(1,is))*(nvirt(isblk1(1,is))+1))/2           7d28s22
        else                                                             7d28s22
         nrow3=nvirt(isblk1(1,is))*nvirt(isblk1(2,is))                   7d28s22
        end if                                                           7d28s22
        call ilimts(irefo(isblk1(3,is)),nvirt(isblk1(4,is)),mynprocg,    7d28s22
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            7d28s22
        ncol=ih+1-il                                                     7d28s22
        if(min(nrow3,ncol).gt.0)then                                     10d18s23
         do i=0,nrow3*ncol-1                                             10d18s23
          t3x=t3x+bc(id3x(is)+i)*bc(i3x(is)+i)
         end do
        end if                                                           10d18s23
       end do                                                            7d28s22
       call dws_gsumf(t3x,1)                                             10d18s23
       write(6,*)('t3x: '),t3x
       ee=ee+t3x+dot4v
       write(6,*)('energy from trace '),ee+shift,ee
      end if                                                            6d20s24
      return
      end
