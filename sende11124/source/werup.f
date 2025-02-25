c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine werup(iumat,itmat,ibmat10,iamat10,jmats,kmats,iqk,ipk,
     $     iden1,ih0mo,numpro,nowpro,igmat,noc,ionex,ioooo,nbasdwsc,    8d31s15
     $     nvirtc,bc,ibc)                                               11d9s22
      implicit real*8 (a-h,o-z)
      include "common.store"
      include "common.hf"
c
c     werner update ...
c
      dimension iumat(1),itmat(1),jmats(1),kmats(1),noc(8),igmat(8,8),
     $     ipk(1),iqk(1),iden1(1),ionex(1),ioooo(1),nbasdwsc(8),        8d31s15
     $     nvirtc(8)                                                    8d31s15
      ibcoffo=ibcoff                                                    8d4s14
c
      iamat=iamat10
      ibmat=ibmat10
      nbmat=0                                                           7d10s13
      ih0u=ih0mo
      nocx=noc(1)                                                       8d1s14
      do isb=2,nsymb                                                    8d1s14
       nocx=max(nocx,noc(isb))                                          8d1s14
      end do                                                            8d1s14
      do isb=1,nsymb
       if(noc(isb)*nbasdwsc(isb).gt.0)then                              8d31s15
        nbmat=nbmat+noc(isb)*nbasdwsc(isb)                              8d31s15
        ibmat1=ibmat                                                     4d19s13
        ibmat2=ibmat1+noc(isb)*noc(isb)                                  9d9s04
        do i=0,noc(isb)*nbasdwsc(isb)-1                                 8d31s15
         bc(ibmat+i)=0d0
        end do
c
c     copy a...
c
        if(nowpro.eq.0)then
         do i=1,noc(isb)                                                  9d9s04
          j1=ibmat1-1+(i-1)*noc(isb)                                      9d9s04
          j2=iamat-1+(i-1)*nbasdwsc(isb)                                8d31s15
          do j=1,noc(isb)                                                 9d9s04
           bc(j1+j)=bc(j2+j)                                              9d9s04
          end do                                                          9d9s04
         end do                                                           9d9s04
         do j=1,noc(isb)
          j1=ibmat2-1+(j-1)*nvirtc(isb)                                    1d6s11
          j2=iamat+noc(isb)-1+(j-1)*nbasdwsc(isb)                       8d31s15
          do i=1,nvirtc(isb)                                               1d6s11
           bc(j1+i)=bc(j2+i)
          end do
         end do
        end if
c
c     this procs part of hTD...
c
        call ilimts(nbasdwsc(isb),1,numpro,nowpro,il,ih,i1,i2,i3,i4)    8d31s15
        nhere=ih+1-il                                                   4d19s13
        if(nhere.gt.0)then                                              4d19s13
         ih=ih0u+nbasdwsc(isb)*(il-1)                                   8d31s15
         it=itmat(isb)+il-1                                             4d19s13
         itmp=ibcoff                                                    4d19s13
         ibcoff=itmp+nbasdwsc(isb)*noc(isb)                             8d31s15
         call enough('werup.  1',bc,ibc)
         if(min(nbasdwsc(isb),noc(isb),nhere).gt.0)then                 6d12s19
         call dgemm('n','n',nbasdwsc(isb),noc(isb),nhere,1d0,bc(ih),    8d31s15
     $        nbasdwsc(isb),bc(it),nbasdwsc(isb),0d0,bc(itmp),          8d31s15
     $        nbasdwsc(isb),                                            8d31s15
     d' werup.  1')
         call dgemm('n','n',noc(isb),noc(isb),noc(isb),1d0,bc(itmp),    4d19s13
     $        nbasdwsc(isb),bc(iden1(isb)),noc(isb),1d0,bc(ibmat1),     8d31s15
     $        noc(isb),                                                 4d19s13
     d' werup.  2')
         end if                                                         6d12s19
         jtmp=itmp+noc(isb)                                             4d19s13
         if(min(nvirtc(isb),noc(isb),nbasdwsc(isb)).gt.0)then           6d12s19
         call dgemm('n','n',nvirtc(isb),noc(isb),noc(isb),1d0,bc(jtmp),    4d19s13
     $        nbasdwsc(isb),bc(iden1(isb)),noc(isb),1d0,bc(ibmat2),     8d31s15
     $        nvirtc(isb),                                               4d19s13
     d' werup.  3')
         end if                                                         6d12s19
         ibcoff=itmp                                                    4d19s13
        end if
c
c     now sum over other symmetry for J*T ...
c
        do iso=1,nsymb                                                  8d5s14
         if(noc(iso).gt.0)then                                          8d5s14
          itmp=ibcoff                                                   8d5s14
          nrow=(noc(iso)*(noc(iso)+1))/2                                8d5s14
          ibcoff=itmp+nrow*noc(isb)*nbasdwsc(isb)                       8d31s15
          call enough('werup.  2',bc,ibc)
          do i=0,nrow*noc(isb)*nbasdwsc(isb)-1                          8d31s15
           bc(itmp+i)=0d0                                               8d5s14
          end do                                                        8d5s14
          do is=1,nsdlk                                                 8d5s14
           if(isblk(1,is).eq.iso.and.isblk(2,is).eq.iso.and.            8d5s14
     $          isblk(3,is).eq.isb)then                                 8d5s14
            call ilimts(noc(isb),noc(isb),numpro,                       8d5s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d5s14
            nhere=ih+1-il                                               8d5s14
            if(nhere.gt.0)then                                          8d5s14
             i10=i1s                                                    8d5s14
             i1n=noc(isb)                                               8d5s14
             joooo=ioooo(is)-1                                          8d5s14
             do i2=i2s,i2e                                              8d5s14
              if(i2.eq.i2e)i1n=i1e                                      8d5s14
              do i1=i10,i1n                                             8d5s14
               do i5=1,noc(isb)                                         8d5s14
                iad1=itmat(isb)+i2-1+nbasdwsc(isb)*(i5-1)               8d31s15
                do i34=1,nrow                                            8d5s14
                 iad2=itmp+i34-1+nrow*(i1-1+nbasdwsc(isb)*(i5-1))       8d31s15
                 iad3=joooo+i34                                         8d5s14
                 bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                    8d5s14
                end do                                                  8d5s14
               end do                                                   8d5s14
               joooo=joooo+nrow                                         8d5s14
              end do                                                    8d5s14
              i10=1                                                     8d5s14
             end do                                                     8d5s14
            end if                                                      8d5s14
           end if                                                       8d5s14
          end do                                                        8d5s14
          do is=1,nsdlk1                                                8d5s14
           if(isblk1(1,is).eq.iso.and.isblk1(2,is).eq.iso.and.          8d5s14
     $          isblk1(3,is).eq.isb)then                                8d5s14
            call ilimts(noc(isb),nvirtc(isb),numpro,nowpro,il,ih,i1s,i1e8d31s15
     $          ,i2s,i2e)                                               8d31s15
            nhere=ih+1-il                                                 7d12s13
            if(nhere.gt.0)then                                          8d5s14
             i10=i1s                                                    8d5s14
             i1n=noc(isb)                                               8d5s14
             jonex=ionex(is)-1                                          8d5s14
             do i2=i2s,i2e                                              8d5s14
              i2p=i2+noc(isb)
              if(i2.eq.i2e)i1n=i1e                                      8d5s14
              do i1=i10,i1n                                             8d5s14
               do i5=1,noc(isb)                                         8d5s14
                iad1=itmat(isb)+i2p-1+nbasdwsc(isb)*(i5-1)              8d31s15
                do i34=1,nrow                                           8d5s14
                 iad2=itmp+i34-1+nrow*(i1-1+nbasdwsc(isb)*(i5-1))       8d31s15
                 iad3=jonex+i34                                         8d5s14
                 bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                    8d5s14
                end do                                                  8d5s14
               end do                                                   8d5s14
               do i5=1,noc(isb)                                         8d5s14
                iad1=itmat(isb)+i1-1+nbasdwsc(isb)*(i5-1)               8d31s15
                do i34=1,nrow                                           8d5s14
                 iad2=itmp+i34-1+nrow*(i2p-1+nbasdwsc(isb)*(i5-1))      8d31s15
                 iad3=jonex+i34                                         8d5s14
                 bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                    8d5s14
                end do                                                  8d5s14
               end do                                                   8d5s14
               jonex=jonex+nrow                                         8d5s14
              end do                                                    8d5s14
              i10=1                                                     8d5s14
             end do                                                     8d5s14
            end if                                                      8d5s14
           end if                                                       8d5s14
          end do                                                        8d5s14
          do is=1,nsdlk                                                 8d5s14
           if(isblk(1,is).eq.iso.and.isblk(2,is).eq.iso.and.            8d5s14
     $          isblk(3,is).eq.isb)then                                 8d5s14
            call ilimts(nvirtc(isb),nvirtc(isb),numpro,nowpro,il,ih,i1s,  8d5s14
     $         i1e,i2s,i2e)                                             8d5s14
            nhere=ih+1-il                                                 8d4s14
            if(nhere.gt.0)then                                            8d4s14
             i10=i1s                                                    8d5s14
             i1n=nvirtc(isb)                                             8d5s14
             jms=jmats(is)-1                                            8d5s14
             do i2=i2s,i2e                                              8d5s14
              i2p=i2+noc(isb)                                           8d5s14
              if(i2.eq.i2e)i1n=i1e                                      8d5s14
              do i1=i10,i1n                                             8d5s14
               i1p=i1+noc(isb)                                          8d5s14
               do i5=1,noc(isb)                                         8d5s14
                iad1=itmat(isb)+i2p-1+nbasdwsc(isb)*(i5-1)              8d31s15
                do i34=1,nrow                                           8d5s14
                 iad2=itmp+i34-1+nrow*(i1p-1+nbasdwsc(isb)*(i5-1))      8d31s15
                 iad3=jms+i34                                           8d5s14
                 bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                    8d5s14
                end do                                                  8d5s14
               end do                                                   8d5s14
               jms=jms+nrow                                             8d5s14
              end do                                                    8d5s14
              i10=1                                                     8d5s14
             end do                                                     8d5s14
            end if                                                      8d5s14
           end if                                                       8d5s14
          end do                                                        8d5s14
          do i2=1,noc(isb)                                              8d5s14
           do i1=1,noc(isb)                                             8d5s14
            iad1=ibmat1+i1-1+noc(isb)*(i2-1)                            8d5s14
            sum=0d0                                                     8d5s14
            do l=1,noc(iso)                                               8d5s14
             iad2=itmp+((l*(l+1))/2)-1+nrow*(i1-1+nbasdwsc(isb)*(i2-1)) 8d31s15
             sum=sum+bc(iad2)                                           8d5s14
            end do                                                      8d5s14
            sum=sum*4d0                                                 8d5s14
            if(iso.eq.isb)then                                          8d5s14
             sum2=0d0                                                   8d5s14
             do l=1,noc(isb)                                            8d5s14
              ix=max(l,i2)                                              8d5s14
              in=min(l,i2)                                              8d5s14
              iad2=itmp+((ix*(ix-1))/2)+in-1+nrow                       8d5s14
     $             *(i1-1+nbasdwsc(isb)*(l-1))                          8d31s15
              sum2=sum2+bc(iad2)                                        8d5s14
             end do                                                     8d5s14
             sum=sum-2d0*sum2                                           8d5s14
            end if                                                      8d5s14
            bc(iad1)=bc(iad1)+sum                                       8d5s14
           end do                                                       8d5s14
           do i1=1,nvirtc(isb)                                             8d5s14
            i1p=i1+noc(isb)                                             8d5s14
            iad1=ibmat2+i1-1+nvirtc(isb)*(i2-1)                            8d5s14
            sum=0d0                                                     8d5s14
            do l=1,noc(iso)                                               8d5s14
             iad2=itmp+((l*(l+1))/2)-1+nrow*(i1p-1+nbasdwsc(isb)*(i2-1))8d31s15
             sum=sum+bc(iad2)                                           8d5s14
            end do                                                      8d5s14
            sum=sum*4d0                                                 8d5s14
            if(iso.eq.isb)then                                          8d5s14
             sum2=0d0                                                   8d5s14
             do l=1,noc(isb)                                            8d5s14
              ix=max(l,i2)                                              8d5s14
              in=min(l,i2)                                              8d5s14
              iad2=itmp+((ix*(ix-1))/2)+in-1+nrow                       8d5s14
     $             *(i1p-1+nbasdwsc(isb)*(l-1))                         8d31s15
              sum2=sum2+bc(iad2)                                        8d5s14
             end do                                                     8d5s14
             sum=sum-2d0*sum2                                           8d5s14
            end if                                                      8d5s14
            bc(iad1)=bc(iad1)+sum                                       8d5s14
           end do                                                       8d5s14
          end do                                                        8d5s14
          ibcoff=itmp                                                   8d5s14
          if(iso.ne.isb)then                                            8d5s14
           itmp=ibcoff                                                  8d5s14
           nrow=noc(iso)*noc(isb)                                       8d5s14
           ibcoff=itmp+nrow*noc(iso)*nbasdwsc(isb)                      8d31s15
           call enough('werup.  3',bc,ibc)
           do i=0,nrow*noc(iso)*nbasdwsc(isb)-1                         8d31s15
            bc(itmp+i)=0d0                                               8d5s14
           end do                                                        8d5s14
           nmatch=0                                                     8d6s14
           do is=1,nsdlk                                                 8d5s14
            if(isblk(1,is).eq.iso.and.isblk(2,is).eq.isb.and.            8d5s14
     $           isblk(3,is).eq.iso)then                                 8d5s14
             nmatch=1
             call ilimts(noc(iso),noc(isb),numpro,                       8d5s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d5s14
             nhere=ih+1-il                                               8d5s14
             if(nhere.gt.0)then                                          8d5s14
              i10=i1s                                                    8d5s14
              i1n=noc(iso)                                               8d5s14
              joooo=ioooo(is)-1                                          8d5s14
              do i2=i2s,i2e                                              8d5s14
               if(i2.eq.i2e)i1n=i1e                                      8d5s14
               do i1=i10,i1n                                             8d5s14
                do i5=1,noc(iso)                                         8d5s14
                 iad1=itmat(iso)+i1-1+nbasdwsc(iso)*(i5-1)              8d31s15
                 do i34=1,nrow                                          8d6s14
                  iad2=itmp+i34-1+nrow*(i2-1+nbasdwsc(isb)*(i5-1))      8d31s15
                  iad3=joooo+i34                                         8d5s14
                  bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                    8d5s14
                 end do                                                  8d5s14
                end do                                                   8d5s14
                joooo=joooo+nrow                                         8d5s14
               end do                                                    8d5s14
               i10=1                                                     8d5s14
              end do                                                     8d5s14
             end if                                                      8d5s14
            else if(isblk(1,is).eq.isb.and.isblk(2,is).eq.iso.and.            8d5s14
     $           isblk(3,is).eq.isb)then                                 8d5s14
             nmatch=1
             call ilimts(noc(isb),noc(iso),numpro,                       8d5s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d5s14
             nhere=ih+1-il                                               8d5s14
             if(nhere.gt.0)then                                          8d5s14
              i10=i1s                                                    8d5s14
              i1n=noc(isb)                                               8d5s14
              joooo=ioooo(is)-1                                          8d5s14
              do i2=i2s,i2e                                              8d5s14
               if(i2.eq.i2e)i1n=i1e                                      8d5s14
               do i1=i10,i1n                                             8d5s14
                do i5=1,noc(iso)                                         8d5s14
                 iad1=itmat(iso)+i2-1+nbasdwsc(iso)*(i5-1)              8d31s15
                 i34=1
                 do i4=1,noc(iso)                                       8d6s14
                  do i3=1,noc(isb)                                      8d6s14
                   iad2=itmp+i4-1+noc(iso)*(i3-1+noc(isb)               8d6s14
     $                  *(i1-1+nbasdwsc(isb)*(i5-1)))                   8d31s15
                   iad3=joooo+i34                                         8d5s14
                   bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                    8d5s14
                   i34=i34+1                                            8d6s14
                  end do                                                8d6s14
                 end do                                                  8d5s14
                end do                                                   8d5s14
                joooo=joooo+nrow                                         8d5s14
               end do                                                    8d5s14
               i10=1                                                     8d5s14
              end do                                                     8d5s14
             end if                                                      8d5s14
            end if                                                       8d5s14
           end do                                                        8d5s14
           if(nmatch.eq.0)then
            write(6,*)('failed to match oooo for J in werup')
            call dws_sync
            call dws_finalize
            stop
           end if
           nmatch=0
           do is=1,nsdlk1                                                8d5s14
            if(isblk1(1,is).eq.iso.and.isblk1(2,is).eq.isb.and.          8d5s14
     $           isblk1(3,is).eq.isb)then                                8d5s14
             nmatch=1
             call ilimts(noc(isb),nvirtc(iso),numpro,nowpro,il,ih,i1s,   8d6s14
     $         i1e,i2s,i2e)                                             8d6s14
             nhere=ih+1-il                                                 7d12s13
             if(nhere.gt.0)then                                          8d5s14
              i10=i1s                                                    8d5s14
              i1n=noc(isb)                                               8d5s14
              jonex=ionex(is)-1                                          8d5s14
              do i2=i2s,i2e                                              8d5s14
               i2p=i2+noc(iso)
               if(i2.eq.i2e)i1n=i1e                                      8d5s14
               do i1=i10,i1n                                             8d5s14
                do i5=1,noc(iso)                                         8d5s14
                 iad1=itmat(iso)+i2p-1+nbasdwsc(iso)*(i5-1)             8d31s15
                 do i34=1,nrow                                           8d5s14
                  iad2=itmp+i34-1+nrow*(i1-1+nbasdwsc(isb)*(i5-1))      8d31s15
                  iad3=jonex+i34                                         8d5s14
                  bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                    8d5s14
                 end do                                                  8d5s14
                end do                                                   8d5s14
                jonex=jonex+nrow                                         8d5s14
               end do                                                    8d5s14
               i10=1                                                     8d5s14
              end do                                                     8d5s14
             end if                                                      8d5s14
            else if(isblk1(1,is).eq.isb.and.isblk1(2,is).eq.iso.and.          8d5s14
     $           isblk1(3,is).eq.isb)then                                8d5s14
             nmatch=1
             call ilimts(noc(isb),nvirtc(iso),numpro,nowpro,il,ih,i1s,   8d6s14
     $         i1e,i2s,i2e)                                             8d6s14
             nhere=ih+1-il                                                 7d12s13
             if(nhere.gt.0)then                                          8d5s14
              i10=i1s                                                    8d5s14
              i1n=noc(isb)                                               8d5s14
              jonex=ionex(is)-1                                          8d5s14
              do i2=i2s,i2e                                              8d5s14
               i2p=i2+noc(iso)
               if(i2.eq.i2e)i1n=i1e                                      8d5s14
               do i1=i10,i1n                                             8d5s14
                do i5=1,noc(iso)                                         8d5s14
                 iad1=itmat(iso)+i2p-1+nbasdwsc(iso)*(i5-1)             8d31s15
                 i34=1
                 do i4=1,noc(iso)
                  do i3=1,noc(isb)
                   iad2=itmp+i4-1+noc(iso)*(i3-1+noc(isb)               8d6s14
     $                  *(i1-1+nbasdwsc(isb)*(i5-1)))                   8d31s15
                   iad3=jonex+i34                                         8d5s14
                   bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                    8d5s14
                   i34=i34+1                                            8d6s14
                  end do                                                8d6s14
                 end do                                                  8d5s14
                end do                                                   8d5s14
                jonex=jonex+nrow                                         8d5s14
               end do                                                    8d5s14
               i10=1                                                     8d5s14
              end do                                                     8d5s14
             end if                                                      8d5s14
            end if                                                       8d5s14
           end do                                                        8d5s14
           if(nmatch.eq.0)then
            write(6,*)('failed to match onex for J in werup')
            call dws_sync
            call dws_finalize
            stop
           end if
           nmatch=0
           do is=1,nsdlk1                                                8d5s14
            if(isblk1(1,is).eq.iso.and.isblk1(2,is).eq.isb.and.          8d5s14
     $           isblk1(3,is).eq.iso)then                                8d5s14
             nmatch=1
             call ilimts(noc(iso),nvirtc(isb),numpro,nowpro,il,ih,i1s,   8d6s14
     $         i1e,i2s,i2e)                                             8d6s14
             nhere=ih+1-il                                                 7d12s13
             if(nhere.gt.0)then                                          8d5s14
              i10=i1s                                                    8d5s14
              i1n=noc(iso)                                               8d5s14
              jonex=ionex(is)-1                                          8d5s14
              do i2=i2s,i2e                                              8d5s14
               i2p=i2+noc(isb)
               if(i2.eq.i2e)i1n=i1e                                      8d5s14
               do i1=i10,i1n                                             8d5s14
                do i5=1,noc(iso)                                         8d5s14
                 iad1=itmat(iso)+i1-1+nbasdwsc(iso)*(i5-1)              8d31s15
                 do i34=1,nrow                                           8d5s14
                  iad2=itmp+i34-1+nrow*(i2p-1+nbasdwsc(isb)*(i5-1))     8d31s15
                  iad3=jonex+i34                                         8d5s14
                  bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                    8d5s14
                 end do                                                  8d5s14
                end do                                                   8d5s14
                jonex=jonex+nrow                                         8d5s14
               end do                                                    8d5s14
               i10=1                                                     8d5s14
              end do                                                     8d5s14
             end if                                                      8d5s14
            else if(isblk1(1,is).eq.isb.and.isblk1(2,is).eq.iso.and.          8d5s14
     $           isblk1(3,is).eq.iso)then                                8d5s14
             nmatch=1
             call ilimts(noc(iso),nvirtc(isb),numpro,nowpro,il,ih,i1s,   8d6s14
     $         i1e,i2s,i2e)                                             8d6s14
             nhere=ih+1-il                                                 7d12s13
             if(nhere.gt.0)then                                          8d5s14
              i10=i1s                                                    8d5s14
              i1n=noc(iso)                                               8d5s14
              jonex=ionex(is)-1                                          8d5s14
              do i2=i2s,i2e                                              8d5s14
               i2p=i2+noc(isb)
               if(i2.eq.i2e)i1n=i1e                                      8d5s14
               do i1=i10,i1n                                             8d5s14
                do i5=1,noc(iso)                                         8d5s14
                 iad1=itmat(iso)+i1-1+nbasdwsc(iso)*(i5-1)              8d31s15
                 i34=1                                                  8d6s14
                 do i4=1,noc(iso)                                       8d6s14
                  do i3=1,noc(isb)                                      8d6s14
                   iad2=itmp+i4-1+noc(iso)*(i3-1+noc(isb)               8d6s14
     $                  *(i2p-1+nbasdwsc(isb)*(i5-1)))                  8d31s15
                   iad3=jonex+i34                                         8d5s14
                   bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                    8d5s14
                   i34=i34+1                                            8d6s14
                  end do                                                8d6s14
                 end do                                                  8d5s14
                end do                                                   8d5s14
                jonex=jonex+nrow                                         8d5s14
               end do                                                    8d5s14
               i10=1                                                     8d5s14
              end do                                                     8d5s14
             end if                                                      8d5s14
            end if                                                       8d5s14
           end do                                                        8d5s14
           if(nmatch.eq.0)then
            write(6,*)('failed to match onexb for J in werup')
            call dws_sync
            call dws_finalize
            stop
           end if
           nmatch=0
           do is=1,nsdlk                                                 8d5s14
            if(isblk(1,is).eq.isb.and.isblk(2,is).eq.iso.and.            8d5s14
     $          isblk(3,is).eq.isb)then                                 8d5s14
             nmatch=1
             call ilimts(nvirtc(isb),nvirtc(iso),numpro,nowpro,il,ih,i1s8d31s15
     $            ,i1e,i2s,i2e)                                         8d31s15
             nhere=ih+1-il                                                 8d4s14
             if(nhere.gt.0)then                                            8d4s14
              i10=i1s                                                    8d5s14
              i1n=nvirtc(isb)                                             8d5s14
              jms=jmats(is)-1                                            8d5s14
              do i2=i2s,i2e                                              8d5s14
               i2p=i2+noc(iso)                                           8d5s14
               if(i2.eq.i2e)i1n=i1e                                      8d5s14
               do i1=i10,i1n                                             8d5s14
                i1p=i1+noc(isb)                                          8d5s14
                do i5=1,noc(iso)                                         8d5s14
                 iad1=itmat(iso)+i2p-1+nbasdwsc(iso)*(i5-1)             8d31s15
                 i34=1
                 do i4=1,noc(iso)
                  do i3=1,noc(isb)
                   iad2=itmp+i4-1+noc(iso)*(i3-1+noc(isb)               8d6s14
     $                  *(i1p-1+nbasdwsc(isb)*(i5-1)))                  8d31s15
                   iad3=jms+i34                                           8d5s14
                   bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                    8d5s14
                   i34=i34+1                                            8d6s14
                  end do                                                8d6s14
                 end do                                                  8d5s14
                end do                                                   8d5s14
                jms=jms+nrow                                             8d5s14
               end do                                                    8d5s14
               i10=1                                                     8d5s14
              end do                                                     8d5s14
             end if                                                      8d5s14
            else if(isblk(1,is).eq.iso.and.isblk(2,is).eq.isb.and.            8d5s14
     $          isblk(3,is).eq.iso)then                                 8d5s14
             nmatch=1
             call ilimts(nvirtc(iso),nvirtc(isb),numpro,nowpro,il,ih,i1s8d31s15
     $            ,i1e,i2s,i2e)                                         8d31s15
             nhere=ih+1-il                                                 8d4s14
             if(nhere.gt.0)then                                            8d4s14
              i10=i1s                                                    8d5s14
              i1n=nvirtc(iso)                                             8d5s14
              jms=jmats(is)-1                                            8d5s14
              do i2=i2s,i2e                                              8d5s14
               i2p=i2+noc(isb)                                           8d5s14
               if(i2.eq.i2e)i1n=i1e                                      8d5s14
               do i1=i10,i1n                                             8d5s14
                i1p=i1+noc(iso)                                          8d5s14
                do i5=1,noc(iso)                                         8d5s14
                 iad1=itmat(iso)+i1p-1+nbasdwsc(iso)*(i5-1)             8d31s15
                 do i34=1,nrow                                           8d5s14
                  iad2=itmp+i34-1+nrow*(i2p-1+nbasdwsc(isb)*(i5-1))     8d31s15
                  iad3=jms+i34                                           8d5s14
                  bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                    8d5s14
                 end do                                                  8d5s14
                end do                                                   8d5s14
                jms=jms+nrow                                             8d5s14
               end do                                                    8d5s14
               i10=1                                                     8d5s14
              end do                                                     8d5s14
             end if                                                      8d5s14
            end if                                                       8d5s14
           end do                                                        8d5s14
           if(nmatch.eq.0)then
            write(6,*)('failed to match jmat for J in werup')
            call dws_sync
            call dws_finalize
            stop
           end if
           do i2=1,noc(isb)
            do i1=1,noc(isb)                                             8d5s14
             iad1=ibmat1+i1-1+noc(isb)*(i2-1)                            8d5s14
             sum=0d0                                                     8d5s14
             do l=1,noc(iso)                                            8d5s14
              iad2=itmp+l-1+noc(iso)*(i2-1+noc(isb)
     $             *(i1-1+nbasdwsc(isb)*(l-1)))                         8d31s15
              sum=sum+bc(iad2)                                          8d5s14
             end do                                                     8d5s14
             bc(iad1)=bc(iad1)-2d0*sum                                  8d6s14
            end do                                                       8d5s14
           end do                                                        8d5s14
           do i2=1,noc(isb)
            do i1=1,nvirtc(isb)                                             8d5s14
             i1p=i1+noc(isb)                                             8d5s14
             iad1=ibmat2+i1-1+nvirtc(isb)*(i2-1)                            8d5s14
             sum=0d0                                                     8d5s14
             do l=1,noc(iso)                                            8d5s14
              iad2=itmp+l-1+noc(iso)*(i2-1+noc(isb)
     $             *(i1p-1+nbasdwsc(isb)*(l-1)))                        8d31s15
              sum=sum+bc(iad2)                                          8d5s14
             end do                                                     8d5s14
             bc(iad1)=bc(iad1)-2d0*sum                                  8d6s14
            end do                                                       8d5s14
           end do                                                        8d5s14
           ibcoff=itmp                                                   8d5s14
          end if                                                        8d5s14
         end if                                                         8d5s14
        end do                                                          8d5s14
c
c     now sum over other symmetry for K*T ...
c
        do iso=1,nsymb                                                  8d5s14
         if(noc(iso).gt.0)then                                          8d5s14
          itmp=ibcoff                                                   8d5s14
          nrow=noc(iso)*noc(iso)
          ibcoff=itmp+nrow*noc(isb)*nbasdwsc(isb)                       8d31s15
          call enough('werup.  4',bc,ibc)
          do i=0,nrow*noc(isb)*nbasdwsc(isb)-1                          8d31s15
           bc(itmp+i)=0d0                                               8d5s14
          end do                                                        8d5s14
          do is=1,nsdlk                                                 8d5s14
           if(isblk(1,is).eq.iso.and.isblk(2,is).eq.iso.and.            8d5s14
     $          isblk(3,is).eq.isb)then                                 8d5s14
            call ilimts(noc(isb),noc(isb),numpro,                       8d5s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d5s14
            nhere=ih+1-il                                               8d5s14
            if(nhere.gt.0)then                                          8d5s14
             i10=i1s                                                    8d5s14
             i1n=noc(isb)                                               8d5s14
             joooo=ioooo(is)-1                                          8d5s14
             nrowo=(noc(iso)*(noc(iso)+1))/2
             do i2=i2s,i2e                                              8d5s14
              if(i2.eq.i2e)i1n=i1e                                      8d5s14
              do i1=i10,i1n                                             8d5s14
               do i5=1,noc(iso)                                         8d5s14
                do i4=1,noc(iso)                                        8d6s14
                 iad1=itmat(iso)+i4-1+nbasdwsc(iso)*(i5-1)              8d31s15
                 do i3=1,noc(iso)                                       8d6s14
                  in=min(i3,i4)
                  ix=max(i3,i4)
                  inx=((ix*(ix-1))/2)+in
                  iad2=itmp+i3-1+noc(iso)*(i5-1+noc(iso)
     $                 *(i1-1+nbasdwsc(isb)*(i2-1)))                    8d31s15
                  iad3=joooo+inx                                        8d6s14
                  bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                    8d5s14
                 end do                                                 8d6s14
                end do                                                  8d5s14
               end do                                                   8d5s14
               joooo=joooo+nrowo                                        8d6s14
              end do                                                    8d5s14
              i10=1                                                     8d5s14
             end do                                                     8d5s14
            end if                                                      8d5s14
           end if                                                       8d5s14
          end do                                                        8d5s14
          nmatch=0
          do is=1,nsdlk1                                                 8d5s14
           if(isblk1(1,is).eq.iso.and.isblk1(2,is).eq.iso.and.            8d5s14
     $          isblk1(3,is).eq.isb)then                                 8d5s14
            call ilimts(noc(isb),nvirtc(isb),numpro,                       8d5s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d5s14
            nhere=ih+1-il                                               8d5s14
            if(nhere.gt.0)then                                          8d5s14
             i10=i1s                                                    8d5s14
             i1n=noc(isb)                                               8d5s14
             jonex=ionex(is)-1                                          8d5s14
             nrowo=(noc(iso)*(noc(iso)+1))/2
             do i2=i2s,i2e                                              8d5s14
              i2p=i2+noc(isb)
              if(i2.eq.i2e)i1n=i1e                                      8d5s14
              do i1=i10,i1n                                             8d5s14
               do i5=1,noc(iso)                                         8d5s14
                do i4=1,noc(iso)                                        8d6s14
                 iad1=itmat(iso)+i4-1+nbasdwsc(iso)*(i5-1)              8d31s15
                 do i3=1,noc(iso)                                       8d6s14
                  in=min(i3,i4)
                  ix=max(i3,i4)
                  inx=((ix*(ix-1))/2)+in
                  iad2=itmp+i3-1+noc(iso)*(i5-1+noc(iso)
     $                 *(i2p-1+nbasdwsc(isb)*(i1-1)))                   8d31s15
                  iad3=jonex+inx                                        8d6s14
                  bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                    8d5s14
                 end do                                                 8d6s14
                end do                                                  8d5s14
               end do                                                   8d5s14
               jonex=jonex+nrowo                                        8d6s14
              end do                                                    8d5s14
              i10=1                                                     8d5s14
             end do                                                     8d5s14
            end if                                                      8d5s14
           end if                                                       8d5s14
           if(isblk1(1,is).eq.isb.and.isblk1(2,is).eq.isb.and.            8d5s14
     $          isblk1(4,is).eq.iso)then                                 8d5s14
            nmatch=1
            call ilimts(noc(iso),nvirtc(iso),numpro,                       8d5s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d5s14
            nhere=ih+1-il                                               8d5s14
            if(nhere.gt.0)then                                          8d5s14
             i10=i1s                                                    8d5s14
             i1n=noc(iso)                                               8d5s14
             jonex=ionex(is)-1                                          8d5s14
             nrowo=(noc(isb)*(noc(isb)+1))/2
             do i2=i2s,i2e                                              8d5s14
              i2p=i2+noc(iso)
              if(i2.eq.i2e)i1n=i1e                                      8d5s14
              do i1=i10,i1n                                             8d5s14
               do i5=1,noc(iso)                                         8d5s14
                iad1=itmat(iso)+i2p-1+nbasdwsc(iso)*(i5-1)              8d31s15
                do i4=1,noc(isb)                                        8d6s14
                 do i3=1,noc(isb)                                       8d6s14
                  in=min(i3,i4)
                  ix=max(i3,i4)
                  inx=((ix*(ix-1))/2)+in
                  iad2=itmp+i1-1+noc(iso)*(i5-1+noc(iso)
     $                 *(i3-1+nbasdwsc(isb)*(i4-1)))                    8d31s15
                  iad3=jonex+inx                                        8d6s14
                  bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                    8d5s14
                 end do                                                 8d6s14
                end do                                                  8d5s14
               end do                                                   8d5s14
               jonex=jonex+nrowo                                        8d6s14
              end do                                                    8d5s14
              i10=1                                                     8d5s14
             end do                                                     8d5s14
            end if                                                      8d5s14
           end if                                                       8d5s14
          end do                                                        8d5s14
          if(nmatch.eq.0)then
           write(6,*)('fail to match onexb in werup '),isb,iso
           call dws_sync
           call dws_finalize
          end if
          imatch=0
          do is=1,nsdlkk                                                 8d6s14
           if(isblkk(1,is).eq.isb.and.isblkk(2,is).eq.iso.and.
     $          isblkk(3,is).eq.iso)then                                 8d6s14
            imatch=1
            call ilimts(nvirtc(iso),nvirtc(isb),numpro,                       8d5s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d5s14
            nhere=ih+1-il                                               8d5s14
            if(nhere.gt.0)then                                          8d5s14
             i10=i1s
             i1n=nvirtc(iso)
             kk=kmats(is)-1
             nrowo=noc(isb)*noc(iso)
             do i2=i2s,i2e
              i2p=i2+noc(isb)
              if(i2.eq.i2e)i1n=i1e
              do i1=i10,i1n
               i1p=i1+noc(iso)
               do i5=1,noc(iso)
                iad1=itmat(iso)+i1p-1+nbasdwsc(iso)*(i5-1)              8d31s15
                do i4=1,noc(iso)
                 do i3=1,noc(isb)
                  iad2=itmp+i4-1+noc(iso)*(i5-1+noc(iso)
     $                 *(i2p-1+nbasdwsc(isb)*(i3-1)))                   8d31s15
                  iad3=kk+i3+noc(isb)*(i4-1)
                  bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)
                 end do
                end do
               end do
               kk=kk+nrowo
              end do
              i10=1
             end do
            end if
           else if(isblkk(1,is).eq.iso.and.isblkk(2,is).eq.isb.and.
     $          isblkk(3,is).eq.isb)then                                 8d6s14
            imatch=1
            call ilimts(nvirtc(isb),nvirtc(iso),numpro,                       8d5s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d5s14
            nhere=ih+1-il                                               8d5s14
            if(nhere.gt.0)then                                          8d5s14
             i10=i1s
             i1n=nvirtc(isb)
             kk=kmats(is)-1
             nrowo=noc(isb)*noc(iso)
             do i2=i2s,i2e
              i2p=i2+noc(iso)
              if(i2.eq.i2e)i1n=i1e
              do i1=i10,i1n
               i1p=i1+noc(isb)
               do i5=1,noc(iso)
                iad1=itmat(iso)+i2p-1+nbasdwsc(iso)*(i5-1)              8d31s15
                do i4=1,noc(isb)
                 do i3=1,noc(iso)
                  iad2=itmp+i3-1+noc(iso)*(i5-1+noc(iso)
     $                 *(i1p-1+nbasdwsc(isb)*(i4-1)))                   8d31s15
                  iad3=kk+i3+noc(iso)*(i4-1)
                  bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)
                 end do
                end do
               end do
               kk=kk+nrowo
              end do
              i10=1
             end do
            end if
           end if                                                       8d6s14
          end do                                                        8d6s14
          if(imatch.eq.0)then
           write(6,*)('missed kmats for '),isb,iso
           write(6,*)('in werup ')
           call dws_sync
           call dws_finalize
           stop
          end if
          do i2=1,noc(isb)
           do i1=1,noc(isb)
            iad1=ibmat1+i1-1+noc(isb)*(i2-1)
            sum=0d0
            do l=1,noc(iso)
             iad2=itmp+l-1+noc(iso)*(l-1+noc(iso)
     $            *(i1-1+nbasdwsc(isb)*(i2-1)))                         8d31s15
             sum=sum+bc(iad2)
            end do
            sum=sum*8d0
            if(iso.eq.isb)then
             sum2=0d0
             do l=1,noc(iso)
              iad2=itmp+i2-1+noc(iso)*(l-1+noc(iso)
     $             *(i1-1+nbasdwsc(isb)*(l-1)))                         8d31s15
              iad3=itmp+l-1+noc(iso)*(i2-1+noc(iso)
     $             *(i1-1+nbasdwsc(isb)*(l-1)))                         8d31s15
              sum2=sum2+bc(iad2)+bc(iad3)
             end do
             sum=sum-2d0*sum2
            end if
            bc(iad1)=bc(iad1)+sum
           end do
           do i1=1,nvirtc(isb)
            iad1=ibmat2+i1-1+nvirtc(isb)*(i2-1)
            i1p=i1+noc(isb)
            sum=0d0
            do l=1,noc(iso)
             iad2=itmp+l-1+noc(iso)*(l-1+noc(iso)
     $            *(i1p-1+nbasdwsc(isb)*(i2-1)))                        8d31s15
             sum=sum+bc(iad2)
            end do
            sum=sum*8d0
            if(iso.eq.isb)then
             sum2=0d0
             do l=1,noc(iso)
              iad2=itmp+i2-1+noc(iso)*(l-1+noc(iso)
     $             *(i1p-1+nbasdwsc(isb)*(l-1)))                        8d31s15
              iad3=itmp+l-1+noc(iso)*(i2-1+noc(iso)
     $             *(i1p-1+nbasdwsc(isb)*(l-1)))                        8d31s15
              sum2=sum2+bc(iad2)+bc(iad3)
             end do
             sum=sum-2d0*sum2
            end if
            bc(iad1)=bc(iad1)+sum
           end do
          end do
          ibcoff=itmp                                                   8d6s14
          if(isb.ne.iso)then                                            8d6s14
c
           nrow=nbasdwsc(isb)*noc(iso)                                  8d31s15
           itmp=ibcoff                                                  8d6s14
           ibcoff=itmp+nrow*noc(iso)*noc(isb)                           8d6s14
           call enough('werup.  5',bc,ibc)
           do i=0,nrow*noc(iso)*noc(isb)-1                              8d6s14
            bc(itmp+i)=0d0                                              8d6s14
           end do                                                       8d6s14
           imatch=0
           do is=1,nsdlk                                                8d6s14
            if(isblk(1,is).eq.isb.and.isblk(2,is).eq.iso.and.           8d6s14
     $           isblk(3,is).eq.isb)then                                8d6s14
             imatch=1                                                   8d6s14
             call ilimts(noc(isb),noc(iso),numpro,                      8d6s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d6s14
             nhere=ih+1-il                                              8d6s14
             if(nhere.gt.0)then                                         8d6s14
              i10=i1s                                                   8d6s14
              i1n=noc(isb)                                              8d6s14
              joooo=ioooo(is)-1                                         8d6s14
              nrowo=noc(isb)*noc(iso)                                   8d6s14
              do i2=i2s,i2e                                             8d6s14
               if(i2.eq.i2e)i1n=i1e                                     8d6s14
               do i1=i10,i1n                                            8d6s14
                do i5=1,noc(iso)                                        8d6s14
                 iad1=itmat(iso)+i2-1+nbasdwsc(iso)*(i5-1)              8d31s15
                 do i4=1,noc(iso)                                       8d6s14
                  do i3=1,noc(isb)                                      8d6s14
                   iad2=itmp+i3-1+nbasdwsc(isb)*(i4-1+noc(iso)*         8d31s15
     $                  (i1-1+noc(isb)*(i5-1)))                         8d6s14
                   iad3=joooo+i3+noc(isb)*(i4-1)                        8d6s14
                   bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                  8d6s14
                  end do                                                8d6s14
                 end do                                                 8d6s14
                end do                                                  8d6s14
                joooo=joooo+nrowo                                       8d6s14
               end do                                                   8d6s14
               i10=1                                                    8d6s14
              end do                                                    8d6s14
             end if                                                     8d6s14
            else if(isblk(1,is).eq.iso.and.isblk(2,is).eq.isb.and.           8d6s14
     $           isblk(3,is).eq.iso)then                                8d6s14
             imatch=1                                                   8d6s14
             call ilimts(noc(iso),noc(isb),numpro,                      8d6s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d6s14
             nhere=ih+1-il                                              8d6s14
             if(nhere.gt.0)then                                         8d6s14
              i10=i1s                                                   8d6s14
              i1n=noc(iso)                                              8d6s14
              joooo=ioooo(is)-1                                         8d6s14
              nrowo=noc(isb)*noc(iso)                                   8d6s14
              do i2=i2s,i2e                                             8d6s14
               if(i2.eq.i2e)i1n=i1e                                     8d6s14
               do i1=i10,i1n                                            8d6s14
                do i5=1,noc(iso)                                        8d6s14
                 iad1=itmat(iso)+i1-1+nbasdwsc(iso)*(i5-1)              8d31s15
                 do i4=1,noc(isb)                                       8d6s14
                  do i3=1,noc(iso)                                      8d6s14
                   iad2=itmp+i4-1+nbasdwsc(isb)*(i3-1+noc(iso)*         8d31s15
     $                  (i2-1+noc(isb)*(i5-1)))                         8d6s14
                   iad3=joooo+i3+noc(iso)*(i4-1)                        8d6s14
                   bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                  8d6s14
                  end do                                                8d6s14
                 end do                                                 8d6s14
                end do                                                  8d6s14
                joooo=joooo+nrowo                                       8d6s14
               end do                                                   8d6s14
               i10=1                                                    8d6s14
              end do                                                    8d6s14
             end if                                                     8d6s14
            end if                                                      8d6s14
           end do                                                       8d6s14
           if(imatch.eq.0)then
            write(6,*)('no match for oooo in werup: '),isb,iso
            call dws_sync
            call dws_finalize
            stop
           end if
           imatch=0                                                     8d6s14
           imatchb=0
           do is=1,nsdlk1                                               8d6s14
            if(isblk1(1,is).eq.isb.and.isblk1(2,is).eq.iso.and.
     $           isblk1(3,is).eq.isb)then
             imatch=1                                                   8d6s14
             call ilimts(noc(isb),nvirtc(iso),numpro,                      8d6s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d6s14
             nhere=ih+1-il                                              8d6s14
             if(nhere.gt.0)then                                         8d6s14
              i10=i1s                                                   8d6s14
              i1n=noc(isb)                                              8d6s14
              jonex=ionex(is)-1                                         8d6s14
              nrowo=noc(isb)*noc(iso)                                   8d6s14
              do i2=i2s,i2e                                             8d6s14
               i2p=i2+noc(iso)                                          8d6s14
               if(i2.eq.i2e)i1n=i1e                                     8d6s14
               do i1=i10,i1n                                            8d6s14
                do i5=1,noc(iso)                                        8d6s14
                 iad1=itmat(iso)+i2p-1+nbasdwsc(iso)*(i5-1)             8d31s15
                 do i4=1,noc(iso)                                       8d6s14
                  do i3=1,noc(isb)                                      8d6s14
                   iad2=itmp+i3-1+nbasdwsc(isb)*(i4-1+noc(iso)          8d31s15
     $                  *(i1-1+noc(isb)*(i5-1)))                        8d6s14
                   iad3=jonex+i3+noc(isb)*(i4-1)                        8d6s14
                   bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                  8d6s14
                  end do                                                8d6s14
                 end do                                                 8d6s14
                end do                                                  8d6s14
                jonex=jonex+nrowo                                       8d6s14
               end do                                                   8d6s14
               i10=1                                                    8d6s14
              end do                                                    8d6s14
             end if                                                     8d6s14
            else if(isblk1(1,is).eq.iso.and.isblk1(2,is).eq.isb.and.
     $           isblk1(3,is).eq.isb)then
             imatch=1                                                   8d6s14
             call ilimts(noc(isb),nvirtc(iso),numpro,                      8d6s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d6s14
             nhere=ih+1-il                                              8d6s14
             if(nhere.gt.0)then                                         8d6s14
              i10=i1s                                                   8d6s14
              i1n=noc(isb)                                              8d6s14
              jonex=ionex(is)-1                                         8d6s14
              nrowo=noc(isb)*noc(iso)                                   8d6s14
              do i2=i2s,i2e                                             8d6s14
               i2p=i2+noc(iso)                                          8d6s14
               if(i2.eq.i2e)i1n=i1e                                     8d6s14
               do i1=i10,i1n                                            8d6s14
                do i5=1,noc(iso)                                        8d6s14
                 iad1=itmat(iso)+i2p-1+nbasdwsc(iso)*(i5-1)             8d31s15
                 do i4=1,noc(isb)                                       8d6s14
                  do i3=1,noc(iso)                                      8d6s14
                   iad2=itmp+i4-1+nbasdwsc(isb)*(i3-1+noc(iso)          8d31s15
     $                  *(i1-1+noc(isb)*(i5-1)))                        8d6s14
                   iad3=jonex+i3+noc(iso)*(i4-1)                        8d6s14
                   bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                  8d6s14
                  end do                                                8d6s14
                 end do                                                 8d6s14
                end do                                                  8d6s14
                jonex=jonex+nrowo                                       8d6s14
               end do                                                   8d6s14
               i10=1                                                    8d6s14
              end do                                                    8d6s14
             end if                                                     8d6s14
            end if
            if(isblk1(1,is).eq.isb.and.isblk1(2,is).eq.iso.and.
     $           isblk1(3,is).eq.iso)then
             imatchb=1                                                   8d6s14
             call ilimts(noc(iso),nvirtc(isb),numpro,                      8d6s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d6s14
             nhere=ih+1-il                                              8d6s14
             if(nhere.gt.0)then                                         8d6s14
              i10=i1s                                                   8d6s14
              i1n=noc(iso)                                              8d6s14
              jonex=ionex(is)-1                                         8d6s14
              nrowo=noc(isb)*noc(iso)                                   8d6s14
              do i2=i2s,i2e                                             8d6s14
               i2p=i2+noc(isb)                                          8d6s14
               if(i2.eq.i2e)i1n=i1e                                     8d6s14
               do i1=i10,i1n                                            8d6s14
                do i5=1,noc(iso)                                        8d6s14
                 iad1=itmat(iso)+i1-1+nbasdwsc(iso)*(i5-1)              8d31s15
                 do i4=1,noc(iso)                                       8d6s14
                  do i3=1,noc(isb)                                      8d6s14
                   iad2=itmp+i2p-1+nbasdwsc(isb)*(i4-1+noc(iso)         8d31s15
     $                  *(i3-1+noc(isb)*(i5-1)))                        8d6s14
                   iad3=jonex+i3+noc(isb)*(i4-1)                        8d6s14
                   bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                  8d6s14
                  end do                                                8d6s14
                 end do                                                 8d6s14
                end do                                                  8d6s14
                jonex=jonex+nrowo                                       8d6s14
               end do                                                   8d6s14
               i10=1                                                    8d6s14
              end do                                                    8d6s14
             end if                                                     8d6s14
            else if(isblk1(1,is).eq.iso.and.isblk1(2,is).eq.isb.and.
     $           isblk1(3,is).eq.iso)then
             imatchb=1                                                   8d6s14
             call ilimts(noc(iso),nvirtc(isb),numpro,                      8d6s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d6s14
             nhere=ih+1-il                                              8d6s14
             if(nhere.gt.0)then                                         8d6s14
              i10=i1s                                                   8d6s14
              i1n=noc(iso)                                              8d6s14
              jonex=ionex(is)-1                                         8d6s14
              nrowo=noc(isb)*noc(iso)                                   8d6s14
              do i2=i2s,i2e                                             8d6s14
               i2p=i2+noc(isb)                                          8d6s14
               if(i2.eq.i2e)i1n=i1e                                     8d6s14
               do i1=i10,i1n                                            8d6s14
                do i5=1,noc(iso)                                        8d6s14
                 iad1=itmat(iso)+i1-1+nbasdwsc(iso)*(i5-1)              8d31s15
                 do i4=1,noc(isb)                                       8d6s14
                  do i3=1,noc(iso)                                      8d6s14
                   iad2=itmp+i2p-1+nbasdwsc(isb)*(i3-1+noc(iso)         8d31s15
     $                  *(i4-1+noc(isb)*(i5-1)))                        8d6s14
                   iad3=jonex+i3+noc(iso)*(i4-1)                        8d6s14
                   bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                  8d6s14
                  end do                                                8d6s14
                 end do                                                 8d6s14
                end do                                                  8d6s14
                jonex=jonex+nrowo                                       8d6s14
               end do                                                   8d6s14
               i10=1                                                    8d6s14
              end do                                                    8d6s14
             end if                                                     8d6s14
            end if
           end do                                                       8d6s14
           if(imatch.eq.0)then
            write(6,*)('no match for onex in werup: '),isb,iso
            call dws_sync
            call dws_finalize
            stop
           end if
           if(imatchb.eq.0)then
            write(6,*)('no match for onexb in werup: '),isb,iso
            call dws_sync
            call dws_finalize
            stop
           end if
           imatch=0                                                     8d7s14
           do is=1,nsdlkk                                               8d7s14
            if(isblkk(1,is).eq.iso.and.isblkk(2,is).eq.isb.and.         8d7s14
     $           isblkk(3,is).eq.iso)then                               8d7s14
             imatch=1                                                   8d7s14
             call ilimts(nvirtc(iso),nvirtc(isb),numpro,                      8d6s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d6s14
             nhere=ih+1-il                                              8d6s14
             if(nhere.gt.0)then                                         8d6s14
              i10=i1s                                                   8d7s14
              i1n=nvirtc(iso)                                            8d7s14
              kk=kmats(is)-1                                            8d7s14
              nrowo=noc(isb)*noc(iso)                                   8d7s14
              do i2=i2s,i2e                                             8d7s14
               if(i2.eq.i2e)i1n=i1e                                     8d7s14
               i2p=i2+noc(isb)                                          8d7s14
               do i1=i10,i1n                                            8d7s14
                i1p=i1+noc(iso)                                         8d7s14
                do i5=1,noc(iso)                                        8d7s14
                 iad1=itmat(iso)+i1p-1+nbasdwsc(iso)*(i5-1)             8d31s15
                 do i4=1,noc(isb)                                       8d7s14
                  do i3=1,noc(iso)                                      8d7s14
                   iad2=itmp+i2p-1+nbasdwsc(isb)*(i3-1+noc(iso)         8d31s15
     $                  *(i4-1+noc(isb)*(i5-1)))                        8d7s14
                   iad3=kk+i3+noc(iso)*(i4-1)                           8d7s14
                   bc(iad2)=bc(iad2)+bc(iad1)*bc(iad3)                  8d7s14
                  end do                                                8d7s14
                 end do                                                 8d7s14
                end do                                                  8d7s14
                kk=kk+nrowo                                             8d7s14
               end do                                                   8d7s14
               i10=1                                                    8d7s14
              end do                                                    8d7s14
             end if                                                     8d7s14
            else if(isblkk(1,is).eq.isb.and.isblkk(2,is).eq.iso.and.         8d7s14
     $           isblkk(3,is).eq.isb)then                               8d7s14
             imatch=1                                                   8d7s14
             call ilimts(nvirtc(isb),nvirtc(iso),numpro,                      8d6s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d6s14
             nhere=ih+1-il                                              8d6s14
             if(nhere.gt.0)then                                         8d6s14
              i10=i1s                                                   8d7s14
              i1n=nvirtc(isb)                                            8d7s14
              kk=kmats(is)-1                                            8d7s14
              nrowo=noc(isb)*noc(iso)                                   8d7s14
              do i2=i2s,i2e                                             8d7s14
               if(i2.eq.i2e)i1n=i1e                                     8d7s14
               i2p=i2+noc(iso)                                          8d7s14
               do i1=i10,i1n                                            8d7s14
                i1p=i1+noc(isb)                                         8d7s14
                do i5=1,noc(iso)                                        8d7s14
                 iad1=itmat(iso)+i2p-1+nbasdwsc(iso)*(i5-1)             8d31s15
                 do i4=1,noc(iso)                                       8d7s14
                  do i3=1,noc(isb)                                      8d7s14
                   iad2=itmp+i1p-1+nbasdwsc(isb)*(i4-1+noc(iso)         8d31s15
     $                  *(i3-1+noc(isb)*(i5-1)))                        8d7s14
                   iad3=kk+i3+noc(isb)*(i4-1)                           8d7s14
                   bc(iad2)=bc(iad2)+bc(iad1)*bc(iad3)                  8d7s14
                  end do                                                8d7s14
                 end do                                                 8d7s14
                end do                                                  8d7s14
                kk=kk+nrowo                                             8d7s14
               end do                                                   8d7s14
               i10=1                                                    8d7s14
              end do                                                    8d7s14
             end if                                                     8d7s14
            end if                                                      8d7s14
           end do                                                       8d7s14
           if(imatch.eq.0)then
            write(6,*)('no match for kmat in werup: '),isb,iso
            call dws_sync
            call dws_finalize
            stop
           end if
           do i2=1,noc(isb)                                             8d7s14
            do i1=1,noc(isb)                                            8d7s14
             iad1=ibmat1+i1-1+noc(isb)*(i2-1)                           8d7s14
             sum=0d0                                                    8d7s14
             do l=1,noc(iso)                                            8d7s14
              iad2=itmp+i1-1+nbasdwsc(isb)*(l-1+noc(iso)                8d31s15
     $             *(i2-1+noc(isb)*(l-1)))                              8d7s14
              sum=sum+bc(iad2)                                          8d7s14
             end do                                                     8d7s14
             bc(iad1)=bc(iad1)-2d0*sum                                  8d7s14
            end do                                                      8d7s14
            do i1=1,nvirtc(isb)                                          8d7s14
             i1p=i1+noc(isb)
             iad1=ibmat2+i1-1+nvirtc(isb)*(i2-1)
             sum=0d0                                                    8d7s14
             do l=1,noc(iso)                                            8d7s14
              iad2=itmp+i1p-1+nbasdwsc(isb)*(l-1+noc(iso)               8d31s15
     $             *(i2-1+noc(isb)*(l-1)))                              8d7s14
              sum=sum+bc(iad2)                                          8d7s14
             end do                                                     8d7s14
             bc(iad1)=bc(iad1)-2d0*sum                                  8d7s14
            end do                                                      8d7s14
           end do                                                       8d7s14
c
           do i=0,nrow*noc(iso)*noc(isb)-1                              8d6s14
            bc(itmp+i)=0d0                                              8d6s14
           end do                                                       8d6s14
           imatch=0
           do is=1,nsdlk                                                8d6s14
            if(isblk(1,is).eq.iso.and.isblk(2,is).eq.isb.and.           8d6s14
     $           isblk(4,is).eq.isb)then                                8d6s14
             imatch=1                                                   8d6s14
             call ilimts(noc(iso),noc(isb),numpro,                      8d6s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d6s14
             nhere=ih+1-il                                              8d6s14
             if(nhere.gt.0)then                                         8d6s14
              i10=i1s                                                   8d6s14
              i1n=noc(iso)                                              8d6s14
              joooo=ioooo(is)-1                                         8d6s14
              nrowo=noc(isb)*noc(iso)                                   8d6s14
              do i2=i2s,i2e                                             8d6s14
               if(i2.eq.i2e)i1n=i1e                                     8d6s14
               do i1=i10,i1n                                            8d6s14
                do i5=1,noc(isb)                                        8d6s14
                 iad1=itmat(isb)+i2-1+nbasdwsc(isb)*(i5-1)              8d31s15
                 do i4=1,noc(isb)                                       8d6s14
                  do i3=1,noc(iso)                                      8d6s14
                   iad2=itmp+i4-1+nbasdwsc(isb)*(i3-1+noc(iso)*         8d31s15
     $                  (i5-1+noc(isb)*(i1-1)))                         8d6s14
                   iad3=joooo+i3+noc(iso)*(i4-1)                        8d6s14
                   bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                  8d6s14
                  end do                                                8d6s14
                 end do                                                 8d6s14
                end do                                                  8d6s14
                joooo=joooo+nrowo                                       8d6s14
               end do                                                   8d6s14
               i10=1                                                    8d6s14
              end do                                                    8d6s14
             end if                                                     8d6s14
            else if(isblk(1,is).eq.isb.and.isblk(2,is).eq.iso.and.           8d6s14
     $           isblk(4,is).eq.iso)then                                8d6s14
             imatch=1                                                   8d6s14
             call ilimts(noc(isb),noc(iso),numpro,                      8d6s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d6s14
             nhere=ih+1-il                                              8d6s14
             if(nhere.gt.0)then                                         8d6s14
              i10=i1s                                                   8d6s14
              i1n=noc(isb)                                              8d6s14
              joooo=ioooo(is)-1                                         8d6s14
              nrowo=noc(isb)*noc(iso)                                   8d6s14
              do i2=i2s,i2e                                             8d6s14
               if(i2.eq.i2e)i1n=i1e                                     8d6s14
               do i1=i10,i1n                                            8d6s14
                do i5=1,noc(isb)                                        8d6s14
                 iad1=itmat(isb)+i1-1+nbasdwsc(isb)*(i5-1)              8d31s15
                 do i4=1,noc(iso)                                       8d6s14
                  do i3=1,noc(isb)                                      8d6s14
                   iad2=itmp+i3-1+nbasdwsc(isb)*(i4-1+noc(iso)*         8d31s15
     $                  (i5-1+noc(isb)*(i2-1)))                         8d6s14
                   iad3=joooo+i3+noc(isb)*(i4-1)                        8d6s14
                   bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                  8d6s14
                  end do                                                8d6s14
                 end do                                                 8d6s14
                end do                                                  8d6s14
                joooo=joooo+nrowo                                       8d6s14
               end do                                                   8d6s14
               i10=1                                                    8d6s14
              end do                                                    8d6s14
             end if                                                     8d6s14
            end if                                                      8d6s14
           end do                                                       8d6s14
           if(imatch.eq.0)then
            write(6,*)('no match for oooo in werup: '),isb,iso
            call dws_sync
            call dws_finalize
            stop
           end if
           imatch=0                                                     8d6s14
           imatchb=0
           do is=1,nsdlk1                                               8d6s14
            if(isblk1(1,is).eq.isb.and.isblk1(2,is).eq.iso.and.
     $           isblk1(4,is).eq.isb)then
             imatch=1                                                   8d6s14
             call ilimts(noc(iso),nvirtc(isb),numpro,                      8d6s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d6s14
             nhere=ih+1-il                                              8d6s14
             if(nhere.gt.0)then                                         8d6s14
              i10=i1s                                                   8d6s14
              i1n=noc(iso)                                              8d6s14
              jonex=ionex(is)-1                                         8d6s14
              nrowo=noc(isb)*noc(iso)                                   8d6s14
              do i2=i2s,i2e                                             8d6s14
               i2p=i2+noc(isb)                                          8d6s14
               if(i2.eq.i2e)i1n=i1e                                     8d6s14
               do i1=i10,i1n                                            8d6s14
                do i5=1,noc(isb)                                        8d6s14
                 iad1=itmat(isb)+i2p-1+nbasdwsc(isb)*(i5-1)             8d31s15
                 do i4=1,noc(iso)                                       8d6s14
                  do i3=1,noc(isb)                                      8d6s14
                   iad2=itmp+i3-1+nbasdwsc(isb)*(i4-1+noc(iso)          8d31s15
     $                  *(i5-1+noc(isb)*(i1-1)))                        8d6s14
                   iad3=jonex+i3+noc(isb)*(i4-1)                        8d6s14
                   bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                  8d6s14
                  end do                                                8d6s14
                 end do                                                 8d6s14
                end do                                                  8d6s14
                jonex=jonex+nrowo                                       8d6s14
               end do                                                   8d6s14
               i10=1                                                    8d6s14
              end do                                                    8d6s14
             end if                                                     8d6s14
            else if(isblk1(1,is).eq.iso.and.isblk1(2,is).eq.isb.and.
     $           isblk1(4,is).eq.isb)then
             imatch=1                                                   8d6s14
             call ilimts(noc(iso),nvirtc(isb),numpro,                      8d6s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d6s14
             nhere=ih+1-il                                              8d6s14
             if(nhere.gt.0)then                                         8d6s14
              i10=i1s                                                   8d6s14
              i1n=noc(iso)                                              8d6s14
              jonex=ionex(is)-1                                         8d6s14
              nrowo=noc(isb)*noc(iso)                                   8d6s14
              do i2=i2s,i2e                                             8d6s14
               i2p=i2+noc(isb)                                          8d6s14
               if(i2.eq.i2e)i1n=i1e                                     8d6s14
               do i1=i10,i1n                                            8d6s14
                do i5=1,noc(isb)                                        8d6s14
                 iad1=itmat(isb)+i2p-1+nbasdwsc(isb)*(i5-1)             8d31s15
                 do i4=1,noc(isb)                                       8d6s14
                  do i3=1,noc(iso)                                      8d6s14
                   iad2=itmp+i4-1+nbasdwsc(isb)*(i3-1+noc(iso)          8d31s15
     $                  *(i5-1+noc(isb)*(i1-1)))                        8d6s14
                   iad3=jonex+i3+noc(iso)*(i4-1)                        8d6s14
                   bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                  8d6s14
                  end do                                                8d6s14
                 end do                                                 8d6s14
                end do                                                  8d6s14
                jonex=jonex+nrowo                                       8d6s14
               end do                                                   8d6s14
               i10=1                                                    8d6s14
              end do                                                    8d6s14
             end if                                                     8d6s14
            end if
            if(isblk1(1,is).eq.isb.and.isblk1(2,is).eq.iso.and.
     $           isblk1(4,is).eq.isb)then
             imatchb=1                                                   8d6s14
             call ilimts(noc(iso),nvirtc(isb),numpro,                      8d6s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d6s14
             nhere=ih+1-il                                              8d6s14
             if(nhere.gt.0)then                                         8d6s14
              i10=i1s                                                   8d6s14
              i1n=noc(iso)                                              8d6s14
              jonex=ionex(is)-1                                         8d6s14
              nrowo=noc(isb)*noc(iso)                                   8d6s14
              do i2=i2s,i2e                                             8d6s14
               i2p=i2+noc(isb)                                          8d6s14
               if(i2.eq.i2e)i1n=i1e                                     8d6s14
               do i1=i10,i1n                                            8d6s14
                do i5=1,noc(isb)                                        8d6s14
                 do i4=1,noc(iso)                                       8d6s14
                  do i3=1,noc(isb)                                      8d6s14
                   iad1=itmat(isb)+i3-1+nbasdwsc(isb)*(i5-1)            8d31s15
                   iad2=itmp+i2p-1+nbasdwsc(isb)*(i4-1+noc(iso)         8d31s15
     $                  *(i5-1+noc(isb)*(i1-1)))                        8d6s14
                   iad3=jonex+i3+noc(isb)*(i4-1)                        8d6s14
                   bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                  8d6s14
                  end do                                                8d6s14
                 end do                                                 8d6s14
                end do                                                  8d6s14
                jonex=jonex+nrowo                                       8d6s14
               end do                                                   8d6s14
               i10=1                                                    8d6s14
              end do                                                    8d6s14
             end if                                                     8d6s14
            else if(isblk1(1,is).eq.iso.and.isblk1(2,is).eq.isb.and.
     $           isblk1(3,is).eq.iso)then
             imatchb=1                                                   8d6s14
             call ilimts(noc(iso),nvirtc(isb),numpro,                      8d6s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d6s14
             nhere=ih+1-il                                              8d6s14
             if(nhere.gt.0)then                                         8d6s14
              i10=i1s                                                   8d6s14
              i1n=noc(iso)                                              8d6s14
              jonex=ionex(is)-1                                         8d6s14
              nrowo=noc(isb)*noc(iso)                                   8d6s14
              do i2=i2s,i2e                                             8d6s14
               i2p=i2+noc(isb)                                          8d6s14
               if(i2.eq.i2e)i1n=i1e                                     8d6s14
               do i1=i10,i1n                                            8d6s14
                do i5=1,noc(isb)                                        8d6s14
                 do i4=1,noc(isb)                                       8d6s14
                  iad1=itmat(isb)+i4-1+nbasdwsc(isb)*(i5-1)             8d31s15
                  do i3=1,noc(iso)                                      8d6s14
                   iad2=itmp+i2p-1+nbasdwsc(isb)*(i3-1+noc(iso)         8d31s15
     $                  *(i5-1+noc(isb)*(i1-1)))                        8d6s14
                   iad3=jonex+i3+noc(iso)*(i4-1)                        8d6s14
                   bc(iad2)=bc(iad2)+bc(iad3)*bc(iad1)                  8d6s14
                  end do                                                8d6s14
                 end do                                                 8d6s14
                end do                                                  8d6s14
                jonex=jonex+nrowo                                       8d6s14
               end do                                                   8d6s14
               i10=1                                                    8d6s14
              end do                                                    8d6s14
             end if                                                     8d6s14
            end if
           end do                                                       8d6s14
           if(imatch.eq.0)then
            write(6,*)('no match for onex in werup: '),isb,iso
            call dws_sync
            call dws_finalize
            stop
           end if
           if(imatchb.eq.0)then
            write(6,*)('no match for onexb in werup: '),isb,iso
            call dws_sync
            call dws_finalize
            stop
           end if
           imatch=0                                                     8d7s14
           do is=1,nsdlkk                                               8d7s14
            if(isblkk(1,is).eq.iso.and.isblkk(2,is).eq.iso.and.         8d7s14
     $           isblkk(3,is).eq.isb)then                               8d7s14
             imatch=1                                                   8d7s14
             call ilimts(nvirtc(isb),nvirtc(isb),numpro,                      8d6s14
     $          nowpro,il,ih,i1s,i1e,i2s,i2e)                           8d6s14
             nhere=ih+1-il                                              8d6s14
             if(nhere.gt.0)then                                         8d6s14
              i10=i1s                                                   8d7s14
              i1n=nvirtc(isb)                                            8d7s14
              kk=kmats(is)-1                                            8d7s14
              nrowo=noc(iso)*noc(iso)                                   8d7s14
              do i2=i2s,i2e                                             8d7s14
               if(i2.eq.i2e)i1n=i1e                                     8d7s14
               i2p=i2+noc(isb)                                          8d7s14
               do i1=i10,i1n                                            8d7s14
                i1p=i1+noc(isb)                                         8d7s14
                do i5=1,noc(isb)                                        8d7s14
                 iad1=itmat(isb)+i1p-1+nbasdwsc(isb)*(i5-1)             8d31s15
                 do i4=1,noc(iso)                                       8d7s14
                  do i3=1,noc(iso)                                      8d7s14
                   iad2=itmp+i2p-1+nbasdwsc(isb)*(i3-1+noc(iso)         8d31s15
     $                  *(i5-1+noc(isb)*(i4-1)))                        8d7s14
                   iad3=kk+i3+noc(iso)*(i4-1)                           8d7s14
                   bc(iad2)=bc(iad2)+bc(iad1)*bc(iad3)                  8d7s14
                  end do                                                8d7s14
                 end do                                                 8d7s14
                end do                                                  8d7s14
                kk=kk+nrowo                                             8d7s14
               end do                                                   8d7s14
               i10=1                                                    8d7s14
              end do                                                    8d7s14
             end if                                                     8d7s14
            end if                                                      8d7s14
           end do                                                       8d7s14
           if(imatch.eq.0)then
            write(6,*)('no match for kmat in werup: '),isb,iso
            call dws_sync
            call dws_finalize
            stop
           end if
           do i2=1,noc(isb)                                             8d7s14
            do i1=1,noc(isb)                                            8d7s14
             iad1=ibmat1+i1-1+noc(isb)*(i2-1)                           8d7s14
             sum=0d0                                                    8d7s14
             do l=1,noc(iso)                                            8d7s14
              iad2=itmp+i1-1+nbasdwsc(isb)*(l-1+noc(iso)                8d31s15
     $             *(i2-1+noc(isb)*(l-1)))                              8d7s14
              sum=sum+bc(iad2)                                          8d7s14
             end do                                                     8d7s14
             bc(iad1)=bc(iad1)-2d0*sum                                  8d7s14
            end do                                                      8d7s14
            do i1=1,nvirtc(isb)                                          8d7s14
             i1p=i1+noc(isb)
             iad1=ibmat2+i1-1+nvirtc(isb)*(i2-1)
             sum=0d0                                                    8d7s14
             do l=1,noc(iso)                                            8d7s14
              iad2=itmp+i1p-1+nbasdwsc(isb)*(l-1+noc(iso)               8d31s15
     $             *(i2-1+noc(isb)*(l-1)))                              8d7s14
              sum=sum+bc(iad2)                                          8d7s14
             end do                                                     8d7s14
             bc(iad1)=bc(iad1)-2d0*sum                                  8d7s14
            end do                                                      8d7s14
           end do                                                       8d7s14
          end if                                                        8d6s14
         end if
        end do
        iamat=iamat+noc(isb)*nbasdwsc(isb)                              8d31s15
        ibmat=ibmat+noc(isb)*nbasdwsc(isb)                              8d31s15
       end if                                                           4d19s13
       ih0u=ih0u+nbasdwsc(isb)*nbasdwsc(isb)                            8d31s15
      end do
      call dws_gsumf(bc(ibmat10),nbmat)                                 7d10s13
        ibmat=ibmat10                                                   8d4s14
        do isb=1,nsymb                                                  8d4s14
         itmp1=ibcoff                                                    8d4s14
         itmp2=itmp1+nbasdwsc(isb)*noc(isb)                             8d31s15
         ibcoff=itmp2+nbasdwsc(isb)*noc(isb)                            8d31s15
         call enough('werup.  6',bc,ibc)
         do i=1,nbasdwsc(isb)                                           8d31s15
          do j=1,nbasdwsc(isb)                                          8d31s15
           iad1=itmat(isb)+j-1+nbasdwsc(isb)*(i-1)                      8d31s15
           iad2=iumat(isb)+i-1+nbasdwsc(isb)*(j-1)                      8d31s15
           bc(iad1)=bc(iad2)                                            8d4s14
          end do                                                        8d4s14
         end do                                                         8d4s14
         ibmatp=ibmat+noc(isb)*noc(isb)                                 8d5s14
         do i=1,noc(isb)                                                8d5s14
          do j=1,noc(isb)                                               8d5s14
           iad1=itmp1+j-1+nbasdwsc(isb)*(i-1)                           8d31s15
           iad2=ibmat+j-1+noc(isb)*(i-1)                                8d5s14
           bc(iad1)=bc(iad2)                                            8d5s14
          end do                                                        8d5s14
          do j=1,nvirtc(isb)                                             8d5s14
           iad1=itmp1+j+noc(isb)-1+nbasdwsc(isb)*(i-1)                  8d31s15
           iad2=ibmatp+j-1+nvirtc(isb)*(i-1)                             8d5s14
           bc(iad1)=bc(iad2)                                            8d5s14
          end do                                                        8d5s14
         end do                                                         8d5s14
         if(min(nbasdwsc(isb),noc(isb)).gt.0)then                       6d12s19
         call dgemm('n','n',nbasdwsc(isb),noc(isb),nbasdwsc(isb),1d0,   8d31s15
     $        bc(itmat(isb)),nbasdwsc(isb),bc(itmp1),nbasdwsc(isb),0d0, 8d31s15
     $        bc(itmp2),nbasdwsc(isb),                                  8d31s15
     d' werup.  4')
         end if                                                         6d12s19
         do i=1,noc(isb)                                                8d5s14
          do j=1,noc(isb)                                               8d5s14
           iad1=itmp2+j-1+nbasdwsc(isb)*(i-1)                           8d31s15
           iad2=ibmat+j-1+noc(isb)*(i-1)                                8d5s14
           bc(iad2)=bc(iad1)                                            8d5s14
          end do                                                        8d5s14
          do j=1,nvirtc(isb)                                             8d5s14
           iad1=itmp2+j+noc(isb)-1+nbasdwsc(isb)*(i-1)                  8d31s15
           iad2=ibmatp+j-1+nvirtc(isb)*(i-1)                             8d5s14
           bc(iad2)=bc(iad1)                                            8d5s14
          end do                                                        8d5s14
         end do                                                         8d5s14
         ibmat=ibmat+noc(isb)*noc(isb)                                  8d4s14
         ibmat=ibmat+noc(isb)*nvirtc(isb)                                  8d4s14
         ibcoff=itmp1                                                   8d4s14
        end do                                                          8d4s14
        ibcoff=ibcoffo
      return
      end
