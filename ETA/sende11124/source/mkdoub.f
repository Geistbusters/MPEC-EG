c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine mkdoub(trans,trans2,nvdim,nrootu,nroot,vinout,hd,      1d23s21
     $     vdold,gdold,ngot,nfdat,nvirt,nsymb,multh,isymmrci,ndoub,     12d28s21
     $     ldavid,maxdsr,bc,ibc)                                        11d10s22
      implicit real*8 (a-h,o-z)                                         1d22s21
      logical ldavid                                                    12d28s21
      dimension trans(nvdim,nroot),vinout(*),hd(*),vdold(*),            1d22s21
     $     gdold(*),nfdat(5,4,*),nvirt(*),multh(8,8),trans2(nvdim,*)    1d23s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      include "common.store"                                            1d22s21
c
c     trans2 takes orthonormalized external only vectors and h*vectors
c     and forms best full eigenvectors and best full h*eigenvectors
c
c     trans takes normalized external only vectors and h*vectors
c     and forms orthonormal external only vectors and h*vectors that
c     span external space of best full eigenvectors.
c
      ivtmp=ibcoff                                                      1d23s21
      igtmp=ivtmp+ndoub*nroot                                           1d23s21
      ibcoff=igtmp+ndoub*nroot                                          1d23s21
      call enough('mkdoub.  1',bc,ibc)
      do i=ivtmp,ibcoff-1                                               1d23s21
       bc(i)=0d0                                                        1d23s21
      end do                                                            1d23s21
      ioff=0                                                            1d22s21
      joff=0                                                            1d24s21
      joffg=0                                                           12d28s21
      jvtmp=ivtmp-1                                                     1d24s21
      jgtmp=igtmp-1                                                     1d24s21
      if(ldavid)then                                                    12d28s21
       if(ngot+nroot.le.maxdsr)then                                     12d28s21
        nrootx=ngot+nroot                                               12d28s21
       else                                                             12d28s21
        nrootx=max(ngot,nroot)                                                    12d28s21
       end if                                                           12d28s21
       ivold=ibcoff                                                     12d28s21
       igold=ivold+ndoub*nrootx                                         12d28s21
       ibcoff=igold+ndoub*nrootx                                        12d28s21
       ndelta=ndoub*nrootx                                              12d28s21
       call enough('mkdoub.  2',bc,ibc)
       jvold=ivold-1                                                    10s28s24
      else                                                              12d28s21
       nrootx=nroot                                                     12d28s21
      end if                                                            12d28s21
      do isb=1,nsymb                                                    1d22s21
       isbv12=multh(isb,isymmrci)                                       1d22s21
       do isbv1=1,nsymb                                                 1d22s21
        isbv2=multh(isbv1,isbv12)                                       1d22s21
        if(isbv2.ge.isbv1)then                                          1d22s21
         if(isbv1.eq.isbv2)then                                         1d22s21
          if(nvirt(isbv1).gt.0)then                                     1d25s21
           nhere2=nvirt(isbv1)*2                                         1d24s21
           nrow=nhere2*nfdat(3,1,isb)                                    1d22s21
           itmp=ibcoff                                                   1d22s21
           itmp1=itmp+nrow*nrootx                                       12d28s21
           itmp2=itmp1+nrow*nrootx                                      12d28s21
           ibcoff=itmp2+nrow*nroot                                       1d23s21
           call enough('mkdoub.  3',bc,ibc)
           jtmp=itmp-1                                                   1d24s21
           jtmp1=itmp1-1                                                  1d24s21
           jtmp2=itmp2-1                                                 1d24s21
           if(min(ngot,nfdat(3,1,isb)).gt.0)then                        11d17s21
            koff=joffg                                                  12d28s21
            do i=0,nfdat(3,1,isb)-1                                      1d22s21
             do ir=0,ngot-1                                              1d22s21
              iad=jtmp1+nhere2*(i+nfdat(3,1,isb)*ir)                     1d22s21
              jad=iad+nvirt(isbv1)                                        1d24s21
              do iv=1,nvirt(isbv1)                                        1d24s21
               bc(iad+iv)=vdold(koff+iv)                                 1d22s21
               bc(jad+iv)=gdold(koff+iv)                                 1d22s21
              end do                                                     1d22s21
              koff=koff+nvirt(isbv1)                                      1d24s21
             end do                                                      1d22s21
            end do                                                       1d22s21
            np=nrootu+1                                                  1d22s21
            if(.not.ldavid)then                                         12d28s21
             call dgemm('n','n',nrow,nroot,ngot,1d0,bc(itmp1),nrow,       1d22s21
     $           trans(np,1),nvdim,0d0,bc(itmp),nrow,                    1d22s21
     d' mkdoub.  1')
            else
             if(ngot+nroot.le.maxdsr)then                               12d28s21
              do i=0,nrow*ngot-1                                        12d28s21
               bc(itmp+i)=bc(itmp1+i)                                   12d28s21
              end do                                                    12d28s21
             end if                                                     12d28s21
            end if                                                      12d28s21
            call dgemm('n','n',nrow,nroot,ngot,1d0,bc(itmp1),nrow,       1d22s21
     $           trans2(np,1),nvdim,0d0,bc(itmp2),nrow,                    1d22s21
     d' mkdoub.  2')
            fact=1d0                                                     1d22s21
           else                                                          1d22s21
            fact=0d0                                                     1d22s21
           end if                                                        1d22s21
           if(nfdat(3,1,isb).gt.0)then                                  11d17s21
            do i=0,nfdat(3,1,isb)-1                                       1d22s21
             do ir=0,nrootu-1                                             1d22s21
              iad=jtmp1+nhere2*(i+nfdat(3,1,isb)*ir)                      1d22s21
              jad=iad+nvirt(isbv1)                                         1d24s21
              do iv=1,nvirt(isbv1)                                         1d24s21
               bc(iad+iv)=vinout(ioff+iv)                                 1d22s21
               bc(jad+iv)=hd(ioff+iv)                                     1d22s21
              end do                                                      1d22s21
              ioff=ioff+nvirt(isbv1)                                      1d22s21
             end do                                                       1d22s21
            end do                                                        1d22s21
            if(.not.ldavid)then                                         12d28s21
             call dgemm('n','n',nrow,nroot,nrootu,1d0,bc(itmp1),nrow,      1d22s21
     $          trans,nvdim,fact,bc(itmp),nrow,                          1d22s21
     d' mkdoub.  3')
             jjoff=joff                                                 12d28s21
             do i=0,nfdat(3,1,isb)-1                                       1d22s21
              do ir=0,nroot-1                                              1d22s21
               iad=jtmp+nhere2*(i+nfdat(3,1,isb)*ir)                       1d22s21
               jad=iad+nvirt(isbv1)                                         1d24s21
               do iv=1,nvirt(isbv1)                                         1d24s21
                vdold(jjoff+iv)=bc(iad+iv)                                  1d22s21
                gdold(jjoff+iv)=bc(jad+iv)                                  1d22s21
               end do                                                      1d22s21
               jjoff=jjoff+nvirt(isbv1)                                 12d28s21
              end do                                                       1d22s21
             end do                                                        1d22s21
            else                                                        12d28s21
             if(ngot+nroot.le.maxdsr)then                               12d28s21
              do i=0,nfdat(3,1,isb)-1                                       1d22s21
               do ir=0,nroot-1                                              1d22s21
                iad=jtmp1+nhere2*(i+nfdat(3,1,isb)*ir)                  12d28s21
                jad=iad+nvirt(isbv1)                                         1d24s21
                kvold=jvold+nvirt(isbv1)*(ir+nrootx*i)                  12d28s21
                kgold=kvold+ndelta                                      12d28s21
                do iv=1,nvirt(isbv1)                                         1d24s21
                 bc(kvold+iv)=bc(iad+iv)                                12d28s21
                 bc(kgold+iv)=bc(jad+iv)                                12d28s21
                end do                                                      1d22s21
               end do                                                       1d22s21
               do ir=0,ngot-1                                           12d28s21
                iad=jtmp+nhere2*(i+nfdat(3,1,isb)*ir)                   12d28s21
                jad=iad+nvirt(isbv1)                                         1d24s21
                kvold=jvold+nvirt(isbv1)*(nroot+ir+nrootx*i)            12d28s21
                kgold=kvold+ndelta                                      12d28s21
                do iv=1,nvirt(isbv1)                                         1d24s21
                 bc(kvold+iv)=bc(iad+iv)                                12d28s21
                 bc(kgold+iv)=bc(jad+iv)                                12d28s21
                end do                                                      1d22s21
               end do                                                   12d28s21
              end do                                                        1d22s21
             end if                                                     12d28s21
            end if                                                      12d28s21
            call dgemm('n','n',nrow,nroot,nrootu,1d0,bc(itmp1),nrow,      1d23s21
     $          trans2,nvdim,fact,bc(itmp2),nrow,                        1d23s21
     d' mkdoub.  4')
            do i=0,nfdat(3,1,isb)-1                                       1d22s21
             do ir=0,nroot-1                                              1d22s21
              iad2=jtmp2+nhere2*(i+nfdat(3,1,isb)*ir)                     1d23s21
              jad2=iad2+nvirt(isbv1)                                       1d24s21
              do iv=1,nvirt(isbv1)                                         1d24s21
               bc(jvtmp+iv)=bc(iad2+iv)                                   1d23s21
               bc(jgtmp+iv)=bc(jad2+iv)                                   1d23s21
              end do                                                      1d22s21
              joff=joff+nvirt(isbv1)                                      1d24s21
              jvtmp=jvtmp+nvirt(isbv1)                                    1d23s21
              jgtmp=jgtmp+nvirt(isbv1)                                    1d23s21
             end do                                                       1d22s21
             joffg=joffg+nvirt(isbv1)*ngot                              12d28s21
            end do                                                        1d22s21
           end if                                                       11d17s21
           ibcoff=itmp                                                   1d22s21
          end if                                                        1d25s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2
          if(ldavid)jvold=jvold+nvirt(isbv1)*nfdat(3,1,isb)*nrootx      10d30s24
         else                                                           1d22s21
          nvv=nvirt(isbv1)*nvirt(isbv2)
         end if                                                         1d22s21
         do l=1,4
          if(min(nvv,nfdat(3,l,isb)).gt.0)then                          1d22s21
           nhere2=nvv*2                                                 1d24s21
           nrow=nhere2*nfdat(3,l,isb)                                   1d22s21
           itmp=ibcoff                                                  1d22s21
           itmp1=itmp+nrow*nrootx                                       12d28s21
           itmp2=itmp1+nrow*nrootx                                      12d28s21
           ibcoff=itmp2+nrow*nroot                                      1d23s21
           call enough('mkdoub.  4',bc,ibc)
           jtmp=itmp-1                                                  1d24s21
           jtmp1=itmp1-1                                                1d24s21
           jtmp2=itmp2-1                                                1d24s21
           koff=joffg                                                   12d28s21
           if(ngot.gt.0)then                                            1d22s21
            do i=0,nfdat(3,l,isb)-1                                     1d22s21
             do ir=0,ngot-1                                             1d22s21
              iad=jtmp1+nhere2*(i+nfdat(3,l,isb)*ir)                    1d22s21
              jad=iad+nvv                                               1d24s21
              do iv=1,nvv                                               1d24s21
               bc(iad+iv)=vdold(koff+iv)                                1d22s21
               bc(jad+iv)=gdold(koff+iv)                                1d22s21
              end do                                                    1d22s21
              koff=koff+nvv                                             1d24s21
             end do                                                     1d22s21
            end do                                                      1d22s21
            np=nrootu+1                                                 1d22s21
            if(.not.ldavid)then                                         12d28s21
             call dgemm('n','n',nrow,nroot,ngot,1d0,bc(itmp1),nrow,      1d22s21
     $           trans(np,1),nvdim,0d0,bc(itmp),nrow,                   1d22s21
     d' mkdoub.  5')
            else
             if(ngot+nroot.le.maxdsr)then                               12d28s21
              do i=0,nrow*ngot-1                                        12d28s21
               bc(itmp+i)=bc(itmp1+i)                                   12d28s21
              end do                                                    12d28s21
             end if                                                     12d28s21
            end if                                                      12d28s21
            call dgemm('n','n',nrow,nroot,ngot,1d0,bc(itmp1),nrow,      1d23s21
     $           trans2(np,1),nvdim,0d0,bc(itmp2),nrow,                 1d23s21
     d' mkdoub.  6')
            fact=1d0                                                    1d22s21
           else                                                         1d22s21
            fact=0d0                                                    1d22s21
           end if                                                       1d22s21
           do i=0,nfdat(3,l,isb)-1                                      1d22s21
            do ir=0,nrootu-1                                            1d22s21
             iad=jtmp1+nhere2*(i+nfdat(3,l,isb)*ir)                     1d22s21
             jad=iad+nvv                                                1d24s21
             do iv=1,nvv                                                1d24s21
              bc(iad+iv)=vinout(ioff+iv)                                1d22s21
              bc(jad+iv)=hd(ioff+iv)                                    1d22s21
             end do                                                     1d22s21
             ioff=ioff+nvv                                              1d24s21
            end do                                                      1d22s21
           end do                                                       1d22s21
           if(.not.ldavid)then                                          12d28s21
            call dgemm('n','n',nrow,nroot,nrootu,1d0,bc(itmp1),nrow,     1d22s21
     $          trans,nvdim,fact,bc(itmp),nrow,                         1d22s21
     d' mkdoub.  7')
            jjoff=joff                                                  12d28s21
            do i=0,nfdat(3,l,isb)-1                                      1d22s21
             do ir=0,nroot-1                                             1d22s21
              iad=jtmp+nhere2*(i+nfdat(3,l,isb)*ir)                      1d24s21
              jad=iad+nvv                                                1d24s21
              do iv=1,nvv                                                1d24s21
               vdold(jjoff+iv)=bc(iad+iv)                               12d28s21
               gdold(jjoff+iv)=bc(jad+iv)                               12d28s21
              end do                                                     1d22s21
              jjoff=jjoff+nvv                                           12d28s21
             end do                                                      1d22s21
            end do                                                       1d22s21
           else                                                         12d28s21
            if(ngot+nroot.le.maxdsr)then                                12d28s21
             do i=0,nfdat(3,l,isb)-1                                    12d28s21
              do ir=0,nroot-1                                              1d22s21
               iad=jtmp1+nhere2*(i+nfdat(3,l,isb)*ir)                   12d28s21
               jad=iad+nvv                                              12d28s21
               kvold=jvold+nvv*(ir+nrootx*i)                            12d28s21
               kgold=kvold+ndelta                                       12d28s21
               do iv=1,nvv                                              12d28s21
                bc(kvold+iv)=bc(iad+iv)                                 12d28s21
                bc(kgold+iv)=bc(jad+iv)                                 12d28s21
               end do                                                      1d22s21
              end do                                                       1d22s21
              do ir=0,ngot-1                                            12d28s21
               iad=jtmp+nhere2*(i+nfdat(3,l,isb)*ir)                    12d28s21
               jad=iad+nvv                                              12d28s21
               kvold=jvold+nvv*(nroot+ir+nrootx*i)                      12d28s21
               kgold=kvold+ndelta                                       12d28s21
               do iv=1,nvv                                              12d28s21
                bc(kvold+iv)=bc(iad+iv)                                 12d28s21
                bc(kgold+iv)=bc(jad+iv)                                 12d28s21
               end do                                                      1d22s21
              end do                                                    12d28s21
             end do                                                        1d22s21
            end if                                                      12d28s21
           end if                                                       12d28s21
           call dgemm('n','n',nrow,nroot,nrootu,1d0,bc(itmp1),nrow,     1d22s21
     $          trans2,nvdim,fact,bc(itmp2),nrow,                       1d23s21
     d' mkdoub.  8')
           do i=0,nfdat(3,l,isb)-1                                      1d22s21
            do ir=0,nroot-1                                             1d22s21
             iad2=jtmp2+nhere2*(i+nfdat(3,l,isb)*ir)                    1d24s21
             jad2=iad2+nvv                                              1d24s21
             do iv=1,nvv                                                1d24s21
              bc(jvtmp+iv)=bc(iad2+iv)                                  1d24s21
              bc(jgtmp+iv)=bc(jad2+iv)                                  1d24s21
             end do                                                     1d22s21
             joff=joff+nvv                                              1d24s21
             jvtmp=jvtmp+nvv                                            1d24s21
             jgtmp=jgtmp+nvv                                            1d24s21
            end do                                                      1d22s21
           end do                                                       1d22s21
           if(ldavid)jvold=jvold+nvv*nfdat(3,l,isb)*nrootx              10d29s24
           joffg=joffg+nvv*nfdat(3,l,isb)*ngot                          12d28s21
           ibcoff=itmp                                                  1d22s21
          end if
         end do
        end if                                                          1d22s21
       end do                                                           1d22s21
      end do                                                            1d22s21
      jvtmp=ivtmp-1                                                     1d23s21
      jgtmp=igtmp-1                                                     1d23s21
      if(mynowprog.eq.0)then                                            1d25s21
       do i=1,ndoub*nroot                                               1d25s21
        vinout(i)=bc(jvtmp+i)                                            1d23s21
        hd(i)=bc(jgtmp+i)                                               1d25s21
       end do                                                            1d23s21
      else                                                              1d25s21
       do i=1,ndoub*nroot                                               1d25s21
        vinout(i)=bc(jvtmp+i)                                            1d23s21
        hd(i)=0d0                                                       1d25s21
       end do                                                            1d23s21
      end if                                                            1d25s21
      call dws_bcast(vinout,ndoub*nroot)                                2d9s21
      if(ldavid)then                                                    12d28s21
       if(nroot+ngot.le.maxdsr)then
        jvold=ivold-1                                                    12d28s21
        jgold=igold-1                                                    12d28s21
        do i=1,ndelta                                                    12d28s21
         vdold(i)=bc(jvold+i)                                            12d28s21
         gdold(i)=bc(jgold+i)                                            12d28s21
        end do                                                           12d28s21
       end if
      end if                                                            12d28s21
      ibcoff=ivtmp                                                      1d23s21
      return                                                            1d22s21
      end                                                               1d22s21
