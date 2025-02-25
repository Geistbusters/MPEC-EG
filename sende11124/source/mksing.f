c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine mksing(trans,trans2,nvdim,nrootu,nroot,ihsdiag,vsold,  1d23s21
     $     gsold,ngot,nff1,ncsf,nvirt,mdon,mdoo,nsymb,multh,isymmrci,   2d18s21
     $     gss,gsso,dotsing,ldavid,maxdsr,nlocals,bc,ibc)               11d10s22
      implicit real*8 (a-h,o-z)
      logical ldavid                                                    12d28s21
      integer*8 ihsdiag(mdoo+1,nsymb,2),i18,i28,i38,i48                 1d22s21
      dimension trans(nvdim,nroot),vsold(*),gsold(*),                   1d23s21
     $     trans2(nvdim,nroot),nff1(mdoo+1,nsymb,*),ncsf(*),nvirt(*),   1d23s21
     $     multh(8,8),gss(*),gsso(*),dotsing(2,*)                       2d18s21
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
c
c     trans2 takes orthonormalized external only vectors and h*vectors
c     and forms best full eigenvectors and best full h*eigenvectors
c
c     trans takes normalized external only vectors and h*vectors
c     and forms orthonormal external only vectors and h*vectors that
c     span external space of best full eigenvectors.
c
      joff=1                                                            1d24s21
      joffg=1                                                           12d28s21
      loff=1                                                            2d18s21
      if(ldavid)then                                                    12d28s21
       if(ngot+nroot.le.maxdsr)then                                     12d28s21
        nrootx=ngot+nroot                                               12d28s21
       else                                                             12d28s21
        nrootx=max(ngot,nroot)                                          12d29s21
       end if                                                           12d28s21
       ivold=ibcoff                                                     12d28s21
       ndelta=(nlocals/nroot)*nrootx                                    12d28s21
       igold=ivold+ndelta                                               12d28s21
       igssold=igold+ndelta                                             1d24s23
       ibcoff=igssold+ndelta                                            1d24s23
       call enough('mksing.  1',bc,ibc)
       jvold=ivold                                                      1d19s23
       jgssold=igssold                                                  1d24s23
      else                                                              12d28s21
       nrootx=nroot                                                     12d28s21
      end if                                                            12d28s21
      do isb=1,nsymb
       isbv=multh(isb,isymmrci)
       do nclo=mdon,mdoo                                                1d21s21
        nclop=nclo+1                                                    1d21s21
        if(min(nvirt(isbv),nff1(nclop,isb,1)).gt.0)then                 1d21s21
         iarg=nclop-mdon                                                 1d21s21
         nxx=nff1(nclop,isb,1)*ncsf(iarg)                               1d21s21
         call ilimts(1,nxx,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)    1d21s21
         nhere=ih+1-il                                                  1d21s21
         if(nhere.gt.0)then                                             1d21s21
          nrow=nvirt(isbv)*nrootu                                       1d21s21
          nrowb=nvirt(isbv)*nroot                                       1d23s21
          ntt=nrowb*nhere                                               1d23s21
          ivec=ibcoff                                                   1d21s21
          ig=ivec+ntt                                                   1d21s21
          ibcoff=ig+ntt                                                 1d21s21
          call enough('mksing.  2',bc,ibc)
          i18=1                                                         1d21s21
          i28=nrow                                                      1d21s21
          i38=il                                                        1d21s21
          i48=ih                                                        1d21s21
          call ddi_get(bc,ibc,ihsdiag(nclop,isb,1),i18,i28,i38,i48,     11d15s22
     $         bc(ivec))                                                11d15s22
          call ddi_get(bc,ibc,ihsdiag(nclop,isb,2),i18,i28,i38,i48,     11d15s22
     $         bc(ig))                                                  11d15s22
          nhere2=nhere*2                                                1d22s21
          nrow=nvirt(isbv)*nhere2
          nrows=nvirt(isbv)*nhere                                       2d18s21
          itmp=ibcoff                                                   1d22s21
          itmp1=itmp+nrow*nrootx                                        12d28s21
          itmp2=itmp1+nrow*nrootx                                       12d28s21
          ibcoff=itmp2+nrow*nroot                                       1d23s21
          istmp1=ibcoff                                                 2d18s21
          istmp2=istmp1+nrows*nrootx                                    1d24s23
          istmp3=istmp2+nrows*nroot                                     2d18s21
          ibcoff=istmp3+nrows*nroot                                     2d18s21
          call enough('mksing.  3',bc,ibc)
          nn=nvirt(isbv)*nhere                                          1d22s21
          if(ngot.gt.0)then                                             1d22s21
           koff=joffg                                                   12d28s21
           do i=0,nhere-1                                               1d22s21
            do iv=0,nvirt(isbv)-1                                       1d22s21
             do ir=0,ngot-1                                             1d22s21
              iad=itmp1+iv+nvirt(isbv)*(i+nhere2*ir)                    1d22s21
              jad=iad+nn                                                1d22s21
              bc(iad)=vsold(koff+ir)                                    1d22s21
              bc(jad)=gsold(koff+ir)                                    1d22s21
              iad=istmp1+iv+nvirt(isbv)*(i+nhere*ir)                    2d18s21
              bc(iad)=gsso(koff+ir)                                     2d18s21
             end do                                                     1d22s21
             koff=koff+ngot                                             1d22s21
            end do                                                      1d22s21
           end do                                                       1d22s21
           np=nrootu+1                                                  1d22s21
           if(ldavid)then                                               12d28s21
            if(nroot+ngot.le.maxdsr)then                                12d28s21
             do i=0,nrow*ngot-1                                          12d28s21
              bc(itmp+i)=bc(itmp1+i)                                     12d28s21
             end do                                                      12d28s21
            end if                                                      12d28s21
           else                                                         12d28s21
            call dgemm('n','n',nrow,nroot,ngot,1d0,bc(itmp1),nrow,       1d22s21
     $          trans(np,1),nvdim,0d0,bc(itmp),nrow,                    1d22s21
     d' mksing.  1')
           end if                                                       12d28s21
           call dgemm('n','n',nrow,nroot,ngot,1d0,bc(itmp1),nrow,       1d22s21
     $          trans2(np,1),nvdim,0d0,bc(itmp2),nrow,                  1d23s21
     d' mksing.  2')
           call dgemm('n','n',nrows,nroot,ngot,1d0,bc(istmp1),nrows,    2d18s21
     $          trans2(np,1),nvdim,0d0,bc(istmp2),nrows,                2d18s21
     d' mksing.  3')
           call dgemm('n','n',nrows,nroot,ngot,1d0,bc(istmp1),nrows,    2d18s21
     $          trans(np,1),nvdim,0d0,bc(istmp3),nrows,                 2d18s21
     d' mksing.  4')
           fact=1d0                                                     1d22s21
          else                                                          1d22s21
           fact=0d0                                                     1d22s21
          end if                                                        1d22s21
          jvec=ivec                                                     1d22s21
          jg=ig                                                         1d22s21
          koff=joff                                                     2d18s21
          do i=0,nhere-1                                                1d22s21
           do iv=0,nvirt(isbv)-1                                        1d22s21
            do ir=0,nrootu-1                                            1d22s21
             iad=itmp1+iv+nvirt(isbv)*(i+nhere2*ir)                     1d22s21
             jad=iad+nn                                                 1d22s21
             bc(iad)=bc(jvec+ir)                                        1d22s21
             bc(jad)=bc(jg+ir)                                          1d22s21
             iad=istmp1+iv+nvirt(isbv)*(i+nhere*ir)                     2d18s21
             bc(iad)=gss(loff+ir)
            end do                                                      1d22s21
            jvec=jvec+nrootu                                            1d22s21
            loff=loff+nrootu                                            2d18s21
            jg=jg+nrootu                                                1d22s21
           end do                                                       1d22s21
          end do                                                        1d22s21
          if(ldavid)then                                                12d28s21
           if(nroot+ngot.le.maxdsr)then                                 12d28s21
            koff=joffg                                                  1d24s23
            do i=0,nhere-1                                                1d22s21
             do iv=0,nvirt(isbv)-1                                        1d22s21
              do ir=0,nroot-1                                             1d22s21
               irp=ir+ndelta                                            12d28s21
               iad=itmp1+iv+nvirt(isbv)*(i+nhere2*ir)                      1d22s21
               jad=iad+nn                                                 1d22s21
               bc(jvold+ir)=bc(iad)                                     12d28s21
               bc(jvold+irp)=bc(jad)                                    12d28s21
               iad=istmp1+iv+nvirt(isbv)*(i+nhere*ir)                   1d24s23
               bc(jgssold+ir)=bc(iad)                                   1d24s23
              end do                                                      1d22s21
              do ir=0,ngot-1                                             1d22s21
               irp=ir+nroot                                             12d28s21
               irpp=irp+ndelta                                          12d28s21
               iad=itmp+iv+nvirt(isbv)*(i+nhere2*ir)                    12d28s21
               jad=iad+nn                                                 1d22s21
               bc(jvold+irp)=bc(iad)                                    12d28s21
               bc(jvold+irpp)=bc(jad)                                   12d28s21
               bc(jgssold+irp)=gsso(koff+ir)                            1d24s23
              end do                                                      1d22s21
              jvold=jvold+nrootx                                        12d28s21
              jgssold=jgssold+nrootx                                    1d24s23
              koff=koff+ngot                                            1d24s23
             end do                                                       1d22s21
            end do                                                        1d22s21
           end if                                                       12d28s21
          else                                                          12d28s21
           call dgemm('n','n',nrow,nroot,nrootu,1d0,bc(itmp1),nrow,      1d22s21
     $         trans,nvdim,fact,bc(itmp),nrow,                          1d22s21
     d' mksing.  5')
           koff=joff                                                     1d22s21
           do i=0,nhere-1                                                1d22s21
            do iv=0,nvirt(isbv)-1                                        1d22s21
             do ir=0,nroot-1                                             1d22s21
              irp=ir+1                                                   2d18s21
              iad=itmp+iv+nvirt(isbv)*(i+nhere2*ir)                      1d22s21
              jad=iad+nn                                                 1d22s21
              vsold(koff+ir)=bc(iad)                                     1d22s21
              gsold(koff+ir)=bc(jad)                                     1d22s21
             end do                                                      1d22s21
             koff=koff+nroot                                             1d22s21
            end do                                                       1d22s21
           end do                                                        1d22s21
          end if                                                        12d28s21
          call dgemm('n','n',nrow,nroot,nrootu,1d0,bc(itmp1),nrow,      1d22s21
     $         trans2,nvdim,fact,bc(itmp2),nrow,                          1d22s21
     d' mksing.  6')
          call dgemm('n','n',nrows,nroot,nrootu,1d0,bc(istmp1),nrows,   2d18s21
     $         trans2,nvdim,fact,bc(istmp2),nrows,                      2d18s21
     d' mksing.  7')
          call dgemm('n','n',nrows,nroot,nrootu,1d0,bc(istmp1),nrows,   2d18s21
     $         trans,nvdim,fact,bc(istmp3),nrows,                       2d18s21
     d' mksing.  8')
          koff=joff                                                     1d22s21
          jvec=ivec                                                     1d23s21
          jg=ig                                                         1d23s21
          if(.not.ldavid)then                                           1d24s23
           do i=0,nhere-1                                                1d22s21
            do iv=0,nvirt(isbv)-1                                        1d22s21
             do ir=0,nroot-1                                             1d22s21
              irp=ir+1                                                   2d18s21
              iad2=itmp2+iv+nvirt(isbv)*(i+nhere2*ir)                      1d22s21
              jad2=iad2+nn                                                 1d22s21
              bc(jvec+ir)=bc(iad2)                                       1d23s21
              dotsing(1,irp)=dotsing(1,irp)+bc(iad2)*bc(iad2)            2d18s21
              bc(jg+ir)=bc(jad2)                                         1d23s21
              iad=istmp2+iv+nvirt(isbv)*(i+nhere*ir)                     2d18s21
              dotsing(2,irp)=dotsing(2,irp)                              2d18s21
     $             +bc(iad2)*bc(iad)                                     2d19s21
              gss(koff+ir)=bc(iad)                                       2d18s21
              iad=istmp3+iv+nvirt(isbv)*(i+nhere*ir)                     2d18s21
              gsso(koff+ir)=bc(iad)                                      2d18s21
             end do                                                      1d22s21
             koff=koff+nroot                                             1d22s21
             jvec=jvec+nroot                                             1d23s21
             jg=jg+nroot                                                 1d23s21
            end do                                                       1d22s21
           end do                                                        1d22s21
          else                                                          1d24s23
           do i=0,nhere-1                                                1d22s21
            do iv=0,nvirt(isbv)-1                                        1d22s21
             do ir=0,nroot-1                                             1d22s21
              irp=ir+1                                                   2d18s21
              iad2=itmp2+iv+nvirt(isbv)*(i+nhere2*ir)                      1d22s21
              jad2=iad2+nn                                                 1d22s21
              bc(jvec+ir)=bc(iad2)                                       1d23s21
              dotsing(1,irp)=dotsing(1,irp)+bc(iad2)*bc(iad2)            2d18s21
              bc(jg+ir)=bc(jad2)                                         1d23s21
              iad=istmp2+iv+nvirt(isbv)*(i+nhere*ir)                     2d18s21
              dotsing(2,irp)=dotsing(2,irp)                              2d18s21
     $             +bc(iad2)*bc(iad)                                     2d19s21
             end do                                                      1d22s21
             koff=koff+nroot                                             1d22s21
             jvec=jvec+nroot                                             1d23s21
             jg=jg+nroot                                                 1d23s21
            end do                                                       1d22s21
           end do                                                        1d22s21
          end if                                                        1d24s23
          joff=koff                                                     1d22s21
          joffg=joffg+ngot*nvirt(isbv)*nhere                            12d28s21
          i18=1                                                         1d21s21
          i28=nrowb                                                      1d21s21
          i38=il                                                        1d21s21
          i48=ih                                                        1d21s21
          call ddi_put(bc,ibc,ihsdiag(nclop,isb,1),i18,i28,i38,i48,     11d15s22
     $         bc(ivec))                                                11d15s22
          call ddi_put(bc,ibc,ihsdiag(nclop,isb,2),i18,i28,i38,i48,     11d15s22
     $         bc(ig))                                                  11d15s22
          ibcoff=ivec                                                   1d22s21
         end if                                                         1d22s21
        end if                                                          1d22s21
       end do                                                           1d22s21
      end do                                                            1d22s21
      if(ldavid)then                                                    12d28s21
       if(nroot+ngot.le.maxdsr)then                                     12d28s21
        jvold=ivold-1                                                   12d28s21
        jgold=igold-1                                                   12d28s21
        jgssold=igssold-1                                               1d24s23
        do i=1,ndelta                                                   12d28s21
         vsold(i)=bc(jvold+i)                                           12d28s21
         gsold(i)=bc(jgold+i)                                           12d28s21
         gsso(i)=bc(jgssold+i)                                          1d24s23
        end do                                                          12d28s21
        ibcoff=ivold                                                    12d28s21
       end if                                                           12d28s21
      end if                                                            12d28s21
      call dws_gsumf(dotsing,nroot*2)                                   2d18s21
      return                                                            1d22s21
      end                                                               1d22s21
