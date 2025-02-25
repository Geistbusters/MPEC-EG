c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine mkdoubuc(trans,trans2,nvdim,nrootu,nroot,ihddiag,vdold,  1d23s21
     $     gdold,ngot,nff2,ncsf,ncsf2,nvirt,mdon,mdoo,nsymb,multh,      7d7s21
     $     isymmrci,ldavid,maxdsr,bc,ibc)                               11d10s22
      implicit real*8 (a-h,o-z)
      logical ldavid                                                    12d28s21
      integer*8 ihddiag(mdoo+1,nsymb,2),i18,i28,i38,i48                 1d22s21
      dimension trans(nvdim,nroot),vdold(*),gdold(*),                   1d23s21
     $     trans2(nvdim,nroot),nff2(mdoo+1,nsymb,*),ncsf(*),nvirt(*),   1d23s21
     $     multh(8,8),ncsf2(4,*)                                        7d7s21
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      if(ldavid)then                                                    1d5s22
       write(6,*)('in mkdoubuc, ldavid = '),ldavid
       write(6,*)('I have no code for that yet!!!!')                    12d28s21
       call dws_synca                                                   12d28s21
       call dws_finalize                                                12d28s21
       stop                                                             12d28s21
      end if                                                            12d28s21
c
c     trans2 takes orthonormalized external only vectors and h*vectors
c     and forms best full eigenvectors and best full h*eigenvectors
c
c     trans takes normalized external only vectors and h*vectors
c     and forms orthonormal external only vectors and h*vectors that
c     span external space of best full eigenvectors.
c
      joff=1                                                            1d24s21
      do isb=1,nsymb
       isbv12=multh(isb,isymmrci)                                       7d7s21
       nvisv=0                                                          7d7s21
       nvnotv=0                                                         7d7s21
       do isbv1=1,nsymb                                                 7d7s21
        isbv2=multh(isbv1,isbv12)                                       7d7s21
        if(isbv2.ge.isbv1)then                                          7d7s21
         if(isbv1.eq.isbv2)then                                         7d7s21
          nvisv=nvisv+nvirt(isbv1)                                      7d7s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         7d7s21
         else                                                           7d7s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 7d7s21
         end if                                                         7d7s21
         nvnotv=nvnotv+nvv                                              7d7s21
        end if                                                          7d7s21
       end do                                                           7d7s21
       do nclo=mdon,mdoo                                                1d21s21
        nclop=nclo+1                                                    1d21s21
        if(min(nvisv+nvnotv,nff2(nclop,isb,1)).gt.0)then                7d7s21
         iarg=nclop-mdon                                                 1d21s21
         nxx=nff2(nclop,isb,1)                                          7d7s21
         call ilimts(1,nxx,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)    1d21s21
         nhere=ih+1-il                                                  1d21s21
         if(nhere.gt.0)then                                             1d21s21
          mvv=(nvisv*ncsf2(1,iarg)+nvnotv*ncsf(iarg))                   7d7s21
          nrow=mvv*nrootu                                               7d7s21
          nrowb=mvv*nroot                                               7d7s21
          ntt=nrowb*nhere                                               1d23s21
          ivec=ibcoff                                                   1d21s21
          ig=ivec+ntt                                                   1d21s21
          ibcoff=ig+ntt                                                 1d21s21
          call enough('mkdoubuc.  1',bc,ibc)
          i18=1                                                         1d21s21
          i28=nrow                                                      1d21s21
          i38=il                                                        1d21s21
          i48=ih                                                        1d21s21
          call ddi_get(bc,ibc,ihddiag(nclop,isb,1),i18,i28,i38,i48,     11d15s22
     $         bc(ivec))                                                11d15s22
          call ddi_get(bc,ibc,ihddiag(nclop,isb,2),i18,i28,i38,i48,     11d15s22
     $         bc(ig))                                                  11d15s22
          nhere2=nhere*2                                                1d22s21
          nrow=mvv*nhere2                                               7d7s21
          nrows=mvv*nhere                                               7d7s21
          itmp=ibcoff                                                   1d22s21
          itmp1=itmp+nrow*nroot                                         1d22s21
          itmp2=itmp1+nrow*nroot                                        1d23s21
          ibcoff=itmp2+nrow*nroot                                       1d23s21
          istmp1=ibcoff                                                 2d18s21
          istmp2=istmp1+nrows*nroot                                     2d18s21
          istmp3=istmp2+nrows*nroot                                     2d18s21
          ibcoff=istmp3+nrows*nroot                                     2d18s21
          call enough('mkdoubuc.  2',bc,ibc)
          nn=mvv*nhere                                                  7d7s21
          if(ngot.gt.0)then                                             1d22s21
           koff=joff                                                    1d22s21
           do i=0,nhere-1                                               1d22s21
            do iv=0,mvv-1                                               7d7s21
             do ir=0,ngot-1                                             1d22s21
              iad=itmp1+iv+mvv*(i+nhere2*ir)                            7d7s21
              jad=iad+nn                                                1d22s21
              bc(iad)=vdold(koff+ir)                                    1d22s21
              bc(jad)=gdold(koff+ir)                                    1d22s21
             end do                                                     1d22s21
             koff=koff+ngot                                             1d22s21
            end do                                                      1d22s21
           end do                                                       1d22s21
           np=nrootu+1                                                  1d22s21
           call dgemm('n','n',nrow,nroot,ngot,1d0,bc(itmp1),nrow,       1d22s21
     $          trans(np,1),nvdim,0d0,bc(itmp),nrow,                    1d22s21
     d' mkdoubuc.  1')
           call dgemm('n','n',nrow,nroot,ngot,1d0,bc(itmp1),nrow,       1d22s21
     $          trans2(np,1),nvdim,0d0,bc(itmp2),nrow,                  1d23s21
     d' mkdoubuc.  2')
           fact=1d0                                                     1d22s21
          else                                                          1d22s21
           fact=0d0                                                     1d22s21
          end if                                                        1d22s21
          jvec=ivec                                                     1d22s21
          jg=ig                                                         1d22s21
          koff=joff                                                     2d18s21
          do i=0,nhere-1                                                1d22s21
           do iv=0,mvv-1                                                7d7s21
            do ir=0,nrootu-1                                            1d22s21
             iad=itmp1+iv+mvv*(i+nhere2*ir)                             7d7s21
             jad=iad+nn                                                 1d22s21
             bc(iad)=bc(jvec+ir)                                        1d22s21
             bc(jad)=bc(jg+ir)                                          1d22s21
            end do                                                      1d22s21
            jvec=jvec+nrootu                                            1d22s21
            jg=jg+nrootu                                                1d22s21
           end do                                                       1d22s21
          end do                                                        1d22s21
          call dgemm('n','n',nrow,nroot,nrootu,1d0,bc(itmp1),nrow,      1d22s21
     $         trans,nvdim,fact,bc(itmp),nrow,                          1d22s21
     d' mkdoubuc.  3')
          call dgemm('n','n',nrow,nroot,nrootu,1d0,bc(itmp1),nrow,      1d22s21
     $         trans2,nvdim,fact,bc(itmp2),nrow,                          1d22s21
     d' mkdoubuc.  4')
          koff=joff                                                     1d22s21
          jvec=ivec                                                     1d23s21
          jg=ig                                                         1d23s21
          do i=0,nhere-1                                                1d22s21
           do iv=0,mvv-1                                                7d7s21
            do ir=0,nroot-1                                             1d22s21
             irp=ir+1                                                   2d18s21
             iad=itmp+iv+mvv*(i+nhere2*ir)                              7d7s21
             jad=iad+nn                                                 1d22s21
             iad2=itmp2+iv+mvv*(i+nhere2*ir)                            7d7s21
             jad2=iad2+nn                                                 1d22s21
             vdold(koff+ir)=bc(iad)                                     1d22s21
             gdold(koff+ir)=bc(jad)                                     1d22s21
             bc(jvec+ir)=bc(iad2)                                       1d23s21
             bc(jg+ir)=bc(jad2)                                         1d23s21
            end do                                                      1d22s21
            koff=koff+nroot                                             1d22s21
            jvec=jvec+nroot                                             1d23s21
            jg=jg+nroot                                                 1d23s21
           end do                                                       1d22s21
          end do                                                        1d22s21
          joff=koff                                                     1d22s21
          i18=1                                                         1d21s21
          i28=nrowb                                                      1d21s21
          i38=il                                                        1d21s21
          i48=ih                                                        1d21s21
          call ddi_put(bc,ibc,ihddiag(nclop,isb,1),i18,i28,i38,i48,     11d15s22
     $         bc(ivec))                                                11d15s22
          call ddi_put(bc,ibc,ihddiag(nclop,isb,2),i18,i28,i38,i48,     11d15s22
     $         bc(ig))                                                  11d15s22
          ibcoff=ivec                                                   1d22s21
         end if                                                         1d22s21
        end if                                                          1d22s21
       end do                                                           1d22s21
      end do                                                            1d22s21
      return                                                            1d22s21
      end                                                               1d22s21
