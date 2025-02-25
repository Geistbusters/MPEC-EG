c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine wrot1(ihsdiag,nff1,ncsf,mdon,mdoop,isymmrci,nvirt,     2d13s22
     $     nsymb,multh,ovr,nroot,vecno,bc,ibc)                          11d9s22
      implicit real*8 (a-h,o-z)
      integer*8 ihsdiag(mdoop,nsymb),i18,i28,i38,i48
      dimension nff1(mdoop,nsymb),ncsf(*),nvirt(*),multh(8,8),
     $     ovr(nroot,*),vecno(*)                                        2d13s22
      include "common.store"                                            8d19s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      i18=1                                                             8d19s21
      do ii=mdon+1,mdoop                                                8d19s21
       iarg=ii-mdon                                                     8d19s21
       do isb=1,nsymb                                                   8d19s21
        isbv=multh(isb,isymmrci)                                        8d19s21
        if(min(nff1(ii,isb),nvirt(isbv)).gt.0)then                      8d19s21
         nrow=nroot*nvirt(isbv)                                         8d19s21
         ncol=ncsf(iarg)*nff1(ii,isb)                                   8d19s21
         call ilimts(1,ncol,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)   8d21s21
         nhere=ih+1-il                                                  8d21s21
         if(nhere.gt.0)then                                             8d21s21
          itmp=ibcoff                                                    8d19s21
          ibcoff=itmp+nrow*nhere                                        8d21s21
          call enough('wrot1.  1',bc,ibc)
          i28=nrow                                                      8d21s21
          i38=il                                                        8d21s21
          i48=ih                                                        8d21s21
          call ddi_get(bc,ibc,ihsdiag(ii,isb),i18,i28,i38,i48,bc(itmp)) 11d15s22
          itmpm=ibcoff                                                  2d13s22
          ibcoff=itmpm+nrow*nhere                                       2d13s22
          call enough('wrot1.  2',bc,ibc)
          nother=nvirt(isbv)*nhere                                      2d13s22
          call dgemm('n','n',nroot,nother,nroot,1d0,vecno,nroot,        2d13s22
     $         bc(itmp),nroot,0d0,bc(itmpm),nroot,                      2d13s22
     d' wrot1.  1')
          do i=0,nrow*nhere-1                                           2d13s22
           bc(itmp+i)=bc(itmpm+i)                                       2d13s22
          end do                                                        2d13s22
          ibcoff=itmpm                                                  2d13s22
          jtmp=itmp-1                                                    8d19s21
          do icol=1,nhere                                               8d21s21
           do iv=1,nvirt(isbv)                                           8d19s21
            do ir=1,nroot                                                8d19s21
             do jr=1,nroot                                               8d19s21
              ovr(jr,ir)=ovr(jr,ir)+bc(jtmp+jr)*bc(jtmp+ir)              8d19s21
             end do                                                      8d19s21
            end do                                                       8d19s21
            jtmp=jtmp+nroot                                              8d19s21
           end do                                                        8d19s21
          end do                                                         8d19s21
          call ddi_put(bc,ibc,ihsdiag(ii,isb),i18,i28,i38,i48,bc(itmp)) 11d15s22
          ibcoff=itmp                                                    8d19s21
         end if                                                         8d21s21
        end if                                                          8d19s21
       end do                                                           8d19s21
      end do                                                            8d19s21
      return                                                            8d18s21
      end                                                               8d18s21
