c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine tofro(veco,vecnono,nroot,nfdat,nvirt,nsymb,multh,      12d21s20
     $     isymmrci,iflag,ndoub,mdoub,bc,ibc)                           11d10s22
      implicit real*8 (a-h,o-z)                                         12d21s20
c
c     iflag =  1: transform veco to vecnono
c     iflag ne 1: transform vecnono to veco
c
      integer*8 i18,i28,i38,i48,i58
      dimension veco(*),vecnono(*),nfdat(5,4,*),nvirt(*),multh(8,8)     12d21s20
      include "common.store"
      equivalence (i58,x58)
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      i18=1                                                             12d21s20
      iflaga=iabs(iflag)                                                1d12s21
      if(iflaga.eq.1)then                                               1d12s21
       do i=1,mdoub*nroot                                               12d21s20
        vecnono(i)=0d0                                                  12d21s20
       end do                                                           12d21s20
       i28=mdoub*nroot                                                  12d21s20
      else                                                              12d21s20
       do i=1,ndoub*nroot                                               12d21s20
        veco(i)=0d0                                                     12d21s20
       end do                                                           12d21s20
       i28=ndoub*nroot                                                  12d21s20
      end if                                                            12d21s20
      ioffo=0                                                           12d21s20
      ioffnono=0                                                        12d21s20
      do isb=1,nsymb                                                    12d21s20
       isbv12=multh(isb,isymmrci)                                       12d21s20
       do isbv1=1,nsymb                                                 12d21s20
        isbv2=multh(isbv1,isbv12)                                       12d21s20
        if(isbv2.ge.isbv1)then                                          12d31s20
         if(isbv1.eq.isbv2)then                                         12d21s20
          nnv=nvirt(isbv1)*nroot                                        12d31s20
          call ilimts(1,nnv,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)   12d21s20
          nhere=ih+1-il                                                 12d21s20
          if(min(nhere,nfdat(3,1,isb)).gt.0)then                        11d17s21
           istarto=ioffo+il                                             12d21s20
           istartn=ioffnono+il                                          12d21s20
           if(iflaga.eq.1)then                                          1d12s21
            itrans=ibcoff                                               12d21s20
            ibcoff=itrans+nfdat(2,1,isb)*nfdat(3,1,isb)                 12d21s20
            call enough('tofro.  1',bc,ibc)
            do i=0,nfdat(2,1,isb)-1
             do j=0,nfdat(3,1,isb)-1
              ji=itrans+j+nfdat(3,1,isb)*i                              12d21s20
              ij=nfdat(4,1,isb)+i+nfdat(2,1,isb)*j                      12d21s20
              bc(ji)=bc(ij)                                             12d21s20
             end do                                                     12d21s20
            end do                                                      12d21s20
            ibcoff=itrans                                               12d21s20
            call dgemm('n','n',nhere,nfdat(2,1,isb),nfdat(3,1,isb),1d0,  12d21s20
     $          veco(istarto),nnv,bc(itrans),nfdat(3,1,isb),0d0,        12d21s20
     $          vecnono(istartn),nnv,                                   12d21s20
     d' tofro.  1')
           else                                                         12d21s20
            call dgemm('n','n',nhere,nfdat(3,1,isb),nfdat(2,1,isb),1d0,  12d21s20
     $       vecnono(istartn),nnv,bc(nfdat(4,1,isb)),nfdat(2,1,isb),0d0,12d21s20
     $          veco(istarto),nnv,                                      12d21s20
     d' tofro.  2')
           end if                                                       12d21s20
          end if                                                        12d21s20
          ioffo=ioffo+nnv*nfdat(3,1,isb)                                12d21s20
          ioffnono=ioffnono+nnv*nfdat(2,1,isb)                          12d21s20
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         12d31s20
         else                                                           12d21s20
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 12d21s20
         end if                                                         12d21s20
         nvv=nvv*nroot                                                  1d6s21
         call ilimts(1,nvv,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)    12d21s20
         nhere=ih+1-il                                                  12d21s20
         do l=1,4                                                       12d21s20
          if(nfdat(3,l,isb).gt.0)then                                   12d21s20
           if(nhere.gt.0)then                                           12d21s20
            istarto=ioffo+il                                             12d21s20
            istartn=ioffnono+il                                          12d21s20
            if(iflaga.eq.1)then                                         1d12s21
             itrans=ibcoff                                               12d21s20
             ibcoff=itrans+nfdat(2,l,isb)*nfdat(3,l,isb)                 12d21s20
             call enough('tofro.  2',bc,ibc)
             do i=0,nfdat(2,l,isb)-1
              do j=0,nfdat(3,l,isb)-1
               ji=itrans+j+nfdat(3,l,isb)*i                              12d21s20
               ij=nfdat(4,l,isb)+i+nfdat(2,l,isb)*j                      12d21s20
               bc(ji)=bc(ij)                                             12d21s20
              end do                                                     12d21s20
             end do                                                      12d21s20
             call dgemm('n','n',nhere,nfdat(2,l,isb),nfdat(3,l,isb),1d0,  12d21s20
     $          veco(istarto),nvv,bc(itrans),nfdat(3,l,isb),0d0,        12d29s20
     $          vecnono(istartn),nvv,                                   12d29s20
     d' tofro.  3')
             ibcoff=itrans                                              12d21s20
            else                                                         12d21s20
             call dgemm('n','n',nhere,nfdat(3,l,isb),nfdat(2,l,isb),1d0,  12d21s20
     $       vecnono(istartn),nvv,bc(nfdat(4,l,isb)),nfdat(2,l,isb),0d0,12d29s20
     $          veco(istarto),nvv,                                      12d29s20
     d' tofro.  4')
            end if                                                       12d21s20
           end if                                                       12d21s20
          end if                                                        12d21s20
          ioffo=ioffo+nvv*nfdat(3,l,isb)                                3d16s21
          ioffnono=ioffnono+nvv*nfdat(2,l,isb)                          3d16s21
         end do                                                         12d21s20
        end if                                                          12d21s20
       end do                                                           12d21s20
      end do                                                            12d21s20
       call dws_synca
      call ddi_create(bc,ibc,i18,i28,i58)                               11d15s22
      call ddi_zero(bc,ibc,i58)                                         11d15s22
      call dws_synca                                                    12d21s20
      if(iflaga.eq.1)then                                               1d12s21
       call ddi_acc(bc,ibc,i58,i18,i18,i18,i28,vecnono)                 11d15s22
      else                                                              12d21s20
       call ddi_acc(bc,ibc,i58,i18,i18,i18,i28,veco)                    11d15s22
      end if                                                            12d21s20
      call dws_synca                                                    12d21s20
      if(iflaga.eq.1)then                                               1d12s21
       call ddi_get(bc,ibc,i58,i18,i18,i18,i28,vecnono)                 11d15s22
      else                                                              12d21s20
       call ddi_get(bc,ibc,i58,i18,i18,i18,i28,veco)                    11d15s22
      end if                                                            12d21s20
      call dws_synca                                                    12d21s20
      call ddi_destroy(i58)                                             12d21s20
      return                                                            12d21s20
      end                                                               12d21s20
