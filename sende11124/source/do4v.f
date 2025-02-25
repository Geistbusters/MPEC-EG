c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine do4v(ga,gb,nwds,vk,nvirt,nsymb,multh,isymmrci,nfdat,   8d13s21
     $     nroot,n2e,isymop,i2eop,ixmt,phase2,sr2,idoubo,irefo,nbasdws, 8d13s21
     $     lpr,bc,ibc)                                                  11d14s22
      implicit real*8 (a-h,o-z)                                         8d12s21
      logical lpr                                                       8d13s21
      dimension gb(*),vk(*),nvirt(*),multh(8,8),nfdat(5,4,*),isymop(*), 8d12s21
     $     i2eop(2,3),ixmt(8,*),idoubo(*),irefo(*),nbasdws(*),ga(*)     8d13s21
      include "common.store"                                            8d12s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
c     4x contributions: Gvv'ir+/-=[(vu|v'u')+/-(vu'|v'u)]Vuu'ir+/-
c     /sqrt((1+delta vv')(1+delta uu'))
      sr2p=sr2*phase2                                                   8d12s21
      ioff=1
      jjoff=1                                                            8d12s21
      idoit=0                                                           8d12s21
      igoal=0
      if(lpr)write(6,*)('Hi, my name is do4v'),lpr
      do isb=1,nsymb                                                    8d12s21
       isbv12=multh(isb,isymmrci)                                       8d12s21
       nfall=nfdat(3,1,isb)+nfdat(3,2,isb)+nfdat(3,3,isb)+nfdat(3,4,isb)8d12s21
       do isbv1=1,nsymb                                                 8d12s21
        isbv2=multh(isbv1,isbv12)                                       8d12s21
        if(isbv1.le.isbv2)then                                          8d12s21
         if(isbv1.eq.isbv2)then                                         8d12s21
          sz=0d0                                                        8d12s21
          joff=jjoff                                                     8d12s21
          do jsbv1=1,nsymb                                               8d12s21
           jsbv2=multh(jsbv1,isbv12)                                    8d12s21
           if(jsbv2.ge.jsbv1)then                                       8d12s21
            if(jsbv1.eq.jsbv2)then                                      8d12s21
             ncol=nfdat(3,1,isb)*nroot                                  8d12s21
c     Gvvir+=[(vu)(vu)+(vu)(vu)]Vuuir+/sqrt(2*2)
c           =(vu)(vu)*Vuuir+
             if(mod(idoit,mynprocg).eq.mynowprog)then                   8d12s21
              itmpi=ibcoff                                              8d12s21
              ibcoff=itmpi+nvirt(isbv1)*nvirt(jsbv1)                    8d12s21
              call enough('do4v.  1',bc,ibc)
              do iz=itmpi,ibcoff-1                                      8d12s21
               bc(iz)=0d0                                               8d12s21
              end do                                                    8d12s21
              nhit=0                                                    8d12s21
              do ii2e=1,n2e                                             8d12s21
               if(multh(isbv1,jsbv1).eq.isymop(i2eop(1,ii2e)))then      8d12s21
                ivu=ixmt(jsbv1,i2eop(1,ii2e))+idoubo(jsbv1)+irefo(jsbv1)8d12s21
     $              +nbasdws(jsbv1)*(idoubo(isbv1)+irefo(isbv1))        8d12s21
                do iv=0,nvirt(isbv1)-1                                  8d12s21
                 ivui=ivu+nbasdws(jsbv1)*iv                             8d12s21
                 jtmpi=itmpi+nvirt(jsbv1)*iv                            8d12s21
                 do jv=0,nvirt(jsbv1)-1                                 8d12s21
                  bc(jtmpi+jv)=bc(jtmpi+jv)+phase2*bc(ivui+jv)          8d12s21
     $                 *bc(ivui+jv)                                     8d12s21
                 end do                                                 8d12s21
                end do                                                  8d12s21
                nhit=1                                                  8d12s21
               end if                                                   8d12s21
              end do                                                    8d12s21
              if(min(nhit,nvirt(jsbv1),ncol,nvirt(isbv1)).gt.0)then     8d26s21
               call dgemm('n','n',nvirt(jsbv1),ncol,nvirt(isbv1),1d0,   8d12s21
     $              bc(itmpi),nvirt(jsbv1),vk(ioff),nvirt(isbv1),1d0,   8d12s21
     $              gb(joff),nvirt(jsbv1),                              8d12s21
     d' do4v.  1')
              end if                                                    8d12s21
              ibcoff=itmpi                                              8d12s21
             end if                                                     8d12s21
             idoit=idoit+1                                              8d12s21
             joff=joff+ncol*nvirt(jsbv1)                                8d12s21
             mvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                      8d12s21
             jsw=0                                                      8d12s21
            else                                                        8d12s21
             mvv=nvirt(jsbv1)*nvirt(jsbv2)                              8d12s21
             jsw=1                                                      8d12s21
            end if                                                      8d12s21
c     Gvv'ir+=[(vu)(v'u)+(vu)(v'u)]Vuuir+
c     /sqrt(2)
            ncol=nroot*nfdat(3,1,isb)                                   8d12s21
            do ii2e=1,n2e
             if(multh(jsbv1,isbv1).eq.isymop(i2eop(1,ii2e)).and.        8d12s21
     $           multh(jsbv2,isbv1).eq.isymop(i2eop(2,ii2e)))then       8d12s21
              ix1=ixmt(jsbv1,i2eop(1,ii2e))+idoubo(jsbv1)+irefo(jsbv1)  8d12s21
     $             +nbasdws(jsbv1)*(idoubo(isbv1)+irefo(isbv1))         8d12s21
              ix2=ixmt(jsbv2,i2eop(2,ii2e))+idoubo(jsbv2)+irefo(jsbv2)  8d12s21
     $             +nbasdws(jsbv2)*(idoubo(isbv1)+irefo(isbv1))         8d12s21
              do iv=0,nvirt(isbv1)-1                                     8d12s21
               if(mod(idoit,mynprocg).eq.mynowprog)then                    8d12s21
                iad=ioff+iv                                             8d12s21
                ix1i=ix1+nbasdws(jsbv1)*iv                              8d12s21
                ix2i=ix2+nbasdws(jsbv2)*iv                              8d12s21
                do jv2=0,nvirt(jsbv2)-1                                 8d12s21
                 ff=sr2p*bc(ix2i+jv2)                                   8d12s21
                 itop=(jv2+jsw*(nvirt(jsbv1)-jv2))-1                    8d12s21
                 do jv1=0,itop                                          8d12s21
                  fff=ff*bc(ix1i+jv1)                                   8d12s21
                  irec=jv1+nvirt(jsbv1)*jv2                             8d12s21
                  itri=((jv2*(jv2-1))/2)+jv1                            8d12s21
                  irow=itri+jsw*(irec-itri)                             8d12s21
                  jad=joff+irow                                         8d12s21
                  do i=0,ncol-1                                         8d12s21
                   gb(jad+mvv*i)=gb(jad+mvv*i)                          8d12s21
     $                  +fff*vk(iad+nvirt(isbv1)*i)                     8d12s21
                  end do                                                8d12s21
                 end do                                                 8d12s21
                end do                                                  8d12s21
               end if                                                   8d12s21
               idoit=idoit+1                                            8d12s21
              end do                                                    8d12s21
             end if                                                     8d12s21
            end do                                                      8d12s21
            joff=joff+mvv*nroot*nfall                                   8d12s21
           end if                                                       8d12s21
          end do                                                        8d12s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         8d12s21
          ioff=ioff+nvirt(isbv1)*nroot*nfdat(3,1,isb)                   8d12s21
          isw=0                                                         8d12s21
         else                                                           8d12s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 8d12s21
          isw=1                                                         8d12s21
         end if                                                         8d12s21
         joff=jjoff                                                     8d12s21
         do jsbv1=1,nsymb                                               8d12s21
          jsbv2=multh(jsbv1,isbv12)                                     8d12s21
          if(jsbv2.ge.jsbv1)then                                        8d12s21
           if(jsbv1.eq.jsbv2)then                                       8d12s21
c      Gvvir+=[(vu)(vu')+(vu')(vu)]Vuu'ir+
c     /sqrt(2)
            ncol=nroot*nfdat(3,1,isb)                                   8d12s21
            do ii2e=1,n2e                                               8d12s21
             if(multh(jsbv1,isbv1).eq.isymop(i2eop(1,ii2e)).and.
     $           multh(jsbv1,isbv2).eq.isymop(i2eop(2,ii2e)))then       8d12s21
              ix1=ixmt(jsbv1,i2eop(1,ii2e))+idoubo(jsbv1)+irefo(jsbv1)  8d12s21
     $             +nbasdws(jsbv1)*(idoubo(isbv1)+irefo(isbv1))         8d12s21
              ix2=ixmt(jsbv1,i2eop(2,ii2e))+idoubo(jsbv1)+irefo(jsbv1)  8d12s21
     $             +nbasdws(jsbv1)*(idoubo(isbv2)+irefo(isbv2))         8d12s21
              do iv2=0,nvirt(isbv2)-1                                    8d12s21
               if(mod(idoit,mynprocg).eq.mynowprog)then                    8d12s21
                itop=(iv2+isw*(nvirt(isbv1)-iv2))-1                       8d12s21
                ix2i=ix2+nbasdws(jsbv1)*iv2                             8d12s21
                do iv1=0,itop                                             8d12s21
                 ix1i=ix1+nbasdws(jsbv1)*iv1                            8d12s21
                 irec=iv1+nvirt(isbv1)*iv2                                8d12s21
                 itri=((iv2*(iv2-1))/2)+iv1                               8d12s21
                 irow=itri+isw*(irec-itri)                                8d12s21
                 iad=ioff+irow                                            8d12s21
                 do jv=0,nvirt(jsbv1)-1                                 8d12s21
                  ff=sr2p*bc(ix1i+jv)*bc(ix2i+jv)                       8d12s21
                  jad=joff+jv                                           8d12s21
                  do i=0,ncol-1                                         8d12s21
                   gb(jad+i*nvirt(jsbv1))=gb(jad+i*nvirt(jsbv1))        8d12s21
     $                  +ff*vk(iad+nvv*i)                               8d12s21
                  end do                                                8d12s21
                 end do                                                 8d12s21
                end do                                                    8d12s21
               end if                                                   8d12s21
               idoit=idoit+1                                            8d12s21
              end do                                                     8d12s21
             end if                                                     8d12s21
            end do                                                      8d12s21
            joff=joff+nvirt(jsbv1)*nroot*nfdat(3,1,isb)                 8d12s21
            mvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                       8d12s21
            jsw=0                                                       8d12s21
           else                                                         8d12s21
            mvv=nvirt(jsbv1)*nvirt(jsbv2)                               8d12s21
            jsw=1                                                       8d12s21
           end if                                                       8d12s21
c     Gvv'ir+/-=[(vu)(v'u')+/-(vu')(v'u)]Vuu'ir+/-
           tf=1d0                                                       8d12s21
           iioff=ioff                                                   8d12s21
           do l=1,4                                                     8d12s21
            if(nfdat(3,l,isb).gt.0)then                                 8d12s21
             ncol=nfdat(3,l,isb)*nroot                                  8d12s21
             do ii2e=1,n2e                                              8d12s21
              if(multh(jsbv1,isbv1).eq.isymop(i2eop(1,ii2e)).and.       8d12s21
     $            multh(jsbv2,isbv2).eq.isymop(i2eop(2,ii2e)))then      8d12s21
               ix1=ixmt(jsbv1,i2eop(1,ii2e))+idoubo(jsbv1)+irefo(jsbv1) 8d12s21
     $              +nbasdws(jsbv1)*(idoubo(isbv1)+irefo(isbv1))        8d12s21
               ix2=ixmt(jsbv2,i2eop(2,ii2e))+idoubo(jsbv2)+irefo(jsbv2) 8d12s21
     $              +nbasdws(jsbv2)*(idoubo(isbv2)+irefo(isbv2))        8d12s21
               do iv2=0,nvirt(isbv2)-1                                  8d12s21
                ix2i=ix2+nbasdws(jsbv2)*iv2                             8d12s21
                itop=(iv2+isw*(nvirt(isbv1)-iv2))-1                     8d12s21
                do iv1=0,itop                                           8d12s21
                 if(mod(idoit,mynprocg).eq.mynowprog)then                   8d12s21
                  ix1i=ix1+nbasdws(jsbv1)*iv1                           8d12s21
                  irec=iv1+nvirt(isbv1)*iv2                             8d12s21
                  itri=((iv2*(iv2-1))/2)+iv1                            8d12s21
                  irow=itri+isw*(irec-itri)                             8d12s21
                  iad=iioff+irow                                        8d12s21
                  do jv2=0,nvirt(jsbv2)-1                               8d12s21
                   jtop=(jv2+jsw*(nvirt(jsbv1)-jv2))-1                  8d12s21
                   ff=phase2*bc(ix2i+jv2)                               8d12s21
                   do jv1=0,jtop                                        8d12s21
                    fff=ff*bc(ix1i+jv1)                                 8d12s21
                    jrec=jv1+nvirt(jsbv1)*jv2                           8d12s21
                    jtri=((jv2*(jv2-1))/2)+jv1                          8d12s21
                    jrow=jtri+jsw*(jrec-jtri)                           8d12s21
                    jad=joff+jrow                                       8d12s21
                    do i=0,ncol-1                                       8d12s21
                     gb(jad+mvv*i)=gb(jad+mvv*i)+fff*vk(iad+nvv*i)      8d12s21
                    end do                                              8d12s21
                   end do                                               8d12s21
                  end do                                                8d12s21
                 end if                                                 8d12s21
                 idoit=idoit+1                                          8d12s21
                end do                                                  8d12s21
               end do                                                   8d12s21
              end if                                                    8d12s21
             end do                                                     8d12s21
             do ii2e=1,n2e                                              8d12s21
              if(multh(jsbv1,isbv2).eq.isymop(i2eop(1,ii2e)).and.       8d12s21
     $            multh(jsbv2,isbv1).eq.isymop(i2eop(2,ii2e)))then      8d12s21
               ix1=ixmt(jsbv1,i2eop(1,ii2e))+idoubo(jsbv1)+irefo(jsbv1) 8d12s21
     $              +nbasdws(jsbv1)*(idoubo(isbv2)+irefo(isbv2))        8d12s21
               ix2=ixmt(jsbv2,i2eop(2,ii2e))+idoubo(jsbv2)+irefo(jsbv2) 8d12s21
     $              +nbasdws(jsbv2)*(idoubo(isbv1)+irefo(isbv1))        8d12s21
               do iv2=0,nvirt(isbv2)-1                                  8d12s21
                ix1i=ix1+nbasdws(jsbv1)*iv2                             8d12s21
                itop=(iv2+isw*(nvirt(isbv1)-iv2))-1                     8d12s21
                do iv1=0,itop                                           8d12s21
                 if(mod(idoit,mynprocg).eq.mynowprog)then                   8d12s21
                  ix2i=ix2+nbasdws(jsbv2)*iv1                           8d12s21
                  irec=iv1+nvirt(isbv1)*iv2                             8d12s21
                  itri=((iv2*(iv2-1))/2)+iv1                            8d12s21
                  irow=itri+isw*(irec-itri)                             8d12s21
                  iad=iioff+irow                                        8d12s21
                  do jv2=0,nvirt(jsbv2)-1                               8d12s21
                   jtop=(jv2+jsw*(nvirt(jsbv1)-jv2))-1                  8d12s21
                   ff=tf*phase2*bc(ix2i+jv2)                               8d12s21
                   do jv1=0,jtop                                        8d12s21
                    fff=ff*bc(ix1i+jv1)                                 8d12s21
                    jrec=jv1+nvirt(jsbv1)*jv2                           8d12s21
                    jtri=((jv2*(jv2-1))/2)+jv1                          8d12s21
                    jrow=jtri+jsw*(jrec-jtri)                           8d12s21
                    jad=joff+jrow                                       8d12s21
                    do i=0,ncol-1                                       8d12s21
                     gb(jad+mvv*i)=gb(jad+mvv*i)+fff*vk(iad+nvv*i)      8d12s21
                    end do                                              8d12s21
                   end do                                               8d12s21
                  end do                                                8d12s21
                 end if                                                 8d12s21
                 idoit=idoit+1                                          8d12s21
                end do                                                  8d12s21
               end do                                                   8d12s21
              end if                                                    8d12s21
             end do                                                     8d12s21
             joff=joff+mvv*ncol                                         8d12s21
             iioff=iioff+nvv*ncol                                       8d12s21
            end if                                                      8d12s21
            tf=-1d0                                                     8d12s21
           end do                                                       8d12s21
          end if                                                        8d12s21
         end do                                                         8d12s21
         ioff=ioff+nvv*nroot*nfall                                      8d12s21
        end if                                                          8d12s21
       end do                                                           8d12s21
       jjoff=joff                                                       8d12s21
      end do                                                            8d12s21
      call dws_gsumf(gb,nwds)                                           8d13s21
      do i=1,nwds                                                       8d13s21
       ga(i)=ga(i)+gb(i)                                                8d13s21
      end do                                                            8d13s21
      if(lpr)then
       write(6,*)('gb in orthogonal basis with 4v: ')
       ioff=0                                                            8d12s21
       do isb=1,nsymb                                                    8d12s21
        isbv12=multh(isb,isymmrci)                                       8d12s21
        do isbv1=1,nsymb                                                 8d12s21
         isbv2=multh(isbv1,isbv12)                                       8d12s21
         if(isbv2.ge.isbv1)then                                          8d12s21
          if(isbv1.eq.isbv2)then                                         8d12s21
           sz=0d0                                                        8d12s21
           nrow=nvirt(isbv1)*nroot                                       8d12s21
           do i=1,nrow*nfdat(3,1,isb)                                    8d12s21
            sz=sz+ga(ioff+i)**2                                          8d12s21
           end do                                                        8d12s21
           sz=sqrt(sz/dfloat(nrow*nfdat(3,1,isb)))                       8d12s21
           if(sz.gt.1d-10)then                                           8d12s21
            write(6,*)('visv for syms '),isb,isbv1,ioff
            call prntm2(ga(ioff+1),nrow,nfdat(3,1,isb),nrow)             8d12s21
           end if                                                        8d12s21
           ioff=ioff+nrow*nfdat(3,1,isb)                                 8d12s21
           nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         8d12s21
          else                                                           8d12s21
           nvv=nvirt(isbv1)*nvirt(isbv2)                                 8d12s21
          end if                                                         8d12s21
          nrow=nvv*nroot                                                 8d12s21
          do l=1,4                                                       8d12s21
           if(nfdat(3,l,isb).gt.0)then                                   8d12s21
            sz=0d0                                                       8d12s21
            do i=1,nrow*nfdat(3,l,isb)                                   8d12s21
             sz=sz+ga(ioff+i)**2                                         8d12s21
            end do                                                       8d12s21
            sz=sqrt(sz/dfloat(nrow*nfdat(3,l,isb)))                      8d12s21
            if(sz.gt.1d-10)then                                          8d12s21
             write(6,*)('vnotv for spins '),l,('symmetries '),isb,isbv1,
     $          isbv2,ioff
             call prntm2(ga(ioff+1),nrow,nfdat(3,l,isb),nrow)            8d12s21
            end if                                                       8d12s21
           end if                                                        8d12s21
           ioff=ioff+nrow*nfdat(3,l,isb)                                 8d12s21
          end do                                                         8d12s21
         end if                                                          8d12s21
        end do                                                           8d12s21
       end do                                                            8d12s21
       write(6,*)('ioff at end: '),ioff
      end if                                                            8d13s21
      return                                                            8d12s21
      end                                                               8d12s21
