c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine chc4v(nsymb,multh,isymmrci,vd,nrootu,nfdat,nvirt,sr2,  11d16s23
     $     bc,ibc,ndoub,dot,i4xb,iprint)                                       11d16s23
      implicit real*8 (a-h,o-z)                                         11d16s23
      dimension multh(8,8),vd(*),nfdat(5,4,*),nvirt(*),i4xb(*)          11d16s23
      include "common.store"                                            12d12s20
      common/fnd2cm/inv(2,8,8,8)                                        4d9s18
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data loopx/1000/
      loop=0
      dot=0d0                                                           11d25s23
c
c     vd is stored for isb ( for isbv1 ( for l vs,root,nfdat))
c     if uc, vd is stored [nroot,(ncsf2(1)*visv+ncsf*vnotv)]if2
      if(iprint.ne.0)then
       igoal=2
       igoal2=2
      else
       igoal=2
       igoal2=2
      end if
      ibcoffo=ibcoff                                                    11d16s23
      xnan=-2d0
      ivtmp=ibcoff                                                      11d16s23
      ibcoff=ivtmp+ndoub*nrootu                                         11d16s23
      call enough('chc4v.vtmp',bc,ibc)                                  11d16s23
      do iz=ivtmp,ibcoff-1
       bc(iz)=xnan
      end do
      ivoff=1                                                           11d16s23
      jvtmp=ivtmp                                                       11d16s23
      sz1=0d0
      sz2=0d0
      do isb=1,nsymb                                                    11d16s23
       isbv12=multh(isb,isymmrci)                                       11d16s23
       nn=nrootu*nfdat(3,1,isb)                                         11d16s23
       nnn=nn+nrootu*(nfdat(3,2,isb)+nfdat(3,3,isb)+nfdat(3,4,isb))     11d16s23
       do isbv1=1,nsymb                                                 11d16s23
        isbv2=multh(isbv1,isbv12)                                       11d16s23
        if(isbv2.ge.isbv1)then                                          11d16s23
         if(isbv1.eq.isbv2)then                                         11d16s23
          do i=0,nfdat(3,1,isb)-1                                       11d16s23
           do ir=0,nrootu-1                                              11d16s23
            iad=jvtmp+ir+nrootu*i                                       11d16s23
            do iv=0,nvirt(isbv1)-1                                      11d16s23
             bc(iad)=vd(ivoff+iv)                                       11d16s23
             sz1=sz1+vd(ivoff+iv)**2
             iad=iad+nn                                                 11d16s23
            end do                                                      11d16s23
            ivoff=ivoff+nvirt(isbv1)                                    11d16s23
           end do                                                       11d16s23
          end do                                                        11d16s23
          jvtmp=jvtmp+nn*nvirt(isbv1)                                   11d16s23
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         11d16s23
         else                                                           11d16s23
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 11d16s23
         end if                                                         11d16s23
         ii=0
         do l=1,4                                                       11d16s23
          do i=0,nfdat(3,l,isb)-1                                       11d16s23
           ip=i+ii                                                      11d16s23
           do ir=0,nrootu-1                                             11d16s23
            iad=jvtmp+ir+nrootu*ip                                      11d16s23
            do ivv=0,nvv-1                                              11d16s23
             bc(iad)=vd(ivoff+ivv)                                      11d16s23
             sz1=sz1+vd(ivoff+ivv)**2
             iad=iad+nnn                                                11d16s23
            end do
            ivoff=ivoff+nvv                                             11d16s23
           end do
          end do                                                        11d16s23
          ii=ii+nfdat(3,l,isb)                                          11d16s23
         end do                                                         11d16s23
         jvtmp=jvtmp+nnn*nvv
        end if                                                          11d16s23
       end do                                                           11d16s23
      end do                                                            11d16s23
      ivoff=ivtmp                                                       11d16s23
      nrows=0
      jjvv=1                                                            11d16s23
      do isb=1,nsymb                                                    11d16s23
       isbv12=multh(isb,isymmrci)                                       11d16s23
       jsbv12=isbv12                                                    11d16s23
       jjvd=ivoff                                                       11d16s23
       nrow=0                                                           11d16s23
       do isbv1=1,nsymb                                                 11d16s23
        isbv2=multh(isbv1,isbv12)                                       11d16s23
        if(isbv2.ge.isbv1)then                                          11d16s23
         if(isbv2.eq.isbv1)then                                         11d16s23
          nrow=nrow+nrootu*nfdat(3,1,isb)*nvirt(isbv1)                  11d16s23
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         11d16s23
         else                                                           11d16s23
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 11d16s23
         end if                                                         11d16s23
         nvv=nvv*nrootu                                                 11d16s23
         do l=1,4                                                       11d16s23
          nrow=nrow+nfdat(3,l,isb)*nvv                                  11d16s23
         end do                                                         11d16s23
        end if                                                          11d16s23
       end do                                                           11d16s23
       nrows=nrows+nrow
       igg=ibcoff                                                       11d16s23
       ivdtmp=igg+nrow                                                  11d16s23
       itmp=ivdtmp+nrow                                                 11d16s23
       ibcoff=itmp+nrow                                                 11d16s23
       call enough('chc4v.vdtmp',bc,ibc)                                11d16s23
       do iz=igg,itmp-1                                                 11d16s23
        bc(iz)=0d0
       end do
       nna=nrootu*(nfdat(3,1,isb)+nfdat(3,2,isb)+nfdat(3,3,isb)         11d16s23
     $      +nfdat(3,4,isb))                                            11d16s23
       nnam=nna-1                                                       11d16s23
       nn=nrootu*nfdat(3,1,isb)                                         11d16s23
       nnm=nn-1                                                         11d16s23
       do jsbv1=1,nsymb                                                 11d16s23
        jsbv2=multh(jsbv1,jsbv12)                                       11d16s23
        if(jsbv2.ge.jsbv1)then                                          11d16s23
         if(jsbv1.eq.jsbv2)then                                         11d16s23
          jjvdp=jjvd+nrootu*nfdat(3,1,isb)*nvirt(jsbv1)                 11d16s23
          nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                         11d16s23
          isw=0                                                         11d16s23
         else                                                           11d16s23
          jjvdp=jjvd                                                    11d16s23
          isw=1                                                         11d16s23
          nvv=nvirt(jsbv1)*nvirt(jsbv2)                                 11d16s23
         end if                                                         11d16s23
         ktmp=ivdtmp                                                    11d16s23
         do ksbv1=1,nsymb                                               11d16s23
          ksbv2=multh(ksbv1,jsbv12)                                     11d16s23
          if(ksbv2.ge.ksbv1)then                                        11d16s23
           i2euj=inv(1,ksbv1,jsbv1,ksbv2)                               11d16s23
           icasej=inv(2,ksbv1,jsbv1,ksbv2)                              11d16s23
           if(ksbv1.eq.jsbv1)then                                       11d16s23
            nnrowj=(nvirt(ksbv1)*(nvirt(ksbv1)+1))/2                    11d16s23
           else                                                         11d16s23
            nnrowj=nvirt(ksbv1)*nvirt(jsbv1)                            11d16s23
           end if                                                       11d16s23
           if(ksbv1.eq.jsbv2)then                                       11d16s23
            nnrowk=(nvirt(ksbv1)*(nvirt(ksbv1)+1))/2                    11d16s23
           else                                                         11d16s23
            nnrowk=nvirt(ksbv1)*nvirt(jsbv2)                            11d16s23
           end if                                                       11d16s23
           if(jsbv12.eq.1)then                                          11d16s23
            mvv=(nvirt(ksbv1)*(nvirt(ksbv1)-1))/2                       11d16s23
            ktmpp=ktmp+nrootu*nvirt(ksbv1)*nfdat(3,1,isb)               11d16s23
            do jv=0,nvirt(jsbv1)-1                                      11d16s23
             jjvdl=jjvd+nn*jv                                           11d16s23
             do kv2=0,nvirt(ksbv2)-1                                    11d16s23
              ktmpl=ktmp+nn*kv2                                         11d16s23
              ntop=(kv2+isw*(nvirt(ksbv1)-kv2))-1                       11d16s23
              xinta=get4int(i4xb,ksbv1,jsbv1,ksbv2,jsbv2,               11d16s23
     $                      kv2,jv,kv2,jv,nvirt,idum,bc,ibc)            11d16s23
              if(abs(xinta).gt.1d-12)then                                      11d16s23
               do ii=0,nnm                                              11d16s23
                orig=bc(igoal2)
                bc(ktmpl+ii)=bc(ktmpl+ii)+xinta*bc(jjvdl+ii)            11d16s23
                if(abs(orig-bc(igoal2)).gt.1d-12.and.iprint.ne.0)
     $               write(6,*)('4goal2a '),orig,xinta,bc(jjvdl+ii),
     $               bc(igoal2),ksbv1,jsbv1,ksbv2,jsbv2,kv2,jv,kv2,jv,
     $                ibc(ibcoff)
               end do                                                   11d16s23
              end if                                                    11d16s23
              do kv1=0,ntop                                             11d16s23
               xinta=get4int(i4xb,ksbv1,jsbv1,ksbv2,jsbv2,              11d16s23
     $                      kv1,jv,kv2,jv,nvirt,idum,bc,ibc)*sr2        11d10s22
               if(abs(xinta).gt.1d-12)then                                     11d16s23
                ktri=((kv2*(kv2-1))/2)+kv1                              11d16s23
                krec=kv1+nvirt(ksbv1)*kv2                               11d16s23
                kcol=ktri+isw*(krec-ktri)                               11d16s23
                ktmpl=ktmpp+nna*kcol                                    11d16s23
                do ii=0,nnm                                             11d16s23
                orig=bc(igoal2)
                 bc(ktmpl+ii)=bc(ktmpl+ii)+xinta*bc(jjvdl+ii)           11d16s23
                if(abs(orig-bc(igoal2)).gt.1d-12.and.iprint.ne.0)
     $               write(6,*)('4goal2b '),orig,xinta,bc(jjvdl+ii),
     $               bc(igoal2),ksbv1,jsbv1,ksbv2,jsbv2,kv1,jv,kv2,jv,
     $                ibc(ibcoff)
                end do                                                  11d16s23
               end if                                                   11d16s23
              end do                                                    11d16s23
             end do                                                     11d16s23
            end do                                                      11d16s23
            do kv=0,nvirt(ksbv1)-1                                      11d16s23
             ktmpl=ktmp+nn*kv                                           11d16s23
             do jv2=0,nvirt(jsbv2)-1                                    11d16s23
              ntop=(jv2+isw*(nvirt(jsbv1)-jv2))-1                       11d16s23
              icol=kv+1+nvirt(ksbv2)*jv2                                11d16s23
              do jv1=0,ntop                                             11d16s23
               xinta=get4int(i4xb,ksbv1,jsbv1,ksbv2,jsbv2,              11d16s23
     $                      kv,jv1,kv,jv2,nvirt,idum,bc,ibc)*sr2        11d10s22
               if(abs(xinta).gt.1d-12)then                              11d16s23
                jtri=((jv2*(jv2-1))/2)+jv1                              11d16s23
                jrec=jv1+nvirt(jsbv1)*jv2                               11d16s23
                jcol=jtri+isw*(jrec-jtri)                               11d16s23
                jjvdl=jjvdp+nna*jcol                                    11d16s23
                do ii=0,nnm                                             11d16s23
                orig=bc(igoal2)
                 bc(ktmpl+ii)=bc(ktmpl+ii)+xinta*bc(jjvdl+ii)           11d16s23
                if(abs(orig-bc(igoal2)).gt.1d-12.and.iprint.ne.0)
     $               write(6,*)('4goal2c '),orig,xinta,bc(jjvdl+ii),
     $               bc(igoal2),ksbv1,jsbv1,ksbv2,jsbv2,kv,jv1,kv,jv2,
     $                ibc(ibcoff)
                end do                                                  11d16s23
               end if                                                   11d16s23
              end do                                                    11d16s23
             end do                                                     11d16s23
            end do                                                      11d16s23
           else                                                         11d16s23
            ktmpp=ktmp                                                  11d16s23
            mvv=nvirt(ksbv1)*nvirt(ksbv2)                               11d16s23
           end if                                                       11d16s23
           do kv2=0,nvirt(ksbv2)-1                                      11d16s23
            ktop=(kv2+isw*(nvirt(ksbv1)-kv2))-1                         11d16s23
            ktri=((kv2*(kv2-1))/2)                                      11d16s23
            krec=nvirt(ksbv1)*kv2                                       11d16s23
            krow=ktri+isw*(krec-ktri)                                   11d16s23
            do jv2=0,nvirt(jsbv2)-1                                     11d16s23
             jtop=(jv2+isw*(nvirt(jsbv1)-jv2))-1                        11d16s23
             jtri=((jv2*(jv2-1))/2)                                     11d16s23
             jrec=nvirt(jsbv1)*jv2                                      11d16s23
             jrow=jtri+isw*(jrec-jtri)                                  11d16s23
             do kv1=0,ktop                                              11d16s23
              ktmpl=ktmpp+nna*(kv1+krow)                                11d16s23
              do jv1=0,jtop                                             11d16s23
               xinta=get4int(i4xb,ksbv1,jsbv1,ksbv2,jsbv2,              11d16s23
     $                      kv1,jv1,kv2,jv2,nvirt,idum,bc,ibc)          11d10s22
               if(abs(xinta).gt.1d-12)then                                     11d16s23
                jjvdl=jjvdp+nna*(jv1+jrow)                              11d16s23
                do ii=0,nnam                                            11d16s23
                orig=bc(igoal2)
                 bc(ktmpl+ii)=bc(ktmpl+ii)+xinta*bc(jjvdl+ii)           11d16s23
                if(abs(orig-bc(igoal2)).gt.1d-12.and.iprint.ne.0)
     $               write(6,*)('4goal2d '),orig,xinta,bc(jjvdl+ii),
     $               bc(igoal2),ksbv1,jsbv1,ksbv2,jsbv2,kv1,jv1,kv2,jv2,
     $                ibc(ibcoff),ibc(ibcoff+1),ibc(ibcoff+2)
                end do                                                  11d16s23
               end if                                                   11d16s23
              end do                                                    11d16s23
             end do                                                     11d16s23
            end do                                                      11d16s23
            do jv1=0,nvirt(jsbv1)-1                                     11d16s23
             jbot=(jv1+1)*(1-isw)                                       11d16s23
             do jv2=jbot,nvirt(jsbv2)-1                                 11d16s23
              jtri=((jv2*(jv2-1))/2)+jv1                                11d16s23
              jrec=jv1+nvirt(jsbv1)*jv2                                 11d16s23
              jcol=jtri+isw*(jrec-jtri)                                 11d16s23
              jjvdl=jjvdp+nna*jcol                                      11d16s23
              do kv1=0,ktop
               xinta=get4int(i4xb,ksbv1,jsbv2,ksbv2,jsbv1,              11d16s23
     $                      kv1,jv2,kv2,jv1,nvirt,idum,bc,ibc)          11d10s22
               if(abs(xinta).gt.1d-12)then                                     11d16s23
                ktmpl=ktmpp+nna*(kv1+krow)                              11d16s23
                do ii=0,nnm                                             11d16s23
                orig=bc(igoal2)
                 bc(ktmpl+ii)=bc(ktmpl+ii)+xinta*bc(jjvdl+ii)           11d16s23
                if(abs(orig-bc(igoal2)).gt.1d-12.and.iprint.ne.0)
     $               write(6,*)('4goal2e '),orig,xinta,bc(jjvdl+ii),
     $               bc(igoal2),ksbv1,jsbv1,ksbv2,jsbv2,kv1,jv2,kv2,jv1,
     $                ibc(ibcoff)
                end do                                                  11d16s23
                do ii=nn,nnam                                           11d16s23
                orig=bc(igoal2)
                 bc(ktmpl+ii)=bc(ktmpl+ii)-xinta*bc(jjvdl+ii)           11d16s23
                if(abs(orig-bc(igoal2)).gt.1d-12.and.iprint.ne.0)
     $               write(6,*)('4goal2f '),orig,xinta,bc(jjvdl+ii),
     $               bc(igoal2),ksbv1,jsbv1,ksbv2,jsbv2,kv1,jv2,kv2,jv1,
     $                ibc(ibcoff)
                end do                                                  11d16s23
               end if                                                   11d16s23
              end do                                                    11d16s23
             end do                                                     11d16s23
            end do                                                      11d16s23
           end do                                                       11d16s23
           ktmp=ktmpp+mvv*nna                                           11d16s23
          end if                                                        11d16s23
         end do                                                         11d16s23
         jjvd=jjvdp+nvv*nna                                             11d16s23
        end if                                                          11d16s23
       end do                                                           11d16s23
       jtmp=ivdtmp                                                      11d16s23
       kgg=igg                                                          11d16s23
       do jsbv1=1,nsymb                                                 11d16s23
        jsbv2=multh(jsbv1,jsbv12)                                       11d16s23
        if(jsbv2.ge.jsbv1)then                                          11d16s23
         if(jsbv1.eq.jsbv2)then                                         11d16s23
          jtmpp=jtmp+nvirt(jsbv1)*nn                                    11d16s23
          kggp=kgg+nvirt(jsbv1)*nn                                      11d16s23
          nvv=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                         11d16s23
          do iv=0,nvirt(jsbv1)-1                                        11d16s23
           do ii=0,nnm                                                  11d16s23
            ij=kgg+iv+nvirt(jsbv1)*ii                                   11d16s23
            ji=jtmp+ii+nn*iv                                            11d16s23
            orig=bc(igoal)
            bc(ij)=bc(ij)+bc(ji)                                        11d16s23
            if(abs(orig-bc(igoal)).gt.1d-12.and.iprint.ne.0)
     $           write(6,*)('4goala: '),orig,bc(ji),bc(igoal),ji
           end do                                                       11d16s23
          end do                                                        11d16s23
         else                                                           11d16s23
          jtmpp=jtmp                                                    11d16s23
          kggp=kgg                                                      11d16s23
          nvv=nvirt(jsbv1)*nvirt(jsbv2)                                 11d16s23
         end if                                                         11d16s23
         nvvm=nvv-1                                                     11d16s23
         do ivv=0,nvvm                                                  11d16s23
          do ii=0,nnam                                                  11d16s23
           ij=kggp+ivv+nvv*ii                                           11d16s23
           ji=jtmpp+ii+nna*ivv                                          11d16s23
           orig=bc(igoal)
           bc(ij)=bc(ij)+bc(ji)                                         11d16s23
           if(abs(orig-bc(igoal)).gt.1d-12.and.iprint.ne.0)write(6,*)
     $          ('4goalb: '),orig,bc(ji),bc(igoal),ji
          end do                                                        11d16s23
         end do                                                         11d16s23
         jtmp=jtmpp+nvv*nna                                             11d16s23
         kgg=kggp+nvv*nna                                               11d16s23
        end if                                                          11d16s23
       end do                                                           11d16s23
       doth=0d0                                                         11d16s23
       do i=0,nrow-1                                                    11d16s23
        orig=doth
        doth=doth+bc(igg+i)*vd(jjvv+i)                                  11d16s23
        term=bc(igg+i)*vd(jjvv+i)                                       11d16s23
       end do                                                           11d16s23
       jjvv=jjvv+nrow                                                   11d16s23
       dot=dot+doth                                                     11d16s23
       ibcoff=ivdtmp                                                    11d16s23
       ivoff=jjvd                                                       11d16s23
      end do                                                            11d16s23
      ibcoff=ibcoffo                                                    11d16s23
      return                                                            11d16s23
      end                                                               11d16s23
